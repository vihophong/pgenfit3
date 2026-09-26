// Simulate a MIXTURE of implanted ion species: each species has its own
// decay network (parmsex) and beam/detector settings (simparms: implant
// rate, implant profile, backgrounds, efficiencies). Each species is run by
// its own simulation object over the same beam time; all implant, decay and
// neutron hits are then merged into one stream and correlated together, so
// decays of one species are also (accidentally) correlated with implants of
// the others, as in data where the species are not separated.
//
// Correlation settings (window, position cut) come from the FIRST simparms.
// Backgrounds of all simparms add up -- normally give them in the first one
// and set betabkgrateg/betabkgrateu/neubkgrate/r2neubkgrate to 0 in the
// others. The beam time and tsoffset should be the same in all.
//
// Output as simulation.cc (trees tree/treebw/..., histograms, window
// parameters), plus: branch "ionspecies" in "tree" (true species of the
// correlated implant) and TParameter<double> nspecies, nimplant_s<i>.
//
// Usage: ./simulation_mix <out.root> <seed> <parmsex_1> <simparms_1> [<parmsex_2> <simparms_2> ...]

#include <iostream>
#include <vector>
#include "simulation.hh"
#include "TFile.h"
#include "TParameter.h"

int main(int argc, char *argv[])
{
    if (argc<5 || (argc-3)%2!=0){
        std::cout<<"Usage: "<<argv[0]<<" <out.root> <seed> <parmsex_1> <simparms_1> [<parmsex_2> <simparms_2> ...]"<<std::endl;
        return 1;
    }
    const int seed=atoi(argv[2]);
    const int nsp=(argc-3)/2;
    std::vector<simulation*> sims;
    for (int i=0;i<nsp;i++){
        std::cout<<"==== species "<<i<<": "<<argv[3+2*i]<<"  "<<argv[4+2*i]<<std::endl;
        simulation* sim=new simulation(argv[3+2*i]);
        sim->readSimulationParameters(argv[4+2*i]);
        // independent random streams per species (seed 0 = random for all)
        sim->setRandomSeed(seed==0 ? 0 : seed+1000003*i);
        sim->setSpeciesIndex(i);
        sims.push_back(sim);
    }
    TFile* fout=new TFile(argv[1],"recreate");
    fout->cd();
    simulation* s0=sims[0];
    s0->BookSimulationTree();
    s0->BookCorrelationTree();
    std::vector<int> nimp(nsp);
    for (int i=0;i<nsp;i++){
        sims[i]->runSimulation();
        nimp[i]=sims[i]->getNImplants();
        std::cout<<"species "<<i<<": "<<nimp[i]<<" implants"<<std::endl;
        if (i>0) s0->mergeHits(sims[i]); // free species i's hits right away
    }
    s0->fillTreeData();
    s0->correlateData();
    fout->cd();
    s0->getIonSimulationTree()->Write();
    s0->getBetaSimulationTree()->Write();
    s0->getNeutronSimulationTree()->Write();
    s0->getCorrelationTree()->Write();
    s0->getMLHTree()->Write();
    s0->getMLHTreeBackward()->Write();
    s0->writeMLHHistos();
    TParameter<double>("nspecies",nsp).Write();
    for (int i=0;i<nsp;i++) TParameter<double>(Form("nimplant_s%d",i),nimp[i]).Write();
    fout->Close();
    return 0;
}
