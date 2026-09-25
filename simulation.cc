//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * Copyright@2019 Vi Ho Phong, email: phong@ribf.riken.jp           *
// *                                                                  *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications.                    *
// ********************************************************************
//
/// \file simulation.cc
/// \brief Implementation of the simulation class
/// Only support 1 isomeric state per isotope for now.

#include "simulation.hh"

#include <iostream>
#include <fstream>


#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TBox.h"

#include "TLine.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TPad.h"
#include "TFrame.h"
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

simulation::simulation(char* inputParm):fneutronData(),fbetaData(),fionData(),fcorrNeutronData_fw(),fcorrNeutronData_bw()
{
    Init(inputParm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

simulation::simulation(decaypath* path):fneutronData(),fbetaData(),fionData(),fcorrNeutronData_fw(),fcorrNeutronData_bw()
{
    Init(path);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

simulation::~simulation()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::Init(char* inputParms)
{
    fdecaypathobj=new decaypath();
    fdecaypathobj->Init(inputParms);
    fdecaypathobj->makePath();

    Int_t index=0;
    for (Int_t i=0;i<fdecaypathobj->getNMember();i++){
        for (Int_t j=0;j<fdecaypathobj->getMember(i)->path.size();j++){
             fripathToId[i][j]=index;
             index++;
        }
    }
    cout<<"index="<<index<<endl;

    fprimImplantEvt=0;
    fprimImplantT=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::BookSimulationTree()
{
    ftreeion=new TTree("ion","tree ion");
    ftreebeta=new TTree("tbeta","tree beta");
    ftreeneu=new TTree("neutron","tree neutron");
    ftreeion->Branch("ion",&fionData,"T/D:Tcorr/D:x/D:y/D:z/D:evt/I:mode/I:id/I:ifl/I");
    ftreeion->Branch("fl_n",&fionData.fl_n,"fl_n/I");
    ftreeion->Branch("fl_i",&fionData.fl_i,"fl_i[fl_n]/I");
    ftreebeta->Branch("beta",&fbetaData,"T/D:Tcorr/D:x/D:y/D:z/D:evt/I:mode/I:id/I:ifl/I");
    ftreebeta->Branch("fl_n",&fbetaData.fl_n,"fl_n/I");
    ftreebeta->Branch("fl_i",&fbetaData.fl_i,"fl_i[fl_n]/I");
    ftreeneu->Branch("neu",&fneutronData,"T/D:Tcorr/D:x/D:y/D:z/D:evt/I:mode/I:id/I:ifl/I");
    ftreeneu->Branch("fl_n",&fneutronData.fl_n,"fl_n/I");
    ftreeneu->Branch("fl_i",&fneutronData.fl_i,"fl_i[fl_n]/I");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::BookCorrelationTree()
{
    fdeltaxy=2.;
    fionbetawindowlow=fsimparms.ionbetawindowlow;
    fionbetawindowup=fsimparms.ionbetawindowup;
    fwindowbetaneutronlow=400000./1e9;
    fwindowbetaneutronup=400000./1e9;

    ftreemlh=new TTree("tree","tree");
    ftreemlh->Branch("x",&fmlh_t,"x/D");
    ftreemlh->Branch("y",&fmlh_mult,"y/I");
    ftreemlh->Branch("ionT",&fmlh_ionT,"ionT/D");
    ftreemlh->Branch("alpha",&fmlh_alpha,"alpha/I");

    ftreemlhbw=new TTree("treebw","treebw");
    ftreemlhbw->Branch("x",&fmlhbw_t,"x/D");
    ftreemlhbw->Branch("y",&fmlhbw_mult,"y/I");

    ftreecorr=new TTree("treecorr","treecorr");

    ftreecorr->Branch("ion",&fionData,"T/D:Tcorr/D:x/D:y/D:z/D:evt/I:mode/I:id/I:ifl/I:fl_n/I");
    ftreecorr->Branch("beta",&fbetaData,"T/D:Tcorr/D:x/D:y/D:z/D:evt/I:mode/I:id/I:ifl/I:fl_n/I");

    ftreecorr->Branch("multf",&fcorrNeutronData_fw.mult,"multf/I");
    ftreecorr->Branch("Tf",fcorrNeutronData_fw.T,"Tf[multf]/D");
    ftreecorr->Branch("Tcorrf",fcorrNeutronData_fw.Tcorr,"Tcorrf[multf]/D");
    ftreecorr->Branch("xf",fcorrNeutronData_fw.x,"xf[multf]/D");
    ftreecorr->Branch("yf",fcorrNeutronData_fw.y,"yf[multf]/D");
    ftreecorr->Branch("zf",fcorrNeutronData_fw.z,"zf[multf]/D");
    ftreecorr->Branch("evtf",fcorrNeutronData_fw.evt,"evtf[multf]/I");
    ftreecorr->Branch("modef",fcorrNeutronData_fw.mode,"modef[multf]/I");
    ftreecorr->Branch("idf",fcorrNeutronData_fw.id,"idf[multf]/I");
    ftreecorr->Branch("iflf",fcorrNeutronData_fw.ifl,"iflf[multf]/I");
    ftreecorr->Branch("fl_nf",fcorrNeutronData_fw.fl_n,"fl_nf[multf]/I");

    ftreecorr->Branch("multb",&fcorrNeutronData_bw.mult,"multb/I");
    ftreecorr->Branch("Tb",fcorrNeutronData_bw.T,"Tb[multb]/D");
    ftreecorr->Branch("Tcorrb",fcorrNeutronData_bw.Tcorr,"Tcorrb[multb]/D");
    ftreecorr->Branch("xb",fcorrNeutronData_bw.x,"xb[multb]/D");
    ftreecorr->Branch("yb",fcorrNeutronData_bw.y,"yb[multb]/D");
    ftreecorr->Branch("zb",fcorrNeutronData_bw.z,"zb[multb]/D");
    ftreecorr->Branch("evtb",fcorrNeutronData_bw.evt,"evtb[multb]/I");
    ftreecorr->Branch("modeb",fcorrNeutronData_bw.mode,"modeb[multb]/I");
    ftreecorr->Branch("idb",fcorrNeutronData_bw.id,"idb[multb]/I");
    ftreecorr->Branch("iflb",fcorrNeutronData_bw.ifl,"iflb[multb]/I");
    ftreecorr->Branch("fl_nb",fcorrNeutronData_bw.fl_n,"fl_nb[multb]/I");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::bookTDiffData(){
    TFile *fintdiff = new TFile("tdiffdataex.root");
    ftreetdiff=0;
    fintdiff->GetObject("tree",ftreetdiff);
    ftreetdiff->SetBranchAddress("tdiff",&ftdiff);
    ientrytdiff=fsimparms.iTdiffEntryBegin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::readSimulationParameters(char *inputfile)
{
    std::ifstream ifscond(inputfile);
    std::string line;

    std::string line_head;
    Double_t  line_val;

    while (!ifscond.eof()){
        std::getline(ifscond,line);
        std::stringstream ss(line);
        ss>>line_head;
        if (line_head.at(0)=='#') continue;
        ss>>line_val;

        if (line_head=="istdifffromfile") fsimparms.isTdiffFromFile=(Int_t)line_val;
        if (line_head=="ntdiffentries") fsimparms.nTdiffEntries=(Int_t)line_val;
        if (line_head=="itdiffentrybegin") fsimparms.iTdiffEntryBegin=(Int_t)line_val;
        if (line_head=="isfiximppos") fsimparms.isFixImplantPosition=(Int_t)line_val;
        if (line_head=="ishistgausbkg") fsimparms.isSpatialDistFromHist=(Int_t)line_val;

        if (line_head=="beamtime") fsimparms.BeamTime=line_val;
        if (line_head=="beamrate") fsimparms.rate=line_val;
        if (line_head=="isoperctg") fsimparms.isoperctg=line_val;
        if (line_head=="betabkgrateg") fsimparms.betabkgrateg=line_val;
        if (line_head=="betabkgrateu") fsimparms.betabkgrateu=line_val;
        if (line_head=="neubkgrate") fsimparms.neurndbkgrate=line_val;
        if (line_head=="r2neubkgrate") fsimparms.r2neurndbkgrate=line_val;
        if (line_head=="randbetaneuperctg") fsimparms.percentage_bkgnb=line_val;
        if (line_head=="randbeta2neuperctg") fsimparms.percentage_bkg2nb=line_val;
        if (line_head=="betaeff") fsimparms.betaeff=line_val;
        if (line_head=="neueff") fsimparms.neutroneff=line_val;
        if (line_head=="betaneutronmodtime") fsimparms.betaneutronmodtime=line_val;
        if (line_head=="beamneutronmodtime") fsimparms.beamneutronmodtime=line_val;

        if (line_head=="xmin") fsimparms.xmin=line_val;
        if (line_head=="xmax") fsimparms.xmax=line_val;
        if (line_head=="ymin") fsimparms.ymin=line_val;
        if (line_head=="ymax") fsimparms.ymax=line_val;

        if (line_head=="ximpmean") fsimparms.ximpmean=line_val;
        if (line_head=="ximpsigma") fsimparms.ximpsigma=line_val;
        if (line_head=="yimpmean") fsimparms.yimpmean=line_val;
        if (line_head=="yimpsigma") fsimparms.yimpsigma=line_val;

        if (line_head=="deltaxylimit") fsimparms.deltaxy=line_val;
        if (line_head=="dxbetamean") fsimparms.dxbetamean=line_val;
        if (line_head=="dxbetasigma") fsimparms.dxbetasigma=line_val;
        if (line_head=="dybetamean") fsimparms.dybetamean=line_val;
        if (line_head=="dybetasigma") fsimparms.dybetasigma=line_val;


        if (line_head=="xbetabkgmean") fsimparms.xbetabkgmean=line_val;
        if (line_head=="xbetabkgsigma") fsimparms.xbetabkgsigma=line_val;
        if (line_head=="ybetabkgmean") fsimparms.ybetabkgmean=line_val;
        if (line_head=="ybetabkgsigma") fsimparms.ybetabkgsigma=line_val;


        if (line_head=="alphaeff") fsimparms.alphaeff=line_val;
        if (line_head=="deltaxylimitalpha") fsimparms.deltaxyalpha=line_val;
        if (line_head=="dxalphamean") fsimparms.dxalphamean=line_val;
        if (line_head=="dxalphasigma") fsimparms.dxalphasigma=line_val;
        if (line_head=="dyalphamean") fsimparms.dyalphamean=line_val;
        if (line_head=="dyalphasigma") fsimparms.dyalphasigma=line_val;

        if (line_head=="tsoffset") fsimparms.tsoffset=line_val;
        if (line_head=="ionbetawindowlow") fsimparms.ionbetawindowlow=line_val;
        if (line_head=="ionbetawindowup") fsimparms.ionbetawindowup=line_val;
        if (line_head=="neuwbeamperctg") fsimparms.neuwbeamperctg=line_val;
    }

    cout<<"*****************\nSimulation parameters:\n"<<endl;
    cout<<"Fix implant position?  "<<fsimparms.isFixImplantPosition<<endl;
    cout<<"Beam time = "<<fsimparms.BeamTime<<endl;
    cout<<"Beam rate = "<<fsimparms.rate<<endl;
    cout<<"Isotope percentage = "<<fsimparms.isoperctg<<endl;
    cout<<"Beta Gaussian background = "<<fsimparms.betabkgrateg<<endl;
    cout<<"Beta Uniform background = "<<fsimparms.betabkgrateu<<endl;
    cout<<"Neutron background = "<<fsimparms.neurndbkgrate<<endl;
    cout<<"2 Neutron background = "<<fsimparms.r2neurndbkgrate<<endl;
    cout<<"Percentage of neutron correlated with random beta background = "<<fsimparms.percentage_bkgnb<<endl;
    cout<<"Percentage of 2 neutron correlated with random beta background = "<<fsimparms.percentage_bkg2nb<<endl;
    cout<<"Beta Efficiency = "<<fsimparms.betaeff<<endl;
    cout<<"Neutron Efficiency = "<<fsimparms.neutroneff<<endl;
    cout<<"Beta neutron moderation time = "<<fsimparms.betaneutronmodtime<<endl;
    cout<<"Beam neutron moderation time = "<<fsimparms.beamneutronmodtime<<endl;

    cout<<"xmin = "<<fsimparms.xmin<<endl;
    cout<<"xmax = "<<fsimparms.xmax<<endl;
    cout<<"ymin = "<<fsimparms.ymin<<endl;
    cout<<"ymax = "<<fsimparms.ymax<<endl;

    cout<<"ximpmean = "<<fsimparms.ximpmean<<endl;
    cout<<"ximpsigma = "<<fsimparms.ximpsigma<<endl;
    cout<<"yimpmean = "<<fsimparms.yimpmean<<endl;
    cout<<"yimpsigma = "<<fsimparms.yimpsigma<<endl;

    cout<<"deltaxylimit = "<<fsimparms.deltaxy<<endl;
    cout<<"dxbetamean = "<<fsimparms.dxbetamean<<endl;
    cout<<"dxbetasigma = "<<fsimparms.dxbetasigma<<endl;
    cout<<"dybetamean = "<<fsimparms.dybetamean<<endl;
    cout<<"dybetasigma = "<<fsimparms.dybetasigma<<endl;


    cout<<"xbetabkgmean = "<<fsimparms.xbetabkgmean<<endl;
    cout<<"xbetabkgsigma = "<<fsimparms.xbetabkgsigma<<endl;
    cout<<"ybetabkgmean = "<<fsimparms.ybetabkgmean<<endl;
    cout<<"ybetabkgsigma = "<<fsimparms.ybetabkgsigma<<endl;

    cout<<"ybetabkgmean = "<<fsimparms.ybetabkgmean<<endl;
    cout<<"ybetabkgsigma = "<<fsimparms.ybetabkgsigma<<endl;

    cout<<"Alpha efficiency = "<<fsimparms.alphaeff<<endl;
    cout<<"deltaxylimitalpha = "<<fsimparms.deltaxyalpha<<endl;
    cout<<"dxalphamean = "<<fsimparms.dxalphamean<<"  dxalphasigma = "<<fsimparms.dxalphasigma<<endl;
    cout<<"dyalphamean = "<<fsimparms.dyalphamean<<"  dyalphasigma = "<<fsimparms.dyalphasigma<<endl;
    cout<<"Offset time stamp = "<<fsimparms.tsoffset<<endl;
    cout<<"Ion-decay correlation window = -"<<fsimparms.ionbetawindowlow<<" s to +"<<fsimparms.ionbetawindowup<<" s"<<endl;
    cout<<"Percentage of neutron correlated with beam = "<<fsimparms.neuwbeamperctg<<endl;

    cout<<"*****************\n"<<endl;


    hhx=NULL;
    hhy=NULL;
    hhimpx=NULL;
    hhimpy=NULL;
    hhdx=NULL;
    hhdy=NULL;

    if (fsimparms.isSpatialDistFromHist!=0){
        cout<<"reading histograms hy.root, hx.root, himpx.root, himpy.root, hdx.root and hdy.root for beta background distribution"<<endl;
        TFile* f1x = TFile::Open("hx.root");
        hhx=(TH1F*)f1x->Get("hx");
        TFile* f1y = TFile::Open("hy.root");
        hhy=(TH1F*)f1y->Get("hy");
        TFile* f1impx = TFile::Open("himpx.root");
        hhimpx=(TH1F*)f1impx->Get("himpx");
        TFile* f1impy = TFile::Open("himpy.root");
        hhimpy=(TH1F*)f1impy->Get("himpy");
        TFile* f1hdx = TFile::Open("hdx.root");
        hhdx=(TH1F*)f1hdx->Get("hdx");
        TFile* f1hdy = TFile::Open("hdy.root");
        hhdy=(TH1F*)f1hdy->Get("hdy");
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::printPathMembers()
{
    fdecaypathobj->printMember();
    fdecaypathobj->printPath();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::registerDecay(MemberDef* decaymember,Int_t pathid,Double_t decaytime)
{
    Int_t uniqueId;
    if (pathid<0)//parent decay
        uniqueId=0;
    else
        uniqueId=fripathToId[decaymember->id][pathid];

    if (decaymember->sim_neumult==100){
        //! alpha decay: its own position spread and detection efficiency,
        //! recorded with the betas (a DSSD decay event), mode 10, no neutrons
        Double_t dxalpha, dyalpha;
        do {
            dxalpha=rseed->Gaus(fsimparms.dxalphamean,fsimparms.dxalphasigma);
            dyalpha=rseed->Gaus(fsimparms.dyalphamean,fsimparms.dyalphasigma);
        } while (!(abs(dxalpha)<fsimparms.deltaxyalpha && abs(dyalpha)<fsimparms.deltaxyalpha));
        simulationdatatype alphahit;
        alphahit.evt=fprimImplantEvt;
        alphahit.Tcorr=decaytime;
        alphahit.T=fprimImplantT+decaytime;
        alphahit.x=dxalpha+fximp;
        alphahit.y=dyalpha+fyimp;
        alphahit.z=0;
        alphahit.mode=10;
        alphahit.id=decaymember->id;
        alphahit.ifl=uniqueId;
        if (pathid>=0){
            alphahit.fl_n=decaymember->path[pathid].size();
            for (Int_t k=0;k<decaymember->path[pathid].size();k++)
                alphahit.fl_i[k]=decaymember->path[pathid][k];
        }else{
            alphahit.fl_n=0;
        }
        if (rseed->Rndm()*100<fsimparms.alphaeff) {
            betaMap.insert(make_pair(alphahit.T,alphahit));
            fbetaEvt++;
        }
        return;
    }

    //! beta decay
    Double_t dxbeta=rseed->Gaus(fsimparms.dxbetamean,fsimparms.dxbetasigma);
    Double_t dybeta=rseed->Gaus(fsimparms.dybetamean,fsimparms.dybetasigma);
    
    while (!(abs(dxbeta)<fsimparms.deltaxy && abs(dybeta)<fsimparms.deltaxy)){
        dxbeta=rseed->Gaus(fsimparms.dxbetamean,fsimparms.dxbetasigma);
        dybeta=rseed->Gaus(fsimparms.dybetamean,fsimparms.dybetasigma);
    }
    if (fsimparms.isSpatialDistFromHist!=0){
        dxbeta=hhdx->GetRandom();
        dybeta=hhdy->GetRandom();
    }

    Double_t xbeta=dxbeta+fximp;
    Double_t ybeta=dybeta+fyimp;

    simulationdatatype betahit;
    betahit.evt=fprimImplantEvt;
    betahit.Tcorr=decaytime;
    betahit.T=fprimImplantT+decaytime;
    betahit.x=xbeta;
    betahit.y=ybeta;
    betahit.z=0;
    betahit.mode=decaymember->sim_neumult;
    betahit.id=decaymember->id;
    betahit.ifl=uniqueId;
    ftsbeta=betahit.T;
    if (pathid>=0){//parent decay
        betahit.fl_n=decaymember->path[pathid].size();
        for (Int_t k=0;k<decaymember->path[pathid].size();k++)// loop all elelemts of the path
            betahit.fl_i[k]=decaymember->path[pathid][k];
    }else{
        betahit.fl_n=0;
    }

    if (rseed->Rndm()*100<fsimparms.betaeff) {
        betaMap.insert(make_pair(betahit.T,betahit));
        fbetaEvt++;
    }

    //! beta decay with neutron
    if (decaymember->sim_neumult==1){//1 neutron decay
        Double_t pneu=rseed->Rndm();
        Double_t dtneu=rseed->Exp(fsimparms.betaneutronmodtime/TMath::Log(2));
        if (pneu<=decaymember->neueff){
            simulationdatatype neuhit;
            neuhit.Tcorr=dtneu;
            neuhit.T=betahit.T+dtneu;
            neuhit.evt=betahit.evt;
            neuhit.mode=betahit.mode;
            neuhit.id=betahit.id;
            neuhit.ifl=uniqueId;
            neuhit.fl_n=betahit.fl_n;
            neuhit.x = betahit.Tcorr;
            memcpy(neuhit.fl_i,betahit.fl_i,sizeof(betahit.fl_i));
            neuMap.insert(make_pair(neuhit.T,neuhit));
            fneuEvt++;
        }


    }else if (decaymember->sim_neumult==2){//2 neutrons decay
        Double_t pneu=rseed->Rndm();
        Double_t dtneu=rseed->Exp(fsimparms.betaneutronmodtime/TMath::Log(2));

        if (pneu<=decaymember->neueff){
            simulationdatatype neuhit1;
            neuhit1.Tcorr=dtneu;
            neuhit1.T=betahit.T+dtneu;
            neuhit1.evt=betahit.evt;
            neuhit1.mode=betahit.mode;
            neuhit1.id=betahit.id;
            neuhit1.ifl=uniqueId;
            neuhit1.fl_n=betahit.fl_n;
            neuhit1.x = betahit.Tcorr;
            memcpy(neuhit1.fl_i,betahit.fl_i,sizeof(betahit.fl_i));
            neuMap.insert(make_pair(neuhit1.T,neuhit1));
            fneuEvt++;
        }

        pneu=rseed->Rndm();
        dtneu=rseed->Exp(fsimparms.betaneutronmodtime/TMath::Log(2));

        if (pneu<=decaymember->neueff){
            simulationdatatype neuhit2;
            neuhit2.Tcorr=dtneu;
            neuhit2.T=betahit.T+dtneu;
            neuhit2.evt=betahit.evt;
            neuhit2.mode=betahit.mode;
            neuhit2.id=betahit.id;
            neuhit2.ifl=uniqueId;
            neuhit2.fl_n=betahit.fl_n;
            neuhit2.x = betahit.Tcorr;
            memcpy(neuhit2.fl_i,betahit.fl_i,sizeof(betahit.fl_i));
            neuMap.insert(make_pair(neuhit2.T,neuhit2));
            fneuEvt++;
        }
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::registerIonImplant(){
    //! beam outside DSSD
    fximp=rseed->Gaus(fsimparms.ximpmean,fsimparms.ximpsigma);
    fyimp=rseed->Gaus(fsimparms.yimpmean,fsimparms.yimpsigma);

    if (fsimparms.isFixImplantPosition!=0){
        fximp=fsimparms.ximpmean;
        fyimp=fsimparms.yimpmean;
    }

    if (fsimparms.isSpatialDistFromHist!=0) {
        fximp=hhimpx->GetRandom();
        fyimp=hhimpy->GetRandom();
    }

    if (fximp<fsimparms.xmax&&fyimp<fsimparms.ymax&&fximp>=fsimparms.xmin&&fyimp>=fsimparms.ymin){
        Double_t dtimp=0;
        if (fsimparms.isTdiffFromFile==1) // f
            dtimp=ftdiff;
        else
            dtimp=rseed->Exp(1/fsimparms.rate);

        fprimImplantT+=dtimp;
        //! neutron associated with beam
        Double_t p1neui=rseed->Rndm()*100;
        if (p1neui<fsimparms.neuwbeamperctg){//note this neuwbeamperctg includes the efficiency
            simulationdatatype neuhit;
            neuhit.Tcorr=rseed->Exp(fsimparms.beamneutronmodtime/TMath::Log(2));
            neuhit.x=0;
            neuhit.y=0;
            neuhit.z=0;
            neuhit.T=neuhit.Tcorr+fprimImplantT;
            neuhit.mode=0;
            neuhit.id=-1;
            neuhit.evt=-fprimImplantEvt;
            neuhit.fl_n=0;
            neuMap.insert(make_pair(neuhit.T,neuhit));
        }
        Double_t pimp=rseed->Rndm()*100;
        if (pimp<fsimparms.isoperctg){
            fprimImplantEvt++;
            simulationdatatype implanthit;
            implanthit.T=fprimImplantT;
            implanthit.Tcorr=dtimp;
            implanthit.x=fximp;
            implanthit.y=fyimp;
            implanthit.z=0;
            implanthit.id=-1;
            implanthit.mode=-1;
            implanthit.evt=fprimImplantEvt;
            implanthit.fl_n=0;
            ionMap.insert(make_pair(implanthit.T,implanthit));
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::doSingle()
{
    //! register ion implantation decay
    Int_t currImptEvt=fprimImplantEvt;
    if (fsimparms.isTdiffFromFile==1){
        ftreetdiff->GetEntry(ientrytdiff); //using tdiff from file
        while (fprimImplantEvt==currImptEvt){
            registerIonImplant();
        }
        ientrytdiff++;
    }else{//normal mode
        while (fprimImplantEvt==currImptEvt){
            registerIonImplant();
        }
    }

    //! register beta decay
    for (Int_t i=0;i<fdecaypathobj->getNMember();i++){
        //! alpha decay with probability AlphaBR (sim_neumult=100, the alpha
        //! link code in the paths); the random number is only drawn for
        //! members that have an alpha branch, so beta-only chains keep the
        //! same random sequence as before
        MemberDef* mi=fdecaypathobj->getMember(i);
        Bool_t isalpha = mi->decay_abr>0 && rseed->Rndm()*100<mi->decay_abr;
        Double_t rneu = isalpha ? 0 : rseed->Rndm();
        if (isalpha){
            mi->sim_neumult=100;
        }else if (rneu<=(fdecaypathobj->getMember(i)->decay_p0n)){//isobaric decay
            fdecaypathobj->getMember(i)->sim_neumult=0;
        }else if (rneu>fdecaypathobj->getMember(i)->decay_p0n&&rneu<=fdecaypathobj->getMember(i)->decay_p0n+fdecaypathobj->getMember(i)->decay_p1n){//decay with 1 delayed neutron
            fdecaypathobj->getMember(i)->sim_neumult=1;
        }else{//decay with 2 delayed neutron
            fdecaypathobj->getMember(i)->sim_neumult=2;
        }

        fdecaypathobj->getMember(i)->sim_ispopulated=true;
        //! isomer "population"
        if (fdecaypathobj->getMember(i)->gspatner>=0){
            Double_t rpop=rseed->Rndm();
            if (rpop<fdecaypathobj->getMember(i)->population_ratio){
                fdecaypathobj->getMember(i)->sim_ispopulated=true;
                fdecaypathobj->getMember(fdecaypathobj->getMember(i)->gspatner)->sim_ispopulated=false;
            }else{
                fdecaypathobj->getMember(i)->sim_ispopulated=false;
                fdecaypathobj->getMember(fdecaypathobj->getMember(i)->gspatner)->sim_ispopulated=true;
            }
        }
        fdecaypathobj->getMember(i)->sim_T=rseed->Exp(fdecaypathobj->getMember(i)->decay_hl/TMath::Log(2));
    }

    //! parent decay
    registerDecay(fdecaypathobj->getMember(0),-1,fdecaypathobj->getMember(0)->sim_T);

    for (Int_t i=0;i<fdecaypathobj->getNMember();i++){
        Int_t npathflow=0;
        Int_t pathid=-9999;
        Double_t decaytime=0;
        Double_t flowtime=0;
        for (Size_t j=0;j<fdecaypathobj->getMember(i)->path.size();j++){//loop all path
            Bool_t isflow=true;
            decaytime=fdecaypathobj->getMember(0)->sim_T;//first member of the path must always be 0
            for (Size_t k=1;k<fdecaypathobj->getMember(i)->path[j].size();k++){// loop all elelemts of the path
                Int_t previd=fdecaypathobj->getMember(i)->path[j][k-1];
                Int_t currid=fdecaypathobj->getMember(i)->path[j][k];
                Int_t currnneu=fdecaypathobj->getMember(i)->nneupath[j][k];
                // A path flows only if every step emitted the neutron number
                // this path assumes AND every member on it was populated (an
                // isomer and its ground state are mutually exclusive per
                // event). pgenfit/pgenfit2 used && here, so the unpopulated
                // partner of an isomer pair still decayed.
                if (currnneu!=fdecaypathobj->getMember(previd)->sim_neumult||!fdecaypathobj->getMember(currid)->sim_ispopulated) isflow=false;
                decaytime+=fdecaypathobj->getMember(currid)->sim_T;
            }
            if (isflow) {
                pathid=j;
                // keep the decay time of the path that actually flowed; pgenfit/
                // pgenfit2 used the time of the LAST path examined, which is
                // wrong for any member reachable by more than one path
                flowtime=decaytime;
                npathflow++;
            }
        }
        if (npathflow>1) {
            cout<<"something wrong"<<endl;
            exit(1);
        }

        if (npathflow==1){
            registerDecay(fdecaypathobj->getMember(i),pathid,flowtime);
            //cout<<"member: "<<fdecaypathobj->getMember(i)->id<<"\t mult="<<fdecaypathobj->getMember(i)->sim_neumult<<endl;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::registerBetaBackground(Int_t opt)
{
    //! generate beta backgroud with gausian spatial distribution
    Double_t xbetabkgg,ybetabkgg;

    if (opt==0){
        xbetabkgg=rseed->Gaus(fsimparms.xbetabkgmean,fsimparms.xbetabkgsigma);
        ybetabkgg=rseed->Gaus(fsimparms.ybetabkgmean,fsimparms.ybetabkgsigma);
        if (fsimparms.isSpatialDistFromHist!=0) {
            xbetabkgg=hhx->GetRandom();
            ybetabkgg=hhy->GetRandom();
        }
    }else{
        xbetabkgg=fsimparms.xmin+rseed->Rndm()*(fsimparms.xmax-fsimparms.xmin);
        ybetabkgg=fsimparms.ymin+rseed->Rndm()*(fsimparms.ymax-fsimparms.ymin);
    }

    if (xbetabkgg<fsimparms.xmax&&ybetabkgg<fsimparms.ymax&&xbetabkgg>=fsimparms.xmin&&ybetabkgg>=fsimparms.ymin){
        Double_t dtbetabkgg=rseed->Exp(1/fsimparms.betabkgrateg);
        ftsbetabkgg+=dtbetabkgg;

        simulationdatatype betahit;
        betahit.x=xbetabkgg;
        betahit.y=ybetabkgg;
        betahit.z=0;
        betahit.T=ftsbetabkgg;
        betahit.Tcorr=-9999;
        betahit.mode=-1;
        if (opt==0) betahit.id=-1;
        else betahit.id=-2;
        betahit.evt=-fbetabkgEvt;
        betahit.fl_n=0;
        fbetabkgEvt++;
        betaMap.insert(make_pair(betahit.T,betahit));

        if (opt==0){//only for beta decay with gaussian distrubution
            //! single neutron associated with background beta
            if (rseed->Rndm()*100<=fsimparms.percentage_bkgnb){
                Double_t dtneu=rseed->Exp(fsimparms.betaneutronmodtime/TMath::Log(2));
                simulationdatatype neuhit;
                neuhit.T=betahit.T+dtneu;
                neuhit.Tcorr=dtneu;
                neuhit.x=0.;
                neuhit.y=0.;
                neuhit.z=0.;
                neuhit.id=-2;
                neuhit.evt=-fneutronbkgEvt;
                neuhit.fl_n=0;
                neuMap.insert(make_pair(neuhit.T,neuhit));
                fneutronbkgEvt++;
            }

            if (rseed->Rndm()*100<=fsimparms.percentage_bkg2nb){
                //! 1st neutron
                Double_t dtneu=rseed->Exp(fsimparms.betaneutronmodtime/TMath::Log(2));
                simulationdatatype neuhit1;
                neuhit1.T=betahit.T+dtneu;
                neuhit1.Tcorr=dtneu;
                neuhit1.x=0.;
                neuhit1.y=0.;
                neuhit1.z=0.;
                neuhit1.id=-3;
                neuhit1.evt=-fneutronbkgEvt;
                neuhit1.fl_n=0;
                neuMap.insert(make_pair(neuhit1.T,neuhit1));
                fneutronbkgEvt++;

                //! 2nd neutron
                dtneu=rseed->Exp(fsimparms.betaneutronmodtime/TMath::Log(2));
                simulationdatatype neuhit2;
                neuhit2.T=betahit.T+dtneu;
                neuhit2.Tcorr=dtneu;
                neuhit2.x=0.;
                neuhit2.y=0.;
                neuhit2.z=0.;
                neuhit2.id=-3;
                neuhit2.evt=-fneutronbkgEvt;
                neuhit2.fl_n=0;
                neuMap.insert(make_pair(neuhit2.T,neuhit2));
                fneutronbkgEvt++;
            }
        }
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::registerNeutronBackground(Int_t opt)
{
    Double_t dtneubkg;
    if (opt==0){
        dtneubkg=rseed->Exp(1/fsimparms.neurndbkgrate);
        ftsneutronbkgg+=dtneubkg;
        simulationdatatype neuhit;
        neuhit.x=0;
        neuhit.y=0;
        neuhit.z=0;
        neuhit.T=ftsneutronbkgg;
        neuhit.Tcorr=-1;
        neuhit.id=-4;
        neuhit.evt=-fneutronbkgEvt;
        neuhit.fl_n=0;
        neuMap.insert(make_pair(neuhit.T,neuhit));
        fneutronbkgEvt++;
    }else{
        //! timming of the unknown source
        dtneubkg=rseed->Exp(1/fsimparms.r2neurndbkgrate);
        ftsneutronbkgg+=dtneubkg;
        //! correlated neutron

        //!1st neutron
        Double_t dtneu=rseed->Exp(fsimparms.beamneutronmodtime/TMath::Log(2));
        simulationdatatype neuhit1;
        neuhit1.x=0;
        neuhit1.y=0;
        neuhit1.z=0;
        neuhit1.T=ftsneutronbkgg+dtneu;
        neuhit1.Tcorr=-1;
        neuhit1.id=-5;
        neuhit1.evt=-fneutronbkgEvt;
        neuhit1.fl_n=0;
        neuMap.insert(make_pair(neuhit1.T,neuhit1));
        fneutronbkgEvt++;

        //!2nd neutron
        dtneu=rseed->Exp(fsimparms.beamneutronmodtime/TMath::Log(2));
        simulationdatatype neuhit2;
        neuhit2.x=0;
        neuhit2.y=0;
        neuhit2.z=0;
        neuhit2.T=ftsneutronbkgg+dtneu;
        neuhit2.Tcorr=-1;
        neuhit2.id=-5;
        neuhit2.evt=-fneutronbkgEvt;
        neuhit2.fl_n=0;
        neuMap.insert(make_pair(neuhit2.T,neuhit2));
        fneutronbkgEvt++;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::runSimulation()
{
    //! begin
    fprimImplantT=fsimparms.tsoffset;

    fprimImplantEvt=0;
    fbetaEvt=0;
    fneuEvt=0;
    cout<<"--------"<<endl;
    if (fsimparms.isTdiffFromFile==1){//if tdiff is from file
        bookTDiffData();
        while (ientrytdiff<fsimparms.nTdiffEntries+fsimparms.iTdiffEntryBegin){
            doSingle();
        }
    }else{//normal mode
        while (fprimImplantT<fsimparms.BeamTime+fsimparms.tsoffset){
            doSingle();
        }
    }
    cout<<fprimImplantEvt<<endl;
    cout<<fbetaEvt<<endl;
    cout<<fneuEvt<<endl;

    //! gaussian distrubution beta bkg
    ftsbetabkgg=fsimparms.tsoffset;
    fbetabkgEvt=0;
    while (ftsbetabkgg<fprimImplantT*1.1){
        registerBetaBackground();
    }
    cout<<fbetabkgEvt<<endl;

    //! uniform distrubution beta bkg
    ftsbetabkgg=fsimparms.tsoffset;
    fbetabkgEvt=0;
    while (ftsbetabkgg<fprimImplantT*1.1){
        registerBetaBackground(1);
    }
    cout<<fbetabkgEvt<<endl;


    cout<<fneutronbkgEvt<<endl;
    //! single neutron bkg
    ftsneutronbkgg=fsimparms.tsoffset;
    fneutronbkgEvt=0;
    while (ftsneutronbkgg<fprimImplantT*1.01){
        registerNeutronBackground();
    }
    cout<<fneutronbkgEvt<<endl;

    //! two neutrons bkg
    ftsneutronbkgg=fsimparms.tsoffset;
    fneutronbkgEvt=0;
    while (ftsneutronbkgg<fprimImplantT*1.01){
        registerNeutronBackground(1);
    }
    cout<<fneutronbkgEvt<<endl;

    cout<<"------------"<<endl;
    cout<<ionMap.size()<<endl;
    cout<<betaMap.size()<<endl;
    cout<<neuMap.size()<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::fillTreeData()
{
    //! ion
    for (ionMap_it=ionMap.begin();ionMap_it!=ionMap.end();ionMap_it++){
        simulationdatatype ionhit = ionMap_it->second;
        copydata(fionData,ionhit);
        ftreeion->Fill();
    }

    //! beta
    for (betaMap_it=betaMap.begin();betaMap_it!=betaMap.end();betaMap_it++){
        simulationdatatype betahit = betaMap_it->second;
        copydata(fbetaData,betahit);
        ftreebeta->Fill();
    }

    //! neutron
    for (neuMap_it=neuMap.begin();neuMap_it!=neuMap.end();neuMap_it++){
        simulationdatatype neuhit = neuMap_it->second;
        copydata(fneutronData,neuhit);
        ftreeneu->Fill();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void simulation::correlateData()
{
    //! setup histograms
    fsim_hdecay=new TH1F("hdecay","hdecay",2000,-fionbetawindowlow,fionbetawindowup);
    fsim_hdecay1nbwd=new TH1F("hdecay1nbwd","hdecay1nbwd",2000,-fionbetawindowlow,fionbetawindowup);
    fsim_hdecaygt0nbwd=new TH1F("hdecaygt0nbwd","hdecaygt0nbwd",2000,-fionbetawindowlow,fionbetawindowup);
    fsim_hdecay2nbwd=new TH1F("hdecay2nbwd","decay2nbwd",2000,-fionbetawindowlow,fionbetawindowup);


    Int_t k=0;
    Int_t ncorrneuf=0;
    Int_t ncorrneub=0;
    Int_t ncorrion=0;

    //! beta
    for (betaMap_it=betaMap.begin();betaMap_it!=betaMap.end();betaMap_it++){
        if (k%10000==0&&k>0) cout<<(Double_t)k/(Double_t)betaMap.size()*100.<<" % completed - "<<ncorrion<<"/"<<ncorrneuf<<"/"<<ncorrneub<<"\r"<<flush;
        double ts=(double)betaMap_it->first;
        simulationdatatype betahit = betaMap_it->second;
        copydata(fbetaData,betahit);

        //! with neutron forward
        fcorrNeutronData_fw.mult=0;
        double ts1 = ts - 0;
        double ts2 = ts + fwindowbetaneutronlow;
        double corrts = 0;
        double check_time = 0;
        neuMap_it = neuMap.lower_bound(ts1);
        while(neuMap_it!=neuMap.end()&&neuMap_it->first<ts2){
            corrts = neuMap_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                simulationdatatype neuhit = neuMap_it->second;
                copydatamult(fcorrNeutronData_fw,neuhit,fcorrNeutronData_fw.mult);
                fcorrNeutronData_fw.mult++;
            }
            neuMap_it++;
        }

        //! with neutron backward
        fcorrNeutronData_bw.mult=0;
        ts1 = ts - fwindowbetaneutronup;
        ts2 = ts + 0;
        corrts = 0;
        check_time = 0;
        neuMap_it = neuMap.lower_bound(ts1);
        while(neuMap_it!=neuMap.end()&&neuMap_it->first<ts2){
            corrts = neuMap_it->first;
            if (corrts!=check_time){
                check_time=corrts;
                simulationdatatype neuhit = neuMap_it->second;
                copydatamult(fcorrNeutronData_bw,neuhit,fcorrNeutronData_bw.mult);
                fcorrNeutronData_bw.mult++;
            }
            neuMap_it++;
        }

        //! with implantation
        int nion=0;
        ts1 = ts - fionbetawindowup;
        ts2 = ts + fionbetawindowlow;
        corrts = 0;
        check_time = 0;
        ionMap_it = ionMap.lower_bound(ts1);
        while(ionMap_it!=ionMap.end()&&ionMap_it->first<ts2){
            corrts = ionMap_it->first;
            simulationdatatype ionhit=ionMap_it->second;

            if (corrts!=check_time&&ionhit.z==betahit.z){
                if  (!((betahit.x-ionhit.x>=-fdeltaxy)&&(betahit.x-ionhit.x<=fdeltaxy)&&(betahit.y-ionhit.y>=-fdeltaxy)&&(betahit.y-ionhit.y<=fdeltaxy))){
                    ionMap_it++;
                    continue;
                }
                copydata(fionData,ionhit);
                //! fill tree here!
                ftreecorr->Fill();

                fmlh_t=ts-corrts;
                fmlh_mult=fcorrNeutronData_fw.mult;
                fmlh_ionT=corrts;
                fmlh_alpha=(betahit.mode==10) ? 1 : 0;
                ftreemlh->Fill();
                fmlhbw_t=fmlh_t;
                fmlhbw_mult=fcorrNeutronData_bw.mult;
                ftreemlhbw->Fill();
                fsim_hdecay->Fill(fmlh_t);
                if (fcorrNeutronData_bw.mult==1) fsim_hdecay1nbwd->Fill(fmlh_t);
                if (fcorrNeutronData_bw.mult>0) fsim_hdecaygt0nbwd->Fill(fmlh_t);
                if (fcorrNeutronData_bw.mult==2) fsim_hdecay2nbwd->Fill(fmlh_t);
                nion++;
            }
            ionMap_it++;
        }
        if (fcorrNeutronData_fw.mult>0) ncorrneuf++;
        if (fcorrNeutronData_bw.mult>0) ncorrneub++;
        if (nion>0) ncorrion++;
        k++;
    }
}

