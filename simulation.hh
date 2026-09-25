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

//! \file simulation.hh
//! \brief Definition of the simulation class

#ifndef simulation_h
#define simulation_h 1

#include "decaypath.hh"
#include "TRandom3.h"


#include <map>
#include <vector>
#include <iostream>

#include "TH1.h"
#include "TTree.h"
#include "TParameter.h"


typedef struct {
    double T; 	 // Calibrated time
    double Tcorr; //correlated time
    double x,y,z;// number of pixel for AIDA, or number of tube for BELEN
    int evt;

    int mode;//decay mode:b0n = 0, b1n = 1, b2n = 2

    int id;//id of decaying iso
    int ifl;//id of flow
    int fl_n;//number of members
    int fl_i[kmaxnri];//index of member
} simulationdatatype;

#define kmaxnmult 100
typedef struct {
    int mult;
    double T[kmaxnmult]; 	 // Calibrated time
    double Tcorr[kmaxnmult]; //correlated time
    double x[kmaxnmult],y[kmaxnmult],z[kmaxnmult];// number of pixel for AIDA, or number of tube for BELEN
    int evt[kmaxnmult];

    int mode[kmaxnmult];//decay mode:b0n = 0, b1n = 1, b2n = 2

    int id[kmaxnmult];//id of decaying iso
    int ifl[kmaxnmult];//id of flow
    int fl_n[kmaxnmult];//number of members
} simulationdatatypemult;

typedef struct{
    Double_t betaeff=50.; //percentage of an isotope
    Double_t neutroneff=68*0.95; //neutron detection efficiency in percentage

    Double_t deltaxy=2.5;
    Double_t dxbetamean=0;Double_t dxbetasigma=deltaxy/4.;
    Double_t dybetamean=0;Double_t dybetasigma=deltaxy/4.;
    Double_t betaneutronmodtime=0.000021; //21 us moderationtime
    Double_t beamneutronmodtime=0.000027; //27 us moderationtime

    Double_t ximpmean=64;Double_t ximpsigma=32;
    Double_t yimpmean=64;Double_t yimpsigma=32;
    Double_t xbetabkgmean=64;Double_t xbetabkgsigma=32;
    Double_t ybetabkgmean=64;Double_t ybetabkgsigma=32;

    Double_t xmin=0;Double_t xmax=128;
    Double_t ymin=0;Double_t ymax=128;
    Double_t tsoffset=3600; //! there is a bug on this offset (don't now why -> to be fixed later)
    //! ion-decay correlation window (s): decays from ionbetawindowlow before
    //! to ionbetawindowup after each implant are correlated (x = t_decay-t_ion)
    Double_t ionbetawindowlow=10;
    Double_t ionbetawindowup=20;
    Double_t neuwbeamperctg=40.;


    Int_t isTdiffFromFile = 0;
    Int_t nTdiffEntries = 20000;
    Int_t iTdiffEntryBegin = 110000;
    Int_t isFixImplantPosition = 0;
    Int_t isSpatialDistFromHist = 0;

    //! some default rate values
    Double_t BeamTime=3600;
    Double_t rate=141.;
    Double_t isoperctg=100.;
    Double_t betabkgrateg=1.;
    Double_t betabkgrateu=1.;
    Double_t neurndbkgrate=1.;//random single neutron background  per second
    Double_t r2neurndbkgrate=1.;//random 2 neutron background per second

    Double_t percentage_bkgnb=7.;
    Double_t percentage_bkg2nb=1.;

    //! alpha decays (pgenfit's alpha extension): detection efficiency in
    //! percent and position spread around the implant, like the beta ones
    Double_t alphaeff=100.;
    Double_t deltaxyalpha=1.;
    Double_t dxalphamean=0;Double_t dxalphasigma=deltaxyalpha/8.;
    Double_t dyalphamean=0;Double_t dyalphasigma=deltaxyalpha/8.;


} simulationparmstype;

class simulation
{
  public:
    simulation(char* inputParm);
    simulation(decaypath* path);
    virtual ~simulation();
    void Init(char* inputParms);
    void Init(decaypath* path){fdecaypathobj=path;}
    void setRandomSeed(Int_t seedno=4357){rseed=new TRandom3(seedno);}

    void BookSimulationTree();
    void BookCorrelationTree();

    void readSimulationParameters(char* inputfile);

    void printPathMembers();

    void registerIonImplant();
    void registerDecay(MemberDef * decaymember, Int_t pathid, Double_t decaytime);
    void registerBetaBackground(Int_t opt=0);
    void registerNeutronBackground(Int_t opt=0);

    void doSingle();

    void runSimulation();

    void fillTreeData();

    void correlateData();

    void bookTDiffData();

    TTree* getIonSimulationTree(){return ftreeion;}
    TTree* getBetaSimulationTree(){return ftreebeta;}
    TTree* getNeutronSimulationTree(){return ftreeneu;}

    TTree* getMLHTree(){return ftreemlh;}
    TTree* getMLHTreeBackward(){return ftreemlhbw;}
    void writeMLHHistos(){
        //! record the correlation window so a fit can check its time range
        TParameter<double>("ionbetawindowlow",fionbetawindowlow).Write();
        TParameter<double>("ionbetawindowup",fionbetawindowup).Write();
        fsim_hdecay->Write();
        fsim_hdecay1nbwd->Write();
        fsim_hdecaygt0nbwd->Write();
        fsim_hdecay2nbwd->Write();
    }
    TTree* getCorrelationTree(){return ftreecorr;}

    void copydata(simulationdatatype& datades,simulationdatatype &datasrc){
        datades.T=datasrc.T;
        datades.Tcorr=datasrc.Tcorr;
        datades.x=datasrc.x;
        datades.y=datasrc.y;
        datades.z=datasrc.z;
        datades.evt=datasrc.evt;

        datades.mode=datasrc.mode;

        datades.id=datasrc.id;
        datades.ifl=datasrc.ifl;
        datades.fl_n=datasrc.fl_n;
        for (Int_t i=0;i<datasrc.fl_n;i++) datades.fl_i[i]=datasrc.fl_i[i];
    }

    void copydatamult(simulationdatatypemult& datades,simulationdatatype &datasrc,Int_t i){
        datades.T[i]=datasrc.T;
        datades.Tcorr[i]=datasrc.Tcorr;
        datades.x[i]=datasrc.x;
        datades.y[i]=datasrc.y;
        datades.z[i]=datasrc.z;
        datades.evt[i]=datasrc.evt;

        datades.mode[i]=datasrc.mode;

        datades.id[i]=datasrc.id;
        datades.ifl[i]=datasrc.ifl;
        datades.fl_n[i]=datasrc.fl_n;
    }

    void resetdata(simulationdatatype& datades){
        datades.T=-9999;
        datades.Tcorr=-9999;
        datades.x=-9999;
        datades.y=-9999;
        datades.z=-9999;
        datades.evt=-9999;

        datades.mode=-9999;

        datades.id=-9999;
        datades.ifl=-9999;
        datades.fl_n=0;

        for (Int_t i=0;i<kmaxnri;i++) datades.fl_i[i]=-9999;
    }
 private:
    decaypath* fdecaypathobj;

    TTree* ftreeion;
    TTree* ftreebeta;
    TTree* ftreeneu;


    //! stuff for simulation
    simulationdatatype fionData;
    simulationdatatype fbetaData;
    simulationdatatype fneutronData;
    simulationdatatypemult fcorrNeutronData_fw;
    simulationdatatypemult fcorrNeutronData_bw;
    TTree* ftreecorr;
    TTree* ftreemlh;
    TTree* ftreemlhbw;
    Double_t fmlh_t;
    Int_t fmlh_mult;
    Double_t fmlh_ionT; //! absolute ion implantation timestamp, for
                         //! analysis-side run-boundary cuts (excludes
                         //! implants too close to the start/end of the
                         //! beam-on window or a beam interruption, where
                         //! the accidental background rate is still
                         //! ramping rather than stationary -- same
                         //! technique as the BRIKEN commissioning paper's
                         //! beam-interruption exclusion).
    Int_t fmlh_alpha;   //! 1 if the correlated decay hit is an alpha (mode 10)
    Double_t fmlhbw_t;
    Int_t fmlhbw_mult;
    TH1F* fsim_hdecay;
    TH1F* fsim_hdecay1nbwd;
    TH1F* fsim_hdecaygt0nbwd;
    TH1F* fsim_hdecay2nbwd;


    Double_t fdeltaxy;
    Double_t fionbetawindowlow;
    Double_t fionbetawindowup;
    Double_t fwindowbetaneutronlow;
    Double_t fwindowbetaneutronup;

    std::multimap < double, simulationdatatype > ionMap;
    std::multimap < double, simulationdatatype >::iterator ionMap_it;
    std::multimap < double, simulationdatatype > betaMap;
    std::multimap < double, simulationdatatype >::iterator betaMap_it;
    std::multimap < double, simulationdatatype > neuMap;
    std::multimap < double, simulationdatatype >::iterator neuMap_it;

    Int_t fripathToId[kmaxnri][kmaxpaths];
    TRandom3 *rseed;

    simulationparmstype fsimparms;
    TH1F* hhx;
    TH1F* hhy;
    TH1F* hhimpx;
    TH1F* hhimpy;
    TH1F* hhdx;
    TH1F* hhdy;

    Int_t fprimImplantEvt;
    Double_t fprimImplantT;
    Int_t fbetaEvt;
    Double_t ftsbeta;

    Int_t fneuEvt;

    Double_t ftsbetabkgg;
    Double_t ftsneutronbkgg;

    Int_t fbetabkgEvt;
    Int_t fneutronbkgEvt;

    Double_t fximp;
    Double_t fyimp;

    Double_t ftdiff;
    TTree* ftreetdiff;
    Long64_t ientrytdiff;


};

#endif
