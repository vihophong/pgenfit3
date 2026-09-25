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
/// \file decaypath.cc
/// \brief Implementation of the decaypath class

#include "decaypath.hh"
#include <iostream>
#include <fstream>


#include "TROOT.h"
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
#include <cmath>
#include <functional>
#include <vector>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

decaypath::decaypath():
    flistofdecaymember()
{
    fdecaypath=new path();
    finputParms=new char[1000];
    flistofdecaymember.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

decaypath::~decaypath()
{
    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        delete *flistofdecaymember_it;
    }
    flistofdecaymember.clear();
    delete finputParms;
    delete fdecaypath;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::Init(char* inputParms)
{
    sprintf(finputParms,"%s",inputParms);
    std::clog<< __PRETTY_FUNCTION__<<" read input files:"<<
               std::endl<<finputParms<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::ProcessMember(MemberDef *obj)
{
    //char temp[100];
    //sprintf(temp,"^{%i}",obj->n+obj->z);
    //obj->name=obj->name.Prepend(temp);

    //! further processing
    if (obj->decay_hl<0){
        obj->decay_hl=-obj->decay_hl;
        obj->is_decay_hl_fix=0;
    }else{//exclude decay half-life =0;
        obj->is_decay_hl_fix=1;
    }
    if (obj->decay_hlerr==0){//if error ==0-> fix amd no contrain on parameter
        obj->is_decay_hl_fix=2;
    }

    // convert half-life into activity
    obj->decay_lamda=log(2)/obj->decay_hl;
    obj->decay_lamdaerr=log(2)/obj->decay_hl/obj->decay_hl*obj->decay_hlerrhi;
    obj->decay_lamdaerrhi=log(2)/obj->decay_hl/obj->decay_hl*obj->decay_hlerr;
    obj->decay_lamdalow=log(2)/obj->decay_hlup;
    obj->decay_lamdaup=log(2)/obj->decay_hllow;
    obj->is_decay_lamda_fix=obj->is_decay_hl_fix;

    if (obj->decay_p1n<0){
        obj->decay_p1n=-obj->decay_p1n;
        obj->is_decay_p1n_fix=0;
    }else if (obj->decay_p1n==0){//exclude decay with pn = 0;
        obj->is_decay_p1n_fix=2;
    }else{
        obj->is_decay_p1n_fix=1;
    }

    if (obj->decay_p2n<0){
        obj->decay_p2n=-obj->decay_p2n;
        obj->is_decay_p2n_fix=0;
    }else if (obj->decay_p2n==0){//exclude decay p2n =0;
        obj->is_decay_p2n_fix=2;
    }else{
        obj->is_decay_p2n_fix=1;
    }

    // alpha branch stays in percent here (as in pgenfit)
    if (obj->decay_abr<0){
        obj->decay_abr=-obj->decay_abr;
        obj->is_decay_abr_fix=0;
    }else if (obj->decay_abr==0){//no alpha decay
        obj->is_decay_abr_fix=2;
    }else{
        obj->is_decay_abr_fix=1;
    }

    if (obj->neueff<0){
        obj->neueff=-obj->neueff;
        obj->is_neueff_fix=0;
    }else{
        obj->is_neueff_fix=1;
    }
    if (obj->neuefferr==0){//if error ==0-> fix amd no contrain on parameter
        obj->is_neueff_fix=2;
    }

    // convert pn in % to pn in 1
    obj->decay_p1n=obj->decay_p1n/100;
    obj->decay_p1nerr=obj->decay_p1nerr/100;
    obj->decay_p1nerrhi=obj->decay_p1nerrhi/100;
    obj->decay_p1nlow=obj->decay_p1nlow/100;
    obj->decay_p1nup=obj->decay_p1nup/100;

    obj->decay_p2n=obj->decay_p2n/100;
    obj->decay_p2nerr=obj->decay_p2nerr/100;
    obj->decay_p2nerrhi=obj->decay_p2nerrhi/100;
    obj->decay_p2nlow=obj->decay_p2nlow/100;
    obj->decay_p2nup=obj->decay_p2nup/100;

    obj->decay_p0n=1-obj->decay_p1n-obj->decay_p2n;
    obj->decay_p0nerr=sqrt(obj->decay_p1nerr*obj->decay_p1nerr+obj->decay_p2nerr*obj->decay_p2nerr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::CopyMember(MemberDef *source, MemberDef *destination)
{
    destination-> id = source->  id;
    destination-> z = source->  z;
    destination-> n = source->  n;
    destination-> name = source->  name;

    destination-> decay_hl = source->  decay_hl;
    destination-> decay_lamda = source->  decay_lamda;

    destination-> decay_p0n = source->  decay_p0n;
    destination-> decay_p1n = source->  decay_p1n;
    destination-> decay_p2n = source->  decay_p2n;

    destination-> decay_hlerr = source->  decay_hlerr;
    destination-> decay_lamdaerr = source->  decay_lamdaerr;
    destination-> decay_p0nerr = source->  decay_p0nerr;
    destination-> decay_p1nerr = source->  decay_p1nerr;
    destination-> decay_p2nerr = source->  decay_p2nerr;

    destination-> decay_hlerrhi = source->  decay_hlerrhi;
    destination-> decay_lamdaerrhi = source->  decay_lamdaerrhi;
    destination-> decay_p0nerrhi = source->  decay_p0nerrhi;
    destination-> decay_p1nerrhi = source->  decay_p1nerrhi;
    destination-> decay_p2nerrhi = source->  decay_p2nerrhi;

    destination-> decay_hlup = source->  decay_hlup;
    destination-> decay_lamdaup = source->  decay_lamdaup;
    destination-> decay_p0nup = source->  decay_p0nup;
    destination-> decay_p1nup = source->  decay_p1nup;
    destination-> decay_p2nup = source->  decay_p2nup;

    destination-> decay_hllow = source->  decay_hllow;
    destination-> decay_lamdalow = source->  decay_lamdalow;
    destination-> decay_p0nlow = source->  decay_p0nlow;
    destination-> decay_p1nlow = source->  decay_p1nlow;
    destination-> decay_p2nlow = source->  decay_p2nlow;

    destination-> population_ratio = source->  population_ratio;
    destination-> population_ratioerr = source->  population_ratioerr;
    destination-> population_ratioup = source->  population_ratioup;
    destination-> population_ratiolow = source->  population_ratiolow;

    destination-> neueff = source->  neueff;
    destination-> neuefferr = source->  neuefferr;
    destination-> neuefferrhi = source->  neuefferrhi;
    destination-> neueffup = source->  neueffup;
    destination-> neuefflow = source->  neuefflow;

    destination-> is_decay_hl_fix = source->  is_decay_hl_fix;
    destination-> is_decay_lamda_fix = source->  is_decay_lamda_fix;
    destination-> is_decay_p0n_fix = source->  is_decay_p0n_fix;
    destination-> is_decay_p1n_fix = source->  is_decay_p1n_fix;
    destination-> is_decay_p2n_fix = source->  is_decay_p2n_fix;
    destination-> is_neueff_fix = source->  is_neueff_fix;

    destination-> gspatner = source->  gspatner;

    destination-> decay_abr = source->  decay_abr;
    destination-> decay_abrerr = source->  decay_abrerr;
    destination-> decay_abrerrhi = source->  decay_abrerrhi;
    destination-> decay_abrlow = source->  decay_abrlow;
    destination-> decay_abrup = source->  decay_abrup;
    destination-> is_decay_abr_fix = source->  is_decay_abr_fix;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::makePath()
{
    std::clog<< __PRETTY_FUNCTION__<<std::endl;
    //! reading input file and put informations to list of members
    Int_t id=0;
    std::string line;
    std::ifstream infile(finputParms);
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        if (line.find_first_not_of(" \t\r")==std::string::npos) continue; // blank line
        if (line[0]=='#') continue;
        Int_t a;
        // decay properies
        MemberDef* obj=new MemberDef();
        obj->id=id;
        //! read info without isomer
//        if (!(iss >> obj->name >> obj->z >> a >> obj->decay_hl >> obj->decay_hlerr >> obj->decay_hllow >> obj->decay_hlup >>
//              obj->decay_p1n >> obj->decay_p1nerr >> obj->decay_p1nlow >> obj->decay_p1nup >>
//              obj->decay_p2n >> obj->decay_p2nerr >> obj->decay_p2nlow >> obj->decay_p2nup)) break;
        obj->population_ratio = 1;
        obj->population_ratioerr = 0;
        obj->population_ratiolow = 0;
        obj->population_ratioup = 2;
        obj->is_population_ratio_fix = 2;
        obj->decay_abr = 0; obj->decay_abrerr = 0; obj->decay_abrerrhi = 0;
        obj->decay_abrlow = 0; obj->decay_abrup = 200;

        //! read info with isomer
        if (!(iss >> obj->name >> obj->z >> a >> obj->decay_hl >> obj->decay_hlerr >> obj->decay_hlerrhi >> obj->decay_hllow >> obj->decay_hlup >>
              obj->decay_p1n >> obj->decay_p1nerr >> obj->decay_p1nerrhi >> obj->decay_p1nlow >> obj->decay_p1nup >>
              obj->decay_p2n >> obj->decay_p2nerr >> obj->decay_p2nerrhi >> obj->decay_p2nlow >> obj->decay_p2nup >>
              obj->neueff >> obj->neuefferr >> obj->neuefferrhi >> obj->neuefflow >> obj->neueffup)) break;
        obj->n=a-obj->z;

        obj->gspatner = -1;

        //! isomeric state
        Bool_t flagisomer = false;
        std::string tempstring(obj->name.Data());
        MemberDef* objisomer;
        if (tempstring.back()=='*') {
            objisomer=new MemberDef();
            CopyMember(obj,objisomer);
            objisomer->id = obj->id + 1;
            cout<<"ISOMER of "<<obj->name<<endl;
            if (!(iss >> objisomer->population_ratio >> objisomer->population_ratioerr >> objisomer->population_ratiolow >> objisomer->population_ratioup >>
                  objisomer->decay_hl >> objisomer->decay_hlerr >> objisomer->decay_hlerrhi >> objisomer->decay_hllow >> objisomer->decay_hlup >>
                  objisomer->decay_p1n >> objisomer->decay_p1nerr >> objisomer->decay_p1nerrhi >> objisomer->decay_p1nlow >> objisomer->decay_p1nup >>
                  objisomer->decay_p2n >> objisomer->decay_p2nerr >> objisomer->decay_p2nerrhi >> objisomer->decay_p2nlow >> objisomer->decay_p2nup >>
                  objisomer->neueff >> objisomer->neuefferr >> objisomer->neuefferrhi >> objisomer->neuefflow >> objisomer->neueffup)) break;
            obj->name = TString(tempstring.substr(0,tempstring.length()-1).data());

            if (objisomer->population_ratio>0){
                obj->is_population_ratio_fix = 1;
                objisomer->is_population_ratio_fix = 1;
            }else{
                objisomer->population_ratio = -objisomer->population_ratio;
                obj->is_population_ratio_fix = 0;
                objisomer->is_population_ratio_fix = 0;
#ifdef ISOMER_SUM_UNITY
                objisomer->is_population_ratio_fix = 1;
#endif
            }

            obj->population_ratio = 1 - objisomer->population_ratio;
            obj->population_ratioerr = objisomer->population_ratioerr;
            obj->population_ratioup = 1 - objisomer->population_ratiolow;
            obj->population_ratiolow = 1 - objisomer->population_ratioup;

            objisomer->gspatner = obj->id;
            flagisomer = true;
        }

        //! optional alpha block (pgenfit format): "AlphaBR errLo errHi low up"
        //! in percent after the base (and isomer) columns; for an isomer row a
        //! second block may follow for the isomer itself. Absent = no alpha.
        auto readAlpha=[&](MemberDef* m){
            std::streampos pos = iss.tellg();
            Double_t abr,err,errhi,lo,up;
            if (iss >> abr >> err >> errhi >> lo >> up){
                m->decay_abr=abr; m->decay_abrerr=err; m->decay_abrerrhi=errhi;
                m->decay_abrlow=lo; m->decay_abrup=up;
            }else{
                iss.clear(); iss.seekg(pos);
            }
        };
        readAlpha(obj);
        if (flagisomer){
            objisomer->decay_abr=0; objisomer->decay_abrerr=0; objisomer->decay_abrerrhi=0;
            objisomer->decay_abrlow=0; objisomer->decay_abrup=200;
            readAlpha(objisomer);
        }

        //! ground state
        ProcessMember(obj);
        flistofdecaymember.emplace(flistofdecaymember.end(),obj);
        id++;

        if (flagisomer){
            //! isomeric state
            ProcessMember(objisomer);
            flistofdecaymember.emplace(flistofdecaymember.end(),objisomer);
            id++;
        }
    }

    //! Paths from the parent (member 0) to every member, by depth-first
    //! search over the decay links:
    //!   beta  (Z+1): N-1 -> 0n, N-2 -> 1n, N-3 -> 2n
    //!   alpha (Z-2, N-2): code 100
    //! pgenfit/pgenfit2 built paths in list order, which needs every
    //! predecessor listed before its daughters; that fails for an alpha
    //! daughter (lower Z than its parent). For beta-only chains listed in
    //! increasing Z the paths, and their order, are the same as before.
    std::vector<MemberDef*> mem(flistofdecaymember.begin(), flistofdecaymember.end());
    const Int_t nmem = mem.size();
    auto linkType=[](MemberDef* from, MemberDef* to)->Int_t{
        Int_t dz=to->z-from->z, dn=to->n-from->n;
        if (dz==1){
            if (dn==-1) return 0;
            if (dn==-2) return 1;
            if (dn==-3) return 2;
        }
        if (dz==-2 && dn==-2) return 100;
        return -1;
    };
    std::vector<Int_t> state(nmem,0); // 0 todo, 1 in progress, 2 done
    std::function<void(Int_t)> build=[&](Int_t i){
        if (state[i]==2) return;
        if (state[i]==1){ std::cerr<<"decay network has a cycle at "<<mem[i]->name<<std::endl; exit(1); }
        state[i]=1;
        if (i!=0){ // the parent is the root: no path leads to it
            for (Int_t j=0;j<nmem;j++){
                if (j==i) continue;
                Int_t t=linkType(mem[j],mem[i]);
                if (t<0) continue;
                build(j);
                if (j==0){
                    mem[i]->path.push_back({0,mem[i]->id});
                    mem[i]->nneupath.push_back({-1,t});
                }else{
                    for (size_t k=0;k<mem[j]->path.size();k++){
                        std::vector<Int_t> row=mem[j]->path[k]; row.push_back(mem[i]->id);
                        std::vector<Int_t> nrow=mem[j]->nneupath[k]; nrow.push_back(t);
                        mem[i]->path.push_back(row);
                        mem[i]->nneupath.push_back(nrow);
                    }
                }
            }
        }
        state[i]=2;
    };
    for (Int_t i=0;i<nmem;i++) build(i);

    //! A path carries no flow if one of its links has a branching fixed at
    //! zero (then its coefficient is identically 0 and it is left out).
    auto linkFlows=[](MemberDef* from, Int_t t)->Bool_t{
        const Bool_t alphaOnly = from->decay_abr>=100 && from->is_decay_abr_fix!=0;
        if (t==100) return from->decay_abr>0;
        if (alphaOnly) return false;
        if (t==0) return !(from->decay_p1n+from->decay_p2n>=1 && from->is_decay_p1n_fix!=0 && from->is_decay_p2n_fix!=0);
        if (t==1) return from->decay_p1n!=0;
        return from->decay_p2n!=0;
    };

    fdecaypath->npaths=0;
    fdecaypath->nri=nmem;
    for (Int_t i=0;i<nmem;i++){
        for (size_t k=0;k<mem[i]->path.size();k++){
            const Int_t ip=fdecaypath->npaths;
            if (ip>=kmaxpaths){ std::cerr<<"too many decay paths (kmaxpaths="<<kmaxpaths<<")"<<std::endl; exit(1); }
            const std::vector<Int_t>& row=mem[i]->path[k];
            const std::vector<Int_t>& nrow=mem[i]->nneupath[k];
            Bool_t isflow=true;
            fdecaypath->ndecay[ip]=row.size();
            for (size_t j=0;j<row.size();j++){
                fdecaypath->decaymap[ip][j]=row[j];
                if (j+1<row.size()){
                    fdecaypath->nneu[ip][j]=nrow[j+1];
                    if (!linkFlows(mem[row[j]],nrow[j+1])) isflow=false;
                }
            }
            fdecaypath->ispathhasflow[ip]=isflow;
            fdecaypath->npaths++;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::writePath()
{
    std::ofstream pathfile("path.txt");
    pathfile<<fdecaypath->nri<<std::endl;
    pathfile<<fdecaypath->npaths<<std::endl;

    for (int i=0;i<fdecaypath->npaths;i++){
        pathfile<<fdecaypath->ndecay[i]<<std::endl;
        pathfile<<fdecaypath->ispathhasflow[i]<<std::endl;
        for (int j=0;j<fdecaypath->ndecay[i];j++){
            pathfile<<fdecaypath->decaymap[i][j]<<"\t"<<fdecaypath->nneu[i][j]<<std::endl;
        }
    }
    //! stuffs for isomers
    Int_t nisomers=0;
    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        if ((*flistofdecaymember_it)->gspatner!=-1){
            nisomers++;
        }
    }
    pathfile<<nisomers<<std::endl;
    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        if ((*flistofdecaymember_it)->gspatner!=-1){
            pathfile<<(*flistofdecaymember_it)->gspatner<<" "<<(*flistofdecaymember_it)->id<<std::endl;
        }
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::printPath()
{
    std::clog<< __PRETTY_FUNCTION__<<std::endl;
    Int_t npaths=0;

    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        cout<<"********* Go for Isotope "<<(*flistofdecaymember_it)->id<<" ("<<(*flistofdecaymember_it)->name<<")"<<endl;
        for (Size_t i=0;i<(*flistofdecaymember_it)->path.size();i++){
            cout<<"row "<<i<<" = ";
            std::list<MemberDef*>::iterator listofdecaymember_it2;
            for (Size_t j=0;j<(*flistofdecaymember_it)->path[i].size();j++){
                cout<<(*flistofdecaymember_it)->path[i][j]<<" ";
                for (listofdecaymember_it2 = flistofdecaymember.begin(); listofdecaymember_it2 != flistofdecaymember.end(); listofdecaymember_it2++)
                {
                    if ((*flistofdecaymember_it)->path[i][j]==(*listofdecaymember_it2)->id){
                        cout<<"("<<(*listofdecaymember_it2)->name<<")\t to \t";
                    }
                }
            }
            cout<<" end"<<endl;
            npaths++;
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::printMember()
{
    std::clog<< __PRETTY_FUNCTION__<<std::endl;
    //! display information of the paths to decay
    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        cout<<(*flistofdecaymember_it)->id<<"\t"<<(*flistofdecaymember_it)->name<<"\t"<<(*flistofdecaymember_it)->z<<"\t"<<
              (*flistofdecaymember_it)->n<<"\t"<<(*flistofdecaymember_it)->n+(*flistofdecaymember_it)->z<<"\t"<<(*flistofdecaymember_it)->decay_hl<<"\t"<<(*flistofdecaymember_it)->decay_lamda<<"\t"<<
              (*flistofdecaymember_it)->decay_p1n<<"\t"<<(*flistofdecaymember_it)->decay_p2n<<"\t"<<(*flistofdecaymember_it)->neueff<<"\t"<<
              (*flistofdecaymember_it)->is_decay_hl_fix<<"\t"<<(*flistofdecaymember_it)->is_decay_p1n_fix<<"\t"<<(*flistofdecaymember_it)->is_decay_p2n_fix<<"\t"<<(*flistofdecaymember_it)->is_neueff_fix<<endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void decaypath::drawPath(char* outputFileName)
{
    std::clog<< __PRETTY_FUNCTION__<<std::endl;
    Double_t minz=100000;
    Double_t minn=100000;
    Double_t maxz=0;
    Double_t maxn=0;

    Double_t expandZ=3;
    Double_t expandN=3;

    for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    {
        if ((*flistofdecaymember_it)->z<minz) minz=(*flistofdecaymember_it)->z;
        if ((*flistofdecaymember_it)->n<minn) minn=(*flistofdecaymember_it)->n;

        if ((*flistofdecaymember_it)->z>maxz) maxz=(*flistofdecaymember_it)->z;
        if ((*flistofdecaymember_it)->n>maxn) maxn=(*flistofdecaymember_it)->n;
    }

    TCanvas* c1=new TCanvas("cc","",900,700) ;
    Double_t xrange[2]={minn-expandN,maxn+expandN};
    Double_t yrange[2]={minz-expandZ,maxz+expandZ};


    // c1->GetFrame()->SetFillColor(21);
    // c1->GetFrame()->SetBorderSize(12);
    // gStyle->SetOptStat(0);


    // Double_t minhalflife=0.0001;//100 ns
    // TH2F *hchart = new TH2F("hist","",185,-0.5,184.5,127,-0.5,126.5);
    // c1->SetLogz(0);
    // hchart->SetTitleSize(0.04);
    // hchart->GetXaxis()->SetTitleOffset(1.0);
    // hchart->GetYaxis()->SetTitleOffset(1.2);
    // hchart->GetYaxis()->CenterTitle();
    // hchart->GetXaxis()->SetLabelSize(0.03);
    // hchart->GetYaxis()->SetLabelSize(0.03);
    // hchart->GetYaxis()->SetTitle("N_{Proton}");
    // hchart->GetXaxis()->SetTitle("N_{Neutron}");
    // hchart->GetXaxis()->SetRangeUser(xrange[0],xrange[1]);
    // hchart->GetYaxis()->SetRangeUser(yrange[0],yrange[1]);
    // hchart->SetMinimum(minhalflife);


    // c1->SetLogz();
    // hchart->SetLineWidth(10);
    // hchart->SetLineColor(1);
    // hchart->Draw("COLZ");

    // //! Drawing magic number
    // Double_t dd = 0.5;
    // TLine a1;
    // //  a1.SetLineWidth(1.5);
    // a1.SetLineWidth(3.0);
    // a1.SetLineColor(7);

    // Int_t magicn[]={8,20,28,50,82,126};

    // for (Int_t i=0;i<6;i++){
    //     a1.DrawLine(magicn[i]-dd,yrange[0]-dd,magicn[i]-dd,yrange[1]+dd); a1.DrawLine(magicn[i]+1-dd,yrange[0]-dd,magicn[i]+1-dd,yrange[1]+dd);
    //     a1.DrawLine(xrange[0]-dd,magicn[i]-dd,xrange[1]+dd,magicn[i]-dd); a1.DrawLine(xrange[0]-dd,magicn[i]+1-dd,xrange[1]+dd,magicn[i]+1-dd);
    // }

    // //! draw paths
    // TLatex latex;
    // latex.SetTextAlign(12);
    // latex.SetTextSize(0.025);

    // TArrow arr;
    // arr.SetLineColor(6);
    // arr.SetFillColor(2);
    // //arr.SetLineStyle(2);
    // // Plot the flow
    // Int_t npaths=0;

    // for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    // {
    //     Int_t isplotiso=0;
    //     for (Size_t i=0;i<(*flistofdecaymember_it)->path.size();i++){
    //         Double_t prevz=0;
    //         Double_t prevn=0;
    //         Double_t previd=0;
    //         Double_t prevp1n=0;
    //         Double_t prevp2n=0;
    //         Bool_t isplot=true;

    //         std::list<MemberDef*>::iterator listofdecaymember_it2;
    //         for (Size_t j=0;j<(*flistofdecaymember_it)->path[i].size();j++){
    //             //! draw arrow
    //             for (listofdecaymember_it2 = flistofdecaymember.begin(); listofdecaymember_it2 != flistofdecaymember.end(); listofdecaymember_it2++)
    //             {
    //                 if ((*flistofdecaymember_it)->path[i][j]==(*listofdecaymember_it2)->id){

    //                     Double_t presz=(*listofdecaymember_it2)->z;
    //                     Double_t presn=(*listofdecaymember_it2)->n;
    //                     Double_t presid=(*listofdecaymember_it2)->id;


    //                     if (prevz+prevn==presz+presn){
    //                         if ((1-prevp1n+prevp2n)==0) {
    //                             isplot=false;
    //                         }
    //                     }else if(presz+presn==prevz+prevn-1){
    //                         if (prevp1n==0) {
    //                             isplot=false;
    //                         }
    //                     }else if((presz+presn==prevz+prevn-2)){
    //                         if (prevp2n==0) {
    //                             isplot=false;
    //                         }
    //                     }

    //                     if (prevz!=0&&isplot){
    //                         arr.DrawArrow(prevn,prevz,presn,presz,0.01,">");
    //                         //! A trick for plotting, draw at the end of each track another arrow
    //                         arr.DrawArrow(presn,presz,presn-1,presz+1,0.01,">");
    //                         isplotiso++;
    //                         //latex.DrawLatex((*listofdecaymember_it2)->n-0.5,(*listofdecaymember_it2)->z,Form("%s",(*listofdecaymember_it2)->name.Data()));
    //                     }

    //                     prevz=presz;
    //                     prevn=presn;
    //                     prevp1n=(*listofdecaymember_it2)->decay_p1n;
    //                     prevp2n=(*listofdecaymember_it2)->decay_p2n;
    //                     previd=presid;
    //                 }
    //             }

    //         }
    //         npaths++;
    //     }
    // }

    // for (flistofdecaymember_it = flistofdecaymember.begin(); flistofdecaymember_it != flistofdecaymember.end(); flistofdecaymember_it++)
    // {
    //     if ((*flistofdecaymember_it)->id==0) latex.DrawLatex((*flistofdecaymember_it)->n-0.5,(*flistofdecaymember_it)->z,Form("%s",(*flistofdecaymember_it)->name.Data()));
    //     else latex.DrawLatex((*flistofdecaymember_it)->n-0.5,(*flistofdecaymember_it)->z,Form("%s",(*flistofdecaymember_it)->name.Data()));
    // }

    //! save to file
    c1->SaveAs(outputFileName);
}

