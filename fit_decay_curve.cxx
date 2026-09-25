// pgenfit3 fit driver: unbinned, extended, SIMULTANEOUS fit of the total
// beta-decay curve (no neutron gates) in forward and backward time.
//
//   pos (x >= startTime): activity of parent + all daughters + bkg (flat)
//   neg (x <  0)        : bkg (flat)
//
// x = t_beta - t_ion from the "tree" of pgenfit's simulation.cc output (or
// real data with the same layout). The neutron multiplicity branch y is
// ignored. A beta cannot precede its own implant, so the backward region is
// a pure accidental (wrong-ion) sample; its rate bkg (counts/s) is the SAME
// parameter as the forward-time background, so it is measured and
// propagated within one likelihood instead of being fixed from a separate
// pre-fit (pgenfit's unbinfit fitOpt=2 fixed nbkg from the backward count).
//
// Usage:
//   ./fit_decay_curve <data.root> <parmsex> [ncpu=1] [effparms|none]
//                     [out_prefix|noplot] [boundaryMarginSec=0]
//                     [startTime=0] [timeRange=10] [linBinFactor=4] [nominos]
//
//   out_prefix        write <out_prefix>_fit.png and <out_prefix>_fit.root;
//                     "noplot" (or "") skips plotting
//   boundaryMarginSec drop implants within this many seconds of the first/
//                     last observed implant (needs the ionT branch); the
//                     accidental rate is not stationary near run edges
//   startTime         forward-time dead-time cut, positive side only
//   timeRange         fit window: pos [startTime,T], neg [-T,0)
//   linBinFactor      linear bin width = parent T1/2 / linBinFactor
//                     (goodness-of-fit and linear plot)
//
// Parameters, from the parmsex file (negative value = float):
//   half-lives of any species; isomer ratio (floats the ground-state
//   feeding py<gs>). P1n/P2n only shape the daughter feeding here, so their
//   float flags are IGNORED (held fixed) unless FLOAT_PN=1 is set.
//   alpha branch pa<k> (parmsex columns 24-28, AlphaBR in percent, as in
//   pgenfit; an isomer row may carry a second block for the isomer).
//   be (parent beta-efficiency factor, effparms column 1) and ab (alpha/beta
//   detection-efficiency ratio, optional effparms columns 10-11, default 1):
//   negative = float freely; positive with nonzero error =
//   Gaussian-constrained; otherwise fixed. Always floating: N0raw, bkg.
//
// Alpha gate (ALPHA_GATE=1; network must contain an alpha branch): the data
// are split, like pgenfit2's neutron gates, into mutually exclusive gates
// fitted simultaneously, each in forward and backward time:
//   beta gate  (bpos/bneg): decays NOT tagged as alpha
//   alpha gate (apos/aneg): alpha-tagged decays, from any species
// with separate accidental rates bkg / bkga. The total "beta-or-alpha minus
// implant" curve is bpos+apos (plotted and gof-checked as "sumpos"); fitting
// the total and the alpha-only curve directly would count every alpha twice.
// The tag is branch ALPHA_BRANCH (default "alpha", written by pgenfit3's
// simulation.cc), tagged if > 0. (pgenfit's own "z" branch holds the member
// id, so its parent alphas, id 0, would look untagged.)
//
// Model-generated toys (GEN_TOYS=N, GEN_SEED=s): after the fit, generate N
// datasets from the FITTED pdf itself (independent events by construction,
// unlike simulation.cc whose decays of one implant are correlated along the
// chain), refit each, and print one "GENTOY PULLDATA ..." line per toy with
// <name>true = the value the toys were generated with. A pull width of 1
// here and <1 with simulation.cc toys isolates the chain correlation.
//
// Environment: FLOAT_PN=1 (see above); ALPHA_GATE, ALPHA_BRANCH (above);
// GEN_TOYS, GEN_SEED (above);
// SCAN_PARAM/SCAN_POINTS/SCAN_MIN/SCAN_MAX for a profile-likelihood scan
// printed as SCANDATA lines.
//
// Output: a "PULLDATA key=value ..." line with every floating parameter
// (plus <name>true for species parameters, from the parmsex file) and the
// per-region Baker-Cousins goodness-of-fit.

#include <decayModel.hh>
#include <RooRealVar.h>
#include <RooArgList.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooMsgService.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooMinimizer.h>
#include <RooRandom.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooCurve.h>
#include <RooBinning.h>

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TTreeFormula.h>
#include <TParameter.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <map>
#include <memory>
#include <chrono>
#include <algorithm>

using namespace RooFit;

struct Species {
    std::string name;
    double l, llow, lup, p1n, p2n;
    bool islvary, isp1nvary, isp2nvary;
    double populationRatio=1.0;
    bool isPopulationRatioVary=false;
    double pa=0;              // alpha branching fraction
    bool ispavary=false;
};

// Same column layout as pgenfit/pgenfit2 (23 columns; an isomer row, name
// ending in '*', carries 4 isomer-ratio columns and a 20-column block for
// the isomer's own decay), then optional alpha blocks exactly as read by
// decaypath.cc. Emits ground state then isomer, matching decaypath.cc's
// member order (and so path.txt's indices).
std::vector<Species> parseParmsex(const char* path)
{
    std::vector<Species> out;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f,line)){
        if (line.empty() || line[0]=='#') continue;
        std::istringstream iss(line);
        std::string name; double z,a,hl,hlerr,hlerrhi,hllow,hlup,p1n,p1nerr,p1nerrhi,p1nlow,p1nup,
                    p2n,p2nerr,p2nerrhi,p2nlow,p2nup,ne,neerr,neerrhi,nelow,neup;
        if (!(iss>>name>>z>>a>>hl>>hlerr>>hlerrhi>>hllow>>hlup
                  >>p1n>>p1nerr>>p1nerrhi>>p1nlow>>p1nup
                  >>p2n>>p2nerr>>p2nerrhi>>p2nlow>>p2nup
                  >>ne>>neerr>>neerrhi>>nelow>>neup)) continue;
        bool isIsomer = !name.empty() && name.back()=='*';
        double isoRatio=0,isoRatioErr=0,isoRatioLow=0,isoRatioUp=0;
        double ihl=0,ihlerr=0,ihlerrhi=0,ihllow=0,ihlup=0,
               ip1n=0,ip1nerr=0,ip1nerrhi=0,ip1nlow=0,ip1nup=0,
               ip2n=0,ip2nerr=0,ip2nerrhi=0,ip2nlow=0,ip2nup=0,
               ine=0,ineerr=0,ineerrhi=0,inelow=0,ineup=0;
        if (isIsomer){
            name.pop_back();
            if (!(iss>>isoRatio>>isoRatioErr>>isoRatioLow>>isoRatioUp
                      >>ihl>>ihlerr>>ihlerrhi>>ihllow>>ihlup
                      >>ip1n>>ip1nerr>>ip1nerrhi>>ip1nlow>>ip1nup
                      >>ip2n>>ip2nerr>>ip2nerrhi>>ip2nlow>>ip2nup
                      >>ine>>ineerr>>ineerrhi>>inelow>>ineup)){
                std::cerr<<"isomer row '"<<name<<"*' is missing columns -- skipped\n";
                continue;
            }
        }
        // half-life bounds are absolute, so the decay-constant limits are
        // their exact inversion (sides swap); see pgenfit2 CAMPAIGN_LOG 1.1
        auto makeSpecies=[](const std::string& nm,double hl_,double hllow_,double hlup_,double p1n_,double p2n_){
            Species s;
            s.name=nm;
            s.islvary=(hl_<0);
            s.l=std::log(2.0)/std::fabs(hl_);
            s.isp1nvary=(p1n_<0); s.p1n=std::fabs(p1n_)/100.0;
            s.isp2nvary=(p2n_<0); s.p2n=std::fabs(p2n_)/100.0;
            s.llow=(hlup_>0) ? std::log(2.0)/hlup_ : 0.0;
            s.lup =(hllow_>0) ? std::log(2.0)/hllow_ : s.l*1e3;
            return s;
        };
        // optional "AlphaBR errLo errHi low up" block (percent)
        auto readAlpha=[&](Species& sp){
            double abr,e1,e2,lo,up;
            if (iss>>abr>>e1>>e2>>lo>>up){ sp.ispavary=(abr<0); sp.pa=std::fabs(abr)/100.0; }
        };
        Species s=makeSpecies(name,hl,hllow,hlup,p1n,p2n);
        if (isIsomer){
            s.isPopulationRatioVary=(isoRatio<0);
            s.populationRatio=1.0-std::fabs(isoRatio);
        }
        readAlpha(s);
        out.push_back(s);
        if (isIsomer){
            Species iso=makeSpecies(name+"m",ihl,ihllow,ihlup,ip1n,ip2n);
            readAlpha(iso);
            out.push_back(iso);
        }
    }
    return out;
}

struct EffFactor { double val=1, err=0; bool vary=false; };
struct Eff { EffFactor be, ab; };

// pgenfit effparms line: be err b1ne err b2ne err n1n2ne errLo errHi
// [ab err]. A decay-curve-only fit uses be (columns 1-2) and the optional
// alpha/beta efficiency ratio ab (columns 10-11); the neutron columns are
// ignored. Negative value = float.
Eff parseEff(const std::string& path)
{
    Eff e;
    if (path.empty() || path=="none") return e;
    std::ifstream f(path);
    if (!f){ std::cerr<<"cannot open effparms "<<path<<" -- using be=1, ab=1 fixed\n"; return e; }
    std::string line;
    while (std::getline(f,line)){
        if (line.empty() || line[0]=='#') continue;
        std::istringstream iss(line);
        std::vector<double> v; double x;
        while (iss>>x) v.push_back(x);
        if (v.size()<2) continue;
        e.be.val=v[0]; e.be.err=v[1];
        if (v.size()>=11){ e.ab.val=v[9]; e.ab.err=v[10]; }
        break;
    }
    for (EffFactor* ef : {&e.be,&e.ab}){ ef->vary = ef->val<0; ef->val = std::fabs(ef->val); }
    return e;
}

// Rescale a log-binned frame from RooFit's per-nominal-bin-width density to
// raw counts per log channel (Schmidt convention); same as pgenfit2.
void toSchmidtCounts(RooPlot* fr, const RooBinning& bins)
{
    double nominalW = 0.0;
    for (int i=0;i<fr->numItems();i++){
        TObject* o = fr->getObject(i);
        if (o && o->InheritsFrom(RooHist::Class())){ nominalW = ((RooHist*)o)->getNominalBinWidth(); break; }
    }
    if (!(nominalW>0)) return;
    const double ratio = (bins.numBins()>0 && bins.lowBound()>0)
        ? std::pow(bins.highBound()/bins.lowBound(), 1.0/bins.numBins()) : 1.0;
    if (!(ratio>1.0)) return;
    double ymax = 0.0;
    for (int i=0;i<fr->numItems();i++){
        TObject* o = fr->getObject(i);
        if (!o) continue;
        if (o->InheritsFrom(RooHist::Class())){
            RooHist* h=(RooHist*)o;
            for (int p=0;p<h->GetN();p++){
                double x,y; h->GetPoint(p,x,y);
                double exl=h->GetErrorXlow(p), exh=h->GetErrorXhigh(p);
                double w = exl+exh;
                if (!(w>0)) continue;
                double s = w/nominalW;
                h->SetPoint(p,x,y*s);
                h->SetPointError(p,exl,exh,h->GetErrorYlow(p)*s,h->GetErrorYhigh(p)*s);
                ymax = std::max(ymax, y*s + h->GetErrorYhigh(p)*s);
            }
        } else if (o->InheritsFrom(RooCurve::Class())){
            RooCurve* c=(RooCurve*)o;
            for (int p=0;p<c->GetN();p++){
                double x,y; c->GetPoint(p,x,y);
                if (!(x>0)) continue;
                double w = x*(std::sqrt(ratio) - 1.0/std::sqrt(ratio));
                c->SetPoint(p,x,y*w/nominalW);
                ymax = std::max(ymax, y*w/nominalW);
            }
        }
    }
    fr->GetYaxis()->SetTitle("counts / log channel");
    if (ymax>0){ fr->SetMinimum(0.0); fr->SetMaximum(ymax*1.15); }
}

RooBinning makeLogBinning(double xlo, double xhi, double rate, double binsPerDoubling)
{
    double lo = std::max(xlo, std::min(1.0/(rate*50.0), (xhi-xlo)*0.5));
    if (!(lo>0) || !(lo<xhi)) lo = (xhi-xlo)*0.02;
    int nbins = (int)std::llround(std::log(xhi/lo)/(std::log(2.0)/binsPerDoubling));
    nbins = std::max(10, std::min(nbins, 2000));
    std::vector<double> edges(nbins+1);
    double r = std::pow(xhi/lo, 1.0/nbins);
    edges[0] = lo;
    for (int i=1;i<=nbins;i++) edges[i] = edges[i-1]*r;
    edges[nbins] = xhi;
    return RooBinning(nbins, edges.data());
}

// Baker-Cousins chi2 = 2*sum[m-d+d*ln(d/m)] between a binned data hist and
// the curve's bin-averaged value (same average pullHist() uses).
double bakerCousins(RooHist* h, RooCurve* c, int& nbins)
{
    double chi2=0; nbins=0;
    for (int ib=0; ib<h->GetN(); ib++){
        double xi, di; h->GetPoint(ib, xi, di);
        double exl=h->GetErrorXlow(ib), exh=h->GetErrorXhigh(ib);
        if (exl<=0) exl=0.5*h->getNominalBinWidth();
        if (exh<=0) exh=0.5*h->getNominalBinWidth();
        double mi = c->average(xi-exl, xi+exh);
        if (!(mi>0)) continue;
        chi2 += 2.0*(mi-di);
        if (di>0) chi2 += 2.0*di*std::log(di/mi);
        nbins++;
    }
    return chi2;
}

int main(int argc, char** argv)
{
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    if (argc<3){
        std::cerr<<"Usage: "<<argv[0]<<" <data.root> <parmsex> [ncpu=1] [effparms|none] [out_prefix|noplot]"
                   " [boundaryMarginSec=0] [startTime=0] [timeRange=10] [linBinFactor=4] [nominos]\n";
        return 1;
    }
    const std::string dataFile = argv[1];
    const std::string parmsexFile = argv[2];
    const int ncpu = (argc>3) ? atoi(argv[3]) : 1;
    const std::string effFile = (argc>4) ? argv[4] : "none";
    const std::string outPrefix = (argc>5) ? argv[5] : "noplot";
    const bool doPlot = !(outPrefix.empty() || outPrefix=="noplot");
    const double boundaryMarginSec = (argc>6) ? atof(argv[6]) : 0.0;
    const double startTime = (argc>7) ? atof(argv[7]) : 0.0;
    const double timeRange = (argc>8) ? atof(argv[8]) : 10.0;
    const double linBinFactor = (argc>9) ? atof(argv[9]) : 4.0;
    const bool forceNoMinos = (argc>10) && (std::string(argv[10])=="nominos");
    const char* floatPnEnv = std::getenv("FLOAT_PN");
    const bool floatPn = floatPnEnv && std::string(floatPnEnv)=="1";
    const char* alphaGateEnv = std::getenv("ALPHA_GATE");
    const bool alphaGate = alphaGateEnv && std::string(alphaGateEnv)=="1";
    const char* alphaBranchEnv = std::getenv("ALPHA_BRANCH");
    const std::string alphaBranch = (alphaBranchEnv && *alphaBranchEnv) ? alphaBranchEnv : "alpha";

    if (!(timeRange>startTime) || startTime<0){
        std::cerr<<"need 0 <= startTime < timeRange\n";
        return 1;
    }

    Int_t nri;
    { std::ifstream f("path.txt"); if (!(f>>nri)){ std::cerr<<"cannot read path.txt (run ./main <parmsex> first)\n"; return 1; } }
    auto species = parseParmsex(parmsexFile.c_str());
    if ((int)species.size()!=nri){
        std::cerr<<"parmsex file has "<<species.size()<<" species but path.txt says nri="<<nri<<"\n";
        return 1;
    }
    std::cout<<"Loaded "<<nri<<" species. Parent="<<species[0].name
             <<" T1/2(input)="<<std::log(2.0)/species[0].l<<" s\n";

    // ---- parameters ----
    RooRealVar x_pos("x_pos","t (s)",startTime,timeRange);
    RooRealVar x_neg("x_neg","t (s)",-timeRange,0.0);
    RooArgList P;
    std::vector<RooRealVar*> l(nri),p1n(nri),p2n(nri),py(nri),pa(nri);
    std::map<std::string,double> trueVal; // parmsex values of species parameters
    for (int k=0;k<nri;k++){
        const Species& s=species[k];
        l[k]=new RooRealVar(Form("l%d",k),"",s.l,s.llow,s.lup);
        p1n[k]=new RooRealVar(Form("p1n%d",k),"",s.p1n,0,1);
        p2n[k]=new RooRealVar(Form("p2n%d",k),"",s.p2n,0,1);
        py[k]=new RooRealVar(Form("py%d",k),"",s.populationRatio,0,1);
        pa[k]=new RooRealVar(Form("pa%d",k),"",s.pa,0,1);
        pa[k]->setConstant(!s.ispavary);
        if (s.pa>0 || s.ispavary)
            std::cout<<"alpha branch: "<<s.name<<" pa="<<s.pa<<(s.ispavary?" (floating)":" (fixed)")<<"\n";
        l[k]->setConstant(!s.islvary);
        p1n[k]->setConstant(!(floatPn && s.isp1nvary));
        p2n[k]->setConstant(!(floatPn && s.isp2nvary));
        py[k]->setConstant(!s.isPopulationRatioVary);
        if (!floatPn && (s.isp1nvary || s.isp2nvary))
            std::cout<<"note: "<<s.name<<" P1n/P2n float flag ignored (set FLOAT_PN=1 to float)\n";
        trueVal[l[k]->GetName()]=s.l; trueVal[p1n[k]->GetName()]=s.p1n;
        trueVal[p2n[k]->GetName()]=s.p2n; trueVal[py[k]->GetName()]=s.populationRatio;
        trueVal[pa[k]->GetName()]=s.pa;
        P.add(*l[k]); P.add(*p1n[k]); P.add(*p2n[k]); P.add(*py[k]); P.add(*pa[k]);
    }
    RooRealVar N0raw("N0raw","parent activity at t=0 (counts/s)",1.0,1e-6,1e12);
    RooRealVar bkg("bkg","accidental rate (counts/s)",1.0,1e-6,1e12);
    RooRealVar bkga("bkga","alpha-gate accidental rate (counts/s)",1.0,1e-6,1e12);
    bkga.setConstant(!alphaGate);
    Eff eff = parseEff(effFile);
    RooRealVar be("be","parent beta-efficiency factor",eff.be.val,0,1);
    RooRealVar ab("ab","alpha/beta detection-efficiency ratio",eff.ab.val,0,100);
    RooArgSet constraintPdfs;
    // negative = free; nonzero error = Gaussian constraint; else fixed
    auto setupFactor=[&](RooRealVar& v, const EffFactor& ef){
        if (ef.vary){
            v.setConstant(false);
        } else if (ef.err>0){
            v.setConstant(false);
            auto* mean=new RooRealVar(Form("%s_mean",v.GetName()),"",ef.val);  mean->setConstant(true);
            auto* sigma=new RooRealVar(Form("%s_sigma",v.GetName()),"",ef.err); sigma->setConstant(true);
            constraintPdfs.add(*new RooGaussian(Form("%s_constraint",v.GetName()),"",v,*mean,*sigma));
        } else {
            v.setConstant(true);
        }
        std::cout<<v.GetName()<<"="<<ef.val<<(ef.vary?" (floating)":(ef.err>0?" (constrained)":" (fixed)"))<<"\n";
    };
    setupFactor(be, eff.be);
    setupFactor(ab, eff.ab);
    // effparms values serve as the truth for toy pulls when these float
    trueVal["be"]=eff.be.val; trueVal["ab"]=eff.ab.val;
    P.add(N0raw); P.add(be); P.add(ab); P.add(bkg); P.add(bkga);

    DecayModel* M = buildDecayModel(P, x_pos, x_neg);
    if (alphaGate && !M->simPdfGated){
        std::cerr<<"ALPHA_GATE=1 but no species in the parmsex file has an alpha branch\n";
        return 1;
    }
    // total mode: cat pos=0/neg=1; gated mode: gcat bpos=0/bneg=1/apos=2/aneg=3
    RooSimultaneous& simPdf = alphaGate ? *M->simPdfGated : *M->simPdf;
    RooCategory& cat = alphaGate ? *M->gcat : *M->cat;
    if (alphaGate) std::cout<<"Alpha-gated fit: beta gate + alpha gate (tag: "<<alphaBranch<<">0)\n";

    // ---- data ----
    TFile* f = TFile::Open(dataFile.c_str());
    if (!f || f->IsZombie()){ std::cerr<<"cannot open "<<dataFile<<"\n"; return 1; }
    TTree* tree = (TTree*)f->Get("tree");
    if (!tree){ std::cerr<<"no 'tree' in "<<dataFile<<"\n"; return 1; }
    // pgenfit3 simulation files record their ion-decay correlation window;
    // a fit window beyond it would silently truncate the data
    {
        auto* wlo=(TParameter<double>*)f->Get("ionbetawindowlow");
        auto* whi=(TParameter<double>*)f->Get("ionbetawindowup");
        if (wlo && whi){
            std::cout<<"Correlation window in file: -"<<wlo->GetVal()<<" s to +"<<whi->GetVal()<<" s\n";
            if (timeRange>whi->GetVal()*(1+1e-9) || timeRange>wlo->GetVal()*(1+1e-9)){
                std::cerr<<"timeRange="<<timeRange<<" s exceeds the file's correlation window (-"<<wlo->GetVal()
                         <<", +"<<whi->GetVal()<<") -- the fit regions would be truncated\n";
                return 1;
            }
        }
    }
    double simx=0, ionT=0;
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("x",1);
    tree->SetBranchAddress("x",&simx);
    const bool haveIonT = (tree->GetBranch("ionT")!=nullptr);
    if (haveIonT){ tree->SetBranchStatus("ionT",1); tree->SetBranchAddress("ionT",&ionT); }
    std::unique_ptr<TTreeFormula> alphaTag;
    if (alphaGate){
        if (!tree->GetBranch(alphaBranch.c_str())){
            std::cerr<<"ALPHA_GATE=1 but the tree has no '"<<alphaBranch<<"' branch (set ALPHA_BRANCH)\n";
            return 1;
        }
        tree->SetBranchStatus(alphaBranch.c_str(),1);
        alphaTag.reset(new TTreeFormula("alphaTag",(alphaBranch+">0").c_str(),tree));
    }
    const Long64_t nEntries = tree->GetEntries();

    double ionTmin=1e300, ionTmax=-1e300;
    const bool doBoundary = boundaryMarginSec>0 && haveIonT;
    if (boundaryMarginSec>0 && !haveIonT)
        std::cout<<"boundaryMarginSec requested but the tree has no ionT branch -- no boundary cut\n";
    if (doBoundary){
        for (Long64_t i=0;i<nEntries;i++){
            tree->GetEntry(i);
            ionTmin=std::min(ionTmin,ionT); ionTmax=std::max(ionTmax,ionT);
        }
        std::cout<<"Run-boundary cut: ionT range ["<<ionTmin<<", "<<ionTmax<<"], excluding "
                 <<boundaryMarginSec<<" s at each edge\n";
    }

    RooDataSet data("data","",RooArgSet(x_pos,x_neg,cat));
    // nPos/nNeg count all decays; nAPos/nANeg the alpha-tagged ones (gated mode)
    Long64_t nPos=0, nNeg=0, nAPos=0, nANeg=0, nOut=0, nDead=0, nBoundary=0;
    for (Long64_t i=0;i<nEntries;i++){
        tree->GetEntry(i);
        if (doBoundary && (ionT<ionTmin+boundaryMarginSec || ionT>ionTmax-boundaryMarginSec)){ nBoundary++; continue; }
        bool isAlpha=false;
        if (alphaTag){ alphaTag->GetNdata(); isAlpha = alphaTag->EvalInstance()!=0; }
        const int gateOffset = isAlpha ? 2 : 0;
        if (simx>=0.0){
            if (simx<startTime){ nDead++; continue; }
            if (simx>timeRange){ nOut++; continue; }
            x_pos.setVal(simx); cat.setIndex(0+gateOffset);
            data.add(RooArgSet(x_pos,x_neg,cat)); nPos++; if (isAlpha) nAPos++;
        } else {
            if (simx<-timeRange){ nOut++; continue; }
            x_neg.setVal(simx); cat.setIndex(1+gateOffset);
            data.add(RooArgSet(x_pos,x_neg,cat)); nNeg++; if (isAlpha) nANeg++;
        }
    }
    std::cout<<"Loaded "<<nPos<<" forward-time + "<<nNeg<<" backward-time events (outside window "<<nOut
             <<", dead time "<<nDead<<", run boundary "<<nBoundary<<")\n";
    if (alphaGate) std::cout<<"  of which alpha-tagged: "<<nAPos<<" forward, "<<nANeg<<" backward\n";
    if (nPos==0 || nNeg==0){ std::cerr<<"empty forward or backward region -- cannot fit\n"; return 1; }

    // ---- data-driven starting values; each bound scaled off its own guess
    // (one shared huge bound breaks Minuit's conditioning, pgenfit2 Zr108) ----
    {
        const double posDur = timeRange-startTime;
        const double bkgGuess = std::max(1e-3, (double)(nNeg-nANeg)/timeRange);
        const double bkgaGuess = std::max(1e-3, (double)nANeg/timeRange);
        const double nSig = std::max(10.0, nPos - (bkgGuess+(alphaGate?bkgaGuess:0.0))*posDur);
        const double l0 = l[0]->getVal();
        // parent + daughters roughly double the parent-only count
        const double n0Guess = 0.5*nSig*l0*std::exp(l0*startTime)/std::max(1e-12,1.0-std::exp(-l0*posDur))/std::max(1e-6,eff.be.val);
        N0raw.setRange(1e-6, std::max(100.0, 50.0*n0Guess)); N0raw.setVal(n0Guess);
        bkg.setRange(1e-6, std::max(100.0, 50.0*bkgGuess));  bkg.setVal(bkgGuess);
        bkga.setRange(1e-6, std::max(100.0, 50.0*bkgaGuess)); bkga.setVal(bkgaGuess);
        std::cout<<"Initial N0raw="<<n0Guess<<" bkg="<<bkgGuess<<" /s"
                 <<(alphaGate ? Form(" bkga=%g /s",bkgaGuess) : "")<<"\n";
    }

    // ---- fit ----
    const Long64_t nTot = nPos+nNeg;
    const bool isLargeN = nTot >= 1000000;
    const bool useMinos = !forceNoMinos && !isLargeN;
    if (isLargeN) std::cout<<nTot<<" events: Minos skipped, Hesse errors only\n";
    std::cout<<std::flush; // flush before NumCPU forks
    RooFitResult* r=nullptr;
    const auto t0 = std::chrono::steady_clock::now();
    {
        std::unique_ptr<RooAbsReal> nll;
        if (ncpu>1) nll.reset(simPdf.createNLL(data, Extended(), EvalBackend("legacy"), NumCPU(ncpu), ExternalConstraints(constraintPdfs)));
        else        nll.reset(simPdf.createNLL(data, Extended(), EvalBackend("legacy"), ExternalConstraints(constraintPdfs)));
        RooMinimizer m(*nll);
        m.setPrintLevel(1);
        m.setOffsetting(true);      // FCN ~1e8 at large N; see pgenfit2
        m.setPrintEvalErrors(-1);   // count, don't store, eval errors
        if (isLargeN) m.setEps(15000);
        m.migrad();
        m.migrad();                 // restart from the first minimum
        m.hesse();
        if (useMinos){
            // N0raw (and bkga, which can sit near 0) excluded: Minos stepping
            // a parameter onto its 1e-6 bound throws
            RooArgSet minosPars;
            std::unique_ptr<RooArgSet> pars(simPdf.getParameters(data));
            for (auto* a : *pars){
                auto* v=dynamic_cast<RooRealVar*>(a);
                const std::string vn = v ? v->GetName() : "";
                if (v && !v->isConstant() && vn!="N0raw" && vn!="bkga") minosPars.add(*v);
            }
            if (!minosPars.empty()) m.minos(minosPars);
        }
        r = m.save();

        const char* scanParam = std::getenv("SCAN_PARAM");
        if (scanParam && *scanParam){
            std::unique_ptr<RooArgSet> pars(simPdf.getParameters(data));
            auto* sv = dynamic_cast<RooRealVar*>(pars->find(scanParam));
            if (!sv) std::cout<<"SCAN: no such parameter '"<<scanParam<<"'\n";
            else {
                const double best=sv->getVal();
                const double err=(sv->getError()>0) ? sv->getError() : std::fabs(best)*0.5+1e-9;
                const char* eN=std::getenv("SCAN_POINTS"); const int npts = eN ? std::max(3,atoi(eN)) : 41;
                const char* eLo=std::getenv("SCAN_MIN"); const char* eHi=std::getenv("SCAN_MAX");
                double xlo = std::max(eLo ? atof(eLo) : best-5*err, sv->getMin());
                double xhi = std::min(eHi ? atof(eHi) : best+5*err, sv->getMax());
                if (best>0 && xlo<=0) xlo=best*1e-3;
                const double nll0=nll->getVal();
                const bool wasConst=sv->isConstant();
                for (int i=0;i<npts;i++){
                    const double xv = xlo+(xhi-xlo)*i/(npts-1);
                    sv->setVal(xv); sv->setConstant(true);
                    RooMinimizer ms(*nll);
                    ms.setPrintLevel(-1); ms.setPrintEvalErrors(-1); ms.setOffsetting(true);
                    ms.migrad();
                    std::unique_ptr<RooFitResult> rs(ms.save());
                    const double v=nll->getVal();
                    std::cout<<"SCANDATA param="<<scanParam<<" x="<<std::setprecision(12)<<xv
                             <<" nll="<<v<<" dnll="<<(v-nll0)<<" status="<<(rs?rs->status():-1)<<"\n"<<std::flush;
                }
                sv->setConstant(wasConst); sv->setVal(best);
            }
        }
    }
    const double fitTimeSec = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    r->Print();

    RooRealVar& l0=*l[0];
    const double hl=std::log(2.0)/l0.getVal();
    std::cout<<"\n===== Result =====\n";
    std::cout<<"parent T1/2 = "<<hl<<" s  (+"<<hl-std::log(2.0)/(l0.getVal()+l0.getErrorHi())
             <<" / -"<<std::log(2.0)/(l0.getVal()+l0.getErrorLo())-hl<<")   input "<<std::log(2.0)/species[0].l<<" s\n";
    std::cout<<"bkg rate    = "<<bkg.getVal()<<" +/- "<<bkg.getError()<<" /s\n";
    std::cout<<"status="<<r->status()<<" covQual="<<r->covQual()<<"  fit time "<<fitTimeSec<<" s\n";
    std::cout<<"\n----- total decays per physically-decaying species (all time) -----\n";
    for (int k=0;k<nri;k++)
        std::cout<<"  "<<species[k].name<<": "<<M->TrueNbetaPerSpecies[k]->getVal()
                 <<" +/- "<<M->TrueNbetaPerSpecies[k]->getPropagatedError(*r)<<"\n";

    // ---- goodness of fit (and plots) ----
    // Forward time: linear bins of width parent T1/2/linBinFactor (fine log
    // bins near t=0 hold <<1 expected count and are poor for a chi2 sum).
    // One entry per fitted region; in gated mode also "sumpos" = bpos+apos,
    // the total beta-or-alpha curve (a derived check, like pgenfit2's sum gate).
    struct Region {
        std::string key, title;
        bool isPos;
        RooAbsPdf* pdf;
        std::unique_ptr<RooDataSet> data;
    };
    auto slice=[&](const std::string& cut){ return std::unique_ptr<RooDataSet>((RooDataSet*)data.reduce(Cut(cut.c_str()))); };
    std::unique_ptr<RooAddPdf> sumGatedPos;
    std::vector<Region> regions;
    if (!alphaGate){
        regions.push_back({"pos","forward time",true,M->decayPos,slice("cat==cat::pos")});
        regions.push_back({"neg","backward time (accidental)",false,M->bkgNeg,slice("cat==cat::neg")});
    } else {
        sumGatedPos.reset(new RooAddPdf("sumGatedPos","",RooArgList(*M->betaPos,*M->alphaPos)));
        regions.push_back({"sumpos","forward time, beta-or-alpha",true,sumGatedPos.get(),slice("gcat==gcat::bpos || gcat==gcat::apos")});
        regions.push_back({"bpos","forward time, beta gate (not alpha-tagged)",true,M->betaPos,slice("gcat==gcat::bpos")});
        regions.push_back({"apos","forward time, alpha gate",true,M->alphaPos,slice("gcat==gcat::apos")});
        regions.push_back({"bneg","backward time, beta gate",false,M->bkgNeg,slice("gcat==gcat::bneg")});
        regions.push_back({"aneg","backward time, alpha gate",false,M->alphaNeg,slice("gcat==gcat::aneg")});
    }
    const int nLin = std::max(1,(int)std::llround((timeRange-startTime)/(hl/linBinFactor)));
    RooBinning linBins(std::min(nLin,5000), startTime, timeRange);
    const int nNegBins = 100;

    auto countFloating=[&](RooAbsPdf* pdf, RooDataSet* d){
        std::unique_ptr<RooArgSet> ps(pdf->getParameters(*d));
        int n=0;
        for (auto* a : *ps){ auto* v=dynamic_cast<RooRealVar*>(a); if (v && !v->isConstant()) n++; }
        return n;
    };
    std::map<std::string,double> gofChi2, gofNbins, gofNdof;
    for (auto& rg : regions){
        if (rg.data->numEntries()==0) continue;
        RooRealVar& obs = rg.isPos ? x_pos : x_neg;
        std::unique_ptr<RooPlot> fr(obs.frame());
        if (rg.isPos) rg.data->plotOn(fr.get(), Binning(linBins), Name("d"));
        else          rg.data->plotOn(fr.get(), Binning(nNegBins), Name("d"));
        rg.pdf->plotOn(fr.get(), Precision(1e-8), Name("c"));
        int nb; double c2=bakerCousins(fr->getHist("d"), fr->getCurve("c"), nb);
        gofChi2[rg.key]=c2; gofNbins[rg.key]=nb; gofNdof[rg.key]=nb-countFloating(rg.pdf,rg.data.get());
    }
    for (auto& kv : gofChi2)
        std::cout<<"GOF "<<kv.first<<": chi2(Baker-Cousins)/ndof = "<<kv.second<<"/"<<gofNdof[kv.first]
                 <<" = "<<kv.second/std::max(1.0,gofNdof[kv.first])<<"\n";

    if (doPlot){
        // One row per fitted region: linear | pull | log-x (counts per log
        // channel) | pull. Forward-time rows first (total, then the gates),
        // then one row of backward-time panels (flat, linear only).
        // components of the total forward-time curve
        RooAbsPdf* bkgPosAll = alphaGate ? M->CombinedBkgPosGated : M->CombinedBkgPos;
        RooAddPdf bkgPlusParent("bkgPlusParent","",RooArgList(*bkgPosAll,*M->ActivityPdfParentPos));
        std::vector<const Region*> posRegions, negRegions;
        for (auto& rg : regions) (rg.isPos ? posRegions : negRegions).push_back(&rg);
        const int nrows = posRegions.size() + 1;
        TCanvas c("cfit","pgenfit3 decay-curve fit",2400,380*nrows);
        c.Divide(4,nrows);
        auto drawPos=[&](int pad, const RooBinning& bins, bool logBins, const Region& rg, bool components){
            std::string title = rg.title + (logBins?" (log x)":" (linear x)") + (components?": total, bkg+parent, bkg":"");
            RooPlot* fr = x_pos.frame(Range(bins.lowBound(),bins.highBound()), Title(title.c_str()));
            rg.data->plotOn(fr, Binning(bins), Name("d"));
            rg.pdf->plotOn(fr, LineColor(kRed), Precision(1e-8), Name("c"));
            if (components){
                bkgPlusParent.plotOn(fr, LineColor(kBlue), LineStyle(kDashed), Precision(1e-8));
                bkgPosAll->plotOn(fr, LineColor(kGreen+2), LineStyle(kDotted), Precision(1e-8));
            }
            RooPlot* pf = x_pos.frame(Range(bins.lowBound(),bins.highBound()), Title("pull"));
            pf->addPlotable(fr->pullHist("d","c"),"P");
            if (logBins) toSchmidtCounts(fr, bins);
            c.cd(pad);   if (logBins) gPad->SetLogx(); else gPad->SetLogy(); fr->Draw();
            c.cd(pad+1); if (logBins) gPad->SetLogx(); pf->Draw();
        };
        auto drawNeg=[&](int pad, const Region& rg){
            RooPlot* fn = x_neg.frame(Title(rg.title.c_str()));
            rg.data->plotOn(fn, Binning(nNegBins), Name("d"));
            rg.pdf->plotOn(fn, LineColor(kRed), Precision(1e-8), Name("c"));
            RooPlot* pn = x_neg.frame(Title("pull"));
            pn->addPlotable(fn->pullHist("d","c"),"P");
            c.cd(pad); fn->Draw();
            c.cd(pad+1); pn->Draw();
        };
        // linear x: the total curve shows the first 20 parent half-lives
        // (bin width T1/2/linBinFactor); the gates, whose shape is set by the
        // daughters, show the whole window in 100 bins. The gof above always
        // uses the whole window.
        const double linHi = std::min(timeRange, startTime+20.0*hl);
        RooBinning linTotalBins(std::max(1,(int)std::llround((linHi-startTime)/(hl/linBinFactor))), startTime, linHi);
        RooBinning linGateBins(100, startTime, timeRange);
        RooBinning logBins = makeLogBinning(startTime>0?startTime:0.0, timeRange, l0.getVal(), 4.0);
        for (size_t i=0;i<posRegions.size();i++){
            const bool isTotal = (i==0);
            drawPos(4*i+1, isTotal ? linTotalBins : linGateBins, false, *posRegions[i], isTotal);
            drawPos(4*i+3, logBins, true, *posRegions[i], isTotal);
        }
        for (size_t i=0;i<negRegions.size() && i<2;i++) drawNeg(4*posRegions.size()+2*i+1, *negRegions[i]);
        c.SaveAs((outPrefix+"_fit.png").c_str());
        TFile fo((outPrefix+"_fit.root").c_str(),"RECREATE");
        c.Write();
        fo.Close();
    }

    // ---- PULLDATA ----
    const std::streamsize prec = std::cout.precision(12);
    std::cout<<"PULLDATA l0="<<l0.getVal()<<" l0err="<<l0.getError()
             <<" l0errhi="<<l0.getErrorHi()<<" l0errlo="<<l0.getErrorLo()<<" l0true="<<species[0].l
             <<" halflife="<<hl
             <<" N0raw="<<N0raw.getVal()<<" N0rawerr="<<N0raw.getError()
             <<" TrueNbeta="<<M->TrueNbetaPerSpecies[0]->getVal()
             <<" TrueNbetaerr="<<M->TrueNbetaPerSpecies[0]->getPropagatedError(*r)
             <<" npos="<<nPos<<" nneg="<<nNeg<<" napos="<<nAPos<<" naneg="<<nANeg
             <<" status="<<r->status()<<" covQual="<<r->covQual()<<" fittimesec="<<fitTimeSec;
    for (auto* a : r->floatParsFinal()){
        auto* v=(RooRealVar*)a;
        std::string vn=v->GetName();
        if (vn=="l0" || vn=="N0raw") continue;
        std::cout<<" "<<vn<<"="<<v->getVal()<<" "<<vn<<"err="<<v->getError()
                 <<" "<<vn<<"errhi="<<v->getErrorHi()<<" "<<vn<<"errlo="<<v->getErrorLo();
        if (trueVal.count(vn)) std::cout<<" "<<vn<<"true="<<trueVal[vn];
    }
    for (auto& kv : gofChi2)
        std::cout<<" gof_"<<kv.first<<"_chi2lambda="<<kv.second<<" gof_"<<kv.first<<"_nbins="<<gofNbins[kv.first]
                 <<" gof_"<<kv.first<<"_ndof="<<gofNdof[kv.first];
    std::cout<<"\n"<<std::flush;
    std::cout.precision(prec);

    // ---- model-generated toys (GEN_TOYS) ----
    const char* genEnv = std::getenv("GEN_TOYS");
    const int nGenToys = genEnv ? atoi(genEnv) : 0;
    if (nGenToys>0){
        const char* seedEnv = std::getenv("GEN_SEED");
        RooRandom::randomGenerator()->SetSeed(seedEnv ? atoi(seedEnv) : 4357);
        std::unique_ptr<RooArgSet> pars(simPdf.getParameters(data));
        std::unique_ptr<RooArgSet> genVals((RooArgSet*)pars->snapshot());
        std::cout<<"Generating "<<nGenToys<<" toys from the fitted pdf\n";
        for (int it=0; it<nGenToys; it++){
            pars->assign(*genVals);
            std::unique_ptr<RooDataSet> toy(simPdf.generate(RooArgSet(x_pos,x_neg,cat), Extended()));
            std::unique_ptr<RooAbsReal> tnll(simPdf.createNLL(*toy, Extended(), EvalBackend("legacy"), ExternalConstraints(constraintPdfs)));
            RooMinimizer tm(*tnll);
            tm.setPrintLevel(-1); tm.setOffsetting(true); tm.setPrintEvalErrors(-1);
            tm.migrad(); tm.migrad(); tm.hesse();
            RooArgSet minosPars;
            for (auto* a : *pars){
                auto* v=dynamic_cast<RooRealVar*>(a);
                const std::string vn = v ? v->GetName() : "";
                if (v && !v->isConstant() && vn!="N0raw" && vn!="bkga") minosPars.add(*v);
            }
            if (!minosPars.empty()) tm.minos(minosPars);
            std::unique_ptr<RooFitResult> tr(tm.save());
            std::cout.precision(12);
            std::cout<<"GENTOY PULLDATA l0="<<l0.getVal()<<" l0err="<<l0.getError()
                     <<" l0errhi="<<l0.getErrorHi()<<" l0errlo="<<l0.getErrorLo()
                     <<" l0true="<<((RooRealVar*)genVals->find("l0"))->getVal()
                     <<" status="<<tr->status()<<" covQual="<<tr->covQual()<<" nevents="<<toy->numEntries();
            for (auto* a : tr->floatParsFinal()){
                auto* v=(RooRealVar*)a;
                std::string vn=v->GetName();
                if (vn=="l0") continue;
                std::cout<<" "<<vn<<"="<<v->getVal()<<" "<<vn<<"err="<<v->getError()
                         <<" "<<vn<<"errhi="<<v->getErrorHi()<<" "<<vn<<"errlo="<<v->getErrorLo()
                         <<" "<<vn<<"true="<<((RooRealVar*)genVals->find(vn.c_str()))->getVal();
            }
            std::cout<<"\n"<<std::flush;
            std::cout.precision(prec);
        }
        pars->assign(*genVals);
    }
    return 0;
}
