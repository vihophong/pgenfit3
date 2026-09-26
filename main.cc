// pgenfit3 model generator.
//
// Reads a parmsex file (same format as pgenfit/pgenfit2, isomer rows
// ending in '*' included), builds the decay network with decaypath, writes
// path.txt, and emits decayModel_cal.cc: a RooFit model of the TOTAL
// beta-decay curve (no neutron gates) fitted simultaneously in forward
// and backward time -- see decayModel.hh.
//
// Usage: ./main <parmsex_file> [prefix]
//
// prefix (e.g. "s0_", for fitting a mixture of several implanted species):
// every parameter and RooFit object name in the generated code gets this
// prefix -- except the shared be, ab, bkg, bkga and the category labels --
// and the builder is named buildDecayModel_<prefix>, so several networks
// can be compiled into one library and fitted together.

#include <decaypath.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>
#include <regex>
#include <set>
using namespace std;

path* readPathFile()
{
    path* fpath=new path;
    std::ifstream pathfile("path.txt");
    pathfile>>fpath->nri;
    pathfile>>fpath->npaths;
    for (int i=0;i<fpath->npaths;i++){
        pathfile>>fpath->ndecay[i];
        pathfile>>fpath->ispathhasflow[i];
        for (int j=0;j<fpath->ndecay[i];j++){
            pathfile>>fpath->decaymap[i][j]>>fpath->nneu[i][j];
        }
    }
    pathfile>>fpath->nisomers;
    for (int i=0;i<fpath->nisomers;i++){
        pathfile>>fpath->isomer_gs_index[i]>>fpath->isomer_ex_index[i];
    }
    pathfile.close();
    return fpath;
}

//! Bateman decomposition of the total activity.
//!
//! For a path 0 -> m1 -> ... -> last, the activity of "last" is
//!   N0raw * Ck * sum_i exp(-l_mi t) / prod_{j!=i} (l_mj - l_mi)
//! with Ck = l_last * prod_steps (py_next * branching_fraction * l_prev),
//! the parent's own lambda dropped since it cancels against N0=N0raw/l0.
//! The parent's own activity is N0raw*be*exp(-l0 t). Terms are grouped by
//! exponential (coeff[m], for the fit pdf) and, separately, by the
//! physically decaying species (actShape/actCoeff/coeffRaw, for the
//! component plots and TrueNbeta counts).
//!
//! Alpha decay (species with hasAlpha[k]; pgenfit's alpha extension):
//! pa<k> is the alpha branching fraction of species k. A beta link from k
//! carries (1-pa<k>)*(p0n|p1n|p2n), an alpha link (code 100) carries pa<k>.
//! A decay of species k is detected with weight (1-pa<k>) + pa<k>*ab
//! (parent: (1-pa0)*be + pa0*ab), ab = alpha/beta detection-efficiency
//! ratio. Species without alpha get no extra factor, so beta-only networks
//! generate exactly the pgenfit2-derived formulas.
void generateDecayCurveModel(path* fpath, const std::vector<bool>& hasAlpha)
{
    std::ofstream ofnc("decayModel_cal.cc");
    ofnc<<"#include <decayModel.hh>"<<std::endl;
    ofnc<<"#include <RooExponential.h>"<<std::endl;
    ofnc<<"#include <RooFormulaVar.h>"<<std::endl;
    ofnc<<"#include <RooPolynomial.h>"<<std::endl;
    ofnc<<"#include <RooRealSumPdf.h>"<<std::endl;
    ofnc<<"#include <RooAddPdf.h>"<<std::endl;
    ofnc<<"#include <RooCategory.h>"<<std::endl;
    ofnc<<"#include <RooSimultaneous.h>"<<std::endl;
    ofnc<<"#include <RooAddition.h>"<<std::endl;
    ofnc<<std::endl;

    const Int_t nri=fpath->nri;

    //! an isomer's feeding fraction is 1-py<gs> (ISOMER_SUM_UNITY)
    std::vector<std::string> pyName(nri);
    for (Int_t k=0;k<nri;k++) pyName[k]="py"+std::to_string(k);
#ifdef ISOMER_SUM_UNITY
    for (Int_t i=0;i<fpath->nisomers;i++){
        pyName[fpath->isomer_ex_index[i]]="pyEff"+std::to_string(fpath->isomer_ex_index[i]);
    }
#endif

    std::vector<std::vector<std::string>> coeff(nri), coeffRaw(nri), actShape(nri), actCoeff(nri);
    //! alpha-gated model: the same Bateman terms split into the beta gate
    //! (decay not tagged as alpha) and the alpha gate (alpha-tagged)
    std::vector<std::vector<std::string>> coeffB(nri), coeffA(nri);
    bool anyAlpha=false;
    for (Int_t k=0;k<nri;k++) anyAlpha = anyAlpha || hasAlpha[k];

    //! detection weight of a decay of species k
    auto weight=[&](Int_t k)->std::string{
        std::string pa="pa"+std::to_string(k);
        if (k==0) return hasAlpha[0] ? "((1-"+pa+")*be+"+pa+"*ab)" : "be";
        return hasAlpha[k] ? "((1-"+pa+")+"+pa+"*ab)" : "";
    };

    //! detection weights of a decay of species k in the beta / alpha gate
    auto weightB=[&](Int_t k)->std::string{
        std::string pa="pa"+std::to_string(k);
        if (k==0) return hasAlpha[0] ? "(1-"+pa+")*be" : "be";
        return hasAlpha[k] ? "(1-"+pa+")" : "";
    };
    auto weightA=[&](Int_t k)->std::string{
        return hasAlpha[k] ? "pa"+std::to_string(k)+"*ab" : "";
    };
    auto times=[](const std::string& a, const std::string& w){ return w.empty() ? a : a+"*"+w; };

    //! parent
    coeffB[0].push_back(times("N0raw",weightB(0)));
    if (hasAlpha[0]) coeffA[0].push_back(times("N0raw",weightA(0)));
    coeff[0].push_back("N0raw*"+weight(0));
    coeffRaw[0].push_back("N0raw*be/l0");
    actShape[0].push_back("exp0");
    actCoeff[0].push_back("N0raw*"+weight(0));

    //! daughters, one Bateman chain per flowing path
    for (Int_t k=0;k<fpath->npaths;k++){
#ifdef PATHFLOW
        if (!fpath->ispathhasflow[k]) continue;
#endif
        Int_t last=fpath->decaymap[k][fpath->ndecay[k]-1];
        std::ostringstream Ck;
        Ck<<"l"<<last;
        for (int i=0;i<fpath->ndecay[k]-1;i++){
            Int_t mi=fpath->decaymap[k][i];
            Int_t mip1=fpath->decaymap[k][i+1];
            const std::string si=std::to_string(mi);
            const std::string beta = hasAlpha[mi] ? "(1-pa"+si+")*" : "";
            std::string frac;
            if (fpath->nneu[k][i]==0) frac=beta+"(1-p1n"+si+"-p2n"+si+")";
            else if (fpath->nneu[k][i]==1) frac=beta+"p1n"+si;
            else if (fpath->nneu[k][i]==2) frac=beta+"p2n"+si;
            else frac="pa"+si; // alpha link
            Ck<<"*"<<pyName[mip1]<<"*"<<frac;
            if (i>0) Ck<<"*l"<<mi;
        }
        for (int i=0;i<fpath->ndecay[k];i++){
            Int_t mi=fpath->decaymap[k][i];
            std::ostringstream denom;
            for (int j=0;j<fpath->ndecay[k];j++){
                if (j!=i) denom<<"(l"<<fpath->decaymap[k][j]<<"-l"<<mi<<")*";
            }
            denom<<"1";
            std::string term="N0raw*"+Ck.str()+"/("+denom.str()+")";
            std::string w=weight(last);
            std::string detected = w.empty() ? term : term+"*"+w;
            coeff[mi].push_back(detected);
            coeffB[mi].push_back(times(term,weightB(last)));
            if (hasAlpha[last]) coeffA[mi].push_back(times(term,weightA(last)));
            coeffRaw[last].push_back("("+term+")/l"+std::to_string(mi));
            actShape[last].push_back("exp"+std::to_string(mi));
            actCoeff[last].push_back(detected);
        }
    }

    ofnc<<"DecayModel* buildDecayModel(RooArgList& P, RooRealVar& x_pos, RooRealVar& x_neg)"<<std::endl;
    ofnc<<"{"<<std::endl;
    ofnc<<"DecayModel* M = new DecayModel();"<<std::endl;

    //! Each RooFormulaVar gets only the parameters its formula references
    //! (a RooArgList with the full parameter set breaks TFormula on large
    //! networks -- see pgenfit2/main.cc). Word-boundary match so "l2" does
    //! not match inside "l20".
    std::vector<std::string> depNamesP={"N0raw","be","ab"};
    for (Int_t k=0;k<nri;k++){
        depNamesP.push_back("l"+std::to_string(k));
        depNamesP.push_back("p1n"+std::to_string(k));
        depNamesP.push_back("p2n"+std::to_string(k));
        depNamesP.push_back("py"+std::to_string(k));
        depNamesP.push_back("pa"+std::to_string(k));
    }
    std::vector<std::string> depNamesLocal;
#ifdef ISOMER_SUM_UNITY
    for (Int_t i=0;i<fpath->nisomers;i++) depNamesLocal.push_back("pyEff"+std::to_string(fpath->isomer_ex_index[i]));
#endif
    auto isWordChar=[](char c){ return std::isalnum((unsigned char)c) || c=='_'; };
    auto findDeps=[&](const std::vector<std::string>& candidates, const std::string& s){
        std::vector<std::string> found;
        for (auto& name : candidates){
            size_t pos=0;
            while ((pos=s.find(name,pos))!=std::string::npos){
                bool leftOk = (pos==0) || !isWordChar(s[pos-1]);
                size_t endPos = pos+name.size();
                bool rightOk = (endPos>=s.size()) || !isWordChar(s[endPos]);
                if (leftOk && rightOk){ found.push_back(name); break; }
                pos+=1;
            }
        }
        return found;
    };
    auto emitDeps=[&](const std::string& varName, const std::string& formula){
        ofnc<<"RooArgList "<<varName<<"deps;"<<std::endl;
        for (auto& d : findDeps(depNamesP, formula))
            ofnc<<varName<<"deps.add(*(RooAbsArg*)P.find(\""<<d<<"\"));"<<std::endl;
        for (auto& d : findDeps(depNamesLocal, formula))
            ofnc<<varName<<"deps.add(*"<<d<<");"<<std::endl;
        return varName+"deps";
    };
    auto emitFormula=[&](const std::string& vn, const std::string& formula){
        std::string deps=emitDeps(vn,formula);
        ofnc<<"RooFormulaVar* "<<vn<<" = new RooFormulaVar(\""<<vn<<"\",\""<<formula<<"\","<<deps<<");"<<std::endl;
    };
    //! sum of small terms via RooAddition (one long formula string crashes TFormula)
    auto emitSum=[&](const std::string& vn, const std::vector<std::string>& terms){
        if (terms.empty()){ emitFormula(vn,"0"); return; }
        if (terms.size()==1){ emitFormula(vn,terms[0]); return; }
        std::ostringstream termList;
        for (size_t t=0;t<terms.size();t++){
            std::string tn=vn+"_t"+std::to_string(t);
            emitFormula(tn,terms[t]);
            termList<<(t?",":"")<<"*"<<tn;
        }
        ofnc<<"RooAddition* "<<vn<<" = new RooAddition(\""<<vn<<"\",\"\",RooArgList("<<termList.str()<<"));"<<std::endl;
    };

#ifdef ISOMER_SUM_UNITY
    for (Int_t i=0;i<fpath->nisomers;i++){
        Int_t ex=fpath->isomer_ex_index[i], gs=fpath->isomer_gs_index[i];
        emitFormula("pyEff"+std::to_string(ex),"1-py"+std::to_string(gs));
    }
#endif

    for (Int_t m=0;m<nri;m++){
        emitSum("C"+std::to_string(m),coeff[m]);
        ofnc<<"RooExponential* exp"<<m<<" = new RooExponential(\"exp"<<m<<"\",\"\",x_pos,(RooRealVar&)*P.find(\"l"<<m<<"\"),true);"<<std::endl;
    }
    for (Int_t m=0;m<nri;m++){
        std::ostringstream shapeArgs, coeffArgs;
        if (actCoeff[m].empty()){
            ofnc<<"RooRealVar* ActCoeff"<<m<<"_0 = new RooRealVar(\"ActCoeff"<<m<<"_0\",\"\",0.0); ActCoeff"<<m<<"_0->setConstant(true);"<<std::endl;
            shapeArgs<<"*exp"<<m;
            coeffArgs<<"*ActCoeff"<<m<<"_0";
        }else{
            for (size_t i=0;i<actCoeff[m].size();i++){
                std::string cn="ActCoeff"+std::to_string(m)+"_"+std::to_string(i);
                emitFormula(cn,actCoeff[m][i]);
                shapeArgs<<(i?",":"")<<"*"<<actShape[m][i];
                coeffArgs<<(i?",":"")<<"*"<<cn;
            }
        }
        ofnc<<"RooRealSumPdf* ActivityPdf"<<m<<" = new RooRealSumPdf(\"ActivityPdf"<<m<<"\",\"\",RooArgList("<<shapeArgs.str()<<"),RooArgList("<<coeffArgs.str()<<"),true);"<<std::endl;
        emitSum("NrawTotal"+std::to_string(m),coeffRaw[m]);
        ofnc<<"M->TrueNbetaPerSpecies.push_back(NrawTotal"<<m<<");"<<std::endl;
    }
    ofnc<<std::endl;

    ofnc<<"RooPolynomial* flatPos = new RooPolynomial(\"flatPos\",\"\",x_pos);"<<std::endl;
    ofnc<<"RooPolynomial* flatNeg = new RooPolynomial(\"flatNeg\",\"\",x_neg);"<<std::endl;
    ofnc<<"RooRealVar& bkg = (RooRealVar&)*P.find(\"bkg\");"<<std::endl;

    std::ostringstream shapes, coefs, totals;
    for (Int_t m=0;m<nri;m++){
        shapes<<"*exp"<<m<<",";
        coefs<<"*C"<<m<<",";
        totals<<(m?",":"")<<"*NrawTotal"<<m;
    }
    ofnc<<"RooRealSumPdf* decayPos = new RooRealSumPdf(\"decayPos\",\"\",RooArgList("<<shapes.str()<<"*flatPos),RooArgList("<<coefs.str()<<"bkg),true);"<<std::endl;
    ofnc<<"RooRealSumPdf* bkgNeg = new RooRealSumPdf(\"bkgNeg\",\"\",RooArgList(*flatNeg),RooArgList(bkg),true);"<<std::endl;
    //! the same forward-time activity without the background term (mixture fits)
    std::string shapesNoFlat=shapes.str(); shapesNoFlat.pop_back();
    std::string coefsNoBkg=coefs.str(); coefsNoBkg.pop_back();
    ofnc<<"RooRealSumPdf* signalPos = new RooRealSumPdf(\"signalPos\",\"\",RooArgList("<<shapesNoFlat<<"),RooArgList("<<coefsNoBkg<<"),true);"<<std::endl;
    ofnc<<"M->signalPos=signalPos;"<<std::endl;
    ofnc<<"RooAddition* TrueNbetaAll = new RooAddition(\"TrueNbetaAll\",\"\",RooArgList("<<totals.str()<<"));"<<std::endl;
    ofnc<<"RooRealSumPdf* CombinedBkgPos = new RooRealSumPdf(\"CombinedBkgPos\",\"\",RooArgList(*flatPos),RooArgList(bkg),true);"<<std::endl;
    if (nri>1){
        ofnc<<"RooAddPdf* ActivityPdfDaughtersPos = new RooAddPdf(\"ActivityPdfDaughtersPos\",\"\",RooArgList(";
        for (Int_t m=1;m<nri;m++) ofnc<<(m>1?",":"")<<"*ActivityPdf"<<m;
        ofnc<<"));"<<std::endl;
    }else{
        ofnc<<"RooRealVar* ActivityDaughtersZero = new RooRealVar(\"ActivityDaughtersZero\",\"\",0.0); ActivityDaughtersZero->setConstant(true);"<<std::endl;
        ofnc<<"RooRealSumPdf* ActivityPdfDaughtersPos = new RooRealSumPdf(\"ActivityPdfDaughtersPos\",\"\",RooArgList(*exp0),RooArgList(*ActivityDaughtersZero),true);"<<std::endl;
    }
    //! alpha-gated model (only for networks with an alpha branch):
    //!   bpos/bneg: decays not tagged as alpha, forward/backward time
    //!   apos/aneg: alpha-tagged decays,        forward/backward time
    //! with separate accidental rates bkg (beta gate) and bkga (alpha gate).
    //! bpos+apos is the total ("beta-or-alpha") curve, apos the alpha-only one.
    if (anyAlpha){
        std::ostringstream cB, cA;
        for (Int_t m=0;m<nri;m++){
            emitSum("CB"+std::to_string(m),coeffB[m]);
            emitSum("CA"+std::to_string(m),coeffA[m]);
            cB<<"*CB"<<m<<",";
            cA<<"*CA"<<m<<",";
        }
        ofnc<<"RooRealVar& bkga = (RooRealVar&)*P.find(\"bkga\");"<<std::endl;
        ofnc<<"RooPolynomial* flatPosA = new RooPolynomial(\"flatPosA\",\"\",x_pos);"<<std::endl;
        ofnc<<"RooPolynomial* flatNegA = new RooPolynomial(\"flatNegA\",\"\",x_neg);"<<std::endl;
        ofnc<<"RooRealSumPdf* betaPos = new RooRealSumPdf(\"betaPos\",\"\",RooArgList("<<shapes.str()<<"*flatPos),RooArgList("<<cB.str()<<"bkg),true);"<<std::endl;
        ofnc<<"RooRealSumPdf* alphaPos = new RooRealSumPdf(\"alphaPos\",\"\",RooArgList("<<shapes.str()<<"*flatPosA),RooArgList("<<cA.str()<<"bkga),true);"<<std::endl;
        ofnc<<"RooRealSumPdf* alphaNeg = new RooRealSumPdf(\"alphaNeg\",\"\",RooArgList(*flatNegA),RooArgList(bkga),true);"<<std::endl;
        ofnc<<"RooFormulaVar* bkgAll = new RooFormulaVar(\"bkgAll\",\"bkg+bkga\",RooArgList(bkg,bkga));"<<std::endl;
        ofnc<<"RooRealSumPdf* CombinedBkgPosGated = new RooRealSumPdf(\"CombinedBkgPosGated\",\"\",RooArgList(*flatPos),RooArgList(*bkgAll),true);"<<std::endl;
        ofnc<<"RooCategory* gcat = new RooCategory(\"gcat\",\"\");"<<std::endl;
        ofnc<<"gcat->defineType(\"bpos\",0); gcat->defineType(\"bneg\",1); gcat->defineType(\"apos\",2); gcat->defineType(\"aneg\",3);"<<std::endl;
        ofnc<<"RooSimultaneous* simPdfGated = new RooSimultaneous(\"simPdfGated\",\"\",*gcat);"<<std::endl;
        ofnc<<"simPdfGated->addPdf(*betaPos,\"bpos\"); simPdfGated->addPdf(*bkgNeg,\"bneg\");"<<std::endl;
        ofnc<<"simPdfGated->addPdf(*alphaPos,\"apos\"); simPdfGated->addPdf(*alphaNeg,\"aneg\");"<<std::endl;
        ofnc<<"M->gcat=gcat; M->simPdfGated=simPdfGated; M->betaPos=betaPos; M->alphaPos=alphaPos; M->alphaNeg=alphaNeg;"<<std::endl;
        ofnc<<"M->CombinedBkgPosGated=CombinedBkgPosGated;"<<std::endl;
        std::string cBs=cB.str(); cBs.pop_back();
        std::string cAs=cA.str(); cAs.pop_back();
        std::string shapesNoFlatG=shapes.str(); shapesNoFlatG.pop_back();
        ofnc<<"M->betaSignalPos = new RooRealSumPdf(\"betaSignalPos\",\"\",RooArgList("<<shapesNoFlatG<<"),RooArgList("<<cBs<<"),true);"<<std::endl;
        ofnc<<"M->alphaSignalPos = new RooRealSumPdf(\"alphaSignalPos\",\"\",RooArgList("<<shapesNoFlatG<<"),RooArgList("<<cAs<<"),true);"<<std::endl;
    }
    ofnc<<"RooCategory* cat = new RooCategory(\"cat\",\"\");"<<std::endl;
    ofnc<<"cat->defineType(\"pos\",0); cat->defineType(\"neg\",1);"<<std::endl;
    ofnc<<"RooSimultaneous* simPdf = new RooSimultaneous(\"simPdf\",\"\",*cat);"<<std::endl;
    ofnc<<"simPdf->addPdf(*decayPos,\"pos\"); simPdf->addPdf(*bkgNeg,\"neg\");"<<std::endl;
    ofnc<<"M->cat=cat; M->simPdf=simPdf; M->decayPos=decayPos; M->bkgNeg=bkgNeg;"<<std::endl;
    ofnc<<"M->TrueNbeta=TrueNbetaAll;"<<std::endl;
    ofnc<<"M->CombinedBkgPos=CombinedBkgPos; M->ActivityPdfParentPos=ActivityPdf0; M->ActivityPdfDaughtersPos=ActivityPdfDaughtersPos;"<<std::endl;
    ofnc<<"return M;"<<std::endl;
    ofnc<<"}"<<std::endl;
}

//! Rename identifiers inside the string literals of the generated file (parameter
//! names in formulas and P.find() keys, RooFit object names) with a prefix, and
//! the builder function, so networks generated with different prefixes can be
//! linked and fitted together. Shared parameters and category labels keep
//! their names.
void applyPrefix(const std::string& file, const std::string& prefix)
{
    std::ifstream in(file);
    std::stringstream buf; buf<<in.rdbuf();
    std::string src=buf.str();
    const std::set<std::string> shared={"be","ab","bkg","bkga","pos","neg","bpos","bneg","apos","aneg"};
    const std::regex strLit("\"([^\"]*)\"");
    const std::regex ident("[A-Za-z_][A-Za-z0-9_]*");
    std::string out;
    auto it=std::sregex_iterator(src.begin(),src.end(),strLit);
    size_t last=0;
    for (; it!=std::sregex_iterator(); ++it){
        const std::smatch& m=*it;
        out+=src.substr(last,m.position(0)-last);
        std::string lit=m[1].str(), renamed;
        size_t l2=0;
        for (auto jt=std::sregex_iterator(lit.begin(),lit.end(),ident); jt!=std::sregex_iterator(); ++jt){
            renamed+=lit.substr(l2,jt->position(0)-l2);
            std::string id=jt->str();
            renamed+= shared.count(id) ? id : prefix+id;
            l2=jt->position(0)+id.size();
        }
        renamed+=lit.substr(l2);
        out+="\""+renamed+"\"";
        last=m.position(0)+m.length(0);
    }
    out+=src.substr(last);
    const std::string fn="DecayModel* buildDecayModel(";
    size_t pos=out.find(fn);
    if (pos!=std::string::npos) out.replace(pos,fn.size(),"DecayModel* buildDecayModel_"+prefix+"(");
    std::ofstream(file)<<out;
}

int main(int argc, char *argv[])
{
    if (argc<2){
        std::cerr<<"Usage: "<<argv[0]<<" <parmsex_file> [prefix]"<<std::endl;
        return 1;
    }
    const std::string prefix = (argc>2) ? argv[2] : "";
    decaypath* fdecaypath = new decaypath;
    fdecaypath->Init(argv[1]);
    fdecaypath->makePath();
    fdecaypath->printMember();
    fdecaypath->printPath();
    fdecaypath->writePath();
    path* fpath = readPathFile();
    std::vector<bool> hasAlpha(fdecaypath->getNMember());
    for (Int_t k=0;k<fdecaypath->getNMember();k++){
        MemberDef* m=fdecaypath->getMember(k);
        hasAlpha[k] = m->decay_abr>0 || m->is_decay_abr_fix==0;
        if (hasAlpha[k]) std::cout<<"alpha branch: "<<m->name<<" "<<m->decay_abr<<" %"<<(m->is_decay_abr_fix==0?" (floating)":"")<<std::endl;
    }
    generateDecayCurveModel(fpath, hasAlpha);
    if (!prefix.empty()) applyPrefix("decayModel_cal.cc", prefix);
    return 0;
}
