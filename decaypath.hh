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

//! \file decaypath.hh
//! \brief Definition of the decaypath class

#ifndef decaypath_h
#define decaypath_h 1

#include "common.hh"


using namespace std;

class decaypath
{
  public:
    decaypath();
    virtual ~decaypath();
    void Init(char* inputParms);

    void ProcessMember(MemberDef* obj);
    void CopyMember(MemberDef* source, MemberDef* destination);

    void makePath();
    path* getDecayPath(){return fdecaypath;}

    void printPath();
    void printMember();
    void drawPath(char* outputFileName);

    void writePath();

    MemberDef* getMember(int rino){flistofdecaymember_it=std::next(flistofdecaymember.begin(), rino);return (MemberDef*) *flistofdecaymember_it;}
    Int_t getNMember(){return (Int_t)flistofdecaymember.size();}
 private:
    char* finputParms;
    path* fdecaypath;
    std::list<MemberDef*> flistofdecaymember;
    std::list<MemberDef*>::iterator flistofdecaymember_it;
};

#endif
