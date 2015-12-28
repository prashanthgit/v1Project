
#ifndef proton_h
#define proton_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <cstdio>

const static Int_t nCuts = 16;
Bool_t flag[nCuts];

class proton {
  public :
    Int_t mFlag;
    TString mOutputFileName;
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    Double_t        mVzT;
    Double_t        mMomentumT;
    Double_t        mPtT;
    Double_t        mNSigmaT;
    Double_t        mMass2T;
    Double_t        mDcaT;
    Double_t        mNHitsT;
    Double_t        mHitsPossT;
    Double_t        mRapidityT;
    Double_t        mPhiPsiT;
    Double_t        mPhiT;
    Double_t        mEtaT;
    UShort_t        mCentralityBin;
    Double_t        mV1Track;

    // List of branches
    TBranch        *b_mVzT;   //!
    TBranch        *b_mMomentumT;   //!
    TBranch        *b_mPtT;   //!
    TBranch        *b_mNSigmaT;   //!
    TBranch        *b_mMass2T;   //!
    TBranch        *b_mDcaT;   //!
    TBranch        *b_mNHitsT;   //!
    TBranch        *b_mHitsPossT;   //!
    TBranch        *b_mRapidityT;   //!
    TBranch        *b_mPhiPsiT;   //!
    TBranch        *b_mPhiT;   //!
    TBranch        *b_mEtaT;   //!
    TBranch        *b_mCentralityBin;   //!
    TBranch        *b_mV1Track;   //!

    proton(Char_t path_file[100]);
    virtual ~proton();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
};

#endif

#ifdef proton_cxx
proton::proton(Char_t path_file[100],Char_t jobid[100]) : fChain(0) {
  cout << jobid << endl;
  TChain tree("proton");
  std::ifstream inputFile(path_file);
  std::string line;
  while(getline(inputFile, line)) {
    if (line.length()){
      cout << line <<endl;
      tree.Add(line.c_str());
    }
  }

  mOutputFileName = Form("FlagP_%s.root",jobid); 

  cout << path_file << endl;
  cout << mOutputFileName <<endl;

  Long64_t nentries = tree.GetEntries();
  cout << "hello 123  " <<nentries <<endl;

  Init(tree);

}

proton::~proton() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t proton::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t proton::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
  }
  return centry;
}

void proton::Init(TTree *tree) {

  Long64_t nentries = tree.GetEntries();
  cout << "hello   " <<nentries <<endl;

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("mVzT", &mVzT, &b_mVzT);
  fChain->SetBranchAddress("mMomentumT", &mMomentumT, &b_mMomentumT);
  fChain->SetBranchAddress("mPtT", &mPtT, &b_mPtT);
  fChain->SetBranchAddress("mNSigmaT", &mNSigmaT, &b_mNSigmaT);
  fChain->SetBranchAddress("mMass2T", &mMass2T, &b_mMass2T);
  fChain->SetBranchAddress("mDcaT", &mDcaT, &b_mDcaT);
  fChain->SetBranchAddress("mNHitsT", &mNHitsT, &b_mNHitsT);
  fChain->SetBranchAddress("mHitsPossT", &mHitsPossT, &b_mHitsPossT);
  fChain->SetBranchAddress("mRapidityT", &mRapidityT, &b_mRapidityT);
  fChain->SetBranchAddress("mPhiPsiT", &mPhiPsiT, &b_mPhiPsiT);
  fChain->SetBranchAddress("mPhiT", &mPhiT, &b_mPhiT);
  fChain->SetBranchAddress("mEtaT", &mEtaT, &b_mEtaT);
  fChain->SetBranchAddress("mCentralityBin", &mCentralityBin, &b_mCentralityBin);
  fChain->SetBranchAddress("mV1Track", &mV1Track, &b_mV1Track);
}


#endif // #ifdef proton_cxx
