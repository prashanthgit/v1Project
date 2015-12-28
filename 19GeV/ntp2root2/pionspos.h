
#ifndef pionspos_h
#define pionspos_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

const static Int_t nCuts = 16;
Bool_t flag[nCuts];

class pionspos {
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

   pionspos(Char_t path_file[100]);
   virtual ~pionspos();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
};

#endif

#ifdef pionspos_cxx
pionspos::pionspos(Char_t path_file[100]) : fChain(0) {
 TString tmpName = path_file;
  TString fileBaseName(gSystem->BaseName(tmpName));
  fileBaseName.ReplaceAll(".root",".1.root");
  mOutputFileName = fileBaseName;

  cout << path_file << endl;
  cout << mOutputFileName <<endl;

  TTree *tree;
  TFile *f = new TFile(path_file); 
  f->GetObject("pionpos",tree);

  Init(tree);

}

pionspos::~pionspos()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pionspos::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t pionspos::LoadTree(Long64_t entry) {
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

void pionspos::Init(TTree *tree) {

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


#endif // #ifdef pionspos_cxx
