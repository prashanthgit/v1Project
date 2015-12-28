
#ifndef alambda_h
#define alambda_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

const static Int_t nCuts = 31;
Bool_t flag[nCuts];

class alambda {
  public :

    Int_t mRapidityBin;
    Int_t mPhiPsiBin;
    Int_t mFlag;
    TString mOutputFileName;

    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        mVzLB;
   Double_t        mDcaV0ToPVLB;
   Double_t        mDecayLengthLB;
   Double_t        mDcaDaughtersLB;
   Double_t        mMomentumLB;
   Double_t        mPtLB;
   Double_t        mV0TypeLB;
   Double_t        mNSigmaPLB;
   Double_t        mNSigmaPiLB;
   Double_t        mMass2PLB;
   Double_t        mMass2PiLB;
   Double_t        mDcaPToPVLB;
   Double_t        mDcaPiToPVLB;
   Double_t        mNHitsPLB;
   Double_t        mNHitsPiLB;
   Double_t        mHitsPossPLB;
   Double_t        mHitsPossPiLB;
   Double_t        mRapidityLB;
   Double_t        mPhiPsiLB;
   Double_t        mPhiLB;
   Double_t        mLBarMass;
   Double_t        mEtaLB;
   UShort_t        mCentralityBin;

   // List of branches
   TBranch        *b_mVzLB;   //!
   TBranch        *b_mDcaV0ToPVLB;   //!
   TBranch        *b_mDecayLengthLB;   //!
   TBranch        *b_mDcaDaughtersLB;   //!
   TBranch        *b_mMomentumLB;   //!
   TBranch        *b_mPtLB;   //!
   TBranch        *b_mV0TypeLB;   //!
   TBranch        *b_mNSigmaPLB;   //!
   TBranch        *b_mNSigmaPiLB;   //!
   TBranch        *b_mMass2PLB;   //!
   TBranch        *b_mMass2PiLB;   //!
   TBranch        *b_mDcaPToPVLB;   //!
   TBranch        *b_mDcaPiToPVLB;   //!
   TBranch        *b_mNHitsPLB;   //!
   TBranch        *b_mNHitsPiLB;   //!
   TBranch        *b_mHitsPossPLB;   //!
   TBranch        *b_mHitsPossPiLB;   //!
   TBranch        *b_mRapidityLB;   //!
   TBranch        *b_mPhiPsiLB;   //!
   TBranch        *b_mPhiLB;   //!
   TBranch        *b_mLBarMass;   //!
   TBranch        *b_mEtaLB;   //!
   TBranch        *b_mCentralityBin;   //!

   alambda(Char_t path_file[100]);
   virtual ~alambda();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
    virtual void     RapidityBin(Double_t);
    virtual void     PhiPsiBin(Double_t);
};

#endif

#ifdef alambda_cxx
alambda::alambda(Char_t path_file[100]) : fChain(0) { 
  TString tmpName = path_file;
  TString fileBaseName(gSystem->BaseName(tmpName));
  fileBaseName.ReplaceAll(".root",".1.root");
  mOutputFileName = fileBaseName;

  cout << path_file << endl;
  cout << mOutputFileName <<endl;

  TTree *tree;
  TFile *f = new TFile(path_file); 
  f->GetObject("alambda",tree);

  Init(tree);

}

alambda::~alambda()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t alambda::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t alambda::LoadTree(Long64_t entry)
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

void alambda::Init(TTree *tree) {

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mVzLB", &mVzLB, &b_mVzLB);
   fChain->SetBranchAddress("mDcaV0ToPVLB", &mDcaV0ToPVLB, &b_mDcaV0ToPVLB);
   fChain->SetBranchAddress("mDecayLengthLB", &mDecayLengthLB, &b_mDecayLengthLB);
   fChain->SetBranchAddress("mDcaDaughtersLB", &mDcaDaughtersLB, &b_mDcaDaughtersLB);
   fChain->SetBranchAddress("mMomentumLB", &mMomentumLB, &b_mMomentumLB);
   fChain->SetBranchAddress("mPtLB", &mPtLB, &b_mPtLB);
   fChain->SetBranchAddress("mV0TypeLB", &mV0TypeLB, &b_mV0TypeLB);
   fChain->SetBranchAddress("mNSigmaPLB", &mNSigmaPLB, &b_mNSigmaPLB);
   fChain->SetBranchAddress("mNSigmaPiLB", &mNSigmaPiLB, &b_mNSigmaPiLB);
   fChain->SetBranchAddress("mMass2PLB", &mMass2PLB, &b_mMass2PLB);
   fChain->SetBranchAddress("mMass2PiLB", &mMass2PiLB, &b_mMass2PiLB);
   fChain->SetBranchAddress("mDcaPToPVLB", &mDcaPToPVLB, &b_mDcaPToPVLB);
   fChain->SetBranchAddress("mDcaPiToPVLB", &mDcaPiToPVLB, &b_mDcaPiToPVLB);
   fChain->SetBranchAddress("mNHitsPLB", &mNHitsPLB, &b_mNHitsPLB);
   fChain->SetBranchAddress("mNHitsPiLB", &mNHitsPiLB, &b_mNHitsPiLB);
   fChain->SetBranchAddress("mHitsPossPLB", &mHitsPossPLB, &b_mHitsPossPLB);
   fChain->SetBranchAddress("mHitsPossPiLB", &mHitsPossPiLB, &b_mHitsPossPiLB);
   fChain->SetBranchAddress("mRapidityLB", &mRapidityLB, &b_mRapidityLB);
   fChain->SetBranchAddress("mPhiPsiLB", &mPhiPsiLB, &b_mPhiPsiLB);
   fChain->SetBranchAddress("mPhiLB", &mPhiLB, &b_mPhiLB);
   fChain->SetBranchAddress("mLBarMass", &mLBarMass, &b_mLBarMass);
   fChain->SetBranchAddress("mEtaLB", &mEtaLB, &b_mEtaLB);
   fChain->SetBranchAddress("mCentralityBin", &mCentralityBin, &b_mCentralityBin);
}

void alambda::RapidityBin(Double_t yy){
  for(int k = -5; k<5; k++){
    Double_t Yrange1 = k*(2./10);
    Double_t Yrange2 = k*(2./10) + (2./10);
    if( (Yrange1 <= yy) && (yy <= Yrange2 )){
      return k+5;
    }
  }
  return 999;

}

void alambda::PhiPsiBin(Double_t phiRp){
  float two_pi = 6.35;
  float phiBinsMixed = 30;
  for(int i=0;i<phiBinsMixed;i++){
    Double_t RPrange1 = i*(two_pi/phiBinsMixed);
    Double_t RPrange2 = i*(two_pi/phiBinsMixed) + (two_pi/phiBinsMixed);
    if( (RPrange1 <= phiRp) && (phiRp <= RPrange2) ){
      return i;

    }
  }
  return 999;
}


#endif // #ifdef alambda_cxx
