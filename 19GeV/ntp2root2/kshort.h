
#ifndef kshort_h
#define kshort_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

const static Int_t nCuts = 19;
Bool_t flag[nCuts];

class kshort {
  public :
    Int_t mRapidityBin;
    Int_t mPhiPsiBin;
    Int_t mFlag;
    TString mOutputFileName;

    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    Double_t        mVzK;
    Double_t        mDcaV0ToPVK;
    Double_t        mDecayLengthK;
    Double_t        mDcaDaughtersK;
    Double_t        mMomentumK;
    Double_t        mPtK;
    Double_t        mV0TypeK;
    Double_t        mNSigmaPiPK;
    Double_t        mNSigmaPiNK;
    Double_t        mMass2PiPK;
    Double_t        mMass2PiNK;
    Double_t        mDcaPiPToPVK;
    Double_t        mDcaPiNToPVK;
    Double_t        mNHitsPiPK;
    Double_t        mNHitsPiNK;
    Double_t        mHitsPossPiPK;
    Double_t        mHitsPossPiNK;
    Double_t        mRapidityK;
    Double_t        mPhiPsiK;
    Double_t        mPhiK;
    Double_t        mKshortMass;
    Double_t        mEtaK;
    UShort_t        mCentralityBin;

    // List of branches
    TBranch        *b_mVzK;   //!
    TBranch        *b_mDcaV0ToPVK;   //!
    TBranch        *b_mDecayLengthK;   //!
    TBranch        *b_mDcaDaughtersK;   //!
    TBranch        *b_mMomentumK;   //!
    TBranch        *b_mPtK;   //!
    TBranch        *b_mV0TypeK;   //!
    TBranch        *b_mNSigmaPiPK;   //!
    TBranch        *b_mNSigmaPiNK;   //!
    TBranch        *b_mMass2PiPK;   //!
    TBranch        *b_mMass2PiNK;   //!
    TBranch        *b_mDcaPiPToPVK;   //!
    TBranch        *b_mDcaPiNToPVK;   //!
    TBranch        *b_mNHitsPiPK;   //!
    TBranch        *b_mNHitsPiNK;   //!
    TBranch        *b_mHitsPossPiPK;   //!
    TBranch        *b_mHitsPossPiNK;   //!
    TBranch        *b_mRapidityK;   //!
    TBranch        *b_mPhiPsiK;   //!
    TBranch        *b_mPhiK;   //!
    TBranch        *b_mKshortMass;   //!
    TBranch        *b_mEtaK;   //!
    TBranch        *b_mCentralityBin;   //!

    kshort(Char_t path_file[100]);
    virtual ~kshort();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual void     RapidityBin(Double_t);
    virtual void     PhiPsiBin(Double_t);
};

#endif

#ifdef kshort_cxx
kshort::kshort(Char_t path_file[100]) : fChain(0) {

  TString tmpName = path_file;
  TString fileBaseName(gSystem->BaseName(tmpName));
  fileBaseName.ReplaceAll(".root",".1.root");
  mOutputFileName = fileBaseName;

  cout << path_file << endl;
  cout << mOutputFileName <<endl;

  TTree *tree;
  TFile *f = new TFile(path_file); 
  f->GetObject("kshort",tree);

  Init(tree);

}

kshort::~kshort()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t kshort::GetEntry(Long64_t entry) {
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t kshort::LoadTree(Long64_t entry) {
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
  }
  return centry;
}

void kshort::Init(TTree *tree) {
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("mVzK", &mVzK, &b_mVzK);
  fChain->SetBranchAddress("mDcaV0ToPVK", &mDcaV0ToPVK, &b_mDcaV0ToPVK);
  fChain->SetBranchAddress("mDecayLengthK", &mDecayLengthK, &b_mDecayLengthK);
  fChain->SetBranchAddress("mDcaDaughtersK", &mDcaDaughtersK, &b_mDcaDaughtersK);
  fChain->SetBranchAddress("mMomentumK", &mMomentumK, &b_mMomentumK);
  fChain->SetBranchAddress("mPtK", &mPtK, &b_mPtK);
  fChain->SetBranchAddress("mV0TypeK", &mV0TypeK, &b_mV0TypeK);
  fChain->SetBranchAddress("mNSigmaPiPK", &mNSigmaPiPK, &b_mNSigmaPiPK);
  fChain->SetBranchAddress("mNSigmaPiNK", &mNSigmaPiNK, &b_mNSigmaPiNK);
  fChain->SetBranchAddress("mMass2PiPK", &mMass2PiPK, &b_mMass2PiPK);
  fChain->SetBranchAddress("mMass2PiNK", &mMass2PiNK, &b_mMass2PiNK);
  fChain->SetBranchAddress("mDcaPiPToPVK", &mDcaPiPToPVK, &b_mDcaPiPToPVK);
  fChain->SetBranchAddress("mDcaPiNToPVK", &mDcaPiNToPVK, &b_mDcaPiNToPVK);
  fChain->SetBranchAddress("mNHitsPiPK", &mNHitsPiPK, &b_mNHitsPiPK);
  fChain->SetBranchAddress("mNHitsPiNK", &mNHitsPiNK, &b_mNHitsPiNK);
  fChain->SetBranchAddress("mHitsPossPiPK", &mHitsPossPiPK, &b_mHitsPossPiPK);
  fChain->SetBranchAddress("mHitsPossPiNK", &mHitsPossPiNK, &b_mHitsPossPiNK);
  fChain->SetBranchAddress("mRapidityK", &mRapidityK, &b_mRapidityK);
  fChain->SetBranchAddress("mPhiPsiK", &mPhiPsiK, &b_mPhiPsiK);
  fChain->SetBranchAddress("mPhiK", &mPhiK, &b_mPhiK);
  fChain->SetBranchAddress("mKshortMass", &mKshortMass, &b_mKshortMass);
  fChain->SetBranchAddress("mEtaK", &mEtaK, &b_mEtaK);
  fChain->SetBranchAddress("mCentralityBin", &mCentralityBin, &b_mCentralityBin);
}

void kshort::RapidityBin(Double_t yy){
  for(int k = -5; k<5; k++){
    Double_t Yrange1 = k*(2./10);
    Double_t Yrange2 = k*(2./10) + (2./10);
    if( (Yrange1 <= yy) && (yy <= Yrange2 )){
      return k+5;
    }
  }
  return 999;

}

void kshort::PhiPsiBin(Double_t phiRp){
  float two_pi = 6.29;
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


#endif // #ifdef kshort_cxx
