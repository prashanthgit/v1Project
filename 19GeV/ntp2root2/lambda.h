
#ifndef lambda_h
#define lambda_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>

const static Int_t nCuts = 31;
Bool_t flag[nCuts];

class lambda {
  public :
    Int_t mRapidityBin;
    Int_t mPhiPsiBin;
    Int_t mFlag;
    TString mOutputFileName;

    TTree          *fChain;   
    Int_t           fCurrent; 

    // Declaration of leaf types
    Double_t        mVzL;
    Double_t        mDcaV0ToPVL;
    Double_t        mDecayLengthL;
    Double_t        mDcaDaughtersL;
    Double_t        mMomentumL;
    Double_t        mPtL;
    Double_t        mV0TypeL;
    Double_t        mNSigmaPL;
    Double_t        mNSigmaPiL;
    Double_t        mMass2PL;
    Double_t        mMass2PiL;
    Double_t        mDcaPToPVL;
    Double_t        mDcaPiToPVL;
    Double_t        mNHitsPL;
    Double_t        mNHitsPiL;
    Double_t        mHitsPossPL;
    Double_t        mHitsPossPiL;
    Double_t        mRapidityL;
    Double_t        mPhiPsiL;
    Double_t        mPhiL;
    Double_t        mLambdaMass;
    Double_t        mEtaL;
    UShort_t        mCentralityBin;

    // List of branches
    TBranch        *b_mVzL;   //!
    TBranch        *b_mDcaV0ToPVL;   //!
    TBranch        *b_mDecayLengthL;   //!
    TBranch        *b_mDcaDaughtersL;   //!
    TBranch        *b_mMomentumL;   //!
    TBranch        *b_mPtL;   //!
    TBranch        *b_mV0TypeL;   //!
    TBranch        *b_mNSigmaPL;   //!
    TBranch        *b_mNSigmaPiL;   //!
    TBranch        *b_mMass2PL;   //!
    TBranch        *b_mMass2PiL;   //!
    TBranch        *b_mDcaPToPVL;   //!
    TBranch        *b_mDcaPiToPVL;   //!
    TBranch        *b_mNHitsPL;   //!
    TBranch        *b_mNHitsPiL;   //!
    TBranch        *b_mHitsPossPL;   //!
    TBranch        *b_mHitsPossPiL;   //!
    TBranch        *b_mRapidityL;   //!
    TBranch        *b_mPhiPsiL;   //!
    TBranch        *b_mPhiL;   //!
    TBranch        *b_mLambdaMass;   //!
    TBranch        *b_mEtaL;   //!
    TBranch        *b_mCentralityBin;   //!

    lambda(Char_t path_file[100]);
    virtual ~lambda();
    virtual Int_t    Cut(Long64_t);
    virtual Long64_t LoadTree(Long64_t);
    virtual void     Init(TTree);
    virtual void     Loop();
    virtual void     Finish();
    virtual void     RapidityBin(Double_t);
    virtual void     PhiPsiBin(Double_t);

};

#endif

#ifdef lambda_cxx
lambda::lambda(Char_t path_file[100]):fChain(0){

  TString tmpName = path_file;
  TString fileBaseName(gSystem->BaseName(tmpName));
  fileBaseName.ReplaceAll(".root",".1.root");
  mOutputFileName = fileBaseName;



  cout << path_file << endl;
  cout << mOutputFileName <<endl;

  TTree *tree;
  TFile *f = new TFile(path_file); 
  f->GetObject("lambda",tree);

  Init(tree);


}

lambda::~lambda() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Long64_t lambda::LoadTree(Long64_t entry) {
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
  }
  return centry;
}

void lambda::Init(TTree *tree) {
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("mVzL", &mVzL, &b_mVzL);
  fChain->SetBranchAddress("mDcaV0ToPVL", &mDcaV0ToPVL, &b_mDcaV0ToPVL);
  fChain->SetBranchAddress("mDecayLengthL", &mDecayLengthL, &b_mDecayLengthL);
  fChain->SetBranchAddress("mDcaDaughtersL", &mDcaDaughtersL, &b_mDcaDaughtersL);
  fChain->SetBranchAddress("mMomentumL", &mMomentumL, &b_mMomentumL);
  fChain->SetBranchAddress("mPtL", &mPtL, &b_mPtL);
  fChain->SetBranchAddress("mV0TypeL", &mV0TypeL, &b_mV0TypeL);
  fChain->SetBranchAddress("mNSigmaPL", &mNSigmaPL, &b_mNSigmaPL);
  fChain->SetBranchAddress("mNSigmaPiL", &mNSigmaPiL, &b_mNSigmaPiL);
  fChain->SetBranchAddress("mMass2PL", &mMass2PL, &b_mMass2PL);
  fChain->SetBranchAddress("mMass2PiL", &mMass2PiL, &b_mMass2PiL);
  fChain->SetBranchAddress("mDcaPToPVL", &mDcaPToPVL, &b_mDcaPToPVL);
  fChain->SetBranchAddress("mDcaPiToPVL", &mDcaPiToPVL, &b_mDcaPiToPVL);
  fChain->SetBranchAddress("mNHitsPL", &mNHitsPL, &b_mNHitsPL);
  fChain->SetBranchAddress("mNHitsPiL", &mNHitsPiL, &b_mNHitsPiL);
  fChain->SetBranchAddress("mHitsPossPL", &mHitsPossPL, &b_mHitsPossPL);
  fChain->SetBranchAddress("mHitsPossPiL", &mHitsPossPiL, &b_mHitsPossPiL);
  fChain->SetBranchAddress("mRapidityL", &mRapidityL, &b_mRapidityL);
  fChain->SetBranchAddress("mPhiPsiL", &mPhiPsiL, &b_mPhiPsiL);
  fChain->SetBranchAddress("mPhiL", &mPhiL, &b_mPhiL);
  fChain->SetBranchAddress("mLambdaMass", &mLambdaMass, &b_mLambdaMass);
  fChain->SetBranchAddress("mEtaL", &mEtaL, &b_mEtaL);
  fChain->SetBranchAddress("mCentralityBin", &mCentralityBin, &b_mCentralityBin);

}

void lambda::Finish(){
}


Int_t lambda::Cut(Long64_t entry) {
}


void lambda::RapidityBin(Double_t yy){
  for(int k = -5; k<5; k++){
    Double_t Yrange1 = k*(2./10);
    Double_t Yrange2 = k*(2./10) + (2./10);
    if( (Yrange1 <= yy) && (yy <= Yrange2 )){
      return k+5;
    }
  }
  return 999;

}

void lambda::PhiPsiBin(Double_t phiRp){
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

#endif // #ifdef lambda_cxx
