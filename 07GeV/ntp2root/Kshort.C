
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "Constants.h"


TChain tree("kshortBg");
const static Int_t nCuts = 19;
Bool_t flag[nCuts];
Int_t mRapidityBin;
Int_t mPhiPsiBin;
Int_t mPtBin;
Int_t mFlag;
TString mOutputFileName;

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

void Kshort(Char_t path_file[100],Char_t jobid[100]){

  mOutputFileName = Form("%s_%s.root",histname[8].c_str(),jobid);

  std::ifstream inputFile(path_file);
  std::string line;
  while(getline(inputFile, line)) {
    cout << line <<endl;
    TFile *f = new TFile(line.c_str());
    tree.Add(line.c_str());
  }

  if (!tree) return;
  Long64_t nentries = tree.GetEntries();
  cout << nentries <<endl;


  tree.SetBranchAddress("mVzK", &mVzK, &b_mVzK);
  tree.SetBranchAddress("mDcaV0ToPVK", &mDcaV0ToPVK, &b_mDcaV0ToPVK);
  tree.SetBranchAddress("mDecayLengthK", &mDecayLengthK, &b_mDecayLengthK);
  tree.SetBranchAddress("mDcaDaughtersK", &mDcaDaughtersK, &b_mDcaDaughtersK);
  tree.SetBranchAddress("mMomentumK", &mMomentumK, &b_mMomentumK);
  tree.SetBranchAddress("mPtK", &mPtK, &b_mPtK);
  tree.SetBranchAddress("mV0TypeK", &mV0TypeK, &b_mV0TypeK);
  tree.SetBranchAddress("mNSigmaPiPK", &mNSigmaPiPK, &b_mNSigmaPiPK);
  tree.SetBranchAddress("mNSigmaPiNK", &mNSigmaPiNK, &b_mNSigmaPiNK);
  tree.SetBranchAddress("mMass2PiPK", &mMass2PiPK, &b_mMass2PiPK);
  tree.SetBranchAddress("mMass2PiNK", &mMass2PiNK, &b_mMass2PiNK);
  tree.SetBranchAddress("mDcaPiPToPVK", &mDcaPiPToPVK, &b_mDcaPiPToPVK);
  tree.SetBranchAddress("mDcaPiNToPVK", &mDcaPiNToPVK, &b_mDcaPiNToPVK);
  tree.SetBranchAddress("mNHitsPiPK", &mNHitsPiPK, &b_mNHitsPiPK);
  tree.SetBranchAddress("mNHitsPiNK", &mNHitsPiNK, &b_mNHitsPiNK);
  tree.SetBranchAddress("mHitsPossPiPK", &mHitsPossPiPK, &b_mHitsPossPiPK);
  tree.SetBranchAddress("mHitsPossPiNK", &mHitsPossPiNK, &b_mHitsPossPiNK);
  tree.SetBranchAddress("mRapidityK", &mRapidityK, &b_mRapidityK);
  tree.SetBranchAddress("mPhiPsiK", &mPhiPsiK, &b_mPhiPsiK);
  tree.SetBranchAddress("mPhiK", &mPhiK, &b_mPhiK);
  tree.SetBranchAddress("mKshortMass", &mKshortMass, &b_mKshortMass);
  tree.SetBranchAddress("mEtaK", &mEtaK, &b_mEtaK);
  tree.SetBranchAddress("mCentralityBin", &mCentralityBin, &b_mCentralityBin);

  fOutputfile = new TFile(mOutputFileName,"recreate");
  TTree* hhh = new TTree("v0","RECREATE");
  hhh -> Branch("mFlag",&mFlag,"mFlag/s");
  hhh -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");
  hhh -> Branch("mRapidityBin",&mRapidityBin,"mRapidityBin/s");
  hhh -> Branch("mPhiPsiBin",&mPhiPsiBin,"mPhiPsiBin/s");
  hhh -> Branch("mPtBin",&mPtBin,"mPtBin/s");
  hhh -> Branch("mKshortMass",&mKshortMass,"mKshortMass/D");




  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = tree.GetEntry(jentry);   nbytes += nb;

    mRapidityBin = RapidityBin(mRapidityK);
    mPhiPsiBin = PhiPsiBin(mPhiPsiK);
    mPtBin = GetPtBin16(mPtK);

    if( mRapidityBin == 999 || mPhiPsiBin == 999 || mPtBin == -999 ) continue;


    Bool_t iCut[9][3]={{0}};

    // Event 
    if(TMath::Abs(mVzK) <= vzDef)								iCut[0][0] =1;
    if(TMath::Abs(mVzK) <= vzMin)								iCut[0][1] =1;
    if(TMath::Abs(mVzK) <= vzMax)								iCut[0][2] =1;

    //  Track
    if(mNHitsPiPK >= nhitsDef && mNHitsPiNK >= nhitsDef)					iCut[1][0] =1;
    if(mNHitsPiPK >= nhitsMin && mNHitsPiNK >= nhitsMin)					iCut[1][1] =1;
    if(mNHitsPiPK >= nhitsMax && mNHitsPiNK >= nhitsMax)					iCut[1][2] =1;

    if(mHitsPossPiPK >= nhitsposDef && mHitsPossPiNK >= nhitsposDef)				iCut[2][0] =1;
    if(mHitsPossPiPK >  nhitsposMin && mHitsPossPiNK >= nhitsposMin)  			        iCut[2][1] =1;
    if(mHitsPossPiPK >= nhitsposMax && mHitsPossPiNK >= nhitsposMax)				iCut[2][2] =1;

    if( (TMath::Abs(mNSigmaPiPK) <= nsigmav0Def) && ( TMath::Abs(mNSigmaPiNK) <= nsigmav0Def) )	iCut[3][0] =1;
    if( (TMath::Abs(mNSigmaPiPK) <= nsigmav0Min) && ( TMath::Abs(mNSigmaPiNK) <= nsigmav0Min) )	iCut[3][1] =1;
    if( (TMath::Abs(mNSigmaPiPK) <= nsigmav0Max) && ( TMath::Abs(mNSigmaPiNK) <= nsigmav0Max) )	iCut[3][2] =1;

    // Momentum
    if( (kshortPtMinDef <= mPtK) && (mPtK <= kshortPtMaxDef))  					iCut[4][0] =1;
    if(mMomentumK <= kshortMomMaxMes && kshortPtMinMes <= mPtK )				iCut[4][1] =1;
    //if( (mMomentumK <= kshortMomMaxBar) && (kshortPtMinBar <= mPtK) && (mPtK <= kshortPtMaxBar))iCut[4][2] =2;

    if( mDcaV0ToPVK < kshortDcaV0PvDef)								iCut[5][0] =1;
    if( mDcaV0ToPVK < kshortDcaV0PvMin)								iCut[5][1] =1;
    if( mDcaV0ToPVK < kshortDcaV0PvMax)								iCut[5][2] =1;

    if( mDecayLengthK > kashortDecaylengthDef )							iCut[6][0] =1;
    if( mDecayLengthK > kashortDecaylengthMin )							iCut[6][1] =1;
    if( mDecayLengthK > kashortDecaylengthMax )							iCut[6][2] =1;

    if( mDcaPiPToPVK > kshortDcaPitoPvDef && mDcaPiNToPVK > kshortDcaPitoPvDef)			iCut[7][0] =1;
    if( mDcaPiPToPVK > kshortDcaPitoPvMin && mDcaPiNToPVK > kshortDcaPitoPvMin)			iCut[7][1] =1;
    if( mDcaPiPToPVK > kshortDcaPitoPvMax && mDcaPiNToPVK > kshortDcaPitoPvMax)			iCut[7][2] =1;


    if( mDcaDaughtersK < kshortDcaDaughtersDef)							iCut[8][0] =1;
    if( mDcaDaughtersK < kshortDcaDaughtersMin)							iCut[8][1] =1;
    if( mDcaDaughtersK < kshortDcaDaughtersMax)							iCut[8][2] =1;

    ////////////// ////////////// ////////////// //////////////
    ////////////// ////////////// ////////////// //////////////


    mFlag = 0;
    flag[mFlag]=1;
    for(int i=0;i<9;i++) {flag[mFlag] =  flag[mFlag] * iCut[i][0];} 
    if(flag[mFlag]) hhh -> Fill();
    mFlag++;


    for(int i=0;i<9;i++){
      for(int j=1;j<=2;j++){
	flag[mFlag] = iCut[i][j];
	if(!flag[mFlag]){mFlag++;continue;}
	for(int k=0;k<9;k++){
	  if(i==k)continue;
	  flag[mFlag] = flag[mFlag] * iCut[k][0];
	}
	if(flag[mFlag]) hhh -> Fill();
	mFlag++;
      }
    }

    ////////////// ////////////// ////////////// //////////////
    ////////////// ////////////// ////////////// //////////////

    if(( jentry %10000) == 0 ) cout << jentry  <<endl;
  }


  hhh -> Write();
  fOutputfile-> Close();




}






void RapidityBin(Double_t yy){
  for(int k = -10; k<10; k++){
    Double_t Yrange1 = k*(1./10);
    Double_t Yrange2 = k*(1./10) + (1./10);
    if( (Yrange1 <= yy) && (yy <= Yrange2 )){
      return k+10;
    }
  }
  return 999;

}

void PhiPsiBin(Double_t phiRp){
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


int GetPtBin16(Float_t pt){
  for(int i=0;i<16;i++) if( (pT_min[i] < pt) && (pt <= pT_max[i] )) return i;
  cout << " some think wrong" <<endl;
  return -999;
}
