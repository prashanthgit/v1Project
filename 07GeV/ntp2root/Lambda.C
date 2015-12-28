
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "Constants.h"

TChain tree("lambdaBg");
const static Int_t nCuts = 31;
Bool_t flag[nCuts];
Int_t mRapidityBin;
Int_t mPhiPsiBin;
Int_t mPtBin;
Int_t mFlag;
TString mOutputFileName;

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


void Lambda(Char_t path_file[100],Char_t jobid[100]){

  mOutputFileName = Form("%s_%s.root",histname[6].c_str(),jobid);

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

  tree.SetBranchAddress("mVzL", &mVzL, &b_mVzL);
  tree.SetBranchAddress("mDcaV0ToPVL", &mDcaV0ToPVL, &b_mDcaV0ToPVL);
  tree.SetBranchAddress("mDecayLengthL", &mDecayLengthL, &b_mDecayLengthL);
  tree.SetBranchAddress("mDcaDaughtersL", &mDcaDaughtersL, &b_mDcaDaughtersL);
  tree.SetBranchAddress("mMomentumL", &mMomentumL, &b_mMomentumL);
  tree.SetBranchAddress("mPtL", &mPtL, &b_mPtL);
  tree.SetBranchAddress("mV0TypeL", &mV0TypeL, &b_mV0TypeL);
  tree.SetBranchAddress("mNSigmaPL", &mNSigmaPL, &b_mNSigmaPL);
  tree.SetBranchAddress("mNSigmaPiL", &mNSigmaPiL, &b_mNSigmaPiL);
  tree.SetBranchAddress("mMass2PL", &mMass2PL, &b_mMass2PL);
  tree.SetBranchAddress("mMass2PiL", &mMass2PiL, &b_mMass2PiL);
  tree.SetBranchAddress("mDcaPToPVL", &mDcaPToPVL, &b_mDcaPToPVL);
  tree.SetBranchAddress("mDcaPiToPVL", &mDcaPiToPVL, &b_mDcaPiToPVL);
  tree.SetBranchAddress("mNHitsPL", &mNHitsPL, &b_mNHitsPL);
  tree.SetBranchAddress("mNHitsPiL", &mNHitsPiL, &b_mNHitsPiL);
  tree.SetBranchAddress("mHitsPossPL", &mHitsPossPL, &b_mHitsPossPL);
  tree.SetBranchAddress("mHitsPossPiL", &mHitsPossPiL, &b_mHitsPossPiL);
  tree.SetBranchAddress("mRapidityL", &mRapidityL, &b_mRapidityL);
  tree.SetBranchAddress("mPhiPsiL", &mPhiPsiL, &b_mPhiPsiL);
  tree.SetBranchAddress("mPhiL", &mPhiL, &b_mPhiL);
  tree.SetBranchAddress("mLambdaMass", &mLambdaMass, &b_mLambdaMass);
  tree.SetBranchAddress("mEtaL", &mEtaL, &b_mEtaL);
  tree.SetBranchAddress("mCentralityBin", &mCentralityBin, &b_mCentralityBin);

  fOutputfile = new TFile(mOutputFileName,"recreate");
  TTree* hhh = new TTree("v0","RECREATE");
  hhh -> Branch("mFlag",&mFlag,"mFlag/s");
  hhh -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");
  hhh -> Branch("mRapidityBin",&mRapidityBin,"mRapidityBin/s");
  hhh -> Branch("mPhiPsiBin",&mPhiPsiBin,"mPhiPsiBin/s");
  hhh -> Branch("mPtBin",&mPtBin,"mPtBin/s");
  hhh -> Branch("mLambdaMass",&mLambdaMass,"mLambdaMass/D");



  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = tree.GetEntry(jentry);   nbytes += nb;

    // cuts 
    if(mLambdaMass < 1.06 || mLambdaMass > 1.24) continue; 

    mRapidityBin = RapidityBin(mRapidityL);
    mPhiPsiBin = PhiPsiBin(mPhiPsiL);
    mPtBin = GetPtBin16(mPtL);

    if( mRapidityBin == 999 || mPhiPsiBin == 999 || mPtBin == -999 ) continue;
    ////////////// ////////////// ////////////// //////////////
    ////////////// ////////////// ////////////// //////////////
    //
    Bool_t iCut[15][3]={{0}};

    // Event 
    if(TMath::Abs(mVzL) <= vzDef)			iCut[0][0] =1;
    if(TMath::Abs(mVzL) <= vzMin)			iCut[0][1] =1;
    if(TMath::Abs(mVzL) <= vzMax)			iCut[0][2] =1;

    // Positive Track
    if(mNHitsPL >= nhitsDef) 				iCut[1][0] =1;
    if(mNHitsPL >= nhitsMin) 				iCut[1][1] =1;
    if(mNHitsPL >= nhitsMax) 				iCut[1][2] =1;

    if(mHitsPossPL >= nhitsposDef)	 		iCut[2][0] =1;
    if(mHitsPossPL >  nhitsposMin)  			iCut[2][1] =1;
    if(mHitsPossPL >= nhitsposMax)	 		iCut[2][2] =1;

    if(TMath::Abs(mNSigmaPL) <= nsigmav0Def) 		iCut[3][0] =1;
    if(TMath::Abs(mNSigmaPL) <= nsigmav0Min) 		iCut[3][1] =1;
    if(TMath::Abs(mNSigmaPL) <= nsigmav0Max)		iCut[3][2] =1;

    if( (mV0TypeL ==2 || mV0TypeL == 3) && (mass2LprotonMinDef < mMass2PL) && (mMass2PL < mass2LprotonMaxDef) ){ iCut[4][0] =1; }else {	iCut[4][0] =1; }
    if( (mV0TypeL ==2 || mV0TypeL == 3) && (mass2LprotonMinMin < mMass2PL) && (mMass2PL < mass2LprotonMaxMin) ){ iCut[4][1] =1; }else { iCut[4][1] =1; } 

    // negative Track
    if(mNHitsPiL >= nhitsDef) 				iCut[5][0] =1;
    if(mNHitsPiL >= nhitsMin) 				iCut[5][1] =1;
    if(mNHitsPiL >= nhitsMax) 				iCut[5][2] =1;

    if(mHitsPossPiL > nhitsposDef)	 		iCut[6][0] =1;
    if(mHitsPossPiL > nhitsposMin) 			iCut[6][1] =1;
    if(mHitsPossPiL > nhitsposMax)	 		iCut[6][2] =1;

    if(TMath::Abs(mNSigmaPiL) <= nsigmav0Def) 		iCut[7][0] =1;
    if(TMath::Abs(mNSigmaPiL) <= nsigmav0Min) 		iCut[7][1] =1;
    if(TMath::Abs(mNSigmaPiL) <= nsigmav0Max) 		iCut[7][2] =1;

    iCut[8][0] = 1;
    if( (mV0TypeL ==1 || mV0TypeL == 3) && (mass2LpionMinMin <= mMass2PiL) && (mMass2PiL <= mass2LpionMaxMin) ){ iCut[8][1] =1; }else {	iCut[8][1] =1; } 

    // Momentum
    if( (lambdaPtMinDef <= mPtL) && (mPtL <= lambdaPtMaxDef)	)				iCut[9][0] =1;
    if(mMomentumL <= lambdaMomMaxBar && lambdaPtMinBar <= mPtL && mPtL <= lambdaPtMaxBar)	iCut[9][1] =1; // proton equvalance
    // if( mMomentumL <= lambdaMomMaxMes && lambdaPtMinMes <= mPtL)				iCut[9][2] =1;


    // V0 
    switch(mV0TypeL){
    case 0: {
	      if( mDcaV0ToPVL < lambdaDcaV0PvDef_0)			iCut[10][0] =1;
	      if( mDcaV0ToPVL < lambdaDcaV0PvMin_0)			iCut[10][1] =1;
	      if( mDcaV0ToPVL < lambdaDcaV0PvMax_0)			iCut[10][2] =1;


	      if( mDecayLengthL > lambdaDecaylengthDef_0 )		iCut[11][0] =1;
	      if( mDecayLengthL > lambdaDecaylengthMin_0 )		iCut[11][1] =1;
	      if( mDecayLengthL > lambdaDecaylengthMax_0 )		iCut[11][2] =1;


	      if( mDcaPToPVL > lambdaDcaPtoPvDef_0)			iCut[12][0] =1;
	      if( mDcaPToPVL > lambdaDcaPtoPvMin_0)			iCut[12][1] =1;
	      if( mDcaPToPVL > lambdaDcaPtoPvMax_0)			iCut[12][2] =1;


	      if( mDcaPiToPVL > lambdaDcaPitoPvDef_0 ) 			iCut[13][0] =1;
	      if( mDcaPiToPVL > lambdaDcaPitoPvMin_0 ) 			iCut[13][1] =1;
	      if( mDcaPiToPVL > lambdaDcaPitoPvMax_0 ) 			iCut[13][2] =1;


	      break;
	    }
    case 1: {
	      if( mDcaV0ToPVL < lambdaDcaV0PvDef_1)			iCut[10][0] =1;
	      if( mDcaV0ToPVL < lambdaDcaV0PvMin_1)			iCut[10][1] =1;
	      if( mDcaV0ToPVL < lambdaDcaV0PvMax_1)			iCut[10][2] =1;


	      if( mDecayLengthL > lambdaDecaylengthDef_1 )		iCut[11][0] =1;
	      if( mDecayLengthL > lambdaDecaylengthMin_1 )		iCut[11][1] =1;
	      if( mDecayLengthL > lambdaDecaylengthMax_1 )		iCut[11][2] =1;


	      if( mDcaPToPVL > lambdaDcaPtoPvDef_1)			iCut[12][0] =1;
	      if( mDcaPToPVL > lambdaDcaPtoPvMin_1)			iCut[12][1] =1;
	      if( mDcaPToPVL > lambdaDcaPtoPvMax_1)			iCut[12][2] =1;


	      if( mDcaPiToPVL > lambdaDcaPitoPvDef_1 ) 			iCut[13][0] =1;
	      if( mDcaPiToPVL > lambdaDcaPitoPvMin_1 ) 			iCut[13][1] =1;
	      if( mDcaPiToPVL > lambdaDcaPitoPvMax_1 ) 			iCut[13][2] =1;


	      break;
	    }
    case 2: {
	      if( mDcaV0ToPVL < lambdaDcaV0PvDef_2)			iCut[10][0] =1;
	      if( mDcaV0ToPVL < lambdaDcaV0PvMin_2)			iCut[10][1] =1;
	      if( mDcaV0ToPVL < lambdaDcaV0PvMax_2)			iCut[10][2] =1;


	      if( mDecayLengthL > lambdaDecaylengthDef_2 )		iCut[11][0] =1;
	      if( mDecayLengthL > lambdaDecaylengthMin_2 )		iCut[11][1] =1;
	      if( mDecayLengthL > lambdaDecaylengthMax_2 )		iCut[11][2] =1;


	      if( mDcaPToPVL > lambdaDcaPtoPvDef_2)			iCut[12][0] =1;
	      if( mDcaPToPVL > lambdaDcaPtoPvMin_2)			iCut[12][1] =1;
	      if( mDcaPToPVL > lambdaDcaPtoPvMax_2)			iCut[12][2] =1;


	      if( mDcaPiToPVL > lambdaDcaPitoPvDef_2 ) 			iCut[13][0] =1;
	      if( mDcaPiToPVL > lambdaDcaPitoPvMin_2 ) 			iCut[13][1] =1;
	      if( mDcaPiToPVL > lambdaDcaPitoPvMax_2 ) 			iCut[13][2] =1;


	      break;
	    }
    case 3: {
	      if( mDcaV0ToPVL < lambdaDcaV0PvDef_3)			iCut[10][0] =1;
	      if( mDcaV0ToPVL < lambdaDcaV0PvMin_3)			iCut[10][1] =1;
	      if( mDcaV0ToPVL < lambdaDcaV0PvMax_3)			iCut[10][2] =1;


	      if( mDecayLengthL > lambdaDecaylengthDef_3 )		iCut[11][0] =1;
	      if( mDecayLengthL > lambdaDecaylengthMin_3 )		iCut[11][1] =1;
	      if( mDecayLengthL > lambdaDecaylengthMax_3 )		iCut[11][2] =1;


	      if( mDcaPToPVL > lambdaDcaPtoPvDef_3)			iCut[12][0] =1;
	      if( mDcaPToPVL > lambdaDcaPtoPvMin_3)			iCut[12][1] =1;
	      if( mDcaPToPVL > lambdaDcaPtoPvMax_3)			iCut[12][2] =1;


	      if( mDcaPiToPVL > lambdaDcaPitoPvDef_3 ) 			iCut[13][0] =1;
	      if( mDcaPiToPVL > lambdaDcaPitoPvMin_3 ) 			iCut[13][1] =1;
	      if( mDcaPiToPVL > lambdaDcaPitoPvMax_3 ) 			iCut[13][2] =1;


	      break;
	    }

    }


    if( mDcaDaughtersL < lambdaDcaDaughtersDef)				iCut[14][0] =1;
    if( mDcaDaughtersL < lambdaDcaDaughtersMin)				iCut[14][1] =1;
    if( mDcaDaughtersL < lambdaDcaDaughtersMax)				iCut[14][2] =1;

    ////////////// ////////////// ////////////// //////////////
    ////////////// ////////////// ////////////// //////////////


    mFlag = 0;
    flag[mFlag]=1;
    for(int i=0;i<15;i++) {flag[mFlag] =  flag[mFlag] * iCut[i][0];} 
    if(flag[mFlag]) hhh -> Fill();
    mFlag++;


    for(int i=0;i<15;i++){
      for(int j=1;j<=2;j++){
	flag[mFlag] = iCut[i][j];
	if(!flag[mFlag]){mFlag++;continue;}
	for(int k=0;k<15;k++){
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
