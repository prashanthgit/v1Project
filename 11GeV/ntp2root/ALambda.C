
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "Constants.h"


TChain tree("alambdaBg");
const static Int_t nCuts = 31;
Bool_t flag[nCuts];
Int_t mRapidityBin;
Int_t mPhiPsiBin;
Int_t mPtBin;
Int_t mFlag;
TString mOutputFileName;

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

void ALambda(Char_t path_file[100],Char_t jobid[100]){

  mOutputFileName = Form("%s_%s.root",histname[7].c_str(),jobid);

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


  tree.SetBranchAddress("mVzLB", &mVzLB, &b_mVzLB);
  tree.SetBranchAddress("mDcaV0ToPVLB", &mDcaV0ToPVLB, &b_mDcaV0ToPVLB);
  tree.SetBranchAddress("mDecayLengthLB", &mDecayLengthLB, &b_mDecayLengthLB);
  tree.SetBranchAddress("mDcaDaughtersLB", &mDcaDaughtersLB, &b_mDcaDaughtersLB);
  tree.SetBranchAddress("mMomentumLB", &mMomentumLB, &b_mMomentumLB);
  tree.SetBranchAddress("mPtLB", &mPtLB, &b_mPtLB);
  tree.SetBranchAddress("mV0TypeLB", &mV0TypeLB, &b_mV0TypeLB);
  tree.SetBranchAddress("mNSigmaPLB", &mNSigmaPLB, &b_mNSigmaPLB);
  tree.SetBranchAddress("mNSigmaPiLB", &mNSigmaPiLB, &b_mNSigmaPiLB);
  tree.SetBranchAddress("mMass2PLB", &mMass2PLB, &b_mMass2PLB);
  tree.SetBranchAddress("mMass2PiLB", &mMass2PiLB, &b_mMass2PiLB);
  tree.SetBranchAddress("mDcaPToPVLB", &mDcaPToPVLB, &b_mDcaPToPVLB);
  tree.SetBranchAddress("mDcaPiToPVLB", &mDcaPiToPVLB, &b_mDcaPiToPVLB);
  tree.SetBranchAddress("mNHitsPLB", &mNHitsPLB, &b_mNHitsPLB);
  tree.SetBranchAddress("mNHitsPiLB", &mNHitsPiLB, &b_mNHitsPiLB);
  tree.SetBranchAddress("mHitsPossPLB", &mHitsPossPLB, &b_mHitsPossPLB);
  tree.SetBranchAddress("mHitsPossPiLB", &mHitsPossPiLB, &b_mHitsPossPiLB);
  tree.SetBranchAddress("mRapidityLB", &mRapidityLB, &b_mRapidityLB);
  tree.SetBranchAddress("mPhiPsiLB", &mPhiPsiLB, &b_mPhiPsiLB);
  tree.SetBranchAddress("mPhiLB", &mPhiLB, &b_mPhiLB);
  tree.SetBranchAddress("mLBarMass", &mLBarMass, &b_mLBarMass);
  tree.SetBranchAddress("mEtaLB", &mEtaLB, &b_mEtaLB);
  tree.SetBranchAddress("mCentralityBin", &mCentralityBin, &b_mCentralityBin);


  fOutputfile = new TFile(mOutputFileName,"recreate");
  TTree* hhh = new TTree("v0","RECREATE");
  hhh -> Branch("mFlag",&mFlag,"mFlag/s");
  hhh -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");
  hhh -> Branch("mRapidityBin",&mRapidityBin,"mRapidityBin/s");
  hhh -> Branch("mPhiPsiBin",&mPhiPsiBin,"mPhiPsiBin/s");
  hhh -> Branch("mPtBin",&mPtBin,"mPtBin/s");
  hhh -> Branch("mLBarMass",&mLBarMass,"mLBarMass/D");




  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = tree.GetEntry(jentry);   nbytes += nb;

    // cuts 
    if(mLBarMass < 1.06 || mLBarMass > 1.24) continue; 

    mRapidityBin = RapidityBin(mRapidityLB);
    mPhiPsiBin = PhiPsiBin(mPhiPsiLB);
    mPtBin = GetPtBin16(mPtLB);


    if( mRapidityBin == 999 || mPhiPsiBin == 999 || mPtBin == -999 ) continue;
    ////////////// ////////////// ////////////// //////////////
    ////////////// ////////////// ////////////// //////////////
    //
    Bool_t iCut[15][3]={{0}};

    // Event 
    if(TMath::Abs(mVzLB) <= vzDef)			iCut[0][0] =1;
    if(TMath::Abs(mVzLB) <= vzMin)			iCut[0][1] =1;
    if(TMath::Abs(mVzLB) <= vzMax)			iCut[0][2] =1;

    // Positive Track
    if(mNHitsPLB >= nhitsDef) 				iCut[1][0] =1;
    if(mNHitsPLB >= nhitsMin) 				iCut[1][1] =1;
    if(mNHitsPLB >= nhitsMax) 				iCut[1][2] =1;

    if(mHitsPossPLB >= nhitsposDef)	 		iCut[2][0] =1;
    if(mHitsPossPLB >  nhitsposMin)  			iCut[2][1] =1;
    if(mHitsPossPLB >= nhitsposMax)	 		iCut[2][2] =1;

    if(TMath::Abs(mNSigmaPLB) <= nsigmav0Def) 		iCut[3][0] =1;
    if(TMath::Abs(mNSigmaPLB) <= nsigmav0Min) 		iCut[3][1] =1;
    if(TMath::Abs(mNSigmaPLB) <= nsigmav0Max) 		iCut[3][2] =1;

    if( (mV0TypeLB ==2 || mV0TypeLB == 3) && (mass2LprotonMinDef < mMass2PLB) && (mMass2PLB < mass2LprotonMaxDef) ){ 	iCut[4][0] =1;}else {  iCut[4][0] =1;}
    if( (mV0TypeLB ==2 || mV0TypeLB == 3) && (mass2LprotonMinMin < mMass2PLB) && (mMass2PLB < mass2LprotonMaxMin) ){ 	iCut[4][1] =1;}else {  iCut[4][1] =1;}

    // negative Track
    if(mNHitsPiLB >= nhitsDef) 				iCut[5][0] =1;
    if(mNHitsPiLB >= nhitsMin) 				iCut[5][1] =1;
    if(mNHitsPiLB >= nhitsMax) 				iCut[5][2] =1;

    if(mHitsPossPiLB > nhitsposDef)	 		iCut[6][0] =1;
    if(mHitsPossPiLB > nhitsposMin) 			iCut[6][1] =1;
    if(mHitsPossPiLB > nhitsposMax)	 		iCut[6][2] =1;

    if(TMath::Abs(mNSigmaPiLB) <= nsigmav0Def) 		iCut[7][0] =1;
    if(TMath::Abs(mNSigmaPiLB) <= nsigmav0Min) 		iCut[7][1] =1;
    if(TMath::Abs(mNSigmaPiLB) <= nsigmav0Max)  	iCut[7][2] =1;

    iCut[8][0] = 1;
    if( (mV0TypeLB ==1 || mV0TypeLB == 3) && (mass2LpionMinMin <= mMass2PiLB) && (mMass2PiLB <= mass2LpionMaxMin) ){ iCut[8][1] =1;}else {  iCut[8][1] =1;}


    // Momentum
    if( (lambdaPtMinDef <= mPtLB) && (mPtLB <= lambdaPtMaxDef)	)	 			iCut[9][0] =1;
    if(mMomentumLB <= lambdaMomMaxBar && lambdaPtMinBar <= mPtLB && mPtLB <= lambdaPtMaxBar)	iCut[9][1] =1;
    //    if(mMomentumLB <= lambdaMomMaxMes && lambdaPtMinMes <= mPtLB)				iCut[9][2] =1;


    // V0 
    switch(mV0TypeLB){
    case 0: {
	      if( mDcaV0ToPVLB < lambdaDcaV0PvDef_0)			iCut[10][0] =1;
	      if( mDcaV0ToPVLB < lambdaDcaV0PvMin_0)			iCut[10][1] =1;
	      if( mDcaV0ToPVLB < lambdaDcaV0PvMax_0)			iCut[10][2] =1;


	      if( mDecayLengthLB > lambdaDecaylengthDef_0 )		iCut[11][0] =1;
	      if( mDecayLengthLB > lambdaDecaylengthMin_0 )		iCut[11][1] =1;
	      if( mDecayLengthLB > lambdaDecaylengthMax_0 )		iCut[11][2] =1;


	      if( mDcaPToPVLB > lambdaDcaPtoPvDef_0)			iCut[12][0] =1;
	      if( mDcaPToPVLB > lambdaDcaPtoPvMin_0)			iCut[12][1] =1;
	      if( mDcaPToPVLB > lambdaDcaPtoPvMax_0)			iCut[12][2] =1;


	      if( mDcaPiToPVLB > lambdaDcaPitoPvDef_0) 			iCut[13][0] =1;
	      if( mDcaPiToPVLB > lambdaDcaPitoPvMin_0) 			iCut[13][1] =1;
	      if( mDcaPiToPVLB > lambdaDcaPitoPvMax_0) 			iCut[13][2] =1;


	      break;
	    }
    case 1: {
	      if( mDcaV0ToPVLB < lambdaDcaV0PvDef_1)			iCut[10][0] =1;
	      if( mDcaV0ToPVLB < lambdaDcaV0PvMin_1)			iCut[10][1] =1;
	      if( mDcaV0ToPVLB < lambdaDcaV0PvMax_1)			iCut[10][2] =1;


	      if( mDecayLengthLB > lambdaDecaylengthDef_1  )		iCut[11][0] =1;
	      if( mDecayLengthLB > lambdaDecaylengthMin_1  )		iCut[11][1] =1;
	      if( mDecayLengthLB > lambdaDecaylengthMax_1  )		iCut[11][2] =1;


	      if( mDcaPToPVLB > lambdaDcaPtoPvDef_1)			iCut[12][0] =1;
	      if( mDcaPToPVLB > lambdaDcaPtoPvMin_1)			iCut[12][1] =1;
	      if( mDcaPToPVLB > lambdaDcaPtoPvMax_1)			iCut[12][2] =1;


	      if( mDcaPiToPVLB > lambdaDcaPitoPvDef_1  ) 			iCut[13][0] =1;
	      if( mDcaPiToPVLB > lambdaDcaPitoPvMin_1  ) 			iCut[13][1] =1;
	      if( mDcaPiToPVLB > lambdaDcaPitoPvMax_1  ) 			iCut[13][2] =1;


	      break;
	    }
    case 2: {
	      if( mDcaV0ToPVLB < lambdaDcaV0PvDef_2)			iCut[10][0] =1;
	      if( mDcaV0ToPVLB < lambdaDcaV0PvMin_2)			iCut[10][1] =1;
	      if( mDcaV0ToPVLB < lambdaDcaV0PvMax_2)			iCut[10][2] =1;


	      if( mDecayLengthLB > lambdaDecaylengthDef_2)		iCut[11][0] =1;
	      if( mDecayLengthLB > lambdaDecaylengthMin_2)		iCut[11][1] =1;
	      if( mDecayLengthLB > lambdaDecaylengthMax_2)		iCut[11][2] =1;


	      if( mDcaPToPVLB > lambdaDcaPtoPvDef_2)			iCut[12][0] =1;
	      if( mDcaPToPVLB > lambdaDcaPtoPvMin_2)			iCut[12][1] =1;
	      if( mDcaPToPVLB > lambdaDcaPtoPvMax_2)			iCut[12][2] =1;


	      if( mDcaPiToPVLB > lambdaDcaPitoPvDef_2) 			iCut[13][0] =1;
	      if( mDcaPiToPVLB > lambdaDcaPitoPvMin_2) 		iCut[13][1] =1;
	      if( mDcaPiToPVLB > lambdaDcaPitoPvMax_2) 		iCut[13][2] =1;


	      break;
	    }
    case 3: {
	      if( mDcaV0ToPVLB < lambdaDcaV0PvDef_3)			iCut[10][0] =1;
	      if( mDcaV0ToPVLB < lambdaDcaV0PvMin_3)			iCut[10][1] =1;
	      if( mDcaV0ToPVLB < lambdaDcaV0PvMax_3)			iCut[10][2] =1;


	      if( mDecayLengthLB > lambdaDecaylengthDef_3  )		iCut[11][0] =1;
	      if( mDecayLengthLB > lambdaDecaylengthMin_3  )		iCut[11][1] =1;
	      if( mDecayLengthLB > lambdaDecaylengthMax_3  )		iCut[11][2] =1;


	      if( mDcaPToPVLB > lambdaDcaPtoPvDef_3)			iCut[12][0] =1;
	      if( mDcaPToPVLB > lambdaDcaPtoPvMin_3)			iCut[12][1] =1;
	      if( mDcaPToPVLB > lambdaDcaPtoPvMax_3)			iCut[12][2] =1;


	      if( mDcaPiToPVLB > lambdaDcaPitoPvDef_3  ) 		iCut[13][0] =1;
	      if( mDcaPiToPVLB > lambdaDcaPitoPvMin_3  ) 		iCut[13][1] =1;
	      if( mDcaPiToPVLB > lambdaDcaPitoPvMax_3  ) 		iCut[13][2] =1;


	      break;
	    }

    }


    if( mDcaDaughtersLB < lambdaDcaDaughtersDef)				iCut[14][0] =1;
    if( mDcaDaughtersLB < lambdaDcaDaughtersMin)				iCut[14][1] =1;
    if( mDcaDaughtersLB < lambdaDcaDaughtersMax)				iCut[14][2] =1;

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
