#define kshort_cxx
#include "kshort.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "Constants.h"

void kshort::Loop() {

  fOutputfile = new TFile(mOutputFileName,"recreate");
  TTree* hhh = new TTree("v0","RECREATE");
  hhh -> Branch("mFlag",&mFlag,"mFlag/s");
  hhh -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");
  hhh -> Branch("mRapidityBin",&mRapidityBin,"mRapidityBin/s");
  hhh -> Branch("mPhiPsiBin",&mPhiPsiBin,"mPhiPsiBin/s");
  hhh -> Branch("mKshortMass",&mKshortMass,"mKshortMass/D");




  if (fChain == 0) return;


  Long64_t nentries = fChain->GetEntriesFast();
  cout << nentries <<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    mRapidityBin = RapidityBin(mRapidityK);
    mPhiPsiBin = PhiPsiBin(mPhiPsiK);
    //if( mRapidityBin == 999 || mPhiPsiBin == 999 ) continue;
    ////////////// ////////////// ////////////// //////////////
    ////////////// ////////////// ////////////// //////////////
//
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
    if( (mMomentumK <= kshortMomMaxBar) && (kshortPtMinBar <= mPtK) && (mPtK <= kshortPtMaxBar))iCut[4][2] =2;

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
