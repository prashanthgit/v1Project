#define proton_cxx
#include "proton.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "Constants.h"

void proton::Loop() {

  fOutputfile = new TFile(mOutputFileName,"recreate");
  TProfile *ppp[nCuts][9];
  char buf[100];
  for(int i=0;i<nCuts;i++){
    for(int c=0;c<9;c++){
      sprintf(buf,"%s_%s_%d_%d",histname[0].c_str(),energyname.c_str(),i,c);
      ppp[i][c] = new TProfile(buf,"",20,-1,1);
      ppp[i][c] -> Sumw2();
    }
  }
  TH1I* hFlag = new TH1I("hFlag","hFlag",16,0,16);



  if (fChain == 0) return;


  Long64_t nentries = fChain->GetEntriesFast();
  cout << nentries <<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    ////////////// ////////////// ////////////// //////////////
    ////////////// ////////////// ////////////// //////////////

    Bool_t iCut[7][3]={{0}};

    // Event 
    if(TMath::Abs(mVzT) <= vzDef)							iCut[0][0] =1;
    if(TMath::Abs(mVzT) <= vzMin)							iCut[0][1] =1;
    if(TMath::Abs(mVzT) <= vzMax)							iCut[0][2] =1;

    if(mNHitsT >= nhitsDef) 								iCut[1][0] =1;
    if(mNHitsT >= nhitsMin) 								iCut[1][1] =1;
    if(mNHitsT >= nhitsMax) 								iCut[1][2] =1;

    if(mHitsPossT > nhitsposDef) 							iCut[2][0] =1;
    if(mHitsPossT > nhitsposMin)	 						iCut[2][1] =1;
    if(mHitsPossT > nhitsposMax)	 						iCut[2][2] =1;

    if(TMath::Abs(mNSigmaT) <= nsigmaDef) 						iCut[3][0] =1;
    if(TMath::Abs(mNSigmaT) <= nsigmaMin) 						iCut[3][1] =1;
    if(TMath::Abs(mNSigmaT) <= nsigmaMax) 						iCut[3][2] =1;

    if( (mass2protonMinDef <= mMass2T) && ( mMass2T <= mass2protonMaxDef) ) 		iCut[4][0] =1;
    if( (mass2protonMinMin <= mMass2T) && ( mMass2T <= mass2protonMaxMin) ) 		iCut[4][1] =1;
    if( (mass2protonMinMax <= mMass2T) && ( mMass2T <= mass2protonMaxMax) ) 		iCut[4][2] =1;

    if( mDcaT <= dcaDef)								iCut[5][0] =1; 
    if( mDcaT <= dcaMin)								iCut[5][1] =1; 
    if( mDcaT <= dcaMax)								iCut[5][2] =1; 


    // Momentum
    if(mMomentumT <= protonPMaxDef && protonPtMinDef <= mPtT && mPtT <= protonPtMaxDef) 	iCut[6][0] =1;
    if(mMomentumT <= protonPMaxMin && protonPtMinMin <= mPtT && mPtT <= protonPtMaxMin) 	iCut[6][1] =1;
    if(mMomentumT <= protonPMaxMax && protonPtMinMax <= mPtT && mPtT <= protonPtMaxMax) 	iCut[6][2] =1;



    ////////////// ////////////// ////////////// //////////////
    ////////////// ////////////// ////////////// //////////////

    mFlag = 0;
    flag[mFlag]=1;
    for(int i=0;i<7;i++) flag[mFlag] =  flag[mFlag] * iCut[i][0]; 
    if(flag[mFlag]) {ppp[mFlag][mCentralityBin] -> Fill(mRapidityT,mV1Track); hFlag->Fill(mFlag);}
    mFlag++;

    for(int i=0;i<7;i++){
      for(int j=1;j<=2;j++){
	flag[mFlag] = iCut[i][j];
	if(!flag[mFlag]){mFlag++;continue;}
	for(int k=0;k<7;k++){
	  if(i==k)continue;
	  flag[mFlag] = flag[mFlag] * iCut[k][0];
	}
	if(flag[mFlag]) {ppp[mFlag][mCentralityBin]-> Fill(mRapidityT,mV1Track); hFlag->Fill(mFlag);}
	mFlag++;
      }
    }


    ////////////// ////////////// ////////////// //////////////
    ////////////// ////////////// ////////////// //////////////

    if(( jentry %10000) == 0 ) cout << jentry  <<endl;
  }


  for(int i=0;i<nCuts;i++){
    for(int c=0;c<9;c++){
      ppp[i][c] -> Write(); 
    }
  }
  hFlag->Write();

  fOutputfile-> Close();


}
