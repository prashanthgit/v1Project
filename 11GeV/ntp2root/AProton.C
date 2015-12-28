
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "Constants.h"

TChain tree("aproton");
const static Int_t nCuts = 16;
Bool_t flag[nCuts];
Int_t mFlag;
TString mOutputFileName;

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



void AProton(Char_t path_file[100],Char_t jobid[100] ){

  mOutputFileName = Form("%s_%s.root",histname[1].c_str(),jobid);

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



  tree.SetBranchAddress("mVzT", &mVzT, &b_mVzT);
  tree.SetBranchAddress("mMomentumT", &mMomentumT, &b_mMomentumT);
  tree.SetBranchAddress("mPtT", &mPtT, &b_mPtT);
  tree.SetBranchAddress("mNSigmaT", &mNSigmaT, &b_mNSigmaT);
  tree.SetBranchAddress("mMass2T", &mMass2T, &b_mMass2T);
  tree.SetBranchAddress("mDcaT", &mDcaT, &b_mDcaT);
  tree.SetBranchAddress("mNHitsT", &mNHitsT, &b_mNHitsT);
  tree.SetBranchAddress("mHitsPossT", &mHitsPossT, &b_mHitsPossT);
  tree.SetBranchAddress("mRapidityT", &mRapidityT, &b_mRapidityT);
  tree.SetBranchAddress("mPhiPsiT", &mPhiPsiT, &b_mPhiPsiT);
  tree.SetBranchAddress("mPhiT", &mPhiT, &b_mPhiT);
  tree.SetBranchAddress("mEtaT", &mEtaT, &b_mEtaT);
  tree.SetBranchAddress("mCentralityBin", &mCentralityBin, &b_mCentralityBin);
  tree.SetBranchAddress("mV1Track", &mV1Track, &b_mV1Track);






  fOutputfile = new TFile(mOutputFileName,"recreate");
  TProfile *ppp[nCuts][9];
  char buf[100];
  for(int i=0;i<nCuts;i++){
    for(int c=0;c<9;c++){
      sprintf(buf,"%s_%s_%d_%d",histname[1].c_str(),energyname.c_str(),i,c);
      ppp[i][c] = new TProfile(buf,"",20,-1,1);
      ppp[i][c] -> Sumw2();
    }
  }
  TH1I* hFlag = new TH1I("hFlag","hFlag",16,0,16);





  Long64_t nentries = tree.GetEntriesFast();
  cout << nentries <<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = tree.GetEntry(jentry);   nbytes += nb;

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


