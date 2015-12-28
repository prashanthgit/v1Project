
{
#include <iostream>
#include <fstream>
#include <cstdio>
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


  //TTree *tree;
  TChain tree("proton");
  std::ifstream inputFile("test.list");
  std::string line;
  while(getline(inputFile, line)) {
    cout << line <<endl;
    TFile *f = new TFile(line.c_str()); 
    tree->Add(line.c_str());
    //f->GetObject("proton",tree);
  }

   if (!tree) return;
  Long64_t nentries = tree->GetEntries();
  cout << nentries <<endl;


  tree->SetBranchAddress("mVzT", &mVzT, &b_mVzT);
  tree->SetBranchAddress("mMomentumT", &mMomentumT, &b_mMomentumT);
  tree->SetBranchAddress("mPtT", &mPtT, &b_mPtT);
  tree->SetBranchAddress("mNSigmaT", &mNSigmaT, &b_mNSigmaT);
  tree->SetBranchAddress("mMass2T", &mMass2T, &b_mMass2T);
  tree->SetBranchAddress("mDcaT", &mDcaT, &b_mDcaT);
  tree->SetBranchAddress("mNHitsT", &mNHitsT, &b_mNHitsT);
  tree->SetBranchAddress("mHitsPossT", &mHitsPossT, &b_mHitsPossT);
  tree->SetBranchAddress("mRapidityT", &mRapidityT, &b_mRapidityT);
  tree->SetBranchAddress("mPhiPsiT", &mPhiPsiT, &b_mPhiPsiT);
  tree->SetBranchAddress("mPhiT", &mPhiT, &b_mPhiT);
  tree->SetBranchAddress("mEtaT", &mEtaT, &b_mEtaT);
  tree->SetBranchAddress("mCentralityBin", &mCentralityBin, &b_mCentralityBin);
  tree->SetBranchAddress("mV1Track", &mV1Track, &b_mV1Track);


  Long64_t nentries = tree->GetEntriesFast();
  cout << nentries <<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
   // Long64_t ientry = LoadTree(jentry);
    //if (ientry < 0) break;
    nb = tree->GetEntry(jentry);   nbytes += nb;
    cout <<mCentralityBin << endl;
  }

}
