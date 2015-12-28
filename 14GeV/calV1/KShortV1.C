
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TProfile.h"
#include <iostream> 

using namespace std; 

////// EP reslution //////
//Double_t RPresolution[9] = {0.331726,0.439187,0.563428,0.634947,0.639296,0.586541,0.458289,0.303581,0.162335};        // 7
//Double_t RPresolution[9] = {0.303036,0.404264,0.515243,0.582435,0.599519,0.561734,0.455123,0.302329,0.161215};        // 11
Double_t RPresolution[9] = {0.196782, 0.267807, 0.338045, 0.39346, 0.413343, 0.396492, 0.32453, 0.225285, 0.120832}; // 14
//Double_t RPresolution[9] = {0.184772,0.260588,0.340869,0.404803,0.438205,0.43153,0.364991,0.254255,0.145082};         // 19
//Double_t RPresolution[9] = {0.092077,0.155365,0.225169,0.28477,0.326373,0.332316,0.290974,0.212434,0.11631};          // 27
//Double_t RPresolution[9] = {0.0305477,0.032294,0.0779149,0.139483,0.187975,0.213022,0.201958,0.155992,0.0976138};     // 39


// Set root files and spectra files
TFile* fS  = new TFile("/star/u/sprastar/v1_u8_res/V0/Kshort.15GeV.root");
TFile* fB  = new TFile("/star/u/sprastar/v1_u8_res/V0/KshortBG.15GeV.root");
TFile* fOut = new TFile("pV1K.15GeV.root","recreate");

char  v1v1[50] =    "pV1K_15GeV";

TProfile* pV1[31][9];

// main variable
Double_t rap;
Double_t phiRp;
Double_t v1Obs;
Double_t v1;

char buf[255];

int bg11 = 10;
int bg12 = 35;
int bg21 = 60;
int bg22 = 85;


void KShortV1(){
  TH1::AddDirectory(false);
  TH2::AddDirectory(false);

  Init();
  v1();
  Finish();
}





void v1(){
  TH1::AddDirectory(false);
  TH2::AddDirectory(false);
  int count =1;

  for(int f=0;f<31;f++){
    for(int c=0;c<9;c++){
      for(int y=0;y<20;y++){
	for(int p=0;p<30;p++){

	  // v1 calculation
	  rap   = (y-9.5)/10;
	  phiRp = (2 *p + 1) * ( 6.29 / (2.0 * 30.0 ));
	  v1Obs = cos(phiRp);
	  v1 = v1Obs / RPresolution[c];

	  // Yield calculation
	  sprintf(buf,"h_%d_%d_%d_%d",f,c,y,p);
	  TH2D* hSMPt = (TH2D*)fS->Get(buf);
	  TH2D* hBMPt = (TH2D*)fB->Get(buf);

	  TH1D* hS = hSMPt->ProjectionX("hS",3,16,"");
	  TH1D* hB = hBMPt->ProjectionX("hB",3,16,"");

	  Int_t sigCount1 = hS -> Integral(bg11,bg12);
	  Int_t bgCount1  = hB -> Integral(bg11,bg12);
	  Int_t sigCount2 = hS -> Integral(bg21,bg22);
	  Int_t bgCount2  = hB -> Integral(bg21,bg22);
	  Int_t sigCount  = (sigCount1 + sigCount2) / 2.0; 
	  Int_t bgCount   = (bgCount1  + bgCount2) / 2.0; 
	  Float_t scal = (bgCount==0)?1:(Float_t)sigCount/bgCount;
	  hB->Scale(scal);

	  TH1D* hhh = new TH1D("hhh","",100,0.38,0.62);
	  hhh->Add(hS,hB,1,-1);

	  int meanBin = hhh->GetMaximumBin();
	  int minBin  = meanBin - 3;
	  int maxBin  = meanBin + 3;
	  Double_t  rawYield = floor( hhh->Integral(minBin,maxBin));


	  cout << count <<endl;

	  // Fill
	  if(  (rawYield > 0) && (45<=meanBin && meanBin <= 53  )){
	    for(int zz=0; zz < rawYield; zz++)  pV1[f][c] -> Fill(rap,v1);
	  }

	  // Delete memory
	  hSMPt -> Delete();
	  hBMPt -> Delete();
	  hS -> Delete();
	  hB -> Delete();
	  hhh -> Delete();

	  count++;

	}
      }
    }
  }


} 


void Init(){
  for(int f=0;f<31;f++){
    for(int c=0;c<9;c++){
      sprintf(buf,"%s_%d_%d",v1v1,f,c);
      pV1[f][c] = new TProfile(buf,"",10,-1,1);
    }
  }
}

void Finish(){
  fOut->cd();
  for(int f=0;f<31;f++){
    for(int i=0;i<9;i++){
      pV1[f][i] -> Write();
    }
  }
  fOut -> Close();
}

