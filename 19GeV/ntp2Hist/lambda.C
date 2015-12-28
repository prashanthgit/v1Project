#define lambda_cxx
#include "lambda.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void lambda::Loop() {



  TH2D* hhh[31][9][20][30];
  char buf[255];
  for(int f=0;f<31;f++){
    for(int c=0;c<9;c++){
      for(int y=0;y<20;y++){
	for(int p=0;p<30;p++){
	  sprintf(buf,"h_%d_%d_%d_%d",f,c,y,p);
	  hhh[f][c][y][p] = new TH2D(buf,"",100,1.06,1.24,16,0,16);
	  //hhh[f][c][y][p] -> Sumw2 (); 
	}
      }
    }
  }


  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntries();
  cout << nentries <<endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    hhh[mFlag][mCentralityBin][mRapidityBin][mPhiPsiBin] -> Fill(mLambdaMass,mPtBin);
    if(( jentry %1000000) == 0 ) cout << jentry  <<endl;
  }


  char buff[255];
  sprintf(buff,"LambdaBG.19GeV.root");
  fOutputfile = new TFile(buff,"recreate");
  fOutputfile>cd();
  for(int f=0;f<31;f++){
    for(int c=0;c<9;c++){
      for(int y=0;y<20;y++){
	for(int p=0;p<30;p++){
	  hhh[f][c][y][p] -> Write(); 
	}
      }
    }
  }

  fOutputfile-> Close();

}

