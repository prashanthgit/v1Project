
{

  TFile* fS  = new TFile("/star/u/sprastar/v1_u8_res/V0/Kshort.19GeV.root");
  TFile* fB  = new TFile("/star/u/sprastar/v1_u8_res/V0/KshortBG.19GeV.root");

  char buf[255];

  int bg11 = 10;
  int bg12 = 35;
  int bg21 = 60;
  int bg22 = 85;

  // Yield calculation
  sprintf(buf,"h_%d_%d_%d_%d",0,4,9,10);
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


  cout <<meanBin <<"\t"<<minBin << "\t"<< maxBin <<endl;


}
