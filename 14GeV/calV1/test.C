
{
  int int1 = 29;
  int int2 = 34;
  int bg1  = 40;
  int bg2  = 80;


  TFile* fSL   = new TFile("/star/u/sprastar/v1_u8_res/V0/Lambda.19GeV.root");
  TFile* fBL   = new TFile("/star/u/sprastar/v1_u8_res/V0/LambdaBG.19GeV.root");

  char buf[100];
  for(int i=10;i<11;i++){
  sprintf(buf,"h_%d_%d_%d_%d",0,4,8,i);
  TH2D* hSMPt = (TH2D*)fSL->Get(buf);
  TH2D* hBMPt = (TH2D*)fBL->Get(buf);
  TH1D* hS = hSMPt->ProjectionX("hS",3,16,"");
  TH1D* hB = hBMPt->ProjectionX("hB",3,16,"");


  Int_t sigCount = hS ->Integral(bg1,bg2);
  Int_t bgCount =  hB ->Integral(bg1,bg2);
  Float_t scal = (bgCount==0)?1:(Float_t)sigCount/bgCount;

  hB->Scale(scal);
  Double_t  rawYield = floor( hS->Integral(int1,int2) - hB->Integral(int1,int2) );

  cout << rawYield <<endl;
  TH1D* hhh = new TH1D("hhh","",100,1.06,1.24);
  hhh->Add(hS,hB,1,-1);


  //hSMPt->Draw();
  new TCanvas();
  hhh->Draw();
  //hB->Draw("SAME");
  }
}
