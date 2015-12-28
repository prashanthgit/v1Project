
int c1_low = 4;
int c1_hi  = 6;


TF1 *PrevFitTMP = new TF1("PrevFitTMP","pol1",-0.8,0.8);
PrevFitTMP->SetParameter(0,0);
PrevFitTMP->SetParError(0,1);


TProfile* pV1[19][9];
TProfile* pV1Data[19];
double yData[19];
double eYData[19];
char buf[255];

TGraphErrors* gFlag;
TGraphErrors* gPoint1;
TGraphErrors* gPoint2;
double xx[19];
double exx[19] = {0};
double sys;


void kshort(){

  GetPlots();
  Addplots();
  Fit();
  GetSys();
  Plot();
}


void Plot(){
  TCanvas *c = new TCanvas("c", "",216,211,1000,500);
  gStyle->SetOptTitle(0);
  c->Range(-21.66667,-0.2083333,61.66667,0.2083333);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetLeftMargin(0.1);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.1);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);


  TH1F *hframe__1 = new TH1F("hframe__1","",1200,-1,20);
  hframe__1->SetMinimum(-0.03);
  hframe__1->SetMaximum(0.03);
  hframe__1->SetDirectory(0);
  hframe__1->SetStats(0);
  hframe__1->SetLineColor(1);
  hframe__1->GetXaxis()->SetTitle("");
  hframe__1->GetXaxis()->CenterTitle(true);
  hframe__1->GetXaxis()->SetMoreLogLabels();
  hframe__1->GetXaxis()->SetNoExponent();
  hframe__1->GetXaxis()->SetNdivisions(040);
  hframe__1->GetXaxis()->SetLabelFont(42);
  hframe__1->GetXaxis()->SetLabelOffset(0.007);
  hframe__1->GetXaxis()->SetLabelSize(0.05);
  hframe__1->GetXaxis()->SetTitleSize(0.05);
  hframe__1->GetXaxis()->SetTickLength(0.02);
  hframe__1->GetXaxis()->SetTitleOffset(1.2);
  hframe__1->GetXaxis()->SetTitleFont(42);
  hframe__1->GetYaxis()->SetTitle("d#font[12]{v}_{1} / d#font[12]{y}");
  hframe__1->GetYaxis()->CenterTitle(true);
  hframe__1->GetYaxis()->SetNdivisions(406);
  hframe__1->GetYaxis()->SetLabelFont(42);
  hframe__1->GetYaxis()->SetLabelOffset(0.006);
  hframe__1->GetYaxis()->SetLabelSize(0.05);
  hframe__1->GetYaxis()->SetTitleSize(0.08);
  hframe__1->GetYaxis()->SetTickLength(0.02);
  hframe__1->GetYaxis()->SetTitleOffset(0.60);
  hframe__1->GetYaxis()->SetTitleFont(42);
  hframe__1->Draw(" ");

  TLine *line = new TLine(-1,0,20,0);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();

  gFlag->Draw("Psame9");
  gPoint1->Draw("Psame9");
  gPoint2->Draw("p[]same9");



 }

void GetSys(){

  double xp[1]; xp[0]=19;
  double yp[1]; yp[0] = yData[0];
  double xpe[1]; xpe[0]=0;
  double ype[1]; ype[0] = eYData[0]; 

  ype[0]=eYData[0];
  gPoint1 = new TGraphErrors(1,xp,yp,xpe,ype);
  gPoint1->SetMarkerStyle(20);
  gPoint1->SetMarkerSize(1.2);
  gPoint1->SetMarkerColor(2);

  double rms2=0;
  for(int f=1;f<19;f++){
    if(f==9 || f==10) continue;
    rms2 += ( (yData[0]-yData[f]) * (yData[0]-yData[f]) );
  }
  ype[0] =  TMath::Sqrt(rms2/2);
  gPoint2 = new TGraphErrors(1,xp,yp,xpe,ype);



 
}


void Fit(){

  TCanvas* forFit = new TCanvas();
  for(int f=0;f<19;f++){
    cout << f<<endl;
    yData[f] = -999;
    eYData[f] = -999;
    pV1Data[f] -> Fit(PrevFitTMP,"EMR"); 
    if(pV1Data[f]->GetEntries() != 0)   yData[f] = PrevFitTMP->GetParameter(1);
    if(pV1Data[f]->GetEntries() != 0)   eYData[f] = PrevFitTMP->GetParError(1);
    xx[f] = f;

  }
  forFit->Close();
  gFlag = new TGraphErrors(19,xx,yData,exx,eYData);
  gFlag->SetMarkerStyle(20);
  gFlag->SetMarkerSize(1.2);

}

void Addplots(){

  char buf[255];
  for(int f=0;f<19;f++){
    sprintf(buf,"%s_%d","p",f);
      pV1Data[f] = new TProfile(buf,buf,10,-1,1);
      for(int c= c1_low; c<=c1_hi;c++){
	pV1Data[f] ->  Add(pV1[f][c]);
      }
  }
}

void GetPlots(){
  TFile* dataRoot = new TFile("./res/v1.19GeV.root");
  for(int f=0;f<19;f++){
    for(int c=0;c<9;c++){
      sprintf(buf,"%s_%d_%d","pV1K_11",f,c);
      pV1[f][c] = (TProfile*)dataRoot->Get(buf);
    }
  }

}
