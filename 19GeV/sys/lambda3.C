

TF1 *PrevFitTMP = new TF1("PrevFitTMP","pol1",-0.8,0.8);
PrevFitTMP->SetParameter(0,0);
PrevFitTMP->SetParError(0,1);


TProfile* pV1[31][9];
double  yData[31][9];
double eYData[31][9];
char buf[255];

TGraphErrors* gFlag;
TGraphErrors* gSys;
double xx[9]={0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5};
double exx[9] = {0};
double yy[9] = {0};
double sys[9];


void lambda3(){

  GetPlots();
  Fit();
  GetSys();
  Plot1();
}


void Plot1(){
  TCanvas *c = new TCanvas("c", "",216,211,1200,600);
  gStyle->SetOptTitle(0);
  c->Range(-21.66667,-0.2083333,61.66667,0.2083333);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetLeftMargin(0.2);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.3);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);


  TH1F *hframe__1 = new TH1F("hframe__1","",9,0,9);
  hframe__1->SetMinimum(-0.15);
  hframe__1->SetMaximum(0.1);
  hframe__1->SetDirectory(0);
  hframe__1->SetStats(0);
  hframe__1->GetXaxis()->CenterTitle(true);
  hframe__1->GetXaxis()->SetNdivisions(-10);
  hframe__1->GetXaxis()->SetLabelFont(42);
  hframe__1->GetXaxis()->SetLabelOffset(0.009);
  hframe__1->GetXaxis()->SetTickLength(0.02);
  hframe__1->GetXaxis()->SetTitleOffset(1.8);
  hframe__1->GetXaxis()->SetTitleSize(0.08);
  hframe__1->GetXaxis()->SetLabelSize(0.08);  
  hframe__1->GetXaxis()->SetLabelColor(1);  
  hframe__1->GetXaxis()->SetTitle("centrality"); 
  hframe__1->GetYaxis()->SetTitle("d#font[12]{v}_{1} / d#font[12]{y}"); 
  hframe__1->GetYaxis()->CenterTitle(true); 
  hframe__1->GetYaxis()->SetNdivisions(705);  
  hframe__1->GetYaxis()->SetLabelFont(42);  
  hframe__1->GetYaxis()->SetLabelSize(0.05); 
  hframe__1->GetYaxis()->SetLabelOffset(0.009);  
  hframe__1->GetYaxis()->SetTitleSize(0.08);
  hframe__1->GetYaxis()->SetTitleOffset(0.7); 
  hframe__1->GetYaxis()->SetTickLength(0.02);  
  hframe__1->GetYaxis()->SetTitleFont(42);  
  hframe__1->GetXaxis()->SetBit(TAxis::kLabelsVert); 
  hframe__1->GetXaxis()->SetBinLabel(1,"80-70%");
  hframe__1->GetXaxis()->SetBinLabel(2,"70-60%");
  hframe__1->GetXaxis()->SetBinLabel(3,"60-50%");
  hframe__1->GetXaxis()->SetBinLabel(4,"50-40%");
  hframe__1->GetXaxis()->SetBinLabel(5,"40-30%");
  hframe__1->GetXaxis()->SetBinLabel(6,"30-20%");
  hframe__1->GetXaxis()->SetBinLabel(7,"20-10%");
  hframe__1->GetXaxis()->SetBinLabel(8,"10-5%");
  hframe__1->GetXaxis()->SetBinLabel(9,"5-0%");
  hframe__1->Draw(" ");

  TLine *line = new TLine(0,0,9,0);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();

  gFlag = new TGraphErrors(9,xx,yData[0],exx,eYData[0]);
  gFlag->SetMarkerSize(1.4);
  gFlag->SetMarkerStyle(20);
  gFlag->Draw("Psame9");

  gSys = new TGraphErrors(9,xx,yy,exx,sys);
  gSys->SetFillColor(1);
  gSys->SetFillStyle(3003);
  gSys->Draw("E3same9");


 }

void GetSys(){

  double rms2[9]={0};
  for(int c=0;c<9;c++){
    for(int f=1;f<31;f++){
      if(f==19 || f==20) continue;
      if(yData[f][c] == -999){
	rms2[c] += ( (yData[0][c]-yData[f-1][c]) * (yData[0][c]-yData[f-1][c]) );
	continue;
      }
      rms2[c] += ( (yData[0][c]-yData[f][c]) * (yData[0][c]-yData[f][c]) );
    }
    sys[c] =  TMath::Sqrt(rms2[c]/2);
  }
}


void Fit(){
  TCanvas* forFit = new TCanvas();
  for(int c=0;c<9;c++){
    for(int f=0;f<31;f++){
      yData[f][c] = -999;
      eYData[f][c] = -999;
      pV1[f][c] -> Fit(PrevFitTMP,"EMR"); 
      if(pV1[f][c]->GetEntries() != 0)   yData[f][c] = PrevFitTMP->GetParameter(1);
      if(pV1[f][c]->GetEntries() != 0)   eYData[f][c] = PrevFitTMP->GetParError(1);
    }
  }
  forFit->Close();
}


void GetPlots(){
  TFile* dataRoot = new TFile("./res/v1.11GeV.root");
  for(int f=0;f<31;f++){
    for(int c=0;c<9;c++){
      sprintf(buf,"%s_%d_%d","pV1L_11",f,c);
      pV1[f][c] = (TProfile*)dataRoot->Get(buf);
    }
  }

}
