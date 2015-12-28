
int c1_low = 4;
int c1_hi  = 6;


TF1 *PrevFitTMP = new TF1("PrevFitTMP","pol1",-0.8,0.8);
PrevFitTMP->SetParameter(0,0);
PrevFitTMP->SetParError(0,1);


TProfile* pV1[6][15][9];
double  yData[6][15][9];
double eYData[6][15][9];
char buf[255];

TGraphErrors* gFlag[6];

TGraphErrors* gSys[6];
double xx[9]={0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5};
double exx[9] = {0};
double yy[9] = {0};
double sys[6][9]={{0}};


void track3(){

  GetPlots();
  Fit();
  GetSys();
  new TCanvas();
  //for(int i=0;i<15;i++){
  //  if(i==9 || i==10)continue;
  //  pV1[0][i][3]->Draw("same");
  //}
  Plot();
}


void Plot(){
  TCanvas *c = new TCanvas("c", "",216,211,1000,800);
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
  c->SetBottomMargin(0.3);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  c->Divide(1,3,0,0);
  int p =0;
  for(int i=0;i<3;i++){
    if(i==0)p=0;
    if(i==1)p=2;
    if(i==2)p=3;
    c->cd(i+1);


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
  hframe__1->GetXaxis()->SetTitleSize(0.1);
  hframe__1->GetXaxis()->SetLabelSize(0.1);  
  hframe__1->GetXaxis()->SetLabelColor(1);  
  hframe__1->GetXaxis()->SetTitle("centrality"); 
  hframe__1->GetYaxis()->SetTitle("d#font[12]{v}_{1} / d#font[12]{y}"); 
  hframe__1->GetYaxis()->CenterTitle(true); 
  hframe__1->GetYaxis()->SetNdivisions(705);  
  hframe__1->GetYaxis()->SetLabelFont(42);  
  hframe__1->GetYaxis()->SetLabelSize(0.08); 
  hframe__1->GetYaxis()->SetLabelOffset(0.009);  
  hframe__1->GetYaxis()->SetTitleSize(0.1);
  hframe__1->GetYaxis()->SetTitleOffset(0.5); 
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

    gFlag[p] = new TGraphErrors(9,xx,yData[p][0],exx,eYData[p][0]);
    gFlag[p]->SetMarkerSize(1.4);
    gFlag[p]->SetMarkerStyle(20);
    gFlag[p]->Draw("Psame9");

    gSys[p] = new TGraphErrors(9,xx,yy,exx,sys[p]);
    gSys[p]->SetFillColor(1);
    gSys[p]->SetFillStyle(3003);
    gSys[p]->Draw("E3same9");
  }

}


void GetSys(){

  for(int p=0;p<6;p++){
    double rms2[9] = {0};
    for(int c=0;c<9;c++){
      for(int f=1;f<13;f++){
	rms2[c] += ( (yData[p][0][c]-yData[p][f][c]) * (yData[p][0][c]-yData[p][f][c]) );
      }
      sys[p][c] = TMath::Sqrt(rms2[c]/2.0);
    }
  }
}



void Fit(){

  TCanvas* forFit = new TCanvas();
  for(int p=0;p<6;p++){
    for(int c=0;c<9;c++){
      for(int f=0;f<15;f++){
	pV1[p][f][c] -> Fit(PrevFitTMP,"EMR"); 
	if(pV1[p][f][c]->GetEntries() != 0)  yData[p][f][c] = PrevFitTMP->GetParameter(1);
	if(pV1[p][f][c]->GetEntries() != 0)  eYData[p][f][c] = PrevFitTMP->GetParError(1);

      }
    }
  }
  forFit->Close();

}


void GetPlots(){
  TFile* dataRoot =  new TFile("./res/v1.11GeV.root");

  for(int f=0;f<15;f++){
    for(int c=0;c<9;c++){
      // Proton
      sprintf(buf,"pV1P_11_%d_%d",f,c);
      pV1[0][f][c] = (TProfile*)dataRoot->Get(buf);

      // Anti-Proton
      sprintf(buf,"pV1AP_11_%d_%d",f,c);
      pV1[1][f][c] = (TProfile*)dataRoot->Get(buf);

      // Pion Pos 
      sprintf(buf,"pV1PiP_11_%d_%d",f,c);
      pV1[2][f][c] = (TProfile*)dataRoot->Get(buf);

      //Pion Neg
      sprintf(buf,"pV1PiN_11_%d_%d",f,c);
      pV1[3][f][c] = (TProfile*)dataRoot->Get(buf);

      // Kaon Pos
      sprintf(buf,"pV1KP_11_%d_%d",f,c);
      pV1[4][f][c] = (TProfile*)dataRoot->Get(buf);

      //Kaon Neg
      sprintf(buf,"pV1KN_11_%d_%d",f,c);
      pV1[5][f][c] = (TProfile*)dataRoot->Get(buf);

    }
  }


}
