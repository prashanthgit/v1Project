
int c1_low = 4;
int c1_hi  = 6;


TF1 *PrevFitTMP = new TF1("PrevFitTMP","pol1",-0.8,0.8);
PrevFitTMP->SetParameter(0,0);
PrevFitTMP->SetParError(0,1);


TProfile* pV1[6][15][9];
TProfile* pV1Data[6][15];
double yData[6][15];
double eYData[6][15];
char buf[255];

TGraphErrors* gFlag[6];
TGraphErrors* gPoint1[6];
TGraphErrors* gPoint2[6];
double xx[15];
double exx[15] = {0};


void track(){

  GetPlots();
  Addplots();
  Fit();
  GetSys();
  Plot();
}


void Plot(){
  TCanvas *c = new TCanvas("c", "",216,211,1200,800);
  gStyle->SetOptTitle(0);
  c->Range(-21.66667,-0.2083333,61.66667,0.2083333);
  c->SetFillColor(0);
  c->SetBorderMode(0);
  c->SetBorderSize(2);
  c->SetTickx(1);
  c->SetTicky(1);
  c->SetLeftMargin(0.5);
  c->SetRightMargin(0.0);
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.2);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  c->Divide(2,3,0.01,0.01);
  c->SetLeftMargin(0.2);
  c_1->SetLeftMargin(0.2);
  c_2->SetLeftMargin(0.2);
  c_3->SetLeftMargin(0.2);
  c_4->SetLeftMargin(0.2);
  c_5->SetLeftMargin(0.2);
  c_6->SetLeftMargin(0.2);
  c_1->SetRightMargin(0.01);
  c_2->SetRightMargin(0.01);
  c_3->SetRightMargin(0.01);
  c_4->SetRightMargin(0.01);
  c_5->SetRightMargin(0.01);
  c_6->SetRightMargin(0.01);

  for(int p=0;p<6;p++){
    c->cd(p+1);


  TH1F *hframe__1 = new TH1F("hframe__1","",1200,-1,17);
  hframe__1->SetMinimum(-0.05);
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
  hframe__1->GetXaxis()->SetLabelSize(0.08);
  hframe__1->GetXaxis()->SetTitleSize(0.05);
  hframe__1->GetXaxis()->SetTickLength(0.02);
  hframe__1->GetXaxis()->SetTitleOffset(1.2);
  hframe__1->GetXaxis()->SetTitleFont(42);
  hframe__1->GetYaxis()->SetTitle("d#font[12]{v}_{1} / d#font[12]{y}");
  hframe__1->GetYaxis()->CenterTitle(true);
  hframe__1->GetYaxis()->SetNdivisions(406);
  hframe__1->GetYaxis()->SetLabelFont(42);
  hframe__1->GetYaxis()->SetLabelOffset(0.006);
  hframe__1->GetYaxis()->SetLabelSize(0.08);
  hframe__1->GetYaxis()->SetTitleSize(0.12);
  hframe__1->GetYaxis()->SetTickLength(0.02);
  hframe__1->GetYaxis()->SetTitleOffset(0.8);
  hframe__1->GetYaxis()->SetTitleFont(42);
  hframe__1->Draw(" ");

  TLine *line = new TLine(-1,0,17,0);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw();

  gFlag[p]->Draw("Psame9");
  gPoint1[p]->Draw("Psame9");
  gPoint2[p]->Draw("p[]same9");
  }

 }


void GetSys(){

  for(int p=0;p<6;p++){
    double xp[1]; xp[0]=15;
    double yp[1]; yp[0] = yData[p][0];
    double xpe[1]; xpe[0]=0;
    double ype[1]; ype[0] = eYData[p][0]; 

    ype[0]=eYData[p][0];
    gPoint1[p] = new TGraphErrors(1,xp,yp,xpe,ype);
    gPoint1[p]->SetMarkerStyle(20);
    gPoint1[p]->SetMarkerSize(1.2);
    gPoint1[p]->SetMarkerColor(2);

    double rms2 = 0;
    for(int f=1;f<13;f++){
      rms2 += ( (yData[p][0]-yData[p][f]) * (yData[p][0]-yData[p][f]) );
    }

    ype[0] =  TMath::Sqrt(rms2/2);
    gPoint2[p] = new TGraphErrors(1,xp,yp,xpe,ype);
  }



}



void Fit(){

  TCanvas* forFit = new TCanvas();
  for(int p=0;p<6;p++){
    for(int f=0;f<15;f++){
      pV1Data[p][f] -> Fit(PrevFitTMP,"EMR"); 
      if(pV1Data[p][f]->GetEntries() != 0)  yData[p][f] = PrevFitTMP->GetParameter(1);
      if(pV1Data[p][f]->GetEntries() != 0)  eYData[p][f] = PrevFitTMP->GetParError(1);
      xx[f] = f;

    }
    gFlag[p] = new TGraphErrors(15,xx,yData[p],exx,eYData[p]);
    gFlag[p]->SetMarkerStyle(20);
    gFlag[p]->SetMarkerSize(1.2);
  }
  forFit->Close();

}

void Addplots(){

  char buf[255];
  for(int p=0;p<6;p++){
    for(int f=0;f<15;f++){
      sprintf(buf,"%d_%d",p,f);
      pV1Data[p][f] = new TProfile(buf,buf,10,-1,1);
      for(int c= c1_low; c<=c1_hi;c++){
	pV1Data[p][f] ->  Add(pV1[p][f][c]);
      }
    }
  }
}

void GetPlots(){
  TFile* dataRoot =  new TFile("./v1.19GeV.root");

  for(int f=0;f<15;f++){
    for(int c=0;c<9;c++){
      // Proton
      sprintf(buf,"pV1P_19GeV_%d_%d",f,c);
      pV1[0][f][c] = (TProfile*)dataRoot->Get(buf);

      // Anti-Proton
      sprintf(buf,"pV1AP_19GeV_%d_%d",f,c);
      pV1[1][f][c] = (TProfile*)dataRoot->Get(buf);

      // Pion Pos 
      sprintf(buf,"pV1PiP_19GeV_%d_%d",f,c);
      pV1[2][f][c] = (TProfile*)dataRoot->Get(buf);

      //Pion Neg
      sprintf(buf,"pV1PiN_19GeV_%d_%d",f,c);
      pV1[3][f][c] = (TProfile*)dataRoot->Get(buf);

      // Kaon Pos
      sprintf(buf,"pV1KP_19GeV_%d_%d",f,c);
      pV1[4][f][c] = (TProfile*)dataRoot->Get(buf);

      //Kaon Neg
      sprintf(buf,"pV1KN_19GeV_%d_%d",f,c);
      pV1[5][f][c] = (TProfile*)dataRoot->Get(buf);

    }
  }



  for(int f=0;f<15;f++){
    for(int c=0;c<9;c++){
      pV1[0][f][c] ->Rebin(2);
      pV1[1][f][c] ->Rebin(2);
      pV1[2][f][c] ->Rebin(2);
      pV1[3][f][c] ->Rebin(2);
      pV1[4][f][c] ->Rebin(2);
      pV1[5][f][c] ->Rebin(2);

    }
  }


}
