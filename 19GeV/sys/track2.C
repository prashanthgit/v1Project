
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

TGraphErrors* gSys[6];
double yy[10]={0};
double eyy[10];
double xx[10] = {-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9};
double exx[10] = {0};


void track2(){

  GetPlots();
  Addplots();
  GetSys();
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
  c->SetLeftMargin(0.2);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.1);
  c->SetBottomMargin(0.2);
  c->SetFrameBorderMode(0);
  c->SetFrameBorderMode(0);
  c->Divide(2,3,0,0);
  for(int p=0;p<6;p++){
    c->cd(p+1);


    TH1F *hframe__1 = new TH1F("hframe__1","",1200,-1,1);
    hframe__1->SetMinimum(-0.03);
    hframe__1->SetMaximum(0.03);
    hframe__1->SetDirectory(0);
    hframe__1->SetStats(0);
    hframe__1->SetLineColor(1);
    hframe__1->GetXaxis()->SetTitle("");
    hframe__1->GetXaxis()->CenterTitle(true);
    hframe__1->GetXaxis()->SetMoreLogLabels();
    hframe__1->GetXaxis()->SetNoExponent();
    hframe__1->GetXaxis()->SetNdivisions(010);
    hframe__1->GetXaxis()->SetLabelFont(42);
    hframe__1->GetXaxis()->SetLabelOffset(0.007);
    hframe__1->GetXaxis()->SetLabelSize(0.08);
    hframe__1->GetXaxis()->SetTitleSize(0.1);
    hframe__1->GetXaxis()->SetTickLength(0.02);
    hframe__1->GetXaxis()->SetTitleOffset(1.0);
    hframe__1->GetXaxis()->SetTitleFont(42);
    hframe__1->GetXaxis()->SetTitle("y");
    hframe__1->GetYaxis()->SetTitle("#font[12]{v}_{1}");
    hframe__1->GetYaxis()->CenterTitle(true);
    hframe__1->GetYaxis()->SetNdivisions(406);
    hframe__1->GetYaxis()->SetLabelFont(42);
    hframe__1->GetYaxis()->SetLabelOffset(0.006);
    hframe__1->GetYaxis()->SetLabelSize(0.08);
    hframe__1->GetYaxis()->SetTitleSize(0.15);
    hframe__1->GetYaxis()->SetTickLength(0.02);
    hframe__1->GetYaxis()->SetTitleOffset(0.5);
    hframe__1->GetYaxis()->SetTitleFont(42);
    hframe__1->Draw(" ");

    TLine *line = new TLine(-1,0,1,0);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();

    TLine *line = new TLine(0,-0.03,0,0.03);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw();

    gStyle->SetErrorX(0.0001);
    pV1Data[p][0] -> SetMarkerStyle(20);
    pV1Data[p][0] -> SetMarkerColor(2);
    pV1Data[p][0] -> SetLineColor(2);
    pV1Data[p][0] -> SetMarkerSize(1.1);
    pV1Data[p][0] -> Draw("same");

    gSys[p]->SetFillColor(1);
    gSys[p]->SetFillStyle(3003);
    gSys[p]->Draw("E3same9");
  }

}


void GetSys(){


  for(int p=0;p<6;p++){
    for(int b=1;b<11;b++){
      double rms = 0;
      for(int f=1;f<13;f++){
	rms += ( ( pV1Data[p][0]->GetBinContent(b) -  pV1Data[p][f]->GetBinContent(b) ) * ( pV1Data[p][0]->GetBinContent(b) -  pV1Data[p][f]->GetBinContent(b) ) ); 
      }
      eyy[b-1] = TMath::Sqrt(rms/2);
    }
    gSys[p] = new TGraphErrors(10,xx,yy,exx,eyy);
  }


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
