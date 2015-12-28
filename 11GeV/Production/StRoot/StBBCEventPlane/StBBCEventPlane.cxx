

#include "TMath.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StBBCEventPlane.h"
#include "TVector2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom.h"
//#define pi = 3.14159;

ClassImp(StBBCEventPlane)


StBBCEventPlane :: ~StBBCEventPlane(){}


StBBCEventPlane :: StBBCEventPlane(StMuEvent* muEvent,int cent){

  //  //7.7 GeV
  //  float egain[16]={0.93224,0.999464,1.25065,0.969271,0.97567,0.872704,0.977678,1.12397,0.957429,1.12241,1.03693,0.916625,0.970594,0.899347,1.08898,1.01173};
  //  float wgain[16]={0.975746,1.11539,1.08283,0.862711,1.04224,0.92109,0.777219,0.896277,0.935381,1.34723,1.04392,1.01167,1.08745,1.00913,1.07499,1.02784};
  //
  //  //11.5 
    float egain[16]={0.93224,0.999464,1.25065,0.969271,0.97567,0.872704,0.977678,1.12397,0.957429,1.12241,1.03693,0.916625,0.970594,0.899347,1.08898,1.01173};
    float wgain[16]={0.975746,1.11539,1.08283,0.862711,1.04224,0.92109,0.777219,0.896277,0.935381,1.34723,1.04392,1.01167,1.08745,1.00913,1.07499,1.02784};
  //
  //  //14.5 GeV
  //  float egain[16]={1.4498,0.969584,1.05667,0.944818,0.734853,0.84428,0.769285,0.98697,0.988923,1.65317,1.0291,0.731249,1.07852,0.981391,1.34576,0.935099};
  //  float wgain[16]={1.05859,0.657341,0.677656,1.17072,1.22709,1.20861,0.731588,1.04389,0.849099,1.16959,0.76936,0.788733,1.02664,1.21351,1.78977,1.0975};
  //
  //  //19.6 GeV
  //float egain[16]={0.951405,1.01101,1.18373,0.940584,0.956778,0.956499,0.996344,1.02918,0.843086,1.09859,1.16699,0.870963,0.992364,0.982628,1.1668,0.985753};
  //float wgain[16]={1.04369,1.16231,1.17224,0.769896,1.0413,0.810565,0.870006,0.970389,0.840919,1.25982,1.00126,0.885717,1.07928,1.07985,1.22874,1.0283};
  //
  //  //27 GeV
  //  float egain[16]={1.02852,0.727557,1.24688,0.796138,1.10474,1.09617,1.14122,1.04763,1.19436,1.17578,1.08493,0.851029,0.823415,0.947453,0.839726,0.90221};
  //  float wgain[16]={0.917759,1.3602,1.15049,0.821606,0.916213,0.833731,0.803203,0.867922,0.751416,1.38092,1.073,0.93318,1.22186,1.00512,1.14631,1.08068};
  //
  //  //39 GeV
  //  float egain[16]={0.972341,0.957723,1.197,0.931167,0.988347,0.953421,0.994432,1.0734,0.960428,1.05807,1.06945,0.922059,0.970653,0.913712,1.10895,1.01235};
  //  float wgain[16]={1.02993,1.14555,1.05731,0.811009,1.05182,0.904392,0.800753,0.921134,0.872874,1.24529,0.998994,1.02481,1.10775,1.05473,1.10359,1.04452};

  mBBC = muEvent->bbcTriggerDetector();
  mCentralityID = cent;
  for(int i=0;i<6;i++) mPsi[i] = -999;
  float pi = TMath::Pi();

  // event loop

  float nHitE[16],nHitW[16];

  for(int i=0;i<16;i++){
    nHitE[i] = 1.0 * mBBC.adc(i)    / egain[i];
    nHitW[i] = 1.0 * mBBC.adc(i+24) / wgain[i];
  }


  //BBC East
  TVector2 mQE;
  Double_t mQEx=0., mQEy=0.;

  Float_t eXsum=0., eYsum=0., eWgt=0.,psi_e=0.,psi_e_s=0.;

  for(int iTile = 0; iTile < 16; iTile++) {

    eXsum += cos(BBC_GetPhi(0,iTile))*nHitE[iTile];
    eYsum += sin(BBC_GetPhi(0,iTile))*nHitE[iTile];
    eWgt += nHitE[iTile];
  }
  mQEx= (eWgt>0.0) ? eXsum/eWgt:0.0;
  mQEy= (eWgt>0.0) ? eYsum/eWgt:0.0;

  mQE.Set(mQEx,mQEy);

  if(mQE.Mod()){

    psi_e=mQE.Phi();

    if(psi_e<0.0) {psi_e +=2*pi;}  //done in BBC_GetPsi()
  }
  //	if(psi_e==0.0) return kFALSE; // done in BBC_GetPsi()

  //BBC West
  TVector2 mQW;
  Double_t mQWx=0., mQWy=0.;

  Float_t wXsum=0.,wYsum=0.,wWgt=0.,psi_w=0.,psi_w_s=0.;

  for(int iTile = 0; iTile < 16; iTile++) {

    wXsum += cos(BBC_GetPhi(1,iTile))*nHitW[iTile];
    wYsum += sin(BBC_GetPhi(1,iTile))*nHitW[iTile];
    wWgt += nHitW[iTile];

  }
  mQWx= (wWgt>0.0) ? wXsum/wWgt:0.0;
  mQWy= (wWgt>0.0) ? wYsum/wWgt:0.0;

  mQW.Set(mQWx,mQWy);

  if(mQW.Mod()){

    psi_w=mQW.Phi();

    if(psi_w<0.0) {psi_w +=2*pi;} //done in BBC_GetPsi()

  }
  //	if(psi_w==0.0) return kFALSE; //done in BBC_GetPsi()


  //BBC Full
  TVector2 mQ;
  Double_t mQx=0., mQy=0.;

  Float_t eXFsum=0., wXFsum=0., eYFsum=0.,
	  wYFsum=0.,eWFgt=0.,wWFgt=0.,psi_f=0.,psi_f_s=0.;

  for(int iTile = 0; iTile < 16; iTile++) {

    eXFsum += cos(BBC_GetPhi(0,iTile))*nHitE[iTile];
    wXFsum += cos(BBC_GetPhi(1,iTile))*nHitW[iTile];

    eYFsum +=sin(BBC_GetPhi(0,iTile))*nHitE[iTile];
    wYFsum +=sin(BBC_GetPhi(1,iTile))*nHitW[iTile];

    eWFgt += nHitE[iTile];
    wWFgt += nHitW[iTile];

  }
  mQx=(eWFgt>0. && wWFgt>0.) ? wXFsum/wWFgt - eXFsum/eWFgt:0.;
  mQy=(eWFgt>0. && wWFgt>0.) ? wYFsum/wWFgt - eYFsum/eWFgt:0.;

  mQ.Set(mQx,mQy);

  if (mQ.Mod()) {
    psi_f= mQ.Phi();
    if (psi_f < 0.)  psi_f += 2*pi;  //done in BBC_GetPsi()
  }

  //	if(psi_f==0.0) return kFALSE; //done in BBC_GetPsi()

  // shift correction


  //Read shift file
  TFile* pShiftFile = new TFile("./shift11.root","READ");

  //BBC psi shift
  TProfile2D* mBBC_shiftWest_c = (TProfile2D*)pShiftFile->Get("mBBC_shiftwest_cos");
  TProfile2D* mBBC_shiftWest_s = (TProfile2D*)pShiftFile->Get("mBBC_shiftwest_sin");
  TProfile2D* mBBC_shiftEast_c = (TProfile2D*)pShiftFile->Get("mBBC_shifteast_cos");
  TProfile2D* mBBC_shiftEast_s = (TProfile2D*)pShiftFile->Get("mBBC_shifteast_sin");
  TProfile2D* mBBC_shiftFull_c = (TProfile2D*)pShiftFile->Get("mBBC_shiftfull_cos");
  TProfile2D* mBBC_shiftFull_s = (TProfile2D*)pShiftFile->Get("mBBC_shiftfull_sin");

  for(int i=1;i<21;i++) {
    mBBCshiftEast_cos[i-1] = mBBC_shiftEast_c->GetBinContent(i, mCentralityID+2);
    mBBCshiftEast_sin[i-1] = mBBC_shiftEast_s->GetBinContent(i, mCentralityID+2);
    mBBCshiftWest_cos[i-1] = mBBC_shiftWest_c->GetBinContent(i, mCentralityID+2);
    mBBCshiftWest_sin[i-1] = mBBC_shiftWest_s->GetBinContent(i, mCentralityID+2);
    mBBCshiftFull_cos[i-1] = mBBC_shiftFull_c->GetBinContent(i, mCentralityID+2);
    mBBCshiftFull_sin[i-1] = mBBC_shiftFull_s->GetBinContent(i, mCentralityID+2);
  }

  pShiftFile->Close();
  //Initilize the value
  psi_e_s=psi_e;
  psi_w_s=psi_w;
  psi_f_s=psi_f;


  for(int j=1;j<21;j++){
    psi_e_s += 2*(-mBBCshiftEast_sin[j-1]*cos(j*psi_e)+mBBCshiftEast_cos[j-1]*sin(j*psi_e))/(float)j ;
    psi_w_s += 2*(-mBBCshiftWest_sin[j-1]*cos(j*psi_w)+mBBCshiftWest_cos[j-1]*sin(j*psi_w))/(float)j ;
    psi_f_s += 2*(-mBBCshiftFull_sin[j-1]*cos(j*psi_f)+mBBCshiftFull_cos[j-1]*sin(j*psi_f))/(float)j ;
  }

  if(psi_e_s<0.0) psi_e_s +=2*pi;
  if(psi_e_s>2*pi) psi_e_s -=2*pi;

  if(psi_w_s<0.0) psi_w_s +=2*pi;
  if(psi_w_s>2*pi) psi_w_s -=2*pi;

  if(psi_f_s<0.0) psi_f_s +=2*pi;
  if(psi_f_s>2*pi) psi_f_s -=2*pi;


  //BBC psi Shift Correlation:
  // ...

  mPsi[0] = psi_e; 
  mPsi[1] = psi_w;
  mPsi[2] = psi_f;
  mPsi[3] = psi_e_s;
  mPsi[4] = psi_w_s;
  mPsi[5] = psi_f_s;


  return;
}



Double_t* StBBCEventPlane::BBC_GetPsi() {
  //for(int i=0;i<6;i++) {if(mPsi[i] == 0) mPsi[i] = kFALSE; }
  return mPsi;

}



//BBC azimuthal diatribution

Double_t StBBCEventPlane::BBC_GetPhi(int e_w,int iTile) {

  //get phi of BBC tile
  float pi = TMath::Pi();
  float phi_div = pi /6.0;
  float bbc_phi=phi_div;
  gRandom->SetSeed(0);
  switch(iTile) {
  case 0: bbc_phi=3*phi_div;//pi/2
	  break;
  case 1: bbc_phi=phi_div;//pi/6
	  break;
  case 2: bbc_phi=-1*phi_div; // -pi/6
	  break;
  case 3: bbc_phi=-3*phi_div; // -pi/2
	  break;
  case 4: bbc_phi=-5*phi_div; // -5pi/6
	  break;
  case 5: bbc_phi=5*phi_div; // 5pi/6
	  break;
  case 6: bbc_phi= (gRandom->Rndm()>0.5) ? 2*phi_div:4*phi_div; // pi/3 or 2pi/3	
	  break;
  case 7: bbc_phi=3*phi_div; // pi/2
	  break;
  case 8: bbc_phi=phi_div;
	  break;
  case 9: bbc_phi=0.;
	  break;
  case 10: bbc_phi=-phi_div;
	   break;
  case 11: bbc_phi=(gRandom->Rndm()>0.5) ? -2*phi_div:-4*phi_div;
	   break;
  case 12: bbc_phi=-3*phi_div;
	   break;
  case 13: bbc_phi=-5*phi_div;
	   break;
  case 14: bbc_phi=pi;
	   break;
  case 15: bbc_phi=5*phi_div;
	   break;
  }

  if(e_w==0){if (bbc_phi > -0.001){ bbc_phi = pi-bbc_phi;}
    else {bbc_phi= -pi-bbc_phi;}
  }

  if(bbc_phi<0.0) bbc_phi +=2*pi;
  if(bbc_phi>2*pi) bbc_phi -=2*pi;

  return bbc_phi;

}
