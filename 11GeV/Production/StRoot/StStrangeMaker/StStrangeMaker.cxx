#include "StStrangeMaker.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuBTofHit.h"
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuBTofPidTraits.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StBBCEventPlane/StBBCEventPlane.h"
#include "StStrangeConstants.h"

#include "StStrangeV0.h"
#include "StStrangeCut.h"
#include "StStrangePIDTracks.h"
#include "StStrangePIDv1Tracks.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TTree.h"

float pi = TMath::Pi();

ClassImp(StStrangeMaker)

StStrangeMaker :: StStrangeMaker(StMuDstMaker* maker) : StMaker(){

  mMuDstMaker = maker;
  mStrangeCut = new StStrangeCut();	
  mEvents = 0;
  mEventsProcessed = 0;

  controlHitos1  = NULL;
  v1StrangeHistosL  = NULL;
  v1StrangeHistosLBG  = NULL;
  v1StrangeHistosAL  = NULL;
  v1StrangeHistosALBG  = NULL;
  v1StrangeHistosK  = NULL;
  v1StrangeHistosKBG  = NULL;
  v1TrackHistosP  = NULL;
  v1TrackHistosAP  = NULL;
  v1TrackHistosPiP  = NULL;
  v1TrackHistosPiN  = NULL;
  v1TrackHistosKP  = NULL;
  v1TrackHistosKN  = NULL;

  mMixedEventFlag = kFALSE;
  mMixedEventNo = -1;

  for(int cent=0;cent<9;cent++){
    for(int vz=0;vz<10;vz++){
      for(int psi=0;psi<phiBinsMixed;psi++){
	mMixedEventCounter[cent][vz][psi]=0;
	for(int ev=0;ev<mMixedEvents;ev++){
	  mMixedEventMagneticField[cent][vz][psi][ev]=0;
	  mMixedEventPrimaryVertex[cent][vz][psi][ev]=0;
	}
      }
    }
  }

}


StStrangeMaker :: ~StStrangeMaker(){}

Int_t StStrangeMaker :: Init(){
  mCentralityBin = -1;
  mPsi = new Double_t[6];
  InitHist();
  return kStOK;
}  



Int_t StStrangeMaker :: Finish(){
  WriteHist();
  return kStOK;
}  


///////////////          Main Maker               /////////////////////////////////
Int_t StStrangeMaker :: Make(){

  mEvents ++;
  hStatsEvents ->Fill("mEvents",1);
  mCentralityBin = -1; 
  mMixedEventFlag = kFALSE;
  mMixedEventNo=-1;
  for(int i=0;i<6;i++) mPsi[i] = -999;  // this array will give all BBC event plane angels E,W,Full and it's shift corrected angles

  StMuEvent* muEvent  =  mMuDstMaker->muDst()->event() ;

  if (isBadRun(muEvent)) return kStOK;
  if(!(mStrangeCut->passEvent(muEvent,mMuDstMaker))) return kStOK; 

  mCentralityBin =  centrality(muEvent); // get centrality bin 
  if( mCentralityBin == -1) return kStOK;

  mRefMult = muEvent->refMult();
  StThreeVectorF pv      = muEvent->primaryVertexPosition();
  float B                = muEvent->magneticField(); 
  mVz = pv.z();  


  ///////////////           BBC event Plane        /////////////////////////////////

  StBBCEventPlane* BbcEP = new StBBCEventPlane(muEvent,mCentralityBin);
  mPsi = BbcEP->BBC_GetPsi();
  for(int i =0; i<6; i++){if (mPsi[i] == -999) return kStOk;}
  mReactionPlaneE  = mPsi[0];
  mReactionPlaneW  = mPsi[1];
  mReactionPlaneF  = mPsi[2];
  mReactionPlaneES = mPsi[3];
  mReactionPlaneWS = mPsi[4];   
  mReactionPlaneFS = mPsi[5];

  // for shift correction	
  for (int nfill=1;nfill<=20;nfill++){
    pBBC_shifteast_sin ->Fill(nfill-1,mCentralityBin+1,sin(nfill*mReactionPlaneE));
    pBBC_shifteast_cos ->Fill(nfill-1,mCentralityBin+1,cos(nfill*mReactionPlaneE));
    pBBC_shiftwest_sin ->Fill(nfill-1,mCentralityBin+1,sin(nfill*mReactionPlaneW));
    pBBC_shiftwest_cos ->Fill(nfill-1,mCentralityBin+1,cos(nfill*mReactionPlaneW));
    pBBC_shiftfull_sin ->Fill(nfill-1,mCentralityBin+1,sin(nfill*mReactionPlaneF));
    pBBC_shiftfull_cos ->Fill(nfill-1,mCentralityBin+1,cos(nfill*mReactionPlaneF));
  }

  // RP Correlation
  mRPCorrelation = mReactionPlaneES - mReactionPlaneWS + pi;
  if(mRPCorrelation > 2* pi) mRPCorrelation -= 2*pi;
  if(mRPCorrelation < 0    ) mRPCorrelation += 2*pi;
  pRPCorrelation -> Fill(mCentralityBin+1,cos(mRPCorrelation));

  hBbcEventPlaneFS[mCentralityBin]-> Fill( mReactionPlaneFS);
  hBbcEventPlaneES[mCentralityBin]-> Fill( mReactionPlaneES);
  hBbcEventPlaneWS[mCentralityBin]-> Fill( mReactionPlaneWS);
  hBbcEventPlaneF[mCentralityBin]-> Fill( mReactionPlaneF);
  hBbcEventPlaneE[mCentralityBin]-> Fill( mReactionPlaneE);
  hBbcEventPlaneW[mCentralityBin]-> Fill( mReactionPlaneW);

  /////////////  End of BBC event plane ////////////////////////



  mMixedEventVzBin = VzBin(mVz); 
  mMixedEventPsiRpBin = PhiRpBinMixed(mReactionPlaneFS); 
  if(mMixedEventVzBin == 999 || mMixedEventPsiRpBin == 999 ) return kStOK;

  mMixedEventNo = mMixedEventCounter[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin];


  ////////////////////////// Track loop////////////////////////
  TObjArray* gTracks = mMuDstMaker->muDst()->globalTracks() ;
  TObjArrayIter GetTracks(gTracks) ;
  StMuTrack* gTrack ;

  while((gTrack = (StMuTrack*)GetTracks.Next())){
    if(!mStrangeCut->passTrack(gTrack)) continue; //Basic track quality cuts,including pt min and max,  no cuts for protons and pions


    // v1 Strange analysis
    StPhysicalHelixD helix_track = gTrack->helix();
    StThreeVectorF dca           = gTrack->dcaGlobal();
    short c                      = gTrack->charge();
    Bool_t tof                   = mStrangePIDTracks->ToFTrack(gTrack);     //1=ToF track, 0=not a ToF track
    Bool_t proton = 0; 
    Bool_t pion   = 0;

    if(tof){ // then use both TPC and ToF for PID
      proton = mStrangePIDTracks->Proton(gTrack);
      pion   = mStrangePIDTracks->Pion(gTrack);
    }else{ // then use only TPC for PID
      proton = mStrangePIDTracks->ProtonTPC(gTrack);
      pion   = mStrangePIDTracks->PionTPC(gTrack);
    }

    //cout <<gTrack->id() <<"\t"<<proton <<"\t"<<pion <<endl;


    track__ track;
    track.helix = helix_track;
    track.dca = dca;
    track.tof = tof; 
    track.mass2   = mStrangePIDTracks->Mass2Track(gTrack);
    track.nhitsfit = gTrack->nHitsFit();
    track.hitsratio = (1.0*gTrack->nHitsFit(kTpcId))/(1.0*gTrack->nHitsPoss(kTpcId));

      if(proton){
	track.nsigma = gTrack->nSigmaProton();
	mMixedEventFlag = kTRUE;
	if(c>0){
	  mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][mMixedEventNo].push_back(track);
	  mPosProton.push_back(track);
	}
	else if(c<0){
	  mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][mMixedEventNo].push_back(track);
	  mNegProton.push_back(track);
	}
      }
      if(pion){
	track.nsigma = gTrack->nSigmaPion();
	mMixedEventFlag = kTRUE;
	if(c>0){
	  mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][mMixedEventNo].push_back(track);
	  mPosPion.push_back(track);
	}
	else if(c<0){
	  mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][mMixedEventNo].push_back(track);
	  mNegPion.push_back(track);
	}
      }

    // For v1 pi,K,p analysis
    Bool_t tofT =  mStrangePIDv1Tracks->ToFTrack(gTrack);     //1=ToF track, 0=not a ToF track
    if(tofT){
      Bool_t  protonT = mStrangePIDv1Tracks->Proton(gTrack);
      Bool_t  pionT   = mStrangePIDv1Tracks->Pion(gTrack);
      Bool_t  kaonT   = mStrangePIDv1Tracks->Kaon(gTrack);
      if(protonT  || pionT || kaonT)v1Track(gTrack,protonT,pionT,kaonT); 
    }

  } // End of track loop 

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  v1Strange(pv,B);
  clearPIDTracks();

  if(mMixedEventFlag){
    mMixedEventMagneticField[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][mMixedEventNo] = B;
    mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][mMixedEventNo] = pv;
    mMixedEventCounter[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin]+= 1;
  }

  if(mMixedEventNo == (mMixedEvents-1)){ 
    mixEvents();
    clearMixedEvents();
  }

  mEventsProcessed++;
  hStatsEvents ->Fill("mEventsProcessed",1);

  return kStOK;
}  //Make()

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


void StStrangeMaker::v1Track(StMuTrack* gTrack,Bool_t proton,Bool_t pion,Bool_t kaon){


 ////// EP reslution //////
 //Double_t RPresolution[9] = {0.331726,0.439187,0.563428,0.634947,0.639296,0.586541,0.458289,0.303581,0.162335};        // 7
 Double_t RPresolution[9] = {0.303036,0.404264,0.515243,0.582435,0.599519,0.561734,0.455123,0.302329,0.161215};        // 11
 //Double_t RPresolution[9] = {0.196782, 0.267807, 0.338045, 0.39346, 0.413343, 0.396492, 0.32453, 0.225285, 0.120832};  // 14
 //Double_t RPresolution[9] = {0.184772,0.260588,0.340869,0.404803,0.438205,0.43153,0.364991,0.254255,0.145082};         // 19
 //Double_t RPresolution[9] = {0.092077,0.155365,0.225169,0.28477,0.326373,0.332316,0.290974,0.212434,0.11631};          // 27
 //Double_t RPresolution[9] = {0.0305477,0.032294,0.0779149,0.139483,0.187975,0.213022,0.201958,0.155992,0.0976138};     // 39



   mChargeT       = gTrack->charge();
   mMomentumT     = gTrack->p().mag();
   mPhiT          = gTrack->phi();  if(mPhiT< 0) mPhiT+= 2*pi;
   mPhiPsiT       = mPhiT- mReactionPlaneFS;  if(mPhiPsiT< 0 )  mPhiPsiT+= 2*pi;
   mEtaT          = gTrack->eta();
   mPtT           = gTrack->pt();
   mMass2T         = mStrangePIDv1Tracks->Mass2Track(gTrack);
   double pz      = TMath::Sqrt(mMomentumT * mMomentumT - mPtT*mPtT);  if(mEtaT<0)pz = -pz;
   double E       = (Float_t)TMath::Sqrt(mMomentumT * mMomentumT + mMass2T);
   mRapidityT     = 0.5 * log ((E+pz) / (E-pz) );
   mVzT           = mVz;
   mDcaT          = gTrack->dcaGlobal().mag();
   mNHitsT        = gTrack->nHitsFit();
   mHitsPossT     = (1.0*gTrack->nHitsFit(kTpcId))/(1.0*gTrack->nHitsPoss(kTpcId));

   mV1Track       = cos(mPhiPsiT) / RPresolution[mCentralityBin];

   if(proton){
     mNSigmaT = gTrack->nSigmaProton();
     if(mChargeT > 0) hProton -> Fill(); 
     if(mChargeT < 0) hAProton -> Fill();
   }

   if(pion){
     mNSigmaT = gTrack->nSigmaPion();
     if(mChargeT > 0) hPionPos -> Fill();
     if(mChargeT < 0) hPionNeg -> Fill();
   }

   if(kaon){
     mNSigmaT = gTrack->nSigmaKaon();
     if(mChargeT > 0) hKaonPos -> Fill();
     if(mChargeT < 0) hKaonNeg -> Fill();
   }

}



void StStrangeMaker::v1Strange(StThreeVectorF pv,float B ){

  // For Lambda Signal
  for(unsigned int i=0;i<mPosProton.size();i++){
    for(unsigned int j=0;j<mNegPion.size();j++){
      StPhysicalHelixD h_pos = mPosProton[i].helix; 
      StPhysicalHelixD h_neg = mNegPion[j].helix;
      StThreeVectorF dca_pos    = mPosProton[i].dca;
      StThreeVectorF dca_neg    = mNegPion[j].dca; 
      Bool_t tof_pos = mPosProton[i].tof;
      Bool_t tof_neg = mNegPion[j].tof;
      mV0Type = ( (tof_pos * 2) + (tof_neg * 1) );
      StStrangeV0 *mV0 = new StStrangeV0(h_pos,h_neg,pv,pv,B,0); 

      if(mStrangeCut->passLambda(mV0,pv,dca_pos,dca_neg,mV0Type)){

	Double_t energy = mV0 -> energy(proton, pion);
	Double_t pZ     = mV0->momentum().z();

	// event
	mVzL = mVz ;

	// V0
	mDcaV0ToPVL = mV0->dcaV02PV();
	mDecayLengthL = mV0->decayLength();
	mDcaDaughtersL = mV0->dcaDaughters();
	mMomentumL = mV0->momentum().mag();
	mPtL = mV0->momentum().perp() ;
	mV0TypeL = mV0Type;

	// Track
	mNSigmaPL = mPosProton[i].nsigma;
	mNSigmaPiL = mNegPion[j].nsigma;
	mMass2PL = mPosProton[i].mass2;
	mMass2PiL = mNegPion[j].mass2;
	mDcaPToPVL = dca_pos.mag();
	mDcaPiToPVL = dca_neg.mag();
	mNHitsPL = mPosProton[i].nhitsfit;
	mNHitsPiL = mNegPion[j].nhitsfit;
	mHitsPossPL = mPosProton[i].hitsratio;
	mHitsPossPiL = mNegPion[j].hitsratio;

	// Lambda
	mRapidityL = 0.5 * log ((energy + pZ) / (energy - pZ) ) ;
	mPhiL = mV0->momentum().phi(); if(mPhiL < 0) mPhiL += 2*pi ;
	mPhiPsiL =  mPhiL - mReactionPlaneFS; if(mPhiPsiL < 0 )  mPhiPsiL += 2*pi;
	mLambdaMass = mV0->mass();
	mEtaL = mV0->momentum().pseudoRapidity();

	hLambda -> Fill();

      }
      delete mV0;
    }
  }

  // For Anit-Lambda Siganl
  for(unsigned int i=0;i<mNegProton.size();i++){
    for(unsigned int j=0;j<mPosPion.size();j++){
      StPhysicalHelixD h_pos = mPosPion[j].helix; 
      StPhysicalHelixD h_neg = mNegProton[i].helix;
      StThreeVectorF dca_pos    = mPosPion[j].dca;
      StThreeVectorF dca_neg    = mNegProton[i].dca; 
      Bool_t tof_pos = mPosPion[j].tof;
      Bool_t tof_neg = mNegProton[i].tof;
      mV0Type = ( (tof_neg*2) + (tof_pos*1) );
      StStrangeV0 *mV0 = new StStrangeV0(h_pos,h_neg,pv,pv,B,0); 

      if(mStrangeCut->passLbar(mV0,pv,dca_pos,dca_neg,mV0Type)) {
	Double_t energy = mV0 -> energy(pion, proton);
	Double_t pZ     = mV0->momentum().z();

	// event
	mVzLB = mVz ;

	// V0
	mDcaV0ToPVLB = mV0->dcaV02PV();
	mDecayLengthLB = mV0->decayLength();
	mDcaDaughtersLB = mV0->dcaDaughters();
	mMomentumLB = mV0->momentum().mag();
	mPtLB = mV0->momentum().perp() ;
	mV0TypeLB = mV0Type;

	// Track
	mNSigmaPLB = mNegProton[i].nsigma;
	mNSigmaPiLB = mPosPion[j].nsigma;
	mMass2PLB = mNegProton[i].mass2;
	mMass2PiLB = mPosPion[j].mass2;
	mDcaPToPVLB = dca_neg.mag();
	mDcaPiToPVLB = dca_pos.mag();
	mNHitsPLB = mNegProton[i].nhitsfit;
	mNHitsPiLB = mPosPion[j].nhitsfit;
	mHitsPossPLB = mNegProton[i].hitsratio;
	mHitsPossPiLB = mPosPion[j].hitsratio;

	// Lambda
	mRapidityLB = 0.5 * log ((energy + pZ) / (energy - pZ) ) ;
	mPhiLB = mV0->momentum().phi(); if(mPhiLB < 0) mPhiLB += 2*pi ;
	mPhiPsiLB =  mPhiLB - mReactionPlaneFS; if(mPhiPsiLB < 0 )  mPhiPsiLB += 2*pi;
	mLBarMass = mV0->mass();
	mEtaLB = mV0->momentum().pseudoRapidity();

	hALambda -> Fill();
      }
      delete mV0;
    }
  }

  // For k-Short Signal
  for(unsigned int i=0;i<mPosPion.size();i++){
    for(unsigned int j=0;j<mNegPion.size();j++){
      StPhysicalHelixD h_pos = mPosPion[i].helix; 
      StPhysicalHelixD h_neg = mNegPion[j].helix;
      StThreeVectorF dca_pos    = mPosPion[i].dca;
      StThreeVectorF dca_neg    = mNegPion[j].dca; 
      StStrangeV0 *mV0 = new StStrangeV0(h_pos,h_neg,pv,pv,B,0); 

      if(mStrangeCut->passKs(mV0,pv,dca_pos,dca_neg)) {
	Double_t energy = mV0 -> energy(pion, pion);
	Double_t pZ     = mV0->momentum().z();

	// event
	mVzK = mVz;

	// V0
	mDcaV0ToPVK = mV0->dcaV02PV();
	mDecayLengthK = mV0->decayLength();
	mDcaDaughtersK = mV0->dcaDaughters();
	mMomentumK = mV0->momentum().mag();
	mPtK = mV0->momentum().perp() ;
	mV0TypeK = mV0Type;

	// Track
	mNSigmaPiPK = mPosPion[i].nsigma;
	mNSigmaPiNK = mNegPion[j].nsigma;
	mMass2PiPK = mPosPion[i].mass2;
	mMass2PiNK = mNegPion[j].mass2;
	mDcaPiPToPVK = dca_pos.mag();
	mDcaPiNToPVK = dca_neg.mag();
	mNHitsPiPK = mPosPion[i].nhitsfit;
	mNHitsPiNK = mNegPion[j].nhitsfit;
	mHitsPossPiPK = mPosPion[i].hitsratio;
	mHitsPossPiNK = mNegPion[j].hitsratio;

	// Kshort
	mRapidityK = 0.5 * log ((energy + pZ) / (energy - pZ) ) ;
	mPhiK = mV0->momentum().phi(); if(mPhiK < 0) mPhiK += 2*pi ;
	mPhiPsiK =  mPhiK - mReactionPlaneFS; if(mPhiPsiK < 0 )  mPhiPsiK += 2*pi;
	mKshortMass = mV0->mass();
	mEtaK = mV0->momentum().pseudoRapidity();

	hKshort-> Fill();
      }
      delete mV0;
    }
  }

} 

void StStrangeMaker::mixEvents(){
  // mixed background for lambda protons mixed with pions in second event


  for(int ev1=0;ev1<mMixedEvents;ev1++){
    for(unsigned int t1=0;t1<mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1].size();t1++){

      for(int ev2=0;ev2<mMixedEvents;ev2++){
	if(ev1 == ev2) continue;
	for(unsigned int t2=0;t2<mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2].size();t2++){

	  StThreeVectorF pv1 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1];
	  StThreeVectorF pv2 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2];
	  float B            =   mMixedEventMagneticField[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1];

	  StPhysicalHelixD h_pos = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].helix;
	  StThreeVectorF dca_pos = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].dca;
	  Bool_t tof_pos         = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].tof;
	  StPhysicalHelixD h_neg =   mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].helix;
	  StThreeVectorF dca_neg =   mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].dca;
	  Bool_t tof_neg         =   mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].tof;

	  mV0Type = ( (tof_pos * 2) + (tof_neg * 1) ); 

	  StStrangeV0 *mV0 = new StStrangeV0(h_pos,h_neg,pv1,pv2,B,1); 
	  if(mStrangeCut->passLambda(mV0,pv1,dca_pos,dca_neg, mV0Type)) {

	    Double_t energy = mV0 -> energy(proton, pion);
	    Double_t pZ     = mV0->momentum().z();

	    // event
	    mVzL = mVz ;

	    // V0
	    mDcaV0ToPVL = mV0->dcaV02PV();
	    mDecayLengthL = mV0->decayLength();
	    mDcaDaughtersL = mV0->dcaDaughters();
	    mMomentumL = mV0->momentum().mag();
	    mPtL = mV0->momentum().perp() ;
	    mV0TypeL = mV0Type;

	    // Track
	    mNSigmaPL = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nsigma;
	    mNSigmaPiL = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nsigma;
	    mMass2PL = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].mass2;
	    mMass2PiL = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].mass2;
	    mDcaPToPVL = dca_pos.mag();
	    mDcaPiToPVL = dca_neg.mag();
	    mNHitsPL = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nhitsfit;
	    mNHitsPiL = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nhitsfit;
	    mHitsPossPL = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].hitsratio;
	    mHitsPossPiL = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].hitsratio;

	    // Lambda
	    mRapidityL = 0.5 * log ((energy + pZ) / (energy - pZ) ) ;
	    mPhiL = mV0->momentum().phi(); if(mPhiL < 0) mPhiL += 2*pi ;
	    mPhiPsiL =  mPhiL - mReactionPlaneFS; if(mPhiPsiL < 0 )  mPhiPsiL += 2*pi;
	    mLambdaMass = mV0->mass();
	    mEtaL = mV0->momentum().pseudoRapidity();

	    hLambdaBG-> Fill();

	  }
	  delete mV0;

	}
      }
    }
  }


  // mixed background for lambda pions mixed with protons in second event
  for(int ev1=0;ev1<mMixedEvents;ev1++){
    for(unsigned int t1=0;t1<mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1].size();t1++){

      for(int ev2=0;ev2<mMixedEvents;ev2++){
	if(ev1 == ev2) continue;
	for(unsigned int t2=0;t2<mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2].size();t2++){

	  StThreeVectorF pv1 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2];
	  StThreeVectorF pv2 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1];
	  float B =  mMixedEventMagneticField[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2];

	  StPhysicalHelixD h_pos = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].helix;
	  StThreeVectorF dca_pos = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].dca;
	  Bool_t tof_pos         = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].tof;
	  StPhysicalHelixD h_neg =   mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].helix;
	  StThreeVectorF dca_neg =   mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].dca;
	  Bool_t tof_neg         =   mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].tof;

	  mV0Type = ( (tof_pos * 2) + (tof_neg * 1) );

	  StStrangeV0 *mV0 = new StStrangeV0(h_pos,h_neg,pv1,pv2,B,1); 
	  if(mStrangeCut->passLambda(mV0,pv1,dca_pos,dca_neg,mV0Type)) {

	    Double_t energy = mV0 -> energy(proton, pion);
	    Double_t pZ     = mV0->momentum().z();

	    // event
	    mVzL = mVz ;

	    // V0
	    mDcaV0ToPVL = mV0->dcaV02PV();
	    mDecayLengthL = mV0->decayLength();
	    mDcaDaughtersL = mV0->dcaDaughters();
	    mMomentumL = mV0->momentum().mag();
	    mPtL = mV0->momentum().perp() ;
	    mV0TypeL = mV0Type;

	    // Track
	    mNSigmaPL = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nsigma;
	    mNSigmaPiL = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nsigma;
	    mMass2PL = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].mass2;
	    mMass2PiL = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].mass2;
	    mDcaPToPVL = dca_pos.mag();
	    mDcaPiToPVL = dca_neg.mag();
	    mNHitsPL = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nhitsfit;
	    mNHitsPiL = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nhitsfit;
	    mHitsPossPL = mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].hitsratio;
	    mHitsPossPiL = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].hitsratio;

	    // Lambda
	    mRapidityL = 0.5 * log ((energy + pZ) / (energy - pZ) ) ;
	    mPhiL = mV0->momentum().phi(); if(mPhiL < 0) mPhiL += 2*pi ;
	    mPhiPsiL =  mPhiL - mReactionPlaneFS; if(mPhiPsiL < 0 )  mPhiPsiL += 2*pi;
	    mLambdaMass = mV0->mass();
	    mEtaL = mV0->momentum().pseudoRapidity();

	    hLambdaBG-> Fill();
	  }
	  delete mV0;

	}
      }
    }
  }







  // mixed background for anti-lambda pions mixed with protons in second event

  for(int ev1=0;ev1<mMixedEvents;ev1++){
    for(unsigned int t1=0;t1<mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1].size();t1++){

      for(int ev2=0;ev2<mMixedEvents;ev2++){
	if(ev1 == ev2) continue;
	for(unsigned int t2=0;t2<mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2].size();t2++){

	  StThreeVectorF pv1 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1];
	  StThreeVectorF pv2 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2];
	  float B =  mMixedEventMagneticField[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1];

	  StPhysicalHelixD h_pos =    mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].helix;
	  StThreeVectorF dca_pos =    mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].dca;
	  Bool_t tof_pos         =    mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].tof;
	  StPhysicalHelixD h_neg =  mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].helix;
	  StThreeVectorF dca_neg =  mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].dca;
	  Bool_t tof_neg         =  mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].tof;

	  mV0Type = ( (tof_neg * 2) + (tof_pos * 1) );

	  StStrangeV0 *mV0 = new StStrangeV0(h_pos,h_neg,pv1,pv2,B,1); 
	  if(mStrangeCut->passLbar(mV0,pv1,dca_pos,dca_neg,mV0Type)) {

	    Double_t energy = mV0 -> energy(pion, proton);
	    Double_t pZ     = mV0->momentum().z();

	    // event
	    mVzLB = mVz ;

	    // V0
	    mDcaV0ToPVLB = mV0->dcaV02PV();
	    mDecayLengthLB = mV0->decayLength();
	    mDcaDaughtersLB = mV0->dcaDaughters();
	    mMomentumLB = mV0->momentum().mag();
	    mPtLB = mV0->momentum().perp() ;
	    mV0TypeLB = mV0Type;

	    // Track
	    mNSigmaPLB = mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nsigma;
	    mNSigmaPiLB = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nsigma;
	    mMass2PLB = mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].mass2;
	    mMass2PiLB = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].mass2;
	    mDcaPToPVLB = dca_neg.mag();
	    mDcaPiToPVLB = dca_pos.mag();
	    mNHitsPLB = mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nhitsfit;
	    mNHitsPiLB = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nhitsfit;
	    mHitsPossPLB = mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].hitsratio;
	    mHitsPossPiLB = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].hitsratio;

	    // Lambda
	    mRapidityLB = 0.5 * log ((energy + pZ) / (energy - pZ) ) ;
	    mPhiLB = mV0->momentum().phi(); if(mPhiLB < 0) mPhiLB += 2*pi ;
	    mPhiPsiLB =  mPhiLB - mReactionPlaneFS; if(mPhiPsiLB < 0 )  mPhiPsiLB += 2*pi;
	    mLBarMass = mV0->mass();
	    mEtaLB = mV0->momentum().pseudoRapidity();

	    hALambdaBG-> Fill();

	  }
	  delete mV0;
	}
      }
    }
  }



  // mixed background for anti-lambda pions mixed with protons in second event
  for(int ev1=0;ev1<mMixedEvents;ev1++){
    for(unsigned int t1=0;t1<mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1].size();t1++){

      for(int ev2=0;ev2<mMixedEvents;ev2++){
	if(ev1 == ev2) continue;
	for(unsigned int t2=0;t2<mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2].size();t2++){

	  StThreeVectorF pv1 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2];
	  StThreeVectorF pv2 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1];
	  float B =  mMixedEventMagneticField[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2];

	  StPhysicalHelixD h_pos    =    mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].helix;
	  StThreeVectorF dca_pos    =    mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].dca;
	  Bool_t tof_pos            =    mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].tof;
	  StPhysicalHelixD h_neg    =  mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].helix;
	  StThreeVectorF dca_neg    =  mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].dca;
	  Bool_t tof_neg            =  mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].tof;

	  mV0Type = ( (tof_neg * 2) + (tof_pos * 1) );

	  StStrangeV0 *mV0 = new StStrangeV0(h_pos,h_neg,pv1,pv2,B,1); 
	  if(mStrangeCut->passLbar(mV0,pv1,dca_pos,dca_neg,mV0Type)) {

	    Double_t energy = mV0 -> energy(pion, proton);
	    Double_t pZ     = mV0->momentum().z();

	    // event
	    mVzLB = mVz ;

	    // V0
	    mDcaV0ToPVLB = mV0->dcaV02PV();
	    mDecayLengthLB = mV0->decayLength();
	    mDcaDaughtersLB = mV0->dcaDaughters();
	    mMomentumLB = mV0->momentum().mag();
	    mPtLB = mV0->momentum().perp() ;
	    mV0TypeLB = mV0Type;

	    // Track
	    mNSigmaPLB = mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nsigma;
	    mNSigmaPiLB = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nsigma;
	    mMass2PLB = mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].mass2;
	    mMass2PiLB = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].mass2;
	    mDcaPToPVLB = dca_neg.mag();
	    mDcaPiToPVLB = dca_pos.mag();
	    mNHitsPLB = mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nhitsfit;
	    mNHitsPiLB = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nhitsfit;
	    mHitsPossPLB = mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].hitsratio;
	    mHitsPossPiLB = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].hitsratio;

	    // Lambda
	    mRapidityLB = 0.5 * log ((energy + pZ) / (energy - pZ) ) ;
	    mPhiLB = mV0->momentum().phi(); if(mPhiLB < 0) mPhiLB += 2*pi ;
	    mPhiPsiLB =  mPhiLB - mReactionPlaneFS; if(mPhiPsiLB < 0 )  mPhiPsiLB += 2*pi;
	    mLBarMass = mV0->mass();
	    mEtaLB = mV0->momentum().pseudoRapidity();



	    hALambdaBG-> Fill();

	  }
	  delete mV0;
	}
      }
    }
  }




  // mixed background for k-Short pions mixed with pions in second event

  for(int ev1=0;ev1<mMixedEvents;ev1++){
    for(unsigned int t1=0;t1<mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1].size();t1++){

      for(int ev2=0;ev2<mMixedEvents;ev2++){
	if(ev1 == ev2) continue;
	for(unsigned int t2=0;t2<mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2].size();t2++){

	  StThreeVectorF pv1 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1];
	  StThreeVectorF pv2 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2];
	  float B =  mMixedEventMagneticField[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1];

	  StPhysicalHelixD h_pos = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].helix;
	  StThreeVectorF dca1    = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].dca;
	  StPhysicalHelixD h_neg =  mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].helix;
	  StThreeVectorF dca2    =  mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].dca;

	  StStrangeV0 *mV0 = new StStrangeV0(h_pos,h_neg,pv1,pv2,B,1); 
	  if(mStrangeCut->passKs(mV0,pv1,dca1,dca2)) {

	    Double_t energy = mV0 -> energy(pion, pion);
	    Double_t pZ     = mV0->momentum().z();

	    // event
	    mVzK = mVz;

	    // V0
	    mDcaV0ToPVK = mV0->dcaV02PV();
	    mDecayLengthK = mV0->decayLength();
	    mDcaDaughtersK = mV0->dcaDaughters();
	    mMomentumK = mV0->momentum().mag();
	    mPtK = mV0->momentum().perp() ;
	    mV0TypeK = mV0Type;

	    // Track
	    mNSigmaPiPK = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nsigma;
	    mNSigmaPiNK = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nsigma;
	    mMass2PiPK = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].mass2;
	    mMass2PiNK = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].mass2;
	    mDcaPiPToPVK = dca1.mag();
	    mDcaPiNToPVK = dca2.mag();
	    mNHitsPiPK = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nhitsfit;
	    mNHitsPiNK = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nhitsfit;
	    mHitsPossPiPK = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].hitsratio;
	    mHitsPossPiNK = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].hitsratio;

	    // Kshort
	    mRapidityK = 0.5 * log ((energy + pZ) / (energy - pZ) ) ;
	    mPhiK = mV0->momentum().phi(); if(mPhiK < 0) mPhiK += 2*pi ;
	    mPhiPsiK =  mPhiK - mReactionPlaneFS; if(mPhiPsiK < 0 )  mPhiPsiK += 2*pi;
	    mKshortMass = mV0->mass();
	    mEtaK = mV0->momentum().pseudoRapidity();

	    hKshortBG -> Fill();

	  }
	  delete mV0;
	}
      }
    }
  }


  // mixed background for k-short  pions mixed with protons in second event
  for(int ev1=0;ev1<mMixedEvents;ev1++){
    for(unsigned int t1=0;t1<mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1].size();t1++){

      for(int ev2=0;ev2<mMixedEvents;ev2++){
	if(ev1 == ev2) continue;
	for(unsigned int t2=0;t2<mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2].size();t2++){

	  StThreeVectorF pv1 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2];
	  StThreeVectorF pv2 =   mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1];
	  float B =  mMixedEventMagneticField[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2];

	  StPhysicalHelixD h_pos = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].helix;
	  StThreeVectorF dca1    = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].dca;
	  StPhysicalHelixD h_neg =  mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].helix;
	  StThreeVectorF dca2    =  mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].dca;

	  StStrangeV0 *mV0 = new StStrangeV0(h_pos,h_neg,pv1,pv2,B,1); 
	  if(mStrangeCut->passKs(mV0,pv1,dca1,dca2)) {

	    Double_t energy = mV0 -> energy(pion, pion);
	    Double_t pZ     = mV0->momentum().z();

	    // event
	    mVzK = mVz;

	    // V0
	    mDcaV0ToPVK = mV0->dcaV02PV();
	    mDecayLengthK = mV0->decayLength();
	    mDcaDaughtersK = mV0->dcaDaughters();
	    mMomentumK = mV0->momentum().mag();
	    mPtK = mV0->momentum().perp() ;
	    mV0TypeK = mV0Type;

	    // Track
	    mNSigmaPiPK = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nsigma;
	    mNSigmaPiNK = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nsigma;
	    mMass2PiPK = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].mass2;
	    mMass2PiNK = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].mass2;
	    mDcaPiPToPVK = dca1.mag();
	    mDcaPiNToPVK = dca2.mag();
	    mNHitsPiPK = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].nhitsfit;
	    mNHitsPiNK = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].nhitsfit;
	    mHitsPossPiPK = mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev2][t2].hitsratio;
	    mHitsPossPiNK = mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev1][t1].hitsratio;

	    // Kshort
	    mRapidityK = 0.5 * log ((energy + pZ) / (energy - pZ) ) ;
	    mPhiK = mV0->momentum().phi(); if(mPhiK < 0) mPhiK += 2*pi ;
	    mPhiPsiK =  mPhiK - mReactionPlaneFS; if(mPhiPsiK < 0 )  mPhiPsiK += 2*pi;
	    mKshortMass = mV0->mass();
	    mEtaK = mV0->momentum().pseudoRapidity();


	    hKshortBG -> Fill();


	  }

	  delete mV0;

	}
      }
    }
  }


}

void StStrangeMaker :: clearPIDTracks(){
  mPosProton.clear();
  mNegProton.clear();
  mPosPion.clear();
  mNegPion.clear();
}

void StStrangeMaker :: clearMixedEvents(){
  for(int ev=0;ev<mMixedEvents;ev++){
    mMixedEventProtonPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev].clear();
    mMixedEventProtonNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev].clear();
    mMixedEventPionPos[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev].clear();
    mMixedEventPionNeg[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev].clear();
    mMixedEventPrimaryVertex[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev] = 0;
    mMixedEventMagneticField[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin][ev] = 0;
  }
  mMixedEventCounter[mCentralityBin][mMixedEventVzBin][mMixedEventPsiRpBin] = 0;
}


Bool_t StStrangeMaker::isBadRun(StMuEvent* ev){
  StRefMultCorr* refmultCorrUtil1 = new StRefMultCorr();
  Int_t RunId      = ev->runNumber();
  refmultCorrUtil1->init(RunId);
  if ( refmultCorrUtil1->isBadRun(RunId) ) {
    delete refmultCorrUtil1;
    return kTRUE;
  }else {
    delete refmultCorrUtil1;
    return kFALSE;
  }
} 



Int_t StStrangeMaker::centrality(StMuEvent* ev){
  StRefMultCorr* refmultCorrUtil = new StRefMultCorr();
  const Double_t zdcCoincidenceRate = 0.0 ; // Hz
  Int_t RunId      = ev->runNumber();
  refmultCorrUtil->init(RunId);
  Int_t refmult = ev->refMult() ;
  Double_t vz      = ev->primaryVertexPosition().z();
  refmultCorrUtil->initEvent(refmult, vz, zdcCoincidenceRate);
  Int_t CentralityID = refmultCorrUtil->getCentralityBin9();
  delete refmultCorrUtil;
  return CentralityID;
}


////////////////////////////// //////////////////////////////

Int_t StStrangeMaker::VzBin(Double_t vz){
  Double_t vzMax = Constant::mVzMax;
  for(int i=-5;i<5;i++){
    Double_t vz1 = i*(vzMax/5.0);
    Double_t vz2 = i*(vzMax/5.0) + (vzMax/5.0) ;
    if( (vz1 <= vz) && (vz <= vz2) ){
      return  i+5;
    }
  }
  return 999;
}

Int_t StStrangeMaker:: PhiRpBinMixed(Double_t phiRp){

  float two_pi = 6.29; 
  for(int i=0;i<phiBinsMixed;i++){
    Double_t RPrange1 = i*(two_pi/phiBinsMixed);
    Double_t RPrange2 = i*(two_pi/phiBinsMixed) + (two_pi/phiBinsMixed);
    if( (RPrange1 <= phiRp) && (phiRp <= RPrange2) ){
      return i;

    }
  }
  return 999;

}


///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

void StStrangeMaker :: InitHist(){

  // Statistics Histograms
  hStatsEvents	= new TH1D("hStatsEvents","",2,0,1);
  hStatsEvents ->Fill("mEvents",0);
  hStatsEvents ->Fill("mEventsProcessed",0);

  // BBC histos

  char buff[255];
  for(int i=0;i<9;i++){
    sprintf(buff,"BBC_F_S_%d",i);
    hBbcEventPlaneFS[i]  =  new TH1D(buff,buff,100,0,6.5);
    hBbcEventPlaneFS[i] -> Sumw2();  

    sprintf(buff,"BBC_E_S_%d",i);
    hBbcEventPlaneES[i]  =  new TH1D(buff,buff,100,0,6.5);
    hBbcEventPlaneES[i] -> Sumw2();  

    sprintf(buff,"BBC_W_S_%d",i);
    hBbcEventPlaneWS[i]  =  new TH1D(buff,buff,100,0,6.5);
    hBbcEventPlaneWS[i] -> Sumw2();  

    sprintf(buff,"BBC_F_%d",i);
    hBbcEventPlaneF[i]  =  new TH1D(buff,buff,100,0,6.5);
    hBbcEventPlaneF[i] -> Sumw2();  

    sprintf(buff,"BBC_E_%d",i);
    hBbcEventPlaneE[i]  =  new TH1D(buff,buff,100,0,6.5);
    hBbcEventPlaneE[i] -> Sumw2();  

    sprintf(buff,"BBC_W_%d",i);
    hBbcEventPlaneW[i]  =  new TH1D(buff,buff,100,0,6.5);
    hBbcEventPlaneW[i] -> Sumw2();  

  }

  pBBC_shifteast_sin = new TProfile2D("mBBC_shifteast_sin","",20,0,20,10,0,10.0,-1.0,1.0);
  pBBC_shifteast_cos = new TProfile2D("mBBC_shifteast_cos","",20,0,20,10,0,10.0,-1.0,1.0);
  pBBC_shiftwest_sin = new TProfile2D("mBBC_shiftwest_sin","",20,0,20,10,0,10.0,-1.0,1.0);
  pBBC_shiftwest_cos = new TProfile2D("mBBC_shiftwest_cos","",20,0,20,10,0,10.0,-1.0,1.0);
  pBBC_shiftfull_sin = new TProfile2D("mBBC_shiftfull_sin","",20,0,20,10,0,10.0,-1.0,1.0);
  pBBC_shiftfull_cos = new TProfile2D("mBBC_shiftfull_cos","",20,0,20,10,0,10.0,-1.0,1.0);
  pRPCorrelation     = new TProfile("pRPcorr","pRPcorr",10,0,10);
  pBBC_shifteast_sin -> Sumw2();  
  pBBC_shifteast_cos-> Sumw2();  
  pBBC_shiftwest_sin-> Sumw2();  
  pBBC_shiftwest_cos-> Sumw2();  
  pBBC_shiftfull_sin-> Sumw2();  
  pBBC_shiftfull_cos-> Sumw2();  
  pRPCorrelation-> Sumw2();  
//////////////////////////////////////////////////////////////////////////////////


  // Track Histos
  hProton  = new TTree("proton","RECREATE");
  hAProton = new TTree("aproton","RECREATE");
  hPionPos = new TTree("pionpos","RECREATE");
  hPionNeg = new TTree("pionneg","RECREATE");
  hKaonPos = new TTree("kaonpos","RECREATE");
  hKaonNeg = new TTree("kaonneg","RECREATE");

  hProton -> Branch("mVzT",&mVzT,"mVzT/D");
  hProton -> Branch("mMomentumT",&mMomentumT,"mMomentumT/D");
  hProton -> Branch("mPtT",&mPtT,"mPtT/D");
  hProton -> Branch("mNSigmaT",&mNSigmaT,"mNSigmaT/D");
  hProton -> Branch("mMass2T",&mMass2T,"mMass2T/D");
  hProton -> Branch("mDcaT",&mDcaT,"mDcaT/D");
  hProton -> Branch("mNHitsT",&mNHitsT,"mNHitsT/D");
  hProton -> Branch("mHitsPossT",&mHitsPossT,"mHitsPossT/D");
  hProton -> Branch("mRapidityT",&mRapidityT,"mRapidityT/D");
  hProton -> Branch("mPhiPsiT",&mPhiPsiT,"mPhiPsiT/D");
  hProton -> Branch("mPhiT",&mPhiT,"mPhiT/D");
  hProton -> Branch("mEtaT",&mEtaT,"mEtaT/D");
  hProton -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");
  hProton -> Branch("mV1Track",&mV1Track,"mV1Track/D");

  hAProton -> Branch("mVzT",&mVzT,"mVzT/D");
  hAProton -> Branch("mMomentumT",&mMomentumT,"mMomentumT/D");
  hAProton -> Branch("mPtT",&mPtT,"mPtT/D");
  hAProton -> Branch("mNSigmaT",&mNSigmaT,"mNSigmaT/D");
  hAProton -> Branch("mMass2T",&mMass2T,"mMass2T/D");
  hAProton -> Branch("mDcaT",&mDcaT,"mDcaT/D");
  hAProton -> Branch("mNHitsT",&mNHitsT,"mNHitsT/D");
  hAProton -> Branch("mHitsPossT",&mHitsPossT,"mHitsPossT/D");
  hAProton -> Branch("mRapidityT",&mRapidityT,"mRapidityT/D");
  hAProton -> Branch("mPhiPsiT",&mPhiPsiT,"mPhiPsiT/D");
  hAProton -> Branch("mPhiT",&mPhiT,"mPhiT/D");
  hAProton -> Branch("mEtaT",&mEtaT,"mEtaT/D");
  hAProton -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");
  hAProton -> Branch("mV1Track",&mV1Track,"mV1Track/D");

  hPionPos -> Branch("mVzT",&mVzT,"mVzT/D");
  hPionPos -> Branch("mMomentumT",&mMomentumT,"mMomentumT/D");
  hPionPos -> Branch("mPtT",&mPtT,"mPtT/D");
  hPionPos -> Branch("mNSigmaT",&mNSigmaT,"mNSigmaT/D");
  hPionPos -> Branch("mMass2T",&mMass2T,"mMass2T/D");
  hPionPos -> Branch("mDcaT",&mDcaT,"mDcaT/D");
  hPionPos -> Branch("mNHitsT",&mNHitsT,"mNHitsT/D");
  hPionPos -> Branch("mHitsPossT",&mHitsPossT,"mHitsPossT/D");
  hPionPos -> Branch("mRapidityT",&mRapidityT,"mRapidityT/D");
  hPionPos -> Branch("mPhiPsiT",&mPhiPsiT,"mPhiPsiT/D");
  hPionPos -> Branch("mPhiT",&mPhiT,"mPhiT/D");
  hPionPos -> Branch("mEtaT",&mEtaT,"mEtaT/D");
  hPionPos -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");
  hPionPos -> Branch("mV1Track",&mV1Track,"mV1Track/D");

  hPionNeg -> Branch("mVzT",&mVzT,"mVzT/D");
  hPionNeg -> Branch("mMomentumT",&mMomentumT,"mMomentumT/D");
  hPionNeg -> Branch("mPtT",&mPtT,"mPtT/D");
  hPionNeg -> Branch("mNSigmaT",&mNSigmaT,"mNSigmaT/D");
  hPionNeg -> Branch("mMass2T",&mMass2T,"mMass2T/D");
  hPionNeg -> Branch("mDcaT",&mDcaT,"mDcaT/D");
  hPionNeg -> Branch("mNHitsT",&mNHitsT,"mNHitsT/D");
  hPionNeg -> Branch("mHitsPossT",&mHitsPossT,"mHitsPossT/D");
  hPionNeg -> Branch("mRapidityT",&mRapidityT,"mRapidityT/D");
  hPionNeg -> Branch("mPhiPsiT",&mPhiPsiT,"mPhiPsiT/D");
  hPionNeg -> Branch("mPhiT",&mPhiT,"mPhiT/D");
  hPionNeg -> Branch("mEtaT",&mEtaT,"mEtaT/D");
  hPionNeg -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");
  hPionNeg -> Branch("mV1Track",&mV1Track,"mV1Track/D");

  hKaonPos -> Branch("mVzT",&mVzT,"mVzT/D");
  hKaonPos -> Branch("mMomentumT",&mMomentumT,"mMomentumT/D");
  hKaonPos -> Branch("mPtT",&mPtT,"mPtT/D");
  hKaonPos -> Branch("mNSigmaT",&mNSigmaT,"mNSigmaT/D");
  hKaonPos -> Branch("mMass2T",&mMass2T,"mMass2T/D");
  hKaonPos -> Branch("mDcaT",&mDcaT,"mDcaT/D");
  hKaonPos -> Branch("mNHitsT",&mNHitsT,"mNHitsT/D");
  hKaonPos -> Branch("mHitsPossT",&mHitsPossT,"mHitsPossT/D");
  hKaonPos -> Branch("mRapidityT",&mRapidityT,"mRapidityT/D");
  hKaonPos -> Branch("mPhiPsiT",&mPhiPsiT,"mPhiPsiT/D");
  hKaonPos -> Branch("mPhiT",&mPhiT,"mPhiT/D");
  hKaonPos -> Branch("mEtaT",&mEtaT,"mEtaT/D");
  hKaonPos -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");
  hKaonPos -> Branch("mV1Track",&mV1Track,"mV1Track/D");

  hKaonNeg -> Branch("mVzT",&mVzT,"mVzT/D");
  hKaonNeg -> Branch("mMomentumT",&mMomentumT,"mMomentumT/D");
  hKaonNeg -> Branch("mPtT",&mPtT,"mPtT/D");
  hKaonNeg -> Branch("mNSigmaT",&mNSigmaT,"mNSigmaT/D");
  hKaonNeg -> Branch("mMass2T",&mMass2T,"mMass2T/D");
  hKaonNeg -> Branch("mDcaT",&mDcaT,"mDcaT/D");
  hKaonNeg -> Branch("mNHitsT",&mNHitsT,"mNHitsT/D");
  hKaonNeg -> Branch("mHitsPossT",&mHitsPossT,"mHitsPossT/D");
  hKaonNeg -> Branch("mRapidityT",&mRapidityT,"mRapidityT/D");
  hKaonNeg -> Branch("mPhiPsiT",&mPhiPsiT,"mPhiPsiT/D");
  hKaonNeg -> Branch("mPhiT",&mPhiT,"mPhiT/D");
  hKaonNeg -> Branch("mEtaT",&mEtaT,"mEtaT/D");
  hKaonNeg -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");
  hKaonNeg -> Branch("mV1Track",&mV1Track,"mV1Track/D");

//////////////////////////////////////////////////////////////////////////////////


  // Strange Histos
  hLambda     = new TTree("lambda","RECREATE");
  hALambda    = new TTree("alambda","RECREATE");
  hKshort     = new TTree("kshort","RECREATE");

  hLambdaBG   = new TTree("lambdaBg","RECREATE");
  hALambdaBG  = new TTree("alambdaBg","RECREATE");
  hKshortBG   = new TTree("kshortBg","RECREATE");

  hLambda -> Branch("mVzL",&mVzL,"mVzL/D");
  hLambda -> Branch("mDcaV0ToPVL",&mDcaV0ToPVL,"mDcaV0ToPVL/D");
  hLambda -> Branch("mDecayLengthL",&mDecayLengthL,"mDecayLengthL/D");
  hLambda -> Branch("mDcaDaughtersL",&mDcaDaughtersL,"mDcaDaughtersL/D");
  hLambda -> Branch("mMomentumL",&mMomentumL,"mMomentumL/D");
  hLambda -> Branch("mPtL",&mPtL,"mPtL/D");
  hLambda -> Branch("mV0TypeL",&mV0TypeL,"mV0TypeL/D");
  hLambda -> Branch("mNSigmaPL",&mNSigmaPL,"mNSigmaPL/D");
  hLambda -> Branch("mNSigmaPiL",&mNSigmaPiL,"mNSigmaPiL/D");
  hLambda -> Branch("mMass2PL",&mMass2PL,"mMass2PL/D");
  hLambda -> Branch("mMass2PiL",&mMass2PiL,"mMass2PiL/D");
  hLambda -> Branch("mDcaPToPVL",&mDcaPToPVL,"mDcaPToPVL/D");
  hLambda -> Branch("mDcaPiToPVL",&mDcaPiToPVL,"mDcaPiToPVL/D");
  hLambda -> Branch("mNHitsPL",&mNHitsPL,"mNHitsPL/D");
  hLambda -> Branch("mNHitsPiL",&mNHitsPiL,"mNHitsPiL/D");
  hLambda -> Branch("mHitsPossPL",&mHitsPossPL,"mHitsPossPL/D");
  hLambda -> Branch("mHitsPossPiL",&mHitsPossPiL,"mHitsPossPiL/D");
  hLambda -> Branch("mRapidityL",&mRapidityL,"mRapidityL/D");
  hLambda -> Branch("mPhiPsiL",&mPhiPsiL,"mPhiPsiL/D");
  hLambda -> Branch("mPhiL",&mPhiL,"mPhiL/D");
  hLambda -> Branch("mLambdaMass",&mLambdaMass,"mLambdaMass/D");
  hLambda -> Branch("mEtaL",&mEtaL,"mEtaL/D");
  hLambda -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");

  hALambda -> Branch("mVzLB",&mVzLB,"mVzLB/D");
  hALambda -> Branch("mDcaV0ToPVLB",&mDcaV0ToPVLB,"mDcaV0ToPVLB/D");
  hALambda -> Branch("mDecayLengthLB",&mDecayLengthLB,"mDecayLengthLB/D");
  hALambda -> Branch("mDcaDaughtersLB",&mDcaDaughtersLB,"mDcaDaughtersLB/D");
  hALambda -> Branch("mMomentumLB",&mMomentumLB,"mMomentumLB/D");
  hALambda -> Branch("mPtLB",&mPtLB,"mPtLB/D");
  hALambda -> Branch("mV0TypeLB",&mV0TypeLB,"mV0TypeLB/D");
  hALambda -> Branch("mNSigmaPLB",&mNSigmaPLB,"mNSigmaPLB/D");
  hALambda -> Branch("mNSigmaPiLB",&mNSigmaPiLB,"mNSigmaPiLB/D");
  hALambda -> Branch("mMass2PLB",&mMass2PLB,"mMass2PLB/D");
  hALambda -> Branch("mMass2PiLB",&mMass2PiLB,"mMass2PiLB/D");
  hALambda -> Branch("mDcaPToPVLB",&mDcaPToPVLB,"mDcaPToPVLB/D");
  hALambda -> Branch("mDcaPiToPVLB",&mDcaPiToPVLB,"mDcaPiToPVLB/D");
  hALambda -> Branch("mNHitsPLB",&mNHitsPLB,"mNHitsPLB/D");
  hALambda -> Branch("mNHitsPiLB",&mNHitsPiLB,"mNHitsPiLB/D");
  hALambda -> Branch("mHitsPossPLB",&mHitsPossPLB,"mHitsPossPLB/D");
  hALambda -> Branch("mHitsPossPiLB",&mHitsPossPiLB,"mHitsPossPiLB/D");
  hALambda -> Branch("mRapidityLB",&mRapidityLB,"mRapidityLB/D");
  hALambda -> Branch("mPhiPsiLB",&mPhiPsiLB,"mPhiPsiLB/D");
  hALambda -> Branch("mPhiLB",&mPhiLB,"mPhiLB/D");
  hALambda -> Branch("mLBarMass",&mLBarMass,"mLBarMass/D");
  hALambda -> Branch("mEtaLB",&mEtaLB,"mEtaLB/D");
  hALambda -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");

  hKshort -> Branch("mVzK",&mVzK,"mVzK/D");
  hKshort -> Branch("mDcaV0ToPVK",&mDcaV0ToPVK,"mDcaV0ToPVK/D");
  hKshort -> Branch("mDecayLengthK",&mDecayLengthK,"mDecayLengthK/D");
  hKshort -> Branch("mDcaDaughtersK",&mDcaDaughtersK,"mDcaDaughtersK/D");
  hKshort -> Branch("mMomentumK",&mMomentumK,"mMomentumK/D");
  hKshort -> Branch("mPtK",&mPtK,"mPtK/D");
  hKshort -> Branch("mV0TypeK",&mV0TypeK,"mV0TypeK/D");
  hKshort -> Branch("mNSigmaPiPK",&mNSigmaPiPK,"mNSigmaPiPK/D");
  hKshort -> Branch("mNSigmaPiNK",&mNSigmaPiNK,"mNSigmaPiNK/D");
  hKshort -> Branch("mMass2PiPK",&mMass2PiPK,"mMass2PiPK/D");
  hKshort -> Branch("mMass2PiNK",&mMass2PiNK,"mMass2PiNK/D");
  hKshort -> Branch("mDcaPiPToPVK",&mDcaPiPToPVK,"mDcaPiPToPVK/D");
  hKshort -> Branch("mDcaPiNToPVK",&mDcaPiNToPVK,"mDcaPiNToPVK/D");
  hKshort -> Branch("mNHitsPiPK",&mNHitsPiPK,"mNHitsPiPK/D");
  hKshort -> Branch("mNHitsPiNK",&mNHitsPiNK,"mNHitsPiNK/D");
  hKshort -> Branch("mHitsPossPiPK",&mHitsPossPiPK,"mHitsPossPiPK/D");
  hKshort -> Branch("mHitsPossPiNK",&mHitsPossPiNK,"mHitsPossPiNK/D");
  hKshort -> Branch("mRapidityK",&mRapidityK,"mRapidityK/D");
  hKshort -> Branch("mPhiPsiK",&mPhiPsiK,"mPhiPsiK/D");
  hKshort -> Branch("mPhiK",&mPhiK,"mPhiK/D");
  hKshort -> Branch("mKshortMass",&mKshortMass,"mKshortMass/D");
  hKshort -> Branch("mEtaK",&mEtaK,"mEtaK/D");
  hKshort -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");

  hLambdaBG -> Branch("mVzL",&mVzL,"mVzL/D");
  hLambdaBG -> Branch("mDcaV0ToPVL",&mDcaV0ToPVL,"mDcaV0ToPVL/D");
  hLambdaBG -> Branch("mDecayLengthL",&mDecayLengthL,"mDecayLengthL/D");
  hLambdaBG -> Branch("mDcaDaughtersL",&mDcaDaughtersL,"mDcaDaughtersL/D");
  hLambdaBG -> Branch("mMomentumL",&mMomentumL,"mMomentumL/D");
  hLambdaBG -> Branch("mPtL",&mPtL,"mPtL/D");
  hLambdaBG -> Branch("mV0TypeL",&mV0TypeL,"mV0TypeL/D");
  hLambdaBG -> Branch("mNSigmaPL",&mNSigmaPL,"mNSigmaPL/D");
  hLambdaBG -> Branch("mNSigmaPiL",&mNSigmaPiL,"mNSigmaPiL/D");
  hLambdaBG -> Branch("mMass2PL",&mMass2PL,"mMass2PL/D");
  hLambdaBG -> Branch("mMass2PiL",&mMass2PiL,"mMass2PiL/D");
  hLambdaBG -> Branch("mDcaPToPVL",&mDcaPToPVL,"mDcaPToPVL/D");
  hLambdaBG -> Branch("mDcaPiToPVL",&mDcaPiToPVL,"mDcaPiToPVL/D");
  hLambdaBG -> Branch("mNHitsPL",&mNHitsPL,"mNHitsPL/D");
  hLambdaBG -> Branch("mNHitsPiL",&mNHitsPiL,"mNHitsPiL/D");
  hLambdaBG -> Branch("mHitsPossPL",&mHitsPossPL,"mHitsPossPL/D");
  hLambdaBG -> Branch("mHitsPossPiL",&mHitsPossPiL,"mHitsPossPiL/D");
  hLambdaBG -> Branch("mRapidityL",&mRapidityL,"mRapidityL/D");
  hLambdaBG -> Branch("mPhiPsiL",&mPhiPsiL,"mPhiPsiL/D");
  hLambdaBG -> Branch("mPhiL",&mPhiL,"mPhiL/D");
  hLambdaBG -> Branch("mLambdaMass",&mLambdaMass,"mLambdaMass/D");
  hLambdaBG -> Branch("mEtaL",&mEtaL,"mEtaL/D");
  hLambdaBG -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");

  hALambdaBG -> Branch("mVzLB",&mVzLB,"mVzLB/D");
  hALambdaBG -> Branch("mDcaV0ToPVLB",&mDcaV0ToPVLB,"mDcaV0ToPVLB/D");
  hALambdaBG -> Branch("mDecayLengthLB",&mDecayLengthLB,"mDecayLengthLB/D");
  hALambdaBG -> Branch("mDcaDaughtersLB",&mDcaDaughtersLB,"mDcaDaughtersLB/D");
  hALambdaBG -> Branch("mMomentumLB",&mMomentumLB,"mMomentumLB/D");
  hALambdaBG -> Branch("mPtLB",&mPtLB,"mPtLB/D");
  hALambdaBG -> Branch("mV0TypeLB",&mV0TypeLB,"mV0TypeLB/D");
  hALambdaBG -> Branch("mNSigmaPLB",&mNSigmaPLB,"mNSigmaPLB/D");
  hALambdaBG -> Branch("mNSigmaPiLB",&mNSigmaPiLB,"mNSigmaPiLB/D");
  hALambdaBG -> Branch("mMass2PLB",&mMass2PLB,"mMass2PLB/D");
  hALambdaBG -> Branch("mMass2PiLB",&mMass2PiLB,"mMass2PiLB/D");
  hALambdaBG -> Branch("mDcaPToPVLB",&mDcaPToPVLB,"mDcaPToPVLB/D");
  hALambdaBG -> Branch("mDcaPiToPVLB",&mDcaPiToPVLB,"mDcaPiToPVLB/D");
  hALambdaBG -> Branch("mNHitsPLB",&mNHitsPLB,"mNHitsPLB/D");
  hALambdaBG -> Branch("mNHitsPiLB",&mNHitsPiLB,"mNHitsPiLB/D");
  hALambdaBG -> Branch("mHitsPossPLB",&mHitsPossPLB,"mHitsPossPLB/D");
  hALambdaBG -> Branch("mHitsPossPiLB",&mHitsPossPiLB,"mHitsPossPiLB/D");
  hALambdaBG -> Branch("mRapidityLB",&mRapidityLB,"mRapidityLB/D");
  hALambdaBG -> Branch("mPhiPsiLB",&mPhiPsiLB,"mPhiPsiLB/D");
  hALambdaBG -> Branch("mPhiLB",&mPhiLB,"mPhiLB/D");
  hALambdaBG -> Branch("mLBarMass",&mLBarMass,"mLBarMass/D");
  hALambdaBG -> Branch("mEtaLB",&mEtaLB,"mEtaLB/D");
  hALambdaBG -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");

  hKshortBG -> Branch("mVzK",&mVzK,"mVzK/D");
  hKshortBG -> Branch("mDcaV0ToPVK",&mDcaV0ToPVK,"mDcaV0ToPVK/D");
  hKshortBG -> Branch("mDecayLengthK",&mDecayLengthK,"mDecayLengthK/D");
  hKshortBG -> Branch("mDcaDaughtersK",&mDcaDaughtersK,"mDcaDaughtersK/D");
  hKshortBG -> Branch("mMomentumK",&mMomentumK,"mMomentumK/D");
  hKshortBG -> Branch("mPtK",&mPtK,"mPtK/D");
  hKshortBG -> Branch("mV0TypeK",&mV0TypeK,"mV0TypeK/D");
  hKshortBG -> Branch("mNSigmaPiPK",&mNSigmaPiPK,"mNSigmaPiPK/D");
  hKshortBG -> Branch("mNSigmaPiNK",&mNSigmaPiNK,"mNSigmaPiNK/D");
  hKshortBG -> Branch("mMass2PiPK",&mMass2PiPK,"mMass2PiPK/D");
  hKshortBG -> Branch("mMass2PiNK",&mMass2PiNK,"mMass2PiNK/D");
  hKshortBG -> Branch("mDcaPiPToPVK",&mDcaPiPToPVK,"mDcaPiPToPVK/D");
  hKshortBG -> Branch("mDcaPiNToPVK",&mDcaPiNToPVK,"mDcaPiNToPVK/D");
  hKshortBG -> Branch("mNHitsPiPK",&mNHitsPiPK,"mNHitsPiPK/D");
  hKshortBG -> Branch("mNHitsPiNK",&mNHitsPiNK,"mNHitsPiNK/D");
  hKshortBG -> Branch("mHitsPossPiPK",&mHitsPossPiPK,"mHitsPossPiPK/D");
  hKshortBG -> Branch("mHitsPossPiNK",&mHitsPossPiNK,"mHitsPossPiNK/D");
  hKshortBG -> Branch("mRapidityK",&mRapidityK,"mRapidityK/D");
  hKshortBG -> Branch("mPhiPsiK",&mPhiPsiK,"mPhiPsiK/D");
  hKshortBG -> Branch("mPhiK",&mPhiK,"mPhiK/D");
  hKshortBG -> Branch("mKshortMass",&mKshortMass,"mKshortMass/D");
  hKshortBG -> Branch("mEtaK",&mEtaK,"mEtaK/D");
  hKshortBG -> Branch("mCentralityBin",&mCentralityBin,"mCentralityBin/s");

}



void StStrangeMaker :: WriteHist(){

  char buff[255];
  sprintf(buff,"ControlHistos.%s.root",mJobIdname.Data());
  controlHitos1 = new TFile(buff,"recreate");
  controlHitos1 -> cd() ; 
  hStatsEvents  -> Write();

  for(int i=0;i<9;i++){
    hBbcEventPlaneES[i]  -> Write();
    hBbcEventPlaneWS[i]  -> Write();
    hBbcEventPlaneFS[i]  -> Write();
    hBbcEventPlaneE[i]   -> Write();
    hBbcEventPlaneW[i]   -> Write();
    hBbcEventPlaneF[i]   -> Write();
  }
  pBBC_shifteast_sin  -> Write();
  pBBC_shifteast_cos  -> Write();
  pBBC_shiftwest_sin  -> Write();
  pBBC_shiftwest_cos  -> Write();
  pBBC_shiftfull_sin  -> Write();
  pBBC_shiftfull_cos  -> Write();
  pRPCorrelation      -> Write();

  controlHitos1    -> Close();


  sprintf(buff,"v1TrackP.%s.root",mJobIdname.Data());
  v1TrackHistosP = new TFile(buff,"recreate");
  v1TrackHistosP->cd();
  hProton -> Write();
  v1TrackHistosP->Close();

  sprintf(buff,"v1TrackAP.%s.root",mJobIdname.Data());
  v1TrackHistosAP = new TFile(buff,"recreate");
  v1TrackHistosAP->cd();
  hAProton -> Write();
  v1TrackHistosAP->Close();

  sprintf(buff,"v1TrackPiP.%s.root",mJobIdname.Data());
  v1TrackHistosPiP = new TFile(buff,"recreate");
  v1TrackHistosPiP->cd();
  hPionPos -> Write();
  v1TrackHistosPiP->Close();

  sprintf(buff,"v1TrackPiN.%s.root",mJobIdname.Data());
  v1TrackHistosPiN = new TFile(buff,"recreate");
  v1TrackHistosPiN->cd();
  hPionNeg -> Write();
  v1TrackHistosPiN->Close();

  sprintf(buff,"v1TrackKP.%s.root",mJobIdname.Data());
  v1TrackHistosKP = new TFile(buff,"recreate");
  v1TrackHistosKP->cd();
  hKaonPos -> Write();
  v1TrackHistosKP->Close();

  sprintf(buff,"v1TrackKN.%s.root",mJobIdname.Data());
  v1TrackHistosKN = new TFile(buff,"recreate");
  v1TrackHistosKN->cd();
  hKaonNeg -> Write();
  v1TrackHistosKN->Close();

  sprintf(buff,"StrangeHistosL.%s.root",mJobIdname.Data());
  v1StrangeHistosL = new TFile(buff,"recreate");
  v1StrangeHistosL->cd();
  hLambda    -> Write(); 
  v1StrangeHistosL->Close();

  sprintf(buff,"StrangeHistosLBG.%s.root",mJobIdname.Data());
  v1StrangeHistosLBG = new TFile(buff,"recreate");
  v1StrangeHistosLBG->cd();
  hLambdaBG  -> Write(); 
  v1StrangeHistosLBG->Close();

  sprintf(buff,"StrangeHistosAL.%s.root",mJobIdname.Data());
  v1StrangeHistosAL = new TFile(buff,"recreate");
  v1StrangeHistosAL->cd();
  hALambda   -> Write(); 
  v1StrangeHistosAL->Close();

  sprintf(buff,"StrangeHistosALBG.%s.root",mJobIdname.Data());
  v1StrangeHistosALBG = new TFile(buff,"recreate");
  v1StrangeHistosALBG->cd();
  hALambdaBG -> Write(); 
  v1StrangeHistosALBG->Close();

  sprintf(buff,"StrangeHistosK.%s.root",mJobIdname.Data());
  v1StrangeHistosK = new TFile(buff,"recreate");
  v1StrangeHistosK->cd();
  hKshort    -> Write(); 
  v1StrangeHistosK->Close();

  sprintf(buff,"StrangeHistosKBG.%s.root",mJobIdname.Data());
  v1StrangeHistosKBG = new TFile(buff,"recreate");
  v1StrangeHistosKBG->cd();
  hKshortBG  -> Write(); 
  v1StrangeHistosKBG->Close();

}







