
#include "StStrangePIDv1Tracks.h"
#include "StStrangeConstants.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "TMath.h"
#include "StLorentzVectorF.hh"
#include "StLorentzVectorD.hh"


ClassImp(StStrangePIDv1Tracks)


StStrangePIDv1Tracks::StStrangePIDv1Tracks() {}
StStrangePIDv1Tracks::~StStrangePIDv1Tracks() {}


Bool_t StStrangePIDv1Tracks::Proton( StMuTrack *t ){

  if(!t) return kFALSE;

  Bool_t protonTPC = kFALSE;
  Bool_t protonToF = kFALSE;
  Bool_t pCut      = kFALSE;
  Bool_t pDCA      = kFALSE;

  protonTPC = TMath::Abs(t->nSigmaProton()) <= Constant::mTrackNSigmaProtonMax;

  Double_t beta = t->btofPidTraits().beta(); 
  Double_t trackP = t->p().mag();
  Double_t trackPt = t->p().perp();
  Double_t massP2  = trackP*trackP*((1.0/(beta*beta)) -1.0) ;
  if( ( Constant::mTrackMass2ProtonMin < massP2)&& (massP2 < Constant::mTrackMass2ProtonMax)) protonToF = kTRUE;

  if( (trackP <= Constant::mProtonMomMax) && (trackPt >= Constant::mProtonPtMin) && (trackPt <= Constant::mProtonPtMax)) pCut = kTRUE;

  if( t->dcaGlobal().mag() <= Constant::mTrackDcaGlobal) pDCA = kTRUE;

  if(!(protonTPC && protonToF && pCut && pDCA)) return kFALSE;
  return kTRUE;
}



Bool_t StStrangePIDv1Tracks::Pion( StMuTrack *t){

  if(!t) return kFALSE;

  Bool_t pionTPC = kFALSE;
  Bool_t pionToF = kFALSE;
  Bool_t pCut     = kFALSE;
  Bool_t pDCA      = kFALSE;

  pionTPC  = TMath::Abs(t->nSigmaPion()) <= Constant::mTrackNSigmaPionMax;

  Double_t beta = t->btofPidTraits().beta(); 
  Double_t trackP = t->p().mag();
  Double_t trackPt = t->p().perp();
  Double_t massP2  = trackP*trackP*((1.0/(beta*beta)) -1.0) ;
  if( ( Constant::mTrackMass2PionMin < massP2) && (massP2 < Constant::mTrackMass2PionMax) ) pionToF = kTRUE;

  if( (trackP <= Constant::mPionMomMax) && (trackPt >= Constant::mPionPtMin)) pCut = kTRUE; 

  if( t->dcaGlobal().mag() <= Constant::mTrackDcaGlobal) pDCA = kTRUE;

  if(!(pionTPC && pionToF && pCut && pDCA)) return kFALSE;
  return kTRUE;
}


Bool_t StStrangePIDv1Tracks::Kaon( StMuTrack *t){

  if(!t) return kFALSE;

  Bool_t kaonTPC = kFALSE;
  Bool_t kaonToF = kFALSE;
  Bool_t pCut    = kFALSE;
  Bool_t pDCA      = kFALSE;

  kaonTPC  = TMath::Abs(t->nSigmaKaon()) <= Constant::mTrackNSigmaKaonMax;

  Double_t beta = t->btofPidTraits().beta(); 
  Double_t trackP = t->p().mag();
  Double_t trackPt = t->p().perp();
  Double_t massP2  = trackP*trackP*((1.0/(beta*beta)) -1.0) ;
  if( ( Constant::mTrackMass2KaonMin < massP2) && (massP2 < Constant::mTrackMass2KaonMax) ) kaonToF = kTRUE;

  if( (trackP <= Constant::mKaonMomMax) && (trackPt >= Constant::mKaonPtMin) ) pCut = kTRUE;

  if( t->dcaGlobal().mag() <= Constant::mTrackDcaGlobal) pDCA = kTRUE;

  if(!(kaonTPC && kaonToF && pCut && pDCA)) return kFALSE;
  return kTRUE;
}


////////////////////////// ////////////////////////// //////////////////////////
////////////////////////// ////////////////////////// //////////////////////////
Double_t StStrangePIDv1Tracks::Mass2Track( StMuTrack *t){
  if(!t) return kFALSE;
  Double_t beta = t->btofPidTraits().beta(); 
  Double_t trackP = t->p().mag();
  Double_t mass  = trackP*trackP*((1.0/(beta*beta)) -1.0) ;
  return mass;
}

Bool_t StStrangePIDv1Tracks::ToFTrack( StMuTrack *t){
  Double_t beta = t->btofPidTraits().beta(); 
  Double_t tof  = t->btofPidTraits().timeOfFlight();
  Double_t matchFlag = t->btofPidTraits().matchFlag();

  if(matchFlag > 0 && tof !=0 && beta !=0 && beta != -999) return 1;
  return 0;
}
