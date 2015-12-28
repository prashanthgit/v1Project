
#include "StStrangePIDTracks.h"
#include "StStrangeConstants.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "TMath.h"
#include "StLorentzVectorF.hh"
#include "StLorentzVectorD.hh"


ClassImp(StStrangePIDTracks)


StStrangePIDTracks::StStrangePIDTracks() {}
StStrangePIDTracks::~StStrangePIDTracks() {}


Bool_t StStrangePIDTracks::Proton( StMuTrack *t ){
  if(!t) return kFALSE;
  Bool_t protonTPC = kFALSE;
  Bool_t protonToF = kFALSE;
  protonTPC = TMath::Abs(t->nSigmaProton()) <= Constant::mV0NSigmaProtonMax;
  Double_t beta = t->btofPidTraits().beta(); 
  Double_t trackP = t->p().mag();
  Double_t massP2  = trackP*trackP*((1.0/(beta*beta)) -1.0) ;
  if( (Constant::mV0Mass2ProtonMin < massP2)&& (massP2 < Constant::mV0Mass2ProtonMax)) protonToF = kTRUE;
  if(!(protonTPC && protonToF)) return kFALSE;
  return kTRUE;
}



Bool_t StStrangePIDTracks::Pion( StMuTrack *t){
  if(!t) return kFALSE;
  Bool_t pionTPC = kFALSE;
  Bool_t pionToF = kFALSE;
  pionTPC  = TMath::Abs(t->nSigmaPion()) <= Constant::mV0NSigmaPionMax;
  Double_t beta = t->btofPidTraits().beta(); 
  Double_t trackP = t->p().mag();
  Double_t massP2  = trackP*trackP*((1.0/(beta*beta)) -1.0) ;
  if( ((0.017-0.013*trackP) < massP2) && (massP2 < Constant::mV0Mass2PionMax) ) pionToF = kTRUE;
  if(!(pionTPC && pionToF)) return kFALSE;
  return kTRUE;
}


////////////////////////// ////////////////////////// //////////////////////////
////////////////////////// ////////////////////////// //////////////////////////

Bool_t StStrangePIDTracks::ProtonTPC( StMuTrack *t ){
  if(!t) return kFALSE;
  Bool_t protonTPC = kFALSE;
  protonTPC = TMath::Abs(t->nSigmaProton()) <= Constant::mV0NSigmaProtonMax;
  if(!protonTPC ) return kFALSE;
  return kTRUE;
}



Bool_t StStrangePIDTracks::PionTPC( StMuTrack *t ){
  if(!t) return kFALSE;
  Bool_t pionTPC = kFALSE;
  pionTPC  = TMath::Abs(t->nSigmaPion()) <= Constant::mV0NSigmaPionMax;
  if(!pionTPC ) return kFALSE;
  return kTRUE;
}




////////////////////////// ////////////////////////// //////////////////////////
////////////////////////// ////////////////////////// //////////////////////////

Bool_t StStrangePIDTracks::ProtonToF( StMuTrack *t ){
  if(!t) return kFALSE;
  Bool_t protonToF = kFALSE;
  Double_t beta = t->btofPidTraits().beta(); 
  Double_t trackP = t->p().mag();
  Double_t massP2  = trackP*trackP*((1.0/(beta*beta)) -1.0) ;
  if( (Constant::mV0Mass2ProtonMin < massP2)&& (massP2 < Constant::mV0Mass2ProtonMax)) protonToF = kTRUE;
  if( !protonToF) return kFALSE;
  return kTRUE;
}



Bool_t StStrangePIDTracks::PionToF( StMuTrack *t ){
  if(!t) return kFALSE;
  Bool_t pionToF = kFALSE;
  Double_t beta = t->btofPidTraits().beta(); 
  Double_t trackP = t->p().mag();
  Double_t massP2  = trackP*trackP*((1.0/(beta*beta)) -1.0) ;
  if( ((0.017-0.013*trackP) < massP2) && (massP2 < Constant::mV0Mass2PionMax) ) pionToF = kTRUE;
  if(!pionToF) return kFALSE;
  return kTRUE;
}


////////////////////////// ////////////////////////// //////////////////////////
////////////////////////// ////////////////////////// //////////////////////////
Double_t StStrangePIDTracks::Mass2Track( StMuTrack *t){
  if(!t) return kFALSE;

  Bool_t isTOF =0;
  Double_t beta = t->btofPidTraits().beta(); 
  Double_t tof  = t->btofPidTraits().timeOfFlight();
  Double_t matchFlag = t->btofPidTraits().matchFlag();
  if(matchFlag > 0 && tof !=0 && beta !=0 && beta != -999) isTOF =1;

  Double_t mass2 = -999;
  if(isTOF){
    Double_t beta = t->btofPidTraits().beta(); 
    Double_t trackP = t->p().mag();
    mass2  = trackP*trackP*((1.0/(beta*beta)) -1.0);
  }
  return mass2;
}

Bool_t StStrangePIDTracks::ToFTrack( StMuTrack *t){
  Double_t beta = t->btofPidTraits().beta(); 
  Double_t tof  = t->btofPidTraits().timeOfFlight();
  Double_t matchFlag = t->btofPidTraits().matchFlag();

  if(matchFlag > 0 && tof !=0 && beta !=0 && beta != -999) return 1;
  return 0;
}
