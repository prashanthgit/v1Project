
#include "StStrangeCut.h"
#include "StStrangeV0.h"
#include "StStrangeConstants.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "TMath.h"
#include "StLorentzVectorF.hh"
#include "StLorentzVectorD.hh"


ClassImp(StStrangeCut)
StStrangeCut::StStrangeCut() {}
StStrangeCut::~StStrangeCut() {}


bool StStrangeCut::passEvent( StMuEvent *ev, StMuDstMaker *muDstMaker ) {
	
	StThreeVectorF pVertex = ev->eventSummary().primaryVertexPosition();
	if(fabs(pVertex.x())<1.e-5 && fabs(pVertex.y())<1.e-5 && fabs(pVertex.z())<1.e-5) return kFALSE;
	if(fabs(pVertex.z())>Constant::mVzMax) return kFALSE;
	
	const Float_t vx = pVertex.x() ;
	const Float_t vy = pVertex.y() ;
	if(TMath::Sqrt(vx*vx+vy*vy)>Constant::mVrMax) return kFALSE ;
	
	bool isTrg = kFALSE;
	for(int i=0;i<nTrigger;i++) {
		if(ev->triggerIdCollection().nominal().isTrigger(Constant::mTriggerId[i])){
			isTrg = kTRUE;
		}
	}
	if(!isTrg) return kFALSE;

	if(ev->refMult()<Constant::mRefMultMin) return kFALSE;
	if(ev->btofTrayMultiplicity()<Constant::mBTofTrayMultMin) return kFALSE;
	if(muDstMaker->muDst()->numberOfPrimaryVertices()>Constant::mNoPrimaryVertMax) return kFALSE;
	if(muDstMaker->muDst()->numberOfPrimaryTracks() < Constant::mNoPrimaryTracksMin || muDstMaker->muDst()->numberOfPrimaryTracks() > Constant::mNoPrimaryTracksMax) return kFALSE; 
	// for BBC  ADC saturation
	mBBC = ev->bbcTriggerDetector();
	if(mBBC.adcSumEast()< Constant::mBBCSumMin || mBBC.adcSumWest()< Constant::mBBCSumMin) return kFALSE;
	return kTRUE;
}

//----------------------------------------------------------------------------------

bool StStrangeCut::passTrack( StMuTrack *t ){ // ToF cut is done in main maker via StrangePIDTrack class
	if(!t) return kFALSE;
	if(TMath::Abs(t -> charge()) != 1) return kFALSE;
	if(t->type()!=global) return kFALSE;
	if(t->flag()<0||t->flag()>1000)	return kFALSE;
	if(t->bad() )kFALSE;
	if(t->p().perp()<Constant::mPtMin) return kFALSE;
	if(t->p().perp()>Constant::mPtMax) return kFALSE;
	if( t->nHitsFit(kTpcId) < Constant::mNHitsFitMin ) return kFALSE;
	if( ((1.0*t->nHitsFit(kTpcId))/(1.0*t->nHitsPoss(kTpcId))) < Constant::mRatioMin ) return kFALSE;
 	if(t->id()<0 || t->id()>=50000) return kFALSE; 	
	
	return kTRUE;
}


//----------------------------------------------------------------------------------


bool StStrangeCut::passKs( StStrangeV0 *v0,StThreeVectorF pVtx, StThreeVectorF dca1, StThreeVectorF dca2 ) {
  if(!v0) return kFALSE;
  v0->setParticleHypothesis(pion, pion);
  StThreeVectorF v0Mom = v0->momentumTrack(pos) + v0->momentumTrack(neg);
  if(v0Mom.dot(v0->v0Position() -pVtx)<=0.) return kFALSE;  // V0 going away from primary vertex or r dot p for v0. cut on it. should be larger than 0,
  if(dca1.mag()<Constant::mKsPionDca2VertexMin) return kFALSE;
  if(dca2.mag()<Constant::mKsPionDca2VertexMin) return kFALSE;
  if(v0->dcaV02PV()>Constant::mKsDca2VertexMax) return kFALSE;
  if(v0->decayLength()<Constant::mKsDecayLengthMin ) return kFALSE;	
  if(v0->dcaDaughters()>Constant::mKsDcaDaughtersMax) return kFALSE;
  if(fabs(v0->mass()-Constant::mMassV0[ks])>Constant::mKsMassWindow) return kFALSE;
  return kTRUE;
}

//----------------------------------------------------------------------------------

bool StStrangeCut::passLambda(StStrangeV0 *v0,StThreeVectorF pVtx, StThreeVectorF dca1, StThreeVectorF dca2, UShort_t v0Type ) {
  if(!v0) return kFALSE;
  v0->setParticleHypothesis(proton, pion);
  StThreeVectorF v0Mom = v0->momentumTrack(pos) + v0->momentumTrack(neg);
  if(v0Mom.dot(v0->v0Position() - pVtx)<=0.) return kFALSE;  // V0 going away from primary vertex or r dot p for v0. cut on it. should be larger than 0,
  if(v0->dcaDaughters()>Constant::mLambdaDcaDaughtersMax) return kFALSE;
  if(fabs(v0->mass()-Constant::mMassV0[lambda])>Constant::mLambdaMassWindow) return kFALSE;

  if(v0Type == 0){
    if(dca1.mag()<Constant::mLambdaProtonDca2VertexMin_0 ) return kFALSE;
    if(dca2.mag()<Constant::mLambdaPionDca2VertexMin_0 ) return kFALSE;  
    if(v0->dcaV02PV()>Constant::mLambdaDca2VertexMax_0 ) return kFALSE; 
    if(v0->decayLength()<Constant::mLambdaDecayLengthMin_0 ) return kFALSE;
  }else if(v0Type == 1){
    if(dca1.mag()<Constant::mLambdaProtonDca2VertexMin_1 ) return kFALSE;
    if(dca2.mag()<Constant::mLambdaPionDca2VertexMin_1 ) return kFALSE;  
    if(v0->dcaV02PV()>Constant::mLambdaDca2VertexMax_1 ) return kFALSE; 
    if(v0->decayLength()<Constant::mLambdaDecayLengthMin_1 ) return kFALSE;
  }else if(v0Type == 2){
    if(dca1.mag()<Constant::mLambdaProtonDca2VertexMin_2 ) return kFALSE;
    if(dca2.mag()<Constant::mLambdaPionDca2VertexMin_2 ) return kFALSE;  
    if(v0->dcaV02PV()>Constant::mLambdaDca2VertexMax_2 ) return kFALSE; 
    if(v0->decayLength()<Constant::mLambdaDecayLengthMin_2 ) return kFALSE;
  }else if(v0Type == 3){
    if(dca1.mag()<Constant::mLambdaProtonDca2VertexMin_3 ) return kFALSE;
    if(dca2.mag()<Constant::mLambdaPionDca2VertexMin_3 ) return kFALSE;  
    if(v0->dcaV02PV()>Constant::mLambdaDca2VertexMax_3 ) return kFALSE; 
    if(v0->decayLength()<Constant::mLambdaDecayLengthMin_3 ) return kFALSE;
  }else{
    return kFALSE;
  }

  return kTRUE;

}

//----------------------------------------------------------------------------------

bool StStrangeCut::passLbar( StStrangeV0 *v0,StThreeVectorF pVtx, StThreeVectorF dca1, StThreeVectorF dca2, UShort_t v0Type ) {
  if(!v0) return kFALSE;
  v0->setParticleHypothesis(pion, proton);
  StThreeVectorF v0Mom = v0->momentumTrack(pos) + v0->momentumTrack(neg);
  if(v0Mom.dot(v0->v0Position() -pVtx)<=0.) return kFALSE;  // V0 going away from primary vertex or r dot p for v0. cut on it. should be larger than 0,
  if(v0->dcaDaughters()>Constant::mLambdaDcaDaughtersMax) return kFALSE;
  if(fabs(v0->mass()-Constant::mMassV0[lambda])>Constant::mLambdaMassWindow) return kFALSE;

  if(v0Type == 0){
    if(dca1.mag()<Constant::mLambdaPionDca2VertexMin_0) return kFALSE;
    if(dca2.mag()<Constant::mLambdaProtonDca2VertexMin_0) return kFALSE;
    if(v0->dcaV02PV()>Constant::mLambdaDca2VertexMax_0) return kFALSE;
    if(v0->decayLength()<Constant::mLambdaDecayLengthMin_0) return kFALSE;
  }else if(v0Type == 1){
    if(dca1.mag()<Constant::mLambdaPionDca2VertexMin_1) return kFALSE;
    if(dca2.mag()<Constant::mLambdaProtonDca2VertexMin_1) return kFALSE;
    if(v0->dcaV02PV()>Constant::mLambdaDca2VertexMax_1) return kFALSE;
    if(v0->decayLength()<Constant::mLambdaDecayLengthMin_1) return kFALSE;
  }else if(v0Type == 2){
    if(dca1.mag()<Constant::mLambdaPionDca2VertexMin_2) return kFALSE;
    if(dca2.mag()<Constant::mLambdaProtonDca2VertexMin_2) return kFALSE;
    if(v0->dcaV02PV()>Constant::mLambdaDca2VertexMax_2) return kFALSE;
    if(v0->decayLength()<Constant::mLambdaDecayLengthMin_2) return kFALSE;
  }else if(v0Type == 3){
    if(dca1.mag()<Constant::mLambdaPionDca2VertexMin_3) return kFALSE;
    if(dca2.mag()<Constant::mLambdaProtonDca2VertexMin_3) return kFALSE;
    if(v0->dcaV02PV()>Constant::mLambdaDca2VertexMax_3) return kFALSE;
    if(v0->decayLength()<Constant::mLambdaDecayLengthMin_3) return kFALSE;
  }else{
    return kFALSE;
  }

  return kTRUE;
}


