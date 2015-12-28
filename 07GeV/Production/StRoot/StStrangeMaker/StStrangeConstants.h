
#ifndef StStrangeConstants_h
#define StStrangeConstants_h

#include "Rtypes.h"  // check what is this, with out don't work

enum {
	nTrigger = 2,
	pion =0,
	proton=2,
	pos=0,
	neg=1,
	ks=0,
	lambda=1

     };


class Constant {
public:
	static Float_t mMass[4];
	static Float_t mMassV0[2];

	
	static Int_t R_TPC;
	static Int_t R_TPC2;
	static Int_t R2_TPC;
	static Int_t Z_TPC;
	static Float_t d_tolerance;
	
	
	// event selectioin
	static Int_t mTriggerId[nTrigger];
	static Float_t mVzMax;
	static Float_t mVrMax;
	static Int_t mRefMultMin;
	static Int_t mRefMultMax;
	static Int_t mBTofTrayMultMin;
//	static Int_t mBTofTrayMultMax;
	static Int_t mBBCSumMin;
	static unsigned int mNoPrimaryVertMax;
	static unsigned int mNoPrimaryTracksMax;
	static unsigned int mNoPrimaryTracksMin;

	
	// track selection
	static Float_t mPtMin;
	static Float_t mPtMax;
	static Int_t mNHitsFitMin;
	static Float_t mRatioMin;
	static Float_t mProtonMassWindow;
	static Float_t mPionMassWindow;
//	static Int_t   mEtaMax; 
//	static Int_t   mEtaMin;

	static Float_t mTrackNSigmaPionMax;
	static Float_t mTrackNSigmaKaonMax;
	static Float_t mTrackNSigmaProtonMax;

	static Float_t mProtonMomMax;
	static Float_t mPionMomMax;
	static Float_t mKaonMomMax;

	static Float_t mProtonPtMax;
//	static Float_t mPionPtMax;
//	static Float_t mKaonPtMax;

	static Float_t mProtonPtMin;
	static Float_t mPionPtMin;
	static Float_t mKaonPtMin;

	static Float_t mV0NSigmaPionMax;
	static Float_t mV0NSigmaProtonMax;

	static Float_t mTrackDcaGlobal;

	static Float_t mTrackMass2ProtonMin;
	static Float_t mTrackMass2ProtonMax;
	static Float_t mTrackMass2PionMin;
	static Float_t mTrackMass2PionMax;
	static Float_t mTrackMass2KaonMin;
	static Float_t mTrackMass2KaonMax;



	// Kshort
	static Float_t mKsDcaDaughtersMax;
	static Float_t mKsPionDca2VertexMin;
	static Float_t mKsDca2VertexMax;
	static Float_t mKsDecayLengthMin;
	static Float_t mKsDecayLengthMax;
	static Float_t mKsMassWindow;
	
	// Lambda
	static Float_t mLambdaDcaDaughtersMax;
	static Float_t mLambdaDecayLengthMax;
	static Float_t mLambdaMassWindow;

	// 0: not tof, 1:pion tof 2:proton tof 3: both tof

	static Float_t mLambdaPionDca2VertexMin_0;
	static Float_t mLambdaProtonDca2VertexMin_0;
	static Float_t mLambdaDca2VertexMax_0;
	static Float_t mLambdaDecayLengthMin_0;

	static Float_t mLambdaPionDca2VertexMin_1;
	static Float_t mLambdaProtonDca2VertexMin_1;
	static Float_t mLambdaDca2VertexMax_1;
	static Float_t mLambdaDecayLengthMin_1;

	static Float_t mLambdaPionDca2VertexMin_2;
	static Float_t mLambdaProtonDca2VertexMin_2;
	static Float_t mLambdaDca2VertexMax_2;
	static Float_t mLambdaDecayLengthMin_2;

	static Float_t mLambdaPionDca2VertexMin_3;
	static Float_t mLambdaProtonDca2VertexMin_3;
	static Float_t mLambdaDca2VertexMax_3;
	static Float_t mLambdaDecayLengthMin_3;

	static Float_t mV0Mass2ProtonMin;
	static Float_t mV0Mass2ProtonMax;
	static Float_t mV0Mass2PionMin;
	static Float_t mV0Mass2PionMax;



ClassDef(Constant, 1)
};              
#endif
