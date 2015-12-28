
#include "StStrangeConstants.h"
#include "TMath.h"

ClassImp(Constant)

Float_t Constant::mMass[4]   = {0.13957, 0.49368, 0.93827, 0.000511};
Float_t Constant::mMassV0[2] = {0.49765, 1.11568};

Int_t Constant :: R_TPC  = 200;
Int_t Constant :: R_TPC2 = 40000;
Int_t Constant :: R2_TPC = 30000;
Int_t Constant :: Z_TPC  = 180;
Float_t Constant :: d_tolerance =0.85;



Int_t Constant::mTriggerId[nTrigger] = {
  310004, // mb
  310014
};  


Float_t Constant::mVzMax             = 60; // ***50
Float_t Constant::mVrMax             = 2;
Int_t Constant::mRefMultMin          = 6;
Int_t Constant::mRefMultMax          = 1000;
Int_t Constant:: mBTofTrayMultMin    = 5;
Int_t Constant:: mBBCSumMin          = 15;
unsigned int Constant:: mNoPrimaryVertMax   = 5;
unsigned int Constant:: mNoPrimaryTracksMin = 0;
unsigned int Constant:: mNoPrimaryTracksMax = 15000;

// track selection
Float_t Constant::mPtMin       = 0.1;
Float_t Constant::mPtMax       = 5.0;
Int_t   Constant::mNHitsFitMin = 12; // ***15               
Float_t Constant::mRatioMin    = 0.0; // *** 0.52
Float_t Constant::mProtonMassWindow = 0.1; // not used
Float_t Constant::mPionMassWindow = 0.1;// not used







// v1 Track cuts
Float_t Constant::mTrackNSigmaPionMax   = 2.4; // *** 2.0
Float_t Constant::mTrackNSigmaKaonMax   = 2.4; // *** 2.0
Float_t Constant::mTrackNSigmaProtonMax = 2.4; // *** 2.0

Float_t Constant::mProtonMomMax  = 3.36;//*** 2.8 
Float_t Constant::mPionMomMax    = 1.92;//*** 1.6
Float_t Constant::mKaonMomMax    = 1.92;//*** 1.6

Float_t Constant::mProtonPtMax    = 2.4;//*** 2.0 
//Float_t Constant::mPionPtMax    = 5.0;//*** 1.6
//Float_t Constant::mKaonPtMax    = 5.0;//*** 1.6

Float_t Constant::mProtonPtMin  = 0.32;//*** 0.4 
Float_t Constant::mPionPtMin    = 0.2;//*** 0.2
Float_t Constant::mKaonPtMin    = 0.2;//*** 0.2

Float_t Constant::mTrackMass2ProtonMin = 0.72;
Float_t Constant::mTrackMass2ProtonMax = 1.1;

Float_t Constant::mTrackMass2PionMin   = -0.15; 
Float_t Constant::mTrackMass2PionMax   = 0.15;

Float_t Constant::mTrackMass2KaonMin   = 0.18;
Float_t Constant::mTrackMass2KaonMax   = 0.40;

Float_t Constant::mTrackDcaGlobal = 3.6;







// v0 cuts
Float_t Constant::mV0NSigmaPionMax       = 3.25; // *** 3.0
Float_t Constant::mV0NSigmaProtonMax     = 3.25; // *** 3.0

// Kshort
Float_t Constant::mKsDcaDaughtersMax   = 1.20; //*** 1.00
Float_t Constant::mKsDecayLengthMax    = 150.0; // not used
Float_t Constant::mKsMassWindow        = 0.1; // to plot
Float_t Constant::mKsPionDca2VertexMin = 0.56; //*** 0.7
Float_t Constant::mKsDca2VertexMax     = 0.96; //*** 0.8 
Float_t Constant::mKsDecayLengthMin    = 2.4; //*** 3.0 



// Lambda
Float_t Constant::mLambdaDcaDaughtersMax     = 1.20; //*** 1.0
Float_t Constant::mLambdaDecayLengthMax      = 150.0;// not used
Float_t Constant::mLambdaMassWindow          = 0.1;

// no tof
Float_t Constant::mLambdaProtonDca2VertexMin_0 = 0.48; //*** 0.6
Float_t Constant::mLambdaPionDca2VertexMin_0   = 1.36; //*** 1.7
Float_t Constant::mLambdaDca2VertexMax_0       = 0.9; //*** 0.75
Float_t Constant::mLambdaDecayLengthMin_0      = 3.2; //*** 4.0
// Pion tof
Float_t Constant::mLambdaProtonDca2VertexMin_1 = 0.4; //*** 0.5
Float_t Constant::mLambdaPionDca2VertexMin_1   = 1.2; //*** 1.5
Float_t Constant::mLambdaDca2VertexMax_1       = 0.9; //*** 0.75
Float_t Constant::mLambdaDecayLengthMin_1      = 2.8; //*** 3.5
// Proton tof
Float_t Constant::mLambdaProtonDca2VertexMin_2 = 0.12; //*** 0.15
Float_t Constant::mLambdaPionDca2VertexMin_2   = 0.64; //*** 0.8
Float_t Constant::mLambdaDca2VertexMax_2       = 1.44; //*** 1.2
Float_t Constant::mLambdaDecayLengthMin_2      = 2.0; //*** 2.5
//both tof
Float_t Constant::mLambdaProtonDca2VertexMin_3 = 0.08; //*** 0.1
Float_t Constant::mLambdaPionDca2VertexMin_3   = 0.56; //*** 0.7
Float_t Constant::mLambdaDca2VertexMax_3       = 1.56; //*** 1.3
Float_t Constant::mLambdaDecayLengthMin_3      = 1.60; //*** 2.0


Float_t Constant::mV0Mass2ProtonMin = 0.5 ;
Float_t Constant::mV0Mass2ProtonMax = 1.5 ;
Float_t Constant::mV0Mass2PionMin =  -999; //defined StRoot/StStrangeMaker/StStrangePIDTracks.cxx
Float_t Constant::mV0Mass2PionMax = 0.1;

