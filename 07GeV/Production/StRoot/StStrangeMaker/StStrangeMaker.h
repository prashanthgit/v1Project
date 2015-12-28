#ifndef STAR_StStrangeMaker
#define STAR_StStrangeMaker

#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "StPhysicalHelixD.hh"
#include "StTrackHelix.hh"
#include "THnSparse.h"
#include "TTree.h"


class StMuDstMaker;
class StMuEvent;
class StMuTrack;

class StStrangeV0;
class StStrangeCut;
class StStrangePIDTracks;
class StStrangePIDv1Tracks;
class TProfile;
class TProfile2D;
class TH2D;


#if !defined(ST_NO_NAMESPACES)
using namespace std;
#endif

const static Int_t mMixedEvents = 5;
const static Int_t phiBinsMixed = 30;

class StStrangeMaker : public StMaker {

  public:
    StStrangeMaker(StMuDstMaker* maker);
    virtual ~StStrangeMaker();
    Int_t Init();
    Int_t Make();
    Int_t Finish();

    void SetEnergy(Float_t energy){mEnergy = energy;}
    void SetJobIdName(TString JobIdName){mJobIdname = JobIdName;}


  protected:

    Int_t  centrality(StMuEvent*);
    Bool_t isBadRun(StMuEvent*);
    void   InitHist();
    void   WriteHist();
    void   clearMixedEvents();
    void   clearPIDTracks();
    void   mixEvents();
    void   v1Strange(StThreeVectorF,float);
    void   v1Track(StMuTrack*,Bool_t,Bool_t,Bool_t);
    Int_t  PhiRpBinMixed(Double_t);
    Int_t  VzBin(Double_t);

  private:

    StMuDstMaker* mMuDstMaker;                      
    StStrangeCut* mStrangeCut;
    StStrangePIDTracks* mStrangePIDTracks;
    StStrangePIDv1Tracks* mStrangePIDv1Tracks;

    struct track__{
      StPhysicalHelixD helix;
      StThreeVectorF dca;
      UShort_t tof; // tof track 1 not 0
      Double_t nsigma;
      Double_t mass2;
      Double_t nhitsfit;
      Double_t hitsratio;
    };
    typedef std::vector<track__> vect_track__;

    vect_track__ mPosProton;
    vect_track__ mNegProton;
    vect_track__ mPosPion;
    vect_track__ mNegPion;
    vect_track__ mPosKaon;
    vect_track__ mNegKaon;

    vect_track__       mMixedEventProtonPos[9][10][30][mMixedEvents];
    vect_track__       mMixedEventProtonNeg[9][10][30][mMixedEvents];
    vect_track__         mMixedEventPionPos[9][10][30][mMixedEvents];
    vect_track__         mMixedEventPionNeg[9][10][30][mMixedEvents];

    StThreeVectorF mMixedEventPrimaryVertex[9][10][30][mMixedEvents];
    float          mMixedEventMagneticField[9][10][30][mMixedEvents];
    Int_t                mMixedEventCounter[9][10][30];


    TFile*  controlHitos1;
    TFile*  v1StrangeHistosL;
    TFile*  v1StrangeHistosLBG;
    TFile*  v1StrangeHistosAL;
    TFile*  v1StrangeHistosALBG;
    TFile*  v1StrangeHistosK;
    TFile*  v1StrangeHistosKBG;
    TFile*  v1TrackHistosP;
    TFile*  v1TrackHistosAP;
    TFile*  v1TrackHistosPiP;
    TFile*  v1TrackHistosPiN;
    TFile*  v1TrackHistosKP;
    TFile*  v1TrackHistosKN;


    Float_t mEnergy;
    TString mJobIdname;
    Int_t   mEvents;
    Int_t   mEventsProcessed;
    Float_t mBField;
    UInt_t  mRefMult;
    double  mVz;
    Int_t   mCentralityBin;
    Int_t   mNoOfGlobalTracks;
    Bool_t  mMixedEventFlag;
    UInt_t  mMixedEventNo;
    Int_t   mMixedEventVzBin;
    Int_t   mMixedEventPsiRpBin;

    Double_t* mPsi;
    Float_t mReactionPlaneE;
    Float_t mReactionPlaneW;
    Float_t mReactionPlaneF;
    Float_t mReactionPlaneES;
    Float_t mReactionPlaneWS;   
    Float_t mReactionPlaneFS;
    Float_t mRPCorrelation;

    // for track variables
    Double_t mV1Track;
    Double_t mVzT;
    Double_t mMomentumT;
    Double_t mPtT;
    Double_t mNSigmaT;
    Double_t mMass2T;
    Double_t mDcaT;
    Double_t mNHitsT;
    Double_t mHitsPossT;
    Double_t mRapidityT;
    Double_t mPhiPsiT;
    Double_t mPhiT;
    Double_t mEtaT;
    Double_t mChargeT;



    UShort_t mV0Type; // 3: Both Tof, 2: Proton Tof, 1: Pion Tof, 0: no tof ;

    //Lambda 


    Double_t mVzL;
    Double_t mDcaV0ToPVL;
    Double_t mDecayLengthL;
    Double_t mDcaDaughtersL;
    Double_t mMomentumL;
    Double_t mPtL;
    Double_t mV0TypeL;
    Double_t mNSigmaPL;
    Double_t mNSigmaPiL;
    Double_t mMass2PL;
    Double_t mMass2PiL;
    Double_t mDcaPToPVL;
    Double_t mDcaPiToPVL;
    Double_t mNHitsPL;
    Double_t mNHitsPiL;
    Double_t mHitsPossPL;
    Double_t mHitsPossPiL;
    Double_t mRapidityL;
    Double_t mPhiPsiL ;
    Double_t mPhiL ;
    Double_t mLambdaMass;
    Double_t mEtaL;

    //Anti-Lambda 
    Double_t mVzLB;
    Double_t mDcaV0ToPVLB;
    Double_t mDecayLengthLB;
    Double_t mDcaDaughtersLB;
    Double_t mMomentumLB;
    Double_t mPtLB;
    Double_t mV0TypeLB;
    Double_t mNSigmaPLB;
    Double_t mNSigmaPiLB;
    Double_t mMass2PLB;
    Double_t mMass2PiLB;
    Double_t mDcaPToPVLB;
    Double_t mDcaPiToPVLB;
    Double_t mNHitsPLB;
    Double_t mNHitsPiLB;
    Double_t mHitsPossPLB;
    Double_t mHitsPossPiLB;
    Double_t mRapidityLB;
    Double_t mPhiPsiLB ;
    Double_t mPhiLB ;
    Double_t mLBarMass;
    Double_t mEtaLB;



    //k-Short 
    Double_t mVzK;
    Double_t mDcaV0ToPVK;
    Double_t mDecayLengthK;
    Double_t mDcaDaughtersK;
    Double_t mMomentumK;
    Double_t mPtK;
    Double_t mV0TypeK;
    Double_t mNSigmaPiPK;
    Double_t mNSigmaPiNK;
    Double_t mMass2PiPK;
    Double_t mMass2PiNK;
    Double_t mDcaPiPToPVK;
    Double_t mDcaPiNToPVK;
    Double_t mNHitsPiPK;
    Double_t mNHitsPiNK;
    Double_t mHitsPossPiPK;
    Double_t mHitsPossPiNK;
    Double_t mRapidityK;
    Double_t mPhiPsiK ;
    Double_t mPhiK ;
    Double_t mKshortMass;
    Double_t mEtaK;



 

    // controle histos
    TH1D* hStatsEvents;
    
    // BBC Histos
    TProfile2D* pBBC_shifteast_sin;
    TProfile2D* pBBC_shifteast_cos;
    TProfile2D* pBBC_shiftwest_sin;
    TProfile2D* pBBC_shiftwest_cos;
    TProfile2D* pBBC_shiftfull_sin;
    TProfile2D* pBBC_shiftfull_cos;
    TProfile*   pRPCorrelation;

    TH1D* hBbcEventPlaneFS[9];
    TH1D* hBbcEventPlaneES[9];
    TH1D* hBbcEventPlaneWS[9];
    TH1D* hBbcEventPlaneF[9];
    TH1D* hBbcEventPlaneE[9];
    TH1D* hBbcEventPlaneW[9];

    // For Strange particles
    TTree* hLambda;
    TTree* hLambdaBG;
    TTree* hALambda;
    TTree* hALambdaBG;
    TTree* hKshort;
    TTree* hKshortBG;
    
    // For Tracks
    TTree* hProton;
    TTree* hAProton;
    TTree* hPionPos;
    TTree* hPionNeg;
    TTree* hKaonPos;
    TTree* hKaonNeg;

    Int_t nP;
    ClassDef(StStrangeMaker,1)
      
};
#endif
