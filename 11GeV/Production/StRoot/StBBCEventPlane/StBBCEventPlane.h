#ifndef StBBCEventPlane_hh
#define StBBCEventPlane_hh

#include "TObject.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StThreeVectorF.hh"


class StMuTrack;
class StMuEvent;


class StBBCEventPlane : public TObject {
	

	public:
		StBBCEventPlane(StMuEvent*, int);
		~StBBCEventPlane();
		Double_t* BBC_GetPsi();//
	
	private:
	    StBbcTriggerDetector mBBC;
	    int mCentralityID;
		Double_t mBbc_phi;
		Double_t mPsi_full;
		Double_t mPsi_e;
		Double_t mPsi_w;	
		Double_t BBC_GetPhi(int,int);// 
	    Double_t mPsi[6];  // return value of get phi
	
		float  mBBCshiftEast_cos[20];
		float  mBBCshiftEast_sin[20];
		float  mBBCshiftWest_cos[20];
		float  mBBCshiftWest_sin[20];
		float  mBBCshiftFull_cos[20];
		float  mBBCshiftFull_sin[20];
	
	
	

ClassDef(StBBCEventPlane,1)
};
#endif
