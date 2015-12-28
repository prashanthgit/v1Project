
#ifndef StStrangeCut_h
#define StStrangeCut_h

#include "TObject.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"

class StMuEvent;
class StMuTrack;
class StStrangeV0;
class StMuDstMaker;

class StStrangeCut : public TObject {
public:
	StStrangeCut();
	~StStrangeCut();
	
	bool passEvent( StMuEvent *,StMuDstMaker*  );
	bool passTrack( StMuTrack * );
	bool     passKs( StStrangeV0 *,StThreeVectorF,StThreeVectorF, StThreeVectorF );
	bool passLambda( StStrangeV0 *,StThreeVectorF,StThreeVectorF, StThreeVectorF, UShort_t);
	bool   passLbar( StStrangeV0 *,StThreeVectorF,StThreeVectorF, StThreeVectorF, UShort_t);

private:
	StBbcTriggerDetector mBBC;
	bool etaSymmetryCut(StMuDstMaker*);	

	ClassDef(StStrangeCut,1)
};

#endif
