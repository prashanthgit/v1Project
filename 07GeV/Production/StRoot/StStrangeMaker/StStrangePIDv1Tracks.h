
#ifndef StStrangePIDv1Tracks_h
#define StStrangePIDv1Tracks_h

#include "TObject.h"

class StMuEvent;
class StMuTrack;
class StMuDstMaker;

class StStrangePIDv1Tracks : public TObject {
public:
	StStrangePIDv1Tracks();
	~StStrangePIDv1Tracks();
	
	Bool_t Proton( StMuTrack *);
	Bool_t Pion( StMuTrack *);
	Bool_t Kaon( StMuTrack *);
	Double_t Mass2Track( StMuTrack *);
	Bool_t ToFTrack(StMuTrack *);
	
private:
	
	ClassDef(StStrangePIDv1Tracks,1)
};

#endif
