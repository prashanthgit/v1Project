
#ifndef StStrangePIDTracks_h
#define StStrangePIDTracks_h

#include "TObject.h"

class StMuEvent;
class StMuTrack;
class StMuDstMaker;

class StStrangePIDTracks : public TObject {
public:
	StStrangePIDTracks();
	~StStrangePIDTracks();
	
	Bool_t Proton( StMuTrack *);
	Bool_t ProtonTPC( StMuTrack *);
	Bool_t ProtonToF( StMuTrack *);
	Bool_t Pion( StMuTrack *);
	Bool_t PionTPC( StMuTrack *);
	Bool_t PionToF( StMuTrack *);
	Double_t Mass2Track( StMuTrack *);
	Bool_t ToFTrack(StMuTrack *);
	
private:
	
	ClassDef(StStrangePIDTracks,1)
};

#endif
