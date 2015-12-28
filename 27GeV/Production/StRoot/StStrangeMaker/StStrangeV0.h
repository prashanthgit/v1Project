#ifndef StStrangeV0_hh
#define StStrangeV0_hh

#include "TObject.h"
#include "StThreeVectorF.hh"
#include "StThreeVectorD.hh"
#include "StPhysicalHelixD.hh"
#include "StTrackHelix.hh"


class StMuTrack;
class StMuEvent;


class StStrangeV0 : public TObject {
public:
	~StStrangeV0();
	StStrangeV0(StPhysicalHelixD, StPhysicalHelixD, StThreeVectorF,StThreeVectorF,float, int); 

	StThreeVectorF momentumTrack(const Int_t i) const;
	StThreeVectorD v0Position()const { return mxV0;}
	Float_t        dcaDaughters() const { return (Float_t)mDcaDaughters; }
	Float_t        dcaV02PV() const { return (Float_t)mDcaV0toPV; }
	Float_t        mass() const { return mM; } // want give you mass unless setParticleHypothesis is defined.
	Float_t        dcaPos() const {return mDcaPos;}
	Float_t        dcaNeg() const {return mDcaNeg;}
	Float_t        decayLength() const {return mV0DecayLength;};
	StThreeVectorF momentum() const;
	Double_t       energy(const Int_t, const Int_t);
	void           setParticleHypothesis(const Int_t ip_pos, const Int_t ip_neg);

	
private:
	
	double closestDistance(StPhysicalHelixD , StPhysicalHelixD , double magnet, const StThreeVectorF&,const StThreeVectorF& , StThreeVectorF& xv0, StThreeVectorF& op1, StThreeVectorF& op2, int);
	void   dcaPToPi(float *fi_p,float *fi_pi, StTrackHelix* proton,StTrackHelix* pion,float *d_root,float *v0,float alfa);

    
protected:
	
	StThreeVectorF mMomentum[2];
	StThreeVectorD mxV0;
	Float_t        mDcaDaughters;
	Float_t        mDcaV0toPV;
	Float_t        mM;
	Float_t        mV0DecayLength;
	Float_t		   mDcaPos;
	Float_t		   mDcaNeg;

	ClassDef(StStrangeV0, 1)
};

#endif
