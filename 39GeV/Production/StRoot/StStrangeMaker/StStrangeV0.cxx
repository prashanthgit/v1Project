#include "StStrangeV0.h"
#include "StStrangeConstants.h"
#include "TMath.h"

#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StPhysicalHelixD.hh"
#include "StLorentzVectorF.hh"

ClassImp(StStrangeV0)

StStrangeV0::~StStrangeV0(){}

StStrangeV0::StStrangeV0(StPhysicalHelixD h_pos, StPhysicalHelixD h_neg, StThreeVectorF pv1,StThreeVectorF pv2,float B, int la) {
			
		StThreeVectorF xv0;
		StThreeVectorF op1;
		StThreeVectorF op2;
		
		double pathLength_pos = h_pos.pathLength(pv1, true); // do scan periods. NOTE: the default is not scan.
		double pathLength_neg = h_neg.pathLength(pv2, true);
		
		mDcaPos = (h_pos.at(pathLength_pos) - pv1).mag();// which is equal to dcaGlobal()
		mDcaNeg = (h_neg.at(pathLength_neg) - pv2).mag(); // which is equal to dcaGlobal() 
		
		mDcaDaughters = closestDistance(h_pos, h_neg, B, pv1,pv2, xv0, op1, op2,la);
		mM = 0; // for in set particlehypothsis function
		mMomentum[pos] = op1; // momentum of postive track
		mMomentum[neg] = op2; // momentum of negative track
		
		mxV0 = xv0;
		StThreeVectorD pv0 = op1 + op2; // momentum of V0 
		StThreeVectorD xv0toPV = xv0 - pv1; // decay length vector
		mV0DecayLength = xv0toPV.mag(); // decay length of v0 
		
		double v0dotPV = xv0toPV.dot(pv0) ;// r dot p for v0. cut on it. should be larger than 0, done in passV0 in StStrangeCut
		double dcav0toPV2 = v0dotPV*v0dotPV/pv0.mag2();
		mDcaV0toPV  = TMath::Sqrt( xv0toPV.mag2() - dcav0toPV2); // dca between V0 and PV
	
}







double StStrangeV0::closestDistance(StPhysicalHelixD  h_pos, StPhysicalHelixD  h_neg, double magnet, const StThreeVectorF& pv1, const StThreeVectorF& pv2, StThreeVectorF& xv0, StThreeVectorF& op1, StThreeVectorF& op2, int la )
{
    double PI = TMath::Pi();
    float x[3],p1[3],p2[3], d;
    StTrackHelix posTrackobject, negTrackobject;
    StTrackHelix* posTrack = &posTrackobject;
    StTrackHelix* negTrack = &negTrackobject;

    //fill StTrackHelix
    //inner helix. good for dca to PV.
    int hrot = h_pos.h();
    float ttheta = h_pos.phase();
    StThreeVectorF origin = h_pos.origin();
    float Xc = h_pos.xcenter();//-pv1.x();
    float Yc = h_pos.ycenter();//-pv1.y();
    float Zc = origin.z();//-pv1.z();
    float r  = 1./h_pos.curvature();
    float Vz = r*tan(h_pos.dipAngle());
    float pt = (h_pos.momentum(magnet/10)* 1.e-13 ).perp();
    float Pz = (h_pos.momentum(magnet/10)* 1.e-13).z();


    posTrack -> Xc = Xc;
    posTrack -> Yc = Yc;
    posTrack -> Zc = Zc;
    posTrack -> r  = r;
    posTrack -> theta = ttheta;
    posTrack -> h  = hrot;
    posTrack -> Vz = Vz;
    posTrack -> pt = pt;
    posTrack -> Pz = Pz;
    posTrack -> Flag = 3;

    hrot = h_neg.h();
    ttheta = h_neg.phase();
    origin = h_neg.origin();
    Xc = h_neg.xcenter();// - pv2.x();
    Yc = h_neg.ycenter();// - pv2.y();
    Zc = origin.z();// - pv2.z();
    r  = 1./h_neg.curvature();
    Vz = r*tan(h_neg.dipAngle());
    pt = (h_neg.momentum(magnet/10) * 1.e-13).perp();
    Pz = (h_neg.momentum(magnet/10) * 1.e-13).z();


    if( la ==1){
      Xc = Xc +( pv1 - pv2).x() ;
      Yc = Yc +( pv1 - pv2).y() ;
      Zc = Zc +( pv1 - pv2).z() ;
    }

    negTrack -> Xc = Xc;
    negTrack -> Yc = Yc;
    negTrack -> Zc = Zc;
    negTrack -> r  = r;
    negTrack -> theta = ttheta;
    negTrack -> h  = hrot;
    negTrack -> Vz = Vz;
    negTrack -> pt = pt;
    negTrack -> Pz = Pz;
    negTrack -> Flag = 3;
	
    //copy Hui Long's 
    float d_root1,d_root2,d_p_pi,theta,r_p,r_pi,zij,rec_zij,rec_zij2;
    float fi_root1_p,fi_root2_p,fi_root1_pi,fi_root2_pi;
    float reci,reci2,recj,recj2;
    float v0_root1[3],v0_root2[3];
    float fi_p,fi_pi,ti,tj;
    int n1,n2,m1,m2,i,j;
    float alfa=0;
	
    //float R2_TPC = 30000.; //hard-coded cut // sprastar this will come from constants
	
    d_root2=500.;
    d_root1=500.;
	
    d_p_pi=TMath::Sqrt((negTrack->Xc-posTrack->Xc)*(negTrack->Xc-posTrack->Xc)+(negTrack->Yc-posTrack->Yc)*(negTrack->Yc-posTrack->Yc));
	
    theta=atan2((negTrack->Yc-posTrack->Yc),(negTrack->Xc-posTrack->Xc));
	
	
    r_p=posTrack->r;   //radium of proton helix
    r_pi=negTrack->r;    //radium of pion helix
	
    //if(d_p_pi > (r_p+r_pi)||d_p_pi<fabs(r_p-r_pi))return 500.;
	
	if(d_p_pi > (r_p+r_pi))return 500.; // sprastar: two helix not crossing each other removed ' = '
	if(r_p > r_pi && d_p_pi < r_p) return 500; // sprastar: One helix falls inside other 
	if(r_pi > r_p && d_p_pi < r_pi) return 500; // ""
	
    if(d_p_pi==(r_p+r_pi)){ //sprastar: >= to ==          //XZHU: will not come here, it will be rare that two tracks with d_p_pi == r_p+r_pi come from a single particle decay. that means their momenta should be in the same direction or in opposite direction.
		fi_root1_p=0.;
		fi_root2_p=0.;
		fi_root1_pi=PI;
		fi_root2_pi=PI;
    }
    else {
		alfa=acos((+r_pi*r_pi-d_p_pi*d_p_pi+r_p*r_p)/2./r_p/r_pi);
		
		fi_root1_p=acos((-r_pi*r_pi+d_p_pi*d_p_pi+r_p*r_p)/2./r_p/d_p_pi);
		fi_root2_p=-fi_root1_p;
		fi_root1_pi=alfa+fi_root1_p;
		fi_root2_pi=- fi_root1_pi;
		
    }
	
    v0_root1[2]=-1000.;
    v0_root2[2]=1000.;
	
    fi_root1_p=fi_root1_p+theta;
    fi_root2_p=fi_root2_p+theta;
    fi_root1_pi=fi_root1_pi+theta;
    fi_root2_pi=fi_root2_pi+theta;
	
    v0_root1[0]=posTrack->Xc+posTrack->r*cos(fi_root1_p);
    v0_root2[0]=posTrack->Xc+posTrack->r*cos(fi_root2_p);
    v0_root1[1]=posTrack->Yc+posTrack->r*sin(fi_root1_p);
    v0_root2[1]=posTrack->Yc+posTrack->r*sin(fi_root2_p);
    fi_root1_p=(fi_root1_p-posTrack->theta)/posTrack->h;
    fi_root2_p=(fi_root2_p-posTrack->theta)/posTrack->h;
    fi_root1_pi=(fi_root1_pi-negTrack->theta)/negTrack->h;
    fi_root2_pi=(fi_root2_pi-negTrack->theta)/negTrack->h;
    fi_root1_p=fi_root1_p-(int)(fi_root1_p/2/PI)*2*PI; //XZHU: try floor() here, instead of (int). since a wide range of period is tried below, the result should be the same. in fact, fi_root1_p should always be less than 0 (0,-2pi for 1 period), since we are looking backward for the position where the helix originates.
    fi_root2_p=fi_root2_p-(int)(fi_root2_p/2/PI)*2*PI;
    fi_root1_pi=fi_root1_pi-(int)(fi_root1_pi/2/PI)*2*PI;
    fi_root2_pi=fi_root2_pi-(int)(fi_root2_pi/2/PI)*2*PI;
	
    if((v0_root1[0]*v0_root1[0]+v0_root1[1]*v0_root1[1])<Constant::R2_TPC){ //XZHU: if the overlapping position is close (25 cm) to the outer layer of TPC. Even if the dca between two tracks is small, it is not interesting to us anymore. We find the right periods here, which is not necessary for particles with pt > 0.15GeV/c. in recHyperon, the period finding is not used.
		rec_zij=1000.;
		n1=-((int)(fabs((posTrack->Zc-pv1.z())/posTrack->Vz/2./PI))+1);
		n2=posTrack->Flag;  //XZHU: it is meanlingless to go forward in a helix to get the dca point where the particle (or helix) should orginates from! so i suggest to decrease this flag to 0, since fi_root1_pi is [-2pi,2pi]
		m1=-((int)(fabs((negTrack->Zc-pv2.z())/negTrack->Vz/2./PI))+1);;
		m2=negTrack->Flag;
		if(n1<-10)n1=-10;
		if(m1<-10)m1=-10;
		
		reci=fi_root1_p;
		recj=fi_root1_pi;
		for(i=n1;i<=n2;i++) {
			for(j=m1;j<=m2;j++) {
				
				ti=fi_root1_p+i*2*PI;
				tj=fi_root1_pi+j*2*PI;
				zij=-posTrack->Zc-posTrack->Vz*ti+negTrack->Zc+negTrack->Vz*tj;
				if(fabs(zij)<fabs(rec_zij)){
					rec_zij=zij;
					reci=ti;
					recj=tj;
				}
			}
		}
		
		fi_root1_p=reci;
		fi_root1_pi=recj;
		d_root1=rec_zij;
		
		dcaPToPi(&fi_root1_p,&fi_root1_pi,posTrack,negTrack,&d_root1,&v0_root1[0],alfa);
		
    }
	
    if((v0_root2[0]*v0_root2[0]+v0_root2[1]*v0_root2[1])<Constant::R2_TPC){
		rec_zij2=1000.;
		n1=-((int)(fabs((posTrack->Zc-pv1.z())/posTrack->Vz/2./PI))+1);
		n2=posTrack->Flag;
		m1=-((int)(fabs((negTrack->Zc-pv2.z())/negTrack->Vz/2./PI))+1);
		m2=negTrack->Flag;
		if(n1<-10)n1=-10;
		if(m1<-10)m1=-10;
		
		reci2=fi_root2_p;
		recj2=fi_root2_pi;
		for(i=n1;i<=n2;i++){
			for(j=m1;j<=m2;j++){
				ti=fi_root2_p+i*2*PI;
				tj=fi_root2_pi+j*2*PI;
				zij=-posTrack->Zc-posTrack->Vz*ti+negTrack->Zc+negTrack->Vz*tj;
				if(fabs(zij)<fabs(rec_zij2)){
					rec_zij2=zij;
					reci2=ti;
					recj2=tj;
				}
			}
		}
		fi_root2_p=reci2;
		fi_root2_pi=recj2;
		d_root2=rec_zij2;
		
		dcaPToPi(&fi_root2_p,&fi_root2_pi,posTrack,negTrack,&d_root2,&v0_root2[0],alfa);
    }
	
	
    //if(fabs(d_root2)>d_tolerance&&fabs(d_root1)>d_tolerance)return 1; //XZHU: apply the dca cut here
    if(fabs(d_root1)<fabs(d_root2)){
		fi_p=fi_root1_p;
		fi_pi=fi_root1_pi;
		*x=v0_root1[0];
		*(x+1)=v0_root1[1];
		*(x+2)=v0_root1[2];
		d=d_root1;
    }
    else{
		fi_p=fi_root2_p;
		fi_pi=fi_root2_pi;
		*x=v0_root2[0];
		*(x+1)=v0_root2[1];
		*(x+2)=v0_root2[2];
		d=d_root2;
    }
	
	
    *p1=posTrack->pt*cos(fi_p+posTrack->h*PI/2.);
    *p2=negTrack->pt*cos(fi_pi+negTrack->h*PI/2.);
	
    *(p1+1)=posTrack->pt*sin(fi_p+posTrack->h*PI/2.);
    *(p2+1)=negTrack->pt*sin(fi_pi+negTrack->h*PI/2.);
    *(p1+2)=posTrack->Pz;
    *(p2+2)=negTrack->Pz;
	
    //set the StThreeVectors
    xv0.setX(x[0]);
    xv0.setY(x[1]);
    xv0.setZ(x[2]);
	
    op1.setX(p1[0]);
    op1.setY(p1[1]);
    op1.setZ(p1[2]);
	
    op2.setX(p2[0]);
    op2.setY(p2[1]);
    op2.setZ(p2[2]);
	
    return d;
}



void  StStrangeV0::dcaPToPi(float *fi_p,float *fi_pi, StTrackHelix* posTrack,StTrackHelix* negTrack,float *d_root,float *v0,float alfa)
{
	
    //float R2_TPC = 30000.; //hard-coded cuts
    //float Z_TPC  = 180.  ;
    float t1,t2;
	
    float Co_a,Co_b,Co_c,Co_e,Co_f,dfi_p,dfi_pi;
    float new_d,x1,y1,z1,x2,y2,z2;
	
	
    Co_a=posTrack->r*posTrack->r+posTrack->Vz*posTrack->Vz;
    Co_b=posTrack->r*negTrack->r*cos(alfa)-posTrack->Vz*negTrack->Vz;  //XZHU: note the assumption here: helicities of two Helixes should be opposite
    Co_c=negTrack->r*negTrack->r+negTrack->Vz*negTrack->Vz;
    Co_e=(*d_root)*posTrack->Vz;
    Co_f=-(*d_root)*negTrack->Vz;
    dfi_p=0.;
    dfi_pi=0.;
    if(fabs(Co_a*Co_c-Co_b*Co_b)>0.){
		dfi_p=(Co_e*Co_c-Co_b*Co_f)/(Co_a*Co_c-Co_b*Co_b);
		
		dfi_pi=(Co_a*Co_f-Co_b*Co_e)/(Co_a*Co_c-Co_b*Co_b);
    }
	
    t1=*fi_p+dfi_p;
    t2=*fi_pi+dfi_pi;
	
    x1=posTrack->Xc+posTrack->r*cos(t1*posTrack->h+posTrack->theta);
    y1=posTrack->Yc+posTrack->r*sin(t1*posTrack->h+posTrack->theta);
    z1=posTrack->Zc+posTrack->Vz*t1;
	
    x2=negTrack->Xc+negTrack->r*cos(t2*negTrack->h+negTrack->theta);
    y2=negTrack->Yc+negTrack->r*sin(t2*negTrack->h+negTrack->theta);
    z2=negTrack->Zc+negTrack->Vz*t2;
    new_d=TMath::Sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    if(fabs(z1)>Constant::Z_TPC||fabs(z2)>Constant::Z_TPC)new_d=500.;
    if((x1*x1+y1*y1)>Constant::R2_TPC||(x2*x2+y2*y2)>Constant::R2_TPC)new_d=500.;
	
	if(new_d<fabs(*d_root)){
		*fi_p=t1*posTrack->h+posTrack->theta;
		*fi_pi=t2*negTrack->h+negTrack->theta;
		
		*d_root=new_d;
		
		*v0=(x2+x1)/2;
		*(v0+1)=(y2+y1)/2;
		*(v0+2)=(z2+z1)/2;
    }
    else{
		*fi_p=(t1-dfi_p)*posTrack->h+posTrack->theta;
		*fi_pi=(t2-dfi_pi)*negTrack->h+negTrack->theta;
		
		*d_root = fabs(*d_root); //XZHU: remove the negative dca bug
		
		*(v0+2)=((negTrack->Vz*(t2-dfi_pi)+negTrack->Zc)+(posTrack->Vz*(t1-dfi_p)+posTrack->Zc))/2.;
    }
    return;
}


StThreeVectorF StStrangeV0::momentumTrack(const Int_t i) const
{
	StThreeVectorF ret(0.,0.,0.);
	if(i==pos||i==neg) return mMomentum[i];
	else return ret;
}

void StStrangeV0::setParticleHypothesis(const Int_t ip_pos, const Int_t ip_neg)
{
	if(ip_pos<0 || ip_neg<0 ) return;
	StLorentzVectorF fourMom_pos(mMomentum[pos], (Float_t)TMath::Sqrt(mMomentum[pos].mag2()+Constant::mMass[ip_pos]*Constant::mMass[ip_pos]));
	StLorentzVectorF fourMom_neg(mMomentum[neg], (Float_t)TMath::Sqrt(mMomentum[neg].mag2()+Constant::mMass[ip_neg]*Constant::mMass[ip_neg]));
	mM = (fourMom_pos + fourMom_neg).m();

/* both same
	double energy1 = (Float_t)TMath::Sqrt(mMomentum[pos].mag2()+Constant::mMass[ip_pos]*Constant::mMass[ip_pos]);
	double energy2 = (Float_t)TMath::Sqrt(mMomentum[neg].mag2()+Constant::mMass[ip_neg]*Constant::mMass[ip_neg]);
	mM = (Float_t)TMath::Sqrt(Constant::mMass[ip_pos]*Constant::mMass[ip_pos] + Constant::mMass[ip_neg]*Constant::mMass[ip_neg] + 2.*energy1*energy2 -2.* mMomentum[pos].dot(mMomentum[neg])  );
*/

}


Double_t StStrangeV0::energy(const Int_t ip_pos, const Int_t ip_neg)
{
	double energy1 = (Float_t)TMath::Sqrt(mMomentum[pos].mag2()+Constant::mMass[ip_pos]*Constant::mMass[ip_pos]);
	double energy2 = (Float_t)TMath::Sqrt(mMomentum[neg].mag2()+Constant::mMass[ip_neg]*Constant::mMass[ip_neg]);
	return energy1 + energy2;
}





StThreeVectorF StStrangeV0::momentum() const
{
	return mMomentum[pos] + mMomentum[neg];
}



