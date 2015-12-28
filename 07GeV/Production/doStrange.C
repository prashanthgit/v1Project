#if !defined(__CINT__) || defined (__MAKECINT__)
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#endif

class StChain;
class StMuDstMaker;
class StStrangeMaker;

void doStrange( Int_t nEvents, Int_t nFiles, TString InputFileList, TString OutputDir, TString JobIdName, Float_t energy){

  gROOT   -> Reset();
  gROOT   -> Macro("loadMuDst.C");
  gSystem -> Load("StStrangeMaker.so") ;
  gSystem -> Load("StBBCEventPlane.so") ;
  gSystem -> Load("StRefMultCorr.so") ;

  StChain* chain = new StChain;
  StMuDstMaker* muDstMaker  =  new StMuDstMaker(0,0,"",InputFileList,"MuDst",nFiles) ;

  muDstMaker -> SetStatus("*",0) ;               // Turn off all branches
  muDstMaker -> SetStatus("MuEvent",1) ;         // Turn on the Event data (esp. Event number)
  muDstMaker -> SetStatus("PrimaryTracks",1) ;   // Turn on the primary track data
  muDstMaker -> SetStatus("GlobalTracks",1) ;
  muDstMaker -> SetStatus("BTofHit",1) ;
  muDstMaker -> SetDebug(0) ;                    // Turn off Debug information

  StStrangeMaker* AnalysisCode = new StStrangeMaker(muDstMaker);
  
  AnalysisCode -> SetJobIdName(JobIdName);
  AnalysisCode -> SetEnergy(energy);
  AnalysisCode -> SetDebug(1);

  chain -> Init() ;
  for(Int_t evNo=0; evNo<nEvents; evNo++){
    chain->Clear();
    Int_t istat = chain->Make();
    if ( istat ) break;
  }
  chain -> Finish() ;
  delete chain ;

}
