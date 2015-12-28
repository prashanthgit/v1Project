#ifndef Constants_h
#define Constants_h

// name of histograme
const static std:string histname[9] = {"pV1P","pV1AP","pV1PiP","pV1PiN","pV1KP","pV1KN","pV1L","pV1A","pV1K"};
const string energyname = "19GeV";

Float_t ptBins[] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };
//                     0    1    2    3    4    5    6    7    8    9   10   11    12   13  14    15
Float_t pT_min[16] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5};
Float_t pT_max[16] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0};


// Vz
const static double vzDef = 50;
const static double vzMin = 40;
const static double vzMax = 60;




// Momentum & pT cuts
// Bar - proton mom equvalanec
// Mes - pi,K mom equvalance
// Protons
const static double protonPMaxDef  = 2.80; 
const static double protonPtMaxDef = 2.00;  
const static double protonPtMinDef = 0.40;  


const static double protonPMaxMin  = 2.24; 
const static double protonPtMaxMin = 1.6;  
const static double protonPtMinMin = 0.48;  

const static double protonPMaxMax  = 3.36; 
const static double protonPtMaxMax = 2.4;  
const static double protonPtMinMax = 0.32;  

// pions
const static double pionsPMaxDef  = 1.6; 
const static double pionsPtMinDef = 0.2;  

const static double pionsPMaxMax  = 1.92; 
const static double pionsPtMinMax = 0.2;  

const static double pionsPMaxMin  = 1.28; 
const static double pionsPtMinMin = 0.4;  

// Kaons
const static double kaonsPtMinDef = 0.2;  
const static double kaonsPMaxDef  = 1.60; 

const static double kaonsPtMinMax = 0.2;  
const static double kaonsPMaxMax  = 1.92; 

const static double kaonsPMaxMin  = 1.28; 
const static double kaonsPtMinMin = 0.4;  

// Lambda & anti Lambda
const static double lambdaPtMinDef  = 0.2;
const static double lambdaPtMaxDef  = 5.0;

const static double lambdaMomMaxBar = 2.8;  
const static double lambdaPtMinBar  = 0.4;
const static double lambdaPtMaxBar  = 2.0; 

const static double lambdaMomMaxMes = 1.6;  
const static double lambdaPtMinMes  = 0.2;


// Kshort
const static double kshortPtMinDef = 0.2; 
const static double kshortPtMaxDef = 5.0;  

const static double kshortMomMaxBar = 2.8;  
const static double kshortPtMinBar  = 0.4;
const static double kshortPtMaxBar  = 2.0; 

const static double kshortMomMaxMes = 1.6; 
const static double kshortPtMinMes  = 0.2;  









// mass cuts 
// for protons
const static double mass2protonMinDef =  0.8; 
const static double mass2protonMaxDef =  1.0; 

const static double mass2protonMinMin =  0.88; 
const static double mass2protonMaxMin =  0.9; 

const static double mass2protonMinMax =  0.72; 
const static double mass2protonMaxMax =  1.1; 

// for pions
const static double mass2pionMinDef =  -0.1; 
const static double mass2pionMaxDef =   0.1; 

const static double mass2pionMinMin =  -0.05; 
const static double mass2pionMaxMin =   0.05; 

const static double mass2pionMinMax =  -0.15; 
const static double mass2pionMaxMax =   0.15; 


// for kaons
const static double mass2kaonsMinDef =  0.20;
const static double mass2kaonsMaxDef =  0.35; 

const static double mass2kaonsMinMin =  0.20; 
const static double mass2kaonsMaxMin =  0.30; 

const static double mass2kaonsMinMax =  0.18; 
const static double mass2kaonsMaxMax =  0.40; 

// Lambda
const static double mass2LprotonMinDef =  0.5; 
const static double mass2LprotonMaxDef =  1.5; 

const static double mass2LprotonMinMin =  0.8; 
const static double mass2LprotonMaxMin =  1.0; 

const static double mass2LpionMinMin = -0.1;
const static double mass2LpionMaxMin =  0.1;

// kShort
// no cuts for K0s










// nSigma
// for v1Tracks (pi, K, p)
const static double nsigmaDef = 2.0; 
const static double nsigmaMin = 1.6; 
const static double nsigmaMax = 2.4; 

// V0 
const static double nsigmav0Def = 3.0; 
const static double nsigmav0Min = 2.0; 
const static double nsigmav0Max = 3.2; 









// all tracks for v1tracks and V0
const static double nhitsDef = 15;  
const static double nhitsMin = 12; 
const static double nhitsMax = 18; 

const static double nhitsposDef = 0.52; 
const static double nhitsposMin = 0.00; 
const static double nhitsposMax = 0.60; 

// for pi,K,p
const static double dcaDef = 3.0; 
const static double dcaMin = 2.4; 
const static double dcaMax = 3.6; 











// V0 cuts for Lambda & anit Lambda
const static double lambdaDcaV0PvDef_0 = 0.75; 
const static double lambdaDcaV0PvMin_0 = 0.60; 
const static double lambdaDcaV0PvMax_0 = 0.90; 

const static double lambdaDecaylengthDef_0 = 4.0; 
const static double lambdaDecaylengthMin_0 = 3.2; 
const static double lambdaDecaylengthMax_0 = 4.8; 

const static double lambdaDcaPtoPvDef_0 = 0.60; 
const static double lambdaDcaPtoPvMin_0 = 0.48; 
const static double lambdaDcaPtoPvMax_0 = 0.72; 

const static double lambdaDcaPitoPvDef_0 = 1.7;
const static double lambdaDcaPitoPvMin_0 = 1.36;
const static double lambdaDcaPitoPvMax_0 = 2.24;



const static double lambdaDcaV0PvDef_1 = 0.75; 
const static double lambdaDcaV0PvMin_1 = 0.60; 
const static double lambdaDcaV0PvMax_1 = 0.90; 

const static double lambdaDecaylengthDef_1 = 3.5; 
const static double lambdaDecaylengthMin_1 = 2.8; 
const static double lambdaDecaylengthMax_1 = 4.2; 

const static double lambdaDcaPtoPvDef_1 = 0.50; 
const static double lambdaDcaPtoPvMin_1 = 0.40; 
const static double lambdaDcaPtoPvMax_1 = 0.60; 

const static double lambdaDcaPitoPvDef_1 = 1.5;
const static double lambdaDcaPitoPvMin_1 = 1.2;
const static double lambdaDcaPitoPvMax_1 = 1.8;




const static double lambdaDcaV0PvDef_2 = 1.2; 
const static double lambdaDcaV0PvMin_2 = 0.96; 
const static double lambdaDcaV0PvMax_2 = 1.44; 

const static double lambdaDecaylengthDef_2 = 2.5; 
const static double lambdaDecaylengthMin_2 = 2.0; 
const static double lambdaDecaylengthMax_2 = 3.0; 

const static double lambdaDcaPtoPvDef_2 = 0.15; 
const static double lambdaDcaPtoPvMin_2 = 0.12; 
const static double lambdaDcaPtoPvMax_2 = 0.18; 

const static double lambdaDcaPitoPvDef_2 = 0.8;
const static double lambdaDcaPitoPvMin_2 = 0.64;
const static double lambdaDcaPitoPvMax_2 = 0.96;




const static double lambdaDcaV0PvDef_3 = 1.3; 
const static double lambdaDcaV0PvMin_3 = 1.04; 
const static double lambdaDcaV0PvMax_3 = 1.56; 

const static double lambdaDecaylengthDef_3 = 2.0; 
const static double lambdaDecaylengthMin_3 = 1.6; 
const static double lambdaDecaylengthMax_3 = 2.4; 

const static double lambdaDcaPtoPvDef_3 = 0.1; 
const static double lambdaDcaPtoPvMin_3 = 0.8; 
const static double lambdaDcaPtoPvMax_3 = 0.12; 

const static double lambdaDcaPitoPvDef_3 = 0.7;
const static double lambdaDcaPitoPvMin_3 = 0.56;
const static double lambdaDcaPitoPvMax_3 = 0.84;

const static double lambdaDcaDaughtersDef = 1.0; 
const static double lambdaDcaDaughtersMin = 0.8;
const static double lambdaDcaDaughtersMax = 1.2;




// V0 cuts for Kshort

const static double kshortDcaV0PvDef = 0.80; 
const static double kshortDcaV0PvMin = 0.64;
const static double kshortDcaV0PvMax = 0.96; 

const static double kashortDecaylengthDef = 3.0; 
const static double kashortDecaylengthMin = 2.4;
const static double kashortDecaylengthMax = 3.6;

const static double kshortDcaPitoPvDef = 0.7;
const static double kshortDcaPitoPvMin = 0.56;
const static double kshortDcaPitoPvMax = 0.84;

const static double kshortDcaDaughtersDef = 1.0;
const static double kshortDcaDaughtersMin = 0.8;
const static double kshortDcaDaughtersMax = 1.2;



#endif

