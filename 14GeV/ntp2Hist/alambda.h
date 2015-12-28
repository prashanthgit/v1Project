
#ifndef alambda_h
#define alambda_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


class alambda {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UShort_t        mFlag;
   UShort_t        mCentralityBin;
   UShort_t        mRapidityBin;
   UShort_t        mPhiPsiBin;
   UShort_t        mPtBin;
   Double_t        mLBarMass;

   // List of branches
   TBranch        *b_mFlag;   //!
   TBranch        *b_mCentralityBin;   //!
   TBranch        *b_mRapidityBin;   //!
   TBranch        *b_mPhiPsiBin;   //!
   TBranch        *b_mPtBin;   //!
   TBranch        *b_mLBarMass;   //!

   alambda(TTree *tree=0);
   virtual ~alambda();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
};

#endif

#ifdef alambda_cxx
alambda::alambda(TTree *tree) : fChain(0) {
  TChain * chain = new TChain("v0","");
  Int_t  nfiles;
  Char_t *file;
  Char_t path_file[250];
  Char_t dirname[]="/star/data05/scratch/keane/14GeV/flag/root/";
  void *dirhandle = gSystem->OpenDirectory(gSystem->ExpandPathName(dirname));

  while (file = gSystem->GetDirEntry(dirhandle)) {
    if (strstr(file,"pV1A_FDC4") != 0) {
      strcpy(path_file,dirname);
      strcat(path_file,file);
      TFile fileRoot(path_file);
      if(fileRoot.IsZombie()) continue;
      chain->Add(path_file);
      cout << nfiles<<" :: file = " << path_file << endl;
      nfiles++;
    }
  }

  tree = chain;
  Init(tree);

}

alambda::~alambda()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t alambda::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t alambda::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

void alambda::Init(TTree *tree) {
 
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("mFlag", &mFlag, &b_mFlag);
   fChain->SetBranchAddress("mCentralityBin", &mCentralityBin, &b_mCentralityBin);
   fChain->SetBranchAddress("mRapidityBin", &mRapidityBin, &b_mRapidityBin);
   fChain->SetBranchAddress("mPhiPsiBin", &mPhiPsiBin, &b_mPhiPsiBin);
   fChain->SetBranchAddress("mPtBin", &mPtBin, &b_mPtBin);
   fChain->SetBranchAddress("mLBarMass", &mLBarMass, &b_mLBarMass);
}

#endif // #ifdef alambda_cxx
