#ifndef kshort_h
#define kshort_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


class kshort {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UShort_t        mFlag;
   UShort_t        mCentralityBin;
   UShort_t        mRapidityBin;
   UShort_t        mPhiPsiBin;
   UShort_t        mPtBin;
   Double_t        mKshortMass;

   // List of branches
   TBranch        *b_mFlag;   //!
   TBranch        *b_mCentralityBin;   //!
   TBranch        *b_mRapidityBin;   //!
   TBranch        *b_mPhiPsiBin;   //!
   TBranch        *b_mPtBin;   //!
   TBranch        *b_mKshortMass;   //!

   kshort(TTree *tree=0);
   virtual ~kshort();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
};

#endif

#ifdef kshort_cxx
kshort::kshort(TTree *tree) : fChain(0) {

  TChain * chain = new TChain("v0","");
  Int_t  nfiles;
  Char_t *file;
  Char_t path_file[250];
  Char_t dirname[]="/star/u/sprastar/data05/19GeV/flag/root/";
  void *dirhandle = gSystem->OpenDirectory(gSystem->ExpandPathName(dirname));

  while (file = gSystem->GetDirEntry(dirhandle)) {
    if (strstr(file,"pV1K_6F84") != 0) {
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

kshort::~kshort() {
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t kshort::GetEntry(Long64_t entry) {
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t kshort::LoadTree(Long64_t entry) {
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}

void kshort::Init(TTree *tree) {
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
   fChain->SetBranchAddress("mKshortMass", &mKshortMass, &b_mKshortMass);
}

#endif // #ifdef kshort_cxx
