//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar  6 16:46:37 2025 by ROOT version 6.26/10
// from TTree evt//cbmout_0
// found on file: 107Ag_implantingAIDA_0001.root
//////////////////////////////////////////////////////////

#ifndef aida_frs_efficiency_h
#define aida_frs_efficiency_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include "TClonesArray.h"

#include <vector>



class aida_frs_efficiency : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain


   // FRS Gated
   TFile* gate = new TFile("/lustre/gamma/gbrunic/G302/analysis/ionbeta/gates/82Nb.root");
   if (!gate){
     std::cerr << "Error: Could not find FRSCUT in file " << std::endl;
     std::exit(1);
   }
   TCutG* zAoQ_cut = (TCutG*)gate->Get("cut_Z_AoQ");
   TCutG* z1z2_cut = (TCutG*)gate->Get("cut_Z_Z2");

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<Int_t> AidaImplantHits_DSSD = {fReader, "AidaImplantHits.DSSD"};
   TTreeReaderArray<Bool_t> AidaImplantHits_Stopped = {fReader, "AidaImplantHits.Stopped"};
   TTreeReaderArray<FrsMultiHitItem> FrsMultiHitItem(fReader, "FrsMultiHitData");

   int totalGatedHits = 0;
   int totalImplantedGatedHits = 0;
   int totalStoppedGatedHits = 0;

   aida_frs_efficiency(TTree * /*tree*/ =0) { }
   virtual ~aida_frs_efficiency() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(aida_frs_efficiency,0);

};

#endif

#ifdef aida_frs_efficiency_cxx
void aida_frs_efficiency::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t aida_frs_efficiency::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef evt_cxx
