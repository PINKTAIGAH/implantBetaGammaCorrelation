#define aida_frs_efficiency_cxx
// The class definition in evt.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("evt.C")
// root> T->Process("evt.C","some options")
// root> T->Process("evt.C+")
//

#include "aida_frs_efficinecy.h"
#include <TH2.h>
#include <TStyle.h>


void evt::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void evt::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t evt::Process(Long64_t entry){
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetLocalEntry(entry);

   int aidaImplantHits = aidaImplantHits_DSSD.GetSize();
    
   // Define set to hold index of altready filled aida events
   std::unordered_set<int> filled_gatedimplanttree{};

   // Loop over frs events
   for (auto const& FrsMultiHitItem : FrsMultiHitItem){
     // Define frs data items
     std::vector<float> const& AoQ = FrsMultiHitItem.Get_ID_AoQ_corr_s2s4_mhtdc();
     std::vector<float> const& Z = FrsMultiHitItem.Get_ID_z41_mhtdc();
     std::vector<float> const& Z2 = FrsMultiHitItem.Get_ID_z42_mhtdc();


     // loop over frs subevents
     for (int i=0; i<AoQ.size(); i++){

       // Check if indexes are out of bounds & skip
       if ( i>=AoQ.size() || i>=Z.size() ){
         continue;
       }

       if ( zAoQ_cut->IsInside(AoQ[i], Z[i]) && z1z2_cut->IsInside(Z[i], Z2[i]) ){
         totalGatedHits++;
         // Check if we have aida hits
         if ( aidaImplantHits>0 ){
           for ( int j=0; j<aidaImplantHits; j++ ){
             if ( implant_dssd[j]==1 && filled_gatedimplanttree.count(j)==0 ){
               totalImplantedGatedHits++;
               filled_gatedimplanttree.insert(j);
             }
           }
         }
       }
     }
   }
   return kTRUE;
}

void evt::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void evt::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   std::cout << "# of counts: " << totalGatedHits << ""
   std::cout << "# of gated counts: " << totalImplantedGatedHits << std::endl;

}
