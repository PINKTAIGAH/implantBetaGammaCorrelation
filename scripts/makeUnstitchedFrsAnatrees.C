#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>

#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "../../../c4Root/c4data/frsData/FrsHitData.h"

// *************************************************************************************
// ****************************** DEFINE SCRIPT CONSTANTS ****************************** 
// *************************************************************************************

// *************************************************************************************
// ****************************** DEFINE CUSTOM STRUCTURES ***************************** 
// *************************************************************************************

struct frs_data_struct{
  long long time;
  double aoq;
  double x2;
  double x4;
  double z41;
  double z42;
  double sci41E;
  double sci42E;
}frs_data;


// *************************************************************************************
// ****************************** START MACRO ******************************************
// *************************************************************************************


void makeUnstitchedFrsAnatrees(const char* input, const char* output) {

  // *************************************************************************************
  // ****************************** OPEN & CREATE FILES **********************************
  // *************************************************************************************
  
  // Open the file
  TFile* file = TFile::Open(input);
  if (!file) {
    std::cerr << "Error: Could not open file " << std::endl;
    return;
  }

  // Get the existing tree
  TTree* old_tree = dynamic_cast<TTree*>(file->Get("evt"));
  if (!old_tree) {
    std::cerr << "Error: Could not find tree evt in file " << std::endl;
    return;
  }


  // *************************************************************************************
  // ****************************** DEFINE TREEREADER ************************************
  // *************************************************************************************

  // Create a tree reader
  TTreeReader reader(old_tree);


  TTreeReaderArray<FrsMultiHitItem> FrsMultiItem(reader, "FrsMultiHitData");
  TTreeReaderArray<FrsHitItem> FrsItem(reader, "FrsHitData");

  // *************************************************************************************
  // ****************************** DEFINE OUTPUT TREES **********************************
  // *************************************************************************************

  // Open the output file
  TFile* outputFile = new TFile(output, "RECREATE");
  if (!outputFile) {
      std::cerr << "Error: Could not create output file " << std::endl;
      return;
  }
  
  TTree* frs_tree = new TTree("frs_tree", "New AIDA Analysis Tree");  
  frs_tree->SetDirectory(0);
  frs_tree->Branch("frs", &frs_data, "time/l:aoq/D:x2/D:x4/D:z41/D:z42/D:sci41E/D:sci42E/D");


  // *************************************************************************************
  // ****************************** LOOP OVER ALL EVENTS *********************************
  // *************************************************************************************

  const char spinner[] = {'-', '\\', '|', '/'};
  int totalEntries = reader.GetEntries(true);

  // Loop over all entries in the old tree
  for (int current_entry = 0; current_entry<totalEntries; ++current_entry) {

    reader.SetEntry(current_entry);


    // *************************************************************************************
    // ****************************** LOOP OVER FRS EVENTS *******************************
    // *************************************************************************************


    for (int mbsIdx=0; mbsIdx<FrsMultiItem.GetSize(); mbsIdx++){

      std::set<int> frs_filledtree = {};

      // Get corresponding FrsHitItem and FrsMultiHitItem
      auto& frs_multi_item = FrsMultiItem[mbsIdx];
      auto& frs_item = FrsItem[mbsIdx];

      // Get values inside the item
      std::vector<float> const& AoQ = frs_multi_item.Get_ID_AoQ_corr_s2s4_mhtdc();
      std::vector<float> const& Z = frs_multi_item.Get_ID_z41_mhtdc();
      std::vector<float> const& Z2 = frs_multi_item.Get_ID_z42_mhtdc();
      std::vector<float> const& X2 = frs_multi_item.Get_ID_s2x_s2s4_mhtdc();
      std::vector<float> const& X4 = frs_multi_item.Get_ID_s4x_mhtdc();
      Long64_t const& frs_wr = frs_item.Get_wr_t();
      float const& sci41E = frs_item.Get_sci_e_41();
      float const& sci42E = frs_item.Get_sci_e_42();

      for (int mhtdcIdx=0; mhtdcIdx<AoQ.size(); mhtdcIdx++){
        
        if ( mhtdcIdx >= Z.size() || mhtdcIdx >= Z2.size() || mhtdcIdx >= X4.size() || mhtdcIdx >= X4.size()) { continue; }

        frs_data.time = frs_wr;
        frs_data.sci41E = sci41E;
        frs_data.sci42E = sci42E;
        frs_data.aoq = AoQ[mhtdcIdx];
        frs_data.x2 = X2[mhtdcIdx];
        frs_data.x4 = X4[mhtdcIdx];
        frs_data.z41 = Z[mhtdcIdx];
        frs_data.z42 = Z2[mhtdcIdx];

        frs_filledtree.insert(mhtdcIdx);

        frs_tree->Fill();
      
      }

    }

    // *************************************************************************************
    // ****************************** PRINT OUT LOOP PROGRESS ******************************
    // *************************************************************************************
    
    // Show the progress of the loop
    if (reader.GetCurrentEntry() % 100000 == 0) {
      int progress = (reader.GetCurrentEntry() * 100) / totalEntries;
      char spin = spinner[reader.GetCurrentEntry() / 1000 % 4];
      std::cout << "\rProcessing the tree " << reader.GetCurrentEntry() << " (" << progress << "%) " << spin << std::flush;
    }
  
  }

  std::cout << std::endl;



  // *************************************************************************************
  // ****************************** WRITE TREES ******************************************
  // *************************************************************************************

  std::cout << "Writing trees ..." << std::endl;
  outputFile->cd();
  
  frs_tree->Write();

  std::cout << "Finished writing trees successfully!" << std::endl << std::endl;

  // *************************************************************************************
  // ****************************** CLEANUP **********************************************
  // *************************************************************************************

  // Close the input file
  file->Close();

  // CLose and write the output file
  outputFile->Close();
  
  // Delete pointer objects
  delete file;
  delete outputFile;
}
