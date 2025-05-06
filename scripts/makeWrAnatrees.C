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

namespace constants{

  const int DSSSD = 1;

  const int TOTAL_GERMANIUM_DETECTORS = 15; // Number of germanium detectors 
  const int TOTAL_GERMANIUM_CRYSTALS = 2; // Number of germanium crystals

}

// *************************************************************************************
// ****************************** DEFINE CUSTOM STRUCTURES ***************************** 
// *************************************************************************************

// Define the data structures for branches
struct implant_data{
  Long64_t time_wr;
  int placeholder;
}aida_implant_data;

struct decay_data{
  Long64_t time_wr;
  Long64_t time_x;
  Long64_t time_y;
  double energy;
  double energy_x;
  double energy_y;
}aida_decay_data;

struct frs_data{
  Long64_t time_wr;
  int placeholder;
}frs_data;

struct frs_mhtdc_data{
  Long64_t time_wr;
  int empty;
}frs_mhtdc_data;
  
struct bplast_data{
  ULong_t time_wr;
  ULong_t time_abs_evt;
}bplast_data;

struct germanium_data{
  ULong_t time_wr;
  Long_t time_abs_evt;
}germanium_data;
  

// *************************************************************************************
// ****************************** DEFINE FUNCTIONS *************************************
// *************************************************************************************


// *************************************************************************************
// ****************************** START MACRO ******************************************
// *************************************************************************************


void makeWrAnatrees(const char* input, const char* output) {

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

  // Create tree reader values for the branches you want to read
  TTreeReaderArray<Long64_t> implant_time(reader, "AidaImplantHits.Time");
  TTreeReaderArray<Int_t> implant_dssd(reader, "AidaImplantHits.DSSD");

  TTreeReaderArray<Long64_t> decay_time(reader, "AidaDecayHits.Time");
  TTreeReaderArray<Long64_t> decay_time_x(reader, "AidaDecayHits.TimeX");
  TTreeReaderArray<Long64_t> decay_time_y(reader, "AidaDecayHits.TimeY");
  TTreeReaderArray<Double_t> decay_energy(reader, "AidaDecayHits.Energy");
  TTreeReaderArray<Double_t> decay_energy_x(reader, "AidaDecayHits.EnergyX");
  TTreeReaderArray<Double_t> decay_energy_y(reader, "AidaDecayHits.EnergyY");
  TTreeReaderArray<Int_t> decay_dssd(reader, "AidaDecayHits.DSSD");

  TTreeReaderArray<ULong_t> bplast_time(reader, "bPlastTwinpeaksCalData.fwr_t");
  TTreeReaderArray<ULong_t> bplast_abs_evt_time(reader, "bPlastTwinpeaksCalData.fabsolute_event_time");
  TTreeReaderArray<UShort_t> bplast_id(reader, "bPlastTwinpeaksCalData.fdetector_id");

  TTreeReaderArray<ULong_t> germanium_time(reader, "GermaniumCalData.fwr_t");
  TTreeReaderArray<Long_t> germanium_abs_evt_time(reader, "GermaniumCalData.fabsolute_event_time");
  TTreeReaderArray<UInt_t> germanium_det(reader, "GermaniumCalData.fdetector_id");
  TTreeReaderArray<UInt_t> germanium_cry(reader, "GermaniumCalData.fcrystal_id");
  TTreeReaderArray<Double_t> germanium_energy(reader, "GermaniumCalData.fchannel_energy");

  TTreeReaderArray<FrsMultiHitItem> FrsMultiItem(reader, "FrsMultiHitData");
  TTreeReaderArray<FrsHitItem> FrsItem(reader, "FrsHitData");

  // *************************************************************************************
  // ****************************** DEFINE OUTPUT FILE ***********************************
  // *************************************************************************************

  // Open the output file
  TFile* outputFile = new TFile(output, "RECREATE");
  if (!outputFile) {
      std::cerr << "Error: Could not create output file " << std::endl;
      return;
  }

  // *************************************************************************************
  // ****************************** DEFINE OUTPUT TREES **********************************
  // *************************************************************************************
  
  // Create implant tree and branches for anatree
  TTree* implant_tree = new TTree("aida_implant_tree", "New AIDA Analysis Tree");
  implant_tree->SetDirectory(0); // Save tree object in memory untill explicitly written to file
  implant_tree->Branch("implant", &aida_implant_data, "time_wr/l:placeholder/I");

  // Create decay tree and branches for anatree
  TTree* decay_tree = new TTree("aida_decay_tree", "New AIDA Analysis Tree");
  decay_tree->SetDirectory(0); // Save tree object im memory untill explicitly written to file
  decay_tree->Branch("decay", &aida_decay_data, "time_wr/l:time_x/l:time_y/l:e/D:ex/D:ey/D");

  // Create frs tree and branches for anatree
  TTree* frs_tree = new TTree("frs_tree", "New AIDA Analysis Tree");
  frs_tree->SetDirectory(0); // Save tree object im memory untill explicitly written to file
  frs_tree->Branch("frs", &frs_data, "time_wr/l:placeholder/I");

  // Create frs mhtdc tree and branches for anatree
  TTree* frs_mhtdc_tree = new TTree("frs_mhtdc_tree", "New AIDA Analysis Tree");
  frs_mhtdc_tree->SetDirectory(0); // Save tree object im memory untill explicitly written to file
  frs_mhtdc_tree->Branch("frs_mhtdc", &frs_mhtdc_data, "time_wr/l:placeholder/I");

  // Create germanium tree and branches for anatree
  TTree* germanium_tree = new TTree("germanium_tree", "New AIDA Analysis Tree");
  germanium_tree->SetDirectory(0); // Save tree object im memory untill explicitly written to file
  germanium_tree->Branch("germanium", &germanium_data, "time_wr/l:time_abs_evt/l");

  // Create bplast tree and branches for anatree
  TTree* bplast_tree = new TTree("bplast_tree", "New AIDA Analysis Tree");
  bplast_tree->SetDirectory(0); // Save tree object im memory untill explicitly written to file
  bplast_tree->Branch("bplast", &bplast_data, "time_wr/l:time_abs_evt/l");

  // *************************************************************************************
  // ****************************** LOOP OVER ALL EVENTS *********************************
  // *************************************************************************************

  const char spinner[] = {'-', '\\', '|', '/'};
  int totalEntries = reader.GetEntries(true);

  // Loop over all entries in the old tree
  while (reader.Next()) {

    // *************************************************************************************
    // ****************************** GET EVENT INFORMATION *********************************
    // *************************************************************************************

    // sizes of events
    int germaniumhits = germanium_time.GetSize();
    int aidadecayhits = decay_time.GetSize();
    int aidaimphits = implant_time.GetSize();
    int bplasthits = bplast_time.GetSize();

    // *************************************************************************************
    // ****************************** LOOP OVER AIDA IMPLANTS  *****************************
    // *************************************************************************************

    // Empty set for each subevents
    std::set<int> aidaimplant_filledtree{};

    if( aidaimphits > 0 ){
      for ( int j = 0; j < aidaimphits; j++ ){

        if ( implant_dssd[j] != constants::DSSSD ){ continue; }
        if ( aidaimplant_filledtree.count(j) != 0 ){ continue; }

        // Fill up implant datastruct
        aida_implant_data.time_wr = implant_time[j];
        aida_implant_data.placeholder = 0;
        implant_tree->Fill();

        aidaimplant_filledtree.insert(j);
      }
    }

    // *************************************************************************************
    // ****************************** LOOP OVER DECAY EVENTS *******************************
    // *************************************************************************************

    // Empty set for each subevents
    std::set<int> aidadecay_filledtree{};

    // Loop over decay subevents
    if( aidadecayhits == 1 && aidaimphits == 0 ){
      for ( int j = 0; j < aidadecayhits; j++ ) {

        if ( decay_dssd[j] != constants::DSSSD ){ continue; }
        if ( TMath::Abs(decay_time[aidadecayhits -1] - decay_time[0]) > 7000 ) { continue; }
        if ( aidadecay_filledtree.count(j) != 0 ){ continue; }

        aida_decay_data.time_wr = decay_time[j];
        aida_decay_data.time_x = decay_time_x[j];
        aida_decay_data.time_y = decay_time_y[j];
        aida_decay_data.energy = decay_energy[j];
        aida_decay_data.energy_x = decay_energy_x[j];
        aida_decay_data.energy_y = decay_energy_y[j];
        decay_tree->Fill();

        aidadecay_filledtree.insert(j);

      } // End of decay loop

    }

    // *************************************************************************************
    // ****************************** LOOP OVER BPLAST  ************************************
    // *************************************************************************************
      
    if ( bplasthits > 0 ){
      for(int j=0; j<bplasthits; j++){

        if ( bplast_id[j] > 64 ) { continue; }

        bplast_data.time_wr = bplast_time[j];
        bplast_data.time_abs_evt = bplast_abs_evt_time[j];
        bplast_tree->Fill();

      }
    }

    // *************************************************************************************
    // ****************************** LOOP OVER GERMANIUM EVENTS *******************************
    // *************************************************************************************

    
    std::set<int> germanium_filledtree{};

    // Loop over decay subevents
    if ( germaniumhits > 0 ){

      for ( int j = 0; j < germaniumhits; j++ ){

        if ( germanium_det[j] > constants::TOTAL_GERMANIUM_DETECTORS ){ continue; }
        if ( germanium_cry[j] > constants::TOTAL_GERMANIUM_CRYSTALS ){ continue; } 
        if ( germanium_energy[j] <= 25. ){ continue; } 
        if ( germanium_filledtree.count(j) != 0 ){ continue; }

        germanium_data.time_wr = germanium_time[j];
        germanium_data.time_abs_evt = germanium_abs_evt_time[j];
        germanium_tree->Fill();

        germanium_filledtree.insert(j);

      } // End of germanium loop

    }


    // *************************************************************************************
    // ****************************** LOOP OVER FRS EVENTS *******************************
    // *************************************************************************************

    for (auto const& frs_item : FrsItem) {

      Long64_t time_wr = frs_item.Get_wr_t();
      Float_t z41 = frs_item.Get_ID_z41();

      if ( z41 == 0.0) { continue; }

      frs_data.time_wr = time_wr;
      frs_data.placeholder = 0;
      frs_tree->Fill();

      // std::vector<Long64_t> const& Time_WR = frs_item.Get_wr_t();

      // for ( const auto time_wr : Time_WR){

      //   frs_data.time_wr = time_wr;
      //   frs_tree->Fill();

      // }
      
    }


    // *************************************************************************************
    // ****************************** LOOP OVER FRS MHTDC EVENTS ***************************
    // *************************************************************************************

    // for (auto const& frs_multi_item : FrsMultiItem) {

    //   std::vector<Long64_t> const& Time_WR = frs_multi_item.;

    //   for ( const auto time_wr : Time_WR){

    //     frs_mhtdc_data.time_wr = time_wr;
    //     frs_mhtdc_data.placeholder = 0;
    //     frs_mhtdc_tree->Fill();

    //   }
      
    // }

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

  std::cout << std::endl << std::endl;

  // *************************************************************************************
  // ****************************** WRITE TREES ******************************************
  // *************************************************************************************

  std::cout << "Writing trees ..." << std::endl;
  outputFile->cd();
  
  implant_tree->Write();
  decay_tree->Write();
  bplast_tree->Write();
  germanium_tree->Write();
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
