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

  const std::map<std::string, std::string> IMPLANT_GATES_INFILE_MAP = {
    {"82nb", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/82nb.root"},
    {"84nb", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/84nb.root"},
    {"84mo", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/84mo.root"},
    {"85mo", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/85mo.root"}
  };

}

// *************************************************************************************
// ****************************** DEFINE CUSTOM STRUCTURES ***************************** 
// *************************************************************************************

// *************************************************************************************
// ****************************** DEFINE FUNCTIONS *************************************
// *************************************************************************************

// *************************************************************************************
// ****************************** START MACRO ******************************************
// *************************************************************************************


void findAidaFrsCoincidence(const char* input) {

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
  // ****************************** CREATE GATE MULTIMAPS ********************************
  // *************************************************************************************

  // Create a map with TCUTG objects of each gate and a implant event counter for each gate defined in the script
  std::map< std::string, std::tuple<TCutG*, TCutG*> > gatedimplant_cuts_map = {};

  // Loop over gate map
  for ( auto itr : constants::IMPLANT_GATES_INFILE_MAP){
    
    // Open gate file
    TFile* gate_file = new TFile(itr.second.c_str());
    if (!gate_file){
      std::cerr << "Error: Could not open gate file " << std::endl;
      std::exit(1);
    }

    std::string gate_name = itr.first; // Get gate name

    // Extract TCutG objects from gates
    TCutG *zaoq_cut = (TCutG*)gate_file->Get("cut_Z_AoQ");
    TCutG *zz2_cut = (TCutG*)gate_file->Get("cut_Z_Z2");
    
    // Add gate information to multimap
    gatedimplant_cuts_map.emplace( gate_name, std::make_tuple(zaoq_cut, zz2_cut) );   
  }

  // *************************************************************************************
  // ****************************** DEFINE TREEREADER ************************************
  // *************************************************************************************

  // Create a tree reader
  TTreeReader reader(old_tree);

  // Create tree reader values for the branches you want to read
  TTreeReaderValue<Bool_t> spill(reader, "EventHeader.fSpillFlag");
  TTreeReaderArray<Long64_t> implant_time(reader, "AidaImplantHits.Time");
  TTreeReaderArray<Int_t> implant_dssd(reader, "AidaImplantHits.DSSD");
  TTreeReaderArray<Bool_t> implant_stopped(reader, "AidaImplantHits.Stopped");
  TTreeReaderArray<Double_t> implant_energy(reader, "AidaImplantHits.Energy");
  TTreeReaderArray<Double_t> implant_energy_x(reader, "AidaImplantHits.EnergyX");
  TTreeReaderArray<Double_t> implant_energy_y(reader, "AidaImplantHits.EnergyY");
  TTreeReaderArray<Double_t> implant_x(reader, "AidaImplantHits.StripX");
  TTreeReaderArray<Double_t> implant_y(reader, "AidaImplantHits.StripY");

  TTreeReaderArray<FrsMultiHitItem> FrsMultiItem(reader, "FrsMultiHitData");
  TTreeReaderArray<FrsHitItem> FrsItem(reader, "FrsHitData");

  // *************************************************************************************
  // ****************************** DEFINE SCRIPT COUNTERS *******************************
  // *************************************************************************************
  
  // std::map<std::string, int> mhtdcFrsAidaCoincidenceCounters = {
  //   {"sci41Hits", 0},
  //   {"aidaImpMulti0", 0},
  //   {"aidaImpMulti1", 0},
  //   {"aidaImpMulti2+", 0}
  // };

  std::map<std::string, int> frsAidaCoincidenceCounters = {
    {"sci41Hits", 0},
    {"aidaImpMulti0", 0},
    {"aidaImpMulti1", 0},
    {"aidaImpMulti2+", 0}
  };

  std::map<std::string, int> nb82FrsAidaCoincidenceCounters = {
    {"sci41Hits", 0},
    {"aidaImpMulti0", 0},
    {"aidaImpMulti1", 0},
    {"aidaImpMulti2+", 0}
  };

  std::map<std::string, int> mo84FrsAidaCoincidenceCounters = {
    {"sci41Hits", 0},
    {"aidaImpMulti0", 0},
    {"aidaImpMulti1", 0},
    {"aidaImpMulti2+", 0}
  };

  std::map<std::string, int> nb84FrsAidaCoincidenceCounters = {
    {"sci41Hits", 0},
    {"aidaImpMulti0", 0},
    {"aidaImpMulti1", 0},
    {"aidaImpMulti2+", 0}
  };

  // Define time variables
  Long64_t wr_start = 0;
  Long64_t wr_end = 0;

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
    int aidaimphits = implant_time.GetSize();

    // Define Event flags
    int spflag = (*spill) ? 1 : 2;

    // Print number of entries for FRS TAC and MHTDC 
    // if (!(FrsItem.GetSize()==0 && FrsMultiItem.GetSize()==0)){
    //   std::cout << "[DEBUG] FRS TAC: " << FrsItem.GetSize() << "   ###   FRS MHTDC: " << FrsMultiItem.GetSize() << std::endl;
    // }

    // *************************************************************************************
    // ****************************** LOOP OVER FRS MHTDC  *********************************
    // *************************************************************************************

    int mhtdc_counter = 0;

    // Loop over all FrsMHTDC dataitems
    for (auto const& frs_multi_item : FrsMultiItem) {
      
      ++mhtdc_counter;
      // std::vector<Float_t> const& AoQ = frs_multi_item.Get_ID_AoQ_corr_s2s4_mhtdc();
      
      // //  Skip frs subevents where implant index wouldbe out of bounds
      // for(int idx=0; idx<AoQ.size(); ++idx){
      //   if (idx >= AoQ.size() || idx >= Z.size()) { continue; }

      //   // Skip if FRS has invalid entires
      //   if ( Z[idx] == -999.) { continue; }

      //   ++mhtdcFrsAidaCoincidenceCounters["sci41Hits"];

      //   // Switch block to handle possibilities
      //   switch (aidaimphits){
      //     case 0: 
      //       ++mhtdcFrsAidaCoincidenceCounters["aidaImpMulti0"];
      //       break;
      //     case 1: 
      //       ++mhtdcFrsAidaCoincidenceCounters["aidaImpMulti1"];
      //       break;
      //     default:
      //      ++mhtdcFrsAidaCoincidenceCounters["aidaImpMulti2+"];
      //      break;
      //   }

        // Apply 82Nb gate

      // } // End of frs subevent loop

    } // End of FRS data item loop


    // *************************************************************************************
    // ****************************** LOOP OVER FRS ****************************************
    // *************************************************************************************
    int tac_counter = 0;
    // Loop over all Frs dataitems
    for (auto const& frs_item : FrsItem) {
      ++tac_counter;
      Float_t const& AoQ = frs_item.Get_ID_AoQ_corr_s2s4();
      Float_t const& Z = frs_item.Get_ID_z41();
      Float_t const& Z2 = frs_item.Get_ID_z42();
      Long64_t time_wr = frs_item.Get_wr_t();
      
      // Get times
      if (wr_start==0 ){ wr_start = time_wr; }
      wr_end = time_wr;

      
      // Skip if FRS has invalid entires
      if ( AoQ == 0. || Z == 0.) { continue; }

      ++frsAidaCoincidenceCounters["sci41Hits"];

      // Switch block to handle possibilities
      switch (aidaimphits){
        case 0: 
          ++frsAidaCoincidenceCounters["aidaImpMulti0"];
          break;
        case 1: 
          ++frsAidaCoincidenceCounters["aidaImpMulti1"];
          break;
        default:
         ++frsAidaCoincidenceCounters["aidaImpMulti2+"];
         break;
      }

      // Extract and Apply 82Nb gate
      auto& [nb82_zaoq_cut, nb82_zz2_cut] = gatedimplant_cuts_map["82nb"];

      if (nb82_zaoq_cut->IsInside(AoQ, Z) && nb82_zz2_cut->IsInside(Z, Z2)){

        // Increase Sci41 counter
        ++nb82FrsAidaCoincidenceCounters["sci41Hits"];

        // Switch block to handle possibilities
        switch (aidaimphits){
          case 0: 
            ++nb82FrsAidaCoincidenceCounters["aidaImpMulti0"];
            break;
          case 1: 
            ++nb82FrsAidaCoincidenceCounters["aidaImpMulti1"];
            break;
          default:
           ++nb82FrsAidaCoincidenceCounters["aidaImpMulti2+"];
           break;
        }

      }

      // Extract and Apply 84Mo gate
      auto& [mo84_zaoq_cut, mo84_zz2_cut] = gatedimplant_cuts_map["84mo"];

      if (mo84_zaoq_cut->IsInside(AoQ, Z) && mo84_zz2_cut->IsInside(Z, Z2)){

        // Increase Sci41 counter
        ++mo84FrsAidaCoincidenceCounters["sci41Hits"];

        // Switch block to handle possibilities
        switch (aidaimphits){
          case 0: 
            ++mo84FrsAidaCoincidenceCounters["aidaImpMulti0"];
            break;
          case 1: 
            ++mo84FrsAidaCoincidenceCounters["aidaImpMulti1"];
            break;
          default:
           ++mo84FrsAidaCoincidenceCounters["aidaImpMulti2+"];
           break;
        }

      }

      // Extract and Apply 84Nb gate
      auto& [nb84_zaoq_cut, nb84_zz2_cut] = gatedimplant_cuts_map["84nb"];

      if (nb84_zaoq_cut->IsInside(AoQ, Z) && nb84_zz2_cut->IsInside(Z, Z2)){

        // Increase Sci41 counter
        ++nb84FrsAidaCoincidenceCounters["sci41Hits"];

        // Switch block to handle possibilities
        switch (aidaimphits){
          case 0: 
            ++nb84FrsAidaCoincidenceCounters["aidaImpMulti0"];
            break;
          case 1: 
            ++nb84FrsAidaCoincidenceCounters["aidaImpMulti1"];
            break;
          default:
           ++nb84FrsAidaCoincidenceCounters["aidaImpMulti2+"];
           break;
        }

      }




    } // End of FRS data item loop

    if (!(tac_counter==0 && mhtdc_counter==0)){
      // std::cout << "[DEBUG_KILL] FRS TAC: " << tac_counter << "   ###   FRS MHTDC: " << mhtdc_counter << std::endl;
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

  std::cout << std::endl << std::endl;

  // *************************************************************************************
  // ****************************** PRINT STATISTICS *************************************
  // *************************************************************************************

  std::cout << "###############################  Run Time  ###############################" << std::endl;
  std::cout << "File run time: " << (wr_end-wr_start)*1e-9/60 << " mins" << std::endl << std::endl;

  // std::cout << "###############################  MHTDC FRS Coincidences  ###############################" << std::endl;
  // std::cout << "Number of valid Sci41 hits: " << mhtdcFrsAidaCoincidenceCounters["sci41Hits"] << std::endl;
  // std::cout << "Number of Sci41 hits coincident w/ no AIDA implant: " << mhtdcFrsAidaCoincidenceCounters["aidaImpMulti0"] 
  //   << " (" << 100*(double)mhtdcFrsAidaCoincidenceCounters["aidaImpMulti0"]/mhtdcFrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  // std::cout << "Number of Sci41 hits coincident w/ 1 AIDA implant: " << mhtdcFrsAidaCoincidenceCounters["aidaImpMulti1"] 
  //   << " (" << 100*(double)mhtdcFrsAidaCoincidenceCounters["aidaImpMulti1"]/mhtdcFrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  // std::cout << "Number of Sci41 hits coincident w/ 2+ AIDA implant: " << mhtdcFrsAidaCoincidenceCounters["aidaImpMulti2+"] 
  //   << " (" << 100*(double)mhtdcFrsAidaCoincidenceCounters["aidaImpMulti2+"]/mhtdcFrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  // std::cout << std::endl;

  
  std::cout << "###############################  FRS Coincidences  ###############################" << std::endl;
  std::cout << "Number of valid Sci41 hits: " << frsAidaCoincidenceCounters["sci41Hits"] << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ no AIDA implant: " << frsAidaCoincidenceCounters["aidaImpMulti0"] 
    << " (" << 100*(double)frsAidaCoincidenceCounters["aidaImpMulti0"]/frsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ 1 AIDA implant: " << frsAidaCoincidenceCounters["aidaImpMulti1"] 
    << " (" << 100*(double)frsAidaCoincidenceCounters["aidaImpMulti1"]/frsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ 2+ AIDA implant: " << frsAidaCoincidenceCounters["aidaImpMulti2+"] 
    << " (" << 100*(double)frsAidaCoincidenceCounters["aidaImpMulti2+"]/frsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << std::endl;

  std::cout << "###############################  82Nb FRS Coincidences  ###############################" << std::endl;
  std::cout << "Number of valid Sci41 hits: " << nb82FrsAidaCoincidenceCounters["sci41Hits"] << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ no AIDA implant: " << nb82FrsAidaCoincidenceCounters["aidaImpMulti0"] 
    << " (" << 100*(double)nb82FrsAidaCoincidenceCounters["aidaImpMulti0"]/nb82FrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ 1 AIDA implant: " << nb82FrsAidaCoincidenceCounters["aidaImpMulti1"] 
    << " (" << 100*(double)nb82FrsAidaCoincidenceCounters["aidaImpMulti1"]/nb82FrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ 2+ AIDA implant: " << nb82FrsAidaCoincidenceCounters["aidaImpMulti2+"] 
    << " (" << 100*(double)nb82FrsAidaCoincidenceCounters["aidaImpMulti2+"]/nb82FrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << std::endl;

  std::cout << "###############################  84Mo FRS Coincidences  ###############################" << std::endl;
  std::cout << "Number of valid Sci41 hits: " << mo84FrsAidaCoincidenceCounters["sci41Hits"] << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ no AIDA implant: " << mo84FrsAidaCoincidenceCounters["aidaImpMulti0"] 
    << " (" << 100*(double)mo84FrsAidaCoincidenceCounters["aidaImpMulti0"]/mo84FrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ 1 AIDA implant: " << mo84FrsAidaCoincidenceCounters["aidaImpMulti1"] 
    << " (" << 100*(double)mo84FrsAidaCoincidenceCounters["aidaImpMulti1"]/mo84FrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ 2+ AIDA implant: " << mo84FrsAidaCoincidenceCounters["aidaImpMulti2+"] 
    << " (" << 100*(double)mo84FrsAidaCoincidenceCounters["aidaImpMulti2+"]/mo84FrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << std::endl;

  std::cout << "###############################  84Nb FRS Coincidences  ###############################" << std::endl;
  std::cout << "Number of valid Sci41 hits: " << nb84FrsAidaCoincidenceCounters["sci41Hits"] << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ no AIDA implant: " << nb84FrsAidaCoincidenceCounters["aidaImpMulti0"] 
    << " (" << 100*(double)nb84FrsAidaCoincidenceCounters["aidaImpMulti0"]/nb84FrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ 1 AIDA implant: " << nb84FrsAidaCoincidenceCounters["aidaImpMulti1"] 
    << " (" << 100*(double)nb84FrsAidaCoincidenceCounters["aidaImpMulti1"]/nb84FrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << "Number of Sci41 hits coincident w/ 2+ AIDA implant: " << nb84FrsAidaCoincidenceCounters["aidaImpMulti2+"] 
    << " (" << 100*(double)nb84FrsAidaCoincidenceCounters["aidaImpMulti2+"]/nb84FrsAidaCoincidenceCounters["sci41Hits"] << "% )" << std::endl;
  std::cout << std::endl;

  // *************************************************************************************
  // ****************************** CLEANUP **********************************************
  // *************************************************************************************

  // Close the input file
  file->Close();

  // Delete pointer objects
  delete file;
}
