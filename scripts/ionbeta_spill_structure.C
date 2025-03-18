#include <iostream>
#include <map>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <string>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

// *************************************************************************************
// ****************************** DEFINE SCRIPT CONSTANTS ****************************** 
// *************************************************************************************


namespace constants{

  const int64_t TIME_SCALE = 1e9; // Timescale of time variables wrt ns
  const int64_t TIME_THRESHOLD = 100* TIME_SCALE; // Time threshold for implant beta correlation
  const int64_t POSITION_THRESHOLD = 1; //  Position window for decay wrt implant pixel as centroid

  const int LIFETIME_BINS = 300;  // Bin # used for lifetime decay plot 
  const int NEIGHBOURING_POSITION_BINS = POSITION_THRESHOLD*2+1; // Bin # used for beta candidate hit pattern histogram

  const std::vector<double> BROKEN_AIDA_X_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_X_STRIPS_DECAY = {128, 191, 192}; 
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_DECAY = {};
}

// *************************************************************************************
// ****************************** DEFINE CUSTOM STRUCTURES ***************************** 
// *************************************************************************************


enum EventType { GATEDIMPLANT, IMPLANT, DECAY }; // Tags for event type 


// *************************************************************************************
// ****************************** DEFINE MAPS FOR EVENTS *******************************
// *************************************************************************************

// Multimaps to hold events from anatrees
std::multimap<int64_t, std::tuple<double, double, int, int, EventType>> implants_map;
std::multimap<int64_t, std::tuple<double, double, int, int, EventType>> good_decays_map;
std::multimap<int64_t, std::tuple<double, int>> germanium_map;

// *************************************************************************************
// ****************************** DEFINE FUNCTIONS *************************************
// *************************************************************************************

bool isNoisyStrip(std::vector<double> noisy_strip_vector, double event_strip){
  // Function which will loop over vector containing noisy strips and check 
  // if current event occured in a noisy strip

  bool isNoisy = false; // Flag telling if event is from noisy strip

  // Loop over vector if not empty
  if (!noisy_strip_vector.empty()){
    // Loop for each noisy strip
    for (auto & noisy_strip : noisy_strip_vector){
      // If event strip corresponds to a noisy strip, change flag and break 
      if ( event_strip == noisy_strip ){ 
        isNoisy = true;
        break; 
      }

    }

  }

  return isNoisy;
}


// *************************************************************************************
// ****************************** START MACRO ******************************************
// *************************************************************************************

void ionbeta_spill_structure(const char* input, const char* output){
  
  // Load in input file 
  TFile* file = TFile::Open(input);
  if (!file){
    std::cerr << "Error: Could not open file " << std::endl;
    std::exit(1);
  }
  std::cout << "File loaded: "<< file->GetName() << std::endl << std::endl;

  
  // Get the trees
  TTree* implant_tree = (TTree*)file->Get("aida_implant_tree");
  TTree* decay_tree = (TTree*)file->Get("aida_decay_tree");
  
  // Open output file
  TFile* outputFile = new TFile(output, "RECREATE");
  if (!outputFile){
    std::cerr << "Error: Could not create output file " << std::endl;
    std::exit(1);
  }

  // Create tree readers objects
  TTreeReader implant_reader(implant_tree);
  TTreeReader decay_reader(decay_tree);

  // Define leaves of variables for implant tree
  TTreeReaderValue<ULong64_t> implant_time(implant_reader, "implant.time");
  TTreeReaderValue<double> implant_x(implant_reader, "implant.x");
  TTreeReaderValue<double> implant_y(implant_reader, "implant.y");
  TTreeReaderValue<int> implant_dssd(implant_reader, "implant.dssd");
  TTreeReaderValue<double> implant_e(implant_reader, "implant.e");
  TTreeReaderValue<double> implant_ex(implant_reader, "implant.ex");
  TTreeReaderValue<double> implant_ey(implant_reader, "implant.ey");
  TTreeReaderValue<Int_t> implant_spill(implant_reader, "implant.sp"); // sp = 1 spill, sp = 2 no spill
  /*TTreeReaderValue<Int_t> implant_bplast(implant_reader, "implant.bp"); // bp = 0 neither fired, bp = 1 only bp1 fired, bp = 2 only bp2 fired, bp = 3 both fired*/

  // Define leaves of variables for decay tree
  TTreeReaderValue<ULong64_t> decay_time(decay_reader, "decay.time");
  TTreeReaderValue<int> decay_dssd(decay_reader, "decay.dssd");
  TTreeReaderValue<double> decay_x(decay_reader, "decay.x");
  TTreeReaderValue<double> decay_y(decay_reader, "decay.y");
  TTreeReaderValue<double> decay_e(decay_reader, "decay.e");
  TTreeReaderValue<double> decay_ex(decay_reader, "decay.ex");
  TTreeReaderValue<double> decay_ey(decay_reader, "decay.ey");
  TTreeReaderValue<Int_t> decay_spill(decay_reader, "decay.sp"); // sp = 1 spill, sp = 2 no spill
  /*TTreeReaderValue<Int_t> decay_bplast(decay_reader, "decay.bp"); // bp = 0 neither fired, bp = 1 only bp1 fired, bp = 2 only bp2 fired, bp = 3 both fired*/

  // *************************************************************************************
  // ****************************** DEFINE HISTOGRAMS ************************************
  // *************************************************************************************

  // Histograms for implant beta matches
  TH1D* h1_implant_decay_time_difference = new TH1D("implant_decay_time_difference", "Implant-Decay Time Difference; dt (s); Counts", 700, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD);

  // *************************************************************************************
  // ****************************** FILL MAPS WITH EVENTS ********************************
  // *************************************************************************************

  // Loop over the implant, decay and germanium trees and fill the maps
  std::cout << "Start filled the maps!" << std::endl << std::endl;

  // Read implant events
  while (implant_reader.Next()){
    /*if(*implant_x > 250 && *implant_x < 350 && *implant_y > 40 && *implant_y < 110){*/
      implants_map.emplace(*implant_time, std::make_tuple(*implant_x, *implant_y, *implant_spill, *implant_dssd, IMPLANT));
    /*}*/
  }
  std::cout << "Finished filling the implant map" << std::endl;
  std::cout << "Number of All implant events cut on region: " << implants_map.size() << std::endl << std::endl;

  // Read decay events
  while (decay_reader.Next()){
    /*if(*decay_x > 250 && *decay_x < 350 && *decay_y > 40 && *decay_y < 110){*/
      good_decays_map.emplace(*decay_time, std::make_tuple(*decay_x, *decay_y, *decay_dssd, *decay_spill/*, decay_bplast*/, DECAY));
    /*}*/
  }
  std::cout << "Finished filling the decay map" << std::endl;
  std::cout << "Number of Decay events: " << good_decays_map.size() << std::endl << std::endl;

  // *************************************************************************************
  // ****************************** LOOP OVER GATED IMPLANTS *****************************
  // *************************************************************************************

  // Define & declare variebles to be used inside loop of gated implant events
  int64_t last_implant_time;
  double implant_pos_x;
  double implant_pos_y;
  int implant_counter = 0;

  // Loop over all gated implant events in map and perform a beta match
  for (auto imp_evt = implants_map.begin(); imp_evt != implants_map.end(); imp_evt++){
    
    // Unpack event variables for current gated implant event
    auto [x, y, spill, dssd, type] = imp_evt->second;

    // Continue loop only if gated implant occured in DSSSD 1 (AIDA)
    if (type == IMPLANT && dssd == 1){

      // Check noisy implant channel strips and skip
      if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_IMPLANT, x) ){ continue; }
      if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_IMPLANT, y) ){ continue; }

      last_implant_time = imp_evt->first; // Unpack white rabbit time of gated implant

      // Set gated implant position
      implant_pos_x = x;
      implant_pos_y = y;

      // *************************************************************************************
      // ****************************** LOOP OVER VALID DECAYS *******************************
      // *************************************************************************************

      // Find the decay event corresponding to the start of our decay loop using our time window
      // The inital decay event will be the one whose time corresponds to our time threshould before the implant occured
      auto decay_start = good_decays_map.lower_bound(last_implant_time - constants::TIME_THRESHOLD);

      // Now loop over decay events starting from our decay start defined above untoll we pass our time threshold
      for(auto decay_evt = decay_start; decay_evt != good_decays_map.end(); decay_evt++){
  
        // Break out of loop if decay events are now outside of time window
        if ( decay_evt->first > last_implant_time + 60*constants::TIME_THRESHOLD ){ break; }

        // Unpack event variables for current decay event
        auto [decay_x, decay_y, decay_dssd, decay_spill,/*, decay_bplast*/ decay_type] = decay_evt->second;

        // Skip if not from correct DSSD
        if ( decay_dssd != 1) { continue; }

        // Check for noisy decay branch strips and skip
        if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_DECAY, decay_x) ){ continue; }
        if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_DECAY, decay_y) ){ continue; }

        // Check if decay event is within position threshold
        if ( decay_type == DECAY && TMath::Abs(decay_x - implant_pos_x) <= constants::POSITION_THRESHOLD && TMath::Abs(decay_y - implant_pos_y) <= constants::POSITION_THRESHOLD ){

          // Find time difference between implant event and decay event
          int64_t time_diff = decay_evt->first - last_implant_time;
          
          // Check if decay event falls within time threshold and decay event occures after implant event 
          if (time_diff > 0 && time_diff < constants::TIME_THRESHOLD) {

            // *************************************************************************************
            // ****************************** FOUND FORWARD BETA CANDIDATE *************************
            // *************************************************************************************
            
            h1_implant_decay_time_difference->Fill(time_diff);
          }

          // Check if decay event falls within time threshold and decay event occures before implant event 
          if (time_diff < 0 && -time_diff < constants::TIME_THRESHOLD) {

            // *************************************************************************************
            // ****************************** FOUND BACKWARD BETA CANDIDATE ************************
            // *************************************************************************************

            h1_implant_decay_time_difference->Fill(time_diff);
          }

        }

      } // End of decay event loop
      
    }

    implant_counter++;

    if (implant_counter % 25 == 0){
      std::cout << "Finished looping over implant " << implant_counter << "/" << implants_map.size() << std::endl;
    }

  } // End of gated implant event loop
  
  // *************************************************************************************
  // ****************************** WRITE HISTOGRAMS TO OUTPUT FILE **********************
  // *************************************************************************************
  
  h1_implant_decay_time_difference->Write();

  std::cout << "Finished writing the histograms" << std::endl;

  // *************************************************************************************
  // ****************************** CLEANUP **********************************************
  // *************************************************************************************

  // Close the file
  file->Close();
  outputFile->Close();

  delete file;
  delete outputFile;

  std::cout << "Finished closing the files" << std::endl;

  std::exit(0);
}
