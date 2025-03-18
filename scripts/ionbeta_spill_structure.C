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
  const int64_t TIME_THRESHOLD = 100 * TIME_SCALE; // Time threshold for implant beta correlation
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
std::multimap<int64_t, std::tuple<double, double, int, int, EventType>> all_implants_map;
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
  TH1F* h1_implant_decay_time_difference = new TH1F("implant_decay_time_difference", "Implant-Decay Time Difference; dt (s); Counts", 120*constants::TIME_SCALE/100, -60*constants::TIME_SCALE, 60*constants::TIME_SCALE);

  // *************************************************************************************
  // ****************************** FILL MAPS WITH EVENTS ********************************
  // *************************************************************************************

  // Loop over the implant, decay and germanium trees and fill the maps
  std::cout << "Start filled the maps!" << std::endl << std::endl;

  // Read gated implant events
  while (gatedimplant_reader.Next()){
    gated_implants_map.emplace(*gatedimplant_time, std::make_tuple(*gatedimplant_x, *gatedimplant_y, *gatedimplant_spill, *gatedimplant_dssd, GATEDIMPLANT));
  }
  std::cout << "Finished filling the gated implant map" << std::endl;
  std::cout << "Number of implant events: " << gated_implants_map.size() << std::endl << std::endl;

  // Read implant events
  while (implant_reader.Next()){
    /*if(*implant_x > 250 && *implant_x < 350 && *implant_y > 40 && *implant_y < 110){*/
      all_implants_map.emplace(*implant_time, std::make_tuple(*implant_x, *implant_y, *implant_spill, *implant_dssd, IMPLANT));
    /*}*/
  }
  std::cout << "Finished filling the implant map" << std::endl;
  std::cout << "Number of All implant events cut on region :" << all_implants_map.size() << std::endl << std::endl;

  // Read decay events
  while (decay_reader.Next()){
    /*if(*decay_x > 250 && *decay_x < 350 && *decay_y > 40 && *decay_y < 110){*/
      good_decays_map.emplace(*decay_time, std::make_tuple(*decay_x, *decay_y, *decay_dssd, *decay_spill/*, decay_bplast*/, DECAY));
    /*}*/
  }
  std::cout << "Finished filling the decay map" << std::endl;
  std::cout << "Number of Decay events :" << good_decays_map.size() << std::endl << std::endl;

  // Read germanium events
  while (germanium_reader.Next()){
    germanium_map.emplace(*germanium_time, std::make_tuple(*germanium_energy, *germanium_spill));
  }
  std::cout << "Finished filling the germanium map" << std::endl;
  std::cout << "Number of Germanium events :" << germanium_map.size() << std::endl << std::endl;


  // *************************************************************************************
  // ****************************** LOOP OVER GATED IMPLANTS *****************************
  // *************************************************************************************

  // Define & declare variebles to be used inside loop of gated implant events
  int64_t last_gatedimplant_time;
  double gatedimplant_pos_x;
  double gatedimplant_pos_y;
  bool found_forward_candidate;
  bool found_backwards_candidate;
  int implantbeta_candidate_counter;
  int matched_implantdecays_counter = 0;
  int matched_backwards_implantdecays_counter = 0;
  int matched_implantbetagamma_counter = 0;

  // Loop over all gated implant events in map and perform a beta match
  for (auto gimp_evt = gated_implants_map.begin(); gimp_evt != gated_implants_map.end(); gimp_evt++){
    
    // Unpack event variables for current gated implant event
    auto [x, y, spill, dssd, type] = gimp_evt->second;

    // Continue loop only if gated implant occured in DSSSD 1 (AIDA)
    if (type == GATEDIMPLANT && dssd == 1){

      // Check noisy implant channel strips and skip
      if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_IMPLANT, x) ){ continue; }
      if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_IMPLANT, y) ){ continue; }

      int implantbeta_candidate_counter = 0; // Reset counter to check all beta candidates
      
      // Reset flags for positive match in new loop
      bool found_forward_candidate = false;
      bool found_backwards_candidate = false;

      last_gatedimplant_time = gimp_evt->first; // Unpack white rabbit time of gated implant

      // Set gated implant position
      gatedimplant_pos_x = x;
      gatedimplant_pos_y = y;

      h2_aida_implant_xy->Fill(x,y); // Fill Histogram with gated implant position
      

      // *************************************************************************************
      // ****************************** LOOP OVER VALID DECAYS *******************************
      // *************************************************************************************

      // Find the decay event corresponding to the start of our decay loop using our time window
      // The inital decay event will be the one whose time corresponds to our time threshould before the implant occured
      auto decay_start = good_decays_map.lower_bound(last_gatedimplant_time - constants::TIME_THRESHOLD);

      // Now loop over decay events starting from our decay start defined above untoll we pass our time threshold
      for(auto decay_evt = decay_start; decay_evt != good_decays_map.end(); decay_evt++){
  
        // Break out of loop if decay events are now outside of time window
        if ( decay_evt->first > last_gatedimplant_time + constants::TIME_THRESHOLD ){ break; }

        // Break out of loop if we have found a forward  & backward match and we are not checking all candidates 
        if ( !constants::CHECK_BETA_CANDITATES && found_forward_candidate && found_backwards_candidate ){ break; }

        // Unpack event variables for current decay event
        auto [decay_x, decay_y, decay_dssd, decay_spill,/*, decay_bplast*/ decay_type] = decay_evt->second;

        // Skip if not from correct DSSD
        if ( decay_dssd == 3) { continue; }

        // Check for noisy decay branch strips and skip
        if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_DECAY, decay_x) ){ continue; }
        if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_DECAY, decay_y) ){ continue; }

        // Check if decay event is onspill and skip if defined by user and is from desired dssd
        if( constants::ONLY_OFFSPILL_DECAY && decay_spill == 1){ continue; }
        
        // Fill histogram
        h2_aida_decay_xy->Fill(decay_x, decay_y);
        /*std::cout << decay_x << " " << gatedimplant_pos_x << " " << decay_y << " " << gatedimplant_pos_y << " " << decay_type << std::endl;*/
        // Check if decay event is within position threshold
        if ( decay_type == DECAY && TMath::Abs(decay_x - gatedimplant_pos_x) <= constants::POSITION_THRESHOLD && TMath::Abs(decay_y - gatedimplant_pos_y) <= constants::POSITION_THRESHOLD ){

          // Find time difference between implant event and decay event
          int64_t time_diff = decay_evt->first - last_gatedimplant_time;
          
          // Check if decay event falls within time threshold and decay event occures after implant event 
          if (time_diff > 0 && time_diff < constants::TIME_THRESHOLD) {

            // *************************************************************************************
            // ****************************** FOUND FORWARD BETA CANDIDATE *************************
            // *************************************************************************************
            
            implantbeta_candidate_counter++; // Increase counter for beta canditade event

            // Check that there has not been a forward beta candidate that has been matched to implant event
            if (!found_forward_candidate){
              
              // Fill histograms
              h1_aida_implant_beta_dt->Fill(time_diff/constants::TIME_SCALE);
              h2_aida_matched_xy->Fill(decay_x,decay_y);
              h1_aida_wr_times->Fill(decay_evt->first);

              matched_implantdecays_counter++; // Increase counter for succesfull matched implant decay
              found_forward_candidate = true; // Change flag for succesfull forward implant decay match
              
              // Fill map with matched decay events
              matched_decays_map.emplace(decay_evt->first, std::make_tuple(decay_x, decay_y, decay_spill, DECAY));

              //************* DEBUG **************
              /*std::cout << decay_x << " " << decay_y << " " << time_diff/constants::TIME_SCALE  << std::endl;*/
              //************* DEBUG **************
            }

          }

          // Check if decay event falls within time threshold and decay event occures before implant event 
          if (time_diff < 0 && -time_diff < constants::TIME_THRESHOLD) {

            // *************************************************************************************
            // ****************************** FOUND BACKWARD BETA CANDIDATE ************************
            // *************************************************************************************

            // Check that there has not been a backwards beta candidate that has been matched to implant event
            if (!found_backwards_candidate){

              h1_aida_implant_beta_dt->Fill(time_diff/constants::TIME_SCALE); // Fill_histogram

              matched_backwards_implantdecays_counter++; // Increase counter for succesfull matched implant decay
              found_backwards_candidate = true; // Change flag for succesfull backward implant decay match
              
              //************* DEBUG **************
              /*std::cout << decay_x << " " << decay_y << " " << time_diff/constants::TIME_SCALE  << std::endl;*/
              //************* DEBUG **************
            }

          }

          // Fill the histograms with decay candidate events
          h2_implantbeta_candidate_hitpattern->Fill(decay_x - gatedimplant_pos_x, decay_y - gatedimplant_pos_y);
          h1_implantbeta_candidate_hitpattern_x->Fill(decay_x - gatedimplant_pos_x);
          h1_implantbeta_candidate_hitpattern_y->Fill(decay_y - gatedimplant_pos_y);

        }

      } // End of decay event loop
      
      // Fill candidate multiplicity histogram
      h1_implantbeta_candidate_multiplicity->Fill(implantbeta_candidate_counter); 

    }

  } // End of gated implant event loop
  
 
  // *************************************************************************************
  // ************************** MATCHED DECAY - DECAY CORRELATION ************************
  // *************************************************************************************
  
  // Loop over all matched decay events
  for( auto matched_decay_evt = matched_decays_map.begin(); matched_decay_evt != matched_decays_map.end(); matched_decay_evt++ ){
    
    // Unpack matched decay event variables
    auto [decay_x, decay_y, decay_spill, decay_type] = matched_decay_evt->second;
    int64_t last_matched_decay_time = matched_decay_evt->first;

    // Find the germanium event starting at the same time as the decay event (50 microsecond grace period)
    auto germanium_start = germanium_map.lower_bound(matched_decay_evt->first - 50e3);

    // *************************************************************************************
    // **************************  LOOP OVER GERMANIUM EVENTS ******************************
    // *************************************************************************************

    // Loop over all germanium events between decay event and end of prompt gamma window
    for( auto germanium_evt = germanium_start; germanium_evt != germanium_map.end(); germanium_evt++ ){

      /*int time_diff = germanium_evt->first - last_matched_decay_time; // Find time difference between decay and germanium event*/
      int time_diff = last_matched_decay_time - germanium_evt->first; // Reverse like in jeroens code

      if ( time_diff > constants::PROMPT_WINDOW_END ){ break; }

      // Check if germanium event is within prompt window
      if ( time_diff > constants::PROMPT_WINDOW_START && time_diff < constants::PROMPT_WINDOW_END ){
          
        // Unpack germanium event items within prompt window
        auto [germanium_energy, germanium_spill] = germanium_evt->second; 

        // Fill gamma energy spectrum and increase counter
        h1_implantbetagamma_spectrum->Fill(germanium_energy); // Fill gamma energy spectrum
        matched_implantbetagamma_counter++;
      }

    } // End of germanium event loop
      
  } // End of matched decay event loop
    
  // *************************************************************************************
  // ****************************** PRINT OUT STATISTICS *********************************
  // *************************************************************************************
   
  // Print results of implant decay correlation algorithm
    std::cout << "Finished processing the data" << std::endl;
    std::cout << "Matched: " << matched_implantdecays_counter<< " out of " << gated_implants_map.size() << " gated implant events" << std::endl;
    std::cout << "Matched: " << matched_implantdecays_counter << " out of " << all_implants_map.size() << " implant events" << std::endl;
    std::cout << "Matched: " << matched_backwards_implantdecays_counter << " backwards gated implant events" << std::endl;
    std::cout << "Matched: " << matched_implantbetagamma_counter << " implant-beta-gamma events" << std::endl << std::endl;
  
  // *************************************************************************************
  // ****************************** WRITE HISTOGRAMS TO OUTPUT FILE **********************
  // *************************************************************************************
  
  h2_aida_implant_xy->Write();
  h1_aida_wr_times->Write();
  h2_aida_matched_xy->Write();
  h2_aida_decay_xy->Write();
  h1_aida_implant_beta_dt->Write();
  if (constants::CHECK_BETA_CANDITATES){ h1_implantbeta_candidate_multiplicity->Write(); }
  if (constants::CHECK_BETA_CANDITATES){ h2_implantbeta_candidate_hitpattern->Write(); }
  if (constants::CHECK_BETA_CANDITATES){ h1_implantbeta_candidate_hitpattern_x->Write(); }
  if (constants::CHECK_BETA_CANDITATES){ h1_implantbeta_candidate_hitpattern_y->Write(); }
  h1_implantbetagamma_spectrum->Write();

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
