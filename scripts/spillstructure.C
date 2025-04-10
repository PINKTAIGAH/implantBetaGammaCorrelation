#include<iostream>
#include<map>
#include<unordered_map>
#include<tuple>
#include<utility>
#include<string>
#include<vector>

#include<TH1F.h>
#include<TH2F.h>
#include<TFile.h>
#include<TTree.h>
#include<THStack.h>
#include<TNtuple.h>
#include<TString.h>
#include<TDirectory.h>
#include<TTreeReader.h>
#include<TTreeReaderArray.h>
#include<TTreeReaderValue.h>

// *************************************************************************************
// ****************************** DEFINE SCRIPT CONSTANTS ****************************** 
// *************************************************************************************


namespace constants{
  const std::string ISOTOPE_TREE = "84nb"; // Name suffix for gatedimplant tree & branch in anatree
  const int DSSD = 1; // Which DSSD will the analysis be run on

  const bool ONLY_OFFSPILL_DECAY = false; // Check for onspill decay matches
  const bool CHECK_BETA_CANDITATES = true; // Check for all beta candidates of an implant
  /*const bool INCLUDE_BACKWARDS_MATCH = true; // Look for reverse time implant beta correlations*/

  const int64_t TIME_SCALE = 1e9; // Timescale of time variables wrt ns
  const int64_t TIME_THRESHOLD = 40 * TIME_SCALE; // Time threshold for implant beta correlation
  const double TIME_PER_BIN = 1e8; // Time per bin in Implant-beta time correlation
  const int64_t SPILLSTRUCTURE_BINS = TIME_THRESHOLD/TIME_PER_BIN;

  const int64_t POSITION_THRESHOLD = 1; //  Position window for decay wrt implant pixel as centroid
  const int NEIGHBOURING_POSITION_BINS = POSITION_THRESHOLD*2+1; // Bin # used for beta candidate hit pattern histograms

  const std::vector<double> BROKEN_AIDA_X_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_X_STRIPS_DECAY = {63, 64, 66, 130, 189, 194, 225, 256, 319, 320}; 
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
std::multimap<int64_t, std::tuple<double, double, int, int, int, EventType>> gated_implants_map;
std::multimap<int64_t, std::tuple<double, double, int, int, int, EventType>> all_implants_map;
std::multimap<int64_t, std::tuple<double, double, int, int, EventType>> good_decays_map;

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

void spillstructure(const char* input, const char* output){
  
  // Load in input file 
  TFile* file = TFile::Open(input);
  if (!file){
    std::cerr << "Error: Could not open file " << std::endl;
    std::exit(1);
  }
  std::cout << "File loaded: "<< file->GetName() << std::endl << std::endl;

  
  // Get the trees
  TTree* implant_tree = (TTree*)file->Get("aida_implant_tree");
  TTree* gatedimplant_tree = (TTree*)file->Get( (std::string("aida_gatedimplant_")+constants::ISOTOPE_TREE+std::string("_tree") ).c_str() );
  TTree* decay_tree = (TTree*)file->Get("aida_decay_tree");
  
  // Open output file
  TFile* outputFile = new TFile(output, "RECREATE");
  if (!outputFile){
    std::cerr << "Error: Could not create output file " << std::endl;
    std::exit(1);
  }

  // Create tree readers objects
  TTreeReader implant_reader(implant_tree);
  TTreeReader gatedimplant_reader(gatedimplant_tree);
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
  TTreeReaderValue<Int_t> implant_bplast(implant_reader, "implant.bp"); // bp = 0 neither fired, bp = 1 only bp1 fired, bp = 2 only bp2 fired, bp = 3 both fired

  // Define leaves of variables for gated implant tree
  TTreeReaderValue<ULong64_t> gatedimplant_time(gatedimplant_reader, ( std::string("gatedimplant_")+constants::ISOTOPE_TREE+std::string(".time") ).c_str());
  TTreeReaderValue<double> gatedimplant_x(gatedimplant_reader, ( std::string("gatedimplant_")+constants::ISOTOPE_TREE+std::string(".x") ).c_str());
  TTreeReaderValue<double> gatedimplant_y(gatedimplant_reader, ( std::string("gatedimplant_")+constants::ISOTOPE_TREE+std::string(".y") ).c_str());
  TTreeReaderValue<int> gatedimplant_dssd(gatedimplant_reader, ( std::string("gatedimplant_")+constants::ISOTOPE_TREE+std::string(".dssd") ).c_str());
  TTreeReaderValue<double> gatedimplant_e(gatedimplant_reader, ( std::string("gatedimplant_")+constants::ISOTOPE_TREE+std::string(".e") ).c_str());
  TTreeReaderValue<double> gatedimplant_ex(gatedimplant_reader, ( std::string("gatedimplant_")+constants::ISOTOPE_TREE+std::string(".ex") ).c_str());
  TTreeReaderValue<double> gatedimplant_ey(gatedimplant_reader, ( std::string("gatedimplant_")+constants::ISOTOPE_TREE+std::string(".ey") ).c_str());
  TTreeReaderValue<Int_t> gatedimplant_spill(gatedimplant_reader, ( std::string("gatedimplant_")+constants::ISOTOPE_TREE+std::string(".sp") ).c_str()); // sp = 1 spill, sp = 2 no spill
  TTreeReaderValue<Int_t> gatedimplant_bplast(gatedimplant_reader, ( std::string("gatedimplant_")+constants::ISOTOPE_TREE+std::string(".bp") ).c_str()); // bp = 0 neither fired, bp = 1 only bp1 fired, bp = 2 only bp2 fired, bp = 3 both fired

  // Define leaves of variables for decay tree
  TTreeReaderValue<ULong64_t> decay_time(decay_reader, "decay.time");
  TTreeReaderValue<ULong64_t> decay_time_x(decay_reader, "decay.time_x");
  TTreeReaderValue<ULong64_t> decay_time_y(decay_reader, "decay.time_y");
  TTreeReaderValue<int> decay_dssd(decay_reader, "decay.dssd");
  TTreeReaderValue<double> decay_x(decay_reader, "decay.x");
  TTreeReaderValue<double> decay_y(decay_reader, "decay.y");
  TTreeReaderValue<double> decay_e(decay_reader, "decay.e");
  TTreeReaderValue<double> decay_ex(decay_reader, "decay.ex");
  TTreeReaderValue<double> decay_ey(decay_reader, "decay.ey");
  TTreeReaderValue<Int_t> decay_spill(decay_reader, "decay.sp"); // sp = 1 spill, sp = 2 no spill

  // *************************************************************************************
  // ****************************** DEFINE HISTOGRAMS ************************************
  // *************************************************************************************

  // Histograms for implant beta matches
  TH2F* h2_aida_implant_xy = new TH2F("aida_implant_xy", "AIDA Implant XY", 384, 0, 384, 128, 0, 128);
  TH2F* h2_aida_decay_xy = new TH2F("aida_decay_xy", "AIDA Decay XY", 384, 0, 384, 128, 0, 128);
  TH2F* h2_aida_matched_onspill_xy = new TH2F("aida_matched_onspill_xy", "AIDA Matched Onspill XY", 384, 0, 384, 128, 0, 128);
  TH2F* h2_aida_matched_offspill_xy = new TH2F("aida_matched_offspill_xy", "AIDA Matched Offspill XY", 384, 0, 384, 128, 0, 128);
  TH1F* h1_aida_implant_beta_spillstructure = new TH1F("aida_implant_beta_spillstructure", "Implant-Decay Spill Structure; dT (ns); Counts", constants::SPILLSTRUCTURE_BINS, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD);
  TH1F* h1_aida_implant_beta_onspillstructure = new TH1F("aida_implant_beta_onspillstructure", "Implant-Decay Spill Structure Onspill; dT (ns); Counts", constants::SPILLSTRUCTURE_BINS, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD);
  TH1F* h1_aida_implant_beta_offspillstructure = new TH1F("aida_implant_beta_offspillstructure", "Implant-Decay Spill Structure Offspill; dT (ns); Counts", constants::SPILLSTRUCTURE_BINS, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD);
  THStack* sh1_onoff_spillstructure = new THStack("onoff_spillstructure", "");


  // *************************************************************************************
  // ****************************** FILL MAPS WITH EVENTS ********************************
  // *************************************************************************************

  // Loop over the implant, decay and germanium trees and fill the maps
  std::cout << "Start filled the maps!" << std::endl << std::endl;

  // Read gated implant events
  while (gatedimplant_reader.Next()){
    if( *gatedimplant_dssd==constants::DSSD && *gatedimplant_bplast==0 ){
      gated_implants_map.emplace(*gatedimplant_time, std::make_tuple(*gatedimplant_x, *gatedimplant_y, *gatedimplant_spill, *gatedimplant_bplast, *gatedimplant_dssd, IMPLANT));
    }
  }
  std::cout << "Finished filling the gated implant map" << std::endl;
  std::cout << "Number of implant events: " << gated_implants_map.size() << std::endl << std::endl;

  // Read implant events
  while (implant_reader.Next()){
    if( *implant_dssd==constants::DSSD && *implant_bplast==0 ){
      all_implants_map.emplace(*implant_time, std::make_tuple(*implant_x, *implant_y, *implant_spill, *implant_bplast, *implant_dssd, IMPLANT));
    }
  }
  std::cout << "Finished filling the implant map" << std::endl;
  std::cout << "Number of All implant events cut on region: " << all_implants_map.size() << std::endl << std::endl;

  // Read decay events
  while (decay_reader.Next()){
    if( *decay_dssd==constants::DSSD && TMath::Abs( (int64_t)(*decay_time_x-*decay_time_y) )<2.2e3 && TMath::Abs(*decay_ex-*decay_ey)<168 && *decay_e>151 && *decay_e<1000 ){
      good_decays_map.emplace(*decay_time, std::make_tuple(*decay_x, *decay_y, *decay_dssd, *decay_spill, DECAY));
    }
  }
  std::cout << "Finished filling the decay map" << std::endl;
  std::cout << "Number of Decay events: " << good_decays_map.size() << std::endl << std::endl;

  // *************************************************************************************
  // ****************************** LOOP OVER GATED IMPLANTS *****************************
  // *************************************************************************************

  // Define & declare variebles to be used inside loop of gated implant events
  int64_t last_gatedimplant_time;
  double gatedimplant_pos_x;
  double gatedimplant_pos_y;
  int implantbeta_candidate_counter;
  int matched_implantdecays_counter = 0;
  int matched_backwards_implantdecays_counter = 0;

  // Loop over all gated implant events in map and perform a beta match
  for (auto gimp_evt = all_implants_map.begin(); gimp_evt != all_implants_map.end(); gimp_evt++){
    
    // Unpack event variables for current gated implant event
    auto [x, y, spill, bplast, dssd, type] = gimp_evt->second;

    // Continue loop only if gated implant occured in DSSSD 1 (AIDA)
    if (type == IMPLANT && dssd == constants::DSSD){

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

        // Unpack event variables for current decay event
        auto [decay_x, decay_y, decay_dssd, decay_spill, decay_type] = decay_evt->second;

        // Skip if not from correct DSSD
        if ( decay_dssd != constants::DSSD) { continue; }

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
            //
            if ( decay_spill == 1 ){
              h1_aida_implant_beta_onspillstructure->Fill(time_diff); // Fill spillstructure histogram
              h2_aida_matched_onspill_xy->Fill(decay_x,decay_y);
            }
            if ( decay_spill == 2 ){
              h1_aida_implant_beta_offspillstructure->Fill(time_diff); // Fill spillstructure histogram
              h2_aida_matched_offspill_xy->Fill(decay_x,decay_y);
            }

            h1_aida_implant_beta_spillstructure->Fill(time_diff);

            // Check that there has not been a forward beta candidate that has been matched to implant event
            if (!found_forward_candidate){

              matched_implantdecays_counter++; // Increase counter for succesfull matched implant decay
              found_forward_candidate = true; // Change flag for succesfull forward implant decay match
            }

          }

          // Check if decay event falls within time threshold and decay event occures before implant event 
          if (time_diff < 0 && -time_diff < constants::TIME_THRESHOLD) {

            // *************************************************************************************
            // ****************************** FOUND BACKWARD BETA CANDIDATE ************************
            // *************************************************************************************
            
            if ( decay_spill == 1 ){
              h1_aida_implant_beta_onspillstructure->Fill(time_diff); // Fill spillstructure histogram
              h2_aida_matched_onspill_xy->Fill(decay_x,decay_y);
            }
            if ( decay_spill == 2 ){
              h1_aida_implant_beta_offspillstructure->Fill(time_diff); // Fill spillstructure histogram
              h2_aida_matched_offspill_xy->Fill(decay_x,decay_y);
            }

            h1_aida_implant_beta_spillstructure->Fill(time_diff);

            // Check that there has not been a backwards beta candidate that has been matched to implant event
            if (!found_backwards_candidate){

              matched_backwards_implantdecays_counter++; // Increase counter for succesfull matched implant decay
              found_backwards_candidate = true; // Change flag for succesfull backward implant decay match
            }

          }

        }

      } // End of decay event loop
      
    }

  } // End of gated implant event loop
  
 
    
  // *************************************************************************************
  // ****************************** PRINT OUT STATISTICS *********************************
  // *************************************************************************************
   
  // Print results of implant decay correlation algorithm
    std::cout << "Finished processing the data" << std::endl;
    std::cout << "Matched: " << matched_implantdecays_counter<< " out of " << gated_implants_map.size() << " gated implant events" << std::endl;
    std::cout << "Matched: " << matched_implantdecays_counter << " out of " << all_implants_map.size() << " implant events" << std::endl;
    std::cout << "Matched: " << matched_backwards_implantdecays_counter << " backwards gated implant events" << std::endl << std::endl;
  
  // *************************************************************************************
  // ****************************** WRITE HISTOGRAMS TO OUTPUT FILE **********************
  // *************************************************************************************
  
  h2_aida_implant_xy->Write();
  h2_aida_decay_xy->Write();
  h2_aida_matched_onspill_xy->Write();
  h2_aida_matched_offspill_xy->Write();
  h1_aida_implant_beta_spillstructure->Write();
  h1_aida_implant_beta_onspillstructure->Write();
  h1_aida_implant_beta_offspillstructure->Write();

  std::cout << "Finished writing the histograms" << std::endl;


  // *************************************************************************************
  // ****************************** WRITE MERGED TO OUTPUT FILE **********************
  // *************************************************************************************
  
  h1_aida_implant_beta_onspillstructure->SetLineColor(kBlue);
  /*h1_aida_implant_beta_onspillstructure->SetLineWidth(3);*/
  sh1_onoff_spillstructure->Add(h1_aida_implant_beta_onspillstructure);

  h1_aida_implant_beta_offspillstructure->SetLineColor(kRed);
  /*h1_aida_implant_beta_offspillstructure->SetLineWidth(3);*/
  sh1_onoff_spillstructure->Add(h1_aida_implant_beta_offspillstructure);

  sh1_onoff_spillstructure->Write();

  std::cout << "Finished writing the canvases" << std::endl;

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
