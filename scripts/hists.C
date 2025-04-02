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
  const std::string ISOTOPE_TREE = "85mo"; // Name suffix for gatedimplant tree & branch in anatree
  const int DSSD = 3; // Which DSSD will the analysis be run on

  const bool ONLY_OFFSPILL_DECAY = false; // Check for onspill decay matches
  const bool CHECK_BETA_CANDITATES = true; // Check for all beta candidates of an implant
  /*const bool INCLUDE_BACKWARDS_MATCH = true; // Look for reverse time implant beta correlations*/

  const int64_t TIME_SCALE = 1e9; // Timescale of time variables wrt ns
  const int64_t TIME_THRESHOLD = 120 * TIME_SCALE; // Time threshold for implant beta correlation
  const int64_t TIME_PER_BIN = 4e6; // Time per bin in Implant-beta time correlation
  const int64_t IMPDECAY_TIME_BINS = TIME_THRESHOLD/TIME_PER_BIN;

  const int64_t POSITION_THRESHOLD = 1; //  Position window for decay wrt implant pixel as centroid
  const int NEIGHBOURING_POSITION_BINS = POSITION_THRESHOLD*2+1; // Bin # used for beta candidate hit pattern histograms

  const std::vector<double> BROKEN_AIDA_X_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_X_STRIPS_DECAY = {128, 191, 192}; 
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_DECAY = {};
}

namespace experimentInfo{
    /*const uint64_t WR_EXPERIMENT_START = 1.7401830e+18; // RUN 21 */
    /*const uint64_t WR_EXPERIMENT_END = 1740190296922399672; // RUN 21 */
    const uint64_t WR_EXPERIMENT_START = 1740223604000299360; // RUN 32 
    const uint64_t WR_EXPERIMENT_END = 1740225294304511440; // RUN 32
    const double SLICES_EVERY = 0.1; // Size of white rabbit histogram bins 
    const int64_t DURATION_IN_SECONDS = (WR_EXPERIMENT_END - WR_EXPERIMENT_START)/1e9; // Duration of experiment
    const int64_t NUMBER_OF_SLICES = DURATION_IN_SECONDS/SLICES_EVERY;  // Number of white rabbit histogram bins
}

// *************************************************************************************
// ****************************** DEFINE CUSTOM STRUCTURES ***************************** 
// *************************************************************************************


enum EventType { GATEDIMPLANT, IMPLANT, DECAY }; // Tags for event type 


// *************************************************************************************
// ****************************** DEFINE MAPS FOR EVENTS *******************************
// *************************************************************************************

// Multimaps to hold events from anatrees
std::multimap<int64_t, std::tuple<double, double, double, double, double, int, int, int, EventType>> all_implants_map;
std::multimap<int64_t, std::tuple<double, double, double, double, double, int, int, EventType>> good_decays_map;

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

void hists(const char* input, const char* output){
  
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
  TTreeReaderValue<Int_t> implant_bplast(implant_reader, "implant.bp"); // bp = 0 neither fired, bp = 1 only bp1 fired, bp = 2 only bp2 fired, bp = 3 both fired

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
  TH1F* h1_implant_eventrate = new TH1F("implant_eventrate", "Implant Event Rate", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1_implant_eventrate_onspill = new TH1F("implant_eventrate_onspill", "Implant Event Rate Onspill", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1_implant_eventrate_offspill = new TH1F("implant_eventrate_offspill", "Implant Event Rate Offspill", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1_decay_eventrate = new TH1F("decay_eventrate", "Decay Event Rate", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1_decay_eventrate_onspill = new TH1F("decay_eventrate_onspill", "Decay Event Rate Onspill", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1_decay_eventrate_offspill = new TH1F("decay_eventrate_offspill", "Decay Event Rate Offspill", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);

  // *************************************************************************************
  // ****************************** FILL MAPS WITH EVENTS ********************************
  // *************************************************************************************

  // Read implant events
  while (implant_reader.Next()){
    if( *implant_dssd==constants::DSSD && *implant_bplast==0 ){
      all_implants_map.emplace(*implant_time, std::make_tuple(*implant_x, *implant_y, *implant_e, *implant_ex, *implant_ey,  *implant_spill, *implant_bplast, *implant_dssd, IMPLANT));
    }
  }
  std::cout << "Finished filling the implant map" << std::endl;
  std::cout << "Number of All implant events cut on region: " << all_implants_map.size() << std::endl << std::endl;

  // Read decay events
  while (decay_reader.Next()){
    if( *decay_dssd==constants::DSSD && TMath::Abs( (int64_t)(*decay_time_x-*decay_time_y) )<5e3 && TMath::Abs(*decay_ex-*decay_ey)<168 && *decay_e>151 && *decay_e<1000 ){
      good_decays_map.emplace(*decay_time, std::make_tuple(*decay_x, *decay_y,*decay_e, *decay_ex, *decay_ey,  *decay_dssd, *decay_spill, DECAY));
    }
  }
  std::cout << "Finished filling the decay map" << std::endl;
  std::cout << "Number of Decay events: " << good_decays_map.size() << std::endl << std::endl;

  // *************************************************************************************
  // ****************************** LOOP OVER DECAYS *************************************
  // *************************************************************************************
  for (auto decay_evt=good_decays_map.begin(); decay_evt!=good_decays_map.end(); decay_evt++){

    auto [x, y, e, ex, ey, dssd, spill, event_type] = decay_evt->second;
    auto time = decay_evt->first;
    if (spill==1){h1_decay_eventrate_onspill->Fill(time);}
    if (spill==2){h1_decay_eventrate_offspill->Fill(time);}
    h1_decay_eventrate->Fill(time);

  }

  // *************************************************************************************
  // ****************************** LOOP OVER IMPLANTS ***********************************
  // *************************************************************************************
  for (auto implant_evt=all_implants_map.begin(); implant_evt!=all_implants_map.end(); implant_evt++){

    auto [x, y, e, ex, ey, spill, bplast, dssd, event_type] = implant_evt->second;
    auto time = implant_evt->first;
    if (spill==1){h1_implant_eventrate_onspill->Fill(time);}
    if (spill==2){h1_implant_eventrate_offspill->Fill(time);}
    h1_implant_eventrate->Fill(time);

  }

  // *************************************************************************************
  // ****************************** MAKE HISTOGRAM STACKS ********************************
  // *************************************************************************************
  THStack* sh1_implant_eventrate_onoffspill = new THStack("implant_eventrate_onoffspill", "");
  THStack* sh1_decay_eventrate_onoffspill = new THStack("decay_eventrate_onoffspill", "");

  h1_implant_eventrate_onspill->SetLineColor(kBlue);
  h1_decay_eventrate_onspill->SetLineColor(kBlue);

  h1_implant_eventrate_offspill->SetLineColor(kRed);
  h1_decay_eventrate_offspill->SetLineColor(kRed);

  sh1_implant_eventrate_onoffspill->Add(h1_implant_eventrate_onspill);
  sh1_implant_eventrate_onoffspill->Add(h1_implant_eventrate_offspill);

  sh1_decay_eventrate_onoffspill->Add(h1_decay_eventrate_onspill);
  sh1_decay_eventrate_onoffspill->Add(h1_decay_eventrate_offspill);

  // *************************************************************************************
  // ****************************** WRITE HISTOGRAMS TO OUTPUT FILE **********************
  // *************************************************************************************

  h1_implant_eventrate->Write();
  h1_implant_eventrate_onspill->Write();
  h1_implant_eventrate_offspill->Write();
  sh1_implant_eventrate_onoffspill->Write();
  h1_decay_eventrate->Write();
  h1_decay_eventrate_onspill->Write();
  h1_decay_eventrate_offspill->Write();
  sh1_decay_eventrate_onoffspill->Write();

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
