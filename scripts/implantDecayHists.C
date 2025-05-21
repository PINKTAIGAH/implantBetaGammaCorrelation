#include<iostream>
#include<map>
#include<unordered_map>
#include<tuple>
#include<utility>
#include<string>
#include<vector>

#include<TH1F.h>
#include<TH2F.h>
#include<THStack.h>
#include<TFile.h>
#include<TTree.h>
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
  const int64_t TIME_THRESHOLD = 50 * TIME_SCALE; // Time threshold for implant beta correlation
  const int64_t TIME_PER_BIN = 1e9; // Time per bin in Implant-beta time correlation
  const int64_t IMPDECAY_TIME_BINS = TIME_THRESHOLD/TIME_PER_BIN;

  const int64_t POSITION_THRESHOLD = 1; //  Position window for decay wrt implant pixel as centroid
  const int NEIGHBOURING_POSITION_BINS = POSITION_THRESHOLD*2+1; // Bin # used for beta candidate hit pattern histograms

  const std::vector<double> BROKEN_AIDA_X_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_X_STRIPS_DECAY = {63, 64, 66, 130, 189, 194, 225, 256, 319, 320}; 
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_DECAY = {};

}

namespace experimentInfo{
    const uint64_t WR_EXPERIMENT_START = 1.7401830e+18; // White rabbit start time of files 
    const uint64_t WR_EXPERIMENT_END = 1.74018310e+18; // White rabbit end time of files,
    /*const uint64_t WR_EXPERIMENT_END = 1.74022529e+18; // White rabbit end time of files,*/
    const double SLICES_EVERY = 700e3; // Size of white rabbit histogram bins 
    const int64_t DURATION_IN_SECONDS = (WR_EXPERIMENT_END - WR_EXPERIMENT_START); // Duration of experiment
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
std::multimap<int64_t, std::tuple<double, double, double, double, double, int, int, int, EventType>> gated_implants_map;
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

void implantDecayHists(const char* input, const char* output){
  
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
  TH1F* h1_implant_event_rate = new TH1F("implant_event_rate", "Implant Event Rate", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1_decay_event_rate = new TH1F("decay_event_rate", "Decay Event Rate", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);

  TH1F* h1_postcut_implant_event_rate = new TH1F("postcut_implant_event_rate", "Implant Event Rate (Post Cut)", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1_postcut_decay_event_rate = new TH1F("postcut_decay_event_rate", "Decay Event Rate (Post Cut)", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1D* h1_postcut_decay_event_rate_onspill = new TH1D("postcut_decay_event_rate_onspill", "Decay Event Rate Onspill (Post Cut)", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1D* h1_postcut_decay_event_rate_offspill = new TH1D("postcut_decay_event_rate_offspill", "Decay Event Rate Offspill (Post Cut)", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  THStack* sh1_decay_event_rate_onoffspill = new THStack("decay_event_rate_onoffspill", "");
  

  TH2F* h2_implant_strip_energy = new TH2F("implant_strip_energy", "Implant Strip vs Energy Matrix", 528, 0, 528, 7000/20, 0, 7000);
  TH2F* h2_decay_strip_energy = new TH2F("decay_strip_energy", "Decay Strip vs Energy Matrix", 528, 0, 528, 1500/20, 0, 1500);

  TH2F* h2_implant_xy_energy = new TH2F("implant_xy_energy", "Implant XY Energy", 7000/20, 0, 7000, 7000/20, 0, 7000);
  TH2F* h2_decay_xy_energy = new TH2F("decay_xy_energy", "Decay XY Energy", 1500/20, 0, 1500, 1500/20, 0, 1500);
  TH2F* h2_implant_xy = new TH2F("implant_xy", "Implant XY", 384, 0, 384, 128, 0, 128);
  TH2F* h2_decay_xy = new TH2F("decay_xy", "Decay XY", 384, 0, 384, 128, 0, 128);

  TH1F* h1_implant_energy_all = new TH1F("implant_energy_all", "Implant Energy (Stopped + Punch-through)", 7000/20, 0, 7000);
  TH1F* h1_implant_energy_stopped = new TH1F("implant_energy_stopped", "Implant Energy (Stopped)", 7000/20, 0, 7000);
  THStack* sh1_implant_energy = new THStack("implant_energy", "Implant Energy (Stopped & Punch-through)");
  TH1F* h1_gatedimplant_energy_all = new TH1F("gatedimplant_energy_all", (constants::ISOTOPE_TREE+std::string(" Implant Energy (Stopped + Punch-through)") ).c_str(), 7000/20, 0, 7000);
  TH1F* h1_gatedimplant_energy_stopped = new TH1F("gatedimplant_energy_stopped", (constants::ISOTOPE_TREE+std::string(" Implant Energy (Stopped + Punch-through)") ).c_str(), 7000/20, 0, 7000);
  THStack* sh1_gatedimplant_energy = new THStack("gatedimplant_energy", (constants::ISOTOPE_TREE+std::string(" Implant Energy (Stopped & Punch-through)") ).c_str());

  TH1F* h1_decay_energy = new TH1F("decay_energy", "Decay Energy", 1500/20, 0, 1500);

  TH2F* h2_matchedimplant_xy_short_dt = new TH2F("matchedimplant_xy_short_dt", "Matched Implant XY (dt #leq 1s)", 384, 0, 384, 128, 0, 128);
  TH2F* h2_matchedimplant_xy_long_dt = new TH2F("matchedimplant_xy_long_dt", "Matched Implant XY (dt #leq 20s)", 384, 0, 384, 128, 0, 128);

  TH1F* h1_decay_xy_dt = new TH1F("decay_xy_dt", "Decay XY Event dt", 88, -4400, 4400);

  TH2F* h2_implant_energy_dt = new TH2F("implant_energy_dt", "Implant Energy vs Implant-Decay dt", 7000/20, 0, 7000, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_decay_energy_dt = new TH2F("decay_energy_dt", "Decay Energy vs Implant-Decay dt", 1500/20, 0, 1500, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_strip_dt = new TH2F("strip_dt", "Strip XY vs Implant-Decay dt", 528, 0, 528, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);

  TH2F* h2_decay_energy_de_dt = new TH2F("decay_energy_de_dt", "Decay XY Energy difference vs Implant-Decay dt", 168/4, 0, 168, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD); 

  TH2F* h2_impdecay_deadtime_hitpattern_leq_500e6 = new TH2F("impdecay_deadtime_hitpattern_leq_500e6", "Implant - Decay dX vs dY (dT <= 500e6 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_leq_500e3 = new TH2F("impdecay_deadtime_hitpattern_leq_500e3", "Implant - Decay dX vs dY (dT <= 500e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_leq_450e3 = new TH2F("impdecay_deadtime_hitpattern_leq_450e3", "Implant - Decay dX vs dY (dT <= 450e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_leq_400e3 = new TH2F("impdecay_deadtime_hitpattern_leq_400e3", "Implant - Decay dX vs dY (dT <= 400e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_leq_350e3 = new TH2F("impdecay_deadtime_hitpattern_leq_350e3", "Implant - Decay dX vs dY (dT <= 350e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_leq_300e3 = new TH2F("impdecay_deadtime_hitpattern_leq_300e3", "Implant - Decay dX vs dY (dT <= 300e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_leq_250e3 = new TH2F("impdecay_deadtime_hitpattern_leq_250e3", "Implant - Decay dX vs dY (dT <= 250e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_leq_200e3 = new TH2F("impdecay_deadtime_hitpattern_leq_200e3", "Implant - Decay dX vs dY (dT <= 200e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_leq_150e3 = new TH2F("impdecay_deadtime_hitpattern_leq_150e3", "Implant - Decay dX vs dY (dT <= 150e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_leq_100e3 = new TH2F("impdecay_deadtime_hitpattern_leq_100e3", "Implant - Decay dX vs dY (dT <= 100e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_leq_50e3 = new TH2F("impdecay_deadtime_hitpattern_leq_50e3", "Implant - Decay dX vs dY (dT <= 50e3 ns)", 768, -384, 384, 256, -128, 128);

  TH2F* h2_impdecay_deadtime_hitpattern_union_500e6 = new TH2F("impdecay_deadtime_hitpattern_union_500e6", "Implant - Decay dX vs dY (450e6 <= dT <= 500e6 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_union_500e3 = new TH2F("impdecay_deadtime_hitpattern_union_500e3", "Implant - Decay dX vs dY (450e3 <= dT <= 500e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_union_450e3 = new TH2F("impdecay_deadtime_hitpattern_union_450e3", "Implant - Decay dX vs dY (400e3 <= dT <= 450e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_union_400e3 = new TH2F("impdecay_deadtime_hitpattern_union_400e3", "Implant - Decay dX vs dY (350e3 <= dT <= 400e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_union_350e3 = new TH2F("impdecay_deadtime_hitpattern_union_350e3", "Implant - Decay dX vs dY (300e3 <= dT <= 350e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_union_300e3 = new TH2F("impdecay_deadtime_hitpattern_union_300e3", "Implant - Decay dX vs dY (250e3 <= dT <= 300e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_union_250e3 = new TH2F("impdecay_deadtime_hitpattern_union_250e3", "Implant - Decay dX vs dY (200e3 <= dT <= 250e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_union_200e3 = new TH2F("impdecay_deadtime_hitpattern_union_200e3", "Implant - Decay dX vs dY (150e3 <= dT <= 200e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_union_150e3 = new TH2F("impdecay_deadtime_hitpattern_union_150e3", "Implant - Decay dX vs dY (100e3 <= dT <= 150e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_union_100e3 = new TH2F("impdecay_deadtime_hitpattern_union_100e3", "Implant - Decay dX vs dY (50e3 <= dT <= 100e3 ns)", 768, -384, 384, 256, -128, 128);
  TH2F* h2_impdecay_deadtime_hitpattern_union_50e3 = new TH2F("impdecay_deadtime_hitpattern_union_50e3", "Implant - Decay dX vs dY (dT <= 50e3 ns)", 768, -384, 384, 256, -128, 128);

  // *************************************************************************************
  // ****************************** FILL MAPS WITH EVENTS ********************************
  // *************************************************************************************

  // Loop over the implant, decay and germanium trees and fill the maps
  std::cout << "Start filled the maps!" << std::endl << std::endl;

  // Read gated implant events
  while (gatedimplant_reader.Next()){
    h1_gatedimplant_energy_all->Fill(*gatedimplant_e);
    if( *gatedimplant_dssd==constants::DSSD && *gatedimplant_bplast==0 ){
      h1_gatedimplant_energy_stopped->Fill(*gatedimplant_e);
      gated_implants_map.emplace(*gatedimplant_time, std::make_tuple(*gatedimplant_x, *gatedimplant_y, *gatedimplant_e, *gatedimplant_ex, *gatedimplant_ey, *gatedimplant_spill, *gatedimplant_bplast, *gatedimplant_dssd, IMPLANT));
    }
  }
  std::cout << "Finished filling the gated implant map" << std::endl;
  std::cout << "Number of implant events: " << gated_implants_map.size() << std::endl << std::endl;

  // Read implant events
  while (implant_reader.Next()){
    h1_implant_event_rate->Fill(*implant_time);
    h1_implant_energy_all->Fill(*implant_e);
    if( *implant_dssd==constants::DSSD && *implant_bplast==0 ){
      h2_implant_xy->Fill(*implant_x, *implant_y);
      h1_postcut_implant_event_rate->Fill(*implant_time);
      h2_implant_strip_energy->Fill(*implant_x, *implant_e);
      h2_implant_strip_energy->Fill(*implant_y+400, *implant_e);
      h2_implant_xy_energy->Fill(*implant_ex, *implant_ey);
      h1_implant_energy_stopped->Fill(*implant_e);
      all_implants_map.emplace(*implant_time, std::make_tuple(*implant_x, *implant_y,*implant_e, *implant_ex, *implant_ey,  *implant_spill, *implant_bplast, *implant_dssd, IMPLANT));
    }
  }
  std::cout << "Finished filling the implant map" << std::endl;
  std::cout << "Number of All implant events cut on region: " << all_implants_map.size() << std::endl << std::endl;

  // Read decay events
  while (decay_reader.Next()){
    h1_decay_event_rate->Fill(*decay_time);
    if( *decay_dssd==constants::DSSD && TMath::Abs((int64_t)*decay_time_x-(int64_t)*decay_time_y)<5e3 && TMath::Abs(*decay_ex-*decay_ey)<168 && *decay_e>151 && *decay_e<3000 ){
      h2_decay_xy->Fill(*decay_x, *decay_y);
      h1_postcut_decay_event_rate->Fill(*decay_time);
      h2_decay_strip_energy->Fill(*decay_x, *decay_e);
      h2_decay_strip_energy->Fill(*decay_y+400, *decay_e);
      h2_decay_xy_energy->Fill(*decay_ex, *decay_ey);
      h1_decay_energy->Fill(*decay_e);
      h1_decay_xy_dt->Fill(*decay_time_x-*decay_time_y);
      good_decays_map.emplace(*decay_time, std::make_tuple(*decay_x, *decay_y,*decay_e, *decay_ex, *decay_ey,  *decay_dssd, *decay_spill, DECAY));

      if (*decay_spill==1){h1_postcut_decay_event_rate_onspill->Fill(*decay_time);}
      if (*decay_spill==2){h1_postcut_decay_event_rate_offspill->Fill(*decay_time);}
    }
  }
  std::cout << "Finished filling the decay map" << std::endl;
  std::cout << "Number of Decay events: " << good_decays_map.size() << std::endl << std::endl;

  // *************************************************************************************
  // ****************************** ANLYSIS **********************************************
  // *************************************************************************************

  // Define & declare variebles to be used inside loop of gated implant events
  int64_t last_gatedimplant_time;
  double gatedimplant_pos_x;
  double gatedimplant_pos_y;
  int implantbeta_candidate_counter;
  int matched_implantdecays_counter = 0;
  int matched_backwards_implantdecays_counter = 0;


  // *************************************************************************************
  // ****************************** LOOP OVER IMPLANTS **********************************
  // *************************************************************************************

  // Loop over all gated implant events in map and perform a beta match
  for (auto gimp_evt = gated_implants_map.begin(); gimp_evt != gated_implants_map.end(); gimp_evt++){
    
    // Unpack event variables for current gated implant event
    auto [x, y, e, ex, ey, spill, bplast, dssd, type] = gimp_evt->second;

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
        auto [decay_x, decay_y, decay_e, decay_ex, decay_ey, decay_dssd, decay_spill, decay_type] = decay_evt->second;

        // Skip if not from correct DSSD
        if ( decay_dssd != constants::DSSD) { continue; }

        // Check for noisy decay branch strips and skip
        if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_DECAY, decay_x) ){ continue; }
        if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_DECAY, decay_y) ){ continue; }

        // Check if decay event is onspill and skip if defined by user and is from desired dssd
        if( constants::ONLY_OFFSPILL_DECAY && decay_spill == 1){ continue; }
        
        // Check if decay event is within position threshold
        if ( decay_type == DECAY && TMath::Abs(decay_x - gatedimplant_pos_x) <= constants::POSITION_THRESHOLD && TMath::Abs(decay_y - gatedimplant_pos_y) <= constants::POSITION_THRESHOLD ){

          // Find time difference between implant event and decay event
          double time_diff = decay_evt->first - last_gatedimplant_time;
          
          // Check if decay event falls within time threshold and decay event occures after implant event 
          if (time_diff > 0 && time_diff < constants::TIME_THRESHOLD) {

            // *************************************************************************************
            // ****************************** FOUND FORWARD BETA CANDIDATE *************************
            // *************************************************************************************
            implantbeta_candidate_counter++; // Increase counter for beta canditade event

            if (time_diff < 1e9) {h2_matchedimplant_xy_short_dt->Fill(decay_x, decay_y);}
            if (time_diff < 20e9) {h2_matchedimplant_xy_long_dt->Fill(decay_x, decay_y);}
            h2_implant_energy_dt->Fill(e, time_diff);
            h2_decay_energy_dt->Fill(decay_e, time_diff);
            h2_strip_dt->Fill(decay_x, time_diff);
            h2_strip_dt->Fill(decay_y+400, time_diff);
            h2_decay_energy_de_dt->Fill(TMath::Abs(decay_ex-decay_ey), time_diff);

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
  // ****************************** IMPLANT INTERUPTIONS *********************************
  // *************************************************************************************

  int64_t interruption_time_scale = 1e9;
  int64_t interruption_time_threshold = 50 * interruption_time_scale;
  int interruption_pos_threshold = 1;
  double interruption_binwidth = 2e6;
  int interruption_counter = 0;
  
  TH2F* h2_gimplant_xy = new TH2F("gimplant_xy", "Gated Implant Hit Pattern", 384, 0, 384, 128, 0, 128);
  TH2F* h2_gimplant_interuption_pixels = new TH2F("gimplant_interuption_pixels", Form("XY of subsequent Gated Implants Occuring within %.1e ns of an implant XY", (double)interruption_time_threshold), 384, 0, 384, 128, 0, 128);
  TH1F* h1_gimplant_interuption_dt = new TH1F("gimplant_interuption_dt", Form("dT of Subsequent GatedImplants Occuring within %.1e ns of an implant;dt;%.1e", (double)interruption_time_threshold, interruption_binwidth), interruption_time_threshold/interruption_binwidth, 0, interruption_time_threshold);

  
  std::cout << "Started interrupted implant loop ..." << std::endl;
  
  for (auto gimp_evt = gated_implants_map.begin(); gimp_evt != gated_implants_map.end(); gimp_evt++){

    last_gatedimplant_time = gimp_evt->first; // Unpack white rabbit time of gated implant
    auto [x, y, e, ex, ey, spill, bplast, dssd, type] = gimp_evt->second;
    gatedimplant_pos_x = x;
    gatedimplant_pos_y = y;
    h2_gimplant_xy->Fill(gatedimplant_pos_x, gatedimplant_pos_y);

    for (auto gimp_second_evt = all_implants_map.begin(); gimp_second_evt!=all_implants_map.end(); gimp_second_evt++){

      int64_t current_gatedimplant_time = gimp_second_evt->first;

      if (last_gatedimplant_time > current_gatedimplant_time || last_gatedimplant_time==current_gatedimplant_time){continue;}
      if (current_gatedimplant_time > last_gatedimplant_time+interruption_time_threshold ){break;}

      auto [x_curr, y_curr, e_curr, ex_curr, ey_curr, spill_curr, bplast_curr, dssd_curr, type_curr] = gimp_second_evt->second;
      int64_t time_diff=current_gatedimplant_time-last_gatedimplant_time;

      if ( TMath::Abs(x_curr - gatedimplant_pos_x) <= interruption_pos_threshold && TMath::Abs(y_curr - gatedimplant_pos_y) <= interruption_pos_threshold ){

        if (time_diff<=interruption_time_threshold){
          interruption_counter++;
          // std::cout << "Interrupted implant WR: " << last_gatedimplant_time << " ##### Time Difference: " << time_diff/(double)interruption_time_scale << " ##### X,Y: " << TMath::Abs(x_curr-gatedimplant_pos_x) << "," << TMath::Abs(y_curr-gatedimplant_pos_y) << std::endl;
          h2_gimplant_interuption_pixels->Fill(gatedimplant_pos_x, gatedimplant_pos_y);
          h1_gimplant_interuption_dt->Fill(time_diff);
          break;
        }
      }

    } 

  }

  std::cout << "Number of interruptions found: " << interruption_counter << std::endl;
  std::cout << "Finished interrupted implant loop!" << std::endl << std::endl;


  // *************************************************************************************
  // ****************************** DEADTIME LOCALISATION ********************************
  // *************************************************************************************

  // Define & declare variebles to be used inside loop of gated implant events
  last_gatedimplant_time = 0;
  gatedimplant_pos_x = 0;
  gatedimplant_pos_y = 0;
  Long64_t max_deadtime_window = 500000;

  // *************************************************************************************
  // ****************************** LOOP OVER IMPLANTS **********************************
  // *************************************************************************************

  // Loop over all gated implant events in map and perform a beta match
  for (auto gimp_evt = all_implants_map.begin(); gimp_evt != all_implants_map.end(); gimp_evt++){
    
    // Unpack event variables for current gated implant event
    auto [x, y, e, ex, ey, spill, bplast, dssd, type] = gimp_evt->second;

    // Continue loop only if gated implant occured in DSSSD 1 (AIDA)
    if (type == IMPLANT && dssd == constants::DSSD){

      // Check noisy implant channel strips and skip
      if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_IMPLANT, x) ){ continue; }
      if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_IMPLANT, y) ){ continue; }

      int implantbeta_candidate_counter = 0; // Reset counter to check all beta candidates
      
      last_gatedimplant_time = gimp_evt->first; // Unpack white rabbit time of gated implant

      // Set gated implant position
      gatedimplant_pos_x = x;
      gatedimplant_pos_y = y;

      // *************************************************************************************
      // ****************************** LOOP OVER VALID DECAYS *******************************
      // *************************************************************************************

      // Find the decay event corresponding to the start of our decay loop using our time window
      // The inital decay event will be the one whose time corresponds to our time threshould before the implant occured
      auto decay_start = good_decays_map.lower_bound(last_gatedimplant_time - 50e3);

      // Now loop over decay events starting from our decay start defined above untoll we pass our time threshold
      for(auto decay_evt = decay_start; decay_evt != good_decays_map.end(); decay_evt++){
  
        // Break out of loop if decay events are now outside of time window
        if ( decay_evt->first > last_gatedimplant_time + max_deadtime_window ){ break; }

        // Unpack event variables for current decay event
        auto [decay_x, decay_y, decay_e, decay_ex, decay_ey, decay_dssd, decay_spill, decay_type] = decay_evt->second;

        // Skip if not from correct DSSD
        if ( decay_dssd != constants::DSSD) { continue; }

        // Check for noisy decay branch strips and skip
        if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_DECAY, decay_x) ){ continue; }
        if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_DECAY, decay_y) ){ continue; }

        // Check if decay event is onspill and skip if defined by user and is from desired dssd
        if( constants::ONLY_OFFSPILL_DECAY && decay_spill == 1){ continue; }
        
        // Find time difference between implant event and decay event
        Long64_t time_diff = decay_evt->first - last_gatedimplant_time;
        double impdecay_dx = decay_x - gatedimplant_pos_x; 
        double impdecay_dy = decay_y - gatedimplant_pos_y; 

        // Draw Hists
        h2_impdecay_deadtime_hitpattern_leq_500e3->Fill(impdecay_dx, impdecay_dy);
        if (time_diff > 450000) { h2_impdecay_deadtime_hitpattern_union_500e3->Fill(impdecay_dx, impdecay_dy);}
        if (time_diff > 450000) {continue;}
        h2_impdecay_deadtime_hitpattern_leq_450e3->Fill(impdecay_dx, impdecay_dy);
        if (time_diff > 400000) { h2_impdecay_deadtime_hitpattern_union_450e3->Fill(impdecay_dx, impdecay_dy);}
        if (time_diff > 400000) {continue;}
        h2_impdecay_deadtime_hitpattern_leq_400e3->Fill(impdecay_dx, impdecay_dy);
        if (time_diff > 350000) { h2_impdecay_deadtime_hitpattern_union_400e3->Fill(impdecay_dx, impdecay_dy);}
        if (time_diff > 350000) {continue;}
        h2_impdecay_deadtime_hitpattern_leq_350e3->Fill(impdecay_dx, impdecay_dy);
        if (time_diff > 300000) { h2_impdecay_deadtime_hitpattern_union_350e3->Fill(impdecay_dx, impdecay_dy);}
        if (time_diff > 300000) {continue;}
        h2_impdecay_deadtime_hitpattern_leq_300e3->Fill(impdecay_dx, impdecay_dy);
        if (time_diff > 250000) { h2_impdecay_deadtime_hitpattern_union_300e3->Fill(impdecay_dx, impdecay_dy);}
        if (time_diff > 250000) {continue;}
        h2_impdecay_deadtime_hitpattern_leq_250e3->Fill(impdecay_dx, impdecay_dy);
        if (time_diff > 200000) { h2_impdecay_deadtime_hitpattern_union_250e3->Fill(impdecay_dx, impdecay_dy);}
        if (time_diff > 200000) {continue;}
        h2_impdecay_deadtime_hitpattern_leq_200e3->Fill(impdecay_dx, impdecay_dy);
        if (time_diff > 150000) { h2_impdecay_deadtime_hitpattern_union_200e3->Fill(impdecay_dx, impdecay_dy);}
        if (time_diff > 150000) {continue;}
        h2_impdecay_deadtime_hitpattern_leq_150e3->Fill(impdecay_dx, impdecay_dy);
        if (time_diff > 100000) { h2_impdecay_deadtime_hitpattern_union_150e3->Fill(impdecay_dx, impdecay_dy);}
        if (time_diff > 100000) {continue;}
        h2_impdecay_deadtime_hitpattern_leq_100e3->Fill(impdecay_dx, impdecay_dy);
        if (time_diff > 50000) { h2_impdecay_deadtime_hitpattern_union_100e3->Fill(impdecay_dx, impdecay_dy);}
        if (time_diff > 50000) {continue;}
        h2_impdecay_deadtime_hitpattern_leq_50e3->Fill(impdecay_dx, impdecay_dy);
        h2_impdecay_deadtime_hitpattern_union_50e3->Fill(impdecay_dx, impdecay_dy);
        
      } // End of decay event loop
      
    }

  } // End of gated implant event loop
  
 
  // *************************************************************************************
  // *************************** DEADTIME LOCALISATION LONG TIMESCALE ********************
  // *************************************************************************************

  // Define & declare variebles to be used inside loop of gated implant events
  last_gatedimplant_time = 0;
  gatedimplant_pos_x = 0;
  gatedimplant_pos_y = 0;
  max_deadtime_window = 500000000;

  // *************************************************************************************
  // ****************************** LOOP OVER IMPLANTS **********************************
  // *************************************************************************************

  // Loop over all gated implant events in map and perform a beta match
  for (auto gimp_evt = all_implants_map.begin(); gimp_evt != all_implants_map.end(); gimp_evt++){
    
    // Unpack event variables for current gated implant event
    auto [x, y, e, ex, ey, spill, bplast, dssd, type] = gimp_evt->second;

    // Continue loop only if gated implant occured in DSSSD 1 (AIDA)
    if (type == IMPLANT && dssd == constants::DSSD){

      // Check noisy implant channel strips and skip
      if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_IMPLANT, x) ){ continue; }
      if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_IMPLANT, y) ){ continue; }

      int implantbeta_candidate_counter = 0; // Reset counter to check all beta candidates
      
      last_gatedimplant_time = gimp_evt->first; // Unpack white rabbit time of gated implant

      // Set gated implant position
      gatedimplant_pos_x = x;
      gatedimplant_pos_y = y;

      // *************************************************************************************
      // ****************************** LOOP OVER VALID DECAYS *******************************
      // *************************************************************************************

      // Find the decay event corresponding to the start of our decay loop using our time window
      // The inital decay event will be the one whose time corresponds to our time threshould before the implant occured
      auto decay_start = good_decays_map.lower_bound(last_gatedimplant_time - 50000);
      auto decay_end = good_decays_map.lower_bound(last_gatedimplant_time + max_deadtime_window);

      // Now loop over decay events starting from our decay start defined above untoll we pass our time threshold
      for(auto decay_evt = decay_start; decay_evt != decay_end; decay_evt++){
  
        // Break out of loop if decay events are now outside of time window
        if ( decay_evt->first > last_gatedimplant_time + max_deadtime_window ){ break; }

        // Unpack event variables for current decay event
        auto [decay_x, decay_y, decay_e, decay_ex, decay_ey, decay_dssd, decay_spill, decay_type] = decay_evt->second;

        // Skip if not from correct DSSD
        if ( decay_dssd != constants::DSSD) { continue; }

        // Check for noisy decay branch strips and skip
        if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_DECAY, decay_x) ){ continue; }
        if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_DECAY, decay_y) ){ continue; }

        // Check if decay event is onspill and skip if defined by user and is from desired dssd
        if( constants::ONLY_OFFSPILL_DECAY && decay_spill == 1){ continue; }
        
        // Find time difference between implant event and decay event
        Long64_t time_diff = decay_evt->first - last_gatedimplant_time;
        double impdecay_dx = decay_x - gatedimplant_pos_x; 
        double impdecay_dy = decay_y - gatedimplant_pos_y; 

        // Draw Hists
        h2_impdecay_deadtime_hitpattern_leq_500e6->Fill(impdecay_dx, impdecay_dy);
        if (time_diff > 490050000) { h2_impdecay_deadtime_hitpattern_union_500e6->Fill(impdecay_dx, impdecay_dy);}
        
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
  // ****************************** MAKE SUBTRACTED HISTOGRAMS ***************************
  // *************************************************************************************

  // implant energy punchthrough
  TH1F* h1_implant_energy_punchthrough = (TH1F*) h1_implant_energy_all->Clone("implant_energy_punchthrough");
  h1_implant_energy_punchthrough->SetTitle("Implant Energy (Punch-though via subtraction)");
  h1_implant_energy_punchthrough->Add(h1_implant_energy_stopped, -1);

  // gatedimplant energy punchthrough
  TH1F* h1_gatedimplant_energy_punchthrough = (TH1F*) h1_gatedimplant_energy_all->Clone("gatedimplant_energy_punchthrough");
  h1_gatedimplant_energy_punchthrough->SetTitle((constants::ISOTOPE_TREE+std::string(" Implant Energy (Stopped & Punch-through)") ).c_str());
  h1_gatedimplant_energy_punchthrough->Add(h1_gatedimplant_energy_stopped, -1);

  
  // *************************************************************************************
  // ****************************** MAKE STACKED HISTOGRAMS ******************************
  // *************************************************************************************

  // Decay Event rates
  h1_postcut_decay_event_rate_onspill->SetLineColor(kBlue);
  h1_postcut_decay_event_rate_offspill->SetLineColor(kRed);
  sh1_decay_event_rate_onoffspill->Add(h1_postcut_decay_event_rate_onspill);
  sh1_decay_event_rate_onoffspill->Add(h1_postcut_decay_event_rate_offspill);

  // Implant Energies
  h1_implant_energy_punchthrough->SetLineColor(kBlue);
  h1_implant_energy_stopped->SetLineColor(kRed);
  sh1_implant_energy->Add(h1_implant_energy_punchthrough);
  sh1_implant_energy->Add(h1_implant_energy_stopped);

  // Gated implant Energies
  h1_gatedimplant_energy_punchthrough->SetLineColor(kBlue);
  h1_gatedimplant_energy_stopped->SetLineColor(kRed);
  sh1_gatedimplant_energy->Add(h1_gatedimplant_energy_punchthrough);
  sh1_gatedimplant_energy->Add(h1_gatedimplant_energy_stopped);

  
  // *************************************************************************************
  // ****************************** MAKE AIDA DX DY PROJECTION HISTOS ********************
  // *************************************************************************************
  TH1F* h1_impdecay_deadtime_50e3_xproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_50e3->ProjectionX("impdecay_deadtime_50e3_xproj", 128-10, 128+10, "Implant - Decay dX vs dY X-Projection (0e3 <= dT <= 50e3 ns)");
  TH1F* h1_impdecay_deadtime_100e3_xproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_100e3->ProjectionX("impdecay_deadtime_100e3_xproj", 128-10, 128+10, "Implant - Decay dX vs dY X-Projection (50e3 <= dT <= 100e3 ns)");
  TH1F* h1_impdecay_deadtime_150e3_xproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_150e3->ProjectionX("impdecay_deadtime_150e3_xproj", 128-10, 128+10, "Implant - Decay dX vs dY X-Projection (100e3 <= dT <= 150e3 ns)");
  TH1F* h1_impdecay_deadtime_200e3_xproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_200e3->ProjectionX("impdecay_deadtime_200e3_xproj", 128-10, 128+10, "Implant - Decay dX vs dY X-Projection (150e3 <= dT <= 200e3 ns)");
  TH1F* h1_impdecay_deadtime_250e3_xproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_250e3->ProjectionX("impdecay_deadtime_250e3_xproj", 128-10, 128+10, "Implant - Decay dX vs dY X-Projection (200e3 <= dT <= 250e3 ns)");
  TH1F* h1_impdecay_deadtime_300e3_xproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_300e3->ProjectionX("impdecay_deadtime_300e3_xproj", 128-10, 128+10, "Implant - Decay dX vs dY X-Projection (250e3 <= dT <= 300e3 ns)");
  TH1F* h1_impdecay_deadtime_350e3_xproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_350e3->ProjectionX("impdecay_deadtime_350e3_xproj", 128-10, 128+10, "Implant - Decay dX vs dY X-Projection (300e3 <= dT <= 350e3 ns)");
  TH1F* h1_impdecay_deadtime_400e3_xproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_400e3->ProjectionX("impdecay_deadtime_400e3_xproj", 128-10, 128+10, "Implant - Decay dX vs dY X-Projection (350e3 <= dT <= 400e3 ns)");
  TH1F* h1_impdecay_deadtime_450e3_xproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_450e3->ProjectionX("impdecay_deadtime_450e3_xproj", 128-10, 128+10, "Implant - Decay dX vs dY X-Projection (400e3 <= dT <= 450e3 ns)");
  TH1F* h1_impdecay_deadtime_500e3_xproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_500e3->ProjectionX("impdecay_deadtime_500e3_xproj", 128-10, 128+10, "Implant - Decay dX vs dY X-Projection (450e3 <= dT <= 500e3 ns)");

  TH1F* h1_impdecay_deadtime_50e3_yproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_50e3->ProjectionY("impdecay_deadtime_50e3_yproj", 384-10, 384+10, "Implant - Decay dX vs dY Y-Projection (0e3 <= dT <= 50e3 ns)");
  TH1F* h1_impdecay_deadtime_100e3_yproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_100e3->ProjectionY("impdecay_deadtime_100e3_yproj", 384-10, 384+10, "Implant - Decay dX vs dY Y-Projection (50e3 <= dT <= 100e3 ns)");
  TH1F* h1_impdecay_deadtime_150e3_yproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_150e3->ProjectionY("impdecay_deadtime_150e3_yproj", 384-10, 384+10, "Implant - Decay dX vs dY Y-Projection (100e3 <= dT <= 150e3 ns)");
  TH1F* h1_impdecay_deadtime_200e3_yproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_200e3->ProjectionY("impdecay_deadtime_200e3_yproj", 384-10, 384+10, "Implant - Decay dX vs dY Y-Projection (150e3 <= dT <= 200e3 ns)");
  TH1F* h1_impdecay_deadtime_250e3_yproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_250e3->ProjectionY("impdecay_deadtime_250e3_yproj", 384-10, 384+10, "Implant - Decay dX vs dY Y-Projection (200e3 <= dT <= 250e3 ns)");
  TH1F* h1_impdecay_deadtime_300e3_yproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_300e3->ProjectionY("impdecay_deadtime_300e3_yproj", 384-10, 384+10, "Implant - Decay dX vs dY Y-Projection (250e3 <= dT <= 300e3 ns)");
  TH1F* h1_impdecay_deadtime_350e3_yproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_350e3->ProjectionY("impdecay_deadtime_350e3_yproj", 384-10, 384+10, "Implant - Decay dX vs dY Y-Projection (300e3 <= dT <= 350e3 ns)");
  TH1F* h1_impdecay_deadtime_400e3_yproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_400e3->ProjectionY("impdecay_deadtime_400e3_yproj", 384-10, 384+10, "Implant - Decay dX vs dY Y-Projection (350e3 <= dT <= 400e3 ns)");
  TH1F* h1_impdecay_deadtime_450e3_yproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_450e3->ProjectionY("impdecay_deadtime_450e3_yproj", 384-10, 384+10, "Implant - Decay dX vs dY Y-Projection (400e3 <= dT <= 450e3 ns)");
  TH1F* h1_impdecay_deadtime_500e3_yproj = (TH1F*)h2_impdecay_deadtime_hitpattern_union_500e3->ProjectionY("impdecay_deadtime_500e3_yproj", 384-10, 384+10, "Implant - Decay dX vs dY Y-Projection (450e3 <= dT <= 500e3 ns)");

  // *************************************************************************************
  // ****************************** WRITE HISTOGRAMS TO OUTPUT FILE **********************
  // *************************************************************************************
  
  h1_implant_event_rate->Write();
  h1_decay_event_rate->Write();

  h1_postcut_implant_event_rate->Write();
  h1_postcut_decay_event_rate->Write();
  sh1_decay_event_rate_onoffspill->Write();

  h2_implant_strip_energy->Write();
  h2_decay_strip_energy->Write();

  h2_implant_xy_energy->Write();
  h2_decay_xy_energy->Write();

  h2_decay_xy->Write();
  h2_implant_xy->Write();

  h1_implant_energy_all->Write();
  h1_implant_energy_punchthrough->Write();
  h1_implant_energy_stopped->Write();
  sh1_implant_energy->Write();

  h1_gatedimplant_energy_all->Write();
  h1_gatedimplant_energy_punchthrough->Write();
  h1_gatedimplant_energy_stopped->Write();
  sh1_gatedimplant_energy->Write();

  h1_decay_energy->Write();

  h2_matchedimplant_xy_short_dt->Write();
  h2_matchedimplant_xy_long_dt->Write();

  h1_decay_xy_dt->Write();

  h2_implant_energy_dt->Write();
  h2_decay_energy_dt->Write();
  h2_strip_dt->Write();

  h2_decay_energy_de_dt->Write();

  h2_gimplant_xy->Write();
  h2_gimplant_interuption_pixels->Write();
  h1_gimplant_interuption_dt->Write();

  TDirectory* leqImpdecayDeadtime = outputFile->mkdir("impdecay_deadtime_leq");
  leqImpdecayDeadtime->cd();
  h2_impdecay_deadtime_hitpattern_leq_500e6->Write();
  h2_impdecay_deadtime_hitpattern_leq_500e3->Write();
  h2_impdecay_deadtime_hitpattern_leq_450e3->Write();
  h2_impdecay_deadtime_hitpattern_leq_400e3->Write();
  h2_impdecay_deadtime_hitpattern_leq_350e3->Write();
  h2_impdecay_deadtime_hitpattern_leq_300e3->Write();
  h2_impdecay_deadtime_hitpattern_leq_250e3->Write();
  h2_impdecay_deadtime_hitpattern_leq_200e3->Write();
  h2_impdecay_deadtime_hitpattern_leq_150e3->Write();
  h2_impdecay_deadtime_hitpattern_leq_100e3->Write();
  h2_impdecay_deadtime_hitpattern_leq_50e3->Write();
  gFile->cd();

  TDirectory* unionImpdecayDeadtime = outputFile->mkdir("impdecay_deadtime_union");
  unionImpdecayDeadtime->cd();
  h2_impdecay_deadtime_hitpattern_union_500e6->Write();
  h2_impdecay_deadtime_hitpattern_union_500e3->Write();
  h2_impdecay_deadtime_hitpattern_union_450e3->Write();
  h2_impdecay_deadtime_hitpattern_union_400e3->Write();
  h2_impdecay_deadtime_hitpattern_union_350e3->Write();
  h2_impdecay_deadtime_hitpattern_union_300e3->Write();
  h2_impdecay_deadtime_hitpattern_union_250e3->Write();
  h2_impdecay_deadtime_hitpattern_union_200e3->Write();
  h2_impdecay_deadtime_hitpattern_union_150e3->Write();
  h2_impdecay_deadtime_hitpattern_union_100e3->Write();
  h2_impdecay_deadtime_hitpattern_union_50e3->Write();
  gFile->cd();


  TDirectory* impdecayDeadtimeXProj = outputFile->mkdir("impdecay_deadtime_xproj");
  impdecayDeadtimeXProj->cd();
  h1_impdecay_deadtime_50e3_xproj->Write();
  h1_impdecay_deadtime_100e3_xproj->Write();
  h1_impdecay_deadtime_150e3_xproj->Write();
  h1_impdecay_deadtime_200e3_xproj->Write();
  h1_impdecay_deadtime_250e3_xproj->Write();
  h1_impdecay_deadtime_300e3_xproj->Write();
  h1_impdecay_deadtime_350e3_xproj->Write();
  h1_impdecay_deadtime_400e3_xproj->Write();
  h1_impdecay_deadtime_450e3_xproj->Write();
  h1_impdecay_deadtime_500e3_xproj->Write();
  gFile->cd();

  TDirectory* impdecayDeadtimeYProj = outputFile->mkdir("impdecay_deadtime_yproj");
  impdecayDeadtimeYProj->cd();
  h1_impdecay_deadtime_50e3_yproj->Write();
  h1_impdecay_deadtime_100e3_yproj->Write();
  h1_impdecay_deadtime_150e3_yproj->Write();
  h1_impdecay_deadtime_200e3_yproj->Write();
  h1_impdecay_deadtime_250e3_yproj->Write();
  h1_impdecay_deadtime_300e3_yproj->Write();
  h1_impdecay_deadtime_350e3_yproj->Write();
  h1_impdecay_deadtime_400e3_yproj->Write();
  h1_impdecay_deadtime_450e3_yproj->Write();
  h1_impdecay_deadtime_500e3_yproj->Write();
  gFile->cd();

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
