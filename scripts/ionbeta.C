#include<iostream>
#include<map>
#include<unordered_map>
#include<cstdint>
#include<tuple>
#include<utility>
#include<string>
#include<vector>
#include<algorithm>

#include<TF1.h>
#include<TH1F.h>
#include<TH2F.h>
#include<TFile.h>
#include<TTree.h>
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
  const std::string ISOTOPE_TREE = "82nb"; // Name suffix for gatedimplant tree & branch in anatree
  const int DSSD = 1; // Which DSSD will the analysis be run on

  const bool ONLY_OFFSPILL_DECAY = false; // Check for onspill decay matches
  const bool CHECK_BETA_CANDITATES = true; // Check for all beta candidates of an implant
  const bool VETO_INTERRUPTED_IMPLANTS = true; // Veto any implant which has a subsequent implant occure within time and position window
  /*const bool INCLUDE_BACKWARDS_MATCH = true; // Look for reverse time implant beta correlations*/

  const int64_t TIME_SCALE = 1e6; // Timescale of time variables wrt ns
  const int64_t TIME_THRESHOLD = 500 * TIME_SCALE; // Time threshold for implant beta correlation
  const double TIME_PER_DT_BIN = 10e6; // Time per bin in Implant-beta time correlation
  const int64_t POSITION_THRESHOLD = 1; //  Position window for decay wrt implant pixel as centroid

  const uint64_t IMPLANT_DEAD_TIME = 300e3; // The deadtime we impose on aida LEC after an implant occures in AIDA
  const int BETA_CANDIDATE_CUT = 1; // Define number of candidate betas a implant must have before plotting

  /*const std::map<, int64_t> PROMPT_GAMMA_WINDOW = { {"start", 14498}, {"final", 16498} };*/
  const int64_t PROMPT_WINDOW_START = 13610; 
  const int64_t PROMPT_WINDOW_END = 16223; 

  const int NEIGHBOURING_POSITION_BINS = POSITION_THRESHOLD*2+1; // Bin # used for beta candidate hit pattern histogram
  const int64_t IMPDECAY_TIME_BINS = 2*TIME_THRESHOLD/TIME_PER_DT_BIN;

  const std::vector<double> BROKEN_AIDA_X_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_X_STRIPS_DECAY = {63, 64, 66, 130, 189, 191, 194, 225, 256, 319, 320, 192, 128, 255}; 
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_DECAY = {};


}

namespace experimentInfo{
    const uint64_t WR_EXPERIMENT_START = 1.7401830e+18; // White rabbit start time of files 
    const uint64_t WR_EXPERIMENT_END = 1.74022529e+18; // White rabbit end time of files,
    const int64_t SLICES_EVERY = 1; // Size of white rabbit histogram bins 
    const int64_t DURATION_IN_SECONDS = (WR_EXPERIMENT_END - WR_EXPERIMENT_START)/1e9; // Duration of experiment
    const int64_t NUMBER_OF_SLICES = DURATION_IN_SECONDS/SLICES_EVERY;  // Number of white rabbit histogram bins
}


// *************************************************************************************
// ****************************** DEFINE CUSTOM STRUCTURES ***************************** 
// *************************************************************************************


enum EventType { GATEDIMPLANT, IMPLANT, DECAY }; // Tags for event type 
enum CorrelationType { FORWARDS, BACKWARDS }; // Tags is ionbeta match is a borwards or backwards time correlations
enum BetaType { CANDIDATE, MATCH }; // Tags for candidate betas (Have they been matched to an implant or not)
enum InterruptedCoincidenceType { NOT_INTERRUPTED, INTERRUPTED }; // Tags to define if an implant has beenn interrupted in a time window


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
// ****************************** DEFINE CLASSES ***************************************
// *************************************************************************************


class TimeRangeManager {
  public:
    struct TimeRange {
      uint64_t start;
      uint64_t end;

      bool contains(uint64_t time) const {
        return time >= start && time <= end;
      }

      bool operator<(const TimeRange& other) const {
        return start < other.start;
      }
    };

    // Add a new time range
    void addRange(uint64_t start, uint64_t end) {
      if (start > end) std::swap(start, end);  // ensure proper order
      ranges.push_back({start, end});
      isMerged = false;
    }

    // Check if a time is within any merged range
    bool contains(uint64_t time) {
      if (!isMerged) mergeRanges();
      for (const auto& range : mergedRanges) {
        if (range.contains(time)) return true;
      }
      return false;
    }

    // Access the merged ranges
    const std::vector<TimeRange>& getMergedRanges() {
      if (!isMerged) mergeRanges();
      return mergedRanges;
    }

  private:
    std::vector<TimeRange> ranges;
    std::vector<TimeRange> mergedRanges;
    bool isMerged = false;

    // Merge overlapping and adjacent ranges
    void mergeRanges() {
      if (ranges.empty()) {
        mergedRanges.clear();
        isMerged = true;
        return;
      }

      std::sort(ranges.begin(), ranges.end());
      mergedRanges.clear();
      mergedRanges.push_back(ranges[0]);

      for (size_t i = 1; i < ranges.size(); ++i) {
        TimeRange& last = mergedRanges.back();
        const TimeRange& current = ranges[i];

        if (current.start <= last.end) {
          last.end = std::max(last.end, current.end);
        } else {
          mergedRanges.push_back(current);
        }
      }

      isMerged = true;
    } 
};


// *************************************************************************************
// ****************************** DEFINE MAPS FOR EVENTS *******************************
// *************************************************************************************

// Multimaps to hold events from anatrees
std::multimap<int64_t, std::tuple<double, double, double, double, double, int, int, int, InterruptedCoincidenceType>> gated_implants_map;
std::multimap<int64_t, std::tuple<double, double, double, double, double, int, int, int, InterruptedCoincidenceType>> all_implants_map;
std::multimap<int64_t, std::tuple<double, double, double, double, double, int, int, int>> good_decays_map;
std::multimap<int64_t, std::pair<std::pair<int, int>, std::vector<std::tuple<CorrelationType, BetaType, int64_t, int64_t, double, double, double, double, double, double> > > > matched_decays_map;
std::multimap<int64_t, std::tuple<double, int>> germanium_map;



// *************************************************************************************
// *************************** DEFINE DEADTIME WINDOW MANAGER **************************
// *************************************************************************************

// Define a TimeRangeManager object to manage the deadtime windows for all implants in DSSD

TimeRangeManager deadtimeWindowManager;

// *************************************************************************************
// ****************************** START MACRO ******************************************
// *************************************************************************************

void ionbeta(const char* input, const char* output){
  
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
  TTree* germanium_tree = (TTree*)file->Get("germanium_tree");
  
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
  TTreeReader germanium_reader(germanium_tree);

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
  TTreeReaderValue<Int_t> decay_bplast(decay_reader, "decay.bp"); // bp = 0 neither fired, bp = 1 only bp1 fired, bp = 2 only bp2 fired, bp = 3 both fired

  // Define leaves of variables for germanium tree
  TTreeReaderValue<ULong64_t> germanium_time(germanium_reader, "germanium.time");
  TTreeReaderValue<Double_t> germanium_energy(germanium_reader, "germanium.energy");
  TTreeReaderValue<Int_t> germanium_spill(germanium_reader, "germanium.sp"); // sp = 1 spill, sp = 2 no spill


  // *************************************************************************************
  // ****************************** DEFINE HISTOGRAMS ************************************
  // *************************************************************************************

  // Histograms for implant beta matches
  TH2F* h2_aida_implant_xy = new TH2F("aida_implant_xy", "AIDA Implant XY", 384, 0, 384, 128, 0, 128);
  TH2F* h2_aida_matched_xy = new TH2F("aida_matched_xy", "AIDA Matched XY", 384, 0, 384, 128, 0, 128);
  TH1F* h1_aida_wr_times = new TH1F("aida_wr_times", "AIDA WR Times", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);

  TH1F* h1_aida_implant_beta_dt = new TH1F("aida_implant_beta_dt", "Implant-Decay #Deltat;Implant-Decay #Deltat (); Counts/", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD/constants::TIME_SCALE, constants::TIME_THRESHOLD/constants::TIME_SCALE);

  // Histograms for beta candidate events
  TH1F* h1_implantbeta_candidate_multiplicity_forwards = new TH1F("implantbeta_candidate_multiplicity_forwards", "Implant-Decay Forwards Candidate Multiplicity; Candidate Multiplicity; Counts", 100, 0, 100);
  TH1F* h1_implantbeta_candidate_multiplicity_backwards = new TH1F("implantbeta_candidate_multiplicity_backwards", "Implant-Decay Backwards Candidate Multiplicity; Candidate Multiplicity; Counts", 100, 0, 100);

  // Histograms for gamma correlated events
  TH1F* h1_implantbetagamma_spectrum_before_ionbeta = new TH1F("implantbetagamma_spectrum_before_ionbeta", "Implant-Beta-Gamma Energy Spectrum (all); Energy (keV); Counts/keV", 2000, 0, 2000);
  TH1F* h1_implantbetagamma_spectrum_after_ionbeta = new TH1F("implantbetagamma_spectrum_after_ionbeta", "Implant-Beta-Gamma Energy Spectrum (ionbeta matched)); Energy (keV); Counts/keV", 2000, 0, 2000);

  // Correlate forward implant-decay dt with other observables
  TH2F* h2_implant_energy_dt_forward = new TH2F("implant_energy_dt_forward", "Implant Energy vs Implant-Decay dt", 7000/20, 0, 7000, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_decay_energy_dt_forward = new TH2F("decay_energy_dt_forward", "Decay Energy vs Implant-Decay dt (Forward)", 1500/20, 0, 1500, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_strip_dt_forward = new TH2F("strip_dt_forward", "Strip XY vs Implant-Decay dt (Forward)", 528, 0, 528, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);

  // Correlate backward implant-decay dt with other observables
  TH2F* h2_implant_energy_dt_backward = new TH2F("implant_energy_dt_backward", "Implant Energy vs Implant-Decay dt (Backward)", 7000/20, 0, 7000, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_decay_energy_dt_backward = new TH2F("decay_energy_dt_backward", "Decay Energy vs Implant-Decay dt (Backward)", 1500/20, 0, 1500, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_strip_dt_backward = new TH2F("strip_dt_backward", "Strip XY vs Implant-Decay dt (Backward)", 528, 0, 528, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);

  // NTuple for plotting in python
  TNtuple* nt_aida_implant_beta_dt = new TNtuple("nt_aida_implant_beta_dt", "Implant Decay dt", "dt");
  TNtuple* nt_all_candidate_ionbeta_dt = new TNtuple("nt_all_candidate_ionbeta_dt", "All beta candidate implant decay dt", "dt");

  // WR times of decays vetoed by implant deadtime 
  TH1F* h1_deadtime_implant_vetoed_decay_candidates_time = new TH1F("deadtime_implant_vetoed_decay_candidates_time", "AIDA WR times of vetoed decay candidates (Implant deadtime)", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1_interrupted_implant_vetoed_decay_candidates_time = new TH1F("interrupted_implant_vetoed_decay_candidates_time", "AIDA WR times of vetoed decay candidates (Implant interruption)", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);

  // Define plots doen with a canddate beta plot
  TH1F* h1_aida_implant_beta_dt_candidate_cut = new TH1F("aida_implant_beta_dt_candidate_cut", "Implant-Decay #Deltat (Beta Candidate Cut);Implant-Decay #deltat (); Counts/", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD/constants::TIME_SCALE, constants::TIME_THRESHOLD/constants::TIME_SCALE);
  TH2F* h2_implant_energy_dt_candidate_cut = new TH2F("implant_energy_dt_candidate_cut", "Implant Energy vs Implant-Decay dt (Beta Candidate Cut)", 7000/20, 0, 7000, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_decay_energy_dt_candidate_cut = new TH2F("decay_energy_dt_candidate_cut", "Decay Energy vs Implant-Decay dt (Beta Candidate Cut)", 1500/20, 0, 1500, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_strip_dt_candidate_cut = new TH2F("strip_dt_candidate_cut", "Strip XY vs Implant-Decay dt (Beta Candidate Cut)", 528, 0, 528, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);

  // Define plots for all candidate implant decay dt
  TH1F* h1_all_candidate_ionbeta_dt = new TH1F("all_candidate_ionbeta_dt", "Implant-Decay dt (All candidates)", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD);

  // Correlate low multiplicity implant-decay dt with other observables
  TH2F* h2_hitpattern_low_multiplicity = new TH2F("hitpattern_low_multiplicity", "AIDA Beta Candidate Hitpattern (Beta Multiplicity < 4)", 384, 0, 384, 128, 0, 128);
  TH2F* h2_implant_energy_dt_low_multiplicity = new TH2F("implant_energy_dt_low_multiplicity", "Implant Energy vs Implant-Decay dt (Beta Multiplicity < 4)", 7000/20, 0, 7000, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_decay_energy_dt_low_multiplicity = new TH2F("decay_energy_dt_low_multiplicity", "Decay Energy vs Implant-Decay dt (Beta Multiplicity < 4)", 1500/20, 0, 1500, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_strip_dt_low_multiplicity = new TH2F("strip_dt_low_multiplicity", "Strip XY vs Implant-Decay dt (Beta Multiplicity < 4)", 528, 0, 528, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH1F* h1_candidate_wr_time_low_multiplicity = new TH1F("candidate_wr_time_low_multiplicity", "AIDA Candidate WR Times (Beta Multiplicity < 4)", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1_candidate_ionbeta_dt_low_multiplicity = new TH1F("candidate_ionbeta_dt_low_multiplicity", "Candidate Implant-Decay dt (Beta Multiplicity < 4)", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD);

  // Correlate high multiplicity implant-decay dt with other observables
  TH2F* h2_hitpattern_high_multiplicity = new TH2F("hitpattern_high_multiplicity", "AIDA Beta Candidate Hitpattern (Beta Multiplicity > 20)", 384, 0, 384, 128, 0, 128);
  TH2F* h2_implant_energy_dt_high_multiplicity = new TH2F("implant_energy_dt_high_multiplicity", "Implant Energy vs Implant-Decay dt (Beta Multiplicity > 20)", 7000/20, 0, 7000, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_decay_energy_dt_high_multiplicity = new TH2F("decay_energy_dt_high_multiplicity", "Decay Energy vs Implant-Decay dt (Beta Multiplicity > 20)", 1500/20, 0, 1500, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_strip_dt_high_multiplicity = new TH2F("strip_dt_high_multiplicity", "Strip XY vs Implant-Decay dt (Beta Multiplicity > 20)", 528, 0, 528, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH1F* h1_candidate_wr_time_high_multiplicity = new TH1F("candidate_wr_time_high_multiplicity", "AIDA Candidate WR Times (Beta Multiplicity > 20)", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1_candidate_ionbeta_dt_high_multiplicity = new TH1F("candidate_ionbeta_dt_high_multiplicity", "Candidate Implant-Decay dt (Beta Multiplicity > 20)", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD);

  // *************************************************************************************
  // ****************************** FILL MAPS WITH EVENTS ********************************
  // *************************************************************************************

  // Loop over the implant, decay and germanium trees and fill the maps
  std::cout << "Start filled the maps!" << std::endl << std::endl; // Read gated implant events
  while (gatedimplant_reader.Next()){
    if( *gatedimplant_dssd==constants::DSSD && *gatedimplant_bplast==0 ){
      // Create new entry for map
      gated_implants_map.emplace(
        *gatedimplant_time,
        std::make_tuple(*gatedimplant_x, *gatedimplant_y, *gatedimplant_e, *gatedimplant_ex, *gatedimplant_ey, *gatedimplant_spill, *gatedimplant_bplast, *gatedimplant_dssd, NOT_INTERRUPTED)
      );
    }
  }
  std::cout << "Finished filling the gated implant map" << std::endl;
  std::cout << "Number of gated implant events: " << gated_implants_map.size() << std::endl << std::endl;

  // Read implant events
  while (implant_reader.Next()){
    if( *implant_dssd==constants::DSSD && *implant_bplast==0 ){
      // Create entry for map
      all_implants_map.emplace(
        *implant_time,
        std::make_tuple(*implant_x, *implant_y, *implant_e, *implant_ex, *implant_ey, *implant_spill, *implant_bplast, *implant_dssd, NOT_INTERRUPTED)
      );

      // Define deadtime range in manager
      deadtimeWindowManager.addRange((uint64_t)*implant_time, (uint64_t)*implant_time+constants::IMPLANT_DEAD_TIME);

    }
  }
  std::cout << "Finished filling the implant map" << std::endl;
  std::cout << "Number of All implant events cut on region: " << all_implants_map.size() << std::endl << std::endl;

  //************* DEBUG **************
  // std::cout << "[DEBUG] Printing deadtime ranges in the manager:\n";
  // for (const auto& r : deadtimeWindowManager.getMergedRanges()) {
  //   std::cout << "[" << r.start << ", " << r.end << "]\n";
  // }
  // for (auto gimp_evt = gated_implants_map.begin(); gimp_evt!=gated_implants_map.end(); gimp_evt++){
  //   uint64_t t = gimp_evt->first;
  //   if(deadtimeWindowManager.contains(t)){
  //     std::cout << "[DEBUG] Gated implant falls within deadtime!" << std::endl;
  //   }
  // }
  // std::exit(0);
  //************* DEBUG **************

  
  // Read decay events
  while (decay_reader.Next()){
    if( *decay_dssd==constants::DSSD && TMath::Abs( (int64_t)(*decay_time_x-*decay_time_y) )<5e3 && TMath::Abs(*decay_ex-*decay_ey)<168 && *decay_e>151 && *decay_e<1000 ){
      good_decays_map.emplace(
        *decay_time,
        std::make_tuple(*decay_x, *decay_y, *decay_e, *decay_ex, *decay_ey, *decay_spill, *decay_bplast, *decay_dssd)
      );
    }
  }
  std::cout << "Finished filling the decay map" << std::endl;
  std::cout << "Number of Decay events: " << good_decays_map.size() << std::endl << std::endl;

  // Read germanium events
  while (germanium_reader.Next()){
    germanium_map.emplace(
      *germanium_time,
      std::make_tuple(*germanium_energy, *germanium_spill)
    );
  }
  std::cout << "Finished filling the germanium map" << std::endl;
  std::cout << "Number of Germanium events: " << germanium_map.size() << std::endl << std::endl;



  // *************************************************************************************
  // ************ INTERRUPTED IMPLANT TAGGING ****** LOOP OVER GATED IMPLANTS ************
  // *************************************************************************************

  // Define varibles used in the loop
  int64_t last_gatedimplant_time = 0;
  int64_t current_gatedimplant_time = 0;
  int interrupted_implant_counter = 0;

  if (constants::VETO_INTERRUPTED_IMPLANTS){

    std::cout << "Started implant interruption tagging..." << std::endl;


    // Loop over all gated implants of interest
    for (auto gimp_evt = gated_implants_map.begin(); gimp_evt != gated_implants_map.end(); gimp_evt++){

      // Unpack event data
      last_gatedimplant_time = gimp_evt->first;
      auto [gimp_x, gimp_y, gimp_e, gimp_ex, gimp_ey, gimp_spill, gimp_bplast, gimp_dssd, gimp_interruption_type] = gimp_evt->second;

      for (auto imp_evt = gated_implants_map.begin(); imp_evt != gated_implants_map.end(); imp_evt++){

        // Unpack time of subsequent implant event
        current_gatedimplant_time = imp_evt->first;
        int64_t time_diff = current_gatedimplant_time - last_gatedimplant_time;
        auto [imp_x, imp_y, imp_e, imp_ex, imp_ey, imp_spill, imp_bplast, imp_dssd, imp_interruption_type] = imp_evt->second;

        // Check if current implant time is outside of time correlation window
        if ( time_diff <= 0 ){continue;}
        if (current_gatedimplant_time > last_gatedimplant_time+constants::TIME_THRESHOLD ){break;}

        // Check if implant is in positional window
        if ( TMath::Abs(gimp_x - imp_x) <= constants::POSITION_THRESHOLD && TMath::Abs(gimp_y - imp_y) <= constants::POSITION_THRESHOLD ){
          
          // Check if implant is within temporal window
          if ( time_diff < constants::TIME_THRESHOLD ){

            //************* DEBUG **************
            // std::cout << "Interrupted implant WR: " << last_gatedimplant_time << " Time Difference: " << time_diff*1e-9 << " ##### X,Y: " << TMath::Abs(gimp_x-imp_x) << "," << TMath::Abs(gimp_y-imp_y) << std::endl;
            //************* DEBUG **************

            // std::cout << "Found an interruption at " << time_diff/(double)constants::TIME_SCALE << std::endl;
            interrupted_implant_counter++;

            // Update interrupted implant tag of gated implant in the map
            imp_evt->second = std::make_tuple( imp_x, imp_y, imp_e, imp_ex, imp_ey, imp_spill, imp_bplast, imp_dssd, INTERRUPTED );
            
            break; // Once there is an interruption, we break from looping through implants 
          }

        }

      } // End of implant loop

    } // End of gated implant loop




    // Loop over all gated implants of interest
    for (auto gimp_evt = gated_implants_map.begin(); gimp_evt != gated_implants_map.end(); gimp_evt++){

      // Unpack event data
      last_gatedimplant_time = gimp_evt->first;
      auto [gimp_x, gimp_y, gimp_e, gimp_ex, gimp_ey, gimp_spill, gimp_bplast, gimp_dssd, gimp_interruption_type] = gimp_evt->second;

      for (auto imp_evt = all_implants_map.begin(); imp_evt != all_implants_map.end(); imp_evt++){

        // Unpack time of subsequent implant event
        current_gatedimplant_time = imp_evt->first;
        int64_t time_diff = current_gatedimplant_time - last_gatedimplant_time;
        auto [imp_x, imp_y, imp_e, imp_ex, imp_ey, imp_spill, imp_bplast, imp_dssd, imp_interruption_type] = imp_evt->second;

        // Check if current implant time is outside of time correlation window
        if ( time_diff <= 0 ){continue;}
        if (current_gatedimplant_time > last_gatedimplant_time+constants::TIME_THRESHOLD ){break;}

        // Check if implant is in positional window
        if ( TMath::Abs(gimp_x - imp_x) <= constants::POSITION_THRESHOLD && TMath::Abs(gimp_y - imp_y) <= constants::POSITION_THRESHOLD ){
          
          // Check if implant is within temporal window
          if ( time_diff < constants::TIME_THRESHOLD ){

            //************* DEBUG **************
            // std::cout << "Interrupted implant WR: " << last_gatedimplant_time << " Time Difference: " << time_diff*1e-9 << " ##### X,Y: " << TMath::Abs(gimp_x-imp_x) << "," << TMath::Abs(gimp_y-imp_y) << std::endl;
            //************* DEBUG **************

            // std::cout << "Found an interruption at " << time_diff/(double)constants::TIME_SCALE << std::endl;
            interrupted_implant_counter++;

            // Update interrupted implant tag of gated implant in the map
            gimp_evt->second = std::make_tuple( gimp_x, gimp_y, gimp_e, gimp_ex, gimp_ey, gimp_spill, gimp_bplast, gimp_dssd, INTERRUPTED );
            
            break; // Once there is an interruption, we break from looping through implants 
          }

        }

      } // End of implant loop

    } // End of gated implant loop




  std::cout << "Number of implant-decay correlations interrupted by subsequent implants: " << interrupted_implant_counter << std::endl << std::endl;

  }

  // *************************************************************************************
  // ****************** ION-BETA MATCHING ****** LOOP OVER GATED IMPLANTS ****************
  // *************************************************************************************

  std::cout << "Started Implant-Decay matching routine..." << std::endl;

  // Define & declare variebles to be used inside loop of gated implant events
  double gatedimplant_pos_x;
  double gatedimplant_pos_y;
  bool found_forward_candidate;
  bool found_backwards_candidate;
  int forwards_implantbeta_candidate_counter;
  int backwards_implantbeta_candidate_counter;
  int matched_implantdecays_counter = 0;
  int matched_backwards_implantdecays_counter = 0;
  int matched_implantbetagamma_counter = 0;
  std::vector<std::tuple<CorrelationType, BetaType, int64_t, int64_t, double, double, double, double, double, double>> beta_candidate_data_vector;

  // Loop over all gated implant events in map and perform a beta match
  for (auto gimp_evt = gated_implants_map.begin(); gimp_evt != gated_implants_map.end(); gimp_evt++){
    
    // Empty vector containing previous beta candidate
    beta_candidate_data_vector.clear();

    // Unpack event variables for current gated implant event
    auto [x, y, e, ex, ey, spill, bplast, dssd, interruption_type] = gimp_evt->second;

    // Check noisy implant channel strips and skip
    if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_IMPLANT, x) ){ continue; }
    if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_IMPLANT, y) ){ continue; }

    // Check if gated implant is interrupted
    if (interruption_type == INTERRUPTED ){
      h1_interrupted_implant_vetoed_decay_candidates_time->Fill(last_gatedimplant_time);
      continue; 
    }


    // ************* DEBUG **************
    // if ( e < 2000 ) { continue; } // Skip any possible throughgoing ion
    // ************* DEBUG **************
    
    // Reset counter to check all beta candidates
    int forwards_implantbeta_candidate_counter = 0;     
    int backwards_implantbeta_candidate_counter = 0;     

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
      auto [decay_x, decay_y, decay_e, decay_ex, decay_ey, decay_spill, decay_bplast, decay_dssd] = decay_evt->second;

      // Check for noisy decay branch strips and skip
      if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_DECAY, decay_x) ){ continue; }
      if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_DECAY, decay_y) ){ continue; }

      // Check if decay event is onspill and skip if defined by user and is from desired dssd
      if( constants::ONLY_OFFSPILL_DECAY && decay_spill == 1){ continue; }
      
      // Check if decay event is within position threshold
      if ( TMath::Abs(decay_x - gatedimplant_pos_x) <= constants::POSITION_THRESHOLD && TMath::Abs(decay_y - gatedimplant_pos_y) <= constants::POSITION_THRESHOLD ){

        // Find time difference between implant event and decay event
        int64_t time_diff = decay_evt->first - last_gatedimplant_time;

        // *************************************************************************************
        // ****************************** FOUND FORWARD BETA CANDIDATE *************************
        // *************************************************************************************

        // Check if decay event falls within time threshold and decay event occures after implant event 
        if (time_diff > 0 && time_diff < constants::TIME_THRESHOLD) {

          // Check decay is outside implant promptflash
          // if ( time_diff < 60e3 ){ continue; }

          // Check if decay occures within a deadtime induced my implant in AIDA
          if ( deadtimeWindowManager.contains(decay_evt->first) ){
            h1_deadtime_implant_vetoed_decay_candidates_time->Fill(decay_evt->first);
            std::cout << "Rejected a forwards ion-beta candidate due to implant deadtime." << std::endl;
            continue; 
          }
          
          forwards_implantbeta_candidate_counter++; // Increase counter for beta canditade event

          BetaType beta_type = CANDIDATE; // Set tag for beta
          CorrelationType correlation_type = FORWARDS; // Set tag for correlation type

          // Check that there has not been a forward beta candidate that has been matched to implant event
          if (!found_forward_candidate){
            
            matched_implantdecays_counter++; // Increase counter for succesfull matched implant decay
            found_forward_candidate = true; // Change flag for succesfull forward implant decay match
            beta_type = MATCH; // Overwrite beta_type
          }


          // Add beta data to the vector
          beta_candidate_data_vector.emplace_back( std::make_tuple(correlation_type, beta_type, time_diff, decay_evt->first, gatedimplant_pos_x, gatedimplant_pos_y, e, decay_x, decay_y, decay_e) );

        }

        // *************************************************************************************
        // ****************************** FOUND BACKWARD BETA CANDIDATE ************************
        // *************************************************************************************

        // Check if decay event falls within time threshold and decay event occures before implant event 
        if (time_diff < 0 && -time_diff < constants::TIME_THRESHOLD) {

          // Check if decay occures within a deadtime induced my implant in AIDA
          if ( deadtimeWindowManager.contains(decay_evt->first) ){
            std::cout << "Rejected backwards ion-beta candidate due to implant deadtime." << std::endl;
            continue; 
          }

          backwards_implantbeta_candidate_counter++; // Increase counter for beta canditade event

          BetaType beta_type = CANDIDATE; // Set tag for beta
          CorrelationType correlation_type = BACKWARDS; // Set tag for correlation type

          // Check that there has not been a backwards beta candidate that has been matched to implant event
          if (!found_backwards_candidate){

            matched_backwards_implantdecays_counter++; // Increase counter for succesfull matched implant decay
            found_backwards_candidate = true; // Change flag for succesfull backward implant decay match
            beta_type = MATCH; // Overwrite beta_type

          }

          // Add beta data to the vector
          beta_candidate_data_vector.emplace_back( std::make_tuple(correlation_type, beta_type, time_diff, decay_evt->first, gatedimplant_pos_x, gatedimplant_pos_y, e, decay_x, decay_y, decay_e) );

        }

      }

    } // End of decay event loop

    // Fill matched candidate map
    matched_decays_map.emplace(
      last_gatedimplant_time, 
      std::make_pair(
        std::make_pair(
          forwards_implantbeta_candidate_counter,
          backwards_implantbeta_candidate_counter
        ), 
        beta_candidate_data_vector
      ) 
    );


  } // End of gated implant event loop

  std::cout << "Finalised Implant-Decay matching routine!" << std::endl << std::endl;
  


  // *************************************************************************************
  // ************************** MATCHED IMPLANT - DECAY PLOTTING *************************
  // *************************************************************************************
  

  std::cout << "Started filling histograms..." << std::endl;
  
  // Loop over each gated implant match 
  for ( auto matched_decay_evt = matched_decays_map.begin(); matched_decay_evt!=matched_decays_map.end(); matched_decay_evt++ ){

    // Unpack the implant time of the gated implant
    int64_t gimp_time = matched_decay_evt->first;
    auto [forwards_candidate_multiplicity, backwards_candidate_multiplicity] = matched_decay_evt->second.first;

    // Plot candidate multiplicity
    if(forwards_candidate_multiplicity>0) {h1_implantbeta_candidate_multiplicity_forwards->Fill(forwards_candidate_multiplicity);}
    if(backwards_candidate_multiplicity>0) {h1_implantbeta_candidate_multiplicity_backwards->Fill(backwards_candidate_multiplicity);}

    // Loop over all matched beta candidate for each gated implant
    for ( auto beta_candidate_data : matched_decay_evt->second.second ){

      // Unpack the data for each beta candidate
      auto [correlation_type, beta_type, time_diff, decay_time, gimp_x, gimp_y, gimp_e, decay_x, decay_y, decay_e] = beta_candidate_data;

      if ( correlation_type==FORWARDS && beta_type==MATCH ){
        h1_aida_implant_beta_dt->Fill(time_diff/constants::TIME_SCALE);
        nt_aida_implant_beta_dt->Fill((double)time_diff);
        h2_aida_matched_xy->Fill(decay_x, decay_y); 
        h1_aida_wr_times->Fill(decay_time);
      }

      if ( correlation_type==BACKWARDS && beta_type==MATCH ){
        h1_aida_implant_beta_dt->Fill(time_diff/constants::TIME_SCALE);
        nt_aida_implant_beta_dt->Fill((double)time_diff);
      }

      if ( correlation_type==FORWARDS ){
        h2_implant_energy_dt_forward->Fill(gimp_e, time_diff);
        h2_decay_energy_dt_forward->Fill(decay_e, time_diff);
        h2_strip_dt_forward->Fill(decay_x, time_diff);
        h2_strip_dt_forward->Fill(decay_y+400, time_diff);

        // All candidate dionbeta dt
        nt_all_candidate_ionbeta_dt->Fill((double)time_diff);
        h1_all_candidate_ionbeta_dt->Fill(time_diff);
      }

      if ( correlation_type==BACKWARDS ){
        h2_implant_energy_dt_backward->Fill(gimp_e, -time_diff);
        h2_decay_energy_dt_backward->Fill(decay_e, -time_diff);
        h2_strip_dt_backward->Fill(decay_x, -time_diff);
        h2_strip_dt_backward->Fill(decay_y+400, -time_diff);

        // All candidate dionbeta dt
        nt_all_candidate_ionbeta_dt->Fill((double)time_diff);
        h1_all_candidate_ionbeta_dt->Fill(time_diff);
      }

      if ( correlation_type==FORWARDS && forwards_candidate_multiplicity>0 && forwards_candidate_multiplicity<=constants::BETA_CANDIDATE_CUT){
        h1_aida_implant_beta_dt_candidate_cut->Fill(time_diff/constants::TIME_SCALE);
        h2_implant_energy_dt_candidate_cut->Fill(gimp_e, time_diff);
        h2_decay_energy_dt_candidate_cut->Fill(decay_e, time_diff);
        h2_strip_dt_candidate_cut->Fill(decay_x, time_diff);
        h2_strip_dt_candidate_cut->Fill(decay_y+400, time_diff);
      }
      
      if ( correlation_type==BACKWARDS && backwards_candidate_multiplicity>0 && backwards_candidate_multiplicity<=constants::BETA_CANDIDATE_CUT){
        h1_aida_implant_beta_dt_candidate_cut->Fill(time_diff/constants::TIME_SCALE);
      }


      if ( correlation_type == FORWARDS && forwards_candidate_multiplicity >= 20){
        std::cout << "[DEBUG] Found high multiplicity forward beta candidate:  Implant WR Time: " << gimp_time << "  #####  Decay WR Time: " << decay_time << "  #####  Implant Energy: " << gimp_e << " MeV  #####  Decay Energy: " << decay_e << " keV  #####  Decay Pos (X, Y): (" << decay_x << "," << decay_y << ")" << std::endl;
        h1_candidate_wr_time_high_multiplicity->Fill(decay_time);
        h1_candidate_ionbeta_dt_high_multiplicity->Fill(time_diff);
        h2_hitpattern_high_multiplicity->Fill(decay_x, decay_y);
        h2_implant_energy_dt_high_multiplicity->Fill(gimp_e, time_diff);
        h2_decay_energy_dt_high_multiplicity->Fill(decay_e, time_diff);
        h2_strip_dt_high_multiplicity->Fill(decay_x, time_diff);
        h2_strip_dt_high_multiplicity->Fill(decay_y+400, time_diff);
      }

      // if ( correlation_type == BACKWARDS && backwards_candidate_multiplicity >= 20){
        
      // }

      if ( correlation_type == FORWARDS && forwards_candidate_multiplicity <= 2 && backwards_candidate_multiplicity > 0){
        h1_candidate_wr_time_low_multiplicity->Fill(decay_time);
        h1_candidate_ionbeta_dt_low_multiplicity->Fill(time_diff);
        h2_hitpattern_low_multiplicity->Fill(decay_x, decay_y);
        h2_implant_energy_dt_low_multiplicity->Fill(gimp_e, time_diff);
        h2_decay_energy_dt_low_multiplicity->Fill(decay_e, time_diff);
        h2_strip_dt_low_multiplicity->Fill(decay_x, time_diff);
        h2_strip_dt_low_multiplicity->Fill(decay_y+400, time_diff);
      }

      // if ( correlation_type == BACKWARDS && backwards_candidate_multiplicity <= 2 && backwards_candidate_multiplicity > 0){
        // 
      // }

    }
    
  }
  
  std::cout << "Finished filling histograms!" << std::endl << std::endl;
 
  // *************************************************************************************
  // ************************** MATCHED DECAY - GAMMA CORRELATION ************************
  // *************************************************************************************

  std::cout << "Started Matched Beta - Gamma matching routine..." << std::endl;
  
  // Loop over all matched decay events
  for( auto matched_decay_evt = matched_decays_map.begin(); matched_decay_evt != matched_decays_map.end(); matched_decay_evt++ ){

    for ( auto matched_beta_data : matched_decay_evt->second.second){

      // Unpack matched decay event variables
      auto [correlation_type, beta_type, useless_1, last_matched_decay_time, useless_2, useless_3, useless_4, useless_5, useless_6, useless_7] = matched_beta_data;

      // Skip any matche dimplant that is correlated backwards
      if ( correlation_type!=FORWARDS ){ continue; }

      // Find the germanium event starting at the same time as the decay event (50 microsecond grace period)
      auto germanium_start = germanium_map.lower_bound(matched_decay_evt->first - 50e3);

      // *************************************************************************************
      // **************************  LOOP OVER GERMANIUM EVENTS ******************************
      // *************************************************************************************

      // Loop over all germanium events between decay event and end of prompt gamma window
      for( auto germanium_evt = germanium_start; germanium_evt != germanium_map.end(); germanium_evt++ ){

        // int time_diff = last_matched_decay_time - germanium_evt->first; // Reverse like in jeroens code
        int time_diff = germanium_evt->first - last_matched_decay_time; // Like regular people

        if ( time_diff > constants::PROMPT_WINDOW_END ){ break; }

        // Check if germanium event is within prompt window
        if ( time_diff > constants::PROMPT_WINDOW_START && time_diff < constants::PROMPT_WINDOW_END ){
            
          // Unpack germanium event items within prompt window
          auto [germanium_energy, germanium_spill] = germanium_evt->second; 

          // Fill gamma energy spectrum and increase counter
          h1_implantbetagamma_spectrum_after_ionbeta->Fill(germanium_energy); // Fill gamma energy spectrum
          matched_implantbetagamma_counter++;
        }

      } // End of germanium event loop

    }
      
  } // End of matched decay event loop

  std::cout << "Finished Matched Beta - Gamma matching routine!" << std::endl << std::endl;
    
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

  if (constants::CHECK_BETA_CANDITATES){
    TDirectory* candidateDir = outputFile->mkdir("candidate_matches");
    candidateDir->cd();
    h1_all_candidate_ionbeta_dt->Write();
    h1_implantbeta_candidate_multiplicity_forwards->Write(); 
    h1_implantbeta_candidate_multiplicity_backwards->Write(); 
    gFile->cd();
  }

  TDirectory* forwardDir = outputFile->mkdir("forward_matches");
  forwardDir->cd();
  h2_implant_energy_dt_forward->Write();
  h2_decay_energy_dt_forward->Write();
  h2_strip_dt_forward->Write();
  h2_strip_dt_forward->Write();
  gFile->cd();

  TDirectory* backwardDir = outputFile->mkdir("backward_matches");
  backwardDir->cd();
  h2_implant_energy_dt_backward->Write();
  h2_decay_energy_dt_backward->Write();
  h2_strip_dt_backward->Write();
  h2_strip_dt_backward->Write();
  gFile->cd();

  TDirectory* candidateCutDir = outputFile->mkdir("candidate_cut_matches");
  candidateCutDir->cd();
  h1_aida_implant_beta_dt_candidate_cut->Write();
  h2_implant_energy_dt_candidate_cut->Write();
  h2_decay_energy_dt_candidate_cut->Write();
  h2_strip_dt_candidate_cut->Write();
  h2_strip_dt_candidate_cut->Write();
  gFile->cd();

  TDirectory* lowCandidateMultiplicity = outputFile->mkdir("low_candidate_multiplicity");
  lowCandidateMultiplicity->cd();
  h1_candidate_wr_time_low_multiplicity->Write();
  h1_candidate_ionbeta_dt_low_multiplicity->Write();
  h2_hitpattern_low_multiplicity->Write();
  h2_implant_energy_dt_low_multiplicity->Write();
  h2_decay_energy_dt_low_multiplicity->Write();
  h2_strip_dt_low_multiplicity->Write();
  gFile->cd();

  TDirectory* highCandidateMultiplicity = outputFile->mkdir("high_candidate_multiplicity");
  highCandidateMultiplicity->cd();
  h1_candidate_wr_time_high_multiplicity->Write();
  h1_candidate_ionbeta_dt_high_multiplicity->Write();
  h2_hitpattern_high_multiplicity->Write();
  h2_implant_energy_dt_high_multiplicity->Write();
  h2_decay_energy_dt_high_multiplicity->Write();
  h2_strip_dt_high_multiplicity->Write();
  gFile->cd();
  
  nt_aida_implant_beta_dt->Write();
  nt_all_candidate_ionbeta_dt->Write();

  h2_aida_implant_xy->Write();
  h1_aida_wr_times->Write();
  h2_aida_matched_xy->Write();
  h1_aida_implant_beta_dt->Write();
  h1_implantbetagamma_spectrum_after_ionbeta->Write();

  h1_deadtime_implant_vetoed_decay_candidates_time->Write();
  h1_interrupted_implant_vetoed_decay_candidates_time->Write();

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
