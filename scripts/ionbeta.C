#include<iostream>
#include<map>
#include<set>
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
  const std::string ISOTOPE_TREE = "84nb"; // Name suffix for gatedimplant tree & branch in anatree
  const int DSSD = 1; // Which DSSD will the analysis be run on

  const bool ONLY_OFFSPILL_DECAY = false; // Check for onspill decay matches
  const bool CHECK_BETA_CANDITATES = true; // Check for all beta candidates of an implant
  const bool VETO_INTERRUPTED_IMPLANTS = false; // Veto any implant which has a subsequent implant occure within time and position window
  /*const bool INCLUDE_BACKWARDS_MATCH = true; // Look for reverse time implant beta correlations*/
  const bool ALL_IMPLANTS_DECAY_MATCH = false; // Look at impdecay matches for all implants aswell

  const int64_t TIME_SCALE = 1e9; // Timescale of time variables wrt ns
  const int64_t TIME_THRESHOLD = 50 * TIME_SCALE; // Time threshold for implant beta correlation
  const double TIME_PER_DT_BIN = 1e9; // Time per bin in Implant-beta time correlation
  const int64_t POSITION_THRESHOLD = 1; //  Position window for decay wrt implant pixel as centroid

  const uint64_t IMPLANT_DEAD_TIME = 350e3; // The deadtime we impose on aida LEC after an implant occures in AIDA
  const double IMPLANT_DEADTIME_LOCAL_RANGE = 5; // How far from the implant is the deadtime applied in x and y
  const int BETA_CANDIDATE_CUT = 5; // Define number of candidate betas a implant must have before plotting
  const int BETA_GAMMA_CANDIDATE_CUT = 5;

  // const std::pair<int64_t, int64_t> PROMPT_GAMMA_WINDOW = { 13229, 17991 };
  const std::pair<int64_t, int64_t> PROMPT_GAMMA_WINDOW = { -18200, -12900 };
  // const std::pair<int64_t, int64_t> PROMPT_GAMMA_WINDOW = { -0, 5000 };

  const int NEIGHBOURING_POSITION_BINS = POSITION_THRESHOLD*2+1; // Bin # used for beta candidate hit pattern histogram
  const int64_t IMPDECAY_TIME_BINS = 2*TIME_THRESHOLD/TIME_PER_DT_BIN;

  const std::vector<double> BROKEN_AIDA_X_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_IMPLANT = {};
  const std::vector<double> BROKEN_AIDA_X_STRIPS_DECAY = {63, 63.5, 64, 66, 130, 189, 191, 194, 225, 256, 319, 320, 192, 128, 255}; 
  const std::vector<double> BROKEN_AIDA_Y_STRIPS_DECAY = {127};

  // Map containing known beta delayed gammas for matched beta 
  const std::map< std::string, std::set<double> > KNOWN_BETA_DELAYED_GAMMAS = {
    { "84nb", {456.2, 540., 579.4, 704., 1036.4, 1119.6, 1426.7} },
    { "70br", {944.6, 656., 1096.3, 782., 963.7, 911.7, 1033.6, 596.0, 690.2, 958., 1604., 348., 1168.8, 1656.2, 1062.0} },
    { "81zr", {113.31, 175.38, 230.1} }
  };

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
enum BetaType { CANDIDATE, MATCH, SECOND_MATCH, THIRD_MATCH, TENTH_MATCH}; // Tags for candidate betas (Have they been matched to an implant or not)
enum InterruptedCoincidenceType { NOT_INTERRUPTED, INTERRUPTED }; // Tags to define if an implant has beenn interrupted in a time window


// *************************************************************************************
// ****************************** DEFINE CUSTOM TYPES **********************************
// *************************************************************************************

typedef std::tuple<CorrelationType, BetaType, int64_t, int64_t, double, double, double, double, double, double, double > BetaCandidateInfoTuple_t;
typedef std::tuple<double> BetaDelatedGammaInfoTuple_t;

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


class TimeRangeManagerGlobal {
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

class TimeRangeManagerLocal {
  public:

    struct TimeRange {
      uint64_t  start;
      uint64_t  end;
      double    posX;
      double    posY;

      bool contains(uint64_t time, double x, double y, double posRange) const {
        return time >= start && time <= end && TMath::Abs(x-posX) <= posRange && TMath::Abs(y-posY) <= posRange;
      }

      bool operator<(const TimeRange& other) const {
        return start < other.start;
      }
    };

    // Constructor
    TimeRangeManagerLocal(double positionThreshold=5){
      posRange = positionThreshold;
    }

    // Add a new time range
    void addRange(uint64_t start, uint64_t end, double posX, double posY) {
      if (start > end) std::swap(start, end);  // ensure proper order
      ranges.push_back({start, end, posX, posY});
    }

    // Check if a time is within any merged range
    bool contains(uint64_t time, double x, double y) {
      for (const auto& range : ranges) {
        if (range.contains(time, x, y, posRange)) return true;
      }
      return false;
    }

    // Access the ranges
    const std::vector<TimeRange>& getRanges() {
      return ranges;
    }

  private:
    double                  posRange;
    std::vector<TimeRange>  ranges;

};

struct BetaCandidateInfo{
  // Struct to hold the information of each beta candidate for an implant decay match

  CorrelationType   correlationType;
  BetaType          betaType;
  int64_t           timeDiff;
  int64_t           decayTime;
  double            gimpPosX;
  double            gimpPosY;
  double            gimpE;
  double            decayPosX;
  double            decayPosY;
  double            decayE;
  int               decaySpill;

  bool              delayedGammaMatch   = false;
  double            delayedGammaEnergy;

  void setBetaParameters( CorrelationType corrType, BetaType bType, int64_t tDiff, int64_t decT, double impX, double impY, double impE, double decX, double decY, double decE, int decSpill){
    // Struct beta param setter
    correlationType = corrType;
    betaType        = bType;
    timeDiff        = tDiff;
    decayTime       = decT;
    gimpPosX        = impX;
    gimpPosY        = impY;
    gimpE           = impE;
    decayPosX       = decX;
    decayPosY       = decY;
    decayE          = decE;
  }

  void setGammaParameters( double gammaE ){
    // Struct gamma param setter
    delayedGammaEnergy  = gammaE;
    delayedGammaMatch   = true; // Flip matched bool
  }

  BetaCandidateInfoTuple_t getBetaParameters(){
    // Struct beta param getter
    return std::make_tuple(correlationType, betaType, timeDiff, decayTime, gimpPosX, gimpPosY, gimpE, decayPosX, decayPosY, decayE, decaySpill);
  } 

  BetaDelatedGammaInfoTuple_t getGammaParameters(){
    // Sruct gamma param getter
    return std::make_tuple(decayE);
  }

  bool hasMatchedGammaInRange(std::set<double> knownGammaEnergySet, double gammaCheckRange = 5. ){
    // Check if matched beta delayed gamma is within range of known gamma
    if (!delayedGammaMatch) { return false; }

    // Check if matched gamma is contained in known gamma energy set 
    return knownGammaEnergySet.lower_bound(delayedGammaEnergy-gammaCheckRange) != knownGammaEnergySet.upper_bound(delayedGammaEnergy+gammaCheckRange);
  }
  

};

struct ImplantDecayInfo{
  // Structure to hold information of implant decay matches per implant

  std::pair<int, int> implantDecayMatchMultiplicities;
  std::vector<BetaCandidateInfo> betaCandidateInfoVector;
};


// *************************************************************************************
// ****************************** DEFINE MAPS FOR EVENTS *******************************
// *************************************************************************************

// Multimaps to hold events from anatrees
std::multimap<int64_t, std::tuple<double, double, double, double, double, int, int, int, InterruptedCoincidenceType>> gated_implants_map;
std::multimap<int64_t, std::tuple<double, double, double, double, double, int, int, int, InterruptedCoincidenceType>> all_implants_map;
std::multimap<int64_t, std::tuple<double, double, double, double, double, int, int, int>> good_decays_map;
// std::multimap<int64_t, std::pair<std::pair<int, int>, std::vector<std::tuple<CorrelationType, BetaType, int64_t, int64_t, double, double, double, double, double, double> > > > matched_decays_map;
std::multimap<int64_t, ImplantDecayInfo> gated_matched_decays_map;
std::multimap<int64_t, ImplantDecayInfo> all_matched_decays_map;
std::multimap<int64_t, std::tuple<double, int>> germanium_map;



// *************************************************************************************
// *************************** DEFINE DEADTIME WINDOW MANAGER **************************
// *************************************************************************************

// Define a TimeRangeManager object to manage the deadtime windows for all implants in DSSD

TimeRangeManagerLocal deadtimeWindowManager(constants::IMPLANT_DEADTIME_LOCAL_RANGE);

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
  TH1F* h1_aida_implant_beta_dt = new TH1F("aida_implant_beta_firstmatch_dt", "Implant-Decay #Deltat (First Match);Implant-Decay #Deltat (); Counts/", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD/constants::TIME_SCALE, constants::TIME_THRESHOLD/constants::TIME_SCALE);
  TH1F* h1_aida_implant_beta_secondmatch_dt = new TH1F("aida_implant_beta_secondmatch_dt", "Implant-Decay #Deltat (Second Match);Implant-Decay #Deltat (); Counts/", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD/constants::TIME_SCALE, constants::TIME_THRESHOLD/constants::TIME_SCALE);
  TH1F* h1_aida_implant_beta_thirdmatch_dt = new TH1F("aida_implant_beta_thirdmatch_dt", "Implant-Decay #Deltat (Third Match);Implant-Decay #Deltat (); Counts/", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD/constants::TIME_SCALE, constants::TIME_THRESHOLD/constants::TIME_SCALE);
  TH1F* h1_aida_implant_beta_tenthmatch_dt = new TH1F("aida_implant_beta_tenthmatch_dt", "Implant-Decay #Deltat (Tenth Match);Implant-Decay #Deltat (); Counts/", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD/constants::TIME_SCALE, constants::TIME_THRESHOLD/constants::TIME_SCALE);

  // Histograms for beta candidate events
  TH1F* h1_implantbeta_candidate_multiplicity_forwards = new TH1F("implantbeta_candidate_multiplicity_forwards", "Implant-Decay Forwards Candidate Multiplicity; Candidate Multiplicity; Counts", 100, 0, 100);
  TH1F* h1_implantbeta_candidate_multiplicity_backwards = new TH1F("implantbeta_candidate_multiplicity_backwards", "Implant-Decay Backwards Candidate Multiplicity; Candidate Multiplicity; Counts", 100, 0, 100);

  // HIST^O
  // Histograms for betaGamma correlations
  TH1F* h1_gamma_spectrum = new TH1F("gamma_spectrum", "Gamma Energy Spectrum ; Energy (keV); Counts/keV", 2000, 0, 2000);
  TH1F* h1_implantbetagamma_spectrum_before_ionbeta = new TH1F("implantbetagamma_spectrum_before_ionbeta", "Beta-Gamma Energy Spectrum (all); Energy (keV); Counts/keV", 2000, 0, 2000);
  TH1F* h1_implantbetagamma_spectrum_before_ionbeta_dt = new TH1F("implantbetagamma_spectrum_dt", "Decay Germanium dt; Time (ns); Counts", 2000, -20000, 2000);
  TH2F* h2_implantbetagamma_spectrum_before_ionbeta_dt_energy = new TH2F("h2_implantbetagamma_spectrum_before_ionbeta_dt_energy", "Time between beta decay and gamma ray vs energy (keV); Time (ns); Energy (keV); Counts", 500, -20e3, 2e3,2000,0,2000);

  // Histograms for implant-beta-gamma backward correlated events
  TH1F* h1_gatedimplantbetagamma_spectrum_after_ionbeta_backwardmatch = new TH1F("gatedimplantbetagamma_spectrum_after_ionbeta_backwardmatch", "Implant-Beta-Gamma Energy Spectrum (backward ionbeta matched)); Energy (keV); Counts/keV", 2000, 0, 2000);
  TH2F* h2_gatedimplantbetagamma_energy_vs_decaystrip_backwardmatch = new TH2F("gatedimplantbetagamma_energy_vs_decaystrip_backwardmatch", "Matched Gamma Energy vs Decay Strip (X then Y);Strip;Energy (keV)", 528, 0, 528, 2000, 0, 2000);
  TH1F* h1_gatedimplantbetagamma_spectrum_after_ionbeta_dt_backwardmatch = new TH1F("gatedimplantbetagamma_spectrum_after_ionbeta_dt_backwardmatch", "Time between beta decay and gamma ray backward match; Time (ns); Counts", 2000, -20000, 2000);
  TH2F* h2_gatedimplantbetagamma_spectrum_after_ionbeta_dt_energy_backwardmatch = new TH2F("gatedimplantbetagamma_spectrum_after_ionbeta_dt_energy_backwardmatch", "Time between beta decay and gamma ray backward Match vs energy (keV); Time (ns); Energy (keV); Counts", 500, -20e3, 2e3,2000,0,2000);
  TH2F* h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_backwardmatch = new TH2F("gatedimplantbetagamma_spectrum_after_ionbeta_square_backwardmatch", "gatedimplantbetagamma_spectrum_after_ionbeta_square_backwardmatch", 2000,0,2000,2000,0,2000);
  TH2F* h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_time_backwardmatch = new TH2F("gatedimplantbetagamma_spectrum_after_ionbeta_square_time_backwardmatch", "gatedimplantbetagamma_spectrum_after_ionbeta_square_time_backwardmatch", 2000,0,2000,200,-100,100);
  // Histograms for implant-beta-gamma firward correlated events
  TH1F* h1_gatedimplantbetagamma_spectrum_after_ionbeta_forwardmatch = new TH1F("gatedimplantbetagamma_spectrum_after_ionbeta_forwardmatch", "Implant-Beta-Gamma Energy Spectrum (forward ionbeta matched)); Energy (keV); Counts/keV", 2000, 0, 2000);
  TH2F* h2_gatedimplantbetagamma_energy_vs_decaystrip_forwardmatch = new TH2F("gatedimplantbetagamma_energy_vs_decaystrip_forwardmatch", "Matched Gamma Energy vs Decay Strip (X then Y);Strip;Energy (keV)", 528, 0, 528, 2000, 0, 2000);
  TH1F* h1_gatedimplantbetagamma_spectrum_after_ionbeta_dt_forwardmatch = new TH1F("gatedimplantbetagamma_spectrum_after_ionbeta_dt_forwardmatch", "Time between beta decay and gamma ray forward match; Time (ns); Counts", 2000, -20000, 2000);
  TH2F* h2_gatedimplantbetagamma_spectrum_after_ionbeta_dt_energy_forwardmatch = new TH2F("gatedimplantbetagamma_spectrum_after_ionbeta_dt_energy_forwardmatch", "Time between beta decay and gamma ray Forward Match vs energy (keV); Time (ns); Energy (keV); Counts", 500, -20e3, 2e3,2000,0,2000);
  TH2F* h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_forwardmatch = new TH2F("gatedimplantbetagamma_spectrum_after_ionbeta_square_forwardmatch", "gatedimplantbetagamma_spectrum_after_ionbeta_square_forwardmatch", 2000,0,2000,2000,0,2000);
  TH2F* h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_time_forwardmatch = new TH2F("gatedimplantbetagamma_spectrum_after_ionbeta_square_time_forwardmatch", "gatedimplantbetagamma_spectrum_after_ionbeta_square_time_forwardmatch", 2000,0,2000,200,-100,100);

  // Histograms for implant-beta-gamma backward correlated events
  TH1F* h1_allimplantbetagamma_spectrum_after_ionbeta_backwardmatch = new TH1F("allimplantbetagamma_spectrum_after_ionbeta_backwardmatch", "Implant-Beta-Gamma Energy Spectrum (backward ionbeta matched)); Energy (keV); Counts/keV", 2000, 0, 2000);
  TH2F* h2_allimplantbetagamma_energy_vs_decaystrip_backwardmatch = new TH2F("allimplantbetagamma_energy_vs_decaystrip_backwardmatch", "Matched Gamma Energy vs Decay Strip (X then Y);Strip;Energy (keV)", 528, 0, 528, 2000, 0, 2000);
  TH1F* h1_allimplantbetagamma_spectrum_after_ionbeta_dt_backwardmatch = new TH1F("allimplantbetagamma_spectrum_after_ionbeta_dt_backwardmatch", "Time between beta decay and gamma ray backward match; Time (ns); Counts", 2000, -20000, 2000);
  TH2F* h2_allimplantbetagamma_spectrum_after_ionbeta_dt_energy_backwardmatch = new TH2F("allimplantbetagamma_spectrum_after_ionbeta_dt_energy_backwardmatch", "Time between beta decay and gamma ray backward Match vs energy (keV); Time (ns); Energy (keV); Counts", 500, -20e3, 2e3,2000,0,2000);
  TH2F* h2_allimplantbetagamma_spectrum_after_ionbeta_square_backwardmatch = new TH2F("allimplantbetagamma_spectrum_after_ionbeta_square_backwardmatch", "allimplantbetagamma_spectrum_after_ionbeta_square_backwardmatch", 2000,0,2000,2000,0,2000);
  TH2F* h2_allimplantbetagamma_spectrum_after_ionbeta_square_time_backwardmatch = new TH2F("allimplantbetagamma_spectrum_after_ionbeta_square_time_backwardmatch", "allimplantbetagamma_spectrum_after_ionbeta_square_time_backwardmatch", 2000,0,2000,200,-100,100);
  // Histograms for implant-beta-gamma firward correlated events
  TH1F* h1_allimplantbetagamma_spectrum_after_ionbeta_forwardmatch = new TH1F("allimplantbetagamma_spectrum_after_ionbeta_forwardmatch", "Implant-Beta-Gamma Energy Spectrum (forward ionbeta matched)); Energy (keV); Counts/keV", 2000, 0, 2000);
  TH2F* h2_allimplantbetagamma_energy_vs_decaystrip_forwardmatch = new TH2F("allimplantbetagamma_energy_vs_decaystrip_forwardmatch", "Matched Gamma Energy vs Decay Strip (X then Y);Strip;Energy (keV)", 528, 0, 528, 2000, 0, 2000);
  TH1F* h1_allimplantbetagamma_spectrum_after_ionbeta_dt_forwardmatch = new TH1F("allimplantbetagamma_spectrum_after_ionbeta_dt_forwardmatch", "Time between beta decay and gamma ray forward match; Time (ns); Counts", 2000, -20000, 2000);
  TH2F* h2_allimplantbetagamma_spectrum_after_ionbeta_dt_energy_forwardmatch = new TH2F("allimplantbetagamma_spectrum_after_ionbeta_dt_energy_forwardmatch", "Time between beta decay and gamma ray Forward Match vs energy (keV); Time (ns); Energy (keV); Counts", 500, -20e3, 2e3,2000,0,2000);
  TH2F* h2_allimplantbetagamma_spectrum_after_ionbeta_square_forwardmatch = new TH2F("allimplantbetagamma_spectrum_after_ionbeta_square_forwardmatch", "allimplantbetagamma_spectrum_after_ionbeta_square_forwardmatch", 2000,0,2000,2000,0,2000);
  TH2F* h2_allimplantbetagamma_spectrum_after_ionbeta_square_time_forwardmatch = new TH2F("allimplantbetagamma_spectrum_after_ionbeta_square_time_forwardmatch", "allimplantbetagamma_spectrum_after_ionbeta_square_time_forwardmatch", 2000,0,2000,200,-100,100);

  // Histograms for forward backwards matche
  TH2F* h2_gatedimplantbetagamma_betagamma_dt_vs_implantbeta_dt = new TH2F("gatedimplantbetagamma_betagamma_dt_vs_implantbeta_dt", "Implant-Beta-Gamma ImplantBeta dt vs BetaGamma dt;betagamma dt (ns); ionbeta dt(ns)", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD, 1000, -20, 2e3);
  TH2F* h2_gatedimplantbetagamma_gammaenergy_vs_implantbeta_dt = new TH2F("gatedimplantbetagamma_gammaenergy_dt_vs_implantbeta_dt", "Implant-Beta-Gamma ImplantBeta dt vs Gamma energy;betagamma dt (ns); Energy (keV)", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD, 2000, 0, 2000);
  TH2F* h2_allimplantbetagamma_betagamma_dt_vs_implantbeta_dt = new TH2F("allimplantbetagamma_betagamma_dt_vs_implantbeta_dt", "Implant-Beta-Gamma ImplantBeta dt vs BetaGamma dt;betagamma dt (ns); ionbeta dt(ns)", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD, 1000, -20, 2e3);
  TH2F* h2_allimplantbetagamma_gammaenergy_vs_implantbeta_dt = new TH2F("allimplantbetagamma_gammaenergy_dt_vs_implantbeta_dt", "Implant-Beta-Gamma ImplantBeta dt vs Gamma energy;betagamma dt (ns); Energy (keV)", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD, constants::TIME_THRESHOLD, 2000, 0, 2000);

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

  // Diagnostic histograms for ionbetagamma mathces to known gammas
  TH1F* h1_aida_implant_beta_dt_knowngamma_match = new TH1F("aida_implant_beta_dt_knowngamma_match", "Implant-Decay #Deltat (Known Gamma Match);Implant-Decay #deltat (); Counts/", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD/constants::TIME_SCALE, constants::TIME_THRESHOLD/constants::TIME_SCALE);
  TH2F* h2_implant_energy_decay_energy_knowngamma_match = new TH2F("implant_energy_decay_energy_knowngamma_match", "Implant Energy vs Decay Energy (Known Gamma Match);Decay Energy (keV);Implant Energy (MeV)", 1500/20, 0, 1500, 7000/20, 0, 7000);
  TH2F* h2_implant_energy_dt_knowngamma_match = new TH2F("implant_energy_dt_knowngamma_match", "Implant Energy vs Implant-Decay dt (Known Gamma Match)", 7000/20, 0, 7000, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_decay_energy_dt_knowngamma_match = new TH2F("decay_energy_dt_knowngamma_match", "Decay Energy vs Implant-Decay dt (Known Gamma Match)", 1500/20, 0, 1500, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_strip_dt_knowngamma_match = new TH2F("strip_dt_knowngamma_match", "Strip XY vs Implant-Decay dt (Known Gamma Match)", 528, 0, 528, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);

  // Diagnostic histograms for ionbetagamma mathces to unknown gammas
  TH1F* h1_aida_implant_beta_dt_unknowngamma_match = new TH1F("aida_implant_beta_dt_unknowngamma_match", "Implant-Decay #Deltat (Unknown Gamma Match);Implant-Decay #deltat (); Counts/", constants::IMPDECAY_TIME_BINS, -constants::TIME_THRESHOLD/constants::TIME_SCALE, constants::TIME_THRESHOLD/constants::TIME_SCALE);
  TH2F* h2_implant_energy_decay_energy_unknowngamma_match = new TH2F("implant_energy_decay_energy_unknowngamma_match", "Implant Energy vs Decay Energy (Unknown Gamma Match);Decay Energy (keV);Implant Energy (MeV)", 1500/20, 0, 1500, 7000/20, 0, 7000);
  TH2F* h2_implant_energy_dt_unknowngamma_match = new TH2F("implant_energy_dt_unknowngamma_match", "Implant Energy vs Implant-Decay dt (Unknown Gamma Match)", 7000/20, 0, 7000, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_decay_energy_dt_unknowngamma_match = new TH2F("decay_energy_dt_unknowngamma_match", "Decay Energy vs Implant-Decay dt (Unknown Gamma Match)", 1500/20, 0, 1500, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);
  TH2F* h2_strip_dt_unknowngamma_match = new TH2F("strip_dt_unknowngamma_match", "Strip XY vs Implant-Decay dt (Unknown Gamma Match)", 528, 0, 528, constants::IMPDECAY_TIME_BINS, 0, constants::TIME_THRESHOLD);

  // Subsequent matches

  // *************************************************************************************
  // ****************************** FILL MAPS WITH EVENTS ********************************
  // *************************************************************************************

  // Loop over the implant, decay and germanium trees and fill the maps
  std::cout << "Start filled the maps!" << std::endl << std::endl; // Read gated implant events
  while (gatedimplant_reader.Next()){
    if( *gatedimplant_dssd==constants::DSSD && *gatedimplant_bplast==0 && *gatedimplant_e>3.5){
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
    if( *implant_dssd==constants::DSSD && *implant_bplast==0  && *implant_e>3.5e3){
      // Create entry for map
      all_implants_map.emplace(
        *implant_time,
        std::make_tuple(*implant_x, *implant_y, *implant_e, *implant_ex, *implant_ey, *implant_spill, *implant_bplast, *implant_dssd, NOT_INTERRUPTED)
      );

      // Define deadtime range in manager
      deadtimeWindowManager.addRange((uint64_t)*implant_time, (uint64_t)*implant_time+constants::IMPLANT_DEAD_TIME, *implant_x, *implant_y);

    }
  }
  std::cout << "Finished filling the implant map" << std::endl;
  std::cout << "Number of All implant events cut on region: " << all_implants_map.size() << std::endl << std::endl;

  //************* DEBUG ************** // std::cout << "[DEBUG] Printing deadtime ranges in the manager:\n";
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
    h1_gamma_spectrum->Fill(*germanium_energy);
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

  std::cout << "Started Gated Implant-Decay matching routine..." << std::endl;

  // Define & declare variebles to be used inside loop of gated implant events
  double gatedimplant_pos_x;
  double gatedimplant_pos_y;
  bool found_forward_match;
  bool found_backward_match;
  int forwards_implantbeta_candidate_counter;
  int backwards_implantbeta_candidate_counter;
  int matched_implantdecays_counter = 0;
  int matched_backwards_implantdecays_counter = 0;
  int matched_implantbetagamma_counter = 0;
  int matched_betagamma_counter = 0;
  std::vector<BetaCandidateInfo> beta_candidate_data_vector;

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
    forwards_implantbeta_candidate_counter = 0;     
    backwards_implantbeta_candidate_counter = 0;     

    // Reset flags for positive match in new loop
    found_forward_match = false;
    found_backward_match = false;

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
    auto decay_start = good_decays_map.lower_bound(last_gatedimplant_time - 5*constants::TIME_THRESHOLD);

    // Now loop over decay events starting from our decay start defined above untoll we pass our time threshold
    for(auto decay_evt = decay_start; decay_evt != good_decays_map.end(); decay_evt++){

      // Break out of loop if decay events are now outside of time window
      // if ( decay_evt->first > last_gatedimplant_time + 5*constants::TIME_THRESHOLD ){ break; }

      // Break out of loop if we have found a forward  & backward match and we are not checking all candidates 
      // if ( !constants::CHECK_BETA_CANDITATES && found_forward_match && found_backward_match ){ break; }

      // Unpack event variables for current decay event
      auto [decay_x, decay_y, decay_e, decay_ex, decay_ey, decay_spill, decay_bplast, decay_dssd] = decay_evt->second;

      // Check for noisy decay branch strips and skip
      if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_DECAY, decay_x) ){ continue; }
      if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_DECAY, decay_y) ){ continue; }

      // Check if decay event is onspill and skip if defined by user 
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
          if ( deadtimeWindowManager.contains(decay_evt->first, decay_x, decay_y) ){
            h1_deadtime_implant_vetoed_decay_candidates_time->Fill(decay_evt->first);
            std::cout << "Rejected a forwards ion-beta candidate due to implant deadtime." << std::endl;
            continue; 
          }
          
          forwards_implantbeta_candidate_counter++; // Increase counter for beta canditade event

          BetaType beta_type = CANDIDATE; // Set tag for beta
          CorrelationType correlation_type = FORWARDS; // Set tag for correlation type

          switch (forwards_implantbeta_candidate_counter){
          case 1: {
            // First Match 
            matched_implantdecays_counter++; // Increase counter for succesfull matched implant decay
            found_forward_match = true; // Change flag for succesfull forward implant decay match
            beta_type = MATCH; // Overwrite beta_type
            break;
          }
          case 2: {
            // Second Match 
            beta_type = SECOND_MATCH; // Overwrite beta_type
            break;
          } 
          case 3: {
            // Third Match 
            beta_type = THIRD_MATCH; // Overwrite beta_type
            break;
          } 
          case 10: {
            // Tenth Match 
            beta_type = TENTH_MATCH; // Overwrite beta_type
            break;
          } 
          default:
            // No longer tracking any matches
            break;
          }

          // Make new struct to contain match info and add parameters
          BetaCandidateInfo matchedBetaInfo;
          matchedBetaInfo.setBetaParameters(correlation_type, beta_type, time_diff, decay_evt->first, gatedimplant_pos_x, gatedimplant_pos_y, e, decay_x, decay_y, decay_e, decay_spill);
          // Add beta data to the vector
          beta_candidate_data_vector.emplace_back( matchedBetaInfo );

        }

        // *************************************************************************************
        // ****************************** FOUND BACKWARD BETA CANDIDATE ************************
        // *************************************************************************************

        // Check if decay event falls within time threshold and decay event occures before implant event 
        if (time_diff < 0 && -time_diff < constants::TIME_THRESHOLD) {

          // Check if decay occures within a deadtime induced my implant in AIDA
          if ( deadtimeWindowManager.contains(decay_evt->first, decay_x, decay_y) ){
            std::cout << "Rejected backwards ion-beta candidate due to implant deadtime." << std::endl;
            continue; 
          }

          backwards_implantbeta_candidate_counter++; // Increase counter for beta canditade event

          BetaType beta_type = CANDIDATE; // Set tag for beta
          CorrelationType correlation_type = BACKWARDS; // Set tag for correlation type

          switch (backwards_implantbeta_candidate_counter){
          case 1: {
            // First
            matched_backwards_implantdecays_counter++; // Increase counter for succesfull matched implant decay
            found_backward_match = true; // Change flag for succesfull forward implant decay match
            beta_type = MATCH; // Overwrite beta_type
            break;
          }
          case 2: {
            // Second Match 
            beta_type = SECOND_MATCH; // Overwrite beta_type
            break;
          } 
          case 3: {
            // Third Match 
            beta_type = THIRD_MATCH; // Overwrite beta_type
            break;
          } 
          case 10: {
            // Tenth Match 
            beta_type = TENTH_MATCH; // Overwrite beta_type
            break;
          } 
          default:
            // No longer tracking any matches
            break;
          }

          // Make new struct to contain match info and add parameters
          BetaCandidateInfo matchedBetaInfo;
          matchedBetaInfo.setBetaParameters(correlation_type, beta_type, time_diff, decay_evt->first, gatedimplant_pos_x, gatedimplant_pos_y, e, decay_x, decay_y, decay_e, decay_spill);
          // Add beta data to the vector
          beta_candidate_data_vector.emplace_back( matchedBetaInfo );

        }

      }

    } // End of decay event loop

    // Make implant-decay infor struct and fill it
    ImplantDecayInfo gatedImplantDecayInfo;
    gatedImplantDecayInfo.implantDecayMatchMultiplicities = std::make_pair(forwards_implantbeta_candidate_counter, backwards_implantbeta_candidate_counter);
    gatedImplantDecayInfo.betaCandidateInfoVector = beta_candidate_data_vector;
    // Fill matched candidate map
    gated_matched_decays_map.emplace(last_gatedimplant_time, gatedImplantDecayInfo);


  } // End of gated implant event loop

  std::cout << "Finalised Gated Implant-Decay matching routine!" << std::endl << std::endl;



  // *************************************************************************************
  // ****************** ION-BETA MATCHING ****** LOOP OVER ALL IMPLANTS ****************
  // *************************************************************************************
  if (constants::ALL_IMPLANTS_DECAY_MATCH){


    std::cout << "Started All Implant-Decay matching routine..." << std::endl;

    // Flush vector for new beta candidates
    beta_candidate_data_vector.clear();

    // Loop over all gated implant events in map and perform a beta match
    for (auto imp_evt = all_implants_map.begin(); imp_evt != all_implants_map.end(); imp_evt++){
      
      // Empty vector containing previous beta candidate
      beta_candidate_data_vector.clear();

      // Unpack event variables for current gated implant event
      auto [x, y, e, ex, ey, spill, bplast, dssd, interruption_type] = imp_evt->second;

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
      forwards_implantbeta_candidate_counter = 0;     
      backwards_implantbeta_candidate_counter = 0;     

      // Reset flags for positive match in new loop
      found_forward_match = false;
      found_backward_match = false;

      last_gatedimplant_time = imp_evt->first; // Unpack white rabbit time of gated implant

      // Set gated implant position
      gatedimplant_pos_x = x;
      gatedimplant_pos_y = y;

      h2_aida_implant_xy->Fill(x,y); // Fill Histogram with gated implant position

      // *************************************************************************************
      // ****************************** LOOP OVER VALID DECAYS *******************************
      // *************************************************************************************

      // Find the decay event corresponding to the start of our decay loop using our time window
      // The inital decay event will be the one whose time corresponds to our time threshould before the implant occured
      auto decay_start = good_decays_map.lower_bound(last_gatedimplant_time - 2*constants::TIME_THRESHOLD);

      // Now loop over decay events starting from our decay start defined above untoll we pass our time threshold
      for(auto decay_evt = decay_start; decay_evt != good_decays_map.end(); decay_evt++){

        // Break out of loop if decay events are now outside of time window
        if ( decay_evt->first > last_gatedimplant_time + 2*constants::TIME_THRESHOLD ){ break; }

        // Break out of loop if we have found a forward  & backward match and we are not checking all candidates 
        if ( !constants::CHECK_BETA_CANDITATES && found_forward_match && found_backward_match ){ break; }

        // Unpack event variables for current decay event
        auto [decay_x, decay_y, decay_e, decay_ex, decay_ey, decay_spill, decay_bplast, decay_dssd] = decay_evt->second;

        // Check for noisy decay branch strips and skip
        if ( isNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_DECAY, decay_x) ){ continue; }
        if ( isNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_DECAY, decay_y) ){ continue; }

        // Check if decay event is onspill and skip if defined by user 
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
            if ( deadtimeWindowManager.contains(decay_evt->first, decay_x, decay_y) ){
              h1_deadtime_implant_vetoed_decay_candidates_time->Fill(decay_evt->first);
              std::cout << "Rejected a forwards ion-beta candidate due to implant deadtime." << std::endl;
              continue; 
            }
            
            forwards_implantbeta_candidate_counter++; // Increase counter for beta canditade event

            BetaType beta_type = CANDIDATE; // Set tag for beta
            CorrelationType correlation_type = FORWARDS; // Set tag for correlation type

            switch (forwards_implantbeta_candidate_counter){
            case 1: {
              // First Match 
              matched_implantdecays_counter++; // Increase counter for succesfull matched implant decay
              found_forward_match = true; // Change flag for succesfull forward implant decay match
              beta_type = MATCH; // Overwrite beta_type
              break;
            }
            case 2: {
              // Second Match 
              beta_type = SECOND_MATCH; // Overwrite beta_type
              break;
            } 
            case 3: {
              // Third Match 
              beta_type = THIRD_MATCH; // Overwrite beta_type
              break;
            } 
            case 10: {
              // Tenth Match 
              beta_type = TENTH_MATCH; // Overwrite beta_type
              break;
            } 
            default:
              // No longer tracking any matches
              break;
            }

            // Make new struct to contain match info and add parameters
            BetaCandidateInfo matchedBetaInfo;
            matchedBetaInfo.setBetaParameters(correlation_type, beta_type, time_diff, decay_evt->first, gatedimplant_pos_x, gatedimplant_pos_y, e, decay_x, decay_y, decay_e, decay_spill);
            // Add beta data to the vector
            beta_candidate_data_vector.emplace_back( matchedBetaInfo );

          }

          // *************************************************************************************
          // ****************************** FOUND BACKWARD BETA CANDIDATE ************************
          // *************************************************************************************

          // Check if decay event falls within time threshold and decay event occures before implant event 
          if (time_diff < 0 && -time_diff < constants::TIME_THRESHOLD) {

            // Check if decay occures within a deadtime induced my implant in AIDA
            if ( deadtimeWindowManager.contains(decay_evt->first, decay_x, decay_y) ){
              std::cout << "Rejected backwards ion-beta candidate due to implant deadtime." << std::endl;
              continue; 
            }

            backwards_implantbeta_candidate_counter++; // Increase counter for beta canditade event

            BetaType beta_type = CANDIDATE; // Set tag for beta
            CorrelationType correlation_type = BACKWARDS; // Set tag for correlation type

            switch (backwards_implantbeta_candidate_counter){
            case 1: {
              // First
              matched_backwards_implantdecays_counter++; // Increase counter for succesfull matched implant decay
              found_backward_match = true; // Change flag for succesfull forward implant decay match
              beta_type = MATCH; // Overwrite beta_type
              break;
            }
            case 2: {
              // Second Match 
              beta_type = SECOND_MATCH; // Overwrite beta_type
              break;
            } 
            case 3: {
              // Third Match 
              beta_type = THIRD_MATCH; // Overwrite beta_type
              break;
            } 
            case 10: {
              // Tenth Match 
              beta_type = TENTH_MATCH; // Overwrite beta_type
              break;
            } 
            default:
              // No longer tracking any matches
              break;
            }

            // Make new struct to contain match info and add parameters
            BetaCandidateInfo matchedBetaInfo;
            matchedBetaInfo.setBetaParameters(correlation_type, beta_type, time_diff, decay_evt->first, gatedimplant_pos_x, gatedimplant_pos_y, e, decay_x, decay_y, decay_e, decay_spill);
            // Add beta data to the vector
            beta_candidate_data_vector.emplace_back( matchedBetaInfo );

          }

        }

      } // End of decay event loop

      // Make implant-decay infor struct and fill it
      ImplantDecayInfo allImplantDecayInfo;
      allImplantDecayInfo.implantDecayMatchMultiplicities = std::make_pair(forwards_implantbeta_candidate_counter, backwards_implantbeta_candidate_counter);
      allImplantDecayInfo.betaCandidateInfoVector = beta_candidate_data_vector;
      // Fill matched candidate map
      all_matched_decays_map.emplace(last_gatedimplant_time, allImplantDecayInfo);


    } // End of gated implant event loop

    std::cout << "Finalised All Implant-Decay matching routine!" << std::endl << std::endl;

  }

  // *************************************************************************************
  // ************************** MATCHED IMPLANT - DECAY PLOTTING *************************
  // *************************************************************************************
  

  std::cout << "Started filling histograms..." << std::endl;
  
  // Loop over each gated implant match 
  for ( auto matched_decay_evt = gated_matched_decays_map.begin(); matched_decay_evt!=gated_matched_decays_map.end(); matched_decay_evt++ ){

    // Unpack the implant time of the gated implant
    int64_t gimp_time = matched_decay_evt->first;
    auto& implantDecayInfo = matched_decay_evt->second;
    auto [forwards_candidate_multiplicity, backwards_candidate_multiplicity] = implantDecayInfo.implantDecayMatchMultiplicities;

    // Plot candidate multiplicity
    if(forwards_candidate_multiplicity>0) {h1_implantbeta_candidate_multiplicity_forwards->Fill(forwards_candidate_multiplicity);}
    if(backwards_candidate_multiplicity>0) {h1_implantbeta_candidate_multiplicity_backwards->Fill(backwards_candidate_multiplicity);}

    // Loop over all matched beta candidate for each gated implant
    for ( auto& beta_candidate_data : implantDecayInfo.betaCandidateInfoVector ){

      // Unpack the data for each beta candidate
      auto [correlation_type, beta_type, time_diff, decay_time, gimp_x, gimp_y, gimp_e, decay_x, decay_y, decay_e, decay_spill] = beta_candidate_data.getBetaParameters();

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

      if ( beta_type == SECOND_MATCH ){
        h1_aida_implant_beta_secondmatch_dt->Fill(time_diff/constants::TIME_SCALE);
      }

      if ( beta_type == THIRD_MATCH ){
        h1_aida_implant_beta_thirdmatch_dt->Fill(time_diff/constants::TIME_SCALE);
      }

      if ( beta_type == TENTH_MATCH ){
        h1_aida_implant_beta_tenthmatch_dt->Fill(time_diff/constants::TIME_SCALE);
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

    }
    
  }
  
  std::cout << "Finished filling histograms!" << std::endl << std::endl;
 
  // *************************************************************************************
  // ************************** GATED MATCHED DECAY - GAMMA CORRELATION ******************
  // *************************************************************************************

  std::cout << "Started Gated Matched Beta - Gamma matching routine..." << std::endl;
  
  // Loop over all matched decay events
  for( auto matched_decay_evt = gated_matched_decays_map.begin(); matched_decay_evt != gated_matched_decays_map.end(); matched_decay_evt++ ){
    // Apply cut to number of candidates
    auto& implantDecayInfo = matched_decay_evt->second;
    auto [forwards_candidate_multiplicity, backwards_candidate_multiplicity] = implantDecayInfo.implantDecayMatchMultiplicities;
    if ( forwards_candidate_multiplicity > constants::BETA_GAMMA_CANDIDATE_CUT) continue;
    if ( backwards_candidate_multiplicity > constants::BETA_GAMMA_CANDIDATE_CUT) continue;

    for ( auto& matched_beta_data : implantDecayInfo.betaCandidateInfoVector){

      // Unpack matched decay event variables
      CorrelationType correlation_type  = matched_beta_data.correlationType;
      BetaType beta_type                = matched_beta_data.betaType;
      int64_t last_matched_decay_time   = matched_beta_data.decayTime;
      int decay_spill                   = matched_beta_data.decaySpill;
      Long64_t ionbetaTimeDiff          = matched_beta_data.timeDiff;
      double decayPosX                  = matched_beta_data.decayPosX;
      double decayPosY                  = matched_beta_data.decayPosY;

      // Check if decay event is onspill and skip if defined by user 
      if( constants::ONLY_OFFSPILL_DECAY && decay_spill == 1){ continue; }

      // Find the germanium event starting at the same time as the decay event (50 microsecond grace period)
      auto germanium_start = germanium_map.lower_bound(last_matched_decay_time - 50e3);

      // *************************************************************************************
      // **************************  LOOP OVER GERMANIUM EVENTS ******************************
      // *************************************************************************************

      // Loop over all germanium events between decay event and end of prompt gamma window
      for( auto germanium_evt = germanium_start; germanium_evt != germanium_map.end(); germanium_evt++ ){

        int64_t time_diff = - ( last_matched_decay_time - germanium_evt->first ); // Reverse like in jeroens code

        if (correlation_type==FORWARDS) h1_gatedimplantbetagamma_spectrum_after_ionbeta_dt_forwardmatch->Fill(time_diff);
        if (correlation_type==BACKWARDS) h1_gatedimplantbetagamma_spectrum_after_ionbeta_dt_backwardmatch->Fill(time_diff);
        
        // Unpack germanium event items outside prompt window
        auto [germanium_energy, germanium_spill] = germanium_evt->second; 

        if ( time_diff > 50e3 ){ break; } // At this point dt is > 50us

        h2_gatedimplantbetagamma_betagamma_dt_vs_implantbeta_dt->Fill(ionbetaTimeDiff, time_diff);
        if (correlation_type==FORWARDS) h2_gatedimplantbetagamma_spectrum_after_ionbeta_dt_energy_forwardmatch->Fill(time_diff,germanium_energy);
        if (correlation_type==BACKWARDS) h2_gatedimplantbetagamma_spectrum_after_ionbeta_dt_energy_backwardmatch->Fill(time_diff,germanium_energy);

        // Check if germanium event is within prompt window
        if ( time_diff > constants::PROMPT_GAMMA_WINDOW.first && time_diff < constants::PROMPT_GAMMA_WINDOW.second ){

          h2_gatedimplantbetagamma_gammaenergy_vs_implantbeta_dt->Fill(ionbetaTimeDiff, germanium_energy);

          // Forward correlations
          if (correlation_type==FORWARDS){
            // Fill gamma energy spectrum and increase counter
            h1_gatedimplantbetagamma_spectrum_after_ionbeta_forwardmatch->Fill(germanium_energy);
            h2_gatedimplantbetagamma_energy_vs_decaystrip_forwardmatch->Fill(decayPosX, germanium_energy);
            h2_gatedimplantbetagamma_energy_vs_decaystrip_forwardmatch->Fill(decayPosY+400, germanium_energy);
            
            matched_implantbetagamma_counter++;

            // Provide gamma energy to info struct
            matched_beta_data.setGammaParameters(germanium_energy);

            for( auto germanium_evt2 = germanium_start; germanium_evt2 != germanium_map.end(); germanium_evt2++ ){
              int64_t time_diff2 = - ( last_matched_decay_time - germanium_evt2->first ); // Reverse like in jeroens code
              if ( time_diff2 > 50e3 ){ break; }
              auto [germanium_energy2, useless] = germanium_evt2->second;
              
              if (TMath::Abs(germanium_evt2->first - germanium_evt->first ) < 500 && germanium_evt2 != germanium_evt ) { h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_forwardmatch->Fill(germanium_energy,germanium_energy2); }
              if (germanium_evt != germanium_evt2) { h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_time_forwardmatch->Fill(germanium_energy,germanium_evt2->first - germanium_evt->first); }
            }
          }

          // Backwards correlations
          if (correlation_type==BACKWARDS){
            // Fill gamma energy spectrum and increase counter
            h1_gatedimplantbetagamma_spectrum_after_ionbeta_backwardmatch->Fill(germanium_energy);
            h2_gatedimplantbetagamma_energy_vs_decaystrip_backwardmatch->Fill(decayPosX, germanium_energy);
            h2_gatedimplantbetagamma_energy_vs_decaystrip_backwardmatch->Fill(decayPosY+400, germanium_energy);
            
            for( auto germanium_evt2 = germanium_start; germanium_evt2 != germanium_map.end(); germanium_evt2++ ){
              int64_t time_diff2 = - ( last_matched_decay_time - germanium_evt2->first ); // Reverse like in jeroens code
              if ( time_diff2 > 50e3 ){ break; }
              auto [germanium_energy2, useless] = germanium_evt2->second;
              
              if (TMath::Abs(germanium_evt2->first - germanium_evt->first ) < 500 && germanium_evt2 != germanium_evt ) { h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_backwardmatch->Fill(germanium_energy,germanium_energy2); }
              if (germanium_evt != germanium_evt2) { h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_time_backwardmatch->Fill(germanium_energy,germanium_evt2->first - germanium_evt->first); }
            }
          }
            
        }

      } // End of germanium event loop

    }
      
  } // End of matched decay event loop

  std::cout << "Finished Gated Matched Beta - Gamma matching routine!" << std::endl << std::endl;


  // *************************************************************************************
  // ************************** ALL MATCHED DECAY - GAMMA CORRELATION ******************
  // *************************************************************************************

  if (constants::ALL_IMPLANTS_DECAY_MATCH){


    std::cout << "Started all Beta - Gamma matching routine..." << std::endl;
    
    // Loop over all matched decay events
    for( auto matched_decay_evt = all_matched_decays_map.begin(); matched_decay_evt != all_matched_decays_map.end(); matched_decay_evt++ ){
      // Apply cut to number of candidates
      auto& implantDecayInfo = matched_decay_evt->second;
      auto [forwards_candidate_multiplicity, backwards_candidate_multiplicity] = implantDecayInfo.implantDecayMatchMultiplicities;
      if ( forwards_candidate_multiplicity > constants::BETA_GAMMA_CANDIDATE_CUT) continue;
      if ( backwards_candidate_multiplicity > constants::BETA_GAMMA_CANDIDATE_CUT) continue;

      for ( auto& matched_beta_data : implantDecayInfo.betaCandidateInfoVector){

        // Unpack matched decay event variables
        CorrelationType correlation_type  = matched_beta_data.correlationType;
        BetaType beta_type                = matched_beta_data.betaType;
        int64_t last_matched_decay_time   = matched_beta_data.decayTime;
        int decay_spill                   = matched_beta_data.decaySpill;
        Long64_t ionbetaTimeDiff          = matched_beta_data.timeDiff;
        double decayPosX                  = matched_beta_data.decayPosX;
        double decayPosY                  = matched_beta_data.decayPosY;

        // Check if decay event is onspill and skip if defined by user 
        if( constants::ONLY_OFFSPILL_DECAY && decay_spill == 1){ continue; }

        // Find the germanium event starting at the same time as the decay event (50 microsecond grace period)
        auto germanium_start = germanium_map.lower_bound(last_matched_decay_time - 50e3);

        // *************************************************************************************
        // **************************  LOOP OVER GERMANIUM EVENTS ******************************
        // *************************************************************************************

        // Loop over all germanium events between decay event and end of prompt gamma window
        for( auto germanium_evt = germanium_start; germanium_evt != germanium_map.end(); germanium_evt++ ){

          int64_t time_diff = - ( last_matched_decay_time - germanium_evt->first ); // Reverse like in jeroens code

          if (correlation_type==FORWARDS) h1_allimplantbetagamma_spectrum_after_ionbeta_dt_forwardmatch->Fill(time_diff);
          if (correlation_type==BACKWARDS) h1_allimplantbetagamma_spectrum_after_ionbeta_dt_backwardmatch->Fill(time_diff);
          
          // Unpack germanium event items outside prompt window
          auto [germanium_energy, germanium_spill] = germanium_evt->second; 

          if ( time_diff > 50e3 ){ break; } // At this point dt is > 50us

          h2_allimplantbetagamma_betagamma_dt_vs_implantbeta_dt->Fill(ionbetaTimeDiff, time_diff);
          if (correlation_type==FORWARDS) h2_allimplantbetagamma_spectrum_after_ionbeta_dt_energy_forwardmatch->Fill(time_diff,germanium_energy);
          if (correlation_type==BACKWARDS) h2_allimplantbetagamma_spectrum_after_ionbeta_dt_energy_backwardmatch->Fill(time_diff,germanium_energy);

          // Check if germanium event is within prompt window
          if ( time_diff > constants::PROMPT_GAMMA_WINDOW.first && time_diff < constants::PROMPT_GAMMA_WINDOW.second ){

            h2_allimplantbetagamma_gammaenergy_vs_implantbeta_dt->Fill(ionbetaTimeDiff, germanium_energy);

            // Forward correlations
            if (correlation_type==FORWARDS){
              // Fill gamma energy spectrum and increase counter
              h1_allimplantbetagamma_spectrum_after_ionbeta_forwardmatch->Fill(germanium_energy);
              h2_allimplantbetagamma_energy_vs_decaystrip_forwardmatch->Fill(decayPosX, germanium_energy);
              h2_allimplantbetagamma_energy_vs_decaystrip_forwardmatch->Fill(decayPosY+400, germanium_energy);
              
              matched_implantbetagamma_counter++;

              // Provide gamma energy to info struct
              matched_beta_data.setGammaParameters(germanium_energy);

              for( auto germanium_evt2 = germanium_start; germanium_evt2 != germanium_map.end(); germanium_evt2++ ){
                int64_t time_diff2 = - ( last_matched_decay_time - germanium_evt2->first ); // Reverse like in jeroens code
                if ( time_diff2 > 50e3 ){ break; }
                auto [germanium_energy2, useless] = germanium_evt2->second;
                
                if (TMath::Abs(germanium_evt2->first - germanium_evt->first ) < 500 && germanium_evt2 != germanium_evt ) { h2_allimplantbetagamma_spectrum_after_ionbeta_square_forwardmatch->Fill(germanium_energy,germanium_energy2); }
                if (germanium_evt != germanium_evt2) { h2_allimplantbetagamma_spectrum_after_ionbeta_square_time_forwardmatch->Fill(germanium_energy,germanium_evt2->first - germanium_evt->first); }
              }
            }

            // Backwards correlations
            if (correlation_type==BACKWARDS){
              // Fill gamma energy spectrum and increase counter
              h1_allimplantbetagamma_spectrum_after_ionbeta_backwardmatch->Fill(germanium_energy);
              h2_allimplantbetagamma_energy_vs_decaystrip_backwardmatch->Fill(decayPosX, germanium_energy);
              h2_allimplantbetagamma_energy_vs_decaystrip_backwardmatch->Fill(decayPosY+400, germanium_energy);
              
              for( auto germanium_evt2 = germanium_start; germanium_evt2 != germanium_map.end(); germanium_evt2++ ){
                int64_t time_diff2 = - ( last_matched_decay_time - germanium_evt2->first ); // Reverse like in jeroens code
                if ( time_diff2 > 50e3 ){ break; }
                auto [germanium_energy2, useless] = germanium_evt2->second;
                
                if (TMath::Abs(germanium_evt2->first - germanium_evt->first ) < 500 && germanium_evt2 != germanium_evt ) { h2_allimplantbetagamma_spectrum_after_ionbeta_square_backwardmatch->Fill(germanium_energy,germanium_energy2); }
                if (germanium_evt != germanium_evt2) { h2_allimplantbetagamma_spectrum_after_ionbeta_square_time_backwardmatch->Fill(germanium_energy,germanium_evt2->first - germanium_evt->first); }
              }
            }
              
          }

        } // End of germanium event loop

      }
        
    } // End of matched decay event loop

    std::cout << "Finished Matched Beta - Gamma matching routine!" << std::endl << std::endl;

  }


  // *************************************************************************************
  // ************************** KNOWN GAMMA - DECAY BACKREFERENCING **********************
  // *************************************************************************************

  // Only run this routine if gated isotope has known gammas defined in script
  auto knownBetaDelayedGammas = constants::KNOWN_BETA_DELAYED_GAMMAS.find(constants::ISOTOPE_TREE); 
  if ( knownBetaDelayedGammas != constants::KNOWN_BETA_DELAYED_GAMMAS.end() ){

    std::cout << "Started Known Gamma - Decay backrefereincing routine..." << std::endl;
    
    auto& knownGammaSet = knownBetaDelayedGammas->second; // Get Known Gamma Set

    // Loop over all matched decay events
    for( auto matched_decay_evt = gated_matched_decays_map.begin(); matched_decay_evt != gated_matched_decays_map.end(); matched_decay_evt++ ){
      // Apply cut to number of candidates
      auto& implantDecayInfo = matched_decay_evt->second;
      auto [forwards_candidate_multiplicity, backwards_candidate_multiplicity] = implantDecayInfo.implantDecayMatchMultiplicities;
      if ( forwards_candidate_multiplicity > constants::BETA_GAMMA_CANDIDATE_CUT) continue;

      for ( auto& matched_beta_data : implantDecayInfo.betaCandidateInfoVector){

        // Define flag for if beta has a known gamma match
        bool hasKnownGammaMatch = matched_beta_data.hasMatchedGammaInRange(knownGammaSet);

        // Unpack matched decay event variables
        auto [correlation_type, beta_type, time_diff, decay_time, gimp_x, gimp_y, gimp_e, decay_x, decay_y, decay_e, decay_spill] = matched_beta_data.getBetaParameters();

        // Check if decay event is onspill and skip if defined by user 
        if( constants::ONLY_OFFSPILL_DECAY && decay_spill == 1){ continue; }

        if ( hasKnownGammaMatch ) { h1_aida_implant_beta_dt_knowngamma_match->Fill(time_diff/constants::TIME_SCALE); }
        if ( !hasKnownGammaMatch ) { h1_aida_implant_beta_dt_unknowngamma_match->Fill(time_diff/constants::TIME_SCALE); }

        if ( hasKnownGammaMatch && correlation_type == FORWARDS ){
          // Fill Diagnostic Histograms
          h2_implant_energy_decay_energy_knowngamma_match->Fill(decay_e, gimp_e);
          h2_implant_energy_dt_knowngamma_match->Fill(gimp_e, time_diff);
          h2_decay_energy_dt_knowngamma_match->Fill(decay_e, time_diff);
          h2_strip_dt_knowngamma_match->Fill(decay_x, time_diff);
          h2_strip_dt_knowngamma_match->Fill(decay_y + 400, time_diff);
        }

        if ( !hasKnownGammaMatch && correlation_type == FORWARDS ){
          // Fill Diagnostic Histograms
          h1_aida_implant_beta_dt_unknowngamma_match->Fill(time_diff);
          h2_implant_energy_decay_energy_unknowngamma_match->Fill(decay_e, gimp_e);
          h2_implant_energy_dt_unknowngamma_match->Fill(gimp_e, time_diff);
          h2_decay_energy_dt_unknowngamma_match->Fill(decay_e, time_diff);
          h2_strip_dt_unknowngamma_match->Fill(decay_x, time_diff);
          h2_strip_dt_unknowngamma_match->Fill(decay_y + 400, time_diff);
        }

      }
        
    } // End of matched decay event loop

    std::cout << "Finished Known Gamma - Decay backrefereincing routine!" << std::endl;

  }


  // *************************************************************************************
  // ************************** ALL DECAY - GAMMA CORRELATION ****************************
  // *************************************************************************************

  std::cout << "Started All Beta - Gamma matching routine..." << std::endl;
  
  // Loop over all matched decay events
  for( auto good_decay_evt = good_decays_map.begin(); good_decay_evt != good_decays_map.end(); good_decay_evt++ ){

    // Unpack decay data 
    int64_t decay_time = good_decay_evt->first;
    auto [decay_x, decay_y, useless_3, useless_4, useless_5, decay_spill, useless_6, useless_7] = good_decay_evt->second;

    // Skip if decay is onspill and have not been selected
    if( constants::ONLY_OFFSPILL_DECAY && decay_spill == 1){ continue; }

    // Skip if decay occures in the deadtime induced by an implant
    if ( deadtimeWindowManager.contains(decay_time, decay_x, decay_y) ){ continue; }

    // Find the germanium event starting at the same time as the decay event (50 microsecond grace period)
    auto germanium_start = germanium_map.lower_bound(decay_time - 50e3);

    // *************************************************************************************
    // **************************  LOOP OVER GERMANIUM EVENTS ******************************
    // *************************************************************************************

    // Loop over all germanium events between decay event and end of prompt gamma window
    for( auto germanium_evt = germanium_start; germanium_evt != germanium_map.end(); germanium_evt++ ){

      int64_t time_diff = -( decay_time - germanium_evt->first ); // Reverse like in jeroens code
      // int64_t time_diff = germanium_evt->first - decay_time; // Like regular people
      h1_implantbetagamma_spectrum_before_ionbeta_dt->Fill(time_diff);
      // Unpack germanium event items outside prompt window
      auto [germanium_energy, germanium_spill] = germanium_evt->second; 
      h2_implantbetagamma_spectrum_before_ionbeta_dt_energy->Fill(time_diff,germanium_energy);

      if ( time_diff > constants::PROMPT_GAMMA_WINDOW.second + 5e3 ){ break; }

      // Check if germanium event is within prompt window
      if ( time_diff > constants::PROMPT_GAMMA_WINDOW.first && time_diff < constants::PROMPT_GAMMA_WINDOW.second ){
          
        // Fill gamma energy spectrum and increase counter
        h1_implantbetagamma_spectrum_before_ionbeta->Fill(germanium_energy); // Fill gamma energy spectrum
        matched_betagamma_counter++;
      }

    } // End of germanium event loop
      
  } // End of matched decay event loop

  std::cout << "Finished All Beta - Gamma matching routine!" << std::endl << std::endl;
    
  // *************************************************************************************
  // ****************************** PRINT OUT STATISTICS *********************************
  // *************************************************************************************
   
  // Print results of implant decay correlation algorithm
    std::cout << "Finished processing the data" << std::endl;
    std::cout << "Matched: " << matched_implantdecays_counter<< " out of " << gated_implants_map.size() << " gated implant events" << std::endl;
    std::cout << "Matched: " << matched_implantdecays_counter << " out of " << all_implants_map.size() << " implant events" << std::endl;
    std::cout << "Matched: " << matched_backwards_implantdecays_counter << " backwards gated implant events" << std::endl;
    std::cout << "Matched: " << matched_implantbetagamma_counter << " implant-beta-gamma events" << std::endl;
    std::cout << "Matched: " << matched_betagamma_counter << " beta-gamma events" << std::endl << std::endl;
  
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

  TDirectory* vetoedDecays = outputFile->mkdir("vetoed_decays");
  vetoedDecays->cd();
  h1_deadtime_implant_vetoed_decay_candidates_time->Write();
  h1_interrupted_implant_vetoed_decay_candidates_time->Write();
  gFile->cd();

  TDirectory* ionbetaMatch = outputFile->mkdir("ionbeta_match");
  ionbetaMatch->cd();
  h2_aida_implant_xy->Write();
  h1_aida_wr_times->Write();
  h2_aida_matched_xy->Write();
  h1_aida_implant_beta_dt->Write();
  h1_aida_implant_beta_secondmatch_dt->Write();
  h1_aida_implant_beta_thirdmatch_dt->Write();
  h1_aida_implant_beta_tenthmatch_dt->Write();
  gFile->cd();

  TDirectory* ionbetagammaForwardCoincidences = outputFile->mkdir("ionbetagamma_forward_coincidences");
  ionbetagammaForwardCoincidences->cd();
  h1_gatedimplantbetagamma_spectrum_after_ionbeta_forwardmatch->Write();
  h1_gatedimplantbetagamma_spectrum_after_ionbeta_dt_forwardmatch->Write();
  h2_gatedimplantbetagamma_spectrum_after_ionbeta_dt_energy_forwardmatch->Write();
  h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_forwardmatch->Write();
  h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_time_forwardmatch->Write();
  h2_gatedimplantbetagamma_energy_vs_decaystrip_forwardmatch->Write();
  if (constants::ALL_IMPLANTS_DECAY_MATCH){
    h1_allimplantbetagamma_spectrum_after_ionbeta_forwardmatch->Write();
    h1_allimplantbetagamma_spectrum_after_ionbeta_dt_forwardmatch->Write();
    h2_allimplantbetagamma_spectrum_after_ionbeta_dt_energy_forwardmatch->Write();
    h2_allimplantbetagamma_spectrum_after_ionbeta_square_forwardmatch->Write();
    h2_allimplantbetagamma_spectrum_after_ionbeta_square_time_forwardmatch->Write();
    h2_allimplantbetagamma_energy_vs_decaystrip_forwardmatch->Write();
  }
  gFile->cd();

  TDirectory* ionbetagammaBackwardCoincidences = outputFile->mkdir("ionbetagamma_backward_coincidences");
  ionbetagammaBackwardCoincidences->cd();
  h1_gatedimplantbetagamma_spectrum_after_ionbeta_backwardmatch->Write();
  h1_gatedimplantbetagamma_spectrum_after_ionbeta_dt_backwardmatch->Write();
  h2_gatedimplantbetagamma_spectrum_after_ionbeta_dt_energy_backwardmatch->Write();
  h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_backwardmatch->Write();
  h2_gatedimplantbetagamma_spectrum_after_ionbeta_square_time_backwardmatch->Write();
  h2_gatedimplantbetagamma_energy_vs_decaystrip_backwardmatch->Write();
  if (constants::ALL_IMPLANTS_DECAY_MATCH){
    h1_allimplantbetagamma_spectrum_after_ionbeta_backwardmatch->Write();
    h1_allimplantbetagamma_spectrum_after_ionbeta_dt_backwardmatch->Write();
    h2_allimplantbetagamma_spectrum_after_ionbeta_dt_energy_backwardmatch->Write();
    h2_allimplantbetagamma_spectrum_after_ionbeta_square_backwardmatch->Write();
    h2_allimplantbetagamma_spectrum_after_ionbeta_square_time_backwardmatch->Write();
    h2_allimplantbetagamma_energy_vs_decaystrip_backwardmatch->Write();
  }
  gFile->cd();

  TDirectory* ionbetagammaAllCoincidence = outputFile->mkdir("ionbetagamma_all_coincidence");
  ionbetagammaAllCoincidence->cd();
  h2_gatedimplantbetagamma_betagamma_dt_vs_implantbeta_dt->Write();
  h2_gatedimplantbetagamma_gammaenergy_vs_implantbeta_dt->Write();
  if (constants::ALL_IMPLANTS_DECAY_MATCH){
    h2_allimplantbetagamma_betagamma_dt_vs_implantbeta_dt->Write();
    h2_allimplantbetagamma_gammaenergy_vs_implantbeta_dt->Write();
  }
  gFile->cd();
  
  TDirectory* betagammaCoincidences = outputFile->mkdir("betagamma_coincidences");
  betagammaCoincidences->cd();
  h1_gamma_spectrum->Write();
  h1_implantbetagamma_spectrum_before_ionbeta->Write();
  h1_implantbetagamma_spectrum_before_ionbeta_dt->Write();
  h2_implantbetagamma_spectrum_before_ionbeta_dt_energy->Write();
  gFile->cd();

  TDirectory* knownGammaCoincidences = outputFile->mkdir("known_gamma_matches");
  knownGammaCoincidences->cd();
  h1_aida_implant_beta_dt_knowngamma_match->Write();
  h2_implant_energy_decay_energy_knowngamma_match->Write();
  h2_implant_energy_dt_knowngamma_match->Write();
  h2_decay_energy_dt_knowngamma_match->Write();
  h2_strip_dt_knowngamma_match->Write();
  h2_strip_dt_knowngamma_match->Write();
  gFile->cd();

  TDirectory* unknownGammaCoincidences = outputFile->mkdir("unknown_gamma_matches");
  unknownGammaCoincidences->cd();
  h1_aida_implant_beta_dt_unknowngamma_match->Write();
  h2_implant_energy_decay_energy_unknowngamma_match->Write();
  h2_implant_energy_dt_unknowngamma_match->Write();
  h2_decay_energy_dt_unknowngamma_match->Write();
  h2_strip_dt_unknowngamma_match->Write();
  h2_strip_dt_unknowngamma_match->Write();
  gFile->cd();
  
  nt_aida_implant_beta_dt->Write();
  nt_all_candidate_ionbeta_dt->Write();

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
