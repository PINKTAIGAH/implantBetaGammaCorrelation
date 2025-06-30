#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <cstdint>
#include <tuple>
#include <utility>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TCutG.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TString.h>
#include <TDirectory.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

// *************************************************************************************
// ****************************** DEFINE SCRIPT CONSTANTS ******************************
// *************************************************************************************

namespace constants
{
  const std::string ISOTOPE_TREE = "zrnbmo"; // Name suffix for gatedimplant tree & branch in anatree
  const int DSSD = 1;                      // Which DSSD will the analysis be run on

  const bool ONLY_OFFSPILL_DECAY = true;      // Check for onspill decay matches
  const std::pair<Long64_t, Long64_t> PROMPT_GAMMA_WINDOW = {-18200, -12900};

  const std::vector<int> BROKEN_AIDA_X_STRIPS_IMPLANT = {};
  const std::vector<int> BROKEN_AIDA_Y_STRIPS_IMPLANT = {};
  const std::vector<int> BROKEN_AIDA_X_STRIPS_DECAY = {63, 64, 66, 130, 189, 191, 194, 225, 256, 319, 320, 322, 383, 381, 192, 128, 255};
  const std::vector<int> BROKEN_AIDA_Y_STRIPS_DECAY = {};

  const std::pair<std::string, std::string> IMPLANT_GATES_INFILE_MAP = {"zrnbmo", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/new_gates/zrnbmo_new.root"};
}

namespace experimentInfo
{
  const uint64_t WR_EXPERIMENT_START = 1740207161207825830;                            // Run 26         
  const uint64_t WR_EXPERIMENT_END = 1740213634452253830;                              // Run 26     
  // const uint64_t WR_EXPERIMENT_START = 1740223603996689830;                               // Run 32   
  // const uint64_t WR_EXPERIMENT_END = 1740225294308847830;                                 // Run 32  
  const Double_t SLICES_EVERY = 0.05;                                                     // Size of white rabbit histogram bins
  const int64_t DURATION_IN_SECONDS = (WR_EXPERIMENT_END - WR_EXPERIMENT_START) / 1e9; // Duration of experiment
  const int64_t NUMBER_OF_SLICES = DURATION_IN_SECONDS / SLICES_EVERY;                 // Number of white rabbit histogram bins

  const int64_t SPILL_PERIOD = 10e9;
  const int64_t NUMBER_OF_SLICES_CUMMULATIVE = (SPILL_PERIOD/1e9) / SLICES_EVERY;                 // Number of white rabbit histogram bins
}


// *************************************************************************************
// ****************************** DEFINE FUNCTIONS *************************************
// *************************************************************************************

bool IsNoisyStrip(std::vector<int> noisyStripVector, std::pair<int, int> decayClusterMinMax){
  // Function which will loop over vector containing noisy strips and check
  // if current event occured in a noisy strip

  bool isNoisy = false; // Flag telling if event is from noisy strip

  // Loop over vector if not empty
  if (!noisyStripVector.empty()){
    // Loop for each noisy strip
    for (auto& noisyStrip : noisyStripVector){
      // If event strip corresponds to a noisy strip, change flag and break
      if ( noisyStrip >= decayClusterMinMax.first && noisyStrip < decayClusterMinMax.second+1 ){
        isNoisy = true;
        break;
      }
    }
  }

  return isNoisy;
}

TH1D* scaleHistogram(TH1D* hist, double xMin, double xMax, double scaleFactor) {
  
  // Make a deep copy of the original histogram
  TH1D* scaledHist =  (TH1D*)hist->Clone();
  scaledHist->SetName(Form("%s_scaled", hist->GetName()));

  // Get the bin numbers corresponding to the x-range
  int binMin = scaledHist->FindBin(xMin);
  int binMax = scaledHist->FindBin(xMax);

  // Loop over the specified bins and scale their content and errors
  for (int bin = binMin; bin <= binMax; ++bin) {
      double content = scaledHist->GetBinContent(bin);
      double error = scaledHist->GetBinError(bin);

      scaledHist->SetBinContent(bin, content * scaleFactor);
      scaledHist->SetBinError(bin, error * scaleFactor);
  }

  return scaledHist;
}



// *************************************************************************************
// ****************************** DEFINE STRUCTS ***************************************
// *************************************************************************************

struct ImplantEvent {
  ULong64_t time;
  ULong64_t timeX; 
  ULong64_t timeY; 
  Double_t x; 
  Double_t y; 
  Double_t e; 
  Double_t eX; 
  Double_t eY;
  Int_t csX;
  Int_t csY;
  Int_t cminX;
  Int_t cmaxX;
  Int_t cminY;
  Int_t cmaxY;
  Int_t dssd;
  Int_t spill;
};

struct DecayEvent {
  ULong64_t time;
  ULong64_t timeX; 
  ULong64_t timeY; 
  Double_t x; 
  Double_t y; 
  Double_t e; 
  Double_t eX; 
  Double_t eY;
  Int_t csX;
  Int_t csY;
  Int_t cminX;
  Int_t cmaxX;
  Int_t cminY;
  Int_t cmaxY;
  Int_t dssd;
  Int_t spill;
};

struct GermaniumEvent {
  ULong64_t time;
  Double_t e;
  Int_t spill;
};

struct FrsEvent {
  ULong64_t time;
  Double_t aoq;
  Double_t x2;
  Double_t x4;
  Double_t z41;
  Double_t z42;
  Double_t sci41E;
  Double_t sci42E;
  Int_t spill;
};


// *************************************************************************************
// ****************************** DEFINE MAPS FOR EVENTS *******************************
// *************************************************************************************


std::multimap<ULong64_t, ImplantEvent>    allImplantsEventMap;
std::multimap<ULong64_t, ImplantEvent>    gatedImplantsEventMap;
std::multimap<ULong64_t, DecayEvent>      decayEventMap;
std::multimap<ULong64_t, GermaniumEvent>  germaniumEventMap;
std::multimap<ULong64_t, FrsEvent>        frsEventMap;
std::multimap<ULong64_t, FrsEvent>        gatedFrsEventMap;


// *************************************************************************************
// *********************************** DEFINE HISTOGRAMS *******************************
// *************************************************************************************

TH1D* h1_frsRates = new TH1D("frs_rates", "Frs Rates", experimentInfo::NUMBER_OF_SLICES, 0, experimentInfo::WR_EXPERIMENT_END-experimentInfo::WR_EXPERIMENT_START);
TH1D* h1_gatedFrsRates = new TH1D("gated_frs_rates", "Gated Frs Rates", experimentInfo::NUMBER_OF_SLICES, 0, experimentInfo::WR_EXPERIMENT_END-experimentInfo::WR_EXPERIMENT_START);
TH1D* h1_implantRates = new TH1D("implant_rates", "Implant Rates", experimentInfo::NUMBER_OF_SLICES, 0, experimentInfo::WR_EXPERIMENT_END-experimentInfo::WR_EXPERIMENT_START);
TH1D* h1_gatedImplantRates = new TH1D("gated_implant_rates", "Gated Implant Rates", experimentInfo::NUMBER_OF_SLICES, 0, experimentInfo::WR_EXPERIMENT_END-experimentInfo::WR_EXPERIMENT_START);
TH1D* h1_decayRates = new TH1D("decay_rates", "Decay Rates", experimentInfo::NUMBER_OF_SLICES, 0, experimentInfo::WR_EXPERIMENT_END-experimentInfo::WR_EXPERIMENT_START);
TH1D* h1_germaniumDecayRates = new TH1D("germanium_decay_rates", "Germanium Decay Rates", experimentInfo::NUMBER_OF_SLICES, 0, experimentInfo::WR_EXPERIMENT_END-experimentInfo::WR_EXPERIMENT_START);

TH1D* h1_frsRatesCummulative = new TH1D("frs_rates_cummulative", "Frs Rates", experimentInfo::NUMBER_OF_SLICES_CUMMULATIVE, 0, experimentInfo::SPILL_PERIOD);
TH1D* h1_gatedFrsRatesCummulative = new TH1D("gated_frs_rates_cummulative", "Gated Frs Rates", experimentInfo::NUMBER_OF_SLICES_CUMMULATIVE, 0, experimentInfo::SPILL_PERIOD);
TH1D* h1_implantRatesCummulative = new TH1D("implant_rates_cummulative", "Implant Rates", experimentInfo::NUMBER_OF_SLICES_CUMMULATIVE, 0, experimentInfo::SPILL_PERIOD);
TH1D* h1_gatedImplantRatesCummulative = new TH1D("gated_implant_rates_cummulative", "Gated Implant Rates", experimentInfo::NUMBER_OF_SLICES_CUMMULATIVE, 0, experimentInfo::SPILL_PERIOD);
TH1D* h1_decayRatesCummulative = new TH1D("decay_rates_cummulative", "Decay Rates", experimentInfo::NUMBER_OF_SLICES_CUMMULATIVE, 0, experimentInfo::SPILL_PERIOD);
TH1D* h1_germaniumDecayRatesCummulative = new TH1D("germanium_decay_rates_cummulative", "Germanium Decay Rates", experimentInfo::NUMBER_OF_SLICES_CUMMULATIVE, 0, experimentInfo::SPILL_PERIOD);

TH2D* h1_decayHitPattern = new TH2D("decay_hit_pattern", "Decay Hit Pattern", 384, 0, 384, 128, 0, 128);
TH2D* h1_decayClusterSizes = new TH2D("decay_cluster_sizes", "Decay Cluster Sizes", 10, 0, 10, 10, 0, 10);
TH2D* h1_implantHitPattern = new TH2D("implant_hit_pattern", "Implant Hit Pattern", 384, 0, 384, 128, 0, 128);
TH2D* h1_implantClusterSizes = new TH2D("implant_cluster_sizes", "Implant Cluster Sizes", 10, 0, 10, 10, 0, 10);


// *************************************************************************************
// ****************************** START MACRO ******************************************
// *************************************************************************************

void makeDavidPlots(const char* input, const char* output){

  // *************************************************************************************
  // ****************************** LOAD IN FILES ****************************************
  // *************************************************************************************

  std::cout << "Started opening files..." << std::endl;

  // Load in input file
  TFile* file = TFile::Open(input);
  if (!file)
  {
    std::cerr << "Error: Could not open file " << std::endl;
    std::exit(1);
  }
  std::cout << "File loaded: " << file->GetName() << std::endl;

  // Get the trees
  TTree* implantTree = (TTree*)file->Get("aida_implant_tree");
  TTree* gatedImplantTree = (TTree*)file->Get((std::string("aida_gatedimplant_") + constants::ISOTOPE_TREE + std::string("_tree")).c_str());
  TTree* decayTree = (TTree*)file->Get("aida_decay_tree");
  TTree* germaniumTree = (TTree*)file->Get("germanium_tree");
  TTree* frsTree = (TTree*)file->Get("frs_tree");

  // Open output file
  TFile* outputFile = new TFile(output, "RECREATE");
  if (!outputFile)
  {
    std::cerr << "Error: Could not create output file " << std::endl;
    std::exit(1);
  }

  // Open gate file
  TFile* gateFile = new TFile(constants::IMPLANT_GATES_INFILE_MAP.second.c_str());
  if (!gateFile){
    std::cerr << "Error: Could not open gate file " << std::endl;
    std::exit(1);
  }

  // Extract TCutG objects from gates
  TCutG* x4aoqCut = (TCutG*)gateFile->Get("cut_x4_AoQ");
  TCutG* zz2Cut = (TCutG*)gateFile->Get("cut_Z_Z2");
  TCutG* zsci42ECut = (TCutG*)gateFile->Get("cut_Z_sci42E");

  std::cout << "Finished opening files!\n" << std::endl;

  // *************************************************************************************
  // ****************************** LOAD IN LEAVES ***************************************
  // *************************************************************************************

  std::cout << "Started reading trees..." << std::endl;

  // Create tree readers objects
  TTreeReader implantReader(implantTree);
  TTreeReader gatedImplantReader(gatedImplantTree);
  TTreeReader decayReader(decayTree);
  TTreeReader germaniumReader(germaniumTree);
  TTreeReader frsReader(frsTree);

  // Define leaves of variables for implant tree
  TTreeReaderValue<ULong64_t> implant_time(implantReader, "implant.time");
  TTreeReaderValue<ULong64_t> implant_time_x(implantReader, "implant.time_x");
  TTreeReaderValue<ULong64_t> implant_time_y(implantReader, "implant.time_y");
  TTreeReaderValue<Double_t> implant_x(implantReader, "implant.x");
  TTreeReaderValue<Double_t> implant_y(implantReader, "implant.y");
  TTreeReaderValue<Double_t> implant_energy(implantReader, "implant.e");
  TTreeReaderValue<Double_t> implant_energy_x(implantReader, "implant.ex");
  TTreeReaderValue<Double_t> implant_energy_y(implantReader, "implant.ey");
  TTreeReaderValue<Int_t> implant_clustersize_x(implantReader, "implant.csx");
  TTreeReaderValue<Int_t> implant_clustersize_y(implantReader, "implant.csy");
  TTreeReaderValue<Int_t> implant_clustermin_x(implantReader, "implant.cminx");
  TTreeReaderValue<Int_t> implant_clustermax_x(implantReader, "implant.cmaxx");
  TTreeReaderValue<Int_t> implant_clustermin_y(implantReader, "implant.cminy");
  TTreeReaderValue<Int_t> implant_clustermax_y(implantReader, "implant.cmaxy");
  TTreeReaderValue<Int_t> implant_dssd(implantReader, "implant.dssd");
  TTreeReaderValue<Int_t> implant_spill(implantReader, "implant.sp");

  // Define leaves of variables for gatedimplant tree
  TTreeReaderValue<ULong64_t> gatedImplant_time(gatedImplantReader,(std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".time")).c_str());
  TTreeReaderValue<ULong64_t> gatedImplant_time_x(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".time_x")).c_str());
  TTreeReaderValue<ULong64_t> gatedImplant_time_y(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".time_y")).c_str());
  TTreeReaderValue<Double_t> gatedImplant_x(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".x")).c_str());
  TTreeReaderValue<Double_t> gatedImplant_y(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".y")).c_str());
  TTreeReaderValue<Double_t> gatedImplant_energy(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".e")).c_str());
  TTreeReaderValue<Double_t> gatedImplant_energy_x(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".ex")).c_str());
  TTreeReaderValue<Double_t> gatedImplant_energy_y(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".ey")).c_str());
  TTreeReaderValue<Int_t> gatedImplant_clustersize_x(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".csx")).c_str());
  TTreeReaderValue<Int_t> gatedImplant_clustersize_y(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".csy")).c_str());
  TTreeReaderValue<Int_t> gatedImplant_clustermin_x(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".cminx")).c_str());
  TTreeReaderValue<Int_t> gatedImplant_clustermax_x(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".cmaxx")).c_str());
  TTreeReaderValue<Int_t> gatedImplant_clustermin_y(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".cminy")).c_str());
  TTreeReaderValue<Int_t> gatedImplant_clustermax_y(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".cmaxy")).c_str());
  TTreeReaderValue<Int_t> gatedImplant_dssd(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".dssd")).c_str());
  TTreeReaderValue<Int_t> gatedImplant_spill(gatedImplantReader, (std::string("gatedimplant_") + constants::ISOTOPE_TREE + std::string(".sp")).c_str());

  // Define leaves of variables for decay tree
  TTreeReaderValue<ULong64_t> decay_time(decayReader, "decay.time");
  TTreeReaderValue<ULong64_t> decay_time_x(decayReader, "decay.time_x");
  TTreeReaderValue<ULong64_t> decay_time_y(decayReader, "decay.time_y");
  TTreeReaderValue<Double_t> decay_x(decayReader, "decay.x");
  TTreeReaderValue<Double_t> decay_y(decayReader, "decay.y");
  TTreeReaderValue<Double_t> decay_energy(decayReader, "decay.e");
  TTreeReaderValue<Double_t> decay_energy_x(decayReader, "decay.ex");
  TTreeReaderValue<Double_t> decay_energy_y(decayReader, "decay.ey");
  TTreeReaderValue<Int_t> decay_clustersize_x(decayReader, "decay.csx");
  TTreeReaderValue<Int_t> decay_clustersize_y(decayReader, "decay.csy");
  TTreeReaderValue<Int_t> decay_clustermin_x(decayReader, "decay.cminx");
  TTreeReaderValue<Int_t> decay_clustermax_x(decayReader, "decay.cmaxx");
  TTreeReaderValue<Int_t> decay_clustermin_y(decayReader, "decay.cminy");
  TTreeReaderValue<Int_t> decay_clustermax_y(decayReader, "decay.cmaxy");
  TTreeReaderValue<Int_t> decay_dssd(decayReader, "decay.dssd");
  TTreeReaderValue<Int_t> decay_spill(decayReader, "decay.sp");

  // Define leaves of variables for germanium tree
  TTreeReaderValue<ULong64_t> germanium_time(germaniumReader, "germanium.time");
  TTreeReaderValue<Double_t> germanium_energy(germaniumReader, "germanium.energy");
  TTreeReaderValue<Int_t> germanium_spill(germaniumReader, "germanium.sp"); // sp = 1 spill, sp = 2 no spill

  // Define leaves of variables for frs tree
  TTreeReaderValue<ULong64_t> frs_time(frsReader, "frs.time");
  TTreeReaderValue<Double_t> frs_aoq(frsReader, "frs.aoq");
  TTreeReaderValue<Double_t> frs_x2(frsReader, "frs.x2");
  TTreeReaderValue<Double_t> frs_x4(frsReader, "frs.x4");
  TTreeReaderValue<Double_t> frs_z41(frsReader, "frs.z41");
  TTreeReaderValue<Double_t> frs_z42(frsReader, "frs.z42");
  TTreeReaderValue<Double_t> frs_sci41E(frsReader, "frs.sci41E");
  TTreeReaderValue<Double_t> frs_sci42E(frsReader, "frs.sci42E");
  TTreeReaderValue<Int_t> frs_spill(frsReader, "frs.sp");


  std::cout << "Finished reading trees!\n" << std::endl;

  // *************************************************************************************
  // ****************************** FILL MAPS WITH EVENTS ********************************
  // *************************************************************************************

  std::cout << "Started filling maps..." << std::endl;

  // Load in all implant maps
  while (implantReader.Next()){

    if ( *implant_dssd == constants::DSSD && TMath::Abs((Long64_t)*implant_time_x-(Long64_t)*implant_time_y) < 5e3 ){

      ImplantEvent tempEvent = {
        *implant_time, *implant_time_x, *implant_time_y, *implant_x, *implant_y, *implant_energy, *implant_energy_x, *implant_energy_y, *implant_clustersize_x, *implant_clustersize_y,
        *implant_clustermin_x, *implant_clustermax_x, *implant_clustermin_y, *implant_clustermax_y, *implant_dssd, *implant_spill 
      };
      allImplantsEventMap.emplace(*implant_time, tempEvent);

    }
  }

  // Load in gated implant maps
  while (gatedImplantReader.Next()){

    if ( *gatedImplant_dssd == constants::DSSD && TMath::Abs((Long64_t)*gatedImplant_time_x-(Long64_t)*gatedImplant_time_y) < 5e3 ){

      ImplantEvent tempEvent = {
        *gatedImplant_time, *gatedImplant_time_x, *gatedImplant_time_y, *gatedImplant_x, *gatedImplant_y, *gatedImplant_energy, *gatedImplant_energy_x, *gatedImplant_energy_y,
        *gatedImplant_clustersize_x, *gatedImplant_clustersize_y, *gatedImplant_clustermin_x, *gatedImplant_clustermax_x, *gatedImplant_clustermin_y, *gatedImplant_clustermax_y, 
        *gatedImplant_dssd, *gatedImplant_spill 
      };
      gatedImplantsEventMap.emplace(*gatedImplant_time, tempEvent);

    }
  }

  // Load in decay maps
  while (decayReader.Next()){

    if ( *decay_dssd == constants::DSSD && *decay_clustersize_x <= 3 && *decay_clustersize_y <= 3 && TMath::Abs((Long64_t)*decay_time_x-(Long64_t)*decay_time_y) < 5e3 && TMath::Abs(*decay_energy_x - *decay_energy_y) < 168. && *decay_energy > 150. && *decay_energy < 1000.){

      DecayEvent tempEvent = {
        *decay_time, *decay_time_x, *decay_time_y, *decay_x, *decay_y, *decay_energy, *decay_energy_x, *decay_energy_y, *decay_clustersize_x, *decay_clustersize_y,
        *decay_clustermin_x, *decay_clustermax_x, *decay_clustermin_y, *decay_clustermax_y,*decay_dssd, *decay_spill 
      };
      decayEventMap.emplace(*decay_time, tempEvent);

    }
  }

  // Load in germanium maps
  while (germaniumReader.Next()){

    if ( *germanium_energy > 20. ){
      GermaniumEvent tempEvent = { *germanium_time, *germanium_energy, *germanium_spill };
      germaniumEventMap.emplace(*germanium_time, tempEvent);

    }
  }

  ULong64_t lastFrsTime = 0;
  // Load in frs maps
  while (frsReader.Next()){

    FrsEvent tempEvent = { *frs_time, *frs_aoq, *frs_x2, *frs_x4, *frs_z41, *frs_z42, *frs_sci41E, *frs_sci42E, *frs_spill };
    if (*frs_time != lastFrsTime) frsEventMap.emplace(*frs_time, tempEvent);

    if ( x4aoqCut->IsInside(*frs_aoq, *frs_x4) && zz2Cut->IsInside(*frs_z41, *frs_z42) && zsci42ECut->IsInside(*frs_sci42E, *frs_z41) ){
      gatedFrsEventMap.emplace(*frs_time, tempEvent);
    }

    lastFrsTime = *frs_time;
  }

  std::cout << "Finished filling maps...\n" << std::endl;


  // *************************************************************************************
  // ****************************** LOOP OVER FRS EVENTS *********************************
  // *************************************************************************************

  bool startCummulative;
  int lastSpill;
  Long64_t lastOnspillTime;
  Long64_t lastTime;
  Double_t offspillThreshold = 0.5e9;
  

  std::cout << "Started filling histograms..." << std::endl;

  startCummulative = false;
  lastOnspillTime = 0;
  lastTime = experimentInfo::WR_EXPERIMENT_START;
  // Loop over frs events
  for (auto& itr : frsEventMap){

    FrsEvent currentFrsEvent = itr.second; // Unpack current event
    h1_frsRates->Fill(currentFrsEvent.time-experimentInfo::WR_EXPERIMENT_START); // Plot time

    if ( (Long64_t)currentFrsEvent.time - lastTime > offspillThreshold ){
      startCummulative = true;
      lastOnspillTime = (Long64_t)currentFrsEvent.time;
    }
    if (startCummulative) h1_frsRatesCummulative->Fill( (Long64_t)currentFrsEvent.time - lastOnspillTime );
    lastTime = currentFrsEvent.time;

  }


  startCummulative = false;
  lastOnspillTime = 0;
  lastTime = experimentInfo::WR_EXPERIMENT_START;
  // Loop over gated frs events
  for (auto& itr : gatedFrsEventMap){

    FrsEvent currentGatedFrsEvent = itr.second; // Unpack current event
    h1_gatedFrsRates->Fill(currentGatedFrsEvent.time-experimentInfo::WR_EXPERIMENT_START); // Plot time

    if ((Long64_t)currentGatedFrsEvent.time - lastTime > offspillThreshold){
      startCummulative = true;
      lastOnspillTime = (Long64_t)currentGatedFrsEvent.time;
    }
    if (startCummulative) h1_gatedFrsRatesCummulative->Fill( (Long64_t)currentGatedFrsEvent.time - lastOnspillTime );
    lastTime = currentGatedFrsEvent.time;

  }


  // *************************************************************************************
  // ****************************** LOOP OVER IMPLANT EVENTS *****************************
  // *************************************************************************************

  startCummulative = false;
  lastOnspillTime = 0;
  lastTime = experimentInfo::WR_EXPERIMENT_START;
  // Loop over all implants
  for (auto& itr : allImplantsEventMap){

    ImplantEvent currentImplantEvent = itr.second; // Unpack current event
    h1_implantRates->Fill(currentImplantEvent.time-experimentInfo::WR_EXPERIMENT_START); // Plot time

    if ( (Long64_t)currentImplantEvent.time - lastTime > offspillThreshold ){
      startCummulative = true;
      lastOnspillTime = (Long64_t)currentImplantEvent.time;
    }
    if (startCummulative) h1_implantRatesCummulative->Fill( (Long64_t)currentImplantEvent.time - lastOnspillTime );
    lastTime = currentImplantEvent.time;

  }

  startCummulative = false;
  lastOnspillTime = 0;
  lastTime = experimentInfo::WR_EXPERIMENT_START;
  // Loop over gated implants
  for (auto& itr : gatedImplantsEventMap){

    ImplantEvent currentGatedImplantEvent = itr.second; // Unpack current event
    h1_gatedImplantRates->Fill(currentGatedImplantEvent.time-experimentInfo::WR_EXPERIMENT_START); // Plot time
    h1_implantHitPattern->Fill(currentGatedImplantEvent.x, currentGatedImplantEvent.y);
    h1_implantClusterSizes->Fill(currentGatedImplantEvent.csX, currentGatedImplantEvent.csY);

    if ( (Long64_t)currentGatedImplantEvent.time - lastTime > offspillThreshold ){
      startCummulative = true;
      lastOnspillTime = (Long64_t)currentGatedImplantEvent.time;
    }
    if (startCummulative) h1_gatedImplantRatesCummulative->Fill( (Long64_t)currentGatedImplantEvent.time - lastOnspillTime );
    lastTime = currentGatedImplantEvent.time;

  }


  // *************************************************************************************
  // ******************************** LOOP OVER DECAY EVENTS *****************************
  // *************************************************************************************

  // Loop over decay events
  int counter = 0;
  int totalEntries = decayEventMap.size();
  startCummulative = false;
  lastSpill = 1;
  lastOnspillTime = 0;

  for (auto& itr : decayEventMap){

    if ( ++counter%100000 == 0 ) std::cout << "\rProgress: " << counter << " / " << totalEntries << " (" <<  counter*100./totalEntries << " %)" << std::flush;

    DecayEvent currentDecayEvent = itr.second; // Unpack current event

    if ( IsNoisyStrip(constants::BROKEN_AIDA_X_STRIPS_DECAY, std::make_pair(currentDecayEvent.cminX, currentDecayEvent.cmaxX)) ) continue;
    if ( IsNoisyStrip(constants::BROKEN_AIDA_Y_STRIPS_DECAY, std::make_pair(currentDecayEvent.cminY, currentDecayEvent.cmaxY)) ) continue;
    
    h1_decayHitPattern->Fill(currentDecayEvent.x, currentDecayEvent.y);
    h1_decayClusterSizes->Fill(currentDecayEvent.csX, currentDecayEvent.csY);

    if (lastSpill == 2 && currentDecayEvent.spill == 1){
      startCummulative = true;
      lastOnspillTime = (Long64_t)currentDecayEvent.time;
    }
    lastSpill = currentDecayEvent.spill;

    if (currentDecayEvent.spill == 1 && constants::ONLY_OFFSPILL_DECAY) continue;

    h1_decayRates->Fill(currentDecayEvent.time-experimentInfo::WR_EXPERIMENT_START); // Plot time

    if (startCummulative) h1_decayRatesCummulative->Fill( (Long64_t)currentDecayEvent.time - lastOnspillTime );

    auto lowerBound = germaniumEventMap.lower_bound(currentDecayEvent.time - 20000); 
    auto upperBound = germaniumEventMap.lower_bound(currentDecayEvent.time - 5000); 

    for (auto germaniumItr=lowerBound; germaniumItr!=upperBound; ++germaniumItr){

      GermaniumEvent currentGermaniumEvent = germaniumItr->second;
      Long64_t timeDiff = (Long64_t)currentGermaniumEvent.time - (Long64_t)currentDecayEvent.time;

      // // Check if within prompt window
      // if (timeDiff < constants::PROMPT_GAMMA_WINDOW.first) continue;
      // if (timeDiff > constants::PROMPT_GAMMA_WINDOW.second + 20e3) break;

      h1_germaniumDecayRates->Fill(currentDecayEvent.time-experimentInfo::WR_EXPERIMENT_START); // Fill Plot
      if (startCummulative) h1_germaniumDecayRatesCummulative->Fill( (Long64_t)currentDecayEvent.time - lastOnspillTime );


    }


  }

  std::cout << "\nFinished filling histograms!\n" << std::endl;


  // *************************************************************************************
  // ********************************** SCALE DECAY HISTOGRAMS ***************************
  // *************************************************************************************

  TH1D* h1_decayRatesCummulativeScaled = scaleHistogram(h1_decayRatesCummulative, 0e9, 3.5e9, 1./7.);
  TH1D* h1_germaniumDecayRatesCummulativeScaled = scaleHistogram(h1_germaniumDecayRatesCummulative, 0e9, 3.5e9, 1./7.);

  // *************************************************************************************
  // ****************************** WRITE HISTOGRAMS TO OUTPUT FILE **********************
  // *************************************************************************************

  std::cout << "Writing histograms to file..." << std::endl;

  outputFile->cd();

  TDirectory* ratesDir = outputFile->mkdir("rates");
  ratesDir->cd();
  h1_frsRates->Write();
  h1_gatedFrsRates->Write();
  h1_implantRates->Write();
  h1_gatedImplantRates->Write();
  h1_decayRates->Write();
  h1_germaniumDecayRates->Write();

  TDirectory* cummulativeRatesDir = outputFile->mkdir("cummulative_rates");
  cummulativeRatesDir->cd();
  h1_frsRatesCummulative->Write();
  h1_gatedFrsRatesCummulative->Write();
  h1_implantRatesCummulative->Write();
  h1_gatedImplantRatesCummulative->Write();
  h1_decayRatesCummulative->Write();
  h1_decayRatesCummulativeScaled->Write();
  h1_germaniumDecayRatesCummulative->Write();
  h1_germaniumDecayRatesCummulativeScaled->Write();

  TDirectory* hitsDir = outputFile->mkdir("hits");
  hitsDir->cd();
  h1_decayHitPattern->Write();
  h1_decayClusterSizes->Write();
  h1_implantHitPattern->Write();
  h1_implantClusterSizes->Write();

  std::cout << "Finished writing histograms to file!\n" << std::endl;

  // *************************************************************************************
  // ****************************** CLEANUP **********************************************
  // *************************************************************************************

  // Close files
  file->Close();
  outputFile->Close();
  gateFile->Close();

  // Clear heap
  delete file;
  delete outputFile;
  delete gateFile;
  delete h1_frsRates;
  delete h1_gatedFrsRates;
  delete h1_implantRates;
  delete h1_gatedImplantRates;
  delete h1_decayRates;
  delete h1_germaniumDecayRates;
  delete h1_frsRatesCummulative;
  delete h1_gatedFrsRatesCummulative;
  delete h1_implantRatesCummulative;
  delete h1_gatedImplantRatesCummulative;
  delete h1_decayRatesCummulative;
  delete h1_germaniumDecayRatesCummulative;
  delete h1_decayHitPattern;
  delete h1_decayClusterSizes;
  delete h1_implantHitPattern;
  delete h1_implantClusterSizes;

  std::cout << "Succesfully ran file!\n" << std::endl;
  std::exit(0);

}