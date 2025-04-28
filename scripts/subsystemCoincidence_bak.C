#include<map>
#include<set>
#include<iostream>
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

  const int DSSD = 1; // Which DSSD will the analysis be run on

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

struct WrSubsystemTime{
  // Struct to hold times for subsystems with only WR times (AIDA, BB7 (MIDAS), FRS)

  ULong64_t time_wr;

  bool operator < (const WrSubsystemTime &other) const { return time_wr < other.time_wr; }
};

struct FebexSubsystemTime{
  // Struct to hold times for TAMEX/FEBEX subsystems with a corrected time (bPlast, Germanium, FATIMA, BGO, BB7 (Mezitek), etc)

  ULong64_t time_wr;
  ULong64_t time_abs_evt;

  // Overloaded < operator that is based on the !!!WHITE RABBIT TIME!!!
  bool operator < (const FebexSubsystemTime &other) const { return time_wr < other.time_wr; }
};

// *************************************************************************************
// ****************************** DEFINE FUNCTIONS *************************************
// *************************************************************************************


// *************************************************************************************
// ****************************** DEFINE CLASSES ***************************************
// *************************************************************************************


// *************************************************************************************
// ****************************** DEFINE SETS FOR EVENTS *******************************
// *************************************************************************************

// Sets to hold times from wr anatrees
std::set<WrSubsystemTime> implantTimeSet;
std::set<WrSubsystemTime> decayTimeSet;
std::set<WrSubsystemTime> frsTimeSet;
std::set<FebexSubsystemTime> bplastTimeSet;
std::set<FebexSubsystemTime> germaniumTimeSet;

// *************************************************************************************
// ****************************** START MACRO ******************************************
// *************************************************************************************

void subsystemCoincidence_bak(const char* input, const char* output){

  // *************************************************************************************
  // ****************************** OPEN & CREATE FILES **********************************
  // *************************************************************************************
  
  // Load in input file 
  TFile* file = TFile::Open(input);
  if (!file){
    std::cerr << "Error: Could not open file " << std::endl;
    std::exit(1);
  }
  std::cout << "File loaded: "<< file->GetName() << std::endl << std::endl;

  // *************************************************************************************
  // ****************************** OPEN TREES *******************************************
  // *************************************************************************************
  
  // Get the trees
  TTree* implantTimeTree = (TTree*)file->Get("aida_implant_tree");
  TTree* decayTimeTree = (TTree*)file->Get("aida_decay_tree");
  TTree* frsTimeTree = (TTree*)file->Get("frs_tree");
  TTree* bplastTimeTree = (TTree*)file->Get("bplast_tree");
  TTree* germaniumTimeTree = (TTree*)file->Get("germanium_tree");
  
  // *************************************************************************************
  // ****************************** DEFINE OUTPUT FILE ***********************************
  // *************************************************************************************

  // Open output file
  TFile* outputFile = new TFile(output, "RECREATE");
  if (!outputFile){
    std::cerr << "Error: Could not create output file " << std::endl;
    std::exit(1);
  }

  // *************************************************************************************
  // ****************************** DEFINE TREEREADER ************************************
  // *************************************************************************************

  // Create tree readers objects
  TTreeReader implantReader(implantTimeTree);
  TTreeReader decayReader(decayTimeTree);
  TTreeReader frsReader(frsTimeTree);
  TTreeReader bplastReader(bplastTimeTree);
  TTreeReader germaniumReader(germaniumTimeTree);

  // Define leaves of variables for implant tree
  TTreeReaderValue<ULong64_t> implantWrTime(implantReader, "implant.time_wr");

  // Define leaves of variables for decay tree
  TTreeReaderValue<ULong64_t> decayWrTime(decayReader, "decay.time_wr");
  TTreeReaderValue<ULong64_t> decay_time_x(decayReader, "decay.time_x");
  TTreeReaderValue<ULong64_t> decay_time_y(decayReader, "decay.time_y");
  TTreeReaderValue<double> decay_e(decayReader, "decay.e");
  TTreeReaderValue<double> decay_ex(decayReader, "decay.ex");
  TTreeReaderValue<double> decay_ey(decayReader, "decay.ey");

  // Define leaves of variables for frs tree
  TTreeReaderValue<ULong64_t> frsWrTime(frsReader, "frs.time_wr");

  // Define leaves of variables for bplast tree
  TTreeReaderValue<ULong64_t> bplastWrTime(bplastReader, "bplast.time_wr");
  TTreeReaderValue<ULong64_t> bplastAbsEvtTime(bplastReader, "bplast.time_abs_evt");

  // Define leaves of variables for germanium tree
  TTreeReaderValue<ULong64_t> germaniumWrTime(germaniumReader, "germanium.time_wr");
  TTreeReaderValue<ULong64_t> germaniumAbsEvtTime(germaniumReader, "germanium.time_abs_evt");


  // *************************************************************************************
  // ****************************** DEFINE HISTOGRAMS ************************************
  // *************************************************************************************

  // Time spectra
  TH1F* h1ImplantWrTimeSpectrum = new TH1F("h1_implant_wrtime_spectrum", "AIDA Implant WR Time Spectrum;time;Counts", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1DecayWrTimeSpectrum = new TH1F("h1_decay_wrtime_spectrum", "AIDA Decay WR Time Spectrum;time;Counts", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1FrsWrTimeSpectrum = new TH1F("h1_frs_wrtime_spectrum", "AIDA FRS WR Time Spectrum;time;Counts", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1BplastWrTimeSpectrum = new TH1F("h1_bplast_wrtime_spectrum", "AIDA Bplast WR Time Spectrum;time;Counts", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1BplastAbsEvtTimeSpectrum = new TH1F("h1_bplast_absevttime_spectrum", "AIDA Bplast Abs Evt Time Spectrum;time;Counts", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1GermaniumWrTimeSpectrum = new TH1F("h1_germanium_wrtime_spectrum", "AIDA Germanium WR Time Spectrum;time;Counts", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);
  TH1F* h1GermaniumAbsEvtTimeSpectrum = new TH1F("h1_germanium_absevttime_spectrum", "AIDA Germanium Abs Evt Time Spectrum;time;Counts", experimentInfo::NUMBER_OF_SLICES, experimentInfo::WR_EXPERIMENT_START, experimentInfo::WR_EXPERIMENT_END);

  // Implant-Decay Coincidence
  TH1F* h1ImplantDecayCoincidenceTime = new TH1F("h1_implant_decay_coincidence", "AIDA Implant Decay Time Spectrum;dt (ns); Counts", 500/20, 0, 500e3);

  // Implant-FRS Coincidence
  TH1F* h1ImplantFrsCoincidenceTime = new TH1F("h1_implant_frs_coincidence", "AIDA Implant FRS Time Spectrum;dt (ns); Counts", 5000, 0, 10e9);
  
  // Implant-Bplast Coinicdence 
  TH1F* h1ImplantBplastCoincidenceWrTime = new TH1F("h1_implant_bplast_wrtime_coincidence", "AIDA Decay bPlast Wr Time Spectrum;dt (ns); Counts", 5000, 0, 100e6);
  TH1F* h1ImplantBplastCoincidenceAbsEvtTime = new TH1F("h1_implant_bplast_absevttime_coincidence", "AIDA Decay Bplast Abs Evt Time Spectrum;dt (ns); Counts", 5000, 0, 100e6);
  
  // Decay-Germanium Coincidence
  TH1F* h1DecayGermaniumCoincidenceWrTime = new TH1F("h1_decay_germanium_wrtime_coincidence", "AIDA Decay Germanium Wr Time Spectrum;dt (ns); Counts", 5000, 0, 100e6);
  TH1F* h1DecayGermaniumCoincidenceAbsEvtTime = new TH1F("h1_decay_germanium_absevttime_coincidence", "AIDA Decay Germanium Abs Evt Time Spectrum;dt (ns); Counts", 5000, 0, 10e9);

  // Decay-Bplast Coinicdence 
  TH1F* h1DecayBplastCoincidenceWrTime = new TH1F("h1_decay_bplast_wrtime_coincidence", "AIDA Decay bPlast Wr Time Spectrum;dt (ns); Counts", 5000, 0, 500e3);
  TH1F* h1DecayBplastCoincidenceAbsEvtTime = new TH1F("h1_decay_bplast_absevttime_coincidence", "AIDA Decay Bplast Abs Evt Time Spectrum;dt (ns); Counts", 5000, 0, 500e3);

  // Wr-AbsEvt Time diff 
  TH1F* h1BplastWrAbsEvtTimeDiff = new TH1F("h1_bplast_wr_absevt_timediff", "bPlast Wr-AbsEvt Time Difference;dt (ns); Counts", 5000, 0, 400);
  TH1F* h1GermaniumWrAbsEvtTimeDiff = new TH1F("h1_germanium_wr_absevt_timediff", "Germanium Wr-AbsEvt Time Difference;dt (ns); Counts", 5000, 5e9, 20e9);

  // *************************************************************************************
  // ****************************** FILL MAPS WITH EVENTS ********************************
  // *************************************************************************************
  
  // Define counter to fee number of loops
  int counter = 0; 

  // Read implant events
  while (implantReader.Next()){
    // Create and fill struct
    WrSubsystemTime timeStruct = {*implantWrTime};
    implantTimeSet.emplace(timeStruct);  
    counter++;
  }
  std::cout << "Finished filling the implant time set!" << std::endl;
  std::cout << "Number of implant time entries: " << implantTimeSet.size() << " / " << counter << std::endl << std::endl;

  // Read decay events
  counter = 0;
  while (decayReader.Next()){
    // Create and fill struct
    if( TMath::Abs( (int64_t)(*decay_time_x-*decay_time_y) )<5e3 && TMath::Abs(*decay_ex-*decay_ey)<168 && *decay_e>151 && *decay_e<1000 ){
      WrSubsystemTime timeStruct = {*decayWrTime};
      decayTimeSet.emplace(timeStruct);  
      counter++;
    }
  }
  std::cout << "Finished filling the decay time set!" << std::endl;
  std::cout << "Number of decay time entries: " << decayTimeSet.size() << " / " << counter << std::endl << std::endl;

  // Read frs events
  counter = 0;
  while (frsReader.Next()){
    // Create and fill struct
    WrSubsystemTime timeStruct = {*frsWrTime};
    frsTimeSet.emplace(timeStruct);  
    counter++;
  }
  std::cout << "Finished filling the frs time set!" << std::endl;
  std::cout << "Number of frs time entries: " << frsTimeSet.size() << " / " << counter  << std::endl << std::endl;

  // Read bplast events 
  counter = 0;
  while (bplastReader.Next()){
    // Create and fill struct
    FebexSubsystemTime timeStruct = {*bplastWrTime, *bplastAbsEvtTime};
    bplastTimeSet.emplace(timeStruct);  
    counter++;
  }
  std::cout << "Finished filling the bplast time set!" << std::endl;
  std::cout << "Number of bplast time entries: " << bplastTimeSet.size() << " / " << counter << std::endl << std::endl;

  // Read germanium events
  counter = 0;
  while (germaniumReader.Next()){
    // Create and fill struct
    FebexSubsystemTime timeStruct = {*germaniumWrTime, *germaniumAbsEvtTime};
    germaniumTimeSet.emplace(timeStruct);  
    counter++;
  }
  std::cout << "Finished filling the germanium time set!" << std::endl;
  std::cout << "Number of germanium time entries: " << germaniumTimeSet.size() << " / " << counter << std::endl << std::endl;



  // *************************************************************************************
  // ****************************** IMPLANT COINCIDENCES *********************************
  // *************************************************************************************

  std::cout << "Started implant coincidences..." << std::endl;
  for ( const auto implantTimeStruct : implantTimeSet ){

    // Fill the time spectra
    h1ImplantWrTimeSpectrum->Fill(implantTimeStruct.time_wr);

    // MAKE IMPLANT DECAY COINCIDENCES
    for (const auto decayTimeStruct : decayTimeSet ){

      // Fill coincidence time spectra
      h1ImplantDecayCoincidenceTime->Fill( decayTimeStruct.time_wr - implantTimeStruct.time_wr );

    }

    // MAKE IMPLANT FRS COINCIDENCES
    for (const auto frsTimeStruct : frsTimeSet ){

      // Fill coincidence time spectra
      h1ImplantFrsCoincidenceTime->Fill( frsTimeStruct.time_wr - implantTimeStruct.time_wr );

    }

    // MAKE IMPLANT BPLAST COINCIDENCES
    for (const auto bplastTimeStruct : bplastTimeSet ){

      // Fill coincidence time spectra
      h1ImplantBplastCoincidenceWrTime->Fill( bplastTimeStruct.time_wr - implantTimeStruct.time_wr );
      h1ImplantBplastCoincidenceAbsEvtTime->Fill( bplastTimeStruct.time_abs_evt - implantTimeStruct.time_wr );

    }

  }
  std::cout << "Finished implant coincidences!" << std::endl << std::endl;


  // *************************************************************************************
  // ****************************** DECAY COINCIDENCES *********************************
  // *************************************************************************************

// std::cout << "Started decay coincidences..." << std::endl;
// for ( const auto decayTimeStruct : decayTimeSet ){

//   // Fill the time spectra
//   h1DecayWrTimeSpectrum->Fill(decayTimeStruct.time_wr);

//   // MAKE DECAY BPLAST COINCIDENCES
//   for (const auto bplastTimeStruct : bplastTimeSet ){

//     // Fill coincidence time spectra
//     h1DecayBplastCoincidenceWrTime->Fill( bplastTimeStruct.time_wr - decayTimeStruct.time_wr );
//     h1DecayBplastCoincidenceAbsEvtTime->Fill( bplastTimeStruct.time_abs_evt - decayTimeStruct.time_wr );

//   }

//   // MAKE DECAY GERMANIUM COINCIDENCES
//   for (const auto germaniumTimeStruct : germaniumTimeSet ){

//     // Fill coincidence time spectra
//     h1DecayGermaniumCoincidenceWrTime->Fill( germaniumTimeStruct.time_wr - decayTimeStruct.time_wr );
//     h1DecayGermaniumCoincidenceAbsEvtTime->Fill( germaniumTimeStruct.time_abs_evt - decayTimeStruct.time_wr );

//   }

// }
// std::cout << "Finished decay coincidences!" << std::endl << std::endl;


  // *************************************************************************************
  // ****************************** FRS COINCIDENCES *********************************
  // *************************************************************************************

  std::cout << "Started frs coincidences..." << std::endl;
  for ( const auto frsTimeStruct : frsTimeSet ){

    // Fill the time spectra
    h1FrsWrTimeSpectrum->Fill(frsTimeStruct.time_wr);

  }
  std::cout << "Finished frs coincidences!" << std::endl << std::endl;


  // *************************************************************************************
  // ****************************** BPLAST COINCIDENCES *********************************
  // *************************************************************************************

  std::cout << "Started bplast coincidences..." << std::endl;
  for ( const auto bplastTimeStruct : bplastTimeSet ){

    // Fill the time spectra
    h1BplastWrTimeSpectrum->Fill(bplastTimeStruct.time_wr);
    h1BplastAbsEvtTimeSpectrum->Fill(bplastTimeStruct.time_abs_evt);

    // Fill the wr - absevt time difference spectra
    h1BplastWrAbsEvtTimeDiff->Fill( bplastTimeStruct.time_wr - bplastTimeStruct.time_abs_evt );

  }
  std::cout << "Finished bplast coincidences!" << std::endl << std::endl;


  // *************************************************************************************
  // ****************************** GERMANIUM COINCIDENCES *********************************
  // *************************************************************************************

  std::cout << "Started germanium coincidences..." << std::endl;
  for ( const auto germaniumTimeStruct : germaniumTimeSet ){

    // Fill the time spectra
    h1GermaniumWrTimeSpectrum->Fill(germaniumTimeStruct.time_wr);
    h1GermaniumAbsEvtTimeSpectrum->Fill(germaniumTimeStruct.time_abs_evt);

    // Fill the wr - absevt time difference spectra
    h1GermaniumWrAbsEvtTimeDiff->Fill( germaniumTimeStruct.time_wr - germaniumTimeStruct.time_abs_evt );

  }
  std::cout << "Finished germanium coincidences!" << std::endl << std::endl;

  
  // *************************************************************************************
  // ****************************** WRITE HISTOGRAMS *************************************
  // *************************************************************************************

  outputFile->cd();

  // Make directory obects
  TDirectory* timeSpectraDir = outputFile->mkdir("subsystem_time_spectra");
  TDirectory* implantSubsystemCoincidenceDir = outputFile->mkdir("implant_subsystem_coincidence");
  TDirectory* decaySubsystemCoincidenceDir = outputFile->mkdir("decay_subsystem_coincidence");
  TDirectory* wrAbsEvtTimeDifferenceDir = outputFile->mkdir("wr_absevt_time_difference");

  // Write timespectra
  timeSpectraDir->cd();
  h1ImplantWrTimeSpectrum->Write();
  h1DecayWrTimeSpectrum->Write();
  h1FrsWrTimeSpectrum->Write();
  h1BplastWrTimeSpectrum->Write();
  h1BplastAbsEvtTimeSpectrum->Write();
  h1GermaniumWrTimeSpectrum->Write();
  h1GermaniumAbsEvtTimeSpectrum->Write();
  gFile->cd();


  // Implant-Decay Coincidence
  implantSubsystemCoincidenceDir->cd();
  h1ImplantDecayCoincidenceTime->Write();

  // Implant-FRS Coincidence
  h1ImplantFrsCoincidenceTime->Write();

  //  Implant-Bplast Coincidence
  h1ImplantBplastCoincidenceWrTime->Write();
  h1ImplantBplastCoincidenceAbsEvtTime->Write();
  gFile->cd();
  
  // Decay-Germanium Coincidence
  decaySubsystemCoincidenceDir->cd();
  h1DecayGermaniumCoincidenceWrTime->Write();
  h1DecayGermaniumCoincidenceAbsEvtTime->Write();

  // Decay-Bplast Coinicdence 
  h1DecayBplastCoincidenceWrTime->Write();
  h1DecayBplastCoincidenceAbsEvtTime->Write();
  gFile->cd();

  // Wr-AbsEvt Time diff 
  wrAbsEvtTimeDifferenceDir->cd();
  h1BplastWrAbsEvtTimeDiff->Write();
  h1GermaniumWrAbsEvtTimeDiff->Write();
  gFile->cd();

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
