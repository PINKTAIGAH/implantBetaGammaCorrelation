#include <iostream>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

double getPercentage(int a, int b){
  // Get Percentage as a double when both parameters are integers
  return 100*(double)a/b;
}

void findAidaCounts(const char* input){

  TChain* tree = new TChain("evt");
  tree->Add(input);

  // Create a tree reader
  TTreeReader reader(tree);

  // Create tree reader values for the branches you want to read
  TTreeReaderValue<bool> spill(reader, "EventHeader.fSpillFlag");

  TTreeReaderArray<int64_t> implant_time(reader, "AidaImplantHits.Time");
  TTreeReaderArray<Int_t> implant_dssd(reader, "AidaImplantHits.DSSD");

  TTreeReaderArray<int64_t> decay_time(reader, "AidaDecayHits.Time");
  TTreeReaderArray<Int_t> decay_dssd(reader, "AidaDecayHits.DSSD");

  const char spinner[] = {'-', '\\', '|', '/'};
  int totalEntries = reader.GetEntries(true);

  // Initialise counters
  std::map< std::string, int > implant_event_counter_map = { {"all", 0}, {"onspill", 0}, {"offspill", 0} };
  std::map< std::string, int > decay_event_counter_map = { {"all", 0}, {"onspill", 0}, {"offspill", 0} };
  std::map< std::string, int > spillflag_counter_map = { {"all", 0}, {"onspill", 0}, {"offspill", 0} };

  // Loop over all the entries in the new tree
  while (reader.Next()){

    // Get subevent mulltiplicities
    int implant_hits = implant_time.GetSize();
    int decay_hits = decay_time.GetSize();

    // Log spillflags
    spillflag_counter_map["all"]++;
    if (*spill) {spillflag_counter_map["onspill"]++;}
    if (!*spill) {spillflag_counter_map["offspill"]++;}

    // Loop over implant subevents
    if (implant_hits>0){
      for (int j=0; j<implant_hits; j++){
        if (implant_dssd[j]==1){
          implant_event_counter_map["all"]++;
          if (*spill) {implant_event_counter_map["onspill"]++;}
          if (!*spill) {implant_event_counter_map["offspill"]++;}
        }
      }
    }

    // Loop over decay subevents
    if (decay_hits>0){
      for (int j=0; j<decay_hits; j++){
        if (decay_dssd[j]==1){
          decay_event_counter_map["all"]++;
          if (*spill) {decay_event_counter_map["onspill"]++;}
          if (!*spill) {decay_event_counter_map["offspill"]++;}
        }
      }
    }

    // Show the progress of the loop
    if (reader.GetCurrentEntry() % 100000 == 0) {
        int progress = (reader.GetCurrentEntry() * 100) / totalEntries;
        char spin = spinner[reader.GetCurrentEntry() / 1000 % 4];
        std::cout << "\rProcessing the tree " << reader.GetCurrentEntry() << " (" << progress << "%) " << spin << std::flush;
    }
  }

  std::cout << std::endl << std::endl;

  std::cout << Form("Total number of events: %lld", reader.GetEntries()) << std::endl << std::endl;

  std::cout << "Number of Implant Hits: " << implant_event_counter_map["all"] << "  ###   Onspill: " << implant_event_counter_map["onspill"] << " (" << Form( "%.1f", getPercentage(implant_event_counter_map["onspill"], implant_event_counter_map["all"]) ) << "%)   ###   Offspill: " << implant_event_counter_map["offspill"] << " (" << Form( "%.1f", getPercentage(implant_event_counter_map["offspill"], implant_event_counter_map["all"]) ) << "%)" << std::endl;

  std::cout << "Number of Decay Hits: " << decay_event_counter_map["all"] << "  ###   Onspill: " << decay_event_counter_map["onspill"] << " (" << Form( "%.1f", getPercentage(decay_event_counter_map["onspill"], decay_event_counter_map["all"]) ) << "%)   ###   Offspill: " << decay_event_counter_map["offspill"] << " (" << Form( "%.1f", getPercentage(decay_event_counter_map["offspill"], decay_event_counter_map["all"]) ) << "%)" << std::endl;
  
  std::cout << "Number of Spill Flags: " << spillflag_counter_map["all"] << "  ###   Onspill: " << spillflag_counter_map["onspill"] << " (" << Form( "%.1f", getPercentage(spillflag_counter_map["onspill"], spillflag_counter_map["all"]) ) << "%)   ###   Offspill: " << spillflag_counter_map["offspill"] << " (" << Form( "%.1f", getPercentage(spillflag_counter_map["offspill"], spillflag_counter_map["all"]) ) << "%)" << std::endl << std::endl;

  std::cout << "ROOT Macro has run succesfully!" << std::endl << std::endl;

  // Close the file
  std::exit(0);

  return;
}
