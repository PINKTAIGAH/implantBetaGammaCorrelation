#include <iostream>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

void findWrTimes(const char* input){

  TChain* tree = new TChain("evt");
  tree->Add(input);

  // Create a tree reader
  TTreeReader reader(tree);

  // Create tree reader values for the branches you want to read
  TTreeReaderArray<Long_t> frs_time(reader, "FrsHitData.fwr_t");

  // Define the counter
  int counter = 0;
  ULong_t wr_start = 0;
  ULong_t wr_end = 0;

  const char spinner[] = {'-', '\\', '|', '/'};
  int totalEntries = reader.GetEntries(true);

  // Loop over all the entries in the new tree
  while (reader.Next()){

    int frshits = frs_time.GetSize();

    if (frshits>0 && wr_start==0 ){ wr_start = frs_time[0]; }
    if ( frshits>0 ){ wr_end = frs_time[0]; }

    // Show the progress of the loop
    if (reader.GetCurrentEntry() % 100000 == 0) {
        int progress = (reader.GetCurrentEntry() * 100) / totalEntries;
        char spin = spinner[reader.GetCurrentEntry() / 1000 % 4];
        std::cout << "\rProcessing the tree " << reader.GetCurrentEntry() << " (" << progress << "%) " << spin << std::flush;
    }
  }

  std::cout << std::endl;


  std::cout << "Wr start time: " << wr_start << std::endl;
  std::cout << "Wr end time: " << wr_end << std::endl;
  std::cout << "File run time: " << (wr_end-wr_start)*1e-9/60 << " mins" << std::endl;

  std::cout << "ROOT Macro has run succesfully!" << std::endl;

  // Close the file
  std::exit(0);

}
