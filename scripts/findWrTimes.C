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
  TTreeReaderArray<long long> time(reader, "bPlastTwinpeaksCalData.fwr_t");

  // Define the counter
  int counter = 0;
  ULong_t wr_start = 0;
  ULong_t wr_end = 0;

  const char spinner[] = {'-', '\\', '|', '/'};
  int totalEntries = reader.GetEntries(true);
  bool jumpedToEnd = false;

  // Loop over all the entries in the new tree
  while (reader.Next()){

    int hits = time.GetSize();

    if (hits>0 && wr_start==0 ){ wr_start = time[0]; }
    if ( hits>0 ){ wr_end = time[0]; }

    // Change entry to last 1000 entries if start time is assigned and we havent aready jumped
    if ( !jumpedToEnd && wr_start!=0 ){
      reader.SetEntry(totalEntries-1000);
      jumpedToEnd = true;
    }

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
