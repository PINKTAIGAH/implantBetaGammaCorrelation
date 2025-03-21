#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

void aida_frs_efficiency(const char* input, const char* gate_dir){

  // Open the TFile
  TFile* file = TFile::Open(input);
  if (!file){
    std::cerr << "Error: Could not open file " << std::endl;
    return;
  }

  // Get the existing tree
  TTree* old_tree = dynamic_cast<TTree*>(file->Get("evt"));
  if (!old_tree){
    std::cerr << "Error: Could not find tree evt in file " << std::endl;
    return;
  }

  // Create a tree reader
  TTreeReader reader(old_tree);

  // Define variables of interest
  TTreeReaderValue<bool> spill(reader, "EventHeader.fSpillFlag");

  TTreeReaderArray<int64_t> implant_time(reader, "AidaImplantHits.Time");
  TTreeReaderArray<Int_t> implant_dssd(reader, "AidaImplantHits.DSSD");
  TTreeReaderArray<int64_t> implant_time_x(reader, "AidaImplantHits.TimeX");
  TTreeReaderArray<int64_t> implant_time_y(reader, "AidaImplantHits.TimeY");
  TTreeReaderArray<Double_t> implant_energy(reader, "AidaImplantHits.Energy");
  TTreeReaderArray<Double_t> implant_energy_x(reader, "AidaImplantHits.EnergyX");
  TTreeReaderArray<Double_t> implant_energy_y(reader, "AidaImplantHits.EnergyY");
  TTreeReaderArray<Double_t> implant_strip_x(reader, "AidaImplantHits.StripX");
  TTreeReaderArray<Double_t> implant_strip_y(reader, "AidaImplantHits.StripY");
  TTreeReaderArray<Int_t> implant_cluster_size_x(reader, "AidaImplantHits.ClusterSizeX");
  TTreeReaderArray<Int_t> implant_cluster_size_y(reader, "AidaImplantHits.ClusterSizeY");

  TTreeReaderArray<FrsMultiHitItem> FrsMultiItem(reader, "FrsMultiHitData");
  /*TTreeReaderArray<Float_t> frs_z(reader, "FrsHitData.fID_z42");*/
  /*TTreeReaderArray<Float_t> frs_AoQ(reader, "FrsHitData.fID_AoQ_s2s4");*/
  /*TTreeReaderArray<int64_t> frs_time(reader, "FrsHitData.fwr_t");*/
  

  /*// Open output file*/
  /*TFile* outputFile = new TFile(output, "RECREATE");*/
  /*if (!outputFile){*/
  /*    std::cerr << "Error: Could not create output file " << std::endl;*/
  /*    return;*/
  /*}*/

  /*// Define histograms*/
  /*TH2D* h2_frs_pid = new TH2D("h2_frs_pid", "FRS PID @ s4 after slits", 800, 1.9, 2.3, 800, 30, 40);*/
  /*TH2D* h2_gated_frs_pid = new TH2D("h2_aida_frs_pid", "Gated FRS PID @ s4 after slits", 800, 1.9, 2.3, 800, 30, 40);*/
  /*TH2D* h2_aida_implanted_hit_pattern = new TH2D("h2_aida_implanted_hit_pattern", "AIDA hit pattern (Implanted)", 384, 0, 384, 128, 0, 128 );*/
  /*TH2D* h2_aida_stopped_hit_pattern = new TH2D("h2_aida_stopped_hit_pattern", "AIDA hit pattern (Stopped)", 384, 0, 384, 128, 0, 128 );*/
  /*TH2D* h2_gated_aida_implanted_hit_pattern = new TH2D("h2_gated_aida_implanted_hit_pattern", "AIDA hit pattern (Gated Implanted)", 384, 0, 384, 128, 0, 128 );*/
  /*TH2D* h2_gated_aida_stopped_hit_pattern = new TH2D("h2_gated_aida_stopped_hit_pattern", "AIDA hit pattern (Gated Stopped)", 384, 0, 384, 128, 0, 128 );*/
  /**/
  // Get total entries
  int totalEntries = reader.GetEntries();

  // Define spinner
  const char spinner[] = {'-', '\\', '|', '/',};

  // FRS Gated
  TFile* gate = new TFile(gate_dir);
  if (!gate){
    std::cerr << "Error: Could not find FRSCUT in file " << std::endl;
    return;
  }
  TCutG* zAoQ_cut = (TCutG*)gate->Get("cut_Z_AoQ");
  TCutG* z1z2_cut = (TCutG*)gate->Get("cut_Z_Z2");

  int totalGatedHits = 0;
  int totalImplantedGatedHits = 0;
  int totalStoppedGatedHits = 0;

  int spflag = 0;

  while (reader.Next()){

    int aidaImplantHits = implant_dssd.GetSize();
    
    if (*spill == true){spflag = 1;}
    if (*spill == false){spflag = 2;}

    // Define set to hold index of altready filled aida events
    std::unordered_set<int> filled_gatedimplanttree{};

    // Loop over frs events
    for (auto const& FrsMultiItem : FrsMultiItem){
      // Define frs data items
      std::vector<float> const& AoQ = FrsMultiItem.Get_ID_AoQ_corr_s2s4_mhtdc();
      std::vector<float> const& Z = FrsMultiItem.Get_ID_z41_mhtdc();
      std::vector<float> const& Z2 = FrsMultiItem.Get_ID_z42_mhtdc();


      // loop over frs subevents
      for (int i=0; i<AoQ.size(); i++){

        // Check if indexes are out of bounds & skip
        if ( i>=AoQ.size() || i>=Z.size() ){
          continue;
        }

        if ( zAoQ_cut->IsInside(AoQ[i], Z[i]) && z1z2_cut->IsInside(Z[i], Z2[i]) ){
          totalGatedHits++;
          // Check if we have aida hits
          if ( aidaImplantHits>0 ){
            for ( int j=0; j<aidaImplantHits; j++ ){
              if ( implant_dssd[j]==1 && filled_gatedimplanttree.count(j)==0 ){
                totalImplantedGatedHits++;
                filled_gatedimplanttree.insert(j);
              }
            }
          }
        }
      }
    }


    // Give update on status of loop
    if (reader.GetCurrentEntry() % (totalEntries/100) == 0){
      int progress = (reader.GetCurrentEntry() * 100) / totalEntries;
      char spin = spinner[reader.GetCurrentEntry() / (totalEntries/100) % 4];
      std::cout << "\rProcessing the tree " << reader.GetCurrentEntry() << " (" << progress << "%) " << spin << std::flush;
    }
  }

  // Compute implant efficiencies
  float implant_efficiency = (float)totalImplantedGatedHits/(float)totalGatedHits*100;

  std::cout << std::endl << std::endl;
  std::cout << "The total number of FRS gated hits @ s4 is: " << totalGatedHits << std::endl; 
  std::cout << std::endl << std::endl;
  std::cout << "The total number of implanted FRS gated hits @ s4 is: " << totalImplantedGatedHits << " (" << implant_efficiency << " % efficiency)"<< std::endl; 

  std::cout << "The macro has ended succesfully!" << std::endl;

  file->Close();
  delete file;

  gate->Close();
  delete gate;

  std::exit(0);
}
