#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>

#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include "../../../c4Root/c4data/frsData/FrsHitData.h"

// *************************************************************************************
// ****************************** DEFINE SCRIPT CONSTANTS ****************************** 
// *************************************************************************************

namespace constants{

  const int DSSSD = 1;

  const int TOTAL_GERMANIUM_DETECTORS = 15; // Number of germanium detectors 
  const int TOTAL_GERMANIUM_CRYSTALS = 2; // Number of germanium crystals

  const std::pair<ULong64_t, ULong64_t> FRS_TIMEWINDOW = {-17083, -10789};

  const std::map<std::string, std::string> IMPLANT_GATES_INFILE_MAP = {
    {"82nb", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/new_gates/82nb_new.root"},
    {"zrnbmo", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/new_gates/zrnbmo_new.root"}
  };

}

// *************************************************************************************
// ****************************** DEFINE CUSTOM STRUCTURES ***************************** 
// *************************************************************************************

// Define the data structures for branches
struct implant_data_struct{
  long long time;
  long long time_x;
  long long time_y;
  double x;
  double y;
  double energy;
  double energy_x;
  double energy_y;
  int cluster_size_x;
  int cluster_size_y;
  int cluster_min_x;
  int cluster_max_x;
  int cluster_min_y;
  int cluster_max_y;
  int stopped;
  int dssd;
  int sp;
  int bp;
}aida_implant_data;
  
struct decay_data_struct{
  long long time;
  long long time_x;
  long long time_y;
  double x;
  double y;
  double energy;
  double energy_x;
  double energy_y;
  int cluster_size_x;
  int cluster_size_y;
  int cluster_min_x;
  int cluster_max_x;
  int cluster_min_y;
  int cluster_max_y;
  int dssd;
  int sp;
  int bp;
}aida_decay_data;

struct bplast_data_struct{
  uint64_t time;
  short id;
  double slow_tot;
  int sp;
  int bp;
}bplast_data;
  
struct germanium_data_struct{
  uint64_t time;
  double energy;
  int sp;
  /*int bp;*/
}germanium_data;

struct frs_data_struct{
  long long time;
  double aoq;
  double x2;
  double x4;
  double z41;
  double z42;
  double sci41E;
  double sci42E;
  int sp;
}frs_data;


// *************************************************************************************
// ****************************** START MACRO ******************************************
// *************************************************************************************


void makeUnstitchedAnatrees(const char* input, const char* output) {

  // *************************************************************************************
  // ****************************** OPEN & CREATE FILES **********************************
  // *************************************************************************************
  
  // Open the file
  TFile* file = TFile::Open(input);
  if (!file) {
    std::cerr << "Error: Could not open file " << std::endl;
    return;
  }

  // Get the existing tree
  TTree* old_tree = dynamic_cast<TTree*>(file->Get("evt"));
  if (!old_tree) {
    std::cerr << "Error: Could not find tree evt in file " << std::endl;
    return;
  }

  // *************************************************************************************
  // ****************************** CREATE GATE MULTIMAPS ********************************
  // *************************************************************************************

  // Create a multimap with TCUTG objects of each gate and a implant event counter for each gate defined in the script
  std::multimap< std::string, std::tuple<TCutG*, TCutG*, TCutG*> > gatedimplant_cuts_map = {};

  // Loop over gate map
  for ( auto itr : constants::IMPLANT_GATES_INFILE_MAP){
    
    // Open gate file
    TFile* gate_file = new TFile(itr.second.c_str());
    if (!gate_file){
      std::cerr << "Error: Could not open gate file " << std::endl;
      std::exit(1);
    }

    std::string gate_name = itr.first; // Get gate name

    // Extract TCutG objects from gates
    TCutG* x4aoq_cut = (TCutG*)gate_file->Get("cut_x4_AoQ");
    TCutG* zz2_cut = (TCutG*)gate_file->Get("cut_Z_Z2");
    TCutG* zsci42E_cut = (TCutG*)gate_file->Get("cut_Z_sci42E");
    
    // Add gate information to multimap
    gatedimplant_cuts_map.emplace( gate_name, std::make_tuple(x4aoq_cut, zz2_cut, zsci42E_cut) );   
  }


  // *************************************************************************************
  // ****************************** DEFINE TREEREADER ************************************
  // *************************************************************************************

  // Create a tree reader
  TTreeReader reader(old_tree);

  // Create tree reader values for the branches you want to read
  TTreeReaderValue<bool> spill(reader, "EventHeader.fSpillFlag");
  TTreeReaderArray<long long> implant_time(reader, "AidaImplantHits.Time");
  TTreeReaderArray<long long> implant_time_x(reader, "AidaImplantHits.TimeX");
  TTreeReaderArray<long long> implant_time_y(reader, "AidaImplantHits.TimeY");
  TTreeReaderArray<Int_t> implant_dssd(reader, "AidaImplantHits.DSSD");
  TTreeReaderArray<Bool_t> implant_stopped(reader, "AidaImplantHits.Stopped");
  TTreeReaderArray<double> implant_energy(reader, "AidaImplantHits.Energy");
  TTreeReaderArray<double> implant_energy_x(reader, "AidaImplantHits.EnergyX");
  TTreeReaderArray<double> implant_energy_y(reader, "AidaImplantHits.EnergyY");
  TTreeReaderArray<Double_t> implant_x(reader, "AidaImplantHits.StripX");
  TTreeReaderArray<Double_t> implant_y(reader, "AidaImplantHits.StripY");
  TTreeReaderArray<Int_t> implant_cluster_size_x(reader, "AidaImplantHits.ClusterSizeX");
  TTreeReaderArray<Int_t> implant_cluster_size_y(reader, "AidaImplantHits.ClusterSizeY");
  TTreeReaderArray<Int_t> implant_cluster_min_x(reader, "AidaImplantHits.StripXMin");
  TTreeReaderArray<Int_t> implant_cluster_max_x(reader, "AidaImplantHits.StripXMax");
  TTreeReaderArray<Int_t> implant_cluster_min_y(reader, "AidaImplantHits.StripYMin");
  TTreeReaderArray<Int_t> implant_cluster_max_y(reader, "AidaImplantHits.StripYMax");

  TTreeReaderArray<Int_t> implant_adc_dssd(reader, "AidaImplantCalAdcData.dssd");
  TTreeReaderArray<Int_t> implant_adc_side(reader, "AidaImplantCalAdcData.side");
  TTreeReaderArray<Int_t> implant_adc_strip(reader, "AidaImplantCalAdcData.strip");
  TTreeReaderArray<Int_t> implant_adc_fee(reader, "AidaImplantCalAdcData.fee");
  TTreeReaderArray<Double_t> implant_adc_energy(reader, "AidaImplantCalAdcData.energy");
  TTreeReaderArray<unsigned long long> implant_adc_time(reader, "AidaImplantCalAdcData.slowTime");

  TTreeReaderArray<u_int64_t> bplast_time(reader, "bPlastTwinpeaksCalData.fabsolute_event_time");
  // TTreeReaderArray<ULong_t> bplast_time(reader, "bPlastTwinpeaksCalData.fwr_t");
  TTreeReaderArray<ushort> bplast_id(reader, "bPlastTwinpeaksCalData.fdetector_id");
  TTreeReaderArray<Double_t> bplast_slow_tot(reader, "bPlastTwinpeaksCalData.fslow_ToT");

  TTreeReaderArray<long long> decay_time(reader, "AidaDecayHits.Time");
  TTreeReaderArray<Int_t> decay_dssd(reader, "AidaDecayHits.DSSD");
  TTreeReaderArray<long long> decay_time_x(reader, "AidaDecayHits.TimeX");
  TTreeReaderArray<long long> decay_time_y(reader, "AidaDecayHits.TimeY");
  TTreeReaderArray<Double_t> decay_energy(reader, "AidaDecayHits.Energy");
  TTreeReaderArray<Double_t> decay_energy_x(reader, "AidaDecayHits.EnergyX");
  TTreeReaderArray<Double_t> decay_energy_y(reader, "AidaDecayHits.EnergyY");
  TTreeReaderArray<Double_t> decay_x(reader, "AidaDecayHits.StripX");
  TTreeReaderArray<Double_t> decay_y(reader, "AidaDecayHits.StripY");
  TTreeReaderArray<Int_t> decay_cluster_size_x(reader, "AidaDecayHits.ClusterSizeX");
  TTreeReaderArray<Int_t> decay_cluster_size_y(reader, "AidaDecayHits.ClusterSizeY");
  TTreeReaderArray<Int_t> decay_cluster_min_x(reader, "AidaDecayHits.StripXMin");
  TTreeReaderArray<Int_t> decay_cluster_max_x(reader, "AidaDecayHits.StripXMax");
  TTreeReaderArray<Int_t> decay_cluster_min_y(reader, "AidaDecayHits.StripYMin");
  TTreeReaderArray<Int_t> decay_cluster_max_y(reader, "AidaDecayHits.StripYMax");

  TTreeReaderArray<Int_t> decay_adc_dssd(reader, "AidaDecayCalAdcData.dssd");
  TTreeReaderArray<Int_t> decay_adc_side(reader, "AidaDecayCalAdcData.side");
  TTreeReaderArray<Int_t> decay_adc_strip(reader, "AidaDecayCalAdcData.strip");
  TTreeReaderArray<Int_t> decay_adc_fee(reader, "AidaDecayCalAdcData.fee");
  TTreeReaderArray<Double_t> decay_adc_energy(reader, "AidaDecayCalAdcData.energy");
  TTreeReaderArray<unsigned long long> decay_adc_time(reader, "AidaDecayCalAdcData.slowTime");

  TTreeReaderArray<uint64_t> germanium_time(reader, "GermaniumCalData.fwr_t");
  TTreeReaderArray<int64_t> germanium_abs_ev_time(reader, "GermaniumCalData.fabsolute_event_time");
  TTreeReaderArray<uint> germanium_det(reader, "GermaniumCalData.fdetector_id");
  TTreeReaderArray<uint> germanium_cry(reader, "GermaniumCalData.fcrystal_id");
  TTreeReaderArray<Double_t> germanium_energy(reader, "GermaniumCalData.fchannel_energy");

  TTreeReaderArray<FrsMultiHitItem> FrsMultiItem(reader, "FrsMultiHitData");
  TTreeReaderArray<FrsHitItem> FrsItem(reader, "FrsHitData");

  // *************************************************************************************
  // ****************************** DEFINE OUTPUT TREES **********************************
  // *************************************************************************************

  // Open the output file
  TFile* outputFile = new TFile(output, "RECREATE");
  if (!outputFile) {
      std::cerr << "Error: Could not create output file " << std::endl;
      return;
  }
  
  TTree* frs_tree = new TTree("frs_tree", "New AIDA Analysis Tree");  
  frs_tree->SetDirectory(0);
  frs_tree->Branch("frs", &frs_data, "time/l:aoq/D:x2/D:x4/D:z41/D:z42/D:sci41E/D:sci42E/D:sp/I");

  // Create implant tree and branches for anatree
  TTree* implant_tree = new TTree("aida_implant_tree", "New AIDA Analysis Tree");
  implant_tree->Branch("implant", &aida_implant_data, "time/l:time_x/l:time_y/l:x/D:y/D:e/D:ex/D:ey/D:csx/I:csy/I:cminx/I:cmaxx/I:cminy/I:cmaxy/I:stopped/I:dssd/I:sp/I:bp/I");

  // Create decay tree and branches for anatree
  TTree* decay_tree = new TTree("aida_decay_tree", "New AIDA Analysis Tree");
  decay_tree->Branch("decay", &aida_decay_data, "time/l:time_x/l:time_y/l:x/D:y/D:e/D:ex/D:ey/D:csx/I:csy/I:cminx/I:cmaxx/I:cminy/I:cmaxy/I:dssd/I:sp/I:bp/I");

  // Create germanium tree and branches for anatree
  TTree* germanium_tree = new TTree("germanium_tree", "New AIDA Analysis Tree");
  germanium_tree->Branch("germanium", &germanium_data, "time/l:energy/D:sp/I");

  // Create bplast tree and branches for anatree*/
  TTree* bplast_tree = new TTree("bplast_tree", "New AIDA Analysis Tree");
  bplast_tree->Branch("bplast", &bplast_data, "time/l:id/S:slow_tot/D:sp/I:bp/I");

  // Define a map to contain all the gated implant trees
  std::map<std::string, TTree*> gatedimplant_trees_map = {};

  // Loop over gate and make trees for each gated implant isotope
  for ( auto itr : constants::IMPLANT_GATES_INFILE_MAP ){
    
    // Get the fist element of the tuple containing the name
    std::string implant_gate_name = itr.first;
    
    // Create name strings for tree and branch
    std::string tree_name = "aida_gatedimplant_" + implant_gate_name + "_tree";
    std::string branch_name = "gatedimplant_" + implant_gate_name;

    // Create tree and define branch structure
    TTree* gatedimplant_tree = new TTree(tree_name.c_str(), "New AIDA Analysis Tree");
    gatedimplant_tree->SetDirectory(0); // Save tree object im memory untill explicitly written to file
    gatedimplant_tree->Branch(branch_name.c_str(), &aida_implant_data, "time/l:time_x/l:time_y/l:x/D:y/D:e/D:ex/D:ey/D:csx/I:csy/I:cminx/I:cmaxx/I:cminy/I:cmaxy/I:stopped/I:dssd/I:sp/I:bp/I");

    // Add tree to set
    gatedimplant_trees_map.emplace(implant_gate_name, gatedimplant_tree);
  }


  // *************************************************************************************
  // **************** MAKE TEMPORARY MAPS FOR FRS-IMPLANT-BPlast CORRELATION *************
  // *************************************************************************************

  std::multimap<uint64_t, frs_data_struct> frs_temporary_map = {};
  std::multimap<uint64_t, implant_data_struct> implant_temporary_map = {};
  std::multimap<uint64_t, bplast_data_struct> bplast_temporary_map = {};

  // *************************************************************************************
  // ****************************** LOOP OVER ALL EVENTS *********************************
  // *************************************************************************************

  const char spinner[] = {'-', '\\', '|', '/'};
  int totalEntries = reader.GetEntries(true);

  // Loop over all entries in the old tree
  for (int current_entry = 0; current_entry<totalEntries; ++current_entry) {

    reader.SetEntry(current_entry);

    // *************************************************************************************
    // ****************************** GET EVENT INFORMATION *********************************
    // *************************************************************************************
     
    // sizes of events
    int germaniumhits = germanium_time.GetSize();
    int aidadecayhits = decay_time.GetSize();
    int aidaimphits = implant_time.GetSize();
    int bplasthits = bplast_id.GetSize();
      
    int spflag = 0;
    int bpflag = 0;
    int bp1flag = 0;
    int bp2flag = 0;

    // Define spill flag for event
    if(*spill == true) spflag = 1;
    if(*spill == false) spflag = 2;
      
    // *************************************************************************************
    // ****************************** LOOP OVER BPLAST  ************************************
    // *************************************************************************************
      
    for (int j = 0; j < bplasthits; j++) {
      if(bplast_id[j] < 64) bp1flag = 1;
      if(bplast_id[j] > 63 && bplast_id[j] < 128) bp2flag = 1;
    }
  
    if (bp1flag == 1 && bp2flag == 0) bpflag = 1;
    if (bp1flag == 0 && bp2flag == 1) bpflag = 2;
    if (bp1flag == 1 && bp2flag == 1) bpflag = 3;
          
    std::set<int> bplast_filledtree = {};
    for(int j=0; j<bplasthits; j++){
      if (bplast_filledtree.count(j) == 0){
        bplast_data.time = bplast_time[j];
        bplast_data.id = bplast_id[j];
        bplast_data.slow_tot = bplast_slow_tot[j];
        bplast_data.sp = spflag;
        bplast_data.bp = bpflag;
        bplast_tree->Fill();
        bplast_filledtree.insert(j);

        bplast_temporary_map.emplace(bplast_time[j], bplast_data);
      }
    }


    // *************************************************************************************
    // ****************************** LOOP OVER AIDA IMPLANTS  *****************************
    // *************************************************************************************

    // Initialise flag for stopped implant
    int stopped = 0;

    std::set<int> aidaimplant_filledtree = {};

    if(aidaimphits > 0){
      for (int j = 0; j < aidaimphits; j++) {

        if (aidaimplant_filledtree.count(j) == 0){
          // Determine if implant subevent has stopped in specific DSSSD
          if (implant_dssd[j] == constants::DSSSD && bpflag == 0){
            stopped = 1;
          } 
          else {
            stopped = 2;
          }
              
          // Fill up implant datastruct
          aida_implant_data.time = implant_time[j];
          aida_implant_data.time_x = implant_time_x[j];
          aida_implant_data.time_y = implant_time_y[j];
          aida_implant_data.stopped = stopped;
          aida_implant_data.dssd = implant_dssd[j];
          aida_implant_data.x = implant_x[j];
          aida_implant_data.y = implant_y[j];
          aida_implant_data.energy = implant_energy[j];
          aida_implant_data.energy_x = implant_energy_x[j];
          aida_implant_data.energy_y = implant_energy_y[j];
          aida_implant_data.cluster_size_x = implant_cluster_size_x[j];
          aida_implant_data.cluster_size_y = implant_cluster_size_y[j];
          aida_implant_data.cluster_min_x = implant_cluster_min_x[j];
          aida_implant_data.cluster_max_x = implant_cluster_max_x[j];
          aida_implant_data.cluster_min_y = implant_cluster_min_y[j];
          aida_implant_data.cluster_max_y = implant_cluster_max_y[j];
          aida_implant_data.sp = spflag;
          aida_implant_data.bp = bpflag;

          // Fill implant tree with subevent
          implant_tree->Fill();

          aidaimplant_filledtree.insert(j);

          implant_temporary_map.emplace(implant_time[j], aida_implant_data);

        }
      }
    }

    // *************************************************************************************
    // ****************************** LOOP OVER DECAY EVENTS *******************************
    // *************************************************************************************

    int timestitchedXStripMultiplicity = 0;
    int timestitchedYStripMultiplicity = 0;

    for (int i=0; i < decay_adc_side.GetSize(); ++i){
      if (decay_adc_side[i] == 1 && decay_adc_dssd[i] == 1) ++timestitchedXStripMultiplicity;
      if (decay_adc_side[i] == -1 && decay_adc_dssd[i] == 1) ++timestitchedYStripMultiplicity;
    }

    // Empty set for each subevents
    std::set<int> aidadecay_filledtree{};

    // Loop over decay subevents
    if(true){
    // if(aidadecayhits == 1 && aidaimphits == 0){
    // if( timestitchedXStripMultiplicity <= 3 && timestitchedYStripMultiplicity <= 3 ){

      for (int i = 0; i < aidadecayhits; i++) {

        if (aidadecay_filledtree.count(i) == 0) {

          aida_decay_data.time = decay_time[i];
          aida_decay_data.time_x = decay_time_x[i];
          aida_decay_data.time_y = decay_time_y[i];
          aida_decay_data.x = decay_x[i];
          aida_decay_data.y = decay_y[i];
          aida_decay_data.energy = decay_energy[i];
          aida_decay_data.energy_x = decay_energy_x[i];
          aida_decay_data.energy_y = decay_energy_y[i];
          aida_decay_data.cluster_size_x = decay_cluster_size_x[i];
          aida_decay_data.cluster_size_y = decay_cluster_size_y[i];
          aida_decay_data.cluster_min_x = decay_cluster_min_x[i];
          aida_decay_data.cluster_max_x = decay_cluster_max_x[i];
          aida_decay_data.cluster_min_y = decay_cluster_min_y[i];
          aida_decay_data.cluster_max_y = decay_cluster_max_y[i];
          aida_decay_data.dssd = decay_dssd[i];
          aida_decay_data.sp = spflag;
          aida_decay_data.bp = bpflag;

          decay_tree->Fill();

          aidadecay_filledtree.insert(i);
        }


      } // End of decay loop

    }

    // *************************************************************************************
    // ****************************** LOOP OVER GERMANIUM EVENTS *******************************
    // *************************************************************************************

    
    std::set<int> germanium_filledtree{};

    // Loop over decay subevents
    for (int j = 0; j < germaniumhits; j++){

      if ( germanium_det[j] <= constants::TOTAL_GERMANIUM_DETECTORS && germanium_cry[j] <= constants::TOTAL_GERMANIUM_CRYSTALS && germanium_energy[j] > 20. ){

        if ( germanium_filledtree.count(j) == 0 ){

          germanium_data.time = germanium_time[j];
          germanium_data.energy = germanium_energy[j];
          germanium_data.sp = spflag;

          germanium_tree->Fill();

        }
        
      germanium_filledtree.insert(j);

      }

    }

    // *************************************************************************************
    // ****************************** LOOP OVER FRS EVENTS *******************************
    // *************************************************************************************


    for (int mbsIdx=0; mbsIdx<FrsMultiItem.GetSize(); mbsIdx++){

      std::set<int> frs_filledtree = {};

      // Get corresponding FrsHitItem and FrsMultiHitItem
      auto& frs_multi_item = FrsMultiItem[mbsIdx];
      auto& frs_item = FrsItem[mbsIdx];

      // Get values inside the item
      std::vector<float> const& AoQ = frs_multi_item.Get_ID_AoQ_corr_s2s4_mhtdc();
      std::vector<float> const& Z = frs_multi_item.Get_ID_z41_mhtdc();
      std::vector<float> const& Z2 = frs_multi_item.Get_ID_z42_mhtdc();
      std::vector<float> const& X2 = frs_multi_item.Get_ID_s2x_s2s4_mhtdc();
      std::vector<float> const& X4 = frs_multi_item.Get_ID_s4x_mhtdc();
      Long64_t const& frs_wr = frs_item.Get_wr_t();
      float const& sci41E = frs_item.Get_sci_e_41();
      float const& sci42E = frs_item.Get_sci_e_42();

      for (int mhtdcIdx=0; mhtdcIdx<AoQ.size(); mhtdcIdx++){
        
        if ( mhtdcIdx >= Z.size() || mhtdcIdx >= Z2.size() || mhtdcIdx >= X4.size() || mhtdcIdx >= X4.size()) { continue; }

        frs_data.time = frs_wr;
        frs_data.sci41E = sci41E;
        frs_data.sci42E = sci42E;
        frs_data.aoq = AoQ[mhtdcIdx];
        frs_data.x2 = X2[mhtdcIdx];
        frs_data.x4 = X4[mhtdcIdx];
        frs_data.z41 = Z[mhtdcIdx];
        frs_data.z42 = Z2[mhtdcIdx];
        frs_data.sp = spflag;

        frs_filledtree.insert(mhtdcIdx);

        frs_tree->Fill();
      
        frs_temporary_map.emplace(frs_wr, frs_data);
      }

    }

    // *************************************************************************************
    // ****************************** PRINT OUT LOOP PROGRESS ******************************
    // *************************************************************************************
    
    // Show the progress of the loop
    if (reader.GetCurrentEntry() % 100000 == 0) {
      int progress = (reader.GetCurrentEntry() * 100) / totalEntries;
      char spin = spinner[reader.GetCurrentEntry() / 1000 % 4];
      std::cout << "\rProcessing the tree " << reader.GetCurrentEntry() << " (" << progress << "%) " << spin << std::flush;
    }
  
  }

  std::cout << std::endl << "Finished looping through TreeReader" << std::endl;

  // *************************************************************************************
  // ********************** MAKE FRS-IMPLANT-BPlast CORRELATIONS **************************
  // *************************************************************************************

  int implantCounter = 0;

  for (auto implantItr=implant_temporary_map.begin(); implantItr!=implant_temporary_map.end(); implantItr++){

    long long implantTime = implantItr->first;

    // THERE IS A HIDDEN BUG HERE. IF WE REMOVE THE +/- 50e3 THERE ARE CASES WHERE THE PROGRAM HANGS.
    // MAYBE EDGE CASE WHEN LOWER_BOUND == UPPER_BOUND??
    const auto& lowerBound = frs_temporary_map.lower_bound(implantTime-constants::FRS_TIMEWINDOW.first - 50e3);
    const auto& upperBound = frs_temporary_map.upper_bound(implantTime+constants::FRS_TIMEWINDOW.second + 50e3);

    for (auto frsItr=lowerBound; frsItr!=upperBound; frsItr++){

      Long64_t timediff = frsItr->first - implantTime ;
      
      // Correlated FRS IMPLANT EVENT 
      if (timediff > constants::FRS_TIMEWINDOW.first && timediff < constants::FRS_TIMEWINDOW.second ){


        // Get data
        implant_data_struct implantData = implantItr->second; 
        frs_data_struct frsData = frsItr->second; 

        // Loop over different gates
        for (auto gimp_itr : gatedimplant_cuts_map){
        
          // Unpack entries inside map 
          std::string gimp_key = gimp_itr.first;
          auto [x4aoq_cut, zz2_cut, zsci42E_cut] = gimp_itr.second;

          if(x4aoq_cut->IsInside(frsData.aoq, frsData.x4) && zz2_cut->IsInside(frsData.z41, frsData.z42) && zsci42E_cut->IsInside(frsData.sci42E, frsData.z41)){
            aida_implant_data = implantData;
            gatedimplant_trees_map[gimp_key]->Fill();
          
          }

        }

      } 

    }
  }


  // *************************************************************************************
  // ****************************** write trees ******************************************
  // *************************************************************************************

  std::cout << "Writing trees ..." << std::endl;
  outputFile->cd();
  
  frs_tree->Write();
  std::cout << "Finished frs tree" << std::endl;
  implant_tree->Write();
  std::cout << "Finished implant tree" << std::endl;
  decay_tree->Write();
  std::cout << "Finished decay tree" << std::endl;
  germanium_tree->Write();
  std::cout << "Finished germanium tree" << std::endl;
  
  // Loop over trees and Write
  for (auto itr : gatedimplant_trees_map){
    
    // Extract tree and write
    itr.second->Write();
    std::cout << "Finished " + itr.first + " tree" << std::endl;
  }

  std::cout << "Finished writing trees successfully!" << std::endl << std::endl;

  // *************************************************************************************
  // ****************************** CLEANUP **********************************************
  // *************************************************************************************

  // Close the input file
  file->Close();

  // CLose and write the output file
  outputFile->Close();
  
  // Delete pointer objects
  delete file;
  delete outputFile;
}
