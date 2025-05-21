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

  const std::string BROKEN_STRIPS_INFILE = "/lustre/gamma/jeroen/S100/analysis/betaion/AIDA_strips.txt";

  const int TOTAL_GERMANIUM_DETECTORS = 15; // Number of germanium detectors 
  const int TOTAL_GERMANIUM_CRYSTALS = 2; // Number of germanium crystals

  const int LOOKAHEAD = 10 ; // Number of entries to look ahead for implant-FRS coincidence
  const std::pair<ULong64_t, ULong64_t> FRS_TIMEWINDOW = {-16890, -10700};

  const std::map<std::string, std::string> IMPLANT_GATES_INFILE_MAP = {
    // {"82nb", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/82nb_bb7.root"},
    // {"81zr", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/81zr_bb7.root"}

    {"82nb", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/82nb.root"},
    {"84nb", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/84nb.root"},
    {"84mo", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/84mo.root"},
    {"85mo", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/85mo.root"},
    {"81zr", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/81zr.root"},
    {"86tc", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/86tc.root"}

    // {"70br", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/70br.root"}
  };

}

// *************************************************************************************
// ****************************** DEFINE CUSTOM STRUCTURES ***************************** 
// *************************************************************************************

// Define the data structures for branches
struct implant_data{
  long long time;
  int stopped;
  int dssd;
  double x;
  double y;
  double energy;
  double energy_x;
  double energy_y;
  int cluster_size_x;
  int cluster_size_y;
  int sp;
  int bp;
}aida_implant_data;
  
struct decay_data{
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
  int dssd;
  int sp;
  int bp;
}aida_decay_data;

struct bplast_data{
  uint64_t time;
  short id;
  double slow_tot;
  int sp;
  int bp;
}bplast_data;
  
struct germanium_data{
  uint64_t time;
  double energy;
  int sp;
  /*int bp;*/
}germanium_data;


// *************************************************************************************
// ****************************** DEFINE FUNCTIONS *************************************
// *************************************************************************************

std::tuple< std::vector<int>, std::vector<int> > loadBrokenStrips(){
  // Read in broken strips for both upstream and downstream DSSSD which are to be skipped 
  
  std::ifstream aida_strips_file(constants::BROKEN_STRIPS_INFILE); // Open file

  // Declare function variables
  std::vector<int> broken_xstrips;
  std::vector<int> broken_ystrips;
  int dssd_number;
  std::string xy;
  int strip_number;
  int threshold;

  std::string line;

  while(std::getline(aida_strips_file, line)){
    // Skip comments
    if(line[0] == '#') continue;

    std::istringstream iss(line);
    iss >> dssd_number >> xy >> strip_number >> threshold;

    // If the threshold is -1, add the strip to the list of strips to skip
    if ( threshold == -1 ){
      if ( xy == "X" ){
        strip_number = strip_number - 1;
        broken_xstrips.push_back(strip_number);
      } 
      else if ( xy == "Y" ){
        strip_number = strip_number - 1;
        broken_ystrips.push_back(strip_number);
      }
    }
  }

  // Print out broken strips 
  for( auto x : broken_xstrips ){
    std::cout << "Broken X strip: " << x << std::endl;
  }
  
  for(auto y : broken_ystrips){
    std::cout << "Broken Y strip: " << y << std::endl;
  }

  return std::make_tuple(broken_xstrips, broken_ystrips);
}


bool isBrokenStrip(std::vector<double> broken_strip_vector, double event_strip){
  // Function which will loop over vector containing broken strips and check 
  // if current event occured in a broken strip

  bool isBroken = false; // Flag telling if event is from broken strip

  // Loop over vector if not empty
  if (!broken_strip_vector.empty()){
    // Loop for each noisy strip
    for (auto & broken_strip : broken_strip_vector){
      // If event strip corresponds to a broken strip, change flag and break 
      if ( event_strip == broken_strip ){ 
        isBroken = true;
        break; 
      }

    }

  }

  return isBroken;
}

// *************************************************************************************
// ****************************** START MACRO ******************************************
// *************************************************************************************


void makeAnatrees(const char* input, const char* output) {

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
  std::multimap< std::string, std::tuple<TCutG*, TCutG*> > gatedimplant_cuts_map = {};

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
    TCutG *x4aoq_cut = (TCutG*)gate_file->Get("cut_x4_AoQ");
    TCutG *zz2_cut = (TCutG*)gate_file->Get("cut_Z_Z2");
    
    // Add gate information to multimap
    gatedimplant_cuts_map.emplace( gate_name, std::make_tuple(x4aoq_cut, zz2_cut) );   
  }

  // Create a map which will hold the counters for each gated implant species
  std::map<std::string, int> gatedimplant_counter_map = {};
  std::map<std::string, int> gatedimplant_stopped_counter_map = {};
  
  // Loop over map and fill with zerod counter for each gate
  for (auto itr : constants::IMPLANT_GATES_INFILE_MAP){

    std::string gate_name = itr.first; // Get gate name

    // Fill map
    gatedimplant_counter_map.emplace( gate_name, 0);
    gatedimplant_stopped_counter_map.emplace( gate_name, 0);
  }

  // *************************************************************************************
  // ****************************** DEFINE TREEREADER ************************************
  // *************************************************************************************

  // Create a tree reader
  TTreeReader reader(old_tree);

  // Create tree reader values for the branches you want to read
  TTreeReaderValue<bool> spill(reader, "EventHeader.fSpillFlag");
  TTreeReaderArray<long long> implant_time(reader, "AidaImplantHits.Time");
  TTreeReaderArray<Int_t> implant_dssd(reader, "AidaImplantHits.DSSD");
  TTreeReaderArray<Bool_t> implant_stopped(reader, "AidaImplantHits.Stopped");
  TTreeReaderArray<double> implant_energy(reader, "AidaImplantHits.Energy");
  TTreeReaderArray<double> implant_energy_x(reader, "AidaImplantHits.EnergyX");
  TTreeReaderArray<double> implant_energy_y(reader, "AidaImplantHits.EnergyY");
  TTreeReaderArray<Double_t> implant_x(reader, "AidaImplantHits.StripX");
  TTreeReaderArray<Double_t> implant_y(reader, "AidaImplantHits.StripY");
  TTreeReaderArray<Int_t> implant_cluster_size_x(reader, "AidaImplantHits.ClusterSizeX");
  TTreeReaderArray<Int_t> implant_cluster_size_y(reader, "AidaImplantHits.ClusterSizeY");

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
  
  // Create implant tree and branches for anatree
  TTree* implant_tree = new TTree("aida_implant_tree", "New AIDA Analysis Tree");
  implant_tree->SetDirectory(0); // Save tree object im memory untill explicitly written to file
  implant_tree->Branch("implant", &aida_implant_data, "time/l:stopped/I:dssd/I:x/D:y/D:e/D:ex/D:ey/D:csx/I:csy/I:sp/I:bp/I");

  // Create decay tree and branches for anatree
  TTree* decay_tree = new TTree("aida_decay_tree", "New AIDA Analysis Tree");
  decay_tree->SetDirectory(0); // Save tree object im memory untill explicitly written to file
  decay_tree->Branch("decay", &aida_decay_data, "time/l:time_x/l:time_y/l:x/D:y/D:e/D:ex/D:ey/D:csx/I:csy/I:dssd/I:sp/I:bp/I");

  // Create germanium tree and branches for anatree
  TTree* germanium_tree = new TTree("germanium_tree", "New AIDA Analysis Tree");
  germanium_tree->SetDirectory(0); // Save tree object im memory untill explicitly written to file
  germanium_tree->Branch("germanium", &germanium_data, "time/l:energy/D:sp/I");

  // Create bplast tree and branches for anatree*/
  // TTree* bplast_tree = new TTree("bplast_tree", "New AIDA Analysis Tree");
  // bplast_tree->Branch("bplast", &bplast_data, "time/l:id/S:slow_tot/D:sp/I:bp/I");

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
    gatedimplant_tree->Branch(branch_name.c_str(), &aida_implant_data, "time/l:stopped/I:dssd/I:x/D:y/D:e/D:ex/D:ey/D:csx/I:csy/I:sp/I:bp/I");

    // Add tree to set
    gatedimplant_trees_map.emplace(implant_gate_name, gatedimplant_tree);
  }

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
    /*int frshits = frs_time.GetSize();*/

    int mult_x = decay_x.GetSize();
    int mult_y = decay_y.GetSize();

      
    int spflag = 0;
    int bpflag = 0;
    int bp1flag = 0;
    int bp2flag = 0;

    /*cout << "# of FRS hits:  " << frshits << endl;*/
    /*cout << "# of implant hits:  " << aidaimphits << endl;*/
    /*cout << "# of decay hits:  " << aidadecayhits << endl;*/
    /*cout << "# of gamma hits:  " << germaniumhits << endl;*/
    /*if (bplasthits > 0) cout << "# of bplast hits:  " << bplasthits << endl;*/

    // Define spill flag for event
    if(*spill == true) spflag = 1;
    if(*spill == false) spflag = 2;
      
    // if (*spill == true) std::cout << "[DEBUG] New Entry" << std::endl;

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
          
    /*cout << bp1flag << "  " << bp2flag << "  " << bpflag << endl;*/
    
    // for(int j=0; j<bplasthits; j++){
    //   bplast_data.time = bplast_time[j];
    //   bplast_data.id = bplast_id[j];
    //   bplast_data.slow_tot = bplast_slow_tot[j];
    //   bplast_data.sp = spflag;
    //   bplast_data.bp = bpflag;
    //   bplast_tree->Fill();
    // }


    // *************************************************************************************
    // ****************************** LOOP OVER AIDA IMPLANTS  *****************************
    // *************************************************************************************

    // Initialise flag for stopped implant
    int stopped = 0;

    if(aidaimphits > 0){
      for (int j = 0; j < aidaimphits; j++) {

        // Determine if implant subevent has stopped in specific DSSSD
        if (implant_dssd[j] == constants::DSSSD && bpflag == 0){
          stopped = 1;
        } 
        else {
          stopped = 2;
        }
              
        // Fill up implant datastruct
        aida_implant_data.time = implant_time[j];
        aida_implant_data.stopped = stopped;
        aida_implant_data.dssd = implant_dssd[j];
        aida_implant_data.x = implant_x[j];
        aida_implant_data.y = implant_y[j];
        aida_implant_data.energy = implant_energy[j];
        aida_implant_data.energy_x = implant_energy_x[j];
        aida_implant_data.energy_y = implant_energy_y[j];
        aida_implant_data.cluster_size_x = implant_cluster_size_x[j];
        aida_implant_data.cluster_size_y = implant_cluster_size_y[j];
        aida_implant_data.sp = spflag;
        aida_implant_data.bp = bpflag;

        // Fill implant tree with subevent
        implant_tree->Fill();
      }
    }

    // *************************************************************************************
    // ****************************** LOOP OVER GATED IMPLANTS  *****************************
    // *************************************************************************************

    // Define a map to hold a set of  subevents which contained FRS and gated implant coincidence
    std::map<std::string, std::set<int>> gatedimplant_filledtree_map = {};
    
    // Loop over all gatedimplant isotopes and define an empty set
    for ( auto itr : constants::IMPLANT_GATES_INFILE_MAP ){
      // Get the fist element of the tuple containing the name
      std::string implant_gate_name = itr.first;

      // Fill map with empty set
      gatedimplant_filledtree_map.emplace(implant_gate_name, std::set<int>());
    }
    

    // std::cout << "[DEBUG] PreCheck" << std::endl;
    // Implants in coincidence with FRS
    if (aidaimphits > 0 && FrsMultiItem.GetSize()>0 && FrsItem.GetSize()>0) {
      // std::cout << "[DEBUG] PostCheck" << std::endl;

      // Apply Lookahead to check frs entrie sof other frs events
      for (int ilook=current_entry; ilook<=current_entry+constants::LOOKAHEAD; ++ilook){

        // Defensive code for end of entry loop
        if (ilook < 0 || ilook > totalEntries - 1) { continue; }

        // Jump to ilook
        reader.SetEntry(ilook);
        
        // Loop over all FrsMHTDC dataitems
        for (auto const& frs_multi_item : FrsMultiItem) {

          for (auto const& frs_item : FrsItem){

            std::vector<float> const& AoQ = frs_multi_item.Get_ID_AoQ_corr_s2s4_mhtdc();
            std::vector<float> const& Z = frs_multi_item.Get_ID_z41_mhtdc();
            std::vector<float> const& Z2 = frs_multi_item.Get_ID_z42_mhtdc();
            std::vector<float> const& X4 = frs_multi_item.Get_ID_s4x_mhtdc();
            Long64_t const& frs_wr = frs_item.Get_wr_t();
            
            //  Skip frs subevents where implant index wouldbe out of bounds
            for(int i =0; i<AoQ.size(); i++){
              if ( i >= Z.size() || i >= Z2.size() || i >= X4.size() ) { continue; }

              // Get FRS variables
              float aoq = AoQ[i];
              float z = Z[i];
              float z2 = Z2[i];
              float x4 = X4[i];

              // Get aida entries of current entry
              reader.SetEntry(current_entry);

              // Loop over all gatedimplant isotopes
              for (auto gimp_itr : gatedimplant_cuts_map){
                
                // Unpack entries inside map 
                std::string gimp_key = gimp_itr.first;
                auto [x4aoq_cut, zz2_cut] = gimp_itr.second;
                
                if(x4aoq_cut->IsInside(aoq, x4) && zz2_cut->IsInside(z, z2) ){

                  for (int j = 0; j < aidaimphits; j++){

                    // Determine if implant subevent has stopped in specific DSSSD
                    if (implant_dssd[j] == constants::DSSSD && bpflag == 0){ stopped = 1; } 
                    else { stopped = 2; }

                    // Define time diff of coincidence
                    Long64_t time_diff = frs_wr - implant_time[j];
                    
                    //Time gate on the frs event and not relying on the timestitching. Idea: require only one implant for the frs event?
                    if (time_diff < constants::FRS_TIMEWINDOW.first) { continue; }
                    if (time_diff > constants::FRS_TIMEWINDOW.second) { break; }

                    if( gatedimplant_filledtree_map[gimp_key].count(j) == 0 ){
                       
                      // std::cout << "[DEBUG] FOUND an FRS-Implant coincidence.  dT = " << frs_wr-implant_time[j]  << " ###  LookAhead = " << ilook << std::endl;
                       
                      aida_implant_data.time = implant_time[j];
                      aida_implant_data.stopped = stopped;
                      aida_implant_data.x = implant_x[j];
                      aida_implant_data.y = implant_y[j];
                      aida_implant_data.energy = implant_energy[j];
                      aida_implant_data.energy_x = implant_energy_x[j];
                      aida_implant_data.energy_y = implant_energy_y[j];
                      aida_implant_data.cluster_size_x = implant_cluster_size_x[j];
                      aida_implant_data.cluster_size_y = implant_cluster_size_y[j];
                      aida_implant_data.sp = spflag;
                      aida_implant_data.bp = bpflag;

                      gatedimplant_trees_map[gimp_key]->Fill();

                      gatedimplant_counter_map[gimp_key]++;
                      if ( stopped == 1 ){ gatedimplant_stopped_counter_map[gimp_key]++; }

                    }

                    gatedimplant_filledtree_map[gimp_key].insert(j);

                  } // End of aida implant loop

                }

              } // End of gated implant map loop

            } // End of frs subevent loop

          } // End of FRS data item loop

        } // End of FRS MHTDC data item loop
        
        // Change back to ilook for next loop
        reader.SetEntry(ilook);

      } // End of lookahead loop

    }

    // Set entry to current for rest of loop
    reader.SetEntry(current_entry);

    // *************************************************************************************
    // ****************************** LOOP OVER DECAY EVENTS *******************************
    // *************************************************************************************

    // Empty set for each subevents
    std::set<int> aidadecay_filledtree{};

    // Loop over decay subevents
    if(aidadecayhits == 1 && aidaimphits == 0){

      for (int i = 0; i < aidadecayhits; i++) {

        if (TMath::Abs(decay_time[aidadecayhits -1] - decay_time[0]) < 7000) {

          if (aidadecay_filledtree.count(i) == 0) {

            // ####### DEBUG ###### 
            // std::cout << "[DEBUG] Decay time: " << decay_time[i] << " ### Decay time X: " << decay_time_x[i] << " ### Decay time Y: " << decay_time_y[i] << 
            //   " ### Decay Strip X: " << decay_x[i] << " ### Decay Y: " << decay_y[i] << " ### Decay Energy: " << decay_energy[i] << " ### Decay Energy X: " << 
            //   decay_energy_x[i] << " ### Decay Energy Y: " << decay_energy_y[i] << " ### Decay DSSD: " << decay_dssd[i] << " ### SP: " << spflag << 
            //   " ### BP: " << spflag << std::endl << std::endl;
            // ####### DEBUG ###### 

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
            aida_decay_data.dssd = decay_dssd[i];
            aida_decay_data.sp = spflag;
            aida_decay_data.bp = bpflag;
            decay_tree->Fill();
          }
        }

        aidadecay_filledtree.insert(i);

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
    // ****************************** PRINT OUT LOOP PROGRESS ******************************
    // *************************************************************************************
    
    // Show the progress of the loop
    if (reader.GetCurrentEntry() % 100000 == 0) {
      int progress = (reader.GetCurrentEntry() * 100) / totalEntries;
      char spin = spinner[reader.GetCurrentEntry() / 1000 % 4];
      std::cout << "\rProcessing the tree " << reader.GetCurrentEntry() << " (" << progress << "%) " << spin << std::flush;
    }
  
  }

  std::cout << std::endl;

  // *************************************************************************************
  // ****************************** PRINT STATISTICS *************************************
  // *************************************************************************************

  // Loop over counter map and print statistics
  for (auto itr : gatedimplant_counter_map){
    std::string key = itr.first;
    std::cout << "Number of " << key << " implants: " << gatedimplant_counter_map[key] << " ### Stopped: " << gatedimplant_stopped_counter_map[key] << " ### Stopping Efc.: " << gatedimplant_stopped_counter_map[key]*100/gatedimplant_counter_map[key] << " %" << std::endl;
  }
  
  // std::cout << std::endl << std::endl;

  // *************************************************************************************
  // ****************************** WRITE TREES ******************************************
  // *************************************************************************************

  std::cout << "Writing trees ..." << std::endl;
  outputFile->cd();
  
  implant_tree->Write();
  decay_tree->Write();
  germanium_tree->Write();
  
  // Loop over trees and Write
  for (auto itr : gatedimplant_trees_map){
    
    // Extract tree and write
    itr.second->Write();
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
