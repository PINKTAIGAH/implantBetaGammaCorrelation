#include<iostream>
#include<fstream>
#include<sstream>
#include<map>
#include<set>
#include<unordered_map>
#include<tuple>
#include<utility>
#include<string>
#include<vector>

#include<TH1F.h>
#include<TH2F.h>
#include <TCutG.h>
#include<TFile.h>
#include<TTree.h>
#include<THStack.h>
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

  const std::string BROKEN_STRIPS_INFILE = "/lustre/gamma/jeroen/S100/analysis/betaion/AIDA_strips.txt";

  const int TOTAL_GERMANIUM_DETECTORS = 15; // Number of germanium detectors 
  const int TOTAL_GERMANIUM_CRYSTALS = 2; // Number of germanium crystals
  const int TOTAL_FEES = 16; // Number of FEES in AIDA DAQ

  const std::map<std::string, std::string> IMPLANT_GATES_INFILE_MAP = {
    {"82nb", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/82Nb.root"},
    {"84nb", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/84Nb.root"},
    {"84mo", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/84Mo.root"},
    {"85mo", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/85Mo.root"}
  };

  const int DSSD = 1;

}

namespace histogram_resolutions{
  const int DATA_ITEM_RATES = 268e6; // 268ms per bin for HEC/LEC dataitems
  const int EVENT_RATES = 262e3; // 262us per bin for Implant/Decay events 
}

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


void eventwiseHists(const char* input, const char* output) {

  // *************************************************************************************
  // ****************************** OPEN & CREATE FILES **********************************
  // *************************************************************************************
  
  std::cout << std::endl << "Opening files ..." << std::endl;
  
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

  std::cout << "Files opened succesfully!" << std::endl << std::endl;

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
    TCutG *zaoq_cut = (TCutG*)gate_file->Get("cut_Z_AoQ");
    TCutG *zz2_cut = (TCutG*)gate_file->Get("cut_Z_Z2");
    
    // Add gate information to multimap
    gatedimplant_cuts_map.emplace( gate_name, std::make_tuple(zaoq_cut, zz2_cut) );   
  }

  // Create a map which will hold the counters for each gated implant species
  std::map<std::string, int> gatedimplant_counter_map = {};
  
  // Loop over map and fill with zerod counter for each gate
  for (auto itr : constants::IMPLANT_GATES_INFILE_MAP){

    std::string gate_name = itr.first; // Get gate name

    // Fill map
    gatedimplant_counter_map.emplace( gate_name, 0);
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

  TTreeReaderArray<Int_t> implant_adc_dssd(reader, "AidaImplantCalAdcData.dssd");
  TTreeReaderArray<Int_t> implant_adc_side(reader, "AidaImplantCalAdcData.side");
  TTreeReaderArray<Int_t> implant_adc_strip(reader, "AidaImplantCalAdcData.strip");
  TTreeReaderArray<Int_t> implant_adc_fee(reader, "AidaImplantCalAdcData.fee");
  TTreeReaderArray<Double_t> implant_adc_energy(reader, "AidaImplantCalAdcData.energy");
  TTreeReaderArray<unsigned long long> implant_adc_time(reader, "AidaImplantCalAdcData.slowTime");

  TTreeReaderArray<ULong_t> bplast_time(reader, "bPlastTwinpeaksCalData.fwr_t");
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
  TTreeReaderArray<Float_t> frs_x4(reader, "FrsHitData.fID_x4");
  TTreeReaderArray<Float_t> frs_x2(reader, "FrsHitData.fID_x2");
  TTreeReaderArray<Float_t> frs_z(reader, "FrsHitData.fID_z41");
  TTreeReaderArray<Float_t> frs_z2(reader, "FrsHitData.fID_z42");
  TTreeReaderArray<Float_t> frs_aoq(reader, "FrsHitData.fID_AoQ_corr_s2s4");
  TTreeReaderArray<long long> frs_time(reader, "FrsHitData.fwr_t");
  TTreeReaderArray<Float_t> frs_dedeg(reader, "FrsHitData.fID_dEdeg_z41");


  // *************************************************************************************
  // ****************************** DEFINE OUTPUT FILE **********************************
  // *************************************************************************************

  // Open the output file
  TFile* outputFile = new TFile(output, "RECREATE");
  if (!outputFile) {
      std::cerr << "Error: Could not create output file " << std::endl;
      return;
  }

  // *************************************************************************************
  // ****************************** FIND WR TIMES OF EXPERIMENT **************************
  // *************************************************************************************

  // Initialise time variables
  ULong_t wr_start = 1.7401830e+18;
  ULong_t wr_end = 1.74022529e+18;

  /*std::cout << "Finding file white rabbit times ..." << std::endl;*/
  /**/
  /*// Loop over all the entries in the new tree*/
  /*while (reader.Next()){*/
  /**/
  /*  int frshits = frs_time.GetSize();*/
  /**/
  /*  if (frshits>0 && wr_start==0 ){ wr_start = frs_time[0]; }*/
  /*  if ( frshits>0 ){ wr_end = frs_time[0]; }*/
  /*}*/
  /**/
  /*reader.Restart(); // Reset reader for next loop*/
  /**/
  /*// Find white rabbit time difference of experiments*/
  ULong_t wr_dt = wr_end - wr_start;

  std::cout << "Found file white rabbit times!" << std::endl;
  std::cout << "WR Start time: " << wr_start << "  #####    WR End time: " << wr_end << std::endl;  
  std::cout << "File run time: " << wr_dt * 1e-9 / 60 << " mins" << std::endl << std::endl; 
  
  // *************************************************************************************
  // ****************************** DEFINE OUTPUT HISTOGRAMS **********************************
  // *************************************************************************************

  // Define binning width for specific plots using file runtime
  int data_item_binwidth = wr_dt / histogram_resolutions::DATA_ITEM_RATES;
  int event_binwidth = wr_dt / histogram_resolutions::EVENT_RATES;
  /*int lec_energy_binwidth = 3000 / histogram_resolutions::LEC_ENERGY; */
  /*int hec_energy_binwidth = 7000 / histogram_resolutions::HEC_ENERGY; */
  
  // Initialise histogram object for 16 histograms for each FEEs of LEC and HEC data items 
  std::map<int, TH1F*> h1_aida_lec_dataitem_all_map = {};
  std::map<int, TH1F*> h1_aida_lec_dataitem_lowbound_map = {};
  std::map<int, TH1F*> h1_aida_lec_dataitem_highbound_map = {};

  std::map<int, TH1F*> h1_aida_hec_dataitem_all_map = {};
  std::map<int, TH1F*> h1_aida_hec_dataitem_lowbound_map = {};
  std::map<int, TH1F*> h1_aida_hec_dataitem_highbound_map = {};

  // Loop over all fee and fill maps
  for(int fee = 1; fee<17; fee++){
    // Create hist object 
    TH1F* h1_aida_lec_dataitem_all = new TH1F( Form("aida_lec_dataitem_all_fee%d", fee), Form("AIDA FEE %d LEC Dataitem Rate (All)", fee), data_item_binwidth, wr_start, wr_end); 
    TH1F* h1_aida_lec_dataitem_lowbound = new TH1F( Form("aida_lec_dataitem_lowbound_fee%d", fee), Form("AIDA FEE %d LEC Dataitem Rate (150 keV #leq E #leq 1500 keV)", fee), data_item_binwidth, wr_start, wr_end); 
    TH1F* h1_aida_lec_dataitem_highbound = new TH1F( Form("aida_lec_dataitem_highbound_fee%d", fee), Form("AIDA FEE %d LEC Dataitem Rate (E #geq 1500 keV)", fee), data_item_binwidth, wr_start, wr_end); 

    TH1F* h1_aida_hec_dataitem_all = new TH1F( Form("aida_hec_dataitem_all_fee%d", fee), Form("AIDA FEE %d HEC Dataitem Rate (All)", fee), data_item_binwidth, wr_start, wr_end); 
    TH1F* h1_aida_hec_dataitem_lowbound = new TH1F( Form("aida_hec_dataitem_lowbound_fee%d", fee), Form("AIDA FEE %d HEC Dataitem Rate (100 MeV #leq E #leq 1000 MeV)", fee), data_item_binwidth, wr_start, wr_end); 
    TH1F* h1_aida_hec_dataitem_highbound = new TH1F( Form("aida_hec_dataitem_highbound_fee%d", fee), Form("AIDA FEE %d HEC Dataitem Rate (E #geq 1000 NeV)", fee), data_item_binwidth, wr_start, wr_end); 

    // Fill Maps
    h1_aida_lec_dataitem_all_map.emplace(fee, h1_aida_lec_dataitem_all);
    h1_aida_lec_dataitem_lowbound_map.emplace(fee, h1_aida_lec_dataitem_lowbound);
    h1_aida_lec_dataitem_highbound_map.emplace(fee, h1_aida_lec_dataitem_highbound);
    h1_aida_hec_dataitem_all_map.emplace(fee, h1_aida_hec_dataitem_all);
    h1_aida_hec_dataitem_lowbound_map.emplace(fee, h1_aida_hec_dataitem_lowbound);
    h1_aida_hec_dataitem_highbound_map.emplace(fee, h1_aida_hec_dataitem_highbound);
  }

  // Make LEC and HEC event and multiplicity histograms
  TH2F* h2_aida_lec_xy_multiplicity = new TH2F("aida_lec_xy_multiplicity", "AIDA LEC XY Strip Multiplicity", 500, 0, 500, 500, 0, 500);
  TH2F* h2_aida_hec_xy_multiplicity = new TH2F("aida_hec_xy_multiplicity", "AIDA HEC XY Strip Multiplicity", 50, 0, 50, 50, 0, 50);

  // Bplast Channel Multiplicity per implant event 
  TH1F* h1_aida_implant_bplast_multiplicity = new TH1F("aida_implant_bplast_multiplicity", "AIDA Implant Event - Bplast Channel Multiplicity (Upstream + Downstream)", 140, 0, 140); 
  TH2F* h2_aida_implant_energy_bplast_multiplicity = new TH2F("aida_implant_energy_bplast_multiplicity", "AIDA Implant Event Energy - Bplast Channel Multiplicity (Upstream + Downstream)", 7000/20, 0, 7000, 140, 0, 140); 

  // Bplast Channel Multiplicity per gatedimplant event 
  std::map<std::string, TH1F*> h1_aida_gatedimplant_bplast_multiplicity_map = {};
  std::map<std::string, TH2F*> h2_aida_gatedimplant_energy_bplast_multiplicity_map = {};
  
  // Loop over all gated implants and fill map
  for (auto itr : constants::IMPLANT_GATES_INFILE_MAP ){
    
    // Make histogram for Bplast channel multiplicity
    std::string gatedimplant_name = itr.first;
    std::string histogram_name = "aida_" + gatedimplant_name + "_implant_bplast_multiplicity";
    std::string histogram_title = "AIDA " + gatedimplant_name + " Implant Event - Bplast Channel Multiplicity (Upstream + Downstream)";
    TH1F* h1_aida_gatedimplant_bplast_multiplicity = new TH1F(histogram_name.c_str(), histogram_title.c_str(), 140, 0, 140);
    h1_aida_gatedimplant_bplast_multiplicity_map.emplace(gatedimplant_name, h1_aida_gatedimplant_bplast_multiplicity);

    // Make histogram for gated impalnt vs bplas multiplicity
    histogram_name = "aida_" + gatedimplant_name + "energy_implant_bplast_multiplicity";
    histogram_title = "AIDA " + gatedimplant_name + " Implant Event Energy - Bplast Channel Multiplicity (Upstream + Downstream)";
    TH2F* h2_aida_gatedimplant_energy_bplast_multiplicity = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 7000/20, 0, 7000, 140, 0, 140);
    h2_aida_gatedimplant_energy_bplast_multiplicity_map.emplace(gatedimplant_name, h2_aida_gatedimplant_energy_bplast_multiplicity);

    // Make histogram of energy distribution for each gated implant
  }

  // FRS histograms
  TH2F* h2_frs_z_z2 = new TH2F("frs_z_z2", "FRS Z vs Z2", 1000, 30, 50, 1000, 30, 50);
  TH2F* h2_frs_z_aoq = new TH2F("frs_z_aoq", "FRS Z vs AoQ", 1000, 1.8,2.5, 1000, 30, 50);
  TH2F* h2_frs_z_aoq_corr = new TH2F("frs_z_aoq_corr", "FRS Z vs AoQ", 1000, 1.8,2.5, 1000, 30, 50);
  TH2F* h2_frs_aoq_sx4 = new TH2F("frs_aoq_sx4", "FRS AoQ vs X4", 1000, -150, 150, 1000, 2.0, 2.5);

  // Decay Hits multiplicity
  TH1F* h1_aida_decay_multiplicity_allspill = new TH1F("aida_decay_multiplicity_allspill", "AIDA decay sub-event multiplicity per event (Allspill)", 50, 0, 50);
  TH1F* h1_aida_decay_multiplicity_onspill = new TH1F("aida_decay_multiplicity_onspill", "AIDA decay sub-event multiplicity per event (Onspill)", 50, 0, 50);
  TH1F* h1_aida_decay_multiplicity_offspill = new TH1F("aida_decay_multiplicity_offspill", "AIDA decay sub-event multiplicity per event (Offspill)", 50, 0, 50);
  THStack* sh1_aida_decay_multiplicity = new THStack("aida_decay_multiplicity", "AIDA decay sub-event multiplicity per event");

  // Plot Decay cluster sizes
  TH1F* h1_aida_decay_cluster_size_x =new TH1F("aida_decay_cluster_size_x", "AIDA Decay X Strip Cluster Size", 10, 0, 10);
  TH1F* h1_aida_decay_cluster_size_y =new TH1F("aida_decay_cluster_size_y", "AIDA Decay Y Strip Cluster Size", 10, 0, 10);
  TH1F* h1_aida_decay_cluster_size =new TH1F("aida_decay_cluster_size", "AIDA Decay Y Strip Cluster Size", 100, 0, 10);

  // Plot Germanium spectra
  TH1F* h1_germanium_spectra_all = new TH1F("germanium_spectra_all", "Gamma spectra (All summed)", 1500, 0, 1500);

  // Plot Decay-Germanium coincidence
  TH1F* h1_decay_germanium_dt = new TH1F("decay_germanium_dt", "AIDA Decay - DEGAS time difference", 50000, -100e3, 100e3);

  // *************************************************************************************
  // ****************************** LOOP OVER ALL EVENTS *********************************
  // *************************************************************************************

  const char spinner[] = {'-', '\\', '|', '/'};
  int totalEntries = reader.GetEntries(true);


  // Loop over all entries in the old tree
  while (reader.Next()) {

    // *************************************************************************************
    // ****************************** GET EVENT INFORMATION *********************************
    // *************************************************************************************

    // sizes of events
    int germaniumhits = germanium_time.GetSize();
    int aidadecayhits = decay_time.GetSize();
    int aidaimphits = implant_time.GetSize();
    int aidaadcimphits = implant_adc_dssd.GetSize();
    int aidaadcdecayhits = decay_adc_dssd.GetSize();
    int bplasthits = bplast_id.GetSize();
    int frshits = frs_time.GetSize();
      
    // Initialisation of flags
    int spflag = 0;
    int bpflag = 0;
    int bp1flag = 0;
    int bp2flag = 0;

    // Initialisation of multiplicities
    int lec_x_multiplicity = 0;
    int lec_y_multiplicity = 0;
    int hec_x_multiplicity = 0;
    int hec_y_multiplicity = 0;

    int bplast_upstream_channel_multiplicity = 0;
    int bplast_downstream_channel_multiplicity = 0;

    // Define spill flag for event
    if(*spill == true) spflag = 1;
    if(*spill == false) spflag = 2;
      

    // *************************************************************************************
    // ****************************** LOOP OVER FRS  ************************************
    // *************************************************************************************

    // if(frshits>0){

    //   // Loop over FRS hits
    //   for(int j=0; j<frshits; j++){
        
    //     h2_frs_z_z2->Fill(frs_z[j], frs_z2[j]);
    //     h2_frs_z_aoq->Fill(frs_z[j], frs_aoq[j]);
    //     h2_frs_aoq_x4->Fill(frs_aoq[j], frs_x4[j]);
    //     h2_frs_dedeg_z->Fill(frs_z[j],frs_dedeg[j]);

    //   } // End of FRS loop

    // }

    for ( const auto& FrsMultiItem : FrsMultiItem ){

      // Define frs item vectors
      std::vector<float> const& AoQ = FrsMultiItem.Get_ID_AoQ_s2s4_mhtdc();
      std::vector<float> const& AoQ_CORR = FrsMultiItem.Get_ID_AoQ_corr_s2s4_mhtdc();
      std::vector<float> const& Z41 = FrsMultiItem.Get_ID_z41_mhtdc();
      std::vector<float> const& Z42 = FrsMultiItem.Get_ID_z42_mhtdc();
      std::vector<float> const& S4X = FrsMultiItem.Get_ID_s4x_mhtdc();

      // Loop over variables in each frs item vector
      for ( int i = 0; i < AoQ.size(); i++ ){

        if ( i >= AoQ.size() || i >= AoQ_CORR.size() || i >= Z41.size() || i >= Z42.size() || i >= S4X.size() ){ continue; }

        // Get frs variables
        double aoq = AoQ[i];
        double aoq_corr = AoQ_CORR[i];
        double z41 = Z41[i];
        double z42 = Z42[i];
        double s4x = S4X[i];

        if ( aoq == -999. || aoq_corr == -999. || z41 == -999. || z42 == -999. || s4x == -000. ){ continue; }

        //************* DEBUG **************
        // std::cout << "[DEBUG] FRS data items --> AoQ: " << aoq << "   #####   AoQ_CORR: " << aoq_corr << "   #####   Z41: " << z41 << "   #####   Z42: " << z42 << "   #####   S4X: " << s4x << std::endl;
        //************* DEBUG **************

        // Fill histograms
        h2_frs_z_z2->Fill(z41, z42);
        h2_frs_z_aoq->Fill(aoq, z42);
        h2_frs_z_aoq_corr->Fill(aoq_corr, z42);
        h2_frs_aoq_sx4->Fill(s4x, aoq_corr);

      }

    }


    // *************************************************************************************
    // ****************************** LOOP OVER BPLAST  ************************************
    // *************************************************************************************
      
    // Loop over bplast subevents to get flags and multiplicities
    for (int j = 0; j < bplasthits; j++) {
      
      // If in upstream bplast
      if(bplast_id[j] < 64){
        bp1flag = 1;
        bplast_upstream_channel_multiplicity++;
      }

      // If in downstream bplast 
      if(bplast_id[j] > 63 && bplast_id[j] < 128){
        bp2flag = 1;
        bplast_downstream_channel_multiplicity++;
      }

    }
  
    // Determine bplast event flags
    if (bp1flag == 1 && bp2flag == 0) bpflag = 1;
    if (bp1flag == 0 && bp2flag == 1) bpflag = 2;
    if (bp1flag == 1 && bp2flag == 1) bpflag = 3;



    // *************************************************************************************
    // ****************************** LOOP OVER AIDA IMPLANTS  *****************************
    // *************************************************************************************

    if(aidaimphits > 0){

      for (int j = 0; j < aidaimphits; j++) {
        
        // Check if subevent is in the correct DSSD
        if(implant_dssd[j] == constants::DSSD){

          // Fill bplast channel multiplicities   
          h1_aida_implant_bplast_multiplicity->Fill(bplast_upstream_channel_multiplicity); 
          h1_aida_implant_bplast_multiplicity->Fill(bplast_downstream_channel_multiplicity + 70); 
          
          // Fill bplast channel multiplicity vs implant energy
          h2_aida_implant_energy_bplast_multiplicity->Fill(implant_energy[j], bplast_upstream_channel_multiplicity); 
          h2_aida_implant_energy_bplast_multiplicity->Fill(implant_energy[j], bplast_downstream_channel_multiplicity + 70); 
      
        }

      } // End of implant loop
    
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
    

    // Implants in coincidence with FRS
    if (aidaimphits > 0) {
      
      // Loop over all FrsMHTDC dataitems
      for (auto const& FrsMultiItem : FrsMultiItem) {

        std::vector<float> const& AoQ = FrsMultiItem.Get_ID_AoQ_corr_s2s4_mhtdc();
        std::vector<float> const& Z = FrsMultiItem.Get_ID_z41_mhtdc();
        std::vector<float> const& Z2 = FrsMultiItem.Get_ID_z42_mhtdc();
        
        //  Skip frs subevents where implant index wouldbe out of bounds
        for(int i =0; i<AoQ.size(); i++){
          if (i >= AoQ.size() || i >= Z.size()) { continue; }

          // Get FRS variables
          double aoq = AoQ[i];
          double z = Z[i];
          double z2 = Z[i];

          // Loop over all gatedimplant isotopes
          for (auto gimp_itr : gatedimplant_cuts_map){
            
            // Unpack entries inside map 
            std::string gimp_key = gimp_itr.first;
            auto [zaoq_cut, zz2_cut] = gimp_itr.second;
            
            // Check if inside gate
            if( zaoq_cut->IsInside(aoq, z) && zz2_cut->IsInside(z, z2) ){
              for (int j = 0; j < aidaimphits; j++){

                // Check if correct DSSD and if same subevent hasnt been filled
                if( implant_dssd[j] == constants::DSSD && gatedimplant_filledtree_map[gimp_key].count(j) == 0 ){
                   
                  // Fill bplast channel multiplicity for each gatedimplant isotope     
                  h1_aida_gatedimplant_bplast_multiplicity_map[gimp_key]->Fill(bplast_upstream_channel_multiplicity); // Upstream
                  h1_aida_gatedimplant_bplast_multiplicity_map[gimp_key]->Fill(bplast_downstream_channel_multiplicity + 70); // Downstream

                  // Fill bplast channel multiplicity vs implant energy
                  h2_aida_gatedimplant_energy_bplast_multiplicity_map[gimp_key]->Fill(implant_energy[j], bplast_upstream_channel_multiplicity); // Upstream
                  h2_aida_gatedimplant_energy_bplast_multiplicity_map[gimp_key]->Fill(implant_energy[j], bplast_downstream_channel_multiplicity + 70); // Downstream

                  gatedimplant_counter_map[gimp_key]++; // Increase counter

                }

                gatedimplant_filledtree_map[gimp_key].insert(j);

              } // End of aida implant loop

            }

          } // End of gated implant map loop

        } // End of frs subevent loop

      } // End of FRS data item loop

    }

    
    // *************************************************************************************
    // ****************************** LOOP OVER AIDA ADC HEC  ******************************
    // *************************************************************************************
    
    // Loop over all HEC data items
    if(aidaadcimphits>0){

      for (int j = 0; j < aidaadcimphits; j++){
          
        // Check if correct DSSSD
        if (implant_adc_dssd[j] == constants::DSSD){

          // Increase multiplicity
          if (implant_adc_side[j] == 1){ hec_x_multiplicity++; }
          if (implant_adc_side[j] == -1){ hec_y_multiplicity++; }

          // Fill Data item rate histograms
          if (implant_adc_energy[j] >= 100 && implant_adc_energy[j] <= 1000){ h1_aida_hec_dataitem_lowbound_map[implant_adc_fee[j]]->Fill(implant_adc_time[j]); }
          if (implant_adc_energy[j] > 1000){ h1_aida_hec_dataitem_highbound_map[implant_adc_fee[j]]->Fill(implant_adc_time[j]); }
          h1_aida_hec_dataitem_all_map[implant_adc_fee[j]]->Fill(implant_adc_time[j]);
        }

      } // End of HEC data items loop
      
      // Fill multiplicity histogram
      h2_aida_hec_xy_multiplicity->Fill(hec_x_multiplicity, hec_y_multiplicity);
      
    }
    

    // *************************************************************************************
    // ****************************** LOOP OVER DECAY EVENTS *******************************
    // *************************************************************************************

    // Empty set for each subevents*/
    std::set<int> aidadecay_filledtree{};

    int decay_multiplicity = 0;
    
    // Loop over decay subevents
    if(aidadecayhits > 0){
    
      for (int i = 0; i < aidadecayhits; i++) {
    
        if (TMath::Abs(decay_time[aidadecayhits -1] - decay_time[0]) < 33000) {
    
          if (aidadecay_filledtree.count(i) == 0) {

            decay_multiplicity++;

            if ( aidadecayhits==1 && TMath::Abs(decay_time_x[i]-decay_time_y[i])<1e3 && TMath::Abs(decay_energy_x[i]-decay_energy_y[i])<168 && decay_energy[i]>150 && decay_energy[i]<1000 ){

              h1_aida_decay_cluster_size_x->Fill(decay_cluster_size_x[i]);
              h1_aida_decay_cluster_size_y->Fill(decay_cluster_size_y[i]);
              h1_aida_decay_cluster_size->Fill(TMath::Sqrt(decay_cluster_size_x[i]+decay_cluster_size_y[i]));

            }
    
          }

        }
    
      } // End of decay loop
    
    }

    // Fill multiplicity histogram 
    if(decay_multiplicity>0) {h1_aida_decay_multiplicity_allspill->Fill(decay_multiplicity);}
    if(decay_multiplicity>0 && spflag==1) {h1_aida_decay_multiplicity_onspill->Fill(decay_multiplicity);}
    if(decay_multiplicity>0 && spflag==2) {h1_aida_decay_multiplicity_offspill->Fill(decay_multiplicity);}


    // *************************************************************************************
    // ****************************** LOOP OVER AIDA ADC LEC  ******************************
    // *************************************************************************************
    
    // Loop over all LEC data items
    if(aidaadcdecayhits>0){

      for (int j = 0; j < aidaadcdecayhits; j++){
          
        // Check if correct DSSSD
        if (decay_adc_dssd[j] == constants::DSSD){

          // Increase multiplicity
          if (decay_adc_side[j] == 1){ lec_x_multiplicity++; }
          if (decay_adc_side[j] == -1){ lec_y_multiplicity++; }

          // Fill Data item rate histograms
          if (decay_adc_energy[j] >= 150 && decay_adc_energy[j] <= 1500){ h1_aida_lec_dataitem_lowbound_map[decay_adc_fee[j]]->Fill(decay_adc_time[j]); }
          if (decay_adc_energy[j] > 1500){ h1_aida_lec_dataitem_highbound_map[decay_adc_fee[j]]->Fill(decay_adc_time[j]); }
          h1_aida_lec_dataitem_all_map[decay_adc_fee[j]]->Fill(decay_adc_time[j]);
        }

      } // End of LEC data items loop
      
      // Fill multiplicity histogram
      h2_aida_lec_xy_multiplicity->Fill(lec_x_multiplicity, lec_y_multiplicity);
      
    }
    

    // *************************************************************************************
    // ****************************** LOOP OVER GERMANIUMS EVENTS **************************
    // *************************************************************************************

    
    std::set<int> germanium_filledtree{};
    
    // Loop over decay subevents*/
    for (int j = 0; j < germaniumhits; j++){
    
      if ( germanium_det[j] <= constants::TOTAL_GERMANIUM_DETECTORS && germanium_cry[j] <= constants::TOTAL_GERMANIUM_CRYSTALS && germanium_energy[j] > 0 ){
    
        if ( germanium_filledtree.count(j) == 0 ){

          h1_germanium_spectra_all->Fill(germanium_energy[j]);

        }
    
      germanium_filledtree.insert(j);
    
      }
    
    }


    // *************************************************************************************
    // ****************************** DEGAS - AIDA DECAY COINCIDENCE ***********************
    // *************************************************************************************

    if (aidadecayhits > 0 && germaniumhits > 0 ){
      
      for (int i = 0; i < aidadecayhits; i++){

        if ( decay_dssd[i]==1 && TMath::Abs(decay_time[aidadecayhits -1] - decay_time[0]) < 33000 && TMath::Abs(decay_time_y[i]-decay_time_x[i])<1e3 && TMath::Abs(decay_energy_x[i]-decay_energy_x[i])<168 && decay_energy[i]>150 & decay_energy[i]<1000 ){

          for (int j = 0; j < germaniumhits; j++){

            if (germanium_det[j] <= 15 && germanium_cry[j] <= 2 && germanium_energy[j] > 0 ){

              //************* DEBUG **************
              // std::cout << "[DEBUG] decay_t - germ_t = " << (decay_time[i] - germanium_time[j])*1e-9 << " seconds  ####### decay_t - germ_abs_evt_t = " << (germanium_time[j] - germanium_abs_ev_time[j])*1e-9 << " seconds" << std::endl;
              //************* DEBUG **************

              h1_decay_germanium_dt->Fill(decay_time[i]-germanium_time[j]);

            }

          }

        }

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
    std::cout << "Number of " << itr.first << " implants: " << itr.second << std::endl;
  }
  

  // *************************************************************************************
  // ****************************** WRITE HISTOGRAMS *************************************
  // *************************************************************************************
  
  // Write FRS histograms
  TDirectory* frsDir = outputFile->mkdir("frs");
  frsDir->cd();
  h2_frs_z_aoq->Write();
  h2_frs_z_aoq_corr->Write();
  h2_frs_z_z2->Write();
  h2_frs_aoq_sx4->Write();
  gFile->cd();

  // Write LEC dataitem All histograms
  TDirectory* lecAllDir = outputFile->mkdir("lec_all");
  lecAllDir->cd();
  for (auto itr : h1_aida_lec_dataitem_all_map){
    itr.second->Write();
  }
  gFile->cd();

  // Write LEC dataitem lowbound histograms
  TDirectory* lecLowboundDir = outputFile->mkdir("lec_lowbound");
  lecLowboundDir->cd();
  for (auto itr : h1_aida_lec_dataitem_lowbound_map){
    itr.second->Write();
  }
  gFile->cd();

  // Write LEC dataitem highbound histograms
  TDirectory* lecHighboundDir = outputFile->mkdir("lec_highbound");
  lecHighboundDir->cd();
  for (auto itr : h1_aida_lec_dataitem_highbound_map){
    itr.second->Write();
  }
  gFile->cd();

  // Write HEC dataitem All histograms
  TDirectory* hecAllDir = outputFile->mkdir("hec_all");
  hecAllDir->cd();
  for (auto itr : h1_aida_hec_dataitem_all_map){
    itr.second->Write();
  }
  gFile->cd();

  // Write HEC dataitem lowbound histograms
  TDirectory* hecLowboundDir = outputFile->mkdir("hec_lowbound");
  hecLowboundDir->cd();
  for (auto itr : h1_aida_hec_dataitem_lowbound_map){
    itr.second->Write();
  }
  gFile->cd();

  // Write HEC dataitem highbound histograms
  TDirectory* hecHighboundDir = outputFile->mkdir("hec_highbound");
  hecHighboundDir->cd();
  for (auto itr : h1_aida_hec_dataitem_highbound_map){
    itr.second->Write();
  }
  gFile->cd();
  
  // Write dataitem multiplicities
  h2_aida_lec_xy_multiplicity->Write();
  h2_aida_hec_xy_multiplicity->Write();

  // Write implant bplast channel multiplicities
  h1_aida_implant_bplast_multiplicity->Write();
  for (auto itr : h1_aida_gatedimplant_bplast_multiplicity_map){
    itr.second->Write();
  }

  // Write implant energy vs bplast channel multiplicities
  h2_aida_implant_energy_bplast_multiplicity->Write();
  for (auto itr : h2_aida_gatedimplant_energy_bplast_multiplicity_map){
    itr.second->Write();
  }

  // Write decay event multiplicity
  h1_aida_decay_multiplicity_allspill->Write();

  h1_aida_decay_multiplicity_onspill->SetLineColor(kBlue);
  h1_aida_decay_multiplicity_offspill->SetLineColor(kRed);
  sh1_aida_decay_multiplicity->Add(h1_aida_decay_multiplicity_onspill);
  sh1_aida_decay_multiplicity->Add(h1_aida_decay_multiplicity_offspill);
  sh1_aida_decay_multiplicity->Write();

  // Write decay cluster sizes
  h1_aida_decay_cluster_size_x->Write();
  h1_aida_decay_cluster_size_y->Write();
  h1_aida_decay_cluster_size->Write();

  // Write gamma plots
  h1_germanium_spectra_all->Write();
  h1_decay_germanium_dt->Write();

  // *************************************************************************************
  // ****************************** CLEANUP **********************************************
  // *************************************************************************************

  // Close the file
  file->Close();
  outputFile->Close();
  delete file;
  delete outputFile;
}
