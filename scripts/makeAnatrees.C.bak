#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

namespace constants{
  const bool DRAW_HISTS = false;
  const bool MAKE_ANATREES = true;

  const int DSSD = 1;
}

void makeAnatrees_old(const char* input, const char* output) {
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

  // Create a tree reader
  TTreeReader reader(old_tree);

  // Create tree reader values for the branches you want to read
  TTreeReaderValue<bool> spill(reader, "EventHeader.fSpillFlag");
  TTreeReaderArray<int64_t> implant_time(reader, "AidaImplantHits.Time");
  TTreeReaderArray<Int_t> implant_dssd(reader, "AidaImplantHits.DSSD");
  TTreeReaderArray<Bool_t> implant_stopped(reader, "AidaImplantHits.Stopped");
  TTreeReaderArray<double> implant_energy(reader, "AidaImplantHits.Energy");
  TTreeReaderArray<double> implant_energy_x(reader, "AidaImplantHits.EnergyX");
  TTreeReaderArray<double> implant_energy_y(reader, "AidaImplantHits.EnergyY");
  TTreeReaderArray<Double_t> implant_x(reader, "AidaImplantHits.StripX");
  TTreeReaderArray<Double_t> implant_y(reader, "AidaImplantHits.StripY");

  // // TTreeReaderArray<u_int64_t> bplast_time(reader, "bPlastTwinpeaksCalData.fabsolute_event_time");
  // // TTreeReaderArray<ULong_t> bplast_time(reader, "bPlastTwinpeaksCalData.fwr_t");
  // TTreeReaderArray<ushort> bplast_id(reader, "bPlastTwinpeaksCalData.fdetector_id");
  // TTreeReaderArray<Double_t> bplast_slow_tot(reader, "bPlastTwinpeaksCalData.fslow_ToT");

  TTreeReaderArray<int64_t> decay_time(reader, "AidaDecayHits.Time");
  TTreeReaderArray<Int_t> decay_dssd(reader, "AidaDecayHits.DSSD");
  TTreeReaderArray<int64_t> decay_time_x(reader, "AidaDecayHits.TimeX");
  TTreeReaderArray<int64_t> decay_time_y(reader, "AidaDecayHits.TimeY");
  TTreeReaderArray<Int_t> decay_adc_dssd(reader, "AidaDecayCalAdcData.dssd");
  TTreeReaderArray<Int_t> decay_adc_side(reader, "AidaDecayCalAdcData.side");
  TTreeReaderArray<Int_t> decay_adc_strip(reader, "AidaDecayCalAdcData.strip");
  TTreeReaderArray<Double_t> decay_energy(reader, "AidaDecayHits.Energy");
  TTreeReaderArray<Double_t> decay_energy_x(reader, "AidaDecayHits.EnergyX");
  TTreeReaderArray<Double_t> decay_energy_y(reader, "AidaDecayHits.EnergyY");
  TTreeReaderArray<Double_t> decay_x(reader, "AidaDecayHits.StripX");
  TTreeReaderArray<Double_t> decay_y(reader, "AidaDecayHits.StripY");
  TTreeReaderArray<Int_t> decay_cluster_size_x(reader, "AidaDecayHits.ClusterSizeX");
  TTreeReaderArray<Int_t> decay_cluster_size_y(reader, "AidaDecayHits.ClusterSizeY");

  TTreeReaderArray<uint64_t> germanium_time(reader, "GermaniumCalData.fwr_t");
  TTreeReaderArray<int64_t> germanium_abs_ev_time(reader, "GermaniumCalData.fabsolute_event_time");
  TTreeReaderArray<uint> germanium_det(reader, "GermaniumCalData.fdetector_id");
  TTreeReaderArray<uint> germanium_cry(reader, "GermaniumCalData.fcrystal_id");
  TTreeReaderArray<Double_t> germanium_energy(reader, "GermaniumCalData.fchannel_energy");

  TTreeReaderArray<FrsMultiHitItem> FrsMultiItem(reader, "FrsMultiHitData");
  // TTreeReaderArray<Float_t> frs_x4(reader, "FrsHitData.fID_x4");
  // TTreeReaderArray<Float_t> frs_x2(reader, "FrsHitData.fID_x2");
  // TTreeReaderArray<Float_t> frs_z(reader, "FrsHitData.fID_z");
  // TTreeReaderArray<Float_t> frs_z2(reader, "FrsHitData.fID_z2");
  // TTreeReaderArray<Float_t> frs_aoq(reader, "FrsHitData.fID_AoQ_corr");
  // TTreeReaderArray<uint64_t> frs_time(reader, "FrsHitData.fwr_t");
  // TTreeReaderArray<Float_t> frs_dedeg(reader, "FrsHitData.fID_dEdeg");

  uint64_t wr_experiment_start = 1.7137464e+18;
  uint64_t wr_experiment_end = 1.7143848e+18;
  int64_t duration_in_seconds = (wr_experiment_end - wr_experiment_start)/1e9;
  int64_t slices_every = 1; //s
  int64_t number_of_slices = duration_in_seconds/slices_every;

  int implanted_82Nb = 0;
  int implanted_84Nb = 0;
  int implanted_84Mo = 0;
  int implanted_85Mo = 0;

  // Load broken strips

  // std::ifstream aida_strips_file("/lustre/gamma/jeroen/S100/analysis/betaion/AIDA_strips.txt");
  
  // std::vector<int> broken_xstrips;
  // std::vector<int> broken_ystrips;
  // int dssd_number;
  // std::string xy;
  // int strip_number;
  // int threshold;

  // std::string line;

  // while(std::getline(aida_strips_file, line))
  // {
  //     if(line[0] == '#') continue;

  //     std::istringstream iss(line);
  //     iss >> dssd_number >> xy >> strip_number >> threshold;
      
  //     // If the threshold is -1, add the strip to the list of strips to skip
  //     if (threshold == -1) {
  //         if (xy == "X") {
  //             strip_number = strip_number - 1;
  //             broken_xstrips.push_back(strip_number);
  //         } else if (xy == "Y") {
  //             strip_number = strip_number - 1;
  //             broken_ystrips.push_back(strip_number);
  //         }
  //     }
  // }

  // // Print out broken strips

  // for(auto x : broken_xstrips) {
  //     std::cout << "Broken X strip: " << x << std::endl;
  // }
  // for(auto y : broken_ystrips) {
  //     std::cout << "Broken Y strip: " << y << std::endl;
  // }

  TFile *Nb82 = new TFile("/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/82Nb.root");
  if (!Nb82) {
      std::cerr << "Error: Could not open gate file " << std::endl;
      return;
  }

  TFile *Nb84 = new TFile("/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/84Nb.root");
  if (!Nb84) {
      std::cerr << "Error: Could not open gate file " << std::endl;
      return;
  }

  TFile *Mo84 = new TFile("/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/84Mo.root");
  if (!Mo84) {
      std::cerr << "Error: Could not open gate file " << std::endl;
      return;
  }

  TFile *Mo85 = new TFile("/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/85Mo.root");
  if (!Mo85) {
      std::cerr << "Error: Could not open gate file " << std::endl;
      return;
  }

  TCutG *nb82_zaoq_cut = (TCutG*)Nb82->Get("cut_Z_AoQ");
  TCutG *nb82_zz2_cut = (TCutG*)Nb82->Get("cut_Z_Z2");

  TCutG *nb84_zaoq_cut = (TCutG*)Nb84->Get("cut_Z_AoQ");
  TCutG *nb84_zz2_cut = (TCutG*)Nb84->Get("cut_Z_Z2");

  TCutG *mo84_zaoq_cut = (TCutG*)Mo84->Get("cut_Z_AoQ");
  TCutG *mo84_zz2_cut = (TCutG*)Mo84->Get("cut_Z_Z2");

  TCutG *mo85_zaoq_cut = (TCutG*)Mo85->Get("cut_Z_AoQ");
  TCutG *mo85_zz2_cut = (TCutG*)Mo85->Get("cut_Z_Z2");

  // Open the output file
  TFile* outputFile = new TFile(output, "RECREATE");
  if (!outputFile) {
      std::cerr << "Error: Could not create output file " << std::endl;
      return;
  }

  // Create a new tree with the same structure as the old tree
  TTree* implant_tree = new TTree("aida_implant_tree", "New AIDA Analysis Tree");
  TTree* gatedimplant_82nb_tree = new TTree("aida_gatedimplant_82nb_tree", "New AIDA Analysis Tree");
  TTree* gatedimplant_84nb_tree = new TTree("aida_gatedimplant_84nb_tree", "New AIDA Analysis Tree");
  TTree* gatedimplant_84mo_tree = new TTree("aida_gatedimplant_84mo_tree", "New AIDA Analysis Tree");
  TTree* gatedimplant_85mo_tree = new TTree("aida_gatedimplant_85mo_tree", "New AIDA Analysis Tree");
  TTree* decay_tree = new TTree("aida_decay_tree", "New AIDA Analysis Tree");

  // TTree* bplast_tree = new TTree("bplast_tree", "New AIDA Analysis Tree");

  TTree* germanium_tree = new TTree("germanium_tree", "New AIDA Analysis Tree");


  // Define the data structures
  struct implant_data
  {
      uint64_t time;
      int stopped;
      int dssd;
      double x;
      double y;
      double energy;
      double energy_x;
      double energy_y;
      int sp;
      // int bp;
  } aida_implant_data;
  
  struct decay_data
  {
      uint64_t time;
      double x;
      double y;
      double energy;
      double energy_x;
      double energy_y;
      int dssd;
      int sp;
      // int bp;
  } aida_decay_data;

  struct bplast_data
  {
      uint64_t time;
      short id;
      double slow_tot;
      int sp;
      // int bp;
  } bplast_data;
  

  struct germanium_data
  {
      uint64_t time;
      double energy;
      int sp;
      // int bp;
  } germanium_data;

  // Create the branches
  implant_tree->Branch("implant", &aida_implant_data, "time/l:stopped/I:dssd/I:x/D:y/D:e/D:ex/D:ey/D:sp/I");
  gatedimplant_82nb_tree->Branch("gatedimplant_82nb", &aida_implant_data, "time/l:stopped/I:dssd/I:x/D:y/D:e/D:ex/D:ey/D:sp/I");
  gatedimplant_84nb_tree->Branch("gatedimplant_84nb", &aida_implant_data, "time/l:stopped/I:dssd/I:x/D:y/D:e/D:ex/D:ey/D:sp/I");
  gatedimplant_84mo_tree->Branch("gatedimplant_84mo", &aida_implant_data, "time/l:stopped/I:dssd/I:x/D:y/D:e/D:ex/D:ey/D:sp/I");
  gatedimplant_85mo_tree->Branch("gatedimplant_85mo", &aida_implant_data, "time/l:stopped/I:dssd/I:x/D:y/D:e/D:ex/D:ey/D:sp/I");
  decay_tree->Branch("decay", &aida_decay_data, "time/l:x/D:y/D:e/D:ex/D:ey/D:dssd/I:sp/I");
  // bplast_tree->Branch("bplast", &bplast_data, "time/l:id/S:slow_tot/D:sp/I:bp/I");
  germanium_tree->Branch("germanium", &germanium_data, "time/l:energy/D:sp/I");

  // Make histogram for drawing germanium data

  TH1F* germanium_energy_hist = new TH1F("germanium_energy_hist", "Germanium Energy", 1.5e3, 0, 1.5e3);
  TH2F* aida_implant_xy = new TH2F("aida_implant_xy", "AIDA Implant XY", 384, 0, 384, 128, 0, 128);
  TH2F* aida_implant_energy_xy = new TH2F("aida_implant_energy_xy", "AIDA Implant Energy XY", 500, 0, 7000, 500, 0, 7000);
  TH2F* aida_implant_gated_nb82_xy = new TH2F("aida_implant_gated_nb82_xy", "AIDA Implant XY (Gated ^{82}Nb)", 384, 0, 384, 128, 0, 128);
  TH2F* aida_implant_gated_nb84_xy = new TH2F("aida_implant_gated_nb84_xy", "AIDA Implant XY (Gated ^{84}Nb)", 384, 0, 384, 128, 0, 128);
  TH2F* aida_implant_gated_mo84_xy = new TH2F("aida_implant_gated_mo84_xy", "AIDA Implant XY (Gated ^{84}Mo)", 384, 0, 384, 128, 0, 128);
  TH2F* aida_implant_gated_mo85_xy = new TH2F("aida_implant_gated_mo85_xy", "AIDA Implant XY (Gated ^{85}Mo)", 384, 0, 384, 128, 0, 128);
  TH2F* aida_decay_xy = new TH2F("aida_decay_xy", "AIDA Decay XY", 384, 0, 384, 128, 0, 128);
  TH2F* aida_decay_energy_xy = new TH2F("aida_decay_energy_xy", "AIDA Decay Energy XY", 500, 0, 5000, 500, 0, 5000);
  TH1F* aida_decay_energy = new TH1F("aida_decay_energy", "AIDA Decay Energy", 300, 0, 3000);
  TH1F* aida_wr_times = new TH1F("aida_wr_times", "AIDA WR Times", number_of_slices, 0, duration_in_seconds);
  TH1F* germanium_decay_energy = new TH1F("germanium_decay_energy", "Germanium Decay Energy", 1.5e3, 0, 1.5e3);

  /*TH2F* frs_z_z2_hist = new TH2F("frs_z_z2_hist", "FRS Z vs Z2", 1000, 58, 68, 1000, 58, 68);*/
  TH2F* frs_z_aoq_hist = new TH2F("frs_z_aoq_hist", "FRS Z vs AoQ", 1000, 1.8,2.5, 1000, 30, 50);
  /*TH2F* frs_aoq_x4_hist = new TH2F("frs_aoq_x4_hist", "FRS AoQ vs X4", 1000, 2.0, 2.5, 1000, -100, 100);*/
  /*TH2F* frs_dedeg_z_hist = new TH2F("frs_dedeg_z_hist", "FRS dEdeg vs Z", 1000, 40, 50, 1000, 58, 68);*/

  TH1F* aida_decay_germanium_dt = new TH1F("aida_decay_germanium_dt", "AIDA Decay-Germanium dt", 5000, -50000, 50000);

  const char spinner[] = {'-', '\\', '|', '/'};
  int totalEntries = reader.GetEntries(true);

  int64_t tom_start = 1714192595223307800;
  int64_t tom_end = 1714207628772594000;

  int64_t last_implant_time = 0;

  // Loop over all entries in the old tree
  while (reader.Next()) {
      // Read the data from the old tree

      // Get the spill flag
      // bool spill_flag = *spill;


    // sizes
    int germaniumhits = germanium_time.GetSize();
    int aidadecayhits = decay_time.GetSize();
    int aidaimphits = implant_time.GetSize();
    // int frshits = frs_time.GetSize();
    // int bplasthits = bplast_id.GetSize();

    int mult_x = decay_x.GetSize();
    int mult_y = decay_y.GetSize();

      
    int spflag = 0;
    int decflag = 0;
      // int bpflag = 0;
      // int bp1flag = 0;
      // int bp2flag = 0;

/*  cout << "# of FRS hits:  " << frshits << endl;
    cout << "# of implant hits:  " << aidaimphits << endl;
    cout << "# of decay hits:  " << aidadecayhits << endl;
    cout << "# of gamma hits:  " << germaniumhits << endl;
    if (bplasthits > 0) cout << "# of bplast hits:  " << bplasthits << endl; */

    if(*spill == true) spflag = 1;
    if(*spill == false) spflag = 2;
      
      
      
    // for (int j = 0; j < bplasthits; j++) {
    //     if(bplast_id[j] < 64) bp1flag = 1;
    //     if(bplast_id[j] > 63 && bplast_id[j] < 128) bp2flag = 1;
    // }
  
    // if (bp1flag == 1 && bp2flag == 0) bpflag = 1;
    // if (bp1flag == 0 && bp2flag == 1) bpflag = 2;
    // if (bp1flag == 1 && bp2flag == 1) bpflag = 3;
          
//  cout << bp1flag << "  " << bp2flag << "  " << bpflag << endl;
    // for(int j=0; j<bplasthits; j++){
    //     bplast_data.time = bplast_time[j];
    //     bplast_data.id = bplast_id[j];
    //     bplast_data.slow_tot = bplast_slow_tot[j];
    //     bplast_data.sp = spflag;
    //     bplast_data.bp = bpflag;
    //     bplast_tree->Fill();
    // }
    int stopped = 0;

    if(aidaimphits > 0){
      for (int j = 0; j < aidaimphits; j++) {

        if (implant_dssd[j] == constants::DSSD && implant_stopped[j] == true){
          stopped = 1;
        } 
        else {
          stopped = 2;
        }
              
        aida_implant_data.time = implant_time[j];
        aida_implant_data.stopped = stopped;
        aida_implant_data.dssd = implant_dssd[j];
        aida_implant_data.x = implant_x[j];
        aida_implant_data.y = implant_y[j];
        aida_implant_data.energy = implant_energy[j];
        aida_implant_data.energy_x = implant_energy_x[j];
        aida_implant_data.energy_y = implant_energy_y[j];
        aida_implant_data.sp = spflag;
        // aida_implant_data.bp = bpflag;

        if (implant_dssd[j]==constants::DSSD && constants::DRAW_HISTS) {aida_implant_xy->Fill(implant_x[j], implant_y[j]);}
        if (implant_dssd[j]==constants::DSSD && constants::DRAW_HISTS) {aida_implant_energy_xy->Fill(implant_energy_x[j], implant_energy_y[j]);}
        implant_tree->Fill();
              // if(implant_x[j] >270 && implant_x[j] < 370 && implant_y[j] > 20 && implant_y[j] < 90){
              //     aida_implant_xy->Fill(implant_x[j], implant_y[j]);
              //     for(int j =0; j<frshits; j++){
              //         frs_z_z2_hist->Fill(frs_z[j], frs_z2[j]);
              //         frs_z_aoq_hist->Fill(frs_z[j], frs_aoq[j]);
              //         frs_aoq_x4_hist->Fill(frs_aoq[j], frs_x4[j]);
              //         frs_dedeg_z_hist->Fill(frs_z[j],frs_dedeg[j]);
              //     }
              // }
      }
    }

    std::set<int> filled_gatedimplanttree_82nb{};
    std::set<int> filled_gatedimplanttree_84nb{};
    std::set<int> filled_gatedimplanttree_84mo{};
    std::set<int> filled_gatedimplanttree_85mo{};

// Implants in coincidence with FRS
    if (aidaimphits >0) {
      for (auto const& FrsMultiItem : FrsMultiItem) {
        // std::cout << "Reached the cut!" << std::endl;
        std::vector<float> const& AoQ = FrsMultiItem.Get_ID_AoQ_corr_s2s4_mhtdc();
        std::vector<float> const& Z = FrsMultiItem.Get_ID_z41_mhtdc();
        std::vector<float> const& Z2 = FrsMultiItem.Get_ID_z42_mhtdc();
        for(int i =0; i<AoQ.size(); i++){
          if (i >= AoQ.size() || i >= Z.size()) {
            // std::cerr << "Error: Index out of bounds!" << std::endl;
            continue; // Skip this iteration if index is out of bounds
          }
          double aoq = AoQ[i];
          double z = Z[i];
          double z2 = Z[i];
          if (constants::DRAW_HISTS){
            frs_z_aoq_hist->Fill(aoq, z);
          }

          if(nb82_zaoq_cut->IsInside(aoq, z) && nb82_zz2_cut->IsInside(z, z2)){
            // std::cout << "Passed the cut!" << std::endl;
            for (int j = 0; j < aidaimphits; j++) {
              if (implant_dssd[j] == constants::DSSD && /*implant_stopped[j] == true &&*/ filled_gatedimplanttree_82nb.count(j) == 0) {
                // std::cout << "Found an implant!" << std::endl;
                aida_implant_data.time = implant_time[j];
                aida_implant_data.x = implant_x[j];
                aida_implant_data.y = implant_y[j];
                aida_implant_data.energy = implant_energy[j];
                aida_implant_data.energy_x = implant_energy_x[j];
                aida_implant_data.energy_y = implant_energy_y[j];
                aida_implant_data.sp = spflag;
                // aida_implant_data.bp = bpflag;
                gatedimplant_82nb_tree->Fill();
                if (constants::DRAW_HISTS){
                  aida_implant_gated_nb82_xy->Fill(implant_x[j], implant_y[j]);
                }
                implanted_82Nb++;
              }
              filled_gatedimplanttree_82nb.insert(j);
            }
          }

          if(nb84_zaoq_cut->IsInside(aoq, z) && nb84_zz2_cut->IsInside(z, z2)){
            // std::cout << "Passed the cut!" << std::endl;
            for (int j = 0; j < aidaimphits; j++) {
              if (implant_dssd[j] == constants::DSSD && /*implant_stopped[j] == true &&*/ filled_gatedimplanttree_84nb.count(j) == 0) {
                // std::cout << "Found an implant!" << std::endl;
                aida_implant_data.time = implant_time[j];
                aida_implant_data.x = implant_x[j];
                aida_implant_data.y = implant_y[j];
                aida_implant_data.energy = implant_energy[j];
                aida_implant_data.energy_x = implant_energy_x[j];
                aida_implant_data.energy_y = implant_energy_y[j];
                aida_implant_data.sp = spflag;
                // aida_implant_data.bp = bpflag;
                gatedimplant_84nb_tree->Fill();
                if (constants::DRAW_HISTS){
                  aida_implant_gated_nb84_xy->Fill(implant_x[j], implant_y[j]);
                }
                implanted_84Nb++;
              }
              filled_gatedimplanttree_84nb.insert(j);
            }
          }

          if(mo84_zaoq_cut->IsInside(aoq, z) && mo84_zz2_cut->IsInside(z, z2)){
            // std::cout << "Passed the cut!" << std::endl;
            for (int j = 0; j < aidaimphits; j++) {
              if (implant_dssd[j] == constants::DSSD && /*implant_stopped[j] == true &&*/ filled_gatedimplanttree_84mo.count(j) == 0) {
                // std::cout << "Found an implant!" << std::endl;
                aida_implant_data.time = implant_time[j];
                aida_implant_data.x = implant_x[j];
                aida_implant_data.y = implant_y[j];
                aida_implant_data.energy = implant_energy[j];
                aida_implant_data.energy_x = implant_energy_x[j];
                aida_implant_data.energy_y = implant_energy_y[j];
                aida_implant_data.sp = spflag;
                // aida_implant_data.bp = bpflag;
                gatedimplant_84mo_tree->Fill();
                if (constants::DRAW_HISTS){
                  aida_implant_gated_mo84_xy->Fill(implant_x[j], implant_y[j]);
                }
                implanted_84Mo++;
              }
              filled_gatedimplanttree_84mo.insert(j);
            }
          }

          if(mo85_zaoq_cut->IsInside(aoq, z) && mo85_zz2_cut->IsInside(z, z2)){
            // std::cout << "Passed the cut!" << std::endl;
            for (int j = 0; j < aidaimphits; j++) {
              if (implant_dssd[j] == constants::DSSD && /*implant_stopped[j] == true &&*/ filled_gatedimplanttree_85mo.count(j) == 0) {
                // std::cout << "Found an implant!" << std::endl;
                aida_implant_data.time = implant_time[j];
                aida_implant_data.x = implant_x[j];
                aida_implant_data.y = implant_y[j];
                aida_implant_data.energy = implant_energy[j];
                aida_implant_data.energy_x = implant_energy_x[j];
                aida_implant_data.energy_y = implant_energy_y[j];
                aida_implant_data.sp = spflag;
                // aida_implant_data.bp = bpflag;
                gatedimplant_85mo_tree->Fill();
                if (constants::DRAW_HISTS){
                  aida_implant_gated_mo85_xy->Fill(implant_x[j], implant_y[j]);
                }
                implanted_85Mo++;
              }
              filled_gatedimplanttree_85mo.insert(j);
            }
          }
                
        }
      }
    }


    std::set<int> filled_germtree{};
    std::set<int> filled_aidadecaytree{};

    // Decays only
    if(aidadecayhits == 1){
      for (int i = 0; i < aidadecayhits; i++) {
        if (TMath::Abs(decay_time_x[i] - decay_time_y[i]) < 5e3 && TMath::Abs(decay_energy_x[i] - decay_energy_y[i]) < 168 && decay_energy[i] > 151 && decay_energy[i] < 3000 && TMath::Abs(decay_time[aidadecayhits -1] - decay_time[0]) < 33000) {
          if (filled_aidadecaytree.count(i) == 0) {
            aida_decay_data.time = decay_time[i];
            aida_decay_data.x = decay_x[i];
            aida_decay_data.y = decay_y[i];
            aida_decay_data.energy = decay_energy[i];
            aida_decay_data.energy_x = decay_energy_x[i];
            aida_decay_data.energy_y = decay_energy_y[i];
            aida_decay_data.dssd = decay_dssd[i];
            aida_decay_data.sp = spflag;
            // aida_decay_data.bp = bpflag;
            decay_tree->Fill();
            if (decay_dssd[i] == constants::DSSD && constants::DRAW_HISTS){
              aida_decay_xy->Fill(decay_x[i], decay_y[i]);
              aida_decay_energy->Fill(decay_energy[i]);
              aida_decay_energy_xy->Fill(decay_energy_x[i], decay_energy_y[i]);
            }
          }
          filled_aidadecaytree.insert(i);
        }
      }
    }

    // Germaniums only
    for (int j = 0; j < germaniumhits; j++){
      if ( germanium_det[j] <= 15 && germanium_cry[j] <= 2 && germanium_energy[j] > 0 ){
        if ( filled_germtree.count(j) == 0 ){
          germanium_data.time = germanium_abs_ev_time[j];
          double energy = germanium_energy[j];
          germanium_data.energy = energy;
          germanium_data.sp = spflag;
          germanium_tree->Fill();
          if (constants::DRAW_HISTS){
            germanium_energy_hist->Fill(germanium_energy[j]);
          }
        }
      filled_germtree.insert(j);
      }
    }
    
    if (constants::DRAW_HISTS){
      // Decays in coincidence with Germanium
      if(aidadecayhits > 0 && germaniumhits > 0){
        for (int i = 0; i < aidadecayhits; i++){
          if (decay_dssd[i] == constants::DSSD){
            for (int j = 0; j < germaniumhits; j++){
              if (germanium_det[j] <= 15 && germanium_cry[j] <= 2 && germanium_energy[j] > 0){
                // if ((decay_time[i] - germanium_abs_ev_time[j]) > 14458 && (decay_time[i] - germanium_abs_ev_time[j]) < 16458) 
                aida_decay_germanium_dt->Fill(decay_time[i]-germanium_abs_ev_time[j]);
              }
            }
          }                
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


  // Write the new tree to the file
  if (constants::DRAW_HISTS){
    aida_implant_xy->Write(); 
    aida_implant_energy_xy->Write(); 
    aida_implant_gated_nb82_xy->Write(); 
    aida_implant_gated_nb84_xy->Write();
    aida_implant_gated_mo84_xy->Write();
    aida_implant_gated_mo85_xy->Write(); 
    aida_decay_xy->Write();
    aida_decay_energy_xy->Write();
    aida_decay_energy->Write();
    frs_z_aoq_hist->Write();
    germanium_energy_hist->Write(); 
    aida_decay_germanium_dt->Write();
  }

  std::cout << std::endl;
  std::cout << "Number of 82Nb implants:  " << implanted_82Nb << std::endl;
  std::cout << "Number of 84Nb implants:  " << implanted_84Nb << std::endl;
  std::cout << "Number of 84Mo implants:  " << implanted_84Mo << std::endl;
  std::cout << "Number of 85Mo implants:  " << implanted_85Mo << std::endl;
  
  if (constants::MAKE_ANATREES){
    implant_tree->Write();
    gatedimplant_82nb_tree->Write();
    gatedimplant_84nb_tree->Write();
    gatedimplant_84mo_tree->Write();
    gatedimplant_85mo_tree->Write();
    decay_tree->Write();
    germanium_tree->Write();
  }

  // Close the file
  file->Close();
  outputFile->Close();
  Nb84->Close();
  Nb82->Close();
  Mo84->Close();
  Mo85->Close();
  delete file;
  delete outputFile;
  delete Nb82;
  delete Nb84;
  delete Mo84;
  delete Mo85;
}
