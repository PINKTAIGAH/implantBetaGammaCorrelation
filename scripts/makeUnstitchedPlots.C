#include<map>
#include<set>
#include<iostream>
#include<unordered_map>
#include<cstdint>
#include<tuple>
#include<utility>
#include<string>
#include<vector>
#include<thread>
#include<mutex>
#include<algorithm>

#include<TF1.h>
#include<TH1F.h>
#include<TH2F.h>
#include<TFile.h>
#include<TTree.h>
#include<TCutG.h>
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

  const std::map <std::string, std::string> IMPLANT_GATES_INFILE_MAP = {
    {"82nb", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/new_gates/82nb_new.root"},
    {"zrnbmo", "/lustre/gamma/gbrunic/G302/analysis/implantBetaGammaCorrelation/gates/new_gates/zrnbmo_new.root"}
  };

  const std::string GATED_IMPLANT_NAME = "zrnbmo";

}


// *************************************************************************************
// ****************************** START MACRO ******************************************
// *************************************************************************************

void makeUnstitchedPlots(const char* input, const char* output){

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
  TTree* frsTree = (TTree*)file->Get("frs_tree");
  TTree* implantTree = (TTree*)file->Get("aida_implant_tree");
  TTree* gatedImplantTree = (TTree*)file->Get((std::string("aida_gatedimplant_") + constants::GATED_IMPLANT_NAME + std::string("_tree")).c_str());
  
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
  // ****************************** CREATE GATE MULTIMAPS ********************************
  // *************************************************************************************

  // Create a multimap with TCUTG objects of each gate and a implant event counter for each gate defined in the script
  std::multimap< std::string, std::tuple<TCutG*, TCutG*, TCutG*, TCutG*, TCutG*, TCutG*> > gatedimplant_cuts_map = {};

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
    TCutG* x2aoq_cut = (TCutG*)gate_file->Get("cut_x2_AoQ");
    TCutG* x4aoq_cut = (TCutG*)gate_file->Get("cut_x4_AoQ");
    TCutG* zz2_cut = (TCutG*)gate_file->Get("cut_Z_Z2");
    TCutG* zaoq_cut = (TCutG*)gate_file->Get("cut_Z_AoQ");
    TCutG* sci41Esci42E_cut = (TCutG*)gate_file->Get("cut_sci41E_sci42E");
    TCutG* zsci42E_cut = (TCutG*)gate_file->Get("cut_Z_sci42E");
    
    // Add gate information to multimap
    gatedimplant_cuts_map.emplace( gate_name, std::make_tuple(x2aoq_cut, x4aoq_cut, zz2_cut, zaoq_cut, sci41Esci42E_cut, zsci42E_cut) );   
  }

  // Make maps for gated histograms
  std::map<std::string, TH2F*> h2_x2_aoq_zz2cut_map = {};
  std::map<std::string, TH2F*> h2_x4_aoq_zz2cut_map = {};
  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_map = {};

  std::map<std::string, TH2F*> h2_x2_aoq_zz2cut_zsci42Ecut_map = {};
  std::map<std::string, TH2F*> h2_x4_aoq_zz2cut_zsci42Ecut_map = {};
  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_zsci42Ecut_map = {};

  std::map<std::string, TH2F*> h2_x2_aoq_zz2cut_sci41Esci42Ecut_map = {};
  std::map<std::string, TH2F*> h2_x4_aoq_zz2cut_sci41Esci42Ecut_map = {};
  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_sci41Esci42Ecut_map = {};

  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_x2aoqcut_map = {};
  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_x4aoqcut_map = {};
  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_zaoqcut_map = {};

  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_x2aoqcut_zsci42Ecut_map = {};
  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_x4aoqcut_zsci42Ecut_map = {};
  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_zaoqcut_zsci42Ecut_map = {};

  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_x2aoqcut_sci41Esci42Ecut_map = {};
  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_x4aoqcut_sci41Esci42Ecut_map = {};
  std::map<std::string, TH2F*> h2_z_aoq_zz2cut_zaoqcut_sci41Esci42Ecut_map = {};
  

  // *************************************************************************************
  // ****************************** DEFINE TREEREADER ************************************
  // *************************************************************************************

  // Create tree readers objects
  TTreeReader frsReader(frsTree);
  TTreeReader implantReader(implantTree);
  TTreeReader gatedImplantReader(gatedImplantTree);

  TTreeReaderValue<ULong64_t> wrTime(frsReader, "frs.time");
  TTreeReaderValue<double> aoq(frsReader, "frs.aoq");
  TTreeReaderValue<double> x2(frsReader, "frs.x2");
  TTreeReaderValue<double> x4(frsReader, "frs.x4");
  TTreeReaderValue<double> z41(frsReader, "frs.z41");
  TTreeReaderValue<double> z42(frsReader, "frs.z42");
  TTreeReaderValue<double> sci41E(frsReader, "frs.sci41E");
  TTreeReaderValue<double> sci42E(frsReader, "frs.sci42E");

  TTreeReaderValue<ULong64_t> impTime(implantReader, "implant.time");
  TTreeReaderValue<Int_t> impDssd(implantReader, "implant.dssd");
  TTreeReaderValue<double> impX(implantReader, "implant.x");
  TTreeReaderValue<double> impY(implantReader, "implant.y");
  TTreeReaderValue<double> impE(implantReader, "implant.e");
  TTreeReaderValue<double> impEx(implantReader, "implant.ex");
  TTreeReaderValue<double> impEy(implantReader, "implant.ey");
  TTreeReaderValue<Int_t> impCsx(implantReader, "implant.csx");
  TTreeReaderValue<Int_t> impCsy(implantReader, "implant.csy");
  TTreeReaderValue<Int_t> impSpill(implantReader, "implant.sp"); // sp = 1 spill, sp = 2 no spill

  TTreeReaderValue<ULong64_t> gimpTime(gatedImplantReader, (std::string("gatedimplant_") + constants::GATED_IMPLANT_NAME + std::string(".time")).c_str());
  TTreeReaderValue<Int_t> gimpStopped(gatedImplantReader, (std::string("gatedimplant_") + constants::GATED_IMPLANT_NAME + std::string(".stopped")).c_str());
  TTreeReaderValue<Int_t> gimpDssd(gatedImplantReader, (std::string("gatedimplant_") + constants::GATED_IMPLANT_NAME + std::string(".dssd")).c_str());
  TTreeReaderValue<double> gimpX(gatedImplantReader, ( std::string("gatedimplant_")+constants::GATED_IMPLANT_NAME+std::string(".x") ).c_str());
  TTreeReaderValue<double> gimpY(gatedImplantReader, ( std::string("gatedimplant_")+constants::GATED_IMPLANT_NAME+std::string(".y") ).c_str());
  TTreeReaderValue<double> gimpE(gatedImplantReader, ( std::string("gatedimplant_")+constants::GATED_IMPLANT_NAME+std::string(".e") ).c_str());
  TTreeReaderValue<double> gimpEx(gatedImplantReader, ( std::string("gatedimplant_")+constants::GATED_IMPLANT_NAME+std::string(".ex") ).c_str());
  TTreeReaderValue<double> gimpEy(gatedImplantReader, ( std::string("gatedimplant_")+constants::GATED_IMPLANT_NAME+std::string(".ey") ).c_str());
  TTreeReaderValue<Int_t> gimpCsx(gatedImplantReader, ( std::string("gatedimplant_")+constants::GATED_IMPLANT_NAME+std::string(".csx") ).c_str());
  TTreeReaderValue<Int_t> gimpCsy(gatedImplantReader, ( std::string("gatedimplant_")+constants::GATED_IMPLANT_NAME+std::string(".csy") ).c_str());
  TTreeReaderValue<Int_t> gimpSpill(gatedImplantReader, ( std::string("gatedimplant_")+constants::GATED_IMPLANT_NAME+std::string(".sp") ).c_str());


  // *************************************************************************************
  // ****************************** DEFINE HISTOGRAMS ************************************
  // *************************************************************************************

  // 1D Histograms
  TH1F* h1_aoq = new TH1F("aoq", "AoQ", 1000, 1.8, 2.5);
  TH1F* h1_z41 = new TH1F("z41", "Z41", 1000, 30, 50);
  TH1F* h1_z42 = new TH1F("z42", "Z42", 1000, 30, 50);
  TH1F* h1_x2 = new TH1F("x2", "X2", 1000, -150, 150);
  TH1F* h1_x4 = new TH1F("x4", "X4", 1000, -150, 150);
  TH1F* h1_sci41E = new TH1F("sci41E", "Sci41 E", 1000, 0, 2000);
  TH1F* h1_sci42E = new TH1F("sci42E", "Sci42 E", 1000, 0, 2000);

  TH1F* h1_implantDt = new TH1F("implantDt", "Implant X-Y Strip dT", 1000, -50e3, 50e3);
  
  // 2D Histograms
  TH2F* h2_z42_aoq = new TH2F("z42_aoq", "Z42 vs AoQ", 1000, 1.8, 2.5, 1000, 30, 50);
  TH2F* h2_z41_z42 = new TH2F("z41_z42", "Z41 vs Z42", 1000, 30, 50, 1000, 30, 50);
  TH2F* h2_x2_aoq = new TH2F("x2_aoq", "X2 vs AoQ", 1000, 1.8, 2.5, 1000, -150, 150);
  TH2F* h2_x4_aoq = new TH2F("x4_aoq", "X4 vs AoQ", 1000, 1.8, 2.5, 1000, -150, 150);
  TH2F* h2_sci41E_sci42E = new TH2F("sci42E_sci42E", "Sci41 E vs Sci42 E", 1000, 0, 2000, 1000, 0, 3000);
  TH2F* h2_z42_sci42E = new TH2F("z42_sci42E", "Z42 vs Sci42 E", 1000, 0, 3000, 1000, 30, 50);
  TH2F* h2_z42_sci41E = new TH2F("z42_sci41E", "Z42 vs Sci41 E", 1000, 0, 2000, 1000, 30, 50);

  TH2F* h2_implantPosition_allImplant = new TH2F("implant_posiiton", "Implant Position", 384, 0, 384, 128, 0, 128);
  TH2F* h2_implantPosition_gatedImplant = new TH2F( (constants::GATED_IMPLANT_NAME+"_implant_position").c_str() , (constants::GATED_IMPLANT_NAME+" Implant Position").c_str(), 384, 0, 384, 128, 0, 128);

  // Loop over all gated implants and fill map
  for (auto itr : constants::IMPLANT_GATES_INFILE_MAP ){

    std::string gatedimplant_name = itr.first;
    std::string histogram_name;
    std::string histogram_title;

    // X2 vs AoQ (w/ Z Z2 cut)
    histogram_name = gatedimplant_name + "_x2_aoq_zz2cut";
    histogram_title = gatedimplant_name + " X2 vs AoQ (w/ Z Z2 cut))";
    TH2F* h2_x2_aoq_zz2cut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, -150, 150);
    h2_x2_aoq_zz2cut_map.emplace(gatedimplant_name, h2_x2_aoq_zz2cut);

    // X4 vs AoQ (w/ Z Z2 cut)
    histogram_name = gatedimplant_name + "_x4_aoq_zz2cut";
    histogram_title = gatedimplant_name + " X4 vs AoQ (w/ Z Z2 cut))";
    TH2F* h2_x4_aoq_zz2cut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, -150, 150);
    h2_x4_aoq_zz2cut_map.emplace(gatedimplant_name, h2_x4_aoq_zz2cut);

    // Z41 vs AoQ (w/ Z Z2 cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut))";
    TH2F* h2_z_aoq_zz2cut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut);

    // X2 vs AoQ (w/ Z Z2 cut and Z Sci42E cut)
    histogram_name = gatedimplant_name + "_x2_aoq_zz2cut_zsci42Ecut";
    histogram_title = gatedimplant_name + " X2 vs AoQ (w/ Z Z2 cut and Z Sci42E cut))";
    TH2F* h2_x2_aoq_zz2cut_zsci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, -150, 150);
    h2_x2_aoq_zz2cut_zsci42Ecut_map.emplace(gatedimplant_name, h2_x2_aoq_zz2cut_zsci42Ecut);

    // X4 vs AoQ (w/ Z Z2 cut and Z Sci42E cut)
    histogram_name = gatedimplant_name + "_x4_aoq_zz2cut_zsci42Ecut";
    histogram_title = gatedimplant_name + " X4 vs AoQ (w/ Z Z2 cut and Z Sci42E cut))";
    TH2F* h2_x4_aoq_zz2cut_zsci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, -150, 150);
    h2_x4_aoq_zz2cut_zsci42Ecut_map.emplace(gatedimplant_name, h2_x4_aoq_zz2cut_zsci42Ecut);

    // Z41 vs AoQ (w/ Z Z2 cut and Z Sci42E cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_zsci42Ecut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and Z Sci42E cut))";
    TH2F* h2_z_aoq_zz2cut_zsci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_zsci42Ecut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_zsci42Ecut);

    // X2 vs AoQ (w/ Z Z2 cut and Sci41E Sci42E cut)
    histogram_name = gatedimplant_name + "_x2_aoq_zz2cut_sci41Esci42Ecut";
    histogram_title = gatedimplant_name + " X2 vs AoQ (w/ Z Z2 cut and Sci41E Sci42E cut))";
    TH2F* h2_x2_aoq_zz2cut_sci41Esci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, -150, 150);
    h2_x2_aoq_zz2cut_sci41Esci42Ecut_map.emplace(gatedimplant_name, h2_x2_aoq_zz2cut_sci41Esci42Ecut);

    // X4 vs AoQ (w/ Z Z2 cut and Sci41E Sci42E cut)
    histogram_name = gatedimplant_name + "_x4_aoq_zz2cut_sci41Esci42Ecut";
    histogram_title = gatedimplant_name + " X4 vs AoQ (w/ Z Z2 cut and Sci41E Sci42E cut))";
    TH2F* h2_x4_aoq_zz2cut_sci41Esci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, -150, 150);
    h2_x4_aoq_zz2cut_sci41Esci42Ecut_map.emplace(gatedimplant_name, h2_x4_aoq_zz2cut_sci41Esci42Ecut);

    // Z41 vs AoQ (w/ Z Z2 cut and Sci41E Sci42E cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_sci41Esci42Ecut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and Sci41E Sci42E cut))";
    TH2F* h2_z_aoq_zz2cut_sci41Esci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_sci41Esci42Ecut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_sci41Esci42Ecut);

    // Z41 vs AoQ (w/ Z Z2 cut and X2 AoQ cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_x2aoqcut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and X2 AoQ cut))";
    TH2F* h2_z_aoq_zz2cut_x2aoqcut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_x2aoqcut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_x2aoqcut);

    // Z41 vs AoQ (w/ Z Z2 cut and X4 AoQ cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_x4aoqcut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and X4 AoQ cut))";
    TH2F* h2_z_aoq_zz2cut_x4aoqcut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_x4aoqcut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_x4aoqcut);

    // Z41 vs AoQ (w/ Z Z2 cut and Z41 AoQ cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_zaoqcut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and Z41 AoQ cut))";
    TH2F* h2_z_aoq_zz2cut_zaoqcut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_zaoqcut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_zaoqcut);

    // Z41 vs AoQ (w/ Z Z2 cut and X2 AoQ cut and Z Sci42E cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_x2aoqcut_zsci42Ecut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and X2 AoQ cut and Z Sci42E cut))";
    TH2F* h2_z_aoq_zz2cut_x2aoqcut_zsci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_x2aoqcut_zsci42Ecut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_x2aoqcut_zsci42Ecut);

    // Z41 vs AoQ (w/ Z Z2 cut and X4 AoQ cut and Z Sci42E cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_x4aoqcut_zsci42Ecut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and X4 AoQ cut and Z Sci42E cut))";
    TH2F* h2_z_aoq_zz2cut_x4aoqcut_zsci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_x4aoqcut_zsci42Ecut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_x4aoqcut_zsci42Ecut);

    // Z41 vs AoQ (w/ Z Z2 cut and Z41 AoQ cut and Z Sci42E cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_zaoqcut_zsci42Ecut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and Z41 AoQ cut and Z Sci42E cut))";
    TH2F* h2_z_aoq_zz2cut_zaoqcut_zsci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_zaoqcut_zsci42Ecut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_zaoqcut_zsci42Ecut);

    // Z41 vs AoQ (w/ Z Z2 cut and X2 AoQ cut and Sci41E Sci42E cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_x2aoqcut_sci41Esci42Ecut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and X2 AoQ cut and Sci41E Sci42E cut))";
    TH2F* h2_z_aoq_zz2cut_x2aoqcut_sci41Esci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_x2aoqcut_sci41Esci42Ecut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_x2aoqcut_sci41Esci42Ecut);

    // Z41 vs AoQ (w/ Z Z2 cut and X4 AoQ cut and Sci41E Sci42E cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_x4aoqcut_sci41Esci42Ecut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and X4 AoQ cut and Sci41E Sci42E cut))";
    TH2F* h2_z_aoq_zz2cut_x4aoqcut_sci41Esci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_x4aoqcut_sci41Esci42Ecut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_x4aoqcut_sci41Esci42Ecut);

    // Z41 vs AoQ (w/ Z Z2 cut and Z41 AoQ cut and Sci41E Sci42E cut)
    histogram_name = gatedimplant_name + "_z_aoq_zz2cut_zaoqcut_sci41Esci42Ecut";
    histogram_title = gatedimplant_name + " Z41 vs AoQ (w/ Z Z2 cut and Z41 AoQ cut and Sci41E Sci42E cut))";
    TH2F* h2_z_aoq_zz2cut_zaoqcut_sci41Esci42Ecut = new TH2F(histogram_name.c_str(), histogram_title.c_str(), 1000, 1.8, 2.5, 1000, 30, 50);
    h2_z_aoq_zz2cut_zaoqcut_sci41Esci42Ecut_map.emplace(gatedimplant_name, h2_z_aoq_zz2cut_zaoqcut_sci41Esci42Ecut);

  }

  // *************************************************************************************
  // ****************************** LOOP OVER FRS ENTRIES ********************************
  // *************************************************************************************

  while ( frsReader.Next()){

    // Fill 1D histograms
    h1_aoq->Fill(*aoq);
    h1_z41->Fill(*z41);
    h1_z42->Fill(*z42);
    h1_x2->Fill(*x2);
    h1_x4->Fill(*x4);
    h1_sci41E->Fill(*sci41E);
    h1_sci42E->Fill(*sci42E);

    // Fill 2D histograms
    h2_z42_aoq->Fill(*aoq, *z42);
    h2_z41_z42->Fill(*z41, *z42);
    h2_x2_aoq->Fill(*aoq, *x2);
    h2_x4_aoq->Fill(*aoq, *x4);
    h2_sci41E_sci42E->Fill(*sci41E, *sci42E);
    h2_z42_sci42E->Fill(*sci42E, *z42);
    h2_z42_sci41E->Fill(*sci41E, *z42);

    for (const auto& itr : gatedimplant_cuts_map){
      
      // Extract all TCutG objects
      const auto [x2aoq_cut, x4aoq_cut, zz2_cut, zaoq_cut, sci41Esci42E_cut, zsci42E_cut] = itr.second;
      std::string implant_key = itr.first;

      // Only Z Z2 cut
      if ( zz2_cut->IsInside(*z41, *z42) ){
        h2_x2_aoq_zz2cut_map[implant_key]->Fill(*aoq, *x2);
        h2_x4_aoq_zz2cut_map[implant_key]->Fill(*aoq, *x4);
        h2_z_aoq_zz2cut_map[implant_key]->Fill(*aoq, *z41);
      }

      // Z Z2 cut && Z SciE42
      if ( zsci42E_cut->IsInside(*sci42E, *z41) && zz2_cut->IsInside(*z41, *z42) ){
        h2_x2_aoq_zz2cut_zsci42Ecut_map[implant_key]->Fill(*aoq, *x2);
        h2_x4_aoq_zz2cut_zsci42Ecut_map[implant_key]->Fill(*aoq, *x4);
        h2_z_aoq_zz2cut_zsci42Ecut_map[implant_key]->Fill(*aoq, *z41);
      }

      // Z Z2 cut && Sci41E SciE42
      if ( sci41Esci42E_cut->IsInside(*sci41E, *sci42E) && zz2_cut->IsInside(*z41, *z42) ){
        h2_x2_aoq_zz2cut_sci41Esci42Ecut_map[implant_key]->Fill(*aoq, *x2);
        h2_x4_aoq_zz2cut_sci41Esci42Ecut_map[implant_key]->Fill(*aoq, *x4);
        h2_z_aoq_zz2cut_sci41Esci42Ecut_map[implant_key]->Fill(*aoq, *z41);
      }

      // Z Z2 cut && X2 AoQ cut
      if ( x2aoq_cut->IsInside(*aoq, *x2) && zz2_cut->IsInside(*z41, *z42) ){
        // W/O x4 reactions gate
        h2_z_aoq_zz2cut_x2aoqcut_map[implant_key]->Fill(*aoq, *z41);
        // W/ Z Sci42 cut 
        if ( zsci42E_cut->IsInside(*sci42E, *z41) ) h2_z_aoq_zz2cut_x2aoqcut_zsci42Ecut_map[implant_key]->Fill(*aoq, *z41);
        // W/ Sci41E Sci42 cut 
        if ( sci41Esci42E_cut->IsInside(*sci41E, *sci42E) ) h2_z_aoq_zz2cut_x2aoqcut_sci41Esci42Ecut_map[implant_key]->Fill(*aoq, *z41);
      }

      // Z Z2 cut && X4 AoQ cut
      if ( x4aoq_cut->IsInside(*aoq, *x4) && zz2_cut->IsInside(*z41, *z42) ){
        // W/O x4 reactions gate
        h2_z_aoq_zz2cut_x4aoqcut_map[implant_key]->Fill(*aoq, *z41);
        // W/ Z Sci42 cut 
        if ( zsci42E_cut->IsInside(*sci42E, *z41) ) h2_z_aoq_zz2cut_x4aoqcut_zsci42Ecut_map[implant_key]->Fill(*aoq, *z41);
        // W/ Sci41E Sci42 cut 
        if ( sci41Esci42E_cut->IsInside(*sci41E, *sci42E) ) h2_z_aoq_zz2cut_x4aoqcut_sci41Esci42Ecut_map[implant_key]->Fill(*aoq, *z41);
      }

      // Z Z2 cut && Z41 AoQ cut
      if ( zaoq_cut->IsInside(*aoq, *z41) && zz2_cut->IsInside(*z41, *z42) ){
        // W/O x4 reactions gate
        h2_z_aoq_zz2cut_zaoqcut_map[implant_key]->Fill(*aoq, *z41);
        // W/ Z Sci42 cut 
        if ( zsci42E_cut->IsInside(*sci42E, *z41) ) h2_z_aoq_zz2cut_zaoqcut_zsci42Ecut_map[implant_key]->Fill(*aoq, *z41);
        // W/ Sci41E Sci42 cut 
        if ( sci41Esci42E_cut->IsInside(*sci41E, *sci42E) ) h2_z_aoq_zz2cut_zaoqcut_sci41Esci42Ecut_map[implant_key]->Fill(*aoq, *z41);
      }

    }

  }

  std::cout << std::endl;


  // *************************************************************************************
  // ************************ LOOP OVER Gated Implant ENTRIES ****************************
  // *************************************************************************************

  while ( gatedImplantReader.Next() ){

    if (*gimpDssd != 1) continue;

    // Fill histograms
    h2_implantPosition_gatedImplant->Fill(*gimpX, *gimpY);

  }

  // *************************************************************************************
  // ****************************** LOOP OVER Implant ENTRIES ****************************
  // *************************************************************************************

  while ( implantReader.Next() ){

    if (*impDssd != 1) continue;

    // Fill histograms
    h2_implantPosition_allImplant->Fill(*impX, *impY);

  }

  // *************************************************************************************
  // ****************************** write trees ******************************************
  // *************************************************************************************

  std::cout << "Writing Histos ..." << std::endl;
  outputFile->cd();

  TDirectory* h1Nogate = outputFile->mkdir("h1_nogate");
  h1Nogate->cd();
  h1_aoq->Write();
  h1_z41->Write();
  h1_z42->Write();
  h1_x2->Write();
  h1_x4->Write();
  h1_sci41E->Write();
  h1_sci42E->Write();
  gFile->cd();

  TDirectory* h2Nogate = outputFile->mkdir("h2_nogate");
  h2Nogate->cd();
  h2_z42_aoq->Write();
  h2_z41_z42->Write();
  h2_x2_aoq->Write();
  h2_x4_aoq->Write();
  h2_sci41E_sci42E->Write();
  h2_z42_sci42E->Write();
  h2_z42_sci41E->Write();
  gFile->cd();
  
  TDirectory* h2ZZ2_gate = outputFile->mkdir("h2_ZZ2_gate");
  h2ZZ2_gate->cd();
  for (const auto& itr : h2_x2_aoq_zz2cut_map) itr.second->Write();
  for (const auto& itr : h2_x4_aoq_zz2cut_map) itr.second->Write();
  for (const auto& itr : h2_z_aoq_zz2cut_map) itr.second->Write();
  gFile->cd();

  TDirectory* h2Zz2Zsci42e_gate = outputFile->mkdir("h2_ZZ2_ZSci42e_gate");
  h2Zz2Zsci42e_gate->cd();
  for (const auto& itr : h2_x2_aoq_zz2cut_zsci42Ecut_map) itr.second->Write();
  for (const auto& itr : h2_x4_aoq_zz2cut_zsci42Ecut_map) itr.second->Write();
  for (const auto& itr : h2_z_aoq_zz2cut_zsci42Ecut_map) itr.second->Write();
  gFile->cd();

  TDirectory* h2Zz2Sci41Esci42e_gate = outputFile->mkdir("h2_ZZ2_Sci41ESci42e_gate");
  h2Zz2Sci41Esci42e_gate->cd();
  for (const auto& itr : h2_x2_aoq_zz2cut_sci41Esci42Ecut_map) itr.second->Write();
  for (const auto& itr : h2_x4_aoq_zz2cut_sci41Esci42Ecut_map) itr.second->Write();
  for (const auto& itr : h2_z_aoq_zz2cut_sci41Esci42Ecut_map) itr.second->Write();
  gFile->cd();

  TDirectory* h2Zz2X2AoQ_gate = outputFile->mkdir("h2_ZZ2_X2AoQ_gate");
  h2Zz2X2AoQ_gate->cd();
  for (const auto& itr : h2_z_aoq_zz2cut_x2aoqcut_map) itr.second->Write();
  for (const auto& itr : h2_z_aoq_zz2cut_x2aoqcut_zsci42Ecut_map) itr.second->Write();
  for (const auto& itr : h2_z_aoq_zz2cut_x2aoqcut_sci41Esci42Ecut_map) itr.second->Write();
  gFile->cd();

  TDirectory* h2Zz2X4AoQ_gate = outputFile->mkdir("h2_ZZ2_X4AoQ_gate");
  h2Zz2X4AoQ_gate->cd();
  for (const auto& itr : h2_z_aoq_zz2cut_x4aoqcut_map) itr.second->Write();
  for (const auto& itr : h2_z_aoq_zz2cut_x4aoqcut_zsci42Ecut_map) itr.second->Write();
  for (const auto& itr : h2_z_aoq_zz2cut_x4aoqcut_sci41Esci42Ecut_map) itr.second->Write();
  gFile->cd();
  
  TDirectory* h2Zz2ZZAoQ_gate = outputFile->mkdir("h2_ZZ2_ZAoQ_gate");
  h2Zz2ZZAoQ_gate->cd();
  for (const auto& itr : h2_z_aoq_zz2cut_zaoqcut_map) itr.second->Write();
  for (const auto& itr : h2_z_aoq_zz2cut_zaoqcut_zsci42Ecut_map) itr.second->Write();
  for (const auto& itr : h2_z_aoq_zz2cut_zaoqcut_sci41Esci42Ecut_map) itr.second->Write();
  gFile->cd();

  TDirectory* gatedImplantPosition = outputFile->mkdir("gated_implant_position");
  gatedImplantPosition->cd();
  h2_implantPosition_allImplant->Write();
  h2_implantPosition_gatedImplant->Write();
  gFile->cd();


  std::cout << "Finished writing histos successfully!" << std::endl << std::endl;

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