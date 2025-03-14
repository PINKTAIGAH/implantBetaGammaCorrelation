#include <iostream>
#include <map>
#include <unordered_map>
#include <tuple>
#include <utility>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

// *************** DEFINE SCRIPT CONSTANTS ***************/
namespace constants {
    const bool ONLY_OFFSPILL_DECAY = true; // Check for onspill decay matches

    const int64_t TIME_SCALE = 1e9; // Timescale of time wrt ns
    const int64_t TIME_THRESHOLD = 100 * TIME_SCALE; // 1000 ms
    const int64_t POSITION_THRESHOLD = 1; // 1 strips in X and Y around the centre

    const double HALF_LIFE = 9.8;
    const double HALF_LIFE_NS = HALF_LIFE * TIME_SCALE;
}


// Tree reader and analyser for implant decay gamma analysis

bool is_overlapping(double imp_x, double imp_y, double imp_dx, double imp_dy,
                    double beta_x, double beta_y, double beta_dx, double beta_dy, double B) {
    return (imp_y + (imp_dy / 2.0) >= (beta_y - ((beta_dy / 2.0) + B)) &&
            imp_y - (imp_dy / 2.0) <= (beta_y + ((beta_dy / 2.0) + B)) &&
            imp_x + (imp_dx / 2.0) >= (beta_x - ((beta_dx / 2.0) + B)) &&
            imp_x - (imp_dx / 2.0) <= (beta_x + ((beta_dx / 2.0) + B)));
}

// Custom hash function for std::pair<int, int>
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& pair) const {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};

void ionbeta(const char* input, const char* output){

    // Load file
    TFile* file = TFile::Open(input);
    if (!file) {
        std::cerr << "Error: Could not open file " << std::endl;
        return;
    }
    std::cout << "File loaded: "<< file->GetName() << std::endl;

    // Get the trees
    TTree* implant_tree = (TTree*)file->Get("aida_implant_tree");
    TTree* gatedimplant_tree = (TTree*)file->Get("aida_gatedimplant_84mo_tree");
    TTree* decay_tree = (TTree*)file->Get("aida_decay_tree");
    TTree* germanium_tree = (TTree*)file->Get("germanium_tree");

    // Open the output file
    TFile* outputFile = new TFile(output, "RECREATE");
    if (!outputFile) {
        std::cerr << "Error: Could not create output file " << std::endl;
        return;
    }

    // Set tree readers

    TTreeReader implant_reader(implant_tree);
    TTreeReader gatedimplant_reader(gatedimplant_tree);
    TTreeReader decay_reader(decay_tree);
    TTreeReader germanium_reader(germanium_tree);

    TTreeReaderValue<ULong64_t> implant_time(implant_reader, "implant.time");
    TTreeReaderValue<double> implant_x(implant_reader, "implant.x");
    TTreeReaderValue<double> implant_y(implant_reader, "implant.y");
    TTreeReaderValue<int> implant_dssd(implant_reader, "implant.dssd");
    TTreeReaderValue<double> implant_e(implant_reader, "implant.e");
    TTreeReaderValue<double> implant_ex(implant_reader, "implant.ex");
    TTreeReaderValue<double> implant_ey(implant_reader, "implant.ey");
    TTreeReaderValue<Int_t> implant_spill(implant_reader, "implant.sp"); // sp = 1 spill, sp = 2 no spill
    // TTreeReaderValue<Int_t> implant_bplast(implant_reader, "implant.bp"); // bp = 0 neither fired, bp = 1 only bp1 fired, bp = 2 only bp2 fired, bp = 3 both fired

    TTreeReaderValue<ULong64_t> gatedimplant_time(gatedimplant_reader, "gatedimplant_84mo.time");
    TTreeReaderValue<double> gatedimplant_x(gatedimplant_reader, "gatedimplant_84mo.x");
    TTreeReaderValue<double> gatedimplant_y(gatedimplant_reader, "gatedimplant_84mo.y");
    TTreeReaderValue<int> gatedimplant_dssd(gatedimplant_reader, "gatedimplant_84mo.dssd");
    TTreeReaderValue<double> gatedimplant_e(gatedimplant_reader, "gatedimplant_84mo.e");
    TTreeReaderValue<double> gatedimplant_ex(gatedimplant_reader, "gatedimplant_84mo.ex");
    TTreeReaderValue<double> gatedimplant_ey(gatedimplant_reader, "gatedimplant_84mo.ey");
    TTreeReaderValue<Int_t> gatedimplant_spill(gatedimplant_reader, "gatedimplant_84mo.sp"); // sp = 1 spill, sp = 2 no spill
    // TTreeReaderValue<Int_t> gatedimplant_bplast(gatedimplant_reader, "gatedimplant.bp"); // bp = 0 neither fired, bp = 1 only bp1 fired, bp = 2 only bp2 fired, bp = 3 both fired

    TTreeReaderValue<ULong64_t> decay_time(decay_reader, "decay.time");
    TTreeReaderValue<double> decay_x(decay_reader, "decay.x");
    TTreeReaderValue<double> decay_y(decay_reader, "decay.y");
    TTreeReaderValue<double> decay_e(decay_reader, "decay.e");
    TTreeReaderValue<double> decay_ex(decay_reader, "decay.ex");
    TTreeReaderValue<double> decay_ey(decay_reader, "decay.ey");
    TTreeReaderValue<Int_t> decay_spill(decay_reader, "decay.sp"); // sp = 1 spill, sp = 2 no spill
    // TTreeReaderValue<Int_t> decay_bplast(decay_reader, "decay.bp"); // bp = 0 neither fired, bp = 1 only bp1 fired, bp = 2 only bp2 fired, bp = 3 both fired

    TTreeReaderValue<ULong64_t> germanium_time(germanium_reader, "germanium.time");
    TTreeReaderValue<Double_t> germanium_energy(germanium_reader, "germanium.energy");
    TTreeReaderValue<Int_t> germanium_spill(germanium_reader, "germanium.sp"); // sp = 1 spill, sp = 2 no spill
    // TTreeReaderValue<Int_t> germanium_bplast(germanium_reader, "germanium.bp"); // bp = 0 neither fired, bp = 1 only bp1 fired, bp = 2 only bp2 fired, bp = 3 both fired
    
    // Experiment information
    uint64_t wr_experiment_start = 1.7401830e+18;
    uint64_t wr_experiment_end = 1.74022529e+18;
    int64_t duration_in_seconds = (wr_experiment_end - wr_experiment_start)/1e9;
    int64_t slices_every = 60; //s
    int64_t number_of_slices = duration_in_seconds/slices_every;

    // half life of 82Nb
    double total_time = constants::HALF_LIFE; // seconds
    int bins_per_half_life = 100;
    int number_of_bins = constants::HALF_LIFE*bins_per_half_life;


    // Histograms
    TH2F* aida_implant_xy = new TH2F("aida_implant_xy", "AIDA Implant XY", 384, 0, 384, 128, 0, 128);
    TH2F* aida_decay_xy = new TH2F("aida_decay_xy", "AIDA Decay XY", 384, 0, 384, 128, 0, 128);
    TH2F* aida_matched_xy = new TH2F("aida_matched_xy", "AIDA Matched XY", 384, 0, 384, 128, 0, 128);
    TH1F* aida_wr_times = new TH1F("aida_wr_times", "AIDA WR Times", number_of_slices, wr_experiment_start, wr_experiment_end);
    TH1F* aida_implant_veto_dt = new TH1F("aida_implant_veto_dt", "Implant-Decay #Deltat;Implant-Decay #Deltat (); Counts/", constants::TIME_THRESHOLD/constants::TIME_SCALE,-constants::TIME_THRESHOLD/constants::TIME_SCALE, constants::TIME_THRESHOLD/constants::TIME_SCALE);


    /*TH1F* germanium_energy_hist = new TH1F("germanium_energy_hist", "Germanium Energy", 1.5e3, 0, 1.5e3);*/
    /*TH1F* aida_implant_decay_time = new TH1F("aida_implant_decay_time", "AIDA Implant Decay Time", number_of_bins, -constants::HALF_LIFE*2, constants::HALF_LIFE*2);*/
    /*TH1F* germanium_decay_energy = new TH1F("germanium_decay_energy", "Germanium Decay Energy", 1.5e3, 0, 1.5e3);*/
    /*TH1F* germanium_aida_wr_dt = new TH1F("germanium_aida_wr_dt", "Germanium AIDA WR Time Difference", 1e3, -1e5, 1e5);*/
    /*TH1F* aida_implant_before_decay_counter = new TH1F("aida_implant_before_decay_counter", "AIDA Implant Before Decay Counter", 10, 0, 10);*/
    /*TH1F* aida_gatedimplant_before_decay_counter = new TH1F("aida_gatedimplant_before_decay_counter", "AIDA Gated Implant Before Decay Counter", 1e3, 0, 1e3);*/
    /*TH2F* aida_posdiff_time = new TH2F("aida_posdiff_time", "AIDA Position Difference vs Time", 1e3, 0, 1e3, 1e3, -1e4, 1e5);*/
    /*TH1F* aida_averageimplant_150s = new TH1F("aida_averageimplant_150s", "AIDA Average Implant 150s", number_of_slices, 0, duration_in_seconds);*/
    /*TH1F* aida_averagegatedimplant_150s = new TH1F("aida_averagegatedimplant_150s", "AIDA Average Gated Implant 150s", number_of_slices, 0, duration_in_seconds);*/
    /*TH2F* aida_implant_dt_vs_pos_x = new TH2F("aida_implant_dt_vs_pos_x", "AIDA Implant-Decay Time Difference vs Position Difference X", 1e4, 0, 3e9, 384, 0, 384);*/
    /*TH2F* aida_implant_dt_vs_pos_y = new TH2F("aida_implant_dt_vs_pos_y", "AIDA Implant-Decay Time Difference vs Position Difference Y", 1e4, 0, 3e9, 384, 0, 384);*/
    /*TH2F* aida_implant_dt_vs_pos_diff = new TH2F("aida_implant_dt_vs_pos_diff", "AIDA Implant-Decay Time Difference vs Position Difference X-Y", 1e4, 0, 3e9, 384, 0, 384);*/
    /**/

    // Make new maps for processing the data

    enum EventType { GATEDIMPLANT, IMPLANT, DECAY };

    std::map<int64_t, std::tuple<double, double, int, int>> implant_map;
    std::map<int64_t, std::tuple<double, double, int, int>> gatedimplant_map;
    std::map<int64_t, std::tuple<double, double, int, int>> decay_map;
    std::map<int64_t, std::tuple<double,int, int>> germanium_map;

    std::multimap<int64_t, std::tuple<double, double, int, int, EventType>> event_map;
    std::multimap<int64_t, std::tuple<double, double, int, int, EventType>> gated_implants;
    std::multimap<int64_t, std::tuple<double, double, int, int, EventType>> all_implants;
    std::multimap<int64_t, std::tuple<double, double, int, int, EventType>> good_impdecays;
    std::multimap<int64_t, std::tuple<double, double, int, int, EventType>> good_implants;
    std::multimap<int64_t, std::tuple<double, double, int, EventType>> good_decays;

    // Loop over the implant and decay trees and fill the maps
    std::cout << "Start filled the maps!" << std::endl;

    // Read gated implant events
    while (gatedimplant_reader.Next()) {
        gated_implants.emplace(*gatedimplant_time, std::make_tuple(*gatedimplant_x, *gatedimplant_y, *gatedimplant_spill, *gatedimplant_dssd, IMPLANT));
    }
    std::cout << "Finished filling the gated implant map" << std::endl;
    std::cout << "Number of implant events: " << gated_implants.size() << std::endl;

    // Read implant events
    while (implant_reader.Next()) {
        // if(*implant_x > 250 && *implant_x < 350 && *implant_y > 40 && *implant_y < 110){
        all_implants.emplace(*implant_time, std::make_tuple(*implant_x, *implant_y, *implant_spill, *implant_dssd, IMPLANT));
        // }
    }
    std::cout << "Finished filling the implant map" << std::endl;
    std::cout << "Number of All implant events cut on region :" << all_implants.size() << std::endl;


    // Read decay events
    while (decay_reader.Next()) {
        // if(*decay_x > 250 && *decay_x < 350 && *decay_y > 40 && *decay_y < 110){
            good_decays.emplace(*decay_time, std::make_tuple(*decay_x, *decay_y, *decay_spill/*, decay_bplast*/, DECAY));
        // }
    }

    std::cout << "Finished filling the decay map" << std::endl;
    std::cout << "Number of Decay events :" << good_decays.size() << std::endl;

    int64_t last_gatedimplant_time;
    double gatedimplant_pos_x;
    double gatedimplant_pos_y;
    int matched_implantdecays = 0;


    // Loop over the gated implant events and for each event find decays in the same position up to 150 seconds. If another implant is found in the same position, stop the correlation.
    for (auto it = gated_implants.begin(); it != gated_implants.end(); it++) {
        auto [x, y, spill, dssd, type] = it->second;

        if (type == IMPLANT && dssd == 1) {
            last_gatedimplant_time = it->first;
            // if(x == 128) continue; // noisy strip
            gatedimplant_pos_x = x;
            gatedimplant_pos_y = y;
            aida_implant_xy->Fill(x,y);
            
            //************* DEGUB **************
            /*std::cout << last_gatedimplant_time -  << std::endl;*/
            //************* DEGUB **************

            // Now loop over the decay events and find the ones that are in the same position as the gated implant up until 150s
            auto decay_start = good_decays.lower_bound(last_gatedimplant_time - 10e3);
            for(auto decay_it = decay_start; decay_it != good_decays.end(); decay_it++){
                auto [decay_x, decay_y, decay_spill/*, decay_bplast*/, decay_type] = decay_it->second;
                //************* DEGUB **************
                /*std::cout << last_gatedimplant_time - decay_it->first << std::endl;*/
                /*std::cout << decay_spill << std::endl;*/
                //************* DEGUB **************
                 
                // std::cout << "shit bricks!" << std::endl;
                if(decay_x == 128 || decay_x == 191 || decay_x == 192) continue;
                if(constants::ONLY_OFFSPILL_DECAY && decay_spill == 1) continue;
                aida_decay_xy->Fill(decay_x,decay_y);
                if (decay_type == DECAY && TMath::Abs(decay_x - gatedimplant_pos_x) <= constants::POSITION_THRESHOLD && TMath::Abs(decay_y - gatedimplant_pos_y) <= constants::POSITION_THRESHOLD){ {
                        // std::cout << "wow, more bricks!" << std::endl;
                        int64_t time_diff = decay_it->first - last_gatedimplant_time;
          
                        //************* DEGUB **************
                        /*std::cout << time_diff << std::endl;*/
                        //************* DEGUB **************
          
                        if (time_diff < constants::TIME_THRESHOLD) {
                            // std::cout << "wow, filling bricks" << std::endl;
                            aida_wr_times->Fill(decay_it->first);
                            aida_matched_xy->Fill(decay_x,decay_y);
                            aida_implant_veto_dt->Fill(time_diff/constants::TIME_SCALE);
                            matched_implantdecays++;
                            break;
                        }
                        else{
                            break;
                        }
                    }
                }
            }

        }
    }

    std::cout << "Finished processing the data" << std::endl;
    std::cout << "Matched: " << matched_implantdecays << " out of " << gated_implants.size() << " gated implant events" << std::endl;
    std::cout << "Matched: " << matched_implantdecays << " out of " << all_implants.size() << " implant events" << std::endl;

    // aida_gatedimplant_before_decay_counter->Write();
    // aida_implant_before_decay_counter->Write();
    // germanium_aida_wr_dt->Write();
    // germanium_energy_hist->Write();
    aida_implant_xy->Write();
    // aida_implant_decay_time->Write();
    aida_wr_times->Write();
    aida_matched_xy->Write();
    // germanium_decay_energy->Write();
    aida_decay_xy->Write();
    // aida_averageimplant_150s->Write();
    // aida_averagegatedimplant_150s->Write();
    aida_implant_veto_dt->Write();
    // aida_implant_dt_vs_pos_x->Write();
    // aida_implant_dt_vs_pos_y->Write();
    // aida_implant_dt_vs_pos_diff->Write();

    std::cout << "Finished writing the histograms" << std::endl;

    // Close the file
    delete file;
    delete outputFile;

    std::cout << "Finished closing the files" << std::endl;

}
