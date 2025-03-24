void plot_onoff_spillstructure(const char* input){
  
  // Open file
  TFile* file = TFile::Open(input, "update");
  if (!file){
    std::cerr << "Error: Could not open file." << std::endl;
    std::exit(1);
  }

  // Check if lifetime histogram exists and fit 
  TH1F* h1_aida_implant_beta_onspillstructure = (TH1F*)file->Get("aida_implant_beta_onspillstructure");
  if (!h1_aida_implant_beta_onspillstructure){
    std::cerr << "Error: Could not find 'aida_implant_beta_onspillstructure' histogram." << std::endl;
    std::exit(1);
  }

  TH1F* h1_aida_implant_beta_offspillstructure = (TH1F*)file->Get("aida_implant_beta_offspillstructure");
  if (!h1_aida_implant_beta_offspillstructure){
    std::cerr << "Error: Could not find 'aida_implant_beta_offspillstructure' histogram." << std::endl;
    std::exit(1);
  }

  TCanvas* c_onoff_spillstructure = new TCanvas("onoff_spillstructure", "", 600,400);

  h1_aida_implant_beta_onspillstructure->SetLineColor(kBlue);
  /*h1_aida_implant_beta_onspillstructure->SetLineWidth(3);*/
  h1_aida_implant_beta_onspillstructure->Draw();

  h1_aida_implant_beta_offspillstructure->SetLineColor(kRed);
  /*h1_aida_implant_beta_offspillstructure->SetLineWidth(3);*/
  h1_aida_implant_beta_offspillstructure->Draw("SAME");

  c_onoff_spillstructure->Write();

  std::cout << "Finished writing the canvases" << std::endl;


  // Clean up and close 
  file->Close();
  delete file;

  std::cout << "Macro has succesfully run!" << std::endl;

  std::exit(0);
}
