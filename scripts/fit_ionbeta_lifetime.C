namespace scriptVariables {
  // Define fit range
  float fitMin = 0.;
  float fitMax = 50.;
  // Define initial parameters for fit 
  float lifetime = 3.2;
  float amplitude = 80.;
  float background = 5.;

  // Define fitting calls
  int Npx = 10000;
}

void fit_ionbeta_lifetime(const char* input){
  
  // Open file
  TFile* file = TFile::Open(input, "update");
  if (!file){
    std::cerr << "Error: Could not open file." << std::endl;
    std::exit(1);
  }

  // Define fitting function
  TF1* fit = new TF1("fit", "[bkg]+[A]*exp(-x*log(2)/[t_1/2])", scriptVariables::fitMin, scriptVariables::fitMax);
  
  // Set fitting calls
  fit->SetNpx(scriptVariables::Npx);
  
  // Set initial parameters
  fit->SetParameters(scriptVariables::background, scriptVariables::amplitude, scriptVariables::lifetime);
  
  // Check if lifetime histogram exists and fit 
  TH1F* histogram = (TH1F*)file->Get("aida_implant_veto_dt");
  if (!histogram){
    std::cerr << "Error: Could not find 'aida_implant_veto_dt' histogram." << std::endl;
    std::exit(1);
  }

  // Check if fitted objects already exist and delete them
  TH1F* old_histogram_fitted = (TH1F*)file->Get("aida_implant_veto_dt_fitted");
  gDirectory->Purge();

  if (old_histogram_fitted){
    gDirectory->Delete("aida_implant_veto_dt_fitted;1");
  }
  delete old_histogram_fitted;

  TF1* old_fit = (TF1*)file->Get("fit");
  if (old_fit){
    gDirectory->Delete("fit;1");
  }
  delete old_fit;

  // Make a clone of the histogram to fit
  TH1F* histogram_fitted = (TH1F*)histogram->Clone("aida_implant_veto_dt_fitted");

  // Fit histogram
  histogram_fitted->Fit("fit", "WWM", "", scriptVariables::fitMin, scriptVariables::fitMax);


  std::cout << std::endl;
  
  // Write new fitted hist and fit function
  fit->Write();
  histogram_fitted->Write();

  // Clean up and close 
  file->Close();
  delete file;

  std::cout << "Macro has succesfully run!" << std::endl;

  std::exit(0);
}
