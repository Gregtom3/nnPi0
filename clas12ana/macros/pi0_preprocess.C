int pi0_preprocess(
		   const char * input_file = "MC_3053_1.root",
		   const char * output_file = "MC_3053_1_preprocess.root"
){

  
  TFile *fOut = new TFile(output_file,"RECREATE");
  TTree *tOut = new TTree("PreProcessedEvents","PreProcessedEvents");

  const int Nmax = 100;
  int _nPart=0;
  float _Etarel[Nmax];
  float _Phirel[Nmax];
  float _P[Nmax];
  int   _pid[Nmax];
  float _eta[Nmax];
  float _phi[Nmax];
  float _vz[Nmax];
  float _chi2[Nmax];
  int   _MCmatch_flag[Nmax];
  float _pcal_energy[Nmax];
  float _ecin_energy[Nmax];
  float _ecout_energy[Nmax];

  tOut->Branch("nPart",&_nPart,"nPart/I");
  tOut->Branch("Etarel",&_Etarel,"Etarel[nPart]/F");
  tOut->Branch("Phirel",&_Phirel,"Phirel[nPart]/F");
  tOut->Branch("P",&_P,"P[nPart]/F");
  tOut->Branch("pid",&_pid,"pid[nPart]/I");
  tOut->Branch("chi2",&_chi2,"chi2[nPart]/F");
  tOut->Branch("eta",&_eta,"eta[nPart]/F");
  tOut->Branch("phi",&_phi,"phi[nPart]/F");
  tOut->Branch("MCmatch_flag",&_MCmatch_flag,"MCmatch_flag[nPart]/I");
  tOut->Branch("pcal_energy",&_pcal_energy,"pcal_energy[nPart]/F");
  tOut->Branch("ecin_energy",&_ecin_energy,"ecin_energy[nPart]/F");
  tOut->Branch("ecout_energy",&_ecout_energy,"ecout_energy[nPart]/F");

  TChain *chain = new TChain("RawEvents");
  chain->Add(input_file);

  TTreeReader tr(chain);

  TTreeReaderArray<Int_t> tr_ievent(tr, "ievent");
  TTreeReaderArray<Int_t> tr_nPart(tr, "nPart");
  TTreeReaderArray<Int_t> tr_nPartMatch(tr, "nPartMatch");

  TTreeReaderArray<Int_t> tr_MCmatch_flag(tr, "MCmatch_flag");
  TTreeReaderArray<Int_t> tr_MCmatch_parent_id(tr, "MCmatch_parent_id");
  TTreeReaderArray<Int_t> tr_MCmatch_parent_pid(tr, "MCmatch_parent_pid");

  TTreeReaderArray<Int_t> tr_pid(tr, "pid");
  TTreeReaderArray<Float_t> tr_px(tr, "px");
  TTreeReaderArray<Float_t> tr_py(tr, "py");
  TTreeReaderArray<Float_t> tr_pz(tr, "pz");
  TTreeReaderArray<Float_t> tr_E(tr, "E");
  TTreeReaderArray<Float_t> tr_theta(tr, "theta");
  TTreeReaderArray<Float_t> tr_eta(tr, "eta");
  TTreeReaderArray<Float_t> tr_phi(tr, "phi");
  TTreeReaderArray<Float_t> tr_vz(tr, "vz");
  TTreeReaderArray<Float_t> tr_chi2(tr, "chi2");
  TTreeReaderArray<Int_t> tr_pcal_sector(tr, "pcal_sector");
  TTreeReaderArray<Int_t> tr_ecin_sector(tr, "ecin_sector");
  TTreeReaderArray<Int_t> tr_ecout_sector(tr, "ecout_sector");
  TTreeReaderArray<Float_t> tr_pcal_energy(tr, "pcal_energy");
  TTreeReaderArray<Float_t> tr_ecin_energy(tr, "ecin_energy");
  TTreeReaderArray<Float_t> tr_ecout_energy(tr, "ecout_energy");
  TTreeReaderArray<Float_t> tr_pcal_lu(tr, "pcal_lu");
  TTreeReaderArray<Float_t> tr_ecin_lu(tr, "ecin_lu");
  TTreeReaderArray<Float_t> tr_ecout_lu(tr, "ecout_lu");
  TTreeReaderArray<Float_t> tr_pcal_lv(tr, "pcal_lv");
  TTreeReaderArray<Float_t> tr_ecin_lv(tr, "ecin_lv");
  TTreeReaderArray<Float_t> tr_ecout_lv(tr, "ecout_lv");
  

  
  // Loop over all RawEvents in TTree
  while(tr.Next()){

    // Get number of particles in event
    _nPart = tr_nPart[0];
    
    
    // Loop over particles in the event
    for(int i = 0 ; i < _nPart; i++){
      
      int pid_center = tr_pid[i];
      float eta_center = tr_eta[i];
      float phi_center = tr_eta[i];

      // Nested loop over photons
      if(pid_center==22){
    
	for(int j = 0 ; j < _nPart; j++){
	  int pid       = tr_pid[j];
	  float eta     = tr_eta[j];
	  float phi     = tr_phi[j];
	  float px      = tr_px[j];
	  float py      = tr_py[j];
	  float pz      = tr_pz[j];
	  float p       = sqrt(px*px + py*py + pz*pz);
	  float vz      = tr_vz[j];
	  float chi2    = tr_chi2[j];
	  int MCmatch_flag = tr_MCmatch_flag[j];
	  float pcal_e  = tr_pcal_energy[j];
	  float ecin_e  = tr_ecin_energy[j];
	  float ecout_e  = tr_ecout_energy[j];

	  float etarel = eta_center-eta;
	  float phirel = min((double)abs(phi_center-phi),(double)(2*M_PI-abs(phi_center-phi))); 

	  _Etarel[j] = etarel;
	  _Phirel[j] = phirel;
	  _P[j] = p;
	  _pid[j] = pid;
	  _eta[j] = eta;
	  _phi[j] = phi;
	  _vz[j] = vz;
	  _chi2[j] = chi2;
	  _MCmatch_flag[j] = MCmatch_flag;
	  _pcal_energy[j] = pcal_e;
	  _ecin_energy[j] = ecin_e;
	  _ecout_energy[j] = ecout_e;
	  
	}
	tOut->Fill();
      }

    }

  }
  
  tOut->Write();
  fOut->Close();
  return 0;
}
