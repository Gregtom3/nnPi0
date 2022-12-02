int pi0_preprocess_particleNet(
		   const char * input_file = "MC_3053_2.root",
		   const char * output_file = "MC_3053_2_preprocess.root"
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
  int   _MCmatch_flag;
  float _pcal_energy[Nmax];
  float _ecin_energy[Nmax];
  float _ecout_energy[Nmax];
  float _pcal_lu[Nmax];
  float _pcal_lv[Nmax];
  float _pcal_lw[Nmax];
  float _pcal_du[Nmax];
  float _pcal_dv[Nmax];
  float _pcal_dw[Nmax];
  float _pcal_m2u[Nmax];
  float _pcal_m2v[Nmax];
  float _pcal_m2w[Nmax];
  float _pcal_m3u[Nmax];
  float _pcal_m3v[Nmax];
  float _pcal_m3w[Nmax];
  float _ecin_lu[Nmax];
  float _ecin_lv[Nmax];
  float _ecin_lw[Nmax];
  float _ecin_du[Nmax];
  float _ecin_dv[Nmax];
  float _ecin_dw[Nmax];
  float _ecin_m2u[Nmax];
  float _ecin_m2v[Nmax];
  float _ecin_m2w[Nmax];
  float _ecin_m3u[Nmax];
  float _ecin_m3v[Nmax];
  float _ecin_m3w[Nmax];
  float _ecout_lu[Nmax];
  float _ecout_lv[Nmax];
  float _ecout_lw[Nmax];
  float _ecout_du[Nmax];
  float _ecout_dv[Nmax];
  float _ecout_dw[Nmax];
  float _ecout_m2u[Nmax];
  float _ecout_m2v[Nmax];
  float _ecout_m2w[Nmax];
  float _ecout_m3u[Nmax];
  float _ecout_m3v[Nmax];
  float _ecout_m3w[Nmax];

  tOut->Branch("nPart",&_nPart,"nPart/I");
  tOut->Branch("Etarel",&_Etarel,"Etarel[nPart]/F");
  tOut->Branch("Phirel",&_Phirel,"Phirel[nPart]/F");
  tOut->Branch("P",&_P,"P[nPart]/F");
  tOut->Branch("pid",&_pid,"pid[nPart]/I");
  tOut->Branch("chi2",&_chi2,"chi2[nPart]/F");
  tOut->Branch("eta",&_eta,"eta[nPart]/F");
  tOut->Branch("phi",&_phi,"phi[nPart]/F");
  tOut->Branch("MCmatch_flag",&_MCmatch_flag,"MCmatch_flag/I");
  tOut->Branch("pcal_energy",&_pcal_energy,"pcal_energy[nPart]/F");
  tOut->Branch("pcal_lu",&_pcal_lu,"pcal_lu[nPart]/F");
  tOut->Branch("pcal_lv",&_pcal_lv,"pcal_lv[nPart]/F");
  tOut->Branch("pcal_lw",&_pcal_lw,"pcal_lw[nPart]/F");
  tOut->Branch("pcal_du",&_pcal_du,"pcal_du[nPart]/F");
  tOut->Branch("pcal_dv",&_pcal_dv,"pcal_dv[nPart]/F");
  tOut->Branch("pcal_dw",&_pcal_dw,"pcal_dw[nPart]/F");
  tOut->Branch("pcal_m2u",&_pcal_m2u,"pcal_m2u[nPart]/F");
  tOut->Branch("pcal_m2v",&_pcal_m2v,"pcal_m2v[nPart]/F");
  tOut->Branch("pcal_m2w",&_pcal_m2w,"pcal_m2w[nPart]/F");
  tOut->Branch("pcal_m3u",&_pcal_m3u,"pcal_m3u[nPart]/F");
  tOut->Branch("pcal_m3v",&_pcal_m3v,"pcal_m3v[nPart]/F");
  tOut->Branch("pcal_m3w",&_pcal_m3w,"pcal_m3w[nPart]/F");
  tOut->Branch("ecin_energy",&_ecin_energy,"ecin_energy[nPart]/F");
  tOut->Branch("ecin_lu",&_ecin_lu,"ecin_lu[nPart]/F");
  tOut->Branch("ecin_lv",&_ecin_lv,"ecin_lv[nPart]/F");
  tOut->Branch("ecin_lw",&_ecin_lw,"ecin_lw[nPart]/F");
  tOut->Branch("ecin_du",&_ecin_du,"ecin_du[nPart]/F");
  tOut->Branch("ecin_dv",&_ecin_dv,"ecin_dv[nPart]/F");
  tOut->Branch("ecin_dw",&_ecin_dw,"ecin_dw[nPart]/F");
  tOut->Branch("ecin_m2u",&_ecin_m2u,"ecin_m2u[nPart]/F");
  tOut->Branch("ecin_m2v",&_ecin_m2v,"ecin_m2v[nPart]/F");
  tOut->Branch("ecin_m2w",&_ecin_m2w,"ecin_m2w[nPart]/F");
  tOut->Branch("ecin_m3u",&_ecin_m3u,"ecin_m3u[nPart]/F");
  tOut->Branch("ecin_m3v",&_ecin_m3v,"ecin_m3v[nPart]/F");
  tOut->Branch("ecin_m3w",&_ecin_m3w,"ecin_m3w[nPart]/F");
  tOut->Branch("ecout_energy",&_ecout_energy,"ecout_energy[nPart]/F");
  tOut->Branch("ecout_lu",&_ecout_lu,"ecout_lu[nPart]/F");
  tOut->Branch("ecout_lv",&_ecout_lv,"ecout_lv[nPart]/F");
  tOut->Branch("ecout_lw",&_ecout_lw,"ecout_lw[nPart]/F");
  tOut->Branch("ecout_du",&_ecout_du,"ecout_du[nPart]/F");
  tOut->Branch("ecout_dv",&_ecout_dv,"ecout_dv[nPart]/F");
  tOut->Branch("ecout_dw",&_ecout_dw,"ecout_dw[nPart]/F");
  tOut->Branch("ecout_m2u",&_ecout_m2u,"ecout_m2u[nPart]/F");
  tOut->Branch("ecout_m2v",&_ecout_m2v,"ecout_m2v[nPart]/F");
  tOut->Branch("ecout_m2w",&_ecout_m2w,"ecout_m2w[nPart]/F");
  tOut->Branch("ecout_m3u",&_ecout_m3u,"ecout_m3u[nPart]/F");
  tOut->Branch("ecout_m3v",&_ecout_m3v,"ecout_m3v[nPart]/F");
  tOut->Branch("ecout_m3w",&_ecout_m3w,"ecout_m3w[nPart]/F");

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
  TTreeReaderArray<Float_t> tr_pcal_lu(tr, "pcal_lu");
  TTreeReaderArray<Float_t> tr_pcal_lv(tr, "pcal_lv");
  TTreeReaderArray<Float_t> tr_pcal_lw(tr, "pcal_lw");
  TTreeReaderArray<Float_t> tr_pcal_du(tr, "pcal_du");
  TTreeReaderArray<Float_t> tr_pcal_dv(tr, "pcal_dv");
  TTreeReaderArray<Float_t> tr_pcal_dw(tr, "pcal_dw");
  TTreeReaderArray<Float_t> tr_pcal_m2u(tr, "pcal_m2u");
  TTreeReaderArray<Float_t> tr_pcal_m2v(tr, "pcal_m2v");
  TTreeReaderArray<Float_t> tr_pcal_m2w(tr, "pcal_m2w");
  TTreeReaderArray<Float_t> tr_pcal_m3u(tr, "pcal_m3u");
  TTreeReaderArray<Float_t> tr_pcal_m3v(tr, "pcal_m3v");
  TTreeReaderArray<Float_t> tr_pcal_m3w(tr, "pcal_m3w");
  TTreeReaderArray<Float_t> tr_ecin_energy(tr, "ecin_energy");
  TTreeReaderArray<Float_t> tr_ecin_lu(tr, "ecin_lu");
  TTreeReaderArray<Float_t> tr_ecin_lv(tr, "ecin_lv");
  TTreeReaderArray<Float_t> tr_ecin_lw(tr, "ecin_lw");
  TTreeReaderArray<Float_t> tr_ecin_du(tr, "ecin_du");
  TTreeReaderArray<Float_t> tr_ecin_dv(tr, "ecin_dv");
  TTreeReaderArray<Float_t> tr_ecin_dw(tr, "ecin_dw");
  TTreeReaderArray<Float_t> tr_ecin_m2u(tr, "ecin_m2u");
  TTreeReaderArray<Float_t> tr_ecin_m2v(tr, "ecin_m2v");
  TTreeReaderArray<Float_t> tr_ecin_m2w(tr, "ecin_m2w");
  TTreeReaderArray<Float_t> tr_ecin_m3u(tr, "ecin_m3u");
  TTreeReaderArray<Float_t> tr_ecin_m3v(tr, "ecin_m3v");
  TTreeReaderArray<Float_t> tr_ecin_m3w(tr, "ecin_m3w");
  TTreeReaderArray<Float_t> tr_ecout_energy(tr, "ecout_energy");
  TTreeReaderArray<Float_t> tr_ecout_lu(tr, "ecout_lu");
  TTreeReaderArray<Float_t> tr_ecout_lv(tr, "ecout_lv");
  TTreeReaderArray<Float_t> tr_ecout_lw(tr, "ecout_lw");
  TTreeReaderArray<Float_t> tr_ecout_du(tr, "ecout_du");
  TTreeReaderArray<Float_t> tr_ecout_dv(tr, "ecout_dv");
  TTreeReaderArray<Float_t> tr_ecout_dw(tr, "ecout_dw");
  TTreeReaderArray<Float_t> tr_ecout_m2u(tr, "ecout_m2u");
  TTreeReaderArray<Float_t> tr_ecout_m2v(tr, "ecout_m2v");
  TTreeReaderArray<Float_t> tr_ecout_m2w(tr, "ecout_m2w");
  TTreeReaderArray<Float_t> tr_ecout_m3u(tr, "ecout_m3u");
  TTreeReaderArray<Float_t> tr_ecout_m3v(tr, "ecout_m3v");
  TTreeReaderArray<Float_t> tr_ecout_m3w(tr, "ecout_m3w");
  

  
  // Loop over all RawEvents in TTree
  while(tr.Next()){

    // Get number of particles in event
    _nPart = tr_nPart[0];
    
    
    // Loop over particles in the event
    for(int i = 0 ; i < _nPart; i++){
      
      int pid_center = tr_pid[i];
      float eta_center = tr_eta[i];
      double phi_center = tr_eta[i];
      
      // Flag if the particle has MC parent
      _MCmatch_flag = tr_MCmatch_flag[i];


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

	  float pcal_e  = tr_pcal_energy[j];
	  float ecin_e  = tr_ecin_energy[j];
	  float ecout_e  = tr_ecout_energy[j];

	  float etarel = eta_center-eta;
	  float phirel = min(abs(phi_center-phi),(2*3.14159265-abs(phi_center-phi))); 

	  _Etarel[j] = etarel;
	  _Phirel[j] = phirel;
	  _P[j] = p;
	  _pid[j] = pid;
	  _eta[j] = eta;
	  _phi[j] = phi;
	  _vz[j] = vz;
	  _chi2[j] = chi2;
	  _pcal_energy[j] = tr_pcal_energy[j];
	  _pcal_lu[j] = tr_pcal_lu[j];
	  _pcal_lv[j] = tr_pcal_lv[j];
	  _pcal_lw[j] = tr_pcal_lw[j];
	  _pcal_du[j] = tr_pcal_du[j];
	  _pcal_dv[j] = tr_pcal_dv[j];
	  _pcal_dw[j] = tr_pcal_dw[j];
	  _pcal_m2u[j] = tr_pcal_m2u[j];
	  _pcal_m2v[j] = tr_pcal_m2v[j];
	  _pcal_m2w[j] = tr_pcal_m2w[j];
	  _pcal_m3u[j] = tr_pcal_m3u[j];
	  _pcal_m3v[j] = tr_pcal_m3v[j];
	  _pcal_m3w[j] = tr_pcal_m3w[j];
	  _ecin_energy[j] = tr_ecin_energy[j];
	  _ecin_lu[j] = tr_ecin_lu[j];
	  _ecin_lv[j] = tr_ecin_lv[j];
	  _ecin_lw[j] = tr_ecin_lw[j];
	  _ecin_du[j] = tr_ecin_du[j];
	  _ecin_dv[j] = tr_ecin_dv[j];
	  _ecin_dw[j] = tr_ecin_dw[j];
	  _ecin_m2u[j] = tr_ecin_m2u[j];
	  _ecin_m2v[j] = tr_ecin_m2v[j];
	  _ecin_m2w[j] = tr_ecin_m2w[j];
	  _ecin_m3u[j] = tr_ecin_m3u[j];
	  _ecin_m3v[j] = tr_ecin_m3v[j];
	  _ecin_m3w[j] = tr_ecin_m3w[j];
	  _ecout_energy[j] = tr_ecout_energy[j];
	  _ecout_lu[j] = tr_ecout_lu[j];
	  _ecout_lv[j] = tr_ecout_lv[j];
	  _ecout_lw[j] = tr_ecout_lw[j];
	  _ecout_du[j] = tr_ecout_du[j];
	  _ecout_dv[j] = tr_ecout_dv[j];
	  _ecout_dw[j] = tr_ecout_dw[j];
	  _ecout_m2u[j] = tr_ecout_m2u[j];
	  _ecout_m2v[j] = tr_ecout_m2v[j];
	  _ecout_m2w[j] = tr_ecout_m2w[j];
	  _ecout_m3u[j] = tr_ecout_m3u[j];
	  _ecout_m3v[j] = tr_ecout_m3v[j];
	  _ecout_m3w[j] = tr_ecout_m3w[j];
	  
	}
	tOut->Fill();
      }

    }

  }
  
  tOut->Write();
  fOut->Close();
  return 0;
}
