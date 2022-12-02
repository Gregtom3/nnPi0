#include "../src/SIDISParticle.h"
#include "../src/SIDISParticlev1.h"
typedef std::map<int, SIDISParticle*> type_map_part;
#include "../src/PID.C"
#include "../src/Constants.h"
#include "../src/Kinematics.h"
#include "../src/FiducialCuts.h"
#include "../src/HipoBankInterface.h"

void DeleteParticlePointers(type_map_part& map){
  for(type_map_part::iterator it = map.begin(); it!=map.end(); ++it){
    delete (it->second);
    map.erase(it);
  }
  return;
}
int pi0_readHipo(const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3053_2.hipo",
		const char * outputFile = "MC_3053_2.root",
		const double _electron_beam_energy = 10.6,
		const int maxEvents = 10000,
		bool hipo_is_mc = true){

 
  // Open TTree and declare branches
  // -------------------------------------
  TFile *fOut = new TFile(outputFile,"RECREATE");
  TTree *tree = new TTree("RawEvents","RawEvents");

  // Initialize important event information
  // --------------------------------------
  double x;  // Monte Carlo
  double y;
  double Q2;
  double nu;
  double W;
  const double Mp = 0.938272;
  const double Me = 0.000511;
  double s = pow(Mp,2)+pow(Me,2)+2*Mp*_electron_beam_energy;

  float reco_x;  // From Reconstructed Electron
  float reco_y;
  float reco_Q2;
  float reco_nu;
  float reco_W;

  // Maximum number of particles
  const int Nmax = 100;

  // Global variables
  int   _ievent=0;
  int   _nPart=0;
  int   _nPartMatch=0;
  
  // Particle kinematics
  float _px[Nmax];
  float _py[Nmax];
  float _pz[Nmax];
  float _E[Nmax];
  float _theta[Nmax];
  float _eta[Nmax];
  float _phi[Nmax];
  float _vz[Nmax];

  // Additional REC::Particle information
  int   _pid[Nmax];
  float _beta[Nmax];
  float _chi2[Nmax];

  // MC matching variables
  int   _MCmatch_flag[Nmax];
  int   _MCmatch_parent_id[Nmax];
  int   _MCmatch_parent_pid[Nmax];
  
  // REC::Calorimeter variables
  int   _pcal_sector[Nmax];
  int   _ecin_sector[Nmax];
  int   _ecout_sector[Nmax];
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




  // TTree Branching
  tree->Branch("x",&reco_x,"reco_x/F");
  tree->Branch("Q2",&reco_Q2,"reco_Q2/F");
  tree->Branch("W",&reco_W,"reco_W/F");
  tree->Branch("y",&reco_y,"reco_y/F");
  tree->Branch("ievent",&_ievent,"ievent/I");
  tree->Branch("nPart",&_nPart,"nPart/I");
  tree->Branch("nPartMatch",&_nPartMatch,"nPartMatch/I");
  tree->Branch("px",&_px,"px[nPart]/F");
  tree->Branch("py",&_py,"py[nPart]/F");
  tree->Branch("pz",&_pz,"pz[nPart]/F");
  tree->Branch("E",&_E,"E[nPart]/F");
  tree->Branch("theta",&_theta,"theta[nPart]/F");
  tree->Branch("eta",&_eta,"eta[nPart]/F");
  tree->Branch("phi",&_phi,"phi[nPart]/F");
  tree->Branch("vz",&_vz,"vz[nPart]/F");
  tree->Branch("pid",&_pid,"pid[nPart]/I");
  tree->Branch("beta",&_beta,"beta[nPart]/F");
  tree->Branch("chi2",&_chi2,"chi2[nPart]/F");
  tree->Branch("MCmatch_flag",&_MCmatch_flag,"MCmatch_flag[nPart]/I");
  tree->Branch("MCmatch_parent_id",&_MCmatch_parent_id,"MCmatch_parent_id[nPart]/I");
  tree->Branch("MCmatch_parent_pid",&_MCmatch_parent_pid,"MCmatch_parent_pid[nPart]/I");
  tree->Branch("pcal_sector",&_pcal_sector,"pcal_sector[nPart]/I");
  tree->Branch("ecin_sector",&_ecin_sector,"ecin_sector[nPart]/I");
  tree->Branch("ecout_sector",&_ecout_sector,"ecout_sector[nPart]/I");
  tree->Branch("pcal_energy",&_pcal_energy,"pcal_energy[nPart]/F");
  tree->Branch("ecin_energy",&_ecin_energy,"ecin_energy[nPart]/F");
  tree->Branch("ecout_energy",&_ecout_energy,"ecout_energy[nPart]/F");
  tree->Branch("pcal_lu",&_pcal_lu,"pcal_lu[nPart]/F");
  tree->Branch("pcal_lv",&_pcal_lv,"pcal_lv[nPart]/F");
  tree->Branch("pcal_lw",&_pcal_lw,"pcal_lw[nPart]/F");
  tree->Branch("pcal_du",&_pcal_du,"pcal_du[nPart]/F");
  tree->Branch("pcal_dv",&_pcal_dv,"pcal_dv[nPart]/F");
  tree->Branch("pcal_dw",&_pcal_dw,"pcal_dw[nPart]/F");
  tree->Branch("pcal_m2u",&_pcal_m2u,"pcal_m2u[nPart]/F");
  tree->Branch("pcal_m2v",&_pcal_m2v,"pcal_m2v[nPart]/F");
  tree->Branch("pcal_m2w",&_pcal_m2w,"pcal_m2w[nPart]/F");
  tree->Branch("pcal_m3u",&_pcal_m3u,"pcal_m3u[nPart]/F");
  tree->Branch("pcal_m3v",&_pcal_m3v,"pcal_m3v[nPart]/F");
  tree->Branch("pcal_m3w",&_pcal_m3w,"pcal_m3w[nPart]/F");

  tree->Branch("ecin_lu",&_ecin_lu,"ecin_lu[nPart]/F");
  tree->Branch("ecin_lv",&_ecin_lv,"ecin_lv[nPart]/F");
  tree->Branch("ecin_lw",&_ecin_lw,"ecin_lw[nPart]/F");
  tree->Branch("ecin_du",&_ecin_du,"ecin_du[nPart]/F");
  tree->Branch("ecin_dv",&_ecin_dv,"ecin_dv[nPart]/F");
  tree->Branch("ecin_dw",&_ecin_dw,"ecin_dw[nPart]/F");
  tree->Branch("ecin_m2u",&_ecin_m2u,"ecin_m2u[nPart]/F");
  tree->Branch("ecin_m2v",&_ecin_m2v,"ecin_m2v[nPart]/F");
  tree->Branch("ecin_m2w",&_ecin_m2w,"ecin_m2w[nPart]/F");
  tree->Branch("ecin_m3u",&_ecin_m3u,"ecin_m3u[nPart]/F");
  tree->Branch("ecin_m3v",&_ecin_m3v,"ecin_m3v[nPart]/F");
  tree->Branch("ecin_m3w",&_ecin_m3w,"ecin_m3w[nPart]/F");

  tree->Branch("ecout_lu",&_ecout_lu,"ecout_lu[nPart]/F");
  tree->Branch("ecout_lv",&_ecout_lv,"ecout_lv[nPart]/F");
  tree->Branch("ecout_lw",&_ecout_lw,"ecout_lw[nPart]/F");
  tree->Branch("ecout_du",&_ecout_du,"ecout_du[nPart]/F");
  tree->Branch("ecout_dv",&_ecout_dv,"ecout_dv[nPart]/F");
  tree->Branch("ecout_dw",&_ecout_dw,"ecout_dw[nPart]/F");
  tree->Branch("ecout_m2u",&_ecout_m2u,"ecout_m2u[nPart]/F");
  tree->Branch("ecout_m2v",&_ecout_m2v,"ecout_m2v[nPart]/F");
  tree->Branch("ecout_m2w",&_ecout_m2w,"ecout_m2w[nPart]/F");
  tree->Branch("ecout_m3u",&_ecout_m3u,"ecout_m3u[nPart]/F");
  tree->Branch("ecout_m3v",&_ecout_m3v,"ecout_m3v[nPart]/F");
  tree->Branch("ecout_m3w",&_ecout_m3w,"ecout_m3w[nPart]/F");

  // Configure CLAS12 Reader and HipoChain
  // -------------------------------------
  clas12root::HipoChain _chain;
  clas12::clas12reader *_config_c12{nullptr};

  _chain.Add(hipoFile);
  _config_c12=_chain.GetC12Reader();
  _config_c12->db()->turnOffQADB(); // Turn off QADB for Monte Carlo Analysis

  // Configure PIDs for final state
  // -------------------------------------
  _config_c12->addExactPid(11,1);     // Exactly 1 electron
  _config_c12->addAtLeastPid(22,2);   // 2 or more photons

  
  // Add RUN::config bank
  // -------------------------------------
  int _idx_RUNconfig = _config_c12->addBank("RUN::config");
  int _irun = _config_c12->getBankOrder(_idx_RUNconfig,"run");
  int _ievnum = _config_c12->getBankOrder(_idx_RUNconfig,"event");
  int _itorus = _config_c12->getBankOrder(_idx_RUNconfig,"torus");
  int _runNumber = 0;

  // Establish CLAS12 event parser
  // -------------------------------------
  auto &_c12=_chain.C12ref();
  
  // Create CLAS12Analysis Objects
  // -------------------------------------
  HipoBankInterface _hipoInterface = HipoBankInterface(_c12);
  FiducialCuts _fiducial = FiducialCuts();
  Kinematics _kin;
  PID _pidhelper;

  // Create particleMap objects
  // -------------------------------------
  type_map_part recoParticleMap;
  type_map_part mcParticleMap;
  
  // Loop over HipoChain
  // -------------------------------------
  int NumNoEle=0;
  int NumNo2g=0;
  
  int whileidx=0;
  while(_chain.Next()==true && (whileidx < maxEvents || maxEvents < 0)){
    if(whileidx%10000==0 && whileidx!=0){
      std::cout << whileidx << " events read | " << _ievent*100.0/whileidx << "% passed event selection" << std::endl;
    }

    whileidx++;
    // Wipe both particleMaps clean
    // -------------------------------------
    DeleteParticlePointers(recoParticleMap);
    DeleteParticlePointers(mcParticleMap);
   
    // Get run specific information
    // -------------------------------------
    _runNumber = _c12->getBank(_idx_RUNconfig)->getInt(_irun,0);
    _runNumber *= _c12->getBank(_idx_RUNconfig)->getFloat(_itorus,0); // Multiply run number by torus bending

    // Loop over reconstructed particles
    // -------------------------------------------------------
    auto particles=_c12->getDetParticles();
    int l=0;
    bool foundEle=false;
    bool found2g=false;
    for(unsigned int idx = 0 ; idx < particles.size() ; idx++){
      // Extract each particle from event one-at-a-time
      // -------------------------------------------------------
      auto particle = particles.at(idx);
      int pid = particle->getPid();      
      float chi2 = particle->getChi2Pid();
      if(chi2>100 || chi2<-100)
	chi2=0;
      float theta = particle->getTheta();
      float eta = _kin.eta(theta);
      float phi = particle->getPhi();
      float p = particle->getP();
      float px = _kin.Px(p,theta,phi);
      float py = _kin.Py(p,theta,phi);
      float pz = _kin.Pz(p,theta,phi);
      float pt = _kin.Pt(px,py);
      float m = 0.0;
      if(pid!=22)
	m = particle->getPdgMass();
      float E  = _kin.E(m,p);
      float beta = particle->getBeta();
      int pindex = particle->getIndex();
      float vx = particle->par()->getVx();
      float vy = particle->par()->getVy();
      float vz = particle->par()->getVz();
      int status = particle->getStatus();
      
      if(pid==11){//only electron is assumed to be the scattered e
	reco_Q2=_kin.Q2(_electron_beam_energy,E,_kin.cth(px,py,pz));
	reco_y=_kin.y(_electron_beam_energy,E);
	reco_nu=_kin.nu(_electron_beam_energy,E);
	reco_W=_kin.W(reco_Q2,0.938272,reco_nu);
	reco_x=_kin.x(reco_Q2,s,reco_y);
      }
    
      // Create new particle structure
      // -------------------------------------------------------
      SIDISParticlev1 *sp = new SIDISParticlev1();
      sp->set_candidate_id( pindex );
    
      sp->set_property( SIDISParticle::part_pid, pid);
      sp->set_property( SIDISParticle::part_px,  px);
      sp->set_property( SIDISParticle::part_py,  py);
      sp->set_property( SIDISParticle::part_pz,  pz);
      sp->set_property( SIDISParticle::part_pt,  pt);
      sp->set_property( SIDISParticle::part_p,  p);
      sp->set_property( SIDISParticle::part_E,   E);
      sp->set_property( SIDISParticle::part_theta,   theta);
      sp->set_property( SIDISParticle::part_eta,   eta);
      sp->set_property( SIDISParticle::part_phi,   phi);
      sp->set_property( SIDISParticle::part_vx,   vx);
      sp->set_property( SIDISParticle::part_vy,   vy);
      sp->set_property( SIDISParticle::part_vz,   vz);
      sp->set_property( SIDISParticle::part_pindex,   pindex);
      sp->set_property( SIDISParticle::part_beta,   beta);
      sp->set_property( SIDISParticle::part_chi2,   chi2);
      sp->set_property( SIDISParticle::part_ID, pindex);
      sp->set_property( SIDISParticle::part_status, status);
      sp->set_property( SIDISParticle::part_parentID, (int)-999);
      sp->set_property( SIDISParticle::part_parentPID, (int)-999);

      // ------------------- MONTE CARLO ----------------------------//
      sp->set_property( SIDISParticle::evtgen_part_pid, -999);
      sp->set_property( SIDISParticle::evtgen_part_px,  (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_py,  (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_pz,  (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_pt,  (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_p,  (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_E,  (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_theta,   (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_eta,   (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_phi,   (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_vx,   (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_vy,   (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_vz,   (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_pindex,   (int)-999);
      sp->set_property( SIDISParticle::evtgen_part_beta,   (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_chi2,   (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_status,   (int)-999);
      sp->set_property( SIDISParticle::evtgen_part_ID, (int)-999);
      sp->set_property( SIDISParticle::evtgen_part_parentID, (int)-999);
      sp->set_property( SIDISParticle::evtgen_part_parentPID, (int)-999);
      sp->set_property( SIDISParticle::evtgen_part_parentparentPID, (int)-999);
    
      // Add detector info to SIDISParticle
      // --------------------------------------------------------------------------
      //    
      if(_hipoInterface.loadBankData(_c12, sp)==false)
	{delete sp; continue;}

      // CUT Fiducial
      // --------------------------------------------------------------------------
      if(_fiducial.FidCutParticle(_c12,11,sp) == false)
	{delete sp; continue;}
      
    
      // Add SIDISParticle to the collection
      // --------------------------------------------------------------------------

      recoParticleMap.insert ( make_pair( l++ , sp) );

    }

    // Apply particle cuts on the map, obtaining a new map
    // --------------------------------------------------------------------------
    recoParticleMap = _pidhelper.applyCuts(recoParticleMap);

    // Loop over recoParticleMap to make sure we have an electron and 2 gammas
    // --------------------------------------------------------------------------
    int numEle = _pidhelper.countPID(recoParticleMap,11);
    int numGamma = _pidhelper.countPID(recoParticleMap,22);
    
    if(numEle==1)
      foundEle=true;
    if(numGamma>=2)
      found2g=true;

    if(foundEle==false){
      NumNoEle++;
      continue;
    }
    if(found2g==false){
      NumNo2g++;
      continue;
    }
      
    
    // Loop over all Monte Carlo particles
    // -------------------------------------
    auto mcparticles=_c12->mcparts();
    for(int idx = 0 ; idx < mcparticles->getRows() && hipo_is_mc ; idx++){
    
      // Create new SIDISParticle
      SIDISParticlev1 *sp = new SIDISParticlev1();
      sp->set_candidate_id( mcParticleMap.size() );

      int pid = mcparticles->getPid(idx);
      float px = mcparticles->getPx(idx);
      float py = mcparticles->getPy(idx);
      float pz = mcparticles->getPz(idx);
      float m = mcparticles->getMass(idx);
    
      float pt = _kin.Pt(px,py);
      float p  = _kin.P(px,py,pz);
      float E  = _kin.E(m,p);

      float theta = _kin.th(pt,pz);
      float eta = _kin.eta(theta);
      float phi   = _kin.phi(px,py);

      float vx = mcparticles->getVx(idx);
      float vy = mcparticles->getVy(idx);
      float vz = mcparticles->getVz(idx);

      int parentID = mcparticles->getParent(idx)-1;
      int parentparentID = mcparticles->getParent(parentID)-1;
    
      int parentPID = mcparticles->getPid(parentID);
      int parentparentPID = mcparticles->getPid(parentparentID);

      if(pid==11 && parentID==1){ // scattered electron
	Q2=_kin.Q2(_electron_beam_energy,E,_kin.cth(px,py,pz));
	y=_kin.y(_electron_beam_energy,E);
	nu=_kin.nu(_electron_beam_energy,E);
	W=_kin.W(Q2,0.938272,nu);
	x=_kin.x(Q2,s,y);
      }
      
      if(mcparticles->getType(idx)!=1) // Reject non-final state
	{continue;} 
	    
      sp->set_property( SIDISParticle::evtgen_part_pid, pid);
      sp->set_property( SIDISParticle::evtgen_part_px,  px);
      sp->set_property( SIDISParticle::evtgen_part_py,  py);
      sp->set_property( SIDISParticle::evtgen_part_pz,  pz);
      sp->set_property( SIDISParticle::evtgen_part_pt,  pt);
      sp->set_property( SIDISParticle::evtgen_part_p,  p);
      sp->set_property( SIDISParticle::evtgen_part_E,   E);
      sp->set_property( SIDISParticle::evtgen_part_theta,   theta);
      sp->set_property( SIDISParticle::evtgen_part_eta,   eta);
      sp->set_property( SIDISParticle::evtgen_part_phi,   phi);
      sp->set_property( SIDISParticle::evtgen_part_vx,   vx);
      sp->set_property( SIDISParticle::evtgen_part_vy,   vy);
      sp->set_property( SIDISParticle::evtgen_part_vz,   vz);
      sp->set_property( SIDISParticle::evtgen_part_pindex,   (int)-999);
      sp->set_property( SIDISParticle::evtgen_part_beta,   (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_chi2,   (float)-999);
      sp->set_property( SIDISParticle::evtgen_part_status,   (int)-999);
      sp->set_property( SIDISParticle::evtgen_part_ID,   idx);
      sp->set_property( SIDISParticle::evtgen_part_parentID,  (int)parentID);
      sp->set_property( SIDISParticle::evtgen_part_parentPID, (int)parentPID);
      sp->set_property( SIDISParticle::evtgen_part_parentparentPID, (int)parentparentPID);
      // Add SIDISParticle to the collection
      mcParticleMap.insert ( make_pair( sp->get_candidate_id() , sp) );   
    }
  
     
    //
    //
    //   PARTICLE MATCHING
    //
    //
 
    // Loop over all reco particles
    /* Loop over all reco particles */
    for (type_map_part::iterator it_reco = recoParticleMap.begin(); it_reco!= recoParticleMap.end() && hipo_is_mc; ++it_reco){

      double reco_theta = (it_reco->second)->get_property_float(SIDISParticle::part_theta); 
      double reco_phi = (it_reco->second)->get_property_float(SIDISParticle::part_phi); 
      double reco_E = (it_reco->second)->get_property_float(SIDISParticle::part_E);
      /* Loop over all MC particles */
      for(type_map_part::iterator it_mc = mcParticleMap.begin(); it_mc!= mcParticleMap.end(); ++it_mc){
	double mc_theta = (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_theta); 
	double mc_phi = (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_phi); 
	double mc_E = (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_E);
	/* Match the *theta* and *phi* of two particles. For details, see https://www.jlab.org/Hall-B/general/thesis/THayward_thesis.pdf */
	double dth = abs(reco_theta-mc_theta);
	double dphi = abs(reco_phi-mc_phi);
	double dE = abs(reco_E - mc_E);
	//int mcpid = (it_mc->second)->get_property_int(SIDISParticle::part_pid);
	//int recopid = (it_reco->second)->get_property_int(SIDISParticle::part_pid);
	if((dE < 0.5) &&
	   (dth < 2*degtorad) && 
	   (dphi < 2*degtorad || abs(dphi - 2*PI) < 4*degtorad)){
	  (it_mc->second)->set_property( SIDISParticle::part_pindex, (it_reco->second)->get_property_int(SIDISParticle::part_pindex));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_pid , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_pid));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_px , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_px));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_py , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_py));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_pz , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_pz));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_pt , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_pt));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_p , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_p));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_vx , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_vx));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_vy , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_vy));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_vz , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_vz));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_E , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_E));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_theta , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_theta));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_eta , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_eta));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_phi , (it_mc->second)->get_property_float(SIDISParticle::evtgen_part_phi));
	  (it_reco->second)->set_property( SIDISParticle::part_parentID , (it_mc->second)->get_property_int(SIDISParticle::evtgen_part_parentID));
	  (it_reco->second)->set_property( SIDISParticle::part_parentPID , (it_mc->second)->get_property_int(SIDISParticle::evtgen_part_parentPID));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_parentID , (it_mc->second)->get_property_int(SIDISParticle::evtgen_part_parentID));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_parentPID , (it_mc->second)->get_property_int(SIDISParticle::evtgen_part_parentPID));
	  (it_reco->second)->set_property( SIDISParticle::evtgen_part_parentparentPID , (it_mc->second)->get_property_int(SIDISParticle::evtgen_part_parentparentPID));
	}
      }
    }
    
    // Find the scattered electron
    TLorentzVector e;
    int idx_e=0;
    for (type_map_part::iterator it_reco = recoParticleMap.begin(); it_reco!= recoParticleMap.end(); ++it_reco){
      int pid = (it_reco->second)->get_property_int(SIDISParticle::part_pid);
      if(pid==11){
	double px =  (it_reco->second)->get_property_float(SIDISParticle::part_px);
	double py =  (it_reco->second)->get_property_float(SIDISParticle::part_py);
	double pz =  (it_reco->second)->get_property_float(SIDISParticle::part_pz);
	double E =  (it_reco->second)->get_property_float(SIDISParticle::part_E);
	e.SetPxPyPzE(px,py,pz,E);
	foundEle=true;
	break;
      }
      idx_e++;
    }
      
      
    // Loop over all recoParticles and fill TTree
    _nPart=0;
    _nPartMatch=0;
    for (type_map_part::iterator it_reco = recoParticleMap.begin(); it_reco!= recoParticleMap.end(); ++it_reco){
      _px[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::part_px); 
      _py[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::part_py); 
      _pz[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::part_pz); 
      _E[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::part_E); 
      _theta[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::part_theta); 
      _eta[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::part_eta); 
      _phi[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::part_phi); 
      _vz[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::part_vz); 
      _pid[_nPart] = (it_reco->second)->get_property_int(SIDISParticle::part_pid); 
      _beta[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::part_beta); 
      _chi2[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::part_chi2); 
      _pcal_sector[_nPart] = (it_reco->second)->get_property_int(SIDISParticle::cal_sector_PCAL); 
      _ecin_sector[_nPart] = (it_reco->second)->get_property_int(SIDISParticle::cal_sector_ECIN); 
      _ecout_sector[_nPart] = (it_reco->second)->get_property_int(SIDISParticle::cal_sector_ECOUT); 
      _pcal_energy[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_energy_PCAL); 
      _ecin_energy[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_energy_ECIN); 
      _ecout_energy[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_energy_ECOUT); 
      _pcal_lu[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_lu_PCAL); 
      _pcal_lv[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_lv_PCAL); 
      _pcal_lw[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_lw_PCAL); 
      _pcal_du[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_du_PCAL); 
      _pcal_dv[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_dv_PCAL); 
      _pcal_dw[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_dw_PCAL); 
      _pcal_m2u[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m2u_PCAL); 
      _pcal_m2v[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m2v_PCAL); 
      _pcal_m2w[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m2w_PCAL); 
      _pcal_m3u[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m3u_PCAL); 
      _pcal_m3v[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m3v_PCAL); 
      _pcal_m3w[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m3w_PCAL); 
      _ecin_lu[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_lu_ECIN); 
      _ecin_lv[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_lv_ECIN); 
      _ecin_lw[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_lw_ECIN); 
      _ecin_du[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_du_ECIN); 
      _ecin_dv[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_dv_ECIN); 
      _ecin_dw[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_dw_ECIN); 
      _ecin_m2u[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m2u_ECIN); 
      _ecin_m2v[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m2v_ECIN); 
      _ecin_m2w[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m2w_ECIN); 
      _ecin_m3u[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m3u_ECIN); 
      _ecin_m3v[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m3v_ECIN); 
      _ecin_m3w[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m3w_ECIN); 
      _ecout_lu[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_lu_ECOUT); 
      _ecout_lv[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_lv_ECOUT); 
      _ecout_lw[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_lw_ECOUT); 
      _ecout_du[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_du_ECOUT); 
      _ecout_dv[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_dv_ECOUT); 
      _ecout_dw[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_dw_ECOUT); 
      _ecout_m2u[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m2u_ECOUT); 
      _ecout_m2v[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m2v_ECOUT); 
      _ecout_m2w[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m2w_ECOUT); 
      _ecout_m3u[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m3u_ECOUT); 
      _ecout_m3v[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m3v_ECOUT); 
      _ecout_m3w[_nPart] = (it_reco->second)->get_property_float(SIDISParticle::cal_m3w_ECOUT); 

      _MCmatch_flag[_nPart] = (it_reco->second)->get_property_int(SIDISParticle::evtgen_part_E)>0? 1:0;       
      _MCmatch_parent_pid[_nPart] =(it_reco->second)->get_property_int(SIDISParticle::evtgen_part_parentPID);       
      _MCmatch_parent_id[_nPart] =(it_reco->second)->get_property_int(SIDISParticle::evtgen_part_parentID);       

      if(_MCmatch_flag[_nPart]==1)
	_nPartMatch++;
      _nPart++;
    }
    tree->Fill();
    _ievent++;
  }  

  cout << _ievent << " events analyzed ... " << endl;
  cout << NumNoEle << " events without an electron passing cuts ... " << endl;
  cout << NumNo2g << " events without 2 gammas passing cuts ... " << endl;

  tree->Write();
  fOut->Close();
  return 0;
}
