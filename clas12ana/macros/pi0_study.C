#include "../Constants.h"
#include "../Kinematics.h"
#include "../SIDISParticle.h"
#include "../SIDISParticlev1.h"
#include "../PostProcess.h"
#include "../PID.h"
#include "../FiducialCuts.h"
#include "../src/HipoBankInterface.h"

typedef std::map<int, SIDISParticle*> type_map_part;

void DeleteParticlePointers(type_map_part& map){
  for(type_map_part::iterator it = map.begin(); it!=map.end(); ++it){
    delete (it->second);
    map.erase(it);
  }
  return;
}
int pi0_study(const char * hipoFile = "/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/45nA_job_3053_1.hipo",
		const char * outputFile = "MC_3053_1.root",
		const double _electron_beam_energy = 10.6,
		const int maxEvents = -1,
		bool hipo_is_mc = true){
  
  // Open TTree and declare branches
  // -------------------------------------
  TFile *fOut = new TFile(outputFile,"RECREATE");
  TTree *tree = new TTree("tree_postprocess","tree_postprocess");

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

  // Global variables
  int   _ievent=0;
  int   _nPart=0;
  int   _nPartMatch=0;
  
  // Particle kinematics
  float _px=0.;
  float _py=0.;
  float _pz=0.;
  float _E=0.;

  // Additional REC::Particle information
  int   _pid=0;
  float _beta=0.0;
  float _chi2=0.0;

  // MC matching variables
  int   _MCmatch_flag=0;
  int   _MCmatch_id=0;
  int   _MCmatch_parentpid=0;
  
  // REC::Calorimeter variables
  int   _pcal_sector=0;
  int   _ecin_sector=0;
  int   _ecout_sector=0;
  float _pcal_energy=0.;
  float _ecin_energy=0.;
  float _ecout_energy=0.;
  float _pcal_lu=0;
  float _pcal_lv=0;
  float _ecin_lu=0;
  float _ecin_lv=0;
  float _ecout_lu=0;
  float _ecout_lv=0;

  // TTree Branching
  tree->Branch("x",&reco_x,"reco_x/F");
  tree->Branch("Q2",&reco_Q2,"reco_Q2/F");
  tree->Branch("W",&reco_W,"reco_W/F");
  tree->Branch("y",&reco_y,"reco_y/F");
  tree->Branch("ievent",&_ievent,"ievent/I");
  tree->Branch("nPart",&_nPart,"nPart/I");
  tree->Branch("nPartMatch",&_nPartMatch,"nPartMatch/I");
  tree->Branch("px",&_px,"px/F");
  tree->Branch("py",&_py,"py/F");
  tree->Branch("pz",&_pz,"pz/F");
  tree->Branch("E",&_E,"E/F");
  tree->Branch("pid",&_pid,"pid/I");
  tree->Branch("beta",&_beta,"beta/F");
  tree->Branch("chi2",&_chi2,"chi2/F");
  tree->Branch("MCmatch_flag",&_MCmatch_flag,"MCmatch_flag/I");
  tree->Branch("MCmatch_id",&_MCmatch_id,"MCmatch_id/I");
  tree->Branch("MCmatch_parentpid",&_MCmatch_parentpid,"MCmatch_parentpid/I");
  tree->Branch("pcal_sector",&_pcal_sector,"pcal_sector/I");
  tree->Branch("ecin_sector",&_ecin_sector,"ecin_sector/I");
  tree->Branch("ecout_sector",&_ecout_sector,"ecout_sector/I");
  tree->Branch("pcal_energy",&_pcal_energy,"pcal_energy/F");
  tree->Branch("ecin_energy",&_ecin_energy,"ecin_energy/F");
  tree->Branch("ecout_energy",&_ecout_energy,"ecout_energy/F");
  tree->Branch("pcal_lu",&_pcal_lu,"pcal_lu/F");
  tree->Branch("ecin_lu",&_ecin_lu,"ecin_lu/F");
  tree->Branch("ecout_lu",&_ecout_lu,"ecout_lu/F");
  tree->Branch("pcal_lv",&_pcal_lv,"pcal_lv/F");
  tree->Branch("ecin_lv",&_ecin_lv,"ecin_lv/F");
  tree->Branch("ecout_lv",&_ecout_lv,"ecout_lv/F");

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
  int NumNoGamma=0;

  while(_chain.Next()==true && (_ievent < maxEvents || maxEvents < 0)){
    if(_ievent%10000==0 && _ievent!=0){
      std::cout << _ievent << " events completed " << std::endl;
    }
    
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
    for(unsigned int idx = 0 ; idx < particles.size() ; idx++){
      // Extract each particle from event one-at-a-time
      // -------------------------------------------------------
      auto particle = particles.at(idx);
      int pid = particle->getPid();      
      float chi2 = particle->getChi2Pid();
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

      // CUT REC::Particle
      // --------------------------------------------------------------------------
      if(_pidhelper.performPIDCuts(sp)==false)
	{delete sp;  continue;}

      // CUT Fiducial
      // --------------------------------------------------------------------------
      if(_fiducial.FidCutParticle(_c12,11,sp) == false)
	{delete sp; continue;}
      
    
      // Add SIDISParticle to the collection
      // --------------------------------------------------------------------------
      if(pid==11){
	foundEle=true;
      }
      recoParticleMap.insert ( make_pair( l++ , sp) );

    }

    if(foundEle==false){
      NumNoEle++;
      continue;
    }
    
    // Loop over all Monte Carlo particles
    // -------------------------------------
    auto mcparticles=_c12->mcparts();
    for(int idx = 0 ; idx < mcparticles->getRows() && hipo_is_MC ; idx++){
    
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
    for (type_map_part::iterator it_reco = recoParticleMap.begin(); it_reco!= recoParticleMap.end() && hipo_is_MC; ++it_reco){

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
    for(unsigned int i = 0; i < recoParticleMap.size() ; i++){
      int pid = recoParticleMap[i]->get_property_int(SIDISParticle::part_pid);
      if(pid==11){
	double px = recoParticleMap[i]->get_property_float(SIDISParticle::part_px);
	double py = recoParticleMap[i]->get_property_float(SIDISParticle::part_py);
	double pz = recoParticleMap[i]->get_property_float(SIDISParticle::part_pz);
	double E = recoParticleMap[i]->get_property_float(SIDISParticle::part_E);
	e.SetPxPyPzE(px,py,pz,E);
	foundEle=true;
	idx_e=i;
	break;
      }
    }
      
    // Loop over all recoParticles and fill TTree
    _nPart=0;
    _nPartMatch=0;
    for (type_map_part::iterator it_reco = recoParticleMap.begin(); it_reco!= recoParticleMap.end() && hipo_is_MC; ++it_reco){
      _px = (it_reco->second)->get_property_float(SIDISParticle::part_px); 
      _py = (it_reco->second)->get_property_float(SIDISParticle::part_py); 
      _pz = (it_reco->second)->get_property_float(SIDISParticle::part_pz); 
      _E = (it_reco->second)->get_property_float(SIDISParticle::part_E); 
      _pid = (it_reco->second)->get_property_int(SIDISParticle::part_pid); 
      _beta = (it_reco->second)->get_property_float(SIDISParticle::part_beta); 
      _chi2 = (it_reco->second)->get_property_float(SIDISParticle::part_chi2); 
      _pcal_sector = (it_reco->second)->get_property_int(SIDISParticle::cal_sector_PCAL); 
      _ecin_sector = (it_reco->second)->get_property_int(SIDISParticle::cal_sector_ECIN); 
      _ecout_sector = (it_reco->second)->get_property_int(SIDISParticle::cal_sector_ECOUT); 
      _pcal_energy = (it_reco->second)->get_property_float(SIDISParticle::cal_energy_PCAL); 
      _ecin_energy = (it_reco->second)->get_property_float(SIDISParticle::cal_energy_ECIN); 
      _ecout_energy = (it_reco->second)->get_property_float(SIDISParticle::cal_energy_ECOUT); 
      _pcal_lu = (it_reco->second)->get_property_float(SIDISParticle::cal_lu_PCAL); 
      _ecin_lu = (it_reco->second)->get_property_float(SIDISParticle::cal_lu_ECIN); 
      _ecout_lu = (it_reco->second)->get_property_float(SIDISParticle::cal_lu_ECOUT); 
      _pcal_lv = (it_reco->second)->get_property_float(SIDISParticle::cal_lv_PCAL); 
      _ecin_lv = (it_reco->second)->get_property_float(SIDISParticle::cal_lv_ECIN); 
      _ecout_lv = (it_reco->second)->get_property_float(SIDISParticle::cal_lv_ECOUT); 

      _MCmatch_flag = (it_reco->second)->get_property_int(SIDISParticle::evtgen_part_E)>0? 1:0;       
      _MCmatch_parent_pid =(it_reco->second)->get_property_int(SIDISParticle::evtgen_part_parentPID);       
      _MCmatch_parent_id =(it_reco->second)->get_property_int(SIDISParticle::evtgen_part_parentID);       

      _nPart++;
      if(_MCmatch_flag==1)
	_nPartMatch++;
    }
    _ievent++;
  }  
  cout << _ievent << " events analyzed ... " << endl;
  cout << NumNoEle << " events without an electron passing cuts ... " << endl;

  tree->Write();
  fOut->Close();
  return 0;
}
