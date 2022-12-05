int index_of_element(std::vector<float> v, float key){
    std::vector<float>::iterator itr = std::find(v.begin(), v.end(), key);
    if (itr != v.cend()) {
        return std::distance(v.begin(), itr);
    }
    else{
        cout << "ERROR: index_of_element() cannot find key...Aborting..." << endl;
        return -1;
    }
}

int pi0_preprocess_catboost(
		   const char * input_file = "MC_3053_2.root",
		   const char * output_file = "MC_3053_2_preprocess.root"
){

  
    TFile *fOut = new TFile(output_file,"RECREATE");
    TTree *tOut = new TTree("PreProcessedEvents","PreProcessedEvents");
    
    int ievent = 0;
    int flag = 0;
    int nPhotons = 0;
    int nHadrons = 0;
    float gE = 0;
    float gTheta = 0;
    float gPhi = 0;
    float g_pcal_e = 0;
    float g_pcal_du = 0;
    float g_pcal_dv = 0;
    float g_pcal_m2u = 0;
    float g_pcal_m2v = 0;
    float g_pcal_m3u = 0;
    float g_pcal_m3v = 0;
    
    float g1R = 0; // R is defined as the angle between the photon of interest and its neighbor
                   // the symbol 'g' signifies we are looking at photon neighbors
                   // the symbol '1' signifies the closest neighbor
    float g2R = 0;
    float g1M = 0; // M is defined as the diphoton mass
    float g2M = 0;
    float g1dE = 0; // dE is defined as the energy difference between the photon of interest and its neighbor
    float g2dE = 0; 
    
    float g1_pcal_e = 0;
    float g1_pcal_du = 0;
    float g1_pcal_dv = 0;
    float g1_pcal_m2u = 0;
    float g1_pcal_m2v = 0;
    float g1_pcal_m3u = 0;
    float g1_pcal_m3v = 0;
    float g2_pcal_e = 0;
    float g2_pcal_du = 0;
    float g2_pcal_dv = 0;
    float g2_pcal_m2u = 0;
    float g2_pcal_m2v = 0;
    float g2_pcal_m3u = 0;
    float g2_pcal_m3v = 0;
    
    
    float h1R = 0; // the symbol 'h' signifies we are looking at hadron neighbors
    float h1M = 0;
    float h1dE = 0;
    float h1q = 0; // q is defined as the charge of the hadron
    float h2R = 0; // the symbol 'h' signifies we are looking at hadron neighbors
    float h2M = 0;
    float h2dE = 0;
    float h2q = 0; // q is defined as the charge of the hadron
    
    float eR = 0; // the symbol 'e' signifies we are looking at the electron 
    float eM = 0;
    float edE = 0; // edE is the fraction of energy of the electron carred by the photon
    
    tOut->Branch("ievent",&ievent,"ievent/I");
    tOut->Branch("flag",&flag,"flag/I");
    tOut->Branch("nPhotons",&nPhotons,"nPhotons/I");
    tOut->Branch("nHadrons",&nHadrons,"nHadrons/I");
    tOut->Branch("gE",&gE,"gE/F");
    tOut->Branch("gTheta",&gTheta,"gTheta/F");
    tOut->Branch("gPhi",&gPhi,"gPhi/F");
    tOut->Branch("g_pcal_e",&g_pcal_e,"g_pcal_e/F");
    tOut->Branch("g_pcal_du",&g_pcal_du,"g_pcal_du/F");
    tOut->Branch("g_pcal_dv",&g_pcal_dv,"g_pcal_dv/F");
    tOut->Branch("g_pcal_m2u",&g_pcal_m2u,"g_pcal_m2u/F");
    tOut->Branch("g_pcal_m2v",&g_pcal_m2v,"g_pcal_m2v/F");
    tOut->Branch("g_pcal_m3u",&g_pcal_m3u,"g_pcal_m3u/F");
    tOut->Branch("g_pcal_m3v",&g_pcal_m3v,"g_pcal_m3v/F");
    tOut->Branch("g1R",&g1R,"g1R/F");
    tOut->Branch("g2R",&g2R,"g2R/F");
    tOut->Branch("g1M",&g1M,"g1M/F");
    tOut->Branch("g2M",&g2M,"g2M/F");
    tOut->Branch("g1dE",&g1dE,"g1dE/F");
    tOut->Branch("g2dE",&g2dE,"g2dE/F");
    tOut->Branch("g1_pcal_e",&g1_pcal_e,"g1_pcal_e/F");
    tOut->Branch("g2_pcal_e",&g2_pcal_e,"g2_pcal_e/F");
    tOut->Branch("g1_pcal_du",&g1_pcal_du,"g1_pcal_du/F");
    tOut->Branch("g1_pcal_dv",&g1_pcal_dv,"g1_pcal_dv/F");
    tOut->Branch("g1_pcal_m2u",&g1_pcal_m2u,"g1_pcal_m2u/F");
    tOut->Branch("g1_pcal_m2v",&g1_pcal_m2v,"g1_pcal_m2v/F");
    tOut->Branch("g1_pcal_m3u",&g1_pcal_m3u,"g1_pcal_m3u/F");
    tOut->Branch("g1_pcal_m3v",&g1_pcal_m3v,"g1_pcal_m3v/F");
    tOut->Branch("g2_pcal_du",&g2_pcal_du,"g2_pcal_du/F");
    tOut->Branch("g2_pcal_dv",&g2_pcal_dv,"g2_pcal_dv/F");
    tOut->Branch("g2_pcal_m2u",&g2_pcal_m2u,"g2_pcal_m2u/F");
    tOut->Branch("g2_pcal_m2v",&g2_pcal_m2v,"g2_pcal_m2v/F");
    tOut->Branch("g2_pcal_m3u",&g2_pcal_m3u,"g2_pcal_m3u/F");
    tOut->Branch("g2_pcal_m3v",&g2_pcal_m3v,"g2_pcal_m3v/F");
    tOut->Branch("h1R",&h1R,"h1R/F");
    tOut->Branch("h2R",&h2R,"h2R/F");
    tOut->Branch("h1M",&h1M,"h1M/F");
    tOut->Branch("h2M",&h2M,"h2M/F");
    tOut->Branch("h1dE",&h1dE,"h1dE/F");
    tOut->Branch("h2dE",&h2dE,"h2dE/F");
    tOut->Branch("h1q",&h1q,"h1q/F");
    tOut->Branch("h2q",&h2q,"h2q/F");
    tOut->Branch("eR",&eR,"eR/F");
    tOut->Branch("eM",&eM,"eM/F");
    tOut->Branch("edE",&edE,"edE/F");
    
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
        int _nPart = tr_nPart[0];


        // Loop over particles in the event
        for(int i = 0 ; i < _nPart; i++){

            int pid_center = tr_pid[i];
            if(pid_center!=22) continue; // Only loop over photons

            gE = tr_E[i];
            gTheta = tr_theta[i];
            gPhi = tr_phi[i];
            g_pcal_e = tr_pcal_energy[i];
            g_pcal_du = tr_pcal_du[i];
            g_pcal_dv = tr_pcal_dv[i];
            g_pcal_m2u = tr_pcal_m2u[i];
            g_pcal_m2v = tr_pcal_m2v[i];
            g_pcal_m3u = tr_pcal_m3u[i];
            g_pcal_m3v = tr_pcal_m3v[i];
            TLorentzVector g(tr_px[i],tr_py[i],tr_pz[i],tr_E[i]);

            // Flag if the particle has MC parent
            flag = tr_MCmatch_flag[i];

            // Loop over electrons in the event
            for(int j = 0 ; j < _nPart ; j++){
                if(tr_pid[j]==11){
		    TLorentzVector e(tr_px[j],tr_py[j],tr_pz[j],tr_E[j]);
                    eR = g.Angle(e.Vect());
                    eM = (e+g).M();
                    edE = g.E()/e.E();
                }
            }
            // Loop over photons in the event
            std::vector<int> ig;
            std::vector<float> gR;
            std::vector<float> gM;
            std::vector<float> gdE;
            std::vector<float> gg_pcal_e;
            std::vector<float> gg_pcal_du;
            std::vector<float> gg_pcal_dv;
            std::vector<float> gg_pcal_m2u;
            std::vector<float> gg_pcal_m2v;
            std::vector<float> gg_pcal_m3u;
            std::vector<float> gg_pcal_m3v;
            for(int j = 0 ; j < _nPart; j++){
		if(i==j) continue;
                int pid = tr_pid[j];
                if(pid!=22) continue;
                TLorentzVector gg(tr_px[j],tr_py[j],tr_pz[j],tr_E[j]);
                ig.push_back(j);
                gR.push_back(g.Angle(gg.Vect()));
                gM.push_back((g+gg).M());
                gdE.push_back(g.E()-gg.E());
                gg_pcal_e.push_back(tr_pcal_energy[j]);
                gg_pcal_du.push_back(tr_pcal_du[j]);
                gg_pcal_dv.push_back(tr_pcal_dv[j]);
                gg_pcal_m2u.push_back(tr_pcal_m2u[j]);
                gg_pcal_m2v.push_back(tr_pcal_m2v[j]);
                gg_pcal_m3u.push_back(tr_pcal_m3u[j]);
                gg_pcal_m3v.push_back(tr_pcal_m3v[j]);
            }
            
            // Loop over hadrons in the event
            std::vector<int> ih;
            std::vector<float> hR;
            std::vector<float> hM;
            std::vector<float> hdE;
            std::vector<float> hq;
            for(int j = 0 ; j < _nPart; j++){
                if(i==j) continue;
                int pid = tr_pid[j];
                if(pid!=211 && pid!=-211 && pid!=2212 && pid!=2112 && pid!=321 && pid!=-321) continue;
                TLorentzVector h(tr_px[j],tr_py[j],tr_pz[j],tr_E[j]);
                ih.push_back(j);
                hR.push_back(g.Angle(h.Vect()));
                hM.push_back((g+h).M());
                hdE.push_back(g.E()-h.E());
                if(pid==211 || pid==2212 || pid==321)
                    hq.push_back(1);
                else if(pid==-211 || pid==-321)
                    hq.push_back(-1);
                else
                    hq.push_back(0);
            }
            
            nPhotons=ig.size()+1; // Add one because of the photon of interest (one not looped over)
            nHadrons=ih.size();
            
            // If there was no other photons or no hadrons, continue
            if(ig.size() == 0 || ih.size()==0){
	      continue;
	    }            
            
            // Sort the R vectors to find closest proximity neighbors to photon
            std::vector<float> gRclone = gR;
            std::vector<float> hRclone = hR;
            
            sort(gR.begin(),gR.end());
            sort(hR.begin(),hR.end());
            
            // Pull the nearest particles by R (if only one particle, duplicate it)
            g1R = gR.at(0);
            if(gR.size()==1) g2R = gR.at(0);
            else g2R = gR.at(1);
            
            h1R = hR.at(0);
            if(hR.size()==1) h2R = hR.at(0);
            else h2R = hR.at(1);
            
            // Get indecies of the closest particles
            int ig1 = index_of_element(gRclone,g1R);
            int ig2 = index_of_element(gRclone,g2R);
            int ih1 = index_of_element(hRclone,h1R);
            int ih2 = index_of_element(hRclone,h2R);
            
            // Pull more information from nearest particles
            g1M = gM[ig1];
            g2M = gM[ig2];
            g1dE = gdE[ig1];
            g2dE = gdE[ig2];
            g1_pcal_e = gg_pcal_e[ig1];
            g2_pcal_e = gg_pcal_e[ig2];
            g1_pcal_du = gg_pcal_du[ig1];
            g2_pcal_du = gg_pcal_du[ig2];
            g1_pcal_dv = gg_pcal_dv[ig1];
            g2_pcal_dv = gg_pcal_dv[ig2];
            g1_pcal_m2u = gg_pcal_m2u[ig1];
            g2_pcal_m2u = gg_pcal_m2u[ig2];
            g1_pcal_m2v = gg_pcal_m2v[ig1];
            g2_pcal_m2v = gg_pcal_m2v[ig2];
            g1_pcal_m3u = gg_pcal_m3u[ig1];
            g2_pcal_m3u = gg_pcal_m3u[ig2];
            g1_pcal_m3v = gg_pcal_m3v[ig1];
            g2_pcal_m3v = gg_pcal_m3v[ig2];
            
            
            h1M = hM[ih1];
            h2M = hM[ih2];
            h1q = hq[ih1];
            h2q = hq[ih2];
            h1dE = hdE[ih2];
            h2dE = hdE[ih1];
            
            // Fill TTree
            tOut->Fill();
        }
        ievent++;
    }
    // Write TTree
    tOut->Write();
    // Close TFile
    fOut->Close();
    return 0;
}
    
