#include "BSM_Analysis.h"

int main (int argc, char *argv[])
{
  
	//TApplication app("App",&argc, argv);
	TFile * MassHisto = new TFile(argv[2], "RECREATE");
	int nDir = 5;
	TDirectory *theDirectory[nDir];
	theDirectory[0] = MassHisto->mkdir("AfterDimuonSelection");
	theDirectory[1] = MassHisto->mkdir("AfterPhotonSelection");
	theDirectory[2] = MassHisto->mkdir("AfterDileptonMassSelection");
	theDirectory[3] = MassHisto->mkdir("AfterRemoveJPsiContamination");
    theDirectory[4] = MassHisto->mkdir("AfterRemoveUpsilonContamination");
	cout <<argv[1]<<endl;
	BSM_Analysis BSM_Analysis_(MassHisto, theDirectory, nDir, argv[1]);
}

BSM_Analysis::BSM_Analysis(TFile* theFile, TDirectory *cdDir[], int nDir, char* fname)
{
	crateHistoMasps(nDir);
	//load PU weights
	TFile file_PUdata("PUdata.root","read");
	TH1F *PUweights = (TH1F*)file_PUdata.Get("analyzeHiMassTau/NVertices_0");
	PUweights->Scale(1/PUweights->Integral());
	TFile file_PUsim("PUsim.root","read");
	TH1F *PUsim = (TH1F*)file_PUsim.Get("analyzeHiMassTau/NVertices_0");
	PUsim->Scale(1/PUsim->Integral());
  
	PUweights->Divide(PUsim);
  
	//configure input file
	TFile *f = TFile::Open(fname);
	f->cd("TNT");
	TTree* BOOM = (TTree*)f->Get("TNT/BOOM");
  
	int nentries = (int) BOOM->GetEntries();
	setBranchAddress(BOOM);
  
	for (int i = 0; i < 10000; ++i)
		//for (int i = 0; i < nentries; ++i)
	{
		BOOM->GetEntry(i);
      
		//define global event weight
		double weight =1.;
		weight=PUweights->GetBinContent(PUweights->FindBin(trueInteractions));
		// TLorentz vector -------------------------------
		TLorentzVector first_muon_vec(0., 0., 0., 0.);
		TLorentzVector Subfirst_muon_vec(0., 0., 0., 0.);
		TLorentzVector dalitz_first_muon_vec(0., 0., 0., 0.);
		TLorentzVector dalitz_Subfirst_muon_vec(0., 0., 0., 0.);
		TLorentzVector Photon0_TL(0., 0., 0., 0.);
		TLorentzVector Photon_medium_TL(0., 0., 0., 0.);
		TLorentzVector Photon_tight_TL(0., 0., 0., 0.);
		TLorentzVector Dalitz_photonless_vec(0., 0., 0., 0.);
		TLorentzVector Dalitz_photonhigh_vec(0., 0., 0., 0.);
		double charge_lead = 0.;
		double charge_slead = 0.;
		int lmuon_counter = 0;
		int smuon_counter = 0;
		int pass_dalitz_id[nDir] = {0};
      
		bool pass_trigger = false;
      
		// For Trigger================
		for (int t = 0 ; t < Trigger_decision->size(); t++){
			string theTriggers = Trigger_names->at(t);
         
			//cout <<theTriggers<<endl;
			string myTrigger   = "HLT_Mu17_Photon30";
			std::size_t found = theTriggers.find(myTrigger);
			if (found!=std::string::npos){
				if (Trigger_decision->at(t) == 1){pass_trigger = true;}
			}
		}
      
		float first_muon_pt = 99999.;
		float dimuon_mass_int = 999999.;
		bool match = false;
		bool umatch = false;
     
		pass_trigger = true; 
		// Select dimuons ================
		for (int m1 = 0; m1 < Muon_pt->size(); m1++)
		{
			if ((Muon_pt->size() > 1) && (abs(Muon_eta->at(m1) <2.4)) && (Muon_pt->at(m1) > 23.0) && (pass_trigger)){
				for (int m2 = 0; m2 < Muon_pt->size(); m2++)
				{
					if ((abs(Muon_eta->at(m2) <2.4)) && (Muon_pt->at(m2) >4.0)){
		  
						if (Muon_charge->at(m1)*Muon_charge->at(m2) < 0){
							first_muon_vec.SetPtEtaPhiE(Muon_pt->at(m1), Muon_eta->at(m1), Muon_phi->at(m1), Muon_energy->at(m1));
							Subfirst_muon_vec.SetPtEtaPhiE(Muon_pt->at(m2), Muon_eta->at(m2), Muon_phi->at(m2), Muon_energy->at(m2));
							float dimuon_mass = (first_muon_vec + Subfirst_muon_vec).M();
		    
							if (dimuon_mass < dimuon_mass_int ){
								first_muon_vec.SetPtEtaPhiE(Muon_pt->at(m1), Muon_eta->at(m1), Muon_phi->at(m1), Muon_energy->at(m1));
								Subfirst_muon_vec.SetPtEtaPhiE(Muon_pt->at(m2), Muon_eta->at(m2), Muon_phi->at(m2), Muon_energy->at(m2));
								dimuon_mass_int = dimuon_mass;
							}
						}
		  
					}
					pass_dalitz_id[0] = 1;
				}
			}
		}
      
		//Photon selection =====================
      
		vector<TLorentzVector> Photon0_vec;
		vector<TLorentzVector> Photon_medium_vec;
		vector<TLorentzVector> Photon_tight_vec;
      
		Photon0_vec.erase( Photon0_vec.begin(), Photon0_vec.end() ) ;
		Photon_medium_vec.erase( Photon_medium_vec.begin(), Photon_medium_vec.end() ) ;
		Photon_tight_vec.erase( Photon_tight_vec.begin(), Photon_tight_vec.end() ) ;
      
		for (int ph = 0; ph < Photon_pt->size(); ph++)
		{
			if ((Photon_pt->size() > 1) && (Photon_pt->at(ph) > 15.0) && (pass_trigger)){
				Photon0_TL.SetPtEtaPhiE(Photon_pt->at(ph), Photon_eta->at(ph), Photon_phi->at(ph), Photon_energy->at(ph));
				Photon0_vec.push_back(Photon0_TL);
				if ((Photon_EleVeto->at(ph) == 1) && (Photon_HoverE->at(ph) <0.05)) {
					if(Photon_eta->at(ph) < 1.44){
						if((Photon_SigmaIEtaIEta->at(ph) < 0.011) && (Photon_PFChIso->at(ph) < 1.5) && (Photon_PFPhoIso->at(ph) < 1.0) && (Photon_PFNeuIso->at(ph) < 0.7)){
							Photon_medium_TL.SetPtEtaPhiE(Photon_pt->at(ph), Photon_eta->at(ph), Photon_phi->at(ph), Photon_energy->at(ph));
							Photon_medium_vec.push_back(Photon_medium_TL);
						}
						if((Photon_SigmaIEtaIEta->at(ph) < 0.011) && (Photon_PFChIso->at(ph) < 0.7) && (Photon_PFPhoIso->at(ph) < 0.4) && (Photon_PFNeuIso->at(ph) < 0.5)){
							Photon_tight_TL.SetPtEtaPhiE(Photon_pt->at(ph), Photon_eta->at(ph), Photon_phi->at(ph), Photon_energy->at(ph));
							Photon_tight_vec.push_back(Photon_tight_TL);
						}
					}
					else {
						if((Photon_SigmaIEtaIEta->at(ph) < 0.033) && (Photon_PFChIso->at(ph) < 1.2) && (Photon_PFPhoIso->at(ph) < 1.5) && (Photon_PFNeuIso->at(ph) < 1.0)){
							Photon_medium_TL.SetPtEtaPhiE(Photon_pt->at(ph), Photon_eta->at(ph), Photon_phi->at(ph), Photon_energy->at(ph));
							Photon_medium_vec.push_back(Photon_medium_TL);
						}
						if((Photon_SigmaIEtaIEta->at(ph) < 0.031) && (Photon_PFChIso->at(ph) < 0.5) && (Photon_PFPhoIso->at(ph) < 1.5) && (Photon_PFNeuIso->at(ph) < 1.0)){
							Photon_tight_TL.SetPtEtaPhiE(Photon_pt->at(ph), Photon_eta->at(ph), Photon_phi->at(ph), Photon_energy->at(ph));
							Photon_tight_vec.push_back(Photon_tight_TL);
						}
					}
				}
			}
		}
      
		//three body mass --------------------------
      
		double m_threebody;
		double  min_photon_pt = 0.0;
		double max_photon_pt = 0.0; 
      
		for(int tb = 0; tb < Photon_medium_vec.size(); tb++)
		{
			if(Photon_medium_vec.size() > 0){
				if (Photon_medium_vec.size() == 1){
					Dalitz_photonless_vec = Photon_medium_vec.at(tb);
				}else {
					if(Photon_medium_vec.at(tb).Pt() > max_photon_pt){
						Dalitz_photonhigh_vec = Photon_medium_vec.at(tb);
						if (tb == 0) {
							Dalitz_photonless_vec = Photon_medium_vec.at(tb);
							min_photon_pt = Photon_medium_vec.at(tb).Pt();
						}
						max_photon_pt = Photon_medium_vec.at(tb).Pt(); 
					} else if (Photon_medium_vec.at(tb).Pt() < min_photon_pt ) {
						Dalitz_photonless_vec = Photon_medium_vec.at(tb);
						min_photon_pt = Photon_medium_vec.at(tb).Pt();
					}
				}
	  
				pass_dalitz_id[1] = 1;	  
			}
		}
      
		// Dilepton mass selection (20GeV)--------------------------------
	
		for(int tb = 0; (tb < Photon_medium_vec.size()) && (tb < Muon_pt->size()); tb++)
		{
			if((Photon_medium_vec.size() > 0) &&(Muon_pt->size() > 1)) {
		
				float dimuon_mass = (first_muon_vec + Subfirst_muon_vec).M();
				if (dimuon_mass < 20.0) // dalitz dilepton mass
				{
					//dalitz_first_muon_vec = first_muon_vec;
					//dalitz_first_muon_vec = Subfirst_muon_vec;
					pass_dalitz_id[2] = 1;
					
					if((dimuon_mass < 2.9) || (dimuon_mass > 3.3)) // remove Jpsi contamination
					{
						//dalitz_first_muon_vec = first_muon_vec;
						//dalitz_first_muon_vec = Subfirst_muon_vec;
						pass_dalitz_id[3] = 1;
						
						if((dimuon_mass < 9.3) || (dimuon_mass > 9.7)) //remove Upsilon contamination
						{
							dalitz_first_muon_vec = first_muon_vec;
							dalitz_first_muon_vec = Subfirst_muon_vec;
							pass_dalitz_id[4] = 1;
						}
					}
				}
				
			}
		}

		
	
	
	
      
  
		for (int i = 0; i < nDir; i++)
		{
			if (pass_dalitz_id[i] == 1){
	    
				_hmap_diMuon_mass[i]->Fill((first_muon_vec + Subfirst_muon_vec).M());
				_hmap_lead_muon_pT[i]->Fill(first_muon_vec.Pt());
				_hmap_slead_muon_pT[i]->Fill(Subfirst_muon_vec.Pt());
				_hmap_lead_muon_eta[i]->Fill(first_muon_vec.Eta());
				_hmap_slead_muon_eta[i]->Fill(Subfirst_muon_vec.Eta());
				_hmap_lead_muon_phi[i]->Fill(first_muon_vec.Phi());
				_hmap_slead_muon_phi[i]->Fill(Subfirst_muon_vec.Phi());
	    
			}
			if ((i > 0) && (pass_dalitz_id[i] == 1)){
				_hmap_threebody_light[i]->Fill((first_muon_vec + Subfirst_muon_vec + Dalitz_photonless_vec).M());
				_hmap_threebody_heavy[i]->Fill((first_muon_vec + Subfirst_muon_vec + Dalitz_photonhigh_vec).M());
			}
			
	
	  
		}
      
	}
	theFile->cd();
	for (int d = 0; d < nDir; d++)
	{
		cdDir[d]->cd();
		_hmap_diMuon_mass[d]->Write();
		_hmap_lead_muon_pT[d]->Write();
		_hmap_slead_muon_pT[d]->Write();
		_hmap_lead_muon_eta[d]->Write();
		_hmap_slead_muon_eta[d]->Write();
		_hmap_lead_muon_phi[d]->Write();
		_hmap_slead_muon_phi[d]->Write();
      
		if(d > 0){
			_hmap_threebody_light[d]->Write();
			_hmap_threebody_heavy[d]->Write();
		}
	  
      
	}
	theFile->Close();
  
}

BSM_Analysis::~BSM_Analysis()
{
	// do anything here that needs to be done at desctruction time
}


void BSM_Analysis::crateHistoMasps (int directories)
{
	for (int i = 0; i < directories; i++)
	{
		// Muon distributions
		_hmap_diMuon_mass[i]      = new TH1F("diMuonMass",      "m_{#mu, #mu}", 300., 0., 300.);
		_hmap_lead_muon_pT[i]     = new TH1F("lead_muon_pT",    "#mu p_{T}",    300, 0., 300.);
		_hmap_slead_muon_pT[i]    = new TH1F("slead_muon_pT",   "#mu p_{T}",    300, 0., 300.);
		_hmap_lead_muon_eta[i]    = new TH1F("lead_muon_eta",   "#mu #eta",     50, -2.5, 2.5);
		_hmap_slead_muon_eta[i]   = new TH1F("slead_muon_eta",  "#mu #eta",     50, -2.5, 2.5);
		_hmap_lead_muon_phi[i]    = new TH1F("lead_muon_phi",   "#mu #phi",     70, -3.5, 3.5);
		_hmap_slead_muon_phi[i]    = new TH1F("slead_muon_phi", "#mu #phi",     70, -3.5, 3.5);
      
		if(i > 0){
			_hmap_threebody_light[i]  = new TH1F("threebody_light_mass", "m_{#mu,#mu,#gamma}",200.,50.,250.);
			_hmap_threebody_heavy[i]  = new TH1F("threebody_heavy_mass", "m_{#mu,#mu,#gamma}",200.,50.,250.);
		}
		
		if(i>1){
			_hmap_diMuon_mass[i]      = new TH1F("diMuonMass",      "m_{#mu, #mu}", 100., 0., 50.);
		}
	
      
	}
}

void BSM_Analysis::setBranchAddress(TTree* BOOM)
{
  
	// Set object pointer
  
	Muon_pt = 0;
	Muon_eta = 0;
	Muon_phi = 0;
	Muon_p = 0;
	Muon_energy = 0;
	Muon_charge = 0;
	Muon_tight = 0;
	Muon_soft = 0;
	Muon_pf = 0;
	Muon_isoCharged = 0;
	Muon_isoSum = 0;
	Muon_isoCharParPt = 0;
	Muon_chi2 = 0;
	Muon_validHits = 0;
	Muon_validHitsInner = 0;
	Muon_matchedStat = 0;
	Muon_dxy_pv = 0;
	//Muon_dxy = 0;
	Muon_TLayers = 0;
	//Muon_dz = 0;/
	Muon_dz_bs = 0;
	Muon_isoNeutralHadron = 0;
	Muon_isoPhoton = 0;
	Muon_isoPU = 0;
	Photon_eta = 0;
	Photon_phi = 0;
	Photon_pt = 0;
	Photon_energy = 0;
	Photon_et = 0;
	Photon_HoverE = 0;
	Photon_phoR9 = 0;
	Photon_SigmaIEtaIEta = 0;
	Photon_SigmaIPhiIPhi = 0;
	Photon_PFChIso = 0;
	Photon_PFPhoIso = 0;
	Photon_PFNeuIso = 0;
	Photon_EleVeto = 0;
	Jet_pt = 0;
	Jet_eta = 0;
	Jet_phi = 0;
	Jet_energy = 0;
	Jet_bDiscriminator = 0;
	Jet_mass = 0;
	Jet_neutralHadEnergyFraction = 0;
	Jet_neutralEmEmEnergyFraction = 0;
	Jet_chargedHadronEnergyFraction = 0;
	Jet_chargedEmEnergyFraction = 0;
	Jet_muonEnergyFraction = 0;
	Jet_electronEnergy = 0;
	Jet_photonEnergy = 0;
	UncorrJet_pt = 0;
	Trigger_names = 0;
	Trigger_decision = 0;
	// Set branch addresses and branch pointers
	if(!BOOM) return;
	BOOM->SetBranchAddress("Trigger_names", &Trigger_names, &b_Trigger_names);
	BOOM->SetBranchAddress("Trigger_decision", &Trigger_decision, &b_Trigger_decision);
	BOOM->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
	BOOM->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
	BOOM->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
	BOOM->SetBranchAddress("Muon_p", &Muon_p, &b_Muon_p);
	BOOM->SetBranchAddress("Muon_energy", &Muon_energy, &b_Muon_energy);
	BOOM->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
	BOOM->SetBranchAddress("Muon_tight", &Muon_tight, &b_Muon_tight);
	BOOM->SetBranchAddress("Muon_soft", &Muon_soft, &b_Muon_soft);
	BOOM->SetBranchAddress("Muon_pf", &Muon_pf, &b_Muon_pf);
	BOOM->SetBranchAddress("Muon_isoCharged", &Muon_isoCharged, &b_Muon_isoCharged);
	BOOM->SetBranchAddress("Muon_isoSum", &Muon_isoSum, &b_Muon_isoSum);
	BOOM->SetBranchAddress("Muon_isoCharParPt", &Muon_isoCharParPt, &b_Muon_isoCharParPt);
	BOOM->SetBranchAddress("Muon_chi2", &Muon_chi2, &b_Muon_chi2);
	BOOM->SetBranchAddress("Muon_validHits", &Muon_validHits, &b_Muon_validHits);
	BOOM->SetBranchAddress("Muon_validHitsInner", &Muon_validHitsInner, &b_Muon_validHitsInner);
	BOOM->SetBranchAddress("Muon_matchedStat", &Muon_matchedStat, &b_Muon_matchedStat);
	// BOOM->SetBranchAddress("Muon_dxy", &Muon_dxy, &b_Muon_dxy); 
	BOOM->SetBranchAddress("Muon_dxy_pv", &Muon_dxy_pv, &b_Muon_dxy_pv);
	BOOM->SetBranchAddress("Muon_TLayers", &Muon_TLayers, &b_Muon_TLayers);
	// BOOM->SetBranchAddress("Muon_dz", &Muon_dz, &b_Muon_dz);
	BOOM->SetBranchAddress("Muon_dz_bs", &Muon_dz_bs, &b_Muon_dz_bs);
	BOOM->SetBranchAddress("Muon_isoNeutralHadron", &Muon_isoNeutralHadron, &b_Muon_isoNeutralHadron);
	BOOM->SetBranchAddress("Muon_isoPhoton", &Muon_isoPhoton, &b_Muon_isoPhoton);
	BOOM->SetBranchAddress("Muon_isoPU", &Muon_isoPU, &b_Muon_isoPU);
	BOOM->SetBranchAddress("Photon_eta", &Photon_eta, &b_Photon_eta);
	BOOM->SetBranchAddress("Photon_phi", &Photon_phi, &b_Photon_phi);
	BOOM->SetBranchAddress("Photon_pt", &Photon_pt, &b_Photon_pt);
	BOOM->SetBranchAddress("Photon_energy", &Photon_energy, &b_Photon_energy);
	BOOM->SetBranchAddress("Photon_et", &Photon_et, &b_Photon_et);
	BOOM->SetBranchAddress("Photon_HoverE", &Photon_HoverE, &b_Photon_HoverE);
	BOOM->SetBranchAddress("Photon_phoR9", &Photon_phoR9, &b_Photon_phoR9);
	BOOM->SetBranchAddress("Photon_SigmaIEtaIEta", &Photon_SigmaIEtaIEta, &b_Photon_SigmaIEtaIEta);
	BOOM->SetBranchAddress("Photon_SigmaIPhiIPhi", &Photon_SigmaIPhiIPhi, &b_Photon_SigmaIPhiIPhi);
	BOOM->SetBranchAddress("Photon_PFChIso", &Photon_PFChIso, &b_Photon_PFChIso);
	BOOM->SetBranchAddress("Photon_PFPhoIso", &Photon_PFPhoIso, &b_Photon_PFPhoIso);
	BOOM->SetBranchAddress("Photon_PFNeuIso", &Photon_PFNeuIso, &b_Photon_PFNeuIso);
	BOOM->SetBranchAddress("Photon_EleVeto", &Photon_EleVeto, &b_Photon_EleVeto);
	BOOM->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
	BOOM->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
	BOOM->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
	BOOM->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
	BOOM->SetBranchAddress("Jet_bDiscriminator", &Jet_bDiscriminator, &b_Jet_bDiscriminator);
	BOOM->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
	BOOM->SetBranchAddress("Jet_neutralHadEnergyFraction", &Jet_neutralHadEnergyFraction, &b_Jet_neutralHadEnergyFraction);
	BOOM->SetBranchAddress("Jet_neutralEmEmEnergyFraction", &Jet_neutralEmEmEnergyFraction, &b_Jet_neutralEmEmEnergyFraction);
	BOOM->SetBranchAddress("Jet_chargedHadronEnergyFraction", &Jet_chargedHadronEnergyFraction, &b_Jet_chargedHadronEnergyFraction);
	BOOM->SetBranchAddress("Jet_chargedEmEnergyFraction", &Jet_chargedEmEnergyFraction, &b_Jet_chargedEmEnergyFraction);
	BOOM->SetBranchAddress("Jet_muonEnergyFraction", &Jet_muonEnergyFraction, &b_Jet_muonEnergyFraction);
	BOOM->SetBranchAddress("Jet_electronEnergy", &Jet_electronEnergy, &b_Jet_electronEnergy);
	BOOM->SetBranchAddress("Jet_photonEnergy", &Jet_photonEnergy, &b_Jet_photonEnergy);
	BOOM->SetBranchAddress("UncorrJet_pt", &UncorrJet_pt, &b_UncorrJet_pt);
	BOOM->SetBranchAddress("npuVertices", &npuVertices, &b_npuVertices);
	BOOM->SetBranchAddress("trueInteractions", &trueInteractions, &b_trueInteractions);
	BOOM->SetBranchAddress("ootnpuVertices", &ootnpuVertices, &b_ootnpuVertices);
	BOOM->SetBranchAddress("npuVerticesp1", &npuVerticesp1, &b_npuVerticesp1);
	BOOM->SetBranchAddress("bestVertices", &bestVertices, &b_bestVertices);
	//BOOM->SetBranchAddress("Met_pt", &Met_pt, &b_Met_pt);
	//BOOM->SetBranchAddress("Met_sumEt", &Met_sumEt, &b_Met_sumEt);
	//BOOM->SetBranchAddress("Met_phi", &Met_phi, &b_Met_phi);
	//BOOM->SetBranchAddress("Met_px", &Met_px, &b_Met_px);
	//BOOM->SetBranchAddress("Met_py", &Met_py, &b_Met_py);
	//BOOM->SetBranchAddress("Met_pz", &Met_pz, &b_Met_pz);
	BOOM->SetBranchAddress("Met_puppi_pt", &Met_puppi_pt, &b_Met_puppi_pt);
	BOOM->SetBranchAddress("Met_puppi_sumEt", &Met_puppi_sumEt, &b_Met_puppi_sumEt);
	BOOM->SetBranchAddress("Met_puppi_phi", &Met_puppi_phi, &b_Met_puppi_phi);
	BOOM->SetBranchAddress("Met_puppi_px", &Met_puppi_px, &b_Met_puppi_px);
	BOOM->SetBranchAddress("Met_puppi_py", &Met_puppi_py, &b_Met_puppi_py);
	BOOM->SetBranchAddress("Met_puppi_pz", &Met_puppi_pz, &b_Met_puppi_pz);
	BOOM->SetBranchAddress("Gen_Met", &Gen_Met, &b_Gen_Met);
	BOOM->SetBranchAddress("Met_type1PF_shiftedPtUp", &Met_type1PF_shiftedPtUp, &b_Met_type1PF_shiftedPtUp);
	BOOM->SetBranchAddress("Met_type1PF_shiftedPtDown", &Met_type1PF_shiftedPtDown, &b_Met_type1PF_shiftedPtDown);
	//BOOM->SetBranchAddress("Met_shiftedPtUp", &Met_shiftedPtUp, &b_Met_shiftedPtUp);
	//BOOM->SetBranchAddress("Met_shiftedPtDown", &Met_shiftedPtDown, &b_Met_shiftedPtDown);
};
