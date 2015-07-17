////////////////////////////////////////////////////////////////////
// Developer: Andres Florez, Universidad de los Andes, Colombia. //
// Developer: Denis Rathjens, University of Hamburg, Germany     //
//////////////////////////////////////////////////////////////////


#ifndef BSM_Analysis_h
#define BSM_Analysis_h

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TBranch.h>
#include <TApplication.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TDirectory.h>
using namespace std;

class BSM_Analysis {
public :
   BSM_Analysis(TFile*, TDirectory* dir[], int nDir, char*);
   ~BSM_Analysis();

   // create Histo maps
   void crateHistoMasps (int);
   // Define maps for histograms
   // For muons
   std::map<unsigned int, TH1*> _hmap_lead_muon_pT;
   std::map<unsigned int, TH1*> _hmap_slead_muon_pT;
   std::map<unsigned int, TH1*> _hmap_lead_muon_eta;
   std::map<unsigned int, TH1*> _hmap_slead_muon_eta;
   std::map<unsigned int, TH1*> _hmap_lead_muon_phi;
   std::map<unsigned int, TH1*> _hmap_slead_muon_phi;
   std::map<unsigned int, TH1*> _hmap_diMuon_mass;
   std::map<unsigned int, TH2*> _h2map_lead_slead_muon_pT;
   // For Photons
   std::map<unsigned int, TH1*> _hmap_threebody_light;
   std::map<unsigned int, TH1*> _hmap_photon_pT;
   std::map<unsigned int, TH1*> _hmap_photon_eta;
   std::map<unsigned int, TH1*> _hmap_photon_phi;
   // Define Branches
   void setBranchAddress(TTree* BOOM);
   vector<string>  *Trigger_names;
   vector<int>     *Trigger_decision;
   vector<double>  *Muon_pt;
   vector<double>  *Muon_eta;
   vector<double>  *Muon_phi;
   vector<double>  *Muon_p;
   vector<double>  *Muon_energy;
   vector<double>  *Muon_charge;
   vector<bool>    *Muon_tight;
   vector<bool>    *Muon_soft;
   vector<bool>    *Muon_pf;
   vector<double>  *Muon_isoCharged;
   vector<double>  *Muon_isoSum;
   vector<double>  *Muon_isoCharParPt;
   vector<double>  *Muon_chi2;
   vector<double>  *Muon_validHits;
   vector<double>  *Muon_validHitsInner;
   vector<double>  *Muon_matchedStat;
   vector<double>  *Muon_dxy;
   vector<double>  *Muon_TLayers;
   vector<double>  *Muon_dz;
   vector<double>  *Muon_isoNeutralHadron;
   vector<double>  *Muon_isoPhoton;
   vector<double>  *Muon_isoPU;
   vector<float>   *Photon_eta;
   vector<float>   *Photon_phi;
   vector<float>   *Photon_pt;
   vector<float>   *Photon_energy;
   vector<float>   *Photon_et;
   vector<float>   *Photon_HoverE;
   vector<float>   *Photon_phoR9;
   vector<float>   *Photon_SigmaIEtaIEta;
   vector<float>   *Photon_SigmaIPhiIPhi;
   vector<float>   *Photon_PFChIso;
   vector<float>   *Photon_PFPhoIso;
   vector<float>   *Photon_PFNeuIso;
   vector<int>     *Photon_EleVeto;
   vector<double>  *Jet_pt;
   vector<double>  *Jet_eta;
   vector<double>  *Jet_phi;
   vector<double>  *Jet_energy;
   vector<double>  *Jet_bDiscriminator;
   vector<double>  *Jet_mass;
   vector<double>  *Jet_neutralHadEnergyFraction;
   vector<double>  *Jet_neutralEmEmEnergyFraction;
   vector<double>  *Jet_chargedHadronEnergyFraction;
   vector<double>  *Jet_chargedEmEnergyFraction;
   vector<double>  *Jet_muonEnergyFraction;
   vector<double>  *Jet_electronEnergy;
   vector<double>  *Jet_photonEnergy;
   vector<double>  *UncorrJet_pt;
   Int_t           npuVertices;
   Float_t         trueInteractions;
   Int_t           ootnpuVertices;
   Int_t           npuVerticesp1;
   Int_t           bestVertices;
   Double_t        Met_pt;
   Double_t        Met_sumEt;
   Double_t        Met_phi;
   Double_t        Met_px;
   Double_t        Met_py;
   Double_t        Met_pz;
   Double_t        Gen_Met;
   Double_t        Met_shiftedPtUp;
   Double_t        Met_shiftedPtDown;

   // List of branches
   TBranch        *b_Trigger_names;
   TBranch        *b_Trigger_decision;
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_p;   //!
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_tight;   //!
   TBranch        *b_Muon_soft;   //!
   TBranch        *b_Muon_pf;   //!
   TBranch        *b_Muon_isoCharged;   //!
   TBranch        *b_Muon_isoSum;   //!
   TBranch        *b_Muon_isoCharParPt;   //!
   TBranch        *b_Muon_chi2;   //!
   TBranch        *b_Muon_validHits;   //!
   TBranch        *b_Muon_validHitsInner;   //!
   TBranch        *b_Muon_matchedStat;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_TLayers;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_isoNeutralHadron;   //!
   TBranch        *b_Muon_isoPhoton;   //!
   TBranch        *b_Muon_isoPU;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_energy;   //!
   TBranch        *b_Photon_et;
   TBranch        *b_Photon_HoverE;
   TBranch        *b_Photon_phoR9;
   TBranch        *b_Photon_SigmaIEtaIEta;
   TBranch        *b_Photon_SigmaIPhiIPhi;
   TBranch        *b_Photon_PFChIso;
   TBranch        *b_Photon_PFPhoIso;
   TBranch        *b_Photon_PFNeuIso;
   TBranch        *b_Photon_EleVeto;
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_bDiscriminator;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_neutralHadEnergyFraction;   //!
   TBranch        *b_Jet_neutralEmEmEnergyFraction;   //!
   TBranch        *b_Jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jet_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_muonEnergyFraction;   //!
   TBranch        *b_Jet_electronEnergy;   //!
   TBranch        *b_Jet_photonEnergy;   //!
   TBranch        *b_UncorrJet_pt;   //!
   TBranch        *b_npuVertices;   //!
   TBranch        *b_trueInteractions;   //!
   TBranch        *b_ootnpuVertices;   //!
   TBranch        *b_npuVerticesp1;   //!
   TBranch        *b_bestVertices;   //!
   TBranch        *b_Met_pt;   //!
   TBranch        *b_Met_sumEt;   //!
   TBranch        *b_Met_phi;   //!
   TBranch        *b_Met_px;   //!
   TBranch        *b_Met_py;   //!
   TBranch        *b_Met_pz;   //!
   TBranch        *b_Gen_Met;   //!
   TBranch        *b_Met_shiftedPtUp;   //!
   TBranch        *b_Met_shiftedPtDown;   //!

};
#endif
