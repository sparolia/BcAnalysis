#ifndef _BcTo3MuAnalyzer_h
#define _BcTo3MuAnalyzer_h

// system include files
#include <memory>

// user include files 
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"




#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"
//
// class declaration
//

class BcTo3MuAnalyzer : public edm::EDAnalyzer  {
  public:
    explicit BcTo3MuAnalyzer(const edm::ParameterSet&);
    ~BcTo3MuAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    double GetLifetime(TLorentzVector b_p4, TVector3 b_vtx, TVector3 jpsi_vtx);
    short isTruthMatch(TLorentzVector sim_p4, TVector3 reco_p3);
    double getDeltaR(TLorentzVector v1_p4, TLorentzVector v2_p4);
    reco::Vertex getPVConstrainedToBS(const edm::Event& iEvent,const edm::EventSetup& iSetup, reco::Vertex pv);


  private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_Label;
    edm::EDGetTokenT<reco::GenParticleCollection> prunedGenParticles_Label;
    edm::EDGetTokenT<reco::GenParticleCollection> genPUProtons_Label;
    edm::EDGetTokenT<pat::PackedGenParticleCollection> packedGenParticles_Label;
    edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
    edm::EDGetTokenT<edm::TriggerResults> triggerResults_Label;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_Label;
    edm::EDGetTokenT<reco::BeamSpot> BS_Label;


    bool OnlyBest_;
    bool isMC_;
    bool OnlyGen_;
    
    TTree* tree_;
    // Event information
    int run, event, lumiblock;

    // Trigger match
    std::vector<short>  *triggerMatchDimuon0, *triggerMatchDimuon20, *triggerMatchDimuon25;
    std::vector<short>  *triggerMatchJpsi, *triggerMatchJpsiTk, *triggerMatchJpsiTkTk;

    std::vector<short> *truthMatchMu1Sim, *truthMatchMu2Sim, *truthMatchMuSim;
    std::vector<short> *truthMatchMu1, *truthMatchMu2, *truthMatchMu;
     
    // Primary vertex
    unsigned short nPrimaryVertices;
    std::vector<double> *primaryVertexChi2;
    std::vector<double> *primaryVertexX, *primaryVertexY, *primaryVertexZ;
    std::vector<double> *allPrimaryVertexX, *allPrimaryVertexY, *allPrimaryVertexZ;
    std::vector<double> *primaryVertexXError, *primaryVertexYError, *primaryVertexZError;
    std::vector<double> *primaryVertexXYError, *primaryVertexXZError, *primaryVertexYZError;
    
    std::vector<double> *primaryVertexBSCChi2;
    std::vector<double> *primaryVertexBSCX, *primaryVertexBSCY, *primaryVertexBSCZ;
    std::vector<double> *primaryVertexBSCXError, *primaryVertexBSCYError, *primaryVertexBSCZError;
    std::vector<double> *primaryVertexBSCXYError, *primaryVertexBSCXZError, *primaryVertexBSCYZError;

    std::vector<double> *jpsiVertexX, *jpsiVertexY, *jpsiVertexZ;
    std::vector<double> *jpsiVertexXError, *jpsiVertexYError, *jpsiVertexZError;
    std::vector<double> *jpsiVertexXYError, *jpsiVertexXZError, *jpsiVertexYZError;
    
    //muon tracks error
    std::vector<double> *mu1XError, *mu1YError, *mu1ZError;
    std::vector<double> *mu1XYError, *mu1XZError, *mu1YZError;

    std::vector<double> *mu2XError, *mu2YError, *mu2ZError;
    std::vector<double> *mu2XYError, *mu2XZError, *mu2YZError;

    std::vector<double> *muXError, *muYError, *muZError;
    std::vector<double> *muXYError, *muXZError, *muYZError;
    // Bc particle
    unsigned int nBc;
    std::vector<double> *Bc_chi2;
    std::vector<double> *Bc_vertexProbability;

    std::vector<double> *Bc_mass, *Bc_pt, *Bc_px, *Bc_py, *Bc_pz, *Bc_eta, *Bc_phi, *Bc_ct;
    std::vector<short> *Bc_charge;

    // J/Psi particles coming from Bc
    std::vector<double> *Bc_jpsi_chi2;
    std::vector<double> *Bc_jpsi_Lxy;
    std::vector<double> *mu1mu2_deltaR;
    std::vector<double> *jpsiTrk_deltaR;
    std::vector<double> *jpsiVertexProbability;

    std::vector<double> *Bc_jpsi_mass,*Bc_jpsi_pt, *Bc_jpsi_px, *Bc_jpsi_py, *Bc_jpsi_pz;
    std::vector<double> *Bc_jpsi_eta, *Bc_jpsi_phi;
    std::vector<double> *bcVertexx, *bcVertexy, *bcVertexz;


    // Muons coming from J/Psi
    std::vector<double> *Bc_jpsi_mu1_pt, *Bc_jpsi_mu1_px, *Bc_jpsi_mu1_py, *Bc_jpsi_mu1_pz;
    std::vector<double> *Bc_jpsi_mu2_pt, *Bc_jpsi_mu2_px, *Bc_jpsi_mu2_py, *Bc_jpsi_mu2_pz;
    std::vector<short> *Bc_jpsi_mu1_charge, *Bc_jpsi_mu2_charge;
    std::vector<double> *Bc_jpsi_mu1_eta, *Bc_jpsi_mu2_eta;
    std::vector<double> *Bc_jpsi_mu1_phi, *Bc_jpsi_mu2_phi;

    // Muon coming from Bc
    unsigned int nMuons;
    std::vector<double> *Bc_mu_pt, *Bc_mu_px, *Bc_mu_py, *Bc_mu_pz;
    std::vector<double> *Bc_mu_px_noFit, *Bc_mu_py_noFit, *Bc_mu_pz_noFit;
    std::vector<short> *Bc_mu_charge;
    std::vector<double> *Bc_mu_eta;
    std::vector<double> *Bc_mu_eta_noFit;
    std::vector<double> *Bc_mu_phi;
    std::vector<double> *Bc_mu_IP3D;
    std::vector<double> *Bc_mu_IP3D_error;
    
    // Muon IDs and other properties
    std::vector<double> *mu_Chi2;
    std::vector<double> *jpsi_mu1_Chi2, *jpsi_mu2_Chi2;
    std::vector<short> *mu_NumHits, *mu_NumPixelHits;
    std::vector<short> *jpsi_mu1_NumHits, *jpsi_mu1_NumPixelHits;
    std::vector<short> *jpsi_mu2_NumHits, *jpsi_mu2_NumPixelHits;
    std::vector<double> *mu_Dxy, *mu_Dz;
    std::vector<double> *jpsi_mu1_Dxy, *jpsi_mu1_Dz;
    std::vector<double> *jpsi_mu2_Dxy, *jpsi_mu2_Dz;
    std::vector<double> *muonDCA;

    std::vector<short> *isMu1Soft, *isMu2Soft;
    std::vector<short> *isMu1Global, *isMu2Global;
    std::vector<short> *isMu1Tracker, *isMu2Tracker;
    std::vector<short> *isMu1Tight, *isMu2Tight;
    std::vector<short> *isMu1PF, *isMu2PF;
    std::vector<short> *isMu1Loose, *isMu2Loose;
    std::vector<short> *isMu1Medium, *isMu2Medium;
    std::vector<short> *isMu1HighPtMuon, *isMu2HighPtMuon;

    std::vector<short> *isMuSoft;
    std::vector<short> *isMuGlobal;
    std::vector<short> *isMuTracker;
    std::vector<short> *isMuTight;
    std::vector<short> *isMuPF;
    std::vector<short> *isMuLoose;
    std::vector<short> *isMuMedium;
    std::vector<short> *isMuHighPtMuon;
    TH1F *hDzTrkPV;
    TH1F *hEventCounter;
    TH1F *hDimuon0TriggerCounter;
    TH1F *hJpsiTkTriggerCounter;
    std::vector<int> *signalDecayPresent;
    std::vector<int> *normalizationDecayPresent;
    std::vector<int> *background1DecayPresent;
    TLorentzVector gen_b_p4, gen_jpsi_p4, gen_jpsi_mu1_p4, gen_jpsi_mu2_p4, gen_mu_p4;
    TLorentzVector gen_munu_p4, gen_tau_p4, gen_taunu1_p4, gen_taunu2_p4, gen_pion_p4;
    TVector3 gen_b_vtx, gen_jpsi_vtx, gen_nutau_vtx;
    double gen_b_ct;

};



#endif
