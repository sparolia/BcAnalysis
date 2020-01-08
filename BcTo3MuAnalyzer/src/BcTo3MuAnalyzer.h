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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"




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
//
// class declaration
//

class BcTo3MuAnalyzer : public edm::EDAnalyzer  {
  public:
    explicit BcTo3MuAnalyzer(const edm::ParameterSet&);
    ~BcTo3MuAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    double GetLifetime(TLorentzVector b_p4, TVector3 b_vtx, TVector3 jpsi_vtx);
    int isTruthMatch(TLorentzVector sim_p4, TVector3 reco_p3);


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
    edm::EDGetTokenT<reco::BeamSpot> BS_Label;

    bool OnlyBest_;
    bool isMC_;
    bool OnlyGen_;
    
    TTree* tree_;
    // Event information
    int run, event, lumiblock;

	  // Trigger match
	  std::vector<int>  *triggerMatchDimuon0, *triggerMatchDimuon20, *triggerMatchDimuon25;
	  std::vector<int>  *triggerMatchJpsi, *triggerMatchJpsiTk, *triggerMatchJpsiTkTk;

    std::vector<int> *truthMatchMuPositiveSim, *truthMatchMuNegativeSim, *truthMatchUnpairedMuSim;
    std::vector<int> *truthMatchMuPositive, *truthMatchMuNegative;
     
    // Primary vertex
	  float primaryVertexChi2;
    unsigned int nPrimaryVertices;
    float primaryVertexX, primaryVertexY, primaryVertexZ;
    float primaryVertexXError, primaryVertexYError, primaryVertexZError;
    float primaryVertexXYError, primaryVertexXZError, primaryVertexYZError;

    // Bc particle
    unsigned int nBc;
    std::vector<float> *Bc_chi2;
    std::vector<float> *Bc_vertexProbability;

    std::vector<float> *Bc_decayVertexX, *Bc_decayVertexY, *Bc_decayVertexZ;
    std::vector<double> *Bc_decayVertexXError, *Bc_decayVertexYError, *Bc_decayVertexZError;
    std::vector<double> *Bc_decayVertexYXError, *Bc_decayVertexZXError, *Bc_decayVertexZYError;
    
    std::vector<float> *Bc_mass, *Bc_px, *Bc_py, *Bc_pz;
    std::vector<int> *Bc_charge;

    // J/Psi particles coming from Bc
    std::vector<float> *Bc_jpsi_chi2;
    std::vector<float> *Bc_jpsi_vertexProbability;

    std::vector<float> *Bc_jpsi_mass, *Bc_jpsi_px, *Bc_jpsi_py, *Bc_jpsi_pz;

    // Muons coming from J/Psi
    std::vector<float> *Bc_jpsi_mu1_pt, *Bc_jpsi_mu1_px, *Bc_jpsi_mu1_py, *Bc_jpsi_mu1_pz;
    std::vector<float> *Bc_jpsi_mu2_pt, *Bc_jpsi_mu2_px, *Bc_jpsi_mu2_py, *Bc_jpsi_mu2_pz;
    std::vector<int> *Bc_jpsi_mu1_charge, *Bc_jpsi_mu2_charge;

    // Muon coming from Bc
    unsigned int nMuons;
    std::vector<float> *Bc_mu_px, *Bc_mu_py, *Bc_mu_pz;
    std::vector<float> *Bc_mu_px_noFit, *Bc_mu_py_noFit, *Bc_mu_pz_noFit;
    std::vector<int> *Bc_mu_charge;

	  // Muon IDs and other properties
	  std::vector<float> *muonPositiveChi2, *muonNegativeChi2;
	  std::vector<int> *muonPositiveNumHits, *muonPositiveNumPixelHits;
	  std::vector<int> *muonNegativeNumHits, *muonNegativeNumPixelHits;
	  std::vector<float> *muonPositiveDxy, *muonPositiveDz;
	  std::vector<float> *muonNegativeDxy, *muonNegativeDz;
	  std::vector<float> *muonDCA;

	  std::vector<bool> *isMuon1Soft, *isMuon2Soft;
	  std::vector<bool> *isMuon1Tight, *isMuon2Tight;
	  std::vector<bool> *isMuon1PF, *isMuon2PF;
	  std::vector<bool> *isMuon1Loose, *isMuon2Loose;

    //std::vector<float> *deltaRMuonPositive;
    //std::vector<float> *deltaRMuonNegative;
    //std::vector<float> *deltaRUnpairedMuon;

    TH1F *hEventCounter;
    TLorentzVector gen_b_p4, gen_jpsi_p4, gen_muonPositive_p4, gen_muonNegative_p4, gen_unpairedMuon_p4;
    TVector3 gen_b_vtx, gen_jpsi_vtx;
    float gen_b_ct;

};



#endif
