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


    // J/Psi particles 
    std::vector<float> *jpsi_chi2;
    std::vector<float> *jpsi_vertexProbability;

    std::vector<float> *jpsi_mass, *jpsi_px, *jpsi_py, *jpsi_pz;

    // Muons coming from J/Psi
    std::vector<float> *jpsi_mu1_pt, *jpsi_mu1_px, *jpsi_mu1_py, *jpsi_mu1_pz;
    std::vector<float> *jpsi_mu2_pt, *jpsi_mu2_px, *jpsi_mu2_py, *jpsi_mu2_pz;
    std::vector<int> *jpsi_mu1_charge, *jpsi_mu2_charge;

	  std::vector<bool> *isMuon1Soft, *isMuon2Soft;
	  std::vector<bool> *isMuon1Tight, *isMuon2Tight;
	  std::vector<bool> *isMuon1PF, *isMuon2PF;
	  std::vector<bool> *isMuon1Loose, *isMuon2Loose;

    TH1F *hEventCounter;
    TLorentzVector gen_b_p4, gen_jpsi_p4, gen_muonPositive_p4, gen_muonNegative_p4, gen_unpairedMuon_p4;
    TVector3 gen_b_vtx, gen_jpsi_vtx;
    float gen_b_ct;

};



#endif
