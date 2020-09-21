// -*- C++ -*-
//
// Package:    RJPsiAnalyzers/BcTo3MuAnalyzer
// Class:      BcTo3MuAnalyzer
//
/**\class BcTo3MuAnalyzer BcTo3MuAnalyzer.cc RJPsiAnalyzers/BcTo3MuAnalyzer/plugins/BcTo3MuAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Gabriel Ramirez Sanchez
//         Created:  Wed, 23 Oct 2019 01:37:59 GMT
//
//

// All includes test
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"



// (Default) user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// user include files
//#include "RJPsiAnalyzers/BcTo3MuAnalyzer/src/BcTo3MuAnalyzer.h"
#include "BcTo3MuAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

// user include files: For kinematic fit
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

//
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "RecoTauTag/ImpactParameter/interface/ImpactParameterAlgorithm.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include <RecoBTag/BTagTools/interface/SignedImpactParameter3D.h>
//#include <RecoBTag/BTagTools/interface/SignedTransverseImpactParameter.h>
//#include "TrackingTools/IPTools/interface/IPTools.h"


#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "TTree.h"
#include "TLorentzVector.h"


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BcTo3MuAnalyzer::BcTo3MuAnalyzer(const edm::ParameterSet& iConfig)
 :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  trakCollection_Label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  prunedGenParticles_Label(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("prunedGenParticles"))),
  genPUProtons_Label(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genPUProtons"))),
  packedGenParticles_Label(consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("packedGenParticles"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
  triggerPrescales_Label(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),
  BS_Label(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),

  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),


  tree_(0),

  // Event information
  run(0), event(0), lumiblock(0),

  // Trigger matching
  triggerMatchDimuon0(0), triggerMatchDimuon20(0), triggerMatchDimuon25(0), 
  triggerMatchJpsi(0), triggerMatchJpsiTk(0), triggerMatchJpsiTkTk(0),

  //Sort of truth match using sim information from pat::muons 
  truthMatchMu1Sim(0), truthMatchMu2Sim(0), truthMatchMuSim(0),
  truthMatchMu1(0), truthMatchMu2(0), truthMatchMu(0),
  // Primary vertex
  nPrimaryVertices(0),
  primaryVertexChi2(0),
  primaryVertexX(0), primaryVertexY(0), primaryVertexZ(0),
  primaryVertexXError(0), primaryVertexYError(0), primaryVertexZError(0),
  primaryVertexXYError(0), primaryVertexXZError(0), primaryVertexYZError(0),
  
  primaryVertexBSCChi2(0),
  primaryVertexBSCX(0), primaryVertexBSCY(0), primaryVertexBSCZ(0),
  primaryVertexBSCXError(0), primaryVertexBSCYError(0), primaryVertexBSCZError(0),
  primaryVertexBSCXYError(0), primaryVertexBSCXZError(0), primaryVertexBSCYZError(0),

  jpsiVertexX(0), jpsiVertexY(0), jpsiVertexZ(0),
  jpsiVertexXError(0), jpsiVertexYError(0), jpsiVertexZError(0),
  jpsiVertexXYError(0), jpsiVertexXZError(0), jpsiVertexYZError(0),

  mu1XError(0), mu1YError(0), mu1ZError(0),
  mu1XYError(0), mu1XZError(0), mu1YZError(0),

  mu2XError(0), mu2YError(0), mu2ZError(0),
  mu2XYError(0), mu2XZError(0), mu2YZError(0),

  muXError(0), muYError(0), muZError(0),
  muXYError(0), muXZError(0), muYZError(0),

  nBc(0),
  Bc_chi2(0),
  Bc_vertexProbability(0),
  Bc_mass(0), Bc_pt(0), Bc_px(0), Bc_py(0), Bc_pz(0), Bc_eta(0), Bc_phi(0), Bc_ct(0), Bc_charge(0),

  Bc_jpsi_chi2(0),
  Bc_jpsi_Lxy(0),
  mu1mu2_deltaR(0),
  jpsiTrk_deltaR(0),
  jpsiVertexProbability(0),
  Bc_jpsi_mass(0), Bc_jpsi_pt(0), Bc_jpsi_px(0), Bc_jpsi_py(0), Bc_jpsi_pz(0),
  Bc_jpsi_eta(0), Bc_jpsi_phi(0), bcVertexx(0), bcVertexy(0), bcVertexz(0),

  // Muons coming from the J/Psi
  Bc_jpsi_mu1_pt(0), Bc_jpsi_mu1_px(0), Bc_jpsi_mu1_py(0), Bc_jpsi_mu1_pz(0),
  Bc_jpsi_mu2_pt(0), Bc_jpsi_mu2_px(0), Bc_jpsi_mu2_py(0), Bc_jpsi_mu2_pz(0),
  Bc_jpsi_mu1_charge(0), Bc_jpsi_mu2_charge(0),
  Bc_jpsi_mu1_eta(0), Bc_jpsi_mu2_eta(0),
  Bc_jpsi_mu1_phi(0), Bc_jpsi_mu2_phi(0),

  // Muon coming from the Bc
  nMuons(0),
  Bc_mu_pt(0), Bc_mu_px(0), Bc_mu_py(0), Bc_mu_pz(0),
  Bc_mu_px_noFit(0), Bc_mu_py_noFit(0), Bc_mu_pz_noFit(0),
  Bc_mu_charge(0),
  Bc_mu_eta(0), Bc_mu_eta_noFit(0),
  Bc_mu_phi(0),
  Bc_mu_IP3D(0),
  Bc_mu_IP3D_error(0),

  // Muon IDs and other properties
  mu_Chi2(0),
  jpsi_mu1_Chi2(0), jpsi_mu2_Chi2(0),
  mu_NumHits(0), mu_NumPixelHits(0),
  jpsi_mu1_NumHits(0), jpsi_mu1_NumPixelHits(0),
  jpsi_mu2_NumHits(0), jpsi_mu2_NumPixelHits(0),
  mu_Dxy(0), mu_Dz(0),
  jpsi_mu1_Dxy(0), jpsi_mu1_Dz(0),
  jpsi_mu2_Dxy(0), jpsi_mu2_Dz(0),
  muonDCA(0),

  
  isMu1Soft(0), isMu2Soft(0),
  isMu1Global(0), isMu2Global(0),
  isMu1Tracker(0), isMu2Tracker(0),
  isMu1Tight(0), isMu2Tight(0),
  isMu1PF(0), isMu2PF(0),
  isMu1Loose(0), isMu2Loose(0),
  isMu1Medium(0), isMu2Medium(0),
  isMu1HighPtMuon(0), isMu2HighPtMuon(0),

  isMuSoft(0),
  isMuGlobal(0),
  isMuTracker(0),
  isMuTight(0),
  isMuPF(0),
  isMuLoose(0),
  isMuMedium(0),
  isMuHighPtMuon(0),

  hDzTrkPV(0),
  hEventCounter(0),
  hDimuon0TriggerCounter(0),
  hJpsiTkTriggerCounter(0),
  signalDecayPresent(0),
  normalizationDecayPresent(0),
  background1DecayPresent(0)
  
{
   //now do what ever initialization is needed
}
  


BcTo3MuAnalyzer::~BcTo3MuAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
BcTo3MuAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  ////////////////////////////////
  // Get event content information
  ////////////////////////////////

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theBuilder);

  edm::Handle<View<pat::PackedCandidate>> thePATTrackHandle;
  iEvent.getByToken(trakCollection_Label, thePATTrackHandle);

  edm::Handle<View<pat::Muon>> thePATMuonHandle;
  iEvent.getByToken(dimuon_Label, thePATMuonHandle);

  edm::Handle<reco::GenParticleCollection> prunedGenParticlesHandle;
  iEvent.getByToken(prunedGenParticles_Label, prunedGenParticlesHandle);

  edm::Handle<reco::GenParticleCollection> genPUProtonsHandle;
  iEvent.getByToken(genPUProtons_Label, genPUProtonsHandle);

  edm::Handle<pat::PackedGenParticleCollection> packedGenParticlesHandle;
  iEvent.getByToken(packedGenParticles_Label, packedGenParticlesHandle);

  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  iEvent.getByToken(triggerResults_Label, triggerResultsHandle);

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescalesHandle;
  iEvent.getByToken(triggerPrescales_Label, triggerPrescalesHandle);
  
  //////////////////////////////
  // Trigger test
  //////////////////////////////
  
  //if(!triggerResultsHandle.isValid())
  //{
  //  LogError("BcTo3MuAnalyzer") << "Missing Trigger collection" << std::endl;
  //  return;
  //}

  //const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerResultsHandle);
  //std::cout << triggerNames.size() << std::endl;
  //for(unsigned int iT = 0; iT != triggerResultsHandle->size(); ++iT)
  //{
  //  if(triggerNames.triggerName(iT)== "HLT_Dimuon0_Jpsi3p5_Mu2_v4" || triggerNames.triggerName(iT) == "HLT_DoubleMu4_JpsiTrk_Displaced_v12")
  //  {
  //    std::cout << "Trigger " << triggerNames.triggerName(iT) << std::endl;
  //    std::cout << "Pass trigger " << triggerResultsHandle->accept(iT) << std::endl;
  //    std::cout << "Trigger prescale " << triggerPrescalesHandle->getPrescaleForIndex(iT) << std::endl;
  //    std::cout << "Trigger number " << iT << std::endl;

  //  }
  //  std::cout << "Trigger " << triggerNames.triggerName(iT) << std::endl;

      
  //  //if(triggerNames.triggerName(iT) == "HLT_DoubleMu4_JpsiTrk_Displaced_v15") hJpsiTkTriggerCounter->Fill(triggerResultsHandle->accept(iT)*1.0);
  //  if(triggerNames.triggerName(iT) == "HLT_Dimuon0_Jpsi3p5_Muon2_v4")
  //  { 
  //    hDimuon0TriggerCounter->Fill(triggerResultsHandle->accept(iT)*1.0);
  //    std::cout << "Pass trigger " << triggerResultsHandle->accept(iT) << std::endl;
  //  }
  //  if(triggerNames.triggerName(iT) == "HLT_DoubleMu4_JpsiTrk_Displaced_v12")
  //  { 
  //    hDimuon0TriggerCounter->Fill(triggerResultsHandle->accept(iT)*1.0);
  //    std::cout << "Pass trigger " << triggerResultsHandle->accept(iT) << std::endl;
  //  }
  //}
  //////////////////////////////
  // Trigger test
  //////////////////////////////
  
  

  //////////////////////////////
  // Getting generated particles pT vectors for the truthMatching
  //////////////////////////////

  gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_mu1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_mu2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_mu_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_munu_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_tau_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_taunu1_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_taunu2_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_pion_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_b_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_nutau_vtx.SetXYZ(0.,0.,0.);
  gen_b_ct = -99;
  

  ////////////////////////////////////////////////////////////////////////////////////////////////////
  // Getting the quadrivectors of the Bc, JPsi, muon positive, muon negative and the unpaired muons.
  // These quadrivector will be used to perform truth match with the reconstructed objects.
  //
  // nParticlesFound counts the number of particle matching the decay chain we are studing. 
  // For example (in the signal sample):
  // - Bc (nParticlesFound = 1) is present. Bc decay product include J/Psi(nParticlesFound = 2) and tau (nParticlesFound = 3).
  //

  bool isJpsi0Present = false, isMuonPositive0Present = false, isMuonNegative0Present = false;
  bool isBcPresent = false;
  bool isJpsi1Present = false, isMuonPositive1Present = false, isMuonNegative1Present = false;
  bool isNuMu1Present = false, isMu1Present = false;
  bool isNuTau1Present = false, isTau1Present = false, isNuMu2Present=false, isMu2Present = false, isNuTau2Present = false;
  int isNormalizationDecayPresent = 0, isSignalDecayPresent = 0, isBackground1DecayPresent = 0;
  if( isMC_ && packedGenParticlesHandle.isValid()){
    for(auto genPruned = prunedGenParticlesHandle->begin(); genPruned != prunedGenParticlesHandle->end(); ++genPruned)
    {
      if(genPruned->pdgId() == 443)
      {
        isJpsi0Present = true;
        for(size_t l = 0; l < genPruned->numberOfDaughters(); l++)
        {
          const reco::Candidate *lDaughter = genPruned->daughter(l);
          if(lDaughter->pdgId() == 13) isMuonNegative0Present = true;
          if(lDaughter->pdgId() == -13) isMuonPositive0Present = true;
        }
      }
      else if(genPruned->pdgId() == 541) 
      {
        isBcPresent = true;
        gen_b_p4.SetPtEtaPhiM(genPruned->pt(),genPruned->eta(), genPruned->phi(), genPruned->mass());
        gen_b_vtx.SetXYZ(genPruned->vx(), genPruned->vy(), genPruned->vz());
        for(size_t i = 0; i < genPruned->numberOfDaughters(); i++)
        {
          const reco::Candidate *iDaughter = genPruned->daughter(i);
          //std::cout << "Id: " << iDaughter->pdgId() << std::endl;
          if(iDaughter->pdgId() == 443) 
          {
            isJpsi1Present = true;
            gen_jpsi_p4.SetPtEtaPhiM(iDaughter->pt(),iDaughter->eta(), iDaughter->phi(), iDaughter->mass());
            gen_jpsi_vtx.SetXYZ(iDaughter->vx(), iDaughter->vy(), iDaughter->vz());
            gen_b_ct = GetLifetime(gen_b_p4, gen_b_vtx, gen_jpsi_vtx);

            for(size_t j = 0; j< iDaughter->numberOfDaughters(); j++)
            {
              const reco::Candidate *jGrandDaughter = iDaughter->daughter(j);
              if(jGrandDaughter->pdgId() == -13)
              {
                gen_jpsi_vtx.SetXYZ(jGrandDaughter->vx(), jGrandDaughter->vy(), jGrandDaughter->vz());
                isMuonPositive1Present = true;
                gen_jpsi_mu1_p4.SetPtEtaPhiM(jGrandDaughter->pt(),jGrandDaughter->eta(), jGrandDaughter->phi(), jGrandDaughter->mass());
              }
              if(jGrandDaughter->pdgId() == 13)
              {
                isMuonNegative1Present = true;
                gen_jpsi_mu2_p4.SetPtEtaPhiM(jGrandDaughter->pt(),jGrandDaughter->eta(), jGrandDaughter->phi(), jGrandDaughter->mass());
              }
            }
          }
          else
          {
            if(abs(iDaughter->pdgId()) == 13 ) 
            {
              gen_mu_p4.SetPtEtaPhiM(iDaughter->pt(),iDaughter->eta(), iDaughter->phi(), iDaughter->mass());
              isMu1Present = true;
            }
            else if(abs(iDaughter->pdgId()) == 14 ) 
            {
              isNuMu1Present = true;
              gen_munu_p4.SetPtEtaPhiM(iDaughter->pt(),iDaughter->eta(), iDaughter->phi(), iDaughter->mass());
            }

            else if(abs(iDaughter->pdgId()) == 15 )
            {
              isTau1Present = true;
              gen_tau_p4.SetPtEtaPhiM(iDaughter->pt(),iDaughter->eta(), iDaughter->phi(), iDaughter->mass());
              for(size_t k = 0; k<iDaughter->numberOfDaughters(); ++k)
              { 
                const reco::Candidate *kGrandDaughter = iDaughter->daughter(k);
                if(abs(kGrandDaughter->pdgId()) == 13)  
                {
                  isMu2Present = true;
                  gen_mu_p4.SetPtEtaPhiM(kGrandDaughter->pt(),kGrandDaughter->eta(), kGrandDaughter->phi(), kGrandDaughter->mass());
                }
                if(abs(kGrandDaughter->pdgId()) == 14) 
                {
                  isNuMu2Present = true;
                  gen_munu_p4.SetPtEtaPhiM(kGrandDaughter->pt(),kGrandDaughter->eta(), kGrandDaughter->phi(), kGrandDaughter->mass());
                }
                if(abs(kGrandDaughter->pdgId()) == 16) 
                {
                  isNuTau2Present = true; 
                  gen_nutau_vtx.SetXYZ(kGrandDaughter->vx(), kGrandDaughter->vy(), kGrandDaughter->vz());
                  gen_taunu2_p4.SetPtEtaPhiM(kGrandDaughter->pt(),kGrandDaughter->eta(), kGrandDaughter->phi(), kGrandDaughter->mass());
                }
              }
            }
            else if(abs(iDaughter->pdgId()) == 16 ) 
            {
              isNuTau1Present = true;
              gen_taunu1_p4.SetPtEtaPhiM(iDaughter->pt(),iDaughter->eta(), iDaughter->phi(), iDaughter->mass());
            }
            else if(abs(iDaughter->pdgId()) == 211 ) 
            {
              gen_pion_p4.SetPtEtaPhiM(iDaughter->pt(),iDaughter->eta(), iDaughter->phi(), iDaughter->mass());
              /*
              for(size_t k = 0; k<iDaughter->numberOfDaughters(); ++k)
              { 

                const reco::Candidate *kGrandDaughter = iDaughter->daughter(k);
                //std::cout << "daughetId: " << kGrandDaughter->pdgId() << std::endl;
                //if(abs(kGrandDaughter->pdgId()) == 13)  
                //{
                //  std::cout << "muon is present" << std::endl;
                //}
              }
              */

            }

          }
        }
      }
      if(isJpsi0Present && isMuonPositive0Present && isMuonNegative0Present) isBackground1DecayPresent = 1;
      if(isJpsi1Present && isMuonPositive1Present && isMuonNegative1Present)
      {
        if(isMu1Present && isNuMu1Present) isNormalizationDecayPresent = 1;
        else if (isTau1Present && isNuTau1Present && isNuTau2Present && isNuMu2Present && isMu2Present) isSignalDecayPresent = 1;
      } 
      isJpsi0Present = false; isMuonPositive0Present = false; isMuonNegative0Present = false;
      isJpsi1Present = false; isMuonPositive1Present = false; isMuonNegative1Present = false;
      isNuMu1Present = false; isMu1Present = false;
      isNuTau1Present = false; isTau1Present = false; isNuMu2Present=false; isMu2Present = false; isNuTau2Present = false;
      if(isSignalDecayPresent || isNormalizationDecayPresent || isBackground1DecayPresent) break;
    }
    // The nEventCounter histogram counts the total number of events with at least one Bc generated.
    if(isBcPresent) hEventCounter->Fill(1.);
    else hEventCounter->Fill(0.);
    isBcPresent = false;
    //if(isSignalDecayPresent) 
    //{
      //std::cout << "Diff vertex: " << (gen_jpsi_vtx - gen_nutau_vtx).Mag() << std::endl;
      //std::cout << "bc vertex: " << gen_b_vtx.Mag() << std::endl;
      //std::cout << "jpsi vertex: " << gen_jpsi_vtx.Mag() << std::endl;
      //std::cout << "tau vertex: " << gen_nutau_vtx.Mag() << std::endl;
    //}
  }
/*
  if( isMC_ && packedGenParticlesHandle.isValid()){
    for(auto genPacked = packedGenParticlesHandle->begin(); genPacked != packedGenParticlesHandle->end(); ++genPacked)
    {
      if(genPacked->pdgId() == 211)
      {
        std::cout << "pdf packed: " << genPacked->pdgId() << std::endl;
        std::cout << "number of daughter: " << genPacked->numberOfDaughters() << std::endl;
        for(size_t l = 0; l < genPacked->numberOfDaughters(); l++)
        {
          const reco::Candidate *lDaughter = genPacked->daughter(l);
          if(lDaughter->pdgId() == 211) std::cout << "pion is present" << std::endl;
        }
      }

    }
  }
 */ 
  //////////////////////////////
  // Get the primary vertex
  //////////////////////////////


  reco::Vertex bestVertex;
  edm::Handle<reco::VertexCollection> thePrimaryVerticesHandle;
  iEvent.getByToken(primaryVertices_Label, thePrimaryVerticesHandle);

  // Getting the first primary vertex of the container

  //for(View<reco::VertexCollection>::const_iterator primVertex = thePrimaryVerticesHandle->begin(); primVertex!= thePrimaryVerticesHandle->end(); primVertex++)
  bestVertex = *(thePrimaryVerticesHandle->begin());

  nPrimaryVertices = thePrimaryVerticesHandle->size();
  nMuons = thePATMuonHandle->size(); 

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

  /////////////////////////////////////////////////////
  // The Bc recontruction will consist first
  // on the pairing of two muons to form the J/Psi
  // and add a third muon to complete the "mesurable"
  // products of the Bc
  /////////////////////////////////////////////////////


  for(View<pat::Muon>::const_iterator patMu1 = thePATMuonHandle->begin(); patMu1 != thePATMuonHandle->end(); ++patMu1)
  {
    for(View<pat::Muon>::const_iterator patMu2 = patMu1+1; patMu2 != thePATMuonHandle->end(); ++patMu2)
    {
      // Skipping the pairing of muons with itself
      if(patMu1 == patMu2) continue;

      // Pairing only opposite signed muons
      //TODO: Save also same sign pairs to study
      if((patMu1->charge())*(patMu2->charge()) != -1) continue;
      
      // Getting tracks from the muons
      reco::TrackRef globalTrackMu1;
      reco::TrackRef globalTrackMu2;
 
      globalTrackMu1 = patMu1->track();
      globalTrackMu2 = patMu2->track();
      /*
      if(patMu1->charge() == 1)
      {
        globalTrackMu1 = patMu1->track();
        globalTrackMu2 = patMu2->track();
      }
      else
      {
        globalTrackMu1 = patMu2->track();
        globalTrackMu2 = patMu1->track();
      }
      */

      // Check for the track reference
      if(globalTrackMu1.isNull() || globalTrackMu2.isNull()) continue;

      if(!(globalTrackMu1->quality(reco::TrackBase::highPurity))) continue;
      if(!(globalTrackMu2->quality(reco::TrackBase::highPurity))) continue;
   
   
      reco::TransientTrack transientTrackMu1((*theBuilder).build(globalTrackMu1));
      reco::TransientTrack transientTrackMu2((*theBuilder).build(globalTrackMu2));


      if(!(transientTrackMu1.impactPointTSCP().isValid())) continue;
      if(!(transientTrackMu2.impactPointTSCP().isValid())) continue;


      FreeTrajectoryState trajectoryStateMu1 = transientTrackMu1.impactPointTSCP().theState();
      FreeTrajectoryState trajectoryStateMu2 = transientTrackMu2.impactPointTSCP().theState();


      // Trajectory state to calculate DCA for the two muons

      
      ClosestApproachInRPhi closestApproachMuons;
      closestApproachMuons.calculate(trajectoryStateMu1, trajectoryStateMu2);

      if(!closestApproachMuons.status()) continue;
      float dca = fabs(closestApproachMuons.distance());

      // The (PDG) mass of the muon and the insignificant mass sigma
      // to avoid singularities in the covarience matrix.
      ParticleMass muonMass = 0.10565837;
      ParticleMass jpsiMass = 3.096916;
      //ParticleMass bcMass = 6.27628;
      float muonMassSigma = muonMass*1.0e-6;

      // Creating a Kinematic particle factory.
      KinematicParticleFactoryFromTransientTrack particleFactory;

      // Inicial chi2 and ndf before the kinematic fits.
      float chi = 0.0;
      float ndf = 0.0;

      std::vector<RefCountedKinematicParticle> muonParticles;

      try
      {
        muonParticles.push_back(particleFactory.particle(transientTrackMu1, muonMass, chi, ndf, muonMassSigma));
        muonParticles.push_back(particleFactory.particle(transientTrackMu2, muonMass, chi, ndf, muonMassSigma));
      }
      catch(...)
      {
        std::cout<<" Exeption cought... continuing 1 "<< std::endl;
        continue;
      }

      KinematicParticleVertexFitter vertexFitter;
      RefCountedKinematicTree jpsiVertexFitTree;
      try
      {
        jpsiVertexFitTree = vertexFitter.fit(muonParticles);
      }
      catch(...)
      {
        std::cout<<" Exeption cought... continuing 2 "<< std::endl;
        continue;
      }
      if(!(jpsiVertexFitTree->isValid())) continue;
      jpsiVertexFitTree->movePointerToTheTop();

      RefCountedKinematicParticle jpsiVertexFit = jpsiVertexFitTree->currentParticle();
      RefCountedKinematicVertex jpsiDecayVertex = jpsiVertexFitTree->currentDecayVertex();
      reco::TransientTrack jpsiTrack =jpsiVertexFit->refittedTransientTrack();

      if(jpsiDecayVertex->chiSquared() <0.0) continue;
      if(jpsiDecayVertex->chiSquared() >50.0) continue;

      if(jpsiVertexFit->currentState().mass()<2.6) continue;
      if(jpsiVertexFit->currentState().mass()>3.6) continue;
      //if(jpsiVertexFit->currentState().globalMomentum().perp() < 8) continue;

      double jpsiProb_tmp = TMath::Prob(jpsiDecayVertex->chiSquared(),(int)jpsiDecayVertex->degreesOfFreedom());

      if(jpsiProb_tmp <0.0) continue;

      // Chosing the closest PV in Z direction to the JPsi trajectory projection.
      double dzMin = 1000000.;
      const reco::VertexCollection* vertices = thePrimaryVerticesHandle.product();
      for(reco::VertexCollection::const_iterator  primVertex = vertices->begin(); primVertex!= vertices->end(); primVertex++)
      {
        //std::cout << "prim vertex z: " << primVertex->z() << std::endl;
        if (abs(dzMin) > abs(jpsiTrack.track().dz(primVertex->position())))
        {
          bestVertex = *(primVertex);
          dzMin = jpsiTrack.track().dz(primVertex->position());
        }
      }   

      //template <typename T> std::string type_name();
      
      //std::cout << typeid(bestVertex).name() << std::endl;
      /*
      std::cout << (*bestVertex).tracksSize() << std::endl;
      for(size_t i = 0; i<bestVertex.tracksSize() ; ++i)
      {
        auto bvTkr = bestVertex.trackRefAt(i);
        std::cout << typeid(bvTkr).name() << std::endl;
        //std::cout << bvTkr.charge() << std::endl;
      }
      reco::Vertex::trackRef_iterator v_iter = bestVertex.tracks_begin();
      reco::Vertex::trackRef_iterator v_iend = bestVertex.tracks_end();
      while (v_iter != v_iend)
      {
        const reco::TrackBaseRef& tkr = *v_iter++;
        const reco::Track* tkp = &(*tkr);
        std::cout << tkp->charge() << std::endl;
      }
      */


      //reco::Vertex bestVertexBSC = getPVConstrainedToBS(iEvent, iSetup, bestVertex);
      reco::Vertex bestVertexBSC = bestVertex;

      TVector3 primaryVertex(bestVertex.x(),bestVertex.y(),0.);
      TVector3 jpsiDecayVertexPosition(jpsiDecayVertex->position().x(),
            jpsiDecayVertex->position().y(),
            0.);
      
      double Lxy_tmp = jpsiDecayVertexPosition.Mag() - primaryVertex.Dot(jpsiDecayVertexPosition) / jpsiDecayVertexPosition.Mag();



        

      //std::cout << boostToJpsiRestFrame.Mag() << std::endl;
      //TLorentzVector jpsiBoosted1 = jpsi4V;
      //jpsiBoosted1.Boost(boostToJpsiRestFrame);

      //TLorentzVector jpsiBoosted2 = jpsi4V;
      //jpsiBoosted2.Boost(-boostToJpsiRestFrame);

      //std::cout << "##########" << std::endl;
      //std::cout << "Positive" << jpsiBoosted1.Pt() << std::endl;
      //std::cout << "Negative" << jpsiBoosted2.Pt() << std::endl;
      //std::cout << "##########" << std::endl;
  
      for(View<pat::Muon>::const_iterator patMuon3 = thePATMuonHandle->begin(); patMuon3 != thePATMuonHandle->end(); ++patMuon3)
      {
        if(patMuon3->charge()==0) continue;
        if(patMu1==patMuon3) continue;
        if(patMu2==patMuon3) continue;
        
        
        //std::cout << "mu1 cov: " << patMu1->track()->covariance() << std::endl;
        reco::TrackRef globalTrackMu;
        globalTrackMu = patMuon3->track();

        if(globalTrackMu.isNull()) continue;
        //if(globalTrackMuExtra->pt()<1.0) continue;

        if(!(globalTrackMu->quality(reco::TrackBase::highPurity))) continue;
        

        //const reco::TransientTrack transientTrackMu((*theBuilder).build(globalTrackMu));
        const reco::TransientTrack transientTrackMu = (theBuilder->build(globalTrackMu));
        //FreeTrajectoryState trajectoryStateMuExtra = transientTrackMuExtra.impactPointTSCP().theState();

        if(!transientTrackMu.impactPointTSCP().isValid()) continue;


        // Now we do the kinematic fit. Constaining the JPsi mass applied to the final Bplos fit.

        std::vector<RefCountedKinematicParticle> vectorFitParticles;
        vectorFitParticles.push_back(particleFactory.particle(transientTrackMu1, muonMass, chi,ndf, muonMassSigma));
        vectorFitParticles.push_back(particleFactory.particle(transientTrackMu2, muonMass, chi,ndf, muonMassSigma));
        vectorFitParticles.push_back(particleFactory.particle(transientTrackMu, muonMass, chi,ndf, muonMassSigma));
        
        MultiTrackKinematicConstraint *jpsiConstraint = new TwoTrackMassKinematicConstraint(jpsiMass);
        KinematicConstrainedVertexFitter constrainedVertexFitter;
        RefCountedKinematicTree vertexFitTree = constrainedVertexFitter.fit(vectorFitParticles,jpsiConstraint);

        auto jpsiGlobalMomentum = jpsiVertexFit->currentState().globalMomentum();
        TLorentzVector jpsi4V;
        jpsi4V.SetPtEtaPhiM(jpsiGlobalMomentum.perp(), jpsiGlobalMomentum.eta(), jpsiGlobalMomentum.phi(), jpsiVertexFit->currentState().mass());

        TLorentzVector mu1_4V(patMu1->px(),patMu1->py(),patMu1->pz(),muonMass);
        TLorentzVector mu2_4V(patMu2->px(),patMu2->py(),patMu2->pz(),muonMass);
        TLorentzVector mu3_4V(patMuon3->px(),patMuon3->py(),patMuon3->pz(),muonMass);

        TLorentzVector bc4V(0.,0.,0.,0.);
        TVector3 bcDecayVertexPosition(0.,0.,0.);
        double bcChi2_tmp = -99;
        double bcProb_tmp = -99;
        double mu1mu2_deltaR_tmp = getDeltaR(mu1_4V, mu2_4V);
        double jpsiTrk_deltaR_tmp = getDeltaR(jpsi4V, mu3_4V);
        //if(mu1mu2_deltaR_tmp > 0.4)continue;

        //if(!vertexFitTree->isValid()) continue;
        if(vertexFitTree->isValid())
        {
          vertexFitTree->movePointerToTheTop();
          RefCountedKinematicParticle bcCandidateParticle = vertexFitTree->currentParticle();
          RefCountedKinematicVertex bcDecayVertex = vertexFitTree->currentDecayVertex();

          //if((bcCandidateParticle->currentState().mass()<0.) || (bcCandidateParticle->currentState().mass()>10.)) continue;
          //if(!bcDecayVertex->vertexIsValid()) continue;
          //if((bcDecayVertex->chiSquared()<0.0) || (bcDecayVertex->chiSquared()>50.0)) continue;
          
          bcDecayVertexPosition.SetXYZ(bcDecayVertex->position().x(),
              bcDecayVertex->position().y(),
              bcDecayVertex->position().z());

          bc4V.SetPtEtaPhiM(bcCandidateParticle->currentState().globalMomentum().perp(),
              bcCandidateParticle->currentState().globalMomentum().eta(),
              bcCandidateParticle->currentState().globalMomentum().phi(),
              bcCandidateParticle->currentState().mass());
          bcChi2_tmp = bcDecayVertex->chiSquared();
          bcProb_tmp = TMath::Prob(bcDecayVertex->chiSquared(),(int)bcDecayVertex->degreesOfFreedom());
          //if(bcProb_tmp < 0.0) continue;

          // Get children from final Bc fit
          vertexFitTree->movePointerToTheFirstChild();
          RefCountedKinematicParticle candidateMu1 = vertexFitTree->currentParticle();

          vertexFitTree->movePointerToTheNextChild();
          RefCountedKinematicParticle candidateMu2 = vertexFitTree->currentParticle();

          vertexFitTree->movePointerToTheNextChild();
          RefCountedKinematicParticle candidateMu = vertexFitTree->currentParticle();

          KinematicParameters kinematicParamMu1 = candidateMu1->currentState().kinematicParameters();
          KinematicParameters kinematicParamMu2 = candidateMu2->currentState().kinematicParameters();
          KinematicParameters kinematicParamMu = candidateMu->currentState().kinematicParameters();

          mu1_4V.SetPtEtaPhiM(kinematicParamMu1.momentum().perp(),
                            kinematicParamMu1.momentum().eta(),
                            kinematicParamMu1.momentum().phi(),
                            muonMass);
          mu2_4V.SetPtEtaPhiM(kinematicParamMu2.momentum().perp(),
                            kinematicParamMu2.momentum().eta(),
                            kinematicParamMu2.momentum().phi(),
                            muonMass);
       
          mu3_4V.SetPtEtaPhiM(kinematicParamMu.momentum().perp(),
                            kinematicParamMu.momentum().eta(),
                            kinematicParamMu.momentum().phi(),
                            muonMass);
          //GlobalVector vectorMu1(candidateMu1->currentState().globalMomentum().x(),
          //    candidateMu1->currentState().globalMomentum().y(),
          //    candidateMu1->currentState().globalMomentum().z());
          //
          //GlobalVector vectorMu2(candidateMu2->currentState().globalMomentum().x(),
          //    candidateMu2->currentState().globalMomentum().y(),
          //    candidateMu2->currentState().globalMomentum().z());

        } else
        {
          //continue;
          if(jpsiTrk_deltaR_tmp > 0.7)continue;
        }

        // Prepare variables for the NN

        //template <typename T> std::string type_name();
        // computing 3d impact parameter

        GlobalVector jpsiDirection(jpsiGlobalMomentum.x(), jpsiGlobalMomentum.y(), jpsiGlobalMomentum.z());
        reco::Vertex::Point jpsiVertexPosition(jpsiDecayVertex->position().x(), jpsiDecayVertex->position().y(), jpsiDecayVertex->position().z());
        const double err00 = jpsiDecayVertex->error().matrix()(0,0);
        const double err11 = jpsiDecayVertex->error().matrix()(1,1);
        const double err22 = jpsiDecayVertex->error().matrix()(2,2);
        const double err01 = jpsiDecayVertex->error().matrix()(0,1);
        const double err02 = jpsiDecayVertex->error().matrix()(0,2);
        const double err12 = jpsiDecayVertex->error().matrix()(1,2);
        reco::Vertex::Error jpsiVertexError;

        jpsiVertexError(0,0) = err00;
        jpsiVertexError(0,1) = err01;
        jpsiVertexError(0,2) = err02;
        jpsiVertexError(1,0) = err01;
        jpsiVertexError(1,1) = err11;
        jpsiVertexError(1,2) = err12;
        jpsiVertexError(2,0) = err02;
        jpsiVertexError(2,1) = err12;
        jpsiVertexError(2,2) = err22;

        const reco::Vertex jpsiVertex(jpsiVertexPosition, jpsiVertexError, jpsiDecayVertex->chiSquared(), jpsiDecayVertex->degreesOfFreedom(), 2);

        SignedImpactParameter3D signed_ip3D;
        Measurement1D ip3D = signed_ip3D.apply(transientTrackMu,jpsiGlobalMomentum,jpsiVertex).second;
        //std::cout << ip3D.value() << std::endl;
        //std::cout << ip3D.error() << std::endl;


        TVector3 boostToJpsiRestFrame = -jpsi4V.BoostVector();

        
        TVector3 properDecayLength = bcDecayVertexPosition - primaryVertex;
        double Bc_ct_tmp = GetLifetime(bc4V, primaryVertex, bcDecayVertexPosition);


        TLorentzVector mu_4V;
        mu_4V.SetPtEtaPhiM(patMuon3->pt(), patMuon3->eta(), patMuon3->phi(), muonMass);

        if((mu_4V + jpsi4V).M()<0. || (mu_4V + jpsi4V).M()>10.) continue;
        

        // Check for trigger matching
  
        int triggerMatchDimuon20_tmp = 0;
        int triggerMatchDimuon25_tmp = 0;
        int triggerMatchJpsiTk_tmp = 0;
        int triggerMatchDimuon0_tmp = 0;
        int triggerMatchJpsi_tmp = 0;
        int triggerMatchJpsiTkTk_tmp = 0;
  
        const pat::Muon* muon1 = &(*patMu1);
        const pat::Muon* muon2 = &(*patMu2);
        const pat::Muon* muon3 = &(*patMuon3);
        
        if(muon1->triggerObjectMatchByPath("HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*")!=nullptr) triggerMatchDimuon20_tmp = 1;
        if(muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) triggerMatchDimuon25_tmp = 1;
        if(muon1->triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v*")!=nullptr && muon3->triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v*")!=nullptr) triggerMatchDimuon0_tmp = 1;
        if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr) triggerMatchJpsiTk_tmp= 1;
        if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr) triggerMatchJpsiTkTk_tmp= 1;
        if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_Jpsi_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_Jpsi_Displaced_v*")!=nullptr) triggerMatchJpsi_tmp= 1;
        
        TVector3 reco_jpsi_mu1_p3, reco_jpsi_mu2_p3;
        reco_jpsi_mu1_p3.SetXYZ(globalTrackMu1->px(),globalTrackMu1->py(),globalTrackMu1->pz());
        reco_jpsi_mu2_p3.SetXYZ(globalTrackMu2->px(),globalTrackMu2->py(),globalTrackMu2->pz());
  
        TVector3 reco_mu_p3;
        reco_mu_p3.SetXYZ(globalTrackMu->px(),globalTrackMu->py(),globalTrackMu->pz());
        
        
        int truthMatchMu1Sim_tmp = 0;
        int truthMatchMu2Sim_tmp = 0;
        int truthMatchMu1_tmp = 0;
        int truthMatchMu2_tmp = 0;
        short truthMatchMuSim_tmp = 0;
        short truthMatchMu_tmp = 0;
        if(isMC_)
        {
          // Truth matching using sim information form pat::muons
          if (patMu1->simMotherPdgId() == 443) truthMatchMu1Sim_tmp =1; 
          if (patMu2->simMotherPdgId() == 443) truthMatchMu2Sim_tmp =1;
  
          // Truth matching using generation information from prunedGenParticles agains the muons after decay reconstruction
          if(isTruthMatch(gen_jpsi_mu1_p4, reco_jpsi_mu1_p3) || isTruthMatch(gen_jpsi_mu2_p4, reco_jpsi_mu1_p3)) truthMatchMu1_tmp = 1;
          if(isTruthMatch(gen_jpsi_mu2_p4, reco_jpsi_mu2_p3) || isTruthMatch(gen_jpsi_mu1_p4, reco_jpsi_mu2_p3)) truthMatchMu2_tmp = 1;


          //std::cout << "pion mass" << gen_pion_p4.M() << std::endl;

          if(abs(gen_pion_p4.M()) > 0.0)
          {
            if(isTruthMatch(gen_pion_p4, reco_mu_p3)) {
              truthMatchMu_tmp = 1;
              /*
              std::cout << "mu1: " << patMu1->simPdgId() << std::endl;
              std::cout << "mu2: " << patMu2->simPdgId() << std::endl;
              std::cout << "mu3: " << patMuon3->simPdgId() << std::endl;
              */
            }
          }
          else{
            if(isTruthMatch(gen_mu_p4, reco_mu_p3)) truthMatchMu_tmp = 1;
          }
          if(isSignalDecayPresent && (abs(patMuon3->simMotherPdgId()) == 15)) truthMatchMuSim_tmp = 1;
          if(isNormalizationDecayPresent && (abs(patMuon3->simMotherPdgId()) == 541)) truthMatchMuSim_tmp = 1;
          if(isBackground1DecayPresent && (abs(patMuon3->simMotherPdgId()) == 541)) truthMatchMuSim_tmp = 1;
          //std::cout << "ID: " << patMuon3->simPdgId() << std::endl;
          //std::cout << "motherID: " << patMuon3->simMotherPdgId() << std::endl;

        }

        truthMatchMuSim->push_back(truthMatchMuSim_tmp);
        truthMatchMu->push_back(truthMatchMu_tmp);

        


        // Filling the decays information
        background1DecayPresent->push_back(isBackground1DecayPresent);
        signalDecayPresent->push_back(isSignalDecayPresent);
        normalizationDecayPresent->push_back(isNormalizationDecayPresent);

        //Filling primary vertex information
        primaryVertexX->push_back(bestVertex.x());
        primaryVertexY->push_back(bestVertex.y());
        primaryVertexZ->push_back(bestVertex.z());

        primaryVertexXError->push_back(bestVertex.covariance(0,0));
        primaryVertexYError->push_back(bestVertex.covariance(1,1));
        primaryVertexZError->push_back(bestVertex.covariance(2,2));
        primaryVertexXYError->push_back(bestVertex.covariance(0,1));
        primaryVertexXZError->push_back(bestVertex.covariance(0,2));
        primaryVertexYZError->push_back(bestVertex.covariance(1,2));

        primaryVertexChi2->push_back(ChiSquaredProbability((double)(bestVertex.chi2()),(double)(bestVertex.ndof())));


        //std::cout << bestVertexBSC.x() << std::endl;
        primaryVertexBSCX->push_back(bestVertexBSC.x());
        primaryVertexBSCY->push_back(bestVertexBSC.y());
        primaryVertexBSCZ->push_back(bestVertexBSC.z());
        primaryVertexBSCXError->push_back(bestVertexBSC.covariance(0,0));
        primaryVertexBSCYError->push_back(bestVertexBSC.covariance(1,1));
        primaryVertexBSCZError->push_back(bestVertexBSC.covariance(2,2));
        primaryVertexBSCXYError->push_back(bestVertexBSC.covariance(0,1));
        primaryVertexBSCXZError->push_back(bestVertexBSC.covariance(0,2));
        primaryVertexBSCYZError->push_back(bestVertexBSC.covariance(1,2));
        primaryVertexBSCChi2->push_back(ChiSquaredProbability((double)(bestVertexBSC.chi2()),(double)(bestVertexBSC.ndof())));

        // Secondary vertex information

        jpsiVertexX->push_back(jpsiDecayVertex->position().x());
        jpsiVertexY->push_back(jpsiDecayVertex->position().y());
        jpsiVertexZ->push_back(jpsiDecayVertex->position().z());

        jpsiVertexXError->push_back(jpsiDecayVertex->error().matrix()(0,0));
        jpsiVertexYError->push_back(jpsiDecayVertex->error().matrix()(1,1));
        jpsiVertexZError->push_back(jpsiDecayVertex->error().matrix()(2,2));
        jpsiVertexXYError->push_back(jpsiDecayVertex->error().matrix()(0,1));
        jpsiVertexXZError->push_back(jpsiDecayVertex->error().matrix()(0,2));
        jpsiVertexYZError->push_back(jpsiDecayVertex->error().matrix()(1,1));

        // muon tracks errors 
        //std::cout << "mu1 cov: " << patMu1->track()->covariance() << std::endl;
        mu1XError->push_back(patMu1->track()->covariance(0,0));
        mu1YError->push_back(patMu1->track()->covariance(1,1));
        mu1ZError->push_back(patMu1->track()->covariance(2,2));
        mu1XYError->push_back(patMu1->track()->covariance(0,1));
        mu1XZError->push_back(patMu1->track()->covariance(0,2));
        mu1YZError->push_back(patMu1->track()->covariance(1,2));

        mu2XError->push_back(patMu2->track()->covariance(0,0));
        mu2YError->push_back(patMu2->track()->covariance(1,1));
        mu2ZError->push_back(patMu2->track()->covariance(2,2));
        mu2XYError->push_back(patMu2->track()->covariance(0,1));
        mu2XZError->push_back(patMu2->track()->covariance(0,2));
        mu2YZError->push_back(patMu2->track()->covariance(1,2));

        muXError->push_back(patMuon3->track()->covariance(0,0));
        muYError->push_back(patMuon3->track()->covariance(1,1));
        muZError->push_back(patMuon3->track()->covariance(2,2));
        muXYError->push_back(patMuon3->track()->covariance(0,1));
        muXZError->push_back(patMuon3->track()->covariance(0,2));
        muYZError->push_back(patMuon3->track()->covariance(1,2));
        // Filling trigger matching for the muons coming from the dimuon
        triggerMatchDimuon20->push_back(triggerMatchDimuon20_tmp);
        triggerMatchDimuon25->push_back(triggerMatchDimuon25_tmp);
        triggerMatchJpsiTk->push_back(triggerMatchJpsiTk_tmp);
        triggerMatchDimuon0->push_back(triggerMatchDimuon0_tmp);
        triggerMatchJpsiTkTk->push_back(triggerMatchJpsiTkTk_tmp);
        triggerMatchJpsi->push_back(triggerMatchJpsi_tmp);
          


        if(isMC_)
        {
          truthMatchMu1Sim->push_back(truthMatchMu1Sim_tmp);
          truthMatchMu2Sim->push_back(truthMatchMu2Sim_tmp);
          truthMatchMu1->push_back(truthMatchMu1_tmp);
          truthMatchMu2->push_back(truthMatchMu2_tmp);
        }
    
        // Muon IDs and other properties
        isMu1Soft->push_back(patMu1->isSoftMuon(bestVertex));
        isMu1Global->push_back(patMu1->isGlobalMuon());
        isMu1Tracker->push_back(patMu1->isTrackerMuon());
        isMu1Tight->push_back(patMu1->isTightMuon(bestVertex));
        isMu1PF->push_back(patMu1->isPFMuon());
        isMu1Loose->push_back(muon::isLooseMuon(*patMu1));
        isMu1Medium->push_back(muon::isMediumMuon(*patMu1));
        isMu1HighPtMuon->push_back(patMu1->isHighPtMuon(bestVertex));
    
        isMu2Soft->push_back(patMu2->isSoftMuon(bestVertex));
        isMu2Global->push_back(patMu2->isGlobalMuon());
        isMu2Tracker->push_back(patMu2->isTrackerMuon());
        isMu2Tight->push_back(patMu2->isTightMuon(bestVertex));
        isMu2PF->push_back(patMu2->isPFMuon());
        isMu2Loose->push_back(muon::isLooseMuon(*patMu2));
        isMu2Medium->push_back(muon::isMediumMuon(*patMu2));
        isMu2HighPtMuon->push_back(patMu2->isHighPtMuon(bestVertex));
        
        mu_Chi2->push_back(globalTrackMu->normalizedChi2());
        mu_NumHits->push_back(globalTrackMu->numberOfValidHits());
        mu_NumPixelHits->push_back(globalTrackMu->hitPattern().numberOfValidPixelHits());
        mu_Dxy->push_back(globalTrackMu->dxy(jpsiVertex.position()));
        mu_Dz->push_back(globalTrackMu->dz(jpsiVertex.position()));
    
        jpsi_mu1_Chi2->push_back(globalTrackMu1->normalizedChi2());
        jpsi_mu1_NumHits->push_back(globalTrackMu1->numberOfValidHits());
        jpsi_mu1_NumPixelHits->push_back(globalTrackMu1->hitPattern().numberOfValidPixelHits());
        jpsi_mu1_Dxy->push_back(globalTrackMu1->dxy(jpsiVertex.position()));
        jpsi_mu1_Dz->push_back(globalTrackMu1->dz(jpsiVertex.position()));
        //std::cout << globalTrackMu1->phi(bestVertex.position()) << std::endl;
    
         
        jpsi_mu2_Chi2->push_back(globalTrackMu2->normalizedChi2());
        jpsi_mu2_NumHits->push_back(globalTrackMu2->numberOfValidHits());
        jpsi_mu2_NumPixelHits->push_back(globalTrackMu2->hitPattern().numberOfValidPixelHits());
        jpsi_mu2_Dxy->push_back(globalTrackMu2->dxy(jpsiVertex.position()));
        jpsi_mu2_Dz->push_back(globalTrackMu2->dz(jpsiVertex.position()));
        muonDCA->push_back(dca);

        // Filling candidates variables now.
        Bc_mass->push_back(bc4V.M());
        Bc_pt->push_back(bc4V.Pt());
        Bc_px->push_back(bc4V.Px());
        Bc_py->push_back(bc4V.Py());
        Bc_pz->push_back(bc4V.Pz());
        Bc_eta->push_back(bc4V.Eta());
        Bc_phi->push_back(bc4V.Phi());
        Bc_ct->push_back(Bc_ct_tmp);
        //Bc_charge->push_back(bcCandidateParticle->currentState().particleCharge());
        Bc_charge->push_back(patMu1->charge() + patMu2->charge() + patMuon3->charge());

        // Filling childen variables
        // First for the muon coming directly from the Bc

        Bc_mu_pt->push_back(mu3_4V.Pt());
        Bc_mu_px->push_back(mu3_4V.Px());
        Bc_mu_py->push_back(mu3_4V.Py());
        Bc_mu_pz->push_back(mu3_4V.Pz());
        Bc_mu_charge->push_back(patMuon3->charge());
        Bc_mu_eta->push_back(mu3_4V.Eta());
        Bc_mu_phi->push_back(mu3_4V.Phi());
        Bc_mu_IP3D->push_back(ip3D.value());
        Bc_mu_IP3D_error->push_back(ip3D.error());
        
        Bc_mu_px_noFit->push_back(patMuon3->px());
        Bc_mu_py_noFit->push_back(patMuon3->py());
        Bc_mu_pz_noFit->push_back(patMuon3->pz());
        Bc_mu_eta_noFit->push_back(patMuon3->eta());

        // For thhe JPsi and the muon from its decay
  
        Bc_jpsi_mass->push_back(jpsiVertexFit->currentState().mass());
        Bc_jpsi_pt->push_back(jpsiGlobalMomentum.perp());
        Bc_jpsi_px->push_back(jpsiGlobalMomentum.x());
        Bc_jpsi_py->push_back(jpsiGlobalMomentum.y());
        Bc_jpsi_pz->push_back(jpsiGlobalMomentum.z());
        Bc_jpsi_eta->push_back(jpsiGlobalMomentum.eta());
        Bc_jpsi_phi->push_back(jpsiGlobalMomentum.phi());
        bcVertexx->push_back(bcDecayVertexPosition.X());
        bcVertexy->push_back(bcDecayVertexPosition.Y());
        bcVertexz->push_back(bcDecayVertexPosition.Z());
 
        Bc_jpsi_mu1_pt->push_back(mu1_4V.Pt());
        Bc_jpsi_mu1_px->push_back(mu1_4V.Px());
        Bc_jpsi_mu1_py->push_back(mu1_4V.Py());
        Bc_jpsi_mu1_pz->push_back(mu1_4V.Pz());
        Bc_jpsi_mu1_charge->push_back(patMu1->charge());
        Bc_jpsi_mu1_eta->push_back(mu1_4V.Eta());
        Bc_jpsi_mu1_phi->push_back(mu1_4V.Phi());
  
        Bc_jpsi_mu2_pt->push_back(mu2_4V.Pt());
        Bc_jpsi_mu2_px->push_back(mu2_4V.Px());
        Bc_jpsi_mu2_py->push_back(mu2_4V.Py());
        Bc_jpsi_mu2_pz->push_back(mu2_4V.Pz());
        Bc_jpsi_mu2_charge->push_back(patMu2->charge());
        Bc_jpsi_mu2_eta->push_back(mu2_4V.Eta());
        Bc_jpsi_mu2_phi->push_back(mu2_4V.Phi());
  
        Bc_jpsi_chi2->push_back(jpsiDecayVertex->chiSquared());
        Bc_jpsi_Lxy->push_back(Lxy_tmp);
        mu1mu2_deltaR->push_back(mu1mu2_deltaR_tmp);
        jpsiTrk_deltaR->push_back(jpsiTrk_deltaR_tmp);
        Bc_chi2->push_back(bcChi2_tmp);
        jpsiVertexProbability->push_back(jpsiProb_tmp);
        Bc_vertexProbability->push_back(bcProb_tmp);
        
        isMuSoft->push_back(patMuon3->isSoftMuon(bestVertex));
        isMuGlobal->push_back(patMuon3->isGlobalMuon());
        isMuTracker->push_back(patMuon3->isTrackerMuon());
        isMuTight->push_back(patMuon3->isTightMuon(bestVertex));
        isMuPF->push_back(patMuon3->isPFMuon());
        isMuLoose->push_back(muon::isLooseMuon(*patMuon3));
        isMuMedium->push_back(muon::isMediumMuon(*patMuon3));
        isMuHighPtMuon->push_back(patMuon3->isHighPtMuon(bestVertex));

        nBc++;
      }
      //if(isMC_) break; // Making sure we only save one Jpsi for each event.
    }
    //if(isMC_ && jpsiFound) break;
  }
  
  if(nBc < 1)
  {
    // For thhe JPsi and the muon from its decay

    triggerMatchDimuon20->push_back(-99);
    triggerMatchDimuon25->push_back(-99);
    triggerMatchJpsiTk->push_back(-99);
    triggerMatchDimuon0->push_back(-99);
    triggerMatchJpsiTkTk->push_back(-99);
    triggerMatchJpsi->push_back(-99);

    if(isMC_)
    {
      truthMatchMu1Sim->push_back(-99);
      truthMatchMu2Sim->push_back(-99);
      truthMatchMu1->push_back(-99);
      truthMatchMu2->push_back(-99);
      truthMatchMuSim->push_back(-99);
      truthMatchMu->push_back(-99);
    }

    background1DecayPresent->push_back(-99);
    signalDecayPresent->push_back(-99);
    normalizationDecayPresent->push_back(-99);
     // Muon IDs and other properties
    primaryVertexChi2->push_back(-99);
    
    primaryVertexX->push_back(-99);
    primaryVertexXError->push_back(-99);
    primaryVertexY->push_back(-99);
    primaryVertexYError->push_back(-99);
    primaryVertexZ->push_back(-99);
    primaryVertexZError->push_back(-99);
    primaryVertexXYError->push_back(-99);
    primaryVertexXZError->push_back(-99);
    primaryVertexYZError->push_back(-99);

    primaryVertexBSCX->push_back(-99);
    primaryVertexBSCXError->push_back(-99);
    primaryVertexBSCY->push_back(-99);
    primaryVertexBSCYError->push_back(-99);
    primaryVertexBSCZ->push_back(-99);
    primaryVertexBSCZError->push_back(-99);
    primaryVertexBSCXYError->push_back(-99);
    primaryVertexBSCXZError->push_back(-99);
    primaryVertexBSCYZError->push_back(-99);
    primaryVertexBSCChi2->push_back(-99);
    
    jpsiVertexX->push_back(-99);
    jpsiVertexXError->push_back(-99);
    jpsiVertexY->push_back(-99);
    jpsiVertexYError->push_back(-99);
    jpsiVertexZ->push_back(-99);
    jpsiVertexZError->push_back(-99);
    jpsiVertexXYError->push_back(-99);
    jpsiVertexXZError->push_back(-99);
    jpsiVertexYZError->push_back(-99);

    mu1XError->push_back(-99);
    mu1YError->push_back(-99);
    mu1ZError->push_back(-99);
    mu1XYError->push_back(-99);
    mu1XZError->push_back(-99);
    mu1YZError->push_back(-99);

    mu2XError->push_back(-99);
    mu2YError->push_back(-99);
    mu2ZError->push_back(-99);
    mu2XYError->push_back(-99);
    mu2XZError->push_back(-99);
    mu2YZError->push_back(-99);
    
    muXError->push_back(-99);
    muYError->push_back(-99);
    muZError->push_back(-99);
    muXYError->push_back(-99);
    muXZError->push_back(-99);
    muYZError->push_back(-99);

    isMu1Soft->push_back(-99);
    isMu1Global->push_back(-99);
    isMu1Tracker->push_back(-99);
    isMu1Tight->push_back(-99);
    isMu1PF->push_back(-99);
    isMu1Loose->push_back(-99);
    isMu1Medium->push_back(-99);
    isMu1HighPtMuon->push_back(-99);

    isMu2Soft->push_back(-99);
    isMu2Global->push_back(-99);
    isMu2Tracker->push_back(-99);
    isMu2Tight->push_back(-99);
    isMu2PF->push_back(-99);
    isMu2Loose->push_back(-99);
    isMu2Medium->push_back(-99);
    isMu2HighPtMuon->push_back(-99);

    mu_Chi2->push_back(-99);
    mu_NumHits->push_back(-99);
    mu_NumPixelHits->push_back(-99);
    mu_Dxy->push_back(-99);
    mu_Dz->push_back(-99);

    jpsi_mu1_Chi2->push_back(-99);
    jpsi_mu1_NumHits->push_back(-99);
    jpsi_mu1_NumPixelHits->push_back(-99);
    jpsi_mu1_Dxy->push_back(-99);
    jpsi_mu1_Dz->push_back(-99);

    jpsi_mu2_Chi2->push_back(-99);
    jpsi_mu2_NumHits->push_back(-99);
    jpsi_mu2_NumPixelHits->push_back(-99);
    jpsi_mu2_Dxy->push_back(-99);
    jpsi_mu2_Dz->push_back(-99);
    muonDCA->push_back(-99);

    Bc_mass->push_back(-99);
    Bc_pt->push_back(-99);
    Bc_px->push_back(-99);
    Bc_py->push_back(-99);
    Bc_pz->push_back(-99);
    Bc_eta->push_back(-99);
    Bc_eta->push_back(-99);
    Bc_ct->push_back(-99);
    Bc_charge->push_back(-99);

    Bc_mu_pt->push_back(-99);
    Bc_mu_px->push_back(-99);
    Bc_mu_py->push_back(-99);
    Bc_mu_pz->push_back(-99);
    Bc_mu_charge->push_back(-99);
    Bc_mu_eta->push_back(-99);
    Bc_mu_phi->push_back(-99);
    Bc_mu_IP3D->push_back(-99);
    Bc_mu_IP3D_error->push_back(-99);
    
    Bc_mu_px_noFit->push_back(-99);
    Bc_mu_py_noFit->push_back(-99);
    Bc_mu_pz_noFit->push_back(-99);
    Bc_mu_eta_noFit->push_back(-99);

    Bc_jpsi_mass->push_back(-99);
    Bc_jpsi_pt->push_back(-99);
    Bc_jpsi_px->push_back(-99);
    Bc_jpsi_py->push_back(-99);
    Bc_jpsi_pz->push_back(-99);
    Bc_jpsi_eta->push_back(-99);
    Bc_jpsi_phi->push_back(-99);
    bcVertexx->push_back(-99);
    bcVertexy->push_back(-99);
    bcVertexz->push_back(-99);
  
    Bc_jpsi_mu1_pt->push_back(-99);
    Bc_jpsi_mu1_px->push_back(-99);
    Bc_jpsi_mu1_py->push_back(-99);
    Bc_jpsi_mu1_pz->push_back(-99);
    Bc_jpsi_mu1_charge->push_back(-99);
    Bc_jpsi_mu1_eta->push_back(-99);
    Bc_jpsi_mu1_phi->push_back(-99);
  
    Bc_jpsi_mu2_pt->push_back(-99);
    Bc_jpsi_mu2_px->push_back(-99);
    Bc_jpsi_mu2_py->push_back(-99);
    Bc_jpsi_mu2_pz->push_back(-99);
    Bc_jpsi_mu2_charge->push_back(-99);
    Bc_jpsi_mu2_eta->push_back(-99);
    Bc_jpsi_mu2_phi->push_back(-99);
    

    Bc_jpsi_chi2->push_back(-99);
    Bc_jpsi_Lxy->push_back(-99);
    mu1mu2_deltaR->push_back(-99);
    jpsiTrk_deltaR->push_back(-99);
    jpsiVertexProbability->push_back(-99);
    Bc_chi2->push_back(-99);
    Bc_vertexProbability->push_back(-99);

    isMuSoft->push_back(-99);
    isMuGlobal->push_back(-99);
    isMuTracker->push_back(-99);
    isMuTight->push_back(-99);
    isMuPF->push_back(-99);
    isMuLoose->push_back(-99);
    isMuMedium->push_back(-99);
    isMuHighPtMuon->push_back(-99);
  }
  bool saveTree = true;
  if(nBc > 0 && saveTree)
  {
    tree_->Fill();
  }
  
  nBc = 0;
  nMuons = 0;
  
  
  nPrimaryVertices = 0;
  primaryVertexChi2->clear();
  
  primaryVertexX->clear();
  primaryVertexXError->clear();
  primaryVertexY->clear();
  primaryVertexYError->clear();
  primaryVertexZ->clear();
  primaryVertexZError->clear();
  primaryVertexXYError->clear();
  primaryVertexXZError->clear();
  primaryVertexYZError->clear();

  primaryVertexBSCX->clear();
  primaryVertexBSCXError->clear();
  primaryVertexBSCY->clear();
  primaryVertexBSCYError->clear();
  primaryVertexBSCZ->clear();
  primaryVertexBSCZError->clear();
  primaryVertexBSCXYError->clear();
  primaryVertexBSCXZError->clear();
  primaryVertexBSCYZError->clear();
  primaryVertexBSCChi2->clear();
  
  jpsiVertexX->clear();
  jpsiVertexXError->clear();
  jpsiVertexY->clear();
  jpsiVertexYError->clear();
  jpsiVertexZ->clear();
  jpsiVertexZError->clear();
  jpsiVertexXYError->clear();
  jpsiVertexXZError->clear();
  jpsiVertexYZError->clear();

  mu1XError->clear();
  mu1YError->clear();
  mu1ZError->clear();
  mu1XYError->clear();
  mu1XZError->clear();
  mu1YZError->clear();
  
  mu2XError->clear();
  mu2YError->clear();
  mu2ZError->clear();
  mu2XYError->clear();
  mu2XZError->clear();
  mu2YZError->clear();
  
  muXError->clear();
  muYError->clear();
  muZError->clear();
  muXYError->clear();
  muXZError->clear();
  muYZError->clear();

  mu_Chi2->clear();
  mu_Dxy->clear();
  mu_Dz->clear();
  mu_NumHits->clear();
  mu_NumPixelHits->clear();
 
  jpsi_mu1_Chi2->clear();
  jpsi_mu1_Dxy->clear();
  jpsi_mu1_Dz->clear();
  jpsi_mu1_NumHits->clear();
  jpsi_mu1_NumPixelHits->clear();
  
  jpsi_mu2_Chi2->clear();
  jpsi_mu2_Dxy->clear();
  jpsi_mu2_Dz->clear();
  jpsi_mu2_NumHits->clear();
  jpsi_mu2_NumPixelHits->clear();
  
  muonDCA->clear();
  
  triggerMatchDimuon0->clear();
  triggerMatchDimuon20->clear();
  triggerMatchDimuon25->clear();
  triggerMatchJpsi->clear();
  triggerMatchJpsiTk->clear();
  triggerMatchJpsiTkTk->clear();

  
  if(isMC_)
  {
    truthMatchMu1Sim->clear();
    truthMatchMu2Sim->clear();
    truthMatchMu1->clear();
    truthMatchMu2->clear();
    truthMatchMuSim->clear();
    truthMatchMu->clear();
  }
  
  isMu1Soft->clear();
  isMu1Global->clear();
  isMu1Tracker->clear();
  isMu1Tight->clear();
  isMu1PF->clear();
  isMu1Loose->clear();
  isMu1Medium->clear();
  isMu1HighPtMuon->clear();
  
  isMu2Soft->clear();
  isMu2Global->clear();
  isMu2Tracker->clear();
  isMu2Tight->clear();
  isMu2PF->clear();
  isMu2Loose->clear();
  isMu2Medium->clear();
  isMu2HighPtMuon->clear();

  isMuSoft->clear();
  isMuGlobal->clear();
  isMuTracker->clear();
  isMuTight->clear();
  isMuPF->clear();
  isMuLoose->clear();
  isMuMedium->clear();
  isMuHighPtMuon->clear();
  signalDecayPresent->clear();
  normalizationDecayPresent->clear();
  background1DecayPresent->clear();

  Bc_charge->clear();
  Bc_mass->clear();
  Bc_pt->clear();
  Bc_px->clear();
  Bc_py->clear();
  Bc_pz->clear();
  Bc_eta->clear();
  Bc_phi->clear();
  Bc_ct->clear();

  Bc_mu_pt->clear();
  Bc_mu_px->clear();
  Bc_mu_py->clear();
  Bc_mu_pz->clear();
  Bc_mu_charge->clear();
  Bc_mu_eta->clear();
  Bc_mu_phi->clear();
  Bc_mu_IP3D->clear();
  Bc_mu_IP3D_error->clear();
  
  Bc_mu_px_noFit->clear();
  Bc_mu_py_noFit->clear();
  Bc_mu_pz_noFit->clear();
  Bc_mu_eta_noFit->clear();

  Bc_jpsi_mass->clear();
  Bc_jpsi_pt->clear();
  Bc_jpsi_px->clear();
  Bc_jpsi_py->clear();
  Bc_jpsi_pz->clear();
  Bc_jpsi_eta->clear();
  Bc_jpsi_phi->clear();
  bcVertexx->clear();
  bcVertexy->clear();
  bcVertexz->clear();
  
  Bc_jpsi_mu1_pt->clear();
  Bc_jpsi_mu1_px->clear();
  Bc_jpsi_mu1_py->clear();
  Bc_jpsi_mu1_pz->clear();
  Bc_jpsi_mu1_charge->clear();
  Bc_jpsi_mu1_eta->clear();
  Bc_jpsi_mu1_phi->clear();
  
  Bc_jpsi_mu2_pt->clear();
  Bc_jpsi_mu2_px->clear();
  Bc_jpsi_mu2_py->clear();
  Bc_jpsi_mu2_pz->clear();
  Bc_jpsi_mu2_charge->clear();
  Bc_jpsi_mu2_eta->clear();
  Bc_jpsi_mu2_phi->clear();

  Bc_jpsi_chi2->clear();
  Bc_jpsi_Lxy->clear();
  mu1mu2_deltaR->clear();
  jpsiTrk_deltaR->clear();
  jpsiVertexProbability->clear();

  Bc_chi2->clear();
  Bc_vertexProbability->clear();
}

double
BcTo3MuAnalyzer::GetLifetime(TLorentzVector b_p4, TVector3 b_vtx, TVector3 jpsi_vtx)
{
  TVector3 deltaVtx = jpsi_vtx - b_vtx;
  TVector3 b_p3 = b_p4.Vect();
  deltaVtx.SetZ(0.);
  b_p3.SetZ(0.);
  Double_t lxy = deltaVtx.Dot(b_p3)/b_p3.Mag();
  return lxy*b_p4.M()/b_p3.Mag();
}
//double
//BcTo3MuAnalyzer::GetDCA(TLorentzVector b_p4, TVector3 vtx1, TVector3 vtx2)
//{
//  TVector3 deltaVtx = jpsi_vtx - b_vtx;
//  TVector3 b_p3 = b_p4.Vect();
//  Double_t lxy = deltaVtx.Dot(b_p3)/b_p3.Mag();
//  return lxy*b_p4.M()/b_p3.Mag();
//}
short
BcTo3MuAnalyzer::isTruthMatch(TLorentzVector sim_p4, TVector3 reco_p3)
{
  TVector3 sim_p3 = sim_p4.Vect();
  Double_t simRecoDeltaR = sim_p3.DeltaR(reco_p3);
  Double_t ptPercent = 1.;
  if(sim_p3.Pt() !=0.) ptPercent = abs(sim_p3.Pt()-reco_p3.Pt())/sim_p3.Pt();
  int match = 0;
  if(simRecoDeltaR < 0.1 && ptPercent < 0.3) match = 1;
  return match;
}
double
BcTo3MuAnalyzer::getDeltaR(TLorentzVector v1_p4, TLorentzVector v2_p4)
{
  TVector3 v1 = v1_p4.Vect();
  TVector3 v2 = v2_p4.Vect();
  return v1.DeltaR(v2);
}

reco::Vertex
BcTo3MuAnalyzer::getPVConstrainedToBS(const edm::Event& iEvent, const edm::EventSetup& iSetup, reco::Vertex pv)
{
  using namespace edm;

  edm::ESHandle<TransientTrackBuilder> theBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theBuilder);

  edm::Handle<edm::View<pat::PackedCandidate>> thePATTrackHandle;
  iEvent.getByToken(trakCollection_Label, thePATTrackHandle);

  edm::Handle<reco::BeamSpot> bs;
  iEvent.getByToken(BS_Label, bs);
  
  std::vector<reco::TransientTrack> tTrksForRefit;
  for(edm::View<pat::PackedCandidate>::const_iterator iTrk = thePATTrackHandle->begin(); iTrk != thePATTrackHandle->end(); ++iTrk)
  //for(unsigned int i = 0; i < thePATTrackHandle->size() ; ++i)
  {
    //const reco::Track* iTrk = thePATTrackHandle.at(i);
    //iTrk->dz(pv);
    double trkDzMax = 0.8;
    if(iTrk->charge() == 0) continue;
    if(!(iTrk->hasTrackDetails())) continue;
    /*
    if(abs(iTrk->pdgId()) == 211) continue;
    if(abs(iTrk->pdgId()) == 11) continue;
    if(abs(iTrk->pdgId()) == 13) continue;
    */
    
    //std::cout << iTrk->pdgId() << std::endl;
    //std::cout << iTrk->vertexRef().id() << std::endl;
    //std::cout << iTrk->dz(pv.position()) << std::endl;
    if(iTrk->dz(pv.position()) < trkDzMax){
      hDzTrkPV->Fill(iTrk->dz(pv.position()));
      //std::cout << iTrk->dz(pv.position()) << std::endl;
      const reco::TransientTrack tTrk((*theBuilder).build(iTrk->pseudoTrack()));
      tTrksForRefit.push_back(tTrk);
    }
  }

  KalmanVertexFitter kvf(true);
  reco::Vertex fVtx;
  //TransientVertex tv = kvf.vertex(tTrksForRefit, *bs);
  TransientVertex tv = kvf.vertex(tTrksForRefit);
  fVtx=tv;

  return fVtx;

}

// ------------ method called once each job just before starting event loop  ------------
void
BcTo3MuAnalyzer::beginJob()
{
  std::cout << "Begin analyzer job" << std::endl;
  
  edm::Service<TFileService> fs;

  hEventCounter = fs->make<TH1F>("nGeneratedEvents", "nGeneratedEvents", 10, 0., 10.);
  hDzTrkPV = fs->make<TH1F>("dzTrkPV", "dzTrkPV", 1000, 0., 10.);
  hDimuon0TriggerCounter = fs->make<TH1F>("dimuon0TriggerCounter", "dimuon0TriggerCounter", 10, 0., 10.);
  hJpsiTkTriggerCounter = fs->make<TH1F>("jpsiTkTriggerCounter", "jpsiTkTriggerCounter", 10, 0., 10.);

  tree_ = fs->make<TTree>("ntuple","Bc+ -> J/Psi mu+ ntuple");

  tree_->Branch("nPrimaryVertices",&nPrimaryVertices);
  tree_->Branch("run",&run, "run/I");
  tree_->Branch("event",&event, "event/I");
  tree_->Branch("lumiblock",&lumiblock, "lumiblock/I");
  
  tree_->Branch("nBc",&nBc,"nBc/i");
  tree_->Branch("nMuons",&nMuons,"nMuons/i");

  tree_->Branch("primaryVertexX",&primaryVertexX);
  tree_->Branch("primaryVertexY",&primaryVertexY);
  tree_->Branch("primaryVertexZ",&primaryVertexZ);
  tree_->Branch("primaryVertexXError",&primaryVertexXError);
  tree_->Branch("primaryVertexYError",&primaryVertexYError);
  tree_->Branch("primaryVertexZError",&primaryVertexZError);
  tree_->Branch("primaryVertexXYError",&primaryVertexXYError);
  tree_->Branch("primaryVertexYZError",&primaryVertexYZError);
  tree_->Branch("primaryVertexXZError",&primaryVertexXZError);
  tree_->Branch("primaryVertexChi2",&primaryVertexChi2);

  tree_->Branch("primaryVertexBSCX",&primaryVertexBSCX);
  tree_->Branch("primaryVertexBSCY",&primaryVertexBSCY);
  tree_->Branch("primaryVertexBSCZ",&primaryVertexBSCZ);
  tree_->Branch("primaryVertexBSCXError",&primaryVertexBSCXError);
  tree_->Branch("primaryVertexBSCYError",&primaryVertexBSCYError);
  tree_->Branch("primaryVertexBSCZError",&primaryVertexBSCZError);
  tree_->Branch("primaryVertexBSCXYError",&primaryVertexBSCXYError);
  tree_->Branch("primaryVertexBSCYZError",&primaryVertexBSCYZError);
  tree_->Branch("primaryVertexBSCXZError",&primaryVertexBSCXZError);
  tree_->Branch("primaryVertexBSCChi2",&primaryVertexBSCChi2);

  tree_->Branch("jpsiVertexX",&jpsiVertexX);
  tree_->Branch("jpsiVertexY",&jpsiVertexY);
  tree_->Branch("jpsiVertexZ",&jpsiVertexZ);
  tree_->Branch("jpsiVertexXError",&jpsiVertexXError);
  tree_->Branch("jpsiVertexYError",&jpsiVertexYError);
  tree_->Branch("jpsiVertexZError",&jpsiVertexZError);
  tree_->Branch("jpsiVertexXYError",&jpsiVertexXYError);
  tree_->Branch("jpsiVertexYZError",&jpsiVertexYZError);
  tree_->Branch("jpsiVertexXZError",&jpsiVertexXZError);
  tree_->Branch("jpsiVertexProbability",&jpsiVertexProbability);

  tree_->Branch("mu1XError",&mu1XError);
  tree_->Branch("mu1YError",&mu1YError);
  tree_->Branch("mu1ZError",&mu1ZError);
  tree_->Branch("mu1XYError",&mu1XYError);
  tree_->Branch("mu1YZError",&mu1YZError);
  tree_->Branch("mu1XZError",&mu1XZError);

  tree_->Branch("mu2XError",&mu2XError);
  tree_->Branch("mu2YError",&mu2YError);
  tree_->Branch("mu2ZError",&mu2ZError);
  tree_->Branch("mu2XYError",&mu2XYError);
  tree_->Branch("mu2YZError",&mu2YZError);
  tree_->Branch("mu2XZError",&mu2XZError);

  tree_->Branch("muXError",&muXError);
  tree_->Branch("muYError",&muYError);
  tree_->Branch("muZError",&muZError);
  tree_->Branch("muXYError",&muXYError);
  tree_->Branch("muYZError",&muYZError);
  tree_->Branch("muXZError",&muXZError);

  tree_->Branch("bcVertexx",&bcVertexx);
  tree_->Branch("bcVertexy",&bcVertexy);
  tree_->Branch("bcVertexz",&bcVertexz);
  tree_->Branch("Bc_vertexProbability",&Bc_vertexProbability);
  
  tree_->Branch("Bc_charge",&Bc_charge);
  tree_->Branch("Bc_mass",&Bc_mass);
  tree_->Branch("Bc_pt",&Bc_pt);
  tree_->Branch("Bc_px",&Bc_px);
  tree_->Branch("Bc_py",&Bc_py);
  tree_->Branch("Bc_pz",&Bc_pz);
  tree_->Branch("Bc_eta",&Bc_eta);
  tree_->Branch("Bc_phi",&Bc_phi);
  tree_->Branch("Bc_ct", &Bc_ct);
  tree_->Branch("Bc_chi2",&Bc_chi2);

  tree_->Branch("Bc_mu_charge",&Bc_mu_charge);
  tree_->Branch("Bc_mu_pt",&Bc_mu_pt);
  tree_->Branch("Bc_mu_px",&Bc_mu_px);
  tree_->Branch("Bc_mu_py",&Bc_mu_py);
  tree_->Branch("Bc_mu_pz",&Bc_mu_pz);
  tree_->Branch("Bc_mu_eta",&Bc_mu_eta);
  tree_->Branch("Bc_mu_phi",&Bc_mu_phi);
  tree_->Branch("Bc_mu_IP3D_error",&Bc_mu_IP3D_error);
  tree_->Branch("Bc_mu_IP3D",&Bc_mu_IP3D);
  tree_->Branch("Bc_mu_px_noFit",&Bc_mu_px_noFit);
  tree_->Branch("Bc_mu_py_noFit",&Bc_mu_py_noFit);
  tree_->Branch("Bc_mu_pz_noFit",&Bc_mu_pz_noFit);
  tree_->Branch("Bc_mu_eta_noFit",&Bc_mu_eta_noFit);


  tree_->Branch("Bc_jpsi_mass",&Bc_jpsi_mass);
  tree_->Branch("Bc_jpsi_pt",&Bc_jpsi_pt);
  tree_->Branch("Bc_jpsi_px",&Bc_jpsi_px);
  tree_->Branch("Bc_jpsi_py",&Bc_jpsi_py);
  tree_->Branch("Bc_jpsi_pz",&Bc_jpsi_pz);
  tree_->Branch("Bc_jpsi_eta",&Bc_jpsi_eta);
  tree_->Branch("Bc_jpsi_phi",&Bc_jpsi_phi);
  tree_->Branch("Bc_jpsi_chi2",&Bc_jpsi_chi2);
  tree_->Branch("Bc_jpsi_Lxy",&Bc_jpsi_Lxy);
  tree_->Branch("mu1mu2_deltaR",&mu1mu2_deltaR);
  tree_->Branch("jpsiTrk_deltaR",&jpsiTrk_deltaR);

  tree_->Branch("Bc_jpsi_mu1_charge",&Bc_jpsi_mu1_charge);
  tree_->Branch("Bc_jpsi_mu1_pt",&Bc_jpsi_mu1_pt);
  tree_->Branch("Bc_jpsi_mu1_px",&Bc_jpsi_mu1_px);
  tree_->Branch("Bc_jpsi_mu1_py",&Bc_jpsi_mu1_py);
  tree_->Branch("Bc_jpsi_mu1_pz",&Bc_jpsi_mu1_pz);
  tree_->Branch("Bc_jpsi_mu1_eta",&Bc_jpsi_mu1_eta);
  tree_->Branch("Bc_jpsi_mu1_phi",&Bc_jpsi_mu1_phi);

  tree_->Branch("Bc_jpsi_mu2_charge",&Bc_jpsi_mu2_charge);
  tree_->Branch("Bc_jpsi_mu2_pt",&Bc_jpsi_mu2_pt);
  tree_->Branch("Bc_jpsi_mu2_px",&Bc_jpsi_mu2_px);
  tree_->Branch("Bc_jpsi_mu2_py",&Bc_jpsi_mu2_py);
  tree_->Branch("Bc_jpsi_mu2_pz",&Bc_jpsi_mu2_pz);
  tree_->Branch("Bc_jpsi_mu2_eta",&Bc_jpsi_mu2_eta);
  tree_->Branch("Bc_jpsi_mu2_phi",&Bc_jpsi_mu2_phi);
  
  tree_->Branch("jpsi_mu2_Chi2",&jpsi_mu2_Chi2);
  tree_->Branch("jpsi_mu2_NumHits",&jpsi_mu2_NumHits);
  tree_->Branch("jpsi_mu2_NumPixelHits",&jpsi_mu2_NumPixelHits);
  tree_->Branch("jpsi_mu2_Dxy",&jpsi_mu2_Dxy);
  tree_->Branch("jpsi_mu2_Dz",&jpsi_mu2_Dz);
  tree_->Branch("jpsi_mu1_Chi2",&jpsi_mu1_Chi2);
  tree_->Branch("jpsi_mu1_NumHits",&jpsi_mu1_NumHits);
  tree_->Branch("jpsi_mu1_NumPixelHits",&jpsi_mu1_NumPixelHits);
  tree_->Branch("jpsi_mu1_Dxy",&jpsi_mu1_Dxy);
  tree_->Branch("jpsi_mu1_Dz",&jpsi_mu1_Dz);
  tree_->Branch("mu_Chi2",&mu_Chi2);
  tree_->Branch("mu_NumHits",&mu_NumHits);
  tree_->Branch("mu_NumPixelHits",&mu_NumPixelHits);
  tree_->Branch("mu_Dxy",&mu_Dxy);
  tree_->Branch("mu_Dz",&mu_Dz);
  tree_->Branch("muonDCA",&muonDCA);
  
  tree_->Branch("triggerMatchJpsi",&triggerMatchJpsi);
  tree_->Branch("triggerMatchJpsiTk",&triggerMatchJpsiTk);
  tree_->Branch("triggerMatchJpsiTkTk",&triggerMatchJpsiTkTk);
  tree_->Branch("triggerMatchDimuon0",&triggerMatchDimuon0);
  tree_->Branch("triggerMatchDimuon20",&triggerMatchDimuon20);
  tree_->Branch("triggerMatchDimuon25",&triggerMatchDimuon25);

  
  tree_->Branch("isMu1Soft",&isMu1Soft);
  tree_->Branch("isMu1Global",&isMu1Global);
  tree_->Branch("isMu1Tracker",&isMu1Tracker);
  tree_->Branch("isMu1Tight",&isMu1Tight);
  tree_->Branch("isMu1PF",&isMu1PF);
  tree_->Branch("isMu1Loose",&isMu1Loose);
  tree_->Branch("isMu1Medium",&isMu1Medium);
  tree_->Branch("isMu1HighPtMuon",&isMu1HighPtMuon);

  tree_->Branch("isMu2Soft",&isMu2Soft);
  tree_->Branch("isMu2Global",&isMu2Global);
  tree_->Branch("isMu2Tracker",&isMu2Tracker);
  tree_->Branch("isMu2Tight",&isMu2Tight);
  tree_->Branch("isMu2PF",&isMu2PF);
  tree_->Branch("isMu2Loose",&isMu2Loose);
  tree_->Branch("isMu2Medium",&isMu2Medium);
  tree_->Branch("isMu2HighPtMuon",&isMu2HighPtMuon);
  
  tree_->Branch("isMuSoft",&isMuSoft);
  tree_->Branch("isMuGlobal",&isMuGlobal);
  tree_->Branch("isMuTracker",&isMuTracker);
  tree_->Branch("isMuTight",&isMuTight);
  tree_->Branch("isMuPF",&isMuPF);
  tree_->Branch("isMuLoose",&isMuLoose);
  tree_->Branch("isMuMedium",&isMuMedium);
  tree_->Branch("isMuHighPtMuon",&isMuHighPtMuon);

  tree_->Branch("signalDecayPresent", &signalDecayPresent);
  tree_->Branch("normalizationDecayPresent", &normalizationDecayPresent);
  tree_->Branch("background1DecayPresent", &background1DecayPresent);

  if(isMC_)
  {

    tree_->Branch("truthMatchMu1Sim",&truthMatchMu1Sim);
    tree_->Branch("truthMatchMu2Sim",&truthMatchMu2Sim);
    tree_->Branch("truthMatchMuSim",&truthMatchMuSim);
    tree_->Branch("truthMatchMu1",&truthMatchMu1);
    tree_->Branch("truthMatchMu2",&truthMatchMu2);
    tree_->Branch("truthMatchMu",&truthMatchMu);
    
    tree_->Branch("gen_b_p4", "TLorentzVector", &gen_b_p4);
    tree_->Branch("gen_jpsi_p4", "TLorentzVector", &gen_jpsi_p4);
    tree_->Branch("gen_jpsi_mu1_p4", "TLorentzVector", &gen_jpsi_mu1_p4);
    tree_->Branch("gen_jpsi_mu2_p4", "TLorentzVector", &gen_jpsi_mu2_p4);
    tree_->Branch("gen_mu_p4", "TLorentzVector", &gen_mu_p4);
    tree_->Branch("gen_munu_p4", "TLorentzVector", &gen_munu_p4);
    tree_->Branch("gen_tau_p4", "TLorentzVector", &gen_tau_p4);
    tree_->Branch("gen_taunu1_p4", "TLorentzVector", &gen_taunu1_p4);
    tree_->Branch("gen_taunu2_p4", "TLorentzVector", &gen_taunu2_p4);
    tree_->Branch("gen_pion_p4", "TLorentzVector", &gen_pion_p4);
    tree_->Branch("gen_b_vtx", "TVector3", &gen_b_vtx);
    tree_->Branch("gen_jpsi_vtx", "TVector3", &gen_jpsi_vtx);
    tree_->Branch("gen_nutau_vtx", "TVector3", &gen_nutau_vtx);
    tree_->Branch("gen_b_ct", &gen_b_ct, "gen_b_ct/D");
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void
BcTo3MuAnalyzer::endJob()
{
  tree_->GetDirectory()->cd();
  tree_->Write();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BcTo3MuAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BcTo3MuAnalyzer);

