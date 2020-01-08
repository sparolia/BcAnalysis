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


// user include files: For kinematic fit
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"



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
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
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
  truthMatchMuPositiveSim(0), truthMatchMuNegativeSim(0),
  truthMatchMuPositive(0), truthMatchMuNegative(0),
  // Primary vertex
  primaryVertexChi2(0),
  nPrimaryVertices(0),
  primaryVertexX(0), primaryVertexY(0), primaryVertexZ(0),
  primaryVertexXError(0), primaryVertexYError(0), primaryVertexZError(0),
  primaryVertexXYError(0), primaryVertexXZError(0), primaryVertexYZError(0),

  // J/Psi particles coming from Bc
  jpsi_chi2(0),
  jpsi_vertexProbability(0),
  jpsi_mass(0), jpsi_pt(0), jpsi_px(0), jpsi_py(0), jpsi_pz(0),

  // Muons coming from the J/Psi
  jpsi_mu1_pt(0), jpsi_mu1_px(0), jpsi_mu1_py(0), jpsi_mu1_pz(0),
  jpsi_mu2_pt(0), jpsi_mu2_px(0), jpsi_mu2_py(0), jpsi_mu2_pz(0),
  jpsi_mu1_charge(0), jpsi_mu2_charge(0),

  // Muon IDs and other properties
  muonPositiveChi2(0), muonNegativeChi2(0),
  muonPositiveNumHits(0), muonPositiveNumPixelHits(0),
  muonNegativeNumHits(0), muonNegativeNumPixelHits(0),
  muonPositiveDxy(0), muonPositiveDz(0),
  muonNegativeDxy(0), muonNegativeDz(0),
  muonDCA(0),

  
  isMuon1Soft(0), isMuon2Soft(0),
  isMuon1Tight(0), isMuon2Tight(0),
  isMuon1PF(0), isMuon2PF(0),
  isMuon1Loose(0), isMuon2Loose(0),

  hEventCounter(0)
  

  
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


  //////////////////////////////
  // Getting generated particles pT vectors for the truthMatching
  //////////////////////////////

  gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muonPositive_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_muonNegative_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_unpairedMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
  gen_b_vtx.SetXYZ(0.,0.,0.);
  gen_jpsi_vtx.SetXYZ(0.,0.,0.);
  gen_b_ct = -99;
  
  
  if( isMC_ && packedGenParticlesHandle.isValid()){
    for(auto genPruned = prunedGenParticlesHandle->begin(); genPruned != prunedGenParticlesHandle->end(); ++genPruned)
    {
      if(genPruned->pdgId() == 541) 
      {
        hEventCounter->Fill(1.);
        break;
      }
    }
  }

  if( isMC_ && packedGenParticlesHandle.isValid()){
    int nParticlesFound = 0;
    for(auto genPruned = prunedGenParticlesHandle->begin(); genPruned != prunedGenParticlesHandle->end(); ++genPruned)
    {
      nParticlesFound = 0;
      if(genPruned->pdgId() == 541) 
      {
        for(size_t i = 0; i < genPruned->numberOfDaughters(); ++i)
        {
          const reco::Candidate *iDaughter = genPruned->daughter(i);
          if(iDaughter->pdgId() == 541)
          {
            nParticlesFound++;
            gen_b_p4.SetPtEtaPhiM(iDaughter->pt(),iDaughter->eta(), iDaughter->phi(), iDaughter->mass());
            gen_b_vtx.SetXYZ(iDaughter->vx(), iDaughter->vy(), iDaughter->vz());
          }
          if(iDaughter->pdgId() == 443) 
          {
            nParticlesFound++;
            gen_jpsi_p4.SetPtEtaPhiM(iDaughter->pt(),iDaughter->eta(), iDaughter->phi(), iDaughter->mass());
            gen_jpsi_vtx.SetXYZ(iDaughter->vx(), iDaughter->vy(), iDaughter->vz());
            gen_b_ct = GetLifetime(gen_b_p4, gen_b_vtx, gen_jpsi_vtx);

            for(size_t j = 0; j< iDaughter->numberOfDaughters(); ++j)
            {
              const reco::Candidate *jGrandDaughter = iDaughter->daughter(j);
              if(jGrandDaughter->pdgId() == -13)
              {
                nParticlesFound++;
                gen_muonPositive_p4.SetPtEtaPhiM(jGrandDaughter->pt(),jGrandDaughter->eta(), jGrandDaughter->phi(), jGrandDaughter->mass());
              }
              if(jGrandDaughter->pdgId() == 13)
              {
                nParticlesFound++;
                gen_muonNegative_p4.SetPtEtaPhiM(jGrandDaughter->pt(),jGrandDaughter->eta(), jGrandDaughter->phi(), jGrandDaughter->mass());
              }
            }
          }
          else
          {
            if(abs(iDaughter->pdgId()) == 13) 
            {
              nParticlesFound++;
              gen_unpairedMuon_p4.SetPtEtaPhiM(iDaughter->pt(),iDaughter->eta(), iDaughter->phi(), iDaughter->mass());

            }
            else
            {
              if(abs(iDaughter->pdgId()) == 15)
              {
                for(size_t k = 0; k<iDaughter->numberOfDaughters(); ++k)
                { 
                  const reco::Candidate *kGrandDaughter = iDaughter->daughter(k);
                  if(abs(kGrandDaughter->pdgId()) == 13)
                  {
                    nParticlesFound++;
                    gen_unpairedMuon_p4.SetPtEtaPhiM(kGrandDaughter->pt(),kGrandDaughter->eta(), kGrandDaughter->phi(), kGrandDaughter->mass());
                  }
                }
              }
            }
          }
        }
        // TODO: Find a intelligen way to avoid the following if:
        if(nParticlesFound == 1) continue;
        if (nParticlesFound != 4)
        {
          gen_b_p4.SetPtEtaPhiM(0.,0.,0.,0.);
          gen_jpsi_p4.SetPtEtaPhiM(0.,0.,0.,0.);
          gen_muonPositive_p4.SetPtEtaPhiM(0.,0.,0.,0.);
          gen_muonNegative_p4.SetPtEtaPhiM(0.,0.,0.,0.);
          gen_unpairedMuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
          gen_b_vtx.SetXYZ(0.,0.,0.);
          gen_jpsi_vtx.SetXYZ(0.,0.,0.);
          gen_b_ct = -99;

        }
      }
    }
  }
  

  //////////////////////////////
  // Get the primary vertex
  //////////////////////////////


  reco::Vertex bestVertex;
  edm::Handle<reco::VertexCollection> thePrimaryVerticesHandle;
  iEvent.getByToken(primaryVertices_Label, thePrimaryVerticesHandle);

  // Getting the first primary vertex of the container

  bestVertex = *(thePrimaryVerticesHandle->begin());

  primaryVertexX = bestVertex.x();
  primaryVertexY = bestVertex.y();
  primaryVertexZ = bestVertex.z();
  
  primaryVertexXError = bestVertex.covariance(0,0);
  primaryVertexYError = bestVertex.covariance(1,1);
  primaryVertexZError = bestVertex.covariance(2,2);
  primaryVertexXYError = bestVertex.covariance(0,1);
  primaryVertexXZError = bestVertex.covariance(0,2);
  primaryVertexYZError = bestVertex.covariance(1,2);


  
  primaryVertexChi2 = ChiSquaredProbability((double)(bestVertex.chi2()),(double)(bestVertex.ndof()));
  nPrimaryVertices = thePrimaryVerticesHandle->size();

  lumiblock = iEvent.id().luminosityBlock();
  run = iEvent.id().run();
  event = iEvent.id().event();

  /////////////////////////////////////////////////////
  // The Bc recontruction will consist first
  // on the pairing of two muons to form the J/Psi
  // and add a third muon to complete the "mesurable"
  // products of the Bc
  /////////////////////////////////////////////////////

  for(View<pat::Muon>::const_iterator patMuon1 = thePATMuonHandle->begin(); patMuon1 != thePATMuonHandle->end(); ++patMuon1)
  {

    for(View<pat::Muon>::const_iterator patMuon2 = patMuon1+1; patMuon2 != thePATMuonHandle->end(); ++patMuon2)
    {
      // Skipping the pairing of muons with itself
      if(patMuon1 == patMuon2) continue;

      // Pairing only opposite signed muons
      //TODO: Save also same sign pairs to study
      if((patMuon1->charge())*(patMuon2->charge()) != -1) continue;
      
      // Getting tracks from the muons
      reco::TrackRef globalTrackMuPositive;
      reco::TrackRef globalTrackMuNegative;
 
      if(patMuon1->charge() == 1)
      {
        globalTrackMuPositive = patMuon1->track();
        globalTrackMuNegative = patMuon2->track();
      }
      else
      {
        globalTrackMuPositive = patMuon2->track();
        globalTrackMuNegative = patMuon1->track();
      }

      // Check for the track reference
      if(globalTrackMuPositive.isNull() || globalTrackMuNegative.isNull()) continue;

      /*G: Preserve low pt muons
      if(globalTrackMuPositive->pt()<4.0) continue;
      if(globalTrackMuNegative->pt()<4.0) continue;
      */

      if(!(globalTrackMuPositive->quality(reco::TrackBase::highPurity))) continue;
      if(!(globalTrackMuNegative->quality(reco::TrackBase::highPurity))) continue;
   
   
      reco::TransientTrack transientTrackMuPositive((*theBuilder).build(globalTrackMuPositive));
      reco::TransientTrack transientTrackMuNegative((*theBuilder).build(globalTrackMuNegative));


      if(!(transientTrackMuPositive.impactPointTSCP().isValid())) continue;
      if(!(transientTrackMuNegative.impactPointTSCP().isValid())) continue;


      FreeTrajectoryState trajectoryStateMuPositive = transientTrackMuPositive.impactPointTSCP().theState();
      FreeTrajectoryState trajectoryStateMuNegative = transientTrackMuNegative.impactPointTSCP().theState();


      // Trajectory state to calculate DCA for the two muons

      ClosestApproachInRPhi closestApproachMuons;
      closestApproachMuons.calculate(trajectoryStateMuPositive, trajectoryStateMuNegative);

      if(!closestApproachMuons.status()) continue;
      float dca = fabs(closestApproachMuons.distance());

      // The (PDG) mass of the muon and the insignificant mass sigma
      // to avoid singularities in the covarience matrix.
      ParticleMass muonMass = 0.10565837;
      ParticleMass jpsiMass = 3.096916;
      float muonMassSigma = muonMass*1.0e-6;

      // Creating a Kinematic particle factory.
      KinematicParticleFactoryFromTransientTrack particleFactory;

      // Inicial chi2 and ndf before the kinematic fits.
      float chi = 0.0;
      float ndf = 0.0;

      std::vector<RefCountedKinematicParticle> muonParticles;

      try
      {
        muonParticles.push_back(particleFactory.particle(transientTrackMuPositive, muonMass, chi, ndf, muonMassSigma));
        muonParticles.push_back(particleFactory.particle(transientTrackMuNegative, muonMass, chi, ndf, muonMassSigma));
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
      RefCountedKinematicVertex jpsiVertexFit_vertex = jpsiVertexFitTree->currentDecayVertex();

      if(jpsiVertexFit_vertex->chiSquared() <0.0) continue;
      if(jpsiVertexFit_vertex->chiSquared() >50.0) continue;

      if(jpsiVertexFit->currentState().mass()<2.9) continue;
      if(jpsiVertexFit->currentState().mass()>3.3) continue;




      double jpsiProb_tmp = TMath::Prob(jpsiVertexFit_vertex->chiSquared(),(int)jpsiVertexFit_vertex->degreesOfFreedom());

      if(jpsiProb_tmp <0.01) continue;

      TLorentzVector jpsi4V;
      auto jpsiGlobalMomentum = jpsiVertexFit->currentState().globalMomentum();
      jpsi4V.SetXYZM(jpsiGlobalMomentum.x(), jpsiGlobalMomentum.y(), jpsiGlobalMomentum.z(), jpsiVertexFit->currentState().mass());

        // For thhe JPsi and the muon from its decay

        jpsi_mass->push_back(jpsiVertexFit->currentState().mass());
        jpsi_pt->push_back(jpsiGlobalMomentum.perp());
        jpsi_px->push_back(jpsiGlobalMomentum.x());
        jpsi_py->push_back(jpsiGlobalMomentum.y());
        jpsi_pz->push_back(jpsiGlobalMomentum.z());

        jpsi_mu1_pt->push_back(globalTrackMuPositive->pt());
        jpsi_mu1_px->push_back(globalTrackMuPositive->px());
        jpsi_mu1_py->push_back(globalTrackMuPositive->py());
        jpsi_mu1_pz->push_back(globalTrackMuPositive->pz());
        jpsi_mu1_charge->push_back(globalTrackMuPositive->charge());

        
        jpsi_mu2_pt->push_back(globalTrackMuNegative->pt());
        jpsi_mu2_px->push_back(globalTrackMuNegative->px());
        jpsi_mu2_py->push_back(globalTrackMuNegative->py());
        jpsi_mu2_pz->push_back(globalTrackMuNegative->pz());
        jpsi_mu2_charge->push_back(globalTrackMuNegative->charge());
        

        jpsi_chi2->push_back(jpsiVertexFit_vertex->chiSquared());


        jpsi_vertexProbability->push_back(jpsiProb_tmp);

        // Check for trigger matching

        int triggerMatchDimuon20_tmp = 0;
        int triggerMatchDimuon25_tmp = 0;
        int triggerMatchJpsiTk_tmp = 0;
        int triggerMatchDimuon0_tmp = 0;
        int triggerMatchJpsi_tmp = 0;
        int triggerMatchJpsiTkTk_tmp = 0;

        const pat::Muon* muon1 = &(*patMuon1);
        const pat::Muon* muon2 = &(*patMuon2);
      
        if(muon1->triggerObjectMatchByPath("HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon20_Jpsi_Barrel_Seagulls_v*")!=nullptr) triggerMatchDimuon20_tmp = 1;
        if(muon1->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon25_Jpsi_v*")!=nullptr) triggerMatchDimuon25_tmp = 1;
        if(muon1->triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_Dimuon0_Jpsi3p5_Muon2_v*")!=nullptr) triggerMatchDimuon0_tmp = 1;
        if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrk_Displaced_v*")!=nullptr) triggerMatchJpsiTk_tmp= 1;
        if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_JpsiTrkTrk_Displaced_v*")!=nullptr) triggerMatchJpsiTkTk_tmp= 1;
        if(muon1->triggerObjectMatchByPath("HLT_DoubleMu4_Jpsi_Displaced_v*")!=nullptr && muon2->triggerObjectMatchByPath("HLT_DoubleMu4_Jpsi_Displaced_v*")!=nullptr) triggerMatchJpsi_tmp= 1;
      
 	      triggerMatchDimuon20->push_back(triggerMatchDimuon20_tmp);
    		triggerMatchDimuon25->push_back(triggerMatchDimuon25_tmp);
    		triggerMatchJpsiTk->push_back(triggerMatchJpsiTk_tmp);
    		triggerMatchDimuon0->push_back(triggerMatchDimuon0_tmp);
    		triggerMatchJpsiTkTk->push_back(triggerMatchJpsiTkTk_tmp);
    		triggerMatchJpsi->push_back(triggerMatchJpsi_tmp);

        // Sort of truth matching using sim information form pat::muons

        int truthMatchMuonPositiveSim = 0;
        int truthMatchMuonNegativeSim = 0;
        
        if (patMuon1->simMotherPdgId() == 443) truthMatchMuonPositiveSim =1; 
        if (patMuon2->simMotherPdgId() == 443) truthMatchMuonNegativeSim =1;

        
        truthMatchMuPositiveSim->push_back(truthMatchMuonPositiveSim);
        truthMatchMuNegativeSim->push_back(truthMatchMuonNegativeSim);

        // Sort of truth matching using generation information from prunedGenParticles agains the muons after decay reconstruction
        TVector3 reco_muonPositive_p3, reco_muonNegative_p3;
        reco_muonPositive_p3.SetXYZ(globalTrackMuPositive->px(),globalTrackMuPositive->py(),globalTrackMuPositive->pz());
        reco_muonNegative_p3.SetXYZ(globalTrackMuNegative->px(),globalTrackMuNegative->py(),globalTrackMuNegative->pz());

        int truthMatchMuonPositive = 0;
        int truthMatchMuonNegative = 0;

        if(isTruthMatch(gen_muonPositive_p4, reco_muonPositive_p3)) truthMatchMuonPositive = 1;
        if(isTruthMatch(gen_muonNegative_p4, reco_muonNegative_p3)) truthMatchMuonNegative = 1;
        
        truthMatchMuPositive->push_back(truthMatchMuonPositive);
        truthMatchMuNegative->push_back(truthMatchMuonNegative);
    		// Muon IDs and other properties

 	      isMuon1Soft->push_back(patMuon1->isSoftMuon(bestVertex));
    		isMuon1Tight->push_back(patMuon1->isTightMuon(bestVertex));
    		isMuon1PF->push_back(patMuon1->isPFMuon());
    		isMuon1Loose->push_back(muon::isLooseMuon(*patMuon1));

    		isMuon2Soft->push_back(patMuon2->isSoftMuon(bestVertex));
    		isMuon2Tight->push_back(patMuon2->isTightMuon(bestVertex));
    		isMuon2PF->push_back(patMuon2->isPFMuon());
    		isMuon2Loose->push_back(muon::isLooseMuon(*patMuon1));

    		muonPositiveChi2->push_back(globalTrackMuPositive->normalizedChi2());
    		muonPositiveNumHits->push_back(globalTrackMuPositive->numberOfValidHits());
    		muonPositiveNumPixelHits->push_back(globalTrackMuPositive->hitPattern().numberOfValidPixelHits());
    		muonPositiveDxy->push_back(globalTrackMuPositive->dxy(bestVertex.position()));
    		muonPositiveDz->push_back(globalTrackMuPositive->dz(bestVertex.position()));

 	 
    		muonNegativeChi2->push_back(globalTrackMuNegative->normalizedChi2());
    		muonNegativeNumHits->push_back(globalTrackMuNegative->numberOfValidHits());
    		muonNegativeNumPixelHits->push_back(globalTrackMuNegative->hitPattern().numberOfValidPixelHits());
    		muonNegativeDxy->push_back(globalTrackMuNegative->dxy(bestVertex.position()));
    		muonNegativeDz->push_back(globalTrackMuNegative->dz(bestVertex.position()));
    		muonDCA->push_back(dca);
    }
  }
  
  if(0)
  {
    // For thhe JPsi and the muon from its decay

    jpsi_mass->push_back(-99);
    jpsi_pt->push_back(-99);
    jpsi_px->push_back(-99);
    jpsi_py->push_back(-99);
    jpsi_pz->push_back(-99);

    jpsi_mu1_pt->push_back(-99);
    jpsi_mu1_px->push_back(-99);
    jpsi_mu1_py->push_back(-99);
    jpsi_mu1_pz->push_back(-99);
    jpsi_mu1_charge->push_back(-99);

    
    jpsi_mu2_pt->push_back(-99);
    jpsi_mu2_px->push_back(-99);
    jpsi_mu2_py->push_back(-99);
    jpsi_mu2_pz->push_back(-99);
    jpsi_mu2_charge->push_back(-99);
    
    jpsi_chi2->push_back(-99);

    jpsi_vertexProbability->push_back(-99);

 	  triggerMatchDimuon20->push_back(-99);
    triggerMatchDimuon25->push_back(-99);
    triggerMatchJpsiTk->push_back(-99);
    triggerMatchDimuon0->push_back(-99);
    triggerMatchJpsiTkTk->push_back(-99);
    triggerMatchJpsi->push_back(-99);

    truthMatchMuPositiveSim->push_back(-99);
    truthMatchMuNegativeSim->push_back(-99);

    truthMatchMuPositive->push_back(-99);
    truthMatchMuNegative->push_back(-99);

     // Muon IDs and other properties

 	  isMuon1Soft->push_back(-99);
    isMuon1Tight->push_back(-99);
    isMuon1PF->push_back(-99);
    isMuon1Loose->push_back(-99);

    isMuon2Soft->push_back(-99);
    isMuon2Tight->push_back(-99);
    isMuon2PF->push_back(-99);
    isMuon2Loose->push_back(-99);

    muonPositiveChi2->push_back(-99);
    muonPositiveNumHits->push_back(-99);
    muonPositiveNumPixelHits->push_back(-99);
    muonPositiveDxy->push_back(-99);
    muonPositiveDz->push_back(-99);

 	 
    muonNegativeChi2->push_back(-99);
    muonNegativeNumHits->push_back(-99);
    muonNegativeNumPixelHits->push_back(-99);
    muonNegativeDxy->push_back(-99);
    muonNegativeDz->push_back(-99);
    muonDCA->push_back(-99);
  
  }
	tree_->Fill();
  
  jpsi_mass->clear();
  jpsi_pt->clear();
  jpsi_px->clear();
  jpsi_py->clear();
  jpsi_pz->clear();
  
  jpsi_mu1_pt->clear();
  jpsi_mu1_px->clear();
  jpsi_mu1_py->clear();
  jpsi_mu1_pz->clear();
  jpsi_mu2_pt->clear();
  jpsi_mu2_px->clear();
  jpsi_mu2_py->clear();
  jpsi_mu2_pz->clear();
 
  jpsi_chi2->clear();
  jpsi_vertexProbability->clear();
  
  
  nPrimaryVertices=0;
  primaryVertexChi2=0;
  
  primaryVertexX=0;
  primaryVertexXError=0;
  primaryVertexY=0;
  primaryVertexYError=0;
  primaryVertexZ=0;
  primaryVertexZError=0;
  primaryVertexXYError=0;
  primaryVertexXZError=0;
  primaryVertexYZError=0;
  
  muonPositiveChi2->clear();
  muonPositiveDxy->clear();
  muonPositiveDz->clear();
  muonPositiveNumHits->clear();
  muonPositiveNumPixelHits->clear();
  
  muonNegativeChi2->clear();
  muonNegativeDxy->clear();
  muonNegativeDz->clear();
  muonNegativeNumHits->clear();
  muonNegativeNumPixelHits->clear();
  
  muonDCA->clear();
  
  triggerMatchDimuon0->clear();
  triggerMatchDimuon20->clear();
  triggerMatchDimuon25->clear();
  triggerMatchJpsi->clear();
  triggerMatchJpsiTk->clear();
  triggerMatchJpsiTkTk->clear();
  
  truthMatchMuPositiveSim->clear();
  truthMatchMuNegativeSim->clear();

  truthMatchMuPositive->clear();
  truthMatchMuNegative->clear();
  
  isMuon1Soft->clear();
  isMuon1Tight->clear();
  isMuon1PF->clear();
  isMuon1Loose->clear();
  
  isMuon2Soft->clear();
  isMuon2Tight->clear();
  isMuon2PF->clear();
  isMuon2Loose->clear();
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
int
BcTo3MuAnalyzer::isTruthMatch(TLorentzVector sim_p4, TVector3 reco_p3)
{
  TVector3 sim_p3 = sim_p4.Vect();
  Double_t simRecoDeltaR = sim_p3.DeltaR(reco_p3);
  int match = 0;
  if(simRecoDeltaR < 0.3) match = 1;
  return match;
}

// ------------ method called once each job just before starting event loop  ------------
void
BcTo3MuAnalyzer::beginJob()
{
	std::cout << "Begin analyzer job" << std::endl;

	edm::Service<TFileService> fs;
  hEventCounter = fs->make<TH1F>("nGeneratedEvents", "nGeneratedEvents", 10, 0., 10.);
	tree_ = fs->make<TTree>("ntuple","Bc+ -> J/Psi mu+ ntuple");


	tree_->Branch("jpsi_mass",&jpsi_mass);
	tree_->Branch("jpsi_pt",&jpsi_pt);
	tree_->Branch("jpsi_px",&jpsi_px);
	tree_->Branch("jpsi_py",&jpsi_py);
	tree_->Branch("jpsi_pz",&jpsi_pz);

	tree_->Branch("jpsi_mu1_pt",&jpsi_mu1_pt);
	tree_->Branch("jpsi_mu1_px",&jpsi_mu1_px);
	tree_->Branch("jpsi_mu1_py",&jpsi_mu1_py);
	tree_->Branch("jpsi_mu1_pz",&jpsi_mu1_pz);
	tree_->Branch("jpsi_mu1_charge",&jpsi_mu1_charge);

	tree_->Branch("jpsi_mu2_pt",&jpsi_mu2_pt);
	tree_->Branch("jpsi_mu2_px",&jpsi_mu2_px);
	tree_->Branch("jpsi_mu2_py",&jpsi_mu2_py);
	tree_->Branch("jpsi_mu2_pz",&jpsi_mu2_pz);
	tree_->Branch("jpsi_mu2_charge",&jpsi_mu2_charge);
	


	tree_->Branch("jpsi_chi2",&jpsi_chi2);
	tree_->Branch("jpsi_vertexProbability",&jpsi_vertexProbability);

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


	tree_->Branch("nPrimaryVertices",&nPrimaryVertices);
	tree_->Branch("run",&run, "run/I");
	tree_->Branch("event",&event, "event/I");
	tree_->Branch("lumiblock",&lumiblock, "lumiblock/I");

	tree_->Branch("muonNegativeChi2",&muonNegativeChi2);
	tree_->Branch("muonNegativeNumHits",&muonNegativeNumHits);
	tree_->Branch("muonNegativeNumPixelHits",&muonNegativeNumPixelHits);
	tree_->Branch("muonNegativeDxy",&muonNegativeDxy);
	tree_->Branch("muonNegativeDz",&muonNegativeDz);
	tree_->Branch("muonPositiveChi2",&muonPositiveChi2);
	tree_->Branch("muonPositiveNumHits",&muonPositiveNumHits);
	tree_->Branch("muonPositiveNumPixelHits",&muonPositiveNumPixelHits);
	tree_->Branch("muonPositiveDxy",&muonPositiveDxy);
	tree_->Branch("muonPositiveDz",&muonPositiveDz);
	tree_->Branch("muonDCA",&muonDCA);

	tree_->Branch("triggerMatchJpsi",&triggerMatchJpsi);
	tree_->Branch("triggerMatchJpsiTk",&triggerMatchJpsiTk);
	tree_->Branch("triggerMatchJpsiTkTk",&triggerMatchJpsiTkTk);
	tree_->Branch("triggerMatchDimuon0",&triggerMatchDimuon0);
	tree_->Branch("triggerMatchDimuon20",&triggerMatchDimuon20);
	tree_->Branch("triggerMatchDimuon25",&triggerMatchDimuon25);

  tree_->Branch("truthMatchMuPositiveSim",&truthMatchMuPositiveSim);
  tree_->Branch("truthMatchMuNegativeSim",&truthMatchMuNegativeSim);

  tree_->Branch("truthMatchMuPositive",&truthMatchMuPositive);
  tree_->Branch("truthMatchMuNegative",&truthMatchMuNegative);
  
	tree_->Branch("isMuon1Soft",&isMuon1Soft);
	tree_->Branch("isMuon1Tight",&isMuon1Tight);
	tree_->Branch("isMuon1PF",&isMuon1PF);
	tree_->Branch("isMuon1Loose",&isMuon1Loose);
	tree_->Branch("isMuon2Soft",&isMuon2Soft);
	tree_->Branch("isMuon2Tight",&isMuon2Tight);
	tree_->Branch("isMuon2PF",&isMuon2PF);
	tree_->Branch("isMuon2Loose",&isMuon2Loose);


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
