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











// (Default) user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// user include files
#include "RJPsiAnalyzers/BcTo3MuAnalyzer/src/BcTo3MuAnalyzer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"


// user include files: For kinematic fit
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"




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
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  triggerResults_Label(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
  BS_Label(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("bslabel"))),

  OnlyBest_(iConfig.getParameter<bool>("OnlyBest")),
  isMC_(iConfig.getParameter<bool>("isMC")),
  OnlyGen_(iConfig.getParameter<bool>("OnlyGen")),
  tree_(0),

  // Event information
  run(0), event(0), lumiblock(0),

  // Primary vertex
  nPrimaryVertices(0),
  primaryVertexX(0), primaryVertexY(0), primaryVertexZ(0),
  primaryVertexXError(0), primaryVertexYError(0), primaryVertexZError(0),
  primaryVertexXYError(0), primaryVertexXZError(0), primaryVertexYZError(0),

  // Bc particle
  nBc(0),
  Bc_chi2(0),
  Bc_vertexProbability(0),
  Bc_decayVertexX(0), Bc_decayVertexY(0), Bc_decayVertexZ(0),
  Bc_decayVertexXError(0), Bc_decayVertexYError(0), Bc_decayVertexZError(0),
  Bc_decayVertexXYError(0), Bc_decayVertexXZError(0), Bc_decayVertexYZError(0),

  Bc_mass(0), Bc_px(0), Bc_py(0), Bc_pz(0),

  // J/Psi particles coming from Bc
  Bc_jpsi_chi2(0),
  Bc_jpsi_vertexProbability(0),
  Bc_jpsi_mass(0), Bc_jpsi_px(0), Bc_jpsi_py(0), Bc_jpsi_pz(0),

  // Muons coming from the J/Psi
  Bc_jpsi_mu1_mass(0), Bc_jpsi_mu1_px(0), Bc_jpsi_mu1_py(0), Bc_jpsi_mu1_pz(0),
  Bc_jpsi_mu2_mass(0), Bc_jpsi_mu2_px(0), Bc_jpsi_mu2_py(0), Bc_jpsi_mu2_pz(0),
  Bc_jpsi_mu1_charge(0), Bc_jpsi_mu2_charge(0),

  // Muon coming from the Bc
  nMuons(0),
  Bc_mu_mass(0), Bc_mu_px(0), Bc_mu_py(0), Bc_mu_pz(0),
  Bc_mu_charge(0)
  

  
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

   for(View<pat::Muon>::const_iterator patMuon1 = thePATMuonHandle->begin(); patMuon1 != thePATMuonHandle->end()-1; ++patMuon1)
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

       if(globalTrackMuPositive->pt()<4.0) continue;
       if(globalTrackMuNegative->pt()<4.0) continue;

       if(!(globalTrackMuPositive->quality(reco::TrackBase::highPurity))) continue;
       if(!(globalTrackMuNegative->quality(reco::TrackBase::highPurity))) continue;
    
    
       reco::TransientTrack transientTrackMuPositive((*theBuilder).build(globalTrackMuPositive));
       reco::TransientTrack transientTrackMuNegative((*theBuilder).build(globalTrackMuNegative));

       if(!(transientTrackMuPositive.impactPointTSCP().isValid())) continue;
       if(!(transientTrackMuNegative.impactPointTSCP().isValid())) continue;


       // Trajectory state to calculate DCA for the two muons
       FreeTrajectoryState trajectoryStateMuPositive = transientTrackMuPositive.impactPointTSCP().theState();
       FreeTrajectoryState trajectoryStateMuNegative = transientTrackMuNegative.impactPointTSCP().theState();

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

       // Now include the extra muon
       for(View<pat::Muon>::const_iterator patMuon3 = thePATMuonHandle->begin(); patMuon3 != thePATMuonHandle->end(); ++patMuon3)
       {
         if(patMuon3->charge()==0) continue;
         if(patMuon1==patMuon3 || patMuon2==patMuon3) continue;
         
         reco::TrackRef globalTrackMuExtra;
         globalTrackMuExtra = patMuon3->track();

         if(globalTrackMuExtra.isNull()) continue;
         //if(globalTrackMuExtra->pt()<1.0) continue;

         if(!(globalTrackMuExtra->quality(reco::TrackBase::highPurity))) continue;
         

         reco::TransientTrack transientTrackMuExtra((*theBuilder).build(globalTrackMuExtra));
         //FreeTrajectoryState trajectoryStateMuExtra = transientTrackMuExtra.impactPointTSCP().theState();

         if(!transientTrackMuExtra.impactPointTSCP().isValid()) continue;

         // JPsi + muon invariant mass (before kinematic vertex fit)

         TLorentzVector muonExtra4V, jpsi4V;
         muonExtra4V.SetXYZM(patMuon3->px(), patMuon3->py(), patMuon->pz(), muonMass);
         auto jpsiGlobalMomentum = jpsiVertexFit->currentState().globalMomentum();
         jpsi4V.SetXYZM(jpsiGlobalMomentum.x(), jpsiGlobalMomentum.y(), jpsiGlobalMomentum.z(), jpsiVertexFit->currentState().mass());

         if((muonExtra4V + jpsi4V).M()<4.2 || (muonExtra4V + jpsi4V).M()>6.8) continue;
         

         // Now we do the kinematic fit. Constaining the JPsi mass applied to the final Bplos fit.

         vector<RefCountedKinematicParticle> vectorFitParticles;
         vectorFitParticles.push_back(particleFactory.particle(transientTrackMuPositive, muonMass, chi,ndf, muonSigma));
         vectorFitParticles.push_back(particleFactory.particle(transientTrackMuNegative, muonMass, chi,ndf, muonSigma));
         vectorFitParticles.push_back(particleFactory.particle(transientTrackMuExtra, muonMass, chi,ndf, muonSigma));
         
         MultiTrackKinematicConstraint *jpsiConstraint = new TwoTrackMassKinematicConstraint(jpsiMass);
         KinematicConstrainedVertexFitter constrainedVertexFitter;
         RefCountedKinematicTree vertexFitTree = constrainedVertexFitter.fit(vectorFitParticles,jpsiConstraint);
         if(!vertexFitTree->isValid()) continue;
         vertexFittree->movePointerToTheTop();

         RefCountedKinematicParticle bcCandidateParticle = vertexFitTree->currentParticle();
         RefCountedKinematicVertex bdDecayVertex = vertexFitTree->currentDecayVertex();

         if(!bcDecayVertex->vertexIsValid()) continue;
         if((bcCandidateParticle->currentState().mass()<5.7) || (bcCandidateParticle->currentState().mass()>6.7)) continue;
         if((bcDecayVertex->chiSquared()<0.0) || (bcDecayVertex.chiSquared()>50.0)) continue;

         double bcProb_tmp = TMath::Prob(bcDecayVertex->chiSquared(),(int)bcDecayVertex->degreesOfFreedom());

         if(bdProb_tmp < 0.01) continue;

         // Get children from final Bc fit
         vertexFitTree->movePointerTotheFirstChild();
         RefCountedKinematicParticle candidateMuPositive = vertexFitTree->currentParticle();

         vertexFitTree->movePointerToTheNextChild();
         RefCountedKinematicParticle candidateMuNegative = vertexFitTree->currentParticle();

         vertexFitTree->movePointerToTheNextChild();
         RefCountedKinematicParticle candidateMuExtra = vertexFitTree->currentParticle();

         KinematicParameters kinematicParamMuPositive = candidateMuPositive->currentState().kinematicParameters();
         KinematicParameters kinematicParamMuNegative = candidateMuNegative->currentState().kinematicParameters();
        
         GlobalVector vectorMuPositive(candidateMuPositive->curretState().globalMomentum().x(),
             candidateMuPositive->currentState().globalMomentum.y(),
             candidateMuPositive->currentState().globalMomentum.z());
         
         GlobalVector vectorMuNegative(candidateMuNegative->curretState().globalMomentum().x(),
             candidateMuNegative->currentState().globalMomentum.y(),
             candidateMuNegative->currentState().globalMomentum.z());

         KinematicParameters kinematicParamMuExtra = candidateMuExtra->currentState().kinematicParameters();

         // Filling candidates variables now.
           

       }

       




     }
   }



}


// ------------ method called once each job just before starting event loop  ------------
void
BcTo3MuAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
BcTo3MuAnalyzer::endJob()
{
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
