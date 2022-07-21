// -*- C++ -*-
//
// Package:    LowPtTauProducerGPU/LowPtTauProducer
// Class:      LowPtTauProducer
// 
/**\class FitTauVertex FitTauVertex.cc FitTauVertexGPU/FitTauVertex/plugins/FitTauVertex.cc

 Description: [one line class summary]

 Implementation:
	 [Notes on implementation]
*/
//
// Original Author:  Marc Huwiler
//         Created:  Tue, 11 Feb 2020 15:22:04 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

// Custom added resources
#include "CUDADataFormats/Track/interface/PixelTrackHeterogeneous.h"
#include "CUDADataFormats/Vertex/interface/ZVertexHeterogeneous.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TLorentzVector.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"   
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <TTree.h>
#include <TFile.h>

#include "RecoPixelVertexing/PixelTrackFitting/interface/BrokenLine.h" // karimaki translate
#include "RecoPixelVertexing/PixelTrackFitting/interface/FitResult.h" //circle fit
#include "CUDADataFormats/Track/interface/TrackSoAHeterogeneousT.h" //state at beam spot
#include <Eigen/Core> //matrices
#include <Eigen/Eigenvalues> //matrices
#include "HeterogeneousCore/CUDAUtilities/interface/eigenSoA.h"

//
// class declaration
//

//using PixelTrackHeterogeneous = HeterogeneousSoA<pixelTrack::TrackSoA>;


class FitTauVertex : public edm::stream::EDProducer<> {
   public:
	  explicit FitTauVertex(const edm::ParameterSet&);
	  ~FitTauVertex();

	  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
	  virtual void beginStream(edm::StreamID) override;
	  virtual void produce(edm::Event&, const edm::EventSetup&) override;
	  virtual void endStream() override;

	  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

	  // ----------member data ---------------------------

	  edm::Handle<PixelTrackHeterogeneous> pixeltracks; 
	  //edm::Handle< std::vector<pat::PackedCandidate> > packedpfcandidates; 
	  edm::Handle<ZVertexHeterogeneous> pixelVertices; 
	  //edm::EDGetToken consumes <pixelTrack::PixelTrackHeterogeneous>(inputTag); 

	  edm::EDGetTokenT<PixelTrackHeterogeneous> pTrkToken; 
	  edm::EDGetTokenT<ZVertexHeterogeneous> vtxToken; 

	  edm::EDGetTokenT<reco::BeamSpot> bsToken;

	  Int_t verbosity; 

	  // For Gen matching 
	  bool isMC; 
	  edm::Handle<reco::GenParticleCollection> genParticles; 
	  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken; 


	  std::vector<int> matchedCandidate; 
  	  std::vector<std::vector<TLorentzVector> > generatedTauPions;
  	  std::vector<ROOT::Math::XYZPoint> genTauDecayVertices; 


	  // extremely basic definition of vertex
	  std::vector<std::vector<uint32_t> > vertices; // We just store a vector of track indiced that are part of the vertex for each event




};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
FitTauVertex::FitTauVertex(const edm::ParameterSet& iConfig)
: pTrkToken(consumes<PixelTrackHeterogeneous>(iConfig.getParameter<edm::InputTag>("tracks"))), 
  vtxToken(consumes<ZVertexHeterogeneous>(iConfig.getParameter<edm::InputTag>("vertices"))), 
  bsToken(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotsrc"))), 
  isMC(iConfig.getParameter<bool>("isMC"))
{
	if (isMC) 
	{
		genParticlesToken = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles")); 
	}

	verbosity = 2; 
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
 
   //if you want to put into the Run
   produces<ExampleData2,InRun>();
*/
   //now do what ever other initialization is needed

  std::cout << "Analyzer running!!!" << std::endl; 

  
  // ========


  //produces<lowPtTauCollectionGPU>().setBranchAlias("hltLowPtTauGPUPixelTracks"); //produces<lowPtTauCollectionGPU>("hltLowPtTauGPUPixelTracks"); 
  
}


FitTauVertex::~FitTauVertex()
{
 
   

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
FitTauVertex::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{   
   using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   iEvent.put(std::make_unique<ExampleData2>(*pIn));
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/

   //std::unique_ptr<lowPtTauCollectionGPU> tauCollection = std::make_unique<lowPtTauCollectionGPU>(); 
   
   iEvent.getByToken(pTrkToken, pixeltracks); 

   iEvent.getByToken(vtxToken, pixelVertices); 

   std::vector<uint32_t> trackCandidates; 

   generatedTauPions.clear(); // Empty the collection for each event 
   

   const auto& vertexsoa = *(iEvent.get(vtxToken).get());

   assert(&vertexsoa);

   auto nv = vertexsoa.nvFinal; 

   if(nv > vertexsoa.MAXVTX)   return;
   //hnVertices->Fill(nv); 

   if (nv <= 0) return; 


   const auto& tracksoa = *(iEvent.get(pTrkToken)); 
   const auto *quality = tracksoa.qualityData();  

   assert(&tracksoa); 

   uint32_t maxNumTracks = tracksoa.stride(); // This is the dimension of the soas, not the actual number of tracks. We iterate over them and break the loop if nHits is 0 (index no longer used). 
  

   //beamSpot
   edm::Handle<reco::BeamSpot> beampspotHandle;
   iEvent.getByToken(bsToken, beampspotHandle);
   float x0 = 0., y0 = 0., z0 = 0., dxdz = 0., dydz = 0.;
   if (!beampspotHandle.isValid()) {
	 edm::LogWarning("SiPixelValidateVerticesFromSoA") << "No beamspot found. returning vertexes with (0,0,Z) ";
   } else {
	 const reco::BeamSpot &bs = *beampspotHandle;
	 x0 = bs.x0();
	 y0 = bs.y0(); 
	 z0 = bs.z0();
	 dxdz = bs.dxdz();
	 dydz = bs.dydz();
   }

   	vertices.clear(); // empty the collection for the new event 

   	
   	/*edm::Handle<edm:: HepMCProduct > genEvtHandle;
	event.getByLabel( "generator", genEvtHandle) ;
	const HepMC::GenEvent* Evt = genEvtHandle->GetEvent() ;
	//
	// this is an example loop over the hierarchy of vertices
	//
	for ( HepMC::GenEvent::vertex_const_iterator
	          itVtx=Evt->vertices_begin(); itVtx!=Evt->vertices_end(); ++itVtx )
	{
	      //
	      // this is an example loop over particles coming out of each vertex in the loop
	      //
	      for ( HepMC::GenVertex::particles_out_const_iterator
	              itPartOut=(*itVtx)->particles_out_const_begin();
	              itPartOut!=(*itVtx)->particles_out_const_end(); ++itPartOut )
	      {
	         itPartOut->Print(); 
	      }
	}*/


   // Gen matching 
	std::vector<Int_t> ppdgId;
	if(isMC)
	{
		iEvent.getByToken(genParticlesToken, genParticles); 

		for( unsigned p=0; p < genParticles->size(); ++p)
		{
			auto genParticle = (*genParticles)[p]; 

			if(TMath::Abs(genParticle.pdgId())!=15) continue;
			if(TMath::Abs(genParticle.status())!=2) continue;

			if (verbosity >= 3) std::cout << "\t Tau found with # of daughters = " << genParticle.numberOfDaughters() << " with mother = " << genParticle.mother(0)->pdgId() << std::endl;

			std::vector<TLorentzVector> gp;
			std::vector<Int_t> matchedTrackIdx; 
			std::vector<ROOT::Math::XYZPoint> genVertices; 
			auto tauVertex = genParticle.vertex(); 
	  //      Bool_t matched = true;
			//std::cout << "Tau vertex z: " << (*genParticles)[p].vz() << std::endl; 

			for( unsigned int d=0; d < (*genParticles)[p].numberOfDaughters(); ++d )
			{
				auto daughter = (*genParticles)[p].daughter(d);

				Int_t taupdgId = daughter->pdgId();

				if (verbosity >= 3) std::cout << "\t\t --> tau decay:" << taupdgId << std::endl;

				if(TMath::Abs(taupdgId)!=211) continue; // Does this break on a nu or a pi0 ? 

				TLorentzVector tlv_gen_pion;

				tlv_gen_pion.SetPtEtaPhiM(daughter->pt(),
					daughter->eta(),
					daughter->phi(),
					daughter->mass());

				gp.push_back(tlv_gen_pion);
				genVertices.push_back(daughter->vertex()); 
				//assert((*genParticles)[p].daughter(d)->vertex() == tauVertex); 
				if ((*genParticles)[p].daughter(d)->vertex() == tauVertex) continue; 

			}

			if(gp.size()==3)
			{
				generatedTauPions.push_back(gp); // This collection will contain the pions (tracks) that are matched to a daughter of a generated tau 
				std::cout << "Tau vertex: " << tauVertex << std::endl; 
				for (auto vertex : genVertices) 
				{
					std::cout << "pi vertex: " << vertex << std::endl; 
				}

				genTauDecayVertices.push_back(genVertices.at(0)); 

				//assert(genVertices.at(0) == genVertices.at(1)); 
			}
		}
		if (verbosity >= 2) std::cout << "\t # of gen. taus with 3prong = " << generatedTauPions.size() << std::endl;
	}



	riemannFit::Matrix3d jacobian;
	double xt, yt; // x translate, y translate, sth small ~0.2 for now but later vertex 
	double R;
	riemannFit::CircleFit karimakiCircle;
	xt = 0.2;
	yt = 0.2;
	const double bField = 3.8;
	const double b = 1/bField;
	riemannFit::LineFit karimakiLine;
	// std::ofstream file("jacobian1.csv");
	for (uint32_t iTk=0; iTk<maxNumTracks; iTk++) // put here vertexsoa.MAXTRACKS ? 
	{
		jacobian << 0,0,0,0,0,0,0,0,0;
		auto nHits = tracksoa.nHits(iTk); 
		if (nHits==0) break; // Since we are looping over the size of the soa, we need to escape at the point where the elements are no longer used. 
	  
		
		// Selection 
	  	auto pt = tracksoa.pt(iTk); 

	  	// Here this is a selection, to be replaced, but I leave it as an example of how to access the quantities 		
		if (pt < 0.5) continue; 

		if (TMath::Abs(tracksoa.eta(iTk)) > 2.3) continue; 

		if (TMath::Abs(tracksoa.charge(iTk)) != 1) continue;  // Instead of pdgId 

		if (tracksoa.chi2(iTk) > 100) continue; 

		if (tracksoa.nHits(iTk) < 3) continue;

		trackCandidates.push_back(iTk); 

		R = tracksoa.charge(iTk) / (TMath::Abs(tracksoa.charge(iTk)) * bField)   ; // F=mv^2/r=qvB => r = pt/qB
		karimakiCircle.par(0) = x0; 
		karimakiCircle.par(1) = y0; 
		karimakiCircle.par(2) = R;

	    /*!< circle covariance matrix: \n
        |cov(X0,X0)|cov(Y0,X0)|cov( R,X0)| \n
        |cov(X0,Y0)|cov(Y0,Y0)|cov( R,Y0)| \n
        |cov(X0, R)|cov(Y0, R)|cov( R, R)|
        */
		karimakiCircle.cov(0,0) = tracksoa.stateAtBS.covariance(iTk)(0); 
		karimakiCircle.cov(0,1) = tracksoa.stateAtBS.covariance(iTk)(1); 
		karimakiCircle.cov(0,2) = tracksoa.stateAtBS.covariance(iTk)(2) / b; 
		karimakiCircle.cov(1,0) = karimakiCircle.cov(0,1);
		karimakiCircle.cov(1,1) = tracksoa.stateAtBS.covariance(iTk)(5);
		karimakiCircle.cov(1,2) = tracksoa.stateAtBS.covariance(iTk)(6) / b;
		karimakiCircle.cov(2,0) = karimakiCircle.cov(0,2);
		karimakiCircle.cov(2,1) = karimakiCircle.cov(1,2);
		karimakiCircle.cov(2,2) = tracksoa.stateAtBS.covariance(iTk)(9) / (b * b);
		
		karimakiLine.par(0) = tracksoa.stateAtBS.state(iTk)(4); //cotan theta
		karimakiLine.par(1) = tracksoa.stateAtBS.state(iTk)(4); //zip

		/*!< line covariance matrix: \n
        |cov(c_t,c_t)|cov(Zip,c_t)| \n
        |cov(c_t,Zip)|cov(Zip,Zip)|
        */
		karimakiLine.cov(0,0) = tracksoa.stateAtBS.covariance(iTk)(12);//cov(cotan(theta),cotan(theta))
		karimakiLine.cov(0,1) = tracksoa.stateAtBS.covariance(iTk)(13);
		karimakiLine.cov(1,0) = karimakiLine.cov(0,1); 
		karimakiLine.cov(1,1) = tracksoa.stateAtBS.covariance(iTk)(14);
		std::cout << " *PRE* jacobian " << iTk << ":" << jacobian << std::endl;

		brokenline::translateKarimaki(karimakiCircle, xt, yt, jacobian);
		//file << jacobian << "\n\n" ;
		std::cout << "jacobian " << iTk << ":" << jacobian << std::endl;
		}
		//file.close();
	 

	// We need at least 3 tracks to form a tauh candidate 
	if (trackCandidates.size() < 3) return; 

	
	unsigned int numTrackCandidates = trackCandidates.size(); 

	Double_t piMassHypothesis = 139.57; // Make a pion hypothesis for each track  This might also get into the config file 

	// Start reconstrucing the 3 body vertex
	uint32_t trkIdxi = 0, trkIdxj = 0, trkIdxk = 0; 
	for (unsigned int i=0; i<numTrackCandidates; i++) 
	{
		for (unsigned int j=i+1; j<numTrackCandidates; j++) 
		{
			for (unsigned int k=j+1; k<numTrackCandidates; k++) 
			{
				trkIdxi = trackCandidates.at(i); 
				trkIdxj = trackCandidates.at(j); 
				trkIdxk = trackCandidates.at(k); 
				Int_t tauCharge = tracksoa.charge(trkIdxi) + tracksoa.charge(trkIdxj) + tracksoa.charge(trkIdxk); 
				//std::cout << "Tau candidate charge: " << tauCharge << std::endl; 
				if (TMath::Abs(tauCharge) != 1) continue; 

				TLorentzVector pi1, pi2, pi3; 
				//TransientVertex transientVertex; 

				pi1.SetPtEtaPhiM(tracksoa.pt(trkIdxi), tracksoa.eta(trkIdxi), tracksoa.phi(trkIdxi), piMassHypothesis); 
				pi2.SetPtEtaPhiM(tracksoa.pt(trkIdxj), tracksoa.eta(trkIdxj), tracksoa.phi(trkIdxj), piMassHypothesis); 
				pi3.SetPtEtaPhiM(tracksoa.pt(trkIdxk), tracksoa.eta(trkIdxk), tracksoa.phi(trkIdxk), piMassHypothesis); 

				TLorentzVector tau = pi1 + pi2 + pi3; 


				// Check gen matching 
				Bool_t isRight = false; 
				Bool_t oneIsRight = false; 
				Bool_t twoAreRight = false; 
				Int_t pid = -999;
				if (isMC) 
				{

					for(unsigned int genParticleIdx=0; genParticleIdx < generatedTauPions.size(); genParticleIdx++)
					{
				    
					    Bool_t isRight1 = false;
					    Bool_t isRight2 = false;
					    Bool_t isRight3 = false;
					    
					    std::vector<TLorentzVector> generatedPions = generatedTauPions[genParticleIdx];
					    
					    for(unsigned int genParticleDaughterIdx=0; genParticleDaughterIdx < generatedPions.size(); genParticleDaughterIdx++)
					    {
					      
					      if(pi1.DeltaR(generatedPions[genParticleDaughterIdx]) < 0.1) isRight1 = true;
					      else if(pi2.DeltaR(generatedPions[genParticleDaughterIdx]) < 0.1) isRight2 = true;
					      else if(pi3.DeltaR(generatedPions[genParticleDaughterIdx]) < 0.1) isRight3 = true;	// TODO: set the criterion in the config file 
					      
					    }
					    
					    isRight = isRight1 && isRight2 && isRight3;

					    twoAreRight = (isRight1 and isRight2) or (isRight2 and isRight3) or (isRight1 and isRight3); 

					    oneIsRight = isRight1 or isRight2 or isRight3; 
					    
					    if(isRight)
					    {
					      pid = ppdgId[genParticleIdx]; 
					      std::cout << "Matching track with pid: " << pid << std::endl; 
					    }
					    // Matching output 
					    if (isRight or oneIsRight or twoAreRight) 
					    {
					    	//std::cout << "Tau with: " << (isRight ? 3 : (twoAreRight ? 2 : (oneIsRight ? 1 : 0))) << " matched pion(s). " << std::endl; 
					    }
					}	
				}

				// if the three are matching we add the track indiced to our collection of vertices
				if (isRight) 
				{
					vertices.push_back({trkIdxi, trkIdxj, trkIdxk}); 
				}
			}

		}
	}
	

    std::cout << "Added vertex collection of size: " << vertices.size() << std::endl; 

    //TODO: Write out relevant information 



	//iEvent.put(std::move(tauCollection)); 


   
 
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
FitTauVertex::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
FitTauVertex::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
FitTauVertex::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
FitTauVertex::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
FitTauVertex::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
FitTauVertex::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FitTauVertex::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FitTauVertex);
