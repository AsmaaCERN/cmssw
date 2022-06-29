// -*- C++ -*-
//
// Package:    ProdTutorial/TrackAndPointsProducer
// Class:      TrackAndPointsProducer
//
/**\class TrackAndPointsProducer TrackAndPointsProducer.cc ProdTutorial/TrackAndPointsProducer/plugins/TrackAndPointsProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Asmaa Sultan F Alsubaei
//         Created:  Thu, 23 Jun 2022 09:00:50 GMT
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
#include <vector>
#include <iostream>
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//
// class declaration
//

class TrackAndPointsProducer : public edm::stream::EDProducer<> {
public:
  explicit TrackAndPointsProducer(const edm::ParameterSet&);
  ~TrackAndPointsProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginStream(edm::StreamID) override;
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endStream() override;

  //void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //void endRun(edm::Run const&, edm::EventSetup const&) override;
  //void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::InputTag src_;
  typedef math::XYZPointD Point;
  typedef std::vector<Point> PointCollection;
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
TrackAndPointsProducer::TrackAndPointsProducer(const edm::ParameterSet& iConfig) {
  //register your products
  /* Examples
  produces<ExampleData2>();
  
  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
  */
  src_ = iConfig.getParameter<edm::InputTag>("src");
  produces<PointCollection>("innerPoint").setBranchAlias("innerPoints");
  produces<PointCollection>("outerPoint").setBranchAlias("outerPoints");
  //now do what ever other initialization is needed
}

TrackAndPointsProducer::~TrackAndPointsProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void TrackAndPointsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;
  using namespace std;
  Handle<TrackCollection> tracks;
  iEvent.getByLabel( src_, tracks );
  // create the vectors. Use auto_ptr, as these pointers will automatically
   // delete when they go out of scope, a very efficient way to reduce memory leaks.
  unique_ptr<PointCollection> innerPoints( new PointCollection );
  unique_ptr<PointCollection> outerPoints( new PointCollection );
  // and already reserve some space for the new data, to control the size
  // of your executible's memory use.
 const int size = tracks->size();
 innerPoints->reserve( size );
 outerPoints->reserve( size );
 // loop over the tracks:
 for( TrackCollection::const_iterator track = tracks->begin(); 
     track != tracks->end(); ++ track ) {
   // fill the points in the vectors
   innerPoints->push_back( track->innerPosition() );
   outerPoints->push_back( track->outerPosition() );
 }
 // and save the vectors
 iEvent.put( innerPoints );
 iEvent.put( outerPoints );

  /* This is an event example
  //Read 'ExampleData' from the Event
  ExampleData const& in = iEvent.get(inToken_);

  //Use the ExampleData to create an ExampleData2 which 
  // is put into the Event
  iEvent.put(std::make_unique<ExampleData2>(in));
  */

  /* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  SetupData& setup = iSetup.getData(setupToken_);
  */
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void TrackAndPointsProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void TrackAndPointsProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
TrackAndPointsProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
TrackAndPointsProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TrackAndPointsProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TrackAndPointsProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TrackAndPointsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackAndPointsProducer);
