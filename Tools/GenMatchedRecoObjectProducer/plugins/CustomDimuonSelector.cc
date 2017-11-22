// -*- C++ -*-
//
// Package:    CustomDimuonSelector
// Class:      CustomDimuonSelector
// 
/**\class CustomDimuonSelector CustomDimuonSelector.cc 
   BoostedTauAnalysis/CustomDimuonSelector/src/CustomDimuonSelector.cc

 Description: create a collection of custom selected muons to put in the event, and stop 
 processing if no muons are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: CustomDimuonSelector.cc,v 1.4 2012/11/08 16:40:22 yohay Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Tools/Common/interface/Common.h"

//
// class declaration
//

class CustomDimuonSelector : public edm::EDFilter {
public:
  explicit CustomDimuonSelector(const edm::ParameterSet&);
  ~CustomDimuonSelector();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();
  virtual bool beginRun(edm::Run&, edm::EventSetup const&);
  virtual bool endRun(edm::Run&, edm::EventSetup const&);
  virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------

  //input tag for base reco muon collection
  edm::EDGetTokenT<edm::View<pat::Muon> > baseMuonTag_;

  //input tag for reco muon collection
  edm::EDGetTokenT<edm::View<pat::Muon> > muonTag_;

  //muon ID to apply
  double isoMax_;

  double isoMin_;

  edm::EDGetTokenT<reco::PFCandidateCollection> particleFlow_;

  //minimum number of objects that must be found to pass the filter
  unsigned int minNumObjsToPassFilter_;
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
CustomDimuonSelector::CustomDimuonSelector(const edm::ParameterSet& iConfig) :
  baseMuonTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("baseMuonTag"))),
  muonTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonTag"))),
  isoMax_(iConfig.getParameter<double>("isoMax")),
  isoMin_(iConfig.getParameter<double>("isoMin")),
  particleFlow_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("particleFlow"))),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{
  produces<std::vector<pat::Muon> >();
}


CustomDimuonSelector::~CustomDimuonSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool CustomDimuonSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<std::vector<pat::Muon> > muonColl(new std::vector<pat::Muon> );

  //get base muons
  edm::Handle<edm::View<pat::Muon> > pBaseMuons;
  iEvent.getByToken(baseMuonTag_, pBaseMuons);

  //get muons
  edm::Handle<edm::View<pat::Muon> > pMuons;
  iEvent.getByToken(muonTag_, pMuons);

  edm::Handle<reco::PFCandidateCollection> pPFCandidates;
  iEvent.getByToken(particleFlow_, pPFCandidates);


  std::vector<pat::Muon> muons;
  muons = Common::getIsolatedPATMuons(pMuons, pBaseMuons, pPFCandidates, isoMax_, isoMin_); 

  //fill output collection
  unsigned int nPassingMuons = 0;
  for (std::vector<pat::Muon>::const_iterator iMuon = muons.begin(); iMuon != muons.end(); 
       ++iMuon) {
      muonColl->push_back(*iMuon);
      ++nPassingMuons;

  }
  iEvent.put(muonColl);

  //if not enough muons passing cuts were found in this event, stop processing
  return (nPassingMuons >= minNumObjsToPassFilter_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
CustomDimuonSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CustomDimuonSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
CustomDimuonSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
CustomDimuonSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
CustomDimuonSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
CustomDimuonSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CustomDimuonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(CustomDimuonSelector);
