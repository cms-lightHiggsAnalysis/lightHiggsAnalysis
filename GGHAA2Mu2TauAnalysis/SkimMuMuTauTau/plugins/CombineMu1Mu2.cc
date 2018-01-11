// -*- C++ -*-
//
// Package:   temp/CombineMu1Mu2
// Class:     CombineMu1Mu2
// 
/**\class CombineMu1Mu2 CombineMu1Mu2.cc temp/CombineMu1Mu2/plugins/CombineMu1Mu2.cc

 Description: [one line class summary]

 Implementation:
    [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//        Created:  Wed, 25 Nov 2015 16:25:51 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Tools/Common/interface/Common.h"
//
//
// class declaration
//

class CombineMu1Mu2 : public edm::EDFilter {
  public:
     explicit CombineMu1Mu2(const edm::ParameterSet&);
     ~CombineMu1Mu2();

     static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  private:
     virtual void beginJob() override;
     virtual bool filter(edm::Event&, const edm::EventSetup&) override;
     virtual void endJob() override;
     
     //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
     //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
     //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
     //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<pat::Muon> > mu1Tag_; 
  edm::EDGetTokenT<edm::View<pat::Muon> > mu2Tag_; 
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
CombineMu1Mu2::CombineMu1Mu2(const edm::ParameterSet& iConfig):
  mu1Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu1Tag"))),
  mu2Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu2Tag")))
{
  //now do what ever initialization is needed
  produces<std::vector<pat::Muon> >();
}


CombineMu1Mu2::~CombineMu1Mu2()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
CombineMu1Mu2::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<std::vector<pat::Muon> > muonColl(new std::vector<pat::Muon> );

  edm::Handle<edm::View<pat::Muon> > pMu1;
  iEvent.getByToken(mu1Tag_, pMu1);
  pat::Muon mu1 = pMu1->at(0);

  edm::Handle<edm::View<pat::Muon> > pMu2;
  iEvent.getByToken(mu2Tag_, pMu2);
  pat::Muon mu2 = pMu2->at(0);


  muonColl->push_back(mu1);
  muonColl->push_back(mu2);
  iEvent.put(muonColl);  
  return true;
}//GetMuOn
// ------------ method called once each job just before starting event loop  ------------
void 
CombineMu1Mu2::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CombineMu1Mu2::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
CombineMu1Mu2::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
CombineMu1Mu2::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
CombineMu1Mu2::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
CombineMu1Mu2::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CombineMu1Mu2::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(CombineMu1Mu2);
