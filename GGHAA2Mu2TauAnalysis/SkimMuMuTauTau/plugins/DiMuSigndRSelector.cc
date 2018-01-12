// -*- C++ -*-
//
// Package:   temp/DiMuSigndRSelector
// Class:     DiMuSigndRSelector
// 
/**\class DiMuSigndRSelector DiMuSigndRSelector.cc temp/DiMuSigndRSelector/plugins/DiMuSigndRSelector.cc

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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
//
// class declaration
//

class DiMuSigndRSelector : public edm::EDFilter {
  public:
     explicit DiMuSigndRSelector(const edm::ParameterSet&);
     ~DiMuSigndRSelector();

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
  edm::EDGetTokenT<edm::View<pat::Muon> > muonsTag_; 
  double dRCut_;
  bool oppositeSign_;
  bool passdR_;
  std::map<std::string, TH1D*> histos1D_;
  
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
DiMuSigndRSelector::DiMuSigndRSelector(const edm::ParameterSet& iConfig):
  mu1Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu1Tag"))),
  muonsTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonsTag"))),
  dRCut_(iConfig.getParameter<double>("dRCut")),
  oppositeSign_(iConfig.getParameter<bool>("oppositeSign")),
  passdR_(iConfig.existsAs<bool>("passdR")? iConfig.getParameter<bool>("passdR"):true),  
  histos1D_()
{
  //now do what ever initialization is needed
  produces<std::vector<pat::Muon> >();
}


DiMuSigndRSelector::~DiMuSigndRSelector()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
DiMuSigndRSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<std::vector<pat::Muon> > muonColl(new std::vector<pat::Muon> );

  edm::Handle<edm::View<pat::Muon> > pMu1;
  iEvent.getByToken(mu1Tag_, pMu1); 
  pat::Muon mu1 = pMu1->at(0); 

  edm::Handle<edm::View<pat::Muon> > pMuons;
  iEvent.getByToken(muonsTag_, pMuons); 
  
  for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end();++iMuon)
  {
    double dR = deltaR(*iMuon, mu1);
    if ( dRCut_==-1 || (dR < dRCut_ && passdR_) || (dR > dRCut_ && !passdR_))
    {
      if (  (oppositeSign_ && mu1.pdgId() == (-1)* iMuon->pdgId() )  ||  (!oppositeSign_ && mu1.pdgId()==iMuon->pdgId() )  )   
        muonColl->push_back(*iMuon);
    }//if dR
  }//for iMuon

  iEvent.put(muonColl);
  return true;
}
// ------------ method called once each job just before starting event loop  ------------
void 
DiMuSigndRSelector::beginJob()
{
  edm::Service<TFileService> fileService;
  histos1D_["WhyIFailed"]=fileService->make<TH1D>("WhyIFailed","Left: fail dR | Right: passdR fail charge",2,-.5,1.5);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiMuSigndRSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
DiMuSigndRSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
DiMuSigndRSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
DiMuSigndRSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
DiMuSigndRSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiMuSigndRSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(DiMuSigndRSelector);
