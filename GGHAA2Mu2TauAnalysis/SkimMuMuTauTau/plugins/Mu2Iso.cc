// -*- C++ -*-
//
// Package:   temp/Mu2Iso
// Class:     Mu2Iso
// 
/**\class Mu2Iso Mu2Iso.cc temp/Mu2Iso/plugins/Mu2Iso.cc

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

class Mu2Iso : public edm::EDFilter {
  public:
     explicit Mu2Iso(const edm::ParameterSet&);
     ~Mu2Iso();

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
  edm::EDGetTokenT<edm::View<pat::Muon> > muonTag_; 
  edm::EDGetTokenT<edm::View<pat::Muon> > mu1Tag_; 
  double relIsoCutVal_;
  bool passRelIso_;
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
Mu2Iso::Mu2Iso(const edm::ParameterSet& iConfig):
  muonTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonTag"))),
  mu1Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu1Tag"))),
  relIsoCutVal_(iConfig.getParameter<double>("relIsoCutVal")),
  passRelIso_(iConfig.getParameter<bool>("passRelIso"))
{
  //now do what ever initialization is needed
  produces<std::vector<pat::Muon> >();
}


Mu2Iso::~Mu2Iso()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
Mu2Iso::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std; 
  std::auto_ptr<std::vector<pat::Muon> > muonColl(new std::vector<pat::Muon> );

  edm::Handle<edm::View<pat::Muon> > pMuons;
  iEvent.getByToken(muonTag_, pMuons);
  if (pMuons->size() < 1)
    return 0;

  edm::Handle<edm::View<pat::Muon> > pMu1;
  iEvent.getByToken(mu1Tag_, pMu1);
  pat::Muon mu1 = pMu1->at(0);

  int CountMuon=0;
  for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end();++iMuon)
  {
    if ( deltaR(*iMuon, mu1) < .1)
      continue;
    reco::MuonPFIsolation iso = iMuon->pfIsolationR04(); 
    double reliso = (iso.sumChargedHadronPt+TMath::Max(0.,iso.sumNeutralHadronEt+iso.sumPhotonEt-0.5*iso.sumPUPt) ) / iMuon->pt();
    if ((reliso < relIsoCutVal_ && passRelIso_) || (reliso > relIsoCutVal_ && !passRelIso_) || relIsoCutVal_ == -1)
    {
      CountMuon+=1;
      muonColl->push_back(*iMuon);
    }//if isLooseMuon
  } //for iMuon

  if (CountMuon >= 1)
  {
    iEvent.put(muonColl);
    return true;
  }
  return false;
}

// ------------ method called once each job just before starting event loop  ------------
void 
Mu2Iso::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Mu2Iso::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
Mu2Iso::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
Mu2Iso::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
Mu2Iso::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
Mu2Iso::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Mu2Iso::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(Mu2Iso);
