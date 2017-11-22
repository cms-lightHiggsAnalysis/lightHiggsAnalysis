// -*- C++ -*-
//
// Package:    CustomMuonSelector
// Class:      CustomMuonSelector
// 
/**\class CustomMuonSelector CustomMuonSelector.cc 
   BoostedTauAnalysis/CustomMuonSelector/src/CustomMuonSelector.cc

 Description: create a collection of custom selected muons to put in the event, and stop 
 processing if no muons are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: CustomMuonSelector.cc,v 1.4 2012/11/08 16:40:22 yohay Exp $
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

class CustomMuonSelector : public edm::EDFilter {
public:
  explicit CustomMuonSelector(const edm::ParameterSet&);
  ~CustomMuonSelector();
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

  //input tag for reco vertex collection
  edm::EDGetTokenT<reco::VertexCollection> vtxTag_;

  //input tag for muons that should not be allowed to pass (i.e. they passed some other selection)
  edm::EDGetTokenT<edm::View<pat::Muon> > vetoMuonTag_;
  edm::EDGetTokenT<edm::View<pat::Muon> > vetoMuonTag2_;
  //muon ID to apply
  std::string muonID_;

  //PF isolation cut
  double PFIsoMax_;

  //detector isolation cut
  double detectorIsoMax_;

  //PU subtraction coefficient for PF isolation
  double PUSubtractionCoeff_;

  //flag indicating whether PF isolation or detector isolation should be used
  bool usePFIso_;

  //flag indicating whether the selected muons should pass the isolation cut
  bool passIso_;

  //|eta| cut
  double etaMax_;
  
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
CustomMuonSelector::CustomMuonSelector(const edm::ParameterSet& iConfig) :
  baseMuonTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("baseMuonTag"))),
  muonTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonTag"))),
  vtxTag_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vtxTag"))),
  vetoMuonTag_(iConfig.existsAs<edm::InputTag>("vetoMuonTag")?
                  consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("vetoMuonTag")):edm::EDGetTokenT<edm::View<pat::Muon> >()),
  vetoMuonTag2_(iConfig.existsAs<edm::InputTag>("vetoMuonTag2")?
                  consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("vetoMuonTag2")):edm::EDGetTokenT<edm::View<pat::Muon> >()),
  muonID_(iConfig.getParameter<std::string>("muonID")),
  PFIsoMax_(iConfig.getParameter<double>("PFIsoMax")),
  detectorIsoMax_(iConfig.getParameter<double>("detectorIsoMax")),
  PUSubtractionCoeff_(iConfig.getParameter<double>("PUSubtractionCoeff")),
  usePFIso_(iConfig.getParameter<bool>("usePFIso")),
  passIso_(iConfig.getParameter<bool>("passIso")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter"))
{
  produces<std::vector<pat::Muon> >();
}


CustomMuonSelector::~CustomMuonSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool CustomMuonSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<std::vector<pat::Muon> > muonColl(new std::vector<pat::Muon> );

  //get base muons
  edm::Handle<edm::View<pat::Muon> > pBaseMuons;
  iEvent.getByToken(baseMuonTag_, pBaseMuons);

  //get muons
  edm::Handle<edm::View<pat::Muon> > pMuons;
  iEvent.getByToken(muonTag_, pMuons);

  //get vertices
  edm::Handle<reco::VertexCollection> pVertices;
  iEvent.getByToken(vtxTag_, pVertices);

  //get veto muons
  edm::Handle<edm::View<pat::Muon> > pVetoMuons;
  if(vetoMuonTag_.isUninitialized()){}
  else iEvent.getByToken(vetoMuonTag_, pVetoMuons);

  edm::Handle<edm::View<pat::Muon> > pVetoMuons2;
  if(vetoMuonTag2_.isUninitialized()){}
  else iEvent.getByToken(vetoMuonTag2_, pVetoMuons2);

  //identify the first good vertex (the "primary" (?))
  reco::Vertex* pPV = Common::getPrimaryVertex(pVertices);

  /*fill STL container with muons passing the 2012 tight selection, PF isolation, and |eta| 
    (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId and 
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1)*/
  std::vector<pat::Muon> muons;
  if (muonID_ == "tight") 
  {
    if(usePFIso_) 
    {
      muons = pMuons.isValid() ? 
        Common::getTightPFIsolatedPATMuons(pMuons, pBaseMuons, pPV, PUSubtractionCoeff_, PFIsoMax_, etaMax_, passIso_) : 
        Common::getTightPFIsolatedPATMuons(pBaseMuons, pPV, PUSubtractionCoeff_, PFIsoMax_, etaMax_, passIso_);
    }//if usePFIso_

    /*fill STL container with muons passing the 2012 tight selection, detector isolation, and 
      |eta| (cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId and 
      https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Muon_Isolation_AN1)*/
    else 
    {
      muons = pMuons.isValid() ? 
      Common::getTightDetectorIsolatedPATMuons(pMuons, pBaseMuons, pPV, detectorIsoMax_,  etaMax_, passIso_) : 
      Common::getTightDetectorIsolatedPATMuons(pBaseMuons, pPV, detectorIsoMax_, etaMax_, passIso_);
    }//else
  }//if muonID_ == tihgt
  else if(muonID_== "tightNew")
  {
    muons=pMuons.isValid()?
      Common::getTightDetectorPATMuons(pMuons, pBaseMuons, pPV, etaMax_):
      Common::getTightDetectorPATMuons(pBaseMuons, pPV, etaMax_);
  }
  else if(muonID_== "loose")
  {
    muons=pMuons.isValid()?
      Common::getLooseDetectorPATMuons(pMuons,pBaseMuons, etaMax_):
      Common::getLooseDetectorPATMuons(pBaseMuons, etaMax_);
  }
  

  /*fill STL container of soft muons (cf. 
    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Soft_Muon)*/
  else if (muonID_ == "soft"){
     muons = pMuons.isValid() ? 
      Common::getSoftPATMuons(pMuons, pBaseMuons, pPV, etaMax_) : 
      Common::getSoftPATMuons(pBaseMuons, pPV, etaMax_);
  }
  
  else if (muonID_ == "tracker"){
     muons =pMuons.isValid() ?
       Common::getTrackerPATMuons(pMuons, pBaseMuons):
       Common::getTrackerPATMuons(pBaseMuons);
  }
  //error: unsupported muon ID
  else throw cms::Exception("CustomMuonSelector") << "Error: unsupported muon ID.\n";

  unsigned int nPassingMuons = 0;
  for (std::vector<pat::Muon>::const_iterator iMuon = muons.begin(); iMuon != muons.end(); ++iMuon) 
  {
    bool checkVeto = false;
    if (pVetoMuons.isValid()) 
    {
      for (edm::View<pat::Muon>::const_iterator iVetoMuon = pVetoMuons->begin(); iVetoMuon != pVetoMuons->end(); ++iVetoMuon) 
      {
        if (deltaR(*iMuon, *iVetoMuon) < .0001 && ( fabs(iMuon->pt()-iVetoMuon->pt()) / iVetoMuon->pt() ) < 0.0001)
          checkVeto = true;
      }//for iVetoMuon
    }//if pVetoMuons.isValid
    if (pVetoMuons2.isValid() && !checkVeto)
    {
      for (edm::View<pat::Muon>::const_iterator iVetoMuon = pVetoMuons2->begin(); iVetoMuon != pVetoMuons2->end(); ++iVetoMuon)
      {
        if (deltaR(*iMuon, *iVetoMuon) < .0001 && ( fabs(iMuon->pt()-iVetoMuon->pt()) / iVetoMuon->pt() ) < 0.0001)
          checkVeto = true;
      }//for iVetoMuon
    }//if pVetoMuons2.isValid
    if (!checkVeto)
    { 
      muonColl->push_back(*iMuon);
      ++nPassingMuons;
    }//if
  }//for iMuon

//  //make an STL container of the veto muon ref keys
//  std::vector<int> vetoMuonRefKeys;
//  if (pVetoMuons.isValid()) 
//  {
//    for (edm::View<pat::Muon>::const_iterator iVetoMuon = pVetoMuons->begin(); iVetoMuon != pVetoMuons->end(); ++iVetoMuon) 
//      vetoMuonRefKeys.push_back(iVetoMuon.key());
//  }// ifpVetoMuons.isValid
//
//  std::vector<int> vetoMuonRefKeys2;
//  if (pVetoMuons2.isValid()) 
//  {
//    for (edm::View<pat::Muon>::const_iterator iVetoMuon2 = pVetoMuons2->begin(); iVetoMuon2 != pVetoMuons2->end(); ++iVetoMuon2)
//      vetoMuonRefKeys2.push_back(iVetoMuon2.key());
//  }//pVetoMuons2.isValid
//  //fill output collection
//  unsigned int nPassingMuons = 0;
//  for (std::vector<pat::Muon>::const_iterator iMuon = muons.begin(); iMuon != muons.end(); ++iMuon) 
//  {
//    if ((std::find(vetoMuonRefKeys.begin(), vetoMuonRefKeys.end(), iMuon.key()) == vetoMuonRefKeys.end())
//        && (std::find(vetoMuonRefKeys2.begin(), vetoMuonRefKeys2.end(),  iMuon.key()) == vetoMuonRefKeys2.end())) 
//    {
//      muonColl->push_back(*iMuon);
//      ++nPassingMuons;
//    }//if
//  }//for iMuon
  iEvent.put(muonColl);

  //if not enough muons passing cuts were found in this event, stop processing
  return (nPassingMuons >= minNumObjsToPassFilter_);
}

// ------------ method called once each job just before starting event loop  ------------
void 
CustomMuonSelector::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CustomMuonSelector::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
CustomMuonSelector::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
CustomMuonSelector::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
CustomMuonSelector::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
CustomMuonSelector::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CustomMuonSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(CustomMuonSelector);
