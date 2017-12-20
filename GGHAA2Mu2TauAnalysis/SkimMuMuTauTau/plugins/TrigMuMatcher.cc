// -*- C++ -*-
//
// Package:    temp/TrigMuMatcher
// Class:      TrigMuMatcher
// 
/**\class TrigMuMatcher TrigMuMatcher.cc temp/TrigMuMatcher/plugins/TrigMuMatcher.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Wed, 25 Nov 2015 16:25:51 GMT
//
//


// system include files
#include <memory>
#include <string>
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
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "Tools/Common/interface/GenTauDecayID.h"
#include "Tools/Common/interface/Common.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

//
//
// class declaration
//

class TrigMuMatcher : public edm::EDFilter {
   public:
      explicit TrigMuMatcher(const edm::ParameterSet&);
      ~TrigMuMatcher();

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
  edm::EDGetTokenT<edm::View<pat::Muon> > mu12Tag_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_; 
  edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjects_;
  std::vector<std::string> trigNames_;
  double dRCut_; 
  bool checkMatchMu1_;

  std::map<std::string, TH1D*> histos1D_;
  std::map<std::string, TH2D*> histos2D_;

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
TrigMuMatcher::TrigMuMatcher(const edm::ParameterSet& iConfig):
  mu12Tag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("mu12Tag"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  triggerObjects_(consumes<std::vector<pat::TriggerObjectStandAlone> >(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  trigNames_(iConfig.getParameter<std::vector<std::string> >("trigNames")),
  dRCut_(iConfig.getParameter<double>("dRCut")),
  checkMatchMu1_(iConfig.getParameter<bool>("checkMatchMu1")),
  histos1D_(),
  histos2D_()
{

   //now do what ever initialization is needed
   produces<std::vector<pat::Muon> >();
}


TrigMuMatcher::~TrigMuMatcher()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TrigMuMatcher::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  std::auto_ptr<std::vector<pat::Muon> > muonColl(new std::vector<pat::Muon> );

  edm::Handle<edm::View<pat::Muon> > pMu12;
  iEvent.getByToken(mu12Tag_, pMu12);
  double highestMuPt = -1;
  for(edm::View<pat::Muon>::const_iterator iMuon=pMu12->begin(); iMuon!=pMu12->end();++iMuon)
  { 
    if (iMuon->pt() > highestMuPt)
      highestMuPt = iMuon->pt();
  }


  edm::Handle<edm::TriggerResults> pTriggerBits;
  edm::Handle<pat::PackedTriggerPrescales> pTriggerPrescales; 
  edm::Handle<std::vector<pat::TriggerObjectStandAlone> > pTriggerObjects;
  iEvent.getByToken(triggerBits_, pTriggerBits);
  iEvent.getByToken(triggerObjects_, pTriggerObjects);
  iEvent.getByToken(triggerPrescales_, pTriggerPrescales);

  const edm::TriggerNames &names = iEvent.triggerNames(*pTriggerBits);
  std::vector<std::string> corrTrigNames;
  std::vector<unsigned int> corrTrigSpot;
  for (unsigned int i = 0, n = pTriggerBits->size(); i < n; ++i) 
  {
   // std::cout << "Trigger " << names.triggerName(i) << ": " << (pTriggerBits->accept(i) ? "PASS" : "fail (or not run)") << std::endl;
    for ( std::string iName : trigNames_ )
    {
      if (names.triggerName(i).find(iName) != std::string::npos )
      {
        corrTrigNames.push_back(names.triggerName(i) );
        corrTrigSpot.push_back(i);
      }//if
    }//foor iName
  }//for i



  pat::TriggerObjectStandAlone TO;
  unsigned int nPassingMuons = 0;
  bool checkPassEvent = false;
  for ( uint iobj = 0; iobj <pTriggerObjects->size(); iobj++ ) // loop through objects
  {
    TO = pTriggerObjects->at(iobj);
    TO.unpackPathNames(names); // get names
    bool checkObjMatch = false; // check that this object passes
    for (unsigned int num : corrTrigSpot) // Loop through the triggers we care about
    {
      const std::string& name = names.triggerName(num); // get name
      if (TO.hasPathName(name, true) && !checkObjMatch) // If obj passes trigger we care about, and that an object didn't pass previous triggers in vector
      {
        double matchedMuPt = -1, dRMatch = -1;
        int nMatches = 0;
        for(edm::View<pat::Muon>::const_iterator iMuon=pMu12->begin(); iMuon!=pMu12->end();++iMuon) // loop through Mu1 and mu2 object
        { 
          double dRCurr = deltaR(*iMuon, TO);
          if (dRCurr < dRCut_ && checkMatchMu1_ && iMuon->pt() == highestMuPt) // if dR bettween obj and iMuon is small enough, and we require highest pt match and that iMuon is highest pt
          {
            checkPassEvent = true;
            nMatches++;
            matchedMuPt = iMuon->pt();
            muonColl->push_back(*iMuon);
            dRMatch = dRCurr; 
          }///if dR match
     
          else if (dRCurr < dRCut_ && !checkMatchMu1_ &&  iMuon->pt() > matchedMuPt) // The last req only matters if 2 muons match the same trigger object
          {
            checkPassEvent = true;
            nMatches++;
            matchedMuPt = iMuon->pt();
            muonColl->push_back(*iMuon);
            dRMatch = dRCurr; 
          }///if dR match
        }//for iMuon
        histos1D_["nMatchesPerTrigObj"]->Fill(nMatches);
        histos1D_["dRTrigObjMatch"]->Fill(dRMatch );
        histos1D_["PtTrigObjMatch"]->Fill(matchedMuPt );
        if (checkObjMatch) // count number of objects passing trigger 
          nPassingMuons++;
      }// if TO 
      histos1D_["nTrigObjects"]->Fill(pTriggerObjects->size() );
      histos1D_["nTrigObjMatched"]->Fill(nPassingMuons );
    }//for num
  }//for iobj
  iEvent.put(muonColl);
  return (checkPassEvent);
}//
   

// ------------ method called once each job just before starting event loop  ------------
void 
TrigMuMatcher::beginJob()
{
  edm::Service<TFileService> fileService;
  histos1D_["nMatchesPerTrigObj"]=fileService->make<TH1D>("nMatchesPerTrigObj","# of Matches per Trigger Object",4,-.5,3.5);
  histos1D_["nTrigObjects"]=fileService->make<TH1D>("nTrigObjects","# of Matches",10,0,10);
  histos1D_["nTrigObjMatched"]=fileService->make<TH1D>("nTrigObjMatched","# of TrigObject Matches per event",10,-.5,9.5);
  histos1D_["dRTrigObjMatch"]=fileService->make<TH1D>("dRTrigObjMatch","Matched #DeltaR(TrigObj,#mu) ",200,-1,1);
  histos1D_["PtTrigObjMatch"]=fileService->make<TH1D>("PtTrigObjMatch","Matched TrigObj p_{T}",501, -1,500);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TrigMuMatcher::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TrigMuMatcher::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TrigMuMatcher::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TrigMuMatcher::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TrigMuMatcher::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrigMuMatcher::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrigMuMatcher);
