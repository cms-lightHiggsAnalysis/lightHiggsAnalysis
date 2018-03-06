// -*- C++ -*-
//
// Package:    tmp/Mu1Mu2MassFilter
// Class:      Mu1Mu2MassFilter
// 
/**\class Mu1Mu2MassFilter Mu1Mu2MassFilter.cc tmp/Mu1Mu2MassFilter/plugins/Mu1Mu2MassFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Mon, 31 Oct 2016 11:43:19 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
//
// class declaration
//

class Mu1Mu2MassFilter : public edm::stream::EDFilter<> {
   public:
      explicit Mu1Mu2MassFilter(const edm::ParameterSet&);
      ~Mu1Mu2MassFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon> > Mu1Mu2_;
      double minMass_;
      double maxMass_;
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
Mu1Mu2MassFilter::Mu1Mu2MassFilter(const edm::ParameterSet& iConfig):
  Mu1Mu2_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("Mu1Mu2"))),
  minMass_(iConfig.getParameter<double>("minMass")),
  maxMass_(iConfig.getParameter<double>("maxMass"))
{
   //now do what ever initialization is needed
    produces<std::vector<pat::Muon> >();
}


Mu1Mu2MassFilter::~Mu1Mu2MassFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
Mu1Mu2MassFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   edm::Handle<edm::View<pat::Muon> > pMu1Mu2;
   iEvent.getByToken(Mu1Mu2_, pMu1Mu2);
   double invMass=0;
   invMass=((*pMu1Mu2)[0].p4()+(*pMu1Mu2)[1].p4()).M();
   if((invMass <= minMass_ && minMass_!=-1)||(maxMass_!=-1 && invMass>=maxMass_))
      return false;
   else
   {
     std::auto_ptr<std::vector<pat::Muon> > muonColl(new std::vector<pat::Muon> );
     for (edm::View<pat::Muon>::const_iterator iMuon = pMu1Mu2->begin(); iMuon != pMu1Mu2->end(); ++iMuon)
       muonColl->push_back(*iMuon);
     iEvent.put(muonColl);
     return true;
   }//else
}//Mu1Mu2MassFilter::filter

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
Mu1Mu2MassFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
Mu1Mu2MassFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
Mu1Mu2MassFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
Mu1Mu2MassFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
Mu1Mu2MassFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
Mu1Mu2MassFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Mu1Mu2MassFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(Mu1Mu2MassFilter);
