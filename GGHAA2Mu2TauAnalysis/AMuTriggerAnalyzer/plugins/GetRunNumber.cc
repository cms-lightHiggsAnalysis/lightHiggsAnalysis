// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer
// Class:      GetRunNumber
// 
/**\class GetRunNumber GetRunNumber.cc GGHAA2Mu2TauAnalysis/AMuTriggerAnalyzer/plugins/GetRunNumber.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Fri, 20 Oct 2017 21:44:21 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include "TH1D.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class GetRunNumber : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GetRunNumber(const edm::ParameterSet&);
      ~GetRunNumber();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      edm::RunNumber_t runNumber;
      std::map<std::string, TH1D*> histos1D_;
      edm::Service<TFileService> fileService;
      // ----------member data ---------------------------
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
GetRunNumber::GetRunNumber(const edm::ParameterSet& iConfig):
 histos1D_()
{
   //now do what ever initialization is needed
   usesResource("TFileService");

}


GetRunNumber::~GetRunNumber()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GetRunNumber::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   runNumber=iEvent.run();
   histos1D_["RunNumber"]->Fill(runNumber);


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
GetRunNumber::beginJob()
{
   histos1D_["RunNumber"]=fileService->make<TH1D>("RunNumber", "RunNumber ",828 , 277981,278808);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GetRunNumber::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GetRunNumber::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GetRunNumber);
