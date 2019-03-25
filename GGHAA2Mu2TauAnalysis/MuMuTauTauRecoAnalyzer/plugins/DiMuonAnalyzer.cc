// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer
// Class:      DiMuonAnalyzer
// 
/**\class DiMuonAnalyzer DiMuonAnalyzer.cc GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/plugins/DiMuonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Redwan Habibullah
//         Created:  Thu, 07 Mar 2019 05:12:54 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

//#include "FWCore/Common/interface/TriggerNames.h"
//#include "DataFormats/Common/interface/TriggerResults.h"

//#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TLatex.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DiMuonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DiMuonAnalyzer(const edm::ParameterSet&);
      ~DiMuonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  
  edm::EDGetTokenT<pat::MuonCollection> MuonL_;
  edm::EDGetTokenT<pat::MuonCollection> MuonT_;
  
  






  TH1D *LeadingPt;
  TH1D *TrailingPt;
  //TH1D *TrigPt;
  TH1D *InvMass;




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
DiMuonAnalyzer::DiMuonAnalyzer(const edm::ParameterSet& iConfig):
  MuonL_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("Muons1"))),
  MuonT_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("Muons2")))
 

{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fl;
 
   LeadingPt = fl->make<TH1D>("LeadingPt" , "Pt of Leading Muon" , 150 , -0.5 , 149.5 );
   TrailingPt=fl->make<TH1D>("TrailingPt","Pt of Trailing or Sub Leading  Muon ",150,-0.5,149.5);
   //TrigPt=fl->make<TH1D>("TriggerPt","Pt of Trigger Object",100,-0.5,99.5);
   InvMass=fl->make<TH1D>("Invariant Mass","Invariant Mass of Dimuon Pair",100,-0.5,99.5);

}


DiMuonAnalyzer::~DiMuonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DiMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace pat;
  using namespace reco;
  using namespace std;
  double InvMassPair;

  Handle<pat::MuonCollection> Muons1;
  iEvent.getByToken(MuonL_,Muons1);
  Handle<pat::MuonCollection> Muons2;
  iEvent.getByToken(MuonT_,Muons2);

  for(pat::MuonCollection::const_iterator iMu = Muons1->begin() ; iMu !=Muons1->end() ; ++iMu)

    {
      LeadingPt->Fill(iMu->pt());
      LeadingPt->SetFillColor(kRed);
      LeadingPt->GetXaxis()->SetTitle("Pt(GeV)");
      LeadingPt->GetYaxis()->SetTitle("# of events");

      for(pat::MuonCollection::const_iterator imu = Muons2->begin() ; imu !=Muons2->end() ; ++imu)

	{
	  TrailingPt->Fill(imu->pt());
	  TrailingPt->SetFillColor(kBlue);
	  TrailingPt->GetXaxis()->SetTitle("Pt(GeV)");
	  TrailingPt->GetYaxis()->SetTitle("# of events");

	  //std::cout<< " Trailing Pt:  " << imu->pt() <<endl;
	  InvMassPair=abs((iMu->p4()+imu->p4()).mass());
	  InvMass->Fill(InvMassPair);
	  InvMass->SetFillColor(kGreen);
          InvMass->GetXaxis()->SetTitle("m_{#mu#mu} (GeV)");
          InvMass->GetYaxis()->SetTitle("# of events");

	}



    }

}


// ------------ method called once each job just before starting event loop  ------------
void 
DiMuonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DiMuonAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiMuonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiMuonAnalyzer);
