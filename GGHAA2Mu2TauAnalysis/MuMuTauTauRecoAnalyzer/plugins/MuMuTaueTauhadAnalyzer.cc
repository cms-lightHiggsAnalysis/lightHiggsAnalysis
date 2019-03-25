// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer
// Class:      MuMuTaueTauhadAnalyzer
// 
/**\class MuMuTaueTauhadAnalyzer MuMuTaueTauhadAnalyzer.cc GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/plugins/MuMuTaueTauhadAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Redwan Habibullah
//         Created:  Tue, 12 Mar 2019 20:43:20 GMT
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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "TLorentzVector.h"

#include "TMath.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TLatex.h"
#include "TAttFill.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"

#include"DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"



#include "Math/VectorUtil.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Math/interface/deltaR.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuMuTaueTauhadAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuMuTaueTauhadAnalyzer(const edm::ParameterSet&);
      ~MuMuTaueTauhadAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<pat::ElectronCollection> Passedelectron_;
  edm::EDGetTokenT<pat::TauCollection> PassedTau_;
  
  edm::EDGetTokenT<pat::ElectronCollection> Cutelectron_;
  edm::EDGetTokenT<pat::TauCollection> CutTau_; 
  
  edm::EDGetTokenT<pat::TauCollection> Tau_;
  edm::EDGetTokenT<pat::ElectronCollection> Ele_;
  
  edm::EDGetTokenT<pat::TauCollection> Tau_pt;
  edm::EDGetTokenT<pat::ElectronCollection> Ele_pt;
  
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> >packedGenToken_;
  
  edm::EDGetTokenT<pat::TauCollection> Tau_dR;
  edm::EDGetTokenT<pat::ElectronCollection> Ele_dR;
  
  edm::EDGetTokenT<pat::MuonCollection> MuonL_;
  edm::EDGetTokenT<pat::MuonCollection> MuonT_;
  //Histograms
  TH1D *Mvisible;
  TH1D *Mvisiblecut;
  TH1D *VMassSelect;
  TH1D *PtSelect;
  TH1D *dRSelect;
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
MuMuTaueTauhadAnalyzer::MuMuTaueTauhadAnalyzer(const edm::ParameterSet& iConfig):
  Passedelectron_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  PassedTau_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("Taus"))),
  Cutelectron_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("Cutelectrons"))),
  CutTau_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("CutTaus"))),
  Tau_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("Tau_mass_select"))),
  Ele_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("Ele_mass_select"))),
  Tau_pt(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("Tau_pt_select"))),
  Ele_pt(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("Ele_pt_select"))),
  Tau_dR(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("Tau_dR_select"))),
  Ele_dR(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("Ele_dR_select"))),
  MuonL_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("Muons1"))),
  MuonT_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("Muons2")))





{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fl;
   Mvisible = fl->make<TH1D>("Mvisible" , "Mass_{visible}" , 100 ,-0.5 , 99.5 );  
   Mvisiblecut =fl->make<TH1D>("Mvisiblecut" , "Mass_{visible} after dR <0.05 and 0.8  cut" , 100 ,-0.5 , 99.5 );
   VMassSelect =fl->make<TH1D>("VMassSelect" , "Mass_{visible} from selected High mass Tau-Ele pair" , 100,-0.5 , 99.5 );
   PtSelect=fl->make<TH1D>("PtSelect","Mass_{visible} from High Pt Tau-Ele pair",100,-0.5,99.5);
   dRSelect=fl->make<TH1D>("dRSelect","Mass_{visible} from lowest dR Tau-Ele pair",100,-0.5,99.5);
   InvMass=fl->make<TH1D>("Invariant Mass","Invariant Mass of Dimuon Pair",100,-0.5,99.5);

}


MuMuTaueTauhadAnalyzer::~MuMuTaueTauhadAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuMuTaueTauhadAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   using namespace reco;
   //passed initial selection
   Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(Passedelectron_,electrons);
   Handle<pat::TauCollection> Taus;
   iEvent.getByToken(PassedTau_,Taus);
   //passed dR Cut
   Handle<pat::ElectronCollection> electrons_c;
   iEvent.getByToken(Cutelectron_,electrons_c);
   Handle<pat::TauCollection> Taus_c;
   iEvent.getByToken(CutTau_,Taus_c);
   //passed High mass selection+ CrossCleaning
   edm::Handle<pat::ElectronCollection> electron_h;
   iEvent.getByToken(Ele_,electron_h);
   edm::Handle<pat::TauCollection> tau_h;
   iEvent.getByToken(Tau_,tau_h);
   //passed High Pt selection+ crosscleaning      
   edm::Handle<pat::ElectronCollection> electron_p;
   iEvent.getByToken(Ele_pt,electron_p);
   edm::Handle<pat::TauCollection> tau_p;
   iEvent.getByToken(Tau_pt,tau_p);
   //passed High mass selection + CrossCleaning      
   edm::Handle<pat::ElectronCollection> electron_d;
   iEvent.getByToken(Ele_dR,electron_d);
   edm::Handle<pat::TauCollection> tau_d;
   iEvent.getByToken(Tau_dR,tau_d);
   //DiMuon Peak
   Handle<pat::MuonCollection> Muons1;
   iEvent.getByToken(MuonL_,Muons1);
   Handle<pat::MuonCollection> Muons2;
   iEvent.getByToken(MuonT_,Muons2);

   double Vmass_select=99999;
   double Pt_select=99999;
   double dR_select=99999;
   double InvMassPair=99999;
   double Vmass=99999;
   double VmassCut=99999;
   //Visible mass before deltaR Cut
   for(pat::ElectronCollection::const_iterator iele = electrons->begin() ; iele !=electrons->end() ; ++iele)
     {
       for(pat::TauCollection::const_iterator itau = Taus->begin() ; itau !=Taus->end() ; ++itau)
	 
	{
	  Vmass=abs((iele->p4() + itau->p4()).mass());
	  Mvisible->Fill(Vmass);
	  Mvisible->SetFillColor(kYellow);

	  Mvisible->GetXaxis()->SetTitle("M_{v} before cut ");
	  Mvisible->GetYaxis()->SetTitle("# of events");


	 }
     }
 
   for(pat::ElectronCollection::const_iterator iele = electrons_c->begin() ; iele !=electrons_c->end();++iele)
     {
       for(pat::TauCollection::const_iterator itau = Taus_c->begin() ; itau !=Taus_c->end() ; ++itau)

	 {
	   VmassCut=abs((iele->p4() + itau->p4()).mass());
	   Mvisiblecut->Fill(VmassCut);
	   Mvisiblecut->SetFillColor(kYellow);

	   Mvisiblecut->GetXaxis()->SetTitle("M_{v} after cut ");
	   Mvisiblecut->GetYaxis()->SetTitle("# of events");


         }
     }
   
   
   //Di-Muon invariant Mass
   for(pat::MuonCollection::const_iterator iMu = Muons1->begin() ; iMu !=Muons1->end() ; ++iMu)
     {
       for(pat::MuonCollection::const_iterator imu = Muons2->begin() ; imu !=Muons2->end() ; ++imu)
	 {
	   InvMassPair=abs((iMu->p4()+imu->p4()).mass());
	   InvMass->Fill(InvMassPair);
	   InvMass->SetFillColor(kGreen);
	   InvMass->GetXaxis()->SetTitle("m_{#mu#mu} (GeV)");
	   InvMass->GetYaxis()->SetTitle("# of events");
	 }
     }
   //High mass selection visible mass
   for(pat::ElectronCollection::const_iterator iele = electron_h->begin() ; iele !=electron_h->end() ; ++iele)
     {

       for(pat::TauCollection::const_iterator itau = tau_h->begin() ; itau !=tau_h->end() ; ++itau)
	 {                                                                                                                                                                                                                                                      
	   Vmass_select=abs((iele->p4() + itau->p4()).mass());
	   VMassSelect->Fill(Vmass_select);
	   VMassSelect->SetFillColor(kOrange);
	   VMassSelect->GetYaxis()->SetTitle(" # of events");
	   VMassSelect->GetXaxis()->SetTitle(" M_{v} of reconstructed High mass Tau-ele Pair");
	                                                                                                                                                                                                                                                        
	 }
     }
   // High Pt selection visible mass
   for(pat::ElectronCollection::const_iterator iele = electron_p->begin() ; iele !=electron_p->end() ; ++iele)
     {
       
       for(pat::TauCollection::const_iterator itau = tau_p->begin() ; itau !=tau_p->end() ; ++itau)

	 {  
	   Pt_select=abs((iele->p4() + itau->p4()).mass());
	   PtSelect->Fill(Pt_select);
	   PtSelect->SetFillColor(kBlue);
	   PtSelect->GetYaxis()->SetTitle(" # of events");
	   PtSelect->GetXaxis()->SetTitle(" M_{v} of reconstructed High Di-Object Pt Tau-ele Pair");

	   
	 }
       
     }
   
 

   for(pat::ElectronCollection::const_iterator iele = electron_d->begin() ; iele !=electron_d->end() ; ++iele)
     {

       for(pat::TauCollection::const_iterator itau = tau_d->begin() ; itau !=tau_d->end() ; ++itau)
	 {
	                                                                                                                                                                                                                                                           
	   dR_select=abs((iele->p4() + itau->p4()).mass());
	                                                                                                                                                                            
	   dRSelect->Fill(dR_select);
	   dRSelect->SetFillColor(kRed);
	   dRSelect->GetYaxis()->SetTitle(" # of events");
	   dRSelect->GetXaxis()->SetTitle(" M_{v} of reconstructed Tau-ele Pair with lowest dR");

					  
	 }					  
     }			  


   





}


// ------------ method called once each job just before starting event loop  ------------
void 
MuMuTaueTauhadAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuTaueTauhadAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuMuTaueTauhadAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuTaueTauhadAnalyzer);
