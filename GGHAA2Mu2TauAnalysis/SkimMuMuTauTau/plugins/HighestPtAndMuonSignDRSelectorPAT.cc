// -*- C++ -*-
//
// Package:    temp/HighestPtAndMuonSignDRSelectorPAT
// Class:      HighestPtAndMuonSignDRSelectorPAT
// 
/**\class HighestPtAndMuonSignDRSelectorPAT HighestPtAndMuonSignDRSelectorPAT.cc temp/HighestPtAndMuonSignDRSelectorPAT/plugins/HighestPtAndMuonSignDRSelectorPAT.cc

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

class HighestPtAndMuonSignDRSelectorPAT : public edm::EDFilter {
   public:
      explicit HighestPtAndMuonSignDRSelectorPAT(const edm::ParameterSet&);
      ~HighestPtAndMuonSignDRSelectorPAT();

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
  double Cut_;
  double Mu1PtCut_;
  double Mu2PtCut_;
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
HighestPtAndMuonSignDRSelectorPAT::HighestPtAndMuonSignDRSelectorPAT(const edm::ParameterSet& iConfig):
  muonTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonTag"))),
  Cut_(iConfig.getParameter<double>("dRCut")),
  Mu1PtCut_(iConfig.getParameter<double>("Mu1PtCut")),
  Mu2PtCut_(iConfig.getParameter<double>("Mu2PtCut")),
  oppositeSign_(iConfig.getParameter<bool>("oppositeSign")),
  passdR_(iConfig.existsAs<bool>("passdR")? iConfig.getParameter<bool>("passdR"):true),  
  histos1D_()
{
   //now do what ever initialization is needed
   produces<std::vector<pat::Muon> >();
}


HighestPtAndMuonSignDRSelectorPAT::~HighestPtAndMuonSignDRSelectorPAT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HighestPtAndMuonSignDRSelectorPAT::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace std;
   edm::Handle<edm::View<pat::Muon> > pMuons;
   iEvent.getByToken(muonTag_, pMuons); 

   if( (pMuons->size()) <= 1)
     return 0;
   std::auto_ptr<std::vector<pat::Muon> > muonColl(new std::vector<pat::Muon> );

   double max_=0.0;
   pat::Muon maxMuon;
   
   pat::Muon secondMuon;
   int count=0;
   for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end();++iMuon)
   {
     count+=1;
     if( (iMuon->pt()) > max_)
     {
       max_=iMuon->pt();
       maxMuon=(*iMuon);
     }
   }
   muonColl->push_back(maxMuon);
    
   //Below is selecting highest pt muons's opposite sign partner. So that's only one pair.
   int CountSecondMuon=0;
   pat::Muon tmpSecondMuon; 
   for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end();++iMuon)
   {
     //cout<<"(*iMuon)->pt() ="<<(*iMuon)->pt()<<"; deltaR(**iMuon, *maxMuon)="<<deltaR(**iMuon, *maxMuon)<<"; sign ="<<((*iMuon)->pdgId() == (1)*((maxMuon)->pdgId()))<<std::endl;
     bool PtRequireMet=false;
     bool dRRequireMet=false;
     bool signRequireMet=false;
     PtRequireMet=iMuon->pt()<(maxMuon.pt()) && iMuon->pt()> Mu2PtCut_ && (maxMuon.pt()>Mu1PtCut_);
     dRRequireMet=(passdR_ && deltaR(*iMuon, maxMuon)< Cut_) || (!passdR_ && deltaR(*iMuon, maxMuon)> Cut_) || Cut_==-1;
     signRequireMet=(oppositeSign_ && (iMuon->pdgId()==(-1)*(maxMuon.pdgId()))) || (!oppositeSign_ && (iMuon->pdgId()==maxMuon.pdgId()));
     if (PtRequireMet && dRRequireMet && signRequireMet)      
     {
       CountSecondMuon+=1; 
       if(CountSecondMuon==1){
         secondMuon=(*iMuon);
         continue;
       }
       tmpSecondMuon=(*iMuon);
       if(tmpSecondMuon.pt()> secondMuon.pt()){
         secondMuon=tmpSecondMuon;

       }
     }//if pt, dR, sign requirement
   }//for iMuon
  
   if(CountSecondMuon == 0)
   {
     return 0;
   }//if
   muonColl->push_back(secondMuon);

   //Start Debugging 
   //Below is selecting all partners.
   double invMassMostClose=0.0;
   for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end()-1; ++iMuon){
     for(edm::View<pat::Muon>::const_iterator iMuon2=iMuon+1; iMuon2!=pMuons->end(); ++iMuon2){
       double tmpInvMass=(iMuon->p4()+iMuon2->p4()).M();
       if(fabs(tmpInvMass-92.0)<fabs((invMassMostClose-92.0))){
         invMassMostClose=tmpInvMass;
       }
     }
   }
   for(edm::View<pat::Muon>::const_iterator iMuon=pMuons->begin(); iMuon!=pMuons->end()-1; ++iMuon){
     for(edm::View<pat::Muon>::const_iterator iMuon2=iMuon+1; iMuon2!=pMuons->end(); ++iMuon2){
       double invMass=0;
       bool containMaxPtMuon=false;
       bool containSecondMaxPtMuon=false;
       bool signRequireMet=false;
       signRequireMet=(oppositeSign_ && (iMuon->pdgId()==(-1)*(iMuon2->pdgId()))) || (!oppositeSign_ && (iMuon->pdgId()==iMuon2->pdgId()));
       if( !signRequireMet){
         continue;
       }
       invMass=(iMuon->p4()+iMuon2->p4()).M();
       if(fabs( iMuon->pt()-maxMuon.pt())<1e-6  ||  (fabs(maxMuon.pt()-iMuon2->pt())<1e-6)){
         containMaxPtMuon=true;
       }
       if(fabs(iMuon->pt()-secondMuon.pt())<1e-6 || (fabs(secondMuon.pt()-iMuon2->pt())<1e-6)){
         containSecondMaxPtMuon=true;
       }
       if ((!containMaxPtMuon|| !containSecondMaxPtMuon )&& fabs(invMass-invMassMostClose)<1e-6 ){
         std::cout<<"PairThatWeMissedFound!"<<std::endl;
       }
     }
   }
   iEvent.put(muonColl);
   return true;
}
// ------------ method called once each job just before starting event loop  ------------
void 
HighestPtAndMuonSignDRSelectorPAT::beginJob()
{
  edm::Service<TFileService> fileService;

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HighestPtAndMuonSignDRSelectorPAT::endJob() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
HighestPtAndMuonSignDRSelectorPAT::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
HighestPtAndMuonSignDRSelectorPAT::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
HighestPtAndMuonSignDRSelectorPAT::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
HighestPtAndMuonSignDRSelectorPAT::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HighestPtAndMuonSignDRSelectorPAT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(HighestPtAndMuonSignDRSelectorPAT);
