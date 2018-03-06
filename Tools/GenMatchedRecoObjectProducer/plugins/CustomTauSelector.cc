// -*- C++ -*-
//
// Package:    GenMatchedRecoObjectProducer
// Class:      CustomTauSelector
// 
/**\class CustomTauSelector CustomTauSelector.cc 
   BoostedTauAnalysis/GenMatchedRecoObjectProducer/src/CustomTauSelector.cc

 Description: create a collection of custom selected hadronic taus to put in the event, and stop 
 processing if no taus are selected

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay,512 1-010,+41227670495,
//         Created:  Fri Aug 24 17:10:12 CEST 2012
// $Id: CustomTauSelector.cc,v 1.6 2012/12/12 16:02:05 yohay Exp $
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

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Tools/Common/interface/Common.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//
// class declaration
//
using namespace edm;
template<class T>
class CustomTauSelector : public edm::EDFilter {
public:
  explicit CustomTauSelector(const edm::ParameterSet&);
  ~CustomTauSelector();
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

  //input tag for reco tau collection
  edm::EDGetTokenT<pat::TauCollection> tauTag_;

  //input tag for base tau collection
  edm::EDGetTokenT<pat::TauCollection> baseTauTag_;

  //input tag for tau isolation energy
  std::string tauHadIsoTag_;

  //input tag for overlap candidate collection
  edm::EDGetTokenT<edm::View<pat::Muon> > overlapMuonTag_;

  //input for checking for nearby muon
  edm::EDGetTokenT<edm::View<pat::Muon> > muons_;
  
  //vector of input tags, 1 for each discriminator the tau should pass
  std::vector<std::string> tauDiscriminatorTags_; 

  //flag indicating whether the selected taus should pass or fail the discriminator
  bool passDiscriminator_;

  //pT cut
  double pTMin_;

  //|eta| cut
  double etaMax_;

  //maximum isolation cut
  double isoMax_;

  //overlap candidate matching cut
  double dROverlapMin_;

  //minimum number of objects that must be found to pass the filter
  unsigned int minNumObjsToPassFilter_;
  
  // dR cone to look if muon was removed or nearby, indicating a di-tau object
  double diTaudRMax_;

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
template<class T>
CustomTauSelector<T>::CustomTauSelector(const edm::ParameterSet& iConfig) :
  tauTag_(iConfig.existsAs<edm::InputTag>("tauTag") ? 
	  consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("tauTag")) : edm::EDGetTokenT<pat::TauCollection>()),
  baseTauTag_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("baseTauTag"))),
  tauHadIsoTag_(iConfig.getParameter<std::string>("tauHadIsoTag")),
  overlapMuonTag_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("overlapMuonTag"))),
  muons_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"))),
  tauDiscriminatorTags_(iConfig.getParameter<std::vector<std::string> >("tauDiscriminatorTags")),
  passDiscriminator_(iConfig.getParameter<bool>("passDiscriminator")),
  pTMin_(iConfig.getParameter<double>("pTMin")),
  etaMax_(iConfig.getParameter<double>("etaMax")),
  isoMax_(iConfig.getParameter<double>("isoMax")),
  dROverlapMin_(iConfig.getParameter<double>("dROverlapMin")),
  minNumObjsToPassFilter_(iConfig.getParameter<unsigned int>("minNumObjsToPassFilter")),
  diTaudRMax_(iConfig.getParameter<double>("diTaudRMax")),
  histos1D_(),
  histos2D_()
{

  produces<std::vector<pat::Tau> >();
}

template<class T>
CustomTauSelector<T>::~CustomTauSelector()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
template<class T>
bool CustomTauSelector<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //create pointer to output collection
  std::auto_ptr<std::vector<pat::Tau> > tauColl(new std::vector<pat::Tau> );

  //get taus
  edm::Handle<pat::TauCollection> pTaus;
  if (tauTag_.isUninitialized()) {}
  else iEvent.getByToken(tauTag_, pTaus);

  //get base tau collection
  edm::Handle<pat::TauCollection> pBaseTaus;
  iEvent.getByToken(baseTauTag_, pBaseTaus);

  //get base muon collection
  edm::Handle<edm::View<pat::Muon> > pOverlapMuons;
  iEvent.getByToken(overlapMuonTag_, pOverlapMuons);

  //get base muon collection
  edm::Handle<edm::View<pat::Muon> > pMuons;
  iEvent.getByToken(muons_, pMuons);

  //fill STL container with taus passing specified discriminators in specified eta and pT range
  std::vector<pat::Tau> taus = pTaus.isValid() ? 
    Common::getPATTaus(pTaus, pBaseTaus, tauDiscriminatorTags_, tauHadIsoTag_, pTMin_, etaMax_, passDiscriminator_, isoMax_) : 
    Common::getPATTaus(pBaseTaus, tauDiscriminatorTags_, tauHadIsoTag_, pTMin_, etaMax_, passDiscriminator_, isoMax_);

  if (pOverlapMuons->size() < 1)
    return (0 >= minNumObjsToPassFilter_);

  //loop over selected taus
  unsigned int nPassingTaus = 0;
  for (std::vector<pat::Tau>::const_iterator iTau = taus.begin(); iTau != taus.end(); ++iTau) 
  {
    //fid the nearest overlap candidate to the tau
    double minDR = 1234567891123.0;
    pat::Muon nearestMuon = pOverlapMuons->at(0);
    for (size_t c = 0; c < pOverlapMuons->size(); ++c)
    {
      const auto obj = pOverlapMuons->at(c);
      const pat::Muon iOverlapCand = obj;
      double dR = reco::deltaR(*iTau, iOverlapCand);
      if (dR  < minDR) 
      {
        nearestMuon = iOverlapCand;
        minDR = dR;
      }// if dR  < minDR
    }//for iOverlapCand

    //if tau doesn't overlap with overlap candidate (or no overlap checking requested)...
    if (reco::deltaR(*iTau, nearestMuon) > dROverlapMin_) 
    {
      double bestdR = 1000000000;
      for (edm::View<pat::Muon>::const_iterator iMu = pMuons->begin(); iMu != pMuons->end(); ++iMu)
      {
        double dR = reco::deltaR(*iTau, *iMu);
        histos1D_["TauHadAllMudR"]->Fill(dR );
        if (dR < diTaudRMax_ && dR < bestdR)
          bestdR = dR;
      }//for iMu
      histos1D_["TauHadTauMudR"]->Fill(bestdR );
      histos1D_["TauHadEta"]->Fill(iTau->eta() );
      histos1D_["TauPt"]->Fill(iTau->pt());
      tauColl->push_back(*iTau);
      ++nPassingTaus;
    }//if pOverlapCands
  }
  iEvent.put(tauColl);

  //if not enough taus passing cuts were found in this event, stop processing
  histos1D_["NPassing"]->Fill((nPassingTaus >= minNumObjsToPassFilter_));
  histos2D_["NTauPassvsNTau"]->Fill(nPassingTaus, taus.size() );
  return (nPassingTaus >= minNumObjsToPassFilter_);
}

// ------------ method called once each job just before starting event loop  ------------
template<class T>
void CustomTauSelector<T>::beginJob()
{
  edm::Service<TFileService> fileService;
  histos1D_["TauPt"]=fileService->make<TH1D>("TauPt", "TauPt",100,0,300);
  histos1D_["NPassing"]=fileService->make<TH1D>("nPassing","nPassing",2,0,2);
  histos1D_["TauHadTauMudR"]=fileService->make<TH1D>("TauHadTauMudR", "TauHadTauMudR", 100, 0, 1.0);
  histos1D_["TauHadAllMudR"]=fileService->make<TH1D>("TauHadAllMudR", "TauHadAllMudR", 100, 0, 1.0);
  histos1D_["TauHadEta"]=fileService->make<TH1D>("TauHadEta", "TauHadEta", 100, 2.5, 2.5);
  histos2D_["NTauPassvsNTau"]=fileService->make<TH2D>("NTauPassvsNTau", "NTauPassvsNTau", 10, -.5, 9.5, 10, -.5, 9.5);
}

// ------------ method called once each job just after ending the event loop  ------------
template<class T>
void CustomTauSelector<T>::endJob() {
}

// ------------ method called when starting to processes a run  ------------
template<class T>
bool CustomTauSelector<T>::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}
 
// ------------ method called when ending the processing of a run  ------------
template<class T>
bool CustomTauSelector<T>::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
template<class T>
bool CustomTauSelector<T>::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
template<class T>
bool CustomTauSelector<T>::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
template<class T>
void CustomTauSelector<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
typedef CustomTauSelector<reco::Muon> CustomTauSepFromMuonSelector;
typedef CustomTauSelector<reco::Photon> CustomTauSepFromPhotonSelector;
DEFINE_FWK_MODULE(CustomTauSepFromMuonSelector);
DEFINE_FWK_MODULE(CustomTauSepFromPhotonSelector);
