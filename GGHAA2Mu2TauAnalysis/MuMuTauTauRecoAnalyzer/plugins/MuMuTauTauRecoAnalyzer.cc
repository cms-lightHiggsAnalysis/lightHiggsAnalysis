// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer
// Class:      MuMuTauTauRecoAnalyzer
// 
/**\class MuMuTauTauRecoAnalyzer MuMuTauTauRecoAnalyzer.cc GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/plugins/MuMuTauTauRecoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Fri, 11 Mar 2016 09:17:50 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "Tools/Common/interface/Common.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TSystem.h"
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Pruner.hh"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//trigger stuff
//#include "DataFormats/HLTReco/interface/TriggerObject.h"
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "DataFormats/HLTReco/interface/TriggerEvent.h"
//#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigData.h"
//#include "FWCore/Common/interface/TriggerNames.h"
using namespace std;
using namespace edm;
using namespace reco;
using namespace fastjet;
using namespace trigger;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MuMuTauTauRecoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuMuTauTauRecoAnalyzer(const edm::ParameterSet&);
      ~MuMuTauTauRecoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      void reset(const bool);
      virtual void endJob() override; 
// ----------member data ---------------------------
  edm::EDGetTokenT<reco::PFJetCollection> corrJetTag_;
  edm::EDGetTokenT<reco::PFTauRefVector> tauTag_;
  edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> Mu1Mu2_;

  edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> IsolatedMuon_;
  edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> Mu3ID_;
  edm::EDGetTokenT<std::vector<reco::PFMET>> Met_;
  edm::EDGetTokenT<std::vector<reco::PFJet>> Jet_;
 // edm::EDGetTokenT<reco::GenParticleCollection> genParticleTag_;
  edm::EDGetTokenT<edm::ValueMap<reco::MuonRefVector>>  jetMuonMapTag_;
  std::vector<double> muHadMassBins_;
  std::vector<double> FourBInvMassBins_;
  TFile *out_;
  TH1F *muHadMass_;
  TH1F *tauHadIso_;
  TH1F *reliso_;
  TH2F *Iso2D_;
  TH1F *invMass_;
  TH1F *HighestPt_;
  TH1F *LowestPt_;
  TH1F *dR_;
  TH1F *dRTauHighestPtMuon_;
  TH1F *MissingEnergy_;
  TH1F *count_jets_;
  TH1F *count_isomu_;
  TH1F *Mu3Pt_;
  std::string outFileName_;
  edm::EDGetTokenT<reco::PFTauDiscriminator> tauHadIsoTag_;
  edm::EDGetTokenT<reco::PFCandidateCollection> particleFlow_;
  std::string inputFile_;
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
MuMuTauTauRecoAnalyzer::MuMuTauTauRecoAnalyzer(const edm::ParameterSet& iConfig):
  corrJetTag_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("corrJetTag"))),
  tauTag_(consumes<reco::PFTauRefVector>(iConfig.getParameter<edm::InputTag>("tauTag"))),
  Mu1Mu2_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2"))),
  IsolatedMuon_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("IsolatedMuon"))),
  Mu3ID_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("Mu3ID"))),
  Met_(consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("Met"))),
  Jet_(consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("Jet"))),
 // genParticleTag_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticleTag"))),
  jetMuonMapTag_(consumes<edm::ValueMap<reco::MuonRefVector> >(iConfig.getParameter<edm::InputTag>("jetMuonMapTag"))),
  muHadMassBins_(iConfig.getParameter<std::vector<double> >("muHadMassBins")),
  FourBInvMassBins_(iConfig.getParameter<std::vector<double>>("FourBInvMassBins")),
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  tauHadIsoTag_(consumes<reco::PFTauDiscriminator>(iConfig.getParameter<edm::InputTag>("tauHadIsoTag"))),
  particleFlow_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("particleFlow"))),
  inputFile_(iConfig.getParameter<std::string>("inputFile"))
{
   //now do what ever initialization is needed
   reset(false);
}


MuMuTauTauRecoAnalyzer::~MuMuTauTauRecoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuMuTauTauRecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
  edm::Handle<reco::PFJetCollection> pPFJets;
  iEvent.getByToken(corrJetTag_, pPFJets);

  edm::Handle<reco::PFTauRefVector> pTaus;
  iEvent.getByToken(tauTag_, pTaus);  
 
  edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pMu1Mu2;
  iEvent.getByToken(Mu1Mu2_, pMu1Mu2);
 
  edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pIsolatedMuon;
  iEvent.getByToken(IsolatedMuon_, pIsolatedMuon);

  edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pMu3ID;
  iEvent.getByToken(Mu3ID_, pMu3ID);
  
  edm::Handle<std::vector<reco::PFMET>> pPFMet;
  iEvent.getByToken(Met_,pPFMet);
  
  edm::Handle<std::vector<reco::PFJet>> pPFJet;
  iEvent.getByToken(Jet_,pPFJet);
  
  //edm::Handle<reco::GenParticleCollection> pGenParticles;
  // iEvent.getByToken(genParticleTag_, pGenParticles);

  edm::Handle<edm::ValueMap<reco::MuonRefVector> > pMuonJetMap;
  iEvent.getByToken(jetMuonMapTag_, pMuonJetMap);

  edm::Handle<reco::PFTauDiscriminator> pTauHadIso;
  iEvent.getByToken(tauHadIsoTag_, pTauHadIso);

  edm::Handle<reco::PFCandidateCollection> pPFCandidates;
  iEvent.getByToken(particleFlow_, pPFCandidates);
  //fill an STL container
  std::vector<reco::Muon*> Mu1Mu2Ptrs;
  double invMass=0;
  double dR=0;
  reco::Muon* HighestPtMu1Mu2;
  reco::Muon* LowestPtMu1Mu2;
  Mu3Pt_->Fill((*pMu3ID)[0]->pt());
  invMass=((*pMu1Mu2)[0]->p4()+(*pMu1Mu2)[1]->p4()).M();
  invMass_->Fill(invMass);
  if((*pMu1Mu2)[0]->pt()> (*pMu1Mu2)[1]->pt()){
    HighestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[0])));
    LowestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[1])));
    }
  else
    {
    HighestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[1])));
    LowestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[0])));
    }
  HighestPt_->Fill(HighestPtMu1Mu2->pt());
  LowestPt_->Fill(LowestPtMu1Mu2->pt());
  
  dR=deltaR(*LowestPtMu1Mu2, *HighestPtMu1Mu2);
  dR_->Fill(dR);
  std::vector<reco::PFCandidate*> PFCandidatePtrs; 
  double reliso=0;
  for (reco::MuonRefVector::const_iterator iMuon = pIsolatedMuon->begin(); iMuon != pIsolatedMuon->end(); ++iMuon) {
    reco::MuonPFIsolation iso = (*iMuon)->pfIsolationR04();
    reliso = (iso.sumChargedHadronPt+TMath::Max(0.,iso.sumNeutralHadronEt+iso.sumPhotonEt-0.5*iso.sumPUPt))/(*iMuon)->pt();
    reliso_->Fill(reliso);
  }
  double count_isomu=pIsolatedMuon->size();
  count_isomu_->Fill(count_isomu);

  for(typename std::vector<reco::PFMET>::const_iterator iPFMet=pPFMet->begin(); iPFMet!=pPFMet->end(); ++iPFMet){
    MissingEnergy_->Fill(iPFMet->energy());
  }
  int count_jets=0;

  for(typename std::vector<reco::PFJet>::const_iterator iPFJet=pPFJet->begin(); iPFJet!=pPFJet->end(); ++iPFJet){
    if(iPFJet->pt()>30){
      count_jets++;
    }
  }
  std::cout<<"count_jets="<<count_jets<<std::endl;
  count_jets_->Fill(count_jets);
 

/*  std::vector<reco::GenParticle*> genObjPtrs; 
   for(typename std::vector<reco::GenParticle>::const_iterator iGenObj=pGenParticles->begin(); iGenObj!=pGenParticles->end(); ++iGenObj)
   {
     const unsigned int absPDGID=fabs(iGenObj->pdgId());
     if(absPDGID==13)
       genObjPtrs.push_back(const_cast<reco::GenParticle*>(&(*iGenObj)));
   }
*/
  double PUWeight = 1.0;
  double HiggsPTWeight = 1.0;
   double tauHadPTWeight = 1.0;
  std::vector<reco::PFTauRef> pTSortedTaus;
  int NJets =pPFJets->size();
  std::cout<<"NJets="<<NJets<<std::endl;
//  for(reco::PFJetCollection::const_iterator iJet=pPFJets->begin(); iJet != pPFJets->end(); ++ iJet)
//    double JetPt=iJet->pt();
  for (reco::PFTauRefVector::const_iterator iTau = pTaus->begin(); iTau != pTaus->end(); 
       ++iTau) { pTSortedTaus.push_back(*iTau); }
  std::vector<reco::PFTauRef> taus = pTSortedTaus;
  std::vector<reco::PFTauRef>::const_iterator iTau = taus.begin();
  std::vector<reco::PFTauRef>::const_iterator endTau = taus.end();
  while (iTau != endTau) {
    const reco::PFJetRef& tauJetRef = (*iTau)->jetRef();
    const reco::MuonRefVector& removedMuons = (*pMuonJetMap)[tauJetRef];
    
    //make a collection of corrected old jets in |eta| < 2.4 not overlapping the W muon or tau
  std::vector<reco::MuonRef> removedMuonRefs;

   for (reco::MuonRefVector::const_iterator iMuon = removedMuons.begin(); 
	 iMuon != removedMuons.end(); ++iMuon) { removedMuonRefs.push_back(*iMuon); }
   const double muHadMass = 
      (removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).M();

   double muHadMassForPlotting = muHadMass;
   muHadMass_->Fill(muHadMassForPlotting, PUWeight*tauHadPTWeight*HiggsPTWeight);
   tauHadIso_->Fill((*pTauHadIso)[*iTau], PUWeight);
   //std::cout<<"Iso=="<<Iso<<"(*pTauHadIso)[*iTau]=="<<((*pTauHadIso)[*iTau])<<std::endl;
   double ditauiso=(*pTauHadIso)[*iTau];
   Iso2D_->Fill(HighestPtMu1Mu2->pt(), ditauiso);
   //ABCD method region defined as:
   //region A Iso<0.5, di-tau iso<10GeV,
   //region B Iso>1.0, di-tau iso<10GeV,
   //region C Iso<0.5, di-tau iso>30GeV,
   //region D Iso>1.0, di-tau iso>30GeV.
   double dRTauHighestPtMuon=0;
   dRTauHighestPtMuon=deltaR(**iTau, *HighestPtMu1Mu2);
   dRTauHighestPtMuon_->Fill(dRTauHighestPtMuon);
   
   //std::vector<reco::Muon*> Mu1Mu2Ptrs;
 // reco::Muon* HighestPtMu1Mu2;
 // reco::Muon* LowestPtMu1Mu2;
 // if((*pMu1Mu2)[0]->pt()> (*pMu1Mu2)[1]->pt()){
  //  HighestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[0])));
  //  LowestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[1])));
   // }
 // else
  //  {
   // HighestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[1])));
   // LowestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[0])));
   // }

       //const double mumuInvMass =
         //(HighestPtMu1Mu2->p4()+LowestPtMu1Mu2->p4()).M();
       //const double FourBInvMass=(HighestPtMu1Mu2->p4()+LowestPtMu1Mu2->p4()+removedMuonRefs[removedMuonRefs.size() - 1]->p4() + (*iTau)->p4()).M();
    ++iTau;
  }
}

void 
MuMuTauTauRecoAnalyzer::reset(const bool doDelete)
{
  if(doDelete && (out_ !=NULL)) delete out_;
  out_=NULL;
  if(doDelete && (muHadMass_!=NULL)) delete muHadMass_;
  muHadMass_=NULL;
  if(doDelete && (tauHadIso_!=NULL)) delete tauHadIso_;
  tauHadIso_=NULL;
  if(doDelete && (reliso_!=NULL)) delete reliso_;
  reliso_=NULL;
  if(doDelete && (Iso2D_!=NULL)) delete Iso2D_;
  Iso2D_=NULL;
  if(doDelete && (invMass_!=NULL)) delete invMass_;
  invMass_=NULL;
  if(doDelete && (HighestPt_!=NULL)) delete HighestPt_;
  HighestPt_=NULL;
  if(doDelete && (LowestPt_!=NULL)) delete LowestPt_;
  LowestPt_=NULL;
  if(doDelete && (dR_!=NULL)) delete dR_;
  dR_=NULL;
  if(doDelete && (dRTauHighestPtMuon_!=NULL)) delete dRTauHighestPtMuon_;
  dRTauHighestPtMuon_=NULL;
  if(doDelete && (MissingEnergy_!=NULL)) delete MissingEnergy_;
  MissingEnergy_=NULL;
  if(doDelete && (count_jets_!=NULL)) delete count_jets_;
  count_jets_=NULL;
  if(doDelete && (count_isomu_!=NULL)) delete count_isomu_;
  count_isomu_=NULL;
  if(doDelete && (Mu3Pt_!=NULL)) delete Mu3Pt_;
  Mu3Pt_=NULL;
}

// ------------ method called once each job just before starting event loop  ------------
void 
MuMuTauTauRecoAnalyzer::beginJob()
{
  out_= new TFile(outFileName_.c_str(),"RECREATE");
  muHadMass_=new TH1F("muHadMass", ";H125a19 m_{#mu+X} (GeV);", muHadMassBins_.size()-1, &muHadMassBins_[0]);  
  tauHadIso_=new TH1F("tauHadIso", ";H125a19 Isolation energy (GeV);", 100, 0.0, 500.0);
  reliso_=new TH1F("reliso",";H125a19 di-mu isolation",50,0.0,5.0);
  Iso2D_=new TH2F("Iso2D", "H750a9 Isolation of di-mu(X-aixs) VS di-tau", 50,0.0,500.0, 50,0.0,500.0);
  Iso2D_->SetXTitle("relative Isolation of di-mu");
  Iso2D_->SetYTitle("Tau Isolation Energy");
  invMass_=new TH1F("invMass", "H125a19 invMass di-mu(GeV)", 20, 0.0, 30.0);
  HighestPt_=new TH1F("HighestPt", "H125a19 HighestPt-Muon's Pt(GeV)", 20,0,500);
  LowestPt_=new TH1F("LowestPt","H125a19 LowestPt-Muon's Pt (GeV)", 20, 0, 400);
  dR_=new TH1F("dR(Mu1,Mu2)", "H125a19 dR(Mu1, Mu2)", 100, 0.0, 1.0);
  dRTauHighestPtMuon_=new TH1F("dR(Mu1, Tau)", "H125a19 dR(Mu1, Tau)", 20, 0.0, 5.0);
  Mu3Pt_=new TH1F("Mu3Pt","Mu3Pt",100,0.0,100.0);
  MissingEnergy_=new TH1F("MissingEnergy","MissingEnergy", 20,0.0,200.0);
  count_jets_=new TH1F("count_jets","count_jets", 20.0, 0.0, 20.0);
  count_isomu_=new TH1F("count_isomu","count_isomu", 10,0.0,10.0);

}


// ------------ method called once each job just after ending the event loop  ------------
void 
MuMuTauTauRecoAnalyzer::endJob() 
{
  out_->cd();
  TCanvas muHadMassCanvas("muHadMassCanvas","",600,600);
  Common::draw1DHistograms(muHadMassCanvas, muHadMass_);
  muHadMassCanvas.Write();
  out_->cd();
  TCanvas tauHadIsoCanvas("tauHadIsoCanvas","",600,600);
  Common::draw1DHistograms(tauHadIsoCanvas, tauHadIso_);
  tauHadIsoCanvas.Write();
  out_->cd();
  TCanvas relisoCanvas("relisoCanvas","",600,600);
  Common::draw1DHistograms(relisoCanvas,reliso_);
  relisoCanvas.Write();
  out_->cd();
  TCanvas Iso2DCanvas("iso2DCanvas","",600,600);
  Common::draw2DHistograms(Iso2DCanvas, Iso2D_);
  Iso2DCanvas.Write();
  out_->cd();
  TCanvas invMassCanvas("invMassCanvas","",600,600);
  Common::draw1DHistograms(invMassCanvas, invMass_);
  invMassCanvas.Write();
  out_->cd();
  TCanvas HighestPtCanvas("HighestPtCanvas","", 600, 600);
  Common::draw1DHistograms(HighestPtCanvas, HighestPt_);
  HighestPtCanvas.Write();
  TCanvas LowestPtCanvas("LowestPtCanvas","",600,600);
  Common::draw1DHistograms(LowestPtCanvas, LowestPt_);
  LowestPtCanvas.Write();
  out_->cd();
  TCanvas dRCanvas("dRCanvas","", 600, 600);
  Common::draw1DHistograms(dRCanvas, dR_);
  dRCanvas.Write();
  out_->cd();
  TCanvas dRTauHighestPtMuonCanvas("dRTauHighestPtMuonCanvas","", 600, 600);
  Common::draw1DHistograms(dRTauHighestPtMuonCanvas, dRTauHighestPtMuon_);
  dRTauHighestPtMuonCanvas.Write();
  out_->cd();
  TCanvas MissingEnergyCanvas("MissingEnergyCanvas","",600,600);
  Common::draw1DHistograms(MissingEnergyCanvas, MissingEnergy_); 
  MissingEnergyCanvas.Write();
  out_->cd();
  TCanvas count_jetsCanvas("count_jetsCanvas","",600,600);
  Common::draw1DHistograms(count_jetsCanvas, count_jets_);
  count_jetsCanvas.Write();
  out_->cd();
  TCanvas count_isomuCanvas("count_isomuCanvas","",600,600);
  Common::draw1DHistograms(count_isomuCanvas, count_isomu_);
  count_isomuCanvas.Write();
  out_->cd();
  TCanvas Mu3PtCanvas("Mu3PtCanvas","",600,600);
  Common::draw1DHistograms(Mu3PtCanvas,Mu3Pt_);
  Mu3PtCanvas.Write();
  out_->Close();
  std::cout<<inputFile_<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuMuTauTauRecoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuMuTauTauRecoAnalyzer);
