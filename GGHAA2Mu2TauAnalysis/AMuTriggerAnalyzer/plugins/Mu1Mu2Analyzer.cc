// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/Mu1Mu2Analyzer
// Class:      Mu1Mu2Analyzer
// 
/**\class Mu1Mu2Analyzer Mu1Mu2Analyzer.cc GGHAA2Mu2TauAnalysis/Mu1Mu2Analyzer/plugins/Mu1Mu2Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mengyao Shi
//         Created:  Mon, 21 Mar 2016 12:00:57 GMT
//
//


// system include files
#include <memory>
#include <cmath>
// user include files
#include "TH1D.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TGraph.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "Tools/Common/interface/GenTauDecayID.h"
#include "Tools/Common/interface/Common.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "DataFormats/PatCandidates/interface/MET.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
using namespace std;
using namespace edm;
using namespace reco;
using namespace trigger;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Mu1Mu2Analyzer : public edm::EDAnalyzer{
   public:
      explicit Mu1Mu2Analyzer(const edm::ParameterSet&);
      ~Mu1Mu2Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<edm::View<pat::Muon>> Mu1Mu2_;
      std::map<std::string, TH1D*> histos1D_;
      std::map<std::string, TH2D*> histos2D_;
      edm::EDGetTokenT<std::vector<pat::MET>>  pfMETsTag_;
      std::vector<double> Mu2PtBins_;
      std::vector<double> invMassBins_;
      bool MC_;
      edm::FileInPath _fp;
      edm::FileInPath _fpIDs;
      edm::FileInPath _fpISOs;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUTag_;
      float EventWeight;
      edm::EDGetTokenT<GenEventInfoProduct> generator_;
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
Mu1Mu2Analyzer::Mu1Mu2Analyzer(const edm::ParameterSet& iConfig):
  Mu1Mu2_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2"))),
  histos1D_(),
  histos2D_(),
  pfMETsTag_  (consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("pfMet"))), 
  Mu2PtBins_(iConfig.getParameter<std::vector<double> >("Mu2PtBins")),
  invMassBins_(iConfig.getParameter<std::vector<double>>("invMassBins")),
  MC_(iConfig.getParameter<bool>("MC")),
  _fp(iConfig.getParameter<edm::FileInPath>("fp")),
  _fpIDs(iConfig.getParameter<edm::FileInPath>("fpIDs")),
  _fpISOs(iConfig.getParameter<edm::FileInPath>("fpISOs")),
  PUTag_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUTag"))),
  generator_(consumes<GenEventInfoProduct>(iConfig.existsAs<edm::InputTag>("Generator") ?
                                           iConfig.getParameter<edm::InputTag>("Generator"):
                                           edm::InputTag()))
{
}


Mu1Mu2Analyzer::~Mu1Mu2Analyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void
Mu1Mu2Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   edm::Handle<edm::View<pat::Muon>> pMu1Mu2;
   iEvent.getByToken(Mu1Mu2_, pMu1Mu2);
   
   edm::Handle<std::vector<pat::MET>> pMets;
   iEvent.getByToken( pfMETsTag_ ,pMets);
  
   edm::Handle<std::vector<PileupSummaryInfo> > pPU;
   if (MC_) iEvent.getByToken(PUTag_, pPU);
   //if MC do Pileup reweighting
   double pu_weight = 1.0;
   double IDs_weight=1.0;//every muon pass through "medium ID"
   double ISOs_weight=1.0; // every muon pass through 0.25 relative isolation
   TFile *_filePU;
   std::string FullFilePath = _fp.fullPath();
   _filePU= TFile::Open(FullFilePath.c_str());
   TH1D *puweight = (TH1D*)_filePU->Get("pileup_scale");

   TFile *_fileIDs;
   TH2F *IDsWeight;
   std::string FullFilePathIDs=_fpIDs.fullPath();
   _fileIDs=TFile::Open(FullFilePathIDs.c_str());
   if (MC_){
      IDsWeight= (TH2F*)_fileIDs->Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/efficienciesMC/pt_abseta_MC");
   }
   else{
      IDsWeight= (TH2F*)_fileIDs->Get("MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta/efficienciesDATA/pt_abseta_DATA");
   }
   
   TFile *_fileISOs;
   TH2F *ISOsWeight;
   std::string FullFilePathISOs=_fpISOs.fullPath();
   _fileISOs=TFile::Open(FullFilePathISOs.c_str());
   if (MC_){
      ISOsWeight=(TH2F*)_fileISOs->Get("LooseISO_MediumID_pt_eta/efficienciesMC/pt_abseta_MC");
   }
   else{
      ISOsWeight=(TH2F*)_fileISOs->Get("LooseISO_MediumID_pt_eta/efficienciesDATA/pt_abseta_DATA");
   }

   float num_PU_vertices = -1;
   //int test=0;
   if (MC_ ) {

      if(pPU.isValid()){
         int count_pu=0;
         for(vector<PileupSummaryInfo>::const_iterator cand=pPU->begin(); cand!=pPU->end();++cand){
            //std::cout << " Pileup Information: bunchXing, nvtx: " << cand->getBunchCrossing() << " " << cand->getPU_NumInteractions() << std::endl;
	    if (cand->getBunchCrossing() == 0)
            { 
               num_PU_vertices=cand->getTrueNumInteractions();
	       //test=cand->getPU_NumInteractions();
               count_pu++;
	       break;
            }
	    //num_PU_vertices=cand->getPU_NumInteractions(); in-time,out-of-time pileup
	    //BX=cand->getBunchCrossing();
         }
         //cout<<"count"<<count_pu<<std::endl;
	 //cout<<"PU_NumInteractions"<<test<<std::endl;
	 //cout<<"NumVertices"<<num_PU_vertices<<std::endl;
      }
      histos1D_["NumVertices"]->Fill(num_PU_vertices);
    
      if (num_PU_vertices!=-1){ 
         int binx = puweight->GetXaxis()->FindBin(num_PU_vertices);
         //cout << " bin x= " << binx << " " << puweight->GetBinContent(binx) << endl;	
         pu_weight=double(puweight->GetBinContent(binx));
      }
      // if Mc, do NLO corrections
      EventWeight = 1.0;
      edm::Handle<GenEventInfoProduct> gen_ev_info;
      iEvent.getByToken(generator_, gen_ev_info);
      if(gen_ev_info.isValid())
      {
         EventWeight = gen_ev_info->weight();
      }
      //cout<<"EventWeight=="<<EventWeight<<std::endl;
 
   }

   _filePU->Close();
   delete(_filePU); 
   double invMass=0;
   pat::Muon HighestPtMu1Mu2;
   pat::Muon LowestPtMu1Mu2;
  
   reco::Candidate::LorentzVector p4DiMu;
   for(edm::View<pat::Muon>::const_iterator iMuon=pMu1Mu2->begin(); iMuon!=pMu1Mu2->end();++iMuon)
   {
      p4DiMu+=iMuon->p4();
      if (iMuon->pt()<120.0 && fabs(iMuon->eta())<2.4){
         float binx=IDsWeight->GetXaxis()->FindBin(iMuon->pt());
         float biny=IDsWeight->GetYaxis()->FindBin(fabs(iMuon->eta()));
         float binxISOs=ISOsWeight->GetXaxis()->FindBin(iMuon->pt());
         float binyISOs=ISOsWeight->GetYaxis()->FindBin(fabs(iMuon->eta()));
         IDs_weight=IDs_weight*IDsWeight->GetBinContent(binx, biny);
         ISOs_weight=ISOs_weight*ISOsWeight->GetBinContent(binxISOs, binyISOs);
      }
   }
   _fileIDs->Close();
   delete(_fileIDs);
   _fileISOs->Close();
   delete(_fileISOs);
   //cout<<"IDsWeight=="<<IDs_weight<<std::endl;
   //cout<<"ISOsWeight=="<<ISOs_weight<<std::endl;
   invMass=p4DiMu.mass();
   const auto Mu1=pMu1Mu2->at(0);
   const auto Mu2=pMu1Mu2->at(1);
   if(Mu1.pt()> Mu2.pt()){
      HighestPtMu1Mu2=Mu1;
      LowestPtMu1Mu2=Mu2;
   }
   else
   {
      HighestPtMu1Mu2=Mu2;
      LowestPtMu1Mu2=Mu1;
   }
   double Mu2Pt=0;
   double dR=0.0;
   double dRMetMu1=0.0;
   double etaOfMu2=0;
   double etaOfMu1=0;
   dR=deltaR(LowestPtMu1Mu2, HighestPtMu1Mu2);
   Mu2Pt=LowestPtMu1Mu2.pt();
   etaOfMu2=LowestPtMu1Mu2.eta();
   etaOfMu1=HighestPtMu1Mu2.eta();
   histos2D_["dRVsMu2Pt"]->Fill(dR, Mu2Pt); 
   histos2D_["Mu1PtMu2Pt"]->Fill(HighestPtMu1Mu2.pt(), Mu2Pt);

   pat::MET Met=(*pMets)[0];
   dRMetMu1=HighestPtMu1Mu2.phi()-Met.phi();
   if(MC_){
      histos1D_["dRMetMu1"]->Fill(dRMetMu1, pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["MetPt"]->Fill(Met.pt(), pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["Mu1Mu2Pt"]->Fill(HighestPtMu1Mu2.pt(), pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["Mu1Mu2Pt"]->Fill(LowestPtMu1Mu2.pt(), pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["Mu1Mu2Eta"]->Fill(HighestPtMu1Mu2.eta(), pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["Mu1Mu2Eta"]->Fill(LowestPtMu1Mu2.eta(), pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["Mu1Pt"]->Fill(HighestPtMu1Mu2.pt(),pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["Mu2Pt"]->Fill(Mu2Pt,pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["dRMu1Mu2"]->Fill(dR,pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["dRMu1Mu2Wider"]->Fill(dR,pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["etaOfMu1"]->Fill(etaOfMu1,pu_weight*IDs_weight*ISOs_weight*EventWeight);
      histos1D_["etaOfMu2"]->Fill(etaOfMu2,pu_weight*IDs_weight*ISOs_weight*EventWeight);
   }
   else{
      histos1D_["dRMetMu1"]->Fill(dRMetMu1*IDs_weight*ISOs_weight);
      histos1D_["MetPt"]->Fill(Met.pt()*IDs_weight*ISOs_weight);
      histos1D_["Mu1Mu2Pt"]->Fill(HighestPtMu1Mu2.pt()*IDs_weight*ISOs_weight);
      histos1D_["Mu1Mu2Pt"]->Fill(LowestPtMu1Mu2.pt()*IDs_weight*ISOs_weight);
      histos1D_["Mu1Mu2Eta"]->Fill(HighestPtMu1Mu2.eta()*IDs_weight*ISOs_weight);
      histos1D_["Mu1Mu2Eta"]->Fill(LowestPtMu1Mu2.eta()*IDs_weight*ISOs_weight);
      histos1D_["Mu1Pt"]->Fill(HighestPtMu1Mu2.pt()*IDs_weight*ISOs_weight);
      histos1D_["Mu2Pt"]->Fill(Mu2Pt*IDs_weight*ISOs_weight);
      histos1D_["dRMu1Mu2"]->Fill(dR*IDs_weight*ISOs_weight);
      histos1D_["dRMu1Mu2Wider"]->Fill(dR*IDs_weight*ISOs_weight);
      histos1D_["etaOfMu1"]->Fill(etaOfMu1*IDs_weight*ISOs_weight);
      histos1D_["etaOfMu2"]->Fill(etaOfMu2*IDs_weight*ISOs_weight);
   }
   for(typename edm::View<pat::Muon>::const_iterator iMu1Mu2=pMu1Mu2->begin();iMu1Mu2!=pMu1Mu2->end(); ++iMu1Mu2)
   {
    
      histos1D_["pt_reco"]->Fill(iMu1Mu2->pt());
   }
     
   if(MC_) 
      histos1D_["invMass"]->Fill(invMass,pu_weight*EventWeight);
   else
      histos1D_["invMass"]->Fill(invMass);
}


// ------------ method called once each job just before starting event loop  ------------
void 
Mu1Mu2Analyzer::beginJob()
{
  edm::Service<TFileService> fileService;
  histos1D_["pt_reco"]=fileService->make<TH1D>("pt of reco muon","pt of Mu1 Mu2 (H750a09)",100, 0.0, 100);
  histos1D_["dRMetMu1"]=fileService->make<TH1D>("dRMetMu1","dRMetMu1",100, 0.0, 5.0);
  histos1D_["MetPt"]=fileService->make<TH1D>("MetPt","MetPt", 100, 0.0, 300.0);
  histos1D_["invMass"]=fileService->make<TH1D>("invMass of Mu1 Mu2","invMass of Mu1 Mu2 (H750a09)",invMassBins_.size()-1, &invMassBins_[0]);
  histos1D_["invMass"]->Sumw2();
  histos1D_["Mu1Mu2Pt"]=fileService->make<TH1D>("Mu1Mu2Pt","Mu1Mu2Pt", 100, 0, 300);
  histos1D_["Mu1Mu2Eta"]=fileService->make<TH1D>("Mu1Mu2Eta","Mu1Mu2Eta", 100, -2.5, 2.5);
  histos2D_["Mu1PtMu2Pt"]=fileService->make<TH2D>("Mu1PtMu2Pt","Mu1Pt vs Mu2Pt", 100, 0, 500, 40, 0, 200);
  histos1D_["Mu2Pt"]=fileService->make<TH1D>("pt of Mu2", "Pt of RecoMu2",Mu2PtBins_.size()-1,&Mu2PtBins_[0]);
  histos1D_["Mu2Pt"]->Sumw2();
  histos1D_["Mu1Pt"]=fileService->make<TH1D>("pt of Mu1", "pt of Mu1", 100, 0, 200);
  histos1D_["dRMu1Mu2"]=fileService->make<TH1D>("dRMu1Mu2", "dRMu1Mu2", 50, 0, 2.0);
  histos1D_["dRMu1Mu2Wider"]=fileService->make<TH1D>("dRMu1Mu2Wider", "dRMu1Mu2Wider", 50, 0, 5.0);
  histos2D_["dRVsMu2Pt"]=fileService->make<TH2D>("dRVsMu2Pt", "dRVsMu2Pt", 50, 0, 5.0, 50, 0, 50.0);
  histos1D_["etaOfMu1"]=fileService->make<TH1D>("Eta of Mu1", "Eta of Mu1", 100, -5.0, 5.0);
  histos1D_["etaOfMu2"]=fileService->make<TH1D>("Eta of Mu2", "Eta of Mu2", 100, -5.0, 5.0);
  histos1D_["NumVertices"]=fileService->make<TH1D>("NumVertices","NumVertices", 70, 0, 70);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Mu1Mu2Analyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Mu1Mu2Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Mu1Mu2Analyzer);
