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
      edm::EDGetTokenT<edm::RefVector<std::vector<reco::Muon>>> Mu1Mu2_;
      std::map<std::string, TH1D*> histos1D_;
      std::map<std::string, TH2D*> histos2D_;
      edm::EDGetTokenT<std::vector<reco::PFMET>>  pfMETsTag_;
      std::vector<double> Mu2PtBins_;
      std::vector<double> invMassBins_;
      bool MC_;
      edm::FileInPath _fp;
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
  Mu1Mu2_(consumes<edm::RefVector<std::vector<reco::Muon>>>(iConfig.getParameter<edm::InputTag>("Mu1Mu2"))),
  histos1D_(),
  histos2D_(),
  pfMETsTag_  (consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("pfMet"))), 
  Mu2PtBins_(iConfig.getParameter<std::vector<double> >("Mu2PtBins")),
  invMassBins_(iConfig.getParameter<std::vector<double>>("invMassBins")),
  MC_(iConfig.getParameter<bool>("MC")),
  _fp(iConfig.getParameter<edm::FileInPath>("fp")),
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
   edm::Handle<edm::RefVector<std::vector<reco::Muon>>> pMu1Mu2;
   iEvent.getByToken(Mu1Mu2_, pMu1Mu2);
   
   edm::Handle<std::vector<reco::PFMET>> pMets;
   iEvent.getByToken( pfMETsTag_ ,pMets);
  
   edm::Handle<std::vector<PileupSummaryInfo> > pPU;
   if (MC_) iEvent.getByToken(PUTag_, pPU);
   //if MC do Pileup reweighting
   double pu_weight = 1.0; 
   TFile *_filePU;
   std::string FullFilePath = _fp.fullPath();
   _filePU= TFile::Open(FullFilePath.c_str());
   TH1D *puweight = (TH1D*)_filePU->Get("pileup_scale");
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
   std::vector<reco::Muon*> Mu1Mu2Ptrs;
   double invMass=0;
   reco::Muon* HighestPtMu1Mu2;
   reco::Muon* LowestPtMu1Mu2;

   invMass=((*pMu1Mu2)[0]->p4()+(*pMu1Mu2)[1]->p4()).M();
   if((*pMu1Mu2)[0]->pt()> (*pMu1Mu2)[1]->pt()){
      HighestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[0])));
      LowestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[1])));
   }
   else
   {
      HighestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[1])));
      LowestPtMu1Mu2=const_cast<reco::Muon*>(&(*((*pMu1Mu2)[0])));
   }
   double Mu2Pt=0;
   double dR=0.0;
   double dRMetMu1=0.0;
   double etaOfMu2=0;
   double etaOfMu1=0;
   dR=deltaR(*LowestPtMu1Mu2, *HighestPtMu1Mu2);
   Mu2Pt=LowestPtMu1Mu2->pt();
   etaOfMu2=LowestPtMu1Mu2->eta();
   etaOfMu1=HighestPtMu1Mu2->eta();
   histos2D_["dRVsMu2Pt"]->Fill(dR, Mu2Pt); 
   histos2D_["Mu1PtMu2Pt"]->Fill(HighestPtMu1Mu2->pt(), Mu2Pt);

   reco::PFMET Met=(*pMets)[0];
   dRMetMu1=deltaR(*HighestPtMu1Mu2,Met);
   if(MC_){
      histos1D_["dRMetMu1"]->Fill(dRMetMu1, pu_weight*EventWeight);
      histos1D_["MetPt"]->Fill(Met.pt(), pu_weight*EventWeight);
      histos1D_["Mu1Mu2Pt"]->Fill(HighestPtMu1Mu2->pt(), pu_weight*EventWeight);
      histos1D_["Mu1Mu2Pt"]->Fill(LowestPtMu1Mu2->pt(), pu_weight*EventWeight);
      histos1D_["Mu1Mu2Eta"]->Fill(HighestPtMu1Mu2->eta(), pu_weight*EventWeight);
      histos1D_["Mu1Mu2Eta"]->Fill(LowestPtMu1Mu2->pt(), pu_weight*EventWeight);
      histos1D_["Mu1Pt"]->Fill(HighestPtMu1Mu2->pt(),pu_weight*EventWeight);
      histos1D_["Mu2Pt"]->Fill(Mu2Pt,pu_weight*EventWeight);
      histos1D_["dRMu1Mu2"]->Fill(dR,pu_weight*EventWeight);
      histos1D_["dRMu1Mu2Wider"]->Fill(dR,pu_weight*EventWeight);
      histos1D_["etaOfMu1"]->Fill(etaOfMu1,pu_weight*EventWeight);
      histos1D_["etaOfMu2"]->Fill(etaOfMu2,pu_weight*EventWeight);
   }
   else{
      histos1D_["dRMetMu1"]->Fill(dRMetMu1);
      histos1D_["MetPt"]->Fill(Met.pt());
      histos1D_["Mu1Mu2Pt"]->Fill(HighestPtMu1Mu2->pt());
      histos1D_["Mu1Mu2Pt"]->Fill(LowestPtMu1Mu2->pt());
      histos1D_["Mu1Mu2Eta"]->Fill(HighestPtMu1Mu2->eta());
      histos1D_["Mu1Mu2Eta"]->Fill(LowestPtMu1Mu2->eta());
      histos1D_["Mu1Pt"]->Fill(HighestPtMu1Mu2->pt());
      histos1D_["Mu2Pt"]->Fill(Mu2Pt);
      histos1D_["dRMu1Mu2"]->Fill(dR);
      histos1D_["dRMu1Mu2Wider"]->Fill(dR);
      histos1D_["etaOfMu1"]->Fill(etaOfMu1);
      histos1D_["etaOfMu2"]->Fill(etaOfMu2);
   }
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMu1Mu2=pMu1Mu2->begin(); iMu1Mu2!=pMu1Mu2->end();++iMu1Mu2)
   {
      Mu1Mu2Ptrs.push_back(const_cast<reco::Muon*>(&(**iMu1Mu2)));
   }
   for(typename edm::RefVector<std::vector<reco::Muon>>::const_iterator iMu1Mu2=pMu1Mu2->begin();iMu1Mu2!=pMu1Mu2->end(); ++iMu1Mu2)
   {
    
      histos1D_["pt_reco"]->Fill((*iMu1Mu2)->pt());
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
