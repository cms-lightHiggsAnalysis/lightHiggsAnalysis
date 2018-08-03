
// -*- C++ -*-
//
// Package:    tester/LooseAnalyzer
// Class:      LooseAnalyzer
// 
/**\class LooseAnalyzer LooseAnalyzer.cc tester/LooseAnalyzer/plugins/LooseAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Redwan Habibullah
//         Created:  Tue, 13 Feb 2018 05:43:04 GMT
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
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "PhysicsTools/SelectorUtils/interface/CutApplicatorBase.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "TLorentzVector.h"


#include "TMath.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "PhysicsTools/SelectorUtils/interface/CutApplicatorWithEventContentBase.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TLatex.h"



using namespace edm;
using namespace std;
using namespace reco;






//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MiniLooseAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
  explicit MiniLooseAnalyzer(const edm::ParameterSet&);
  ~MiniLooseAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  float GsfEleEInverseMinusPInverse(pat::ElectronCollection::const_iterator ele);
  float dEtaInSeed(pat::ElectronCollection::const_iterator ele);
  double GsfEleEffAreaPFIsoCut(pat::ElectronCollection::const_iterator ele,const edm::Event& iEvent);
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  // ----------member data ---------------------------
  

  float effA;
  float rhos;




  edm::EDGetTokenT<pat::ElectronCollection> electronSrc_;
  edm::EDGetTokenT<double>rho_;
  
  //float effA;
  //float rhos;
  std::vector<double> Etamin = {0.0000,1.000,1.4790,2.000,2.2000,2.3000,2.4000};
  std::vector<double> Etamax ={1.0000,1.4790,2.0000,2.2000,2.3000,2.4000,5.000};
  std::vector<double> EA = {0.1703,0.1715,0.1213,0.1635,0.1937,0.2393};
  

  //Declare the histograms for the variables
  TH1D *sieiehisto;
  TH1D *HoverEhisto;
  TH1D *dPhiIn;
  TH1D *Ooemoop;
  TH1D *delEtaInSeed;
  TH1D*PFIso;
  TH1D*Pt1;
  TH1D *sieiehisto2;
  TH1D *HoverEhisto2;
  TH1D *dPhiIn2;
  TH1D *Ooemoop2;
  TH1D *delEtaInSeed2;
  TH1D*PFIso2;
  TH1D*Pt2;
  //cout << "1" << endl;
  
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
MiniLooseAnalyzer::MiniLooseAnalyzer(const edm::ParameterSet& iConfig):
  electronSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  rho_(consumes<double>(iConfig.getParameter<edm::InputTag>("Rho")))

{
   //now do what ever initialization is needed
  edm::Service<TFileService> f;
  sieiehisto = f->make<TH1D>("Sigmaietaieta" , "Sigmaietaieta_{EB}" , 100 , 0 , 0.05 );
  HoverEhisto = f->make<TH1D>("H/E" , "H/E_(EB}" , 100 , 0 , 0.4 );
  dPhiIn = f->make<TH1D>("dPhIn","#Delta_{PhiSuperClusterTrackAtVtx}EB",100,0,0.4);
  Ooemoop =f->make<TH1D>("Ooemoop","#frac{1}{E} - #frac{1}{p}_{EB}",100,-.25,0.25);
  delEtaInSeed =f->make<TH1D>("dEtaInSeed","#delta_#etainseed_{EB}",100,-.01,.01);
  PFIso =f->make<TH1D>("PFIso","EffecAreaPFIso_{EB}",100,0,0.1);
  Pt1 =f->make<TH1D>("Pt1","Pt_{EB}",200,0,200);
  sieiehisto2 = f->make<TH1D>("Sigmaietaieta2" , "Sigmaietaieta_{EE}" , 100 , 0 , 0.05 );
  HoverEhisto2 = f->make<TH1D>("H/E" , "H/E_{EE}" , 100 , 0 , 0.4 );
  dPhiIn2 = f->make<TH1D>("dPhIn2","#Delta_{PhiSuperClusterTrackAtVtx}EE",100,0,0.4);
  Ooemoop2 =f->make<TH1D>("Ooemoop2","#frac{1}{E} - #frac{1}{p}_{EE}",100,-0.25,0.25);
  delEtaInSeed2 =f->make<TH1D>("dEtaInSeed2","#delta_#etainseed_{EE}",100,-.01,.01);
  PFIso2 =f->make<TH1D>("PFIso2","EffecAreaPFIso_{EE}",100,0,0.20);
  Pt2=f->make<TH1D>("Pt2","Pt_{EE}",200,0,200);
  //cout<< "2" <<endl;
}


MiniLooseAnalyzer::~MiniLooseAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//
  float MiniLooseAnalyzer::GsfEleEInverseMinusPInverse (pat::ElectronCollection::const_iterator ele)
  {
    const float ecal_energy_inverse = 1.0/ele->ecalEnergy();
    const float eSCoverP = ele->eSuperClusterOverP();
    return std::abs(1.0 - eSCoverP)*ecal_energy_inverse;
  }


float MiniLooseAnalyzer::dEtaInSeed (pat::ElectronCollection::const_iterator ele){
  return ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull() ?
    ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
}

double MiniLooseAnalyzer::GsfEleEffAreaPFIsoCut(pat::ElectronCollection::const_iterator ele,const edm::Event& iEvent)
{
  //Compute the combined isolation with effective area correction


  for(unsigned i=0;i<EA.size();i++)
    {
      if(ele->eta() > Etamin[i] && ele->eta() <Etamax[i])
        {
          effA=EA[i];
        }
    }
  const pat::Electron::PflowIsolationVariables& pfIso  = ele->pfIsolationVariables();
  const float chad = pfIso.sumChargedHadronPt;
  const float nhad = pfIso.sumNeutralHadronEt;
  const float pho = pfIso.sumPhotonEt;
  // float  eA = effA;
  //float rho = (float)(*_rhoHandle); // std::max likes float arguments
  edm::Handle<double>_rhoHandle;
  iEvent.getByToken(rho_,_rhoHandle);
  rhos = *(_rhoHandle.product());
  float iso = chad + std::max(0.0f, nhad + pho-effA*rhos);

  // Divide by pT if the relative isolation is requested
  //if( _isRelativeIso )
  iso /= ele->pt();
  
  // Apply the cut and return the result
  return iso;
  
}















// ------------ method called for each event  ------------
void
MiniLooseAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // using namespace edm;
  Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronSrc_,electrons);
  
  for(pat::ElectronCollection::const_iterator iele = electrons->begin() ; iele !=electrons->end() ; ++iele) 
    {
      if(iele->isEB())
	{
	  // cout<< "3" <<endl;
	  sieiehisto->Fill(iele->full5x5_sigmaIetaIeta());
	   HoverEhisto->Fill(iele->hadronicOverEm());
	  dPhiIn->Fill(iele->deltaPhiSuperClusterTrackAtVtx());
	  Ooemoop->Fill((double)GsfEleEInverseMinusPInverse(iele)); 
	  delEtaInSeed->Fill((double)dEtaInSeed(iele));
	  PFIso->Fill((GsfEleEffAreaPFIsoCut(iele,iEvent)));
	  Pt1->Fill(iele->pt());
	}
      
      if(iele->isEE())
	{
	  sieiehisto2->Fill(iele->full5x5_sigmaIetaIeta());
	  HoverEhisto2->Fill(iele->hadronicOverEm());
	  dPhiIn2->Fill(iele->deltaPhiSuperClusterTrackAtVtx());
	  Ooemoop2->Fill((double)GsfEleEInverseMinusPInverse(iele));
	  delEtaInSeed2->Fill((double)dEtaInSeed(iele));
	  PFIso2->Fill((GsfEleEffAreaPFIsoCut(iele,iEvent)));
	  Pt2->Fill(iele->pt());
	  
	  
	}

      
      
      
    }
  
  /*#ifdef THIS_IS_AN_EVENT_EXAMPLE
    Handle<ExampleData> pIn;
    iEvent.getByLabel("example",pIn);
    #endif
    
    #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
    #endif*/
  
  
  

}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniLooseAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniLooseAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniLooseAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniLooseAnalyzer);
