// -*- C++ -*-  
//
// Package:    GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer
// Class:      PrimaryAnalyzer
// 
/**\class PrimaryAnalyzer PrimaryAnalyzer.cc GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/plugins/PrimaryAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Redwan Habibullah
//         Created:  Wed, 18 Apr 2018 20:01:49 GMT
//
//


// system include files
#include <memory>
#include <algorithm>
#include <map>

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




typedef math::XYZPoint Point;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class PrimaryAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit PrimaryAnalyzer(const edm::ParameterSet&);
      ~PrimaryAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  // restrictions for the electron to be
  // considered
  // double minPt_;
  //double maxEta_;
 

  //Load input collection of Rho's and Taus
  edm::EDGetTokenT<pat::ElectronCollection> electronSrc_;
  edm::EDGetTokenT<pat::TauCollection> TauSrc_;
  edm::EDGetTokenT<reco::VertexCollection> vtx_;
  edm::EDGetTokenT<reco::BeamSpot> thebs_;
  edm::EDGetTokenT<reco::GenParticleCollection> particleSrcToken_;
  edm::EDGetTokenT<pat::TauCollection> Tau_;
  edm::EDGetTokenT<pat::ElectronCollection> Ele_;
  edm::EDGetTokenT<pat::TauCollection> Tau_pt;
  edm::EDGetTokenT<pat::ElectronCollection> Ele_pt;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> >packedGenToken_;
  edm::EDGetTokenT<pat::TauCollection> Tau_dR;
  edm::EDGetTokenT<pat::ElectronCollection> Ele_dR;
  edm::EDGetTokenT<pat::MuonCollection> MuonL_;
  edm::EDGetTokenT<pat::MuonCollection> MuonT_;

 
 // internal variables for gen matching
  double maxDeltaR_;


//Initialize Histograms
  TH1D *DeltaR;
  TH1D *DeltaMatch;
  TH1D *Mvisible;
  TH1F *DeltaRele;
  TH2F *Profile;
  TH1D *DeltaRcut;
  TH1D *Mvisiblecut;
  TH1D *VMassSelect;
  TH1D *PtSelect;
  TH1D *dRSelect;
  TH1D *Pt;
  TH2F *Profile2;
  TH2F *Profile3;
  TH2F *Profile4;
  TH1D *InvMass;



  //int match =0;
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
PrimaryAnalyzer::PrimaryAnalyzer(const edm::ParameterSet& iConfig):
  //declare a tag
  electronSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  TauSrc_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("Taus"))),
  vtx_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
  thebs_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BM"))),
  //minPt_ (iConfig.getParameter<double>("minPt")),
  //maxEta_ (iConfig.getParameter<double>("maxEta")),
  //maxDeltaR_ (iConfig.getParameter<double>("maxDeltaR")),
  particleSrcToken_ (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("particleSrc"))),
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
   DeltaR = fl->make<TH1D>("DeltaR" , "#Delta_{R}" , 1000 , -0.005 , 4.995 );
   DeltaMatch =fl->make<TH1D>("DeltaMatch","#Delta_{R} High Pair ",1000,-0.05,0.95); 
   DeltaRcut = fl->make<TH1D>("DeltaRcut" , "#Delta_{R} after cut" , 1000 , -0.005 , 4.995 );
   Mvisible = fl->make<TH1D>("Mvisible" , "Mass_{visible}" , 100 ,-0.5 , 99.5 );
   //Mvisible = fl->make<TH1D>("Mvisible" , "Mass_{visible}" , 30 ,0 ,2.0 );

   Mvisiblecut =fl->make<TH1D>("Mvisiblecut" , "Mass_{visible}_after cut" , 100 ,-0.5 , 99.5 );
   VMassSelect =fl->make<TH1D>("VMassSelect" , "Mass_{visible} from selected High mass Tau-Ele pair" , 100,-0.5 , 99.5 );
   //DeltaRele = fl->make<TH1F>("DeltaRele","DeltaR_{gen/ele}",100,-.005,0.045);
   Profile =fl->make<TH2F>("Comparison","Plot of Mass_{visible} vs #DeltaR",1000,-.005,4.995,100,-0.5,499.5);
   PtSelect=fl->make<TH1D>("PtSelect","Mass_{visible} from High Pt Tau-Ele pair",100,-0.5,99.5);
   dRSelect=fl->make<TH1D>("dRSelect","Mass_{visible} from lowest dR Tau-Ele pair",100,-0.5,99.5); 
   Pt=fl->make<TH1D>("Pt","Di-Object_{Pt} Ele-Tau",500,-0.5,499.5);
   InvMass=fl->make<TH1D>("Invariant Mass","Invariant Mass of Dimuon Pair",100,-0.5,99.5);
   Profile2 =fl->make<TH2F>("Comparison2","Plot of Di-Object_{Pt} vs #DeltaR",1000,-.005,4.995,100,-0.5,499.5);
   Profile3 =fl->make<TH2F>("Comparison3","Plot of Di-Object_{Pt} vs Mass_{visibe}",100,-0.5,499.5,100,-0.5,499.5);
   Profile4 =fl->make<TH2F>("Comparison4","Plot of Di-Object_{Pt} vs Mass_{visibe} post #DletaR cut",100,-0.5,499.5,100,-0.5,499.5);
   

}


PrimaryAnalyzer::~PrimaryAnalyzer()
{
  
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
PrimaryAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronSrc_,electrons);
  Handle<pat::TauCollection> Taus;
  iEvent.getByToken(TauSrc_,Taus);
  Handle<reco::VertexCollection> Vertex;
  iEvent.getByToken(vtx_,Vertex);
  edm::Handle<reco::BeamSpot> thebs;
  iEvent.getByToken(thebs_,thebs);
  edm::Handle<reco::GenParticleCollection> particles;
  iEvent.getByToken(particleSrcToken_, particles);
  edm::Handle<pat::ElectronCollection> electron_h;
  iEvent.getByToken(Ele_,electron_h);
  edm::Handle<pat::TauCollection> tau_h;
  iEvent.getByToken(Tau_,tau_h);
  edm::Handle<pat::ElectronCollection> electron_p;
  iEvent.getByToken(Ele_pt,electron_p);
  edm::Handle<pat::TauCollection> tau_p;
  iEvent.getByToken(Tau_pt,tau_p);
  edm::Handle<pat::ElectronCollection> electron_d;
  iEvent.getByToken(Ele_dR,electron_d);
  edm::Handle<pat::TauCollection> tau_d;
  iEvent.getByToken(Tau_dR,tau_d);
  Handle<pat::MuonCollection> Muons1;
  iEvent.getByToken(MuonL_,Muons1);
  
  Handle<pat::MuonCollection> Muons2;
  iEvent.getByToken(MuonT_,Muons2);



  
  const reco::Vertex& pv=*Vertex->begin();
  Point p;
  p=pv.position();
  double dR=999999;
  double Vmass=999999;
  double Vmass_select=999999;
  double Pt_select=999999;
  double dR_select=999999;
  double Pt_di=999999;
  double Vmax;
  double dz_tau;
  double dxy_tau;
  double dz_ele;
  double dxy_ele;
  int nDoublect= 0;
  double InvMassPair=99999;

  vector <unsigned int> TausAlreadyInPlot;
  vector <unsigned int> ElectronsAlreadyInPlot;
  map<int,int> TauElemap;
  // int tau_idx=0;
  int ele_idx =0;
  //double minPhi= 0.1;
  //double mindR=0.05;
  


  for(pat::MuonCollection::const_iterator iMu = Muons1->begin() ; iMu !=Muons1->end() ; ++iMu)

    {
     

      for(pat::MuonCollection::const_iterator imu = Muons2->begin() ; imu !=Muons2->end() ; ++imu)

        {
         

          //std::cout<< " Trailing Pt:  " << imu->pt() <<endl;                                                                                                                                                                                                                 
          InvMassPair=abs((iMu->p4()+imu->p4()).mass());
          InvMass->Fill(InvMassPair);
          InvMass->SetFillColor(kGreen);
          InvMass->GetXaxis()->SetTitle("m_{#mu#mu} (GeV)");
          InvMass->GetYaxis()->SetTitle("# of events");
	  
        }



    }


  for(pat::ElectronCollection::const_iterator iele = electron_h->begin() ; iele !=electron_h->end() ; ++iele)
    {
      
	      for(pat::TauCollection::const_iterator itau = tau_h->begin() ; itau !=tau_h->end() ; ++itau)
		{
		  //double deltaphi=deltaPhi(iele->phi(),itau->phi());
		  double delta_R=deltaR(*iele,*itau);
		  DeltaMatch->Fill(delta_R);                                                                                                                                                                                                                                     
		  DeltaMatch->SetFillColor(kSpring);                                                                                                                                                                                                                                 
		  DeltaMatch->GetXaxis()->SetTitle("#Delta_{R} High Mass Pair");                                                                                                                                                                                                               
		  DeltaMatch->GetYaxis()->SetTitle("# of events"); 
		  //if (delta_R > mindR)
		  //{
		  Vmass_select=abs((iele->p4() + itau->p4()).mass());
		  //cout<<"Vmass: "<<Vmass_select<<endl;
		  //cout<<"Vmass Value: "<<Vmass_select<<endl;
		  VMassSelect->Fill(Vmass_select);
		  VMassSelect->SetFillColor(kOrange);
		  VMassSelect->GetYaxis()->SetTitle(" # of events");
		  VMassSelect->GetXaxis()->SetTitle(" M_{v} of reconstructed High mass Tau-ele Pair");
		  //  }
		}
     	
    }


  for(pat::ElectronCollection::const_iterator iele = electron_d->begin() ; iele !=electron_d->end() ; ++iele)
    {
      
	      for(pat::TauCollection::const_iterator itau = tau_d->begin() ; itau !=tau_d->end() ; ++itau)
		{
		  //double deltaphi=deltaPhi(iele->phi(),itau->phi());
		  //double delta_R=deltaR(*iele,*itau);
		  //DeltaMatch->Fill(delta_R);                                                                                                                                                                                                                                     
		  //DeltaMatch->SetFillColor(kSpring);                                                                                                                                                                                                                                 
		  //DeltaMatch->GetXaxis()->SetTitle("#Delta_{R} High Mass Pair");                                                                                                                                                                                                               
		  //DeltaMatch->GetYaxis()->SetTitle("# of events"); 
		  //if (delta_R > mindR)
		  //{
		  dR_select=abs((iele->p4() + itau->p4()).mass());
		  //cout<<"Vmass: "<<Vmass_select<<endl;
		  //cout<<"Vmass Value: "<<Vmass_select<<endl;
		  dRSelect->Fill(dR_select);
		  dRSelect->SetFillColor(kRed);
		  dRSelect->GetYaxis()->SetTitle(" # of events");
		  dRSelect->GetXaxis()->SetTitle(" M_{v} of reconstructed Tau-ele Pair with lowest dR");
		  //  }
		}
     	
    }











  


   for(pat::ElectronCollection::const_iterator iele = electron_p->begin() ; iele !=electron_p->end() ; ++iele)
    {
      //++counter;

      //double delta_Rp=(*iele,*itau)


      for(pat::TauCollection::const_iterator itau = tau_p->begin() ; itau !=tau_p->end() ; ++itau)

	{

	  //double delta_Rp=deltaR(*iele,*itau);
	  //if (delta_Rp > mindR)
	  //{  
	  Pt_select=abs((iele->p4() + itau->p4()).mass());
	  PtSelect->Fill(Pt_select);
	  PtSelect->SetFillColor(kBlue);
	  PtSelect->GetYaxis()->SetTitle(" # of events");
	  PtSelect->GetXaxis()->SetTitle(" M_{v} of reconstructed High Di-Object Pt Tau-ele Pair");
  
	  //}
	}

    }






  for(pat::ElectronCollection::const_iterator iele = electrons->begin() ; iele !=electrons->end() ; ++iele,++ele_idx)
    {
      //size_ele=(electrons.product())->size();
      /*if( iele->genLepton() ){
	     float DeltaRel = ROOT::Math::VectorUtil::DeltaR(iele->genLepton()->p4(), iele->p4());
       
           DeltaRele->Fill(DeltaRel);
	   DeltaRele->SetFillColor(kBlue);
      }*/
      //countele++;
      //double dRMax=5000;
      dz_ele=iele->gsfTrack()->dz(p);
      dxy_ele=iele->gsfTrack()->dxy(p);
      
      //pat::Tau outputTau(*itau);
      pat::ElectronRef inputElectronRef(electrons, ele_idx);

      int tau_idx = 0;
      //vector <unsigned int> TausAlreadyInPlot;
      //int nDoublect= 0;  
      for(pat::TauCollection::const_iterator itau = Taus->begin() ; itau !=Taus->end() ; ++itau,++tau_idx)
	
	
	{ 

	  
	  pat::Tau outputTau(*itau);
	  pat::TauRef inputTauRef(Taus, tau_idx);
	  
	  

	  
	  pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(itau->leadChargedHadrCand().get());
	  dz_tau=packedLeadTauCand->dz(p);
	  dxy_tau=packedLeadTauCand->dxy(p);
	  //dz=itau->leadChargedHadrCand()->dz(p);
	  
	  dR = reco::deltaR(*iele, *itau);
	  Vmass=abs((iele->p4() + itau->p4()).mass());
	  Pt_di=((iele->p4() + itau->p4()).pt());
	  if( (itau->pt()>10) && abs((itau->eta())<2.3) && (itau->tauID("decayModeFinding")) && (itau->tauID("byIsolationMVArun2v1DBoldDMwLTraw") >-0.5) && (dxy_tau<0.2) && (dz_tau < 0.5) && (iele->isEB()) && (dxy_ele<0.05) &&(dxy_ele <0.10) )
	    {


	     
	      //dR = reco::deltaR(*iele, *itau);
	      DeltaR->Fill(dR);
	      //Vmass=abs((iele->p4() + itau->p4()).mass());
	      TausAlreadyInPlot.push_back(inputTauRef.key());
	      ElectronsAlreadyInPlot.push_back(inputElectronRef.key());
	    
	      Mvisible->Fill(Vmass);
	      Mvisible->GetXaxis()->SetTitle("M_{v} before cut ");
              Mvisible->GetYaxis()->SetTitle("# of events");

	      Pt->Fill(Pt_di);
	      
	      
	      Profile->Fill(dR,Vmass);
	      Profile2->Fill(dR,Pt_di);
	      Profile3->Fill(Vmass,Pt_di);
	      if((dR < 0.8) && (dR>0.05))
		{
		  Mvisiblecut->Fill(Vmass);
		  if (Vmass>20)

                    {
                      cout<< " High Mass Point" << endl;
                      cout<< "================"<<endl;
                    }

                  cout<<" Visible mass_EB : " << Vmass <<endl;

		  DeltaRcut->Fill(dR);
		  Profile4->Fill(Vmass,Pt_di);
		}
	      if (Vmax<Vmass)
		Vmax=Vmass;
	      
	    }
	  if((itau->pt()>10) && abs((itau->eta())<2.3) && (itau->tauID("decayModeFinding")) && (itau->tauID("byIsolationMVArun2v1DBoldDMwLTraw") >-0.5) && (dxy_tau<0.2) && (dz_tau < 0.5) && (iele->isEE()) && (dxy_ele<0.10) &&(dz_ele <0.20) )
            {
	      
              DeltaR->Fill(dR);
	      DeltaR->SetFillColor(kRed);
	      DeltaR->GetXaxis()->SetTitle("#DeltaR");
	      DeltaR->GetYaxis()->SetTitle("# of events");
	      TausAlreadyInPlot.push_back(inputTauRef.key());
	      ElectronsAlreadyInPlot.push_back(inputElectronRef.key());
             
              Mvisible->Fill(Vmass);
	      Mvisible->SetFillColor(kYellow);
	      
	      
	      Pt->Fill(Pt_di);
	      Pt->SetFillColor(kSpring);
	      Pt->GetXaxis()->SetTitle("Pt(Gev)");
	      Pt->GetYaxis()->SetTitle("# of events");
	      
	      Profile->Fill(dR,Vmass);
	      Profile->GetXaxis()->SetTitle("#DeltaR");
	      Profile->GetYaxis()->SetTitle("Mass_{visble}");
	      Profile->SetOption("COLZ");
	      

	      Profile2->GetYaxis()->SetTitle("Di-Object_{Pt}");
	      Profile2->GetXaxis()->SetTitle("#DeltaR");
	      Profile2->SetOption("COLZ");
	      
	      Profile3->GetXaxis()->SetTitle("M_{v}");
	      Profile3->GetYaxis()->SetTitle("Di-Object_{Pt}");
	      Profile3->SetOption("COLZ");

	      if((dR < 0.8) && (dR >0.05) )
		{
		  Mvisiblecut->Fill(Vmass);
		  if (Vmass>20)

		    {		    
		      cout<< " High Mass Point" << endl;
		      cout<< "================"<<endl;
		    } 
  
		  cout<<" Visible mass_EE : " << Vmass <<endl;
	   
		  Mvisiblecut->SetFillColor(kYellow);
		  Mvisiblecut->GetXaxis()->SetTitle("M_{v} (Gev)");
		  Mvisiblecut->GetYaxis()->SetTitle("# of events");
                  
		  
		  DeltaRcut->Fill(dR);
		  DeltaRcut->SetFillColor(kRed);
		  // TausAlreadyinPlot.Push_back(inputTauRef->key());
		  
		  
		  Profile4->Fill(Vmass,Pt_di);
		  Profile4->GetXaxis()->SetTitle("Mass_{visible} (Gev)");
		  Profile4->GetYaxis()->SetTitle("Di-Object_{Pt} (Gev)");
		  Profile4->SetOption("COLZ");
		}

             
            }
	  
	  
	   
	 
	  if((std::find(TausAlreadyInPlot.begin(), TausAlreadyInPlot.end(),(inputTauRef.key()))!= TausAlreadyInPlot.end()) && (std::find(ElectronsAlreadyInPlot.begin(), ElectronsAlreadyInPlot.end(),(inputElectronRef.key()))!= ElectronsAlreadyInPlot.end()) )
	{nDoublect++;
	  	}
	  
	}
    }
  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
PrimaryAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PrimaryAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PrimaryAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrimaryAnalyzer);
