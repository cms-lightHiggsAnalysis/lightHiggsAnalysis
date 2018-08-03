
// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer
// Class:      GenMatcher
// 
/**\class PrimaryAnalyzer GenMatcher.cc GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/plugins/GenMatcher.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Redwan Habibullah
//         Created:  Wed, 19 June 2018 20:01:49 GMT
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

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"

#include"DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "TAttFill.h"
#include "Math/VectorUtil.h"

#include "Math/LorentzVector.h"

typedef math::XYZPoint Point;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class GenMatcher : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GenMatcher(const edm::ParameterSet&);
      ~GenMatcher();
  bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);
  bool isDaughter(const reco::Candidate * Daughter, const reco::Candidate * particle);
  bool isPdgMatch(const reco::Candidate * Ancestor, const reco::Candidate * particle, int Id);
  
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
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> >packedGenToken_;
  // internal variables for gen matching
  double maxDeltaR_;
  int matchcount_ele=0;
  //Initialize Histograms
  TH1D *DeltaR;
  TH1D *Mvisible;
  TH1F *DeltaRele;
  TH2F *Profile;
  TH1D *DeltaRcut;
  TH1D *Mvisiblecut;
  TH1D *MvisibleGen;
  TH1D *MinvariantGen;
  TH1D *MvisibleGenTau;

  
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
GenMatcher::GenMatcher(const edm::ParameterSet& iConfig):
  //declare a tag
  electronSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  TauSrc_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("Taus"))),
  vtx_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
  thebs_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BM"))),
  //minPt_ (iConfig.getParameter<double>("minPt")),
  //maxEta_ (iConfig.getParameter<double>("maxEta")),
  //maxDeltaR_ (iConfig.getParameter<double>("maxDeltaR")),
  prunedGenToken_ (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_ (consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
  
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fl;
  DeltaR = fl->make<TH1D>("DeltaR" , "#Delta_{R}" , 1000 , +0.005 , 5.005 );
  DeltaRcut = fl->make<TH1D>("DeltaRcut" , "#Delta_{R} after cut" , 1000 , +0.005 , 5.005 );
  Mvisible = fl->make<TH1D>("Mvisible" , "Mass_{visible}" , 100 ,0.5 , 100.5 );
  Mvisiblecut =fl->make<TH1D>("Mvisiblecut" , "Mass_{visible} after cut" , 100 ,+0.5 , 50.5 );
  DeltaRele = fl->make<TH1F>("DeltaRele","DeltaR_{gen/ele}",100,+.005,0.055);
  Profile =fl->make<TH2F>("Comparison","Profile plot of Mass_{visible} vs #DeltaR",1000,+.005,5.005,100,+0.5,500.5);
  //pdgId=fl->make<TH1F>("PdgId","PdgIds for matched electrons",100,+0.5,50.5);
  MvisibleGen =fl->make<TH1D>("MvisibleGen" , "Mass_{visible} of Gen Particles",100,+0.5,50.5);
  MinvariantGen=fl->make<TH1D>("MinvariantGen","Mass_{invariant} of Gen Particles",100,+0.5,50.5);
  MvisibleGenTau=fl->make<TH1D>("MvisibleGen","Mass_{visible} of Gen Particles involving Tau pathway",100,+0.5,50.5);

}


GenMatcher::~GenMatcher()
{
  
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}
bool GenMatcher::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
{
  //particle is already the ancestor
  if(ancestor == particle ) return true;
  
  //otherwise loop on mothers, if any and return true if the ancestor is found
  for(size_t i=0;i< particle->numberOfMothers();i++)
    {
      if(isAncestor(ancestor,particle->mother(i))) return true;
    }
  //if we did not return yet, then particle and ancestor are not relatives
  return false;
}


bool GenMatcher::isDaughter(const reco::Candidate* Daughter, const reco::Candidate * particle)
{
  //particle is already the ancestor                                                                                                                                                                                                                                           
  if(Daughter== particle ) return true;
  
  //otherwise loop on Daughters, if any and return true if the ancestor is found                                                                                                                                                                                                 
  for(size_t i=0;i< particle->numberOfDaughters();i++)
    {
      if(isDaughter(Daughter,particle->daughter(i))) return true;
    }
  //if we did not return yet, then particle and ancestor are not relatives                                                                                                                                                                                                     t
  return false;
}


bool GenMatcher::isPdgMatch(const reco::Candidate * ancestor, const reco::Candidate * particle, int Id)
{
  
  
  if((fabs(ancestor->pdgId())==Id) && (fabs(particle->pdgId())==Id) && (ancestor==particle)) return true;
  for(size_t i=0;i<particle->numberOfMothers();i++)
    {
      if(isPdgMatch(ancestor,particle->mother(i),Id)) return true;
    }
  return false;
}




//
// member functions
//

// ------------ method called for each event  ------------
void
GenMatcher::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  Handle<reco::BeamSpot> thebs;
  iEvent.getByToken(thebs_,thebs);
  //Handle<reco::GenParticleCollection> pruned;
  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_, pruned);
  //Handle<pat::PackedGenParticle> packed;
  Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_, packed);


  
  const reco::Vertex& pv=*Vertex->begin();
  Point p;
  p=pv.position();
  double dR=999999;
  double Vmass=999999;
  double Vmax=0;
  // int countele=0;
  //int counttau=0;
  double dz_tau;
  double dxy_tau;
  double dz_ele;
  double dxy_ele;
  // int size_ele= 0;
  //int size_tau= 0;
  

//Invariant mass plot
  for(size_t i=0; i<(pruned.product())->size();i++)
    {
      if(fabs((*pruned)[i].pdgId())==36)
	{
	  TLorentzVector inv_v4;
	  const Candidate * PsuedoScalar = &(*pruned)[i];
	  //cout << "PdgID: " << PsuedoScalar->pdgId() <<endl;
	  //std::cout << "  found daugthers: " << std::endl;
	  for(size_t j=0; j<(packed.product())->size() ;j++)
	    {

      
	      const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
	      const GenParticleRef motherInPrunedCollectionRef=(*packed)[j].motherRef();
	      if((motherInPrunedCollection != nullptr) && (isAncestor( PsuedoScalar , motherInPrunedCollection)))
		{
		  cout << "     PdgID: " << (*packed)[j].pdgId() <<endl;
		  cout<< " key: " << motherInPrunedCollectionRef.key()<<endl;
		  //matchcount_ele++;
	
		  
		      //cout<<"Photon event"<<endl;
		      //cout<<"*******************"<<endl;
		  cout<<"mother pdgID: " <<(*packed)[j].mother(0)->pdgId()<<endl;
		  cout<<"mother key: " <<(*packed)[j].motherRef().key()<<endl;
		  cout<<"Grandmother pdgID: "<<(*packed)[j].mother(0)->mother(0)->pdgId()<<endl;
		  cout<<"Grandmother key: "<<(*packed)[j].motherRef()->motherRef().key()<<endl;
		  
		  inv_v4 +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
		}
	      //else if((motherInPrunedCollection != nullptr)){
	      //cout << "Not a pseudoscalar descendant" << endl;
	      //cout<<"mother pdgID: " <<(*packed)[j].mother(0)->pdgId()<<endl;
	      //cout<<"mother key: " <<(*packed)[j].motherRef().key()<<endl;
		//cout<<"Grandmother pdgID: "<<(*packed)[j].mother(0)->mother(0)->pdgId()<<endl;
		//cout<<"Grandmother key: "<<(*packed)[j].motherRef()->motherRef().key()<<endl;
	    }

	  //Checking the low mass tail
	  /* if(abs(visible_v4.M())<19)
	    {
	      for(size_t j=0; j<(packed.product())->size() ;j++)
		{


		  const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
		  const GenParticleRef motherInPrunedCollectionRef=(*packed)[j].motherRef();
		  if((motherInPrunedCollection != nullptr) && (isAncestor( PsuedoScalar , motherInPrunedCollection)))
		    {
		      cout << "     PdgID: " << (*packed)[j].pdgId() <<endl;
		      //cout<< " key: " << motherInPrunedCollectionRef.key()<<endl;                                                                                                                                                                                                
		      //matchcount_ele++;                                                                                                                                                                                                                                          

		      //if ((*packed)[j].pdgId()==22)                                                                                                                                                                                                                              
		      //{                                                                                                                                                                                                                                                          
                      //cout<<"Photon event"<<endl;                                                                                                                                                                                                                            
                      //cout<<"*******************"<<endl;                                                                                                                                                                                                                     
                      cout<<"mother pdgID: " <<(*packed)[j].mother(0)->pdgId()<<endl;
                      cout<<"mother key: " <<(*packed)[j].motherRef().key()<<endl;
                      cout<<"Grandmother pdgID: "<<(*packed)[j].mother(0)->mother(0)->pdgId()<<endl;
                      cout<<"Grandmother key: "<<(*packed)[j].motherRef()->motherRef().key()<<endl;
		    





		    }

		  else if((motherInPrunedCollection != nullptr)){                                                                                                                                                                                                                
		    cout << "Not a pseudoscalar descendant" << endl;
		    cout << "     PdgID: " << (*packed)[j].pdgId() <<endl;
		    
		    //cout << "Not a pseudoscalar descendant" << endl;                                                                                                                                                                                                               
		    cout<<"mother pdgID: " <<(*packed)[j].mother(0)->pdgId()<<endl;                                                                                                                                                                                                
		  //cout<<"mother key: " <<(*packed)[j].motherRef().key()<<endl;                                                                                                                                                                                                   
		  //cout<<"Grandmother pdgID: "<<(*packed)[j].mother(0)->mother(0)->pdgId()<<endl;                                                                                                                                                                               
		  //cout<<"Grandmother key: "<<(*packed)[j].motherRef()->motherRef().key()<<endl;  
		  }
		
		}
	    } */
	  
	  MinvariantGen->Fill(abs(inv_v4.M()));
	  MinvariantGen->SetFillColor(kMagenta+2);
	}
	

    }

  

   for(size_t i=0, k=0; i<(pruned.product())->size() && k<(pruned.product())->size();i++,k++)
    {


      //cout<<"PdgID: "<<(*pruned)[i].pdgId()<<endl;
      //cout<<" index: "<<i<<endl;


      if(fabs((*pruned)[i].pdgId())==36)
        {
          TLorentzVector invariant_v4_Tau_e;
	  TLorentzVector invariant_v4_Tau_had;
          const Candidate * PsuedoScalar = &(*pruned)[i];
          //cout << "PdgID: " << PsuedoScalar->pdgId() <<endl;
	  //cout << "found daugthers: "<< endl;
	  //cout<<"PdgID: "<<(*pruned)[i].pdgId()<<endl;
 
          for(size_t j=0; j<(packed.product())->size() ;j++)
            {


              const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
              const GenParticleRef motherInPrunedCollectionRef=(*packed)[j].motherRef();
	      //Adding for Tau_e
              if((motherInPrunedCollection != nullptr) && (isAncestor( PsuedoScalar , motherInPrunedCollection)) &&  (fabs((*packed)[j].pdgId())==11) && (fabs((motherInPrunedCollection->pdgId()))==15) && (motherInPrunedCollection->status()!=1) )
                {
                  cout << "Tau_e_PdgID: " << (*packed)[j].pdgId() <<endl;
                  //cout<< " key: " << motherInPrunedCollectionRef.key()<<endl;
                  //matchcount_ele++;

                
                  invariant_v4_Tau_e +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
		  //cout<<"value_e: "<< invariant_v4_Tau_e 
		}

	      //Adding for Tau_had
	      
	      
	      if((motherInPrunedCollection != nullptr) && !(isPdgMatch(&(*pruned)[k],motherInPrunedCollection,11)) && !(isPdgMatch(&(*pruned)[k],motherInPrunedCollection,13)) && (isAncestor( PsuedoScalar , motherInPrunedCollection)) && (fabs((*packed)[j].pdgId())!=11) && (fabs((*packed)[j].pdgId())!=13) && (motherInPrunedCollection->status()!=1) && ((*packed)[j].statusFlags().isFirstCopy()) && (isPdgMatch(&(*pruned)[k],motherInPrunedCollection,15)) && (fabs((*packed)[j].pdgId())!=12) && (fabs((*packed)[j].pdgId())!=16)  && (fabs((*packed)[j].pdgId())!=14)  )
	         
                                                                                                                                                       




		{

		  cout<<"Tau_had_pdgID: "<<(*packed)[j].pdgId() <<endl;
		  invariant_v4_Tau_had +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
		  //cout<<"value_had: "<< invariant_v4_Tau_had<<endl;





		}

	      MvisibleGenTau->Fill(fabs((invariant_v4_Tau_e+invariant_v4_Tau_had).M()));
	      //MinvariantGen->Fill((invariant_v4_Tau_had).M());
	      MvisibleGenTau->SetFillColor(kPink+2);

		

	    }

	}
    }

 	      
 //visiblemass

  for(size_t i=0; i<(pruned.product())->size();i++)
    {
      if(fabs((*pruned)[i].pdgId())==36)
        {
          TLorentzVector visible_v4;
          const Candidate * PsuedoScalar = &(*pruned)[i];
          //cout << "PdgID: " << PsuedoScalar->pdgId() <<endl;                                                                                                                                                                                                                 
          //std::cout << "  found daugthers: " << std::endl;                                                                                                                                                                                                                   
          for(size_t j=0; j<(packed.product())->size() ;j++)
            {


              const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
              const GenParticleRef motherInPrunedCollectionRef=(*packed)[j].motherRef();
              if((motherInPrunedCollection != nullptr) && (isAncestor( PsuedoScalar , motherInPrunedCollection)) && (fabs((*packed)[j].pdgId())!=12) && (fabs((*packed)[j].pdgId())!=16)  && (fabs((*packed)[j].pdgId())!=14) )
                {
                  cout << "     PdgID: " << (*packed)[j].pdgId() <<endl;
                  cout<< " key: " << motherInPrunedCollectionRef.key()<<endl;
                  //matchcount_ele++;                                                                                                                                                                                                                                          


		  //cout<<"Photon event"<<endl;                                                                                                                                                                                                                            
		  //cout<<"*******************"<<endl;                                                                                                                                                                                                                     
                  //cout<<"mother pdgID: " <<(*packed)[j].mother(0)->pdgId()<<endl;
                  //cout<<"mother key: " <<(*packed)[j].motherRef().key()<<endl;
                  //cout<<"Grandmother pdgID: "<<(*packed)[j].mother(0)->mother(0)->pdgId()<<endl;
                  //cout<<"Grandmother key: "<<(*packed)[j].motherRef()->motherRef().key()<<endl;

                  visible_v4 +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
                } 
                                                                                                                                                                                                                                        

	      MvisibleGen->Fill(abs(visible_v4.M()));
	      MvisibleGen->SetFillColor(kGreen+2);
	    }
	}
    }


    
		      
 	   
   
                                                                                             

	       
  
 








for(pat::ElectronCollection::const_iterator iele = electrons->begin() ; iele !=electrons->end() ; ++iele)
    {
      //size_ele=(electrons.product())->size();
      /*  if(iele->genLepton()->motherRef().isNonnull())
      {
	cout<<"1"<<endl;
x	if(fabs(iele->genLepton()->motherRef()->pdgId())==15)
   {
     cout<<"2"<<endl;
	
	  matchcount_ele++;
	  pdgId->Fill((float)(iele->genLepton()->motherRef()->pdgId()));

   }
	    
   }*/

      


      //The electron matching loop
      if( iele->genLepton() ){
	     float DeltaRel = ROOT::Math::VectorUtil::DeltaR(iele->genLepton()->p4(), iele->p4());
       
           DeltaRele->Fill(DeltaRel);
	   DeltaRele->SetFillColor(kBlue);
      }
      //countele++;
      //double dRMax=5000;
      dz_ele=iele->gsfTrack()->dz(p);
      dxy_ele=iele->gsfTrack()->dxy(p);

      for(pat::TauCollection::const_iterator itau = Taus->begin() ; itau !=Taus->end() ; ++itau)
	
	
	{  // counttau++;	  
	  // p=pv.position();
	  // size_tau=(Taus.product())->size();
	  
	  pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(itau->leadChargedHadrCand().get());
	  dz_tau=packedLeadTauCand->dz(p);
	  dxy_tau=packedLeadTauCand->dxy(p);
	  //dz=itau->leadChargedHadrCand()->dz(p);
	  
	  dR = reco::deltaR(*iele, *itau);
	  Vmass=abs((iele->p4() + itau->p4()).mass());
	  if( (itau->pt()>10) && abs((itau->eta())<2.3) && (itau->tauID("decayModeFinding")) && (itau->tauID("byIsolationMVArun2v1DBoldDMwLTraw") >-0.5) && (dxy_tau<0.2) && (dz_tau < 0.5) && (iele->isEB()) && (dxy_ele<0.05) &&(dxy_ele <0.10) )
	    {
	      //dR = reco::deltaR(*iele, *itau);
	      DeltaR->Fill(dR);
	      //Vmass=abs((iele->p4() + itau->p4()).mass());
	      Mvisible->Fill(Vmass);
	      Profile->Fill(dR,Vmass);
	      if(dR < 0.8)
		{
		  Mvisiblecut->Fill(Vmass);
		  DeltaRcut->Fill(dR);
		}
	      if (Vmax<Vmass)
		Vmax=Vmass;
	      //if (size_ele==1 && (size_ele==size_tau))
	      //continue;
	    }
	  if((itau->pt()>10) && abs((itau->eta())<2.3) && (itau->tauID("decayModeFinding")) && (itau->tauID("byIsolationMVArun2v1DBoldDMwLTraw") >-0.5) && (dxy_tau<0.2) && (dz_tau < 0.5) && (iele->isEE()) && (dxy_ele<0.10) &&(dz_ele <0.20) )
            {
              //dR = reco::deltaR(*iele, *itau);
              DeltaR->Fill(dR);
	      DeltaR->SetFillColor(kRed);
              //Vmass=abs((iele->p4() + itau->p4()).mass());
              Mvisible->Fill(Vmass);
	      Mvisible->SetFillColor(kYellow);
	      Profile->Fill(dR,Vmass);
	      Profile->GetXaxis()->SetTitle("#DeltaR");
	      Profile->GetYaxis()->SetTitle("Mass_{visble}");
	      Profile->SetOption("COLZ");
	      if(dR < 0.8)
		{
		  Mvisiblecut->Fill(Vmass);
		  Mvisiblecut->SetFillColor(kYellow);
                  DeltaRcut->Fill(dR);
		  DeltaRcut->SetFillColor(kRed);
		}
              if (Vmax<Vmass)
                Vmax=Vmass;
	      //if (size_ele==1 && (size_ele==size_tau))
	      //continue;
            }
	  
	  
	  
	}
    }
  //cout<<" "<<Vmax<<endl;
  //cout<<"Size_ele "<<size_ele<<endl;
  
  //cout<<"Size_tau "<<size_tau<<endl;

  //cout<<"ele:"<<countele<<"tau:"<<counttau<<endl;
  


}
// ------------ method called once each job just before starting event loop  ------------
void 
GenMatcher::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenMatcher::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenMatcher::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenMatcher);
