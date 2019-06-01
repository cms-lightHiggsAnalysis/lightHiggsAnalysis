// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer
// Class:      GenMassPlotter
// 
/**\class PrimaryAnalyzer GenMassPlotter.cc GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/plugins/GenMassPlotter.cc

 Description: [one line class summary]
Takes the gen particle collection and plots the invariant mass,visible mass of all the products coming from a Pseudoscalar and the  also does the visible mass of products in Tau Pathway
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

class GenMassPlotter : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GenMassPlotter(const edm::ParameterSet&);
      ~GenMassPlotter();
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
  //edm::EDGetTokenT<pat::ElectronCollection> electronSrc_;
  //edm::EDGetTokenT<pat::TauCollection> TauSrc_;
  //edm::EDGetTokenT<reco::VertexCollection> vtx_;
  //edm::EDGetTokenT<reco::BeamSpot> thebs_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> >packedGenToken_;
  // internal variables for gen matching
  double maxDeltaR_;
  int matchcount_ele=0;
  //Initialize Histograms
  //TH1D *DeltaR;
  //TH1D *Mvisible;
  //TH1F *DeltaRele;
  //TH2F *Profile;
  //TH1D *DeltaRcut;
  //TH1D *Mvisiblecut;
  TH1D *MvisibleGen;
  TH1D *MinvariantGen;
  TH1D *MvisibleGenTau;
  TH1D *MvisibleGenTau2;
  TH1D *MvisibleGenTau3;
  int Tau_e_loop=0;
  int Tau_had_loop=0;
  int tetefill=0;
  int tethfill=0;
  int ththfill=0;

  
  
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
GenMassPlotter::GenMassPlotter(const edm::ParameterSet& iConfig):
 
  prunedGenToken_ (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_ (consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
  
{
  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fl;
  
  MvisibleGen =fl->make<TH1D>("MvisibleGen" , "Mass_{visible} of Gen Particles",100,-0.5,49.5);
  MinvariantGen=fl->make<TH1D>("MinvariantGen","Mass_{invariant} of Gen Particles",100,-0.5,49.5);
  MvisibleGenTau=fl->make<TH1D>("MvisibleGenTau","Mass_{visible} of Gen Particles for Tau_{e},Tau_{had}",500,-0.5,49.5);
  MvisibleGenTau2=fl->make<TH1D>("MvisibleGenTau","Mass_{visible} of Gen Particles involving Tau_{e},Tau_{e}",500,-0.5,49.5);
  MvisibleGenTau3=fl->make<TH1D>("MvisibleGenTau","Mass_{visible} of Gen Particles involving Tau_{had},Tau_{had}",500,-0.5,49.5);



}


GenMassPlotter::~GenMassPlotter()
{
  
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}
bool GenMassPlotter::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
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


bool GenMassPlotter::isDaughter(const reco::Candidate* Daughter, const reco::Candidate * particle)
{
  //particle is already the ancestor                                                                                                                                                                                                                                           
  if(Daughter== particle ) return true;
  
  //otherwise loop on Daughters, if any and return true if the ancestor is found                                                                                                                                                                                                 
  for(size_t i=0;i< particle->numberOfDaughters();i++)
    {
      if(isDaughter(Daughter,particle->daughter(i))) return true;
    }
  //if we did not return yet, then particle and ancestor are not relatives                                                                                                                                                                                                    
  return false;
}


bool GenMassPlotter::isPdgMatch(const reco::Candidate * ancestor, const reco::Candidate * particle, int Id)
{
  
  //Check the mother in the packed collection and a candidate in the packed collection against a common pdgId->The pdgId can be changed to check for various ancestors
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
GenMassPlotter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  
  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_, pruned);
  //Handle<pat::PackedGenParticle> packed;
  Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_, packed);


  //vector <unsigned int> CandInPlot;
  vector<pat::PackedGenParticle*>Packed;
    

//Invariant mass plot
  for(size_t i=0; i<(pruned.product())->size() ;i++)
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
		  //cout << "     PdgID: " << (*packed)[j].pdgId() <<endl;
		  //cout<< " key: " << motherInPrunedCollectionRef.key()<<endl;
		  //matchcount_ele++
		  
		      //cout<<"Photon event"<<endl;
		      //cout<<"*******************"<<endl;
		  //cout<<"mother pdgID: " <<(*packed)[j].mother(0)->pdgId()<<endl;
		  //cout<<"mother key: " <<(*packed)[j].motherRef().key()<<endl;
		  //cout<<"Grandmother pdgID: "<<(*packed)[j].mother(0)->mother(0)->pdgId()<<endl;
		  //cout<<"Grandmother key: "<<(*packed)[j].motherRef()->motherRef().key()<<endl;
		  
		  inv_v4 +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
		}
	      
	      
	    }  
	      MinvariantGen->Fill(abs(inv_v4.M()));
	      MinvariantGen->SetFillColor(kMagenta+2);
	      MinvariantGen->GetXaxis()->SetTitle("Mass_{invariant} of Generator Particles (Gev)");
	      MinvariantGen->GetYaxis()->SetTitle("# of events");

	    
	  
	  
	}
      
    }
  //VIsible mass for Tau Pathway
  for(size_t i=0; i<(pruned.product())->size();i++)
    {
      
      
      //cout<<"PdgID: "<<(*pruned)[i].pdgId()<<endl;
      //cout<<" index: "<<i<<endl;
      
      
      if(fabs((*pruned)[i].pdgId())==36)
        {
	  TLorentzVector invariant_v4_Tau_e;
	  TLorentzVector invariant_v4_Tau_had;
	  TLorentzVector invariant_v4_Tau_mu;
          const Candidate * PsuedoScalar = &(*pruned)[i];
          //cout << "PdgID: " << PsuedoScalar->pdgId() <<endl;
	  //cout << "found daugthers: "<< endl;
	  //cout<<"PdgID: "<<(*pruned)[i].pdgId()<<endl;
	  //int jid=0;
          for(size_t j=0; j<(packed.product())->size() ;j++)
            {
	      
	      
              const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
              const GenParticleRef motherInPrunedCollectionRef=(*packed)[j].motherRef();
	      //pat::PackedGenParticleRef PackedRef(packed,jid);

	      //Adding for Tau_e
	      
		  
		

              if((motherInPrunedCollection != nullptr) && (isAncestor( PsuedoScalar , motherInPrunedCollection)) &&  /*(fabs((*packed)[j].pdgId())==11) && */ (fabs((motherInPrunedCollection->pdgId()))==15) && (motherInPrunedCollection->status()!=1) )
                {
                  //cout << "Tau_e_PdgID: " << (*packed)[j].pdgId() <<endl;
		  //cout<<"mother pdgID: " <<(*packed)[j].mother(0)->pdgId()<<endl;
		  //cout<<"mother key: " <<(*packed)[j].motherRef().key()<<endl;
		  //cout<<"Grandmother pdgID: "<<(*packed)[j].mother(0)->mother(0)->pdgId()<<endl;
		  //cout<<"Grandmother key: "<<(*packed)[j].motherRef()->motherRef().key()<<endl;
		  //cout<<"Great-Grandmother pdgID: "<<(*packed)[j].mother(0)->mother(0)->mother(0)->pdgId()<<endl;
		  //cout<<"Great-Grandmother key: "<<(*packed)[j].motherRef()->motherRef()->motherRef().key()<<endl;
		  if( (fabs((*packed)[j].pdgId())==13) )
		
		  
		    {
		      
                      invariant_v4_Tau_mu +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());


		    
		    }
		 
		  if((fabs((*packed)[j].pdgId())==11) && ((*packed)[j].statusFlags().isFirstCopy()) )
                 
		    {   
		      // cout << "Tau_e_PdgID: " << (*packed)[j].pdgId() <<endl;
		      
		      
		      invariant_v4_Tau_e +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
		      ++Tau_e_loop;
		      //cout << "Tau_e_PdgID: " << (*packed)[j].pdgId() <<endl;
		      //cout<<"X: "<<(*packed)[j].p4().X()<<endl;
		      //cout<<"Y: "<<(*packed)[j].p4().Y()<<endl;
		      //cout<<"Z: "<<(*packed)[j].p4().Z()<<endl;
		      //cout<<"T: "<<(*packed)[j].p4().T()<<endl;

		    }
		   
		    
		}
	      
	      //Adding for Tau_had
	      /*Veto against ancestry of Muons and electrons
		Veto against Muon and electron final products
		Check for first copy of product
		Veto against any Neutrino in final product
		Check for Tau somewhere in ancestry
		Check for Pseudoscalars in start of ancestry
	      */
	      for( size_t k=0; k<(pruned.product())->size();k++)
		{
		  
		  if((motherInPrunedCollection != nullptr) && !(isPdgMatch(&(*pruned)[k],motherInPrunedCollection,11)) && !(isPdgMatch(&(*pruned)[k],motherInPrunedCollection,13)) && (isAncestor( PsuedoScalar , motherInPrunedCollection)) && (fabs((*packed)[j].pdgId())!=11) && (fabs((*packed)[j].pdgId())!=13) && (motherInPrunedCollection->status()!=1) && ((*packed)[j].statusFlags().isFirstCopy()) && (isPdgMatch(&(*pruned)[k],motherInPrunedCollection,15)) && (fabs((*packed)[j].pdgId())!=12) && (fabs((*packed)[j].pdgId())!=16)  && (fabs((*packed)[j].pdgId())!=14)  )
	        
		    
		    
		    
		    
		    
		    {
		      
		      cout<<"Tau_had_pdgID: "<<(*packed)[j].pdgId() <<endl;
		      cout<<"mother pdgID: " <<(*packed)[j].mother(0)->pdgId()<<endl;
		      cout<<"mother key: " <<(*packed)[j].motherRef().key()<<endl;
		      //cout<<"Grandmother pdgID: "<<(*packed)[j].mother(0)->mother(0)->pdgId()<<endl;
		      //cout<<"Grandmother key: "<<(*packed)[j].motherRef()->motherRef().key()<<endl;
		      //cout<<"Great-Grandmother pdgID: "<<(*packed)[j].mother(0)->mother(0)->mother(0)->pdgId()<<endl;
		      //cout<<"Great-Grandmother key: "<<(*packed)[j].motherRef()->motherRef()->motherRef().key()<<endl;
		      //CandInPlot.push_back(PackedRef.key());
		      
		      // if(std::find(Packed.begin(),Packed.end(),((*packed)[j]))!=Packed.end())
		      std::vector<pat::PackedGenParticle*>::iterator check;
		      check=std::find(Packed.begin(),Packed.end(),const_cast<pat::PackedGenParticle*>((*packed).ptrAt(j).get()));
		      if(check==Packed.end())
		      {
			Packed.push_back(const_cast<pat::PackedGenParticle*>((*packed).ptrAt(j).get()));
			invariant_v4_Tau_had +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
		      
			++Tau_had_loop;
		      //invariant_v4_Tau_e +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
		      cout<<"X: "<<(*packed)[j].p4().X()<<endl;
		      cout<<"Y: "<<(*packed)[j].p4().Y()<<endl;
		      cout<<"Z: "<<(*packed)[j].p4().Z()<<endl;
		      cout<<"T: "<<(*packed)[j].p4().T()<<endl;
		      }
 
		      
		      		      
		      
		    }
		  
		}
	      
	    
	    }
	  if ( invariant_v4_Tau_e.X()!=0 &&  invariant_v4_Tau_e.Y()!=0 &&  (invariant_v4_Tau_e.Z()!=0) &&  (invariant_v4_Tau_e.T()!=0) && (invariant_v4_Tau_e.Mag()>0) && (invariant_v4_Tau_had.X()!=0) && (invariant_v4_Tau_had.Y()!=0) && (invariant_v4_Tau_had.Z()!=0) && (invariant_v4_Tau_had.T()!=0) && (invariant_v4_Tau_had.Mag() >0) && (invariant_v4_Tau_mu.X()==0) && (invariant_v4_Tau_mu.Y()==0) && (invariant_v4_Tau_mu.Z()==0) && (invariant_v4_Tau_mu.T()==0) && (invariant_v4_Tau_mu.Mag()==0) )   	    
	    {MvisibleGenTau->Fill((invariant_v4_Tau_e+invariant_v4_Tau_had).Mag());
	    
	  //MinvariantGenTau->Fill((invariant_v4_Tau_had).M());
	  //MvisibleGenTau->Fill(abs((invariant_v4_Tau_e).M()));
      
	  MvisibleGenTau->SetFillColor(kPink+2);
	  MvisibleGenTau->GetXaxis()->SetTitle("Mass_{visible} of Generator Particles in Tau_{e},Tau_{had} final state (Gev)");
	  MvisibleGenTau->GetYaxis()->SetTitle("# of events");
	  ++tethfill;
  
	    }
	  
	  if ( invariant_v4_Tau_e.X()!=0 &&  invariant_v4_Tau_e.Y()!=0 &&  (invariant_v4_Tau_e.Z()!=0) &&  (invariant_v4_Tau_e.T()!=0) && (invariant_v4_Tau_e.Mag()>0) && (invariant_v4_Tau_had.X()==0) && (invariant_v4_Tau_had.Y()==0) && (invariant_v4_Tau_had.Z()==0) && (invariant_v4_Tau_had.T()==0) && (invariant_v4_Tau_had.Mag()==0) && (invariant_v4_Tau_mu.X()==0) && (invariant_v4_Tau_mu.Y()==0) && (invariant_v4_Tau_mu.Z()==0) && (invariant_v4_Tau_mu.T()==0) && (invariant_v4_Tau_mu.Mag()==0) )

	    
	    {
	      MvisibleGenTau2->Fill((invariant_v4_Tau_e+invariant_v4_Tau_had).Mag());
	      MvisibleGenTau2->SetFillColor(kAzure+6);
	      MvisibleGenTau2->GetXaxis()->SetTitle("Mass_{visible} of Generator Particles in Tau_{e},Tau_{e} final state (Gev)");
	      MvisibleGenTau2->GetYaxis()->SetTitle("# of events");
	      ++tetefill;

	    }
	  

	  if ( invariant_v4_Tau_e.X()==0 &&  invariant_v4_Tau_e.Y()==0 &&  (invariant_v4_Tau_e.Z()==0) &&  (invariant_v4_Tau_e.T()==0) && (invariant_v4_Tau_e.Mag()==0) && (invariant_v4_Tau_had.X()!=0) && (invariant_v4_Tau_had.Y()!=0) && (invariant_v4_Tau_had.Z()!=0) && (invariant_v4_Tau_had.T()!=0) && (invariant_v4_Tau_had.Mag()>0) && (invariant_v4_Tau_mu.X()==0) && (invariant_v4_Tau_mu.Y()==0) && (invariant_v4_Tau_mu.Z()==0) && (invariant_v4_Tau_mu.T()==0) && (invariant_v4_Tau_mu.Mag()==0) )

	    {
	      MvisibleGenTau3->Fill((/*invariant_v4_Tau_e+*/invariant_v4_Tau_had).Mag());
              MvisibleGenTau3->SetFillColor(kBlue+3);
              MvisibleGenTau3->GetXaxis()->SetTitle("Mass_{visible} of Generator Particles in Tau_{had},Tau_{had} final state (Gev)");
              MvisibleGenTau3->GetYaxis()->SetTitle("# of events");
	      ++ththfill;
	      





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
                  //cout << "     PdgID: " << (*packed)[j].pdgId() <<endl;
                  //cout<< " key: " << motherInPrunedCollectionRef.key()<<endl;
                  //matchcount_ele++;                                                                                                                                                                                                                                          
		  
		  
		  visible_v4 +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
                } 
              
	    }
	      MvisibleGen->Fill(abs(visible_v4.M()));
	      MvisibleGen->SetFillColor(kGreen+2);
	      MvisibleGen->GetXaxis()->SetTitle("Mass_{visible} of Generator Particles (Gev)");
	      MvisibleGen->GetYaxis()->SetTitle("# of events");

	    
	}
    }
    
   cout<<"te_l: "<<Tau_e_loop<<" th_l: "<<Tau_had_loop<<" Total: "<<Tau_e_loop+Tau_had_loop <<endl;
   cout<<" tete: "<<tetefill<<" teth: "<<tethfill<<" thth :"<<ththfill<<" Tot: "<<tetefill+ththfill+tethfill<<endl;

}


// ------------ method called once each job just before starting event loop  ------------
void 
GenMassPlotter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenMassPlotter::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenMassPlotter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenMassPlotter);

