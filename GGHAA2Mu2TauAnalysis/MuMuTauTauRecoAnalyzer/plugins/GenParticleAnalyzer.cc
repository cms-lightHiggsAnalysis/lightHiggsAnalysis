// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer
// Class:      GenParticleAnalyzer
// 
/**\class GenParticleAnalyzer GenParticleAnalyzer.cc GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/plugins/GenParticleAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Redwan Habibullah
//         Created:  Mon, 25 Feb 2019 21:02:55 GMT
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
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include"DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "TLorentzVector.h"
#include "TMath.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class GenParticleAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GenParticleAnalyzer(const edm::ParameterSet&);
      ~GenParticleAnalyzer();
  bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);
  bool isDaughter(const reco::Candidate * Daughter, const reco::Candidate * particle);
  bool isPdgMatch(const reco::Candidate * Ancestor, const reco::Candidate * particle, int Id);


      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> >packedGenToken_;
  
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
GenParticleAnalyzer::GenParticleAnalyzer(const edm::ParameterSet& iConfig):
  
  prunedGenToken_ (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_ (consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
{
   //now do what ever initialization is needed
  //  usesResource("TFileService");

}


GenParticleAnalyzer::~GenParticleAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}
bool GenParticleAnalyzer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
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


bool GenParticleAnalyzer::isDaughter(const reco::Candidate* Daughter, const reco::Candidate * particle)
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

bool GenParticleAnalyzer::isPdgMatch(const reco::Candidate * ancestor, const reco::Candidate * particle, int Id)
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
GenParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  int TauCount=0;
  int Tau_eCount=0;
  int Tau_muCount=0;
  int Tau_hadCount=0;

  using namespace edm;
  //using namespace edm;
  using namespace std;
  using namespace reco;

  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_, pruned);
  //Handle<pat::PackedGenParticle> packed;                                                                                                                                                                                                                                     
  Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_, packed);
  

  for(size_t i=0; i<(pruned.product())->size() ;i++)
    {
      /* if(fabs((*pruned)[i].pdgId())==36)
        {
	  const Candidate * PsuedoScalar = &(*pruned)[i];
	  
	  }*/
      
      if(fabs((*pruned)[i].pdgId())==15)
	{
	  cout<<" Tau found " <<endl;
	  cout << " pdgId: " << (*pruned)[i].pdgId() <<  " pt: " << (*pruned)[i].pt() << " eta: " << (*pruned)[i].eta() <<  " phi: " << (*pruned)[i].phi()<<endl;   
	  const Candidate * Tau = &(*pruned)[i];
	  //check for a pseudscalar mother, then check for a electron daughter.
	  //Look for an electron in the daughters, continue if there are any muons.
	  //increase a counter if there are multipleTaus with electron daughters and continue if this is about 2.
	  //have to find out hot to get rid of two taus decaying hardonically.
	  if (fabs(Tau->mother()->pdgId())==36)
	    {
	      ++TauCount;
	      cout << " Tau from PseudoScalar " <<endl;
	      unsigned  n=Tau->numberOfDaughters();
	      for ( size_t j =0; j < n ; j++)
		
		{	      const Candidate * Daughter=Tau->daughter(j);
		  cout<< "Daughter No:" << j << " pdgId " << Daughter->pdgId() <<endl;
		  if (fabs(Daughter->pdgId())==13)
		    {  ++Tau_muCount;
		  //continue;
		  cout<<" Tau_mu found "<<endl;
		    }
		  //continue;
		  if (fabs(Daughter->pdgId())==11)
		    {
		    ++Tau_eCount;
		  cout<<" Tau_e found "<<endl;
		    }
		  if(fabs(Daughter->pdgId()==111) || fabs(Daughter->pdgId()==211))
		    {
		    ++Tau_hadCount;
		  cout<<"Tau_had found "<<endl;
		  continue;
		    }
		}
	      //if (Tau_eCount == 2)
	      //     if (fabs(Daughter->pdgId())==13)
	      //continue;
	      
	      
	    }
	}
    }
  
  
  cout<< "Taus from Pseudoscalar :" << TauCount<<endl;
  if (Tau_eCount==2)
    cout<< "Tau_e Tau_e event" <<endl;
  // #ifdef THIS_IS_AN_EVENT_EXAMPLE
  //    Handle<ExampleData> pIn;
  //    iEvent.getByLabel("example",pIn);
  // #endif
  
  // #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  //    ESHandle<SetupData> pSetup;
  //    iSetup.get<SetupRecord>().get(pSetup);
  // #endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
GenParticleAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenParticleAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenParticleAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(GenParticleAnalyzer);
