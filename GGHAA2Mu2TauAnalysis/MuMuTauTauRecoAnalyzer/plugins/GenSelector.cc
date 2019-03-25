
// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer
// Class:      GenSelector
// 
/**\class GenSelector GenSelector.cc GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/plugins/GenSelector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Redwan Habibullah
//         Created:  Mon, 23 Jul 2018 02:20:47 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

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


#include "DataFormats/Math/interface/deltaR.h"
#include "TLatex.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"

#include"DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TAttFill.h"
#include "Math/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

//
// class declaration
//

class GenSelector : public edm::stream::EDFilter<> {
   public:
      explicit GenSelector(const edm::ParameterSet&);
      ~GenSelector();
  bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);
  bool isDaughter(const reco::Candidate * Daughter, const reco::Candidate * particle);
  bool isPdgMatch(const reco::Candidate * Ancestor, const reco::Candidate * particle, int Id);


      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

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
GenSelector::GenSelector(const edm::ParameterSet& iConfig):


  prunedGenToken_ (consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  packedGenToken_ (consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed")))
{

   //now do what ever initialization is needed
  // produces<edm::View<pat::PackedGenParticle> >("TaueleGenParticle");
  produces<pat::PackedGenParticleCollection>("TaueleGenParticle");
  
}


GenSelector::~GenSelector()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

bool GenSelector::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle)
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


bool GenSelector::isDaughter(const reco::Candidate* Daughter, const reco::Candidate * particle)
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


bool GenSelector::isPdgMatch(const reco::Candidate * ancestor, const reco::Candidate * particle, int Id)
{


  if((fabs(ancestor->pdgId())==Id) && (fabs(particle->pdgId())==Id) && (ancestor==particle)) return true;
  for(size_t i=0;i<particle->numberOfMothers();i++)
    {
      if(isPdgMatch(ancestor,particle->mother(i),Id)) return true;
    }
  return false;
}





















// ------------ method called on each new Event  ------------
bool
GenSelector::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace std;
  using namespace reco;

  Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_, pruned);
  //Handle<pat::PackedGenParticle> packed;                                                                                                                                                                                                         
  Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_, packed);
  //unique_ptr<edm::View<pat::PackedGenParticle> > MatchPackedGens(new View<pat::PackedGenParticle>);
  // std::unique_ptr<pat::PackedGenParticleCollection>   MatchPackedGens(new pat::PackedGenParticleCollection);
auto out=std::unique_ptr<pat::PackedGenParticleCollection> (new pat::PackedGenParticleCollection);

  int Tau_e_count=0;
  //int Tau_had_count=0;
  int Passcount=0;

  for(size_t i=0; i<(pruned.product())->size();i++)
    {
      
      
      //cout<<"PdgID: "<<(*pruned)[i].pdgId()<<endl;                                                                                                                                                                                                                           
      //cout<<" index: "<<i<<endl;                                                                                                                                                                                                                                             
      
      //for(size_t k=0;k<(pruned.product())->size();k++)
	//{ 
	  if(fabs((*pruned)[i].pdgId())==36)
	    {
	      //TLorentzVector invariant_v4_Tau_e;
	      //TLorentzVector invariant_v4_Tau_had;
	      TLorentzVector invariant_v4_Tau_mu;
	      const Candidate * PsuedoScalar = &(*pruned)[i];
	      //cout << "PdgID: " << PsuedoScalar->pdgId() <<endl;                                                                                                                                                                                                                 
	      //cout << "found daugthers: "<< endl;                                                                                                                                                                                                                                
	      //cout<<"PdgID: "<<(*pruned)[i].pdgId()<<endl;                                                                                                                                                                                                                       
	      
	      for(size_t j=0; j<(packed.product())->size() ;j++)
		{
		  
		  
		  const Candidate * motherInPrunedCollection = (*packed)[j].mother(0) ;
		  const GenParticleRef motherInPrunedCollectionRef=(*packed)[j].motherRef();
		  //Adding for Tau_e                                                                                                                                                                                                                                               
		  if((motherInPrunedCollection != nullptr) && (isAncestor( PsuedoScalar , motherInPrunedCollection)) /* (fabs((*packed)[j].pdgId())==11)*/ && (fabs((motherInPrunedCollection->pdgId()))==15) && (motherInPrunedCollection->status()!=1) )
		    {
		      //cout << "Tau_e_PdgID: " << (*packed)[j].pdgId() <<endl;
		      //cout << "Tau_e_Pt: " << (*packed)[j].pt() <<endl;
		      //cout<<"Tau_e_eta: " <<(*packed)[j].eta()<<endl;
		      //cout<< " key: " << motherInPrunedCollectionRef.key()<<endl;                                                                                                                                                                                                		      //matchcount_ele++;            
		      //MatchPackedGens->push_back((*packed)[j]);
		      //different way to have the whoel object
		      if( (fabs((*packed)[j].pdgId())==13) )
			
			
			{
			  
			  invariant_v4_Tau_mu +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
			  
			  
			  
			}
		      
		      
		      if((fabs((*packed)[j].pdgId())==11) && ((*packed)[j].statusFlags().isFirstCopy()) && (invariant_v4_Tau_mu.X()==0) && (invariant_v4_Tau_mu.Y()==0) && (invariant_v4_Tau_mu.Z()==0) && (invariant_v4_Tau_mu.T()==0) && (invariant_v4_Tau_mu.Mag()==0) )

			{

		      const auto obj=packed->at(j);
		      pat::PackedGenParticle newObj=obj; 
		      //pat::PackedGenParticleCollection Obj = std::<vector<pat::PackedGenParticle>>;
		      //MatchPackedGens->push_back(newObj);
		      out->push_back(newObj);

		      ++Tau_e_count;
		      cout<<"Tau_e_loop"<<endl;
			}
		      //invariant_v4_Tau_e +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
		      //cout<<"value_e: "<< invariant_v4_Tau_e                                                                                                                                                                                                                     
		    }
		  
		  //Adding for Tau_had                                                                                                                                                                                                                                      
		  
		  /*for(size_t k=0;k<(pruned.product())->size();k++)  
		    {
		  if((motherInPrunedCollection != nullptr) && !(isPdgMatch(&(*pruned)[k],motherInPrunedCollection,11)) && !(isPdgMatch(&(*pruned)[k],motherInPrunedCollection,13)) && (isAncestor( PsuedoScalar , motherInPrunedCollection)) && (fabs((*packed)[j].pdgId())!=11) && (fabs((*packed)[j].pdgId())!=13) && (motherInPrunedCollection->status()!=1) && ((*packed)[j].statusFlags().isFirstCopy()) && (isPdgMatch(&(*pruned)[k],motherInPrunedCollection,15)) && (fabs((*packed)[j].pdgId())!=12) && (fabs((*packed)[j].pdgId())!=16)  && (fabs((*packed)[j].pdgId())!=14)  )
		    {
		      
		      //cout<<"Tau_had_pdgID: "<<(*packed)[j].pdgId() <<endl;
		      //invariant_v4_Tau_had +=TLorentzVector((*packed)[j].p4().X(),(*packed)[j].p4().Y(),(*packed)[j].p4().Z(),(*packed)[j].p4().T());
		      //cout<<"value_had: "<< invariant_v4_Tau_had<<endl;                                                                                                                                                                                                          
		      //++Tau_had_count;
		      cout<<"Tau_had_loop"<<endl;
		      
		      
		      
		    }
		  
		  //MvisibleGenTau->Fill(fabs((invariant_v4_Tau_e+invariant_v4_Tau_had).M()));
		  //MinvariantGen->Fill((invariant_v4_Tau_had).M());                                                                                                                                                                                                               
		  //MvisibleGenTau->SetFillColor(kPink+2);
		  
		  } */
		  
		}
	      
	    }
	  //}
	  
    }
  
  
  
  
  
















  /*using namespace edm;
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
  */


  //iEvent.put(move(MatchPackedGens),"TaueleGenParticle");
  iEvent.put(move(out),"TaueleGenParticle");
  
  
  //if (Passcount<Tau_e_count+Tau_had_count)
  if(Passcount<Tau_e_count)
    return true;
  else return false;
  
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
GenSelector::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
GenSelector::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
GenSelector::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
GenSelector::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GenSelector::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GenSelector::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenSelector::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenSelector);
