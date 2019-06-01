// -*- C++ -*-
//
// Package:    GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer
// Class:      TauEleFilter
// 
/**\class TauEleFilter TauEleFilter.cc GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/plugins/TauEleFilter.cc

 Description: [one line class summary]
Selects Taus and electrons that pass all the cuts and then finds the best pair(based on maximum visible mass or the maximum DI object Pt and puts them to new collections
     [Notes on implementation]
*/
//
// Original Author:  Redwan Habibullah
//         Created:  Mon, 16 Jul 2018 19:57:49 GMT
//
//


// system include files
#include <memory>
#include <algorithm>
#include <map>

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

typedef math::XYZPoint Point;


//
// class declaration
//

class TauEleFilter : public edm::stream::EDFilter<> {
   public:
      explicit TauEleFilter(const edm::ParameterSet&);
      ~TauEleFilter();

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

  edm::EDGetTokenT<pat::ElectronCollection> electronSrc_;
  edm::EDGetTokenT<pat::TauCollection> TauSrc_;
  edm::EDGetTokenT<reco::VertexCollection> vtx_;
  edm::EDGetTokenT<reco::BeamSpot> thebs_;
  edm::EDGetTokenT<reco::GenParticleCollection> particleSrcToken_;
  unsigned int TauReg=0;
  unsigned int TauMode=0;
  unsigned int TauMed=0;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
///
// constructors and destructor
//
TauEleFilter::TauEleFilter(const edm::ParameterSet& iConfig):
  
  electronSrc_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  TauSrc_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("Taus"))),
  vtx_ (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertex"))),
  thebs_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BM"))),
  particleSrcToken_ (consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("particleSrc")))
{
  //now do what ever initialization is needed
  produces<pat::ElectronCollection>("PassedElectron");
  produces<pat::TauCollection>("PassedTau");
  produces<pat::ElectronCollection>("CutElectron");
  produces<pat::TauCollection>("CutTaus");
  produces<pat::ElectronCollection>("ReuseElectron");
  produces<pat::TauCollection>("ReuseTaus");
  
  produces<std::vector<pat::Electron>>("HighVmassElectron");
  produces<std::vector<pat::Tau>>("HighVmassTaus");
  produces<std::vector<pat::Electron>>("HighPtElectrons");
  produces<std::vector<pat::Tau>>("HighPtTaus");
  produces<std::vector<pat::Electron>>("dRElectrons");
  produces<std::vector<pat::Tau>>("dRTaus");
}


TauEleFilter::~TauEleFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TauEleFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electronSrc_,electrons);
  unique_ptr<pat::ElectronCollection> MatchElectrons(new pat::ElectronCollection);
  unique_ptr<pat::ElectronCollection> TransElectrons(new pat::ElectronCollection);
  unique_ptr<pat::ElectronCollection> RepeatElectrons(new pat::ElectronCollection);
  unique_ptr <vector<pat::Electron>>MassEle(new vector<pat::Electron>);                                                                                                                                                                                        
  unique_ptr <vector<pat::Electron>>PtEle(new vector<pat::Electron>);
  unique_ptr <vector<pat::Electron>>dREle(new vector<pat::Electron>);

  Handle<pat::TauCollection> Taus;
  iEvent.getByToken(TauSrc_,Taus);
  unique_ptr<pat::TauCollection> MatchTaus(new pat::TauCollection);
  unique_ptr<pat::TauCollection> TransTaus(new pat::TauCollection);
  unique_ptr<pat::TauCollection> RepeatTaus(new pat::TauCollection);

  unique_ptr <vector<pat::Tau>>MassTau(new vector<pat::Tau>);
  unique_ptr <vector<pat::Tau>>PtTau(new vector<pat::Tau>);
  unique_ptr <vector<pat::Tau>>dRTau(new vector<pat::Tau>);

  Handle<reco::VertexCollection> Vertex;
  iEvent.getByToken(vtx_,Vertex);
  edm::Handle<reco::BeamSpot> thebs;
  iEvent.getByToken(thebs_,thebs);
  edm::Handle<reco::GenParticleCollection> particles;
  iEvent.getByToken(particleSrcToken_, particles);

  const reco::Vertex& pv=*Vertex->begin();
  Point p;
  p=pv.position();
  double dR=999999;
  double dR_m=999999;
  double dR_m_min=999999;
  // double dR_p=999999;
  double Vmass=999999;                                                                                                                                                                                                                                                       
  double Vmax=0;
  double Ptmax=0;
  double Pt_event=999999;
  // int countele=0;                                                                                                                                                                                                                                                           
  //int counttau=0;                                                                                                                                                                                                                                                            
  double dz_tau;
  double dxy_tau;
  double dz_ele;
  double dxy_ele;
  //int nDoublect= 0;
  
  vector <unsigned int> TausAlreadyInPlot;
  vector <unsigned int> ElectronsAlreadyInPlot;
  vector<pat::Electron> MatchedElectrons;
  vector<pat::Tau>MatchedTaus;
  vector<pat::Electron>PtElectrons;
  vector<pat::Tau>PtTaus;

  unsigned int EBcount = 0;
  unsigned int EEcount = 0;
  unsigned int Passcount= 0;
  unsigned int Taucount = 0;
  unsigned int Modecount =0;
  unsigned int Medcount =0;
  unsigned int int_ele=0;
  unsigned int int_tau=0;
  unsigned int int_ele_pt=0;
  unsigned int int_tau_pt=0;

  //General Selection loop

  int ele_idx =0;
  for(pat::ElectronCollection::const_iterator iele = electrons->begin() ; iele !=electrons->end() ; ++iele,++ele_idx)
    {

      dz_ele=iele->gsfTrack()->dz(p);
      dxy_ele=iele->gsfTrack()->dxy(p);
      pat::ElectronRef inputElectronRef(electrons,ele_idx);
     
      int tau_idx = 0;
      
      for(pat::TauCollection::const_iterator itau = Taus->begin() ; itau !=Taus->end() ; ++itau,++tau_idx)


	{
	  ++Taucount;
	  pat::Tau outputTau(*itau);
	  pat::TauRef inputTauRef(Taus, tau_idx);

	  pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(itau->leadChargedHadrCand().get());
	  dz_tau=packedLeadTauCand->dz(p);
	  dxy_tau=packedLeadTauCand->dxy(p);
	  
	  dR = reco::deltaR(*iele, *itau);
          //Vmass=abs((iele->p4() + itau->p4()).mass());
	  if(itau->tauID("decayModeFinding"))
	    {
	      ++Modecount;
	    }
	  if(itau->tauID("byIsolationMVArun2v1DBoldDMwLTraw") > -0.5)
	    {
	      ++Medcount;
	    }
          if( (itau->pt()>10) && abs((itau->eta())<2.3) && (itau->tauID("decayModeFinding")) && (itau->tauID("byIsolationMVArun2v1DBoldDMwLTraw") >-0.5) && (dxy_tau<0.2) && (dz_tau < 0.5) && (iele->isEB()) && (dxy_ele<0.05) &&(dxy_ele <0.10) )
            {


	      TausAlreadyInPlot.push_back(inputTauRef.key());
              ElectronsAlreadyInPlot.push_back(inputElectronRef.key());
              //MatchedElectrons.push_back(*iele);
	      MatchElectrons->push_back(*iele);
	      //MatchedTaus.push_back(*itau);
	      MatchTaus->push_back(*itau);
	      
	      
	      
              //cout<<"KeyEB_Tau:"<<inputTauRef.key()<<" KeyEB_ele:"<<inputElectronRef.key()<<endl;
	      ++EBcount;
	      if((dR < 0.8) && (dR> 0.05))
		{
		  cout<<"Fill Histo_EB"<<endl;
		  MatchedElectrons.push_back(*iele);
		  MatchedTaus.push_back(*itau);
		  PtElectrons.push_back(*iele);
		  PtTaus.push_back(*itau);
		  TransElectrons->push_back(*iele);
                  TransTaus->push_back(*itau);
		  RepeatElectrons->push_back(*iele);
                  RepeatTaus->push_back(*itau);

		  
		  
		}
	      
	      
	    }
	  
	  if((itau->pt()>10) && abs((itau->eta())<2.3) && (itau->tauID("decayModeFinding")) && (itau->tauID("byIsolationMVArun2v1DBoldDMwLTraw") >-0.5) && (dxy_tau<0.2) && (dz_tau < 0.5) && (iele->isEE()) && (dxy_ele<0.10) &&(dz_ele <0.20) )
            {
              //dR = reco::deltaR(*iele, *itau);                                                                                                                                                                                                                               
              //DeltaR->Fill(dR);                                                                                                                                                                                                                                                            //DeltaR->SetFillColor(kRed);                                                                                                                                                                                                                                    
              TausAlreadyInPlot.push_back(inputTauRef.key());
              ElectronsAlreadyInPlot.push_back(inputElectronRef.key());
	      MatchElectrons->push_back(*iele);
              MatchTaus->push_back(*itau);
	      
	      ++EEcount;
	      if((dR < 0.8) && (dR > 0.05))
                {
                  cout<<"Fill Histo_EE"<<endl;
		  MatchedElectrons.push_back(*iele);
		  MatchedTaus.push_back(*itau);
		  PtElectrons.push_back(*iele);
		  PtTaus.push_back(*itau);		
		  TransElectrons->push_back(*iele);
		  TransTaus->push_back(*itau);
		  RepeatElectrons->push_back(*iele);
                  RepeatTaus->push_back(*itau);

		}

	      
	    }
	  
	}
      
      
      
    }
  
  
  //dR  based Selection loop for Virtual mass
  for(unsigned i=0;i<MatchedElectrons.size();i++)
    {
      for(unsigned j=0;j<MatchedTaus.size();j++)
	{
	  
	  // Vmass=abs(((MatchedElectrons[i]).p4() +(MatchedTaus[j]).p4()).mass());
	  dR_m=reco::deltaR((MatchedElectrons[i]).eta(),(MatchedElectrons[i]).phi(),(MatchedTaus[i]).eta(),(MatchedTaus[i]).phi());
	  //if ((Vmass > Vmax))-> Need to ask
	  if ((dR_m < dR_m_min))
	    {
	      //Vmax=Vmass;
	      dR_m_min=dR_m;
	      
	      
	      int_ele=i;
	      int_tau=i;
	    }
	  
	  
	}
      
      
    }

  if((MatchedElectrons.size()>0) && (MatchedTaus.size()>0))
    {
      dREle->push_back(MatchedElectrons[int_ele]);
      dRTau->push_back(MatchedTaus[int_tau]);
    }









  for(unsigned i=0;i<MatchedElectrons.size() && i<MatchedTaus.size();i++)                                                                                                                                                                                                                            
    {                                                                                                                                                                                                                                                                          
                                                                                                                                                                                                                                                                               
          Vmass=abs(((MatchedElectrons[i]).p4() +(MatchedTaus[i]).p4()).mass());                                                                                                                                                                                            
       
          if ((Vmass > Vmax))                                                                                                                                                                                                                                 
          
            {                                                                                                                                                                                                                                                                  
              Vmax=Vmass;                                                                                                                                                                                                                                                    
	     
                                                                                                                                                                                                                                                                               
              int_ele=i;                                                                                                                                                                                                                                                       
              int_tau=i;
	    }
	  
    }
  
  
  
  




  if((MatchedElectrons.size()>0) && (MatchedTaus.size()>0))
    {
      MassEle->push_back(MatchedElectrons[int_ele]);                                                                                                                                                                                                                       
      MassTau->push_back(MatchedTaus[int_tau]);
    }
  
  
  //Pt based selection loop for Virtual mass
  for(unsigned i=0;i<PtElectrons.size();i++)
    {
      for(unsigned j=0;j<PtTaus.size();j++)
	{
	  Pt_event=(((PtElectrons[i]).p4()+(PtTaus[i]).p4()).pt());
	 
	  if(Ptmax<Pt_event)
	    {
	      Ptmax=Pt_event;
	      
	      int_ele_pt=i;
	      int_tau_pt=i;
	    }
	}
      
    }
  
  if((PtElectrons.size()>0) && (PtTaus.size()>0))
    {
      PtEle->push_back(PtElectrons[int_ele_pt]);
      PtTau->push_back(PtTaus[int_tau_pt]);
    }
  
  
  
  
  
  iEvent.put(move(RepeatElectrons), "ReuseElectron");
  iEvent.put(move(RepeatTaus), "ReuseTaus");
  iEvent.put(move(TransElectrons), "CutElectron");
  iEvent.put(move(TransTaus), "CutTaus");
  iEvent.put(move(MatchElectrons), "PassedElectron");
  iEvent.put(move(MatchTaus), "PassedTau");
  iEvent.put(move(MassEle),"HighVmassElectron");
  iEvent.put(move(MassTau),"HighVmassTaus");
  iEvent.put(move(PtEle),"HighPtElectrons");
  iEvent.put(move(PtTau),"HighPtTaus");
  iEvent.put(move(dREle),"dRElectrons");
  iEvent.put(move(dRTau),"dRTaus");




 
  if (Taucount!=0)
    ++TauReg;
  if (Modecount!=0)
    ++TauMode;
  if (Medcount !=0)
    ++TauMed;
  cout<<"Event with Cleaned Taus " << TauReg <<endl;
  cout<<"Event with DecayMode Passing Taus " << TauMode <<endl;
  cout<<"Event with Med isolation Taus " << TauMed <<endl;

  if(Passcount < EBcount + EEcount)
    return true;
  else return false;
  
}


// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
TauEleFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
TauEleFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
TauEleFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
TauEleFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TauEleFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TauEleFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TauEleFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TauEleFilter);
