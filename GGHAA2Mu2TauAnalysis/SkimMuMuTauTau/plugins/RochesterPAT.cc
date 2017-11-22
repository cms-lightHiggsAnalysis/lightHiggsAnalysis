/**\class RochesterPAT.cc
 *
 * Original Author:  Nicola De Filippis 
 *
 */

#include "GGHAA2Mu2TauAnalysis/SkimMuMuTauTau/plugins/RochesterPAT.h"

#include "FWCore/Framework/interface/ESHandle.h"

// Muons:
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

// Candidate handling
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/Common/interface/AssociationVector.h"
#include "DataFormats/Common/interface/ValueMap.h"

// system include files
#include <Math/VectorUtil.h>
#include <memory>
#include <vector>


using namespace edm;
using namespace std;
using namespace math;

// constructor
RochesterPAT::RochesterPAT(const edm::ParameterSet& pset) {
  isData           = pset.getParameter<bool>("isData");
  muonLabel        = consumes<edm::View<pat::Muon> >(pset.getParameter<edm::InputTag>("muonCollection"));
  iName = "RochesterMuPAT";
  produces<std::vector<pat::Muon> >(iName); 
  RochesterDir_=pset.getParameter<edm::FileInPath>("fp");
  std::string rochCorrDataDirPath=RochesterDir_.fullPath();
  rochCorrDataDirPath.erase(rochCorrDataDirPath.length()-10);
  rc=RoccoR(rochCorrDataDirPath);
}


RochesterPAT::~RochesterPAT() {
}


//
// member functions
//
void RochesterPAT::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {



  // muons
  auto_ptr<std::vector<pat::Muon> > Gmuon( new std::vector<pat::Muon> );
  edm::Handle<edm::View<pat::Muon> > muons;
  edm::View<pat::Muon>::const_iterator mIter;    
  iEvent.getByToken(muonLabel, muons);
  CLHEP::HepRandomEngine& engine= rng->getEngine(iEvent.streamID());

		 


  double corrPt=0.;
  double factor=1.0;
    
  // Loop over muons
  for (mIter = muons->begin(); mIter != muons->end(); ++mIter ) 
  {
    pat::Muon* calibmu = mIter->clone(); 
    reco::Candidate::PolarLorentzVector p4Polar_;
    if (isData)
      factor = rc.kScaleDT(mIter->charge(), mIter->pt(), mIter->eta(), mIter->phi(), 0, 0);
    else
    {
        double u1 = engine.flat();
        double u2 = engine.flat();
        factor= rc.kScaleAndSmearMC(mIter->charge(), mIter->pt(), mIter->eta(), mIter->phi(), mIter->bestTrack()->hitPattern().trackerLayersWithMeasurement(), u1, u2, 0, 0);
    }//else
    corrPt=mIter->pt()*factor; 
    p4Polar_ = reco::Candidate::PolarLorentzVector(corrPt, mIter->eta(), mIter->phi(), mIter->mass());
    calibmu->setP4(p4Polar_);
    Gmuon->push_back( *calibmu );

  }

  
  const string iName = "RochesterMuPAT";
  iEvent.put( Gmuon, iName );

}
DEFINE_FWK_MODULE(RochesterPAT);

