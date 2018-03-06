#ifndef  Rochester_h
#define  Rochester_h

/**\class Rochester
 *
 *
 * Original Author:  Nicola De Filippis
 *
 *
 */

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "Rochester/RoccoR.cc"
#include "CLHEP/Random/RandomEngine.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

using namespace edm;
using namespace std;
using namespace reco;
using namespace math;

class Rochester : public edm::EDProducer {
 public:
  explicit Rochester(const edm::ParameterSet& );
  ~Rochester();

 private:
  virtual void produce(edm::Event&, const edm::EventSetup&);

  bool isData;	
  edm::EDGetTokenT<edm::View<pat::Muon> > muonLabel;
  string iName;
  RoccoR rc;
  edm::FileInPath RochesterDir_;
  edm::Service<edm::RandomNumberGenerator> rng;
};

#endif
