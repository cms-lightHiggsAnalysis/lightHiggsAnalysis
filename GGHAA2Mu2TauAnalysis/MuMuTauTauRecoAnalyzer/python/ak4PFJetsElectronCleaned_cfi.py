import FWCore.ParameterSet.Config as cms

recoElectronsForJetCleaning = cms.EDFilter('ElectronFilter',
                                           vertex = cms.InputTag("offlinePrimaryVerticesWithBS"),
                                           Rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                           electrons = cms.InputTag("gedGsfElectrons"),
                                           conv = cms.InputTag("conversions"),
                                           BM = cms.InputTag("offlineBeamSpot"),
                                           Tracks = cms.InputTag("electronGsfTracks"),
                                           Passcount =cms.uint32(1),
                                           )

ak4PFJetsElectronCleaned = cms.EDProducer(
    'ElectronCleanedJetProducer',
    jetSrc = cms.InputTag("ak4PFJets"),
    electronSrc = cms.InputTag("recoElectronsForJetCleaning"),
    pfCandSrc = cms.InputTag("particleFlow"),
    )
