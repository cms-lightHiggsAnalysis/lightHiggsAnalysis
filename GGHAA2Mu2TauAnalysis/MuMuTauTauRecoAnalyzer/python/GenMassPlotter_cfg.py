import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Demo')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
    )

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/r/rhabibul/BoostedDiTau/CMSSW_8_0_30/src/GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/test/fifth.root'
        )
                            )

process.demo = cms.EDAnalyzer("GenMassPlotter",
                              #electrons = cms.InputTag("Loose","MiniLooseElectron","USER"),
                              #Taus=cms.InputTag("slimmedTausElectronCleaned"),
                              #vertex=cms.InputTag("offlineSlimmedPrimaryVertices"),
                              #BM = cms.InputTag("offlineBeamSpot"),
                              pruned  = cms.InputTag("prunedGenParticles"),
                              packed =cms.InputTag("packedGenParticles"),
                              
                              
                              )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('GenVisinvDmode2.root')
                                   )


process.p = cms.Path(process.demo)

