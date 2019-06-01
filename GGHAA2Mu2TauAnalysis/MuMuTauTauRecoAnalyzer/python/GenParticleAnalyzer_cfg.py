import FWCore.ParameterSet.Config as cms

process = cms.Process("select")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('select')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
    )

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(250) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/user/r/rhabibul/BoostedDiTau/CMSSW_8_0_30/src/GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/python/VMass_Pt_dR_selection_update.root'
        'file:/afs/cern.ch/work/r/rhabibul/RootFiles/Mass_input_11GeV_new.root'
        )
                            )

process.demo = cms.EDAnalyzer("GenParticleAnalyzer",
                              #electrons = cms.InputTag("Loose","MiniLooseElectron","USER"),
                              #Taus=cms.InputTag("slimmedTausElectronCleaned"),
                              #vertex=cms.InputTag("offlineSlimmedPrimaryVertices"),
                              #BM = cms.InputTag("offlineBeamSpot"),
                              pruned  = cms.InputTag("prunedGenParticles"),
                              packed =cms.InputTag("packedGenParticles"),
                              
                              
                              )


#process.TFileService = cms.Service("TFileService",
 #                                  fileName = cms.string('GenTaueSelect.root')
  #                                 )

process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("file:/afs/cern.ch/work/r/rhabibul/Prospectus/TestSelect.root")
                               )


process.p = cms.Path(process.demo)
process.ep=cms.EndPath(process.out)
