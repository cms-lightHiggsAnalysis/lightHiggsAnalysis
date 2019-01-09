import FWCore.ParameterSet.Config as cms

process = cms.Process("Select")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('Select')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
    )

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.options=cms.untracked.PSet(SkipEvent = cms.untracked.vstring('ProductNotFound'))  

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
       #'file:/afs/cern.ch/user/r/rhabibul/BoostedDiTau/CMSSW_8_0_30/src/GGHAA2Mu2TauAnalysis/MuMuTauTauRecoAnalyzer/python/VMass_Pt_selection.root'
       #'file:/afs/cern.ch/work/r/rhabibul/Prospectus/mA_11_filtered.root'
        "file:/afs/cern.ch/work/r/rhabibul/Prospectus/mA_19_dR_Test.root"
        #'file:/afs/cern.ch/work/r/rhabibul/RootFiles/VMass_Pt_dR_selection.root'
        )
                            
                            )



process.select = cms.EDAnalyzer("PrimaryAnalyzer",
                                electrons = cms.InputTag("Loose","MiniLooseElectron","USER"),
                                Taus=cms.InputTag("slimmedTausElectronCleaned"),
                                vertex=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                BM = cms.InputTag("offlineBeamSpot"),
                                particleSrc  = cms.InputTag("prunedGenParticles"),
                                Tau_mass_select =cms.InputTag("demo","HighVmassTaus","Demo"),
                                Ele_mass_select=cms.InputTag("demo","HighVmassElectron","Demo"),
                                Tau_pt_select =cms.InputTag("demo","HighPtTaus","Demo"),
                                Ele_pt_select =cms.InputTag("demo","HighPtElectrons","Demo"),
                                #packed=cms.InputTag("demo","TaueleGenParticle","Demo"),
                              )


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('file:/afs/cern.ch/work/r/rhabibul/RootFiles/mA_19_dR_result_a.root')
                                   )


process.p = cms.Path(process.select)

