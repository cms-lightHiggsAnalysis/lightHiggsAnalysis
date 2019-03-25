# Import CMS python class definitions such as Process, Source, and EDProducer
import FWCore.ParameterSet.Config as cms



# Set up a process, named RECO in this case
process = cms.Process("complete")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('complete')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
    )

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


#Include the number of events to run over
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )











# Configure the object that reads the input file
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/r/rhabibul/RootFiles/Mass_input_11GeV_new.root')
                            )

# Configure an object that produces a new data objec
process.Loose = cms.EDFilter("MiniElectronFilter",
                               vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               Rho = cms.InputTag("fixedGridRhoFastjetAll"),
                               electrons = cms.InputTag("slimmedElectrons"),
                               conv = cms.InputTag("reducedConversions"),
                               BM = cms.InputTag("offlineBeamSpot")
                               #Tracks = cms.InputTag("electronGsfTracks"),
                             )

process.High = cms.EDFilter("TauEleFilter",
                              electrons = cms.InputTag("Loose","MiniLooseElectron","complete"),
                              Taus=cms.InputTag("slimmedTausElectronCleaned"),
                              #Taus=cms.InputTag("slimmedTaus"),
                              vertex=cms.InputTag("offlineSlimmedPrimaryVertices"),
                              BM = cms.InputTag("offlineBeamSpot"),
                              particleSrc  = cms.InputTag("prunedGenParticles"),


                              )


process.DiMuon= cms.EDFilter("DiMuonFilter",  
                            vertex=cms.InputTag("offlineSlimmedPrimaryVertices"), 
                            muonSrc =cms.InputTag("slimmedMuons"),
                            bits = cms.InputTag("TriggerResults","","HLT"),
                            prescales = cms.InputTag("patTrigger"),
                            objects = cms.InputTag("selectedPatTrigger")

                              )

#process.CrossClean=cms.EDFilter("CrossCleaner",
#                               electrons = cms.InputTag("High","PassedElectron","complete"),
#                               Taus=cms.InputTag("High","PassedTau","complete"),
 #                              Muons1 = cms.InputTag("DiMuon","LeadingMuon","complete"),
  #                             Muons2=cms.InputTag("DiMuon","TrailingMuon","complete"),
   #                            Mu1TaudRCut=cms.double(0.8),
    #                           Mu2TaudRCut=cms.double(0.8),
     #                          Mu1EledRCut=cms.double(0.4),
      #                         Mu2EledRCut=cms.double(0.4)


#)
process.Analyze = cms.EDAnalyzer("PrimaryAnalyzer",
                                 electrons = cms.InputTag("Loose","MiniLooseElectron","complete"),
                                 Taus=cms.InputTag("slimmedTausElectronCleaned"),
                                 #Taus=cms.InputTag("slimmedTaus"), 
                                 vertex=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 BM = cms.InputTag("offlineBeamSpot"),
                                 particleSrc  = cms.InputTag("prunedGenParticles"),
                                 #Tau_mass_select =cms.InputTag("CrossClean","CrossCleanedHighMassTaus","complete"),
                                 Tau_mass_select =cms.InputTag("High","HighVmassTaus","complete"),
                                 Ele_mass_select=cms.InputTag("High","HighVmassElectron","complete"),
                                 Tau_pt_select =cms.InputTag("High","HighPtTaus","complete"),
                                 Ele_pt_select =cms.InputTag("High","HighPtElectrons","complete"),
                                 Ele_dR_select=cms.InputTag("High","dRElectrons","complete"),
                                 Tau_dR_select=cms.InputTag("High","dRTaus","complete"),
                                 #packed=cms.InputTag("demo","TaueleGenParticle","Demo"),
                                 Muons1 = cms.InputTag("DiMuon","LeadingMuon","complete"),
                                 Muons2=cms.InputTag("DiMuon","TrailingMuon","complete") 
                                 )


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


# Configure the object that writes an output file
#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string("file:/afs/cern.ch/work/r/rhabibul/Prospectus/Combined.root")
#                               )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('file:/afs/cern.ch/work/r/rhabibul/RootFiles/Combined_11GeV_Cleaned_test_PreCleaned.root')
                                   )



# Configure a path and endpath to run the producer and output modules
process.p = cms.Path(process.Loose*process.High*process.DiMuon*process.Analyze)
