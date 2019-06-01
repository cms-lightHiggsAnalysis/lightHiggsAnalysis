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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))










# Configure the object that reads the input file
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/r/rhabibul/RootFiles/Mass_input_11GeV_new_lowered.root')
                            )
process.GenModeFilter = cms.EDFilter("GenModeFilter",
                              
                                     pruned  = cms.InputTag("prunedGenParticles"),
                                     packed =cms.InputTag("packedGenParticles"),
)
process.PrimaryMuonAnalyzer=cms.EDAnalyzer("PrimaryMuonAnalyzer",
                               muonSrc =cms.InputTag("slimmedMuons"),
                                                            )
                                     
                                     
# Configure an object that produces a new data objec
process.LooseFilter = cms.EDFilter("MiniElectronFilter",
                               vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               Rho = cms.InputTag("fixedGridRhoFastjetAll"),
                               electrons = cms.InputTag("slimmedElectrons"),
                               conv = cms.InputTag("reducedConversions"),
                               BM = cms.InputTag("offlineBeamSpot")
                               #Tracks = cms.InputTag("electronGsfTracks"),
                             )
process.ElectronAnalyzer = cms.EDAnalyzer("MiniLooseAnalyzer",
                                          electrons = cms.InputTag("LooseFilter","MiniLooseElectron","complete"),
                                          Rho = cms.InputTag("fixedGridRhoFastjetAll"),
                                          #pruned  = cms.InputTag("prunedGenParticles"),
                                          #packed =cms.InputTag("packedGenParticles"),
                                          
)
process.PrimaryElectronAnalyzer = cms.EDAnalyzer("PrimaryElectronAnalyzer",
                                                 electrons = cms.InputTag("slimmedElectrons"),
                                                 pruned  = cms.InputTag("prunedGenParticles"),
                                                 packed =cms.InputTag("packedGenParticles"),


)




process.PairFilter = cms.EDFilter("TauEleFilter",
                              electrons = cms.InputTag("LooseFilter","MiniLooseElectron","complete"),
                              Taus=cms.InputTag("slimmedTausElectronCleaned"),
                              #Taus=cms.InputTag("slimmedTaus"),
                              vertex=cms.InputTag("offlineSlimmedPrimaryVertices"),
                              BM = cms.InputTag("offlineBeamSpot"),
                              particleSrc  = cms.InputTag("prunedGenParticles"),


                              )
process.PairAnalyzer = cms.EDAnalyzer("PostPrimarySelectionAnalyzer",
                                      electrons = cms.InputTag("PairFilter","PassedElectron","complete"), 
                                      Taus=cms.InputTag("PairFilter","PassedTau","complete"),
                                      Cutelectrons = cms.InputTag("PairFilter","CutElectron","complete"),
                                      CutTaus=cms.InputTag("PairFilter","CutTaus","complete"),
                                      pruned  = cms.InputTag("prunedGenParticles"),
                                      packed =cms.InputTag("packedGenParticles"),

                                 )
process.PrimaryTauAnalyzer=cms.EDAnalyzer("PrimaryTauAnalyzer",
                                          Taus=cms.InputTag("slimmedTausElectronCleaned"),
                                          vertex=cms.InputTag("offlineSlimmedPrimaryVertices"),
                                          
                                          )
process.DiMuFilter= cms.EDFilter("DiMuonFilter",  
                             vertex=cms.InputTag("offlineSlimmedPrimaryVertices"), 
                             muonSrc =cms.InputTag("slimmedMuons"),
                             bits = cms.InputTag("TriggerResults","","HLT"),
                             prescales = cms.InputTag("patTrigger"),
                             objects = cms.InputTag("selectedPatTrigger"),
                             #electrons = cms.InputTag("High","CutElectron","complete"),
                             #Taus=cms.InputTag("High","CutTaus","complete")

                              )


process.DiMuAnalyzer = cms.EDAnalyzer("DiMuonAnalyzer",
                                     Muons1 = cms.InputTag("DiMuFilter","LeadingMuon","complete"),
                                     Muons2=cms.InputTag("DiMuFilter","TrailingMuon","complete"),
                                     #electrons = cms.InputTag("High","CutElectron","complete"),
                                     #Taus=cms.InputTag("High","CutTaus","complete"),
                                     
                                     #Taus=cms.InputTag("slimmedTaus"),                                                                                                                                                                                                            
                                 )

process.CrossCleanFilter=cms.EDFilter("CrossCleaner",
                               electrons = cms.InputTag("PairFilter","CutElectron","complete"),
                               Taus=cms.InputTag("PairFilter","CutTaus","complete"),
                               Muons1 = cms.InputTag("DiMuFilter","LeadingMuon","complete"),
                               Muons2=cms.InputTag("DiMuFilter","TrailingMuon","complete"),
                               Mu1TaudRCut=cms.double(0.8),
                               Mu2TaudRCut=cms.double(0.8),
                               Mu1EledRCut=cms.double(0.4),
                               Mu2EledRCut=cms.double(0.4)


                                )
process.FinalAnalyzer = cms.EDAnalyzer("MuMuTaueTauhadAnalyzer",
                                       electrons = cms.InputTag("PairFilter","PassedElectron","complete"),
                                       Taus=cms.InputTag("PairFilter","PassedTau","complete"),
                                       Cutelectrons = cms.InputTag("PairFilter","ReuseElectron","complete"),
                                       CutTaus=cms.InputTag("PairFilter","ReuseTaus","complete"),
                                       Tau_mass_select =cms.InputTag("CrossCleanFilter","CrossCleanedHighMassTaus","complete"),
                                       Ele_mass_select=cms.InputTag("CrossCleanFilter","CrossCleanedHighMassElectron","complete"),
                                       Tau_pt_select =cms.InputTag("PairFilter","HighPtTaus","complete"),
                                       Ele_pt_select =cms.InputTag("PairFilter","HighPtElectrons","complete"),
                                       Ele_dR_select=cms.InputTag("PairFilter","dRElectrons","complete"),
                                       Tau_dR_select=cms.InputTag("PairFilter","dRTaus","complete"),
                                       Muons1 = cms.InputTag("CrossCleanFilter","CrossCleanedLeadingMuon","complete"),
                                       Muons2=cms.InputTag("CrossCleanFilter","CrossCleanedTrailingMuon","complete"),
                                       pruned  = cms.InputTag("prunedGenParticles"),
                                       packed=cms.InputTag("demo","TaueleGenParticle","Demo")
                                 )


process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )




process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('file:/afs/cern.ch/work/r/rhabibul/RootFiles/Combined_11GeV_test_steps_swap_Gen_Test.root')
                                   )



# Configure a path and endpath to run the producer and output modules
#process.p = cms.Path(process.Loose*process.PairFilter*process.Primary*process.DiMuon*process.DiMuAnalyze*process.CrossClean*process.Analyze)
process.p = cms.Path(process.GenModeFilter*process.PrimaryMuonAnalyzer*process.DiMuFilter*process.DiMuAnalyzer*process.PrimaryElectronAnalyzer*process.LooseFilter*process.ElectronAnalyzer*process.PrimaryTauAnalyzer*process.PairFilter*process.PairAnalyzer*process.CrossCleanFilter*process.FinalAnalyzer)
