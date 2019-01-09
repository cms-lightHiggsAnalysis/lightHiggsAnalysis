# Import CMS python class definitions such as Process, Source, and EDProducer
import FWCore.ParameterSet.Config as cms



# Set up a process, named RECO in this case
process = cms.Process("USER")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.categories.append('USER')
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    limit = cms.untracked.int32(-1)
    )

#process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


#Include the number of events to run over
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )











# Configure the object that reads the input file
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/r/rhabibul/Prospectus/mA_11.root')
                            )

# Configure an object that produces a new data objec
process.Loose = cms.EDFilter("MiniElectronFilter",
                               vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               Rho = cms.InputTag("fixedGridRhoFastjetAll"),
                               electrons = cms.InputTag("slimmedElectrons"),
                               conv = cms.InputTag("reducedConversions"),
                               BM = cms.InputTag("offlineBeamSpot"),
                               #Tracks = cms.InputTag("electronGsfTracks"),
                               
                               
                               
                              
                               
                               
                               
                               )

# Configure the object that writes an output file
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string("file:/afs/cern.ch/work/r/rhabibul/Prospectus/mA_11_cleaned.root")
                               )



# Configure a path and endpath to run the producer and output modules
process.p = cms.Path(process.Loose)
process.ep = cms.EndPath(process.out)
#Needs to be there
