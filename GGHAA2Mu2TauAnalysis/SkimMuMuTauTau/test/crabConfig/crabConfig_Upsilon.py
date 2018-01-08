from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'MiniAOD_Upsilon_AntiIsoMu2_TauDMMedIso_JAN1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../SkimSequence/NewMiniAOD_Skim.py'

config.Data.inputDataset = '/UpsilonMuMu_UpsilonPt6_TuneCUEP8M1_13TeV-pythia8-evtgen/dntaylor-2017-11-20_Skim_MuMuTauTau_80XSprin16_v1_MINIAODoutput-a6b69b8e27a72c8d16ba62a9c88868dd/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20
config.Data.outLFNDirBase = '/store/group/phys_higgs/HiggsExo/ktos'
config.Data.publication = True
config.Data.outputDatasetTag = 'MiniAOD_Upsilon_AntiIsoMu2_TauDMMedIso_JAN1'

config.Site.storageSite = 'T2_CH_CERN'
