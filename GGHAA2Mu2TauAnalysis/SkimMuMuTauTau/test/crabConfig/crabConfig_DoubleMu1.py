from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'MiniAOD_DoubleMu1_DiMu_TauDM_NOV16'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../SkimSequence/NewMiniAOD_Skim.py' 

config.Data.inputDataset = '/DoubleMuon/dntaylor-2017-11-03_Skim_MuMuTauTau_v4_MINIAODoutput-58fa585557bb6dc4894eaeabc5791e7a/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.outLFNDirBase = '/store/group/phys_higgs/HiggsExo/ktos/' #'/store/user/ktos/' #'/store/group/phys_higgs/HiggsExo/ktos/' #'/store/user/ktos/'
config.Data.publication = True
config.Data.outputDatasetTag = 'MiniAOD_DoubleMu1_DiMu_TauDM_NOV16'

config.Site.storageSite = 'T2_CH_CERN'
