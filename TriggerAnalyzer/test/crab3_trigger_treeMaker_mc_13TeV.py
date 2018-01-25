from WMCore.Configuration import Configuration 
config = Configuration() 
config.section_("General") 
config.General.requestName = 'ttjets_LJ_23rdJan' ## change


config.section_("JobType") 
config.JobType.pluginName = 'Analysis' 
config.JobType.psetName = 'trigger_analyzer_mc_25ns_cfg.py' 
config.JobType.allowUndistributedCMSSW = False
#config.JobType.inputFiles = ['Fall15_25nsV2_DATA.db','Fall15_25nsV2_MC.db']


config.section_("Data") 
config.Data.inputDataset = '/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM'
#'/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
#'/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext4-v1/MINIAODSIM'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/' 
config.Data.splitting = 'FileBased' 
config.Data.unitsPerJob = 1
config.Data.totalUnits = 200  ## only for ttbar
config.Data.publication = True 
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/' 
config.Data.outputDatasetTag = '23rdJan_trigger_csvRWT_13TeV'
### change user Space 
config.Data.outLFNDirBase = '/store/user/lwming/' 
config.Data.ignoreLocality = True

config.section_("Site") 
config.Site.storageSite = 'T3_US_FNALLPC' 
