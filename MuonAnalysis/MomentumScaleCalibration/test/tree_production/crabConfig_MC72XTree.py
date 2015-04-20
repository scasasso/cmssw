from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
#config.General.requestName = 'Run2_DYToMuMu_PU40bx25_tsg_castor_PHYS14_25_V1-v2_test1'
config.General.requestName = 'DYToMuMu_M-50_Tune4C_13TeV-pythia8_Phys14DR-PU20bx25_castor_PHYS14_25_V1-v3'
config.General.workArea = 'crab_projects'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAODVBTFMuonsOnlyTreeMC_CRAB_cfg.py'
config.JobType.outputFiles = ['zmumuTree.root','0_zmumuHisto.root']

config.section_("Data")
#config.Data.inputDataset = '/DYToMuMu_M-50_Tune4C_13TeV-pythia8/Phys14DR-PU40bx25_tsg_castor_PHYS14_25_V1-v2/MINIAODSIM'
config.Data.inputDataset = '/DYToMuMu_M-50_Tune4C_13TeV-pythia8/Phys14DR-PU20bx25_castor_PHYS14_25_V1-v3/MINIAODSIM'
config.Data.inputDBS = 'global'
#config.Data.splitting = 'EventBased'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.outLFN = "/store/user/scasasso/MuScleFit_Phys14Tree"
#config.Data.outLFN = '/lustre/cms/lustre/store/user/emiglior/TEST/'
#config.Data.outLFN = '/store/caf/user/emiglior/Alignment/MuScleFit/Run2/tree/'

config.section_("Site")
#config.Site.storageSite = 'T2_IT_Bari'
config.Site.storageSite = 'T2_CH_CERN'
