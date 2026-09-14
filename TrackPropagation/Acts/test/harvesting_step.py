import FWCore.ParameterSet.Config as cms

process = cms.Process("HARVESTING")

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("DQMServices.Core.DQMStore_cfi")
process.load("DQMServices.Components.DQMFileSaver_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring("file:step3_inDQM.root")
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    wantSummary = cms.untracked.bool(False)
)

process.DQMStore.verbose = 0

from DQMServices.Core.DQMEDHarvester import DQMEDHarvester

process.dqmSaver = cms.EDAnalyzer(
    "DQMFileSaver",
    convention = cms.untracked.string("Offline"),
    workflow = cms.untracked.string("/Global/CMSSW_X_Y_Z/RECO"),
    dirName = cms.untracked.string("."),
    saveByRun = cms.untracked.int32(-1),
    saveAtJobEnd = cms.untracked.bool(True),
    forceRunNumber = cms.untracked.int32(1)
)

process.dqmsave_step = cms.EndPath(process.dqmSaver)

process.schedule = cms.Schedule(process.dqmsave_step)
