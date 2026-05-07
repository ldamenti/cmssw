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

process.actsEfficiencyClient = DQMEDHarvester("DQMGenericClient",

    subDirs = cms.untracked.vstring("Tracking/ACTS/*"),

    efficiency = cms.vstring(
        "effic_vs_pT 'ACTS efficiency vs p_{T};p_{T} [GeV];efficiency' num_assoc(simToReco)_pT num_simul_pT",
        "effic_vs_eta 'ACTS efficiency vs #eta;#eta;efficiency' num_assoc(simToReco)_eta num_simul_eta",
        "effic_vs_phi 'ACTS efficiency vs #phi;#phi;efficiency' num_assoc(simToReco)_phi num_simul_phi",
        "effic_vs_vertpos 'ACTS efficiency vs production radius;r [cm];efficiency' num_assoc(simToReco)_vertpos num_simul_vertpos",
        "effic_vs_zpos 'ACTS efficiency vs production z;z [cm];efficiency' num_assoc(simToReco)_zpos num_simul_zpos",
        "effic_vs_hit 'ACTS efficiency vs hits;hits;efficiency' num_assoc(simToReco)_hit num_simul_hit",
        "effic_vs_layer 'ACTS efficiency vs layers;layers;efficiency' num_assoc(simToReco)_layer num_simul_layer",

        "fakerate_vs_eta 'ACTS fake rate vs #eta;#eta;fake rate' num_assoc(recoToSim)_eta num_reco_eta fake",
        "fakerate_vs_pT 'ACTS fake rate vs p_{T};p_{T} [GeV];fake rate' num_assoc(recoToSim)_pT num_reco_pT fake",
    ),

    resolution = cms.vstring(),
    outputFileName = cms.untracked.string("")
)

process.dqmSaver = cms.EDAnalyzer(
    "DQMFileSaver",
    convention = cms.untracked.string("Offline"),
    workflow = cms.untracked.string("/Global/CMSSW_X_Y_Z/RECO"),
    dirName = cms.untracked.string("."),
    saveByRun = cms.untracked.int32(-1),
    saveAtJobEnd = cms.untracked.bool(True),
    forceRunNumber = cms.untracked.int32(1)
)

process.trackingHarvesting_step = cms.Path(process.actsEfficiencyClient)

process.dqmsave_step = cms.EndPath(process.dqmSaver)

process.schedule = cms.Schedule(
    process.trackingHarvesting_step,
    process.dqmsave_step
)
