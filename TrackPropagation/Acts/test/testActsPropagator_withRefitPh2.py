import FWCore.ParameterSet.Config as cms

from Configuration.ProcessModifiers.dd4hep_cff import dd4hep
from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9

process = cms.Process("RECO", Phase2C22I13M9, dd4hep)

# ============================================================
# Input / output
# ============================================================

globalPath = "/eos/user/l/ldamenti/DatasetFarm_files/Phase2/"
filename = "step2Ph2_Np20k_Pt0p5to5GeV_Pions_etaNeg4toPos4_PhiNegPitoPosPi.root"
outputPath = "/eos/user/l/ldamenti/TempFiles/RefStudy/Pions/Phase2/"
nEvents = -1

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(nEvents)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        "root://eosuser.cern.ch/" + globalPath + filename
    ),
    secondaryFileNames = cms.untracked.vstring()
)

# process.source.skipEvents = cms.untracked.uint32(60000)

# ============================================================
# Core CMSSW services / geometry / reconstruction
# ============================================================

process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")

process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.StandardSequences.L1Reco_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.RecoSim_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

# NOTE: D121 defines a different Detectors scenario (details here: https://github.com/cms-sw/cmssw/blob/master/Configuration/Geometry/python/dictRun4Geometry.py)
process.load('Configuration.Geometry.GeometryDD4hepExtendedRun4D121Reco_cff')

process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerTopology_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cff")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

process.load("Validation.RecoTrack.TrackValidation_cff")
from Validation.RecoTrack.MultiTrackValidator_cfi import *
process.load("SimGeneral.TrackingAnalysis.trackingParticleNumberOfLayersProducer_cff")

process.load("DQMServices.Core.DQMStoreNonLegacy_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic", "")

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    wantSummary = cms.untracked.bool(False),
    Rethrow = cms.untracked.vstring()
)

# ============================================================
# ACTS tracking geometry
# ============================================================

process.trackinGeoProducer = cms.ESProducer("TrackerGeomBuilderPh2WithActsESProducer",

    saveJsonfile = cms.untracked.bool(False),
    saveObjfile = cms.untracked.bool(False),
    outputObjFile = cms.untracked.string("testSlice.obj"),
    rangeZ = cms.untracked.vdouble(-1000, 1000),
    rangeR = cms.untracked.vdouble(0, 1200),

    saveSvgfile = cms.untracked.bool(False),
    outputSvgFile = cms.untracked.string("CMSPhase1Blueprint.svg"),

    mapMaterial = cms.untracked.bool(True),
    MaterialMaps = cms.untracked.string("/eos/user/l/ldamenti/MaterialMaps/PhaseII/MaterialMapsPh2_2e6Tracks_100binsUniform.json"),

    ActsLogLevel = cms.untracked.string("info")
)

# ============================================================
# ACTS refit
# ============================================================

from RecoTracker.TrackProducer.TrackRefitters_cff import TrackRefitter

process.load("TrackPropagation.Acts.actsRefit_withPCA_cff")

process.ActsSmoother.errorRescaling = cms.double(10)
process.ActsTrackRefitter.usePropagatorForPCA = cms.bool(False)
process.ActsTrackRefitter.GeometricInnerState = cms.bool(False)
process.ActsPropagator.LoggerLevel = cms.string("info")

# ============================================================
# CMSSW refit
# ============================================================

process.StandardTrackRefitter = TrackRefitter.clone(
    src = cms.InputTag("generalTracks"),
    NavigationSchool = cms.string(""),
    Fitter = cms.string("FlexibleKFFittingSmoother"),
    Propagator = cms.string("RungeKuttaTrackerPropagator")
)

# ============================================================
# TrackingParticle association / efficiency
# ============================================================

from SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi import quickTrackAssociatorByHits

process.StandardTrackAssociation = quickTrackAssociatorByHits.clone(
    label_tp = cms.InputTag("mix", "MergedTrackTruth"),
    label_tr = cms.InputTag("StandardTrackRefitter")
)

process.ActsTrackAssociation = quickTrackAssociatorByHits.clone(
    label_tp = cms.InputTag("mix", "MergedTrackTruth"),
    label_tr = cms.InputTag("ActsTrackRefitter")
)

process.ActsTrackValidator = multiTrackValidator.clone(
    associators = cms.untracked.VInputTag(
        cms.InputTag("ActsTrackAssociation")
    ),
    UseAssociators = cms.bool(True),

    dirName = cms.string("Tracking/ACTS/Tracking"),

    label = cms.VInputTag(
        cms.InputTag("ActsTrackRefitter")
    ),

    label_tp_effic = cms.InputTag("mix", "MergedTrackTruth"),
    label_tp_fake  = cms.InputTag("mix", "MergedTrackTruth"),

    doSummaryPlots = cms.untracked.bool(True),
    doSimPlots = cms.untracked.bool(True),
    doSimTrackPlots = cms.untracked.bool(True),
    doRecoTrackPlots = cms.untracked.bool(True),

    trackCollectionForDrCalculation = cms.InputTag("ActsTrackRefitter")
)

process.StandardTrackValidator = multiTrackValidator.clone(
    associators = cms.untracked.VInputTag(
        cms.InputTag("StandardTrackAssociation")
    ),
    UseAssociators = cms.bool(True),

    dirName = cms.string("Tracking/CMSSW/MuonTracking"),

    label = cms.VInputTag(
        cms.InputTag("StandardTrackRefitter")
    ),

    label_tp_effic = cms.InputTag("mix", "MergedTrackTruth"),
    label_tp_fake  = cms.InputTag("mix", "MergedTrackTruth"),

    doSummaryPlots = cms.untracked.bool(True),
    doSimPlots = cms.untracked.bool(True),
    doSimTrackPlots = cms.untracked.bool(True),
    doRecoTrackPlots = cms.untracked.bool(True),

    trackCollectionForDrCalculation = cms.InputTag("StandardTrackRefitter")
)

# ============================================================
# DQM output
# ============================================================

process.DQMoutput = cms.OutputModule(
    "DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string("DQMIO"),
        filterName = cms.untracked.string("")
    ),
    fileName = cms.untracked.string(outputPath + "step3_inDQM.root"),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# ============================================================
# Misc required for premixing/playback input
# ============================================================

process.mix.playback = True
process.mix.digitizers = cms.PSet()

for a in process.aliases:
    delattr(process, a)

process.RandomNumberGeneratorService.restoreStateLabel = cms.untracked.string("randomEngineStateProducer")

# ============================================================
# Sequences / paths / schedule
# ============================================================
process.Customval = cms.Sequence(
    process.tpClusterProducer *
    process.trackingParticleNumberOfLayersProducer *
    process.StandardTrackAssociation *
    process.ActsTrackAssociation *
    process.ActsTrackValidator  *
    process.StandardTrackValidator
)

process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)

process.refitAndAnalysis_step = cms.Path(
    process.actsTrackRefit *
    process.StandardTrackRefitter
)

process.CustomValidationPath = cms.Path(process.Customval)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.L1Reco_step,
    process.reconstruction_step,
    process.recosim_step,
    process.refitAndAnalysis_step,
    process.CustomValidationPath,
    process.DQMoutput_step
)

# Early deletion only
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)


process.MessageLogger.suppressWarning = cms.untracked.vstring(
    "earlyGeneralTracks"
)
