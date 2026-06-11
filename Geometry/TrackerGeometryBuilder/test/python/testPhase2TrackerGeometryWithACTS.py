import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.dd4hep_cff import dd4hep
from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9

process = cms.Process("DumpGeometry", Phase2C22I13M9, dd4hep)

process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer",
  usePhase2Stacks = cms.bool(True)
)

# Load DD4hep geometry
# NOTE: D121 defines a different Detectors scenario (details here: https://github.com/cms-sw/cmssw/blob/master/Configuration/Geometry/python/dictRun4Geometry.py)
process.load('Configuration.Geometry.GeometryDD4hepExtendedRun4D121Reco_cff')

# Tracker Geometry and Tracker Topology 
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerTopology_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.ProcessModifiers.dd4hep_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase2_realistic", "") 

# ====== Logging ======
process.MessageLogger.cerr.threshold = 'DEBUG'
process.MessageLogger.cerr.FwkReport.reportEvery = 1

# ====== Input ======
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

# ===== Construct the ACTS Tracking Geometry =====
process.trackinGeoProducer = cms.ESProducer("TrackerGeomBuilderPh2WithActsESProducer", 
    # Option to save the sensitive surfaces in a JSON file
    saveJsonfile   = cms.untracked.bool(False),
    # Options to save the detector elements in an OBJ file
    saveObjfile    = cms.untracked.bool(False),
    outputObjFile  = cms.untracked.string("TBPS_l1_d3.obj"),
    rangeZ         = cms.untracked.vdouble(246, 294),  # Min, Max (mm)
    rangeR         = cms.untracked.vdouble(213, 300),      # Min, Max (mm) 
    # Options to save the Tracker blueprint on an SVG file
    saveSvgfile    = cms.untracked.bool(True),
    outputSvgFile  = cms.untracked.string("CMSPhase2Blueprint.svg"),
    # Option to map the material from a JSON file
    mapMaterial    = cms.untracked.bool(False),
    MaterialMaps   = cms.untracked.string("MaterialMaps.json"),

    ActsLogLevel    = cms.untracked.string("verbose")
)

process.get = cms.EDAnalyzer("EventSetupRecordDataGetter",
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('ACTSTrackerGeometryRecord'),
        data = cms.vstring('TrackingGeometryWithDetEls'),
    )),
    verbose = cms.untracked.bool(True))

process.p = cms.Path(process.get)