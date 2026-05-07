import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.dd4hep_cff import dd4hep

process = cms.Process("DumpGeometry", dd4hep)

process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer",
  usePhase2Stacks = cms.bool(False)
)

# Load DD4hep geometry
process.load("Configuration.Geometry.GeometryDD4hepExtended2023_cff")     # PHASE 1

# Tracker Geometry and Tracker Topology 
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerTopology_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cff")

process.TrackerAdditionalParametersPerDetESModule = cms.ESProducer("TrackerAdditionalParametersPerDetESModule")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "auto:phase1_2023_realistic", "") 

# ====== Logging ======
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

# ===== Construct the ACTS Tracking Geometry =====
process.trackinGeoProducer = cms.ESProducer("TrackerGeomBuilderWithActsESProducer", 
    # Option to save the sensitive surfaces in a JSON file
    saveJsonfile   = cms.untracked.bool(False),
    # Options to save the detector elements in an OBJ file
    saveObjfile    = cms.untracked.bool(False),
    outputObjFile  = cms.untracked.string("testSlice.obj"),
    rangeZ         = cms.untracked.vdouble(-1000, 1000),  # Min, Max (mm)
    rangeR         = cms.untracked.vdouble(0, 1200),      # Min, Max (mm) 
    # Options to save the Tracker blueprint on an SVG file
    saveSvgfile    = cms.untracked.bool(False),
    outputSvgFile  = cms.untracked.string("CMSPhase1Blueprint.svg"),
    # Option to map the material from a JSON file
    mapMaterial    = cms.untracked.bool(False), # NOTE: needs to be false while producing the material maps
    MaterialMaps   = cms.untracked.string(""),

    ActsLogLevel    = cms.untracked.string("info")
)

# ===== Create the Material Maps =====
process.createMaterialFile = cms.EDProducer("ActsJsonMaterialMapProducer",
    G4InputFile = cms.untracked.string("/eos/user/l/ldamenti/G4MaterialFiles/geant4MaterialFile_1e6Tracks_ActsUnits_WithMyFormula.root"),
    OutputFile  = cms.untracked.string("/eos/user/l/ldamenti/MaterialMaps/MaterialMaps_1e6Tracks_NewBP_FullLayers.json"),
    Nevents = cms.untracked.int32(1000000),    # NOTE: MAX value = number of tracks of the G4 file
    ActsLogLevel    = cms.untracked.string("info")
)

# ====== Paths ======
process.p = cms.Path(process.createMaterialFile)



