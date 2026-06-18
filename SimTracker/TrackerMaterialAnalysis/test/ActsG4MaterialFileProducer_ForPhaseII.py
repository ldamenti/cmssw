#!/usr/bin/env cmsRun

import FWCore.ParameterSet.Config as cms

from Configuration.ProcessModifiers.dd4hep_cff import dd4hep
from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9

# process = cms.Process('Geometry', dd4hep)
process = cms.Process('Geometry', Phase2C22I13M9, dd4hep)

## Define the geometry to be used
# NOTE: D121 defines a different Detectors scenario (details here: https://github.com/cms-sw/cmssw/blob/master/Configuration/Geometry/python/dictRun4Geometry.py)
process.load('Configuration.Geometry.GeometryDD4hepExtendedRun4D121Reco_cff')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')

## MC Related stuff
process.load('Configuration.StandardSequences.Generator_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

### Loading 10GeV neutrino gun generator
process.load("SimTracker.TrackerMaterialAnalysis.single10GeVNeutrino_cfi")

### Load vertex generator w/o smearing
from Configuration.StandardSequences.VtxSmeared import VtxSmeared
process.load(VtxSmeared['NoSmear'])

# detector simulation (Geant4-based) with tracking material accounting 
process.load("SimTracker.TrackerMaterialAnalysis.ActsTrackingMaterialProducer_cff")
#For some reason now neutrino are no longer tracked, so we need to force it.
process.trackingMaterialProducer.StackingAction.TrackNeutrino = True
process.trackingMaterialProducer.Generator.HepMCProductLabel = cms.InputTag("generatorSmeared")
process.trackingMaterialProducer.Watchers[0].ActsTrackingMaterialProducer.SelectedVolumes = cms.vstring('BEAM_1', 'BEAM_2', 'Tracker_1')
process.trackingMaterialProducer.Watchers[0].ActsTrackingMaterialProducer.ActsOutputFileName = cms.string("/eos/user/l/ldamenti/G4MaterialFiles/PhaseII/ph2_geant4MaterialFile_4e6Tracks.root")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(2000000)
)

# Input source
process.source = cms.Source("EmptySource")

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring(
        'drop *',                                                       # drop all objects
        'keep MaterialAccountingTracks_trackingMaterialProducer_*_*'),  # but the material accounting informations
    fileName = cms.untracked.string('file:material.root')
)

process.path = cms.Path(process.generator
                        * process.VtxSmeared
                        * process.generatorSmeared
                        * process.trackingMaterialProducer)
