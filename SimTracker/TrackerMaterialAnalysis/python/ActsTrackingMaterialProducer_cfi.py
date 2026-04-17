import FWCore.ParameterSet.Config as cms
import copy
from SimG4Core.Application.g4SimHits_cfi import *

trackingMaterialProducer = copy.deepcopy(g4SimHits)
trackingMaterialProducer.Generator.HepMCProductLabel = 'generatorSmeared'
trackingMaterialProducer.Physics.type = 'SimG4Core/Physics/DummyPhysics'
trackingMaterialProducer.Physics.DummyEMPhysics = True
trackingMaterialProducer.Physics.CutsPerRegion = False
trackingMaterialProducer.UseMagneticField = False
trackingMaterialProducer.StackingAction.TrackNeutrino = True  

trackingMaterialProducer.Watchers = cms.VPSet(cms.PSet(
    ActsTrackingMaterialProducer = cms.PSet( PrimaryTracksOnly = cms.bool(True),
        txtOutFile = cms.untracked.string('VolumesZPosition.txt'),
        hgcalzfront = cms.double(3190.5),
        SelectedVolumes = cms.vstring('BEAM', 'Tracker')
    ),
    type = cms.string('ActsTrackingMaterialProducer') 
))
