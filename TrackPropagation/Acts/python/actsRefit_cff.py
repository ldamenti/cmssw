import FWCore.ParameterSet.Config as cms

## Load propagator
from TrackPropagation.Acts.ActsPropagators_cfi import *

from TrackingTools.TrackRefitter.TracksToTrajectories_cff import *

from SimG4Core.Application.g4SimHits_cfi import g4SimHits as _g4SimHits

# load this to do a track refit
from RecoTracker.TrackProducer.TrackRefitters_cff import *
from RecoVertex.V0Producer.generalV0Candidates_cff import *

ActsFitter = RKTrajectoryFitter.clone ( 
    ComponentName = cms.string('ActsFitter'),
    Propagator = cms.string('ActsPropagator') 
)

ActsSmoother = RKTrajectorySmoother.clone (
     ComponentName = cms.string('ActsSmoother'),
     Propagator = cms.string('ActsPropagator'),
 
     ## modify rescaling to have a more stable fit during the backward propagation
     # errorRescaling = cms.double(2.0) # use the default
)

ActsFitterSmoother = KFFittingSmootherWithOutliersRejectionAndRK.clone(
    ComponentName = cms.string('ActsFitterSmoother'),
    Fitter = cms.string('ActsFitter'),
    Smoother = cms.string('ActsSmoother'),

    # use the default value in KFFittingSmoother
    BreakTrajWith2ConsecutiveMissing = cms.bool(True),
    EstimateCut = cms.double(-1),
    HighEtaSwitch = cms.double(5),
    LogPixelProbabilityCut = cms.double(0),
    MaxFractionOutliers = cms.double(0.3),
    MaxNumberOfOutliers = cms.int32(3),
    MinDof = cms.int32(2),
    MinNumberOfHits = cms.int32(5),
    MinNumberOfHitsHighEta = cms.int32(5),
    NoInvalidHitsBeginEnd = cms.bool(True),
    NoOutliersBeginEnd = cms.bool(False),
    RejectTracks = cms.bool(False),
    appendToDataLabel = cms.string('')
)

## Different versions have different refitter cffs
from RecoTracker.TrackProducer.TrackRefitters_cff import *

# configure the refitter with the Acts propagator
# automatically uses the generalTracks collection as input
# ActsTrackRefitter = TrackActsRefitter.clone() # USES THE CUSTOM ALGORTHM 
ActsTrackRefitter = TrackRefitter.clone()
ActsTrackRefitter.src = cms.InputTag("generalTracks")
ActsTrackRefitter.Fitter = cms.string('ActsFitterSmoother')
ActsTrackRefitter.Propagator = cms.string('ActsPropagator')
ActsTrackRefitter.TrajectoryInEvent = cms.bool(True)

actsTrackRefit = cms.Sequence(ActsTrackRefitter)

