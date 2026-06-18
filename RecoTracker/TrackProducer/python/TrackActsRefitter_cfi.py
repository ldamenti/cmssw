import FWCore.ParameterSet.Config as cms

TrackActsRefitter = cms.EDProducer("TrackActsRefitter",

    src = cms.InputTag("generalTracks"),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    Fitter = cms.string("KFFittingSmootherWithOutliersRejectionAndRK"),
    TTRHBuilder = cms.string("WithAngleAndTemplate"),
    AlgorithmName = cms.string("undefAlgorithm"),

    constraint = cms.string(""),
    srcConstr = cms.InputTag(""),

    useHitsSplitting = cms.bool(False),
    TrajectoryInEvent = cms.bool(True),
    GeometricInnerState = cms.bool(False),

    NavigationSchool = cms.string("SimpleNavigationSchool"),
    MeasurementTracker = cms.string(""),
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
)