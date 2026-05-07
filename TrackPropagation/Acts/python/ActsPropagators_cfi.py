import FWCore.ParameterSet.Config as cms

# Default Acts propagator setup
ActsPropagator = cms.ESProducer("ActsPropagatorESProducer",
                                ComponentName = cms.string("ActsPropagator"),
                                PropagationDirection=cms.string("alongMomentum"),
                                LoggerLevel = cms.string("verbose"),
                                # Remove?
                                MaxDPhi = cms.double(1.6),
                                ptMin = cms.double(0.1),
                                useRungeKutta = cms.bool(False),
                                useOldAnalPropLogic = cms.bool(False),
                                Mass = cms.double(0.105658)
                                )

ActsPropagatorOppositeMomentum = cms.ESProducer("ActsPropagatorESProducer",
                                                ComponentName = cms.string("ActsPropagatorOppositeMomentum"),
                                                PropagationDirection=cms.string("oppositeMomentum"),
                                                LoggerLevel = cms.string("verbose"),
                                                # Remove?
                                                MaxDPhi = cms.double(1.6),
                                                ptMin = cms.double(0.1),
                                                useRungeKutta = cms.bool(False),
                                                useOldAnalPropLogic = cms.bool(False),
                                                Mass = cms.double(0.105658)
                                                )