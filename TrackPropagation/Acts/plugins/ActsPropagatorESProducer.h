#ifndef TrackPropagators_ESProducers_ActsPropagatorESProducer_h
#define TrackPropagators_ESProducers_ActsPropagatorESProducer_h

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include <memory>

#include "ActsDataFormats/GeometrySurface/interface/CMSDetectorElement.h"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Geometry/Records/interface/ACTSTrackerGeometryRecord.h"

#include "TrackPropagation/Acts/interface/ActsPropagator.h"

/*
 * ActsPropagatorESProducer
 *
 * Produces an ActsPropagator for track propagation
 *
 */

class ActsPropagatorESProducer : public edm::ESProducer {
public:
  ActsPropagatorESProducer(const edm::ParameterSet &p);
  ~ActsPropagatorESProducer() override;

  std::unique_ptr<Propagator> produce(const TrackingComponentsRecord &);

private:

  edm::ESGetToken<TrackingGeometryWithDetEls, ACTSTrackerGeometryRecord> ACTStrkGeomInfoToken_;  
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;

  double mass_;
  double maxDPhi_;
  double ptMin_;
  PropagationDirection dir_;
  bool useRK_;
  bool useOldAnalPropLogic_;
  Acts::Logging::Level logLevel_;
};

#endif
