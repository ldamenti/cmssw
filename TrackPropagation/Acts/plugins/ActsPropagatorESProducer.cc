#include "ActsPropagatorESProducer.h"

#include "TrackPropagation/Acts/interface/ActsPropagator.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include <memory>
#include <string>

using namespace edm;

namespace {
  PropagationDirection stringToDirection(std::string const& iName) {
    PropagationDirection dir = alongMomentum;

    if (iName == "oppositeToMomentum")
      dir = oppositeToMomentum;
    if (iName == "alongMomentum")
      dir = alongMomentum;
    if (iName == "anyDirection")
      dir = anyDirection;
    return dir;
  }

  Acts::Logging::Level stringToLogLevel(std::string const& iName) {
    Acts::Logging::Level logLevel = Acts::Logging::Level::INFO;

    if (iName == "info")
      logLevel = Acts::Logging::Level::INFO;
    if (iName == "verbose")
      logLevel = Acts::Logging::Level::VERBOSE;
    return logLevel;
  }

}  // namespace

ActsPropagatorESProducer::ActsPropagatorESProducer(const edm::ParameterSet &p)
    : mass_(p.getParameter<double>("Mass")),
      maxDPhi_(p.getParameter<double>("MaxDPhi")),
      ptMin_(p.getParameter<double>("ptMin")),
      dir_(stringToDirection(p.getParameter<std::string>("PropagationDirection"))),
      useRK_(p.getParameter<bool>("useRungeKutta")),
      useOldAnalPropLogic_(p.getParameter<bool>("useOldAnalPropLogic")),
      logLevel_(stringToLogLevel(p.getParameter<std::string>("LoggerLevel")))
{
  auto cc = setWhatProduced(this, p.getParameter<std::string>("ComponentName"));
  ACTStrkGeomInfoToken_ = cc.consumes();
  trackerGeomToken_ = cc.consumes();
  magFieldToken_ = cc.consumes();
}

ActsPropagatorESProducer::~ActsPropagatorESProducer() {}

std::unique_ptr<Propagator> ActsPropagatorESProducer::produce(const TrackingComponentsRecord &iRecord) {

  std::cout << "USING ACTS PROPAGATOR" << std::endl;

  return std::make_unique<ActsPropagator>(
      dir_, mass_,  &iRecord.get(magFieldToken_), maxDPhi_, useRK_, ptMin_, useOldAnalPropLogic_, 
      std::make_shared<TrackingGeometryWithDetEls>(iRecord.get(ACTStrkGeomInfoToken_)),
      &iRecord.get(trackerGeomToken_),
      logLevel_);
}
