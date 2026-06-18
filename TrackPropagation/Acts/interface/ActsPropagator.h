#ifndef TrackPropagation_ActsPropagator_h
#define TrackPropagation_ActsPropagator_h

#include <memory>
#include "DataFormats/GeometryCommonDetAlgo/interface/DeepCopyPointerByClone.h"

// CMS includes
// - Propagator
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "DataFormats/GeometrySurface/interface/PlaneBuilder.h"

#include "TrackPropagation/Acts/interface/ActsConcretePropagator.h"
#include "TrackPropagation/Acts/interface/ComputeFreeJacobian.h"

// - Acts
#include "ActsDataFormats/GeometrySurface/interface/CMSMagneticFieldProvider.hpp"
#include "ActsDataFormats/GeometrySurface//interface/CMSDetectorElement.h"
#include "ActsDataFormats/GeometrySurface/interface/TrackingGeometryWithDetEls.h"
#include "Geometry/TrackerGeometryBuilder/interface/SurfaceConverters.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "Acts/Definitions/Units.hpp"

using DetElVect = std::vector<std::shared_ptr<Acts::CMSDetectorElement>>;
using MyConcretePropagator = ConcretePropagator<Acts::Propagator<Acts::EigenStepper<>, Acts::Navigator>>;
using TsosWP = std::pair<TrajectoryStateOnSurface, double>;

class MagneticField;

class ActsPropagator : public Propagator {
public:

  ActsPropagator(PropagationDirection dir,
                  const float mass,
                  const MagneticField* mf = nullptr,
                  const float maxDPhi = 1.6,
                  bool useRungeKutta = false,
                  float ptMin = -1.,
                  bool useOldGeoPropLogic = true,
                  //
                  std::shared_ptr<TrackingGeometryWithDetEls> TrkandDetEls = nullptr,
                  const TrackerGeometry* trkGeo_cmssw = nullptr,
                  const Acts::Logging::Level& lvl = Acts::Logging::Level::INFO);

  ~ActsPropagator() override;


  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const FreeTrajectoryState &,
                                                                const Plane &) const override;

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const FreeTrajectoryState &,
                                                                const Cylinder &) const override;

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const TrajectoryStateOnSurface &,
                                                                const Plane &) const override;

  std::pair<TrajectoryStateOnSurface, double> propagateWithPath(const TrajectoryStateOnSurface &,
                                                                const Cylinder &) const override;

  FreeTrajectoryState propagateWithPathToPerigeeInsideBP(const TrajectoryStateOnSurface& tsos) const;

  ActsPropagator *clone() const override { return new ActsPropagator(*this); }

  const MagneticField *magneticField() const override { return field; }


  private:

  // helpfull methods to be used during the propagation:
  std::shared_ptr<const Acts::Surface> inflateStartPlaneIfRectOrTrap(const Acts::Surface& targetSurf, const Acts::GeometryContext& gctx, double marginMm) const;
  Acts::BoundMatrix convertCovCMSSWtoACTS(const TrajectoryStateOnSurface& tsos, AlgebraicSymMatrix55 cov_cmssw) const;
  AlgebraicSymMatrix55 convertCovACTStoCMSSW(const Surface& surf, double phi, double theta, const Acts::BoundMatrix& cov_acts) const;
  LocalTrajectoryParameters GetLocalTrajectoryParameters(const Acts::BoundTrackParameters& pActs, const Surface& surf) const;
  std::shared_ptr<const Acts::Surface> buildTargetSurf(const Surface& cmssw_surf) const;
  TsosWP tsosWithPathFromActs(const Acts::BoundTrackParameters& start_param, const Surface& surfTarget, Acts::Direction actsDir) const;
  TrajectoryStateOnSurface boundFreeTrajectoryState(const FreeTrajectoryState &fts) const;



  const MagneticField *field;
  PropagationDirection propDir_;
  std::shared_ptr<TrackingGeometryWithDetEls> trkGeo_and_DetEls_;
  const TrackerGeometry* trkGeo_cmssw_;
  DeepCopyPointerByClone<MyConcretePropagator> ConcProp;
  Acts::Direction actsPropDir_ = Acts::Direction::Forward();

  // Logger
  Acts::Logging::Level m_Level;
  std::shared_ptr<const Acts::Logger> m_logger;
  const Acts::Logger& logger() const { return *m_logger; }

  mutable int nTot  = 0, nNoSurf = 0, nNoProp = 0, nNoCov = 0;

};

#endif
