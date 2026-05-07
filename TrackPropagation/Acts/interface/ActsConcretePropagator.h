#ifndef _ACTSCONCRETEPROPAGATOR_H_
#define _ACTSCONCRETEPROPAGATOR_H_

#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/MaterialInteractor.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"

#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "Acts/Visualization/PlyVisualization3D.hpp"

#include <Acts/Surfaces/RectangleBounds.hpp>
#include <Acts/Surfaces/TrapezoidBounds.hpp>

#include <iostream>
#include <iostream>
#include <limits>
#include <cmath>
#include <optional>
#include <iterator>

// DEBUG
static std::string surfaceType(const Acts::Surface& s) {
  using Acts::Surface;
  switch (s.type()) {
    case Surface::Plane:    return "Plane";
    case Surface::Cylinder: return "Cylinder";
    case Surface::Disc:     return "Disc";
    case Surface::Cone:     return "Cone";
    case Surface::Perigee:  return "Perigee";
    case Surface::Straw:    return "Straw";
    default:                return "Other";
  }
}

template <typename ParametersT>
inline Acts::Direction chooseDirectionTowardTarget(const Acts::GeometryContext& gctx,
                                                   const Acts::Surface& targetSurf,
                                                   const ParametersT& startParameters) {
  const auto pos = startParameters.position(gctx);
  const auto u   = startParameters.direction();

  // Usa l'overload con extend=true se nella tua build è disponibile
  auto interF = targetSurf.intersect(gctx, pos,  u);  // avanti lungo +u
  auto interB = targetSurf.intersect(gctx, pos, -u);  // avanti lungo -u

  auto firstPath = [](auto& inter) -> std::optional<double> {
    auto it = std::begin(inter);
    auto ie = std::end(inter);
    if (it == ie) return std::nullopt;
    return it->pathLength();
  };

  const auto sF = firstPath(interF);  // pathLength lungo +u
  const auto sB = firstPath(interB);  // pathLength lungo -u

  const bool okF = sF && (*sF > 0.0);
  const bool okB = sB && (*sB > 0.0);

  // Preferisci il verso che rende la target "davanti"
  if (okB && !okF) {
    std::cout << "Based on pathlenth, direction chosen: BACKWARD" << std::endl;
    return Acts::Direction::Backward(); // vai lungo -u
  }
  if (okF && !okB) {
    std::cout << "Based on pathlenth, direction chosen: FORWARD" << std::endl;
    return Acts::Direction::Forward();  // vai lungo +u
  }

  // Se entrambi sono "davanti" (raro), scegli quello più vicino
  if (okF && okB) {
    return (*sB < *sF) ? Acts::Direction::Backward() : Acts::Direction::Forward();
  }

  // Nessuna intersezione "davanti": parallelo o target non raggiungibile da qui
  return Acts::Direction::Forward();
}

template <typename ParametersT>
inline void dumpTargetGeometryPre(const Acts::GeometryContext& gctx,
                           const Acts::Surface& target,
                           const ParametersT& startPars,
                           const char* tag) {
  const auto x = startPars.position(gctx);
  const auto u = startPars.direction(); // unit direction

  const auto n = target.normal(gctx, x, u);

  const double d = n.dot(x - target.center(gctx));
  const double ndotu = n.dot(u);

  std::cout << ("ActsPropDbg")
    << tag
    << " x=(" << x.x() << "," << x.y() << "," << x.z() << ")"
    << " u=(" << u.x() << "," << u.y() << "," << u.z() << ")"
    << " n=(" << n.x() << "," << n.y() << "," << n.z() << ")"
    << " signedDist=" << d
    << " ndotu=" << ndotu
    << " s~" << ((std::abs(ndotu) > 1e-12) ? (-d/ndotu) : 1e99) << std::endl;
}

inline void dumpBoundsCheck(const Acts::GeometryContext& gctx,
                     const Acts::Surface& target,
                     const Acts::Vector3& gpos,
                     const char* tag) {
  auto lp = target.globalToLocal(gctx, gpos, Acts::Vector3::Zero());

  std::cout << "ActsPropDbg"
    << tag << " gpos=(" << gpos.x() << "," << gpos.y() << "," << gpos.z() << ")" << std::endl;

  if (!lp.ok()) {
    std::cout << "ActsPropDbg" << tag << " globalToLocal FAILED" << std::endl;
    return;
  }

  const auto l = lp.value();
  const bool inside = target.bounds().inside(l); 

  std::cout << "ActsPropDbg"
    << tag
    << " local=(" << l[0] << "," << l[1] << ")"
    << " insideBounds=" << inside << std::endl;
}

// End of DEBUG

// DEBUG:
static std::shared_ptr<const Acts::Surface>
inflateTargetPlaneIfRectOrTrap(const Acts::Surface& targetSurf,
                               const Acts::GeometryContext& gctx,
                               double marginMm, 
                               double trapFactor) {
  // 1) Deve essere un piano
  if (targetSurf.type() != Acts::Surface::SurfaceType::Plane) {
    return nullptr;
  }

  // 2) Stesso transform del target originale
  const auto trf = targetSurf.transform(gctx);

  // 3) Discrimina per tipo bounds
  const auto bType = targetSurf.bounds().type();

  // RectangleBounds
  if (bType == Acts::SurfaceBounds::BoundsType::eRectangle) {
    const auto* rect = dynamic_cast<const Acts::RectangleBounds*>(&targetSurf.bounds());
    if (!rect) return nullptr;

    const double hx = rect->halfLengthX() + marginMm;
    const double hy = rect->halfLengthY() + marginMm;

    auto newBounds = std::make_shared<Acts::RectangleBounds>(hx, hy);
    return std::static_pointer_cast<const Acts::Surface>(Acts::Surface::makeShared<Acts::PlaneSurface>(trf, newBounds));

  }

  // TrapezoidBounds
  if (bType == Acts::SurfaceBounds::BoundsType::eTrapezoid) {
    const auto* trap = dynamic_cast<const Acts::TrapezoidBounds*>(&targetSurf.bounds());
    if (!trap) return nullptr; // safety

    // In Acts TrapezoidBounds è descritto da:
    // halfXnegY, halfXposY, halfY  (lato corto a y negativo). :contentReference[oaicite:1]{index=1}
    const double halfXnegY =
        trap->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY)* trapFactor;
    const double halfXposY =
        trap->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY)* trapFactor;
    const double halfY =
        trap->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthY)* trapFactor;
    const double ang =
        trap->get(Acts::TrapezoidBounds::BoundValues::eRotationAngle);

    auto newBounds = std::make_shared<Acts::TrapezoidBounds>(halfXnegY, halfXposY, halfY, ang);

    return std::static_pointer_cast<const Acts::Surface>(Acts::Surface::makeShared<Acts::PlaneSurface>(trf, newBounds));
  }

  // Altri planar bounds (EllipseBounds, InfiniteBounds, ecc.) -> non gestiti
  return nullptr;
}
// End DEBUG

struct PlanePrecheckResult {
  double pathLength = 0.;
  Acts::Vector2 local{0., 0.};
  Acts::Vector3 global{0., 0., 0.};
  bool insideBounds = false;
};

template <typename ParametersT>
std::optional<PlanePrecheckResult>
intersectInfinitePlanePrecheck(const Acts::GeometryContext& gctx,
                               const Acts::Surface& targetSurf,
                               const ParametersT& startParameters,
                               Acts::Direction propDir) {
  Acts::Vector3 pos = startParameters.position(gctx);
  Acts::Vector3 dir = startParameters.direction();

  // If propagating backward, flip the direction for the geometric test
  if (propDir == Acts::Direction::Backward()) {
    dir = -dir;
  }

  auto intersections = targetSurf.intersect(
      gctx,
      pos,
      dir,
      Acts::BoundaryTolerance::Infinite());

  for (const auto& isec : intersections) {
    if (isec.status() != Acts::IntersectionStatus::reachable) {
      continue;
    }

    if (isec.pathLength() < 0.) {
      continue;
    }

    const Acts::Vector3 ipos = isec.position();

    auto lpos = targetSurf.globalToLocal(gctx, ipos, dir);
    if (!lpos.ok()) {
      continue;
    }

    const bool inside = targetSurf.insideBounds(
        *lpos,
        Acts::BoundaryTolerance::None());

    return PlanePrecheckResult{
        isec.pathLength(),
        *lpos,
        ipos,
        inside};
  }

  return std::nullopt;
}

using BoundParameters = Acts::GenericBoundTrackParameters<Acts::ParticleHypothesis>;

// using PropagationOutput = std::pair<PropagationSummary, Acts::RecordedMaterial>;
struct PropagationAlgorithm_Config {
  /// Switch the logger to sterile - for timing measurements
  bool sterileLogger = false;
  /// Modify the behavior of the material interaction: energy loss
  bool energyLoss = true;
  /// Modify the behavior of the material interaction: scattering
  bool multipleScattering = true;
  /// Modify the behavior of the material interaction: record
  bool recordMaterialInteractions = true;
  /// looper protection
  double ptLoopers = 1 * Acts::UnitConstants::MeV;
  /// Max step size steering
  double maxStepSize = 1 * Acts::UnitConstants::m; 
  /// Max path limit
  double pathLimit = 30 * Acts::UnitConstants::m;
  /// Propagation direction
  Acts::Direction propDir = Acts::Direction::Forward();
  /// Transport Covariance 
  bool covarianceTransport = true;
};

template <typename propagator_t>
class ConcretePropagator {
 public:
 
  struct ForceCovTransport {
    template <typename propagator_state_t, typename stepper_t, typename navigator_t>
    void act(propagator_state_t& state,
            const stepper_t&,
            const navigator_t&,
            const Acts::Logger&) const {
      // std::cout << "Forcing Stepper to transport the covariance" << std::endl;
      state.stepping.covTransport = true;
      // std::cout << "Covariance value for this step: " << state.stepping.cov << std::endl;
    }

    // opzionale: se vuoi evitare ulteriori concept check su abort
    template <typename propagator_state_t, typename stepper_t, typename navigator_t>
    bool checkAbort(propagator_state_t&,
                    const stepper_t&,
                    const navigator_t&,
                    const Acts::Logger&) const {
      return false;
    }
  };

  explicit ConcretePropagator(propagator_t propagator)
      : m_propagator{std::move(propagator)} {}

    ~ConcretePropagator() noexcept {}

    ConcretePropagator* clone() const {
            return new ConcretePropagator(*this);
        }

  template <typename StartParameters>
  std::optional<std::pair<BoundParameters, double>> execute(
      const PropagationAlgorithm_Config& cfg,
      const Acts::Logger& logger,
      const StartParameters& startParameters,
      const Acts::Surface& targetSurf) const {
    ACTS_DEBUG("Test propagation/extrapolation starts");

    ACTS_VERBOSE("Starting propagation with these initial parameters: \n" << 
                 "Position: " << startParameters.position(Acts::GeometryContext{}).transpose() << "\n" <<
                 "Direction: " << startParameters.direction().transpose()); 
    ACTS_VERBOSE("Propagating to target surface of type " << targetSurf.type() << 
                 " and position " <<  targetSurf.center(Acts::GeometryContext{}).transpose());

    auto gctx = Acts::GeometryContext{};
    Acts::RecordedMaterial recordedMaterial; 

    // The step length logger for testing & end of world aborter
    using MaterialInteractor = Acts::MaterialInteractor;
    using SteppingLogger = Acts::detail::SteppingLogger;
    using EndOfWorld = Acts::EndOfWorldReached;

    // 

    // Actor list
    using TargetAborter = Acts::ForcedSurfaceReached;
    using ActorList = Acts::ActorList<SteppingLogger, MaterialInteractor, TargetAborter, EndOfWorld>;
    // using ActorList = Acts::ActorList<SteppingLogger, MaterialInteractor, EndOfWorld>;
    using PropagatorOptions = typename propagator_t::template Options<ActorList>;

    PropagatorOptions options(Acts::GeometryContext{}, Acts::MagneticFieldContext{});
    // Activate loop protection at some pt value
    options.loopProtection = startParameters.transverseMomentum() < cfg.ptLoopers;
    // Switch the material interaction on/off & eventually into logging mode
    auto& mInteractor = options.actorList.template get<MaterialInteractor>();
    mInteractor.multipleScattering = cfg.multipleScattering;
    mInteractor.energyLoss = cfg.energyLoss;
    mInteractor.recordInteractions = cfg.recordMaterialInteractions;
    mInteractor.noiseUpdateMode = Acts::NoiseUpdateMode::addNoise;

    options.direction = cfg.propDir;
    // options.navigation.nearLimit = -50 * Acts::UnitConstants::um;//-1e-4; 
    options.navigation.nearLimit = 0.;//-1e-4; 
    // options.navigation.surfaceTolerance = 1e-3; 
    options.stepping.maxStepSize = cfg.maxStepSize;

    auto& sLogger = options.actorList.template get<SteppingLogger>();
    sLogger.sterile = cfg.sterileLogger;
    options.pathLimit = cfg.pathLimit;

    const Acts::Surface* targetToUse = &targetSurf;
    std::shared_ptr<const Acts::Surface> inflated = nullptr;
    inflated = inflateTargetPlaneIfRectOrTrap(targetSurf, Acts::GeometryContext{}, 9 /*mm*/, 1.4);
    if (inflated) {
      targetToUse = inflated.get();
      ACTS_VERBOSE("Using inflated target bounds (+5 mm)");
      ACTS_VERBOSE("Using surface with bounds:");
      auto b_values = targetToUse->bounds().values();
      for (size_t i = 0; i < b_values.size(); ++i) {
          ACTS_VERBOSE(" " << b_values[i]);
      }
    }

    // TEST
    //options.direction = Acts::Direction::Forward();
 
    // auto x0 = startParameters.position(gctx);
    // auto u0 = startParameters.direction().normalized();

    // auto hasPositiveIntersection = [&](const Acts::Vector3& u) {
    //   auto res = targetToUse->intersect(gctx, x0, u);
    //   for (auto const& is : res) {
    //     if (is.pathLength() >= 0.0) return true;
    //   }
    //   return false;
    // };

    // bool okF = hasPositiveIntersection( u0);
    // bool okB = hasPositiveIntersection(-u0);

    // if(okF) ACTS_VERBOSE("hasPositiveIntersection");
    // if(okB) ACTS_VERBOSE("hasNegativeIntersection");

    // ACTS_VERBOSE(okF << " " << okB);

    // // Se l'utente chiede Backward ma Backward non ha intersezione positiva, flippo a Forward
    // if (cfg.propDir == Acts::Direction::Backward() && !okB && okF) {
    //   ACTS_VERBOSE("Requested BACKWARD but target is not reachable in backward direction. Switching to FORWARD.");
    //   options.direction = Acts::Direction::Forward();
    // } else if (cfg.propDir == Acts::Direction::Forward() && !okF && okB) {
    //   ACTS_VERBOSE("Requested FORWARD but target is not reachable in forward direction. Switching to BACKWARD.");
    //   options.direction = Acts::Direction::Backward();
    // } else {
    //   options.direction = cfg.propDir;
    // }
    // End of TEST

    // printer.target = targetToUse;

    ACTS_VERBOSE("START cov diag = " << startParameters.covariance()->diagonal().transpose());

    auto precheck = intersectInfinitePlanePrecheck(gctx, *targetToUse, startParameters, options.direction);

    if (!precheck) {
      ACTS_VERBOSE("[PRECHECK] No reachable intersection with infinite target plane.");
      ACTS_VERBOSE("[PRECHECK] Inverting propagation direction and check again");
      options.direction = (options.direction == Acts::Direction::Forward()) ? Acts::Direction::Backward() : Acts::Direction::Forward();
      precheck = intersectInfinitePlanePrecheck(gctx, *targetToUse, startParameters, options.direction);
      if (!precheck) {
        ACTS_VERBOSE("[PRECHECK] Still no reachable intersection with infinite target plane. Proceeding with original direction.");
        options.direction = cfg.propDir; // restore original direction
      } else {
        ACTS_VERBOSE("[PRECHECK] After inverting direction, found reachable intersection with infinite target plane at: " << precheck->global.transpose());
        ACTS_VERBOSE("[PRECHECK] Local coords on target: " << precheck->local.transpose());
        ACTS_VERBOSE("[PRECHECK] Path length: " << precheck->pathLength);
        ACTS_VERBOSE("[PRECHECK] Inside bounds? " << precheck->insideBounds);
      }

    } else {
      ACTS_VERBOSE("[PRECHECK] Infinite-plane intersection at: " << precheck->global.transpose());
      ACTS_VERBOSE("[PRECHECK] Local coords on target: " << precheck->local.transpose());
      ACTS_VERBOSE("[PRECHECK] Path length: " << precheck->pathLength);
      ACTS_VERBOSE("[PRECHECK] Inside bounds? " << precheck->insideBounds);
    }

    // DEBUG
    // if (options.direction == Acts::Direction::Backward()) {
    //   dumpTargetGeometryPre(gctx, *targetToUse, startParameters, "[SMOOTH PRE]");
    // }
    // End DEBUG

    // DEBUG
    // options.direction = chooseDirectionTowardTarget(gctx, *targetToUse, startParameters);
    // End DEBUG

    //auto propRes = m_propagator.propagate(startParameters, *targetToUse, options);

    auto& targetAborter = options.actorList.template get<TargetAborter>();
    targetAborter.surface = targetToUse;
    targetAborter.nearLimit = -200 * Acts::UnitConstants::um;
    targetAborter.boundaryTolerance = Acts::BoundaryTolerance::Infinite();
    

    //if(precheck->insideBounds) {
    // ACTS_VERBOSE("[PROP Type] Use standard propagation to target surface.");
    // cit: Due to the geometry of the perigee surface the overstepping tolerance is sometimes not met.
    // from: https://github.com/acts-project/acts/blob/52492222b450e23064ccf9a977b61f5781916675/Core/include/Acts/Propagator/Propagator.ipp#L413-L416
    // defined here: https://github.com/acts-project/acts/blob/52492222b450e23064ccf9a977b61f5781916675/Core/include/Acts/Propagator/StandardAborters.hpp#L159-L165
    
    
    auto propRes = m_propagator.template propagate<StartParameters,
                                                  PropagatorOptions>(startParameters, *targetToUse, options);
    // auto propRes = m_propagator.template propagate<StartParameters, PropagatorOptions>(startParameters, options);
    
    if (!propRes.ok() || !propRes->endParameters) {
      ACTS_VERBOSE("Propagator Failed! Returning nullopt");
      return std::nullopt;
    }

    if (!propRes->endParameters->covariance().has_value()) {
      ACTS_VERBOSE("END parameters have NO covariance");
    } else {
      ACTS_VERBOSE("END cov diag = " << propRes->endParameters->covariance()->diagonal().transpose());
    }
  
    BoundParameters finalBound = *propRes->endParameters;
    double length = propRes->pathLength;
    return std::make_optional(std::make_pair(finalBound, length));
                                           
    // } else {
    //   ACTS_VERBOSE("[PROP Type] Use custom fallback propagation to infinite plane.");

    //   std::cout << "Target surface is not reachable within bounds from the start point. Fall back: performing infinite-plane propagation." << std::endl;
    //   options.pathLimit = std::abs(precheck->pathLength);
    //   auto propRes = m_propagator.template propagate<StartParameters, PropagatorOptions>(startParameters, options);
    //   if (!propRes.ok() || !propRes->endParameters) {
    //     ACTS_VERBOSE("Propagator Failed! Returning nullopt");
    //     std::cout << "Propagator failed during fallback infinite-plane propagation. Returning nullopt." << std::endl;
    //     return std::nullopt;
    //   }

    //   // Check if we have a valid covariance at the end of the fallback propagation. If not, return nullopt to signal failure.
    //   const auto& endPars = *propRes->endParameters;
    //   if (!endPars.covariance().has_value()) {
    //     std::cout << "Fallback propagation reached the infinite plane but did not produce a valid covariance. Returning nullopt." << std::endl;
    //     return std::nullopt;
    //   }

    //   // Now check if the final position is actually on the target surface. If not, return nullopt.
    //   /// NOTE: here we're checking if the final position is ON the infinite plane, not if we're inside the boounds (tolerance in this case is 1e-4 (0.1um), i.e. the default value)
    //   auto lposRes = targetSurf.globalToLocal(gctx, endPars.position(gctx), endPars.direction());
    //   if (!lposRes.ok()) {
    //     return std::nullopt;
    //   }
    //   // Acts::Vector2 lpos = *lposRes;

    //   // Get the local covariance and check if it's positive definite. If not, return nullopt.
    //   const auto& cov = *endPars.covariance();
    //   Acts::SquareMatrix2 localCov;
    //   localCov << cov(Acts::eBoundLoc0, Acts::eBoundLoc0),
    //               cov(Acts::eBoundLoc0, Acts::eBoundLoc1),
    //               cov(Acts::eBoundLoc1, Acts::eBoundLoc0),
    //               cov(Acts::eBoundLoc1, Acts::eBoundLoc1);

    //   if (!(localCov.determinant() > 0.)) {
    //     std::cout << "Fallback propagation reached the infinite plane but produced a non-positive-definite local covariance. Returning nullopt." << std::endl;
    //     return std::nullopt;
    //   }

    //   if (!propRes->endParameters->covariance().has_value()) {
    //     ACTS_VERBOSE("END parameters have NO covariance");
    //   } else {
    //     ACTS_VERBOSE("END cov diag = " << propRes->endParameters->covariance()->diagonal().transpose());
    //   }
    
    //   BoundParameters finalBound = *propRes->endParameters;
    //   double length = propRes->pathLength;
    //   return std::make_optional(std::make_pair(finalBound, length));
    // }

  }

 private:
  propagator_t m_propagator;
};

#endif
