#include <sstream>

#include "DataFormats/Common/interface/OrphanHandle.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackerRecHit2D/interface/TrackingRecHitLessFromGlobalPosition.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/thread_safety_macros.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoTracker/TrackProducer/interface/TrackProducerActsAlgorithm.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TRecHit1DMomConstraint.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TRecHit2DPosConstraint.h"
#include "RecoTracker/TransientTrackingRecHit/interface/TkTransientTrackingRecHitBuilder.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderWithPropagator.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/TrackFitters/interface/RecHitSorter.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/GsfTools/interface/GetComponents.h"

// ======================================================================
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "TrackPropagation/Acts/interface/ActsPropagator.h"
#include "TrackingTools/PatternTools/interface/TrajectoryExtrapolatorToLine.h"
// ======================================================================

// #define VI_DEBUG
// #define STAT_TSB

#ifdef VI_DEBUG
#define DPRINT(x) std::cout << x << ": "
#else
#define DPRINT(x) LogTrace(x)
#endif

namespace {
#ifdef STAT_TSB
  struct StatCount {
    long long totTrack = 0;
    long long totLoop = 0;
    long long totGsfTrack = 0;
    long long totFound = 0;
    long long totLost = 0;
    long long totAlgo[15];
    void track(int l) {
      if (l > 0)
        ++totLoop;
      else
        ++totTrack;
    }
    void hits(int f, int l) {
      totFound += f;
      totLost += l;
    }
    void gsf() { ++totGsfTrack; }
    void algo(int a) {
      if (a >= 0 && a < 15)
        ++totAlgo[a];
    }

    void print() const {
      std::cout << "TrackProducer stat\nTrack/Loop/Gsf/FoundHits/LostHits//algos " << totTrack << '/' << totLoop << '/'
                << totGsfTrack << '/' << totFound << '/' << totLost << '/';
      for (auto a : totAlgo)
        std::cout << '/' << a;
      std::cout << std::endl;
    }
    StatCount() {}
    ~StatCount() { print(); }
  };
  StatCount statCount;

#else
  struct StatCount {
    void track(int) {}
    void hits(int, int) {}
    void gsf() {}
    void algo(int) {}
  };
  CMS_THREAD_SAFE StatCount statCount;
#endif

}  // namespace

template <>
bool TrackProducerActsAlgorithm<reco::Track>::buildTrack(const TrajectoryFitter* theFitter,
                                                     const Propagator* thePropagator,
                                                     AlgoProductCollection& algoResults,
                                                     TransientTrackingRecHit::RecHitContainer& hits,
                                                     TrajectoryStateOnSurface& theTSOS,
                                                     const TrajectorySeed& seed,
                                                     float ndof,
                                                     const reco::BeamSpot& bs,
                                                     SeedRef seedRef,
                                                     int qualityMask,
                                                     signed char nLoops) {
  //variable declarations

  PropagationDirection seedDir = seed.direction();

  //perform the fit: the result's size is 1 if it succeded, 0 if fails
  Trajectory&& trajTmp =
      theFitter->fitOne(seed, hits, theTSOS, (nLoops > 0) ? TrajectoryFitter::looper : TrajectoryFitter::standard);
  if UNLIKELY (!trajTmp.isValid()) {
    DPRINT("TrackFitters") << "fit failed " << algo_ << ": " << hits.size() << '|' << int(nLoops) << ' ' << std::endl;
    return false;
  }

  auto theTraj = new Trajectory(std::move(trajTmp));
  theTraj->setSeedRef(seedRef);

  statCount.hits(theTraj->foundHits(), theTraj->lostHits());
  statCount.algo(int(algo_));

  // TrajectoryStateOnSurface innertsos;
  // if (theTraj->direction() == alongMomentum) {
  //  innertsos = theTraj->firstMeasurement().updatedState();
  // } else {
  //  innertsos = theTraj->lastMeasurement().updatedState();
  // }

  ndof = 0;
  for (auto const& tm : theTraj->measurements()) {
    auto const& h = tm.recHitR();
    if (h.isValid())
      ndof = ndof + float(h.dimension()) * h.weight();  // two virtual calls!
  }

  ndof -= 5.f;
  if UNLIKELY (theTSOS.magneticField()->nominalValue() == 0)
    ++ndof;  // same as -4

#if defined(VI_DEBUG) || defined(EDM_ML_DEBUG)
  int chit[7] = {};
  int kk = 0;
  for (auto const& tm : theTraj->measurements()) {
    ++kk;
    auto const& hit = tm.recHitR();
    if (!hit.isValid())
      ++chit[0];
    if (hit.det() == nullptr)
      ++chit[1];
    if (trackerHitRTTI::isUndef(hit))
      continue;
    if (0)
      std::cout << "h " << kk << ": " << hit.localPosition() << ' ' << hit.localPositionError() << ' ' << tm.estimate()
                << std::endl;
    if (hit.dimension() != 2) {
      ++chit[2];
    } else if (trackerHitRTTI::isFromDet(hit)) {
      auto const& thit = static_cast<BaseTrackerRecHit const&>(hit);
      auto const& clus = thit.firstClusterRef();
      if (clus.isPixel())
        ++chit[3];
      else if (thit.isMatched()) {
        ++chit[4];
      } else if (thit.isProjected()) {
        ++chit[5];
      } else {
        ++chit[6];
      }
    }
  }

  std::ostringstream ss;
  ss << algo_ << ": " << hits.size() << '|' << theTraj->measurements().size() << '|' << int(nLoops) << ' ';
  for (auto c : chit)
    ss << c << '/';
  ss << std::endl;
  DPRINT("TrackProducer") << ss.str();

#endif

  //if geometricInnerState_ is false the state for projection to beam line is the state attached to the first hit: to be used for loopers
  //if geometricInnerState_ is true the state for projection to beam line is the one from the (geometrically) closest measurement to the beam line: to be sued for non-collision tracks
  //the two shouuld give the same result for collision tracks that are NOT loopers
  TrajectoryStateOnSurface stateForProjectionToBeamLineOnSurface;
  if (geometricInnerState_) {
    stateForProjectionToBeamLineOnSurface =
        theTraj->closestMeasurement(GlobalPoint(bs.x0(), bs.y0(), bs.z0())).updatedState();
  } else {
    if (theTraj->direction() == alongMomentum) {
      stateForProjectionToBeamLineOnSurface = theTraj->firstMeasurement().updatedState();
    } else {
      stateForProjectionToBeamLineOnSurface = theTraj->lastMeasurement().updatedState();
    }
  }

  if UNLIKELY (!stateForProjectionToBeamLineOnSurface.isValid()) {
    edm::LogError("CannotPropagateToBeamLine") << "the state on the closest measurement isnot valid. skipping track.";
    delete theTraj;
    return false;
  }
  const FreeTrajectoryState& stateForProjectionToBeamLine = *stateForProjectionToBeamLineOnSurface.freeState();

  LogDebug("TrackProducer") << "stateForProjectionToBeamLine=" << stateForProjectionToBeamLine;

  //  TSCBLBuilderNoMaterial tscblBuilder;
  //  TrajectoryStateClosestToBeamLine tscbl = tscblBuilder(stateForProjectionToBeamLine,bs);

  // =============================================
  // FreeTrajectoryState NewStateForPCA = stateForProjectionToBeamLine;
  FreeTrajectoryState fts_acts = stateForProjectionToBeamLine;
  TrajectoryStateClosestToBeamLine tscbl;

  // ++++++++ NEW IDEA (Propagate to PCA using ACTS)++++++++
  const auto* actsProp = dynamic_cast<const ActsPropagator*>(thePropagator);
  bool propOk = true;
  if (actsProp != nullptr) {
    auto stateInsideBP = actsProp->propagateWithPathToPerigeeInsideBP(stateForProjectionToBeamLineOnSurface);
    if (stateInsideBP.position().basicVector().mag2() > 0.) {
      fts_acts = stateInsideBP;
    } else {
      std::cout << "ACTS propagation to perigee inside BP failed; " << std::endl;
      propOk = false;
    }
  }

  /// NOTE: code from https://github.com/cms-sw/cmssw/blob/1f8b4e6f3e9f5bb865a5b8c45a7f6b4fb0b86d94/TrackingTools/PatternTools/src/TSCBLBuilderWithPropagator.cc#L29
  if (propOk) {
    GlobalPoint tp(fts_acts.position().x(), fts_acts.position().y(), fts_acts.position().z());
    GlobalPoint bspos(bs.position().x(), bs.position().y(), bs.position().z());
    GlobalVector bsvec(bs.dxdz(), bs.dydz(), 1.);
    Line bsline(bspos, bsvec);
    GlobalVector hyp(tp.x() - bspos.x(), tp.y() - bspos.y(), tp.z() - bspos.z()); 
    double l = bsline.direction().dot(hyp);
    GlobalPoint closepoint = bspos + l * bsvec;
    tscbl = TrajectoryStateClosestToBeamLine(fts_acts, closepoint, bs);
  }
  else {
    // Fallback: if ACTS propagation failed, use the old TSCBLBuilder with the original state
    TSCBLBuilderNoMaterial tscblBuilder;
    tscbl = tscblBuilder(fts_acts, bs);
  }
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++


  // std::cout << "Starting propagation to the PCA" << std::endl;
  // /// NOTE: this is the only difference wrt the cmssw code:
  // /// We give to the tscblBuilder not the state provided by cmssw, but we first propagate it INTO the bp

  // const Surface& oldSurf = stateForProjectionToBeamLineOnSurface.surface();
  // const GlobalVector oldX = oldSurf.toGlobal(LocalVector(1., 0., 0.)).unit();
  // const GlobalVector oldY = oldSurf.toGlobal(LocalVector(0., 1., 0.)).unit();
  // const GlobalVector oldZ = oldSurf.toGlobal(LocalVector(0., 0., 1.)).unit();

  // const GlobalPoint pos0 = stateForProjectionToBeamLine.position();

  // // Take the transverse radial direction of the starting point,
  // // move along that direction to R
  // // and keep the same z coordinate as the starting point.
  // GlobalVector radial(pos0.x(), pos0.y(), 0.);
  // if (radial.mag2() < 1e-12) {
  //   const GlobalVector mom0 = stateForProjectionToBeamLine.momentum();
  //   radial = GlobalVector(mom0.x(), mom0.y(), 0.);
  // }
  // radial = radial.unit();

  // const double targetR = 2.3; // cm
  // const GlobalPoint planeCenter(targetR * radial.x(), targetR * radial.y(), pos0.z());

  // Surface::RotationType planeRot(oldX, oldY, oldZ);

  // auto planeInsideBP = Plane::build(Surface::PositionType(planeCenter), planeRot, new RectangularPlaneBounds(5., 300., 0.));

  // DEBUG
  // {
  //   const GlobalVector newY_ = planeInsideBP->toGlobal(LocalVector(0., 1., 0.)).unit();
  //   std::cout << "[BP_PLANE_SIZE_CHECK]"
  //             << " newYDotBeam " << newY_.dot(GlobalVector(0., 0., 1.))
  //             << "\n";

  //   const Surface& oldSurf = stateForProjectionToBeamLineOnSurface.surface();

  //   auto printVec = [](const char* name, const GlobalVector& v) {
  //     std::cout
  //       << name << "=("
  //       << v.x() << "," << v.y() << "," << v.z() << ")";
  //   };

  //   auto dot = [](const GlobalVector& a, const GlobalVector& b) {
  //     return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
  //   };

  //   const GlobalPoint oldC = oldSurf.position();
  //   const GlobalPoint newC = planeInsideBP->position();

  //   const GlobalVector oldX = oldSurf.toGlobal(LocalVector(1., 0., 0.)).unit();
  //   const GlobalVector oldY = oldSurf.toGlobal(LocalVector(0., 1., 0.)).unit();
  //   const GlobalVector oldZ = oldSurf.toGlobal(LocalVector(0., 0., 1.)).unit();

  //   const GlobalVector newX = planeInsideBP->toGlobal(LocalVector(1., 0., 0.)).unit();
  //   const GlobalVector newY = planeInsideBP->toGlobal(LocalVector(0., 1., 0.)).unit();
  //   const GlobalVector newZ = planeInsideBP->toGlobal(LocalVector(0., 0., 1.)).unit();

  //   const double oldR = std::sqrt(oldC.x() * oldC.x() + oldC.y() * oldC.y());
  //   const double newR = std::sqrt(newC.x() * newC.x() + newC.y() * newC.y());

  //   const GlobalVector dC = newC - oldC;

  //   std::cout
  //     << "[SURF_TRANSLATION_CHECK]"
  //     << " oldCenter=(" << oldC.x() << "," << oldC.y() << "," << oldC.z() << ")"
  //     << " newCenter=(" << newC.x() << "," << newC.y() << "," << newC.z() << ")"
  //     << " oldR " << oldR
  //     << " newR " << newR
  //     << " dCenter=(" << dC.x() << "," << dC.y() << "," << dC.z() << ")"

  //     << " oldX_newX " << dot(oldX, newX)
  //     << " oldY_newY " << dot(oldY, newY)
  //     << " oldZ_newZ " << dot(oldZ, newZ)

  //     << " oldX_newY " << dot(oldX, newY)
  //     << " oldX_newZ " << dot(oldX, newZ)
  //     << " oldY_newX " << dot(oldY, newX)
  //     << " oldY_newZ " << dot(oldY, newZ)
  //     << " oldZ_newX " << dot(oldZ, newX)
  //     << " oldZ_newY " << dot(oldZ, newY)
  //     << "\n";

  //   std::cout << "[SURF_AXES] ";
  //     printVec(" oldX", oldX);
  //     std::cout << " ";
  //     printVec(" oldY", oldY);
  //     std::cout << " ";
  //     printVec(" oldZ", oldZ);
  //     std::cout << " ";
  //     printVec(" newX", newX);
  //     std::cout << " ";
  //     printVec(" newY", newY);
  //     std::cout << " ";
  //     printVec(" newZ", newZ);
  //     std::cout << "\n";

  //   if (std::abs(dot(oldX, newX) - 1.) > 1e-6 ||
  //       std::abs(dot(oldY, newY) - 1.) > 1e-6 ||
  //       std::abs(dot(oldZ, newZ) - 1.) > 1e-6 ||
  //       std::abs(dot(oldX, newY)) > 1e-6 ||
  //       std::abs(dot(oldX, newZ)) > 1e-6 ||
  //       std::abs(dot(oldY, newX)) > 1e-6 ||
  //       std::abs(dot(oldY, newZ)) > 1e-6 ||
  //       std::abs(dot(oldZ, newX)) > 1e-6 ||
  //       std::abs(dot(oldZ, newY)) > 1e-6) {
  //     std::cout << "[SURF_TRANSLATION_MISMATCH] local axes are not identical\n";
  //   } else {
  //     std::cout << "[SURF_TRANSLATION_OK] local axes are identical; only center changed\n";
  //   }
  // }
  // END DEBUG

  // -> BEGIN OLD 
  // // Read old TSOS surface frame
  // const Surface& oldSurf = stateForProjectionToBeamLineOnSurface.surface();
  // const GlobalPoint oldSurfPos = oldSurf.position();
  // const GlobalVector oldNormal = oldSurf.toGlobal(LocalVector(0., 0., 1.)).unit();
  // const GlobalVector oldYAxis = oldSurf.toGlobal(LocalVector(0., 1., 0.)).unit();

  // const GlobalVector beamAxis(0., 0., 1.);

  // // Normal orientation wrt radial direction
  // double normalSign = +1.0;
  // GlobalVector oldRadial(oldSurfPos.x(), oldSurfPos.y(), 0.);
  // if (oldRadial.mag2() > 1e-12) {
  //   oldRadial = oldRadial.unit();

  //   const double nDotRadial = oldNormal.dot(oldRadial);
  //   normalSign = (nDotRadial >= 0.) ? +1.0 : -1.0;
  // }

  // // Local Y orientation wrt beam pipe
  // const double yDotBeam = oldYAxis.dot(beamAxis);
  // const bool oldYAlongBeam = std::abs(yDotBeam) > 0.9;
  // const double yBeamSign = (yDotBeam >= 0.) ? +1.0 : -1.0;

  // // Build detector-like plane inside beam pipe
  // const double targetR = 1.5; // cm

  // const GlobalPoint pos0 = stateForProjectionToBeamLine.position();
  // const GlobalVector mom0 = stateForProjectionToBeamLine.momentum();

  // const double phi0 = std::atan2(mom0.y(), mom0.x());

  // const GlobalPoint planeCenter(targetR * std::cos(phi0), targetR * std::sin(phi0), pos0.z());

  // GlobalVector radialNew(planeCenter.x(), planeCenter.y(), 0.);
  // if (radialNew.mag2() < 1e-12) {
  //   radialNew = GlobalVector(std::cos(phi0), std::sin(phi0), 0.);
  // }
  // radialNew = radialNew.unit();

  // // local z: radial normal, with same inward/outward sign as old surface
  // GlobalVector zAxis = normalSign * radialNew;

  // // local y: if old local Y is along beam pipe, preserve its beam direction
  // GlobalVector yAxis;
  // if (oldYAlongBeam) {
  //   yAxis = GlobalVector(0., 0., yBeamSign);
  // } else {
  //   // Fallback: detector-like barrel convention
  //   yAxis = GlobalVector(0., 0., 1.);
  // }

  // // Make y perpendicular to z, robustly
  // yAxis = yAxis - yAxis.dot(zAxis) * zAxis;
  // if (yAxis.mag2() < 1e-12) {
  //   yAxis = GlobalVector(0., 0., 1.);
  // }
  // yAxis = yAxis.unit();

  // // Build right-handed frame: x cross y = z
  // GlobalVector xAxis = yAxis.cross(zAxis).unit();
  // yAxis = zAxis.cross(xAxis).unit();

  // Surface::RotationType planeRot(xAxis, yAxis, zAxis);

  // auto planeInsideBP = Plane::build(Surface::PositionType(planeCenter), planeRot, new RectangularPlaneBounds(5., 300., 0.));

  // Final debug for new surface DEBUG
  // {
  //   const GlobalVector newX =
  //       planeInsideBP->toGlobal(LocalVector(1., 0., 0.)).unit();
  //   const GlobalVector newY =
  //       planeInsideBP->toGlobal(LocalVector(0., 1., 0.)).unit();
  //   const GlobalVector newZ =
  //       planeInsideBP->toGlobal(LocalVector(0., 0., 1.)).unit();

  //   const double newNDotRadial = newZ.dot(radialNew);
  //   const double newYDotBeam = newY.dot(beamAxis);

  //   std::cout
  //     << "[NEW_BP_SURF_FRAME]"
  //     << " center=(" << planeCenter.x() << "," << planeCenter.y() << "," << planeCenter.z() << ")"
  //     << " newX=(" << newX.x() << "," << newX.y() << "," << newX.z() << ")"
  //     << " newY=(" << newY.x() << "," << newY.y() << "," << newY.z() << ")"
  //     << " newZ=(" << newZ.x() << "," << newZ.y() << "," << newZ.z() << ")"
  //     << " newNDotRadial=" << newNDotRadial
  //     << " newYDotBeam=" << newYDotBeam
  //     << "\n";

  //   std::cout
  //   << "[OLD_TO_NEW_FRAME_MATCH]"
  //   << " oldYDotBeam " << yDotBeam
  //   << " oldYAlongBeam " << oldYAlongBeam
  //   << " normalSign " << normalSign
  //   << "\n";
  // }
  //END DEBUG
  // -> END OLD

  // // Propagate to the new surface
  // auto stateInsideBP_WP = thePropagator->propagateWithPath(stateForProjectionToBeamLineOnSurface, *planeInsideBP);

  // if (stateInsideBP_WP.first.isValid()) {
  //   NewStateForPCA = *stateInsideBP_WP.first.freeState();
  // } else {
  //   std::cout << "State not propagated properly to detector-like BP plane\n";
  // }

  // TrajectoryStateClosestToBeamLine tscbl;
  // if (usePropagatorForPCA_) {
  //   //std::cout << "PROPAGATOR FOR PCA" << std::endl;
  //   // std::cout << "Performing the actual propagation to the PCA" << std::endl;
  //   TSCBLBuilderWithPropagator tscblBuilder(*thePropagator);
  //   tscbl = tscblBuilder(NewStateForPCA, bs);
  // } else {
  //   TSCBLBuilderNoMaterial tscblBuilder;
  //   tscbl = tscblBuilder(NewStateForPCA, bs);
  // }
  // std::cout << "Propagation to the PCA concluded" << std::endl;
  // =============================================

  if UNLIKELY (!tscbl.isValid()) {
    delete theTraj;
    return false;
  }

  GlobalPoint v = tscbl.trackStateAtPCA().position();
  math::XYZPoint pos(v.x(), v.y(), v.z());
  GlobalVector p = tscbl.trackStateAtPCA().momentum();
  math::XYZVector mom(p.x(), p.y(), p.z());

  LogDebug("TrackProducer") << "pos=" << v << " mom=" << p << " pt=" << p.perp() << " mag=" << p.mag();

  auto theTrack = new reco::Track(theTraj->chiSquared(),
                                  int(ndof),  //FIXME fix weight() in TrackingRecHit
                                  pos,
                                  mom,
                                  tscbl.trackStateAtPCA().charge(),
                                  tscbl.trackStateAtPCA().curvilinearError(),
                                  algo_);

  if (originalAlgo_ != reco::TrackBase::undefAlgorithm)
    theTrack->setOriginalAlgorithm(originalAlgo_);
  if (algoMask_.any())
    theTrack->setAlgoMask(algoMask_);
  theTrack->setQualityMask(qualityMask);
  theTrack->setNLoops(nLoops);
  theTrack->setStopReason(stopReason_);

  LogDebug("TrackProducer") << "theTrack->pt()=" << theTrack->pt();

  LogDebug("TrackProducer") << "track done\n";

  AlgoProduct aProduct{theTraj, theTrack, seedDir, 0};
  algoResults.push_back(aProduct);

  statCount.track(nLoops);

  return true;
}

template <>
bool TrackProducerActsAlgorithm<reco::GsfTrack>::buildTrack(const TrajectoryFitter* theFitter,
                                                        const Propagator* thePropagator,
                                                        AlgoProductCollection& algoResults,
                                                        TransientTrackingRecHit::RecHitContainer& hits,
                                                        TrajectoryStateOnSurface& theTSOS,
                                                        const TrajectorySeed& seed,
                                                        float ndof,
                                                        const reco::BeamSpot& bs,
                                                        SeedRef seedRef,
                                                        int qualityMask,
                                                        signed char nLoops) {
  PropagationDirection seedDir = seed.direction();

  std::cout << "Second method called" << std::endl;

  Trajectory&& trajTmp =
      theFitter->fitOne(seed, hits, theTSOS, (nLoops > 0) ? TrajectoryFitter::looper : TrajectoryFitter::standard);
  if UNLIKELY (!trajTmp.isValid())
    return false;

  auto theTraj = new Trajectory(std::move(trajTmp));
  theTraj->setSeedRef(seedRef);

#ifdef EDM_ML_DEBUG
  TrajectoryStateOnSurface innertsos;
  TrajectoryStateOnSurface outertsos;

  if (theTraj->direction() == alongMomentum) {
    innertsos = theTraj->firstMeasurement().updatedState();
    outertsos = theTraj->lastMeasurement().updatedState();
  } else {
    innertsos = theTraj->lastMeasurement().updatedState();
    outertsos = theTraj->firstMeasurement().updatedState();
  }
  std::ostringstream ss;
  auto dc = [&](TrajectoryStateOnSurface const& tsos) {
    GetComponents comps(tsos);
    auto const& components = comps();
    auto sinTheta = std::sin(tsos.globalMomentum().theta());
    for (auto const& ic : components)
      ss << ic.weight() << "/";
    ss << "\n";
    for (auto const& ic : components)
      ss << ic.localParameters().vector()[0] / sinTheta << "/";
    ss << "\n";
    for (auto const& ic : components)
      ss << std::sqrt(ic.localError().matrix()(0, 0)) / sinTheta << "/";
  };
  ss << "\ninner comps\n";
  dc(innertsos);
  GetComponents icomps(innertsos);
  auto const& tsosComponentsInner = icomps();

  ss << "\nouter comps\n";
  dc(outertsos);
  GetComponents ocomps(outertsos);
  auto const& tsosComponentsOuter = ocomps();

  LogDebug("TrackProducer") << "Nr. of first / last states = " << tsosComponentsInner.size() << " "
                            << tsosComponentsOuter.size() << ss.str();
#endif

  ndof = 0;
  for (auto const& tm : theTraj->measurements()) {
    auto const& h = tm.recHitR();
    if (h.isValid())
      ndof = ndof + h.dimension() * h.weight();
  }

  ndof = ndof - 5;
  if UNLIKELY (theTSOS.magneticField()->nominalValue() == 0)
    ++ndof;  // same as -4

  //if geometricInnerState_ is false the state for projection to beam line is the state attached to the first hit: to be used for loopers
  //if geometricInnerState_ is true the state for projection to beam line is the one from the (geometrically) closest measurement to the beam line: to be sued for non-collision tracks
  //the two shouuld give the same result for collision tracks that are NOT loopers
  TrajectoryStateOnSurface stateForProjectionToBeamLineOnSurface;
  if (geometricInnerState_) {
    stateForProjectionToBeamLineOnSurface =
        theTraj->closestMeasurement(GlobalPoint(bs.x0(), bs.y0(), bs.z0())).updatedState();
  } else {
    if (theTraj->direction() == alongMomentum) {
      stateForProjectionToBeamLineOnSurface = theTraj->firstMeasurement().updatedState();
    } else {
      stateForProjectionToBeamLineOnSurface = theTraj->lastMeasurement().updatedState();
    }
  }

  if UNLIKELY (!stateForProjectionToBeamLineOnSurface.isValid()) {
    edm::LogError("CannotPropagateToBeamLine") << "the state on the closest measurement isnot valid. skipping track.";
    delete theTraj;
    return false;
  }

  const FreeTrajectoryState& stateForProjectionToBeamLine = *stateForProjectionToBeamLineOnSurface.freeState();

  LogDebug("GsfTrackProducer") << "stateForProjectionToBeamLine=" << stateForProjectionToBeamLine;

  //  TSCBLBuilderNoMaterial tscblBuilder;
  //  TrajectoryStateClosestToBeamLine tscbl = tscblBuilder(stateForProjectionToBeamLine,bs);

  TrajectoryStateClosestToBeamLine tscbl;
  if (usePropagatorForPCA_) {
    TSCBLBuilderWithPropagator tscblBuilder(*thePropagator);
    tscbl = tscblBuilder(stateForProjectionToBeamLine, bs);
  } else {
    TSCBLBuilderNoMaterial tscblBuilder;
    tscbl = tscblBuilder(stateForProjectionToBeamLine, bs);
  }

  if UNLIKELY (tscbl.isValid() == false) {
    delete theTraj;
    return false;
  }

  GlobalPoint v = tscbl.trackStateAtPCA().position();
  math::XYZPoint pos(v.x(), v.y(), v.z());
  GlobalVector p = tscbl.trackStateAtPCA().momentum();
  math::XYZVector mom(p.x(), p.y(), p.z());

  LogDebug("GsfTrackProducer") << "pos=" << v << " mom=" << p << " pt=" << p.perp() << " mag=" << p.mag();

  auto theTrack =
      new reco::GsfTrack(theTraj->chiSquared(),
                         int(ndof),  //FIXME fix weight() in TrackingRecHit
                         //			       theTraj->foundHits(),//FIXME to be fixed in Trajectory.h
                         //			       0, //FIXME no corresponding method in trajectory.h
                         //			       theTraj->lostHits(),//FIXME to be fixed in Trajectory.h
                         pos,
                         mom,
                         tscbl.trackStateAtPCA().charge(),
                         tscbl.trackStateAtPCA().curvilinearError());
  theTrack->setAlgorithm(algo_);
  if (originalAlgo_ != reco::TrackBase::undefAlgorithm)
    theTrack->setOriginalAlgorithm(originalAlgo_);
  if (algoMask_.any())
    theTrack->setAlgoMask(algoMask_);

  theTrack->setStopReason(stopReason_);

  LogDebug("GsfTrackProducer") << "track done\n";

  AlgoProduct aProduct{theTraj, theTrack, seedDir, 0};

  LogDebug("GsfTrackProducer") << "track done1\n";
  algoResults.push_back(aProduct);
  LogDebug("GsfTrackProducer") << "track done2\n";

  statCount.gsf();
  return true;
}
