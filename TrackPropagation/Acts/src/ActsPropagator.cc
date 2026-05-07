#include <sstream>

#include "TrackPropagation/Acts/interface/ActsPropagator.h"

// CMSSW
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TrajectoryState/interface/SurfaceSideDefinition.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TrackingTools/AnalyticalJacobians/interface/AnalyticalCurvilinearJacobian.h"

#include "TrackPropagation/Acts/interface/ComputeLocalBoundJacobian.h"

// CLHEP
#include <CLHEP/Units/SystemOfUnits.h>

namespace {
  void printMatrix(const AlgebraicSymMatrix55& m) {
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      std::cout << m(i,j) << " ";
    }
    std::cout << std::endl;
  }
}
}

Acts::SquareMatrix3 localToGlobalRotation(const Surface& surf) {
  const GlobalPoint  g0  = surf.toGlobal(LocalPoint(0., 0.));
  const GlobalPoint  g1x = surf.toGlobal(LocalPoint(1., 0.));
  const GlobalPoint  g1y = surf.toGlobal(LocalPoint(0., 1.));

  GlobalVector gx = (g1x - g0).unit();
  GlobalVector gy = (g1y - g0).unit();

  // z da prodotto vettoriale (può avere verso ambiguo)
  GlobalVector gz = gx.cross(gy).unit();
  gy = gz.cross(gx).unit(); // ri-ortogonalizza

  // *** FIX: forza gz ad avere lo stesso verso del "vero" z locale CMSSW ***
  // Se gz è il +z locale, allora surf.toLocal(gz) dovrebbe avere componente z positiva.
  // LocalVector zLoc = surf.toLocal(GlobalVector(gz.x(), gz.y(), gz.z()));
  // if (zLoc.z() < 0) {
  //   gz = GlobalVector(-gz.x(), -gz.y(), -gz.z());
  //   gy = GlobalVector(-gy.x(), -gy.y(), -gy.z());  // per mantenere terna destrorsa con gx
  // }
  // *** END FIX ***
  LocalVector zLoc = surf.toLocal(GlobalVector(gz.x(), gz.y(), gz.z()));
  if (zLoc.z() < 0) {
    gz = -gz;
    gy = -gy;
  }

  Acts::SquareMatrix3 R;
  R(0,0)=gx.x(); R(1,0)=gx.y(); R(2,0)=gx.z();
  R(0,1)=gy.x(); R(1,1)=gy.y(); R(2,1)=gy.z();
  R(0,2)=gz.x(); R(1,2)=gz.y(); R(2,2)=gz.z();

  return R;
}

/** Constructor.
 */
ActsPropagator::ActsPropagator(PropagationDirection dir,
                                const float mass,
                                const MagneticField* mf,
                                const float maxDPhi,
                                bool useRungeKutta,
                                float ptMin,
                                bool useOldAnalPropLogic,
                                //
                                std::shared_ptr<TrackingGeometryWithDetEls> TrkandDetEls,
                                const TrackerGeometry* trkGeo_cmssw,
                                const Acts::Logging::Level& lvl)
    : field(mf),
      trkGeo_and_DetEls_(TrkandDetEls),
      trkGeo_cmssw_(trkGeo_cmssw),
      actsPropDir_((dir == oppositeToMomentum) ? Acts::Direction::Backward() : Acts::Direction::Forward()),
      m_Level(lvl), m_logger(Acts::getDefaultLogger("ActsPropagator", lvl)) {

  // ===== Define the ACTS propagator to be used in this class =====
  // I) EigenStepper
  std::shared_ptr<const CMSMagneticFieldProvider>  magFieldPtr;
  if(field != nullptr){
    magFieldPtr = std::make_shared<const CMSMagneticFieldProvider>(*field, m_Level);
  } else {
    throw cms::Exception("ACTSMagneticFieldProvider") << " B field ptr is null";
  }
  Acts::EigenStepper<> es(magFieldPtr); 
  // II) Navigator
  Acts::Navigator::Config navi_cfg;
  navi_cfg.trackingGeometry = trkGeo_and_DetEls_->trackingGeometry;
  navi_cfg.resolveSensitive = true;
  navi_cfg.resolveMaterial = true;
  navi_cfg.resolvePassive = false; // Does this make the navigation more "complex"?
  std::shared_ptr<const Acts::Logger> navi_logger = Acts::getDefaultLogger("Navigator", m_Level);
  Acts::Navigator navi(navi_cfg, std::move(navi_logger));
  std::shared_ptr<const Acts::Logger> prop_logger = Acts::getDefaultLogger("Propagator", m_Level);
  // III) Propagator
  Acts::Propagator prop(es, navi, std::move(prop_logger));

  // IV) Define the concrete propagator
  ConcProp = DeepCopyPointerByClone<MyConcretePropagator>(new MyConcretePropagator(prop));
}

/** Destructor.
 */
ActsPropagator::~ActsPropagator() {}

//
////////////////////////////////////////////////////////////////////////////
//

// Define the helpfull methods to be used during the propagation:

Acts::BoundMatrix ActsPropagator::convertCovCMSSWtoACTS(const TrajectoryStateOnSurface& tsos, AlgebraicSymMatrix55 cov_cmssw) const {

  // Define the Jacobian:
  ComputeLocalBoundJacobian computeJ;
  Eigen::Matrix<double,5,5> C5e;
  for(int i=0;i<5;++i)
    for(int j=0;j<5;++j)
      C5e(i,j) = cov_cmssw(i,j);

  Acts::ActsMatrix<6,5> J = computeJ.FromCMSSWtoACTS(tsos);

  Acts::BoundMatrix C6out = J * C5e * J.transpose(); // 6x6

  C6out(Acts::eBoundTime, Acts::eBoundTime) = 1e-6;

  return C6out;
}

AlgebraicSymMatrix55 ActsPropagator::convertCovACTStoCMSSW(const Surface& surf, double phi, double theta, const Acts::BoundMatrix& cov_acts) const {
  ComputeLocalBoundJacobian computeJ;
  Acts::ActsMatrix<5,6> K = computeJ.FromACTStoCMSSW(surf, phi, theta);

  Acts::ActsMatrix<5,5> Cc =  K * cov_acts * K.transpose();
  AlgebraicSymMatrix55 M;
  for (int i=0;i<5;++i) {
    for (int j=0;j<=i;++j) {
      M(i,j) = Cc(i,j);
    }
  }

  return M;
}

LocalTrajectoryParameters ActsPropagator::GetLocalTrajectoryParameters(const Acts::BoundTrackParameters& pActs, const Surface& surf) const {

  Acts::GeometryContext gctx;

  const auto& pars = pActs.parameters();

  const double loc0 = pars[Acts::eBoundLoc0];
  const double loc1 = pars[Acts::eBoundLoc1];
  const double qop  = pars[Acts::eBoundQOverP];

  const Acts::Vector2 actsLoc(loc0, loc1);
  const Acts::Vector3 actsDir = pActs.direction();

  const Acts::Vector3 actsGlobal = pActs.referenceSurface().localToGlobal(gctx, actsLoc, actsDir);

  GlobalPoint gp(actsGlobal[0] * 0.1, actsGlobal[1] * 0.1, actsGlobal[2] * 0.1);

  LocalPoint lp = surf.toLocal(gp);

  GlobalVector dirG(actsDir[0], actsDir[1], actsDir[2]);
  LocalVector dirL = surf.toLocal(dirG);

  const double uz = dirL.z();

  // if (std::abs(uz) < 1e-2) {
  //   std::cout << "[BAD LOCAL PARAMS] uz=" << uz
  //             << " dirL=(" << dirL.x() << "," << dirL.y() << "," << dirL.z() << ")"
  //             << " surf pos=" << surf.position()
  //             << "\n";
  // }

  const double uz_safe = (std::abs(uz) < 1e-12) ? std::copysign(1e-12, uz == 0. ? 1. : uz) : uz;

  const double dxdz = dirL.x() / uz_safe;
  const double dydz = dirL.y() / uz_safe;

  const float pzSign = (uz >= 0.) ? +1.f : -1.f;

  return LocalTrajectoryParameters(qop, dxdz, dydz, lp.x(), lp.y(), pzSign, true);
}


// LocalTrajectoryParameters ActsPropagator::GetLocalTrajectoryParameters(
//     const Acts::BoundTrackParameters& pActs,
//     const Surface& surf) const {

//   const auto& pars = pActs.parameters();

//   const double loc0 = pars[Acts::eBoundLoc0];
//   const double loc1 = pars[Acts::eBoundLoc1];
//   const double qop  = pars[Acts::eBoundQOverP];

//   // ACTS local position [mm] -> CMSSW [cm]
//   LocalPoint lp(loc0 * 0.1, loc1 * 0.1);

//   // Usa la direzione ACTS finale direttamente
//   const auto dirA = pActs.direction();
//   GlobalVector dirG(dirA[0], dirA[1], dirA[2]);

//   // Porta la direzione nel frame locale CMSSW della surface target
//   LocalVector dirL = surf.toLocal(dirG);

//   // Protezione da uz ~ 0
//   const double uz = dirL.z();
//   const double uz_safe = (std::abs(uz) < 1e-12) ? std::copysign(1e-12, (uz == 0. ? 1. : uz)) : uz;

//   const double dxdz = dirL.x() / uz_safe;
//   const double dydz = dirL.y() / uz_safe;

//   const float pzSign = (uz >= 0.) ? +1.f : -1.f;

//   // std::cout << "\n[DBG LTP]"
//   //         << " loc=(" << lp.x() << "," << lp.y() << ")"
//   //         << " dirG=(" << dirG.x() << "," << dirG.y() << "," << dirG.z() << ")"
//   //         << " dirL=(" << dirL.x() << "," << dirL.y() << "," << dirL.z() << ")"
//   //         << " dxdz=" << dxdz
//   //         << " dydz=" << dydz
//   //         << " pzSign=" << pzSign
//   //         << "\n";

//   // const GlobalPoint g0  = surf.toGlobal(LocalPoint(0., 0.));
//   // const GlobalPoint g1x = surf.toGlobal(LocalPoint(1., 0.));
//   // const GlobalPoint g1y = surf.toGlobal(LocalPoint(0., 1.));

//   // GlobalVector exCms = (g1x - g0).unit();
//   // GlobalVector eyCms = (g1y - g0).unit();
//   // GlobalVector nCms  = exCms.cross(eyCms).unit();

//   // std::cout << "\n[CMSSW normal] "
//   //           << nCms.x() << " "
//   //           << nCms.y() << " "
//   //           << nCms.z() << "\n";

//   // const auto& refSurfACTS = pActs.referenceSurface();
//   // auto nActs = refSurfACTS.normal(Acts::GeometryContext{}, Acts::Vector3{loc0, loc1, 0.}, pActs.direction());
//   // std::cout << "\n[ACTS normal] "
//   //           << nActs[0] << " "
//   //           << nActs[1] << " "
//   //           << nActs[2] << "\n";

//   return LocalTrajectoryParameters(qop, dxdz, dydz, lp.x(), lp.y(), pzSign, true);
// }

const Local2DPoint center(0.,0.); 
const Local3DPoint locz(0.,0.,1.);
const Local3DPoint locx(1.,0.,0.);
const Local3DPoint locy(0.,1.,0.);
const GlobalPoint origin(0.,0.,0.);
std::shared_ptr<const Acts::Surface> ActsPropagator::buildTargetSurf(const Surface& cmssw_surf) const {
  // ===== Define the transformation (i.e. Rotation and Translation) of the CMSSW surface =====
  Acts::Transform3 t = Acts::Transform3::Identity();
  Acts::RotationMatrix3 R;

  GlobalPoint position = cmssw_surf.toGlobal(center);
  GlobalPoint zpos = cmssw_surf.toGlobal(locz);
  GlobalPoint xpos = cmssw_surf.toGlobal(locx);
  GlobalPoint ypos = cmssw_surf.toGlobal(locy);
  GlobalVector dz = zpos - position;
  GlobalVector dx = xpos - position;
  GlobalVector dy = ypos - position;

  Eigen::Vector3d dxV(dx.x(), dx.y(), dx.z());
  Eigen::Vector3d dyV(dy.x(), dy.y(), dy.z());
  Eigen::Vector3d dVz(dz.x(), dz.y(), dz.z());

  dxV.normalize();
  dyV.normalize();
  dxV = dxV - dxV.dot(dVz) * dVz; // make dxV ortogonal to dz
  dVz.normalize();

  dyV = dVz.cross(dxV);  // make dy ortogonal to dx e dz

  Eigen::Matrix3d Rot;
  Rot.col(0) = dxV;
  Rot.col(1) = dyV;
  Rot.col(2) = dVz;

  t.prerotate(Rot);
  t.pretranslate(Acts::Vector3(position.x()*10, position.y()*10, position.z()*10)); // from cm to mm


  // ===== Define the ACTS surface considering two tipes of bounds (i.e. Rectangle and Trapezoid) =====
  std::shared_ptr<Acts::Surface> acts_surf = nullptr;
  bool boundFound = false;
  auto rect = dynamic_cast<const RectangularPlaneBounds*>(&cmssw_surf.bounds());
  if (rect){
    boundFound = true;
    ACTS_VERBOSE("Building a Rectangular surface");
    const std::size_t kValues = Acts::RectangleBounds::BoundValues::eSize;
    std::array<double, kValues> bValues{};
    std::vector<double> bVector = {-rect->width()  / 2 * 10,  // cm → mm
                                    -rect->length() / 2 * 10,
                                    rect->width()  / 2 * 10,
                                    rect->length() / 2 * 10};

    std::copy_n(bVector.begin(), kValues, bValues.begin());
    acts_surf = Acts::Surface::makeShared<Acts::PlaneSurface>(t, std::move(std::make_shared<const Acts::RectangleBounds>(bValues)));
  }
  auto trap = dynamic_cast<const TrapezoidalPlaneBounds*>(&cmssw_surf.bounds());
  if (trap){
    boundFound = true;
    ACTS_VERBOSE("Building a Trapezoidal surface");
    auto makeActsTransformFromCmsswSurfaceMm = [&](const Surface& s) -> Acts::Transform3 {
      // Build an ACTS transform whose local frame is *exactly* the CMSSW LocalPoint frame.
      // We use 1 cm steps in CMSSW local coordinates to define the basis.
      GlobalPoint g0 = s.toGlobal(LocalPoint(0., 0., 0.));
      GlobalPoint gx = s.toGlobal(LocalPoint(1., 0., 0.));
      GlobalPoint gy = s.toGlobal(LocalPoint(0., 1., 0.));
      GlobalPoint gz = s.toGlobal(LocalPoint(0., 0., 1.));

      Acts::Vector3 O(g0.x() * 10., g0.y() * 10., g0.z() * 10.); // cm->mm
      Acts::Vector3 X((gx.x() - g0.x()) * 10., (gx.y() - g0.y()) * 10., (gx.z() - g0.z()) * 10.);
      Acts::Vector3 Y((gy.x() - g0.x()) * 10., (gy.y() - g0.y()) * 10., (gy.z() - g0.z()) * 10.);
      Acts::Vector3 Z((gz.x() - g0.x()) * 10., (gz.y() - g0.y()) * 10., (gz.z() - g0.z()) * 10.);

      // Normalize and enforce orthonormal, right-handed frame (numerically robust)
      X.normalize();
      // Make Y orthogonal to X
      Y = (Y - (Y.dot(X)) * X);
      Y.normalize();
      // Recompute Z from X x Y to guarantee handedness
      Z = X.cross(Y);
      Z.normalize();

      Acts::Transform3 t = Acts::Transform3::Identity();
      t.linear().col(0) = X;
      t.linear().col(1) = Y;
      t.linear().col(2) = Z;
      t.translation()   = O;
      return t;
    };

    // 1) Read CMSSW trapezoid parameters (they are half-lengths in cm)
    auto params = trap->parameters();
    std::cout << "7" << std::endl;
    double halfBottom = params[0] * 10.;  // mm  (BOTTOM edge at y=-halfY in CMSSW local)
    double halfTop    = params[1] * 10.;  // mm  (TOP    edge at y=+halfY in CMSSW local)
    double halfY      = params[3] * 10.;  // mm

    if (halfY <= 0.) {
      throw std::runtime_error("Invalid trapezoid halfY (apothem) <= 0");
    }

    // 2) Build ACTS transform from CMSSW surface frame (this removes all X/Y sign ambiguities)
    Acts::Transform3 tFixed = makeActsTransformFromCmsswSurfaceMm(cmssw_surf);

    // 3) Bounds: map bottom/top to negY/posY
    std::array<double, Acts::TrapezoidBounds::BoundValues::eSize> b{};
    b[Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY] = halfBottom; // y = -halfY
    b[Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY] = halfTop;    // y = +halfY
    b[Acts::TrapezoidBounds::BoundValues::eHalfLengthY]     = halfY;
    b[Acts::TrapezoidBounds::BoundValues::eRotationAngle]   = 0.0;

    auto bounds = std::make_shared<const Acts::TrapezoidBounds>(b);
    acts_surf = Acts::Surface::makeShared<Acts::PlaneSurface>(tFixed, bounds);
  } 

  auto cylinder = dynamic_cast<const SimpleCylinderBounds*>(&cmssw_surf.bounds());
  if (cylinder) {
    boundFound = true;
    ACTS_VERBOSE("Building a Cylinder surface");
    // CMSSW: cm
    const double rMax_cm = 0.5 * cylinder->width();
    const double rMin_cm = rMax_cm - cylinder->thickness();
    const double halfZ_cm = 0.5 * cylinder->length();
    // ACTS: mm
    const double targetR_mm = 10.0 * rMin_cm;
    const double halfZ_mm = 10.0 * halfZ_cm;

    const auto& pos = cmssw_surf.position();

    Acts::Transform3 trf = Acts::Transform3::Identity();
    trf.translation() = Acts::Vector3(10.0 * pos.x(), 10.0 * pos.y(), 10.0 * pos.z());

    acts_surf = Acts::Surface::makeShared<Acts::CylinderSurface>(trf, targetR_mm, halfZ_mm);
  }

  if(!boundFound){
    ACTS_VERBOSE("Surface bounds not found! Building a surface with custom rectangular bounds.");
    const std::size_t kValues = Acts::RectangleBounds::BoundValues::eSize;
    std::array<double, kValues> bValues{};
    std::vector<double> bVector = {-50 * Acts::UnitConstants::mm,  // cm → mm
                                   -50 * Acts::UnitConstants::mm,
                                    50 * Acts::UnitConstants::mm,
                                    50 * Acts::UnitConstants::mm};

    std::copy_n(bVector.begin(), kValues, bValues.begin());
    acts_surf = Acts::Surface::makeShared<Acts::PlaneSurface>(t, std::move(std::make_shared<const Acts::RectangleBounds>(bValues)));   
  }

  auto t_final = acts_surf->transform(Acts::GeometryContext{});
  auto b_values = acts_surf->bounds().values();
  ACTS_VERBOSE("Parameters of the Target surface:");
  ACTS_VERBOSE("Surface type: " << acts_surf->type());
  ACTS_VERBOSE("Surface Position in acts [mm]: " << t_final.translation().transpose());
  ACTS_VERBOSE("Bounds: ");
  for (size_t i = 0; i < b_values.size(); ++i) {
      ACTS_VERBOSE(" " << b_values[i]);
  }

  return acts_surf;
}

TsosWP ActsPropagator::tsosWithPathFromActs(const Acts::BoundTrackParameters& start_param, const Surface& surfTarget, Acts::Direction actsDir) const {
  auto makeInvalid = []() -> TsosWP {
    return TsosWP{TrajectoryStateOnSurface(), -1.0};
  };

  // ===== Convert the target surface from CMSSW into an ACTS's surface =====
  ACTS_VERBOSE("Converting Target Surface from CMSSW to ACTS...");
  SurfaceConverters surfConv(trkGeo_cmssw_, *trkGeo_and_DetEls_, m_Level);
  std::shared_ptr<const Acts::Surface> surf_ACTS_target = nullptr;
  surf_ACTS_target = surfConv.fromCMSSWtoACTS(surfTarget);
  if(!surf_ACTS_target) {
    // ===== If the match is not found, build a surface using CMSSW parameters =====
    ACTS_VERBOSE("[WARNING] Failed to build TSOS, Target Surface not found. Building a new one from cmssw parameters and continuing with that.");
    surf_ACTS_target = buildTargetSurf(surfTarget);
    nNoSurf += 1;
  }

  // ===== Propagate to the target surface =====
  PropagationAlgorithm_Config cfg;
  cfg.propDir = actsDir;
  cfg.covarianceTransport = true;
  std::shared_ptr<const Acts::Logger> prop_logger = Acts::getDefaultLogger("Concrete Propagator", m_Level);
  ACTS_VERBOSE("Executing the actual propagation...");
  auto propResOpt = ConcProp->execute(cfg, *prop_logger, start_param, *surf_ACTS_target); 

  if (!propResOpt) {
    nNoProp += 1;
    ACTS_VERBOSE("ACTS propagation failed to this surface; returning invalid TSOS.");
    return makeInvalid();
  }

  const auto& [bParam, length] = *propResOpt;

  const Surface* outputSurface = &surfTarget;
  Plane::PlanePointer tangentPlane;

  if (dynamic_cast<const Cylinder*>(&surfTarget) != nullptr) {
    ACTS_VERBOSE("Target CMSSW surface is a Cylinder: building tangent Plane for output TSOS.");
    Acts::GeometryContext gctx;
    const auto& pars = bParam.parameters();

    Acts::Vector2 actsLoc(pars[Acts::eBoundLoc0], pars[Acts::eBoundLoc1]);
    Acts::Vector3 actsDir = bParam.direction();
    Acts::Vector3 actsGlobal = bParam.referenceSurface().localToGlobal(gctx, actsLoc, actsDir);

    GlobalPoint gp(actsGlobal[0] * 0.1, actsGlobal[1] * 0.1, actsGlobal[2] * 0.1);

    const double gpR = std::sqrt(gp.x() * gp.x() + gp.y() * gp.y());
    const double gpPhi = std::atan2(gp.y(), gp.x());

    ACTS_VERBOSE("Cylinder output ACTS local [mm]: loc0=" << actsLoc[0] << " loc1=" << actsLoc[1]);
    ACTS_VERBOSE("Cylinder output ACTS global [mm]: " << actsGlobal[0] << " " << actsGlobal[1] << " " << actsGlobal[2]);
    ACTS_VERBOSE("Cylinder output CMSSW global [cm]: " << gp.x() << " " << gp.y() << " " << gp.z() << " R=" << gpR << " phi=" << gpPhi);

    GlobalVector zAxis(gp.x(), gp.y(), 0.);
    if (zAxis.mag2() < 1e-12) {
      ACTS_VERBOSE("Cylinder tangent normal is ill-defined near beamline; using fallback normal (1,0,0).");
      zAxis = GlobalVector(1., 0., 0.);
    } else {
      zAxis = zAxis.unit();
    }

    GlobalVector yAxis(0., 0., 1.);
    GlobalVector xAxis = yAxis.cross(zAxis).unit();
    yAxis = zAxis.cross(xAxis).unit();

    ACTS_VERBOSE("Tangent plane basis:" << " xAxis=(" << xAxis.x() << "," << xAxis.y() << "," << xAxis.z() << ")" << " yAxis=(" << yAxis.x() << "," << yAxis.y() << "," << yAxis.z() << ")" << " zAxis=(" << zAxis.x() << "," << zAxis.y() << "," << zAxis.z() << ")");

    Surface::RotationType rot(
        xAxis.x(), yAxis.x(), zAxis.x(),
        xAxis.y(), yAxis.y(), zAxis.y(),
        xAxis.z(), yAxis.z(), zAxis.z()
    );

    tangentPlane = Plane::build(gp, rot, new RectangularPlaneBounds(50., 50., 0.));

    ACTS_VERBOSE("Built tangent Plane at: " << tangentPlane->position());
    ACTS_VERBOSE("Tangent Plane check toGlobal(0,0): " << tangentPlane->toGlobal(LocalPoint(0., 0.)));
    ACTS_VERBOSE("Tangent Plane local position of ACTS final point: " << tangentPlane->toLocal(gp));

    outputSurface = tangentPlane.get();
  } else {
    ACTS_VERBOSE("Target CMSSW surface is not a Cylinder: using original target surface for output TSOS.");
  }

  auto ltp = GetLocalTrajectoryParameters(bParam, *outputSurface);

  const double phi   = bParam.parameters()[Acts::eBoundPhi];
  const double theta = bParam.parameters()[Acts::eBoundTheta];
  auto covOpt = bParam.covariance();
  if (!covOpt) {
    nNoCov += 1;
    throw cms::Exception("ACTSPropagation") << "No covariance in ACTS parameters";
  }
  const Acts::BoundMatrix& cov = *covOpt;

  LocalTrajectoryError lte(convertCovACTStoCMSSW(*outputSurface, phi, theta, cov));

  SurfaceSideDefinition::SurfaceSide side = actsDir == Acts::Direction::Forward() ? SurfaceSideDefinition::beforeSurface : SurfaceSideDefinition::afterSurface;

  TrajectoryStateOnSurface tsos_fromACTS(ltp, lte, *outputSurface, field, side);
  TsosWP newTsosWP_FromACTS(tsos_fromACTS, length * 0.1);

  auto pos = tsos_fromACTS.globalPosition();
  double x = pos.x();
  double y = pos.y();
  double z = pos.z();
  double r   = std::sqrt(x*x + y*y);
  double phi_dbg = std::atan2(y, x);
  auto pDir = "Forward";
  if (actsDir == Acts::Direction::Backward()) 
    pDir = "Backward";

  ACTS_VERBOSE("F_tsos parameters-> Pos: X " << x << " Y " << y << " Z " << z << " R " << r << " phi " << phi_dbg << " q: " << tsos_fromACTS.globalParameters().charge() << 
               " propDir: " << pDir << 
               " sigedInvMom: " << tsos_fromACTS.globalParameters().signedInverseMomentum() <<
               " signedInvTransvMom: " << tsos_fromACTS.globalParameters().signedInverseTransverseMomentum());
  // ACTS_VERBOSE("Local Error: ");
  // printMatrix(tsos_fromACTS.localError().matrix());

  ACTS_VERBOSE("Final tsos position: " << tsos_fromACTS.globalPosition());

  return newTsosWP_FromACTS;
}

TrajectoryStateOnSurface ActsPropagator::boundFreeTrajectoryState(const FreeTrajectoryState &fts) const {

  ACTS_VERBOSE("Bound Input FreeTrajectoryState to a planar start surface");
  // Get momentum and position of the fts
  GlobalPoint  pos = fts.position();

  // Generate the ortogonal local frame system of the plane
  /// NOTE: as defined in https://github.com/cms-sw/cmssw/blob/7a8df40ec95f56efec1c0596201f0827d0fa49ff/TrackingTools/PatternTools/src/TransverseImpactPointExtrapolator.cc#L98
  //x is perpendicular to the momentum
  GlobalVector xLocal = GlobalVector(-fts.momentum().y(), fts.momentum().x(), 0).unit();
  //y along global Z
  GlobalVector yLocal(0., 0., 1.);
  //z accordingly
  GlobalVector zLocal(xLocal.cross(yLocal));
  const Surface::PositionType& surfPos(pos);
  Surface::RotationType rotation(xLocal, yLocal, zLocal);
  ReferenceCountingPointer<Plane> surface = PlaneBuilder().plane(surfPos, rotation);

  // Bound the fts to the plane
  TrajectoryStateOnSurface newTsos(fts, *surface);

  // {
  //   LocalVector dL = newTsos.surface().toLocal(newTsos.globalDirection().unit());
  //   std::cout << "\n[BOUND FTS CHECK]\n"
  //             << "  global dir=" << newTsos.globalDirection() << "\n"
  //             << "  dirL=(" << dL.x() << "," << dL.y() << "," << dL.z() << ")\n"
  //             << "  tx=" << newTsos.localParameters().dxdz()
  //             << " ty=" << newTsos.localParameters().dydz() << "\n";
  // }

  // Sanity checks
  if (!newTsos.isValid()) {
    throw cms::Exception("InvalidTSOS") << "TrajectoryStateOnSurface is not valid";
  }

  if (!newTsos.hasError()) {
    throw cms::Exception("MissingError")<< "TrajectoryStateOnSurface has no error";
  }

  return newTsos;
}

//
////////////////////////////////////////////////////////////////////////////
//

/** The methods propagateWithPath() are identical to the corresponding
 *  methods propagate() in what concerns the resulting
 *  TrajectoryStateOnSurface, but they provide in addition the
 *  exact path length along the trajectory.
 */

std::shared_ptr<const Acts::Surface> ActsPropagator::inflateStartPlaneIfRectOrTrap(const Acts::Surface& targetSurf,
                                                                                   const Acts::GeometryContext& gctx,
                                                                                   double marginMm) const {
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
      return Acts::Surface::makeShared<Acts::PlaneSurface>(trf, newBounds);
    }

    // TrapezoidBounds
    if (bType == Acts::SurfaceBounds::BoundsType::eTrapezoid) {
      const auto* trap = dynamic_cast<const Acts::TrapezoidBounds*>(&targetSurf.bounds());
      if (!trap) return nullptr; // safety

      // In Acts TrapezoidBounds è descritto da:
      // halfXnegY, halfXposY, halfY  (lato corto a y negativo). :contentReference[oaicite:1]{index=1}
      const double halfXnegY =
          trap->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY)* 1.17;
      const double halfXposY =
          trap->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY)* 1.17;
      const double halfY =
          trap->get(Acts::TrapezoidBounds::BoundValues::eHalfLengthY)* 1.17;
      const double ang =
          trap->get(Acts::TrapezoidBounds::BoundValues::eRotationAngle);

      auto newBounds = std::make_shared<Acts::TrapezoidBounds>(halfXnegY, halfXposY, halfY, ang);

      return Acts::Surface::makeShared<Acts::PlaneSurface>(trf, newBounds);
    }

    return nullptr;
  }

std::pair<TrajectoryStateOnSurface, double> ActsPropagator::propagateWithPath(const FreeTrajectoryState &fts,
                                                                              const Plane &pDest) const {                                                                         
  // ===== Get parameters from FreeTrajectoryState =====
  ACTS_VERBOSE("Method called: 1");

  // ===== Convert the fts in a tsos =====
  auto tsos = boundFreeTrajectoryState(fts);

  // ===== Call the method which uses tsos and plane =====
  auto newTsosWP_FromACTS = propagateWithPath(tsos, pDest);

  return newTsosWP_FromACTS;
}

std::pair<TrajectoryStateOnSurface, double> ActsPropagator::propagateWithPath(const FreeTrajectoryState &fts,
                                                                              const Cylinder &cDest) const {
  // ===== Get parameters from FreeTrajectoryState =====
  ACTS_VERBOSE("Method called: 2");

  // ===== Convert the fts in a tsos =====
  auto tsos = boundFreeTrajectoryState(fts);

  // ===== Call the method which uses tsos and plane =====
  auto newTsosWP_FromACTS = propagateWithPath(tsos, cDest);

  return newTsosWP_FromACTS;


  // // ################## TO BE CHANGED WHEN INCLUDING OTHER RECO STEPS ##################                                                                              
  // // ===== Get parameters from FreeTrajectoryState =====
  // ACTS_VERBOSE("Method called: 2");
  // GlobalPoint gPoint = fts.position();
  // GlobalVector gMomentum = fts.momentum();
  // double qOverP = fts.signedInverseMomentum();
  // // const auto& curvErr = fts.curvilinearError(); // TEMPORARY
  // // AlgebraicSymMatrix55 covMat_cmssw_init = curvErr.matrix(); // TEMPORARY  

  // Acts::Vector4 pos4 = {gPoint.x()*10, gPoint.y()*10, gPoint.z()*10, 0};
  // Acts::Vector3 dir = {gMomentum.x(), gMomentum.y(), gMomentum.z()};
  // dir.normalize();
  // // Acts::BoundMatrix covMat_acts_init = convertCovCMSSWtoACTS(covMat_cmssw_init);
  // Acts::BoundMatrix covMat_acts_init = Acts::BoundMatrix::Identity() * 1e-3; // TEMPORARY

  // auto start_param = Acts::BoundTrackParameters::createCurvilinear(pos4, dir, qOverP, covMat_acts_init, Acts::ParticleHypothesis::muon()); // TEMPORARY

  // return tsosWithPathFromActs(start_param, cDest, actsPropDir_);
}

std::pair<TrajectoryStateOnSurface, double> ActsPropagator::propagateWithPath(const TrajectoryStateOnSurface &tsos, 
                                                                              const Plane &pDest) const {


  /// NOTE: this is the only method called during refit 
  // ===== Get parameters from FreeTrajectoryState =====
  nTot += 1;
  ACTS_VERBOSE("Method called: 3");
  ACTS_VERBOSE("Initial tsos position: " << tsos.globalPosition());
  ACTS_VERBOSE("Initial tsos direction: " << tsos.globalDirection());
  ACTS_VERBOSE("Target surface center: " << pDest.position());
  ACTS_VERBOSE("The initial distance to target surface is " << pDest.localZ(tsos.globalPosition()));
  ACTS_VERBOSE("Initial q/P: " << tsos.signedInverseMomentum());
  ACTS_VERBOSE("Initial q/P (from globalParameters): " << tsos.globalParameters().signedInverseMomentum());
  ACTS_VERBOSE("Initial q/Pt (from globalParameters): " << tsos.globalParameters().signedInverseTransverseMomentum());
  ACTS_VERBOSE("Initial charge (from globalParameters): " << tsos.globalParameters().charge());

  // std::cout << tsos.localParameters().pzSign() << std::endl;

  GlobalPoint gPoint = tsos.globalPosition();
  GlobalVector gDir = tsos.globalDirection();
  double qOverP = tsos.signedInverseMomentum();
  const auto& localErr = tsos.localError();
  AlgebraicSymMatrix55 covMat_cmssw_init = localErr.matrix(); 
  
  // DEBUG
  {
    const auto& lp = tsos.localParameters();
    GlobalVector dG = tsos.globalDirection().unit();
    LocalVector dL = tsos.surface().toLocal(dG);

    const double uz = dL.z();
    const double tx = lp.dxdz();
    const double ty = lp.dydz();

    // if (std::abs(uz) < 5e-2 || std::abs(tx) > 10. || std::abs(ty) > 10.) {
    //   std::cout << "\n[BAD INPUT TSOS LOCAL PARAMS]\n"
    //             << "  tx=" << tx << " ty=" << ty << "\n"
    //             << "  uz=" << uz
    //             << " dirL=(" << dL.x() << "," << dL.y() << "," << dL.z() << ")\n"
    //             << "  global pos=" << tsos.globalPosition() << "\n"
    //             << "  global dir=" << tsos.globalDirection() << "\n"
    //             << "  surface pos=" << tsos.surface().position() << "\n"
    //             << "  target pos=" << pDest.position() << "\n";
    // }
  }
  // END DEBUG


  // ===== Use the CMSSW parameters to define the ACTS ones =====
  Acts::Vector4 pos4 = {gPoint.x()*10, gPoint.y()*10, gPoint.z()*10, 0};
  Acts::Vector3 dir = {gDir.x(), gDir.y(), gDir.z()};
  dir.normalize();
  // Acts::BoundMatrix covMat_acts_init = convertCovCMSSWtoACTS(tsos, covMat_cmssw_init);

  // ===== Create the Intial Parameters for ACTS =====
  ACTS_VERBOSE("Creating ACTS initial parameters...");
  SurfaceConverters surfConv(trkGeo_cmssw_, *trkGeo_and_DetEls_, m_Level);
  ACTS_VERBOSE("Converting Initial Surface from CMSSW to ACTS...");                                                                             
  auto startToUse = surfConv.fromCMSSWtoACTS(tsos.surface());

  if(!startToUse) {
    // ===== If the match is not found, build a surface using CMSSW parameters =====
    /// NOTE: This method is used essentially in the case a fts is given 
    startToUse = buildTargetSurf(tsos.surface());
  }

  // DEBUG
  // // ===== CMSSW local frame =====
  // const auto& surf = tsos.surface();
  // GlobalPoint g0_cm   = surf.toGlobal(LocalPoint(0.,0.));
  // GlobalPoint g1x_cm  = surf.toGlobal(LocalPoint(1.,0.));
  // GlobalPoint g1y_cm  = surf.toGlobal(LocalPoint(0.,1.));

  // GlobalVector ex_cm = (g1x_cm - g0_cm).unit();
  // GlobalVector ey_cm = (g1y_cm - g0_cm).unit();
  // GlobalVector ez_cm = ex_cm.cross(ey_cm).unit();

  // std::cout << "\n[DBG SURF 1] CMSSW FRAME\n";
  // std::cout << "origin_cm = " << g0_cm.x() << " " << g0_cm.y() << " " << g0_cm.z() << "\n";
  // std::cout << "ex_cm     = " << ex_cm.x() << " " << ex_cm.y() << " " << ex_cm.z() << "\n";
  // std::cout << "ey_cm     = " << ey_cm.x() << " " << ey_cm.y() << " " << ey_cm.z() << "\n";
  // std::cout << "ez_cm     = " << ez_cm.x() << " " << ez_cm.y() << " " << ez_cm.z() << "\n";

  // // ===== ACTS local frame =====
  // auto T = startToUse->transform(Acts::GeometryContext{});
  // auto Oa = T.translation();
  // auto ex_a = T.rotation().col(0);
  // auto ey_a = T.rotation().col(1);
  // auto ez_a = T.rotation().col(2);

  // std::cout << "\n[DBG SURF 2] ACTS FRAME\n";
  // std::cout << "origin_a [mm] = " << Oa.transpose() << "\n";
  // std::cout << "ex_a          = " << ex_a.transpose() << "\n";
  // std::cout << "ey_a          = " << ey_a.transpose() << "\n";
  // std::cout << "ez_a          = " << ez_a.transpose() << "\n";

  // std::cout << "\n[DBG SURF 3] FRAME OVERLAPS\n";
  // std::cout << "ex_cm·ex_a = " << (ex_cm.x()*ex_a(0) + ex_cm.y()*ex_a(1) + ex_cm.z()*ex_a(2)) << "\n";
  // std::cout << "ey_cm·ey_a = " << (ey_cm.x()*ey_a(0) + ey_cm.y()*ey_a(1) + ey_cm.z()*ey_a(2)) << "\n";
  // std::cout << "ez_cm·ez_a = " << (ez_cm.x()*ez_a(0) + ez_cm.y()*ez_a(1) + ez_cm.z()*ez_a(2)) << "\n";

  // std::cout << "ex_cm·ey_a = " << (ex_cm.x()*ey_a(0) + ex_cm.y()*ey_a(1) + ex_cm.z()*ey_a(2)) << "\n";
  // std::cout << "ex_cm·ez_a = " << (ex_cm.x()*ez_a(0) + ex_cm.y()*ez_a(1) + ex_cm.z()*ez_a(2)) << "\n";
  // std::cout << "ey_cm·ex_a = " << (ey_cm.x()*ex_a(0) + ey_cm.y()*ex_a(1) + ey_cm.z()*ex_a(2)) << "\n";
  // END DEBUG

  Acts::BoundMatrix covMat_acts_init = convertCovCMSSWtoACTS(tsos, covMat_cmssw_init);

  //DEBUG
  // std::cout << "\n[DBG PHYS COV]\n";

  // std::cout << "Var loc0 = "
  // << covMat_acts_init(Acts::eBoundLoc0,Acts::eBoundLoc0) << "\n";

  // std::cout << "Var loc1 = "
  // << covMat_acts_init(Acts::eBoundLoc1,Acts::eBoundLoc1) << "\n";

  // std::cout << "Var phi = "
  // << covMat_acts_init(Acts::eBoundPhi,Acts::eBoundPhi) << "\n";

  // std::cout << "Var theta = "
  // << covMat_acts_init(Acts::eBoundTheta,Acts::eBoundTheta) << "\n";

  // std::cout << "Cov(phi,theta)= "
  // << covMat_acts_init(Acts::eBoundPhi,Acts::eBoundTheta) << "\n";                                                                              
  //END DEBUG

  Acts::GeometryContext gctx;
  auto t_initial = startToUse->transform(Acts::GeometryContext{});
  auto b_values = startToUse->bounds().values();
  ACTS_VERBOSE("Position of the initial surface [mm]: " << t_initial.translation().transpose());
  ACTS_VERBOSE("Bounds of initial surface:");
  for (size_t i = 0; i < b_values.size(); ++i) {
    ACTS_VERBOSE(" " << b_values[i]);
  }
  const auto* startplane = dynamic_cast<const Plane*>(&tsos.surface());
  ACTS_VERBOSE("Distance from initial surface and initial position " << startplane->localZ(tsos.globalPosition()));
  if(!startToUse) {
    ACTS_VERBOSE("Failed to find initial surface");
  }

  auto inflated = inflateStartPlaneIfRectOrTrap(*startToUse, Acts::GeometryContext{}, 9 /*mm*/);
  if (inflated) {
    startToUse = inflated;
    ACTS_VERBOSE("Using inflated start bounds");
  }

  // ===== Define the propagation direction =====
  if(actsPropDir_ == Acts::Direction::Forward()){
    ACTS_VERBOSE("Propagation Direction FORWARD");
  } else if(actsPropDir_ == Acts::Direction::Backward()) {
    ACTS_VERBOSE("Propagation Direction BACKWARD");
  } else {
    throw cms::Exception("ACTSPropagation") << "Propagation direction not found";
  }
  const auto gpos = pos4.segment<3>(0);
  ACTS_VERBOSE("Initial position for precheck: " << gpos.transpose());
  ACTS_VERBOSE("Initial direction for precheck: " << dir.transpose());
  ACTS_VERBOSE("Initial surface center " << startToUse->center(Acts::GeometryContext{})[0] << " " << startToUse->center(Acts::GeometryContext{})[1] << " " << startToUse->center(Acts::GeometryContext{})[2]);
  ACTS_VERBOSE("Bounds of initial surface:");
  auto b_valuesInflated = startToUse->bounds().values();
  for (size_t i = 0; i < b_valuesInflated.size(); ++i) {
    ACTS_VERBOSE(" " << b_valuesInflated[i]);
  }
  const auto& Tr = startToUse->transform(gctx);
  Acts::Vector3 pLocal3 = Tr.inverse() * gpos;
  ACTS_VERBOSE("manual local x,y,z = " << pLocal3.transpose());
  

  auto res_StartParam = Acts::BoundTrackParameters::create(Acts::GeometryContext{}, 
                                                           startToUse, 
                                                           pos4, dir, qOverP,
                                                           covMat_acts_init,
                                                           Acts::ParticleHypothesis::muon(),
                                                           1e-2);

  if (!res_StartParam.ok()) {
    // FallBack
    ACTS_VERBOSE("Failed to build ACTS initial parameters! Trying to further inflate the surface bounds");
    auto more_inflated = inflateStartPlaneIfRectOrTrap(*startToUse, Acts::GeometryContext{}, 90000 /*mm*/);
    auto b_values = more_inflated->bounds().values();
    ACTS_VERBOSE("New Bounds:");
    for (size_t i = 0; i < b_values.size(); ++i) {
      ACTS_VERBOSE(" " << b_values[i]);
    }
    res_StartParam = Acts::BoundTrackParameters::create(Acts::GeometryContext{}, 
                                                        more_inflated, 
                                                        pos4, dir, qOverP,
                                                        covMat_acts_init,
                                                        Acts::ParticleHypothesis::muon(),
                                                        0.1);
    if(!res_StartParam.ok()) {
      ACTS_VERBOSE("Even with super inflated bounds, failed to build ACTS initial parameters. I quit.");
      throw cms::Exception("ACTSPropagator") << " failed to build ACTS initial parameters: " << res_StartParam.error();
    }
  }
  // DEBUG
  // {
  //   ComputeLocalBoundJacobian computeJ;
  //   const auto J = computeJ.FromCMSSWtoACTS(tsos);

  //   const auto& surf = tsos.surface();
  //   const auto& lp = tsos.localParameters();

  //   const double tx0 = lp.dxdz();
  //   const double ty0 = lp.dydz();

  //   GlobalVector dG0 = tsos.globalDirection().unit();
  //   LocalVector dL0 = surf.toLocal(dG0);
  //   const double signUz = (dL0.z() >= 0.) ? +1. : -1.;

  //   auto wrapToPi = [](double a) {
  //     while (a > M_PI) a -= 2. * M_PI;
  //     while (a < -M_PI) a += 2. * M_PI;
  //     return a;
  //   };

  //   auto makeDirFromTxTy = [&](double tx, double ty) -> Acts::Vector3 {
  //     const double s = std::sqrt(1. + tx * tx + ty * ty);

  //     LocalVector dL(signUz * tx / s,
  //                   signUz * ty / s,
  //                   signUz / s);

  //     GlobalVector dG = surf.toGlobal(dL).unit();

  //     Acts::Vector3 dirA(dG.x(), dG.y(), dG.z());
  //     dirA.normalize();
  //     return dirA;
  //   };

  //   auto makeActsParsFromTxTy = [&](double tx, double ty) {
  //     Acts::Vector3 dirPert = makeDirFromTxTy(tx, ty);

  //     return Acts::BoundTrackParameters::create(
  //         Acts::GeometryContext{},
  //         startToUse,
  //         pos4,
  //         dirPert,
  //         qOverP,
  //         covMat_acts_init,
  //         Acts::ParticleHypothesis::muon(),
  //         1e-2);
  //   };

  //   const auto& p0 = res_StartParam.value().parameters();
  //   const double phi0 = p0[Acts::eBoundPhi];
  //   const double theta0 = p0[Acts::eBoundTheta];

  //   const double h = 1e-6;

  //   auto p_tx = makeActsParsFromTxTy(tx0 + h, ty0);
  //   auto p_ty = makeActsParsFromTxTy(tx0, ty0 + h);

  //   if (!p_tx.ok() || !p_ty.ok()) {
  //     std::cout << "\n[DBG ACTS FD] Failed to build perturbed ACTS parameters\n";
  //   } else {
  //     const auto& par_tx = p_tx.value().parameters();
  //     const auto& par_ty = p_ty.value().parameters();

  //     const double dphi_dtx_fd =
  //         wrapToPi(par_tx[Acts::eBoundPhi] - phi0) / h;
  //     const double dphi_dty_fd =
  //         wrapToPi(par_ty[Acts::eBoundPhi] - phi0) / h;

  //     const double dtheta_dtx_fd =
  //         (par_tx[Acts::eBoundTheta] - theta0) / h;
  //     const double dtheta_dty_fd =
  //         (par_ty[Acts::eBoundTheta] - theta0) / h;

  //     std::cout << "\n================ ACTS FD JACOBIAN TEST ================\n";
  //     std::cout << "eta-like theta = " << theta0
  //               << " phi = " << phi0
  //               << " tx = " << tx0
  //               << " ty = " << ty0
  //               << " signUz = " << signUz << "\n";

  //     std::cout << "dphi/dtx   analytic = " << J(Acts::eBoundPhi,   1)
  //               << "   ACTS FD = " << dphi_dtx_fd
  //               << "   diff = " << (J(Acts::eBoundPhi, 1) - dphi_dtx_fd) << "\n";

  //     std::cout << "dphi/dty   analytic = " << J(Acts::eBoundPhi,   2)
  //               << "   ACTS FD = " << dphi_dty_fd
  //               << "   diff = " << (J(Acts::eBoundPhi, 2) - dphi_dty_fd) << "\n";

  //     std::cout << "dtheta/dtx analytic = " << J(Acts::eBoundTheta, 1)
  //               << "   ACTS FD = " << dtheta_dtx_fd
  //               << "   diff = " << (J(Acts::eBoundTheta, 1) - dtheta_dtx_fd) << "\n";

  //     std::cout << "dtheta/dty analytic = " << J(Acts::eBoundTheta, 2)
  //               << "   ACTS FD = " << dtheta_dty_fd
  //               << "   diff = " << (J(Acts::eBoundTheta, 2) - dtheta_dty_fd) << "\n";

  //     std::cout << "========================================================\n";
  //   }
  // }
  // END DEBUG
  

  return tsosWithPathFromActs(res_StartParam.value(), pDest, actsPropDir_);
}

std::pair<TrajectoryStateOnSurface, double> ActsPropagator::propagateWithPath(const TrajectoryStateOnSurface &tsos, 
                                                                              const Cylinder &cDest) const {
  // ################## TO BE CHANGED WHEN INCLUDING OTHER RECO STEPS ################## 
  // ===== Get parameters from FreeTrajectoryState =====
  ACTS_VERBOSE("Method called: 4");
  ACTS_VERBOSE("Initial tsos position: " << tsos.globalPosition());
  ACTS_VERBOSE("Initial tsos direction: " << tsos.globalDirection());
  ACTS_VERBOSE("Target surface center: " << cDest.position());
  ACTS_VERBOSE("Initial q/P: " << tsos.signedInverseMomentum());
  ACTS_VERBOSE("Initial q/P (from globalParameters): " << tsos.globalParameters().signedInverseMomentum());
  ACTS_VERBOSE("Initial q/Pt (from globalParameters): " << tsos.globalParameters().signedInverseTransverseMomentum());
  ACTS_VERBOSE("Initial charge (from globalParameters): " << tsos.globalParameters().charge());

  // std::cout << tsos.localParameters().pzSign() << std::endl;

  GlobalPoint gPoint = tsos.globalPosition();
  GlobalVector gDir = tsos.globalDirection();
  double qOverP = tsos.signedInverseMomentum();
  const auto& localErr = tsos.localError();
  AlgebraicSymMatrix55 covMat_cmssw_init = localErr.matrix();  

  // ===== Use the CMSSW parameters to define the ACTS ones =====
  Acts::Vector4 pos4 = {gPoint.x()*10, gPoint.y()*10, gPoint.z()*10, 0};
  Acts::Vector3 dir = {gDir.x(), gDir.y(), gDir.z()};
  dir.normalize();
  // Acts::BoundMatrix covMat_acts_init = convertCovCMSSWtoACTS(tsos, covMat_cmssw_init);

  // ===== Create the Intial Parameters for ACTS =====
  ACTS_VERBOSE("Creating ACTS initial parameters...");
  SurfaceConverters surfConv(trkGeo_cmssw_, *trkGeo_and_DetEls_, m_Level);
  ACTS_VERBOSE("Converting Initial Surface from CMSSW to ACTS...");                                                                             
  auto startToUse = surfConv.fromCMSSWtoACTS(tsos.surface());

  if(!startToUse) {
    // ===== If the match is not found, build a surface using CMSSW parameters =====
    /// NOTE: This method is used essentially in the case a fts is given 
    startToUse = buildTargetSurf(tsos.surface());
  }

  Acts::BoundMatrix covMat_acts_init = convertCovCMSSWtoACTS(tsos, covMat_cmssw_init);

  Acts::GeometryContext gctx;
  auto t_initial = startToUse->transform(Acts::GeometryContext{});
  auto b_values = startToUse->bounds().values();
  ACTS_VERBOSE("Position of the initial surface [mm]: " << t_initial.translation().transpose());
  ACTS_VERBOSE("Bounds of initial surface:");
  for (size_t i = 0; i < b_values.size(); ++i) {
    ACTS_VERBOSE(" " << b_values[i]);
  }
  const auto* startplane = dynamic_cast<const Plane*>(&tsos.surface());
  ACTS_VERBOSE("Distance from initial surface and initial position " << startplane->localZ(tsos.globalPosition()));
  if(!startToUse) {
    ACTS_VERBOSE("Failed to find initial surface");
  }

  auto inflated = inflateStartPlaneIfRectOrTrap(*startToUse, Acts::GeometryContext{}, 9 /*mm*/);
  if (inflated) {
    startToUse = inflated;
    ACTS_VERBOSE("Using inflated start bounds");
  }

  // ===== Define the propagation direction =====
  if(actsPropDir_ == Acts::Direction::Forward()){
    ACTS_VERBOSE("Propagation Direction FORWARD");
  } else if(actsPropDir_ == Acts::Direction::Backward()) {
    ACTS_VERBOSE("Propagation Direction BACKWARD");
  } else {
    throw cms::Exception("ACTSPropagation") << "Propagation direction not found";
  }
  const auto gpos = pos4.segment<3>(0);
  ACTS_VERBOSE("Initial position for precheck: " << gpos.transpose());
  ACTS_VERBOSE("Initial direction for precheck: " << dir.transpose());
  ACTS_VERBOSE("Initial surface center " << startToUse->center(Acts::GeometryContext{})[0] << " " << startToUse->center(Acts::GeometryContext{})[1] << " " << startToUse->center(Acts::GeometryContext{})[2]);
  ACTS_VERBOSE("Bounds of initial surface:");
  auto b_valuesInflated = startToUse->bounds().values();
  for (size_t i = 0; i < b_valuesInflated.size(); ++i) {
    ACTS_VERBOSE(" " << b_valuesInflated[i]);
  }
  const auto& Tr = startToUse->transform(gctx);
  Acts::Vector3 pLocal3 = Tr.inverse() * gpos;
  ACTS_VERBOSE("manual local x,y,z = " << pLocal3.transpose());
  

  auto res_StartParam = Acts::BoundTrackParameters::create(Acts::GeometryContext{}, 
                                                           startToUse, 
                                                           pos4, dir, qOverP,
                                                           covMat_acts_init,
                                                           Acts::ParticleHypothesis::muon(),
                                                           1e-2);

  if (!res_StartParam.ok())
    throw cms::Exception("ACTSPropagator") << " failed to build ACTS initial parameters: " << res_StartParam.error();

  return tsosWithPathFromActs(res_StartParam.value(), cDest, actsPropDir_);
}

