#include "TrackPropagation/Acts/interface/ComputeLocalBoundJacobian.h"

#include <cmath>
#include <algorithm>

bool debugOn_J = false;

#include <iostream>
#include <iomanip>

// COORDINATE FRAME DEBUG
namespace {

template <typename TVec3>
void printVec3(const std::string& name, const TVec3& v) {
  std::cout << std::setw(20) << name << " = ("
            << std::setprecision(12)
            << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
}

void debugCompareLocalFrames(const Surface& cmsSurf,
                             const Acts::Surface& actsSurf,
                             const GlobalPoint& cmsGlobalPoint) {
  using std::cout;

  cout << "\n====================================================\n";
  cout << "[DBG FRAME] CMSSW local basis vs ACTS local basis\n";
  cout << "====================================================\n";

  // ====================================================
  // 1) Base locale CMSSW espressa in coordinate globali
  // ====================================================
  const GlobalPoint g00 = cmsSurf.toGlobal(LocalPoint(0., 0.));
  const GlobalPoint g10 = cmsSurf.toGlobal(LocalPoint(1., 0.));
  const GlobalPoint g01 = cmsSurf.toGlobal(LocalPoint(0., 1.));

  GlobalVector ex_cms = g10 - g00;  // asse x locale CMSSW, passo 1 cm
  GlobalVector ey_cms = g01 - g00;  // asse y locale CMSSW, passo 1 cm

  const double norm_ex_cms = ex_cms.mag();
  const double norm_ey_cms = ey_cms.mag();

  ex_cms *= 1.0 / norm_ex_cms;
  ey_cms *= 1.0 / norm_ey_cms;

  GlobalVector ez_cms = ex_cms.cross(ey_cms);
  ez_cms *= 1.0 / ez_cms.mag();

  printVec3("ex_cms(global)", std::array<double,3>{ex_cms.x(), ex_cms.y(), ex_cms.z()});
  printVec3("ey_cms(global)", std::array<double,3>{ey_cms.x(), ey_cms.y(), ey_cms.z()});
  printVec3("ez_cms(global)", std::array<double,3>{ez_cms.x(), ez_cms.y(), ez_cms.z()});

  cout << "norm ex_cms raw      = " << norm_ex_cms << " [should be ~1 cm]\n";
  cout << "norm ey_cms raw      = " << norm_ey_cms << " [should be ~1 cm]\n";

  // ====================================================
  // 2) Punto globale usato nel test
  // ====================================================
  Acts::Vector3 actsGlobalPoint(cmsGlobalPoint.x(),
                                cmsGlobalPoint.y(),
                                cmsGlobalPoint.z() );

  printVec3("cmsGlobalPoint", std::array<double,3>{
      cmsGlobalPoint.x(), cmsGlobalPoint.y(), cmsGlobalPoint.z()});
  printVec3("actsGlobalPoint", std::array<double,3>{
      actsGlobalPoint[0], actsGlobalPoint[1], actsGlobalPoint[2]});

  // ====================================================
  // 3) Geometria della surface ACTS
  // ====================================================
  Acts::GeometryContext gctx;
  const auto& trf = actsSurf.transform(gctx);

  const Acts::Vector3 actsCenter_mm = trf.translation();
  const Acts::Vector3 actsCenter = {actsCenter_mm[0] / 10.0, actsCenter_mm[1] / 10.0, actsCenter_mm[2] / 10.0}; // mm -> cm 
  const Acts::Vector3 actsUx = trf.rotation().col(0);
  const Acts::Vector3 actsUy = trf.rotation().col(1);
  const Acts::Vector3 actsUz = trf.rotation().col(2);

  printVec3("acts center", std::array<double,3>{
      actsCenter[0], actsCenter[1], actsCenter[2]});
  printVec3("acts u(global)", std::array<double,3>{
      actsUx[0], actsUx[1], actsUx[2]});
  printVec3("acts v(global)", std::array<double,3>{
      actsUy[0], actsUy[1], actsUy[2]});
  printVec3("acts n(global)", std::array<double,3>{
      actsUz[0], actsUz[1], actsUz[2]});

  const Acts::Vector3 d = actsGlobalPoint - actsCenter;

  printVec3("point-center", std::array<double,3>{d[0], d[1], d[2]});
  cout << "dist along acts n    = " << d.dot(actsUz) << "\n";
  cout << "proj along acts u    = " << d.dot(actsUx) << "\n";
  cout << "proj along acts v    = " << d.dot(actsUy) << "\n";

  // ====================================================
  // 4) Confronto diretto tra le basi locali
  //    A_ij = eActs_i dot eCms_j
  // ====================================================
  const Acts::Vector3 exCmsActs(ex_cms.x(), ex_cms.y(), ex_cms.z());
  const Acts::Vector3 eyCmsActs(ey_cms.x(), ey_cms.y(), ey_cms.z());
  const Acts::Vector3 ezCmsActs(ez_cms.x(), ez_cms.y(), ez_cms.z());

  const double u_ex = actsUx.dot(exCmsActs);
  const double u_ey = actsUx.dot(eyCmsActs);
  const double v_ex = actsUy.dot(exCmsActs);
  const double v_ey = actsUy.dot(eyCmsActs);
  const double n_ez = actsUz.dot(ezCmsActs);

  cout << "\n[DBG FRAME] basis overlaps\n";
  cout << "u·ex_cms            = " << u_ex << "\n";
  cout << "u·ey_cms            = " << u_ey << "\n";
  cout << "v·ex_cms            = " << v_ex << "\n";
  cout << "v·ey_cms            = " << v_ey << "\n";
  cout << "n·ez_cms            = " << n_ez << "\n";

  cout << "\n[DBG FRAME] expected interpretation\n";
  cout << " - u·ex ~ ±1, v·ey ~ ±1, mixed ~ 0 => same basis up to flips\n";
  cout << " - mixed terms not ~0              => axes rotated/swapped in plane\n";
  cout << " - n·ez not ~±1                    => normals not aligned\n";

  // ====================================================
  // 5) Tentativo globalToLocal con due direzioni
  // ====================================================
  const Acts::Vector3 dirActsN = actsUz;
  const Acts::Vector3 dirCmsN(ez_cms.x(), ez_cms.y(), ez_cms.z());

  printVec3("dir acts normal", std::array<double,3>{
      dirActsN[0], dirActsN[1], dirActsN[2]});
  printVec3("dir cms normal", std::array<double,3>{
      dirCmsN[0], dirCmsN[1], dirCmsN[2]});

  const auto resActsN = actsSurf.globalToLocal(gctx, actsGlobalPoint, dirActsN, 1e-2);
  const auto resCmsN  = actsSurf.globalToLocal(gctx, actsGlobalPoint, dirCmsN, 1e-2);

  cout << "\n[DBG FRAME] globalToLocal status\n";
  cout << "with acts normal ok? " << resActsN.ok() << "\n";
  if (resActsN.ok()) {
    const auto l = *resActsN;
    cout << "  local(actsN)      = (" << l[0] << ", " << l[1] << ")\n";
  }

  cout << "with cms normal ok?  " << resCmsN.ok() << "\n";
  if (resCmsN.ok()) {
    const auto l = *resCmsN;
    cout << "  local(cmsN)       = (" << l[0] << ", " << l[1] << ")\n";
  }

  // Se falliscono entrambe, fermati qui: il punto non viene localizzato su ACTS.
  if (!resActsN.ok() && !resCmsN.ok()) {
    cout << "\n[DBG FRAME] globalToLocal failed for both directions.\n";
    cout << "Likely causes:\n";
    cout << " - point is not on / close enough to the Acts surface\n";
    cout << " - wrong Acts surface passed in\n";
    cout << " - direction choice incompatible with this surface/API usage\n";
    return;
  }

  // Usa il primo risultato valido
  const Acts::Vector2 actsLoc = resActsN.ok() ? *resActsN : *resCmsN;

  // ====================================================
  // 6) Base locale ACTS via differenze finite attorno al punto trovato
  // ====================================================
  constexpr double h = 1e-3;  // mm

  Acts::Vector2 lp0 = actsLoc;
  Acts::Vector2 lp_loc0_p = actsLoc;
  Acts::Vector2 lp_loc1_p = actsLoc;
  lp_loc0_p[0] += h;
  lp_loc1_p[1] += h;

  const Acts::Vector3 g0_acts = actsSurf.localToGlobal(gctx, lp0,       dirActsN);
  const Acts::Vector3 g1_acts = actsSurf.localToGlobal(gctx, lp_loc0_p, dirActsN);
  const Acts::Vector3 g2_acts = actsSurf.localToGlobal(gctx, lp_loc1_p, dirActsN);

  Acts::Vector3 ex_acts = (g1_acts - g0_acts) / h;  // global / mm
  Acts::Vector3 ey_acts = (g2_acts - g0_acts) / h;  // global / mm

  ex_acts /= ex_acts.norm();
  ey_acts /= ey_acts.norm();

  const Acts::Vector3 ez_acts = ex_acts.cross(ey_acts).normalized();

  printVec3("ex_acts(global)", std::array<double,3>{
      ex_acts[0], ex_acts[1], ex_acts[2]});
  printVec3("ey_acts(global)", std::array<double,3>{
      ey_acts[0], ey_acts[1], ey_acts[2]});
  printVec3("ez_acts(global)", std::array<double,3>{
      ez_acts[0], ez_acts[1], ez_acts[2]});

  // ====================================================
  // 7) Matrice 2x2 di mapping:
  //    [dloc]_mm = 10 * A * [dx,dy]_cm
  //    A_ij = eActs_i dot eCms_j
  // ====================================================
  const double A00 = ex_acts.dot(exCmsActs);
  const double A01 = ex_acts.dot(eyCmsActs);
  const double A10 = ey_acts.dot(exCmsActs);
  const double A11 = ey_acts.dot(eyCmsActs);

  cout << "\n[DBG FRAME] 2x2 mapping CMSSW(x,y) -> ACTS(loc0,loc1)\n";
  cout << "A = [ [" << A00 << ", " << A01 << "],\n";
  cout << "      [" << A10 << ", " << A11 << "] ]\n";

  cout << "\n[DBG FRAME] position Jacobian block (units included)\n";
  cout << "Jpos = 10 * A = [ [" << 10. * A00 << ", " << 10. * A01 << "],\n";
  cout << "                  [" << 10. * A10 << ", " << 10. * A11 << "] ]\n";

  cout << "\n[DBG FRAME] sanity checks\n";
  cout << "ex_acts·ex_cms      = " << A00 << "\n";
  cout << "ex_acts·ey_cms      = " << A01 << "\n";
  cout << "ey_acts·ex_cms      = " << A10 << "\n";
  cout << "ey_acts·ey_cms      = " << A11 << "\n";

  cout << "\nInterpretation:\n";
  cout << " - A ~ identity         => loc0=+x, loc1=+y\n";
  cout << " - A ~ diag(-1,1)       => loc0=-x, loc1=+y\n";
  cout << " - A ~ [[0,1],[1,0]]    => axes swapped\n";
  cout << " - off-diagonal != 0    => in-plane rotation / mixed basis\n";
}

}  // namespace

// END COORDINATE DEBUG

namespace {
  // Unit conversion
  constexpr double cm_to_mm = 10.0;
  constexpr double mm_to_cm = 0.1;

  // Numerical epsilons
  constexpr double eps_denom = 1e-9;   // for dx^2+dy^2, 1-dz^2, uz^2 safeguards

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
}  // namespace 

ComputeLocalBoundJacobian::ComputeLocalBoundJacobian() {}

Acts::ActsMatrix<6,5>
ComputeLocalBoundJacobian::FromCMSSWtoACTS(const TrajectoryStateOnSurface& tsos, double eps_xy) const {
  
  // If caller passes eps_xy, keep it but don't let it be smaller than eps_denom
  const double eps = std::max(eps_xy, eps_denom);

  // Extract tx, ty from CMSSW local params.
  const auto& lp = tsos.localParameters();
  const double tx = lp.dxdz();
  const double ty = lp.dydz();

  // ===== DEBUG B: slopes consistency =====
  if(debugOn_J) {
    GlobalVector dG_debug = tsos.globalDirection().unit();
    LocalVector dL_debug  = tsos.surface().toLocal(dG_debug);

    double tx_dir = dL_debug.x() / dL_debug.z();
    double ty_dir = dL_debug.y() / dL_debug.z();

    std::cout << "\n[DBG B] SLOPES CONSISTENCY\n";
    std::cout << "tx(lp)      = " << tx << "   ty(lp)      = " << ty << "\n";
    std::cout << "tx(fromDir) = " << tx_dir << "   ty(fromDir) = " << ty_dir << "\n";
    std::cout << "dL.z        = " << dL_debug.z() << "\n";
    std::cout << "diff tx     = " << (tx - tx_dir) << "\n";
    std::cout << "diff ty     = " << (ty - ty_dir) << "\n";
  }

  // DEBUG: tx,ty dei local params sono coerenti con la globalDirection                                         
  // GlobalVector dG_debug = tsos.globalDirection().unit();
  // LocalVector  dL_debug = tsos.surface().toLocal(dG_debug);

  // double tx_dir = dL_debug.x()/dL_debug.z();
  // double ty_dir = dL_debug.y()/dL_debug.z();

  // std::cout << "\n[CHK3] tx(lp)=" << tx << " ty(lp)=" << ty
  //           << " | tx(fromDir)=" << tx_dir << " ty(fromDir)=" << ty_dir
  //           << " | dL.z=" << dL_debug.z()
  //           << "\n";
  // End DEBUG

  // Build local direction v = (tx, ty, +1) and then choose the correct sign of z
  const double s2 = 1.0 + tx*tx + ty*ty;
  const double s  = std::sqrt(s2);
  const double invs  = 1.0 / s;
  const double invs3 = 1.0 / (s2 * s);

  const auto& surf = tsos.surface();

  // --- IMPORTANT FIX ---
  // Determine the sign of local z using CMSSW's definition of the local frame:
  GlobalVector dG(tsos.globalDirection().x(),
                  tsos.globalDirection().y(),
                  tsos.globalDirection().z());
  dG = dG.unit();
  LocalVector dL = surf.toLocal(dG);

  const double signUz = (dL.z() >= 0.0) ? 1.0 : -1.0;
  // ---------------------

  // local unit direction u(tx,ty) with correct z sign
  const Acts::Vector3 v(tx, ty, 1.0);
  const Acts::Vector3 u = (signUz * invs) * v;

  // du/dtx, du/dty in local frame
  const Acts::Vector3 ex(1.0, 0.0, 0.0);
  const Acts::Vector3 ey(0.0, 1.0, 0.0);
  const Acts::Vector3 du_dtx = signUz * (invs * ex - (tx * invs3) * v);
  const Acts::Vector3 du_dty = signUz * (invs * ey - (ty * invs3) * v);

  // local->global rotation 
  const Acts::SquareMatrix3 R = localToGlobalRotation(surf);

  // global direction and derivatives
  const Acts::Vector3 d = R * u;
  const Acts::Vector3 a = R * du_dtx;
  const Acts::Vector3 b = R * du_dty;

  const double dx = d(0), dy = d(1), dz = d(2);

  // phi derivatives
  // const double denomPhi = std::max(dx*dx + dy*dy, eps);
  const double denomPhi = dx*dx + dy*dy;
  const double dphi_dtx = (dx * a(1) - dy * a(0)) / denomPhi;
  const double dphi_dty = (dx * b(1) - dy * b(0)) / denomPhi;

  // theta derivatives
  const double denomThetaSq = std::max(1.0 - dz*dz, eps);
  const double denomTheta   = std::sqrt(denomThetaSq);
  const double dtheta_dtx = -a(2) / denomTheta;
  const double dtheta_dty = -b(2) / denomTheta;

  // ===== DEBUG C: FD check for CMSSW -> ACTS angular derivatives =====
  if(debugOn_J) {
    auto wrapToPi = [](double a) {
      while (a > M_PI) a -= 2. * M_PI;
      while (a < -M_PI) a += 2. * M_PI;
      return a;
    };

    auto phiThetaFromTxTy = [&](double txp, double typ) {
      const double s2p = 1.0 + txp * txp + typ * typ;
      const double sp = std::sqrt(s2p);
      const double invsp = 1.0 / sp;

      const Acts::Vector3 vp(txp, typ, 1.0);
      const Acts::Vector3 up = (signUz * invsp) * vp;
      const Acts::Vector3 dp = R * up;

      const double phip = std::atan2(dp(1), dp(0));
      const double thetap = std::acos(std::max(-1.0, std::min(1.0, dp(2))));
      return std::pair<double,double>(phip, thetap);
    };

    const double h = 1e-6;

    auto [phi0, theta0] = phiThetaFromTxTy(tx, ty);
    auto [phi_tx, theta_tx] = phiThetaFromTxTy(tx + h, ty);
    auto [phi_ty, theta_ty] = phiThetaFromTxTy(tx, ty + h);

    const double dphi_dtx_fd   = wrapToPi(phi_tx - phi0) / h;
    const double dphi_dty_fd   = wrapToPi(phi_ty - phi0) / h;
    const double dtheta_dtx_fd = (theta_tx - theta0) / h;
    const double dtheta_dty_fd = (theta_ty - theta0) / h;

    std::cout << "\n[DBG C] FD CHECK CMSSW->ACTS\n";
    std::cout << "dphi/dtx   ana=" << dphi_dtx   << "   fd=" << dphi_dtx_fd
              << "   diff=" << (dphi_dtx - dphi_dtx_fd) << "\n";
    std::cout << "dphi/dty   ana=" << dphi_dty   << "   fd=" << dphi_dty_fd
              << "   diff=" << (dphi_dty - dphi_dty_fd) << "\n";
    std::cout << "dtheta/dtx ana=" << dtheta_dtx << "   fd=" << dtheta_dtx_fd
              << "   diff=" << (dtheta_dtx - dtheta_dtx_fd) << "\n";
    std::cout << "dtheta/dty ana=" << dtheta_dty << "   fd=" << dtheta_dty_fd
              << "   diff=" << (dtheta_dty - dtheta_dty_fd) << "\n";
  }


  // DEBUG
  // auto wrapToPi = [](double a) {
  //   while (a > M_PI) a -= 2*M_PI;
  //   while (a < -M_PI) a += 2*M_PI;
  //   return a;
  // };

  // auto phiThetaFromTxTy = [&](double txp, double typ) {
  //   // rebuild u from slopes with SAME signUz logic
  //   const double s2p = 1.0 + txp*txp + typ*typ;
  //   const double sp  = std::sqrt(s2p);
  //   const double invsp = 1.0 / sp;

  //   const Acts::Vector3 vp(txp, typ, 1.0);
  //   const Acts::Vector3 up = (signUz * invsp) * vp;

  //   const Acts::Vector3 dp = R * up;
  //   const double dxp = dp(0), dyp = dp(1), dzp = dp(2);

  //   const double phip   = std::atan2(dyp, dxp);
  //   const double thetap = std::acos(std::max(-1.0, std::min(1.0, dzp)));
  //   return std::pair<double,double>(phip, thetap);
  // };

  // const double h = 1e-6; // step for FD

  // auto [phi0, theta0] = phiThetaFromTxTy(tx, ty);

  // // tx perturb
  // auto [phi_tx, theta_tx] = phiThetaFromTxTy(tx + h, ty);
  // const double dphi_dtx_fd   = wrapToPi(phi_tx - phi0) / h;
  // const double dtheta_dtx_fd = (theta_tx - theta0) / h;

  // // ty perturb
  // auto [phi_ty, theta_ty] = phiThetaFromTxTy(tx, ty + h);
  // const double dphi_dty_fd   = wrapToPi(phi_ty - phi0) / h;
  // const double dtheta_dty_fd = (theta_ty - theta0) / h;

  // std::cout << "\n[FD] phi0=" << phi0 << " theta0=" << theta0
  //           << " signUz=" << signUz
  //           << "\n[FD] dphi/dtx   ana=" << dphi_dtx   << " fd=" << dphi_dtx_fd
  //           << "\n[FD] dphi/dty   ana=" << dphi_dty   << " fd=" << dphi_dty_fd
  //           << "\n[FD] dtheta/dtx ana=" << dtheta_dtx << " fd=" << dtheta_dtx_fd
  //           << "\n[FD] dtheta/dty ana=" << dtheta_dty << " fd=" << dtheta_dty_fd
  //           << "\n";
  // -------- END FD CHECK --------
  // End DEBUG

  
  // debugCompareLocalFrames(tsos.surface(), actsSurface, tsos.globalPosition());
  

  // indices
  constexpr int rLoc0  = 0, rLoc1  = 1, rPhi = 2, rTheta = 3, rQop = 4, rTime = 5;
  constexpr int cQop   = 0, cTx    = 1, cTy  = 2, cX     = 3, cY    = 4;

  Acts::ActsMatrix<6,5> J = Acts::ActsMatrix<6,5>::Zero();

  // UNIT FIX: ACTS loc0/loc1 in mm, CMSSW x/y in cm
  J(rLoc0, cX) = cm_to_mm;  // 10
  J(rLoc1, cY) = cm_to_mm;  // 10

  // angles from slopes
  J(rPhi,   cTx) = dphi_dtx;
  J(rPhi,   cTy) = dphi_dty;
  J(rTheta, cTx) = dtheta_dtx;
  J(rTheta, cTy) = dtheta_dty;

  // q/p direct
  J(rQop, cQop) = 1.0;

  (void)rTime;
  return J;
}

Acts::ActsMatrix<5,6>
ComputeLocalBoundJacobian::FromACTStoCMSSW(const Surface& surf,
                                          double phi, double theta) const {
  // Rotation matrices
  const Acts::SquareMatrix3 R  = localToGlobalRotation(surf);
  const Acts::SquareMatrix3 Rt = R.transpose();

  // global direction d(phi,theta)
  const double cphi = std::cos(phi), sphi = std::sin(phi);
  const double cth  = std::cos(theta), sth = std::sin(theta);

  const Acts::Vector3 d(cphi*sth, sphi*sth, cth);
  const Acts::Vector3 dd_dphi(-sphi*sth, cphi*sth, 0.0);
  const Acts::Vector3 dd_dtheta(cphi*cth, sphi*cth, -sth);

  // local direction u and derivatives
  const Acts::Vector3 u = Rt * d;
  const Acts::Vector3 du_dphi   = Rt * dd_dphi;
  const Acts::Vector3 du_dtheta = Rt * dd_dtheta;

  // ===== DEBUG A: frame consistency =====
  if(debugOn_J) {
    GlobalVector dG(d(0), d(1), d(2));
    LocalVector u_cms = surf.toLocal(dG);

    GlobalVector ddphiG(dd_dphi(0), dd_dphi(1), dd_dphi(2));
    GlobalVector ddthetaG(dd_dtheta(0), dd_dtheta(1), dd_dtheta(2));

    LocalVector duphi_cms   = surf.toLocal(ddphiG);
    LocalVector dutheta_cms = surf.toLocal(ddthetaG);

    std::cout << "\n[DBG A] FRAME CONSISTENCY IN FromACTStoCMSSW\n";
    std::cout << "u(R^T d)        = " << u.transpose() << "\n";
    std::cout << "u(surf.toLocal) = "
              << u_cms.x() << " " << u_cms.y() << " " << u_cms.z() << "\n";

    std::cout << "du/dphi   (R^T) = " << du_dphi.transpose() << "\n";
    std::cout << "du/dphi   (CMS) = "
              << duphi_cms.x() << " " << duphi_cms.y() << " " << duphi_cms.z() << "\n";

    std::cout << "du/dtheta (R^T) = " << du_dtheta.transpose() << "\n";
    std::cout << "du/dtheta (CMS) = "
              << dutheta_cms.x() << " " << dutheta_cms.y() << " " << dutheta_cms.z() << "\n";
  }


  // DEBUG: Am I using the right rotation matrix?
  // std::cout << "u from ACTS: " << u.transpose() << "\n"
  //           << "u from CMSSW: " << surf.toLocal(GlobalVector(d.x(),d.y(),d.z())) << std::endl;
  // const auto u_phi_cmssw   = surf.toLocal(GlobalVector(dd_dphi.x(), dd_dphi.y(), dd_dphi.z()));
  // const auto u_theta_cmssw = surf.toLocal(GlobalVector(dd_dtheta.x(), dd_dtheta.y(), dd_dtheta.z()));
  // std::cout << "du/dphi  ACTS: " << du_dphi.transpose()   << "\n"
  //         << "du/dphi  CMSSW: " << u_phi_cmssw << "\n"
  //         << "du/dtheta ACTS: " << du_dtheta.transpose() << "\n"
  //         << "du/dtheta CMSSW: " << u_theta_cmssw << "\n";
  // End DEBUG

  const double ux = u(0), uy = u(1), uz = u(2);

  // protect uz^2
  // const double uz2 = std::max(uz*uz, eps_denom);
  const double uz2 = uz*uz;
  
  // const double sgnUz = (uz >= 0.) ? 1.0 : -1.0;
  // auto dtx_dalpha = [&](const Acts::Vector3& du_dalpha) {
  //   return sgnUz * (uz * du_dalpha(0) - ux * du_dalpha(2)) / uz2;
  // };
  // auto dty_dalpha = [&](const Acts::Vector3& du_dalpha) {
  //   return sgnUz * (uz * du_dalpha(1) - uy * du_dalpha(2)) / uz2;
  // };

  auto dtx_dalpha = [&](const Acts::Vector3& du_dalpha) {
    return (uz * du_dalpha(0) - ux * du_dalpha(2)) / uz2;
  };
  auto dty_dalpha = [&](const Acts::Vector3& du_dalpha) {
    return (uz * du_dalpha(1) - uy * du_dalpha(2)) / uz2;
  };

  const double dtx_dphi   = dtx_dalpha(du_dphi);
  const double dtx_dtheta = dtx_dalpha(du_dtheta);
  const double dty_dphi   = dty_dalpha(du_dphi);
  const double dty_dtheta = dty_dalpha(du_dtheta);

  // ===== DEBUG D: FD check for ACTS -> CMSSW derivatives =====
  if(debugOn_J) {
    auto txTyFromPhiTheta = [&](double phip, double thetap) {
      const double cphip = std::cos(phip), sphip = std::sin(phip);
      const double cthp  = std::cos(thetap), sthp = std::sin(thetap);

      Acts::Vector3 dp(cphip * sthp, sphip * sthp, cthp);
      Acts::Vector3 up = Rt * dp;

      double uzp = up(2);
      double uzp_safe = (std::abs(uzp) < 1e-12) ? std::copysign(1e-12, uzp) : uzp;

      return std::pair<double,double>(up(0)/uzp_safe, up(1)/uzp_safe);
    };

    const double h = 1e-6;

    auto [tx0, ty0] = txTyFromPhiTheta(phi, theta);
    auto [tx_phi, ty_phi] = txTyFromPhiTheta(phi + h, theta);
    auto [tx_th,  ty_th ] = txTyFromPhiTheta(phi, theta + h);

    const double dtx_dphi_fd   = (tx_phi - tx0) / h;
    const double dty_dphi_fd   = (ty_phi - ty0) / h;
    const double dtx_dtheta_fd = (tx_th  - tx0) / h;
    const double dty_dtheta_fd = (ty_th  - ty0) / h;

    std::cout << "\n[DBG D] FD CHECK ACTS->CMSSW\n";
    std::cout << "dtx/dphi   ana=" << dtx_dphi   << "   fd=" << dtx_dphi_fd
              << "   diff=" << (dtx_dphi - dtx_dphi_fd) << "\n";
    std::cout << "dtx/dtheta ana=" << dtx_dtheta << "   fd=" << dtx_dtheta_fd
              << "   diff=" << (dtx_dtheta - dtx_dtheta_fd) << "\n";
    std::cout << "dty/dphi   ana=" << dty_dphi   << "   fd=" << dty_dphi_fd
              << "   diff=" << (dty_dphi - dty_dphi_fd) << "\n";
    std::cout << "dty/dtheta ana=" << dty_dtheta << "   fd=" << dty_dtheta_fd
              << "   diff=" << (dty_dtheta - dty_dtheta_fd) << "\n";
  }

  // indices
  constexpr int rQop = 0, rTx = 1, rTy = 2, rX = 3, rY = 4;
  constexpr int cLoc0 = 0, cLoc1 = 1, cPhi = 2, cTheta = 3, cQop = 4, cTime = 5;

  Acts::ActsMatrix<5,6> K = Acts::ActsMatrix<5,6>::Zero();

  // q/p
  K(rQop, cQop) = 1.0;

  // UNIT FIX: CMSSW x/y in cm, ACTS loc0/loc1 in mm
  K(rX, cLoc0) = mm_to_cm;  // 0.1
  K(rY, cLoc1) = mm_to_cm;  // 0.1

  // slopes from angles
  K(rTx, cPhi)   = dtx_dphi;
  K(rTx, cTheta) = dtx_dtheta;
  K(rTy, cPhi)   = dty_dphi;
  K(rTy, cTheta) = dty_dtheta;

  (void)cTime;
  return K;
}


