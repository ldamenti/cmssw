#include "TrackPropagation/Acts/interface/ComputeFreeJacobian.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

ComputeFreeJacobian::ComputeFreeJacobian() {}

Eigen::Matrix<double, 8, 6> 
ComputeFreeJacobian::FromCMSSWtoACTS(const FreeTrajectoryState& fts) const {
  const auto& cmsCov = fts.cartesianError().matrix();
  const auto& pCMS = fts.momentum();

  Acts::Vector3 p{pCMS.x(), pCMS.y(), pCMS.z()};

  const double q = fts.charge();
  const double pMag = p.norm();

  Acts::Vector3 dir = p / pMag;

  // Jacobiana: ACTS free 8 params vs CMSSW cartesian 6 params
  Eigen::Matrix<double, 8, 6> J;
  J.setZero();

  // x,y,z: CMSSW cm -> ACTS mm
  J(Acts::eFreePos0, 0) = 10;
  J(Acts::eFreePos1, 1) = 10;
  J(Acts::eFreePos2, 2) = 10;

  // direction = p / |p|
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d dDir_dP = (I - dir * dir.transpose()) / pMag;

  J.block<3,3>(Acts::eFreeDir0, 3) = dDir_dP;

  // q/p
  // d(q/p)/dP_i = -q * p_i / |p|^3
  Eigen::RowVector3d dQoP_dP = -q * p.transpose() / std::pow(pMag, 3);

  J.block<1,3>(Acts::eFreeQOverP, 3) = dQoP_dP;

  return J;
}

Eigen::Matrix<double, 6, 8>
ComputeFreeJacobian::FromACTStoCMSSW(const Acts::Vector3& direction, double qOverP, double charge) const {

  Eigen::Matrix<double, 6, 8> K;
  K.setZero();

  // ACTS position mm -> CMSSW position cm
  K(0, Acts::eFreePos0) = 0.1;
  K(1, Acts::eFreePos1) = 0.1;
  K(2, Acts::eFreePos2) = 0.1;

  const double absQOverP = std::abs(qOverP);
  const double pMag = 1.0 / absQOverP;

  // p_vec = |p| * direction
  //
  // d p_vec / d direction = |p| I
  K.block<3, 3>(3, Acts::eFreeDir0) = pMag * Eigen::Matrix3d::Identity();

  // |p| = 1 / |q/p|
  //
  // d|p|/d(q/p) = -sign(q/p) / (q/p)^2
  //
  // p_vec = |p| * direction
  const double dPmag_dQOverP = -std::copysign(1.0, qOverP) / (qOverP * qOverP);

  K.block<3, 1>(3, Acts::eFreeQOverP) = direction * dPmag_dQOverP;


  return K;
}


