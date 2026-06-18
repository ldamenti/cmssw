#ifndef _COMPUTEFREEJACOBIAN_H_
#define _COMPUTEFREEJACOBIAN_H_

#include <cmath>
#include <algorithm>

// ACTS
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"

// CMSSW
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

class ComputeFreeJacobian {

  public:

    ComputeFreeJacobian();

    Eigen::Matrix<double, 8, 6> FromCMSSWtoACTS(const FreeTrajectoryState& fts) const;
    Eigen::Matrix<double, 6, 8> FromACTStoCMSSW(const Acts::Vector3& direction, double qOverP, double charge) const;
    
};


#endif