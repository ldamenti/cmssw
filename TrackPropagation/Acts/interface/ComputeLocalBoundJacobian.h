#ifndef _COMPUTELOCALBOUNDJACOBIAN_H_
#define _COMPUTELOCALBOUNDJACOBIAN_H_

#include <cmath>
#include <algorithm>

// ACTS
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

// CMSSW
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"

class ComputeLocalBoundJacobian {

  public:

    ComputeLocalBoundJacobian();

    Acts::ActsMatrix<6,5> FromCMSSWtoACTS(const TrajectoryStateOnSurface& tsos, double eps_xy = 1e-9) const;
    Acts::ActsMatrix<5,6> FromACTStoCMSSW(const Surface& surf, double phi, double theta) const;

};


#endif