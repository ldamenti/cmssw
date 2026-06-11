#ifndef Geometry_TrackerGeometryBuilder_ActsGeoBuilderUtils_H
#define Geometry_TrackerGeometryBuilder_ActsGeoBuilderUtils_H

#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Detector/KdtSurfacesProvider.hpp"

#include <memory>
#include <string>
#include <vector>

using KdtSurfacesDim2Bin100 = Acts::Experimental::KdtSurfaces<2u, 100u>;

class ActsGeoBuilderUtils {
public:
    ActsGeoBuilderUtils() = default;
    ~ActsGeoBuilderUtils() = default;

    std::vector<std::shared_ptr<Acts::Surface>> SelectActiveSurfaces_PhaseII(const KdtSurfacesDim2Bin100& surfaces, const std::string Layer_name);
 
private:

};

#endif
