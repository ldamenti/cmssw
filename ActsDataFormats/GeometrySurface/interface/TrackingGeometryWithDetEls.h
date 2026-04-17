#pragma once
#include <memory>
#include "ActsDataFormats/GeometrySurface/interface/CMSDetectorElement.h"
#include "Acts/Geometry/TrackingGeometry.hpp"

using DetElVect = std::vector<std::shared_ptr<Acts::CMSDetectorElement>>;
struct TrackingGeometryWithDetEls {
  DetElVect detElements;
  std::shared_ptr<Acts::TrackingGeometry> trackingGeometry;
};
