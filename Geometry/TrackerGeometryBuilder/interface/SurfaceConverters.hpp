#ifndef SURFACE_CONVERTERS_HPP
#define SURFACE_CONVERTERS_HPP

// ACTS:
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Utilities/Logger.hpp"
// CMSSW:
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/SimpleCylinderBounds.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
// Custom:
#include "ActsDataFormats/GeometrySurface/interface/CMSDetectorElement.h"
#include "ActsDataFormats/GeometrySurface/interface/TrackingGeometryWithDetEls.h"
#include <unordered_map>
#include <typeinfo>
#include <iostream>

// DEBUG
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"

using Key3D = std::tuple<double, double, double>;
using DetIdRaw = uint32_t;
using CmsSurfaceToDetId = std::unordered_map<const Surface*, DetIdRaw>;
using DetIdToDetEl = std::unordered_map<DetIdRaw, std::shared_ptr<const Acts::CMSDetectorElement>>;

struct CmsSurfaceEntry {
  uint32_t detId;
  const GeomDet* geomDet;
  const Surface* surface;
  GlobalPoint pos;
  GlobalVector normal;
  std::string boundsName;
  bool isLeaf;
};

class SurfaceConverters {

    public:
    SurfaceConverters(const TrackerGeometry* trkGeo,
                    const TrackingGeometryWithDetEls& actsGeo,
                    const Acts::Logging::Level& lvl);
    ~SurfaceConverters();

    std::shared_ptr<const Acts::Surface> fromCMSSWtoACTS(const Surface& surf_cmssw) const;

    std::shared_ptr<Surface> fromACTStoCMSSW(const Acts::Surface& surf_acts);

    private:
    Acts::Logging::Level m_Level;
    std::unique_ptr<const Acts::Logger> m_logger;
    const TrackerGeometry* m_trackerGeometry;
    TrackingGeometryWithDetEls m_TrkgGeoActs;

    std::vector<CmsSurfaceEntry> cmsSurfaces_;

    const Acts::Logger& logger() const { return *m_logger; }

    Acts::ObjVisualization3D Surf_obj;
    CmsSurfaceToDetId cmsSurfaceToDetId_;
    DetIdToDetEl detIdToDetEl_;

};
 


#endif