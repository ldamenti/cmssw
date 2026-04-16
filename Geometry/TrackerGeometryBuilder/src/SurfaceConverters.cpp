#include "Geometry/TrackerGeometryBuilder/interface/SurfaceConverters.hpp"

SurfaceConverters::SurfaceConverters(const TrackerGeometry* trkGeo,
                                     const TrackingGeometryWithDetEls& actsGeo,
                                     const Acts::Logging::Level& lvl)
  :m_Level(lvl), m_logger(Acts::getDefaultLogger("SurfaceConverters", lvl)),
   m_trackerGeometry(trkGeo), m_TrkgGeoActs(actsGeo) {

    for (auto const* det : m_trackerGeometry->dets()) {
      cmsSurfaceToDetId_[&det->surface()] = det->geographicalId().rawId();
    }

    for (auto const& detEl : m_TrkgGeoActs.detElements) {
      detIdToDetEl_[detEl->detID()] = detEl;   // keep ownership here
    }
}

SurfaceConverters::~SurfaceConverters(){}

std::shared_ptr<const Acts::Surface> SurfaceConverters::fromCMSSWtoACTS(const Surface& surf_cmssw) const {

  uint32_t detId = 0;

  auto itFast = cmsSurfaceToDetId_.find(&surf_cmssw);
  if (itFast != cmsSurfaceToDetId_.end()) {
    // Fast pointer comparison succeeded
    ACTS_VERBOSE("Successfully matched CMSSW surface to DetId using fast detID comparison.");
    detId = itFast->second;
  } else {
    ACTS_VERBOSE("Failed to match CMSSW surface to any TrackerGeometry surface.");
    auto& p = surf_cmssw.position();
    ACTS_VERBOSE("Surface position: R = "
                  << std::sqrt(p.x()*p.x() + p.y()*p.y())
                  << ", Z = " << p.z() << " cm");
    ACTS_VERBOSE("Bounds type: " << typeid(surf_cmssw.bounds()).name());
    return {};
  }

  auto it2 = detIdToDetEl_.find(detId);
  if (it2 == detIdToDetEl_.end()) {
    ACTS_VERBOSE("Failed to find the corresponding ACTS surface associated with CMSSW DetId; cannot convert to ACTS surface.");
    return {};
  }

  // Print CMSSW ans ACTS surface properties if needed:
  auto& pos = surf_cmssw.position();
  auto suft_acts = &it2->second->surface();
  auto subdet = it2->second->subDetector();
  auto t_final = suft_acts->transform(Acts::GeometryContext{});
  auto b_values = suft_acts->bounds().values();
  ACTS_VERBOSE("Converting CMSSW surface with bounds: " << typeid(surf_cmssw.bounds()).name());
  ACTS_VERBOSE("Surface Position in cmssw [mm]: " << pos.x()*10 << ", " <<  pos.y()*10 << ", " << pos.z()*10);
  ACTS_VERBOSE("Result of conversion. Surface with type: " << suft_acts->type() << "; Subdetector: " << subdet);
  ACTS_VERBOSE("Surface Position in acts [mm]: " << t_final.translation().transpose());
  ACTS_VERBOSE("Bounds: ");
  for (size_t i = 0; i < b_values.size(); ++i) {
      ACTS_VERBOSE(" " << b_values[i]);
  }
  ACTS_VERBOSE("Consistency check (localToGlobal of the local point 0,0,0):");
  ACTS_VERBOSE("ACTS: " << suft_acts->localToGlobal(Acts::GeometryContext{}, Acts::Vector2{.0,.0}, Acts::Vector3{0,0,0}).transpose());
  ACTS_VERBOSE("CMSSW: " << surf_cmssw.toGlobal(LocalPoint(.0,.0)));
       

  auto& surf_acts = it2->second->surface();               // Surface&

  GlobalPoint cmsPos = surf_cmssw.position();
  GlobalVector cmsN = surf_cmssw.toGlobal(LocalVector(0.f, 0.f, 1.f));
  Acts::Vector3 cmsNormal(cmsN.x(), cmsN.y(), cmsN.z());
  cmsNormal.normalize();
  Acts::Vector3 actsPos(cmsPos.x() * 10., cmsPos.y() * 10., cmsPos.z() * 10.);
  Acts::Vector3 actsNormal = surf_acts.normal(Acts::GeometryContext{}, actsPos, Acts::Vector3{0,0,1});
  actsNormal.normalize();
  double cosAngle = cmsNormal.dot(actsNormal);
  ACTS_VERBOSE("CMSSW normal: " << cmsNormal.transpose());
  ACTS_VERBOSE("ACTS  normal: " << actsNormal.transpose());
  ACTS_VERBOSE("cos(angle)   : " << cosAngle);

  {
  GlobalPoint g0  = surf_cmssw.toGlobal(LocalPoint(0., 0.));
  GlobalPoint gx1 = surf_cmssw.toGlobal(LocalPoint(1., 0.));
  GlobalPoint gy1 = surf_cmssw.toGlobal(LocalPoint(0., 1.));
  GlobalPoint gz1 = surf_cmssw.toGlobal(LocalPoint(0., 0., 1.));

  Acts::Vector3 exCms(gx1.x()-g0.x(), gx1.y()-g0.y(), gx1.z()-g0.z());
  Acts::Vector3 eyCms(gy1.x()-g0.x(), gy1.y()-g0.y(), gy1.z()-g0.z());
  Acts::Vector3 ezCms(gz1.x()-g0.x(), gz1.y()-g0.y(), gz1.z()-g0.z());
  exCms.normalize();
  eyCms.normalize();
  ezCms.normalize();

  const auto T = surf_acts.transform(Acts::GeometryContext{});
  Acts::Vector3 exA = T.linear().col(0).normalized();
  Acts::Vector3 eyA = T.linear().col(1).normalized();
  Acts::Vector3 ezA = T.linear().col(2).normalized();

  ACTS_VERBOSE("[DBG FRAME FULL] exCms·exA = " << exCms.dot(exA));
  ACTS_VERBOSE("[DBG FRAME FULL] exCms·eyA = " << exCms.dot(eyA));
  ACTS_VERBOSE("[DBG FRAME FULL] exCms·ezA = " << exCms.dot(ezA));
  ACTS_VERBOSE("[DBG FRAME FULL] eyCms·exA = " << eyCms.dot(exA));
  ACTS_VERBOSE("[DBG FRAME FULL] eyCms·eyA = " << eyCms.dot(eyA));
  ACTS_VERBOSE("[DBG FRAME FULL] eyCms·ezA = " << eyCms.dot(ezA));
  ACTS_VERBOSE("[DBG FRAME FULL] ezCms·exA = " << ezCms.dot(exA));
  ACTS_VERBOSE("[DBG FRAME FULL] ezCms·eyA = " << ezCms.dot(eyA));
  ACTS_VERBOSE("[DBG FRAME FULL] ezCms·ezA = " << ezCms.dot(ezA));
}



  if (std::abs(std::abs(cosAngle) - 1.) > 1e-3) {
    ACTS_WARNING("Surface normals not parallel! Possible transform mismatch.");
  } else if (cosAngle < 0.) {
    ACTS_WARNING("Surface normals are antiparallel (flip)! Same plane, opposite normal.");
  }

  return surf_acts.getSharedPtr();   
}

std::shared_ptr<Surface> SurfaceConverters::fromACTStoCMSSW(const Acts::Surface& surf_acts) {
    
  // Get the detector elemeent associated to the surface to obtain the detId
  auto detEl = dynamic_cast<const Acts::CMSDetectorElement*>(surf_acts.associatedDetectorElement());
  auto detID = detEl->detID();

  // Get the GeoDet unit from the detId
  const GeomDet* geomDet = m_trackerGeometry->idToDet(detID);    

  return std::make_shared<Plane>(geomDet->specificSurface()); 
}
