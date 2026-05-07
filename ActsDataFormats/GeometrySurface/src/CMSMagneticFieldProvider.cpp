#include "ActsDataFormats/GeometrySurface/interface/CMSMagneticFieldProvider.hpp"

CMSMagneticFieldProvider::CMSMagneticFieldProvider(const MagneticField& bField, const Acts::Logging::Level& lvl): 
                            bField_(bField),
                            vbField_(dynamic_cast<const VolumeBasedMagneticField*>(&bField)),
                            m_level(lvl),
                            m_logger(Acts::getDefaultLogger("CMSMagneticFieldProvider", lvl)) {
    if (!vbField_) {
      ACTS_VERBOSE("MagneticField is not a VolumeBasedMagneticField; using MagneticField::inTesla fallback");
    } else {
      ACTS_VERBOSE("Successfully casted MagneticField into VolumeBasedMagneticField");
    }
};

Acts::MagneticFieldProvider::Cache CMSMagneticFieldProvider::makeCache(const Acts::MagneticFieldContext&) const {
  return Acts::MagneticFieldProvider::Cache(CMSCache{});
}

Acts::Result<Acts::Vector3> CMSMagneticFieldProvider::getField(const Acts::Vector3& position, Cache& cache) const {
  
  // Convert the acts global point into an acts one
  GlobalPoint gPoint(position[0] / 10., position[1] / 10., position[2] / 10.); // mm -> cm

  // Fallback in case the B filed is not volume based:
  if (!vbField_) {
    GlobalVector field = bField_.inTesla(gPoint) * Acts::UnitConstants::T;
    ACTS_VERBOSE("B field from generic MagneticField in " << gPoint << " [cm] is " << field / Acts::UnitConstants::T << " [T]");
    return Acts::Result<Acts::Vector3>::success(Acts::Vector3(field.x(), field.y(), field.z()));
  }

  // If it's a volume base, first take the volume stored in the cache
  auto& cmsCache = cache.as<CMSCache>();
  const MagVolume* mVol = cmsCache.volume;

  // If we have no volume in the cache or the gPoint is not 
  // anymore in that volume, we recompute a new magnetic volume
  if (mVol == nullptr || !mVol->inside(gPoint)) {
    ACTS_VERBOSE("Volume stored in the cache is nullptr or point outside the cached volume.");
    ACTS_VERBOSE("Computing a new Magnetic Volume");
    mVol = vbField_->findVolume(gPoint);
    cmsCache.volume = mVol;
  } else {
    ACTS_VERBOSE("Global point is inside cached volumes, continuing with that");
  }

  if (mVol != nullptr) {
    // Get the B field
    GlobalVector B = mVol->fieldInTesla(gPoint) * Acts::UnitConstants::T;
    ACTS_VERBOSE("B field in " << gPoint << " [cm] is " << B / Acts::UnitConstants::T << " [T]");
    return Acts::Result<Acts::Vector3>::success(Acts::Vector3(B.x(), B.y(), B.z()));
  } 
  
  ACTS_VERBOSE("No MagVolume found for point " << gPoint << " [cm]");
  GlobalVector field = bField_.inTesla(gPoint) * Acts::UnitConstants::T;
  ACTS_VERBOSE("Proceeding with value stored in MagneticField: " << field / Acts::UnitConstants::T << " [T]");
  return Acts::Result<Acts::Vector3>::success(Acts::Vector3(field.x(), field.y(), field.z()));

}


