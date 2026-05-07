#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Any.hpp"
#include "Acts/Utilities/Result.hpp"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "MagneticField/VolumeBasedEngine/interface/VolumeBasedMagneticField.h"
#include "MagneticField/VolumeGeometry/interface/MagVolume.h"

class CMSMagneticFieldProvider : public Acts::MagneticFieldProvider {
public:
  struct CMSCache {
    const MagVolume* volume = nullptr;
  };

  CMSMagneticFieldProvider(const MagneticField& bField, const Acts::Logging::Level& lvl);

  virtual Acts::MagneticFieldProvider::Cache makeCache(const Acts::MagneticFieldContext& mctx) const override;

  virtual Acts::Result<Acts::Vector3> getField(const Acts::Vector3& position, Cache& cache) const override;

  ~CMSMagneticFieldProvider() = default;

private:
  const MagneticField& bField_;
  const VolumeBasedMagneticField* vbField_;

  Acts::Logging::Level m_level;
  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }
    
};