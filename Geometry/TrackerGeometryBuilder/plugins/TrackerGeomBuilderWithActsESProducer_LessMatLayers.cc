#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Utilities/interface/typelookup.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/ACTSTrackerGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "CondFormats/AlignmentRecord/interface/TrackerAlignmentRcd.h"
#include "CondFormats/Alignment/interface/Alignments.h"
#include "CondFormats/Alignment/interface/AlignTransform.h"
#include "Alignment/CommonAlignment/interface/Alignable.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"

#include "DataFormats/TrackerCommon/interface/PixelBarrelName.h"
#include "DataFormats/TrackerCommon/interface/PixelEndcapName.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DetectorDescription/Core/interface/DDRotationMatrix.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "Math/RotationZ.h"
#include "FWCore/Framework/interface/ESConsumesCollector.h"

#include "Math/Rotation3D.h"
#include "Math/AxisAngle.h"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/KdtSurfacesProvider.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometryVisitor.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/ObjVisualization3D.hpp"
#include "ActsPlugins/Json/JsonDetectorElement.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp" 
#include "Acts/Geometry/LayerBlueprintNode.hpp" 
#include "Acts/Geometry/MaterialDesignatorBlueprintNode.hpp" 
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Detector/KdtSurfacesProvider.hpp"
#include "Acts/Detector/ProtoDetector.hpp"
#include "Acts/Detector/detail/ReferenceGenerators.hpp"
#include "Acts/Detector/interface/ISurfacesProvider.hpp"
#include "Acts/Detector/LayerStructureBuilder.hpp"
#include "Acts/Detector/detail/BlueprintHelper.hpp"
#include "Acts/Detector/detail/BlueprintDrawer.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "ActsPlugins/ActSVG/DetectorSvgConverter.hpp"
#include "ActsPlugins/Json/JsonSurfacesReader.hpp"
#include "ActsPlugins/Json/JsonMaterialDecorator.hpp"
#include "ActsPlugins/Json/SurfaceJsonConverter.hpp"
#include "ActsPlugins/ActSVG/TrackingGeometrySvgConverter.hpp"
#include "ActsPlugins/ActSVG/SurfaceArraySvgConverter.hpp"
#include "ActsPlugins/Root/RootMaterialTrackIo.hpp"

#include "Acts/Material/PropagatorMaterialAssigner.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulater.hpp"
#include "Acts/Material/BinnedSurfaceMaterial.hpp"
#include "Acts/Material/AccumulatedMaterialSlab.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/MaterialValidater.hpp"
#include "Acts/Material/MaterialMapper.hpp"

#include <fstream>
#include <iomanip>
#include <random>
#include <nlohmann/json.hpp> 
#include "TChain.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <vector>

#include "ActsDataFormats/GeometrySurface/interface/CMSDetectorElement.h"
#include "ActsDataFormats/GeometrySurface/interface/JsonMaterialWriter.hpp"
#include "ActsDataFormats/GeometrySurface/interface/TrackingGeometryWithDetEls.h"

using json = nlohmann::json;
using KdtSurfacesDim2Bin100 = Acts::Experimental::KdtSurfaces<2u, 100u>;
const std::array<Acts::AxisDirection, 2UL> casts{Acts::AxisDirection::AxisZ, Acts::AxisDirection::AxisR};

using DetElVect = std::vector<std::shared_ptr<Acts::CMSDetectorElement>>;

namespace {
  // ===== Helper function to visualize the sufaces =====
  void writeSurfacesObj(const std::vector<std::shared_ptr<Acts::Surface>>& surfaces, 
                        const std::string& fileName)
  {
      Acts::ObjVisualization3D obj;

      for (size_t i = 0; i < surfaces.size(); ++i) {
          if (!surfaces[i]) {
              std::cerr << "[ERROR] Surface " << i << " is nullptr!" << std::endl;
              continue;
          }

          Acts::GeometryView3D::drawSurface(
              obj,
              *surfaces[i],
              Acts::GeometryContext{},
              Acts::Transform3::Identity(),
              Acts::ViewConfig{}
          );
      }

      obj.write(fileName);
      obj.clear();
  };

  Acts::Transform3 GenerateTranslation(double dx, double dy, double dz) {
      return Acts::Transform3::Identity() * Acts::Translation3{Acts::Vector3{dx, dy, dz}};
  };

  std::vector<std::shared_ptr<Acts::Surface>> SelectActiveSurfaces_PhaseI(const KdtSurfacesDim2Bin100& surfaces, std::string Layer_name){

    Acts::RangeXD<2, double, std::array> Range({0, 0}, {0, 0});

    if(Layer_name=="PixL0"){
      Range = Acts::RangeXD<2, double, std::array>({-290, 15}, {280, 35});
    }
    else if(Layer_name=="PixL1"){
      Range = Acts::RangeXD<2, double, std::array>({-290,40}, {280,80});
    }
    else if(Layer_name=="PixL2"){
      Range = Acts::RangeXD<2, double, std::array>({-290,85}, {280,150});
    }
    else if(Layer_name=="PixL3"){
      Range = Acts::RangeXD<2, double, std::array>({-290,150}, {280,200});
    }
    else if(Layer_name=="PixelPos2"){
      Range = Acts::RangeXD<2, double, std::array>({450,25}, {1160,220}); //550
    }
    else if(Layer_name=="PixelPos1"){
      Range = Acts::RangeXD<2, double, std::array>({350,25}, {445,220});
    }
    else if(Layer_name=="PixelPos0"){
      Range = Acts::RangeXD<2, double, std::array>({285,25}, {345,220});
    }
    else if(Layer_name=="PixelNeg2"){
      Range = Acts::RangeXD<2, double, std::array>({-1160,25}, {-450,220});
    }
    else if(Layer_name=="PixelNeg1"){
      Range = Acts::RangeXD<2, double, std::array>({-445,25}, {-350,220});
    }
    else if(Layer_name=="PixelNeg0"){
      Range = Acts::RangeXD<2, double, std::array>({-345,25}, {-285,220});
    }
    // Default TID:
    else if(Layer_name=="TIDNeg2"){
      Range = Acts::RangeXD<2, double, std::array>({-1160,220}, {-1000,530});
    }
    else if(Layer_name=="TIDNeg1"){
      Range = Acts::RangeXD<2, double, std::array>({-1000,220}, {-850,530});
    }
    else if(Layer_name=="TIDNeg0"){
      Range = Acts::RangeXD<2, double, std::array>({-850,220}, {-700,530});
    }
    // Extra TID:
    else if(Layer_name=="TIDNeg2_2"){
      Range = Acts::RangeXD<2, double, std::array>({-1130,220}, {-1050,520});
    }
    else if(Layer_name=="TIDNeg2_1"){
      Range = Acts::RangeXD<2, double, std::array>({-1050,220}, {-980,520});
    }
    else if(Layer_name=="TIDNeg1_2"){
      Range = Acts::RangeXD<2, double, std::array>({-980,220}, {-920,520});
    }
    else if(Layer_name=="TIDNeg1_1"){
      Range = Acts::RangeXD<2, double, std::array>({-920,220}, {-860,520});
    }
    else if(Layer_name=="TIDNeg0_2"){
      Range = Acts::RangeXD<2, double, std::array>({-860,220}, {-790,520});
    }
    else if(Layer_name=="TIDNeg0_1"){
      Range = Acts::RangeXD<2, double, std::array>({-790,220}, {-730,520});
    }
    // Default:
    // else if(Layer_name=="TIB3"){
    //   Range = Acts::RangeXD<2, double, std::array>({-730,480}, {730,560});
    // }
    // else if(Layer_name=="TIB2"){
    //   Range = Acts::RangeXD<2, double, std::array>({-730,400}, {730,480});
    // }
    // else if(Layer_name=="TIB1"){
    //   Range = Acts::RangeXD<2, double, std::array>({-730,300}, {730,400});
    // }
    // else if(Layer_name=="TIB0"){
    //   Range = Acts::RangeXD<2, double, std::array>({-730,200}, {730,300});
    // }
    // Extra TIB:
    else if(Layer_name=="TIB7"){
      Range = Acts::RangeXD<2, double, std::array>({-730,500}, {730,535});
    }
    else if(Layer_name=="TIB6"){
      Range = Acts::RangeXD<2, double, std::array>({-730,460}, {730,500});
    }
    else if(Layer_name=="TIB5"){
      Range = Acts::RangeXD<2, double, std::array>({-730,420}, {730,460});
    }
    else if(Layer_name=="TIB4"){
      Range = Acts::RangeXD<2, double, std::array>({-730,380}, {730,420});
    }
    else if(Layer_name=="TIB3"){
      Range = Acts::RangeXD<2, double, std::array>({-730,340}, {730,380});
    }
    else if(Layer_name=="TIB2"){
      Range = Acts::RangeXD<2, double, std::array>({-730,300}, {730,340});
    }
    else if(Layer_name=="TIB1"){
      Range = Acts::RangeXD<2, double, std::array>({-730,255}, {730,300});
    }
    else if(Layer_name=="TIB0"){
      Range = Acts::RangeXD<2, double, std::array>({-730,230}, {730,255});
    }
      // Default TID:
    else if(Layer_name=="TIDPos2"){
      Range = Acts::RangeXD<2, double, std::array>({1000,220}, {1160,530});
    }
    else if(Layer_name=="TIDPos1"){
      Range = Acts::RangeXD<2, double, std::array>({850,220}, {1000,530});
    }
    else if(Layer_name=="TIDPos0"){
      Range = Acts::RangeXD<2, double, std::array>({700,220}, {850,530});
    }
    // Extra TID:
    else if(Layer_name=="TIDPos0_1"){
      Range = Acts::RangeXD<2, double, std::array>({1050,220}, {1130,520});
    }
    else if(Layer_name=="TIDPos0_2"){
      Range = Acts::RangeXD<2, double, std::array>({980,220}, {1050,520});
    }
    else if(Layer_name=="TIDPos1_1"){
      Range = Acts::RangeXD<2, double, std::array>({920,220}, {980,520});
    }
    else if(Layer_name=="TIDPos1_2"){
      Range = Acts::RangeXD<2, double, std::array>({860,220}, {920,520});
    }
    else if(Layer_name=="TIDPos2_1"){
      Range = Acts::RangeXD<2, double, std::array>({790,220}, {860,520});
    }
    else if(Layer_name=="TIDPos2_2"){
      Range = Acts::RangeXD<2, double, std::array>({730,220}, {790,520});
    }
    // else if(Layer_name=="NegTEC8"){
    //   Range = Acts::RangeXD<2, double, std::array>({  -2800,220}, {   -2550,1200});
    // }
      else if(Layer_name=="NegTEC8_2"){
      Range = Acts::RangeXD<2, double, std::array>({  -2790, 220}, { -2667,1200});
    }
    else if(Layer_name=="NegTEC8_1"){
      Range = Acts::RangeXD<2, double, std::array>({  -2667, 220}, { -2554,1200});
    }
    // else if(Layer_name=="NegTEC7"){
    //   Range = Acts::RangeXD<2, double, std::array>({  -2550,220}, {   -2350,1200});
    // }
      else if(Layer_name=="NegTEC7_2"){
      Range = Acts::RangeXD<2, double, std::array>({  -2554, 220}, { -2453,1200});
    }
    else if(Layer_name=="NegTEC7_1"){
      Range = Acts::RangeXD<2, double, std::array>({  -2453, 220}, { -2347,1200});
    }
    // else if(Layer_name=="NegTEC6"){
    //   Range = Acts::RangeXD<2, double, std::array>({  -2350,220}, {   -2150,1200});
    // }
      else if(Layer_name=="NegTEC6_2"){
      Range = Acts::RangeXD<2, double, std::array>({  -2347, 220}, { -2247,1200});
    }
    else if(Layer_name=="NegTEC6_1"){
      Range = Acts::RangeXD<2, double, std::array>({  -2247, 220}, { -2150,1200});
    }
    // else if(Layer_name=="NegTEC5"){
    //   Range = Acts::RangeXD<2, double, std::array>({  -2150,220}, {   -1950,1200});
    // }
      else if(Layer_name=="NegTEC5_2"){
      Range = Acts::RangeXD<2, double, std::array>({  -2150, 220}, { -2058,1200});
    }
    else if(Layer_name=="NegTEC5_1"){
      Range = Acts::RangeXD<2, double, std::array>({  -2058, 220}, { -1950,1200});
    }
    // else if(Layer_name=="NegTEC4"){
    //   Range = Acts::RangeXD<2, double, std::array>({  -1950,220}, {   -1810,1200});
    // }
      else if(Layer_name=="NegTEC4_2"){
      Range = Acts::RangeXD<2, double, std::array>({  -1950, 220}, { -1880,1200});
    }
    else if(Layer_name=="NegTEC4_1"){
      Range = Acts::RangeXD<2, double, std::array>({ -1880, 220}, { -1810,1200});
    }
    // else if(Layer_name=="NegTEC3"){
    //   Range = Acts::RangeXD<2, double, std::array>({  -1810,220}, {   -1675,1200});
    // }
      else if(Layer_name=="NegTEC3_2"){
      Range = Acts::RangeXD<2, double, std::array>({ -1810, 220}, { -1740,1200});
    }
    else if(Layer_name=="NegTEC3_1"){
      Range = Acts::RangeXD<2, double, std::array>({ -1740, 220}, { -1670,1200});
    }
    // else if(Layer_name=="NegTEC2"){
    //   Range = Acts::RangeXD<2, double, std::array>({  -1675,220}, {   -1525,1200});
    // }
      else if(Layer_name=="NegTEC2_2"){
      Range = Acts::RangeXD<2, double, std::array>({ -1670, 220}, { -1600,1200});
    }
    else if(Layer_name=="NegTEC2_1"){
      Range = Acts::RangeXD<2, double, std::array>({ -1600, 220}, { -1530,1200});
    }
    // else if(Layer_name=="NegTEC1"){
    //   Range = Acts::RangeXD<2, double, std::array>({  -1525,220}, {   -1400,1200});
    // }
    else if(Layer_name=="NegTEC1_2"){
      Range = Acts::RangeXD<2, double, std::array>({ -1530, 220}, { -1460,1200});
    }
    else if(Layer_name=="NegTEC1_1"){
      Range = Acts::RangeXD<2, double, std::array>({ -1460, 220}, { -1390,1200});
    }
    // else if(Layer_name=="NegTEC0"){
    //   Range = Acts::RangeXD<2, double, std::array>({  -1400,220}, {   -1200,1200});
    // }
    else if(Layer_name=="NegTEC0_2"){
      Range = Acts::RangeXD<2, double, std::array>({ -1390, 220}, { -1320,1200});
    }
    else if(Layer_name=="NegTEC0_1"){
      Range = Acts::RangeXD<2, double, std::array>({ -1320, 220}, { -1250,1200});
    }
    // Default TOB:
    // else if(Layer_name=="TOB0"){
    //   Range = Acts::RangeXD<2, double, std::array>({  -1160,560}, {    1160, 650});
    // }
    // else if(Layer_name=="TOB1"){
    //   Range = Acts::RangeXD<2, double, std::array>({-1160, 650}, {1160, 750});
    // }
    // else if(Layer_name=="TOB2"){
    //   Range = Acts::RangeXD<2, double, std::array>({-1160, 750}, {1160, 840});
    // }
    // else if(Layer_name=="TOB3"){
    //   Range = Acts::RangeXD<2, double, std::array>({-1160, 840}, {1160, 940});
    // }
    // else if(Layer_name=="TOB4"){
    //   Range = Acts::RangeXD<2, double, std::array>({-1160, 940}, {1160, 1040});
    // }
    // else if(Layer_name=="TOB5"){
    //   Range = Acts::RangeXD<2, double, std::array>({-1160, 1040}, {1160, 1200});
    // }
    // Extra TOB:
    else if(Layer_name=="TOB0"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 560}, {1160, 610});
    }
    else if(Layer_name=="TOB1"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 610}, {1160, 650});
    }
    else if(Layer_name=="TOB2"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 650}, {1160, 695});
    }
    else if(Layer_name=="TOB3"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 695}, {1160, 740});
    }
    else if(Layer_name=="TOB4"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 740}, {1160, 782});
    }
    else if(Layer_name=="TOB5"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 782}, {1160, 828});
    }
    else if(Layer_name=="TOB6"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 828}, {1160, 870});
    }
    else if(Layer_name=="TOB7"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 870}, {1160, 920});
    }
    else if(Layer_name=="TOB8"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 920}, {1160, 967});
    }
    else if(Layer_name=="TOB9"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 967}, {1160, 1030});
    }
    else if(Layer_name=="TOB10"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 1030}, {1160, 1080});
    }
    else if(Layer_name=="TOB11"){
      Range = Acts::RangeXD<2, double, std::array>({-1160, 1080}, {1160, 1135});
    }





    // else if(Layer_name=="PosTEC8"){
    //   Range = Acts::RangeXD<2, double, std::array>({  2550, 220}, {   2800,1200});
    // }
    else if(Layer_name=="PosTEC8_2"){
      Range = Acts::RangeXD<2, double, std::array>({  2667, 220}, {   2790,1200});
    }
    else if(Layer_name=="PosTEC8_1"){
      Range = Acts::RangeXD<2, double, std::array>({  2554, 220}, {   2667,1200});
    }
    // else if(Layer_name=="PosTEC7"){
    //   Range = Acts::RangeXD<2, double, std::array>({  2350, 220}, {   2550,1200});
    // }
    else if(Layer_name=="PosTEC7_2"){
      Range = Acts::RangeXD<2, double, std::array>({  2453, 220}, {   2554,1200});
    }
    else if(Layer_name=="PosTEC7_1"){
      Range = Acts::RangeXD<2, double, std::array>({  2347, 220}, {   2453,1200});
    }
    // else if(Layer_name=="PosTEC6"){
    //   Range = Acts::RangeXD<2, double, std::array>({  2150, 220}, {   2350,1200});
    // }
    else if(Layer_name=="PosTEC6_2"){
      Range = Acts::RangeXD<2, double, std::array>({  2247, 220}, {   2347,1200});
    }
    else if(Layer_name=="PosTEC6_1"){
      Range = Acts::RangeXD<2, double, std::array>({  2150, 220}, {   2247,1200});
    }
    // else if(Layer_name=="PosTEC5"){
    //   Range = Acts::RangeXD<2, double, std::array>({  1950, 220}, {   2150,1200});
    // }
    else if(Layer_name=="PosTEC5_2"){
      Range = Acts::RangeXD<2, double, std::array>({  2058, 220}, {   2150,1200});
    }
    else if(Layer_name=="PosTEC5_1"){
      Range = Acts::RangeXD<2, double, std::array>({  1950, 220}, {   2058,1200});
    }
    // else if(Layer_name=="PosTEC4"){
    //   Range = Acts::RangeXD<2, double, std::array>({  1810, 220}, {   1950,1200});
    // }
    else if(Layer_name=="PosTEC4_2"){
      Range = Acts::RangeXD<2, double, std::array>({  1880, 220}, {   1950,1200});
    }
    else if(Layer_name=="PosTEC4_1"){
      Range = Acts::RangeXD<2, double, std::array>({  1810, 220}, {   1880,1200});
    }
    // else if(Layer_name=="PosTEC3"){
    //   Range = Acts::RangeXD<2, double, std::array>({  1675, 220}, {   1810,1200});
    // }
    else if(Layer_name=="PosTEC3_2"){
      Range = Acts::RangeXD<2, double, std::array>({  1740, 220}, {   1810,1200});
    }
    else if(Layer_name=="PosTEC3_1"){
      Range = Acts::RangeXD<2, double, std::array>({  1670, 220}, {   1740,1200});
    }
    // else if(Layer_name=="PosTEC2"){
    //   Range = Acts::RangeXD<2, double, std::array>({  1525, 220}, {   1675,1200});
    // }
    else if(Layer_name=="PosTEC2_2"){
      Range = Acts::RangeXD<2, double, std::array>({  1600, 220}, {   1670,1200});
    }
    else if(Layer_name=="PosTEC2_1"){
      Range = Acts::RangeXD<2, double, std::array>({  1530, 220}, {   1600,1200});
    }
    // else if(Layer_name=="PosTEC1"){
    //   Range = Acts::RangeXD<2, double, std::array>({  1400, 220}, {   1525,1200});
    // }
    else if(Layer_name=="PosTEC1_2"){
      Range = Acts::RangeXD<2, double, std::array>({  1460, 220}, {   1530,1200});
    }
    else if(Layer_name=="PosTEC1_1"){
      Range = Acts::RangeXD<2, double, std::array>({  1390, 220}, {   1460,1200});
    }
    // else if(Layer_name=="PosTEC0"){
    //   Range = Acts::RangeXD<2, double, std::array>({  1200, 220}, {   1400,1200});
    // }
    else if(Layer_name=="PosTEC0_2"){
      Range = Acts::RangeXD<2, double, std::array>({  1320, 220}, {   1390,1200});
    }
    else if(Layer_name=="PosTEC0_1"){
      Range = Acts::RangeXD<2, double, std::array>({  1250, 220}, {   1320,1200});
    }

    std::vector<std::shared_ptr<Acts::Surface>> ActiveSurfaces = surfaces.surfaces(Range);

    return ActiveSurfaces;
  }

  void makeBinning(auto& layer, Acts::SurfaceArrayNavigationPolicy::LayerType layer_type, double Bin0, double Bin1){
    layer.setNavigationPolicyFactory(Acts::NavigationPolicyFactory()
                              .add<Acts::SurfaceArrayNavigationPolicy>(
                                  Acts::SurfaceArrayNavigationPolicy::Config{
                                      .layerType = layer_type,
                                      .bins = {Bin0, Bin1}}) 
                              .add<Acts::TryAllNavigationPolicy>(
                                  Acts::TryAllNavigationPolicy::Config{.sensitives = true})                                          
                              .asUniquePtr()); 
  }

  void AddExtraLayer(std::string gapName, 
                    bool isBarrel, 
                    Acts::Experimental::CylinderContainerBlueprintNode* cont,
                    Acts::Transform3 transform,
                    std::shared_ptr<Acts::CylinderVolumeBounds> bounds){

    Acts::Transform3 base{Acts::Transform3::Identity()};  
    Acts::AxisDirection DirectionBuild;
    if(!isBarrel){
      DirectionBuild =  Acts::AxisDirection::AxisZ;
    } else {
      DirectionBuild =  Acts::AxisDirection::AxisR;
    }              
    cont->addCylinderContainer(gapName.c_str(), DirectionBuild, [&](auto& gap) {
      gap.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
        .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

      gap.addMaterial(gapName.c_str(), [&](Acts::Experimental::MaterialDesignatorBlueprintNode& mat){
        if(!isBarrel){
          mat.configureFace(Acts::CylinderVolumeBounds::Face::NegativeDisc, {Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
          mat.configureFace(Acts::CylinderVolumeBounds::Face::PositiveDisc, {Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
        }
        else{
          mat.configureFace(Acts::CylinderVolumeBounds::Face::OuterCylinder, {Acts::AxisDirection::AxisRPhi, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, 200});                  
        }

        mat.addCylinderContainer(gapName.c_str(), Acts::AxisDirection::AxisZ, [&](auto& L) {
          L.addStaticVolume(base * transform, bounds, gapName.c_str());
        });
      });       
    });

  };

  struct MatSurfaceSelector : public Acts::TrackingGeometryMutableVisitor {
    std::vector<Acts::Surface*> surfaces;
    void visitSurface(Acts::Surface& surface) override {
      if(surface.surfaceMaterial() != nullptr && !rangeContainsValue(surfaces, &surface)) {
          surfaces.push_back(&surface);
      }
    }
  };

  template <typename MakeLayerFn>
  void AddDiskLayer_and_Material(Acts::Experimental::CylinderContainerBlueprintNode* cont,
                                std::string LayerName,
                                MakeLayerFn&& makeLayerFunc,
                                double z_pos,
                                double bin0,
                                double bin1){
    
      std::string matName = LayerName + "_Material";
      Acts::Transform3 base{Acts::Transform3::Identity()};

      cont->addMaterial(matName.c_str(), [&](Acts::Experimental::MaterialDesignatorBlueprintNode& mat) {
        //mat.configureFace(Acts::CylinderVolumeBounds::Face::OuterCylinder, {Acts::AxisDirection::AxisRPhi, Acts::AxisBoundaryType::Bound, 20}, {Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, 20});
        mat.configureFace(Acts::CylinderVolumeBounds::Face::NegativeDisc, {Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
        mat.configureFace(Acts::CylinderVolumeBounds::Face::PositiveDisc, {Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
        //mat.configureFace(Acts::CylinderVolumeBounds::Face::InnerCylinder, {Acts::AxisDirection::AxisRPhi, Acts::AxisBoundaryType::Bound, 20}, {Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, 20});
        

        mat.addCylinderContainer(LayerName, Acts::AxisDirection::AxisZ, [&](auto& L) {
          L.addLayer(LayerName,[&](auto& layer) {
                        makeLayerFunc(base * Acts::Translation3{Acts::Vector3{0, 0, z_pos}}, layer, LayerName, bin0, bin1);
                    });
        });
      });
  };

  void AddCylinderLayer_and_Material(Acts::Experimental::CylinderContainerBlueprintNode* cont,
                                    std::string LayerName,
                                    KdtSurfacesDim2Bin100& KdtSurfaces, 
                                    double bin0,
                                    double bin1){

    std::string matName = LayerName + "_Material";
    Acts::Transform3 base{Acts::Transform3::Identity()};

    cont->addMaterial(matName.c_str(), [&](Acts::Experimental::MaterialDesignatorBlueprintNode& mat) {
                      mat.configureFace(Acts::CylinderVolumeBounds::Face::OuterCylinder, {Acts::AxisDirection::AxisRPhi, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, 200});
                      //mat.configureFace(Acts::CylinderVolumeBounds::Face::NegativeDisc, {Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 15}, {Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 25});
                      //mat.configureFace(Acts::CylinderVolumeBounds::Face::PositiveDisc, {Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 15}, {Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 25});
                      mat.configureFace(Acts::CylinderVolumeBounds::Face::InnerCylinder, {Acts::AxisDirection::AxisRPhi, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, 200});
                      
                      mat.addCylinderContainer(LayerName, Acts::AxisDirection::AxisR, [&](auto& L) {
                        L.addLayer(LayerName, [&](auto& layer) {
                          std::vector<std::shared_ptr<Acts::Surface>> surfaces = SelectActiveSurfaces_PhaseI(KdtSurfaces, LayerName);
                          // Binning:
                          makeBinning(layer, Acts::SurfaceArrayNavigationPolicy::LayerType::Cylinder, bin0, bin1);

                          if(LayerName.find("TOB") == std::string::npos){
                            layer.setSurfaces(surfaces)
                                .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Cylinder)
                                .setEnvelope(Acts::ExtentEnvelope{{
                                    .z = {1*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm},
                                    .r = {1*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm},
                                }})
                                .setTransform(base);
                          }
                          else{
                            layer.setSurfaces(surfaces)
                                .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Cylinder)
                                .setEnvelope(Acts::ExtentEnvelope{{
                                    .z = {5*Acts::UnitConstants::mm, 5*Acts::UnitConstants::mm},
                                    .r = {1*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm},
                                }})
                                .setTransform(base)
                                .setUseCenterOfGravity(false, false, false);
                          }
                              
                        })
                        .setUseCenterOfGravity(false, false, false); 
                      });
                    });
  };

  class TrackerGeomBuilderWithActsESProducer_LessMatLayers : public edm::ESProducer {
  public:
    explicit TrackerGeomBuilderWithActsESProducer_LessMatLayers(const edm::ParameterSet& ps);
    ~TrackerGeomBuilderWithActsESProducer_LessMatLayers() override = default;
    //std::unique_ptr<Acts::TrackingGeometry> produce(const ACTSTrackerGeometryRecord& iRecord);  // shared to the vector of detector element
    std::shared_ptr<TrackingGeometryWithDetEls> produce(const ACTSTrackerGeometryRecord& iRecord);

  private:
    edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken_;
    edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopoToken_;
    edm::ESGetToken<Alignments, TrackerAlignmentRcd> trackerAlignToken_;

    bool saveObjfile_, saveSvgfile_, mapMaterial_, saveJsonfile_;
    std::string outputObjFile_, outputSvgFile_, materialFile_, ActsLogLevel_;
    std::vector<double> rangeZ_;
    std::vector<double> rangeR_;
  };
}

TrackerGeomBuilderWithActsESProducer_LessMatLayers::TrackerGeomBuilderWithActsESProducer_LessMatLayers(const edm::ParameterSet& ps)
    : saveObjfile_(ps.getUntrackedParameter<bool>("saveObjfile")),
      saveSvgfile_(ps.getUntrackedParameter<bool>("saveSvgfile")), 
      mapMaterial_(ps.getUntrackedParameter<bool>("mapMaterial")), 
      saveJsonfile_(ps.getUntrackedParameter<bool>("saveJsonfile")),
      outputObjFile_(ps.getUntrackedParameter<std::string>("outputObjFile")),
      outputSvgFile_(ps.getUntrackedParameter<std::string>("outputSvgFile")),
      materialFile_(ps.getUntrackedParameter<std::string>("MaterialMaps")),
      ActsLogLevel_(ps.getUntrackedParameter<std::string>("ActsLogLevel")),
      rangeZ_(ps.getUntrackedParameter<std::vector<double>>("rangeZ")),
      rangeR_(ps.getUntrackedParameter<std::vector<double>>("rangeR")) {

    auto cc = setWhatProduced(this);
    trackerGeomToken_ = cc.consumes();
    trackerTopoToken_ = cc.consumes();
    trackerAlignToken_ = cc.consumes();
  }

//std::unique_ptr<Acts::TrackingGeometry> TrackerGeomBuilderWithActsESProducer_LessMatLayers::produce(const ACTSTrackerGeometryRecord& iRecord) { 
std::shared_ptr<TrackingGeometryWithDetEls> TrackerGeomBuilderWithActsESProducer_LessMatLayers::produce(const ACTSTrackerGeometryRecord& iRecord) { 

  const TrackerGeometry& trackerGeom = iRecord.get(trackerGeomToken_);
  const TrackerTopology& trackerTopo = iRecord.get(trackerTopoToken_);
  const Alignments& trackerAlign = iRecord.get(trackerAlignToken_);
  
  DetElVect DetEl_vector;

  std::map<unsigned int, HepGeom::Transform3D> detID_to_alignInfo;

  const Local2DPoint center(0.,0.); 
  const Local3DPoint locz(0.,0.,1.);
  const Local3DPoint locx(1.,0.,0.);
  const Local3DPoint locy(0.,1.,0.);
  const GlobalPoint origin(0.,0.,0.);
  
  // ===== Loop over all the CMSSW detector elments =====
  for (const auto& det : trackerGeom.dets()) { // <-- DetElements might contain NOT physical surfaces (is an assemply of DetElementUnits) that's why I loop on DetElementUnits (sensitive surfaces ONLY)
  // for (const auto& det : trackerGeom.dets()) {
  // (const auto& det : trackerGeom.detUnits()) {

    // ===== Get the surface described with CMSSW =====
    auto cmssw_surf = det->surface();

    // ===== Define the transformation (i.e. Rotation and Translation) of the CMSSW surface =====
    auto& pos = cmssw_surf.position();
    auto& rot = cmssw_surf.rotation();
    Acts::Transform3 t = Acts::Transform3::Identity();
    Acts::RotationMatrix3 R;

    // ===== Get the Alignment INFO ===== (not needed as the aligment info are already stored in the transformation)
    DetId ID = det->geographicalId();
    R << rot.xx(), rot.yx(), rot.zx(), rot.xy(), rot.yy(), rot.zy(), rot.xz(), rot.yz(), rot.zz(); // OLD
    
    // ===== Try a different method to get R and t =====
    GlobalPoint position = trackerGeom.idToDet(ID)->toGlobal(center);
    GlobalPoint zpos = trackerGeom.idToDet(ID)->toGlobal(locz);
    GlobalPoint xpos = trackerGeom.idToDet(ID)->toGlobal(locx);
    GlobalPoint ypos = trackerGeom.idToDet(ID)->toGlobal(locy);
    GlobalVector dz = zpos - position;
    GlobalVector dx = xpos - position;
    GlobalVector dy = ypos - position;

    Eigen::Vector3d dxV(dx.x(), dx.y(), dx.z());
    Eigen::Vector3d dyV(dy.x(), dy.y(), dy.z());
    Eigen::Vector3d dVz(dz.x(), dz.y(), dz.z());

    dxV.normalize();
    dyV.normalize();
    dxV = dxV - dxV.dot(dVz) * dVz; // make dxV ortogonal to dz
    dVz.normalize();

    dyV = dVz.cross(dxV);  // make dy ortogonal to dx e dz

    Eigen::Matrix3d Rot;
    Rot.col(0) = dxV;
    Rot.col(1) = dyV;
    Rot.col(2) = dVz;

    std::string subDet;
    if(ID.subdetId() == PixelSubdetector::PixelBarrel) {
      subDet = "PixelBarrel";
    } 
    else if (ID.subdetId() == PixelSubdetector::PixelEndcap) {
      subDet = "PixelEndcap";
    }     
    else if (ID.subdetId() == StripSubdetector::TID) {
      subDet = "TID";
    }     
    else if (ID.subdetId() == StripSubdetector::TOB) {
      subDet = "TOB";
    }     
    else if (ID.subdetId() == StripSubdetector::TIB) {
      subDet = "TIB";
    }     
    else if (ID.subdetId() == StripSubdetector::TEC)  {
      subDet = "TEC";
    }    

    t.prerotate(Rot);
    t.pretranslate(Acts::Vector3(position.x()*10, position.y()*10, position.z()*10)); // from cm to mm

    // ===== Define the ACTS surface considering two tipes of bounds (i.e. Rectangle and Trapezoid) =====
    std::shared_ptr<Acts::Surface> acts_surf = nullptr;
    auto bounds = dynamic_cast<const RectangularPlaneBounds*>(&cmssw_surf.bounds());
    if (bounds){
      const std::size_t kValues = Acts::RectangleBounds::BoundValues::eSize;
      std::array<double, kValues> bValues{};
      std::vector<double> bVector = {-bounds->width()  / 2 * 10,  // cm → mm
                                     -bounds->length() / 2 * 10,
                                      bounds->width()  / 2 * 10,
                                      bounds->length() / 2 * 10};

      std::copy_n(bVector.begin(), kValues, bValues.begin());
      acts_surf = Acts::Surface::makeShared<Acts::PlaneSurface>(t, std::move(std::make_shared<const Acts::RectangleBounds>(bValues)));
    }
    else{
      auto trap = dynamic_cast<const TrapezoidalPlaneBounds*>(&cmssw_surf.bounds());

      // --- Replace your trapezoid construction block with this (CMSSW -> ACTS, unambiguous) ---
      // Assumes you have:
      //   - const Surface& cmssw_surf;              // the CMSSW surface for this module
      //   - const TrapezoidalPlaneBounds* trap;     // the CMSSW trapezoid bounds (or equivalent ptr/ref)
      //   - Acts::Surface::SharedPtr acts_surf;     // output
      //   - debugCheckTrapezoidConsistency(...)     // your debug helper (optional)

      auto makeActsTransformFromCmsswSurfaceMm = [&](const Surface& s) -> Acts::Transform3 {
        // Build an ACTS transform whose local frame is *exactly* the CMSSW LocalPoint frame.
        // We use 1 cm steps in CMSSW local coordinates to define the basis.
        GlobalPoint g0 = s.toGlobal(LocalPoint(0., 0., 0.));
        GlobalPoint gx = s.toGlobal(LocalPoint(1., 0., 0.));
        GlobalPoint gy = s.toGlobal(LocalPoint(0., 1., 0.));
        GlobalPoint gz = s.toGlobal(LocalPoint(0., 0., 1.));

        Acts::Vector3 O(g0.x() * 10., g0.y() * 10., g0.z() * 10.); // cm->mm
        Acts::Vector3 X((gx.x() - g0.x()) * 10., (gx.y() - g0.y()) * 10., (gx.z() - g0.z()) * 10.);
        Acts::Vector3 Y((gy.x() - g0.x()) * 10., (gy.y() - g0.y()) * 10., (gy.z() - g0.z()) * 10.);
        Acts::Vector3 Z((gz.x() - g0.x()) * 10., (gz.y() - g0.y()) * 10., (gz.z() - g0.z()) * 10.);

        // Normalize and enforce orthonormal, right-handed frame (numerically robust)
        X.normalize();
        // Make Y orthogonal to X
        Y = (Y - (Y.dot(X)) * X);
        Y.normalize();
        // Recompute Z from X x Y to guarantee handedness
        Z = X.cross(Y);
        Z.normalize();

        Acts::Transform3 t = Acts::Transform3::Identity();
        t.linear().col(0) = X;
        t.linear().col(1) = Y;
        t.linear().col(2) = Z;
        t.translation()   = O;
        return t;
      };

      // 1) Read CMSSW trapezoid parameters (they are half-lengths in cm)
      auto params = trap->parameters();
      double halfBottom = params[0] * 10.;  // mm  (BOTTOM edge at y=-halfY in CMSSW local)
      double halfTop    = params[1] * 10.;  // mm  (TOP    edge at y=+halfY in CMSSW local)
      double halfY      = params[3] * 10.;  // mm

      if (halfY <= 0.) {
        throw std::runtime_error("Invalid trapezoid halfY (apothem) <= 0");
      }

      // 2) Build ACTS transform from CMSSW surface frame (this removes all X/Y sign ambiguities)
      Acts::Transform3 tFixed = makeActsTransformFromCmsswSurfaceMm(cmssw_surf);

      // 3) Bounds: map bottom/top to negY/posY
      std::array<double, Acts::TrapezoidBounds::BoundValues::eSize> b{};
      b[Acts::TrapezoidBounds::BoundValues::eHalfLengthXnegY] = halfBottom; // y = -halfY
      b[Acts::TrapezoidBounds::BoundValues::eHalfLengthXposY] = halfTop;    // y = +halfY
      b[Acts::TrapezoidBounds::BoundValues::eHalfLengthY]     = halfY;
      b[Acts::TrapezoidBounds::BoundValues::eRotationAngle]   = 0.0;

      auto bounds = std::make_shared<const Acts::TrapezoidBounds>(b);
      acts_surf = Acts::Surface::makeShared<Acts::PlaneSurface>(tFixed, bounds);

    }

    // ===== Define the CMSDetectorElementData (i.e. the struct needed to define the CMS Detector Element in ACTS) =====
    Acts::CMSDetectorElementData cmsDetData;
    cmsDetData.surf_ = acts_surf;
    cmsDetData.trans_ = t;

    DetId detid = det->geographicalId();
    cmsDetData.detID_ = detid.rawId();
    
    if (detid.subdetId() == PixelSubdetector::PixelBarrel) {
      cmsDetData.subDetector_ = std::string("PixelBarrel");
      cmsDetData.thickness_ = 0.285*Acts::UnitConstants::mm;
    } else if (detid.subdetId() == PixelSubdetector::PixelEndcap) {
      cmsDetData.subDetector_ = std::string("PixelEndcap");
      cmsDetData.thickness_ = 0.3*Acts::UnitConstants::mm;
    } else if (detid.subdetId() == StripSubdetector::TIB) {
      cmsDetData.subDetector_ = std::string("TIB");
      cmsDetData.thickness_ = 0.32*Acts::UnitConstants::mm;
    } else if (detid.subdetId() == StripSubdetector::TID) {
      cmsDetData.subDetector_ = std::string("TID");
      cmsDetData.thickness_ = 0.32*Acts::UnitConstants::mm;
    } else if (detid.subdetId() == StripSubdetector::TOB) {
      cmsDetData.subDetector_ = std::string("TOB");
      cmsDetData.thickness_ = 0.5*Acts::UnitConstants::mm;
    } else if (detid.subdetId() == StripSubdetector::TEC) {
      cmsDetData.subDetector_ = std::string("TEC");
      cmsDetData.thickness_ = 0.5*Acts::UnitConstants::mm;
    }

    // ===== Save the converted detector element in a vector =====
    DetEl_vector.push_back(std::make_shared<Acts::CMSDetectorElement>(cmsDetData));

  }

  // Collect all the surfaces of the detector elements:
  std::vector<std::shared_ptr<Acts::Surface>> surface_vector;
  for(auto DetEl : DetEl_vector){
    surface_vector.push_back((*DetEl).surface().getSharedPtr());
  }

  // Save all the surfaces in a Kdt vector:
  const Acts::GeometryContext viewContext;
  KdtSurfacesDim2Bin100 Kdtsurfaces(viewContext, surface_vector, casts);


  // ===== OPTIONAL: Save the surfaces associate to the detector elements in an obj file =====
  if(saveObjfile_) {

    std::cout << ">>> Flag saveObjfile True <<<" << std::endl;
    std::cout << "Saving detector elements as obj file in the range: " << std::endl;
    std::cout << "R = [" << rangeR_[0] << ", " << rangeR_[1] << "]" << std::endl;
    std::cout << "Z = [" << rangeZ_[0] << ", " << rangeZ_[1] << "]" << std::endl;

    // Select a subset of surfaces:
    std::array<double, 2> min = {rangeZ_[0],rangeR_[0]};
    std::array<double, 2> max = {rangeZ_[1],rangeR_[1]};
    Acts::RangeXD<2, double, std::array> RangeClass(min, max);
    std::vector<std::shared_ptr<Acts::Surface>> surfaceQuery = Kdtsurfaces.surfaces(RangeClass);

    writeSurfacesObj(surfaceQuery, outputObjFile_); 

  }

  // ===== Save all the sensitive surfaces in a JSON file =====
  Acts::GeometryContext geoC;
  std::vector<nlohmann::json> jSurfaces;
  for (const auto& sSurface : surface_vector) {
      jSurfaces.push_back(Acts::SurfaceJsonConverter::toJson(geoC, *sSurface));
  }

  nlohmann::json jSurfacesAll;
  jSurfacesAll["surfaces"] = jSurfaces;

  if(saveJsonfile_) {
    std::cout << ">>> Storing all the surfaces into a json file <<<" << std::endl;
    std::ofstream file("CMSPhaseI_Sensitive_All.json");
    file << jSurfacesAll.dump(4) << '\n';
  }

  // ===== Make the blueprint =====
  Acts::Transform3 base{Acts::Transform3::Identity()};

  Acts::Experimental::Blueprint::Config cfg;
  cfg.envelope[Acts::AxisDirection::AxisZ] = {5*Acts::UnitConstants::mm, 5*Acts::UnitConstants::mm};
  cfg.envelope[Acts::AxisDirection::AxisR] = {10*Acts::UnitConstants::mm, 20*Acts::UnitConstants::mm};
  auto root = std::make_unique<Acts::Experimental::Blueprint>(cfg);

  root->addMaterial("GlobalMaterial", [&](Acts::Experimental::MaterialDesignatorBlueprintNode& mat) {
    using enum Acts::AxisDirection;
    using enum Acts::AxisBoundaryType;
    using enum Acts::CylinderVolumeBounds::Face;

    // Configure cylinder faces with proper binning
    mat.configureFace(Acts::CylinderVolumeBounds::Face::OuterCylinder, {Acts::AxisDirection::AxisRPhi, Acts::AxisBoundaryType::Bound, 20}, {Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, 20});
    mat.configureFace(Acts::CylinderVolumeBounds::Face::NegativeDisc, {Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 15}, {Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 25});
    mat.configureFace(Acts::CylinderVolumeBounds::Face::PositiveDisc, {Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 15}, {Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 25});

    mat.addCylinderContainer("CMS", Acts::AxisDirection::AxisR, [&](auto& CMS) {
      CMS.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

      CMS.addCylinderContainer("SubDetectors", Acts::AxisDirection::AxisZ, [&](auto& det) {
        det.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

        det.addCylinderContainer("NegativeTEC", Acts::AxisDirection::AxisZ, [&](auto& ec) {
              ec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
              //.setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

              auto makeLayer = [&](const Acts::Transform3& trf, Acts::Experimental::LayerBlueprintNode& layer, std::string disk_name, double Bin_Phi, double Bin_R) {
                std::vector<std::shared_ptr<Acts::Surface>> surfaces_disk = SelectActiveSurfaces_PhaseI(Kdtsurfaces, disk_name);
                // Binning:
                makeBinning(layer, Acts::SurfaceArrayNavigationPolicy::LayerType::Disc, Bin_Phi, Bin_R);

                layer.setSurfaces(surfaces_disk)
                    .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Disc)
                    .setEnvelope(Acts::ExtentEnvelope{{
                        .z = {5*Acts::UnitConstants::mm, 5*Acts::UnitConstants::mm},
                        .r = {1*Acts::UnitConstants::mm, 20*Acts::UnitConstants::mm},
                    }})
                    .setTransform(base)
                    .setUseCenterOfGravity(false, false, true); // To fix the transaltion on x-y
              };

              AddDiskLayer_and_Material(&ec, "NegTEC0_1", makeLayer, -1320*Acts::UnitConstants::mm, 80, 4);
              AddDiskLayer_and_Material(&ec, "NegTEC0_2", makeLayer, -1320*Acts::UnitConstants::mm, 80, 4);
              AddDiskLayer_and_Material(&ec, "NegTEC1_1", makeLayer, -1463*Acts::UnitConstants::mm, 80, 5);
              AddDiskLayer_and_Material(&ec, "NegTEC1_2", makeLayer, -1463*Acts::UnitConstants::mm, 80, 5);
              AddDiskLayer_and_Material(&ec, "NegTEC2_1", makeLayer, -1603*Acts::UnitConstants::mm, 80, 5);
              AddDiskLayer_and_Material(&ec, "NegTEC2_2", makeLayer, -1603*Acts::UnitConstants::mm, 80, 5);
              AddDiskLayer_and_Material(&ec, "NegTEC3_1", makeLayer, -1744*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "NegTEC3_2", makeLayer, -1744*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "NegTEC4_1", makeLayer, -1883*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "NegTEC4_2", makeLayer, -1883*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "NegTEC5_1", makeLayer, -2058*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "NegTEC5_2", makeLayer, -2058*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "NegTEC6_1", makeLayer, -2248*Acts::UnitConstants::mm, 80, 7);
              AddDiskLayer_and_Material(&ec, "NegTEC6_2", makeLayer, -2248*Acts::UnitConstants::mm, 80, 7);
              AddDiskLayer_and_Material(&ec, "NegTEC7_1", makeLayer, -2454*Acts::UnitConstants::mm, 80, 7);
              AddDiskLayer_and_Material(&ec, "NegTEC7_2", makeLayer, -2454*Acts::UnitConstants::mm, 80, 7);
              
              AddDiskLayer_and_Material(&ec, "NegTEC8_1", makeLayer, -2666*Acts::UnitConstants::mm, 80, 7);
              ec.addLayer("NegTEC8_2", [&](auto& layer) {
                makeLayer(base * Acts::Translation3{Acts::Vector3{0, 0, -2666*Acts::UnitConstants::mm}}, layer, "NegTEC8_2", 80, 7);
              });


            });
        
        AddExtraLayer("extra_ecneg2", false, &det, GenerateTranslation(0, 0, -1185*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 1200*Acts::UnitConstants::mm, 30*Acts::UnitConstants::mm));
                    
        det.addCylinderContainer("Pixel_TIB_TID_TOB", Acts::AxisDirection::AxisR, [&](auto& barr) { 
          barr.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

          barr.addCylinderContainer("Pixel", Acts::AxisDirection::AxisZ, [&](auto& cyl) { 
            cyl.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

            AddExtraLayer("extra_after_frwPixelNeg1", false, &cyl, GenerateTranslation(0, 0, -800*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 170*Acts::UnitConstants::mm, 10*Acts::UnitConstants::mm));
            AddExtraLayer("extra_after_frwPixelNeg2", false, &cyl, GenerateTranslation(0, 0, -700*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 170*Acts::UnitConstants::mm, 10*Acts::UnitConstants::mm));
            AddExtraLayer("extra_after_frwPixelNeg3", false, &cyl, GenerateTranslation(0, 0, -600*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 170*Acts::UnitConstants::mm, 10*Acts::UnitConstants::mm));
             
            cyl.addCylinderContainer(
                "PixelNegativeEndcap", Acts::AxisDirection::AxisZ, [&](auto& ec) {
                  ec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

                  auto makeLayer = [&](const Acts::Transform3& trf, auto& layer, std::string disk_name, double Bin_Phi, double Bin_R) {
                    std::vector<std::shared_ptr<Acts::Surface>> surfaces_disk = SelectActiveSurfaces_PhaseI(Kdtsurfaces, disk_name);
                    // Binning:
                    makeBinning(layer, Acts::SurfaceArrayNavigationPolicy::LayerType::Disc, Bin_Phi, Bin_R);

                    layer.setSurfaces(surfaces_disk)
                        .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Disc)
                        .setEnvelope(Acts::ExtentEnvelope{{
                            .z = {5*Acts::UnitConstants::mm, 5*Acts::UnitConstants::mm},
                            .r = {10*Acts::UnitConstants::mm, 20*Acts::UnitConstants::mm},
                        }})
                        .setTransform(base)
                        .setUseCenterOfGravity(false, false, true); // To fix the transaltion on x-y
                  };

                  AddDiskLayer_and_Material(&ec, "PixelNeg0", makeLayer, -322*Acts::UnitConstants::mm, 36, 2);
                  AddDiskLayer_and_Material(&ec, "PixelNeg1", makeLayer, -395*Acts::UnitConstants::mm, 36, 2);
                  AddDiskLayer_and_Material(&ec, "PixelNeg2", makeLayer, -493*Acts::UnitConstants::mm, 36, 2);

                });
            
            cyl.addCylinderContainer(
                "PixelBarrel", Acts::AxisDirection::AxisR, [&](auto& brl) {
                  brl.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                      .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);


                  AddCylinderLayer_and_Material(&brl, "PixL0", Kdtsurfaces, 12, 8);
                  AddCylinderLayer_and_Material(&brl, "PixL1", Kdtsurfaces, 32, 8);
                  AddCylinderLayer_and_Material(&brl, "PixL2", Kdtsurfaces, 48, 8);
                  AddCylinderLayer_and_Material(&brl, "PixL3", Kdtsurfaces, 64, 8);
                });
            
            cyl.addCylinderContainer(
                "PixelPositiveEndcap", Acts::AxisDirection::AxisZ, [&](auto& ec) {
                  ec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

                  auto makeLayer = [&](const Acts::Transform3& trf, auto& layer, std::string disk_name, double Bin_Phi, double Bin_R) {
                    std::vector<std::shared_ptr<Acts::Surface>> surfaces_disk = SelectActiveSurfaces_PhaseI(Kdtsurfaces, disk_name);
                    // Binning:
                    makeBinning(layer, Acts::SurfaceArrayNavigationPolicy::LayerType::Disc, Bin_Phi, Bin_R);

                    layer.setSurfaces(surfaces_disk)
                        .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Disc)
                        .setEnvelope(Acts::ExtentEnvelope{{
                            .z = {5*Acts::UnitConstants::mm, 5*Acts::UnitConstants::mm},
                            .r = {10*Acts::UnitConstants::mm, 20*Acts::UnitConstants::mm},
                        }})
                        .setTransform(base)
                        .setUseCenterOfGravity(false, false, true); // To fix the transaltion on x-y
                  };

                  AddDiskLayer_and_Material(&ec, "PixelPos0", makeLayer, 322*Acts::UnitConstants::mm, 36, 2);
                  AddDiskLayer_and_Material(&ec, "PixelPos1", makeLayer, 395*Acts::UnitConstants::mm, 36, 2);
                  AddDiskLayer_and_Material(&ec, "PixelPos2", makeLayer, 493*Acts::UnitConstants::mm, 36, 2);
                });

            AddExtraLayer("extra_after_frwPixelPos1", false, &cyl, GenerateTranslation(0, 0, 600*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 170*Acts::UnitConstants::mm, 10*Acts::UnitConstants::mm));
            AddExtraLayer("extra_after_frwPixelPos2", false, &cyl, GenerateTranslation(0, 0, 700*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 170*Acts::UnitConstants::mm, 10*Acts::UnitConstants::mm));
            AddExtraLayer("extra_after_frwPixelPos3", false, &cyl, GenerateTranslation(0, 0, 800*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 170*Acts::UnitConstants::mm, 10*Acts::UnitConstants::mm));
              
          });

          AddExtraLayer("ExtraPix_TIDMatLayer1", true, &barr, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(189*Acts::UnitConstants::mm, 190*Acts::UnitConstants::mm, 1115*Acts::UnitConstants::mm));
          AddExtraLayer("ExtraPix_TIDMatLayer3", true, &barr, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(202*Acts::UnitConstants::mm, 203*Acts::UnitConstants::mm, 1115*Acts::UnitConstants::mm));
          
          barr.addCylinderContainer("TIB_TID", Acts::AxisDirection::AxisZ, [&](auto& cyl) {
            cyl.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

            AddExtraLayer("fixingHeight", false, &cyl, GenerateTranslation(0, 0, -1110*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(220*Acts::UnitConstants::mm, 500*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));
        
            
            cyl.addCylinderContainer(
                "NegativeTID", Acts::AxisDirection::AxisZ, [&](auto& ec) {
                  ec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

                  auto makeLayer = [&](const Acts::Transform3& trf, auto& layer, std::string disk_name, double Bin_Phi, double Bin_R) {
                    std::vector<std::shared_ptr<Acts::Surface>> surfaces_disk = SelectActiveSurfaces_PhaseI(Kdtsurfaces, disk_name);
                    // Binning:
                    makeBinning(layer, Acts::SurfaceArrayNavigationPolicy::LayerType::Disc, Bin_Phi, Bin_R);

                    layer.setSurfaces(surfaces_disk)
                        .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Disc)
                        .setEnvelope(Acts::ExtentEnvelope{{
                            .z = {5*Acts::UnitConstants::mm, 5*Acts::UnitConstants::mm},
                            .r = {10*Acts::UnitConstants::mm, 20*Acts::UnitConstants::mm},
                        }})
                        .setTransform(base)
                        .setUseCenterOfGravity(false, false, true); // To fix the transaltion on x-y
                  };

                  // AddDiskLayer_and_Material(&ec, "TIDNeg0", makeLayer, -790*Acts::UnitConstants::mm, 40, 3);
                  // AddDiskLayer_and_Material(&ec, "TIDNeg0_1", makeLayer, -790*Acts::UnitConstants::mm, 40, 3);
                  // AddDiskLayer_and_Material(&ec, "TIDNeg0_2", makeLayer, -790*Acts::UnitConstants::mm, 40, 3);

                  // AddDiskLayer_and_Material(&ec, "TIDNeg1", makeLayer, -920*Acts::UnitConstants::mm, 40, 3);
                  // AddDiskLayer_and_Material(&ec, "TIDNeg1_1", makeLayer, -920*Acts::UnitConstants::mm, 40, 3);
                  // AddDiskLayer_and_Material(&ec, "TIDNeg1_2", makeLayer, -920*Acts::UnitConstants::mm, 40, 3);

                  // AddDiskLayer_and_Material(&ec, "TIDNeg2", makeLayer, -1050*Acts::UnitConstants::mm, 40, 3);
                  // ec.addLayer("TIDNeg2", [&](auto& layer) {
                  //   makeLayer(base * Acts::Translation3{Acts::Vector3{0, 0, -1050*Acts::UnitConstants::mm}}, layer, "TIDNeg2", 40, 3);
                  // });

                  // ec.addLayer("TIDNeg2_2", [&](auto& layer) {
                  //   makeLayer(base * Acts::Translation3{Acts::Vector3{0, 0, -1050*Acts::UnitConstants::mm}}, layer, "TIDNeg2_2", 40, 3);
                  // });

                  AddExtraLayer("ExtraNegTIDMatLayer2", false, &ec, GenerateTranslation(0, 0, -722*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(250*Acts::UnitConstants::mm, 500*Acts::UnitConstants::mm, 8*Acts::UnitConstants::mm));

                  AddDiskLayer_and_Material(&ec, "TIDNeg0_1", makeLayer, -790*Acts::UnitConstants::mm, 40, 3);
                  AddDiskLayer_and_Material(&ec, "TIDNeg0_2", makeLayer, -790*Acts::UnitConstants::mm, 40, 3);
                  AddDiskLayer_and_Material(&ec, "TIDNeg1_1", makeLayer, -920*Acts::UnitConstants::mm, 40, 3);
                  AddDiskLayer_and_Material(&ec, "TIDNeg1_2", makeLayer, -920*Acts::UnitConstants::mm, 40, 3);
                  AddDiskLayer_and_Material(&ec, "TIDNeg2_1", makeLayer, -1050*Acts::UnitConstants::mm, 40, 3);
                  AddDiskLayer_and_Material(&ec, "TIDNeg2_2", makeLayer, -1050*Acts::UnitConstants::mm, 40, 3);
                });
              
            cyl.addCylinderContainer(
                "TIB", Acts::AxisDirection::AxisR, [&](auto& brl) {
                  brl.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                      //.setResizeStrategy(Acts::VolumeResizeStrategy::Gap);
                      .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

                  // auto makeLayer = [&](const std::string& name, double Bin_Phi, double Bin_Z) {
                  //   brl.addLayer(name, [&](auto& layer) {
                  //     std::vector<std::shared_ptr<Acts::Surface>> surfaces = SelectActiveSurfaces_PhaseI(Kdtsurfaces, name);
                  //     // Binning:
                  //     makeBinning(layer, Acts::SurfaceArrayNavigationPolicy::LayerType::Cylinder, Bin_Phi, Bin_Z);

                  //     layer.setSurfaces(surfaces)
                  //         .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Cylinder)
                  //         .setEnvelope(Acts::ExtentEnvelope{{
                  //             .z = {5*Acts::UnitConstants::mm, 5*Acts::UnitConstants::mm},
                  //             .r = {1*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm},
                  //         }})
                  //         .setTransform(base)
                  //         .setUseCenterOfGravity(false, false, false); // To fix the transaltion on x-y
                  //   });
                  // };
                  
                  // Default:
                  // AddCylinderLayer_and_Material(&brl, "TIB0", Kdtsurfaces, 26, 10);
                  // AddExtraLayer("ExtraTIBMatLayer1", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(300*Acts::UnitConstants::mm, 301*Acts::UnitConstants::mm, 700*Acts::UnitConstants::mm));   
                  // AddCylinderLayer_and_Material(&brl, "TIB1", Kdtsurfaces, 38, 10);
                  // AddExtraLayer("ExtraTIBMatLayer2", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(380*Acts::UnitConstants::mm, 381*Acts::UnitConstants::mm, 700*Acts::UnitConstants::mm));   
                  // AddCylinderLayer_and_Material(&brl, "TIB2", Kdtsurfaces, 44, 10);
                  // AddExtraLayer("ExtraTIBMatLayer3", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(460*Acts::UnitConstants::mm, 461*Acts::UnitConstants::mm, 700*Acts::UnitConstants::mm));   
                  // AddCylinderLayer_and_Material(&brl, "TIB3", Kdtsurfaces, 56, 10);
                  
                  // Extra: They are not symmetric in Z, but the extra layers force the overall volume to be symmetric
                  AddCylinderLayer_and_Material(&brl, "TIB0", Kdtsurfaces, 26, 10);
                  AddCylinderLayer_and_Material(&brl, "TIB1", Kdtsurfaces, 26, 10);
                  AddExtraLayer("ExtraTIBMatLayer1", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(290*Acts::UnitConstants::mm, 291*Acts::UnitConstants::mm, 700*Acts::UnitConstants::mm));   
                  AddExtraLayer("ExtraTIBMatLayer2", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(300*Acts::UnitConstants::mm, 301*Acts::UnitConstants::mm, 700*Acts::UnitConstants::mm));   
                  AddCylinderLayer_and_Material(&brl, "TIB2", Kdtsurfaces, 38, 10);
                  AddCylinderLayer_and_Material(&brl, "TIB3", Kdtsurfaces, 38, 10);
                  AddExtraLayer("ExtraTIBMatLayer3", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(375*Acts::UnitConstants::mm, 376*Acts::UnitConstants::mm, 700*Acts::UnitConstants::mm));   
                  AddExtraLayer("ExtraTIBMatLayer4", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(385*Acts::UnitConstants::mm, 386*Acts::UnitConstants::mm, 700*Acts::UnitConstants::mm));   
                  AddCylinderLayer_and_Material(&brl, "TIB4", Kdtsurfaces, 44, 10);
                  AddCylinderLayer_and_Material(&brl, "TIB5", Kdtsurfaces, 44, 10);
                  AddExtraLayer("ExtraTIBMatLayer5", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(455*Acts::UnitConstants::mm, 456*Acts::UnitConstants::mm, 700*Acts::UnitConstants::mm));   
                  AddExtraLayer("ExtraTIBMatLayer6", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(465*Acts::UnitConstants::mm, 466*Acts::UnitConstants::mm, 700*Acts::UnitConstants::mm));   
                  AddCylinderLayer_and_Material(&brl, "TIB6", Kdtsurfaces, 56, 10);
                  AddCylinderLayer_and_Material(&brl, "TIB7", Kdtsurfaces, 56, 10);

                  //makeLayer("TIB3", 56,10);
                });
            
            cyl.addCylinderContainer(
                "PositiveTID", Acts::AxisDirection::AxisZ, [&](auto& ec) {
                  ec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

                  auto makeLayer = [&](const Acts::Transform3& trf, auto& layer, std::string disk_name, double Bin_Phi, double Bin_R) {
                    std::vector<std::shared_ptr<Acts::Surface>> surfaces_disk = SelectActiveSurfaces_PhaseI(Kdtsurfaces, disk_name);
                    // Binning:
                    makeBinning(layer, Acts::SurfaceArrayNavigationPolicy::LayerType::Disc, Bin_Phi, Bin_R);

                    layer.setSurfaces(surfaces_disk)
                        .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Disc)
                        .setEnvelope(Acts::ExtentEnvelope{{
                            .z = {1*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm},
                            .r = {1*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm},
                        }})
                        .setTransform(base)
                        .setUseCenterOfGravity(false, false, true); // To fix the transaltion on x-y
                  };
                  //makeLayer(translation, layer_to_move, "Layer_Name", Bin_Phi, Bin_R)
                  // AddDiskLayer_and_Material(&ec, "TIDPos0", makeLayer, 790*Acts::UnitConstants::mm, 40, 3);
                  // AddDiskLayer_and_Material(&ec, "TIDPos1", makeLayer, 920*Acts::UnitConstants::mm, 40, 3);
                  // ec.addLayer("TIDPos2", [&](auto& layer) {
                  //   makeLayer(base * Acts::Translation3{Acts::Vector3{0, 0, 1050*Acts::UnitConstants::mm}}, layer, "TIDPos2", 40, 3);
                  // });

                  AddExtraLayer("ExtraPosTIDMatLayer2", false, &ec, GenerateTranslation(0, 0, 722*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(250*Acts::UnitConstants::mm, 500*Acts::UnitConstants::mm, 8*Acts::UnitConstants::mm));

                  AddDiskLayer_and_Material(&ec, "TIDPos0_1", makeLayer, 790*Acts::UnitConstants::mm, 40, 3);
                  AddDiskLayer_and_Material(&ec, "TIDPos0_2", makeLayer, 790*Acts::UnitConstants::mm, 40, 3);
                  AddDiskLayer_and_Material(&ec, "TIDPos1_1", makeLayer, 920*Acts::UnitConstants::mm, 40, 3);
                  AddDiskLayer_and_Material(&ec, "TIDPos1_2", makeLayer, 920*Acts::UnitConstants::mm, 40, 3);
                  AddDiskLayer_and_Material(&ec, "TIDPos2_1", makeLayer, 1050*Acts::UnitConstants::mm, 40, 3);
                  AddDiskLayer_and_Material(&ec, "TIDPos2_2", makeLayer, 1050*Acts::UnitConstants::mm, 40, 3);
            });
            
          });
          
          AddExtraLayer("ExtraTIB_TOBMatLayer2", true, &barr, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(545*Acts::UnitConstants::mm, 546*Acts::UnitConstants::mm, 1115*Acts::UnitConstants::mm));
          
          barr.addCylinderContainer(
              "TOB", Acts::AxisDirection::AxisR, [&](auto& brl) {
                brl.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                    .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

                // auto makeLayer = [&](const std::string& name, double Bin_Phi, double Bin_Z) {
                //     brl.addLayer(name, [&](auto& layer) {
                //       std::vector<std::shared_ptr<Acts::Surface>> surfaces = SelectActiveSurfaces_PhaseI(Kdtsurfaces, name);
                //       // Binning:
                //       makeBinning(layer, Acts::SurfaceArrayNavigationPolicy::LayerType::Cylinder, Bin_Phi, Bin_Z);

                //       layer.setSurfaces(surfaces)
                //           .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Cylinder)
                //           .setEnvelope(Acts::ExtentEnvelope{{
                //               .z = {5*Acts::UnitConstants::mm, 5*Acts::UnitConstants::mm},
                //               .r = {1*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm},
                //           }})
                //           .setTransform(base)
                //           .setUseCenterOfGravity(false, false, false); // To fix the transaltion on x-y
                //     });
                //   };

                // makeLayer("Layer_name", Bin_Phi, Bin_Z)
                // AddCylinderLayer_and_Material(&brl, "TOB0", Kdtsurfaces, 42, 24);
                // AddExtraLayer("ExtraTOBMatLayer1", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(650*Acts::UnitConstants::mm, 655*Acts::UnitConstants::mm, 1000*Acts::UnitConstants::mm));
                // AddCylinderLayer_and_Material(&brl, "TOB1", Kdtsurfaces, 48, 24);
                // AddExtraLayer("ExtraTOBMatLayer2", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(740*Acts::UnitConstants::mm, 750*Acts::UnitConstants::mm, 1000*Acts::UnitConstants::mm));
                // AddCylinderLayer_and_Material(&brl, "TOB2", Kdtsurfaces, 54, 24);
                // AddExtraLayer("ExtraTOBMatLayer3", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(805*Acts::UnitConstants::mm, 810*Acts::UnitConstants::mm, 1000*Acts::UnitConstants::mm));
                // AddCylinderLayer_and_Material(&brl, "TOB3", Kdtsurfaces, 60, 24);
                // AddExtraLayer("ExtraTOBMatLayer4", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(940*Acts::UnitConstants::mm, 943*Acts::UnitConstants::mm, 1000*Acts::UnitConstants::mm));
                // AddCylinderLayer_and_Material(&brl, "TOB4", Kdtsurfaces, 66, 24);
                // AddExtraLayer("ExtraTOBMatLayer5", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(1000*Acts::UnitConstants::mm, 1050*Acts::UnitConstants::mm, 1000*Acts::UnitConstants::mm));
                // AddCylinderLayer_and_Material(&brl, "TOB5", Kdtsurfaces, 74, 24);

                
                AddCylinderLayer_and_Material(&brl, "TOB0", Kdtsurfaces, 42, 24);
                AddCylinderLayer_and_Material(&brl, "TOB1", Kdtsurfaces, 42, 24);
                AddExtraLayer("ExtraTOBMatLayer1", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(650*Acts::UnitConstants::mm, 655*Acts::UnitConstants::mm, 1115*Acts::UnitConstants::mm));   
                AddCylinderLayer_and_Material(&brl, "TOB2", Kdtsurfaces, 48, 24);
                AddCylinderLayer_and_Material(&brl, "TOB3", Kdtsurfaces, 48, 24);
                AddExtraLayer("ExtraTOBMatLayer2", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(740*Acts::UnitConstants::mm, 750*Acts::UnitConstants::mm, 1115*Acts::UnitConstants::mm));   
                AddCylinderLayer_and_Material(&brl, "TOB4", Kdtsurfaces, 54, 24);
                AddCylinderLayer_and_Material(&brl, "TOB5", Kdtsurfaces, 54, 24);
                AddExtraLayer("ExtraTOBMatLayer3", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(825*Acts::UnitConstants::mm, 826*Acts::UnitConstants::mm, 1115*Acts::UnitConstants::mm));   
                AddCylinderLayer_and_Material(&brl, "TOB6", Kdtsurfaces, 60, 24);
                AddCylinderLayer_and_Material(&brl, "TOB7", Kdtsurfaces, 60, 24);
                AddExtraLayer("ExtraTOBMatLayer4", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(920*Acts::UnitConstants::mm, 921*Acts::UnitConstants::mm, 1115*Acts::UnitConstants::mm));   
                AddCylinderLayer_and_Material(&brl, "TOB8", Kdtsurfaces, 66, 24);
                AddCylinderLayer_and_Material(&brl, "TOB9", Kdtsurfaces, 66, 24);
                AddExtraLayer("ExtraTOBMatLayer5", true, &brl, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(1025*Acts::UnitConstants::mm, 1026*Acts::UnitConstants::mm, 1115*Acts::UnitConstants::mm));   
                AddCylinderLayer_and_Material(&brl, "TOB10", Kdtsurfaces, 74, 24);
                AddCylinderLayer_and_Material(&brl, "TOB11", Kdtsurfaces, 74, 24);
                


                //makeLayer("TOB5", 74,10);
              });
        });
        
        det.addCylinderContainer("PositiveTEC", Acts::AxisDirection::AxisZ, [&](auto& ec) {
            ec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
              //.setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

              auto makeLayer = [&](const Acts::Transform3& trf, auto& layer, std::string disk_name, double Bin_Phi, double Bin_R) {
                std::vector<std::shared_ptr<Acts::Surface>> surfaces_disk = SelectActiveSurfaces_PhaseI(Kdtsurfaces, disk_name);
                // Binning:
                makeBinning(layer, Acts::SurfaceArrayNavigationPolicy::LayerType::Disc, Bin_Phi, Bin_R);

                layer.setSurfaces(surfaces_disk)
                    .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Disc)
                    .setEnvelope(Acts::ExtentEnvelope{{
                        .z = {5*Acts::UnitConstants::mm, 5*Acts::UnitConstants::mm},
                        .r = {1*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm},
                    }})
                    .setTransform(base)
                    .setUseCenterOfGravity(false, false, true); // To fix the transaltion on x-y
              };

              AddDiskLayer_and_Material(&ec, "PosTEC0_1", makeLayer, 1320*Acts::UnitConstants::mm, 80, 4);
              AddDiskLayer_and_Material(&ec, "PosTEC0_2", makeLayer, 1320*Acts::UnitConstants::mm, 80, 4);
              AddDiskLayer_and_Material(&ec, "PosTEC1_1", makeLayer, 1463*Acts::UnitConstants::mm, 80, 5);
              AddDiskLayer_and_Material(&ec, "PosTEC1_2", makeLayer, 1463*Acts::UnitConstants::mm, 80, 5);
              AddDiskLayer_and_Material(&ec, "PosTEC2_1", makeLayer, 1603*Acts::UnitConstants::mm, 80, 5);
              AddDiskLayer_and_Material(&ec, "PosTEC2_2", makeLayer, 1603*Acts::UnitConstants::mm, 80, 5);
              AddDiskLayer_and_Material(&ec, "PosTEC3_1", makeLayer, 1744*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "PosTEC3_2", makeLayer, 1744*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "PosTEC4_1", makeLayer, 1883*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "PosTEC4_2", makeLayer, 1883*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "PosTEC5_1", makeLayer, 2058*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "PosTEC5_2", makeLayer, 2058*Acts::UnitConstants::mm, 80, 6);
              AddDiskLayer_and_Material(&ec, "PosTEC6_1", makeLayer, 2248*Acts::UnitConstants::mm, 80, 7);
              AddDiskLayer_and_Material(&ec, "PosTEC6_2", makeLayer, 2248*Acts::UnitConstants::mm, 80, 7);
              AddDiskLayer_and_Material(&ec, "PosTEC7_1", makeLayer, 2454*Acts::UnitConstants::mm, 80, 7);
              AddDiskLayer_and_Material(&ec, "PosTEC7_2", makeLayer, 2454*Acts::UnitConstants::mm, 80, 7);
              AddDiskLayer_and_Material(&ec, "PosTEC8_1", makeLayer, 2666*Acts::UnitConstants::mm, 80, 7);
              ec.addLayer("PosTEC8_2", [&](auto& layer) {
                makeLayer(base * Acts::Translation3{Acts::Vector3{0, 0, 2666*Acts::UnitConstants::mm}}, layer, "PosTEC8_2", 80, 7);
              });

            });


        AddExtraLayer("extra_ecpos2", false, &det, GenerateTranslation(0, 0, 1185*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 1200*Acts::UnitConstants::mm, 30*Acts::UnitConstants::mm));
   
      });
 
      // CMS.addStaticVolume(base, std::make_shared<Acts::CylinderVolumeBounds>(0*Acts::UnitConstants::mm, 10*Acts::UnitConstants::mm, 3000*Acts::UnitConstants::mm), "BeamPipe");

      CMS.addCylinderContainer("BeampipeVolume", Acts::AxisDirection::AxisR, [&](auto& bp) {
        AddExtraLayer("BeamPipe", true, &bp, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(0*Acts::UnitConstants::mm, 24*Acts::UnitConstants::mm, 3000*Acts::UnitConstants::mm));          
      });

    });
  });

  // ===== Construct the TrackingGeometry from the blueprint =====
  Acts::GeometryContext gctx;
  auto logLevel = (ActsLogLevel_ == "verbose") ? Acts::Logging::VERBOSE : Acts::Logging::INFO;
  auto logger = Acts::getDefaultLogger("UnitTests", logLevel);
  Acts::Experimental::BlueprintOptions BluePrint_otp;
  std::shared_ptr<Acts::TrackingGeometry> trackingGeometry = std::move(root->construct(BluePrint_otp, gctx, *logger));

  // ===== OPTIONAL: Save the surfaces associate to the detector elements in an SVG file =====
  if(saveSvgfile_) {
    std::cout << ">>> Flag saveSvgfile True <<<" << std::endl;
    std::cout << "Saving tracking geometry blueprint in " << outputSvgFile_ << std::endl;
    actsvg::views::z_r view_zr; 
    auto svg_objs_rz = ActsPlugins::Svg::drawTrackingGeometry(gctx, *trackingGeometry, view_zr, true, true);
    ActsPlugins::Svg::toFile(svg_objs_rz, outputSvgFile_);
  }

  // ===== Decorate the material using a JSON material file =====
  if (mapMaterial_) {
    std::cout << ">>> Flag mapMaterial_ True <<<" << std::endl;
    std::cout << "Mapping material on sensitive surfaces, from " << materialFile_ << std::endl;
    Acts::MaterialMapJsonConverter::Config dec_cfg;
    Acts::JsonMaterialDecorator jsonMatDec(dec_cfg, materialFile_, Acts::Logging::Level::VERBOSE);

    MatSurfaceSelector sel;
    trackingGeometry->apply(sel);
    std::vector<Acts::Surface*> surfVec = sel.surfaces;

    // Loop over the surfaces and decorate them:
    for(auto& surf : surfVec){
      jsonMatDec.decorate(*surf);
    }
  }

  auto result = std::make_shared<TrackingGeometryWithDetEls>();
    result->detElements = DetEl_vector;
    result->trackingGeometry = trackingGeometry;
  
  return result;
}

DEFINE_FWK_EVENTSETUP_MODULE(TrackerGeomBuilderWithActsESProducer_LessMatLayers);