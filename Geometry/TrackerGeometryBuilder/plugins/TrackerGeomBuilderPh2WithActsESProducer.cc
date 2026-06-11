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

#include "Geometry/TrackerGeometryBuilder/interface/ActsGeoBuilderUtils.h"

using json = nlohmann::json;
const std::array<Acts::AxisDirection, 2UL> casts{Acts::AxisDirection::AxisZ, Acts::AxisDirection::AxisR};

using DetElVect = std::vector<std::shared_ptr<Acts::CMSDetectorElement>>;

namespace {
  struct myContext{
    /// Magnetic and Geometry contrext
    Acts::GeometryContext geoContext;
    Acts::MagneticFieldContext magFieldContext;
    /// Number of event
    std::size_t EvNumber; 
  };

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

  void AddExtraLayer_noMat(std::string gapName, 
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

      gap.addStaticVolume(base * transform, bounds, gapName.c_str());    
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
                                double bin0,
                                double bin1){
    
      std::string matName = LayerName + "_Material";

      const bool isTilted = LayerName.find("Tilted") != std::string::npos;
      const bool isPosTilted = LayerName.find("PosTilted") != std::string::npos;
      const bool isNegTilted = LayerName.find("NegTilted") != std::string::npos;

      cont->addMaterial(matName.c_str(), [&](Acts::Experimental::MaterialDesignatorBlueprintNode& mat) {
        // if (isTilted) {       
        //   if (isPosTilted) {
        //     // Positive tilted disk: material only on negative-z face
        //     mat.configureFace(Acts::CylinderVolumeBounds::Face::PositiveDisc,{Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200},{Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
        //     mat.configureFace(Acts::CylinderVolumeBounds::Face::NegativeDisc,{Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200},{Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
        //   } else if (isNegTilted) {
        //     // Negative tilted disk: material only on positive-z face
        //     mat.configureFace(Acts::CylinderVolumeBounds::Face::PositiveDisc,{Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200},{Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
        //   }
        // } else {
        //   // Normal disk layer
        //   mat.configureFace(Acts::CylinderVolumeBounds::Face::NegativeDisc,{Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200},{Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
        //   mat.configureFace(Acts::CylinderVolumeBounds::Face::PositiveDisc,{Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200},{Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
        // }

        // Normal disk layer
        mat.configureFace(Acts::CylinderVolumeBounds::Face::NegativeDisc,{Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200},{Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
        mat.configureFace(Acts::CylinderVolumeBounds::Face::PositiveDisc,{Acts::AxisDirection::AxisR, Acts::AxisBoundaryType::Bound, 200},{Acts::AxisDirection::AxisPhi, Acts::AxisBoundaryType::Bound, 200});
        

        mat.addCylinderContainer(LayerName, Acts::AxisDirection::AxisZ, [&](auto& L) {
          L.addLayer(LayerName,[&](auto& layer) {
                        makeLayerFunc(layer, LayerName, bin0, bin1);
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
                      // const bool test = LayerName.find("TBPS") != std::string::npos;
                      // if(!test) {
                      //   mat.configureFace(Acts::CylinderVolumeBounds::Face::OuterCylinder, {Acts::AxisDirection::AxisRPhi, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, 200});
                      //   mat.configureFace(Acts::CylinderVolumeBounds::Face::InnerCylinder, {Acts::AxisDirection::AxisRPhi, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, 200});  
                      // }
                      mat.configureFace(Acts::CylinderVolumeBounds::Face::OuterCylinder, {Acts::AxisDirection::AxisRPhi, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, 200});
                      mat.configureFace(Acts::CylinderVolumeBounds::Face::InnerCylinder, {Acts::AxisDirection::AxisRPhi, Acts::AxisBoundaryType::Bound, 200}, {Acts::AxisDirection::AxisZ, Acts::AxisBoundaryType::Bound, 200});  
                      
                      
                      mat.addCylinderContainer(LayerName, Acts::AxisDirection::AxisR, [&](auto& L) {
                        L.addLayer(LayerName, [&](auto& layer) {
                          ActsGeoBuilderUtils geoUtils;
                          std::vector<std::shared_ptr<Acts::Surface>> surfaces = geoUtils.SelectActiveSurfaces_PhaseII(KdtSurfaces, LayerName);
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
                                    .z = {1*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm},
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

}

class TrackerGeomBuilderPh2WithActsESProducer : public edm::ESProducer {
public:
  explicit TrackerGeomBuilderPh2WithActsESProducer(const edm::ParameterSet& ps);
  ~TrackerGeomBuilderPh2WithActsESProducer() override = default;
  //std::unique_ptr<Acts::TrackingGeometry> produce(const ACTSTrackerGeometryRecord& iRecord);  // shared to the vector of detector element
  std::shared_ptr<TrackingGeometryWithDetEls> produce(const ACTSTrackerGeometryRecord& iRecord);

private:
  template <typename KdtType, typename LayerType>
  void makeLayer(const KdtType& kdt, LayerType& layer, const std::string& diskName, double binPhi, double binR);

  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken_;
  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopoToken_;
  edm::ESGetToken<Alignments, TrackerAlignmentRcd> trackerAlignToken_;

  bool saveObjfile_, saveSvgfile_, mapMaterial_, saveJsonfile_;
  std::string outputObjFile_, outputSvgFile_, materialFile_, ActsLogLevel_;
  std::vector<double> rangeZ_;
  std::vector<double> rangeR_;
};


TrackerGeomBuilderPh2WithActsESProducer::TrackerGeomBuilderPh2WithActsESProducer(const edm::ParameterSet& ps)
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

template <typename KdtType, typename LayerType>
void TrackerGeomBuilderPh2WithActsESProducer::makeLayer(const KdtType& Kdtsurf, LayerType& layer, const std::string& diskName, double binPhi, double binR) {
  ActsGeoBuilderUtils geoUtils;
  auto surfaces_disk = geoUtils.SelectActiveSurfaces_PhaseII(Kdtsurf, diskName);

  makeBinning(layer, Acts::SurfaceArrayNavigationPolicy::LayerType::Disc, binPhi, binR);

  layer.setSurfaces(surfaces_disk)
      .setLayerType(Acts::Experimental::LayerBlueprintNode::LayerType::Disc)
      .setEnvelope(Acts::ExtentEnvelope{{
          .z = {0.1 * Acts::UnitConstants::mm, 0.1 * Acts::UnitConstants::mm},
          .r = {1 * Acts::UnitConstants::mm, 1 * Acts::UnitConstants::mm},
      }})
      .setTransform(Acts::Transform3::Identity())
      .setUseCenterOfGravity(false, false, true);
}

//std::unique_ptr<Acts::TrackingGeometry> TrackerGeomBuilderPh2WithActsESProducer::produce(const ACTSTrackerGeometryRecord& iRecord) { 
std::shared_ptr<TrackingGeometryWithDetEls> TrackerGeomBuilderPh2WithActsESProducer::produce(const ACTSTrackerGeometryRecord& iRecord) { 

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

    // DEBUG pront the normal to the surface:
    // std::cout << "[CMSSW] Normal to module " << det->geographicalId().rawId() << ": " << cmssw_surf.normalVector() << std::endl;

    // ===== Define the transformation (i.e. Rotation and Translation) of the CMSSW surface =====
    auto& pos = cmssw_surf.position();
    auto& rot = cmssw_surf.rotation();
    Acts::Transform3 t = Acts::Transform3::Identity();
    Acts::RotationMatrix3 R;

    // ===== Get the Alignment INFO ===== (not needed as the aligment info are already stored in the transformation)
    DetId ID = det->geographicalId();
    // CLHEP::HepRotation Align_rot;
    // CLHEP::Hep3Vector Align_pos;
    // auto it = detID_to_alignInfo.find(ID.rawId());
    // if( it != detID_to_alignInfo.end()){
    //   auto this_tr = it->second;
    //   Align_rot = this_tr.getRotation();
    //   Align_pos = this_tr.getTranslation();
    //   DEBUG: compare rot e trs from surface and from alignm info
    //   std::cout << "-----------------\n"
    //             << "Module ID: " <<  ID.rawId() << "\n"
    //             << "Rotation and Translation from surface -> Rotation: " << rot << "; translation: " << pos << "\n"
    //             << "Rotation and Translation from alignment -> Rotation: " << Align_rot << "; translation: " << Align_pos << "\n" << std::endl;
    // }
    // else {
    //   std::cout << "[ERROR] No alignment info found for module " << ID.rawId() << std::endl;
    // }
    R << rot.xx(), rot.yx(), rot.zx(), rot.xy(), rot.yy(), rot.zy(), rot.xz(), rot.yz(), rot.zz(); // OLD
    // R << rot.xx(), rot.xy(), rot.xz(), rot.yx(), rot.yy(), rot.yz(), rot.zx(), rot.zy(), rot.zz(); // Transpose (Overlap in r)

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
    /// NOTE: Using only subdetId(), it is not possible to fully distinguish
    /// between all Phase-2 tracker subdetectors. For that, TrackerTopology
    /// information would also be needed (not essential for this use case).
    if(ID.subdetId() == PixelSubdetector::PixelBarrel) {
      subDet = "TBPX";
    } 
    else if (ID.subdetId() == PixelSubdetector::PixelEndcap) {
      subDet = "TEPX/TFPX";
    }     
    else if (ID.subdetId() == StripSubdetector::TID) {
      subDet = "TEDD";
    }     
    else if (ID.subdetId() == StripSubdetector::TOB) {
      subDet = "TB2S/TBPS";
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

  auto makeDisk = [this, &Kdtsurfaces]( auto& layer, const std::string& diskName, double binPhi, double binR) {
    this->makeLayer(Kdtsurfaces, layer, diskName, binPhi, binR);
  };

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

      CMS.addCylinderContainer("SubDetectors", Acts::AxisDirection::AxisZ, [&](auto& det){
        det.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

        det.addCylinderContainer("NegEC", Acts::AxisDirection::AxisZ, [&](auto& Negec){
          Negec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

          Negec.addCylinderContainer("NegTEDD_1to2_NegTFPX_8", Acts::AxisDirection::AxisR, [&](auto& ec) {
            ec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
              .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

            ec.addCylinderContainer("NegTEDD_1to2", Acts::AxisDirection::AxisZ, [&](auto& tedd){
              tedd.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                  .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

              AddExtraLayer_noMat("ExtraNegTEDD3", false, &tedd, GenerateTranslation(0, 0, -1250*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(200*Acts::UnitConstants::mm, 1200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

              AddDiskLayer_and_Material(&tedd, "TEDD_D1_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tedd, "TEDD_D2_neg", makeDisk, 80, 4);

              AddExtraLayer_noMat("ExtraNegTEDD4", false, &tedd, GenerateTranslation(0, 0, -1600*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(200*Acts::UnitConstants::mm, 1200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

            });
              
            ec.addCylinderContainer("NegTFPX_8", Acts::AxisDirection::AxisZ, [&](auto& tfpx){
              tfpx.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                  .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

              AddExtraLayer_noMat("ExtraNegTFPX1", false, &tfpx, GenerateTranslation(0, 0, -1250*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));
              
              AddDiskLayer_and_Material(&tfpx, "TFPX_D8_neg", makeDisk, 80, 4);
              
              AddExtraLayer_noMat("ExtraNegTFPX2", false, &tfpx, GenerateTranslation(0, 0, -1500*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

            });

          });

          Negec.addCylinderContainer("NegTEDD_3to5_NegTEPX", Acts::AxisDirection::AxisR, [&](auto& ec) {
            ec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

            ec.addCylinderContainer("NegTEDD_3to5", Acts::AxisDirection::AxisZ, [&](auto& tedd){
              tedd.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

              AddExtraLayer_noMat("ExtraNegTEDD1", false, &tedd, GenerateTranslation(0, 0, -1780*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(300*Acts::UnitConstants::mm, 1200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

              AddDiskLayer_and_Material(&tedd, "TEDD_D3_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tedd, "TEDD_D4_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tedd, "TEDD_D5_neg", makeDisk, 80, 4);

              AddExtraLayer_noMat("ExtraNegTEDD2", false, &tedd, GenerateTranslation(0, 0, -2700*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(300*Acts::UnitConstants::mm, 1200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

            });

            ec.addCylinderContainer("NegTEPX", Acts::AxisDirection::AxisZ, [&](auto& tepx){
              tepx.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

              AddExtraLayer_noMat("ExtraNegTEPX1", false, &tepx, GenerateTranslation(0, 0, -1700*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(60*Acts::UnitConstants::mm, 260*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

              AddDiskLayer_and_Material(&tepx, "TEPX_D1_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tepx, "TEPX_D2_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tepx, "TEPX_D3_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tepx, "TEPX_D4_neg", makeDisk, 80, 4);

              AddExtraLayer_noMat("ExtraNegTEPX2", false, &tepx, GenerateTranslation(0, 0, -2700*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(60*Acts::UnitConstants::mm, 260*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));
            });

          });
        });
        
        det.addCylinderContainer("Barrel", Acts::AxisDirection::AxisR, [&](auto& barr) {
          barr.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
              .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);


          barr.addCylinderContainer("TBPX_TFPX_1to7", Acts::AxisDirection::AxisZ, [&](auto& pix){
            pix.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

            pix.addCylinderContainer("NegTFPX_1to7", Acts::AxisDirection::AxisZ, [&](auto& tfpx) {
              tfpx.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

              auto makeDisk = [this, &Kdtsurfaces](auto& layer, const std::string& diskName, double binPhi, double binR) {
                this->makeLayer(Kdtsurfaces, layer, diskName, binPhi, binR);
              };

              AddDiskLayer_and_Material(&tfpx, "TFPX_D1_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D2_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D3_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D4_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D5_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D6_neg", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D7_neg", makeDisk, 80, 4);
              
            });

            pix.addCylinderContainer("TBPX", Acts::AxisDirection::AxisR, [&](auto& tbpx) {
              tbpx.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                  .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

              AddCylinderLayer_and_Material(&tbpx, "TBPX_L1", Kdtsurfaces, 12, 8);
              AddCylinderLayer_and_Material(&tbpx, "TBPX_L2", Kdtsurfaces, 12, 8);
              AddCylinderLayer_and_Material(&tbpx, "TBPX_L3", Kdtsurfaces, 12, 8);
              AddCylinderLayer_and_Material(&tbpx, "TBPX_L4", Kdtsurfaces, 12, 8);
            });

            pix.addCylinderContainer("PosTFPX_1to7", Acts::AxisDirection::AxisZ, [&](auto& tfpx) {
              tfpx.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

              AddDiskLayer_and_Material(&tfpx, "TFPX_D1_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D2_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D3_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D4_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D5_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D6_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tfpx, "TFPX_D7_pos", makeDisk, 80, 4);
            });

          });

          barr.addCylinderContainer("TBPS", Acts::AxisDirection::AxisR, [&](auto& tbps){
            tbps.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

            AddExtraLayer("InterTBPS_1", true, &tbps, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(175*Acts::UnitConstants::mm, 190*Acts::UnitConstants::mm, 1204*Acts::UnitConstants::mm));

            tbps.addCylinderContainer("TBPS_L1", Acts::AxisDirection::AxisZ, [&](auto& l1){
              l1.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

              l1.addCylinderContainer("PosTilted_L1", Acts::AxisDirection::AxisZ, [&](auto& pT_L1){
                pT_L1.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D1", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D2", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D3", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D4", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D5", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D6", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D7", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D8", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D9", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D10", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D11", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D12", makeDisk, 18, 1);

                // AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1_All", makeDisk, 18, 1);

                // Hybdrid:
                AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1_1to3", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1_3to6", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D7", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D8", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D9", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D10", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D11", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&pT_L1, "TBPS_PosTilted_L1D12", makeDisk, 18, 1);

                // AddExtraLayer_noMat("ExtraPosTilted_L1", false, &pT_L1, GenerateTranslation(0, 0, 1225*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(213*Acts::UnitConstants::mm, 300*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));
              });

              l1.addCylinderContainer("Barrel_L1", Acts::AxisDirection::AxisR, [&](auto& b_L1){
                b_L1.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                    .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

                AddExtraLayer_noMat("Barrel_L1_InnerGap", true, &b_L1, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(206 * Acts::UnitConstants::mm, 207 * Acts::UnitConstants::mm, 145 * Acts::UnitConstants::mm));
                AddCylinderLayer_and_Material(&b_L1, "TBPS_L1", Kdtsurfaces, 12, 8);
                AddCylinderLayer_and_Material(&b_L1, "TBPS_L2", Kdtsurfaces, 12, 8);

              });

              l1.addCylinderContainer("NegTilted_L1", Acts::AxisDirection::AxisZ, [&](auto& nT_L1){
                nT_L1.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D1", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D2", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D3", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D4", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D5", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D6", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D7", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D8", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D9", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D10", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D11", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D12", makeDisk, 18, 1);

                //AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1_All", makeDisk, 18, 1);

                // Hybdrid:
                AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1_1to3", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1_3to6", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D7", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D8", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D9", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D10", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D11", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&nT_L1, "TBPS_NegTilted_L1D12", makeDisk, 18, 1);

                // AddExtraLayer_noMat("ExtraNegTilted_L1", false, &nT_L1, GenerateTranslation(0, 0, -1225*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(213*Acts::UnitConstants::mm, 300*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));
              });

            });

            AddExtraLayer("InterTBPS_2", true, &tbps, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(310*Acts::UnitConstants::mm, 330*Acts::UnitConstants::mm, 1204*Acts::UnitConstants::mm));

            tbps.addCylinderContainer("TBPS_L2", Acts::AxisDirection::AxisZ, [&](auto& l2){
              l2.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

              l2.addCylinderContainer("PosTilted_L2", Acts::AxisDirection::AxisZ, [&](auto& pT_L2){
                pT_L2.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D1", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D2", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D3", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D4", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D5", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D6", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D7", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D8", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D9", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D10", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D11", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2D12", makeDisk, 18, 1);

                //AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2_All", makeDisk, 18, 1);

                // Hybrid:
                AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2_1to4", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2_5to8", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&pT_L2, "TBPS_PosTilted_L2_9to12", makeDisk, 18, 1);

                AddExtraLayer_noMat("ExtraPosTilted_L2", false, &pT_L2, GenerateTranslation(0, 0, 1225*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(335*Acts::UnitConstants::mm, 415*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));
              });


              l2.addCylinderContainer("Barrel_L2", Acts::AxisDirection::AxisR, [&](auto& b_L2){
                b_L2.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                    .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

                AddExtraLayer_noMat("Barrel_L2_InnerGap", true, &b_L2, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(335 * Acts::UnitConstants::mm, 336 * Acts::UnitConstants::mm, 215 * Acts::UnitConstants::mm));
                AddCylinderLayer_and_Material(&b_L2, "TBPS_L3", Kdtsurfaces, 12, 8);
                AddCylinderLayer_and_Material(&b_L2, "TBPS_L4", Kdtsurfaces, 12, 8);
              });

              l2.addCylinderContainer("NegTilted_L2", Acts::AxisDirection::AxisZ, [&](auto& nT_L2){
                nT_L2.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D1", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D2", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D3", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D4", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D5", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D6", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D7", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D8", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D9", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D10", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D11", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2D12", makeDisk, 18, 1);

                // AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2_All", makeDisk, 18, 1);

                // Hybrid:
                AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2_1to4", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2_5to8", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&nT_L2, "TBPS_NegTilted_L2_9to12", makeDisk, 18, 1);

                AddExtraLayer_noMat("ExtraNegTilted_L2", false, &nT_L2, GenerateTranslation(0, 0, -1225*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(335*Acts::UnitConstants::mm, 415*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));
              });

            });

            AddExtraLayer("InterTBPS_3", true, &tbps, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(420*Acts::UnitConstants::mm, 480*Acts::UnitConstants::mm, 1204*Acts::UnitConstants::mm));

            tbps.addCylinderContainer("TBPS_L3", Acts::AxisDirection::AxisZ, [&](auto& l3){
              l3.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

              l3.addCylinderContainer("PosTilted_L3", Acts::AxisDirection::AxisZ, [&](auto& pT_L3){
                pT_L3.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D1", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D2", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D3", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D4", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D5", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D6", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D7", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D8", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D9", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D10", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D11", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3D12", makeDisk, 18, 1);

                // AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3_All", makeDisk, 18, 1);

                // Hybdrid:
                AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3_1to6", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&pT_L3, "TBPS_PosTilted_L3_6to12", makeDisk, 18, 1);
                AddExtraLayer_noMat("ExtraPosTilted_L3", false, &pT_L3, GenerateTranslation(0, 0, 1225*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(490*Acts::UnitConstants::mm, 566*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));
              
              });

              l3.addCylinderContainer("Barrel_L3", Acts::AxisDirection::AxisR, [&](auto& b_L3){
                b_L3.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                    .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);

                AddExtraLayer_noMat("Barrel_L3_InnerGap", true, &b_L3, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(490 * Acts::UnitConstants::mm, 491 * Acts::UnitConstants::mm, 300 * Acts::UnitConstants::mm));    
                AddCylinderLayer_and_Material(&b_L3, "TBPS_L5", Kdtsurfaces, 12, 8);
                AddCylinderLayer_and_Material(&b_L3, "TBPS_L6", Kdtsurfaces, 12, 8);
              });

              l3.addCylinderContainer("NegTilted_L3", Acts::AxisDirection::AxisZ, [&](auto& nT_L3){
                nT_L3.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D1", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D2", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D3", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D4", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D5", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D6", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D7", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D8", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D9", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D10", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D11", makeDisk, 18, 1);
                // AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3D12", makeDisk, 18, 1);

                //AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3_All", makeDisk, 18, 1);

                // Hybdrid:
                AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3_1to6", makeDisk, 18, 1);
                AddDiskLayer_and_Material(&nT_L3, "TBPS_NegTilted_L3_6to12", makeDisk, 18, 1);
                AddExtraLayer_noMat("ExtraNegTilted_L3", false, &nT_L3, GenerateTranslation(0, 0, -1225*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(490*Acts::UnitConstants::mm, 566*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));
            
              });

            });

            AddExtraLayer("InterTBPS_4", true, &tbps, GenerateTranslation(0, 0, 0), std::make_shared<Acts::CylinderVolumeBounds>(580*Acts::UnitConstants::mm, 590*Acts::UnitConstants::mm, 1204*Acts::UnitConstants::mm));

          });


          barr.addCylinderContainer("TB2S", Acts::AxisDirection::AxisR, [&](auto& tb2s){
            tb2s.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap)
                .setResizeStrategy(Acts::VolumeResizeStrategy::Gap);
            
            AddCylinderLayer_and_Material(&tb2s, "TB2S_L1", Kdtsurfaces, 12, 8);
            AddCylinderLayer_and_Material(&tb2s, "TB2S_L2", Kdtsurfaces, 12, 8);
            AddCylinderLayer_and_Material(&tb2s, "TB2S_L3", Kdtsurfaces, 12, 8);
            AddCylinderLayer_and_Material(&tb2s, "TB2S_L4", Kdtsurfaces, 12, 8);
            AddCylinderLayer_and_Material(&tb2s, "TB2S_L5", Kdtsurfaces, 12, 8);
            AddCylinderLayer_and_Material(&tb2s, "TB2S_L6", Kdtsurfaces, 12, 8);
          });

          
        });

        det.addCylinderContainer("PosEC", Acts::AxisDirection::AxisZ, [&](auto& Posec){
          Posec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);
          
          Posec.addCylinderContainer("PosTEDD_3to5_PosTEPX", Acts::AxisDirection::AxisR, [&](auto& ec) {
            ec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

            ec.addCylinderContainer("PosTEDD_3to5", Acts::AxisDirection::AxisZ, [&](auto& tedd){
              tedd.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

              AddExtraLayer_noMat("ExtraPosTEDD1", false, &tedd, GenerateTranslation(0, 0, 1780*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(300*Acts::UnitConstants::mm, 1200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

              AddDiskLayer_and_Material(&tedd, "TEDD_D3_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tedd, "TEDD_D4_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tedd, "TEDD_D5_pos", makeDisk, 80, 4);

              AddExtraLayer_noMat("ExtraPosTEDD2", false, &tedd, GenerateTranslation(0, 0, 2700*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(300*Acts::UnitConstants::mm, 1200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

            });

            ec.addCylinderContainer("PosTEPX", Acts::AxisDirection::AxisZ, [&](auto& tepx){
              tepx.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

              AddExtraLayer_noMat("ExtraPosTEPX1", false, &tepx, GenerateTranslation(0, 0, 1700*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(60*Acts::UnitConstants::mm, 260*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

              AddDiskLayer_and_Material(&tepx, "TEPX_D1_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tepx, "TEPX_D2_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tepx, "TEPX_D3_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tepx, "TEPX_D4_pos", makeDisk, 80, 4);

              AddExtraLayer_noMat("ExtraPosTEPX2", false, &tepx, GenerateTranslation(0, 0, 2700*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(60*Acts::UnitConstants::mm, 260*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

            });

          });

          Posec.addCylinderContainer("PosTEDD_1to2_PosTFPX_8", Acts::AxisDirection::AxisR, [&](auto& ec) {
            ec.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

            ec.addCylinderContainer("PosTEDD_1to2", Acts::AxisDirection::AxisZ, [&](auto& tedd){
              tedd.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

              AddExtraLayer_noMat("ExtraPosTEDD3", false, &tedd, GenerateTranslation(0, 0, 1250*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(200*Acts::UnitConstants::mm, 1200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

              AddDiskLayer_and_Material(&tedd, "TEDD_D1_pos", makeDisk, 80, 4);
              AddDiskLayer_and_Material(&tedd, "TEDD_D2_pos", makeDisk, 80, 4);

              AddExtraLayer_noMat("ExtraPosTEDD4", false, &tedd, GenerateTranslation(0, 0, 1600*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(200*Acts::UnitConstants::mm, 1200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

            });

            ec.addCylinderContainer("PosTFPX_8", Acts::AxisDirection::AxisZ, [&](auto& tfpx){
              tfpx.setAttachmentStrategy(Acts::VolumeAttachmentStrategy::Gap);

              AddExtraLayer_noMat("ExtraPosTFPX1", false, &tfpx, GenerateTranslation(0, 0, 1250*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));
              AddDiskLayer_and_Material(&tfpx, "TFPX_D8_pos", makeDisk, 80, 4);
              AddExtraLayer_noMat("ExtraPosTFPX2", false, &tfpx, GenerateTranslation(0, 0, 1500*Acts::UnitConstants::mm), std::make_shared<Acts::CylinderVolumeBounds>(24*Acts::UnitConstants::mm, 200*Acts::UnitConstants::mm, 1*Acts::UnitConstants::mm));

            });
          });

        });

      });
      
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
    Acts::JsonMaterialDecorator jsonMatDec(dec_cfg, materialFile_, Acts::Logging::Level::INFO);

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

DEFINE_FWK_EVENTSETUP_MODULE(TrackerGeomBuilderPh2WithActsESProducer);