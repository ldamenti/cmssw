#ifndef ActsTrackingMaterialProducer_h
#define ActsTrackingMaterialProducer_h
#include <cstdint>
#include <string>
#include <vector>
#include <fstream>
#include <G4VTouchable.hh>

#include "SimG4Core/Watcher/interface/SimProducer.h"
#include "SimG4Core/Notification/interface/Observer.h"

#include "G4LogicalVolume.hh"

#include "SimDataFormats/ValidationFormats/interface/MaterialAccountingTrack.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"

class BeginOfJob;
class EndOfJob;
class BeginOfEvent;
class BeginOfTrack;
class EndOfTrack;
class G4Step;

class G4StepPoint;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4TouchableHistory;
namespace edm {
  class ParameterSet;
}

class ActsTrackingMaterialProducer : public SimProducer,
                                 public Observer<const BeginOfJob*>,
                                 public Observer<const EndOfJob*>,
                                 public Observer<const BeginOfEvent*>,
                                 public Observer<const BeginOfTrack*>,
                                 public Observer<const G4Step*>,
                                 public Observer<const EndOfTrack*>,
                                 public Observer<const EndOfRun*> {
public:
  ActsTrackingMaterialProducer(const edm::ParameterSet&);
  ~ActsTrackingMaterialProducer() override;

private:
  void update(const BeginOfJob*) override;
  void update(const BeginOfEvent*) override;
  void update(const BeginOfTrack*) override;
  void update(const G4Step*) override;
  void update(const EndOfTrack*) override;
  void update(const EndOfJob*) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void update(const EndOfRun*) override;

  bool isSelected(const G4VTouchable* touch);
  bool isSelectedFast(const G4TouchableHistory* touch);

private:
  bool m_primaryTracks;
  std::vector<std::string> m_selectedNames;
  std::vector<const G4LogicalVolume*> m_selectedVolumes;
  std::string m_txtOutFile;
  double m_hgcalzfront;
  MaterialAccountingTrack m_track;
  const G4VPhysicalVolume* m_track_volume;
  std::string m_ActsOutputFileName;
  std::vector<MaterialAccountingTrack>* m_tracks;
  TFile* output_file_;
  TProfile* radLen_vs_eta_;
  bool isHGCal;
  bool isHFNose;
  static constexpr float innerHGCalEta = 2.4;
  static constexpr float outerHGCalEta = 2.0;
  static constexpr float innerHFnoseEta = 4.;
  static constexpr float outerHFnoseEta = 3.3;
  std::ofstream outVolumeZpositionTxt;

  // #########################################################
  TTree* tTree_ = nullptr;
  std::vector<float> mat_step_length;
  std::vector<float> mat_X0;
  std::vector<float> mat_L0;
  std::vector<float> mat_x;
  std::vector<float> mat_y;
  std::vector<float> mat_z;
  // Totals (summary) expected by ActsPlugins::RootMaterialTrackIo
  float t_X0 = 0.f;
  float t_L0 = 0.f;

  std::uint32_t event_id = 0;

  float v_x = 0.f, v_y = 0.f, v_z = 0.f;          // vertex position
  float v_px = 0.f, v_py = 0.f, v_pz = 0.f;       // momentum at vertex
  float v_phi = 0.f, v_eta = 0.f;

  bool have_vertex_ = false;
  std::vector<float> mat_dx, mat_dy, mat_dz;
  std::vector<float> mat_Z, mat_A, mat_rho;

  std::vector<float> mat_sx, mat_sy, mat_sz; // start (pre-step)
  std::vector<float> mat_ex, mat_ey, mat_ez; // end   (post-step)
  // #########################################################
  

};

#endif  // ActsTrackingMaterialProducer_h
