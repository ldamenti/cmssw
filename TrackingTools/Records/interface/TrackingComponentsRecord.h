#ifndef RecoTracker_Record_TrackingComponentsRecord_h
#define RecoTracker_Record_TrackingComponentsRecord_h

#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
//#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/ACTSTrackerGeometryRecord.h"

#include "FWCore/Utilities/interface/mplVector.h"

// ORIGINAL VERSION
// class TrackingComponentsRecord : public edm::eventsetup::DependentRecordImplementation<
//                                      TrackingComponentsRecord,
//                                      edm::mpl::Vector<IdealMagneticFieldRecord, GlobalTrackingGeometryRecord> > {};

// MODIFIED VERSION:
class TrackingComponentsRecord : public edm::eventsetup::DependentRecordImplementation<TrackingComponentsRecord,
                                                                        edm::mpl::Vector<IdealMagneticFieldRecord, 
                                                                                         GlobalTrackingGeometryRecord,
                                                                                         TrackerDigiGeometryRecord,
                                                                                         ACTSTrackerGeometryRecord> > {};

#endif