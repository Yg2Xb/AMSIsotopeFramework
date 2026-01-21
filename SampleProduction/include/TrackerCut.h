/***********************************************************
 * File: TrackerCut.h
 *
 * Modern C++ header file for AMS Tracker detector cuts.
 *
 * History:
 * 20241029 - created by ZX.Yan
 * 20241123 - Modified for optimization and caching
 ***********************************************************/

#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include "tracker_var.h"
#include "Tool.h"

// Forward declaration
class selectdata;

namespace AMS_Iso {

// Tracker Status Cache
struct TrackerStatus {
    // Hit Maps
    std::array<bool, Tracker::LAYER_COUNT> hasXYHit;    // Layer-wise XY hit
    std::array<bool, Tracker::LAYER_COUNT> hasYHit;     // Layer-wise Y hit
    int innerLayerHits;                                 // Count of inner hits
    
    // Secondary Hits (Cached for cutBackground)
    int secondaryHitCountX;
    int secondaryHitCountY;
    
    // L2 Status (Cached for Templates)
    bool hasL2XY;
    bool hasL2QStatusGood;

    // Charge Values (Cached)
    double innerQ;
    double innerQRMS;
    double L38InnerAveQ;
    double L1Q_Unbiased;
    double L1Q_Normal;
    double L2Q;
    
    // Charge Quality Status (Cached)
    int L1QStatus_Unbiased;
    int L1QStatus_Normal;

    // Rigidity
    double rigidity;

    TrackerStatus() 
        : hasXYHit{}, hasYHit{}, innerLayerHits(0)
        , secondaryHitCountX(0), secondaryHitCountY(0)
        , hasL2XY(false), hasL2QStatusGood(false)
        , innerQ(0.0), innerQRMS(0.0), L38InnerAveQ(0.0)
        , L1Q_Unbiased(0.0), L1Q_Normal(0.0)
        , L1QStatus_Unbiased(-1), L1QStatus_Normal(-1)
        , rigidity(0.0) 
    {}
};

class TrackerCut {
public:
    explicit TrackerCut(selectdata* event = nullptr);
    ~TrackerCut() = default;

    // Disable copy/assign
TrackerCut(const TrackerCut&) = delete;
TrackerCut& operator=(const TrackerCut&) = delete;

TrackerCut(TrackerCut&&) = default;
TrackerCut& operator=(TrackerCut&&) = default;


    // --- Accessors for Cached Variables (Optimization) ---
    double getInnerQ() const { return status_.innerQ; }
    double getInnerQRMS() const { return status_.innerQRMS; }
    double getL1Q_Unbiased() const { return status_.L1Q_Unbiased; }
    double getL1Q_Normal() const { return status_.L1Q_Normal; }
    double getL2Q() const { return status_.L2Q; }
    double getL38InnerAveQ() const { return status_.L38InnerAveQ; }
    // Basic property access (Legacy, if needed)
    double getRigidity(int algorithm = Tracker::Algorithm::DEFAULT,
                       int alignment = Tracker::Alignment::DEFAULT,
                       int span = Tracker::Span::DEFAULT) const;

    // --- Basic Cuts ---
    CutResult<4> cutBasicAndFiducial(bool isISS, bool isUnbiased = false, bool forInTrkEffNum = false) const;
    CutResult<1> cutPhysTrigger(bool isISS) const;

    // --- Charge Cuts ---
    CutResult<5> cutL1Unbiased(int charge, bool isISS = true,
                         bool forEfficiency = false, bool forBackground = false, float coe = 1.) const;
    CutResult<5> cutL1Norm(int charge, bool isISS = true, float coe = 1.) const;
    CutResult<3> cutUTOFQ(int charge, bool isISS = true,
                          bool forEfficiency = false, bool forBackground = false, float coe = 1.) const;
    CutResult<3> cutInnerQ(int charge, bool isISS = true,
                           bool forEfficiency = false, bool forBackground = false, float coe = 1.) const;
    CutResult<4> cutInnerTracker(int charge = 0, bool isISS = true,
                                bool forEfficiency = false, bool forBackground = false) const;
    CutResult<3> cutBackground(int charge = 0, bool isISS = true,
                             bool forEfficiency = false, bool forBackground = false) const;

    // --- Complex Cuts ---
    CutResult<10> cutTracker(int charge, bool isISS = true) const;
    CutResult<2> TwoAccTrackerCut(int charge, bool isISS, bool forBackground = false, double coe = 1.) const;
    
    CutResult<12> chargeTempCut(int charge, int fragZ, bool isISS, bool forBackground) const;

    bool AccUndepCut(int charge, bool isISS, bool forBackground = false, double coe = 1.) const;
    bool Q_L1_BkgIndependCut(int charge, bool isISS) const;
    
    bool FragSampleSel(int charge, int fragZ, int c, int selector, bool isISS, bool forBackground) const;
    
    //eff
    CutResult<2> getEfficiencyTrigger(int charge, bool isISS, bool forBackground = false) const;
    CutResult<2> getEfficiencynAcc(int charge, bool isISS, bool forBackground = false) const;
    CutResult<2> getEfficiencyL1QLowLimit(int charge, bool isISS, bool forBackground = false) const;
    CutResult<2> getEfficiencyL1PickUp(int charge, bool isISS, bool forBackground = false) const;
    CutResult<2> getEfficiencyInnerTracking(int charge, bool isISS) const;
    CutResult<2> getEfficiencyInnerTrackerQ(int charge, bool isISS, bool forBackground = false) const;
    CutResult<2> getEfficiencyUTOFQ(int charge, bool isISS, bool forBackground = false) const;
    CutResult<2> getEfficiencyBkgReduction(int charge, bool isISS) const;
    CutResult<2> TwoAccTrackerCut_OneTrk(int charge, bool isISS) const;

    // Helpers
    double getRadius(bool isUnphysical, int layer) const;
    bool isInFiducial(bool isUnphysical, int layer) const;
    int getSecondaryHitCount(int direction) const;  // 0: XY hit, 1: Y hit (Now returns cached value)

private:
    selectdata* event_;       
    TrackerStatus status_;     
    bool isISS_;               

    // Internal Helpers
    void initializeStatus();
    bool validateLayer(int layer) const;
    std::array<double, 2> getLayerPosition(bool isUnphysical, int layer) const;
    bool checkFiducialCut(int layer, 
                          const std::array<double, 2>& position,
                          double radius) const;
};

} // namespace AMS_Iso