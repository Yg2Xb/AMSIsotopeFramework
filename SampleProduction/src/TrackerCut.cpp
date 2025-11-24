/***********************************************************
 * File: TrackerCut.cpp
 *
 * Modern C++ implementation file for AMS Tracker detector cuts.
 * Optimized for caching and performance.
 *
 * History:
 * 20241029 - created by ZX.Yan
 * 20241123 - Optimized
 ***********************************************************/
#include <cmath>
#include <bitset>
#include "TrackerCut.h"
#include "selectdata.h"
#include "Tool.h"

namespace AMS_Iso {

TrackerCut::TrackerCut(selectdata* event)
    : event_(event)
    , status_()
    , isISS_(event ? event->isreal : false)
{
    if (event_) {
        initializeStatus();
    }
}

void TrackerCut::initializeStatus() {
    // 1. Hit Map & Secondary Hits Calculation
    status_.innerLayerHits = 0;
    status_.secondaryHitCountX = 0;
    status_.secondaryHitCountY = 0;
    
    // Assuming tk_hitb is bitset representing layers
    std::bitset<32> XYHits(event_->tk_hitb[0]); // 0: XY
    std::bitset<32> YHits(event_->tk_hitb[1]);  // 1: Y
      
    double L38InnerHits = 0;
    status_.L38InnerAveQ = 0.0;

    for (int layer = 0; layer < Tracker::LAYER_COUNT; ++layer) {
        // Cache hit status
        status_.hasXYHit[layer] = XYHits.test(layer);
        status_.hasYHit[layer] = YHits.test(layer);
        
        // Calculate Inner Hits (L2-L8, index 1-7)
        if (layer > 0 && layer < Tracker::LAYER_COUNT - 1 && status_.hasYHit[layer]) {
            ++status_.innerLayerHits;
            // Calculate L3-L8 Avg Q
            if(layer > 1) {
                status_.L38InnerAveQ += event_->tk_qln[Tracker::ChargeReco::DEFAULT][layer]
                                                        [Tracker::Direction::DEFAULT]; 
                L38InnerHits++;
            }
        }

        // Calculate Secondary Hits (bitset logic from getSecondaryHitCount moved here)
        // Assuming event_->betah2hb is also bitset for secondary hits
        if (layer >= 1 && layer < Tracker::LAYER_COUNT - 1) {
            if (std::bitset<32>(event_->betah2hb[0]).test(layer)) ++status_.secondaryHitCountX;
            if (std::bitset<32>(event_->betah2hb[1]).test(layer)) ++status_.secondaryHitCountY;
        }
    }

    // 2. Charge Values Caching
    status_.L38InnerAveQ = L38InnerHits > 0 ? status_.L38InnerAveQ/L38InnerHits : 0.0;
    
    status_.innerQ = event_->tk_qin[Tracker::ChargeReco::DEFAULT][Tracker::Direction::DEFAULT];
    status_.innerQRMS = event_->tk_qrmn[Tracker::ChargeReco::DEFAULT][Tracker::Direction::DEFAULT];
    
    status_.L1Q_Unbiased = event_->tk_exqln[Tracker::ChargeReco::DEFAULT][0][Tracker::Direction::DEFAULT];
    status_.L1Q_Normal   = event_->tk_qln   [Tracker::ChargeReco::DEFAULT][0][Tracker::Direction::DEFAULT];
    status_.L2Q   = event_->tk_qln   [Tracker::ChargeReco::DEFAULT][1][Tracker::Direction::DEFAULT];

    // 3. Charge Status Caching (Masked)
    constexpr int GOOD_CHARGE_MASK = 0x10013D;
    status_.L1QStatus_Unbiased = event_->tk_exqls[0] & GOOD_CHARGE_MASK;
    status_.L1QStatus_Normal   = event_->tk_qls[0]   & GOOD_CHARGE_MASK;

    // 4. L2 Status Caching
    status_.hasL2XY = status_.hasXYHit[1]; // L2 is index 1
    status_.hasL2QStatusGood = ((event_->tk_qls[1] & GOOD_CHARGE_MASK) == 0);
}

// Accessor implementations
double TrackerCut::getRigidity(int algorithm, int alignment, int span) const {
    if (!event_) return -999999.9;
    return event_->tk_rigidity1[algorithm][alignment][span];
}

int TrackerCut::getSecondaryHitCount(int direction) const {
    if (!event_ || direction < 0 || direction > 1) return -1;
    return (direction == 0) ? status_.secondaryHitCountX : status_.secondaryHitCountY;
}

// ... (getRadius, isInFiducial, validateLayer, getLayerPosition, checkFiducialCut: KEEP AS IS) ...
double TrackerCut::getRadius(bool isUnphysical, int layer) const {
    if (!validateLayer(layer)) return -1.0;
    auto position = getLayerPosition(isUnphysical, layer);
    return std::hypot(position[0], position[1]);
}
bool TrackerCut::isInFiducial(bool isUnphysical, int layer) const {
    if (!validateLayer(layer)) return false;
    auto position = getLayerPosition(isUnphysical, layer);
    return checkFiducialCut(layer, position, std::hypot(position[0], position[1]));
}
bool TrackerCut::validateLayer(int layer) const { return layer >= 0 && layer < Tracker::LAYER_COUNT && event_; }
std::array<double, 2> TrackerCut::getLayerPosition(bool isUnphysical, int layer) const {
    if (isUnphysical) return {event_->tk_pos[layer][0], event_->tk_pos[layer][1]};
    return {event_->tk_pos[layer][0], event_->tk_pos[layer][1]};
}
bool TrackerCut::checkFiducialCut(int layer, const std::array<double, 2>& position, double radius) const {
    return radius < Tracker::FiducialCuts::R_POS[layer] && std::abs(position[1]) < Tracker::FiducialCuts::Y_POS[layer];
}

// ===================== Cuts Implementation (Using Cached Values) =====================

CutResult<4> TrackerCut::cutBasicAndFiducial(bool isISS) const {
    if (!event_) return CutResult<4>();

    std::array<bool, Tracker::LAYER_COUNT> layerFV;
    int innerFVCount = 0;
    for (int layer = 0; layer < Tracker::LAYER_COUNT; ++layer) {
        layerFV[layer] = isInFiducial(false, layer);
        if (layerFV[layer] && layer > 0 && layer < 8) ++innerFVCount;
    }

    double tofBeta = event_->tof_betah; 
    std::array<bool, 4> cuts{
        event_->itrtrack >= 0 && event_->ibetah >= 0,
        event_->tof_btype < 10 && tofBeta > 0.4,
        innerFVCount >= 5 && layerFV[0] && layerFV[1] && (layerFV[2] || layerFV[3]) && (layerFV[4] || layerFV[5]) && (layerFV[6] || layerFV[7]),
        true
    };
    return CutResult<4>(cuts);
}

CutResult<5> TrackerCut::cutL1Unbiased(int charge, bool isISS, 
                                bool forEfficiency, bool forBackground, float coe) const {
    if (!event_) return CutResult<5>();

    double l1Q = status_.L1Q_Unbiased; // CACHED
    int qStatus = status_.L1QStatus_Unbiased; // CACHED
    
    bool hasUnbXYSignal = (event_->tk_exqln[0][0][0] > 1.5) && (event_->tk_exqln[0][0][1] > 1.5);

    double low = 0.46 + (charge - 3) * 0.16;
    double up = (charge <= 5 ) ? 0.65 : 0.65+(charge-5)*0.03;

    std::array<bool, 5> cuts{
        (l1Q < charge + coe*up),
        (l1Q > charge - coe*low),
        qStatus == 0,
        hasUnbXYSignal,
        true
    };
    return CutResult<5>(cuts, true);
}

CutResult<5> TrackerCut::cutL1Norm(int charge, bool isISS, float coe) const {
    if (!event_) return CutResult<5>();

    bool hasL1XY = status_.hasXYHit[0]; // CACHED
    
    // Chi2 calculation still needs array access unless cached
    double l1ChiY = (status_.innerLayerHits - 2) * event_->tk_chis1[0][0][2][1] - 
                    (status_.innerLayerHits - 3) * event_->tk_chis1[0][0][1][1];

    double l1Q = status_.L1Q_Normal; // CACHED
    int qStatus = status_.L1QStatus_Normal; // CACHED
    
    double low = 0.46 + (charge - 3) * 0.16;
    double up = (charge <= 5) ? 0.65 : 0.65+(charge-5)*0.03;
    
    std::array<bool, 5> cuts{
        (l1Q < charge + coe*up),
        (l1Q > charge - coe*low),
        qStatus == 0, 
        hasL1XY,
        l1ChiY < 10.0
    };
    return CutResult<5>(cuts, true);
}

CutResult<3> TrackerCut::cutUTOFQ(int charge, bool isISS,
                                 bool forEfficiency, bool forBackground, float coe) const {
    if (!event_) return CutResult<3>();
    double UpperTOFmeanQ = (event_->tof_ql[0] + event_->tof_ql[1]) * 0.5; // Mul instead of div
    if(event_->tof_ql[0]==0 || event_->tof_ql[1]==0) UpperTOFmeanQ = UpperTOFmeanQ*2;
    std::array<bool, 3> cuts{
        UpperTOFmeanQ > (charge - coe*0.6) && UpperTOFmeanQ < (charge + coe*1.5),
        UpperTOFmeanQ > (charge - coe*0.5) && UpperTOFmeanQ < (charge + coe*0.5),
        UpperTOFmeanQ > (charge - coe*0.7) && UpperTOFmeanQ < (charge + coe*0.7)
    };
    return CutResult<3>(cuts,false);
}

CutResult<3> TrackerCut::cutInnerQ(int charge, bool isISS,
                                  bool forEfficiency, bool forBackground, float coe) const {
    if (!event_) return CutResult<3>();
    double Q = status_.innerQ; // CACHED
    double RMS = status_.innerQRMS; // CACHED
    std::array<bool, 3> cuts{
        Q > (charge - coe*0.45) && Q < (charge + coe*0.45),
        RMS < 0.55,
        true
    };
    return CutResult<3>(cuts);
}

CutResult<4> TrackerCut::cutInnerTracker(int charge, bool isISS,
                                         bool forEfficiency, bool forBackground) const {
    if (!event_) return CutResult<4>();
    
    // Use CACHED map
    bool hasInnerHits = (status_.innerLayerHits >= 5) && 
                       status_.hasYHit[1] && 
                       (status_.hasYHit[2] || status_.hasYHit[3]) && 
                       (status_.hasYHit[4] || status_.hasYHit[5]) && 
                       (status_.hasYHit[6] || status_.hasYHit[7]);

    double innerNormChisY = event_->tk_chis1[Tracker::Algorithm::DEFAULT]
                                          [Tracker::Alignment::DEFAULT]
                                          [Tracker::Span::INNER][1];

    std::array<bool, 4> cuts{
        hasInnerHits,
        innerNormChisY < 10,
        !forEfficiency || cutBasicAndFiducial(isISS).total,
        true
    };
    return CutResult<4>(cuts);
}

CutResult<3> TrackerCut::cutBackground(int charge, bool isISS,
                                     bool forEfficiency, bool forBackground) const {
    if (!event_) return CutResult<3>();

    bool single = event_->ntrack == 1;
    bool low2nd = std::abs(event_->betah2r) < 0.5;
    int yH = status_.secondaryHitCountY; // CACHED
    int xyH = status_.secondaryHitCountX; // CACHED

    std::array<bool, 3> cuts{
        forBackground || (single || (yH < 5 || xyH < 3) || low2nd),  
        (single || (yH < 5 || xyH < 3) || low2nd),                  
        single                                                                
    };
    return CutResult<3>(cuts, false);
}

CutResult<1> TrackerCut::cutPhysTrigger(bool isISS) const {
    if (!event_) return CutResult<1>();
    int ptrig = isISS ? event_->physbpatt2 : 0x3EL;
    return CutResult<1>({ (ptrig & 0x3EL) != 0 });
}

CutResult<10> TrackerCut::cutTracker(int charge, bool isISS) const {
    if (!event_) return CutResult<10>();
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    auto inQ = cutInnerQ(charge, isISS);
    auto l1u = cutL1Unbiased(charge, isISS);
    auto l1n = cutL1Norm(charge, isISS);
    auto tof = cutUTOFQ(charge, isISS);
    auto bg = cutBackground(charge, isISS);

    return CutResult<10>({
        phys.total && basic.total && inTrk.total && inQ.total && l1u.total && tof.total && bg.total,
        phys.total, basic.total, inTrk.total, inQ.total, l1u.total, tof.total, bg.total, l1n.total, 
        phys.total && basic.total && inTrk.total && inQ.total && l1n.total && tof.total && bg.total
    }, false);
}

CutResult<2> TrackerCut::cutUnphysical(int charge, bool isISS) const {
    if (!event_) return CutResult<2>();
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    auto inQ = cutInnerQ(charge, isISS);
    auto l1u = cutL1Unbiased(charge, isISS);
    auto tof = cutUTOFQ(charge, isISS);
    auto bg = cutBackground(charge, isISS);

    return CutResult<2>({
        !phys.total && basic.total && inTrk.total && inQ.total && l1u.total && tof.total && bg.total,
        !phys.total && basic.total && inTrk.total && inQ.total && l1u.total && tof.total && bg.details[2]
    }, false);
}

CutResult<2> TrackerCut::getDenominatorL1PickUp(int charge, bool isISS) const {
    if (!event_) return CutResult<2>();
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    auto inQ = cutInnerQ(charge, isISS);
    auto l1u = cutL1Unbiased(charge, isISS, false, false, 0.4);
    auto tof = cutUTOFQ(charge, isISS);
    auto bg = cutBackground(charge, isISS);

    std::array<bool, 2> cuts{
        phys.total && basic.total && inTrk.total && inQ.total && l1u.total && tof.total && bg.details[0],
        phys.total && basic.total && inTrk.total && inQ.total && tof.total && bg.details[2]
    };
    return CutResult<2>(cuts, false);
}

bool TrackerCut::AccUndepCut(int charge, bool isISS, bool forBackground) const {
    if (!event_) return false;
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    auto inQ = cutInnerQ(charge, isISS);
    auto tof = cutUTOFQ(charge, isISS);
    auto bg = cutBackground(charge, isISS, false, forBackground);
    return phys.total && basic.total && inTrk.total && inQ.total && tof.total && bg.details[0];
}

bool TrackerCut::Q_L1_BkgIndependCut(int charge, bool isISS) const {
    if (!event_) return false;
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    return phys.total && basic.total && inTrk.total;
}

CutResult<2> TrackerCut::TwoAccTrackerCut(int charge, bool isISS, bool forBackground) const {
    if (!event_) return CutResult<2>();
    bool accunb = AccUndepCut(charge, isISS, forBackground);
    return CutResult<2>({
        accunb && cutL1Unbiased(charge, isISS).total,
        accunb && cutL1Norm(charge, isISS).total 
    }, false);
}

CutResult<14> TrackerCut::chargeTempCut(int charge, int fragZ, bool isISS, bool forBackground) const {
    std::array<bool, 14> cuts; cuts.fill(false);
    if (!event_) return CutResult<14>(cuts, false);
    
    // 1. Base
    if (!cutPhysTrigger(isISS).total) return CutResult<14>(cuts, false);
    if (!cutBasicAndFiducial(isISS).total) return CutResult<14>(cuts, false);
    if (!cutInnerTracker(charge, isISS).total) return CutResult<14>(cuts, false);

    // 2. Variables (Cached)
    auto bg = cutBackground(charge, isISS, false, forBackground);
    double innerQ = status_.innerQ;
    double L38 = status_.L38InnerAveQ;

    // 3. L1 Selections
    auto l1n = cutL1Norm(charge, isISS);
    auto l1u = cutL1Unbiased(charge, isISS);
    bool L1N_Qual = l1n.details[2] && l1n.details[3] && l1n.details[4];
    bool L1U_Qual = l1u.details[2] && l1u.details[3];

    // 4. Strict
    double bkg_coe = (charge <= 3) ? 0.4 : (charge==2 ? 0.1 : 0.8); // Fix charge==2 logic order
    if (charge==2) bkg_coe = 0.1; // Ensure override

    auto s_InQ  = cutInnerQ(charge, isISS, false, false, bkg_coe);
    auto s_TOF  = cutUTOFQ(charge, isISS, false, false, bkg_coe);
    auto s_L1N  = cutL1Norm(charge, isISS, bkg_coe);
    auto s_L1U  = cutL1Unbiased(charge, isISS, false, false, bkg_coe);

    // --- Group 1: L1 Signal Study ---
    double Q_Low  = charge - 0.55;
    double Q_High = charge + 0.45;
    // Logic check: cutInnerQ check RMS (details[1])
    bool rms_ok = cutInnerQ(charge, isISS).details[1]; 

    if (rms_ok && bg.details[0]) {
        cuts[1] = L1N_Qual; cuts[0] = L1U_Qual; // Any
        if (innerQ > Q_Low && innerQ < Q_High) { cuts[3] = L1N_Qual; cuts[2] = L1U_Qual; } // Pass
        if (innerQ < Q_Low) { cuts[5] = L1N_Qual; cuts[4] = L1U_Qual; } // Frag
    }

    // --- Group 2: Templates ---
    if (s_InQ.total && s_TOF.details[1] && bg.details[1]) {
        cuts[7] = L1N_Qual; cuts[6] = L1U_Qual;
    }

    // L2 Template (Cached vars)
    bool L38_ok = std::abs(L38 - charge) < 0.45 * bkg_coe;
    if (status_.hasL2XY && status_.hasL2QStatusGood && L38_ok && s_TOF.details[2] && bg.details[1]) {
        cuts[9] = s_L1N.total; cuts[8] = s_L1U.total;
    }

    // --- Group 3: Inner Signal Study ---
    // InnerSig: L1 selects X (Normal/Unb Q window check manually)
    if(rms_ok && bg.details[0]){
        double q_l1_u = status_.L1Q_Unbiased;
        double q_l1_n = status_.L1Q_Normal;
        cuts[11] = L1N_Qual && q_l1_n > (charge - 0.2) && q_l1_n < (charge + 0.4);
        cuts[10] = L1U_Qual && q_l1_u > (charge - 0.2) && q_l1_u < (charge + 0.4);
    }

    // InnerTemp: Pure X
    if (s_TOF.details[2] && bg.details[1]) {
        cuts[13] = s_L1N.total; cuts[12] = s_L1U.total;
    }

    return CutResult<14>(cuts, false);
}

// Selector for Fragmentation
std::array<bool, 2> TrackerCut::FragSampleSel(int charge, int fragZ, int selector, bool isISS, bool forBackground) const {
    std::array<bool, 2> pass{false, false};
    if (!event_) return pass;

    // 1. Common Pre-checks (Trigger, Geo, L1Geo, Bkg, InnerRMS)
    if (!cutPhysTrigger(isISS).total) return pass;
    if (!cutBasicAndFiducial(isISS).total) return pass;
    if (!cutInnerTracker(charge, isISS).total) return pass; // Inner Geo
    if (!cutBackground(charge, isISS, false, forBackground).details[0]) return pass; // Bkg Condition
    if (status_.innerQRMS >= 0.55) return pass; // RMS Cut

    // 2. Calculate Inner Status Candidates (Sacrifice CPU for Code Size)
    double q = status_.innerQ;
    double X_lo = charge - 0.55, X_hi = charge + 0.45;
    double Y_lo = fragZ  - 0.55, Y_hi = fragZ  + 0.45;

    // [0]Num:TargetY, [1]Den1:Any, [2]Den2:PassX, [3]Den3:Frag
    bool inner_opts[4] = {
        (q > Y_lo && q < Y_hi),         // 0: X -> Y
        true,                           // 1: X -> Any (RMS checked)
        (q > X_lo && q < X_hi),         // 2: X -> X
        (q < X_lo)          // 3: X -> Frag
    };

    if (!inner_opts[selector & 3]) return pass; // Safety mask & 3

    // 3. L1 Checks (Only calculated if Inner passed)
    auto l1n = cutL1Norm(charge, isISS);
    auto l1u = cutL1Unbiased(charge, isISS);
    double q_l1_u = status_.L1Q_Unbiased;
    double q_l1_n = status_.L1Q_Normal;
    bool L1N_X = l1n.details[2] && l1n.details[3] && l1n.details[4] && q_l1_n > (charge - 0.2) && q_l1_n < (charge + 0.4);
    bool L1U_X = l1u.details[2] && l1u.details[3] && q_l1_u > (charge - 0.2) && q_l1_u < (charge + 0.4);
    
    pass[0] = L1U_X; // Unbiased
    pass[1] = L1N_X; // Normal
    
    return pass;
}

} // namespace AMS_Iso