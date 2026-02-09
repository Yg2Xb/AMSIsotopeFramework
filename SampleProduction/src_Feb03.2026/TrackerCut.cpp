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
    if (isUnphysical) return {event_->tk_pos1s[layer][0], event_->tk_pos1s[layer][1]};
    return {event_->tk_pos[layer][0], event_->tk_pos[layer][1]};
}

bool TrackerCut::checkFiducialCut(int layer, const std::array<double, 2>& position, double radius) const {
    return radius < Tracker::FiducialCuts::R_POS[layer] && std::abs(position[1]) < Tracker::FiducialCuts::Y_POS[layer];
}

// ===================== Cuts Implementation (Using Cached Values) =====================

CutResult<4> TrackerCut::cutBasicAndFiducial(bool isISS, bool isUnbiased, bool ForInTrkEffNum) const {
    if (!event_) return CutResult<4>();

    std::array<bool, Tracker::LAYER_COUNT> layerFV;
    int innerFVCount = 0;
    for (int layer = 0; layer < Tracker::LAYER_COUNT; ++layer) {
        layerFV[layer] = isInFiducial(isUnbiased, layer);
        if (layerFV[layer] && layer > 0 && layer < 8) ++innerFVCount;
    }

    double tofBeta = event_->tof_betah; 
    std::array<bool, 4> cuts{
        isUnbiased || (event_->itrtrack >= 0 && event_->ibetah >= 0),
        isUnbiased || (event_->tof_btype < 10 && tofBeta > 0.4),
        innerFVCount >= 5 && (ForInTrkEffNum || layerFV[0]) && layerFV[1] && (layerFV[2] || layerFV[3]) && (layerFV[4] || layerFV[5]) && (layerFV[6] || layerFV[7]),
        true
    };
    return CutResult<4>(cuts);
}

CutResult<5> TrackerCut::cutL1Unbiased(int charge, bool isISS, 
                                bool forEfficiency, bool forBackground, float coe) const {
    if (!event_) return CutResult<5>();

    double l1Q = status_.L1Q_Unbiased; // CACHED
    int qStatus = status_.L1QStatus_Unbiased; // CACHED

    bool hasUnbXYSignal = (event_->tk_exqln[0][0][0] > 0) && (event_->tk_exqln[0][0][1] > 0);

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
    double l1ChiY = (status_.innerLayerHits - 2) * 
                event_->tk_chis1[Tracker::Algorithm::DEFAULT]
                                [Tracker::Alignment::DEFAULT]
                                [Tracker::Span::INNER_L1]  // 包含L1
                                [Tracker::Direction::Y] - 
                    (status_.innerLayerHits - 3) * 
                event_->tk_chis1[Tracker::Algorithm::DEFAULT]
                                [Tracker::Alignment::DEFAULT]
                                [Tracker::Span::INNER]     
                                [Tracker::Direction::Y];

    double l1Q = status_.L1Q_Normal; // CACHED
    int qStatus = status_.L1QStatus_Normal; // CACHED
    
    double low = 0.46 + (charge - 3) * 0.16;
    double up = (charge <= 5) ? 0.65 : 0.65+(charge-5)*0.03;
    
    std::array<bool, 5> cuts{
        (l1Q < charge + coe*up),
        (l1Q > charge - coe*low),
        qStatus == 0, 
        hasL1XY,
        true //l1ChiY < 10.0 Jan.25.2026 remove all L1 and L1Inner Chi2 as 
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
        !forEfficiency || cutBasicAndFiducial(isISS, false, true).total,
        true
    };
    return CutResult<4>(cuts);
}

CutResult<3> TrackerCut::cutBackground(int charge, bool isISS,
                                     bool forEfficiency, bool forBackground) const {
    if (!event_) return CutResult<3>();

    bool single = event_->ntrack == 1;
    bool low2nd = event_->betah2r < 0.5;
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
    //int ptrig = isISS ? event_->physbpatt2 : 0x3EL;
    int ptrig = isISS ? event_->physbpatt2 : event_->physbpatt1; //2025.12.07, apply mc phystrigger
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
        phys.total && basic.total && inTrk.total && inQ.total && l1n.total  && tof.total && bg.total
    }, false);
}

bool TrackerCut::AccUndepCut(int charge, bool isISS, bool forBackground, double coe) const {
    if (!event_) return false;
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    auto inQ = cutInnerQ(charge, isISS, false, forBackground, coe);
    auto tof = cutUTOFQ(charge, isISS, false, forBackground, coe);
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

CutResult<2> TrackerCut::TwoAccTrackerCut(int charge, bool isISS, bool forBackground, double coe) const {
    if (!event_) return CutResult<2>();
    bool accunb = AccUndepCut(charge, isISS, forBackground, coe);
    return CutResult<2>({
        accunb && cutL1Unbiased(charge, isISS, false, false, coe).total,
        accunb && cutL1Norm(charge, isISS, coe).total 
    }, false);
}

CutResult<12> TrackerCut::chargeTempCut(int zsrc, int fragZ, bool isISS, bool forBackground) const {
    std::array<bool, 12> cuts; cuts.fill(false);
    if (!event_) return CutResult<12>(cuts, false);
    
    // 1. Base
    if (!cutPhysTrigger(isISS).total) return CutResult<12>(cuts, false);
    if (!cutBasicAndFiducial(isISS).total) return CutResult<12>(cuts, false);
    if (!cutInnerTracker(zsrc, isISS).total) return CutResult<12>(cuts, false);

    // 2. Variables (Cached)
    auto bg = cutBackground(zsrc, isISS, false, forBackground);
    double innerQ = status_.innerQ;
    double L38 = status_.L38InnerAveQ;

    // 3. L1 Selections
    auto l1n = cutL1Norm(zsrc, isISS);
    auto l1u = cutL1Unbiased(zsrc, isISS);
    bool L1N_Qual = l1n.details[2] && l1n.details[3] && l1n.details[4] ;
    bool L1U_Qual = l1u.details[2] && l1u.details[3];

    // 4. Strict
    double bkg_coe = 0.6;
    if(zsrc==2) bkg_coe = 0.1;

    auto s_InQ  = cutInnerQ(zsrc, isISS, false, false, bkg_coe);
    auto s_TOF  = cutUTOFQ(zsrc, isISS, false, false, bkg_coe);
    auto s_L1N  = cutL1Norm(zsrc, isISS, bkg_coe);
    auto s_L1U  = cutL1Unbiased(zsrc, isISS, false, false, bkg_coe);

    // Indices: 0 L1Sig any, 1 L1Sig PassFullSel for BL1, 2 L1Temp, 3 L2Temp, 4 InnerTemp
    // --- Group 1: L1 Signal Study ---
    bool rms_ok = cutInnerQ(zsrc, isISS).details[1]; 

    if (bg.details[0]) {
        
        if (rms_ok || event_->tk_qrms[1] < 0.55) { cuts[0] = L1U_Qual; cuts[1] = L1N_Qual; }  // Any
        
        double innerLowLimit[7] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        double innerUpLimit[7] =  {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
        
        if (rms_ok && innerQ > zsrc - innerLowLimit[zsrc-2] && innerQ < zsrc + innerUpLimit[zsrc-2]) { cuts[2] = L1U_Qual; cuts[3] = L1N_Qual; } // Pass Loose TrackerCut except L1 
        
        if (AccUndepCut(zsrc, isISS, forBackground, 1.0)) { cuts[4] = L1U_Qual; cuts[5] = L1N_Qual; } // Pass Full TrackerCut except L1 
    }

    // --- Group 2: Templates ---
    // L1 Template
    if (s_InQ.total && s_TOF.details[1] && bg.details[0]) {
         cuts[6] = L1U_Qual; cuts[7] = L1N_Qual;
    }

    // L2 Template
    bool L38_ok = std::abs(L38 - zsrc) < 0.45 * bkg_coe;
    if (status_.hasL2XY && status_.hasL2QStatusGood && L38_ok && s_TOF.details[2] && bg.details[0]) {
        cuts[8] = s_L1U.total; 
        cuts[9] = s_L1N.total ; 
    }

    // InnerTemp: Pure X
    double LowerTOFmeanQ = (event_->tof_ql[2] + event_->tof_ql[3]) * 0.5; 
    //if(event_->tof_ql[2]==0 || event_->tof_ql[3]==0) LowerTOFmeanQ = LowerTOFmeanQ*2;
    if (rms_ok && cutUTOFQ(zsrc, isISS, false, false, 1.).details[1] && LowerTOFmeanQ > zsrc - 0.5 && LowerTOFmeanQ < zsrc + 0.5 && bg.details[0]) {
        cuts[10] = cutL1Norm(zsrc, isISS, 0.6).total ; 
        cuts[11] = cutL1Unbiased(zsrc, isISS, false, false, 0.6).total;
    }

    return CutResult<12>(cuts, false);
}

// Selector for Fragmentation
bool TrackerCut::FragSampleSel(int zsrc, int fragZ, int c, int selector, bool isISS, bool forBackground) const {
    bool pass = false;
    if (!event_) return pass;

    auto chargetempcut_beam = chargeTempCut(zsrc, fragZ, isISS, forBackground);
    auto chargetempcut_frag = chargeTempCut(fragZ,  fragZ, isISS, forBackground);

    bool inner_opts[5] = {
        chargetempcut_beam.details[0+c],         // 0: L1Beam, L2Any
        chargetempcut_beam.details[2+c],         // 1: L1Beam, L2Beam, loose cut
        chargetempcut_beam.details[4+c],         // 2: L1Beam, L2Beam, full cut
        chargetempcut_frag.details[2+c],         // 3: L1Beam, L2Frag, loose cut
        chargetempcut_frag.details[4+c],         // 4: L1Beam, L2Frag, full cut
    };

    if (selector < 0 || selector > 4) return pass; 
    if (!inner_opts[selector]) return pass;

    // 3. L1 Checks (Only calculated if Inner passed)
    double q_l1_u = status_.L1Q_Unbiased;
    double q_l1_n = status_.L1Q_Normal;
    bool L1N_Beam =  q_l1_n > (zsrc - 0.5) && q_l1_n < (zsrc + 0.5);
    bool L1U_Beam =  q_l1_u > (zsrc - 0.5) && q_l1_u < (zsrc + 0.5);
    
    pass = c == 0 ? L1U_Beam : L1N_Beam; 
    
    return pass;
}

CutResult<2> TrackerCut::getEfficiencyTrigger(int charge, bool isISS, bool forBackground) const {
    if (!event_) return CutResult<2>();
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    auto inQ = cutInnerQ(charge, isISS);
    auto l1u = cutL1Unbiased(charge, isISS);
    auto tof = cutUTOFQ(charge, isISS);
    auto bg = cutBackground(charge, isISS, false, forBackground);
    bool common = basic.total && inTrk.total && inQ.total && l1u.total && tof.total && bg.total;
    
    bool num = common && phys.total;
    bool den = common && !phys.total;
    return CutResult<2>({num , den});
}

CutResult<2> TrackerCut::getEfficiencynAcc(int charge, bool isISS, bool forBackground) const {
    if (!event_) return CutResult<2>();
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    auto inQ = cutInnerQ(charge, isISS);
    auto l1u = cutL1Unbiased(charge, isISS);
    auto tof = cutUTOFQ(charge, isISS);
    auto bg = cutBackground(charge, isISS, false, forBackground);

    bool tofbztrig = (event_->tofflag[1]==0);
    int NACC = (((event_->antipatt&1) != 0)+((event_->antipatt&2) != 0)+((event_->antipatt&4) != 0)+((event_->antipatt&8) != 0)+((event_->antipatt&16) != 0)+((event_->antipatt&32) != 0)+((event_->antipatt&64) != 0)+((event_->antipatt&128) != 0));
    bool trigflag_phynacc5 = ( (tofbztrig && NACC < 5) || ((event_->physbpatt1&0x3A) != 0) );
    bool trigflag_phynacc8 = ( (tofbztrig && NACC < 8) || ((event_->physbpatt1&0x3A) != 0) );

    bool passnAcc5 = isISS ? phys.total && NACC < 5 : trigflag_phynacc5; 
    bool passnAcc8 = isISS ? phys.total && NACC < 8 : trigflag_phynacc8; 
    
    bool common = basic.total && inTrk.total && inQ.total && l1u.total && tof.total && bg.total;
    return CutResult<2>({ 
        common && (isISS && event_->run > 1456503197) && passnAcc5, 
        common && (isISS && event_->run > 1456503197) && passnAcc8
    });
}

CutResult<2> TrackerCut::getEfficiencyL1QLowLimit(int charge, bool isISS, bool forBackground) const {
    if (!event_) return CutResult<2>();
    
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    auto inQ = cutInnerQ(charge, isISS, false, false, 1.0);
    auto tof = cutUTOFQ(charge, isISS, false, false, 1.0);
    //double UpperTOFmeanQ = (event_->tof_ql[0] + event_->tof_ql[1]) * 0.5; // Mul instead of div
    //if(event_->tof_ql[0]==0 || event_->tof_ql[1]==0) UpperTOFmeanQ = UpperTOFmeanQ*2;
    //bool tofq = UpperTOFmeanQ > charge - 0.6 && UpperTOFmeanQ < charge + 1.5;
    auto bg = cutBackground(charge, isISS, false, forBackground);
    auto l1u = cutL1Unbiased(charge, isISS);

    bool den = phys.total && basic.total && inTrk.total && inQ.total && tof.total && bg.total;
    bool num = den && l1u.details[1] && l1u.details[2] && l1u.details[3];

    return CutResult<2>({ num, den });
}

CutResult<2> TrackerCut::getEfficiencyL1PickUp(int charge, bool isISS, bool forBackground) const {
    if (!event_) return CutResult<2>();

    // 1. Base Selection (Common Denominator Base)
    auto tof = cutUTOFQ(charge, isISS, 1.0);
    auto l1u = cutL1Unbiased(charge, isISS);
    bool common = getEfficiencyL1QLowLimit(charge, isISS, forBackground).details[0] && tof.details[1] && l1u.total;

    //bool tightQ = (status_.L1Q_Unbiased > charge - 0.3 && status_.L1Q_Unbiased < charge + 0.3);
    bool den = common;// && tightQ;

    auto l1n = cutL1Norm(charge, isISS);
    bool num = den && l1n.details[1] && l1n.details[2] && l1n.details[3] && l1n.details[4] ;

    return CutResult<2>({ num, den });
}

CutResult<2> TrackerCut::getEfficiencyInnerTracking(int charge, bool isISS) const {
    if (!event_) return CutResult<2>();

    auto phys = cutPhysTrigger(isISS);
    auto basicForEff = cutBasicAndFiducial(isISS, true, false); //full unbiased fiducial volume cut
    bool extraBasicCuts = event_->itrdtracks >= 0 && event_->ibetahs >= 0 && event_->betahs > 0.4;
    //auto l1u = cutL1Unbiased(charge, isISS, 1.0); 
    bool externalL1 = event_->tk_l1qxy[0] > 0 && event_->tk_l1qxy[1] > 0; //l1 xy signal
    //extra basic cuts
    //unbiased tofq
    bool StricTOFQCut = event_->tof_qls[0] > charge - 0.4 && event_->tof_qls[0] < charge + 0.5 &&
                        event_->tof_qls[1] > charge - 0.4 && event_->tof_qls[1] < charge + 0.5 &&
                        event_->tof_qls[2] > charge - 0.4 && event_->tof_qls[2] < charge + 0.5 &&
                        event_->tof_qls[3] > charge - 0.4 && event_->tof_qls[3] < charge + 0.5;
    //bool StricUnbL1QCut = (status_.L1Q_Unbiased > charge - 0.6 && status_.L1Q_Unbiased < charge + 0.7);
    //unbiased external q
    bool StrictExternalQCut = event_->tk_exqvn[0][0] > charge - 0.6 && event_->tk_exqvn[0][0] < charge + 0.7;
    bool TOFRecCut = event_->tof_chiscs < 20 && event_->tof_chists < 20;//unbiased tof rec
    bool den = extraBasicCuts && phys.total && basicForEff.total && externalL1 && StricTOFQCut && StrictExternalQCut && TOFRecCut;
    
    auto inTrk = cutInnerTracker(charge, isISS, true);
    bool num = den && inTrk.total; // cutInnerTracker contain basic cut,haha!

    return CutResult<2>({ num, den });
}

CutResult<2> TrackerCut::getEfficiencyInnerTrackerQ(int charge, bool isISS, bool forBackground) const {
    if (!event_) return CutResult<2>();
    
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    auto l1n = cutL1Norm(charge, isISS, 1.0);
    auto tof = cutUTOFQ(charge, isISS);
    auto bg = cutBackground(charge, isISS, false, forBackground);
    bool StricL1QCut = (status_.L1Q_Normal > charge - 0.5 && status_.L1Q_Normal < charge + 0.5);
    double UpperTOFmeanQ = (event_->tof_ql[0] + event_->tof_ql[1]) * 0.5;
    double LowerTOFmeanQ = (event_->tof_ql[2] + event_->tof_ql[3]) * 0.5;
    bool StricTOFQCut = UpperTOFmeanQ > charge - 0.5 && UpperTOFmeanQ < charge + 0.5 &&
                        LowerTOFmeanQ > charge - 0.5 && LowerTOFmeanQ < charge + 0.5;

    bool den = phys.total && basic.total && inTrk.total && l1n.total && tof.details[0] && StricL1QCut && StricTOFQCut && bg.total;

    auto inQ = cutInnerQ(charge, isISS);
    bool num = den && inQ.total;

    return CutResult<2>({ num, den });
}

CutResult<2> TrackerCut::getEfficiencyUTOFQ(int charge, bool isISS, bool forBackground) const {
    if (!event_) return CutResult<2>();
    
    auto phys = cutPhysTrigger(isISS);
    auto basic = cutBasicAndFiducial(isISS);
    auto inTrk = cutInnerTracker(charge, isISS);
    auto inQ = cutInnerQ(charge, isISS, false, false, 1.0);
    auto l1n = cutL1Norm(charge, isISS, 1.0);
    auto bg = cutBackground(charge, isISS, false, forBackground);
    bool StricL1QCut = (status_.L1Q_Normal > charge - 0.5 && status_.L1Q_Normal < charge + 0.5);

    bool den = phys.total && basic.total && inTrk.total && inQ.total && l1n.total && bg.total  && StricL1QCut;

    auto tof = cutUTOFQ(charge, isISS);
    bool num = den && tof.total;

    return CutResult<2>({ num, den });
}

CutResult<2> TrackerCut::getEfficiencyBkgReduction(int charge, bool isISS) const {
    if (!event_) return CutResult<2>();
    auto num = TwoAccTrackerCut(charge, isISS, false, 0.4);
    auto den = TwoAccTrackerCut(charge, isISS, true, 0.4);
    return CutResult<2>({ num.details[0], den.details[0] });
}

CutResult<2> TrackerCut::TwoAccTrackerCut_OneTrk(int charge, bool isISS) const {
    if (!event_) return CutResult<2>();
    bool accunb = AccUndepCut(charge, isISS) && cutBackground(charge, isISS).details[2];
    return CutResult<2>({
        accunb && cutL1Unbiased(charge, isISS, false, false, 1.0).total,
        accunb && cutL1Norm(charge, isISS, 1.0).total 
    }, false);
}

} // namespace AMS_Iso