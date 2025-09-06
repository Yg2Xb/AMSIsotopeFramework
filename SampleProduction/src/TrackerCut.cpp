/***********************************************************
 *  File: TrackerCut.cpp
 *
 *  Modern C++ implementation file for AMS Tracker detector cuts.
 *
 *  History:
 *    20241029 - created by ZX.Yan
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
    status_.innerLayerHits = 0;
    status_.L38InnerAveQ = 0;
    double L38InnerHits = 0;

    // y : bending direction
    std::bitset<32> XYHits(event_->tk_hitb[0]); //0: Tracker x and y
    std::bitset<32> YHits(event_->tk_hitb[1]);  //1: Only y
     
    for (int layer = 0; layer < Tracker::LAYER_COUNT; ++layer) {
        status_.hasXYHit[layer] = XYHits.test(layer);
        status_.hasYHit[layer] = YHits.test(layer);
        
        if (layer > 0 && layer < Tracker::LAYER_COUNT - 1 && status_.hasYHit[layer]) {
            ++status_.innerLayerHits;
            if(layer > 1)
            {
                status_.L38InnerAveQ += event_->tk_qln[Tracker::ChargeReco::DEFAULT][layer]
                                                        [Tracker::Direction::DEFAULT]; 
                L38InnerHits++;
            }
        }
    }

    status_.L38InnerAveQ = L38InnerHits > 0 ? status_.L38InnerAveQ/L38InnerHits : 0.0;

    status_.rigidity = event_->tk_rigidity1
        [Tracker::Algorithm::DEFAULT]
        [Tracker::Alignment::DEFAULT]
        [Tracker::Span::DEFAULT];
}

double TrackerCut::getRigidity(int algorithm, int alignment, int span) const {
    if (!event_) return -1.0;
    return event_->tk_rigidity1[algorithm][alignment][span];
}

int TrackerCut::getSecondaryHitCount(int direction) const {
    if (!event_ || direction < 0 || direction > 1) return -1;
    
    int hitCount = 0;
    std::bitset<32> hits(event_->betah2hb[direction]);
    
    for (int layer = 1; layer < Tracker::LAYER_COUNT - 1; ++layer) {
        if (hits.test(layer)) {
            ++hitCount;
        }
    }
    return hitCount;
}

double TrackerCut::getRadius(bool isUnphysical, int layer) const {
    if (!validateLayer(layer)) return -1.0;
    
    auto position = getLayerPosition(isUnphysical, layer);
    return std::hypot(position[0], position[1]);
}

bool TrackerCut::isInFiducial(bool isUnphysical, int layer) const {
    if (!validateLayer(layer)) return false;

    auto position = getLayerPosition(isUnphysical, layer);
    double radius = std::hypot(position[0], position[1]);
    
    return checkFiducialCut(layer, position, radius);
}

// 内部辅助函数实现
bool TrackerCut::validateLayer(int layer) const {
    return layer >= 0 && layer < Tracker::LAYER_COUNT && event_;
}

std::array<double, 2> TrackerCut::getLayerPosition(bool isUnphysical, int layer) const {
    if (isUnphysical) {
        return {event_->tk_pos1s[layer][0], event_->tk_pos1s[layer][1]};
    }
    return {event_->tk_pos[layer][0], event_->tk_pos[layer][1]};
}

bool TrackerCut::checkFiducialCut(int layer, 
                                 const std::array<double, 2>& position,
                                 double radius) const {
    return radius < Tracker::FiducialCuts::R_POS[layer] && 
           std::abs(position[1]) < Tracker::FiducialCuts::Y_POS[layer];
}

// cut函数实现
CutResult<4> TrackerCut::cutBasicAndFiducial(bool isISS) const {
    if (!event_) return CutResult<4>();

    std::array<bool, Tracker::LAYER_COUNT> layerFV;
    int innerFVCount = 0;
    
    for (int layer = 0; layer < Tracker::LAYER_COUNT; ++layer) {
        layerFV[layer] = isInFiducial(false, layer);
        if (layerFV[layer] && layer > 0 && layer < 8) {
            ++innerFVCount;
        }
    }

    double tofBetaForBasic = event_->tof_betah; 

    std::array<bool, 4> cuts{
        //event_->nlevel1 > 0,  // Level-1触发器
        //ihep other using itrtracks, but i use itrtrack, diff choose method, see dst
        event_->itrtrack >= 0 && event_->ibetah >= 0,// && event_->trd_nhitk >= 0,  // TkInner matched TOF
        event_->tof_btype < 10 && tofBetaForBasic > 0.4,  // TOF基本cut
        innerFVCount >= 5 && layerFV[0] && layerFV[1] && (layerFV[2] || layerFV[3]) && (layerFV[4] || layerFV[5]) && (layerFV[6] || layerFV[7]),  // Fiducial
        true
    };

    return CutResult<4>(cuts);
}

CutResult<5> TrackerCut::cutL1Unbiased(int charge, bool isISS, 
                                bool forEfficiency, bool forBackground, float coe) const {
    if (!event_) return CutResult<5>();

    double l1UnbiasedQ = event_->tk_exqln[Tracker::ChargeReco::DEFAULT][0]
                                       [Tracker::Direction::DEFAULT];
    constexpr int GOOD_CHARGE_MASK = 0x10013D;
    int l1QStatus = event_->tk_exqls[0] bitand GOOD_CHARGE_MASK;  // 使用bitand
    
    bool hasXYSignal = (event_->tk_exqln[Tracker::ChargeReco::DEFAULT][0][0] > 0) && 
                      (event_->tk_exqln[Tracker::ChargeReco::DEFAULT][0][1] > 0);

    double upperLimit = (charge <= 5 ) ? 0.65 : 0.65+(charge-5)*0.03;
    double lowerLimit = 0.46 + (charge - 3) * 0.16;

    std::array<bool, 5> cuts{
        (!isISS || l1UnbiasedQ < charge + coe*upperLimit),
        (l1UnbiasedQ > charge - coe*lowerLimit),
        l1QStatus == 0,
        hasXYSignal,
        true
    };

    return CutResult<5>(cuts, true);
}

CutResult<5> TrackerCut::cutL1Norm(int charge, bool isISS, float coe) const {
    if (!event_) return CutResult<5>();

    // 检查第一层是否有XY方向的击中
    bool hasLayer1XY = status_.hasXYHit[0];
    
    // 计算L1层的Y方向卡方值
    double l1ChisqY = (status_.innerLayerHits + 1 - 3) * 
                      event_->tk_chis1[Tracker::Algorithm::DEFAULT][Tracker::Alignment::DEFAULT][Tracker::Span::INNER_L1][1] -
                      (status_.innerLayerHits - 3) * 
                      event_->tk_chis1[Tracker::Algorithm::DEFAULT][Tracker::Alignment::DEFAULT][Tracker::Span::INNER][1];

    // 检查L1层电荷重建质量
    constexpr int CHARGE_QUALITY_MASK = 0x10013D;
    int l1QStatus = event_->tk_qls[0] bitand CHARGE_QUALITY_MASK;
    double chargelowLimit = 0.46 + (charge - 3) * 0.16;
    double chargeupLimit = (charge <= 5) ? 0.65 : 0.65+(charge-5)*0.03;
    double L1Q = event_->tk_qln[Tracker::ChargeReco::DEFAULT][0][Tracker::Direction::DEFAULT];
    bool goodL1Charge = (L1Q > charge - chargelowLimit) && (!isISS || L1Q < charge + chargeupLimit) && (l1QStatus == 0);
    bool goodL1Charge_Strict = (L1Q > charge - 0.4*chargelowLimit) && (!isISS || L1Q < charge + 0.4*chargeupLimit) && (l1QStatus == 0);;

    std::array<bool, 5> cuts{
        (!isISS || L1Q < charge + coe*chargeupLimit),
        (L1Q > charge - coe*chargelowLimit),
        l1QStatus == 0, 
        hasLayer1XY,
        l1ChisqY < 10.0
    };

    return CutResult<5>(cuts, true);
}

CutResult<1> TrackerCut::cutUTOFQ(int charge, bool isISS,
                                 bool forEfficiency, bool forBackground, float coe) const {
    if (!event_) return CutResult<1>();

    std::array<double, 2> tofUpperQ{event_->tof_ql[0], event_->tof_ql[1]};
    double meanUpperQ = Tools::calculateAverage(tofUpperQ.data(), 2, 0);

    std::array<bool, 1> cuts{
        meanUpperQ > (charge - coe*0.6) && meanUpperQ < (charge + coe*1.5),
    };

    return CutResult<1>(cuts,true);
}

CutResult<3> TrackerCut::cutInnerQ(int charge, bool isISS,
                                  bool forEfficiency, bool forBackground, float coe) const {
    if (!event_) return CutResult<3>();

    double innerQ = event_->tk_qin[Tracker::ChargeReco::DEFAULT]
                                [Tracker::Direction::DEFAULT];
    double innerQRMS = event_->tk_qrmn[Tracker::ChargeReco::DEFAULT]
                                [Tracker::Direction::DEFAULT];

    std::array<bool, 3> cuts{
        innerQ > (charge - coe*0.45) && innerQ < (charge + coe*0.45),
        innerQRMS < 0.55,
        true
    };

    return CutResult<3>(cuts);
}

CutResult<4> TrackerCut::cutInnerTracker(int charge, bool isISS,
                                        bool forEfficiency, bool forBackground) const {
    if (!event_) return CutResult<4>();
    
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

    // 优化：预先计算常用条件
    bool singleTrack = event_->ntrack == 1;
    bool Low2ndRig = std::abs(event_->betah2r) < 0.5; // 25.6.12 add abs!!!
    int yHits = getSecondaryHitCount(1);
    int xyHits = getSecondaryHitCount(0);

    std::array<bool, 3> cuts{
        singleTrack || (yHits < 5 || xyHits < 3) || Low2ndRig,  // TWIKI背景cut
        true,                                                   // no背景cut
        singleTrack                                         // strict
    };

    return CutResult<3>(cuts, false);
}

CutResult<1> TrackerCut::cutPhysTrigger(bool isISS) const {
    if (!event_) return CutResult<1>();

    //int physbpattp=(isISS)?event_->physbpatt2:event_->physbpatt1;//Trigger ISS/Trigger MC
    int physbpattp=(isISS)?event_->physbpatt2:0x3EL;//Trigger ISS/no MC cut

    std::array<bool, 1> cuts{
        (physbpattp&0x3EL)!=0
    };

    return CutResult<1>(cuts);
}

CutResult<9> TrackerCut::cutTracker(int charge, bool isISS) const {
    if (!event_) return CutResult<9>();

    // 获取各个cut结果
    auto physTrig = cutPhysTrigger(isISS);
    auto basicFid = cutBasicAndFiducial(isISS);
    auto innerTrk = cutInnerTracker(charge, isISS);
    auto innerQ = cutInnerQ(charge, isISS);
    auto unbiasedL1Cut = cutL1Unbiased(charge, isISS);
    auto L1Cut = cutL1Norm(charge, isISS);
    auto utofQ = cutUTOFQ(charge, isISS);
    auto bg = cutBackground(charge, isISS);

    // 直接返回CutResult
    return CutResult<9>({
        // 总结果
        physTrig.total && basicFid.total && innerTrk.total && 
        innerQ.total && unbiasedL1Cut.total && utofQ.total && bg.total,
        
        // 各个独立cut结果
        physTrig.total,
        basicFid.total,
        innerTrk.total,
        innerQ.total,
        unbiasedL1Cut.total,
        utofQ.total,
        bg.total,
        
        // 无背景reduction
        physTrig.total && basicFid.total && innerTrk.total && 
        innerQ.total && unbiasedL1Cut.total && utofQ.total
    }, false);
}

CutResult<2> TrackerCut::cutUnphysical(int charge, bool isISS) const {
    if (!event_) return CutResult<2>();

    // 获取各个cut结果
    auto physTrig = cutPhysTrigger(isISS);
    auto basicFid = cutBasicAndFiducial(isISS);
    auto innerTrk = cutInnerTracker(charge, isISS);
    auto innerQ = cutInnerQ(charge, isISS);
    auto unbiasedL1Cut = cutL1Unbiased(charge, isISS);
    auto L1Cut = cutL1Norm(charge, isISS);
    auto utofQ = cutUTOFQ(charge, isISS);
    auto bg = cutBackground(charge, isISS);

    return CutResult<2>({
        // 总结果（注意物理触发取反）
        !(physTrig.total) && basicFid.total && innerTrk.total && 
        innerQ.total && unbiasedL1Cut.total && utofQ.total && bg.total,
    
        // 严格
        !physTrig.total && basicFid.total && innerTrk.total && 
        innerQ.total && unbiasedL1Cut.total && utofQ.total && bg.details[2]
    }, false);
}

CutResult<2> TrackerCut::getDenominatorL1PickUp(int charge, bool isISS) const {
    if (!event_) return CutResult<2>();

    // 获取各个cut结果
    auto physTrig = cutPhysTrigger(isISS);
    auto basicFid = cutBasicAndFiducial(isISS);
    auto innerTrk = cutInnerTracker(charge, isISS);
    auto innerQ = cutInnerQ(charge, isISS);
    auto unbiasedL1Cut = cutL1Unbiased(charge, isISS, false, false, 0.4);
    auto utofQ = cutUTOFQ(charge, isISS);
    auto bg = cutBackground(charge, isISS);

    std::array<bool, 2> cuts{
        // Wei Cut
        physTrig.total && basicFid.total && innerTrk.total && 
        innerQ.total && unbiasedL1Cut.total && utofQ.total && bg.details[0],
        
        // 
        physTrig.total && basicFid.total && innerTrk.total && 
        innerQ.total && utofQ.total && bg.details[2]
    };

    return CutResult<2>(cuts, false);
}

CutResult<6> TrackerCut::chargeTempFitCut(int charge, bool isISS, bool isNormalL1) const {
    if (!event_) return CutResult<6>();

    bool basicCuts = cutPhysTrigger(isISS).total && 
                    cutBasicAndFiducial(isISS).total && 
                    cutInnerTracker(charge, isISS).total;

    const std::vector<float> coes = {0.2, 0.4, 1.0};
    std::array<bool, 6> cuts;

    for (int i = 0; i < 3; ++i) {
        float coe = coes[i];
        
        // L1Q信号: 使用InnerQ cut
        cuts[i] = basicCuts && 
                  cutInnerQ(charge, isISS, false, false, coe).total && 
                  cutUTOFQ(charge, isISS, false, false, coe).total;

        // L2Q模板: 使用L38InnerAveQ cut
        bool innerL38Cut = std::abs(status_.L38InnerAveQ - charge) < coe*0.4; // 25.0605, 0.45 to 0.4
        bool l2QStatus = (event_->tk_qls[1]&0x10013D)==0;
        cuts[i+3] = basicCuts && innerL38Cut && 
                    cutUTOFQ(charge, isISS, false, false, coe).total && l2QStatus &&
                    (isNormalL1 ? 
                        cutL1Norm(charge, isISS, coe).total :
                        cutL1Unbiased(charge, isISS, false, false, coe).total);
    }

    return CutResult<6>(cuts, false);
}

} // namespace AMS_Iso