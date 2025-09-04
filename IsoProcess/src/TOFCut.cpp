/***********************************************************
 *  File: TOFCut.cpp
 *
 *  Modern C++ implementation file for AMS TOF detector cuts.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/
#include <cmath>
#include <algorithm>
#include "TOFCut.h"
#include "selectdata.h"

namespace AMS_Iso {

TOFCut::TOFCut(selectdata* event)
    : event_(event)
    , betaTOF_(0.0)
    , isISS_(event ? event->isreal : false)
{
    if (event_) {
        calculateCharges();
        betaTOF_ = event_->tof_betah;
    }
}

void TOFCut::calculateCharges() {
    // 计算上层电荷
    std::array<double, 2> upperQ{event_->tof_ql[0], event_->tof_ql[1]};
    charge_.upperCharge = Tools::calculateAverage(upperQ.data(), 2, 0);

    // 计算下层电荷
    std::array<double, 2> lowerQ{event_->tof_ql[2], event_->tof_ql[3]};
    charge_.lowerCharge = Tools::calculateAverage(lowerQ.data(), 2, 0);
}

double TOFCut::getBeta() const {
    return event_ ? betaTOF_ : -1.;
}

bool TOFCut::cutBasic() const {
    return betaTOF_ > 0.4 && event_->tof_btype < 10;
}

bool TOFCut::cutBasicCharge(int charge) const {
    return charge_.upperCharge > (charge - 0.6) && 
           charge_.upperCharge < (charge + 1.5);
}

CutResult<2> TOFCut::cutBetaQuality(int charge) const {
    std::array<bool, 2> cuts{
        true,
        event_->tof_chisc < 5 && event_->tof_chist < 10
    };

    return CutResult<2>(cuts);
}

bool TOFCut::isInBar(int layer) const {
    int barDigit = static_cast<int>(event_->tof_pass * std::pow(10, layer - 3)) % 10;
    return barDigit > 0;
}

bool TOFCut::isLayerEdge(int layer) const {
    if (!event_) return true;

    int barId = event_->tof_barid[layer];
    bool isInBarArea = isInBar(layer);
    
    if (layer == 2) {  // 梯形层
        return barId == TOFConstants::TRAP_EDGE_BARS[0] || 
               barId == TOFConstants::TRAP_EDGE_BARS[1] || 
               !isInBarArea;
    } else {
        return barId == TOFConstants::EDGE_BARS[0] || 
               barId == TOFConstants::EDGE_BARS[1] || 
               !isInBarArea;
    }
}

bool TOFCut::isTrapezoidEdge(int layer) const {
    if (!event_ || layer < 2 || layer > 3) return true;

    bool isInBarArea = isInBar(layer);
    int barId = event_->tof_barid[layer];
    int edgeIndex = (layer == 2) ? TOFConstants::TRAP_EDGE_BARS[1] : 
                                  TOFConstants::EDGE_BARS[1];
    double cutValue = (layer == 2) ? 60.0 : 25.0;
    
    bool isAtEdge = (barId == 0 || barId == edgeIndex || !isInBarArea);
    double position = event_->tof_pos[layer][0];
    
    return isAtEdge && std::abs(position) > cutValue;
}

CutResult<4> TOFCut::cutEdges() const {
    std::array<bool, 4> cuts;
    for (int i = 0; i < TOFConstants::LAYER_COUNT; ++i) {
        cuts[i] = !isLayerEdge(i);
    }
    return CutResult<4>(cuts);
}

CutResult<2> TOFCut::cutTrapezoidEdges() const {
    std::array<bool, 2> cuts{
        !isTrapezoidEdge(2),
        !isTrapezoidEdge(3)
    };
    return CutResult<2>(cuts);
}

CutResult<2> TOFCut::cutTOF(int charge, bool isISS) const {
    std::array<bool, 2> cuts{
        //cutTrapezoidEdges().total,
        isISS ? event_->tof_goodgeo[0] == 1 : cutTrapezoidEdges().total,
        cutBetaQuality(charge).total
    };
    return CutResult<2>(cuts);
}

CutResult<2> TOFCut::cutTOFExcludeLayer4(int charge, bool isISS) const {
    std::array<bool, 2> cuts{
        //cutEdges().total,
        isISS ? event_->tof_goodgeo[1] == 1 : cutEdges().total,
        cutBetaQuality(charge).total
    };
    return CutResult<2>(cuts);
}

bool TOFCut::cutNoGeometry(int charge) const {
    return cutBetaQuality(charge).total;
}

} // namespace AMS_Iso