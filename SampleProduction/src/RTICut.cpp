/***********************************************************
 *  File: RTICut.cpp
 *
 *  Modern C++ implementation file for AMS RTI (Run Time Information)
 *  selections and exposure calculations.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/
#include <cmath>
#include "RTICut.h"
#include "selectdata.h"

namespace AMS_Iso {

RTICut::RTICut(selectdata* event)
    : event_(event)
{}

double RTICut::getCutoffRigidity(int degree) const {
    if (!event_) return 0.0;
    
    // 计算数组索引：(degree / 5 - 5)
    int index = degree / 5 - 5;
    return std::abs(event_->mcutoffi[index][1]);
}

bool RTICut::isPhotonTriggerRun() const {
    return (event_->run >= PHOTON_START && event_->run < PHOTON_END);
}

CutResult<10> RTICut::cutRTI() const {
    if (!event_ || event_->irti != 0) {
        return CutResult<10>();
    }

    bool isPhotonRun = isPhotonTriggerRun();
    
    std::array<bool, 10> cuts{
        // 触发率检查
        event_->rtintrig / event_->rtinev > 0.98,
        
        // 粒子率检查
        [this, isPhotonRun]() {
            double ratio = event_->rtinpar / event_->rtintrig;
            double minThreshold = isPhotonRun ? 0.02 : 
                               (0.07 / 1600.0 * event_->rtintrig);
            return ratio > minThreshold && ratio < 0.25;
        }(),
        
        // 错误率检查
        event_->rtinerr >= 0 && 
        event_->rtinerr / event_->rtinev < 0.1,
        
        // 事例数检查
        event_->rtinpar > 0 && event_->rtinev < 1800,
        
        // RTI状态检查
        (event_->rtigood & 0x3F) == 0,
        
        // 存活时间比例检查
        event_->rtilf > (isPhotonRun ? 0.35 : 0.5),
        
        // 天顶角检查
        event_->zenith < 40,
        
        // SAA区域检查
        !event_->issaa,
        
        // 曝光检查
        event_->rtinexl[0][0] > 700 && event_->rtinexl[1][0] > 500,
        
        // 坏运行检查
        !event_->isbadrun
    };

    return CutResult<10>(cuts);
}

RTIExposure RTICut::calculateExposure() const {
    if (!event_) return RTIExposure();

    double efficiency = event_->rtinev / 
                      static_cast<double>(event_->rtinev + event_->rtinerr);
    
    return RTIExposure(
        event_->rtilf * efficiency,
        efficiency,
        event_->rtilf
    );
}

} // namespace AMS_Iso