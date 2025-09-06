/***********************************************************
 *  File: RTICut.h
 *
 *  Modern C++ header file for AMS RTI (Run Time Information) 
 *  selections and exposure calculations.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/

#pragma once

#include <array>
#include <algorithm>  // 为 std::all_of
#include <cstdint>   // 为 uint32_t
#include "Tool.h"

// 前向声明
class selectdata;

namespace AMS_Iso {

// RTI曝光时间结构
struct RTIExposure {
    double value{0.0};           // 曝光时间值
    double efficiency{1.0};      // 效率比例
    double liveFraction{1.0};    // 存活时间比例

    RTIExposure() = default;
    RTIExposure(double v, double e, double l) 
        : value(v), efficiency(e), liveFraction(l) {}
};

class RTICut {
public:
    explicit RTICut(selectdata* event = nullptr);
    ~RTICut() = default;

    RTICut(const RTICut&) = delete;
    RTICut& operator=(const RTICut&) = delete;

    // 主要功能函数
    double getCutoffRigidity(int degree = 30) const;  // 获取截断刚度
    CutResult<10> cutRTI() const;                     // 执行RTI切割
    RTIExposure calculateExposure() const;            // 计算曝光时间

private:
    selectdata* event_;    // 不负责内存管理
    
    // 内部辅助函数
    bool isPhotonTriggerRun() const;

    // 光子触发运行时间常量
    static constexpr uint32_t PHOTON_START = 1620025528;
    static constexpr uint32_t PHOTON_END = 1635856717;
};

} // namespace AMS_Iso