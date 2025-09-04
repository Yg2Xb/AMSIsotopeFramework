/***********************************************************
 *  File: TrackerCut.h
 *
 *  Modern C++ header file for AMS Tracker detector cuts.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/

#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include "tracker_var.h"
#include "Tool.h"

// 前向声明
class selectdata;

namespace AMS_Iso {

// Tracker状态信息
struct TrackerStatus {
    std::array<bool, Tracker::LAYER_COUNT> hasHit;     // 是否有击中
    std::array<bool, Tracker::LAYER_COUNT> hasXHit;    // 是否有X方向击中
    int innerLayerHits;                                // 内层击中数
    double rigidity;                                   // 刚度值
    double L38InnerAveQ;

    TrackerStatus() 
        : hasHit{}
        , hasXHit{}
        , innerLayerHits(0)
        , rigidity(0.0) 
        , L38InnerAveQ(0.0) 
    {}
};

class TrackerCut {
public:
    explicit TrackerCut(selectdata* event = nullptr);
    ~TrackerCut() = default;

    // 禁止拷贝和赋值
    TrackerCut(const TrackerCut&) = delete;
    TrackerCut& operator=(const TrackerCut&) = delete;

    // 基本属性访问
    double getRigidity(int algorithm = Tracker::Algorithm::DEFAULT,
                      int alignment = Tracker::Alignment::DEFAULT,
                      int span = Tracker::Span::DEFAULT) const;

    // 基本切割
    CutResult<4> cutBasicAndFiducial(bool isISS) const;
    CutResult<1> cutPhysTrigger(bool isISS) const;

    // 电荷相关切割
    CutResult<5> cutL1Unbiased(int charge, bool isISS = true,
                        bool forEfficiency = false, bool forBackground = false, float coe = 1.) const;
    CutResult<5> cutL1Norm(int charge, bool isISS = true, float coe = 1.) const;
    CutResult<1> cutUTOFQ(int charge, bool isISS = true,
                         bool forEfficiency = false, bool forBackground = false, float coe = 1.) const;
    CutResult<3> cutInnerQ(int charge, bool isISS = true,
                          bool forEfficiency = false, bool forBackground = false, float coe = 1.) const;
    CutResult<4> cutInnerTracker(int charge = 0, bool isISS = true,
                                bool forEfficiency = false, bool forBackground = false) const;
    CutResult<3> cutBackground(int charge = 0, bool isISS = true,
                             bool forEfficiency = false, bool forBackground = false) const;

    // 综合切割
    CutResult<9> cutTracker(int charge, bool isISS = true) const;
    CutResult<2> cutUnphysical(int charge, bool isISS = true) const;
    CutResult<2> getDenominatorL1PickUp(int charge, bool isISS) const;
    CutResult<6> chargeTempFitCut(int charge, bool isISS, bool isNormalL1 = false) const;

    // 辅助函数
    double getRadius(bool isUnphysical, int layer) const;
    bool isInFiducial(bool isUnphysical, int layer) const;
    int getSecondaryHitCount(int direction) const;  // 0: XY hit, 1: Y hit

private:
    selectdata* event_;        // 数据指针（不负责内存管理）
    TrackerStatus status_;     // Tracker状态信息
    bool isISS_;              // 是否ISS数据

    // 内部辅助函数
    void initializeStatus();
    bool validateLayer(int layer) const;
    std::array<double, 2> getLayerPosition(bool isUnphysical, int layer) const;
    bool checkFiducialCut(int layer, 
                         const std::array<double, 2>& position,
                         double radius) const;
};

} // namespace AMS_Iso