/***********************************************************
 *  File: TOFCut.h
 *
 *  Modern C++ header file for AMS TOF detector cuts.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/

#pragma once

#include <array>
#include <vector>
#include "Tool.h"

// 前向声明
class selectdata;

namespace AMS_Iso {

// TOF层结构常量
struct TOFConstants {
    static constexpr int LAYER_COUNT = 4;
    static constexpr int UPPER_LAYERS = 2;
    static constexpr int LOWER_LAYERS = 2;
    static constexpr std::array<int, 2> EDGE_BARS{{0, 7}};  // 非梯形层边缘条
    static constexpr std::array<int, 2> TRAP_EDGE_BARS{{0, 9}};  // 梯形层边缘条
};

// TOF质量识别结构
struct TOFChargeResult {
    double upperCharge;    // 上层电荷
    double lowerCharge;    // 下层电荷
    
    TOFChargeResult() : upperCharge(0), lowerCharge(0) {}
};

class TOFCut {
public:
    explicit TOFCut(selectdata* event = nullptr);
    ~TOFCut() = default;

    // 禁止拷贝和赋值
    TOFCut(const TOFCut&) = delete;
    TOFCut& operator=(const TOFCut&) = delete;

    // 基本属性访问
    double getBeta() const;
    TOFChargeResult getCharge() const;

    // 基本切割
    bool cutBasic() const;
    bool cutBasicCharge(int charge) const;
    CutResult<2> cutBetaQuality(int charge) const;

    // 边缘切割
    bool isLayerEdge(int layer) const;
    bool isTrapezoidEdge(int layer) const;
    CutResult<4> cutEdges() const;
    CutResult<2> cutTrapezoidEdges() const;

    // 综合切割
    CutResult<2> cutTOF(int charge, bool isISS) const;
    CutResult<2> cutTOFExcludeLayer4(int charge, bool isISS) const;
    bool cutNoGeometry(int charge) const;

private:
    selectdata* event_;    // 数据指针（不负责内存管理）
    double betaTOF_;       // TOF beta值
    TOFChargeResult charge_;  // 电荷测量结果
    bool isISS_;          // 是否ISS数据

    // 内部辅助函数
    void calculateCharges();
    bool isInBar(int layer) const;
    std::array<double, 2> getLayerCharges(const std::array<int, 2>& layers) const;
};

} // namespace AMS_Iso