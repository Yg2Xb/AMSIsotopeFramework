/***********************************************************
 *  File: tracker_var.h
 *
 *  Modern C++ header file for AMS Tracker constants.
 *
 *  History:
 *    20241029 - Modernized by ZX.Yan
 *    20231200 - Original version
 ***********************************************************/

#pragma once

#include <array>

namespace AMS_Iso {
namespace Tracker {

// 基本参数
inline constexpr int LAYER_COUNT = 9;

// 算法选择
namespace Algorithm {
    inline constexpr int Choutko = 0;          // Choutko
    inline constexpr int GBL = 1;              // Generalized Broken Lines
    inline constexpr int KALMAN = 2;           // Kalman Filter
    
    inline constexpr int DEFAULT = GBL;        // 默认使用GBL
}

// 对准方式
namespace Alignment {
    inline constexpr int INNER_PG = 0;         // inner/PG
    inline constexpr int EXTERNAL_CIEMAT = 1;  // external layers/PG + CIEMAT
    inline constexpr int V6_ALL = 2;           // V6 (Inner, external layers)
    
    inline constexpr int DEFAULT = V6_ALL;     // 默认使用V6
}

// 跨度类型
namespace Span {
    inline constexpr int FULL_SPAN = 0;
    inline constexpr int INNER = 1;
    inline constexpr int INNER_L1 = 2;
    inline constexpr int INNER_L9 = 3;
    inline constexpr int FORCE_FULL_SPAN = 4;
    inline constexpr int INNER_UP = 5;
    inline constexpr int INNER_LOW = 6;
    
    inline constexpr int DEFAULT = INNER;      // 默认使用INNER
}

// 电荷重建
namespace ChargeReco {
    inline constexpr int QYJ = 0;
    inline constexpr int QH = 1;
    
    inline constexpr int DEFAULT = QYJ;        // 默认使用QYJ
}

// 方向
namespace Direction {
    inline constexpr int X = 0;
    inline constexpr int Y = 1;
    inline constexpr int XY = 2;
    
    inline constexpr int DEFAULT = XY;         // 默认使用XY
}

// 位置切割
namespace FiducialCuts {
    inline constexpr std::array<double, LAYER_COUNT> Y_POS{
        47.0, 40.0, 44.0, 44.0, 36.0, 36.0, 44.0, 44.0, 29.0
    };

    inline constexpr std::array<double, LAYER_COUNT> R_POS{
        62.0, 62.0, 46.0, 46.0, 46.0, 46.0, 46.0, 46.0, 43.0
    };
}

} // namespace Tracker
} // namespace AMS_Iso