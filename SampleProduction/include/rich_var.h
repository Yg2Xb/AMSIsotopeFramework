/***********************************************************
 *  File: rich_var.h
 *
 *  Modern C++ header file for AMS RICH basic values.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/

#pragma once

#include <array>

namespace AMS_Iso {
namespace RICH {

// RICH Beta类型
// 0: getBeta() - 默认方法
// 1: Beta - 仅使用未带电荷的hits重建速度
// 2: BetaRefit - 考虑每个hit的光电子数的Beta估计
inline constexpr int iRichBeta{0};

// RICH Z位置
inline constexpr std::array<double, 2> rich_z{-74.55, -75.55};

// NaF和Aerogel几何限制 (x,y)
inline constexpr std::array<double, 2> naf_limits{17.0, 17.0};
inline constexpr std::array<double, 2> aer_limits{19.0, 19.0};

// Aerogel几何参数
inline constexpr double disr_limit{58.5};                // 距离限制
inline constexpr std::array<double, 2> aer_corner{40.5, 40.5};  // 角落位置
inline constexpr std::array<double, 2> aer_border{28.5, 29.5};  // 边界位置

// Aerogel坏位置定义
inline constexpr int nbadpos{9};
inline constexpr std::array<double, 9> badpos_x{-2.0, 20.0, -29.0, 30.0, 16.0, 30.0, -30.0, -10.0, 12.0};
inline constexpr std::array<double, 9> badpos_wx{2.0, 2.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
inline constexpr std::array<double, 9> badpos_y{-60.0, -42.0, -18.0, 16.0, 22.0, 38.0, 40.0, 56.0, 56.0};
inline constexpr std::array<double, 9> badpos_wy{2.0, 2.0, 2.0, 2.0, 2.0, 1.5, 2.0, 1.5, 2.0};

// Aerogel坏瓦片定义
inline constexpr int kNBadTiles{7};
inline constexpr std::array<int, 7> kBadTile{3, 7, 12, 20, 87, 100, 108};

// PMT和光电子比例切割
inline constexpr std::array<int, 2> cut_pmt{2, 2};        // [0]:Aero [1]:NaF
inline constexpr std::array<double, 2> cut_per{0.4, 0.4}; // [0]:Aero [1]:NaF // 25.6.12 NaF 0.45->0.4

} // namespace RICH
} // namespace AMS_Iso