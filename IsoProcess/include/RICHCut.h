/***********************************************************
 *  File: RICHCut.h
 *
 *  Modern C++ header file for AMS RICH detector cuts.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/


 #pragma once

 #include <vector>
 #include <array>
 #include "rich_var.h"
 #include "Tool.h"
 
 // 前向声明
 class selectdata;
 
 namespace AMS_Iso {
 
 class RICHCut {
 public:
     explicit RICHCut(selectdata* event = nullptr);
     ~RICHCut() = default;
 
     RICHCut(const RICHCut&) = delete;
     RICHCut& operator=(const RICHCut&) = delete;
 
     // 基本属性访问
     double getBeta(int betaType = RICH::iRichBeta) const;
     bool isNaF() const { return richRegion == 1; }
     bool isAerogel() const { return richRegion == 0; }
 
     // 主要切割函数
     CutResult<5> cutBasic() const;          // 5个基本切割
     CutResult<3> cutCharge(int charge) const;  // 3个电荷切割
     
     std::array<double, 3> getModifiedPosition(bool interpolate) const;
     bool isBadAerogelPosition(double x, double y) const;
     bool isBadAerogelTile() const;
     CutResult<6> cutGeometry(bool interpolate = true) const;  // 6个几何切割
     
     CutResult<3> cutRICH(int charge, bool isISS, bool interpolate = true) const;  // 3个RICH切割
     CutResult<2> cutNoGeometry(int charge) const;  // 2个非几何切割
     CutResult<2> cutRICHforBkg(int charge, bool isISS, bool interpolate = true) const;  
 
 
 private:
     selectdata* event_;    // 不负责内存管理
     int richRegion;        // 0-Aerogel, 1-NaF
 
     // 内部辅助函数
 };
 
 } // namespace AMS_Iso