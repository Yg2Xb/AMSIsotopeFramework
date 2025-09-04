/***********************************************************
 *  File: RICHCut.cpp
 *
 *  Modern C++ implementation file for AMS RICH detector cuts.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/
 #include <cmath>
 #include "RICHCut.h"
 #include "selectdata.h"
 
 namespace AMS_Iso {
 
 RICHCut::RICHCut(selectdata* event)
     : event_(event)
     , richRegion(event ? event->rich_NaF : 0)
 {}
 
 double RICHCut::getBeta(int betaType) const {
     if (!event_) return -1.0;
     return event_->rich_beta[betaType];  // 直接返回
 }
 
 CutResult<5> RICHCut::cutBasic() const {
     if (!event_) return CutResult<5>();
 
     std::array<bool, 5> cuts{
         event_->rich_good && event_->rich_clean,
         event_->rich_pb > 0.01,
         event_->rich_pmt > RICH::cut_pmt[richRegion],
         (event_->rich_npe[0] / event_->rich_npe[2]) > RICH::cut_per[richRegion],
         event_->rich_beta[0] > 0
     };
 
     return CutResult<5>(cuts);
 }
 
 CutResult<3> RICHCut::cutCharge(int charge) const {
     std::array<double, 2> lowerQ{event_->tof_ql[2], event_->tof_ql[3]};
     double tofLQ = Tools::calculateAverage(lowerQ.data(), 2, 0);
 
     double richQ = std::sqrt(event_->rich_q[0]);
     std::array<bool, 3> cuts{
         tofLQ > (charge - 0.6),
         richQ > (charge - 1),
         richQ < (charge + 2)
     };
 
     return CutResult<3>(cuts);
 }
 
 std::array<double, 3> RICHCut::getModifiedPosition(bool interpolate) const {
     std::array<double, 3> position;
     std::copy_n(event_->rich_pos, 3, position.begin());
     
     if (interpolate) {
         auto interpolateResult = Tools::calculateXYAtZ(event_->tk_pos,
                                 event_->tk_dir,
                                 RICH::rich_z[richRegion],
                                 event_->itrtrack);
         
         if (interpolateResult) {
             // 如果插值成功
             position[0] = interpolateResult->x;
             position[1] = interpolateResult->y;
             position[2] = RICH::rich_z[richRegion];
         } else {
             // 如果插值失败
             //std::cout << "Interpolation failed for event = " 
             //         << event_->event  << ", itrtrack = " << event_->itrtrack << std::endl;
             return { -999, -999, -999 };
         }
         
         /*
         try {
             Tools::modifyPositionByZ(RICH::rich_z[richRegion], position,
                                    event_->rich_theta, event_->rich_phi);
         } catch (const std::invalid_argument& e) {
             std::cout << "There are illegal rich angles: events = " 
                      << event_->event 
                      << ", theta = " << event_->rich_theta 
                      << ", phi = " << event_->rich_phi << std::endl;
             return { -999, -999, -999 };
         }
         */
     }
     
     return position;
 }
 
 bool RICHCut::isBadAerogelPosition(double x, double y) const {
     for (int i = 0; i < RICH::nbadpos; ++i) {
         if (x > RICH::badpos_x[i] && x < RICH::badpos_x[i] + RICH::badpos_wx[i] &&
             y > RICH::badpos_y[i] && y < RICH::badpos_y[i] + RICH::badpos_wy[i]) {
             return true;
         }
     }
     return false;
 }
 
 bool RICHCut::isBadAerogelTile() const {
     return std::find(RICH::kBadTile.begin(), RICH::kBadTile.end(), 
                     event_->rich_tile) != RICH::kBadTile.end();
 }
 
 CutResult<6> RICHCut::cutGeometry(bool interpolate) const {
     auto pos = getModifiedPosition(interpolate);
     
     if (isNaF()) {
         bool naFCut = (std::abs(pos[0]) < RICH::naf_limits[0]) && 
                      (std::abs(pos[1]) < RICH::naf_limits[1]);
         return CutResult<6>({naFCut, true, true, true, true, true});
     }
 
     std::array<bool, 6> cuts{
         std::abs(pos[0]) > RICH::aer_limits[0] || std::abs(pos[1]) > RICH::aer_limits[1],
         std::hypot(pos[0], pos[1]) < RICH::disr_limit,
         std::abs(pos[0]) < RICH::aer_corner[0] || std::abs(pos[1]) < RICH::aer_corner[1],
         std::max(std::abs(pos[0]), std::abs(pos[1])) <= RICH::aer_border[0] || std::max(std::abs(pos[0]), std::abs(pos[1])) >= RICH::aer_border[1],
         !isBadAerogelPosition(pos[0], pos[1]),
         !isBadAerogelTile()
     };
 
     return CutResult<6>(cuts);
 }
 
 CutResult<3> RICHCut::cutRICH(int charge, bool isISS, bool interpolate) const {
     std::array<bool, 3> cuts{
         //cutGeometry(interpolate).total,
         isISS ? event_->rich_goodgeo : cutGeometry(interpolate).total,
         cutBasic().total,
         cutCharge(charge).total
     };
 
     return CutResult<3>(cuts);
 }
 
 CutResult<2> RICHCut::cutNoGeometry(int charge) const {
     std::array<bool, 2> cuts{
         cutBasic().total,
         cutCharge(charge).total
     };
 
     return CutResult<2>(cuts);
 }
 
 CutResult<2> RICHCut::cutRICHforBkg(int charge, bool isISS, bool interpolate) const {
     std::array<bool, 2> cuts{
         //cutGeometry(interpolate).total,
         isISS ? event_->rich_goodgeo : cutGeometry(interpolate).total,
         event_->rich_NaF ? true : (event_->rich_npe[0] / event_->rich_npe[2]) > RICH::cut_per[richRegion]
     };
 
     return CutResult<2>(cuts);
 }
 
 
 } // namespace AMS_Iso