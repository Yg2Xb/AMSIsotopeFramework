#ifndef EVENT_PROCESSOR_HH
#define EVENT_PROCESSOR_HH

#include <vector>
#include <string>
#include <bitset>
#include <cmath>
#include <TTree.h>
#include "Tool.h"
#include "basic_var.h"

namespace AMS_Iso {

// 最小质量元素：Be-7, B-10, C-12, N-14, O-16
struct ElementInfo {
    int charge;
    int mass;
    std::string name;
};
static const std::vector<ElementInfo> elements = {
    {4,  7, "Beryllium"},
    {5, 10, "Boron"},
    {6, 12, "Carbon"},
    {7, 14, "Nitrogen"},
    {8, 16, "Oxygen"}
};
static const std::vector<std::string> detNames = {"TOF", "NaF", "AGL"};
static constexpr int NUM_ELEMENTS = 5;
static constexpr int NUM_DETECTORS = 3;

// bitmask（与旧代码一致）
constexpr unsigned int BKG_BASIC_MASK     = 0x870780;
constexpr unsigned int BKG_UNB_BASIC_MASK = 0xE00780;
constexpr unsigned int BKG_TOF_MASK       = 0x0C000000;
constexpr unsigned int BKG_RICH_MASK      = 0x30000000;

// cutoff 安全因子
constexpr double SAFETY_FACTORS[3] = {1.06, 1.005, 1.0005};

class EventProcessor {
public:
    // branches
    Float_t trk_ql1{}, trk_ql1_unbias{}, trk_ql2{};
    Float_t QTrkInner[3]{}, rich_q[2]{}, tof_ql[4]{}, tk_qrmn[2][3]{}, tk_qln[2][9][3]{};
    int tk_qls[9]{}, tk_hitb[2]{};
    bool ChargeCutsSelectJudge[5][3]{};
    bool betacutoffcut[5]{};
    double usedEk{}, TOFEk{}, NaFEk{}, AGLEk{};
    double TOFBeta{}, NaFBeta{}, AGLBeta{};
    double L1InnerRig{}, InnerRig{}, cutOffRig{};
    unsigned int cutStatus = 0;
    bool isMC = false;

private:
    TTree* fTree = nullptr;

    static inline bool in_cut(double val, double z_val, double coe, double low, double high) {
        // val ∈ (z - coe*low, z + coe*high)
        return val > z_val - coe * low && val < z_val + coe * high;
    }

public:
    void setBranchAddresses(TTree* tree) {
        fTree = tree;
        if (!fTree) return;
        fTree->SetBranchAddress("trk_ql1", &trk_ql1);
        fTree->SetBranchAddress("trk_ql1_unbias", &trk_ql1_unbias);
        fTree->SetBranchAddress("trk_ql2", &trk_ql2);
        fTree->SetBranchAddress("L1InnerRig", &L1InnerRig);
        fTree->SetBranchAddress("InnerRig",   &InnerRig);
        fTree->SetBranchAddress("usedEk", &usedEk);
        fTree->SetBranchAddress("TOFEk", &TOFEk);
        fTree->SetBranchAddress("NaFEk", &NaFEk);
        fTree->SetBranchAddress("AGLEk", &AGLEk);
        fTree->SetBranchAddress("TOFBeta", &TOFBeta);
        fTree->SetBranchAddress("NaFBeta", &NaFBeta);
        fTree->SetBranchAddress("AGLBeta", &AGLBeta);
        fTree->SetBranchAddress("cutOffRig", &cutOffRig);
        fTree->SetBranchAddress("cutStatus", &cutStatus);
        fTree->SetBranchAddress("QTrkInner", QTrkInner);
        fTree->SetBranchAddress("rich_q", rich_q);
        fTree->SetBranchAddress("tof_ql", tof_ql);
        fTree->SetBranchAddress("tk_qls", tk_qls);
        fTree->SetBranchAddress("tk_hitb", tk_hitb);
        fTree->SetBranchAddress("tk_qrmn", tk_qrmn);
        // fTree->SetBranchAddress("tk_qln", tk_qln);
        fTree->SetBranchAddress("ChargeCutsSelectJudge", ChargeCutsSelectJudge);
        fTree->SetBranchAddress("betacutoffcut", betacutoffcut);
    }

    inline void getEntry(Long64_t i) { if (fTree) fTree->GetEntry(i); }

    // masks
    inline bool passBasicCut() const    { return (cutStatus & BKG_BASIC_MASK)     == BKG_BASIC_MASK; }
    inline bool passUnbiasedCut() const { return (cutStatus & BKG_UNB_BASIC_MASK) == BKG_UNB_BASIC_MASK; }
    inline bool passTofCut() const      { return (cutStatus & BKG_TOF_MASK)       == BKG_TOF_MASK; }
    inline bool passRichCut() const     { return (cutStatus & BKG_RICH_MASK)      == BKG_RICH_MASK; }

    // tracker 状态
    inline bool getNormalL1XY() const   { return std::bitset<32>(tk_hitb[0]).test(0); }
    inline bool getL2XY() const         { return std::bitset<32>(tk_hitb[0]).test(1); }
    inline bool getQl2StatusCut() const { return ((tk_qls[1] & 0x10013D) == 0); }

    // TOF 合并 Q
    inline double getTofQUp() const {
        double val = 0.5 * (tof_ql[0] + tof_ql[1]);
        return (tof_ql[0] * tof_ql[1] == 0) ? val * 2 : val;
    }
    inline double getTofQLow() const {
        double val = 0.5 * (tof_ql[2] + tof_ql[3]);
        return (tof_ql[2] * tof_ql[3] == 0) ? val * 2 : val;
    }
    inline double getRichQ() const { return std::sqrt(std::max(0.f, rich_q[0])); }

    // 探测器通过（TOF 不需 RICH；NaF/AGL 需）
    inline bool passDetectorCut(int idet) const {
        if (idet == 0) return isValidBeta(TOFBeta) && passTofCut();
        if (idet == 1) return isValidBeta(NaFBeta) && passRichCut();
        if (idet == 2) return isValidBeta(AGLBeta) && passRichCut();
        return false;
    }

    // cutoff 判定：用 beta 边界
    inline bool passBeyondCutoff(int idet, int z, int mass, const std::vector<double>& rbins_beta) const {
        const double beta = (idet == 0) ? TOFBeta : ((idet == 1) ? NaFBeta : AGLBeta);
        const int betaBin = findBin(rbins_beta, beta);
        return (betaBin >= 0) ? isBeyondCutoff(rbins_beta[betaBin], cutOffRig, SAFETY_FACTORS[idet], z, mass, isMC) : false;
    }

    inline bool passInnerRMSCut() const { return tk_qrmn[0][2] < 0.55; }

    // 旧的 L1Signal/unbiased 判定保留（主程序里已替换为动态 L1Q + innerQ 固定窗）
    inline bool passQTrkInner0UnbiasedL1Cut(int z) const { return QTrkInner[0] > 3.45 && QTrkInner[0] < z + 0.45; }
    inline bool passQTrkInner0L1Cut(int z) const         { return QTrkInner[0] > 3.50 && QTrkInner[0] < z + 0.50; }

    // 模板用 innerQ（保持旧逻辑 ±0.45*coe）
    inline bool passQTrkInnerCut(int index, int z, double coe = 0.8, double low = 0.45, double high = 0.45) const {
        return in_cut(QTrkInner[index], z, coe, low, high);
    }

    // TOF 上层窗：围绕 Z 非对称 Z-0.6*coe ~ Z+1.5*coe
    inline bool passTofQUpCut(int z, double coe = 0.8, double low = 0.6, double high = 1.5) const {
        return in_cut(getTofQUp(), z, coe, low, high);
    }

    // L1 模板窗：调用处传入动态 L1Q 窗 low/high
    inline bool passTrkQL1Cut_L2Template(int z, double coe = 0.8, double low = 0.46, double high = 0.65) const {
        return in_cut(trk_ql1, z, coe, low, high);
    }
    inline bool passTrkQL1UnbiasCut_L2Template(int z, double coe = 0.8, double low = 0.46, double high = 0.65) const {
        return in_cut(trk_ql1_unbias, z, coe, low, high);
    }

    // RICH Q cut（旧逻辑）
    inline bool passRichQCut(int z) const {
        return getRichQ() > z - 1 && getRichQ() < z + 2 && getTofQLow() > z - 0.6;
    }
};

} // namespace AMS_Iso
#endif