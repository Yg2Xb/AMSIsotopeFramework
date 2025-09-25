// MakeChargeHists.cpp
// Build: g++ -O2 -std=c++17 MakeChargeHists.cpp $(root-config --cflags --libs) -o MakeChargeHists
// Run:   ./MakeChargeHists

#include <iostream>
#include <ctime>
#include <vector>
#include <string>
#include <cmath>
#include <cassert>

#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TString.h>

#include "EventProcessor.hh"
#include "basic_var.h"
#include "Tool.h"

using namespace AMS_Iso;

enum HistTypeIdx {
    L1Signal_Idx = 0,
    L1Temp_Idx,
    L2Temp_Idx,
    L1Signal_unb_Idx,
    L1Temp_unb_Idx,
    L2Temp_unb_Idx
};

static const int NUM_BASE_TYPES = 3;     // L1Signal, L1Temp, L2Temp
static const int NUM_ALL_TYPES  = 6;     // 含 unb 版本

struct HistInfo {
    std::string name;
    std::string title;
    std::string xLabel;
    std::string yLabel;
    bool isUnb;     // 使用 unb 变量/判据
    HistTypeIdx idx;
};

static const std::vector<HistInfo> histTypes = {
    {"L1Signal",   "L1 Q Signal",            "L1 Q",       "E_{k}/n [GeV/n]", false, L1Signal_Idx},
    {"L1Temp",     "L1 Q Template",          "L1 Q",       "E_{k}/n [GeV/n]", false, L1Temp_Idx},
    {"L2Temp",     "L2 Q Template",          "L2 Q",       "E_{k}/n [GeV/n]", false, L2Temp_Idx},
    {"L1Signal_unb","L1 Q Signal (unb)",     "L1 Q (unb)", "E_{k}/n [GeV/n]", true,  L1Signal_unb_Idx},
    {"L1Temp_unb", "L1 Q Template (unb)",    "L1 Q (unb)", "E_{k}/n [GeV/n]", true,  L1Temp_unb_Idx},
    {"L2Temp_unb", "L2 Q Template (unb)",    "L2 Q",       "E_{k}/n [GeV/n]", true,  L2Temp_unb_Idx}
};

// 从 RigidityBins 推 Ek/Beta 分箱（A/Z = 2/4）
static inline void buildEkBetaBins(std::vector<double>& EkBins, std::vector<double>& BetaBins) {
    const auto& Rbins = Binning::RigidityBins;
    EkBins.clear(); EkBins.reserve(Rbins.size());
    BetaBins.clear(); BetaBins.reserve(Rbins.size());
    for (double R : Rbins) {
        double beta = rigidityToBeta(R, 2, 4, false);
        double ekn  = rigidityToKineticEnergy(R, 2, 4);
        BetaBins.push_back(beta);
        EkBins.push_back(ekn);
    }
}

// L1Q 动态窗（随 Z），供模板使用
static inline void getL1ChargeLowHigh(int z, double& low, double& high) {
    low  = 0.46 + (z - 3) * 0.16;
    high = (z <= 5) ? 0.65 : (0.65 + (z - 5) * 0.03);
}

// 质量：calculateMass(beta, alpha, Rig, charge) → invMass
static inline double invMass_from(double beta, double Rig, int Z) {
    const double alpha = 0.0;
    MassResult mr = calculateMass(beta, alpha, Rig, Z);
    return mr.isValid() ? mr.invMass : -1.0;
}

// 是否通过落在 cutoff 之外（按 BetaBins）
static inline bool passCutoffByBeta(const EventProcessor& ep, int det, int Z, int A, const std::vector<double>& BetaBins) {
    // 交由 ep 的现有逻辑，如果你已有 passBeyondCutoff 可直接用；否则用 BetaBins 判定
    // 这里示例直接调用 ep.betacutoffcut[det?] 若存在；否则总是 true
    // 为兼容性，这里保留一个通用接口：
    if constexpr (true) {
        // 若 EventProcessor 提供 passBeyondCutoff(det, Z, A, BetaBins)
        return ep.passBeyondCutoff(det, Z, A, BetaBins);
    }
    return true;
}

// 单条 chain 的完整生成
static void ProduceOneChain(const char* outPath,
                            const char* chainTag,    // "_L1InnerRig" 或 "_InnerRig"
                            const EventProcessor& epProto, // 用于获取常量（不读 entry）
                            double (*getRig)(const EventProcessor&)) {
    const double coe = 0.8;

    // 分箱
    std::vector<double> EkBins, BetaBins;
    buildEkBetaBins(EkBins, BetaBins);

    // 直方图容器
    TH2D* hQMaps[NUM_ELEMENTS][NUM_DETECTORS][NUM_ALL_TYPES] = {};
    TH1F* hL1Source[NUM_ELEMENTS][NUM_DETECTORS] = {};
    TH1F* hL2FragTotal[NUM_ELEMENTS][NUM_DETECTORS] = {};
    TH2F* hL2FragMass[NUM_ELEMENTS][NUM_DETECTORS] = {};

    // 创建所有 hist
    for (int ie = 0; ie < NUM_ELEMENTS; ++ie) {
        const auto& e = elements[ie];
        for (int id = 0; id < NUM_DETECTORS; ++id) {
            for (int it = 0; it < NUM_ALL_TYPES; ++it) {
                const auto& hi = histTypes[it];
                int nQbins; double qMin; double qMax;
                if (hi.idx == L1Signal_Idx || hi.idx == L1Signal_unb_Idx) {
                    nQbins = 600; qMin = 3.0; qMax = 9.0; // 宽些
                } else if (hi.idx == L2Temp_Idx || hi.idx == L2Temp_unb_Idx) {
                    nQbins = 400; qMin = e.charge - 2; qMax = e.charge + 2;
                } else {
                    nQbins = 400; qMin = e.charge - 2; qMax = e.charge + 2;
                }
                TString name = Form("%s_%s_%s%s", hi.name.c_str(), e.name.c_str(), detNames[id].c_str(), chainTag);
                TString title = Form("%s for %s @%s;%s;%s",
                                     hi.title.c_str(), e.name.c_str(), detNames[id].c_str(),
                                     hi.xLabel.c_str(), hi.yLabel.c_str());
                hQMaps[ie][id][it] = new TH2D(name, title,
                                              nQbins, qMin, qMax,
                                              (int)EkBins.size()-1, EkBins.data());
            }
            // 补充
            {
                TString n1 = Form("L1source_%s_%s%s", elements[ie].name.c_str(), detNames[id].c_str(), chainTag);
                TString t1 = Form("L1 source for %s @%s;E_{k}/n [GeV/n];Counts", elements[ie].name.c_str(), detNames[id].c_str());
                hL1Source[ie][id] = new TH1F(n1, t1, (int)EkBins.size()-1, EkBins.data());

                TString n2 = Form("L2frag_%s_%s%s", elements[ie].name.c_str(), detNames[id].c_str(), chainTag);
                TString t2 = Form("L2 frag total for %s @%s;E_{k}/n [GeV/n];Counts", elements[ie].name.c_str(), detNames[id].c_str());
                hL2FragTotal[ie][id] = new TH1F(n2, t2, (int)EkBins.size()-1, EkBins.data());

                TString n3 = Form("L2fragMass_%s_%s%s", elements[ie].name.c_str(), detNames[id].c_str(), chainTag);
                TString t3 = Form("L2 frag invMass for %s @%s;1/M [1/GeV];E_{k}/n [GeV/n]", elements[ie].name.c_str(), detNames[id].c_str());
                hL2FragMass[ie][id] = new TH2F(n3, t3, 240, 0.0, 3.0, (int)EkBins.size()-1, EkBins.data());
            }
        }
    }

    // 输入树
    TChain chain("saveTree");
    chain.Add("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/Bor_BeToC.root");
    chain.Add("/eos/ams/user/z/zuhao/yanzx/Isotope/NewData/7_15_tree.root");

    EventProcessor ep;
    ep.setBranchAddresses(&chain);

    Long64_t nEntries = chain.GetEntries();
    std::cout << "[Info] Entries: " << nEntries << " " << chainTag << std::endl;

    // 常量窗
    const double inner_low_min = 3.45; // 你旧版的下限
    // 上限为 Z + 0.45（事件内按 Z 计算）

    std::clock_t c0 = std::clock();
    time_t t0 = time(nullptr);

    for (Long64_t i = 0; i < nEntries; ++i) {
        ep.getEntry(i);

        // 允许两套基本/无偏选择并存：作图时分别用
        bool passBasic = ep.passBasicCut();
        bool passUnb   = ep.passUnbiasedCut();
        if (!passBasic && !passUnb) continue;

        for (int ie = 0; ie < NUM_ELEMENTS; ++ie) {
            const auto& E = elements[ie];
            int Z = E.charge, A = E.mass;

            // 动态 L1 窗（仅模板用）
            double l1w_low, l1w_high;
            getL1ChargeLowHigh(Z, l1w_low, l1w_high);

            for (int id = 0; id < NUM_DETECTORS; ++id) {
                if (!ep.passDetectorCut(id)) continue;
                if (!passCutoffByBeta(ep, id, Z, A, BetaBins)) continue;

                bool needRich = (id != 0); // TOF不需要；NaF/AGL需要
                double Ek   = (id == 0) ? ep.TOFEk   : (id == 1 ? ep.NaFEk   : ep.AGLEk);
                double beta = (id == 0) ? ep.TOFBeta : (id == 1 ? ep.NaFBeta : ep.AGLBeta);

                // Q 值
                double qL1     = ep.trk_ql1;
                double qL1_unb = ep.trk_ql1_unbias;
                double qL2     = ep.trk_ql2;
                double qInner  = ep.QTrkInner[0];

                // 共同 cut
                bool inner_ok_signal = (qInner > inner_low_min && qInner < Z + 0.45); // 仅 Signal 使用这个语义
                bool inner_ok_temp   = std::abs(qInner - Z) < 0.45 * coe;             // 模板用 ±0.45*coe
                bool tof_ok          = (ep.tofQ() > (Z - 0.6*coe) && ep.tofQ() < (Z + 1.5*coe));

                // rich
                bool rich_ok = (!needRich) || ep.passRichQCut(Z);

                // L1 动态窗函数（模板）
                auto inL1Dyn = [&](double q){ return (q > (Z - l1w_low*coe) && q < (Z + l1w_high*coe)); };

                // ---------- L1Signal / L1Signal_unb ----------
                if (passBasic && inner_ok_signal && ep.getNormalL1XY() && ep.passInnerRMSCut())
                    hQMaps[ie][id][L1Signal_Idx]->Fill(qL1, Ek);
                if (passUnb   && inner_ok_signal && ep.passInnerRMSCut())
                    hQMaps[ie][id][L1Signal_unb_Idx]->Fill(qL1_unb, Ek);

                // ---------- L1Temp / L1Temp_unb ----------
                if (inner_ok_temp && tof_ok && rich_ok) {
                    if (passBasic && ep.getNormalL1XY() && inL1Dyn(qL1))
                        hQMaps[ie][id][L1Temp_Idx]->Fill(qL1, Ek);
                    if (passUnb && inL1Dyn(qL1_unb))
                        hQMaps[ie][id][L1Temp_unb_Idx]->Fill(qL1_unb, Ek);
                }

                // ---------- L2Temp / L2Temp_unb ----------
                // 需要 L2 状态良好且 L1Q 在动态窗内（各自对应），并沿用 inner_ok_temp、tof_ok、rich_ok
                if (ep.getL2XY() && ep.getQl2StatusCut() && std::abs(ep.QTrkInner[2] - Z) < 0.45 * coe && tof_ok && rich_ok) {
                    if (passBasic && ep.getNormalL1XY() && inL1Dyn(qL1))
                        hQMaps[ie][id][L2Temp_Idx]->Fill(qL2, Ek);
                    if (passUnb && inL1Dyn(qL1_unb))
                        hQMaps[ie][id][L2Temp_unb_Idx]->Fill(qL2, Ek);
                }

                // ---------- 补充直方图 ----------
                // L1source：按 L1Signal 的 cut（basic、inner_ok_signal、L1XY、RMS），对 Ek 计数
                if (passBasic && ep.getNormalL1XY() && ep.passInnerRMSCut() && inner_ok_signal)
                    hL1Source[ie][id]->Fill(Ek);

                // L2frag：满足（inner_ok_signal）且 L1 任一版本在动态窗内；探测器 cut 已过；RICH: TOF 不要，NaF/AGL 要
                {
                    bool l1_any_dyn = inL1Dyn(qL1) || inL1Dyn(qL1_unb);
                    if (inner_ok_signal && l1_any_dyn && rich_ok)
                        hL2FragTotal[ie][id]->Fill(Ek);
                }

                // L2fragMass：在 L2frag 基础上再计算质量；三探测器都需要 rich_ok
                {
                    bool l1_any_dyn = inL1Dyn(qL1) || inL1Dyn(qL1_unb);
                    if (inner_ok_signal && l1_any_dyn && ep.passRichQCut(Z)) {
                        double Rig = getRig(ep); // 该 chain 的刚度
                        double invM = invMass_from(beta, Rig, Z);
                        if (invM > 0 && std::isfinite(invM))
                            hL2FragMass[ie][id]->Fill(invM, Ek);
                    }
                }

            } // det
        } // element

        if ((i % 1000000) == 0) {
            double cpu = double(std::clock() - c0) / CLOCKS_PER_SEC;
            double wall = difftime(time(nullptr), t0);
            double rss = getCurrentRSS_MB();
            std::cout << "Entry " << i << "/" << nEntries
                      << " CPU " << cpu << "s Wall " << wall << "s RSS " << rss << " MB " << chainTag << std::endl;
        }
    } // loop

    // 写文件
    std::unique_ptr<TFile> fout(TFile::Open(outPath, "RECREATE"));
    for (int ie = 0; ie < NUM_ELEMENTS; ++ie) {
        for (int id = 0; id < NUM_DETECTORS; ++id) {
            for (int it = 0; it < NUM_ALL_TYPES; ++it) {
                if (hQMaps[ie][id][it]) hQMaps[ie][id][it]->Write();
            }
            if (hL1Source[ie][id])    hL1Source[ie][id]->Write();
            if (hL2FragTotal[ie][id]) hL2FragTotal[ie][id]->Write();
            if (hL2FragMass[ie][id])  hL2FragMass[ie][id]->Write();
        }
    }
    fout->Close();
    std::cout << "[Info] Wrote " << outPath << std::endl;
}

int main() {
    // 输入链在 ProduceOneChain 内部设置；如需修改，请改那里
    // 为两个 chain 分别产出
    {
        // chain = L1InnerRig（偏置链）
        EventProcessor epDummy;
        ProduceOneChain("/eos/user/z/zixuan/Isotope/chargeTempFit/ChargeTemp_Hist_Z4to8_L1InnerRig.root",
                        "_L1InnerRig", epDummy,
                        [](const EventProcessor& ep){ return ep.L1InnerRig; });
    }
    {
        // chain = InnerRig（即 unbias，用 InnerRig）
        EventProcessor epDummy;
        ProduceOneChain("/eos/user/z/zixuan/Isotope/chargeTempFit/ChargeTemp_Hist_Z4to8_InnerRig.root",
                        "_InnerRig", epDummy,
                        [](const EventProcessor& ep){ return ep.InnerRig; });
    }
    return 0;
}