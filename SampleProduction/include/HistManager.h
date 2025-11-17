#ifndef HISTMANAGER_H
#define HISTMANAGER_H

#include <memory>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TChain.h"
#include "basic_var.h"

// 命名空间应与你的项目保持一致，这里假设是 AMS_Iso
namespace AMS_Iso {

// 使用 using 别名来简化智能指针的类型声明
using H1Ptr = std::unique_ptr<TH1D>;
using H2Ptr = std::unique_ptr<TH2D>;

class HistManager {
public:
    // 构造函数，使用 const 引用以提高效率
    HistManager(const std::string& output_filename,
                bool isISS,
                const std::vector<std::string>& chains,
                int charge,
                const IsotopeVar* iso,
                int UseMass,
                int FragmentZ); 

    // 新增：准备输出树（在事例循环之前调用）
    void PrepareFilteredTree(TChain* dataChain);
    
    // 新增：获取输出树供 selectdata 填充
    TTree* GetFilteredTree() const { return m_filteredTree; }
    
    // 修改 Save() 声明，添加可选参数控制是否保存 TTree
    void Save(bool saveTree = true);

    // ======== ID 区域 ========
    // ISS
    std::vector<std::vector<std::vector<H1Ptr>>> ISS_IDH1; // [chain][det][iso]
    // MC: 扩展为 Niso 维度（每个同位素一份 1/Mass vs Ek/n 模板）
    std::vector<std::vector<std::vector<H2Ptr>>> MC_IDH1; // [chain][det][iso]

    std::vector<std::vector<std::vector<H2Ptr>>> IDH2; // [chain][det][iso]
    std::vector<std::vector<H2Ptr>> IDH3; // [chain][det]

    // beta study - 增加Z维度
    std::vector<std::vector<H1Ptr>> IDH4a; // [z][chain]
    std::vector<std::vector<H1Ptr>> IDH4b; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5a; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5a2; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5a3; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5b; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5b2; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH5b3; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH6a; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH6b; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH7a; // [z][chain]
    std::vector<std::vector<H2Ptr>> IDH7b; // [z][chain]

    // ======== BKG 区域 ========
    // 统一研究碎裂产物为 Be 的三种同位素（Mass7/9/10），替代原来依据输入 iso->getIsotopeCount() 的维度

    // ISS
    std::vector<std::vector<std::vector<H1Ptr>>> ISS_BKGH1; // [chain][source][det]
    std::vector<std::vector<std::vector<std::vector<H2Ptr>>>> ISS_BKGH2; // [chain][source][det][charge_type]
    std::vector<std::vector<std::vector<H1Ptr>>> ISS_BKGH3; // [chain][source][det]
    // 将 ISS_BKGH4 扩展同位素维度为3（Be7/Be9/Be10）
    std::vector<std::vector<std::vector<std::vector<H2Ptr>>>> ISS_BKGH4; // [chain][source][det][BeIso=3]

    // MC
    std::vector<std::vector<H1Ptr>> MC_BKGH1; // [chain][det]
    std::vector<std::vector<H1Ptr>> MC_BKGH1b; // [chain][det], without Q, L1, Bkg cut
    // 以下三个均固定 iso 维度为 3（Be7,Be9,Be10）
    std::vector<std::vector<std::vector<H1Ptr>>> MC_BKGH2;  // [chain][det][BeIso=3]
    std::vector<std::vector<std::vector<H1Ptr>>> MC_BKGH2b;  // [chain][det][BeIso=3]
    std::vector<std::vector<std::vector<H1Ptr>>> MC_BKGH2c;  // [chain][det][BeIso=3]
    std::vector<std::vector<std::vector<H1Ptr>>> MC_BKGH3a; // [chain][det][BeIso=3]
    std::vector<std::vector<std::vector<H1Ptr>>> MC_BKGH3b; // [chain][det][BeIso=3]
    std::vector<std::vector<std::vector<H1Ptr>>> MC_BKGH3c; // [chain][det][BeIso=3]

    // ======== FLUX 区域 ========
    std::vector<std::vector<std::vector<std::vector<std::vector<H1Ptr>>>>> FLUXH1; // [chain][cut_group][num_den][det][iso/1]

    // iss
    std::vector<H1Ptr> ISS_FLUXH2; // rig expoT, only 1
    std::vector<std::vector<H1Ptr>> ISS_FLUXH3; // ek expoT, [det][iso]

    // MC 独有
    std::vector<std::vector<std::vector<std::vector<H1Ptr>>>> MC_FLUXH2; // acceptance num, [chain][cut_group][det][gen or rec]
    std::vector<H1Ptr> MC_FLUXH3; // only 1

private:
    std::unique_ptr<TFile> m_outputFile;
    TTree* m_filteredTree = nullptr;  // 新增：输出树（由 TFile 管理）

    // 新增：定义需要保存的分支
    static std::vector<std::string> GetActiveBranches();

    // 辅助函数，用于简化直方图的创建
    template <typename HistType, typename... Args>
    std::unique_ptr<HistType> createHist(Args&&... args) {
        // 1. 使用 new 创建直方图的裸指针
        auto hist = new HistType(std::forward<Args>(args)...);
        
        // 2. 关键步骤：阻止 ROOT 的 TFile/TDirectory 获取所有权
        // 这使得 std::unique_ptr 成为内存的唯一管理者
        hist->SetDirectory(nullptr);
        
        // 3. 将裸指针包装在 unique_ptr 中并返回
        return std::unique_ptr<HistType>(hist);
    }
};

} // namespace AMS_Iso

#endif