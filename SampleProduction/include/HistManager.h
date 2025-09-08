#ifndef HISTMANAGER_H
#define HISTMANAGER_H

#include <memory>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "basic_var.h"

// 命名空间应与你的项目保持一致，这里假设是 AMS_Iso
namespace AMS_Iso {

// 使用 using 别名来简化智能指针的类型声明
using H1Ptr = std::unique_ptr<TH1F>;
using H2Ptr = std::unique_ptr<TH2F>;

class HistManager {
public:
    // 构造函数，使用 const 引用以提高效率
    HistManager(const std::string& output_filename,
                bool isISS,
                const std::vector<std::string>& chains,
                int charge,
                const IsotopeVar* iso,
                int UseMass);

    // 析构函数，智能指针会自动清理，这里无需手动操作
    ~HistManager() = default;

    // 保存所有直方图到文件
    void Save();

    // ======== ID 区域 ========
    std::vector<std::vector<std::vector<H1Ptr>>> ISS_IDH1; // [chain][det][iso]
    std::vector<std::vector<H2Ptr>> MC_IDH1; // [chain][det]

    std::vector<std::vector<std::vector<H2Ptr>>> IDH2; // [chain][det][iso/1]
    std::vector<std::vector<H2Ptr>> IDH3; // [chain][det]
    std::vector<H1Ptr> IDH4a; // [chain]
    std::vector<H1Ptr> IDH4b; // [chain]
    std::vector<H2Ptr> IDH5a; // [chain]
    std::vector<H2Ptr> IDH5b; // [chain]
    std::vector<H2Ptr> IDH6a; // [chain]
    std::vector<H2Ptr> IDH6b; // [chain]
    std::vector<H2Ptr> IDH7a; // [chain]
    std::vector<H2Ptr> IDH7b; // [chain]

    // ======== BKG 区域 ========
    // ISS
    std::vector<std::vector<std::vector<H1Ptr>>> ISS_BKGH1; // [chain][source][det]
    std::vector<std::vector<std::vector<std::vector<H2Ptr>>>> ISS_BKGH2; // [chain][source][charge_type][det]
    std::vector<std::vector<std::vector<H1Ptr>>> ISS_BKGH3; // [chain][source][det]
    std::vector<std::vector<std::vector<H2Ptr>>> ISS_BKGH4; // [chain][source][det]
    // MC
    std::vector<std::vector<H1Ptr>> MC_BKGH1; // [chain][det]
    std::vector<std::vector<std::vector<H1Ptr>>> MC_BKGH2; // [chain][det][iso]
    std::vector<std::vector<std::vector<H1Ptr>>> MC_BKGH3a; // [chain][det][iso]
    std::vector<std::vector<std::vector<H1Ptr>>> MC_BKGH3b; // [chain][det][iso]

    // ======== FLUX 区域 ========
    std::vector<std::vector<std::vector<std::vector<std::vector<H1Ptr>>>>> FLUXH1; // [chain][cut_group][num_den][det][iso/1]

    // iss
    std::vector<H1Ptr> ISS_FLUXH2; // rig expoT, only 1
    std::vector<std::vector<std::vector<H1Ptr>>> ISS_FLUXH3; // ek expoT, [det][iso]

    // MC 独有
    std::vector<std::vector<std::vector<std::vector<H1Ptr>>>> MC_FLUXH2; // acceptance num, [chain][cut_group][det][gen or rec]
    std::vector<H1Ptr> MC_FLUXH3; // only 1

private:
    std::unique_ptr<TFile> m_outputFile;

    // 辅助函数，用于简化直方图的创建
    template <typename HistType, typename... Args>
    std::unique_ptr<HistType> createHist(Args&&... args) {
        return std::make_unique<HistType>(std::forward<Args>(args)...);
    }
};

} // namespace AMS_Iso

#endif
