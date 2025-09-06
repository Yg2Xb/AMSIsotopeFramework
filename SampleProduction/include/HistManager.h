#ifndef HIST_MANAGER_H
#define HIST_MANAGER_H

#include <string>
#include <map>
#include <vector>
#include <memory>

#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>

// 转发 ProductRegistry 的声明以避免循环依赖
namespace AMS_Iso {
    struct HistogramBlueprint;
}


namespace AMS_Iso {

class HistManager {
public:
    explicit HistManager(const std::string& output_filename);
    ~HistManager();

    void CreateProducts(const std::vector<HistogramBlueprint>& blueprints);
    void Save();

    TH1* GetHistogram(const std::string& name);

    // 为方便使用提供类型安全的版本
    TH1F* GetHist1F(const std::string& name);
    TH2F* GetHist2F(const std::string& name);

private:
    std::unique_ptr<TFile> m_outputFile;
    std::map<std::string, TH1*> m_histograms;
};

} // namespace AMS_Iso

#endif // HIST_MANAGER_H

