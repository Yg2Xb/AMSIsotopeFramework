#include "IsoToolbox/Logger.h"
#include "IsoToolbox/ConfigManager.h"
#include "DataModel/AMSDstTreeA.h" // <-- 1. 包含我们的数据模型

#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <vector>

using Logger = IsoToolbox::Logger;

int main(int argc, char** argv) {
    
    // --- 1. 配置和日志初始化 (这部分我们已经完成了) ---
    // (为了简化，我们暂时硬编码配置，稍后会从YAML加载)
    Logger::GetInstance().SetLevel(Logger::Level::INFO);
    Logger::Info("Starting IsoProcess application...");

    // --- 2. 准备数据输入 ---
    Logger::Info("Setting up input data chain...");
    TChain* chain = new TChain("amstreea");
    
    // 我们暂时硬编码一个文件名作为测试。
    // 未来我们会从YAML配置文件中读取一个文件列表。
    const char* inputFile = "root://eosams.cern.ch//eos/ams/user/z/zuhao/yanzx/ams_data/amsd69n_IHEPQcutDST/1305853512_10.root";
    chain->Add(inputFile);

    Long64_t nEntries = chain->GetEntries();
    if (nEntries == 0) {
        Logger::Error("No entries found in the input file: {}", inputFile);
        return 1;
    }
    Logger::Info("Chain created with {} entries.", nEntries);

    // --- 3. 连接数据模型 ---
    // 这是我们新框架的核心步骤
    AMSDstTreeA eventReader;       // 创建数据读取器对象
    eventReader.Init(chain);       // 调用Init()来完成所有SetBranchAddress

    // --- 4. 事件循环 ---
    Logger::Info("Starting event loop...");
    for (Long64_t i = 0; i < nEntries; ++i) {
        
        chain->GetEntry(i); // 将数据读入到eventReader的成员变量中

        // --- 验证步骤：打印前10个事件的信息 ---
        if (i < 10) {
            std::cout << "==> Processing Entry " << i 
                      << ", Run: " << eventReader.run 
                      << ", Event: " << eventReader.event 
                      << std::endl;
        }

        // ------------------------------------------------------------------
        //  未来，这里就是我们调用所有物理模块(Cuts, Hists)的地方
        //
        //  if (!cutManager.pass(eventReader)) continue;
        //  histManager.fill(eventReader);
        //
        // ------------------------------------------------------------------
    }
    Logger::Info("Event loop finished.");

    delete chain;
    return 0;
}