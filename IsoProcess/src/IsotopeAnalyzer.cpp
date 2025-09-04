/***********************************************************
 *  File: IsotopeAnalyzer.cpp
 *
 *  Implementation file for AMS Isotopes analysis framework.
 *
 *  History:
 *    20241029 - created by ZX.Yan
 ***********************************************************/
#include "IsotopeAnalyzer.h"
#include <iostream>
#include <stdexcept>

namespace AMS_Iso {

void IsotopeAnalyzer::setConfig(const TString& outDir, const TString& outName,
                               const TString& inData, const TString& inOptions,
                               int UseMass) {
    outDir_ = outDir;
    outName_ = outName;
    inData_ = inData;
    inOptions_ = inOptions;
    UseMass_ = UseMass;
    isISS_ = !inOptions.Contains("|MC");
}

void IsotopeAnalyzer::readDataFrom(TChain* chain, const TString& filename) {
    if (!chain) return;
    
    // 添加文件到链中
    chain->Add(filename);
    
    // 验证是否成功添加
    if (chain->GetEntries() == 0) {
        std::cerr << "Error: No entries found in " << filename << std::endl;
        return;
    }
}

void IsotopeAnalyzer::initialize() {
    try {
        // 创建数据链
        dataChain = std::make_unique<TChain>("amstreea");
        if (!dataChain) {
            throw std::runtime_error("Failed to create TChain");
        }
        readDataFrom(dataChain.get(), inData_);

        // 创建输出文件
        fOutRoot = std::make_unique<TFile>(
            Form("%s/%s", outDir_.Data(), outName_.Data()), "recreate");
        if (!fOutRoot || fOutRoot->IsZombie()) {
            throw std::runtime_error("Failed to create output file");
        }
        std::cout << "Output file: " << fOutRoot->GetName() << std::endl;
        fOutRoot->cd();

        // 设置Tree分支
        setupBranches();

        // 获取同位素信息
        int charge = std::atoi(inOptions_.Data());
        isotope = &AMS_Iso::getIsotopeVar(charge);
        if (charge != isotope->getCharge()) {
            throw std::runtime_error("Wrong Isotope Charge Read in!");
        }

        std::cout << "Isotope name: " << isotope->getName() << std::endl;
        std::cout << "Isotope number: " << isotope->getIsotopeCount() << std::endl;

        // 初始化直方图管理器
        histManager_ = std::make_unique<HistogramManager>();
        histManager_->createHistograms(isotope, UseMass_);

        std::cout << "Initialization completed successfully" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error in initialization: " << e.what() << std::endl;
        cleanup();
        throw;
    }
}

void IsotopeAnalyzer::setupBranches() {
    dataChain->SetBranchStatus("*", 0);
    const std::vector<const char*> activeBranches = {
        
        //"run", "mmom", "rich_NaF", "rich_pos", "rich_theta", "rich_phi", "rich_usedm", "rich_hit", "rich_good", "rich_clean", "rich_pb", "rich_pmt", "rich_npe", "rich_q", "tof_pos", "tof_chisc", "tof_chist", "tof_ql", "tk_pos", "tk_dir"
        
        //"run", "event", "rich_NaF", "ntrack", "mpar", "mmom", "mch", "mevcoo", "mevdir", "mfcoo", "mfdir", "mpare", "mparp", "mparc", "msume", "msumc", "ntrmccl", "mtrmom", "mtrpar", "mtrcoo", "mtrcoo1", "mtrpri", "mtrz", "mevcoo1", "mevdir1", "mevmom1", "tk_qin","tof_ql", "rich_q", "tk_exqln", "tk_qln"
        
        //"run", "event", "rich_NaF","rich_pos","rich_pos_tkItr","rich_theta","rich_phi", "rich_theta_tkItr", "rich_phi_tkItr", "tof_goodgeo","iso_ll","iso_da1","iso_da2", "iso_nda1", "iso_nda2", "btstat_new","cutoffpi"

        //"run", "event", "rich_NaF","tk_qin","tof_ql", "rich_q", "iso_ll", "tk_exqln", "tk_qln", "btstat_new","cutoffpi","mcutoffi"
        
        "run", "event", "rich_NaF", "tk_qrmn", "tk_qln", "tk_qls", "tk_hitb", "tof_ql", "tk_exqln", "rich_q", "btstat_new", "cutoffpi", "mcutoffi", "mpar", "mmom", "mch", "mtrmom", "mparc","mtrpar"
    };
    for (const auto* branch : activeBranches) {
        dataChain->SetBranchStatus(branch, 1);
    }
    //dataChain->SetBranchStatus("*", 1);

    saveTree = std::unique_ptr<TTree>(dataChain->CloneTree(0));
    if (!saveTree) {
        throw std::runtime_error("Failed to clone tree");
    }
    saveTree->SetName("saveTree");
    dataChain->SetBranchStatus("*", 1);
}

void IsotopeAnalyzer::write() {
    if (!fOutRoot || fOutRoot->IsZombie()) {
        throw std::runtime_error("Invalid output file for writing");
    }

    try {
        fOutRoot->cd();
        
        // 使用直方图管理器写入直方图
        histManager_->write(fOutRoot.get());

        // 写入Tree
        if (saveTree) {
            saveTree->Write("saveTree", TObject::kOverwrite);
        }

        fOutRoot->Flush();

    } catch (const std::exception& e) {
        std::cerr << "Error writing histograms: " << e.what() << std::endl;
        throw;
    }
}

void IsotopeAnalyzer::cleanup() {
    // 清理直方图管理器
    if (histManager_) {
        histManager_.reset();
    }

    // 清理文件和树
    if (saveTree) {
        saveTree->SetDirectory(nullptr);
        saveTree.reset();
    }
    
    if (fOutRoot) {
        if (!fOutRoot->IsZombie()) {
            fOutRoot->Close();
        }
        fOutRoot.reset();
    }
    
    if (dataChain) {
        dataChain.reset();
    }

    // 重置其他成员
    isotope = nullptr;
    outDir_.Clear();
    outName_.Clear();
    inData_.Clear();
    inOptions_.Clear();
    UseMass_ = -1;
    isISS_ = true;
}

} // namespace AMS_Iso