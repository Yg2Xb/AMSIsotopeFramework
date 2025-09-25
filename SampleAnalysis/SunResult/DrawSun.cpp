#include "../Tool.h"
#include <memory>
#include <vector>

using namespace AMS_Iso;

// 简单的同位素信息结构体
struct IsotopeInfo {
    const char* name;
    int Z{4}, A;
    std::unique_ptr<TH1D> flux;
    std::unique_ptr<TF1> fluxFit;
    std::unique_ptr<TF1> rigidityFlux;
    Color_t color;
    double fracYmin;
    double fracYmax;
    
    // 构造函数
    IsotopeInfo(const char* n, int z, int a, Color_t c, double ymin, double ymax)
        : name(n), Z(z), A(a), color(c), fracYmin(ymin), fracYmax(ymax) {}
    
    // 移动构造函数
    IsotopeInfo(IsotopeInfo&& other) noexcept
        : name(other.name), Z(other.Z), A(other.A),
          flux(std::move(other.flux)),
          fluxFit(std::move(other.fluxFit)),
          rigidityFlux(std::move(other.rigidityFlux)),
          color(other.color),
          fracYmin(other.fracYmin),
          fracYmax(other.fracYmax) {}
};

void DrawSun() {
    // 通用画图样式设置
    auto SetStyle = [](auto* h, const char* ytitle, double ymin = -99, double ymax = -99) {
        h->SetTitle("");
        h->SetMarkerStyle(20);
        h->SetMarkerSize(1);
        h->GetXaxis()->SetTitle("Ek/n [GeV/n]");
        for (auto axis : {h->GetXaxis(), h->GetYaxis()}) {
            axis->SetTitleSize(0.065);
            axis->SetLabelSize(0.055);
            axis->SetTitleFont(62);
            axis->SetLabelFont(62);
            axis->SetTickLength(0.03);
        }
        h->GetXaxis()->SetTitleOffset(1.1);
        h->GetYaxis()->SetTitleOffset(1.15);
        h->GetYaxis()->SetTitle(ytitle);
        h->GetYaxis()->SetNdivisions(510);
        if (ymin >= -90) h->GetYaxis()->SetRangeUser(ymin, ymax);
    };

    // 初始化
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    const std::string outDir = "/eos/ams/user/z/zuhao/yanzx/Isotope/SunResult/";
    gSystem->Exec(("mkdir -p " + outDir).c_str());

    // 初始化同位素信息
    std::vector<IsotopeInfo> isotopes;
    isotopes.emplace_back("Be7", 4, 7, kBlue, 0.4, 0.8);
    isotopes.emplace_back("Be9", 4, 9, kGreen+2, 0.1, 0.6);
    isotopes.emplace_back("Be10", 4, 10, kOrange+1, 0.0, 0.16);

    // 打开文件并获取直方图
    std::unique_ptr<TFile> f(TFile::Open("/eos/ams/user/z/zetong/ISOTOPE/Results/Sun_results.root"));
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // 获取直方图并进行拟合
    for (auto& iso : isotopes) {
        iso.flux.reset((TH1D*)f->Get(Form("szt_fluxbe%d", iso.A)));
        if (!iso.flux) {
            std::cerr << "Error: isotope flux is null" << std::endl;
            return;
        }
        /*
        iso.fluxFit.reset(SplineFit(iso.flux.get(), 
                                   FluxAnalys::xpointsEkBe.data(),
                                   FluxAnalys::xpointsEkBe.size(),
                                   0x38, "b2e2",
                                   Form("f_%s_Ek", iso.name),
                                   FluxAnalys::xpointsEkBe.front(),
                                   FluxAnalys::xpointsEkBe.back()));
        */
    }

    // 创建画布并绘制flux
    std::unique_ptr<TCanvas> c(new TCanvas("c", "c", 800, 600));
    c->SetLogx();
    c->SetLogy(0);
    c->SetTopMargin(0.05);
    c->SetRightMargin(0.05);

    // 分别绘制每个同位素的flux
    for (const auto& iso : isotopes) {
        if (!iso.flux) continue;
        
        c->Clear();
        SetStyle(iso.flux.get(), Form("#Phi^{%d}Be [m^{-2} sr^{-1} s^{-1} GeV/n^{-1}]", iso.A));
        iso.flux->SetMarkerColor(iso.color);
        iso.flux->SetLineColor(iso.color);
        iso.flux->SetMarkerSize(1.);
        iso.flux->SetLineWidth(2);
        iso.flux->Draw("PE");
        //iso.fluxFit->SetLineWidth(2);
        //iso.fluxFit->SetLineColor(iso.color);
        c->SaveAs(Form("%sflux_%s.png", outDir.c_str(), iso.name));
    }

    // 计算总flux（保持不变）
    auto totalFlux = std::unique_ptr<TH1D>((TH1D*)isotopes[0].flux->Clone("totalFlux"));
    totalFlux->Reset();

    for (const auto& iso : isotopes) {
        if (!iso.flux) continue;
        totalFlux->Add(iso.flux.get());
    }

    // 分别绘制每个同位素的fraction
    c->SetLogy(0);
    for (const auto& iso : isotopes) {
        if (!iso.flux) continue;
        
        c->Clear();
        auto fraction = std::unique_ptr<TH1D>((TH1D*)iso.flux->Clone(Form("frac_%s", iso.name)));
        fraction->GetListOfFunctions()->Clear();  
        fraction->Divide(totalFlux.get());
        
        SetStyle(fraction.get(), Form("#Phi^{%d}Be Fraction", iso.A), iso.fracYmin, iso.fracYmax);
        fraction->SetMarkerColor(iso.color);
        fraction->SetLineColor(iso.color);
        fraction->Draw("PE");
        
        c->SaveAs(Form("%sfraction_%s.png", outDir.c_str(), iso.name));
    }

    // 创建rigidity flux
    for (auto& iso : isotopes) {
        if (!iso.fluxFit) continue;
        
        iso.rigidityFlux.reset(new TF1(Form("f_%s_R", iso.name),
            [&iso](double* x, double* p) {
                double R = x[0];
                double Ek = rigidityToKineticEnergy(R, iso.Z, iso.A);
                double flux_ek = iso.fluxFit->Eval(Ek);
                return flux_ek * dEk_dR(R, iso.Z, iso.A);
            }, 1, 30, 0));
    }

    // 创建总Be flux函数
    auto totalBeFlux = std::make_unique<TF1>("TotalBeFlux_Rig",
        [&isotopes](double* x, double*) -> double {
            double sum = 0;
            for (const auto& iso : isotopes) {
                if (iso.rigidityFlux) sum += iso.rigidityFlux->Eval(x[0]);
            }
            return sum;
        }, 1, 30, 0);

    // 读取AMS数据并比较
    std::unique_ptr<TFile> amsFile(TFile::Open("/eos/ams/user/z/zuhao/yanzx/Isotope/Bkg/AMS2011to2018PhysReport_BerylliumFlux.root"));
    if (amsFile && !amsFile->IsZombie()) {
        TGraphAsymmErrors* amsData = (TGraphAsymmErrors*)amsFile->Get("graph1");
        if (amsData) {
            c->Clear();
            c->SetLogx(1);
            c->SetLogy(0);
            SetStyle(amsData, "#PhiBe [m^{-2} sr^{-1} s^{-1} GV^{-1}]"); 
            amsData->SetMarkerStyle(20);
            amsData->SetMarkerSize(1);
            amsData->GetXaxis()->SetTitle("Rigidity[GV]");
            amsData->GetXaxis()->SetRangeUser(1.9,20);
            amsData->GetYaxis()->SetRangeUser(2e-9,0.2);
            amsData->Draw("AP");
            totalBeFlux->SetLineColor(kBlue);
            totalBeFlux->SetLineWidth(3);
            totalBeFlux->Draw("SAME");
            
            auto legend = std::make_unique<TLegend>(0.65, 0.65, 0.86, 0.86);
            legend->AddEntry(amsData, "PhysicsReport", "P");
            legend->AddEntry(totalBeFlux.get(), "This Work", "L");
            legend->Draw();
            
            c->SaveAs((outDir + "rigFlux_comparison.png").c_str());
        }
    }
}