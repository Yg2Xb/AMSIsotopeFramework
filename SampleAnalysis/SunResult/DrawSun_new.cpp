#include "../Tool.h"
#include <memory>
#include <vector>

using namespace AMS_Iso;

struct IsotopeInfo {
    const char* name;
    int Z{4}, A;
    std::unique_ptr<TH1D> flux;
    Color_t color;
    double fracYmin;
    double fracYmax;

    IsotopeInfo(const char* n, int z, int a, Color_t c, double ymin, double ymax)
        : name(n), Z(z), A(a), color(c), fracYmin(ymin), fracYmax(ymax) {}

    IsotopeInfo(IsotopeInfo&& other) noexcept
        : name(other.name), Z(other.Z), A(other.A),
          flux(std::move(other.flux)),
          color(other.color),
          fracYmin(other.fracYmin),
          fracYmax(other.fracYmax) {}

    IsotopeInfo(const IsotopeInfo&) = delete;
    IsotopeInfo& operator=(const IsotopeInfo&) = delete;
};

static void SetHistStyle(TH1* h, const char* ytitle, bool linearX = false,
                         double ymin = std::numeric_limits<double>::quiet_NaN(),
                         double ymax = std::numeric_limits<double>::quiet_NaN()) {
    h->SetTitle("");
    h->SetMarkerStyle(20);
    h->SetMarkerSize(1.0);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitle("Ek/n [GeV/n]");
    h->GetYaxis()->SetTitle(ytitle);
    for (auto axis : {h->GetXaxis(), h->GetYaxis()}) {
        axis->SetTitleSize(0.065);
        axis->SetLabelSize(0.055);
        axis->SetTitleFont(62);
        axis->SetLabelFont(62);
        axis->SetTickLength(0.03);
    }
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleOffset(1.2);
    h->GetYaxis()->SetNdivisions(508);
    if (ymin == ymin && ymax == ymax) { // check not-NaN
        h->GetYaxis()->SetRangeUser(ymin, ymax);
    }
    // 注：是否 logx 由画布设置控制，这里只负责样式
}

void DrawSun_new() {
    // 全局样式
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // 输出目录与文件
    const std::string outDir = "/eos/user/z/zixuan/Isotope/Sun/";
    gSystem->Exec(Form("mkdir -p %s", outDir.c_str()));
    std::unique_ptr<TFile> fout(TFile::Open((outDir + "Sun.root").c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) {
        std::cerr << "[ERROR] Cannot create output file at " << outDir << "Sun.root\n";
        return;
    }

    // 输入 ROOT 文件（xrootd）
    const char* inPath = "root://eosams.cern.ch//eos/ams/user/z/zetong/ISOTOPE/Results/Sun_results.root";
    std::unique_ptr<TFile> f(TFile::Open(inPath));
    if (!f || f->IsZombie()) {
        std::cerr << "[ERROR] Opening input file failed: " << inPath << "\n";
        return;
    }

    // 初始化同位素
    std::vector<IsotopeInfo> isotopes;
    isotopes.emplace_back("Be7", 4, 7, kBlue,       0.40, 0.80);
    isotopes.emplace_back("Be9", 4, 9, kGreen + 2,  0.10, 0.60);
    isotopes.emplace_back("Be10",4,10, kOrange + 1, 0.00, 0.16);

    // 读取 flux 直方图
    for (auto& iso : isotopes) {
        auto key = Form("szt_fluxbe%d", iso.A);
        iso.flux.reset(dynamic_cast<TH1D*>(f->Get(key)));
        if (!iso.flux) {
            std::cerr << "[ERROR] Missing histogram: " << key << "\n";
            return;
        }
        iso.flux->SetDirectory(nullptr); // 脱离输入文件
    }

    // 画布
    std::unique_ptr<TCanvas> c(new TCanvas("c", "c", 900, 700));
    c->SetTopMargin(0.05);
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.165);

    // 1) 各自的 flux（对数或线性由你决定；这里沿用你之前：logx=true, logy=false）
    for (const auto& iso : isotopes) {
        c->Clear();
        c->SetLogx(true);
        c->SetLogy(false);
        SetHistStyle(iso.flux.get(),
                     Form("#Phi^{%d}Be [m^{-2}sr^{-1}s^{-1}GeV/n^{-1}]", iso.A));
        iso.flux->SetMarkerColor(iso.color);
        iso.flux->SetLineColor(iso.color);
        iso.flux->Draw("PE");
        c->SaveAs(Form("%sflux_%s.png", outDir.c_str(), iso.name));
        // 保存到 root
        fout->cd();
        auto hsave = static_cast<TH1D*>(iso.flux->Clone(Form("flux_%s", iso.name)));
        hsave->Write();
    }

    // 2) flux × (E/n)（你原来写的是 pow(cen,1)，这里维持）
    for (const auto& iso : isotopes) {
        auto* h3 = static_cast<TH1D*>(iso.flux->Clone(Form("fluxEn_%s", iso.name)));
        h3->Reset();
        for (int i = 1; i <= iso.flux->GetNbinsX(); ++i) {
            const double val = iso.flux->GetBinContent(i);
            const double err = iso.flux->GetBinError(i);
            const double cen = iso.flux->GetBinCenter(i);
            h3->SetBinContent(i, val * cen);
            h3->SetBinError(i, err * cen);
        }

        c->Clear();
        c->SetLogx(true);
        c->SetLogy(false);
        SetHistStyle(h3, Form("#Phi^{%d}Be #times (E_{k}/n) [m^{-2}sr^{-1}s^{-1}]", iso.A));
        h3->SetMarkerColor(iso.color);
        h3->SetLineColor(iso.color);

        // 自动 y 轴范围
        double ymin = +1e300, ymax = -1e300;
        for (int i = 1; i <= h3->GetNbinsX(); ++i) {
            double v = h3->GetBinContent(i);
            if (v <= 0) continue;
            ymin = std::min(ymin, v);
            ymax = std::max(ymax, v);
        }
        if (ymax > 0 && ymin < 1e290) h3->GetYaxis()->SetRangeUser(ymin * 0.8, ymax * 1.5);

        h3->Draw("PE");
        c->SaveAs(Form("%sfluxEn_%s.png", outDir.c_str(), iso.name));

        // 保存到 root
        fout->cd();
        auto hsave = static_cast<TH1D*>(h3->Clone(Form("fluxEn_%s", iso.name)));
        hsave->Write();

        delete h3;
    }

    // 3) 总 Be flux
    auto totalFlux = std::unique_ptr<TH1D>(static_cast<TH1D*>(isotopes[0].flux->Clone("totalFlux")));
    totalFlux->Reset();
    for (const auto& iso : isotopes) totalFlux->Add(iso.flux.get());
    fout->cd();
    totalFlux->Write();

    // 4) 各自 fraction
    for (const auto& iso : isotopes) {
        auto fraction = std::unique_ptr<TH1D>(static_cast<TH1D*>(iso.flux->Clone(Form("fraction_%s", iso.name))));
        fraction->GetListOfFunctions()->Clear();
        fraction->Divide(totalFlux.get());

        c->Clear();
        c->SetLogx(true);
        c->SetLogy(false);
        SetHistStyle(fraction.get(), Form("#Phi^{%d}Be Fraction", iso.A), false, iso.fracYmin, iso.fracYmax);
        fraction->SetMarkerColor(iso.color);
        fraction->SetLineColor(iso.color);
        fraction->Draw("PE");
        c->SaveAs(Form("%sfraction_%s.png", outDir.c_str(), iso.name));

        // 保存到 root
        fout->cd();
        auto hsave = static_cast<TH1D*>(fraction->Clone());
        hsave->Write();
    }

    // 5) 计算并绘制 ratio: Be10 / Be9（线性 x，magenta）
    TH1D* hBe10 = isotopes[2].flux.get(); // A=10
    TH1D* hBe9  = isotopes[1].flux.get(); // A=9
    auto ratio = std::unique_ptr<TH1D>(static_cast<TH1D*>(hBe10->Clone("Be10Be9_ratio")));
    ratio->GetListOfFunctions()->Clear();
    ratio->Divide(hBe9);

    c->Clear();
    c->SetLogx(false);  // 线性 x
    c->SetLogy(false);
    // 参照示意图设置 y 轴范围，默认 [0.1, 0.5]，你可调整
    SetHistStyle(ratio.get(), "#Phi^{10}Be / #Phi^{9}Be", true, 0.10, 0.50);
    ratio->SetMarkerColor(kMagenta + 1);
    ratio->SetLineColor(kMagenta + 1);
    ratio->SetMarkerStyle(20);
    ratio->SetMarkerSize(1.1);
    ratio->Draw("PE");
    c->SaveAs(Form("%sBe10Be9_ratio.png", outDir.c_str()));

    // 保存 ratio 到 root
    fout->cd();
    auto hsaveRatio = static_cast<TH1D*>(ratio->Clone());
    hsaveRatio->Write();

    // 结束
    fout->Write();
    fout->Close();
    std::cout << "[DONE] All figures saved to " << outDir << " and histograms written to " << outDir << "Sun.root\n";
}