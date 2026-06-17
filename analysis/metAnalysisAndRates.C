// To execute: root -b -l -q 'metAnalysisAndRates.C+'

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TSystem.h"
#include "TROOT.h"
#include "analysisHelperFunctions.h"

const int   nColors = 6;
const Color_t cols[nColors] = { kRed, kBlue, kGreen+2, kMagenta+1, kOrange+1, kCyan+2 };

// Skip background events with passHSTP == false (HSTP filter removes high-energy
// pileup transients that are not modelled correctly in dijet MC).
const bool applyHSTPFilter = true;

// Variable bin edges for turn-on histograms [GeV]:
//   5 GeV bins [0, 100],  10 GeV bins [100, 200],  20 GeV bins [200, 400],  50 GeV bins [400, 600]
const int nTurnOnBins = 44;
const Double_t turnOnBinEdges[45] = {
     0,   5,  10,  15,  20,  25,  30,  35,  40,  45,  50,
    55,  60,  65,  70,  75,  80,  85,  90,  95, 100,
    110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
    220, 240, 260, 280, 300, 320, 340, 360, 380, 400,
    450, 500, 550, 600
};

// Variable bin edges for MET magnitude distributions [GeV]:
//   5 GeV bins  [0, 200],  10 GeV bins  [200, 400],  50 GeV bins  [400, 600]
const int nMETBins = 64;
const Double_t metBinEdges[65] = {
      0,   5,  10,  15,  20,  25,  30,  35,  40,  45,  50,
     55,  60,  65,  70,  75,  80,  85,  90,  95, 100,
    105, 110, 115, 120, 125, 130, 135, 140, 145, 150,
    155, 160, 165, 170, 175, 180, 185, 190, 195, 200,
    210, 220, 230, 240, 250, 260, 270, 280, 290, 300,
    310, 320, 330, 340, 350, 360, 370, 380, 390, 400,
    450, 500, 550, 600
};

// -----------------------------------------------------------------------
// Normalize histogram to unit area
void normalizeHist(TH1F* h) {
    if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
}

// -----------------------------------------------------------------------
// Overlay signal (red) vs background (blue) normalized to unity, logy, min 1e-8
void drawOverlay(TH1F* sig, TH1F* back, const std::string& title,
                 const std::string& xLabel, const std::string& outputPath) {
    normalizeHist(sig);
    normalizeHist(back);

    sig->SetLineColor(kRed);   sig->SetLineWidth(2);
    back->SetLineColor(kBlue); back->SetLineWidth(2);
    sig->SetTitle(title.c_str());
    sig->GetXaxis()->SetTitle(xLabel.c_str());
    sig->GetYaxis()->SetTitle(Form("Fraction of Events / %.4g GeV", sig->GetBinWidth(1)));

    double ymax = std::max(sig->GetMaximum(), back->GetMaximum()) * 5.0; // headroom for logy
    sig->SetMaximum(ymax);
    sig->SetMinimum(1e-8);

    TCanvas c("c", title.c_str(), 700, 600);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
    gPad->SetLogy();

    sig->Draw("HIST");
    back->Draw("HIST SAME");

    TLegend leg(0.55, 0.75, 0.88, 0.88);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.04);
    leg.AddEntry(sig,  "Signal",     "l");
    leg.AddEntry(back, "Background", "l");
    leg.Draw();

    c.SaveAs(outputPath.c_str());
}

// -----------------------------------------------------------------------
// Compute median of a TH1F
double getMedian(TH1F* h) {
    double median = 0.0, prob = 0.5;
    h->GetQuantiles(1, &median, &prob);
    return median;
}

// -----------------------------------------------------------------------
// |DeltaPhi| folded into [0, pi]
double absDeltaPhi(double phi1, double phi2) {
    double d = std::abs(phi1 - phi2);
    if (d > M_PI) d = 2.0*M_PI - d;
    return d;
}

// -----------------------------------------------------------------------
// Overlay signal vs background for MET X/Y components with mean & median in legend
void drawComponentOverlay(TH1F* sig, TH1F* back, const std::string& title,
                          const std::string& xLabel, const std::string& outputPath) {
    double sigMean   = sig->GetMean(),  sigMedian   = getMedian(sig);
    double backMean  = back->GetMean(), backMedian  = getMedian(back);

    normalizeHist(sig);
    normalizeHist(back);

    sig->SetLineColor(kRed);   sig->SetLineWidth(2);
    back->SetLineColor(kBlue); back->SetLineWidth(2);
    sig->SetTitle(title.c_str());
    sig->GetXaxis()->SetTitle(xLabel.c_str());
    sig->GetYaxis()->SetTitle(Form("Fraction of Events / %.4g GeV", sig->GetBinWidth(1)));

    double ymax = std::max(sig->GetMaximum(), back->GetMaximum()) * 5.0;
    sig->SetMaximum(ymax);
    sig->SetMinimum(1e-8);

    TCanvas c("c", title.c_str(), 700, 600);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetTicks(0,1);
    gPad->SetLogy();

    sig->Draw("HIST");
    back->Draw("HIST SAME");

    TLegend leg(0.14, 0.88, 0.54, 0.96);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.024);
    leg.AddEntry(sig,  Form("Sig. (mean=%+.1f,med=%+.1f GeV)", sigMean,  sigMedian),  "l");
    leg.AddEntry(back, Form("Back. (mean=%+.1f,med=%+.1f GeV)", backMean, backMedian), "l");
    leg.Draw();

    c.SaveAs(outputPath.c_str());
}

// -----------------------------------------------------------------------
// Overlay two algorithms (e.g. GEP vs gFEX) for signal and background on one canvas.
// sig1/back1 drawn solid, sig2/back2 drawn dashed; red=algo1, blue=algo2.
void drawAlgoComparison(TH1F* sig1, TH1F* back1, TH1F* sig2, TH1F* back2,
                        const std::string& label1, const std::string& label2,
                        const std::string& xLabel, const std::string& outputPath,
                        const std::string& signalName = "") {
    normalizeHist(sig1); normalizeHist(back1);
    normalizeHist(sig2); normalizeHist(back2);

    double ymax = std::max({sig1->GetMaximum(), back1->GetMaximum(),
                            sig2->GetMaximum(), back2->GetMaximum()}) * 5.0;
    std::string yTitle = Form("Fraction of Events / %.4g GeV", sig1->GetBinWidth(1));

    auto style = [&](TH1F* h, Color_t col, int ls) {
        h->SetLineColor(col); h->SetLineWidth(2); h->SetLineStyle(ls);
        h->GetXaxis()->SetTitle(xLabel.c_str());
        h->GetYaxis()->SetTitle(yTitle.c_str());
        h->SetMaximum(ymax); h->SetMinimum(1e-8);
    };
    style(sig1,  kRed,  1); style(back1, kRed,  2);
    style(sig2,  kBlue, 1); style(back2, kBlue, 2);

    TCanvas c("c", "", 700, 600);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
    gPad->SetLogy();

    sig1->Draw("HIST");
    back1->Draw("HIST SAME");
    sig2->Draw("HIST SAME");
    back2->Draw("HIST SAME");

    double legTop = 0.88, legH = 0.06 * (!signalName.empty() + 4);
    TLegend leg(0.38, legTop - legH, 0.88, legTop);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);
    if (!signalName.empty())
        leg.AddEntry((TObject*)nullptr, signalName.c_str(), "");
    leg.AddEntry(sig1,  (label1 + " (sig)").c_str(), "l");
    leg.AddEntry(back1, (label1 + " (bkg)").c_str(), "l");
    leg.AddEntry(sig2,  (label2 + " (sig)").c_str(), "l");
    leg.AddEntry(back2, (label2 + " (bkg)").c_str(), "l");
    leg.Draw();

    c.SaveAs(outputPath.c_str());
}

// -----------------------------------------------------------------------
// Overlay multiple signal (or background) histograms from different algorithm configs
void drawOverlayMulti(std::vector<TH1F*>& sigs, std::vector<TH1F*>& backs,
                      const std::vector<std::string>& labels,
                      const std::string& title, const std::string& xLabel,
                      const std::string& outputPath, const std::string& signalName = "") {
    if (sigs.empty() && backs.empty()) return;
    for (auto* h : sigs)  normalizeHist(h);
    for (auto* h : backs) normalizeHist(h);

    double ymax = 0;
    for (auto* h : sigs)  ymax = std::max(ymax, h->GetMaximum());
    for (auto* h : backs) ymax = std::max(ymax, h->GetMaximum());
    ymax *= 5.0;
    TH1F* refH = !sigs.empty() ? sigs[0] : backs[0];
    std::string yTitle = Form("Fraction of Events / %.4g GeV", refH->GetBinWidth(1));

    TCanvas c("c", title.c_str(), 700, 600);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
    gPad->SetLogy();

    // header + 2 entries (sig+bkg) per config
    int nConfigs = (int)std::max(sigs.size(), backs.size());
    double legTop = 0.88, legH = 0.06 * (!signalName.empty() + 2 * nConfigs);
    TLegend leg(0.38, legTop - legH, 0.88, legTop);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);
    if (!signalName.empty())
        leg.AddEntry((TObject*)nullptr, signalName.c_str(), "");

    bool first = true;
    for (unsigned int i = 0; i < sigs.size(); i++) {
        sigs[i]->SetLineColor(cols[i % nColors]);
        sigs[i]->SetLineWidth(2);
        sigs[i]->SetLineStyle(1);
        sigs[i]->GetXaxis()->SetTitle(xLabel.c_str());
        sigs[i]->GetYaxis()->SetTitle(yTitle.c_str());
        sigs[i]->SetMaximum(ymax);
        sigs[i]->SetMinimum(1e-8);
        sigs[i]->Draw(first ? "HIST" : "HIST SAME");
        first = false;
        leg.AddEntry(sigs[i], (labels[i] + " (sig)").c_str(), "l");
    }
    for (unsigned int i = 0; i < backs.size(); i++) {
        backs[i]->SetLineColor(cols[i % nColors]);
        backs[i]->SetLineWidth(2);
        backs[i]->SetLineStyle(2);
        backs[i]->GetXaxis()->SetTitle(xLabel.c_str());
        backs[i]->GetYaxis()->SetTitle(yTitle.c_str());
        backs[i]->SetMaximum(ymax);
        backs[i]->SetMinimum(1e-8);
        backs[i]->Draw(first ? "HIST" : "HIST SAME");
        first = false;
        leg.AddEntry(backs[i], (labels[i] + " (bkg)").c_str(), "l");
    }
    leg.Draw();
    c.SaveAs(outputPath.c_str());
}

// -----------------------------------------------------------------------
// Build a TGraph of signal efficiency vs threshold from an unweighted signal histogram.
// Efficiency = integral from threshold bin upward / total integral.
TGraph* makeEffVsThresholdGraph(TH1F* h, Color_t col) {
    int nBins = h->GetNbinsX();
    double total = h->Integral(1, nBins);
    if (total <= 0) return new TGraph();
    std::vector<double> thresholds, effs;
    for (int iBin = 1; iBin <= nBins; iBin++) {
        thresholds.push_back(h->GetBinLowEdge(iBin));
        effs.push_back(h->Integral(iBin, nBins) / total);
    }
    TGraph* g = new TGraph((int)thresholds.size(), thresholds.data(), effs.data());
    g->SetLineColor(col);   g->SetLineWidth(2);
    g->SetMarkerColor(col); g->SetMarkerStyle(20); g->SetMarkerSize(0.7);
    return g;
}

// -----------------------------------------------------------------------
// Signal efficiency vs threshold overlay for multiple algorithm configs
void drawEffVsThresholdMulti(std::vector<TH1F*>& sigs,
                              const std::vector<std::string>& labels,
                              const std::string& title, const std::string& xLabel,
                              const std::string& outputPath, const std::string& signalName = "") {
    if (sigs.empty()) return;
    TCanvas c("c", title.c_str(), 700, 600);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);

    int nConfigs = (int)sigs.size();
    double legTop = 0.88, legH = 0.06 * (!signalName.empty() + nConfigs);
    TLegend leg(0.38, legTop - legH, 0.88, legTop);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);
    if (!signalName.empty())
        leg.AddEntry((TObject*)nullptr, signalName.c_str(), "");

    const Color_t mcols[] = { kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kOrange+1, kCyan+2, kViolet+1 };
    const int nMcols = 8;
    std::vector<TGraph*> graphs;
    for (unsigned int i = 0; i < sigs.size(); i++) {
        TGraph* g = makeEffVsThresholdGraph(sigs[i], mcols[i % nMcols]);
        g->SetTitle((title + ";" + xLabel + ";Signal Efficiency").c_str());
        graphs.push_back(g);
        g->Draw(i == 0 ? "AL" : "L SAME");
        leg.AddEntry(g, labels[i].c_str(), "l");
    }
    // Fix y-axis range after drawing
    if (!graphs.empty()) {
        graphs[0]->GetYaxis()->SetRangeUser(0.0, 1.05);
    }
    leg.Draw();
    c.SaveAs(outputPath.c_str());
    for (auto* g : graphs) delete g;
}

// -----------------------------------------------------------------------
// Build a TGraphErrors of rate vs threshold from a weighted background histogram.
// Rate error = sqrt(sum of squared bin errors) from threshold bin upward.
TGraphErrors* makeRateGraph(TH1F* h, Color_t col) {
    int nBins = h->GetNbinsX();
    std::vector<double> thresholds, rates, xErrs, rateErrs;
    for (int iBin = 1; iBin <= nBins; iBin++) {
        double err = 0.0;
        double rate = h->IntegralAndError(iBin, nBins, err);
        thresholds.push_back(h->GetBinLowEdge(iBin));
        rates.push_back(rate);
        xErrs.push_back(0.0);
        rateErrs.push_back(err);
    }
    TGraphErrors* g = new TGraphErrors(thresholds.size(),
                                       thresholds.data(), rates.data(),
                                       xErrs.data(),      rateErrs.data());
    g->SetLineColor(col);   g->SetLineWidth(2);
    g->SetMarkerColor(col); g->SetMarkerStyle(20); g->SetMarkerSize(0.7);
    return g;
}

// -----------------------------------------------------------------------
// Rate vs threshold from a weighted background histogram (rate in Hz)
void drawRateVsThreshold(TH1F* back_weighted, const std::string& title,
                         const std::string& xLabel, const std::string& outputPath) {
    TGraphErrors* g = makeRateGraph(back_weighted, kBlue);
    g->SetTitle((title + ";" + xLabel + ";Rate [Hz]").c_str());

    TCanvas c("c", title.c_str(), 700, 600);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
    gPad->SetLogy();
    g->Draw("AP");
    c.SaveAs(outputPath.c_str());
    delete g;
}

// -----------------------------------------------------------------------
// Rate vs threshold overlay for multiple algorithm configs
void drawRateVsThresholdMulti(const std::vector<TH1F*>& backs_weighted,
                              const std::vector<std::string>& labels,
                              const std::string& title, const std::string& xLabel,
                              const std::string& outputPath, const std::string& signalName = "",
                              const std::string& yLabel = "Rate [Hz]") {
    if (backs_weighted.empty()) return;
    TCanvas c("c", title.c_str(), 700, 600);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
    gPad->SetLogy();

    int nConfigs = (int)backs_weighted.size();
    double legTop = 0.88, legH = 0.06 * (!signalName.empty() + nConfigs);
    TLegend leg(0.38, legTop - legH, 0.88, legTop);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);
    if (!signalName.empty())
        leg.AddEntry((TObject*)nullptr, signalName.c_str(), "");

    const Color_t mcols[] = { kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kOrange+1, kCyan+2, kViolet+1 };
    const int nMcols = 8;
    std::vector<TGraphErrors*> graphs;
    for (unsigned int i = 0; i < backs_weighted.size(); i++) {
        TGraphErrors* g = makeRateGraph(backs_weighted[i], mcols[i % nMcols]);
        g->SetTitle((title + ";" + xLabel + ";" + yLabel).c_str());
        graphs.push_back(g);
        g->Draw(i == 0 ? "AP" : "P SAME");
        leg.AddEntry(g, labels[i].c_str(), "lp");
    }
    leg.Draw();
    c.SaveAs(outputPath.c_str());
    for (auto* g : graphs) delete g;
}

// -----------------------------------------------------------------------
// Overlay N normalized distributions on one canvas (no sig/bkg distinction)
void drawMultiDist(std::vector<TH1F*> hists, const std::vector<std::string>& labels,
                   const std::string& title, const std::string& xLabel,
                   const std::string& outputPath, bool logy = true,
                   const std::string& units = "GeV") {
    if (hists.empty()) return;
    const Color_t mcols[] = { kBlack, kRed, kBlue, kGreen+2, kMagenta+1,
                               kOrange+1, kCyan+2, kViolet+2, kGray+2, kPink+1 };
    for (auto* h : hists) normalizeHist(h);
    double ymax = 0;
    for (auto* h : hists) ymax = std::max(ymax, h->GetMaximum());
    ymax *= logy ? 5.0 : 1.3;
    std::string yTitle = units.empty()
        ? Form("Fraction of Events / %.4g", hists[0]->GetBinWidth(1))
        : Form("Fraction of Events / %.4g %s", hists[0]->GetBinWidth(1), units.c_str());

    TCanvas c("c", title.c_str(), 800, 600);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
    if (logy) gPad->SetLogy();

    double legTop = 0.88, legH = 0.055 * hists.size();
    TLegend leg(0.40, legTop - legH, 0.88, legTop);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);

    for (unsigned int i = 0; i < hists.size(); i++) {
        hists[i]->SetLineColor(mcols[i % 10]); hists[i]->SetLineWidth(2);
        hists[i]->SetTitle(title.c_str());
        hists[i]->GetXaxis()->SetTitle(xLabel.c_str());
        hists[i]->GetYaxis()->SetTitle(yTitle.c_str());
        hists[i]->SetMaximum(ymax); hists[i]->SetMinimum(logy ? 1e-8 : 0.0);
        hists[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
        leg.AddEntry(hists[i], labels[i].c_str(), "l");
    }
    leg.Draw();
    c.SaveAs(outputPath.c_str());
}

// -----------------------------------------------------------------------
// Find the MET threshold that gives a target rate (Hz) from a weighted background histogram
double findThreshold(TH1F* back_hw, double targetRateHz) {
    int nBins = back_hw->GetNbinsX();
    for (int iBin = 1; iBin <= nBins; iBin++) {
        double rate = back_hw->Integral(iBin, nBins);
        if (rate <= targetRateHz) return back_hw->GetBinLowEdge(iBin);
    }
    return back_hw->GetXaxis()->GetXmax(); // above range — threshold not reached
}

// -----------------------------------------------------------------------
// Draw turn-on curves (efficiency vs truth NonInt MET) for multiple algorithms,
// overlaid on one canvas at a single rate point
void drawTurnOnOverlay(std::vector<TH1F*> effs, const std::vector<std::string>& labels,
                       const std::string& title, const std::string& outputPath,
                       const std::vector<double>& thresholds = {},
                       const std::string& rateLabel = "",
                       TH1F* truthDist = nullptr,
                       double legX1 = 0.45, double legY1 = 0.15) {
    if (effs.empty()) return;
    const Color_t  mcols[]    = { kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kOrange+1 };
    const Style_t  mkstyles[] = { 20, 21, 22, 23, 29, 33 };

    TCanvas c("c", title.c_str(), 700, 600);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);

    int nLegRows = (int)effs.size() + (!rateLabel.empty() ? 1 : 0);
    double legH = 0.058 * nLegRows;
    TLegend leg(legX1, legY1, legX1 + 0.33, legY1 + legH);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);
    if (!rateLabel.empty())
        leg.AddEntry((TObject*)nullptr, rateLabel.c_str(), "");

    for (unsigned int i = 0; i < effs.size(); i++) {
        effs[i]->SetLineColor(mcols[i % 6]);
        effs[i]->SetMarkerColor(mcols[i % 6]);
        effs[i]->SetMarkerStyle(mkstyles[i % 6]);
        effs[i]->SetMarkerSize(0.8);
        effs[i]->SetLineWidth(2);
        effs[i]->SetTitle(title.c_str());
        effs[i]->GetXaxis()->SetTitle("Truth MET_{NonInt} [GeV]");
        effs[i]->GetYaxis()->SetTitle("Signal Efficiency");
        effs[i]->SetMaximum(1.10); effs[i]->SetMinimum(0.0);
        effs[i]->Draw(i == 0 ? "EP" : "EP SAME");
        std::string lbl = labels[i];
        if (!thresholds.empty() && i < (int)thresholds.size())
            lbl += Form(" (thr=%.1f GeV)", thresholds[i]);
        leg.AddEntry(effs[i], lbl.c_str(), "lp");
    }

    // Draw normalized signal truth MET distribution as gray shaded histogram on top
    if (truthDist) {
        TH1F* hTruth = (TH1F*)truthDist->Clone("hTruthDist_turnOn");
        hTruth->SetDirectory(0);
        if (hTruth->Integral() > 0) hTruth->Scale(1.0 / hTruth->Integral());
        hTruth->SetFillColorAlpha(kGray+1, 0.35);
        hTruth->SetLineColor(kGray+2);
        hTruth->SetLineWidth(1);
        hTruth->Draw("HIST SAME");
    }

    leg.Draw();
    c.SaveAs(outputPath.c_str());
}

// -----------------------------------------------------------------------
// Overlay multiple Rate-vs-Efficiency TGraph* on one canvas
void drawRateVsEffOverlay(std::vector<TGraph*> graphs,
                          const std::vector<std::string>& labels,
                          const std::string& outputPath,
                          const std::string& signalName = "") {
    if (graphs.empty()) return;
    const Color_t  mcols[]    = { kBlack, kRed, kBlue, kGreen+2, kMagenta+1, kOrange+1, kCyan+2, kViolet+1 };
    const Style_t  mkstyles[] = { 20, 21, 22, 23, 29, 33, 20, 21 };
    const int nStyles = 8;

    TCanvas c("c", "", 700, 600);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
    gPad->SetLogy(); gPad->SetLogx();

    int nLeg = (int)(!signalName.empty()) + (int)graphs.size();
    TLegend leg(0.30, 0.19, 0.83, 0.19 + 0.06 * nLeg);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.025);
    if (!signalName.empty())
        leg.AddEntry((TObject*)nullptr, signalName.c_str(), "");

    for (unsigned int i = 0; i < graphs.size(); i++) {
        graphs[i]->SetLineColor(mcols[i % nStyles]);
        graphs[i]->SetMarkerColor(mcols[i % nStyles]);
        graphs[i]->SetMarkerStyle(mkstyles[i % nStyles]);
        graphs[i]->SetMarkerSize(0.8);
        graphs[i]->SetLineWidth(2);
        graphs[i]->GetXaxis()->SetTitle("Signal Efficiency");
        graphs[i]->GetYaxis()->SetTitle("Estimated Background Rate [Hz]");
        graphs[i]->Draw(i == 0 ? "AP" : "P SAME");
        leg.AddEntry(graphs[i], labels[i].c_str(), "lp");
    }
    leg.Draw();
    c.SaveAs(outputPath.c_str());
}

// -----------------------------------------------------------------------
// Combined gFEX + GEP selection: best (gFEX threshold, GEP threshold) at a target rate
struct BestThresh2D {
    double t1;   // gFEX MET threshold [GeV]
    double t2;   // GEP MET threshold [GeV]
    double eff;  // signal efficiency at best point
    double rate; // background rate [Hz] at best point
};

BestThresh2D findBestThresholds2D(const RateEff2DOut& out2D, double targetRateHz) {
    TH2* hEff  = out2D.hEff_vsThr_vsR;
    TH2* hRate = out2D.hRate_vsThr_vsR;
    int nbX = hEff->GetNbinsX(), nbY = hEff->GetNbinsY();
    BestThresh2D best = {hEff->GetXaxis()->GetXmax(), hEff->GetYaxis()->GetXmax(), 0.0, 0.0};
    for (int ix = 1; ix <= nbX; ++ix) {
        for (int iy = 1; iy <= nbY; ++iy) {
            double rate = hRate->GetBinContent(ix, iy);
            double eff  = hEff ->GetBinContent(ix, iy);
            if (rate <= targetRateHz && eff > best.eff) {
                best.t1   = hEff->GetXaxis()->GetBinLowEdge(ix);
                best.t2   = hEff->GetYaxis()->GetBinLowEdge(iy);
                best.eff  = eff;
                best.rate = rate;
            }
        }
    }
    return best;
}

// Extract Pareto-optimal rate-vs-eff frontier from a 2D threshold scan
TGraph* extractFrontier2D(const RateEff2DOut& out2D) {
    TH2* hEff  = out2D.hEff_vsThr_vsR;
    TH2* hRate = out2D.hRate_vsThr_vsR;
    int nbX = hEff->GetNbinsX(), nbY = hEff->GetNbinsY();
    std::vector<std::pair<double,double>> pts;
    for (int ix = 1; ix <= nbX; ++ix)
        for (int iy = 1; iy <= nbY; ++iy)
            pts.emplace_back(hEff->GetBinContent(ix,iy), hRate->GetBinContent(ix,iy));
    std::sort(pts.begin(), pts.end()); // sort ascending by eff
    // Pareto frontier: traverse from high eff to low, keep only improving (lower) rate points
    std::vector<std::pair<double,double>> frontier;
    double minRate = 1e18;
    for (int i = (int)pts.size()-1; i >= 0; --i) {
        if (pts[i].second < minRate) {
            minRate = pts[i].second;
            frontier.emplace_back(pts[i].first, pts[i].second);
        }
    }
    std::reverse(frontier.begin(), frontier.end());
    TGraph* g = new TGraph((int)frontier.size());
    for (int i = 0; i < (int)frontier.size(); ++i)
        g->SetPoint(i, frontier[i].first, frontier[i].second);
    return g;
}

// -----------------------------------------------------------------------
// Each signal/background entry is a pair: .first = HERNTupler input ntuple,
// .second = MET emulator output (contains metTree + emulEventInfoTree only).
void analyze_files(std::vector<std::pair<std::string, std::string>> signalFiles,
                   std::vector<std::pair<std::string, std::string>> backgroundFiles,
                   std::vector<std::string> labels,
                   std::string outputDir,
                   std::string signalName = "") {

    gSystem->mkdir(outputDir.c_str(), true);

    // Per-file histogram vectors for multi-file overlays
    std::vector<TH1F*> sig_h_TotalMET_vec,  back_h_TotalMET_vec;
    std::vector<TH1F*> sig_h_TotalMETX_vec, back_h_TotalMETX_vec;
    std::vector<TH1F*> sig_h_TotalMETY_vec, back_h_TotalMETY_vec;
    std::vector<TH1F*> sig_h_TowerMet_vec,  back_h_TowerMet_vec;
    std::vector<TH1F*> sig_h_JetMet_vec,    back_h_JetMet_vec;
    std::vector<TH1F*> sig_h_SumET_vec,     back_h_SumET_vec;
    std::vector<TH1F*> sig_h_SumJetET_vec,  back_h_SumJetET_vec;
    std::vector<TH1F*> sig_h_SumTowerET_vec,back_h_SumTowerET_vec;
    std::vector<TH1F*> sig_h_gMET_vec,      back_h_gMET_vec;
    std::vector<TH1F*> sig_h_gMET_NC_vec,   back_h_gMET_NC_vec;
    std::vector<TH1F*> sig_h_gMET_Rms_vec,  back_h_gMET_Rms_vec;
    std::vector<TH1F*> sig_h_jMET_vec,      back_h_jMET_vec;
    std::vector<TH1F*> sig_h_metTruthNonInt_vec, back_h_metTruthNonInt_vec;

    // Weighted background histograms for rate vs threshold
    std::vector<TH1F*> back_hw_TotalMET_vec;
    std::vector<TH1F*> back_hw_gMET_vec;
    std::vector<TH1F*> back_hw_gMET_NC_vec;
    std::vector<TH1F*> back_hw_gMET_Rms_vec;
    std::vector<TH1F*> back_hw_jMET_vec;
    std::vector<TH1F*> back_hw_JetMET_vec;
    std::vector<TH1F*> back_hw_TowerMET_vec;

    // 80 kHz efficiency histograms and thresholds per file for multi-file turn-on comparison
    std::vector<TH1F*> eff_gMET_80kHz_vec;
    std::vector<TH1F*> eff_gMET_NC_80kHz_vec;
    std::vector<TH1F*> eff_gMET_Rms_80kHz_vec;
    std::vector<TH1F*> eff_jMET_80kHz_vec;
    std::vector<TH1F*> eff_JetMET_80kHz_vec;
    std::vector<TH1F*> eff_TowerMET_80kHz_vec;
    std::vector<TH1F*> eff_TotalMET_80kHz_vec;
    std::vector<double> thr_gMET_80kHz_vec;
    std::vector<double> thr_gMET_NC_80kHz_vec;
    std::vector<double> thr_gMET_Rms_80kHz_vec;
    std::vector<double> thr_jMET_80kHz_vec;
    std::vector<double> thr_JetMET_80kHz_vec;
    std::vector<double> thr_TowerMET_80kHz_vec;
    std::vector<double> thr_TotalMET_80kHz_vec;
    std::vector<TH1F*> sig_h_metTruthNonInt_unscaled_vec; // unscaled clone for truth overlay

    // Resimulated gFEX MET histogram vectors (populated only when hasGFexSimMET)
    std::vector<TH1F*> sig_h_gMET_JwoJSim_vec, back_h_gMET_JwoJSim_vec;
    std::vector<TH1F*> sig_h_gMET_NCSim_vec,   back_h_gMET_NCSim_vec;
    std::vector<TH1F*> sig_h_gMET_RmsSim_vec,  back_h_gMET_RmsSim_vec;
    std::vector<TH1F*> back_hw_gMET_JwoJSim_vec;
    std::vector<TH1F*> back_hw_gMET_NCSim_vec;
    std::vector<TH1F*> back_hw_gMET_RmsSim_vec;
    std::vector<TH1F*> eff_gMET_JwoJSim_80kHz_vec;
    std::vector<TH1F*> eff_gMET_NCSim_80kHz_vec;
    std::vector<TH1F*> eff_gMET_RmsSim_80kHz_vec;
    std::vector<double> thr_gMET_JwoJSim_80kHz_vec;
    std::vector<double> thr_gMET_NCSim_80kHz_vec;
    std::vector<double> thr_gMET_RmsSim_80kHz_vec;

    bool hasSumJetET   = false;
    bool hasSumTowerET = false;

    for (unsigned int fileIt = 0; fileIt < signalFiles.size(); fileIt++) {
        std::cout << "Processing file " << fileIt << ": " << labels[fileIt] << "\n";
        // HERNTupler input ntuples: provide all upstream trees (gFEX MET, truth, gepCellsTowers jets, etc.)
        TFile* sigF  = TFile::Open(signalFiles[fileIt].first.c_str(),     "READ");
        TFile* backF = TFile::Open(backgroundFiles[fileIt].first.c_str(), "READ");
        // MET emulator outputs: provide metTree (and emulEventInfoTree for ordering validation)
        TFile* sigEmu  = TFile::Open(signalFiles[fileIt].second.c_str(),     "READ");
        TFile* backEmu = TFile::Open(backgroundFiles[fileIt].second.c_str(), "READ");
        if (!sigF  || sigF->IsZombie())  { std::cerr << "Cannot open " << signalFiles[fileIt].first      << "\n"; continue; }
        if (!backF || backF->IsZombie()) { std::cerr << "Cannot open " << backgroundFiles[fileIt].first  << "\n"; continue; }
        if (!sigEmu  || sigEmu->IsZombie())  { std::cerr << "Cannot open " << signalFiles[fileIt].second     << "\n"; continue; }
        if (!backEmu || backEmu->IsZombie()) { std::cerr << "Cannot open " << backgroundFiles[fileIt].second << "\n"; continue; }
        // --- TTrees ---
        TTree* metTreeSig              = (TTree*)sigEmu->Get("metTree");
        TTree* metTreeBack             = (TTree*)backEmu->Get("metTree");
        TTree* gFexMETTreeSig          = (TTree*)sigF->Get("gFexMETTree");
        TTree* gFexMETTreeBack         = (TTree*)backF->Get("gFexMETTree");
        TTree* gFexMETNoiseCutTreeSig  = (TTree*)sigF->Get("gFexMETNoiseCutTree");
        TTree* gFexMETNoiseCutTreeBack = (TTree*)backF->Get("gFexMETNoiseCutTree");
        TTree* gFexMETRmsTreeSig       = (TTree*)sigF->Get("gFexMETRmsTree");
        TTree* gFexMETRmsTreeBack      = (TTree*)backF->Get("gFexMETRmsTree");
        TTree* jFexMETTreeSig          = (TTree*)sigF->Get("jFexMETTree");
        TTree* jFexMETTreeBack         = (TTree*)backF->Get("jFexMETTree");
        TTree* metTruthTreeSig         = (TTree*)sigF->Get("metTruthTree");
        TTree* metTruthTreeBack        = (TTree*)backF->Get("metTruthTree");
        TTree* coreEMTopoTreeSig       = (TTree*)sigF->Get("metCoreAntiKt4EMTopoTree");
        TTree* coreEMTopoTreeBack      = (TTree*)backF->Get("metCoreAntiKt4EMTopoTree");
        TTree* coreEMPFlowTreeSig      = (TTree*)sigF->Get("metCoreAntiKt4EMPFlowTree");
        TTree* coreEMPFlowTreeBack     = (TTree*)backF->Get("metCoreAntiKt4EMPFlowTree");
        TTree* eventInfoTreeBack       = (TTree*)backF->Get("eventInfoTree");

        // --- Resimulated gFEX MET trees (only present in resim ntuples) ---
        TTree* gFexMETJwoJSimTreeSig      = (TTree*)sigF->Get("gFexMETJwoJSimTree");
        TTree* gFexMETJwoJSimTreeBack     = (TTree*)backF->Get("gFexMETJwoJSimTree");
        TTree* gFexMETNoiseCutSimTreeSig  = (TTree*)sigF->Get("gFexMETNoiseCutSimTree");
        TTree* gFexMETNoiseCutSimTreeBack = (TTree*)backF->Get("gFexMETNoiseCutSimTree");
        TTree* gFexMETRmsSimTreeSig       = (TTree*)sigF->Get("gFexMETRmsSimTree");
        TTree* gFexMETRmsSimTreeBack      = (TTree*)backF->Get("gFexMETRmsSimTree");
        bool hasGFexSimMET = (gFexMETJwoJSimTreeSig      != nullptr) &&
                             (gFexMETJwoJSimTreeBack     != nullptr) &&
                             (gFexMETNoiseCutSimTreeSig  != nullptr) &&
                             (gFexMETNoiseCutSimTreeBack != nullptr) &&
                             (gFexMETRmsSimTreeSig       != nullptr) &&
                             (gFexMETRmsSimTreeBack      != nullptr);
        if (hasGFexSimMET)
            std::cout << "  Resimulated gFEX MET trees found\n";
        else
            std::cout << "  Resimulated gFEX MET trees not found — Sim MET plots skipped\n";

        // --- WTA-cone GEP jet tree (gepCellsTowers, SK or non-SK based on config label) ---
        bool isSKConfig = (labels[fileIt].find("NoSK") == std::string::npos) &&
                          (labels[fileIt].find("SK")   != std::string::npos);
        // Tower MET is meaningless when Overlap Removal is on — drop it from GEP algo-comparison plots.
        bool hasOverlapRemoval = signalFiles[fileIt].second.find("_OR_") != std::string::npos;
        TTree* gepWTAConeCellsTowersJetsTreeSig  = isSKConfig
            ? (TTree*)sigF->Get("gepWTAConeCellsTowersSKJetsTree")
            : (TTree*)sigF->Get("gepWTAConeCellsTowersJetsTree");
        TTree* gepWTAConeCellsTowersJetsTreeBack = isSKConfig
            ? (TTree*)backF->Get("gepWTAConeCellsTowersSKJetsTree")
            : (TTree*)backF->Get("gepWTAConeCellsTowersJetsTree");
        bool hasWTAConeJets = (gepWTAConeCellsTowersJetsTreeSig  != nullptr) &&
                              (gepWTAConeCellsTowersJetsTreeBack != nullptr);
        if (hasWTAConeJets)
            std::cout << "  Using jet tree: "
                      << (isSKConfig ? "gepWTAConeCellsTowersSKJetsTree" : "gepWTAConeCellsTowersJetsTree") << "\n";
        else
            std::cout << "  WTA-cone jet tree not found — jet-level event-property plots skipped\n";

        // --- Truth AntiKt4 WZ-dressed jet tree (optional) ---
        TTree* truthAntiKt4WZDressedJetsTreeSig  = (TTree*)sigF->Get("truthAntiKt4TruthDressedWZJets");
        TTree* truthAntiKt4WZDressedJetsTreeBack = (TTree*)backF->Get("truthAntiKt4TruthDressedWZJets");
        bool hasTruthAntiKt4WZDressed = (truthAntiKt4WZDressedJetsTreeSig  != nullptr) &&
                                        (truthAntiKt4WZDressedJetsTreeBack != nullptr);

        // --- In-time pileup truth jet tree (optional) ---
        TTree* inTimeAntiKt4TruthJetsTreeSig  = (TTree*)sigF->Get("inTimeAntiKt4TruthJetsTree");
        TTree* inTimeAntiKt4TruthJetsTreeBack = (TTree*)backF->Get("inTimeAntiKt4TruthJetsTree");
        bool hasInTimeAntiKt4TruthJets = (inTimeAntiKt4TruthJetsTreeSig  != nullptr) &&
                                         (inTimeAntiKt4TruthJetsTreeBack != nullptr);

        // --- Branch variables: emulated MET ---
        double sig_TotalMET,  sig_TotalMETX,  sig_TotalMETY;
        double sig_TowerMet,  sig_JetMet,     sig_SumET;
        double back_TotalMET, back_TotalMETX, back_TotalMETY;
        double back_TowerMet, back_JetMet,    back_SumET;

        // --- Branch variables: GEP scalar sums (optional — guard if absent) ---
        double sig_SumJetET  = 0.0, sig_SumTowerET  = 0.0;
        double back_SumJetET = 0.0, back_SumTowerET = 0.0;
        hasSumJetET   = metTreeSig->FindBranch("SumJetET")   != nullptr;
        hasSumTowerET = metTreeSig->FindBranch("SumTowerET") != nullptr;
        if (hasSumJetET) {
            metTreeSig->SetBranchAddress("SumJetET",  &sig_SumJetET);
            metTreeBack->SetBranchAddress("SumJetET", &back_SumJetET);
        }
        if (hasSumTowerET) {
            metTreeSig->SetBranchAddress("SumTowerET",  &sig_SumTowerET);
            metTreeBack->SetBranchAddress("SumTowerET", &back_SumTowerET);
        }

        metTreeSig->SetBranchAddress("TotalMET",  &sig_TotalMET);
        metTreeSig->SetBranchAddress("TotalMETX", &sig_TotalMETX);
        metTreeSig->SetBranchAddress("TotalMETY", &sig_TotalMETY);
        metTreeSig->SetBranchAddress("TowerMet",  &sig_TowerMet);
        metTreeSig->SetBranchAddress("JetMet",    &sig_JetMet);
        metTreeSig->SetBranchAddress("SumET",     &sig_SumET);

        metTreeBack->SetBranchAddress("TotalMET",  &back_TotalMET);
        metTreeBack->SetBranchAddress("TotalMETX", &back_TotalMETX);
        metTreeBack->SetBranchAddress("TotalMETY", &back_TotalMETY);
        metTreeBack->SetBranchAddress("TowerMet",  &back_TowerMet);
        metTreeBack->SetBranchAddress("JetMet",    &back_JetMet);
        metTreeBack->SetBranchAddress("SumET",     &back_SumET);

        // --- Branch variables: GEP MET X/Y components (optional — guard if absent) ---
        double sig_JetMetX  = 0.0, sig_JetMetY  = 0.0;
        double back_JetMetX = 0.0, back_JetMetY = 0.0;
        double sig_TowerMetX  = 0.0, sig_TowerMetY  = 0.0;
        double back_TowerMetX = 0.0, back_TowerMetY = 0.0;
        bool hasJetMetXY   = metTreeSig->FindBranch("JetMetX")   != nullptr;
        bool hasTowerMetXY = metTreeSig->FindBranch("TowerMetX") != nullptr;
        if (hasJetMetXY) {
            metTreeSig->SetBranchAddress("JetMetX",  &sig_JetMetX);
            metTreeSig->SetBranchAddress("JetMetY",  &sig_JetMetY);
            metTreeBack->SetBranchAddress("JetMetX", &back_JetMetX);
            metTreeBack->SetBranchAddress("JetMetY", &back_JetMetY);
        }
        if (hasTowerMetXY) {
            metTreeSig->SetBranchAddress("TowerMetX",  &sig_TowerMetX);
            metTreeSig->SetBranchAddress("TowerMetY",  &sig_TowerMetY);
            metTreeBack->SetBranchAddress("TowerMetX", &back_TowerMetX);
            metTreeBack->SetBranchAddress("TowerMetY", &back_TowerMetY);
        }

        // --- Branch variables: gFEX MET ---
        double sig_gMET, sig_gMETX, sig_gMETY, sig_gSumET;
        double back_gMET, back_gMETX, back_gMETY, back_gSumET;

        gFexMETTreeSig->SetBranchAddress("gMET",   &sig_gMET);
        gFexMETTreeSig->SetBranchAddress("gMETX",  &sig_gMETX);
        gFexMETTreeSig->SetBranchAddress("gMETY",  &sig_gMETY);
        gFexMETTreeSig->SetBranchAddress("gSumET", &sig_gSumET);

        gFexMETTreeBack->SetBranchAddress("gMET",   &back_gMET);
        gFexMETTreeBack->SetBranchAddress("gMETX",  &back_gMETX);
        gFexMETTreeBack->SetBranchAddress("gMETY",  &back_gMETY);
        gFexMETTreeBack->SetBranchAddress("gSumET", &back_gSumET);

        // --- Branch variables: gFEX MET NoiseCut ---
        double sig_gMET_NC = 0.0, sig_gMETX_NC = 0.0, sig_gMETY_NC = 0.0, sig_gSumET_NC = 0.0;
        double back_gMET_NC = 0.0, back_gMETX_NC = 0.0, back_gMETY_NC = 0.0, back_gSumET_NC = 0.0;

        gFexMETNoiseCutTreeSig->SetBranchAddress("gMET",   &sig_gMET_NC);
        gFexMETNoiseCutTreeSig->SetBranchAddress("gMETX",  &sig_gMETX_NC);
        gFexMETNoiseCutTreeSig->SetBranchAddress("gMETY",  &sig_gMETY_NC);
        gFexMETNoiseCutTreeSig->SetBranchAddress("gSumET", &sig_gSumET_NC);
        gFexMETNoiseCutTreeBack->SetBranchAddress("gMET",   &back_gMET_NC);
        gFexMETNoiseCutTreeBack->SetBranchAddress("gMETX",  &back_gMETX_NC);
        gFexMETNoiseCutTreeBack->SetBranchAddress("gMETY",  &back_gMETY_NC);
        gFexMETNoiseCutTreeBack->SetBranchAddress("gSumET", &back_gSumET_NC);

        // --- Branch variables: gFEX MET Rms ---
        double sig_gMET_Rms = 0.0, sig_gMETX_Rms = 0.0, sig_gMETY_Rms = 0.0, sig_gSumET_Rms = 0.0;
        double back_gMET_Rms = 0.0, back_gMETX_Rms = 0.0, back_gMETY_Rms = 0.0, back_gSumET_Rms = 0.0;

        gFexMETRmsTreeSig->SetBranchAddress("gMET",   &sig_gMET_Rms);
        gFexMETRmsTreeSig->SetBranchAddress("gMETX",  &sig_gMETX_Rms);
        gFexMETRmsTreeSig->SetBranchAddress("gMETY",  &sig_gMETY_Rms);
        gFexMETRmsTreeSig->SetBranchAddress("gSumET", &sig_gSumET_Rms);
        gFexMETRmsTreeBack->SetBranchAddress("gMET",   &back_gMET_Rms);
        gFexMETRmsTreeBack->SetBranchAddress("gMETX",  &back_gMETX_Rms);
        gFexMETRmsTreeBack->SetBranchAddress("gMETY",  &back_gMETY_Rms);
        gFexMETRmsTreeBack->SetBranchAddress("gSumET", &back_gSumET_Rms);

        // --- Branch variables: jFEX MET (magnitude only) ---
        double sig_jMET = 0.0, back_jMET = 0.0;
        jFexMETTreeSig->SetBranchAddress("jMET",  &sig_jMET);
        jFexMETTreeBack->SetBranchAddress("jMET", &back_jMET);

        // --- Branch variables: resimulated gFEX MET (no X/Y components) ---
        double sig_gMET_JwoJSim = 0.0, sig_gSumET_JwoJSim = 0.0;
        double back_gMET_JwoJSim = 0.0, back_gSumET_JwoJSim = 0.0;
        double sig_gMET_NCSim = 0.0, sig_gSumET_NCSim = 0.0;
        double back_gMET_NCSim = 0.0, back_gSumET_NCSim = 0.0;
        double sig_gMET_RmsSim = 0.0, sig_gSumET_RmsSim = 0.0;
        double back_gMET_RmsSim = 0.0, back_gSumET_RmsSim = 0.0;
        if (hasGFexSimMET) {
            gFexMETJwoJSimTreeSig->SetBranchAddress("gMET",   &sig_gMET_JwoJSim);
            gFexMETJwoJSimTreeSig->SetBranchAddress("gSumET", &sig_gSumET_JwoJSim);
            gFexMETJwoJSimTreeBack->SetBranchAddress("gMET",   &back_gMET_JwoJSim);
            gFexMETJwoJSimTreeBack->SetBranchAddress("gSumET", &back_gSumET_JwoJSim);
            gFexMETNoiseCutSimTreeSig->SetBranchAddress("gMET",   &sig_gMET_NCSim);
            gFexMETNoiseCutSimTreeSig->SetBranchAddress("gSumET", &sig_gSumET_NCSim);
            gFexMETNoiseCutSimTreeBack->SetBranchAddress("gMET",   &back_gMET_NCSim);
            gFexMETNoiseCutSimTreeBack->SetBranchAddress("gSumET", &back_gSumET_NCSim);
            gFexMETRmsSimTreeSig->SetBranchAddress("gMET",   &sig_gMET_RmsSim);
            gFexMETRmsSimTreeSig->SetBranchAddress("gSumET", &sig_gSumET_RmsSim);
            gFexMETRmsSimTreeBack->SetBranchAddress("gMET",   &back_gMET_RmsSim);
            gFexMETRmsSimTreeBack->SetBranchAddress("gSumET", &back_gSumET_RmsSim);
        }

        // --- Branch variables: truth MET ---
        double sig_metTruthNonInt, sig_metTruthInt, sig_metTruthIntOut;
        double back_metTruthNonInt, back_metTruthInt, back_metTruthIntOut;

        metTruthTreeSig->SetBranchAddress("metTruthNonInt", &sig_metTruthNonInt);
        metTruthTreeSig->SetBranchAddress("metTruthInt",    &sig_metTruthInt);
        metTruthTreeSig->SetBranchAddress("metTruthIntOut", &sig_metTruthIntOut);
        metTruthTreeBack->SetBranchAddress("metTruthNonInt", &back_metTruthNonInt);
        metTruthTreeBack->SetBranchAddress("metTruthInt",    &back_metTruthInt);
        metTruthTreeBack->SetBranchAddress("metTruthIntOut", &back_metTruthIntOut);

        // Truth NonInt MET X/Y components (optional — guard if absent)
        double sig_metTruthNonIntX  = 0.0, sig_metTruthNonIntY  = 0.0;
        double back_metTruthNonIntX = 0.0, back_metTruthNonIntY = 0.0;
        bool hasTruthNonIntXY = metTruthTreeSig->FindBranch("metTruthNonIntX") != nullptr;
        if (hasTruthNonIntXY) {
            metTruthTreeSig->SetBranchAddress("metTruthNonIntX",  &sig_metTruthNonIntX);
            metTruthTreeSig->SetBranchAddress("metTruthNonIntY",  &sig_metTruthNonIntY);
            metTruthTreeBack->SetBranchAddress("metTruthNonIntX", &back_metTruthNonIntX);
            metTruthTreeBack->SetBranchAddress("metTruthNonIntY", &back_metTruthNonIntY);
        }

        // --- Branch variables: core MET EMTopo ---
        double sig_coreEMTopo_SoftClus_MET   = 0.0, sig_coreEMTopo_PVSoftTrk_MET  = 0.0, sig_coreEMTopo_SoftClusEM_MET = 0.0;
        double back_coreEMTopo_SoftClus_MET  = 0.0, back_coreEMTopo_PVSoftTrk_MET = 0.0, back_coreEMTopo_SoftClusEM_MET = 0.0;

        coreEMTopoTreeSig->SetBranchAddress("SoftClus_MET",   &sig_coreEMTopo_SoftClus_MET);
        coreEMTopoTreeSig->SetBranchAddress("PVSoftTrk_MET",  &sig_coreEMTopo_PVSoftTrk_MET);
        coreEMTopoTreeSig->SetBranchAddress("SoftClusEM_MET", &sig_coreEMTopo_SoftClusEM_MET);
        coreEMTopoTreeBack->SetBranchAddress("SoftClus_MET",   &back_coreEMTopo_SoftClus_MET);
        coreEMTopoTreeBack->SetBranchAddress("PVSoftTrk_MET",  &back_coreEMTopo_PVSoftTrk_MET);
        coreEMTopoTreeBack->SetBranchAddress("SoftClusEM_MET", &back_coreEMTopo_SoftClusEM_MET);

        // --- Branch variables: core MET EMPFlow ---
        double sig_coreEMPFlow_SoftClus_MET  = 0.0, sig_coreEMPFlow_PVSoftTrk_MET  = 0.0;
        double back_coreEMPFlow_SoftClus_MET = 0.0, back_coreEMPFlow_PVSoftTrk_MET = 0.0;

        coreEMPFlowTreeSig->SetBranchAddress("SoftClus_MET",  &sig_coreEMPFlow_SoftClus_MET);
        coreEMPFlowTreeSig->SetBranchAddress("PVSoftTrk_MET", &sig_coreEMPFlow_PVSoftTrk_MET);
        coreEMPFlowTreeBack->SetBranchAddress("SoftClus_MET",  &back_coreEMPFlow_SoftClus_MET);
        coreEMPFlowTreeBack->SetBranchAddress("PVSoftTrk_MET", &back_coreEMPFlow_PVSoftTrk_MET);

        // --- Background event weights + JZ slice index ---
        std::vector<double>* eventWeightsValuesBack = nullptr;
        int sampleJZSliceBack = -1;
        bool passHSTPValuesBack = true;
        eventInfoTreeBack->SetBranchAddress("eventWeights",  &eventWeightsValuesBack);
        eventInfoTreeBack->SetBranchAddress("sampleJZSlice", &sampleJZSliceBack);
        eventInfoTreeBack->SetBranchAddress("passHSTP",      &passHSTPValuesBack);

        // --- Branch variables: WTA-cone GEP jets (optional) ---
        std::vector<double>* gepWTAConeCellsTowersJetsEtValuesSig          = nullptr;
        std::vector<double>* gepWTAConeCellsTowersJetsEtaValuesSig         = nullptr;
        std::vector<double>* gepWTAConeCellsTowersJetsPhiValuesSig         = nullptr;
        std::vector<unsigned int>* gepWTAConeCellsTowersJetsNConstsSig     = nullptr;
        std::vector<double>* gepWTAConeCellsTowersJetsEtValuesBack         = nullptr;
        std::vector<double>* gepWTAConeCellsTowersJetsEtaValuesBack        = nullptr;
        std::vector<double>* gepWTAConeCellsTowersJetsPhiValuesBack        = nullptr;
        std::vector<unsigned int>* gepWTAConeCellsTowersJetsNConstsBack    = nullptr;
        if (hasWTAConeJets) {
            gepWTAConeCellsTowersJetsTreeSig->SetBranchAddress("Et",           &gepWTAConeCellsTowersJetsEtValuesSig);
            gepWTAConeCellsTowersJetsTreeSig->SetBranchAddress("Eta",          &gepWTAConeCellsTowersJetsEtaValuesSig);
            gepWTAConeCellsTowersJetsTreeSig->SetBranchAddress("Phi",          &gepWTAConeCellsTowersJetsPhiValuesSig);
            gepWTAConeCellsTowersJetsTreeSig->SetBranchAddress("NConstituents",&gepWTAConeCellsTowersJetsNConstsSig);
            gepWTAConeCellsTowersJetsTreeBack->SetBranchAddress("Et",           &gepWTAConeCellsTowersJetsEtValuesBack);
            gepWTAConeCellsTowersJetsTreeBack->SetBranchAddress("Eta",          &gepWTAConeCellsTowersJetsEtaValuesBack);
            gepWTAConeCellsTowersJetsTreeBack->SetBranchAddress("Phi",          &gepWTAConeCellsTowersJetsPhiValuesBack);
            gepWTAConeCellsTowersJetsTreeBack->SetBranchAddress("NConstituents",&gepWTAConeCellsTowersJetsNConstsBack);
        }

        // --- Branch variables: truth AntiKt4 WZ-dressed jets (optional) ---
        std::vector<double>* truthAntiKt4WZDressedJetsEtValuesSig  = nullptr;
        std::vector<double>* truthAntiKt4WZDressedJetsEtValuesBack = nullptr;
        if (hasTruthAntiKt4WZDressed) {
            truthAntiKt4WZDressedJetsTreeSig->SetBranchAddress("Et",  &truthAntiKt4WZDressedJetsEtValuesSig);
            truthAntiKt4WZDressedJetsTreeBack->SetBranchAddress("Et", &truthAntiKt4WZDressedJetsEtValuesBack);
        }

        // --- Branch variables: in-time pileup truth jets (optional) ---
        std::vector<double>* inTimeAntiKt4TruthJetsEtValuesSig  = nullptr;
        std::vector<double>* inTimeAntiKt4TruthJetsEtValuesBack = nullptr;
        if (hasInTimeAntiKt4TruthJets) {
            inTimeAntiKt4TruthJetsTreeSig->SetBranchAddress("Et",  &inTimeAntiKt4TruthJetsEtValuesSig);
            inTimeAntiKt4TruthJetsTreeBack->SetBranchAddress("Et", &inTimeAntiKt4TruthJetsEtValuesBack);
        }

        // --- Histograms ---
        std::string tag = std::to_string(fileIt);
        TH1F* sig_h_TotalMET   = new TH1F(("sig_h_TotalMET_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_TotalMET  = new TH1F(("back_h_TotalMET_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_TotalMETX  = new TH1F(("sig_h_TotalMETX_" +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_TotalMETX = new TH1F(("back_h_TotalMETX_"+tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_TotalMETY  = new TH1F(("sig_h_TotalMETY_" +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_TotalMETY = new TH1F(("back_h_TotalMETY_"+tag).c_str(), "", 80, -400, 400);
        // gFEX MET X/Y components
        TH1F* sig_h_gMETX      = new TH1F(("sig_h_gMETX_"     +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_gMETX     = new TH1F(("back_h_gMETX_"    +tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_gMETY      = new TH1F(("sig_h_gMETY_"     +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_gMETY     = new TH1F(("back_h_gMETY_"    +tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_gMETX_NC   = new TH1F(("sig_h_gMETX_NC_"  +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_gMETX_NC  = new TH1F(("back_h_gMETX_NC_" +tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_gMETY_NC   = new TH1F(("sig_h_gMETY_NC_"  +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_gMETY_NC  = new TH1F(("back_h_gMETY_NC_" +tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_gMETX_Rms  = new TH1F(("sig_h_gMETX_Rms_" +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_gMETX_Rms = new TH1F(("back_h_gMETX_Rms_"+tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_gMETY_Rms  = new TH1F(("sig_h_gMETY_Rms_" +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_gMETY_Rms = new TH1F(("back_h_gMETY_Rms_"+tag).c_str(), "", 80, -400, 400);
        // Truth NonInt MET X/Y components
        TH1F* sig_h_metTruthNonIntX  = new TH1F(("sig_h_TruthNonIntX_" +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_metTruthNonIntX = new TH1F(("back_h_TruthNonIntX_"+tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_metTruthNonIntY  = new TH1F(("sig_h_TruthNonIntY_" +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_metTruthNonIntY = new TH1F(("back_h_TruthNonIntY_"+tag).c_str(), "", 80, -400, 400);
        // GEP Jet and Tower MET X/Y components
        TH1F* sig_h_JetMetX    = new TH1F(("sig_h_JetMetX_"   +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_JetMetX   = new TH1F(("back_h_JetMetX_"  +tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_JetMetY    = new TH1F(("sig_h_JetMetY_"   +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_JetMetY   = new TH1F(("back_h_JetMetY_"  +tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_TowerMetX  = new TH1F(("sig_h_TowerMetX_" +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_TowerMetX = new TH1F(("back_h_TowerMetX_"+tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_TowerMetY  = new TH1F(("sig_h_TowerMetY_" +tag).c_str(), "", 80, -400, 400);
        TH1F* back_h_TowerMetY = new TH1F(("back_h_TowerMetY_"+tag).c_str(), "", 80, -400, 400);
        TH1F* sig_h_TowerMet   = new TH1F(("sig_h_TowerMet_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_TowerMet  = new TH1F(("back_h_TowerMet_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_JetMet     = new TH1F(("sig_h_JetMet_"    +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_JetMet    = new TH1F(("back_h_JetMet_"   +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_SumET      = new TH1F(("sig_h_SumET_"     +tag).c_str(), "", 52, 0,  1040);
        TH1F* back_h_SumET     = new TH1F(("back_h_SumET_"    +tag).c_str(), "", 52, 0,  1040);
        TH1F* sig_h_SumJetET   = new TH1F(("sig_h_SumJetET_"  +tag).c_str(), "", 52, 0,  1040);
        TH1F* back_h_SumJetET  = new TH1F(("back_h_SumJetET_" +tag).c_str(), "", 52, 0,  1040);
        TH1F* sig_h_SumTowerET  = new TH1F(("sig_h_SumTowerET_" +tag).c_str(), "", 52, 0,  1040);
        TH1F* back_h_SumTowerET = new TH1F(("back_h_SumTowerET_"+tag).c_str(), "", 52, 0,  1040);
        back_h_SumJetET->Sumw2(); back_h_SumTowerET->Sumw2();
        TH1F* sig_h_gMET       = new TH1F(("sig_h_gMET_"      +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_gMET      = new TH1F(("back_h_gMET_"     +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_gMET_NC    = new TH1F(("sig_h_gMET_NC_"   +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_gMET_NC   = new TH1F(("back_h_gMET_NC_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_gMET_Rms   = new TH1F(("sig_h_gMET_Rms_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_gMET_Rms  = new TH1F(("back_h_gMET_Rms_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_jMET       = new TH1F(("sig_h_jMET_"      +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_jMET      = new TH1F(("back_h_jMET_"     +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_metTruthNonInt       = new TH1F(("sig_h_TruthNonInt_"      +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_metTruthNonInt_coarse = new TH1F(("sig_h_TruthNonInt_coarse_"+tag).c_str(), "", 30, 0, 600); // 20 GeV bins for turn-on overlay
        TH1F* back_h_metTruthNonInt = new TH1F(("back_h_TruthNonInt_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_metTruthInt     = new TH1F(("sig_h_TruthInt_"    +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_metTruthInt    = new TH1F(("back_h_TruthInt_"   +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_metTruthIntOut  = new TH1F(("sig_h_TruthIntOut_" +tag).c_str(), "", 25, 0, 100);
        TH1F* back_h_metTruthIntOut = new TH1F(("back_h_TruthIntOut_"+tag).c_str(), "", 25, 0, 100);

        // Core MET term histograms
        TH1F* sig_h_coreEMTopo_SoftClus_MET    = new TH1F(("sig_coreEMTopo_SoftClus_"  +tag).c_str(), "", 70, 0, 350);
        TH1F* back_h_coreEMTopo_SoftClus_MET   = new TH1F(("back_coreEMTopo_SoftClus_" +tag).c_str(), "", 70, 0, 350);
        TH1F* sig_h_coreEMTopo_PVSoftTrk_MET   = new TH1F(("sig_coreEMTopo_PVSoftTrk_" +tag).c_str(), "", 20, 0, 100);
        TH1F* back_h_coreEMTopo_PVSoftTrk_MET  = new TH1F(("back_coreEMTopo_PVSoftTrk_"+tag).c_str(), "", 20, 0, 100);
        TH1F* sig_h_coreEMTopo_SoftClusEM_MET  = new TH1F(("sig_coreEMTopo_SoftClusEM_" +tag).c_str(), "", 40, 0, 200);
        TH1F* back_h_coreEMTopo_SoftClusEM_MET = new TH1F(("back_coreEMTopo_SoftClusEM_"+tag).c_str(), "", 40, 0, 200);
        TH1F* sig_h_coreEMPFlow_SoftClus_MET   = new TH1F(("sig_coreEMPFlow_SoftClus_"  +tag).c_str(), "", 20, 0, 100);
        TH1F* back_h_coreEMPFlow_SoftClus_MET  = new TH1F(("back_coreEMPFlow_SoftClus_" +tag).c_str(), "", 20, 0, 100);
        TH1F* sig_h_coreEMPFlow_PVSoftTrk_MET  = new TH1F(("sig_coreEMPFlow_PVSoftTrk_" +tag).c_str(), "", 10, 0, 50);
        TH1F* back_h_coreEMPFlow_PVSoftTrk_MET = new TH1F(("back_coreEMPFlow_PVSoftTrk_"+tag).c_str(), "", 10, 0, 50);

        // Weighted background histograms for rate plots
        TH1F* back_hw_TotalMET = new TH1F(("back_hw_TotalMET_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET     = new TH1F(("back_hw_gMET_"    +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_NC  = new TH1F(("back_hw_gMET_NC_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_Rms = new TH1F(("back_hw_gMET_Rms_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_jMET     = new TH1F(("back_hw_jMET_"    +tag).c_str(), "", nMETBins, metBinEdges);
        back_hw_jMET->Sumw2();
        TH1F* back_hw_JetMET   = new TH1F(("back_hw_JetMET_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_TowerMET = new TH1F(("back_hw_TowerMET_"+tag).c_str(), "", nMETBins, metBinEdges);

        // Resimulated gFEX MET distribution histograms and weighted background histograms
        TH1F* sig_h_gMET_JwoJSim  = new TH1F(("sig_h_gMET_JwoJSim_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_gMET_JwoJSim = new TH1F(("back_h_gMET_JwoJSim_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_gMET_NCSim    = new TH1F(("sig_h_gMET_NCSim_"   +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_gMET_NCSim   = new TH1F(("back_h_gMET_NCSim_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_gMET_RmsSim   = new TH1F(("sig_h_gMET_RmsSim_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_gMET_RmsSim  = new TH1F(("back_h_gMET_RmsSim_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_JwoJSim = new TH1F(("back_hw_gMET_JwoJSim_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_NCSim   = new TH1F(("back_hw_gMET_NCSim_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_RmsSim  = new TH1F(("back_hw_gMET_RmsSim_" +tag).c_str(), "", nMETBins, metBinEdges);
        back_hw_gMET_JwoJSim->Sumw2(); back_hw_gMET_NCSim->Sumw2(); back_hw_gMET_RmsSim->Sumw2();

        // 2D histograms for combined gFEX + GEP selection (x = gFEX MET, y = GEP MET), 10 GeV bins 0-600
        const int n2D = 60; const double lo2D = 0.0, hi2D = 600.0;
        TH2F* sig_h2_JwoJ_Jet   = new TH2F(("sig_h2_JwoJ_Jet_"  +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* sig_h2_JwoJ_Tower = new TH2F(("sig_h2_JwoJ_Tower_"+tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* sig_h2_NC_Jet     = new TH2F(("sig_h2_NC_Jet_"    +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* sig_h2_NC_Tower   = new TH2F(("sig_h2_NC_Tower_"  +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* sig_h2_Rms_Jet    = new TH2F(("sig_h2_Rms_Jet_"   +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* sig_h2_Rms_Tower  = new TH2F(("sig_h2_Rms_Tower_" +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* back_hw2_JwoJ_Jet   = new TH2F(("back_hw2_JwoJ_Jet_"  +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* back_hw2_JwoJ_Tower = new TH2F(("back_hw2_JwoJ_Tower_"+tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* back_hw2_NC_Jet     = new TH2F(("back_hw2_NC_Jet_"    +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* back_hw2_NC_Tower   = new TH2F(("back_hw2_NC_Tower_"  +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* back_hw2_Rms_Jet    = new TH2F(("back_hw2_Rms_Jet_"   +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* back_hw2_Rms_Tower  = new TH2F(("back_hw2_Rms_Tower_" +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        // Sumw2 so errors are correctly propagated through the 2D cumulative sum
        back_hw2_JwoJ_Jet->Sumw2();   back_hw2_JwoJ_Tower->Sumw2();
        back_hw2_NC_Jet->Sumw2();     back_hw2_NC_Tower->Sumw2();
        back_hw2_Rms_Jet->Sumw2();    back_hw2_Rms_Tower->Sumw2();

        // TOB MET vs truth NonInt MET calibration (x = truth, y = TOB) — signal and background
        TH2F* sig_h2_gJwoJ_TOBMet_vs_truthMET     = new TH2F(("sig_h2_gJwoJ_TOBMet_vs_truthMET_"    +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* sig_h2_gNC_TOBMet_vs_truthMET       = new TH2F(("sig_h2_gNC_TOBMet_vs_truthMET_"      +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* sig_h2_gRms_TOBMet_vs_truthMET      = new TH2F(("sig_h2_gRms_TOBMet_vs_truthMET_"     +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* sig_h2_JetMET_TOBMet_vs_truthMET    = new TH2F(("sig_h2_JetMET_TOBMet_vs_truthMET_"   +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        TH2F* sig_h2_TowerMET_TOBMet_vs_truthMET  = new TH2F(("sig_h2_TowerMET_TOBMet_vs_truthMET_" +tag).c_str(), "", n2D, lo2D, hi2D, n2D, lo2D, hi2D);
        // Variable bin edges for background: fine in region of interest, coarse outside
        // x (truth MET): 10 GeV bins [0,150], 50 GeV bins [150,600]
        const int nBkTruthX = 24;
        const Double_t bkTruthXedges[25] = {
              0,  10,  20,  30,  40,  50,  60,  70,  80,  90, 100,
            110, 120, 130, 140, 150,
            200, 250, 300, 350, 400, 450, 500, 550, 600
        };
        // y (TOB MET): 10 GeV bins [0,200], 50 GeV bins [200,600]
        const int nBkTOBY = 28;
        const Double_t bkTOBYedges[29] = {
              0,  10,  20,  30,  40,  50,  60,  70,  80,  90, 100,
            110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
            250, 300, 350, 400, 450, 500, 550, 600
        };
        TH2F* back_h2_gJwoJ_TOBMet_vs_truthMET    = new TH2F(("back_h2_gJwoJ_TOBMet_vs_truthMET_"   +tag).c_str(), "", nBkTruthX, bkTruthXedges, nBkTOBY, bkTOBYedges);
        TH2F* back_h2_gNC_TOBMet_vs_truthMET      = new TH2F(("back_h2_gNC_TOBMet_vs_truthMET_"     +tag).c_str(), "", nBkTruthX, bkTruthXedges, nBkTOBY, bkTOBYedges);
        TH2F* back_h2_gRms_TOBMet_vs_truthMET     = new TH2F(("back_h2_gRms_TOBMet_vs_truthMET_"    +tag).c_str(), "", nBkTruthX, bkTruthXedges, nBkTOBY, bkTOBYedges);
        TH2F* back_h2_JetMET_TOBMet_vs_truthMET   = new TH2F(("back_h2_JetMET_TOBMet_vs_truthMET_"  +tag).c_str(), "", nBkTruthX, bkTruthXedges, nBkTOBY, bkTOBYedges);
        TH2F* back_h2_TowerMET_TOBMet_vs_truthMET = new TH2F(("back_h2_TowerMET_TOBMet_vs_truthMET_"+tag).c_str(), "", nBkTruthX, bkTruthXedges, nBkTOBY, bkTOBYedges);
        back_h2_gJwoJ_TOBMet_vs_truthMET->Sumw2(); back_h2_gNC_TOBMet_vs_truthMET->Sumw2();
        back_h2_gRms_TOBMet_vs_truthMET->Sumw2();  back_h2_JetMET_TOBMet_vs_truthMET->Sumw2();
        back_h2_TowerMET_TOBMet_vs_truthMET->Sumw2();

        // Residual (truth - TOB) / truth  1D distributions — signal and background
        TH1F* sig_h1_gJwoJ_relResidual     = new TH1F(("sig_h1_gJwoJ_relResidual_"    +tag).c_str(), "", 100, -3.0, 3.0);
        TH1F* sig_h1_gNC_relResidual       = new TH1F(("sig_h1_gNC_relResidual_"      +tag).c_str(), "", 100, -3.0, 3.0);
        TH1F* sig_h1_gRms_relResidual      = new TH1F(("sig_h1_gRms_relResidual_"     +tag).c_str(), "", 100, -3.0, 3.0);
        TH1F* sig_h1_JetMET_relResidual    = new TH1F(("sig_h1_JetMET_relResidual_"   +tag).c_str(), "", 100, -3.0, 3.0);
        TH1F* sig_h1_TowerMET_relResidual  = new TH1F(("sig_h1_TowerMET_relResidual_" +tag).c_str(), "", 100, -3.0, 3.0);
        TH1F* back_h1_gJwoJ_relResidual    = new TH1F(("back_h1_gJwoJ_relResidual_"   +tag).c_str(), "", 100, -3.0, 3.0);
        TH1F* back_h1_gNC_relResidual      = new TH1F(("back_h1_gNC_relResidual_"     +tag).c_str(), "", 100, -3.0, 3.0);
        TH1F* back_h1_gRms_relResidual     = new TH1F(("back_h1_gRms_relResidual_"    +tag).c_str(), "", 100, -3.0, 3.0);
        TH1F* back_h1_JetMET_relResidual   = new TH1F(("back_h1_JetMET_relResidual_"  +tag).c_str(), "", 100, -3.0, 3.0);
        TH1F* back_h1_TowerMET_relResidual = new TH1F(("back_h1_TowerMET_relResidual_"+tag).c_str(), "", 100, -3.0, 3.0);
        back_h1_gJwoJ_relResidual->Sumw2(); back_h1_gNC_relResidual->Sumw2();
        back_h1_gRms_relResidual->Sumw2();  back_h1_JetMET_relResidual->Sumw2();
        back_h1_TowerMET_relResidual->Sumw2();

        // Residual vs truth NonInt MET 2D — signal and background
        const int nSumETbins2D = 20; const double hiSumET2D = 1000.0;
        TH2F* sig_h2_gJwoJ_relResidual_vs_truthMET     = new TH2F(("sig_h2_gJwoJ_relResidual_vs_truthMET_"    +tag).c_str(), "", n2D, lo2D, hi2D, 100, -3.0, 2.0);
        TH2F* sig_h2_gNC_relResidual_vs_truthMET       = new TH2F(("sig_h2_gNC_relResidual_vs_truthMET_"      +tag).c_str(), "", n2D, lo2D, hi2D, 100, -3.0, 2.0);
        TH2F* sig_h2_gRms_relResidual_vs_truthMET      = new TH2F(("sig_h2_gRms_relResidual_vs_truthMET_"     +tag).c_str(), "", n2D, lo2D, hi2D, 100, -3.0, 2.0);
        TH2F* sig_h2_JetMET_relResidual_vs_truthMET    = new TH2F(("sig_h2_JetMET_relResidual_vs_truthMET_"   +tag).c_str(), "", n2D, lo2D, hi2D, 100, -3.0, 2.0);
        TH2F* sig_h2_TowerMET_relResidual_vs_truthMET  = new TH2F(("sig_h2_TowerMET_relResidual_vs_truthMET_" +tag).c_str(), "", n2D, lo2D, hi2D, 100, -3.0, 2.0);
        TH2F* back_h2_gJwoJ_relResidual_vs_truthMET    = new TH2F(("back_h2_gJwoJ_relResidual_vs_truthMET_"   +tag).c_str(), "", nBkTruthX, bkTruthXedges, 100, -3.0, 2.0);
        TH2F* back_h2_gNC_relResidual_vs_truthMET      = new TH2F(("back_h2_gNC_relResidual_vs_truthMET_"     +tag).c_str(), "", nBkTruthX, bkTruthXedges, 100, -3.0, 2.0);
        TH2F* back_h2_gRms_relResidual_vs_truthMET     = new TH2F(("back_h2_gRms_relResidual_vs_truthMET_"    +tag).c_str(), "", nBkTruthX, bkTruthXedges, 100, -3.0, 2.0);
        TH2F* back_h2_JetMET_relResidual_vs_truthMET   = new TH2F(("back_h2_JetMET_relResidual_vs_truthMET_"  +tag).c_str(), "", nBkTruthX, bkTruthXedges, 100, -3.0, 2.0);
        TH2F* back_h2_TowerMET_relResidual_vs_truthMET = new TH2F(("back_h2_TowerMET_relResidual_vs_truthMET_"+tag).c_str(), "", nBkTruthX, bkTruthXedges, 100, -3.0, 2.0);
        back_h2_gJwoJ_relResidual_vs_truthMET->Sumw2(); back_h2_gNC_relResidual_vs_truthMET->Sumw2();
        back_h2_gRms_relResidual_vs_truthMET->Sumw2();  back_h2_JetMET_relResidual_vs_truthMET->Sumw2();
        back_h2_TowerMET_relResidual_vs_truthMET->Sumw2();

        // Residual vs TOB SumET 2D — signal and background
        TH2F* sig_h2_gJwoJ_relResidual_vs_sumET     = new TH2F(("sig_h2_gJwoJ_relResidual_vs_sumET_"    +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 50, -3.0, 2.0);
        TH2F* sig_h2_gNC_relResidual_vs_sumET       = new TH2F(("sig_h2_gNC_relResidual_vs_sumET_"      +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 50, -3.0, 2.0);
        TH2F* sig_h2_gRms_relResidual_vs_sumET      = new TH2F(("sig_h2_gRms_relResidual_vs_sumET_"     +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 50, -3.0, 2.0);
        TH2F* sig_h2_JetMET_relResidual_vs_sumET    = new TH2F(("sig_h2_JetMET_relResidual_vs_sumET_"   +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 50, -3.0, 2.0);
        TH2F* sig_h2_TowerMET_relResidual_vs_sumET  = new TH2F(("sig_h2_TowerMET_relResidual_vs_sumET_" +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 50, -3.0, 2.0);
        TH2F* back_h2_gJwoJ_relResidual_vs_sumET    = new TH2F(("back_h2_gJwoJ_relResidual_vs_sumET_"   +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 50, -3.0, 2.0);
        TH2F* back_h2_gNC_relResidual_vs_sumET      = new TH2F(("back_h2_gNC_relResidual_vs_sumET_"     +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 50, -3.0, 2.0);
        TH2F* back_h2_gRms_relResidual_vs_sumET     = new TH2F(("back_h2_gRms_relResidual_vs_sumET_"    +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 50, -3.0, 2.0);
        TH2F* back_h2_JetMET_relResidual_vs_sumET   = new TH2F(("back_h2_JetMET_relResidual_vs_sumET_"  +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 50, -3.0, 2.0);
        TH2F* back_h2_TowerMET_relResidual_vs_sumET = new TH2F(("back_h2_TowerMET_relResidual_vs_sumET_"+tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 50, -3.0, 2.0);
        back_h2_gJwoJ_relResidual_vs_sumET->Sumw2(); back_h2_gNC_relResidual_vs_sumET->Sumw2();
        back_h2_gRms_relResidual_vs_sumET->Sumw2();  back_h2_JetMET_relResidual_vs_sumET->Sumw2();
        back_h2_TowerMET_relResidual_vs_sumET->Sumw2();

        // Absolute residual (truth - TOB) [GeV] 1D — signal and background
        const double hiAbsRes = 300.0; // clamp range [GeV]
        TH1F* sig_h1_gJwoJ_absResidual     = new TH1F(("sig_h1_gJwoJ_absResidual_"    +tag).c_str(), "", 60, -hiAbsRes, hiAbsRes);
        TH1F* sig_h1_gNC_absResidual       = new TH1F(("sig_h1_gNC_absResidual_"      +tag).c_str(), "", 60, -hiAbsRes, hiAbsRes);
        TH1F* sig_h1_gRms_absResidual      = new TH1F(("sig_h1_gRms_absResidual_"     +tag).c_str(), "", 60, -hiAbsRes, hiAbsRes);
        TH1F* sig_h1_JetMET_absResidual    = new TH1F(("sig_h1_JetMET_absResidual_"   +tag).c_str(), "", 60, -hiAbsRes, hiAbsRes);
        TH1F* sig_h1_TowerMET_absResidual  = new TH1F(("sig_h1_TowerMET_absResidual_" +tag).c_str(), "", 60, -hiAbsRes, hiAbsRes);
        TH1F* back_h1_gJwoJ_absResidual    = new TH1F(("back_h1_gJwoJ_absResidual_"   +tag).c_str(), "", 60, -hiAbsRes, hiAbsRes);
        TH1F* back_h1_gNC_absResidual      = new TH1F(("back_h1_gNC_absResidual_"     +tag).c_str(), "", 60, -hiAbsRes, hiAbsRes);
        TH1F* back_h1_gRms_absResidual     = new TH1F(("back_h1_gRms_absResidual_"    +tag).c_str(), "", 60, -hiAbsRes, hiAbsRes);
        TH1F* back_h1_JetMET_absResidual   = new TH1F(("back_h1_JetMET_absResidual_"  +tag).c_str(), "", 60, -hiAbsRes, hiAbsRes);
        TH1F* back_h1_TowerMET_absResidual = new TH1F(("back_h1_TowerMET_absResidual_"+tag).c_str(), "", 60, -hiAbsRes, hiAbsRes);
        back_h1_gJwoJ_absResidual->Sumw2(); back_h1_gNC_absResidual->Sumw2();
        back_h1_gRms_absResidual->Sumw2();  back_h1_JetMET_absResidual->Sumw2();
        back_h1_TowerMET_absResidual->Sumw2();
        // Absolute residual vs truth NonInt MET 2D
        TH2F* sig_h2_gJwoJ_absResidual_vs_truthMET     = new TH2F(("sig_h2_gJwoJ_absResidual_vs_truthMET_"    +tag).c_str(), "", n2D, lo2D, hi2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* sig_h2_gNC_absResidual_vs_truthMET       = new TH2F(("sig_h2_gNC_absResidual_vs_truthMET_"      +tag).c_str(), "", n2D, lo2D, hi2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* sig_h2_gRms_absResidual_vs_truthMET      = new TH2F(("sig_h2_gRms_absResidual_vs_truthMET_"     +tag).c_str(), "", n2D, lo2D, hi2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* sig_h2_JetMET_absResidual_vs_truthMET    = new TH2F(("sig_h2_JetMET_absResidual_vs_truthMET_"   +tag).c_str(), "", n2D, lo2D, hi2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* sig_h2_TowerMET_absResidual_vs_truthMET  = new TH2F(("sig_h2_TowerMET_absResidual_vs_truthMET_" +tag).c_str(), "", n2D, lo2D, hi2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* back_h2_gJwoJ_absResidual_vs_truthMET    = new TH2F(("back_h2_gJwoJ_absResidual_vs_truthMET_"   +tag).c_str(), "", nBkTruthX, bkTruthXedges, 60, -hiAbsRes, hiAbsRes);
        TH2F* back_h2_gNC_absResidual_vs_truthMET      = new TH2F(("back_h2_gNC_absResidual_vs_truthMET_"     +tag).c_str(), "", nBkTruthX, bkTruthXedges, 60, -hiAbsRes, hiAbsRes);
        TH2F* back_h2_gRms_absResidual_vs_truthMET     = new TH2F(("back_h2_gRms_absResidual_vs_truthMET_"    +tag).c_str(), "", nBkTruthX, bkTruthXedges, 60, -hiAbsRes, hiAbsRes);
        TH2F* back_h2_JetMET_absResidual_vs_truthMET   = new TH2F(("back_h2_JetMET_absResidual_vs_truthMET_"  +tag).c_str(), "", nBkTruthX, bkTruthXedges, 60, -hiAbsRes, hiAbsRes);
        TH2F* back_h2_TowerMET_absResidual_vs_truthMET = new TH2F(("back_h2_TowerMET_absResidual_vs_truthMET_"+tag).c_str(), "", nBkTruthX, bkTruthXedges, 60, -hiAbsRes, hiAbsRes);
        back_h2_gJwoJ_absResidual_vs_truthMET->Sumw2(); back_h2_gNC_absResidual_vs_truthMET->Sumw2();
        back_h2_gRms_absResidual_vs_truthMET->Sumw2();  back_h2_JetMET_absResidual_vs_truthMET->Sumw2();
        back_h2_TowerMET_absResidual_vs_truthMET->Sumw2();
        // Absolute residual vs TOB SumET 2D
        TH2F* sig_h2_gJwoJ_absResidual_vs_sumET     = new TH2F(("sig_h2_gJwoJ_absResidual_vs_sumET_"    +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* sig_h2_gNC_absResidual_vs_sumET       = new TH2F(("sig_h2_gNC_absResidual_vs_sumET_"      +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* sig_h2_gRms_absResidual_vs_sumET      = new TH2F(("sig_h2_gRms_absResidual_vs_sumET_"     +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* sig_h2_JetMET_absResidual_vs_sumET    = new TH2F(("sig_h2_JetMET_absResidual_vs_sumET_"   +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* sig_h2_TowerMET_absResidual_vs_sumET  = new TH2F(("sig_h2_TowerMET_absResidual_vs_sumET_" +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* back_h2_gJwoJ_absResidual_vs_sumET    = new TH2F(("back_h2_gJwoJ_absResidual_vs_sumET_"   +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* back_h2_gNC_absResidual_vs_sumET      = new TH2F(("back_h2_gNC_absResidual_vs_sumET_"     +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* back_h2_gRms_absResidual_vs_sumET     = new TH2F(("back_h2_gRms_absResidual_vs_sumET_"    +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* back_h2_JetMET_absResidual_vs_sumET   = new TH2F(("back_h2_JetMET_absResidual_vs_sumET_"  +tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 60, -hiAbsRes, hiAbsRes);
        TH2F* back_h2_TowerMET_absResidual_vs_sumET = new TH2F(("back_h2_TowerMET_absResidual_vs_sumET_"+tag).c_str(), "", nSumETbins2D, 0.0, hiSumET2D, 60, -hiAbsRes, hiAbsRes);
        back_h2_gJwoJ_absResidual_vs_sumET->Sumw2(); back_h2_gNC_absResidual_vs_sumET->Sumw2();
        back_h2_gRms_absResidual_vs_sumET->Sumw2();  back_h2_JetMET_absResidual_vs_sumET->Sumw2();
        back_h2_TowerMET_absResidual_vs_sumET->Sumw2();

        // Turn-on histograms: denom (all signal) + 9 numerators (3 algos × 3 rates)
        TH1F* h_turnOn_denom             = new TH1F(("h_turnOn_denom_"             +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_20kHz    = new TH1F(("h_turnOn_num_gMET_20kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_40kHz    = new TH1F(("h_turnOn_num_gMET_40kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_80kHz    = new TH1F(("h_turnOn_num_gMET_80kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_JetMET_20kHz  = new TH1F(("h_turnOn_num_JetMET_20kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_JetMET_40kHz  = new TH1F(("h_turnOn_num_JetMET_40kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_JetMET_80kHz  = new TH1F(("h_turnOn_num_JetMET_80kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TowerMET_20kHz = new TH1F(("h_turnOn_num_TowerMET_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TowerMET_40kHz = new TH1F(("h_turnOn_num_TowerMET_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TowerMET_80kHz = new TH1F(("h_turnOn_num_TowerMET_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TotalMET_20kHz = new TH1F(("h_turnOn_num_TotalMET_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TotalMET_40kHz = new TH1F(("h_turnOn_num_TotalMET_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TotalMET_80kHz = new TH1F(("h_turnOn_num_TotalMET_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        // gFEX NoiseCut and Rms individual turn-on numerators
        TH1F* h_turnOn_num_gMET_NC_20kHz  = new TH1F(("h_turnOn_num_gMET_NC_20kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NC_40kHz  = new TH1F(("h_turnOn_num_gMET_NC_40kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NC_80kHz  = new TH1F(("h_turnOn_num_gMET_NC_80kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_Rms_20kHz = new TH1F(("h_turnOn_num_gMET_Rms_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_Rms_40kHz = new TH1F(("h_turnOn_num_gMET_Rms_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_Rms_80kHz = new TH1F(("h_turnOn_num_gMET_Rms_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        // jFEX MET turn-on numerators
        TH1F* h_turnOn_num_jMET_20kHz    = new TH1F(("h_turnOn_num_jMET_20kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_jMET_40kHz    = new TH1F(("h_turnOn_num_jMET_40kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_jMET_80kHz    = new TH1F(("h_turnOn_num_jMET_80kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        // Resimulated gFEX MET turn-on numerators
        TH1F* h_turnOn_num_gMET_JwoJSim_20kHz = new TH1F(("h_turnOn_num_gMET_JwoJSim_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_JwoJSim_40kHz = new TH1F(("h_turnOn_num_gMET_JwoJSim_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_JwoJSim_80kHz = new TH1F(("h_turnOn_num_gMET_JwoJSim_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NCSim_20kHz   = new TH1F(("h_turnOn_num_gMET_NCSim_20kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NCSim_40kHz   = new TH1F(("h_turnOn_num_gMET_NCSim_40kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NCSim_80kHz   = new TH1F(("h_turnOn_num_gMET_NCSim_80kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_RmsSim_20kHz  = new TH1F(("h_turnOn_num_gMET_RmsSim_20kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_RmsSim_40kHz  = new TH1F(("h_turnOn_num_gMET_RmsSim_40kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_RmsSim_80kHz  = new TH1F(("h_turnOn_num_gMET_RmsSim_80kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        // Combined gFEX+GEP turn-on numerators (6 combos × 3 rate targets = 18)
        // Declared here as nullptr; created after the 2D scan finds the best thresholds
        TH1F* h_turnOn_num_combo_JwoJ_Jet_20kHz   = new TH1F(("h_turnOn_num_combo_JwoJ_Jet_20kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Jet_40kHz   = new TH1F(("h_turnOn_num_combo_JwoJ_Jet_40kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Jet_80kHz   = new TH1F(("h_turnOn_num_combo_JwoJ_Jet_80kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Tower_20kHz = new TH1F(("h_turnOn_num_combo_JwoJ_Tower_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Tower_40kHz = new TH1F(("h_turnOn_num_combo_JwoJ_Tower_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Tower_80kHz = new TH1F(("h_turnOn_num_combo_JwoJ_Tower_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Jet_20kHz     = new TH1F(("h_turnOn_num_combo_NC_Jet_20kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Jet_40kHz     = new TH1F(("h_turnOn_num_combo_NC_Jet_40kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Jet_80kHz     = new TH1F(("h_turnOn_num_combo_NC_Jet_80kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Tower_20kHz   = new TH1F(("h_turnOn_num_combo_NC_Tower_20kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Tower_40kHz   = new TH1F(("h_turnOn_num_combo_NC_Tower_40kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Tower_80kHz   = new TH1F(("h_turnOn_num_combo_NC_Tower_80kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Jet_20kHz    = new TH1F(("h_turnOn_num_combo_Rms_Jet_20kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Jet_40kHz    = new TH1F(("h_turnOn_num_combo_Rms_Jet_40kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Jet_80kHz    = new TH1F(("h_turnOn_num_combo_Rms_Jet_80kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Tower_20kHz  = new TH1F(("h_turnOn_num_combo_Rms_Tower_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Tower_40kHz  = new TH1F(("h_turnOn_num_combo_Rms_Tower_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Tower_80kHz  = new TH1F(("h_turnOn_num_combo_Rms_Tower_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);

        // Per-JZ-slice MET histograms for background (gFEX MET, GEP Jet MET, GEP Tower MET)
        TH1F* back_h_gMET_jz[nJZSlices_];
        TH1F* back_h_gMET_NC_jz[nJZSlices_];
        TH1F* back_h_gMET_Rms_jz[nJZSlices_];
        TH1F* back_h_jMET_jz[nJZSlices_];
        TH1F* back_h_JetMET_jz[nJZSlices_];
        TH1F* back_h_TowerMET_jz[nJZSlices_];
        for (unsigned int jz = 0; jz < nJZSlices_; ++jz) {
            std::string jztag = tag + "_jz" + std::to_string(jz);
            back_h_gMET_jz[jz]     = new TH1F(("back_h_gMET_jz_"    +jztag).c_str(), ("gFEX JwoJ MET, JZ"        +std::to_string(jz)+";gFEX JwoJ MET [GeV];Events").c_str(),        nMETBins, metBinEdges);
            back_h_gMET_NC_jz[jz]  = new TH1F(("back_h_gMET_NC_jz_" +jztag).c_str(), ("gFEX NoiseCut MET, JZ"    +std::to_string(jz)+";gFEX NoiseCut MET [GeV];Events").c_str(),    nMETBins, metBinEdges);
            back_h_gMET_Rms_jz[jz] = new TH1F(("back_h_gMET_Rms_jz_"+jztag).c_str(), ("gFEX RMS MET, JZ"         +std::to_string(jz)+";gFEX RMS MET [GeV];Events").c_str(),         nMETBins, metBinEdges);
            back_h_jMET_jz[jz]     = new TH1F(("back_h_jMET_jz_"    +jztag).c_str(), ("jFEX MET, JZ"             +std::to_string(jz)+";jFEX MET [GeV];Events").c_str(),             nMETBins, metBinEdges);
            back_h_JetMET_jz[jz]   = new TH1F(("back_h_JetMET_jz_"  +jztag).c_str(), ("GEP Jet MET, JZ"          +std::to_string(jz)+";GEP Jet MET [GeV];Events").c_str(),          nMETBins, metBinEdges);
            back_h_TowerMET_jz[jz] = new TH1F(("back_h_TowerMET_jz_"+jztag).c_str(), ("GEP Tower MET, JZ"        +std::to_string(jz)+";GEP Tower MET [GeV];Events").c_str(),        nMETBins, metBinEdges);
        }

        // --- Event-property histograms per 80 kHz selection (signal & weighted background) ---
        // Index: 0 = Inclusive, 1–5 = individual 80 kHz trigger selections
        const int nSel80 = 6;
        const char* selShortNames[nSel80] = {
            "Incl", "JwoJ_80kHz", "NC_80kHz", "Rms_80kHz", "Jet_80kHz", "Tower_80kHz"
        };
        TH1F* sig_sel_truthMET[nSel80],    *back_sel_truthMET[nSel80];
        TH1F* sig_sel_SumET[nSel80],       *back_sel_SumET[nSel80];
        TH1F* sig_sel_gSumET[nSel80],      *back_sel_gSumET[nSel80];
        TH1F* sig_sel_METsig[nSel80],      *back_sel_METsig[nSel80];    // GEP MET / sqrt(SumET)
        TH1F* sig_sel_gMETsig[nSel80],     *back_sel_gMETsig[nSel80];   // gFEX MET / sqrt(gSumET)
        TH1F* sig_sel_dPhi_GEP_gFEX[nSel80],   *back_sel_dPhi_GEP_gFEX[nSel80];
        TH1F* sig_sel_dPhi_GEP_truth[nSel80],  *back_sel_dPhi_GEP_truth[nSel80];
        TH1F* sig_sel_dPhi_truth_TOB[nSel80],  *back_sel_dPhi_truth_TOB[nSel80];
        TH1F* sig_sel_nJets[nSel80],         *back_sel_nJets[nSel80];
        TH1F* sig_sel_nTruthJets[nSel80],   *back_sel_nTruthJets[nSel80];
        TH1F* sig_sel_nPileupJets[nSel80],  *back_sel_nPileupJets[nSel80];
        TH1F* sig_sel_jet1pt[nSel80],      *back_sel_jet1pt[nSel80];
        TH1F* sig_sel_jet2pt[nSel80],      *back_sel_jet2pt[nSel80];
        TH1F* sig_sel_jet1eta[nSel80],     *back_sel_jet1eta[nSel80];
        TH1F* sig_sel_dPhi_jet1_MET[nSel80],  *back_sel_dPhi_jet1_MET[nSel80];
        TH2F* sig_sel_phi2D_TOB_truth[nSel80], *back_sel_phi2D_TOB_truth[nSel80];
        for (int iSel = 0; iSel < nSel80; iSel++) {
            std::string st = tag + "_" + selShortNames[iSel];
            sig_sel_truthMET[iSel]         = new TH1F(("sig_sel_truthMET_"       +st).c_str(), "", nMETBins, metBinEdges);
            back_sel_truthMET[iSel]        = new TH1F(("back_sel_truthMET_"      +st).c_str(), "", nMETBins, metBinEdges);
            sig_sel_SumET[iSel]            = new TH1F(("sig_sel_SumET_"          +st).c_str(), "", 52, 0, 1040);
            back_sel_SumET[iSel]           = new TH1F(("back_sel_SumET_"         +st).c_str(), "", 52, 0, 1040);
            sig_sel_gSumET[iSel]           = new TH1F(("sig_sel_gSumET_"         +st).c_str(), "", 52, 0, 1040);
            back_sel_gSumET[iSel]          = new TH1F(("back_sel_gSumET_"        +st).c_str(), "", 52, 0, 1040);
            sig_sel_METsig[iSel]           = new TH1F(("sig_sel_METsig_"         +st).c_str(), "", 50, 0, 25);
            back_sel_METsig[iSel]          = new TH1F(("back_sel_METsig_"        +st).c_str(), "", 50, 0, 25);
            sig_sel_gMETsig[iSel]          = new TH1F(("sig_sel_gMETsig_"        +st).c_str(), "", 50, 0, 25);
            back_sel_gMETsig[iSel]         = new TH1F(("back_sel_gMETsig_"       +st).c_str(), "", 50, 0, 25);
            sig_sel_dPhi_GEP_gFEX[iSel]   = new TH1F(("sig_sel_dPhi_GEP_gFEX_" +st).c_str(), "", 32, 0, M_PI);
            back_sel_dPhi_GEP_gFEX[iSel]  = new TH1F(("back_sel_dPhi_GEP_gFEX_"+st).c_str(), "", 32, 0, M_PI);
            sig_sel_dPhi_GEP_truth[iSel]   = new TH1F(("sig_sel_dPhi_GEP_truth_" +st).c_str(), "", 32, 0, M_PI);
            back_sel_dPhi_GEP_truth[iSel]  = new TH1F(("back_sel_dPhi_GEP_truth_"+st).c_str(), "", 32, 0, M_PI);
            sig_sel_dPhi_truth_TOB[iSel]   = new TH1F(("sig_sel_dPhi_truth_TOB_" +st).c_str(), "", 32, 0, M_PI);
            back_sel_dPhi_truth_TOB[iSel]  = new TH1F(("back_sel_dPhi_truth_TOB_"+st).c_str(), "", 32, 0, M_PI);
            sig_sel_nJets[iSel]            = new TH1F(("sig_sel_nJets_"          +st).c_str(), "", 10, 0, 10);
            back_sel_nJets[iSel]           = new TH1F(("back_sel_nJets_"         +st).c_str(), "", 10, 0, 10);
            sig_sel_nTruthJets[iSel]       = new TH1F(("sig_sel_nTruthJets_"     +st).c_str(), "", 10, 0, 10);
            back_sel_nTruthJets[iSel]      = new TH1F(("back_sel_nTruthJets_"    +st).c_str(), "", 10, 0, 10);
            sig_sel_nPileupJets[iSel]      = new TH1F(("sig_sel_nPileupJets_"    +st).c_str(), "", 10, 0, 10);
            back_sel_nPileupJets[iSel]     = new TH1F(("back_sel_nPileupJets_"   +st).c_str(), "", 10, 0, 10);
            sig_sel_jet1pt[iSel]           = new TH1F(("sig_sel_jet1pt_"         +st).c_str(), "", 50, 0, 500);
            back_sel_jet1pt[iSel]          = new TH1F(("back_sel_jet1pt_"        +st).c_str(), "", 50, 0, 500);
            sig_sel_jet2pt[iSel]           = new TH1F(("sig_sel_jet2pt_"         +st).c_str(), "", 50, 0, 500);
            back_sel_jet2pt[iSel]          = new TH1F(("back_sel_jet2pt_"        +st).c_str(), "", 50, 0, 500);
            sig_sel_jet1eta[iSel]          = new TH1F(("sig_sel_jet1eta_"        +st).c_str(), "", 50, -5, 5);
            back_sel_jet1eta[iSel]         = new TH1F(("back_sel_jet1eta_"       +st).c_str(), "", 50, -5, 5);
            sig_sel_dPhi_jet1_MET[iSel]    = new TH1F(("sig_sel_dPhi_jet1_MET_" +st).c_str(), "", 32, 0, M_PI);
            back_sel_dPhi_jet1_MET[iSel]   = new TH1F(("back_sel_dPhi_jet1_MET_"+st).c_str(), "", 32, 0, M_PI);
            sig_sel_phi2D_TOB_truth[iSel]  = new TH2F(("sig_sel_phi2D_TOB_truth_" +st).c_str(), ";Truth MET #phi [rad];TOB MET #phi [rad]", 32, -M_PI, M_PI, 32, -M_PI, M_PI);
            back_sel_phi2D_TOB_truth[iSel] = new TH2F(("back_sel_phi2D_TOB_truth_"+st).c_str(), ";Truth MET #phi [rad];TOB MET #phi [rad]", 32, -M_PI, M_PI, 32, -M_PI, M_PI);
            back_sel_phi2D_TOB_truth[iSel]->Sumw2();
            back_sel_truthMET[iSel]->Sumw2();       back_sel_SumET[iSel]->Sumw2();
            back_sel_gSumET[iSel]->Sumw2();         back_sel_METsig[iSel]->Sumw2();
            back_sel_gMETsig[iSel]->Sumw2();        back_sel_dPhi_GEP_gFEX[iSel]->Sumw2();
            back_sel_dPhi_GEP_truth[iSel]->Sumw2(); back_sel_dPhi_truth_TOB[iSel]->Sumw2();
            back_sel_nJets[iSel]->Sumw2();
            back_sel_nTruthJets[iSel]->Sumw2();
            back_sel_nPileupJets[iSel]->Sumw2();
            back_sel_jet1pt[iSel]->Sumw2();         back_sel_jet2pt[iSel]->Sumw2();
            back_sel_jet1eta[iSel]->Sumw2();        back_sel_dPhi_jet1_MET[iSel]->Sumw2();
        }

        // Clamp a value to the last bin of h (avoids overflow for variable-bin histograms)
        auto clampVal = [](TH1F* h, double v) {
            double xmax = h->GetXaxis()->GetXmax();
            return v >= xmax ? xmax - 1e-9 : v;
        };

        // --- Signal event loop ---
        unsigned int nSig = metTreeSig->GetEntries();
        std::cout << "  Signal events: " << nSig << "\n";
        for (unsigned int iEvt = 0; iEvt < nSig; iEvt++) {
            metTreeSig->GetEntry(iEvt);
            gFexMETTreeSig->GetEntry(iEvt);
            gFexMETNoiseCutTreeSig->GetEntry(iEvt);
            gFexMETRmsTreeSig->GetEntry(iEvt);
            metTruthTreeSig->GetEntry(iEvt);
            // metTruthNonIntX/Y already points in the standard MET direction (sum of NonInt particle
            // momenta = direction of neutrinos). The DAOD note's constraint is consistent with
            // momentum conservation, with the Int/IntOut/IntMuon terms (not NonInt) carrying the sign flip.
            coreEMTopoTreeSig->GetEntry(iEvt);
            coreEMPFlowTreeSig->GetEntry(iEvt);

            sig_h_TotalMET->Fill(clampVal(sig_h_TotalMET, sig_TotalMET));
            sig_h_TotalMETX->Fill(sig_TotalMETX);
            sig_h_TotalMETY->Fill(sig_TotalMETY);
            sig_h_TowerMet->Fill(clampVal(sig_h_TowerMet, sig_TowerMet));
            sig_h_JetMet->Fill(clampVal(sig_h_JetMet, sig_JetMet));
            sig_h_SumET->Fill(sig_SumET);
            if (hasSumJetET)   sig_h_SumJetET->Fill(sig_SumJetET);
            if (hasSumTowerET) sig_h_SumTowerET->Fill(sig_SumTowerET);
            sig_h_gMET->Fill(clampVal(sig_h_gMET, sig_gMET));
            sig_h_gMET_NC->Fill(clampVal(sig_h_gMET_NC, sig_gMET_NC));
            sig_h_gMET_Rms->Fill(clampVal(sig_h_gMET_Rms, sig_gMET_Rms));
            jFexMETTreeSig->GetEntry(iEvt);
            sig_h_jMET->Fill(clampVal(sig_h_jMET, sig_jMET));
            sig_h_gMETX->Fill(sig_gMETX); sig_h_gMETY->Fill(sig_gMETY);
            sig_h_gMETX_NC->Fill(sig_gMETX_NC); sig_h_gMETY_NC->Fill(sig_gMETY_NC);
            sig_h_gMETX_Rms->Fill(sig_gMETX_Rms); sig_h_gMETY_Rms->Fill(sig_gMETY_Rms);
            if (hasGFexSimMET) {
                gFexMETJwoJSimTreeSig->GetEntry(iEvt);
                gFexMETNoiseCutSimTreeSig->GetEntry(iEvt);
                gFexMETRmsSimTreeSig->GetEntry(iEvt);
                sig_h_gMET_JwoJSim->Fill(clampVal(sig_h_gMET_JwoJSim, sig_gMET_JwoJSim));
                sig_h_gMET_NCSim->Fill(clampVal(sig_h_gMET_NCSim,     sig_gMET_NCSim));
                sig_h_gMET_RmsSim->Fill(clampVal(sig_h_gMET_RmsSim,   sig_gMET_RmsSim));
            }
            sig_h_metTruthNonIntX->Fill(sig_metTruthNonIntX);
            sig_h_metTruthNonIntY->Fill(sig_metTruthNonIntY);
            sig_h_JetMetX->Fill(sig_JetMetX); sig_h_JetMetY->Fill(sig_JetMetY);
            sig_h_TowerMetX->Fill(sig_TowerMetX); sig_h_TowerMetY->Fill(sig_TowerMetY);
            // 2D combined: x=gFEX MET, y=GEP MET (clamped to [0, hi2D])
            auto clamp2D = [&](double v) { return std::min(v, hi2D - 1e-9); };
            sig_h2_JwoJ_Jet->Fill(clamp2D(sig_gMET),     clamp2D(sig_JetMet));
            sig_h2_JwoJ_Tower->Fill(clamp2D(sig_gMET),   clamp2D(sig_TowerMet));
            sig_h2_NC_Jet->Fill(clamp2D(sig_gMET_NC),    clamp2D(sig_JetMet));
            sig_h2_NC_Tower->Fill(clamp2D(sig_gMET_NC),  clamp2D(sig_TowerMet));
            sig_h2_Rms_Jet->Fill(clamp2D(sig_gMET_Rms),  clamp2D(sig_JetMet));
            sig_h2_Rms_Tower->Fill(clamp2D(sig_gMET_Rms),clamp2D(sig_TowerMet));
            sig_h_metTruthNonInt->Fill(clampVal(sig_h_metTruthNonInt, sig_metTruthNonInt));
            sig_h_metTruthNonInt_coarse->Fill(std::min(sig_metTruthNonInt, 599.9));
            sig_h_metTruthInt->Fill(clampVal(sig_h_metTruthInt, sig_metTruthInt));
            sig_h_metTruthIntOut->Fill(clampVal(sig_h_metTruthIntOut, sig_metTruthIntOut));
            sig_h_coreEMTopo_SoftClus_MET->Fill(clampVal(sig_h_coreEMTopo_SoftClus_MET, sig_coreEMTopo_SoftClus_MET));
            sig_h_coreEMTopo_PVSoftTrk_MET->Fill(clampVal(sig_h_coreEMTopo_PVSoftTrk_MET, sig_coreEMTopo_PVSoftTrk_MET));
            sig_h_coreEMTopo_SoftClusEM_MET->Fill(clampVal(sig_h_coreEMTopo_SoftClusEM_MET, sig_coreEMTopo_SoftClusEM_MET));
            sig_h_coreEMPFlow_SoftClus_MET->Fill(clampVal(sig_h_coreEMPFlow_SoftClus_MET, sig_coreEMPFlow_SoftClus_MET));
            sig_h_coreEMPFlow_PVSoftTrk_MET->Fill(clampVal(sig_h_coreEMPFlow_PVSoftTrk_MET, sig_coreEMPFlow_PVSoftTrk_MET));
            h_turnOn_denom->Fill(clampVal(h_turnOn_denom, sig_metTruthNonInt));

            // Calibration / resolution fills (guard against truth = 0 to avoid division by zero)
            if (sig_metTruthNonInt > 0.0) {
                double truthCl = std::min(sig_metTruthNonInt, hi2D - 1e-9);
                auto fillCalib = [&](TH2F* h2corr,
                                     TH1F* h1relRes, TH2F* h2relResVsTruth, TH2F* h2relResVsSumET,
                                     TH1F* h1absRes, TH2F* h2absResVsTruth, TH2F* h2absResVsSumET,
                                     double tobMET, double tobSumET) {
                    h2corr->Fill(truthCl, std::min(tobMET, hi2D - 1e-9));
                    double relRes = (sig_metTruthNonInt - tobMET) / sig_metTruthNonInt;
                    double relResClamp = std::max(-3.0 + 1e-9, std::min(2.0 - 1e-9, relRes));
                    h1relRes->Fill(relResClamp);
                    h2relResVsTruth->Fill(truthCl, relResClamp);
                    h2relResVsSumET->Fill(std::min(tobSumET, hiSumET2D - 1e-9), relResClamp);
                    double absRes = sig_metTruthNonInt - tobMET;
                    double absResClamp = std::max(-hiAbsRes + 1e-9, std::min(hiAbsRes - 1e-9, absRes));
                    h1absRes->Fill(absResClamp);
                    h2absResVsTruth->Fill(truthCl, absResClamp);
                    h2absResVsSumET->Fill(std::min(tobSumET, hiSumET2D - 1e-9), absResClamp);
                };
                fillCalib(sig_h2_gJwoJ_TOBMet_vs_truthMET,
                          sig_h1_gJwoJ_relResidual,    sig_h2_gJwoJ_relResidual_vs_truthMET,    sig_h2_gJwoJ_relResidual_vs_sumET,
                          sig_h1_gJwoJ_absResidual,    sig_h2_gJwoJ_absResidual_vs_truthMET,    sig_h2_gJwoJ_absResidual_vs_sumET,    sig_gMET,     sig_gSumET);
                fillCalib(sig_h2_gNC_TOBMet_vs_truthMET,
                          sig_h1_gNC_relResidual,      sig_h2_gNC_relResidual_vs_truthMET,      sig_h2_gNC_relResidual_vs_sumET,
                          sig_h1_gNC_absResidual,      sig_h2_gNC_absResidual_vs_truthMET,      sig_h2_gNC_absResidual_vs_sumET,      sig_gMET_NC,  sig_gSumET_NC);
                fillCalib(sig_h2_gRms_TOBMet_vs_truthMET,
                          sig_h1_gRms_relResidual,     sig_h2_gRms_relResidual_vs_truthMET,     sig_h2_gRms_relResidual_vs_sumET,
                          sig_h1_gRms_absResidual,     sig_h2_gRms_absResidual_vs_truthMET,     sig_h2_gRms_absResidual_vs_sumET,     sig_gMET_Rms, sig_gSumET_Rms);
                fillCalib(sig_h2_JetMET_TOBMet_vs_truthMET,
                          sig_h1_JetMET_relResidual,   sig_h2_JetMET_relResidual_vs_truthMET,   sig_h2_JetMET_relResidual_vs_sumET,
                          sig_h1_JetMET_absResidual,   sig_h2_JetMET_absResidual_vs_truthMET,   sig_h2_JetMET_absResidual_vs_sumET,   sig_JetMet,   sig_SumET);
                fillCalib(sig_h2_TowerMET_TOBMet_vs_truthMET,
                          sig_h1_TowerMET_relResidual, sig_h2_TowerMET_relResidual_vs_truthMET, sig_h2_TowerMET_relResidual_vs_sumET,
                          sig_h1_TowerMET_absResidual, sig_h2_TowerMET_absResidual_vs_truthMET, sig_h2_TowerMET_absResidual_vs_sumET, sig_TowerMet, sig_SumET);
            }
        }

        // --- Background event loop ---
        unsigned int nBack = metTreeBack->GetEntries();
        std::cout << "  Background events: " << nBack << "\n";
        for (unsigned int iEvt = 0; iEvt < nBack; iEvt++) {
            metTreeBack->GetEntry(iEvt);
            gFexMETTreeBack->GetEntry(iEvt);
            gFexMETNoiseCutTreeBack->GetEntry(iEvt);
            gFexMETRmsTreeBack->GetEntry(iEvt);
            metTruthTreeBack->GetEntry(iEvt);
            coreEMTopoTreeBack->GetEntry(iEvt);
            coreEMPFlowTreeBack->GetEntry(iEvt);
            eventInfoTreeBack->GetEntry(iEvt);
            if (applyHSTPFilter && !passHSTPValuesBack) continue;
            double w = eventWeightsValuesBack->at(0);
            back_h_TotalMET->Fill(clampVal(back_h_TotalMET, back_TotalMET), w);
            back_h_TotalMETX->Fill(back_TotalMETX, w);
            back_h_TotalMETY->Fill(back_TotalMETY, w);
            back_h_TowerMet->Fill(clampVal(back_h_TowerMet, back_TowerMet), w);
            back_h_JetMet->Fill(clampVal(back_h_JetMet, back_JetMet), w);
            back_h_SumET->Fill(back_SumET, w);
            if (hasSumJetET)   back_h_SumJetET->Fill(back_SumJetET, w);
            if (hasSumTowerET) back_h_SumTowerET->Fill(back_SumTowerET, w);
            back_h_gMET->Fill(clampVal(back_h_gMET, back_gMET), w);
            back_h_gMET_NC->Fill(clampVal(back_h_gMET_NC, back_gMET_NC), w);
            back_h_gMET_Rms->Fill(clampVal(back_h_gMET_Rms, back_gMET_Rms), w);
            jFexMETTreeBack->GetEntry(iEvt);
            back_h_jMET->Fill(clampVal(back_h_jMET, back_jMET), w);
            back_h_gMETX->Fill(back_gMETX, w); back_h_gMETY->Fill(back_gMETY, w);
            back_h_gMETX_NC->Fill(back_gMETX_NC, w); back_h_gMETY_NC->Fill(back_gMETY_NC, w);
            back_h_gMETX_Rms->Fill(back_gMETX_Rms, w); back_h_gMETY_Rms->Fill(back_gMETY_Rms, w);
            back_h_metTruthNonIntX->Fill(back_metTruthNonIntX, w);
            back_h_metTruthNonIntY->Fill(back_metTruthNonIntY, w);
            back_h_JetMetX->Fill(back_JetMetX, w); back_h_JetMetY->Fill(back_JetMetY, w);
            back_h_TowerMetX->Fill(back_TowerMetX, w); back_h_TowerMetY->Fill(back_TowerMetY, w);
            back_h_metTruthNonInt->Fill(clampVal(back_h_metTruthNonInt, back_metTruthNonInt), w);
            back_h_metTruthInt->Fill(clampVal(back_h_metTruthInt, back_metTruthInt), w);
            back_h_metTruthIntOut->Fill(clampVal(back_h_metTruthIntOut, back_metTruthIntOut), w);
            back_h_coreEMTopo_SoftClus_MET->Fill(clampVal(back_h_coreEMTopo_SoftClus_MET, back_coreEMTopo_SoftClus_MET), w);
            back_h_coreEMTopo_PVSoftTrk_MET->Fill(clampVal(back_h_coreEMTopo_PVSoftTrk_MET, back_coreEMTopo_PVSoftTrk_MET), w);
            back_h_coreEMTopo_SoftClusEM_MET->Fill(clampVal(back_h_coreEMTopo_SoftClusEM_MET, back_coreEMTopo_SoftClusEM_MET), w);
            back_h_coreEMPFlow_SoftClus_MET->Fill(clampVal(back_h_coreEMPFlow_SoftClus_MET, back_coreEMPFlow_SoftClus_MET), w);
            back_h_coreEMPFlow_PVSoftTrk_MET->Fill(clampVal(back_h_coreEMPFlow_PVSoftTrk_MET, back_coreEMPFlow_PVSoftTrk_MET), w);

            back_hw_TotalMET->Fill(clampVal(back_hw_TotalMET, back_TotalMET), w);
            back_hw_gMET->Fill(clampVal(back_hw_gMET, back_gMET), w);
            back_hw_gMET_NC->Fill(clampVal(back_hw_gMET_NC, back_gMET_NC), w);
            back_hw_gMET_Rms->Fill(clampVal(back_hw_gMET_Rms, back_gMET_Rms), w);
            back_hw_jMET->Fill(clampVal(back_hw_jMET, back_jMET), w);
            back_hw_JetMET->Fill(clampVal(back_hw_JetMET, back_JetMet), w);
            back_hw_TowerMET->Fill(clampVal(back_hw_TowerMET, back_TowerMet), w);
            if (hasGFexSimMET) {
                gFexMETJwoJSimTreeBack->GetEntry(iEvt);
                gFexMETNoiseCutSimTreeBack->GetEntry(iEvt);
                gFexMETRmsSimTreeBack->GetEntry(iEvt);
                back_h_gMET_JwoJSim->Fill(clampVal(back_h_gMET_JwoJSim, back_gMET_JwoJSim), w);
                back_h_gMET_NCSim->Fill(clampVal(back_h_gMET_NCSim,     back_gMET_NCSim),   w);
                back_h_gMET_RmsSim->Fill(clampVal(back_h_gMET_RmsSim,   back_gMET_RmsSim),  w);
                back_hw_gMET_JwoJSim->Fill(clampVal(back_hw_gMET_JwoJSim, back_gMET_JwoJSim), w);
                back_hw_gMET_NCSim->Fill(clampVal(back_hw_gMET_NCSim,     back_gMET_NCSim),   w);
                back_hw_gMET_RmsSim->Fill(clampVal(back_hw_gMET_RmsSim,   back_gMET_RmsSim),  w);
            }
            // 2D combined weighted histograms
            auto clamp2D_b = [&](double v) { return std::min(v, hi2D - 1e-9); };
            back_hw2_JwoJ_Jet->Fill(clamp2D_b(back_gMET),     clamp2D_b(back_JetMet), w);
            back_hw2_JwoJ_Tower->Fill(clamp2D_b(back_gMET),   clamp2D_b(back_TowerMet), w);
            back_hw2_NC_Jet->Fill(clamp2D_b(back_gMET_NC),    clamp2D_b(back_JetMet), w);
            back_hw2_NC_Tower->Fill(clamp2D_b(back_gMET_NC),  clamp2D_b(back_TowerMet), w);
            back_hw2_Rms_Jet->Fill(clamp2D_b(back_gMET_Rms),  clamp2D_b(back_JetMet), w);
            back_hw2_Rms_Tower->Fill(clamp2D_b(back_gMET_Rms),clamp2D_b(back_TowerMet), w);

            if (sampleJZSliceBack >= 0 && sampleJZSliceBack < (int)nJZSlices_) {
                back_h_gMET_jz[sampleJZSliceBack]->Fill(clampVal(back_h_gMET_jz[sampleJZSliceBack],         back_gMET),     w);
                back_h_gMET_NC_jz[sampleJZSliceBack]->Fill(clampVal(back_h_gMET_NC_jz[sampleJZSliceBack],   back_gMET_NC),  w);
                back_h_gMET_Rms_jz[sampleJZSliceBack]->Fill(clampVal(back_h_gMET_Rms_jz[sampleJZSliceBack], back_gMET_Rms), w);
                back_h_jMET_jz[sampleJZSliceBack]->Fill(clampVal(back_h_jMET_jz[sampleJZSliceBack],         back_jMET),     w);
                back_h_JetMET_jz[sampleJZSliceBack]->Fill(clampVal(back_h_JetMET_jz[sampleJZSliceBack],     back_JetMet),   w);
                back_h_TowerMET_jz[sampleJZSliceBack]->Fill(clampVal(back_h_TowerMET_jz[sampleJZSliceBack], back_TowerMet), w);
            }

            // Calibration / resolution fills — background (weighted, guard against truth = 0)
            if (back_metTruthNonInt > 0.0) {
                double truthCl_b = std::min(back_metTruthNonInt, hi2D - 1e-9);
                auto fillCalibBack = [&](TH2F* h2corr,
                                         TH1F* h1relRes, TH2F* h2relResVsTruth, TH2F* h2relResVsSumET,
                                         TH1F* h1absRes, TH2F* h2absResVsTruth, TH2F* h2absResVsSumET,
                                         double tobMET, double tobSumET) {
                    h2corr->Fill(truthCl_b, std::min(tobMET, hi2D - 1e-9), w);
                    double relRes = (back_metTruthNonInt - tobMET) / back_metTruthNonInt;
                    double relResClamp = std::max(-3.0 + 1e-9, std::min(2.0 - 1e-9, relRes));
                    h1relRes->Fill(relResClamp, w);
                    h2relResVsTruth->Fill(truthCl_b, relResClamp, w);
                    h2relResVsSumET->Fill(std::min(tobSumET, hiSumET2D - 1e-9), relResClamp, w);
                    double absRes = back_metTruthNonInt - tobMET;
                    double absResClamp = std::max(-hiAbsRes + 1e-9, std::min(hiAbsRes - 1e-9, absRes));
                    h1absRes->Fill(absResClamp, w);
                    h2absResVsTruth->Fill(truthCl_b, absResClamp, w);
                    h2absResVsSumET->Fill(std::min(tobSumET, hiSumET2D - 1e-9), absResClamp, w);
                };
                fillCalibBack(back_h2_gJwoJ_TOBMet_vs_truthMET,
                              back_h1_gJwoJ_relResidual,    back_h2_gJwoJ_relResidual_vs_truthMET,    back_h2_gJwoJ_relResidual_vs_sumET,
                              back_h1_gJwoJ_absResidual,    back_h2_gJwoJ_absResidual_vs_truthMET,    back_h2_gJwoJ_absResidual_vs_sumET,    back_gMET,     back_gSumET);
                fillCalibBack(back_h2_gNC_TOBMet_vs_truthMET,
                              back_h1_gNC_relResidual,      back_h2_gNC_relResidual_vs_truthMET,      back_h2_gNC_relResidual_vs_sumET,
                              back_h1_gNC_absResidual,      back_h2_gNC_absResidual_vs_truthMET,      back_h2_gNC_absResidual_vs_sumET,      back_gMET_NC,  back_gSumET_NC);
                fillCalibBack(back_h2_gRms_TOBMet_vs_truthMET,
                              back_h1_gRms_relResidual,     back_h2_gRms_relResidual_vs_truthMET,     back_h2_gRms_relResidual_vs_sumET,
                              back_h1_gRms_absResidual,     back_h2_gRms_absResidual_vs_truthMET,     back_h2_gRms_absResidual_vs_sumET,     back_gMET_Rms, back_gSumET_Rms);
                fillCalibBack(back_h2_JetMET_TOBMet_vs_truthMET,
                              back_h1_JetMET_relResidual,   back_h2_JetMET_relResidual_vs_truthMET,   back_h2_JetMET_relResidual_vs_sumET,
                              back_h1_JetMET_absResidual,   back_h2_JetMET_absResidual_vs_truthMET,   back_h2_JetMET_absResidual_vs_sumET,   back_JetMet,   back_SumET);
                fillCalibBack(back_h2_TowerMET_TOBMet_vs_truthMET,
                              back_h1_TowerMET_relResidual, back_h2_TowerMET_relResidual_vs_truthMET, back_h2_TowerMET_relResidual_vs_sumET,
                              back_h1_TowerMET_absResidual, back_h2_TowerMET_absResidual_vs_truthMET, back_h2_TowerMET_absResidual_vs_sumET, back_TowerMet, back_SumET);
            }
        }

        // --- Threshold computation at 20 / 40 / 80 kHz ---
        double thr_gMET_20kHz      = findThreshold(back_hw_gMET,     20e3);
        double thr_gMET_40kHz      = findThreshold(back_hw_gMET,     40e3);
        double thr_gMET_80kHz      = findThreshold(back_hw_gMET,     80e3);
        double thr_gMET_NC_20kHz   = findThreshold(back_hw_gMET_NC,  20e3);
        double thr_gMET_NC_40kHz   = findThreshold(back_hw_gMET_NC,  40e3);
        double thr_gMET_NC_80kHz   = findThreshold(back_hw_gMET_NC,  80e3);
        double thr_gMET_Rms_20kHz  = findThreshold(back_hw_gMET_Rms, 20e3);
        double thr_gMET_Rms_40kHz  = findThreshold(back_hw_gMET_Rms, 40e3);
        double thr_gMET_Rms_80kHz  = findThreshold(back_hw_gMET_Rms, 80e3);
        double thr_jMET_20kHz      = findThreshold(back_hw_jMET,     20e3);
        double thr_jMET_40kHz      = findThreshold(back_hw_jMET,     40e3);
        double thr_jMET_80kHz      = findThreshold(back_hw_jMET,     80e3);
        double thr_JetMET_20kHz    = findThreshold(back_hw_JetMET,   20e3);
        double thr_JetMET_40kHz    = findThreshold(back_hw_JetMET,   40e3);
        double thr_JetMET_80kHz    = findThreshold(back_hw_JetMET,   80e3);
        double thr_TowerMET_20kHz  = findThreshold(back_hw_TowerMET,  20e3);
        double thr_TowerMET_40kHz  = findThreshold(back_hw_TowerMET,  40e3);
        double thr_TowerMET_80kHz  = findThreshold(back_hw_TowerMET,  80e3);
        double thr_TotalMET_20kHz  = findThreshold(back_hw_TotalMET,  20e3);
        double thr_TotalMET_40kHz  = findThreshold(back_hw_TotalMET,  40e3);
        double thr_TotalMET_80kHz  = findThreshold(back_hw_TotalMET,  80e3);
        double thr_gMET_JwoJSim_20kHz = hasGFexSimMET ? findThreshold(back_hw_gMET_JwoJSim, 20e3) : 0.0;
        double thr_gMET_JwoJSim_40kHz = hasGFexSimMET ? findThreshold(back_hw_gMET_JwoJSim, 40e3) : 0.0;
        double thr_gMET_JwoJSim_80kHz = hasGFexSimMET ? findThreshold(back_hw_gMET_JwoJSim, 80e3) : 0.0;
        double thr_gMET_NCSim_20kHz   = hasGFexSimMET ? findThreshold(back_hw_gMET_NCSim,   20e3) : 0.0;
        double thr_gMET_NCSim_40kHz   = hasGFexSimMET ? findThreshold(back_hw_gMET_NCSim,   40e3) : 0.0;
        double thr_gMET_NCSim_80kHz   = hasGFexSimMET ? findThreshold(back_hw_gMET_NCSim,   80e3) : 0.0;
        double thr_gMET_RmsSim_20kHz  = hasGFexSimMET ? findThreshold(back_hw_gMET_RmsSim,  20e3) : 0.0;
        double thr_gMET_RmsSim_40kHz  = hasGFexSimMET ? findThreshold(back_hw_gMET_RmsSim,  40e3) : 0.0;
        double thr_gMET_RmsSim_80kHz  = hasGFexSimMET ? findThreshold(back_hw_gMET_RmsSim,  80e3) : 0.0;
        if (hasGFexSimMET)
            std::cout << "  Thresholds [GeV] — gMET JwoJSim 20/40/80 kHz: "
                      << thr_gMET_JwoJSim_20kHz << " / " << thr_gMET_JwoJSim_40kHz << " / " << thr_gMET_JwoJSim_80kHz << "\n"
                      << "                    gMET NCSim   20/40/80 kHz: "
                      << thr_gMET_NCSim_20kHz   << " / " << thr_gMET_NCSim_40kHz   << " / " << thr_gMET_NCSim_80kHz   << "\n"
                      << "                    gMET RmsSim  20/40/80 kHz: "
                      << thr_gMET_RmsSim_20kHz  << " / " << thr_gMET_RmsSim_40kHz  << " / " << thr_gMET_RmsSim_80kHz  << "\n";
        std::cout << "  Thresholds [GeV] — gMET JwoJ 20/40/80 kHz: "
                  << thr_gMET_20kHz     << " / " << thr_gMET_40kHz     << " / " << thr_gMET_80kHz     << "\n"
                  << "                    gMET NC   20/40/80 kHz: "
                  << thr_gMET_NC_20kHz  << " / " << thr_gMET_NC_40kHz  << " / " << thr_gMET_NC_80kHz  << "\n"
                  << "                    gMET Rms  20/40/80 kHz: "
                  << thr_gMET_Rms_20kHz << " / " << thr_gMET_Rms_40kHz << " / " << thr_gMET_Rms_80kHz << "\n"
                  << "                    jMET      20/40/80 kHz: "
                  << thr_jMET_20kHz     << " / " << thr_jMET_40kHz     << " / " << thr_jMET_80kHz     << "\n"
                  << "                    JetMET    20/40/80 kHz: "
                  << thr_JetMET_20kHz   << " / " << thr_JetMET_40kHz   << " / " << thr_JetMET_80kHz   << "\n"
                  << "                    TowerMET  20/40/80 kHz: "
                  << thr_TowerMET_20kHz << " / " << thr_TowerMET_40kHz << " / " << thr_TowerMET_80kHz << "\n"
                  << "                    TotalMET  20/40/80 kHz: "
                  << thr_TotalMET_20kHz << " / " << thr_TotalMET_40kHz << " / " << thr_TotalMET_80kHz << "\n";

        // --- 2D combined threshold scan ---
        auto out2D_JwoJ_Jet   = MakeRateVsEff_ScanRMin(sig_h2_JwoJ_Jet,   back_hw2_JwoJ_Jet);
        auto out2D_JwoJ_Tower = MakeRateVsEff_ScanRMin(sig_h2_JwoJ_Tower, back_hw2_JwoJ_Tower);
        auto out2D_NC_Jet     = MakeRateVsEff_ScanRMin(sig_h2_NC_Jet,     back_hw2_NC_Jet);
        auto out2D_NC_Tower   = MakeRateVsEff_ScanRMin(sig_h2_NC_Tower,   back_hw2_NC_Tower);
        auto out2D_Rms_Jet    = MakeRateVsEff_ScanRMin(sig_h2_Rms_Jet,    back_hw2_Rms_Jet);
        auto out2D_Rms_Tower  = MakeRateVsEff_ScanRMin(sig_h2_Rms_Tower,  back_hw2_Rms_Tower);

        // Best (t1=gFEX, t2=GEP) at each rate target for each combination
        auto best_JwoJ_Jet_20   = findBestThresholds2D(out2D_JwoJ_Jet,   20e3);
        auto best_JwoJ_Jet_40   = findBestThresholds2D(out2D_JwoJ_Jet,   40e3);
        auto best_JwoJ_Jet_80   = findBestThresholds2D(out2D_JwoJ_Jet,   80e3);
        auto best_JwoJ_Tower_20 = findBestThresholds2D(out2D_JwoJ_Tower, 20e3);
        auto best_JwoJ_Tower_40 = findBestThresholds2D(out2D_JwoJ_Tower, 40e3);
        auto best_JwoJ_Tower_80 = findBestThresholds2D(out2D_JwoJ_Tower, 80e3);
        auto best_NC_Jet_20     = findBestThresholds2D(out2D_NC_Jet,     20e3);
        auto best_NC_Jet_40     = findBestThresholds2D(out2D_NC_Jet,     40e3);
        auto best_NC_Jet_80     = findBestThresholds2D(out2D_NC_Jet,     80e3);
        auto best_NC_Tower_20   = findBestThresholds2D(out2D_NC_Tower,   20e3);
        auto best_NC_Tower_40   = findBestThresholds2D(out2D_NC_Tower,   40e3);
        auto best_NC_Tower_80   = findBestThresholds2D(out2D_NC_Tower,   80e3);
        auto best_Rms_Jet_20    = findBestThresholds2D(out2D_Rms_Jet,    20e3);
        auto best_Rms_Jet_40    = findBestThresholds2D(out2D_Rms_Jet,    40e3);
        auto best_Rms_Jet_80    = findBestThresholds2D(out2D_Rms_Jet,    80e3);
        auto best_Rms_Tower_20  = findBestThresholds2D(out2D_Rms_Tower,  20e3);
        auto best_Rms_Tower_40  = findBestThresholds2D(out2D_Rms_Tower,  40e3);
        auto best_Rms_Tower_80  = findBestThresholds2D(out2D_Rms_Tower,  80e3);
        std::cout << "  Combined best thresholds [gFEX GeV | GEP GeV | eff | rate Hz]\n"
                  << "    JwoJ+Jet  20kHz: " << best_JwoJ_Jet_20.t1 << " | " << best_JwoJ_Jet_20.t2 << " | " << best_JwoJ_Jet_20.eff << " | " << best_JwoJ_Jet_20.rate << "\n"
                  << "    JwoJ+Jet  40kHz: " << best_JwoJ_Jet_40.t1 << " | " << best_JwoJ_Jet_40.t2 << " | " << best_JwoJ_Jet_40.eff << " | " << best_JwoJ_Jet_40.rate << "\n"
                  << "    JwoJ+Jet  80kHz: " << best_JwoJ_Jet_80.t1 << " | " << best_JwoJ_Jet_80.t2 << " | " << best_JwoJ_Jet_80.eff << " | " << best_JwoJ_Jet_80.rate << "\n"
                  << "    JwoJ+Tower 20kHz: " << best_JwoJ_Tower_20.t1 << " | " << best_JwoJ_Tower_20.t2 << " | " << best_JwoJ_Tower_20.eff << " | " << best_JwoJ_Tower_20.rate << "\n"
                  << "    JwoJ+Tower 40kHz: " << best_JwoJ_Tower_40.t1 << " | " << best_JwoJ_Tower_40.t2 << " | " << best_JwoJ_Tower_40.eff << " | " << best_JwoJ_Tower_40.rate << "\n"
                  << "    JwoJ+Tower 80kHz: " << best_JwoJ_Tower_80.t1 << " | " << best_JwoJ_Tower_80.t2 << " | " << best_JwoJ_Tower_80.eff << " | " << best_JwoJ_Tower_80.rate << "\n"
                  << "    NC+Jet  20kHz: " << best_NC_Jet_20.t1 << " | " << best_NC_Jet_20.t2 << " | " << best_NC_Jet_20.eff << " | " << best_NC_Jet_20.rate << "\n"
                  << "    NC+Jet  40kHz: " << best_NC_Jet_40.t1 << " | " << best_NC_Jet_40.t2 << " | " << best_NC_Jet_40.eff << " | " << best_NC_Jet_40.rate << "\n"
                  << "    NC+Jet  80kHz: " << best_NC_Jet_80.t1 << " | " << best_NC_Jet_80.t2 << " | " << best_NC_Jet_80.eff << " | " << best_NC_Jet_80.rate << "\n"
                  << "    NC+Tower 20kHz: " << best_NC_Tower_20.t1 << " | " << best_NC_Tower_20.t2 << " | " << best_NC_Tower_20.eff << " | " << best_NC_Tower_20.rate << "\n"
                  << "    NC+Tower 40kHz: " << best_NC_Tower_40.t1 << " | " << best_NC_Tower_40.t2 << " | " << best_NC_Tower_40.eff << " | " << best_NC_Tower_40.rate << "\n"
                  << "    NC+Tower 80kHz: " << best_NC_Tower_80.t1 << " | " << best_NC_Tower_80.t2 << " | " << best_NC_Tower_80.eff << " | " << best_NC_Tower_80.rate << "\n"
                  << "    Rms+Jet  20kHz: " << best_Rms_Jet_20.t1 << " | " << best_Rms_Jet_20.t2 << " | " << best_Rms_Jet_20.eff << " | " << best_Rms_Jet_20.rate << "\n"
                  << "    Rms+Jet  40kHz: " << best_Rms_Jet_40.t1 << " | " << best_Rms_Jet_40.t2 << " | " << best_Rms_Jet_40.eff << " | " << best_Rms_Jet_40.rate << "\n"
                  << "    Rms+Jet  80kHz: " << best_Rms_Jet_80.t1 << " | " << best_Rms_Jet_80.t2 << " | " << best_Rms_Jet_80.eff << " | " << best_Rms_Jet_80.rate << "\n"
                  << "    Rms+Tower 20kHz: " << best_Rms_Tower_20.t1 << " | " << best_Rms_Tower_20.t2 << " | " << best_Rms_Tower_20.eff << " | " << best_Rms_Tower_20.rate << "\n"
                  << "    Rms+Tower 40kHz: " << best_Rms_Tower_40.t1 << " | " << best_Rms_Tower_40.t2 << " | " << best_Rms_Tower_40.eff << " | " << best_Rms_Tower_40.rate << "\n"
                  << "    Rms+Tower 80kHz: " << best_Rms_Tower_80.t1 << " | " << best_Rms_Tower_80.t2 << " | " << best_Rms_Tower_80.eff << " | " << best_Rms_Tower_80.rate << "\n";

        // --- Second signal pass: fill turn-on numerators ---
        for (unsigned int iEvt = 0; iEvt < nSig; iEvt++) {
            metTreeSig->GetEntry(iEvt);
            gFexMETTreeSig->GetEntry(iEvt);
            gFexMETNoiseCutTreeSig->GetEntry(iEvt);
            gFexMETRmsTreeSig->GetEntry(iEvt);
            metTruthTreeSig->GetEntry(iEvt);
            double truthMET = clampVal(h_turnOn_denom, sig_metTruthNonInt);
            // Individual gFEX and GEP turn-ons
            if (sig_gMET     > thr_gMET_20kHz)     h_turnOn_num_gMET_20kHz->Fill(truthMET);
            if (sig_gMET     > thr_gMET_40kHz)     h_turnOn_num_gMET_40kHz->Fill(truthMET);
            if (sig_gMET     > thr_gMET_80kHz)     h_turnOn_num_gMET_80kHz->Fill(truthMET);
            if (sig_gMET_NC  > thr_gMET_NC_20kHz)  h_turnOn_num_gMET_NC_20kHz->Fill(truthMET);
            if (sig_gMET_NC  > thr_gMET_NC_40kHz)  h_turnOn_num_gMET_NC_40kHz->Fill(truthMET);
            if (sig_gMET_NC  > thr_gMET_NC_80kHz)  h_turnOn_num_gMET_NC_80kHz->Fill(truthMET);
            if (sig_gMET_Rms > thr_gMET_Rms_20kHz) h_turnOn_num_gMET_Rms_20kHz->Fill(truthMET);
            if (sig_gMET_Rms > thr_gMET_Rms_40kHz) h_turnOn_num_gMET_Rms_40kHz->Fill(truthMET);
            if (sig_gMET_Rms > thr_gMET_Rms_80kHz) h_turnOn_num_gMET_Rms_80kHz->Fill(truthMET);
            jFexMETTreeSig->GetEntry(iEvt);
            if (sig_jMET     > thr_jMET_20kHz)     h_turnOn_num_jMET_20kHz->Fill(truthMET);
            if (sig_jMET     > thr_jMET_40kHz)     h_turnOn_num_jMET_40kHz->Fill(truthMET);
            if (sig_jMET     > thr_jMET_80kHz)     h_turnOn_num_jMET_80kHz->Fill(truthMET);
            if (sig_JetMet   > thr_JetMET_20kHz)   h_turnOn_num_JetMET_20kHz->Fill(truthMET);
            if (sig_JetMet   > thr_JetMET_40kHz)   h_turnOn_num_JetMET_40kHz->Fill(truthMET);
            if (sig_JetMet   > thr_JetMET_80kHz)   h_turnOn_num_JetMET_80kHz->Fill(truthMET);
            if (sig_TowerMet  > thr_TowerMET_20kHz)  h_turnOn_num_TowerMET_20kHz->Fill(truthMET);
            if (sig_TowerMet  > thr_TowerMET_40kHz)  h_turnOn_num_TowerMET_40kHz->Fill(truthMET);
            if (sig_TowerMet  > thr_TowerMET_80kHz)  h_turnOn_num_TowerMET_80kHz->Fill(truthMET);
            if (sig_TotalMET  > thr_TotalMET_20kHz)  h_turnOn_num_TotalMET_20kHz->Fill(truthMET);
            if (sig_TotalMET  > thr_TotalMET_40kHz)  h_turnOn_num_TotalMET_40kHz->Fill(truthMET);
            if (sig_TotalMET  > thr_TotalMET_80kHz)  h_turnOn_num_TotalMET_80kHz->Fill(truthMET);
            if (hasGFexSimMET) {
                gFexMETJwoJSimTreeSig->GetEntry(iEvt);
                gFexMETNoiseCutSimTreeSig->GetEntry(iEvt);
                gFexMETRmsSimTreeSig->GetEntry(iEvt);
                if (sig_gMET_JwoJSim > thr_gMET_JwoJSim_20kHz) h_turnOn_num_gMET_JwoJSim_20kHz->Fill(truthMET);
                if (sig_gMET_JwoJSim > thr_gMET_JwoJSim_40kHz) h_turnOn_num_gMET_JwoJSim_40kHz->Fill(truthMET);
                if (sig_gMET_JwoJSim > thr_gMET_JwoJSim_80kHz) h_turnOn_num_gMET_JwoJSim_80kHz->Fill(truthMET);
                if (sig_gMET_NCSim   > thr_gMET_NCSim_20kHz)   h_turnOn_num_gMET_NCSim_20kHz->Fill(truthMET);
                if (sig_gMET_NCSim   > thr_gMET_NCSim_40kHz)   h_turnOn_num_gMET_NCSim_40kHz->Fill(truthMET);
                if (sig_gMET_NCSim   > thr_gMET_NCSim_80kHz)   h_turnOn_num_gMET_NCSim_80kHz->Fill(truthMET);
                if (sig_gMET_RmsSim  > thr_gMET_RmsSim_20kHz)  h_turnOn_num_gMET_RmsSim_20kHz->Fill(truthMET);
                if (sig_gMET_RmsSim  > thr_gMET_RmsSim_40kHz)  h_turnOn_num_gMET_RmsSim_40kHz->Fill(truthMET);
                if (sig_gMET_RmsSim  > thr_gMET_RmsSim_80kHz)  h_turnOn_num_gMET_RmsSim_80kHz->Fill(truthMET);
            }
            // Combined gFEX + GEP turn-ons at best thresholds
            if (sig_gMET     > best_JwoJ_Jet_20.t1   && sig_JetMet   > best_JwoJ_Jet_20.t2)   h_turnOn_num_combo_JwoJ_Jet_20kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Jet_40.t1   && sig_JetMet   > best_JwoJ_Jet_40.t2)   h_turnOn_num_combo_JwoJ_Jet_40kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Jet_80.t1   && sig_JetMet   > best_JwoJ_Jet_80.t2)   h_turnOn_num_combo_JwoJ_Jet_80kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Tower_20.t1 && sig_TowerMet > best_JwoJ_Tower_20.t2) h_turnOn_num_combo_JwoJ_Tower_20kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Tower_40.t1 && sig_TowerMet > best_JwoJ_Tower_40.t2) h_turnOn_num_combo_JwoJ_Tower_40kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Tower_80.t1 && sig_TowerMet > best_JwoJ_Tower_80.t2) h_turnOn_num_combo_JwoJ_Tower_80kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Jet_20.t1     && sig_JetMet   > best_NC_Jet_20.t2)     h_turnOn_num_combo_NC_Jet_20kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Jet_40.t1     && sig_JetMet   > best_NC_Jet_40.t2)     h_turnOn_num_combo_NC_Jet_40kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Jet_80.t1     && sig_JetMet   > best_NC_Jet_80.t2)     h_turnOn_num_combo_NC_Jet_80kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Tower_20.t1   && sig_TowerMet > best_NC_Tower_20.t2)   h_turnOn_num_combo_NC_Tower_20kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Tower_40.t1   && sig_TowerMet > best_NC_Tower_40.t2)   h_turnOn_num_combo_NC_Tower_40kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Tower_80.t1   && sig_TowerMet > best_NC_Tower_80.t2)   h_turnOn_num_combo_NC_Tower_80kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Jet_20.t1    && sig_JetMet   > best_Rms_Jet_20.t2)    h_turnOn_num_combo_Rms_Jet_20kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Jet_40.t1    && sig_JetMet   > best_Rms_Jet_40.t2)    h_turnOn_num_combo_Rms_Jet_40kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Jet_80.t1    && sig_JetMet   > best_Rms_Jet_80.t2)    h_turnOn_num_combo_Rms_Jet_80kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Tower_20.t1  && sig_TowerMet > best_Rms_Tower_20.t2)  h_turnOn_num_combo_Rms_Tower_20kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Tower_40.t1  && sig_TowerMet > best_Rms_Tower_40.t2)  h_turnOn_num_combo_Rms_Tower_40kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Tower_80.t1  && sig_TowerMet > best_Rms_Tower_80.t2)  h_turnOn_num_combo_Rms_Tower_80kHz->Fill(truthMET);

            // --- Event properties at 80 kHz selections ---
            if (hasWTAConeJets) gepWTAConeCellsTowersJetsTreeSig->GetEntry(iEvt);
            if (hasTruthAntiKt4WZDressed) truthAntiKt4WZDressedJetsTreeSig->GetEntry(iEvt);
            if (hasInTimeAntiKt4TruthJets) inTimeAntiKt4TruthJetsTreeSig->GetEntry(iEvt);
            {
                double phi_GEP  = std::atan2(sig_TotalMETY, sig_TotalMETX);
                double phi_gFEX = std::atan2(sig_gMETY, sig_gMETX);
                double METsig   = (sig_SumET  > 0) ? sig_TotalMET / std::sqrt(sig_SumET)  : 0.0;
                double gMETsig  = (sig_gSumET > 0) ? sig_gMET    / std::sqrt(sig_gSumET) : 0.0;
                int    nJets    = 0;
                int    nWTAConeJets = 0;
                double jet1pt = 0, jet1eta = 0, jet1phi = 0, jet2pt = 0;
                if (hasWTAConeJets && gepWTAConeCellsTowersJetsEtValuesSig &&
                    !gepWTAConeCellsTowersJetsEtValuesSig->empty()) {
                    nJets   = (int)gepWTAConeCellsTowersJetsEtValuesSig->size();
                    jet1pt  = gepWTAConeCellsTowersJetsEtValuesSig->at(0);
                    jet1eta = gepWTAConeCellsTowersJetsEtaValuesSig->at(0);
                    jet1phi = gepWTAConeCellsTowersJetsPhiValuesSig->at(0);
                    if (nJets >= 2) jet2pt = gepWTAConeCellsTowersJetsEtValuesSig->at(1);
                    for (int iJ = 0; iJ < nJets; iJ++)
                        if (gepWTAConeCellsTowersJetsEtValuesSig->at(iJ) > 25.0) nWTAConeJets++;
                }
                int nTruthJets = 0;
                if (hasTruthAntiKt4WZDressed && truthAntiKt4WZDressedJetsEtValuesSig)
                    for (double et : *truthAntiKt4WZDressedJetsEtValuesSig)
                        if (et > 15.0) nTruthJets++;
                int nPileupJets = 0;
                if (hasInTimeAntiKt4TruthJets && inTimeAntiKt4TruthJetsEtValuesSig)
                    for (double et : *inTimeAntiKt4TruthJetsEtValuesSig)
                        if (et > 15.0) nPileupJets++;
                // Per-selection TOB MET phi: Incl/JwoJ=gFEX JwoJ, NC, Rms, Jet, Tower
                double phi_TOBsig[nSel80] = {
                    phi_gFEX,
                    phi_gFEX,
                    std::atan2(sig_gMETY_NC,  sig_gMETX_NC),
                    std::atan2(sig_gMETY_Rms, sig_gMETX_Rms),
                    hasJetMetXY   ? std::atan2(sig_JetMetY,   sig_JetMetX)   : 0.0,
                    hasTowerMetXY ? std::atan2(sig_TowerMetY, sig_TowerMetX) : 0.0
                };
                bool passes80[nSel80] = {
                    true,
                    sig_gMET     > thr_gMET_80kHz,
                    sig_gMET_NC  > thr_gMET_NC_80kHz,
                    sig_gMET_Rms > thr_gMET_Rms_80kHz,
                    sig_JetMet   > thr_JetMET_80kHz,
                    sig_TowerMet > thr_TowerMET_80kHz
                };
                for (int iSel = 0; iSel < nSel80; iSel++) {
                    if (!passes80[iSel]) continue;
                    sig_sel_truthMET[iSel]->Fill(clampVal(sig_sel_truthMET[iSel], sig_metTruthNonInt));
                    sig_sel_SumET[iSel]->Fill(std::min(sig_SumET, 1039.9));
                    sig_sel_gSumET[iSel]->Fill(std::min(sig_gSumET, 1039.9));
                    if (METsig  < 25.0) sig_sel_METsig[iSel]->Fill(METsig);
                    if (gMETsig < 25.0) sig_sel_gMETsig[iSel]->Fill(gMETsig);
                    sig_sel_dPhi_GEP_gFEX[iSel]->Fill(absDeltaPhi(phi_GEP, phi_gFEX));
                    if (hasTruthNonIntXY && (sig_metTruthNonIntX != 0.0 || sig_metTruthNonIntY != 0.0)) {
                        double phi_truth = std::atan2(sig_metTruthNonIntY, sig_metTruthNonIntX);
                        sig_sel_dPhi_GEP_truth[iSel]->Fill(absDeltaPhi(phi_GEP, phi_truth));
                        bool hasTOBphi = (iSel < 4) || (iSel == 4 && hasJetMetXY) || (iSel == 5 && hasTowerMetXY);
                        if (hasTOBphi) {
                            sig_sel_dPhi_truth_TOB[iSel]->Fill(absDeltaPhi(phi_truth, phi_TOBsig[iSel]));
                            sig_sel_phi2D_TOB_truth[iSel]->Fill(phi_truth, phi_TOBsig[iSel]);
                        }
                    }
                    if (hasWTAConeJets) {
                        sig_sel_nJets[iSel]->Fill(nWTAConeJets);
                        if (nJets >= 1) {
                            sig_sel_jet1pt[iSel]->Fill(jet1pt);
                            sig_sel_jet1eta[iSel]->Fill(jet1eta);
                            sig_sel_dPhi_jet1_MET[iSel]->Fill(absDeltaPhi(jet1phi, phi_GEP));
                        }
                        if (nJets >= 2) sig_sel_jet2pt[iSel]->Fill(jet2pt);
                    }
                    if (hasTruthAntiKt4WZDressed)
                        sig_sel_nTruthJets[iSel]->Fill(nTruthJets);
                    if (hasInTimeAntiKt4TruthJets)
                        sig_sel_nPileupJets[iSel]->Fill(nPileupJets);
                }
            }
        }

        // --- Background event-property pass: fill selection histograms at 80 kHz ---
        for (unsigned int iEvt = 0; iEvt < nBack; iEvt++) {
            metTreeBack->GetEntry(iEvt);
            gFexMETTreeBack->GetEntry(iEvt);
            gFexMETNoiseCutTreeBack->GetEntry(iEvt);
            gFexMETRmsTreeBack->GetEntry(iEvt);
            metTruthTreeBack->GetEntry(iEvt);
            eventInfoTreeBack->GetEntry(iEvt);
            if (applyHSTPFilter && !passHSTPValuesBack) continue;
            if (hasWTAConeJets) gepWTAConeCellsTowersJetsTreeBack->GetEntry(iEvt);
            if (hasTruthAntiKt4WZDressed) truthAntiKt4WZDressedJetsTreeBack->GetEntry(iEvt);
            if (hasInTimeAntiKt4TruthJets) inTimeAntiKt4TruthJetsTreeBack->GetEntry(iEvt);
            double w = eventWeightsValuesBack->at(0);

            double phi_GEP  = std::atan2(back_TotalMETY, back_TotalMETX);
            double phi_gFEX = std::atan2(back_gMETY, back_gMETX);
            double METsig   = (back_SumET  > 0) ? back_TotalMET / std::sqrt(back_SumET)  : 0.0;
            double gMETsig  = (back_gSumET > 0) ? back_gMET    / std::sqrt(back_gSumET) : 0.0;
            int    nJets    = 0;
            int    nWTAConeJets = 0;
            double jet1pt = 0, jet1eta = 0, jet1phi = 0, jet2pt = 0;
            if (hasWTAConeJets && gepWTAConeCellsTowersJetsEtValuesBack &&
                !gepWTAConeCellsTowersJetsEtValuesBack->empty()) {
                nJets   = (int)gepWTAConeCellsTowersJetsEtValuesBack->size();
                jet1pt  = gepWTAConeCellsTowersJetsEtValuesBack->at(0);
                jet1eta = gepWTAConeCellsTowersJetsEtaValuesBack->at(0);
                jet1phi = gepWTAConeCellsTowersJetsPhiValuesBack->at(0);
                if (nJets >= 2) jet2pt = gepWTAConeCellsTowersJetsEtValuesBack->at(1);
                for (int iJ = 0; iJ < nJets; iJ++)
                    if (gepWTAConeCellsTowersJetsEtValuesBack->at(iJ) > 25.0) nWTAConeJets++;
            }
            int nTruthJets = 0;
            if (hasTruthAntiKt4WZDressed && truthAntiKt4WZDressedJetsEtValuesBack)
                for (double et : *truthAntiKt4WZDressedJetsEtValuesBack)
                    if (et > 15.0) nTruthJets++;
            int nPileupJets = 0;
            if (hasInTimeAntiKt4TruthJets && inTimeAntiKt4TruthJetsEtValuesBack)
                for (double et : *inTimeAntiKt4TruthJetsEtValuesBack)
                    if (et > 15.0) nPileupJets++;
            // Per-selection TOB MET phi: Incl/JwoJ=gFEX JwoJ, NC, Rms, Jet, Tower
            double phi_TOBback[nSel80] = {
                phi_gFEX,
                phi_gFEX,
                std::atan2(back_gMETY_NC,  back_gMETX_NC),
                std::atan2(back_gMETY_Rms, back_gMETX_Rms),
                hasJetMetXY   ? std::atan2(back_JetMetY,   back_JetMetX)   : 0.0,
                hasTowerMetXY ? std::atan2(back_TowerMetY, back_TowerMetX) : 0.0
            };
            bool passesB[nSel80] = {
                true,
                back_gMET     > thr_gMET_80kHz,
                back_gMET_NC  > thr_gMET_NC_80kHz,
                back_gMET_Rms > thr_gMET_Rms_80kHz,
                back_JetMet   > thr_JetMET_80kHz,
                back_TowerMet > thr_TowerMET_80kHz
            };
            for (int iSel = 0; iSel < nSel80; iSel++) {
                if (!passesB[iSel]) continue;
                back_sel_truthMET[iSel]->Fill(clampVal(back_sel_truthMET[iSel], back_metTruthNonInt), w);
                back_sel_SumET[iSel]->Fill(std::min(back_SumET, 1039.9), w);
                back_sel_gSumET[iSel]->Fill(std::min(back_gSumET, 1039.9), w);
                if (METsig  < 25.0) back_sel_METsig[iSel]->Fill(METsig,  w);
                if (gMETsig < 25.0) back_sel_gMETsig[iSel]->Fill(gMETsig, w);
                back_sel_dPhi_GEP_gFEX[iSel]->Fill(absDeltaPhi(phi_GEP, phi_gFEX), w);
                if (hasTruthNonIntXY && (back_metTruthNonIntX != 0.0 || back_metTruthNonIntY != 0.0)) {
                    double phi_truth = std::atan2(back_metTruthNonIntY, back_metTruthNonIntX);
                    back_sel_dPhi_GEP_truth[iSel]->Fill(absDeltaPhi(phi_GEP, phi_truth), w);
                    bool hasTOBphi = (iSel < 4) || (iSel == 4 && hasJetMetXY) || (iSel == 5 && hasTowerMetXY);
                    if (hasTOBphi) {
                        back_sel_dPhi_truth_TOB[iSel]->Fill(absDeltaPhi(phi_truth, phi_TOBback[iSel]), w);
                        back_sel_phi2D_TOB_truth[iSel]->Fill(phi_truth, phi_TOBback[iSel], w);
                    }
                }
                if (hasWTAConeJets) {
                    back_sel_nJets[iSel]->Fill(nWTAConeJets, w);
                    if (nJets >= 1) {
                        back_sel_jet1pt[iSel]->Fill(jet1pt, w);
                        back_sel_jet1eta[iSel]->Fill(jet1eta, w);
                        back_sel_dPhi_jet1_MET[iSel]->Fill(absDeltaPhi(jet1phi, phi_GEP), w);
                    }
                    if (nJets >= 2) back_sel_jet2pt[iSel]->Fill(jet2pt, w);
                }
                if (hasTruthAntiKt4WZDressed)
                    back_sel_nTruthJets[iSel]->Fill(nTruthJets, w);
                if (hasInTimeAntiKt4TruthJets)
                    back_sel_nPileupJets[iSel]->Fill(nPileupJets, w);
            }
        }

        // --- Compute efficiencies via binomial division ---
        auto makeEff = [&](TH1F* num, const std::string& name) -> TH1F* {
            TH1F* eff = (TH1F*)num->Clone(name.c_str());
            eff->SetDirectory(0);
            eff->Divide(num, h_turnOn_denom, 1.0, 1.0, "B");
            return eff;
        };
        TH1F* eff_gMET_20kHz      = makeEff(h_turnOn_num_gMET_20kHz,      "eff_gMET_20kHz_"     +tag);
        TH1F* eff_gMET_40kHz      = makeEff(h_turnOn_num_gMET_40kHz,      "eff_gMET_40kHz_"     +tag);
        TH1F* eff_gMET_80kHz      = makeEff(h_turnOn_num_gMET_80kHz,      "eff_gMET_80kHz_"     +tag);
        TH1F* eff_gMET_NC_20kHz   = makeEff(h_turnOn_num_gMET_NC_20kHz,   "eff_gMET_NC_20kHz_"  +tag);
        TH1F* eff_gMET_NC_40kHz   = makeEff(h_turnOn_num_gMET_NC_40kHz,   "eff_gMET_NC_40kHz_"  +tag);
        TH1F* eff_gMET_NC_80kHz   = makeEff(h_turnOn_num_gMET_NC_80kHz,   "eff_gMET_NC_80kHz_"  +tag);
        TH1F* eff_gMET_Rms_20kHz  = makeEff(h_turnOn_num_gMET_Rms_20kHz,  "eff_gMET_Rms_20kHz_" +tag);
        TH1F* eff_gMET_Rms_40kHz  = makeEff(h_turnOn_num_gMET_Rms_40kHz,  "eff_gMET_Rms_40kHz_" +tag);
        TH1F* eff_gMET_Rms_80kHz  = makeEff(h_turnOn_num_gMET_Rms_80kHz,  "eff_gMET_Rms_80kHz_" +tag);
        TH1F* eff_jMET_20kHz      = makeEff(h_turnOn_num_jMET_20kHz,      "eff_jMET_20kHz_"     +tag);
        TH1F* eff_jMET_40kHz      = makeEff(h_turnOn_num_jMET_40kHz,      "eff_jMET_40kHz_"     +tag);
        TH1F* eff_jMET_80kHz      = makeEff(h_turnOn_num_jMET_80kHz,      "eff_jMET_80kHz_"     +tag);
        TH1F* eff_JetMET_20kHz    = makeEff(h_turnOn_num_JetMET_20kHz,    "eff_JetMET_20kHz_"   +tag);
        TH1F* eff_JetMET_40kHz    = makeEff(h_turnOn_num_JetMET_40kHz,    "eff_JetMET_40kHz_"   +tag);
        TH1F* eff_JetMET_80kHz    = makeEff(h_turnOn_num_JetMET_80kHz,    "eff_JetMET_80kHz_"   +tag);
        TH1F* eff_TowerMET_20kHz  = makeEff(h_turnOn_num_TowerMET_20kHz,  "eff_TowerMET_20kHz_" +tag);
        TH1F* eff_TowerMET_40kHz  = makeEff(h_turnOn_num_TowerMET_40kHz,  "eff_TowerMET_40kHz_" +tag);
        TH1F* eff_TowerMET_80kHz  = makeEff(h_turnOn_num_TowerMET_80kHz,  "eff_TowerMET_80kHz_" +tag);
        TH1F* eff_TotalMET_20kHz  = makeEff(h_turnOn_num_TotalMET_20kHz,  "eff_TotalMET_20kHz_" +tag);
        TH1F* eff_TotalMET_40kHz  = makeEff(h_turnOn_num_TotalMET_40kHz,  "eff_TotalMET_40kHz_" +tag);
        TH1F* eff_TotalMET_80kHz  = makeEff(h_turnOn_num_TotalMET_80kHz,  "eff_TotalMET_80kHz_" +tag);
        // Combined efficiencies
        TH1F* eff_combo_JwoJ_Jet_20kHz   = makeEff(h_turnOn_num_combo_JwoJ_Jet_20kHz,   "eff_combo_JwoJ_Jet_20kHz_"  +tag);
        TH1F* eff_combo_JwoJ_Jet_40kHz   = makeEff(h_turnOn_num_combo_JwoJ_Jet_40kHz,   "eff_combo_JwoJ_Jet_40kHz_"  +tag);
        TH1F* eff_combo_JwoJ_Jet_80kHz   = makeEff(h_turnOn_num_combo_JwoJ_Jet_80kHz,   "eff_combo_JwoJ_Jet_80kHz_"  +tag);
        TH1F* eff_combo_JwoJ_Tower_20kHz = makeEff(h_turnOn_num_combo_JwoJ_Tower_20kHz, "eff_combo_JwoJ_Tower_20kHz_"+tag);
        TH1F* eff_combo_JwoJ_Tower_40kHz = makeEff(h_turnOn_num_combo_JwoJ_Tower_40kHz, "eff_combo_JwoJ_Tower_40kHz_"+tag);
        TH1F* eff_combo_JwoJ_Tower_80kHz = makeEff(h_turnOn_num_combo_JwoJ_Tower_80kHz, "eff_combo_JwoJ_Tower_80kHz_"+tag);
        TH1F* eff_combo_NC_Jet_20kHz     = makeEff(h_turnOn_num_combo_NC_Jet_20kHz,     "eff_combo_NC_Jet_20kHz_"    +tag);
        TH1F* eff_combo_NC_Jet_40kHz     = makeEff(h_turnOn_num_combo_NC_Jet_40kHz,     "eff_combo_NC_Jet_40kHz_"    +tag);
        TH1F* eff_combo_NC_Jet_80kHz     = makeEff(h_turnOn_num_combo_NC_Jet_80kHz,     "eff_combo_NC_Jet_80kHz_"    +tag);
        TH1F* eff_combo_NC_Tower_20kHz   = makeEff(h_turnOn_num_combo_NC_Tower_20kHz,   "eff_combo_NC_Tower_20kHz_"  +tag);
        TH1F* eff_combo_NC_Tower_40kHz   = makeEff(h_turnOn_num_combo_NC_Tower_40kHz,   "eff_combo_NC_Tower_40kHz_"  +tag);
        TH1F* eff_combo_NC_Tower_80kHz   = makeEff(h_turnOn_num_combo_NC_Tower_80kHz,   "eff_combo_NC_Tower_80kHz_"  +tag);
        TH1F* eff_combo_Rms_Jet_20kHz    = makeEff(h_turnOn_num_combo_Rms_Jet_20kHz,    "eff_combo_Rms_Jet_20kHz_"   +tag);
        TH1F* eff_combo_Rms_Jet_40kHz    = makeEff(h_turnOn_num_combo_Rms_Jet_40kHz,    "eff_combo_Rms_Jet_40kHz_"   +tag);
        TH1F* eff_combo_Rms_Jet_80kHz    = makeEff(h_turnOn_num_combo_Rms_Jet_80kHz,    "eff_combo_Rms_Jet_80kHz_"   +tag);
        TH1F* eff_combo_Rms_Tower_20kHz  = makeEff(h_turnOn_num_combo_Rms_Tower_20kHz,  "eff_combo_Rms_Tower_20kHz_" +tag);
        TH1F* eff_combo_Rms_Tower_40kHz  = makeEff(h_turnOn_num_combo_Rms_Tower_40kHz,  "eff_combo_Rms_Tower_40kHz_" +tag);
        TH1F* eff_combo_Rms_Tower_80kHz  = makeEff(h_turnOn_num_combo_Rms_Tower_80kHz,  "eff_combo_Rms_Tower_80kHz_" +tag);
        TH1F* eff_gMET_JwoJSim_20kHz = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_JwoJSim_20kHz, "eff_gMET_JwoJSim_20kHz_"+tag) : nullptr;
        TH1F* eff_gMET_JwoJSim_40kHz = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_JwoJSim_40kHz, "eff_gMET_JwoJSim_40kHz_"+tag) : nullptr;
        TH1F* eff_gMET_JwoJSim_80kHz = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_JwoJSim_80kHz, "eff_gMET_JwoJSim_80kHz_"+tag) : nullptr;
        TH1F* eff_gMET_NCSim_20kHz   = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_NCSim_20kHz,   "eff_gMET_NCSim_20kHz_"+tag)   : nullptr;
        TH1F* eff_gMET_NCSim_40kHz   = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_NCSim_40kHz,   "eff_gMET_NCSim_40kHz_"+tag)   : nullptr;
        TH1F* eff_gMET_NCSim_80kHz   = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_NCSim_80kHz,   "eff_gMET_NCSim_80kHz_"+tag)   : nullptr;
        TH1F* eff_gMET_RmsSim_20kHz  = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_RmsSim_20kHz,  "eff_gMET_RmsSim_20kHz_"+tag)  : nullptr;
        TH1F* eff_gMET_RmsSim_40kHz  = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_RmsSim_40kHz,  "eff_gMET_RmsSim_40kHz_"+tag)  : nullptr;
        TH1F* eff_gMET_RmsSim_80kHz  = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_RmsSim_80kHz,  "eff_gMET_RmsSim_80kHz_"+tag)  : nullptr;

        // --- Per-file signal vs background overlays ---
        // Extract signal tag from the emulation-output basename (includes config tags: N_Towers, jetEt, SK/OR)
        std::string sigPath = signalFiles[fileIt].second;
        std::string sigBase = sigPath.substr(sigPath.rfind('/') + 1);
        std::string sigTag  = sigBase.substr(0, sigBase.rfind('.'));  // strip .root
        std::string fDir = outputDir + "metPlots/" + sigTag + "_" + labels[fileIt] + "/";
        gSystem->mkdir(fDir.c_str(), true);

        drawOverlay(sig_h_TotalMET,        back_h_TotalMET,        "Total MET (GEP)",      "Total MET (GEP) [GeV]",          fDir + "TotalMET.pdf");
        drawComponentOverlay(sig_h_TotalMETX, back_h_TotalMETX, "Total MET_{x} (GEP)", "Total MET_{x} (GEP) [GeV]",      fDir + "TotalMETx.pdf");
        drawComponentOverlay(sig_h_TotalMETY, back_h_TotalMETY, "Total MET_{y} (GEP)", "Total MET_{y} (GEP) [GeV]",      fDir + "TotalMETy.pdf");
        drawOverlay(sig_h_TowerMet,        back_h_TowerMet,        "Tower MET (GEP)",       "Tower MET (GEP) [GeV]",          fDir + "TowerMET.pdf");
        drawOverlay(sig_h_JetMet,          back_h_JetMet,          "Jet MET (GEP)",         "Jet MET (GEP) [GeV]",          fDir + "JetMET.pdf");
        drawOverlay(sig_h_SumET,           back_h_SumET,           "GEP TOB #Sigma E_{T}",         "GEP TOB #Sigma E_{T} [GeV]",          fDir + "SumET.pdf");
        if (hasSumJetET)   drawOverlay(sig_h_SumJetET,   back_h_SumJetET,   "GEP H_{T} (Sum Jet E_{T})",    "GEP H_{T} [GeV]",                     fDir + "SumJetET.pdf");
        if (hasSumTowerET) drawOverlay(sig_h_SumTowerET, back_h_SumTowerET, "GEP Tower #Sigma E_{T}",       "GEP Tower #Sigma E_{T} [GeV]",        fDir + "SumTowerET.pdf");
        drawOverlay(sig_h_gMET,            back_h_gMET,            "gFEX MET (JwoJ)",             "MET (JwoJ) [GeV]",          fDir + "gFEX_MET_JwoJ.pdf");
        drawOverlay(sig_h_gMET_NC,         back_h_gMET_NC,         "gFEX MET (NoiseCut)",         "MET (NoiseCut) [GeV]",          fDir + "gFEX_MET_NoiseCut.pdf");
        drawOverlay(sig_h_gMET_Rms,        back_h_gMET_Rms,        "gFEX MET (Rms)",              "MET (Rms) [GeV]",          fDir + "gFEX_MET_Rms.pdf");
        drawOverlay(sig_h_jMET,            back_h_jMET,            "jFEX MET",                    "jFEX MET [GeV]",                fDir + "jFEX_MET.pdf");
        drawOverlay(sig_h_metTruthNonInt,  back_h_metTruthNonInt,  "Truth MET (NonInt)",          "Truth MET (NonInt) [GeV]",          fDir + "TruthMET_NonInt.pdf");
        drawOverlay(sig_h_metTruthInt,     back_h_metTruthInt,     "Truth MET (Int)",             "Truth MET (Int) [GeV]",          fDir + "TruthMET_Int.pdf");
        drawOverlay(sig_h_metTruthIntOut,  back_h_metTruthIntOut,  "Truth MET (IntOut)",          "Truth MET (IntOut) [GeV]",          fDir + "TruthMET_IntOut.pdf");
        // --- MET X/Y symmetry plots (mean & median in legend) ---
        drawComponentOverlay(sig_h_gMETX,     back_h_gMETX,     "gFEX MET_{x} (JwoJ)",     "gFEX MET_{x} (JwoJ) [GeV]", fDir + "gFEX_METX_JwoJ.pdf");
        drawComponentOverlay(sig_h_gMETY,     back_h_gMETY,     "gFEX MET_{y} (JwoJ)",     "gFEX MET_{y} (JwoJ) [GeV]", fDir + "gFEX_METY_JwoJ.pdf");
        drawComponentOverlay(sig_h_gMETX_NC,  back_h_gMETX_NC,  "gFEX MET_{x} (NoiseCut)", "gFEX MET_{x} (NoiseCut) [GeV]", fDir + "gFEX_METX_NoiseCut.pdf");
        drawComponentOverlay(sig_h_gMETY_NC,  back_h_gMETY_NC,  "gFEX MET_{y} (NoiseCut)", "gFEX MET_{y} (NoiseCut) [GeV]", fDir + "gFEX_METY_NoiseCut.pdf");
        drawComponentOverlay(sig_h_gMETX_Rms, back_h_gMETX_Rms, "gFEX MET_{x} (Rms)",      "gFEX MET_{x} (rms) [GeV]", fDir + "gFEX_METX_Rms.pdf");
        drawComponentOverlay(sig_h_gMETY_Rms, back_h_gMETY_Rms, "gFEX MET_{y} (Rms)",      "gFEX MET_{y} (rms) [GeV]", fDir + "gFEX_METY_Rms.pdf");
        if (hasTruthNonIntXY) {
            drawComponentOverlay(sig_h_metTruthNonIntX, back_h_metTruthNonIntX, "Truth MET_{x} (NonInt)", "Truth MET_{x} (NonInt) [GeV]", fDir + "TruthMET_NonIntX.pdf");
            drawComponentOverlay(sig_h_metTruthNonIntY, back_h_metTruthNonIntY, "Truth MET_{y} (NonInt)", "Truth MET_{y} (NonInt) [GeV]", fDir + "TruthMET_NonIntY.pdf");
        }
        if (hasJetMetXY) {
            drawComponentOverlay(sig_h_JetMetX, back_h_JetMetX, "GEP Jet MET_{x}", "GEP Jet MET_{x} [GeV]", fDir + "JetMETx.pdf");
            drawComponentOverlay(sig_h_JetMetY, back_h_JetMetY, "GEP Jet MET_{y}", "GEP Jet MET_{y} [GeV]", fDir + "JetMETy.pdf");
        }
        if (hasTowerMetXY) {
            drawComponentOverlay(sig_h_TowerMetX, back_h_TowerMetX, "GEP Tower MET_{x}", "GEP Tower MET_{x} [GeV]", fDir + "TowerMETx.pdf");
            drawComponentOverlay(sig_h_TowerMetY, back_h_TowerMetY, "GEP Tower MET_{y}", "GEP Tower MET_{y} [GeV]", fDir + "TowerMETy.pdf");
        }
        // --- Per-file GEP vs gFEX MET comparison ---
        drawAlgoComparison(sig_h_TotalMET, back_h_TotalMET, sig_h_gMET, back_h_gMET,
                           "Jet Tagger (GEP)", "gFEX",
                           "MET [GeV]", fDir + "GEP_vs_gFEX_MET.pdf", signalName);

        // --- Per-file core term sig vs bkg overlays ---
        drawOverlay(sig_h_coreEMTopo_SoftClus_MET,   back_h_coreEMTopo_SoftClus_MET,   "Core EMTopo SoftClus MET",   "Core EMTopo SoftClus MET [GeV]", fDir + "CoreEMTopo_SoftClus_MET.pdf");
        drawOverlay(sig_h_coreEMTopo_PVSoftTrk_MET,  back_h_coreEMTopo_PVSoftTrk_MET,  "Core EMTopo PVSoftTrk MET",  "Core EMTopo PVSoftTrk MET [GeV]", fDir + "CoreEMTopo_PVSoftTrk_MET.pdf");
        drawOverlay(sig_h_coreEMTopo_SoftClusEM_MET, back_h_coreEMTopo_SoftClusEM_MET, "Core EMTopo SoftClusEM MET", "Core EMTopo SoftClusEM MET [GeV]", fDir + "CoreEMTopo_SoftClusEM_MET.pdf");
        drawOverlay(sig_h_coreEMPFlow_SoftClus_MET,  back_h_coreEMPFlow_SoftClus_MET,  "Core EMPFlow SoftClus MET",  "Core EMPFlow SoftClus MET [GeV]", fDir + "CoreEMPFlow_SoftClus_MET.pdf");
        drawOverlay(sig_h_coreEMPFlow_PVSoftTrk_MET, back_h_coreEMPFlow_PVSoftTrk_MET, "Core EMPFlow PVSoftTrk MET", "Core EMPFlow PVSoftTrk MET [GeV]", fDir + "CoreEMPFlow_PVSoftTrk_MET.pdf");

        // --- Multi-MET comparison: all types on one canvas, signal ---
        drawMultiDist({sig_h_TotalMET, sig_h_gMET, sig_h_metTruthNonInt,
                       sig_h_coreEMTopo_SoftClus_MET, sig_h_coreEMTopo_PVSoftTrk_MET, sig_h_coreEMTopo_SoftClusEM_MET,
                       sig_h_coreEMPFlow_SoftClus_MET, sig_h_coreEMPFlow_PVSoftTrk_MET},
                      {"GEP TotalMET", "gFEX MET", "Truth NonInt",
                       "EMTopo SoftClus", "EMTopo PVSoftTrk", "EMTopo SoftClusEM",
                       "EMPFlow SoftClus", "EMPFlow PVSoftTrk"},
                      "MET comparison — signal", "MET [GeV]", fDir + "METComparison_sig.pdf");

        // --- Multi-MET comparison: all types on one canvas, background ---
        drawMultiDist({back_h_TotalMET, back_h_gMET, back_h_metTruthNonInt,
                       back_h_coreEMTopo_SoftClus_MET, back_h_coreEMTopo_PVSoftTrk_MET, back_h_coreEMTopo_SoftClusEM_MET,
                       back_h_coreEMPFlow_SoftClus_MET, back_h_coreEMPFlow_PVSoftTrk_MET},
                      {"GEP TotalMET", "gFEX MET", "Truth NonInt",
                       "EMTopo SoftClus", "EMTopo PVSoftTrk", "EMTopo SoftClusEM",
                       "EMPFlow SoftClus", "EMPFlow PVSoftTrk"},
                      "MET comparison — background", "MET [GeV]", fDir + "METComparison_bkg.pdf");

        // --- Event properties per 80 kHz selection ---
        {
            std::string selDir = fDir + "SelectionProperties/";
            gSystem->mkdir(selDir.c_str(), true);
            std::vector<std::string> selLegLabels = {
                "Before Selection",
                Form("gFEX JwoJ MET > %.0f GeV",  thr_gMET_80kHz),
                Form("gFEX NC MET > %.0f GeV",     thr_gMET_NC_80kHz),
                Form("gFEX Rms MET > %.0f GeV",    thr_gMET_Rms_80kHz),
                Form("GEP Jet MET > %.0f GeV",     thr_JetMET_80kHz),
                Form("GEP Tower MET > %.0f GeV",   thr_TowerMET_80kHz)
            };
            auto toVec = [nSel80](TH1F** a) { return std::vector<TH1F*>(a, a + nSel80); };
            // Signal — shape comparison across selections
            drawMultiDist(toVec(sig_sel_truthMET),       selLegLabels, "Truth MET_{NonInt} — signal at 80 kHz sel.",                    "Truth MET_{NonInt} [GeV]",                                           selDir + "sig_sel_truthMET.pdf");
            drawMultiDist(toVec(sig_sel_SumET),          selLegLabels, "GEP TOB #Sigma E_{T} — signal at 80 kHz sel.",                  "GEP TOB #Sigma E_{T} [GeV]",                                         selDir + "sig_sel_SumET.pdf");
            drawMultiDist(toVec(sig_sel_gSumET),         selLegLabels, "gFEX TOB #Sigma E_{T} — signal at 80 kHz sel.",                 "gFEX TOB #Sigma E_{T} [GeV]",                                        selDir + "sig_sel_gSumET.pdf");
            drawMultiDist(toVec(sig_sel_METsig),         selLegLabels, "GEP MET significance — signal at 80 kHz sel.",                  "GEP MET / #sqrt{GEP TOB #Sigma E_{T}} [#sqrt{GeV}]",                 selDir + "sig_sel_METsig.pdf",         true,  "#sqrt{GeV}");
            drawMultiDist(toVec(sig_sel_gMETsig),        selLegLabels, "gFEX MET significance — signal at 80 kHz sel.",                 "gFEX MET / #sqrt{gFEX TOB #Sigma E_{T}} [#sqrt{GeV}]",               selDir + "sig_sel_gMETsig.pdf",        true,  "#sqrt{GeV}");
            drawMultiDist(toVec(sig_sel_dPhi_GEP_gFEX), selLegLabels, "#Delta#phi(GEP MET, gFEX MET) — signal at 80 kHz sel.",         "|#Delta#phi(GEP MET, gFEX MET)| [rad]",                              selDir + "sig_sel_dPhi_GEP_gFEX.pdf",  false, "rad");
            if (hasTruthNonIntXY) {
                drawMultiDist(toVec(sig_sel_dPhi_GEP_truth), selLegLabels, "#Delta#phi(GEP MET, truth MET) — signal at 80 kHz sel.",    "|#Delta#phi(GEP MET, truth MET)| [rad]",                             selDir + "sig_sel_dPhi_GEP_truth.pdf",  false, "rad");
                drawMultiDist(toVec(sig_sel_dPhi_truth_TOB), selLegLabels, "#Delta#phi(truth MET, TOB MET) — signal at 80 kHz sel.",    "|#Delta#phi(truth MET, TOB MET)| [rad]",                             selDir + "sig_sel_dPhi_truth_TOB.pdf",  false, "rad");
                for (int iSel = 0; iSel < nSel80; iSel++) {
                    TCanvas c2D(("c_sig_phi2D_"+std::string(selShortNames[iSel])).c_str(), "", 700, 600);
                    c2D.SetRightMargin(0.15);
                    sig_sel_phi2D_TOB_truth[iSel]->SetTitle(Form("Signal — %s;Truth MET #phi [rad];TOB MET #phi [rad]", selLegLabels[iSel].c_str()));
                    sig_sel_phi2D_TOB_truth[iSel]->Draw("COLZ");
                    gPad->SetLogz(sig_sel_phi2D_TOB_truth[iSel]->GetMaximum() > 0);
                    c2D.SaveAs((selDir + "sig_sel_phi2D_TOB_truth_" + selShortNames[iSel] + ".pdf").c_str());
                }
            }
            if (hasTruthAntiKt4WZDressed)
                drawMultiDist(toVec(sig_sel_nTruthJets), selLegLabels, "Truth AntiKt4 WZ-dressed jet multiplicity (E_{T}>15 GeV) — signal at 80 kHz sel.", "N_{truth AntiKt4 WZ-dressed jets}",  selDir + "sig_sel_nTruthJets.pdf", false, "");
            if (hasInTimeAntiKt4TruthJets)
                drawMultiDist(toVec(sig_sel_nPileupJets), selLegLabels, "Truth Pileup Jets with E_{T} #geq 15 GeV — signal at 80 kHz sel.", "N_{truth pileup jets} (E_{T} #geq 15 GeV)", selDir + "sig_sel_nPileupJets.pdf", false, "");
            if (hasWTAConeJets) {
                drawMultiDist(toVec(sig_sel_nJets),         selLegLabels, "WTA-cone jet multiplicity (p_{T}>25 GeV) — signal at 80 kHz sel.",              "N_{WTA-cone jets}",                              selDir + "sig_sel_nJets.pdf",         false, "");
                drawMultiDist(toVec(sig_sel_jet1pt),        selLegLabels, "Leading WTA-cone jet p_{T} — signal at 80 kHz sel.",                            "Leading WTA-cone jet p_{T} [GeV]",               selDir + "sig_sel_jet1pt.pdf");
                drawMultiDist(toVec(sig_sel_jet2pt),        selLegLabels, "Subleading WTA-cone jet p_{T} — signal at 80 kHz sel.",                         "Subleading WTA-cone jet p_{T} [GeV]",            selDir + "sig_sel_jet2pt.pdf");
                drawMultiDist(toVec(sig_sel_jet1eta),       selLegLabels, "Leading WTA-cone jet #eta — signal at 80 kHz sel.",                             "Leading WTA-cone jet #eta",                      selDir + "sig_sel_jet1eta.pdf",       false, "");
                drawMultiDist(toVec(sig_sel_dPhi_jet1_MET), selLegLabels, "#Delta#phi(leading WTA-cone jet, GEP MET) — signal at 80 kHz sel.",             "|#Delta#phi(leading WTA-cone jet, GEP MET)| [rad]", selDir + "sig_sel_dPhi_jet1_MET.pdf", false, "rad");
            }
            // Background — fake-trigger topology
            drawMultiDist(toVec(back_sel_truthMET),      selLegLabels, "Truth MET_{NonInt} — bkg at 80 kHz sel.",                       "Truth MET_{NonInt} [GeV]",                                           selDir + "back_sel_truthMET.pdf");
            drawMultiDist(toVec(back_sel_SumET),         selLegLabels, "GEP TOB #Sigma E_{T} — bkg at 80 kHz sel.",                     "GEP TOB #Sigma E_{T} [GeV]",                                         selDir + "back_sel_SumET.pdf");
            drawMultiDist(toVec(back_sel_gSumET),        selLegLabels, "gFEX TOB #Sigma E_{T} — bkg at 80 kHz sel.",                    "gFEX TOB #Sigma E_{T} [GeV]",                                        selDir + "back_sel_gSumET.pdf");
            drawMultiDist(toVec(back_sel_METsig),        selLegLabels, "GEP MET significance — bkg at 80 kHz sel.",                     "GEP MET / #sqrt{GEP TOB #Sigma E_{T}} [#sqrt{GeV}]",                 selDir + "back_sel_METsig.pdf",        true,  "#sqrt{GeV}");
            drawMultiDist(toVec(back_sel_gMETsig),       selLegLabels, "gFEX MET significance — bkg at 80 kHz sel.",                    "gFEX MET / #sqrt{gFEX TOB #Sigma E_{T}} [#sqrt{GeV}]",               selDir + "back_sel_gMETsig.pdf",       true,  "#sqrt{GeV}");
            drawMultiDist(toVec(back_sel_dPhi_GEP_gFEX),selLegLabels, "#Delta#phi(GEP MET, gFEX MET) — bkg at 80 kHz sel.",            "|#Delta#phi(GEP MET, gFEX MET)| [rad]",                              selDir + "back_sel_dPhi_GEP_gFEX.pdf", false, "rad");
            if (hasTruthNonIntXY) {
                drawMultiDist(toVec(back_sel_dPhi_GEP_truth), selLegLabels, "#Delta#phi(GEP MET, truth MET) — bkg at 80 kHz sel.",      "|#Delta#phi(GEP MET, truth MET)| [rad]",                             selDir + "back_sel_dPhi_GEP_truth.pdf", false, "rad");
                drawMultiDist(toVec(back_sel_dPhi_truth_TOB), selLegLabels, "#Delta#phi(truth MET, TOB MET) — bkg at 80 kHz sel.",      "|#Delta#phi(truth MET, TOB MET)| [rad]",                             selDir + "back_sel_dPhi_truth_TOB.pdf", false, "rad");
                for (int iSel = 0; iSel < nSel80; iSel++) {
                    TCanvas c2D(("c_back_phi2D_"+std::string(selShortNames[iSel])).c_str(), "", 700, 600);
                    c2D.SetRightMargin(0.15);
                    back_sel_phi2D_TOB_truth[iSel]->SetTitle(Form("Background — %s;Truth MET #phi [rad];TOB MET #phi [rad]", selLegLabels[iSel].c_str()));
                    back_sel_phi2D_TOB_truth[iSel]->Draw("COLZ");
                    gPad->SetLogz(back_sel_phi2D_TOB_truth[iSel]->GetMaximum() > 0);
                    c2D.SaveAs((selDir + "back_sel_phi2D_TOB_truth_" + selShortNames[iSel] + ".pdf").c_str());
                }
            }
            if (hasTruthAntiKt4WZDressed)
                drawMultiDist(toVec(back_sel_nTruthJets), selLegLabels, "Truth AntiKt4 WZ-dressed jet multiplicity (E_{T}>15 GeV) — bkg at 80 kHz sel.",  "N_{truth AntiKt4 WZ-dressed jets}",  selDir + "back_sel_nTruthJets.pdf", false, "");
            if (hasInTimeAntiKt4TruthJets)
                drawMultiDist(toVec(back_sel_nPileupJets), selLegLabels, "Truth Pileup Jets with E_{T} #geq 15 GeV — bkg at 80 kHz sel.", "N_{truth pileup jets} (E_{T} #geq 15 GeV)", selDir + "back_sel_nPileupJets.pdf", false, "");
            if (hasWTAConeJets) {
                drawMultiDist(toVec(back_sel_nJets),         selLegLabels, "WTA-cone jet multiplicity (p_{T}>25 GeV) — bkg at 80 kHz sel.",               "N_{WTA-cone jets}",                              selDir + "back_sel_nJets.pdf",         false, "");
                drawMultiDist(toVec(back_sel_jet1pt),        selLegLabels, "Leading WTA-cone jet p_{T} — bkg at 80 kHz sel.",                             "Leading WTA-cone jet p_{T} [GeV]",               selDir + "back_sel_jet1pt.pdf");
                drawMultiDist(toVec(back_sel_jet2pt),        selLegLabels, "Subleading WTA-cone jet p_{T} — bkg at 80 kHz sel.",                          "Subleading WTA-cone jet p_{T} [GeV]",            selDir + "back_sel_jet2pt.pdf");
                drawMultiDist(toVec(back_sel_jet1eta),       selLegLabels, "Leading WTA-cone jet #eta — bkg at 80 kHz sel.",                              "Leading WTA-cone jet #eta",                      selDir + "back_sel_jet1eta.pdf",       false, "");
                drawMultiDist(toVec(back_sel_dPhi_jet1_MET), selLegLabels, "#Delta#phi(leading WTA-cone jet, GEP MET) — bkg at 80 kHz sel.",              "|#Delta#phi(leading WTA-cone jet, GEP MET)| [rad]", selDir + "back_sel_dPhi_jet1_MET.pdf", false, "rad");
            }
        }

        // --- Per-JZ-slice MET distributions ---
        TString jzDir = (fDir + "JZSlices/").c_str();
        gSystem->mkdir(jzDir);
        OverlayAndSave(back_h_gMET_jz,     nJZSlices_, "c_gMET_jz",     jzDir + "gFEX_MET_JZSlices.pdf",        0);
        OverlayAndSave(back_h_gMET_NC_jz,  nJZSlices_, "c_gMET_NC_jz",  jzDir + "gFEX_NoiseCut_MET_JZSlices.pdf", 0);
        OverlayAndSave(back_h_gMET_Rms_jz, nJZSlices_, "c_gMET_Rms_jz", jzDir + "gFEX_RMS_MET_JZSlices.pdf",     0);
        OverlayAndSave(back_h_jMET_jz,     nJZSlices_, "c_jMET_jz",     jzDir + "jFEX_MET_JZSlices.pdf",        0);
        OverlayAndSave(back_h_JetMET_jz,   nJZSlices_, "c_JetMET_jz",   jzDir + "GEP_JetMET_JZSlices.pdf",     0);
        OverlayAndSave(back_h_TowerMET_jz, nJZSlices_, "c_TowerMET_jz", jzDir + "GEP_TowerMET_JZSlices.pdf",   0);

        // --- Per-file rate vs threshold plots ---
        drawRateVsThreshold(back_hw_TotalMET,  "Rate vs Emulated MET threshold",          "MET threshold [GeV]", fDir + "Rate_TotalMET.pdf");
        drawRateVsThreshold(back_hw_gMET,      "Rate vs gFEX MET threshold (JwoJ)",       "MET threshold [GeV]", fDir + "Rate_gFEX_MET_JwoJ.pdf");
        drawRateVsThreshold(back_hw_gMET_NC,   "Rate vs gFEX MET threshold (NoiseCut)",   "MET threshold [GeV]", fDir + "Rate_gFEX_MET_NoiseCut.pdf");
        drawRateVsThreshold(back_hw_gMET_Rms,  "Rate vs gFEX MET threshold (Rms)",        "MET threshold [GeV]", fDir + "Rate_gFEX_MET_Rms.pdf");
        drawRateVsThreshold(back_hw_jMET,      "Rate vs jFEX MET threshold",              "MET threshold [GeV]", fDir + "Rate_jFEX_MET.pdf");
        drawRateVsThreshold(back_hw_JetMET,    "Rate vs GEP Jet MET threshold",           "MET threshold [GeV]", fDir + "Rate_JetMET.pdf");
        drawRateVsThreshold(back_hw_TowerMET,  "Rate vs GEP Tower MET threshold",         "MET threshold [GeV]", fDir + "Rate_TowerMET.pdf");

        // --- Per-file signal efficiency vs threshold plots ---
        {
            std::vector<TH1F*> effHists = hasOverlapRemoval
                ? std::vector<TH1F*>{sig_h_gMET, sig_h_gMET_NC, sig_h_gMET_Rms, sig_h_jMET, sig_h_JetMet, sig_h_TotalMET}
                : std::vector<TH1F*>{sig_h_gMET, sig_h_gMET_NC, sig_h_gMET_Rms, sig_h_jMET, sig_h_JetMet, sig_h_TowerMet, sig_h_TotalMET};
            std::vector<std::string> effLabels = hasOverlapRemoval
                ? std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Total MET"}
                : std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
            drawEffVsThresholdMulti(effHists, effLabels,
                                    "Signal Efficiency vs MET Threshold", "MET threshold [GeV]",
                                    fDir + "SigEff_vs_Threshold_AlgoComparison.pdf", signalName);
            // gFEX-only version
            std::vector<TH1F*> gfexHists  = {sig_h_gMET, sig_h_gMET_NC, sig_h_gMET_Rms};
            std::vector<std::string> gfexL = {"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms"};
            drawEffVsThresholdMulti(gfexHists, gfexL,
                                    "Signal Efficiency vs gFEX MET Threshold", "MET threshold [GeV]",
                                    fDir + "SigEff_vs_Threshold_gFEX_AlgoComparison.pdf", signalName);
            // GEP-only version
            std::vector<TH1F*> gepHists = hasOverlapRemoval
                ? std::vector<TH1F*>{sig_h_JetMet, sig_h_TotalMET}
                : std::vector<TH1F*>{sig_h_JetMet, sig_h_TowerMet, sig_h_TotalMET};
            std::vector<std::string> gepL = hasOverlapRemoval
                ? std::vector<std::string>{"GEP Jet MET", "GEP Total MET"}
                : std::vector<std::string>{"GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
            drawEffVsThresholdMulti(gepHists, gepL,
                                    "Signal Efficiency vs GEP MET Threshold", "MET threshold [GeV]",
                                    fDir + "SigEff_vs_Threshold_GEP_AlgoComparison.pdf", signalName);
        }

        // --- TOB MET vs truth NonInt MET calibration and resolution (signal and background) ---
        {
            std::string calDir = fDir + "Calibration/";
            gSystem->mkdir(calDir.c_str(), true);

            // 2D correlation: x = truth NonInt MET, y = TOB MET
            // zmin: colour scale minimum; xmax_cap: x-axis display cap (0 = use histogram range)
            auto drawTOBvsTruth = [&](TH2F* h, const std::string& path, const std::string& tobLabel,
                                      double zmin = 1e-5, double xmax_cap = 0.0) {
                if (h->Integral() <= 0) return;
                h->Scale(1.0 / h->Integral());
                h->SetMinimum(zmin);
                h->GetXaxis()->SetTitle("Truth MET_{NonInt} [GeV]");
                h->GetYaxis()->SetTitle((tobLabel + " MET [GeV]").c_str());
                if (xmax_cap > 0) h->GetXaxis()->SetRangeUser(0, xmax_cap);
                TCanvas cTmp(("cCalCorr_"+std::string(h->GetName())).c_str(), "", 700, 600);
                cTmp.SetLogz();
                cTmp.SetRightMargin(0.15); cTmp.SetLeftMargin(0.14);
                cTmp.SetBottomMargin(0.14); cTmp.SetTicks(1, 1);
                h->Draw("COLZ");
                double xmax = (xmax_cap > 0) ? xmax_cap : h->GetXaxis()->GetXmax();
                double ymax = h->GetYaxis()->GetXmax();
                double r = h->GetCorrelationFactor(1, 2);
                TProfile* prof = h->ProfileX((std::string(h->GetName())+"_pfx").c_str(), 1, -1, "s");
                prof->SetMarkerStyle(20); prof->SetMarkerSize(0.5);
                prof->SetMarkerColor(kBlack); prof->SetLineColor(kBlack);
                TF1* fitFn = new TF1((std::string("lf_")+h->GetName()).c_str(), "pol1", 0, xmax);
                prof->Fit(fitFn, "QN");
                double slope = fitFn->GetParameter(1), intercept = fitFn->GetParameter(0);
                fitFn->SetLineColor(kBlue); fitFn->SetLineWidth(2);
                fitFn->Draw("SAME"); prof->Draw("SAME");
                double rng = std::min(xmax, ymax);
                TLine* diag = new TLine(0, 0, rng, rng);
                diag->SetLineColor(kRed); diag->SetLineStyle(2); diag->SetLineWidth(2);
                diag->Draw("SAME");
                TLegend leg(0.16, 0.62, 0.61, 0.88);
                leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);
                leg.AddEntry(diag, "TOB MET = Truth MET", "l");
                leg.AddEntry(fitFn, Form("Fit: y = %.3f x + %.1f GeV", slope, intercept), "l");
                leg.AddEntry((TObject*)nullptr, Form("r = %.4f", r), "");
                leg.Draw();
                cTmp.SaveAs(path.c_str());
                delete prof; delete fitFn;
            };
            // Signal
            drawTOBvsTruth(sig_h2_gJwoJ_TOBMet_vs_truthMET,    calDir + "sig_gFEX_JwoJ_TOBMet_vs_truthMET.pdf",    "gFEX JwoJ");
            drawTOBvsTruth(sig_h2_gNC_TOBMet_vs_truthMET,      calDir + "sig_gFEX_NoiseCut_TOBMet_vs_truthMET.pdf","gFEX NoiseCut");
            drawTOBvsTruth(sig_h2_gRms_TOBMet_vs_truthMET,     calDir + "sig_gFEX_Rms_TOBMet_vs_truthMET.pdf",     "gFEX Rms");
            drawTOBvsTruth(sig_h2_JetMET_TOBMet_vs_truthMET,   calDir + "sig_GEP_JetMET_TOBMet_vs_truthMET.pdf",   "GEP Jet");
            drawTOBvsTruth(sig_h2_TowerMET_TOBMet_vs_truthMET, calDir + "sig_GEP_TowerMET_TOBMet_vs_truthMET.pdf", "GEP Tower");
            // Background (lower z-floor to show high-MET tails; x capped at 300 GeV)
            drawTOBvsTruth(back_h2_gJwoJ_TOBMet_vs_truthMET,    calDir + "back_gFEX_JwoJ_TOBMet_vs_truthMET.pdf",    "gFEX JwoJ",    1e-14, 300.0);
            drawTOBvsTruth(back_h2_gNC_TOBMet_vs_truthMET,      calDir + "back_gFEX_NoiseCut_TOBMet_vs_truthMET.pdf","gFEX NoiseCut",1e-14, 300.0);
            drawTOBvsTruth(back_h2_gRms_TOBMet_vs_truthMET,     calDir + "back_gFEX_Rms_TOBMet_vs_truthMET.pdf",     "gFEX Rms",     1e-14, 300.0);
            drawTOBvsTruth(back_h2_JetMET_TOBMet_vs_truthMET,   calDir + "back_GEP_JetMET_TOBMet_vs_truthMET.pdf",   "GEP Jet",      1e-14, 300.0);
            drawTOBvsTruth(back_h2_TowerMET_TOBMet_vs_truthMET, calDir + "back_GEP_TowerMET_TOBMet_vs_truthMET.pdf", "GEP Tower",    1e-14, 300.0);

            // 1D residual distributions, individual and overlaid
            // xLbl: x-axis label; yUnit: appended to y-axis bin-size string (e.g., "" or "GeV")
            auto drawRelResidual1D = [&](TH1F* h, const std::string& path, const std::string& tobLabel,
                                         const std::string& xLbl = "(Truth - TOB) / Truth MET_{NonInt}",
                                         const std::string& yUnit = "") {
                if (h->Integral() <= 0) return;
                h->Scale(1.0 / h->Integral());
                TCanvas cTmp(("cRes1D_"+std::string(h->GetName())).c_str(), "", 700, 600);
                cTmp.SetLeftMargin(0.14); cTmp.SetBottomMargin(0.14); cTmp.SetTicks(1, 1);
                cTmp.SetLogy();
                h->SetLineColor(kBlue); h->SetLineWidth(2);
                std::string yT = Form("Fraction of Events / %.4g", h->GetBinWidth(1));
                if (!yUnit.empty()) yT += " " + yUnit;
                h->SetTitle((tobLabel + ";" + xLbl + ";" + yT).c_str());
                h->SetMinimum(1e-5);
                h->Draw("HIST");
                TLine* zero = new TLine(0, h->GetMinimum(), 0, h->GetMaximum() * 1.5);
                zero->SetLineColor(kRed); zero->SetLineStyle(2); zero->SetLineWidth(2);
                zero->Draw("SAME");
                cTmp.SaveAs(path.c_str());
            };
            // Helper to draw the 5-algorithm overlay for a given set of residual histograms
            auto drawRelResidualOverlay = [&](std::vector<TH1F*> resH, const std::string& path,
                                              const std::string& xLbl = "(Truth - TOB) / Truth MET_{NonInt}",
                                              const std::string& yUnit = "") {
                const Color_t rcols[] = {kBlack, kRed, kOrange+1, kBlue, kGreen+2};
                std::vector<std::string> resL = {"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms",
                                                 "GEP Jet MET", "GEP Tower MET"};
                double ymaxR = 0;
                for (auto* h : resH) ymaxR = std::max(ymaxR, h->GetMaximum());
                ymaxR *= 5.0;
                std::string resYTitle = Form("Fraction of Events / %.4g", resH[0]->GetBinWidth(1));
                if (!yUnit.empty()) resYTitle += " " + yUnit;
                TCanvas cOvl(("cResOvl_"+path).c_str(), "", 800, 600);
                cOvl.SetLeftMargin(0.13); cOvl.SetBottomMargin(0.14); cOvl.SetTicks(1, 1);
                cOvl.SetLogy();
                TLegend legR(0.6, 0.62, 0.88, 0.88);
                legR.SetBorderSize(0); legR.SetFillStyle(0); legR.SetTextSize(0.030);
                for (unsigned int i = 0; i < resH.size(); i++) {
                    resH[i]->SetLineColor(rcols[i]); resH[i]->SetLineWidth(2);
                    resH[i]->SetMaximum(ymaxR); resH[i]->SetMinimum(1e-5);
                    resH[i]->GetXaxis()->SetTitle(xLbl.c_str());
                    resH[i]->GetYaxis()->SetTitle(resYTitle.c_str());
                    resH[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
                    legR.AddEntry(resH[i], resL[i].c_str(), "l");
                }
                TLine* zeroR = new TLine(0, 1e-5, 0, ymaxR);
                zeroR->SetLineColor(kGray+2); zeroR->SetLineStyle(2); zeroR->SetLineWidth(2);
                zeroR->Draw("SAME");
                legR.Draw();
                cOvl.SaveAs(path.c_str());
            };
            // Signal relative residuals
            drawRelResidual1D(sig_h1_gJwoJ_relResidual,    calDir + "sig_gFEX_JwoJ_relResidual.pdf",    "gFEX JwoJ");
            drawRelResidual1D(sig_h1_gNC_relResidual,      calDir + "sig_gFEX_NoiseCut_relResidual.pdf","gFEX NoiseCut");
            drawRelResidual1D(sig_h1_gRms_relResidual,     calDir + "sig_gFEX_Rms_relResidual.pdf",     "gFEX Rms");
            drawRelResidual1D(sig_h1_JetMET_relResidual,   calDir + "sig_GEP_JetMET_relResidual.pdf",   "GEP Jet MET");
            drawRelResidual1D(sig_h1_TowerMET_relResidual, calDir + "sig_GEP_TowerMET_relResidual.pdf", "GEP Tower MET");
            drawRelResidualOverlay({sig_h1_gJwoJ_relResidual, sig_h1_gNC_relResidual, sig_h1_gRms_relResidual,
                                    sig_h1_JetMET_relResidual, sig_h1_TowerMET_relResidual},
                                   calDir + "sig_relResidual_overlay.pdf");
            // Signal absolute residuals
            drawRelResidual1D(sig_h1_gJwoJ_absResidual,    calDir + "sig_gFEX_JwoJ_absResidual.pdf",    "gFEX JwoJ",    "Truth - TOB MET [GeV]", "GeV");
            drawRelResidual1D(sig_h1_gNC_absResidual,      calDir + "sig_gFEX_NoiseCut_absResidual.pdf","gFEX NoiseCut","Truth - TOB MET [GeV]", "GeV");
            drawRelResidual1D(sig_h1_gRms_absResidual,     calDir + "sig_gFEX_Rms_absResidual.pdf",     "gFEX Rms",     "Truth - TOB MET [GeV]", "GeV");
            drawRelResidual1D(sig_h1_JetMET_absResidual,   calDir + "sig_GEP_JetMET_absResidual.pdf",   "GEP Jet MET",  "Truth - TOB MET [GeV]", "GeV");
            drawRelResidual1D(sig_h1_TowerMET_absResidual, calDir + "sig_GEP_TowerMET_absResidual.pdf", "GEP Tower MET","Truth - TOB MET [GeV]", "GeV");
            drawRelResidualOverlay({sig_h1_gJwoJ_absResidual, sig_h1_gNC_absResidual, sig_h1_gRms_absResidual,
                                    sig_h1_JetMET_absResidual, sig_h1_TowerMET_absResidual},
                                   calDir + "sig_absResidual_overlay.pdf",
                                   "Truth - TOB MET [GeV]", "GeV");
            // Background relative residuals
            drawRelResidual1D(back_h1_gJwoJ_relResidual,    calDir + "back_gFEX_JwoJ_relResidual.pdf",    "gFEX JwoJ");
            drawRelResidual1D(back_h1_gNC_relResidual,      calDir + "back_gFEX_NoiseCut_relResidual.pdf","gFEX NoiseCut");
            drawRelResidual1D(back_h1_gRms_relResidual,     calDir + "back_gFEX_Rms_relResidual.pdf",     "gFEX Rms");
            drawRelResidual1D(back_h1_JetMET_relResidual,   calDir + "back_GEP_JetMET_relResidual.pdf",   "GEP Jet MET");
            drawRelResidual1D(back_h1_TowerMET_relResidual, calDir + "back_GEP_TowerMET_relResidual.pdf", "GEP Tower MET");
            drawRelResidualOverlay({back_h1_gJwoJ_relResidual, back_h1_gNC_relResidual, back_h1_gRms_relResidual,
                                    back_h1_JetMET_relResidual, back_h1_TowerMET_relResidual},
                                   calDir + "back_relResidual_overlay.pdf");
            // Background absolute residuals
            drawRelResidual1D(back_h1_gJwoJ_absResidual,    calDir + "back_gFEX_JwoJ_absResidual.pdf",    "gFEX JwoJ",    "Truth - TOB MET [GeV]", "GeV");
            drawRelResidual1D(back_h1_gNC_absResidual,      calDir + "back_gFEX_NoiseCut_absResidual.pdf","gFEX NoiseCut","Truth - TOB MET [GeV]", "GeV");
            drawRelResidual1D(back_h1_gRms_absResidual,     calDir + "back_gFEX_Rms_absResidual.pdf",     "gFEX Rms",     "Truth - TOB MET [GeV]", "GeV");
            drawRelResidual1D(back_h1_JetMET_absResidual,   calDir + "back_GEP_JetMET_absResidual.pdf",   "GEP Jet MET",  "Truth - TOB MET [GeV]", "GeV");
            drawRelResidual1D(back_h1_TowerMET_absResidual, calDir + "back_GEP_TowerMET_absResidual.pdf", "GEP Tower MET","Truth - TOB MET [GeV]", "GeV");
            drawRelResidualOverlay({back_h1_gJwoJ_absResidual, back_h1_gNC_absResidual, back_h1_gRms_absResidual,
                                    back_h1_JetMET_absResidual, back_h1_TowerMET_absResidual},
                                   calDir + "back_absResidual_overlay.pdf",
                                   "Truth - TOB MET [GeV]", "GeV");

            // 2D residual vs x variable (truth MET or SumET), with mean profile overlay
            // zmin: z-axis minimum; xmax_cap: x-axis display cap (0 = use histogram range);
            // yLbl: y-axis label (defaults to relative residual label)
            auto drawRelResidual2D = [&](TH2F* h, const std::string& path, const std::string& xLabel,
                                         double zmin = 1e-5, double xmax_cap = 0.0,
                                         const std::string& yLbl = "(Truth - TOB) / Truth MET_{NonInt}") {
                if (h->Integral() <= 0) return;
                h->Scale(1.0 / h->Integral());
                h->SetMinimum(zmin);
                h->GetXaxis()->SetTitle(xLabel.c_str());
                h->GetYaxis()->SetTitle(yLbl.c_str());
                h->GetYaxis()->SetTitleOffset(1.6);
                h->GetXaxis()->SetTitleOffset(1.2);
                if (xmax_cap > 0) h->GetXaxis()->SetRangeUser(0, xmax_cap);
                TCanvas cTmp(("cRes2D_"+std::string(h->GetName())).c_str(), "", 700, 600);
                cTmp.SetLogz();
                cTmp.SetRightMargin(0.15); cTmp.SetLeftMargin(0.18);
                cTmp.SetBottomMargin(0.16); cTmp.SetTicks(1, 1);
                h->Draw("COLZ");
                TProfile* prof = h->ProfileX((std::string(h->GetName())+"_pfx").c_str(), 1, -1, "s");
                prof->SetMarkerStyle(20); prof->SetMarkerSize(0.5);
                prof->SetMarkerColor(kBlack); prof->SetLineColor(kBlack);
                prof->Draw("SAME");
                double xlo = h->GetXaxis()->GetXmin();
                double xhi = (xmax_cap > 0) ? xmax_cap : h->GetXaxis()->GetXmax();
                TLine* zero = new TLine(xlo, 0, xhi, 0);
                zero->SetLineColor(kRed); zero->SetLineStyle(2); zero->SetLineWidth(2);
                zero->Draw("SAME");
                TLegend leg(0.20, 0.85, 0.58, 0.93);
                leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);
                leg.AddEntry(zero, "Perfect resolution", "l");
                leg.AddEntry(prof, "Mean profile", "lp");
                leg.Draw();
                cTmp.SaveAs(path.c_str());
                delete prof;
            };
            const std::string absYLbl = "Truth - TOB MET [GeV]";
            // Signal: relative residual vs truth NonInt MET
            drawRelResidual2D(sig_h2_gJwoJ_relResidual_vs_truthMET,    calDir + "sig_gFEX_JwoJ_relResidual_vs_truthMET.pdf",    "Truth MET_{NonInt} [GeV]");
            drawRelResidual2D(sig_h2_gNC_relResidual_vs_truthMET,      calDir + "sig_gFEX_NoiseCut_relResidual_vs_truthMET.pdf","Truth MET_{NonInt} [GeV]");
            drawRelResidual2D(sig_h2_gRms_relResidual_vs_truthMET,     calDir + "sig_gFEX_Rms_relResidual_vs_truthMET.pdf",     "Truth MET_{NonInt} [GeV]");
            drawRelResidual2D(sig_h2_JetMET_relResidual_vs_truthMET,   calDir + "sig_GEP_JetMET_relResidual_vs_truthMET.pdf",   "Truth MET_{NonInt} [GeV]");
            drawRelResidual2D(sig_h2_TowerMET_relResidual_vs_truthMET, calDir + "sig_GEP_TowerMET_relResidual_vs_truthMET.pdf", "Truth MET_{NonInt} [GeV]");
            // Signal: relative residual vs TOB SumET
            drawRelResidual2D(sig_h2_gJwoJ_relResidual_vs_sumET,    calDir + "sig_gFEX_JwoJ_relResidual_vs_sumET.pdf",    "gFEX JwoJ #Sigma E_{T} [GeV]");
            drawRelResidual2D(sig_h2_gNC_relResidual_vs_sumET,      calDir + "sig_gFEX_NoiseCut_relResidual_vs_sumET.pdf","gFEX NoiseCut #Sigma E_{T} [GeV]");
            drawRelResidual2D(sig_h2_gRms_relResidual_vs_sumET,     calDir + "sig_gFEX_Rms_relResidual_vs_sumET.pdf",     "gFEX Rms #Sigma E_{T} [GeV]");
            drawRelResidual2D(sig_h2_JetMET_relResidual_vs_sumET,   calDir + "sig_GEP_JetMET_relResidual_vs_sumET.pdf",   "GEP #Sigma E_{T} [GeV]");
            drawRelResidual2D(sig_h2_TowerMET_relResidual_vs_sumET, calDir + "sig_GEP_TowerMET_relResidual_vs_sumET.pdf", "GEP #Sigma E_{T} [GeV]");
            // Signal: absolute residual vs truth NonInt MET
            drawRelResidual2D(sig_h2_gJwoJ_absResidual_vs_truthMET,    calDir + "sig_gFEX_JwoJ_absResidual_vs_truthMET.pdf",    "Truth MET_{NonInt} [GeV]", 1e-5, 0.0, absYLbl);
            drawRelResidual2D(sig_h2_gNC_absResidual_vs_truthMET,      calDir + "sig_gFEX_NoiseCut_absResidual_vs_truthMET.pdf","Truth MET_{NonInt} [GeV]", 1e-5, 0.0, absYLbl);
            drawRelResidual2D(sig_h2_gRms_absResidual_vs_truthMET,     calDir + "sig_gFEX_Rms_absResidual_vs_truthMET.pdf",     "Truth MET_{NonInt} [GeV]", 1e-5, 0.0, absYLbl);
            drawRelResidual2D(sig_h2_JetMET_absResidual_vs_truthMET,   calDir + "sig_GEP_JetMET_absResidual_vs_truthMET.pdf",   "Truth MET_{NonInt} [GeV]", 1e-5, 0.0, absYLbl);
            drawRelResidual2D(sig_h2_TowerMET_absResidual_vs_truthMET, calDir + "sig_GEP_TowerMET_absResidual_vs_truthMET.pdf", "Truth MET_{NonInt} [GeV]", 1e-5, 0.0, absYLbl);
            // Signal: absolute residual vs TOB SumET
            drawRelResidual2D(sig_h2_gJwoJ_absResidual_vs_sumET,    calDir + "sig_gFEX_JwoJ_absResidual_vs_sumET.pdf",    "gFEX JwoJ #Sigma E_{T} [GeV]",    1e-5, 0.0, absYLbl);
            drawRelResidual2D(sig_h2_gNC_absResidual_vs_sumET,      calDir + "sig_gFEX_NoiseCut_absResidual_vs_sumET.pdf","gFEX NoiseCut #Sigma E_{T} [GeV]", 1e-5, 0.0, absYLbl);
            drawRelResidual2D(sig_h2_gRms_absResidual_vs_sumET,     calDir + "sig_gFEX_Rms_absResidual_vs_sumET.pdf",     "gFEX Rms #Sigma E_{T} [GeV]",     1e-5, 0.0, absYLbl);
            drawRelResidual2D(sig_h2_JetMET_absResidual_vs_sumET,   calDir + "sig_GEP_JetMET_absResidual_vs_sumET.pdf",   "GEP #Sigma E_{T} [GeV]",          1e-5, 0.0, absYLbl);
            drawRelResidual2D(sig_h2_TowerMET_absResidual_vs_sumET, calDir + "sig_GEP_TowerMET_absResidual_vs_sumET.pdf", "GEP #Sigma E_{T} [GeV]",          1e-5, 0.0, absYLbl);
            // Background: relative residual vs truth NonInt MET (lower z-floor; x capped at 300 GeV)
            drawRelResidual2D(back_h2_gJwoJ_relResidual_vs_truthMET,    calDir + "back_gFEX_JwoJ_relResidual_vs_truthMET.pdf",    "Truth MET_{NonInt} [GeV]", 1e-14, 300.0);
            drawRelResidual2D(back_h2_gNC_relResidual_vs_truthMET,      calDir + "back_gFEX_NoiseCut_relResidual_vs_truthMET.pdf","Truth MET_{NonInt} [GeV]", 1e-14, 300.0);
            drawRelResidual2D(back_h2_gRms_relResidual_vs_truthMET,     calDir + "back_gFEX_Rms_relResidual_vs_truthMET.pdf",     "Truth MET_{NonInt} [GeV]", 1e-14, 300.0);
            drawRelResidual2D(back_h2_JetMET_relResidual_vs_truthMET,   calDir + "back_GEP_JetMET_relResidual_vs_truthMET.pdf",   "Truth MET_{NonInt} [GeV]", 1e-14, 300.0);
            drawRelResidual2D(back_h2_TowerMET_relResidual_vs_truthMET, calDir + "back_GEP_TowerMET_relResidual_vs_truthMET.pdf", "Truth MET_{NonInt} [GeV]", 1e-14, 300.0);
            // Background: relative residual vs TOB SumET (lower z-floor; full SumET range)
            drawRelResidual2D(back_h2_gJwoJ_relResidual_vs_sumET,    calDir + "back_gFEX_JwoJ_relResidual_vs_sumET.pdf",    "gFEX JwoJ #Sigma E_{T} [GeV]",    1e-14);
            drawRelResidual2D(back_h2_gNC_relResidual_vs_sumET,      calDir + "back_gFEX_NoiseCut_relResidual_vs_sumET.pdf","gFEX NoiseCut #Sigma E_{T} [GeV]", 1e-14);
            drawRelResidual2D(back_h2_gRms_relResidual_vs_sumET,     calDir + "back_gFEX_Rms_relResidual_vs_sumET.pdf",     "gFEX Rms #Sigma E_{T} [GeV]",     1e-14);
            drawRelResidual2D(back_h2_JetMET_relResidual_vs_sumET,   calDir + "back_GEP_JetMET_relResidual_vs_sumET.pdf",   "GEP #Sigma E_{T} [GeV]",          1e-14);
            drawRelResidual2D(back_h2_TowerMET_relResidual_vs_sumET, calDir + "back_GEP_TowerMET_relResidual_vs_sumET.pdf", "GEP #Sigma E_{T} [GeV]",          1e-14);
            // Background: absolute residual vs truth NonInt MET (lower z-floor; x capped at 300 GeV)
            drawRelResidual2D(back_h2_gJwoJ_absResidual_vs_truthMET,    calDir + "back_gFEX_JwoJ_absResidual_vs_truthMET.pdf",    "Truth MET_{NonInt} [GeV]", 1e-14, 300.0, absYLbl);
            drawRelResidual2D(back_h2_gNC_absResidual_vs_truthMET,      calDir + "back_gFEX_NoiseCut_absResidual_vs_truthMET.pdf","Truth MET_{NonInt} [GeV]", 1e-14, 300.0, absYLbl);
            drawRelResidual2D(back_h2_gRms_absResidual_vs_truthMET,     calDir + "back_gFEX_Rms_absResidual_vs_truthMET.pdf",     "Truth MET_{NonInt} [GeV]", 1e-14, 300.0, absYLbl);
            drawRelResidual2D(back_h2_JetMET_absResidual_vs_truthMET,   calDir + "back_GEP_JetMET_absResidual_vs_truthMET.pdf",   "Truth MET_{NonInt} [GeV]", 1e-14, 300.0, absYLbl);
            drawRelResidual2D(back_h2_TowerMET_absResidual_vs_truthMET, calDir + "back_GEP_TowerMET_absResidual_vs_truthMET.pdf", "Truth MET_{NonInt} [GeV]", 1e-14, 300.0, absYLbl);
            // Background: absolute residual vs TOB SumET (lower z-floor; full SumET range)
            drawRelResidual2D(back_h2_gJwoJ_absResidual_vs_sumET,    calDir + "back_gFEX_JwoJ_absResidual_vs_sumET.pdf",    "gFEX JwoJ #Sigma E_{T} [GeV]",    1e-14, 0.0, absYLbl);
            drawRelResidual2D(back_h2_gNC_absResidual_vs_sumET,      calDir + "back_gFEX_NoiseCut_absResidual_vs_sumET.pdf","gFEX NoiseCut #Sigma E_{T} [GeV]", 1e-14, 0.0, absYLbl);
            drawRelResidual2D(back_h2_gRms_absResidual_vs_sumET,     calDir + "back_gFEX_Rms_absResidual_vs_sumET.pdf",     "gFEX Rms #Sigma E_{T} [GeV]",     1e-14, 0.0, absYLbl);
            drawRelResidual2D(back_h2_JetMET_absResidual_vs_sumET,   calDir + "back_GEP_JetMET_absResidual_vs_sumET.pdf",   "GEP #Sigma E_{T} [GeV]",          1e-14, 0.0, absYLbl);
            drawRelResidual2D(back_h2_TowerMET_absResidual_vs_sumET, calDir + "back_GEP_TowerMET_absResidual_vs_sumET.pdf", "GEP #Sigma E_{T} [GeV]",          1e-14, 0.0, absYLbl);
        }

        // --- gFEX algorithm comparison (JwoJ vs NoiseCut vs Rms), jFEX, and GEP types overlaid ---
        {
            std::vector<TH1F*> rateHists = hasOverlapRemoval
                ? std::vector<TH1F*>{back_hw_gMET, back_hw_gMET_NC, back_hw_gMET_Rms, back_hw_jMET, back_hw_JetMET, back_hw_TotalMET}
                : std::vector<TH1F*>{back_hw_gMET, back_hw_gMET_NC, back_hw_gMET_Rms, back_hw_jMET, back_hw_JetMET, back_hw_TowerMET, back_hw_TotalMET};
            std::vector<std::string> rateLabels = hasOverlapRemoval
                ? std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Total MET"}
                : std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
            drawRateVsThresholdMulti(rateHists, rateLabels,
                                     "Rate vs MET threshold", "MET threshold [GeV]",
                                     fDir + "Rate_AlgoComparison.pdf", "",
                                     "Estimated Background Rate [Hz]");
        }
        // --- gFEX algorithm-only comparison ---
        {
            std::vector<TH1F*> rateHists = {back_hw_gMET, back_hw_gMET_NC, back_hw_gMET_Rms};
            std::vector<std::string> rateLabels = {"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms"};
            drawRateVsThresholdMulti(rateHists, rateLabels,
                                     "Rate vs gFEX MET threshold", "MET threshold [GeV]",
                                     fDir + "Rate_gFEX_AlgoComparison.pdf", "",
                                     "Estimated Background Rate [Hz]");
        }
        // --- GEP algorithm-only comparison ---
        {
            std::vector<TH1F*> rateHists = hasOverlapRemoval
                ? std::vector<TH1F*>{back_hw_JetMET, back_hw_TotalMET}
                : std::vector<TH1F*>{back_hw_JetMET, back_hw_TowerMET, back_hw_TotalMET};
            std::vector<std::string> rateLabels = hasOverlapRemoval
                ? std::vector<std::string>{"GEP Jet MET", "GEP Total MET"}
                : std::vector<std::string>{"GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
            drawRateVsThresholdMulti(rateHists, rateLabels,
                                     "Rate vs GEP MET threshold", "MET threshold [GeV]",
                                     fDir + "Rate_GEP_AlgoComparison.pdf", "",
                                     "Estimated Background Rate [Hz]");
        }

        // --- Per-file turn-on curves: all individual algorithms ---
        std::vector<std::string> allTOLabels = hasOverlapRemoval
            ? std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Total MET"}
            : std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
        auto allTOEffs = [&](TH1F* g, TH1F* gNC, TH1F* gRms, TH1F* j, TH1F* jet, TH1F* tower, TH1F* total) {
            return hasOverlapRemoval
                ? std::vector<TH1F*>{g, gNC, gRms, j, jet, total}
                : std::vector<TH1F*>{g, gNC, gRms, j, jet, tower, total};
        };
        auto allTOThrs = [&](double g, double gNC, double gRms, double j, double jet, double tower, double total) {
            return hasOverlapRemoval
                ? std::vector<double>{g, gNC, gRms, j, jet, total}
                : std::vector<double>{g, gNC, gRms, j, jet, tower, total};
        };
        drawTurnOnOverlay(
            allTOEffs(eff_gMET_20kHz, eff_gMET_NC_20kHz, eff_gMET_Rms_20kHz, eff_jMET_20kHz, eff_JetMET_20kHz, eff_TowerMET_20kHz, eff_TotalMET_20kHz),
            allTOLabels,
            "Turn-on at 20 kHz", fDir + "TurnOn_20kHz.pdf",
            allTOThrs(thr_gMET_20kHz, thr_gMET_NC_20kHz, thr_gMET_Rms_20kHz, thr_jMET_20kHz, thr_JetMET_20kHz, thr_TowerMET_20kHz, thr_TotalMET_20kHz),
            "Rate = 20 kHz", sig_h_metTruthNonInt_coarse);
        drawTurnOnOverlay(
            allTOEffs(eff_gMET_40kHz, eff_gMET_NC_40kHz, eff_gMET_Rms_40kHz, eff_jMET_40kHz, eff_JetMET_40kHz, eff_TowerMET_40kHz, eff_TotalMET_40kHz),
            allTOLabels,
            "Turn-on at 40 kHz", fDir + "TurnOn_40kHz.pdf",
            allTOThrs(thr_gMET_40kHz, thr_gMET_NC_40kHz, thr_gMET_Rms_40kHz, thr_jMET_40kHz, thr_JetMET_40kHz, thr_TowerMET_40kHz, thr_TotalMET_40kHz),
            "Rate = 40 kHz", sig_h_metTruthNonInt_coarse);
        drawTurnOnOverlay(
            allTOEffs(eff_gMET_80kHz, eff_gMET_NC_80kHz, eff_gMET_Rms_80kHz, eff_jMET_80kHz, eff_JetMET_80kHz, eff_TowerMET_80kHz, eff_TotalMET_80kHz),
            allTOLabels,
            "Turn-on at 80 kHz", fDir + "TurnOn_80kHz.pdf",
            allTOThrs(thr_gMET_80kHz, thr_gMET_NC_80kHz, thr_gMET_Rms_80kHz, thr_jMET_80kHz, thr_JetMET_80kHz, thr_TowerMET_80kHz, thr_TotalMET_80kHz),
            "Rate = 80 kHz", sig_h_metTruthNonInt_coarse);
        // gFEX-only algorithm comparison turn-ons
        drawTurnOnOverlay(
            {eff_gMET_20kHz, eff_gMET_NC_20kHz, eff_gMET_Rms_20kHz},
            {"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms"},
            "gFEX algorithm comparison — Turn-on at 20 kHz", fDir + "TurnOn_gFEX_AlgoComparison_20kHz.pdf",
            {thr_gMET_20kHz, thr_gMET_NC_20kHz, thr_gMET_Rms_20kHz}, "Rate = 20 kHz",
            sig_h_metTruthNonInt_coarse);
        drawTurnOnOverlay(
            {eff_gMET_40kHz, eff_gMET_NC_40kHz, eff_gMET_Rms_40kHz},
            {"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms"},
            "gFEX algorithm comparison — Turn-on at 40 kHz", fDir + "TurnOn_gFEX_AlgoComparison_40kHz.pdf",
            {thr_gMET_40kHz, thr_gMET_NC_40kHz, thr_gMET_Rms_40kHz}, "Rate = 40 kHz",
            sig_h_metTruthNonInt_coarse);
        drawTurnOnOverlay(
            {eff_gMET_80kHz, eff_gMET_NC_80kHz, eff_gMET_Rms_80kHz},
            {"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms"},
            "gFEX algorithm comparison — Turn-on at 80 kHz", fDir + "TurnOn_gFEX_AlgoComparison_80kHz.pdf",
            {thr_gMET_80kHz, thr_gMET_NC_80kHz, thr_gMET_Rms_80kHz}, "Rate = 80 kHz",
            sig_h_metTruthNonInt_coarse);
        // GEP-only algorithm comparison turn-ons
        std::vector<std::string> gepTOLabels = hasOverlapRemoval
            ? std::vector<std::string>{"GEP Jet MET", "GEP Total MET"}
            : std::vector<std::string>{"GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
        auto gepTOEffs = [&](TH1F* jet, TH1F* tower, TH1F* total) {
            return hasOverlapRemoval
                ? std::vector<TH1F*>{jet, total}
                : std::vector<TH1F*>{jet, tower, total};
        };
        auto gepTOThrs = [&](double jet, double tower, double total) {
            return hasOverlapRemoval
                ? std::vector<double>{jet, total}
                : std::vector<double>{jet, tower, total};
        };
        drawTurnOnOverlay(
            gepTOEffs(eff_JetMET_20kHz, eff_TowerMET_20kHz, eff_TotalMET_20kHz),
            gepTOLabels,
            "GEP algorithm comparison — Turn-on at 20 kHz", fDir + "TurnOn_GEP_AlgoComparison_20kHz.pdf",
            gepTOThrs(thr_JetMET_20kHz, thr_TowerMET_20kHz, thr_TotalMET_20kHz), "Rate = 20 kHz",
            sig_h_metTruthNonInt_coarse);
        drawTurnOnOverlay(
            gepTOEffs(eff_JetMET_40kHz, eff_TowerMET_40kHz, eff_TotalMET_40kHz),
            gepTOLabels,
            "GEP algorithm comparison — Turn-on at 40 kHz", fDir + "TurnOn_GEP_AlgoComparison_40kHz.pdf",
            gepTOThrs(thr_JetMET_40kHz, thr_TowerMET_40kHz, thr_TotalMET_40kHz), "Rate = 40 kHz",
            sig_h_metTruthNonInt_coarse);
        drawTurnOnOverlay(
            gepTOEffs(eff_JetMET_80kHz, eff_TowerMET_80kHz, eff_TotalMET_80kHz),
            gepTOLabels,
            "GEP algorithm comparison — Turn-on at 80 kHz", fDir + "TurnOn_GEP_AlgoComparison_80kHz.pdf",
            gepTOThrs(thr_JetMET_80kHz, thr_TowerMET_80kHz, thr_TotalMET_80kHz), "Rate = 80 kHz",
            sig_h_metTruthNonInt_coarse);

        // --- Combined gFEX+GEP turn-on curves at best thresholds ---
        std::string comboDir = fDir + "Combined/";
        gSystem->mkdir(comboDir.c_str(), true);
        // Helper lambda for combined threshold label
        auto comboThrLabel = [](const BestThresh2D& b) {
            return Form("gFEX>%.0f GeV & GEP>%.0f GeV", b.t1, b.t2);
        };
        // 20 kHz combined turn-ons — Jet and Tower in separate plots
        drawTurnOnOverlay(
            {eff_JetMET_20kHz,
             eff_combo_JwoJ_Jet_20kHz, eff_combo_NC_Jet_20kHz, eff_combo_Rms_Jet_20kHz},
            {"GEP Jet MET",
             std::string("JwoJ+Jet (")+comboThrLabel(best_JwoJ_Jet_20)+")",
             std::string("NC+Jet (")  +comboThrLabel(best_NC_Jet_20)+")",
             std::string("Rms+Jet (") +comboThrLabel(best_Rms_Jet_20)+")"},
            "GEP Jet MET + Combined Turn-on at 20 kHz", comboDir + "TurnOn_combined_Jet_20kHz.pdf",
            {}, "Rate = 20 kHz", sig_h_metTruthNonInt_coarse, 0.38, 0.15);
        drawTurnOnOverlay(
            {eff_TowerMET_20kHz,
             eff_combo_JwoJ_Tower_20kHz, eff_combo_NC_Tower_20kHz, eff_combo_Rms_Tower_20kHz},
            {"GEP Tower MET",
             std::string("JwoJ+Tower (")+comboThrLabel(best_JwoJ_Tower_20)+")",
             std::string("NC+Tower (")  +comboThrLabel(best_NC_Tower_20)+")",
             std::string("Rms+Tower (") +comboThrLabel(best_Rms_Tower_20)+")"},
            "GEP Tower MET + Combined Turn-on at 20 kHz", comboDir + "TurnOn_combined_Tower_20kHz.pdf",
            {}, "Rate = 20 kHz", sig_h_metTruthNonInt_coarse, 0.38, 0.15);
        // 40 kHz combined turn-ons
        drawTurnOnOverlay(
            {eff_JetMET_40kHz,
             eff_combo_JwoJ_Jet_40kHz, eff_combo_NC_Jet_40kHz, eff_combo_Rms_Jet_40kHz},
            {"GEP Jet MET",
             std::string("JwoJ+Jet (")+comboThrLabel(best_JwoJ_Jet_40)+")",
             std::string("NC+Jet (")  +comboThrLabel(best_NC_Jet_40)+")",
             std::string("Rms+Jet (") +comboThrLabel(best_Rms_Jet_40)+")"},
            "GEP Jet MET + Combined Turn-on at 40 kHz", comboDir + "TurnOn_combined_Jet_40kHz.pdf",
            {}, "Rate = 40 kHz", sig_h_metTruthNonInt_coarse, 0.38, 0.15);
        drawTurnOnOverlay(
            {eff_TowerMET_40kHz,
             eff_combo_JwoJ_Tower_40kHz, eff_combo_NC_Tower_40kHz, eff_combo_Rms_Tower_40kHz},
            {"GEP Tower MET",
             std::string("JwoJ+Tower (")+comboThrLabel(best_JwoJ_Tower_40)+")",
             std::string("NC+Tower (")  +comboThrLabel(best_NC_Tower_40)+")",
             std::string("Rms+Tower (") +comboThrLabel(best_Rms_Tower_40)+")"},
            "GEP Tower MET + Combined Turn-on at 40 kHz", comboDir + "TurnOn_combined_Tower_40kHz.pdf",
            {}, "Rate = 40 kHz", sig_h_metTruthNonInt_coarse, 0.38, 0.15);
        // 80 kHz combined turn-ons
        drawTurnOnOverlay(
            {eff_JetMET_80kHz,
             eff_combo_JwoJ_Jet_80kHz, eff_combo_NC_Jet_80kHz, eff_combo_Rms_Jet_80kHz},
            {"GEP Jet MET",
             std::string("JwoJ+Jet (")+comboThrLabel(best_JwoJ_Jet_80)+")",
             std::string("NC+Jet (")  +comboThrLabel(best_NC_Jet_80)+")",
             std::string("Rms+Jet (") +comboThrLabel(best_Rms_Jet_80)+")"},
            "GEP Jet MET + Combined Turn-on at 80 kHz", comboDir + "TurnOn_combined_Jet_80kHz.pdf",
            {}, "Rate = 80 kHz", sig_h_metTruthNonInt_coarse, 0.38, 0.15);
        drawTurnOnOverlay(
            {eff_TowerMET_80kHz,
             eff_combo_JwoJ_Tower_80kHz, eff_combo_NC_Tower_80kHz, eff_combo_Rms_Tower_80kHz},
            {"GEP Tower MET",
             std::string("JwoJ+Tower (")+comboThrLabel(best_JwoJ_Tower_80)+")",
             std::string("NC+Tower (")  +comboThrLabel(best_NC_Tower_80)+")",
             std::string("Rms+Tower (") +comboThrLabel(best_Rms_Tower_80)+")"},
            "GEP Tower MET + Combined Turn-on at 80 kHz", comboDir + "TurnOn_combined_Tower_80kHz.pdf",
            {}, "Rate = 80 kHz", sig_h_metTruthNonInt_coarse, 0.38, 0.15);

        // --- Rate vs Efficiency curves (background rate vs signal efficiency) ---
        auto out_RateVsEff_gMET     = MakeRateVsEff(sig_h_gMET,     back_hw_gMET);
        auto out_RateVsEff_gMET_NC  = MakeRateVsEff(sig_h_gMET_NC,  back_hw_gMET_NC);
        auto out_RateVsEff_gMET_Rms = MakeRateVsEff(sig_h_gMET_Rms, back_hw_gMET_Rms);
        auto out_RateVsEff_jMET     = MakeRateVsEff(sig_h_jMET,     back_hw_jMET);
        auto out_RateVsEff_JetMET   = MakeRateVsEff(sig_h_JetMet,   back_hw_JetMET);
        auto out_RateVsEff_TowerMET = MakeRateVsEff(sig_h_TowerMet, back_hw_TowerMET);
        auto out_RateVsEff_TotalMET = MakeRateVsEff(sig_h_TotalMET, back_hw_TotalMET);

        // Style helper
        auto styleRVE = [](TGraph* g, Color_t col) {
            g->SetMarkerColor(col); g->SetLineColor(col);
            g->SetMarkerStyle(20);  g->SetMarkerSize(0.8); g->SetLineWidth(2);
            g->GetXaxis()->SetTitle("Signal Efficiency");
            g->GetYaxis()->SetTitle("Estimated Background Rate [Hz]");
        };
        styleRVE(out_RateVsEff_gMET.gRate_vsEff,     kBlack);
        styleRVE(out_RateVsEff_gMET_NC.gRate_vsEff,  kRed);
        styleRVE(out_RateVsEff_gMET_Rms.gRate_vsEff, kOrange+1);
        styleRVE(out_RateVsEff_jMET.gRate_vsEff,     kCyan+2);
        styleRVE(out_RateVsEff_JetMET.gRate_vsEff,   kBlue);
        styleRVE(out_RateVsEff_TowerMET.gRate_vsEff, kGreen+2);
        styleRVE(out_RateVsEff_TotalMET.gRate_vsEff, kMagenta+1);

        // Individual PDFs
        auto drawRVEsingle = [&](TGraph* g, const std::string& title, const std::string& path) {
            TCanvas c("c", "", 700, 600);
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
            gPad->SetLogy();
            g->SetTitle((title+";Signal Efficiency;Estimated Background Rate [Hz]").c_str());
            g->Draw("AP");
            c.SaveAs(path.c_str());
        };
        drawRVEsingle(out_RateVsEff_gMET.gRate_vsEff,     "gFEX MET (JwoJ)",    fDir + "RateVsEff_gFEX_MET_JwoJ.pdf");
        drawRVEsingle(out_RateVsEff_gMET_NC.gRate_vsEff,  "gFEX MET (NoiseCut)",fDir + "RateVsEff_gFEX_MET_NoiseCut.pdf");
        drawRVEsingle(out_RateVsEff_gMET_Rms.gRate_vsEff, "gFEX MET (Rms)",     fDir + "RateVsEff_gFEX_MET_Rms.pdf");
        drawRVEsingle(out_RateVsEff_jMET.gRate_vsEff,     "jFEX MET",           fDir + "RateVsEff_jFEX_MET.pdf");
        drawRVEsingle(out_RateVsEff_JetMET.gRate_vsEff,   "GEP Jet MET",        fDir + "RateVsEff_JetMET.pdf");
        drawRVEsingle(out_RateVsEff_TowerMET.gRate_vsEff, "GEP Tower MET",      fDir + "RateVsEff_TowerMET.pdf");
        drawRVEsingle(out_RateVsEff_TotalMET.gRate_vsEff, "GEP Total MET",      fDir + "RateVsEff_TotalMET.pdf");

        // Overlay: all individual algorithms — gFEX differentiated by JwoJ/NC/Rms, jFEX, and GEP types
        {
            std::vector<TGraph*> allRVE = hasOverlapRemoval
                ? std::vector<TGraph*>{out_RateVsEff_gMET.gRate_vsEff, out_RateVsEff_gMET_NC.gRate_vsEff, out_RateVsEff_gMET_Rms.gRate_vsEff, out_RateVsEff_jMET.gRate_vsEff,
                                       out_RateVsEff_JetMET.gRate_vsEff, out_RateVsEff_TotalMET.gRate_vsEff}
                : std::vector<TGraph*>{out_RateVsEff_gMET.gRate_vsEff, out_RateVsEff_gMET_NC.gRate_vsEff, out_RateVsEff_gMET_Rms.gRate_vsEff, out_RateVsEff_jMET.gRate_vsEff,
                                       out_RateVsEff_JetMET.gRate_vsEff, out_RateVsEff_TowerMET.gRate_vsEff, out_RateVsEff_TotalMET.gRate_vsEff};
            std::vector<std::string> allRVELabels = hasOverlapRemoval
                ? std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Total MET"}
                : std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
            drawRateVsEffOverlay(allRVE, allRVELabels, (fDir + "RateVsEff_overlay.pdf"), signalName);
        }
        // gFEX-only algorithm comparison
        drawRateVsEffOverlay(
            {out_RateVsEff_gMET.gRate_vsEff, out_RateVsEff_gMET_NC.gRate_vsEff, out_RateVsEff_gMET_Rms.gRate_vsEff},
            {"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms"},
            (fDir + "RateVsEff_gFEX_AlgoComparison.pdf"), signalName);
        // GEP-only algorithm comparison
        {
            std::vector<TGraph*> gepRVE = hasOverlapRemoval
                ? std::vector<TGraph*>{out_RateVsEff_JetMET.gRate_vsEff, out_RateVsEff_TotalMET.gRate_vsEff}
                : std::vector<TGraph*>{out_RateVsEff_JetMET.gRate_vsEff, out_RateVsEff_TowerMET.gRate_vsEff, out_RateVsEff_TotalMET.gRate_vsEff};
            std::vector<std::string> gepRVELabels = hasOverlapRemoval
                ? std::vector<std::string>{"GEP Jet MET", "GEP Total MET"}
                : std::vector<std::string>{"GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
            drawRateVsEffOverlay(gepRVE, gepRVELabels, (fDir + "RateVsEff_GEP_AlgoComparison.pdf"), signalName);
        }

        // --- Combined gFEX+GEP rate-vs-eff frontiers ---
        TGraph* frontier_JwoJ_Jet   = extractFrontier2D(out2D_JwoJ_Jet);
        TGraph* frontier_JwoJ_Tower = extractFrontier2D(out2D_JwoJ_Tower);
        TGraph* frontier_NC_Jet     = extractFrontier2D(out2D_NC_Jet);
        TGraph* frontier_NC_Tower   = extractFrontier2D(out2D_NC_Tower);
        TGraph* frontier_Rms_Jet    = extractFrontier2D(out2D_Rms_Jet);
        TGraph* frontier_Rms_Tower  = extractFrontier2D(out2D_Rms_Tower);
        drawRateVsEffOverlay(
            {frontier_JwoJ_Jet, frontier_JwoJ_Tower, frontier_NC_Jet,
             frontier_NC_Tower, frontier_Rms_Jet, frontier_Rms_Tower},
            {"JwoJ + GEP Jet", "JwoJ + GEP Tower", "NoiseCut + GEP Jet",
             "NoiseCut + GEP Tower", "Rms + GEP Jet", "Rms + GEP Tower"},
            (comboDir + "RateVsEff_combined_frontiers.pdf"), signalName);
        // Combined vs individual overlay (best combined per gFEX type vs each individual)
        {
            std::vector<TGraph*> cviRVE = hasOverlapRemoval
                ? std::vector<TGraph*>{out_RateVsEff_gMET.gRate_vsEff,   out_RateVsEff_gMET_NC.gRate_vsEff, out_RateVsEff_gMET_Rms.gRate_vsEff,
                                       out_RateVsEff_JetMET.gRate_vsEff, out_RateVsEff_TotalMET.gRate_vsEff,
                                       frontier_JwoJ_Jet, frontier_NC_Jet, frontier_Rms_Jet}
                : std::vector<TGraph*>{out_RateVsEff_gMET.gRate_vsEff,   out_RateVsEff_gMET_NC.gRate_vsEff, out_RateVsEff_gMET_Rms.gRate_vsEff,
                                       out_RateVsEff_JetMET.gRate_vsEff, out_RateVsEff_TowerMET.gRate_vsEff, out_RateVsEff_TotalMET.gRate_vsEff,
                                       frontier_JwoJ_Jet, frontier_NC_Jet, frontier_Rms_Jet};
            std::vector<std::string> cviLabels = hasOverlapRemoval
                ? std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms",
                                           "GEP Jet MET", "GEP Total MET",
                                           "JwoJ + GEP Jet (combined)", "NC + GEP Jet (combined)", "Rms + GEP Jet (combined)"}
                : std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms",
                                           "GEP Jet MET", "GEP Tower MET", "GEP Total MET",
                                           "JwoJ + GEP Jet (combined)", "NC + GEP Jet (combined)", "Rms + GEP Jet (combined)"};
            drawRateVsEffOverlay(cviRVE, cviLabels, (comboDir + "RateVsEff_combined_vs_individual.pdf"), signalName);
        }

        // --- AOD vs Sim gFEX MET comparison (only for resim ntuples) ---
        if (hasGFexSimMET) {
            // Distribution overlays: AOD vs Sim
            drawAlgoComparison(sig_h_gMET, back_h_gMET, sig_h_gMET_JwoJSim, back_h_gMET_JwoJSim,
                               "gFEX JwoJ (AOD)", "gFEX JwoJ (Sim)",
                               "MET [GeV]", fDir + "gFEX_JwoJ_AOD_vs_Sim.pdf", signalName);
            drawAlgoComparison(sig_h_gMET_NC, back_h_gMET_NC, sig_h_gMET_NCSim, back_h_gMET_NCSim,
                               "gFEX NoiseCut (AOD)", "gFEX NoiseCut (Sim)",
                               "MET [GeV]", fDir + "gFEX_NC_AOD_vs_Sim.pdf", signalName);
            drawAlgoComparison(sig_h_gMET_Rms, back_h_gMET_Rms, sig_h_gMET_RmsSim, back_h_gMET_RmsSim,
                               "gFEX Rms (AOD)", "gFEX Rms (Sim)",
                               "MET [GeV]", fDir + "gFEX_Rms_AOD_vs_Sim.pdf", signalName);
            // Rate vs threshold: AOD vs Sim
            drawRateVsThresholdMulti({back_hw_gMET,     back_hw_gMET_JwoJSim},
                                     {"gFEX JwoJ (AOD)", "gFEX JwoJ (Sim)"},
                                     "Rate vs gFEX JwoJ MET (AOD vs Sim)", "MET threshold [GeV]",
                                     fDir + "Rate_gFEX_JwoJ_AOD_vs_Sim.pdf", "", "Estimated Background Rate [Hz]");
            drawRateVsThresholdMulti({back_hw_gMET_NC,   back_hw_gMET_NCSim},
                                     {"gFEX NoiseCut (AOD)", "gFEX NoiseCut (Sim)"},
                                     "Rate vs gFEX NoiseCut MET (AOD vs Sim)", "MET threshold [GeV]",
                                     fDir + "Rate_gFEX_NC_AOD_vs_Sim.pdf", "", "Estimated Background Rate [Hz]");
            drawRateVsThresholdMulti({back_hw_gMET_Rms,  back_hw_gMET_RmsSim},
                                     {"gFEX Rms (AOD)", "gFEX Rms (Sim)"},
                                     "Rate vs gFEX Rms MET (AOD vs Sim)", "MET threshold [GeV]",
                                     fDir + "Rate_gFEX_Rms_AOD_vs_Sim.pdf", "", "Estimated Background Rate [Hz]");
            // Turn-on curves: AOD vs Sim at 20/40/80 kHz
            drawTurnOnOverlay(
                {eff_gMET_80kHz, eff_gMET_JwoJSim_80kHz},
                {"gFEX JwoJ (AOD)", "gFEX JwoJ (Sim)"},
                "gFEX JwoJ Turn-on at 80 kHz (AOD vs Sim)", fDir + "TurnOn_gFEX_JwoJ_AOD_vs_Sim_80kHz.pdf",
                {thr_gMET_80kHz, thr_gMET_JwoJSim_80kHz}, "Rate = 80 kHz", sig_h_metTruthNonInt_coarse);
            drawTurnOnOverlay(
                {eff_gMET_NC_80kHz, eff_gMET_NCSim_80kHz},
                {"gFEX NoiseCut (AOD)", "gFEX NoiseCut (Sim)"},
                "gFEX NoiseCut Turn-on at 80 kHz (AOD vs Sim)", fDir + "TurnOn_gFEX_NC_AOD_vs_Sim_80kHz.pdf",
                {thr_gMET_NC_80kHz, thr_gMET_NCSim_80kHz}, "Rate = 80 kHz", sig_h_metTruthNonInt_coarse);
            drawTurnOnOverlay(
                {eff_gMET_Rms_80kHz, eff_gMET_RmsSim_80kHz},
                {"gFEX Rms (AOD)", "gFEX Rms (Sim)"},
                "gFEX Rms Turn-on at 80 kHz (AOD vs Sim)", fDir + "TurnOn_gFEX_Rms_AOD_vs_Sim_80kHz.pdf",
                {thr_gMET_Rms_80kHz, thr_gMET_RmsSim_80kHz}, "Rate = 80 kHz", sig_h_metTruthNonInt_coarse);
            // Rate vs efficiency: AOD vs Sim
            auto out_RVE_JwoJSim = MakeRateVsEff(sig_h_gMET_JwoJSim, back_hw_gMET_JwoJSim);
            auto out_RVE_NCSim   = MakeRateVsEff(sig_h_gMET_NCSim,   back_hw_gMET_NCSim);
            auto out_RVE_RmsSim  = MakeRateVsEff(sig_h_gMET_RmsSim,  back_hw_gMET_RmsSim);
            styleRVE(out_RVE_JwoJSim.gRate_vsEff, kBlack);
            styleRVE(out_RVE_NCSim.gRate_vsEff,   kRed);
            styleRVE(out_RVE_RmsSim.gRate_vsEff,  kOrange+1);
            drawRateVsEffOverlay(
                {out_RateVsEff_gMET.gRate_vsEff,    out_RVE_JwoJSim.gRate_vsEff},
                {"gFEX JwoJ (AOD)", "gFEX JwoJ (Sim)"},
                (fDir + "RateVsEff_gFEX_JwoJ_AOD_vs_Sim.pdf"), signalName);
            drawRateVsEffOverlay(
                {out_RateVsEff_gMET_NC.gRate_vsEff, out_RVE_NCSim.gRate_vsEff},
                {"gFEX NoiseCut (AOD)", "gFEX NoiseCut (Sim)"},
                (fDir + "RateVsEff_gFEX_NC_AOD_vs_Sim.pdf"), signalName);
            drawRateVsEffOverlay(
                {out_RateVsEff_gMET_Rms.gRate_vsEff, out_RVE_RmsSim.gRate_vsEff},
                {"gFEX Rms (AOD)", "gFEX Rms (Sim)"},
                (fDir + "RateVsEff_gFEX_Rms_AOD_vs_Sim.pdf"), signalName);
            drawRateVsEffOverlay(
                {out_RateVsEff_gMET.gRate_vsEff,    out_RateVsEff_gMET_NC.gRate_vsEff,  out_RateVsEff_gMET_Rms.gRate_vsEff,
                 out_RVE_JwoJSim.gRate_vsEff,        out_RVE_NCSim.gRate_vsEff,           out_RVE_RmsSim.gRate_vsEff},
                {"gFEX JwoJ (AOD)", "gFEX NoiseCut (AOD)", "gFEX Rms (AOD)",
                 "gFEX JwoJ (Sim)", "gFEX NoiseCut (Sim)", "gFEX Rms (Sim)"},
                (fDir + "RateVsEff_gFEX_AOD_vs_Sim_all.pdf"), signalName);
        }

        // Store clones for multi-file overlay — SetDirectory(0) detaches from
        // the file so ROOT doesn't delete them when the file is closed.
        auto cloneDetached = [](TH1F* h) -> TH1F* {
            TH1F* c = (TH1F*)h->Clone();
            c->SetDirectory(0);
            return c;
        };
        sig_h_TotalMET_vec.push_back(cloneDetached(sig_h_TotalMET));
        back_h_TotalMET_vec.push_back(cloneDetached(back_h_TotalMET));
        sig_h_TotalMETX_vec.push_back(cloneDetached(sig_h_TotalMETX));
        back_h_TotalMETX_vec.push_back(cloneDetached(back_h_TotalMETX));
        sig_h_TotalMETY_vec.push_back(cloneDetached(sig_h_TotalMETY));
        back_h_TotalMETY_vec.push_back(cloneDetached(back_h_TotalMETY));
        sig_h_TowerMet_vec.push_back(cloneDetached(sig_h_TowerMet));
        back_h_TowerMet_vec.push_back(cloneDetached(back_h_TowerMet));
        sig_h_JetMet_vec.push_back(cloneDetached(sig_h_JetMet));
        back_h_JetMet_vec.push_back(cloneDetached(back_h_JetMet));
        sig_h_SumET_vec.push_back(cloneDetached(sig_h_SumET));
        back_h_SumET_vec.push_back(cloneDetached(back_h_SumET));
        sig_h_SumJetET_vec.push_back(cloneDetached(sig_h_SumJetET));
        back_h_SumJetET_vec.push_back(cloneDetached(back_h_SumJetET));
        sig_h_SumTowerET_vec.push_back(cloneDetached(sig_h_SumTowerET));
        back_h_SumTowerET_vec.push_back(cloneDetached(back_h_SumTowerET));
        sig_h_gMET_vec.push_back(cloneDetached(sig_h_gMET));
        back_h_gMET_vec.push_back(cloneDetached(back_h_gMET));
        sig_h_gMET_NC_vec.push_back(cloneDetached(sig_h_gMET_NC));
        back_h_gMET_NC_vec.push_back(cloneDetached(back_h_gMET_NC));
        sig_h_gMET_Rms_vec.push_back(cloneDetached(sig_h_gMET_Rms));
        back_h_gMET_Rms_vec.push_back(cloneDetached(back_h_gMET_Rms));
        sig_h_jMET_vec.push_back(cloneDetached(sig_h_jMET));
        back_h_jMET_vec.push_back(cloneDetached(back_h_jMET));
        sig_h_metTruthNonInt_vec.push_back(cloneDetached(sig_h_metTruthNonInt));
        back_h_metTruthNonInt_vec.push_back(cloneDetached(back_h_metTruthNonInt));
        back_hw_TotalMET_vec.push_back(cloneDetached(back_hw_TotalMET));
        back_hw_gMET_vec.push_back(cloneDetached(back_hw_gMET));
        back_hw_gMET_NC_vec.push_back(cloneDetached(back_hw_gMET_NC));
        back_hw_gMET_Rms_vec.push_back(cloneDetached(back_hw_gMET_Rms));
        back_hw_jMET_vec.push_back(cloneDetached(back_hw_jMET));
        back_hw_JetMET_vec.push_back(cloneDetached(back_hw_JetMET));
        back_hw_TowerMET_vec.push_back(cloneDetached(back_hw_TowerMET));
        eff_gMET_80kHz_vec.push_back(cloneDetached(eff_gMET_80kHz));
        eff_gMET_NC_80kHz_vec.push_back(cloneDetached(eff_gMET_NC_80kHz));
        eff_gMET_Rms_80kHz_vec.push_back(cloneDetached(eff_gMET_Rms_80kHz));
        eff_jMET_80kHz_vec.push_back(cloneDetached(eff_jMET_80kHz));
        eff_JetMET_80kHz_vec.push_back(cloneDetached(eff_JetMET_80kHz));
        eff_TowerMET_80kHz_vec.push_back(cloneDetached(eff_TowerMET_80kHz));
        eff_TotalMET_80kHz_vec.push_back(cloneDetached(eff_TotalMET_80kHz));
        thr_gMET_80kHz_vec.push_back(thr_gMET_80kHz);
        thr_gMET_NC_80kHz_vec.push_back(thr_gMET_NC_80kHz);
        thr_gMET_Rms_80kHz_vec.push_back(thr_gMET_Rms_80kHz);
        thr_jMET_80kHz_vec.push_back(thr_jMET_80kHz);
        thr_JetMET_80kHz_vec.push_back(thr_JetMET_80kHz);
        thr_TowerMET_80kHz_vec.push_back(thr_TowerMET_80kHz);
        thr_TotalMET_80kHz_vec.push_back(thr_TotalMET_80kHz);
        sig_h_metTruthNonInt_unscaled_vec.push_back(cloneDetached(sig_h_metTruthNonInt_coarse));
        if (hasGFexSimMET) {
            sig_h_gMET_JwoJSim_vec.push_back(cloneDetached(sig_h_gMET_JwoJSim));
            back_h_gMET_JwoJSim_vec.push_back(cloneDetached(back_h_gMET_JwoJSim));
            sig_h_gMET_NCSim_vec.push_back(cloneDetached(sig_h_gMET_NCSim));
            back_h_gMET_NCSim_vec.push_back(cloneDetached(back_h_gMET_NCSim));
            sig_h_gMET_RmsSim_vec.push_back(cloneDetached(sig_h_gMET_RmsSim));
            back_h_gMET_RmsSim_vec.push_back(cloneDetached(back_h_gMET_RmsSim));
            back_hw_gMET_JwoJSim_vec.push_back(cloneDetached(back_hw_gMET_JwoJSim));
            back_hw_gMET_NCSim_vec.push_back(cloneDetached(back_hw_gMET_NCSim));
            back_hw_gMET_RmsSim_vec.push_back(cloneDetached(back_hw_gMET_RmsSim));
            eff_gMET_JwoJSim_80kHz_vec.push_back(cloneDetached(eff_gMET_JwoJSim_80kHz));
            eff_gMET_NCSim_80kHz_vec.push_back(cloneDetached(eff_gMET_NCSim_80kHz));
            eff_gMET_RmsSim_80kHz_vec.push_back(cloneDetached(eff_gMET_RmsSim_80kHz));
            thr_gMET_JwoJSim_80kHz_vec.push_back(thr_gMET_JwoJSim_80kHz);
            thr_gMET_NCSim_80kHz_vec.push_back(thr_gMET_NCSim_80kHz);
            thr_gMET_RmsSim_80kHz_vec.push_back(thr_gMET_RmsSim_80kHz);
        }
        sigF->Close();
        backF->Close();
    } // file loop

    // --- Multi-file overlays (sig + bkg overlaid, dashed=bkg) ---
    if (signalFiles.size() > 1) {
        std::string mDir = "multiFileOverlay_MET/";
        gSystem->mkdir(mDir.c_str(), true);
        drawOverlayMulti(sig_h_TotalMET_vec,      back_h_TotalMET_vec,      labels, "Total MET (GEP)",     "MET [GeV]",          mDir + "TotalMET.pdf",        signalName);
        drawOverlayMulti(sig_h_TowerMet_vec,      back_h_TowerMet_vec,      labels, "Tower MET (GEP)",     "MET [GeV]",          mDir + "TowerMET.pdf",        signalName);
        drawOverlayMulti(sig_h_JetMet_vec,        back_h_JetMet_vec,        labels, "Jet MET (GEP)",       "MET [GeV]",          mDir + "JetMET.pdf",          signalName);
        drawOverlayMulti(sig_h_SumET_vec,         back_h_SumET_vec,         labels, "GEP TOB #Sigma E_{T}",         "GEP TOB #Sigma E_{T} [GeV]",   mDir + "SumET.pdf",      signalName);
        if (hasSumJetET)   drawOverlayMulti(sig_h_SumJetET_vec,   back_h_SumJetET_vec,   labels, "GEP H_{T} (Sum Jet E_{T})",    "GEP H_{T} [GeV]",              mDir + "SumJetET.pdf",   signalName);
        if (hasSumTowerET) drawOverlayMulti(sig_h_SumTowerET_vec, back_h_SumTowerET_vec, labels, "GEP Tower #Sigma E_{T}",       "GEP Tower #Sigma E_{T} [GeV]", mDir + "SumTowerET.pdf", signalName);
        drawOverlayMulti(sig_h_gMET_vec,          back_h_gMET_vec,          labels, "gFEX MET (JwoJ)",     "MET [GeV]",          mDir + "gFEX_MET_JwoJ.pdf",     signalName);
        drawOverlayMulti(sig_h_gMET_NC_vec,       back_h_gMET_NC_vec,       labels, "gFEX MET (NoiseCut)", "MET [GeV]",          mDir + "gFEX_MET_NoiseCut.pdf", signalName);
        drawOverlayMulti(sig_h_gMET_Rms_vec,      back_h_gMET_Rms_vec,      labels, "gFEX MET (Rms)",      "MET [GeV]",          mDir + "gFEX_MET_Rms.pdf",      signalName);
        drawOverlayMulti(sig_h_jMET_vec,          back_h_jMET_vec,          labels, "jFEX MET",            "MET [GeV]",          mDir + "jFEX_MET.pdf",          signalName);
        drawOverlayMulti(sig_h_metTruthNonInt_vec,back_h_metTruthNonInt_vec,labels, "Truth MET (NonInt)",  "MET [GeV]",          mDir + "TruthMET_NonInt.pdf",   signalName);

        drawRateVsThresholdMulti(back_hw_TotalMET_vec,  labels, "Rate vs Emulated MET threshold",        "MET threshold [GeV]", mDir + "Rate_TotalMET.pdf",        signalName);
        drawRateVsThresholdMulti(back_hw_gMET_vec,      labels, "Rate vs gFEX MET threshold (JwoJ)",     "MET threshold [GeV]", mDir + "Rate_gFEX_MET_JwoJ.pdf",   signalName);
        drawRateVsThresholdMulti(back_hw_gMET_NC_vec,   labels, "Rate vs gFEX MET threshold (NoiseCut)", "MET threshold [GeV]", mDir + "Rate_gFEX_MET_NoiseCut.pdf",signalName);
        drawRateVsThresholdMulti(back_hw_gMET_Rms_vec,  labels, "Rate vs gFEX MET threshold (Rms)",      "MET threshold [GeV]", mDir + "Rate_gFEX_MET_Rms.pdf",    signalName);
        drawRateVsThresholdMulti(back_hw_jMET_vec,      labels, "Rate vs jFEX MET threshold",            "MET threshold [GeV]", mDir + "Rate_jFEX_MET.pdf",        signalName);
        drawRateVsThresholdMulti(back_hw_JetMET_vec,    labels, "Rate vs GEP Jet MET threshold",         "MET threshold [GeV]", mDir + "Rate_JetMET.pdf",           signalName);
        drawRateVsThresholdMulti(back_hw_TowerMET_vec,  labels, "Rate vs GEP Tower MET threshold",       "MET threshold [GeV]", mDir + "Rate_TowerMET.pdf",         signalName);

        // --- Multi-file signal efficiency vs threshold ---
        drawEffVsThresholdMulti(sig_h_TotalMET_vec,  labels, "Signal Efficiency vs GEP Total MET Threshold",    "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_TotalMET.pdf",   signalName);
        drawEffVsThresholdMulti(sig_h_gMET_vec,      labels, "Signal Efficiency vs gFEX MET Threshold (JwoJ)",  "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_gFEX_JwoJ.pdf", signalName);
        drawEffVsThresholdMulti(sig_h_gMET_NC_vec,   labels, "Signal Efficiency vs gFEX MET Threshold (NC)",    "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_gFEX_NC.pdf",   signalName);
        drawEffVsThresholdMulti(sig_h_gMET_Rms_vec,  labels, "Signal Efficiency vs gFEX MET Threshold (Rms)",   "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_gFEX_Rms.pdf",  signalName);
        drawEffVsThresholdMulti(sig_h_jMET_vec,      labels, "Signal Efficiency vs jFEX MET Threshold",         "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_jFEX_MET.pdf",  signalName);
        drawEffVsThresholdMulti(sig_h_JetMet_vec,    labels, "Signal Efficiency vs GEP Jet MET Threshold",      "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_JetMET.pdf",    signalName);
        drawEffVsThresholdMulti(sig_h_TowerMet_vec,  labels, "Signal Efficiency vs GEP Tower MET Threshold",    "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_TowerMET.pdf",  signalName);

        // --- Rate vs Efficiency multi-file overlays ---
        {
            std::vector<TGraph*> rveJet, rveTower, rveTotal;
            for (unsigned int i = 0; i < labels.size(); i++) {
                auto outJ = MakeRateVsEff(sig_h_JetMet_vec[i],   back_hw_JetMET_vec[i]);
                auto outT = MakeRateVsEff(sig_h_TowerMet_vec[i], back_hw_TowerMET_vec[i]);
                auto outTot = MakeRateVsEff(sig_h_TotalMET_vec[i], back_hw_TotalMET_vec[i]);
                rveJet.push_back((TGraph*)outJ.gRate_vsEff->Clone(Form("rve_JetMET_%u",   i)));
                rveTower.push_back((TGraph*)outT.gRate_vsEff->Clone(Form("rve_TowerMET_%u", i)));
                rveTotal.push_back((TGraph*)outTot.gRate_vsEff->Clone(Form("rve_TotalMET_%u", i)));
            }
            drawRateVsEffOverlay(rveJet,   labels, mDir + "RateVsEff_JetMET_multi.pdf",   signalName);
            drawRateVsEffOverlay(rveTower, labels, mDir + "RateVsEff_TowerMET_multi.pdf",  signalName);
            drawRateVsEffOverlay(rveTotal, labels, mDir + "RateVsEff_TotalMET_multi.pdf",  signalName);
            for (auto* g : rveJet)   delete g;
            for (auto* g : rveTower) delete g;
            for (auto* g : rveTotal) delete g;
        }

        // --- GEP vs gFEX comparison, one canvas per config ---
        for (unsigned int i = 0; i < labels.size(); i++) {
            drawAlgoComparison(sig_h_TotalMET_vec[i], back_h_TotalMET_vec[i],
                               sig_h_gMET_vec[i],    back_h_gMET_vec[i],
                               "GEP Total MET", "gFEX MET",
                               "MET [GeV]",
                               mDir + "GEP_vs_gFEX_MET_" + labels[i] + ".pdf",
                               signalName);
        }
        // --- 80 kHz turn-on comparison across configs (one canvas per algorithm) ---
        auto makeConfigLabels = [&](const std::string& algoName) {
            std::vector<std::string> v;
            for (const auto& lbl : labels) v.push_back(algoName + " (" + lbl + ")");
            return v;
        };
        TH1F* truthOverlay = sig_h_metTruthNonInt_unscaled_vec.empty() ? nullptr : sig_h_metTruthNonInt_unscaled_vec[0];
        drawTurnOnOverlay(eff_gMET_80kHz_vec,
                          makeConfigLabels("gFEX JwoJ"),
                          "gFEX JwoJ Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_gFEX_JwoJ_multi.pdf",
                          thr_gMET_80kHz_vec, "Rate = 80 kHz", truthOverlay);
        drawTurnOnOverlay(eff_gMET_NC_80kHz_vec,
                          makeConfigLabels("gFEX NoiseCut"),
                          "gFEX NoiseCut Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_gFEX_NoiseCut_multi.pdf",
                          thr_gMET_NC_80kHz_vec, "Rate = 80 kHz", truthOverlay);
        drawTurnOnOverlay(eff_gMET_Rms_80kHz_vec,
                          makeConfigLabels("gFEX Rms"),
                          "gFEX Rms Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_gFEX_Rms_multi.pdf",
                          thr_gMET_Rms_80kHz_vec, "Rate = 80 kHz", truthOverlay);
        drawTurnOnOverlay(eff_jMET_80kHz_vec,
                          makeConfigLabels("jFEX"),
                          "jFEX MET Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_jFEX_MET_multi.pdf",
                          thr_jMET_80kHz_vec, "Rate = 80 kHz", truthOverlay);
        drawTurnOnOverlay(eff_JetMET_80kHz_vec,
                          makeConfigLabels("GEP Jet MET"),
                          "GEP Jet MET Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_JetMET_multi.pdf",
                          thr_JetMET_80kHz_vec, "Rate = 80 kHz", truthOverlay);
        drawTurnOnOverlay(eff_TowerMET_80kHz_vec,
                          makeConfigLabels("GEP Tower MET"),
                          "GEP Tower MET Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_TowerMET_multi.pdf",
                          thr_TowerMET_80kHz_vec, "Rate = 80 kHz", truthOverlay,
                          0.45, 0.20);
        drawTurnOnOverlay(eff_TotalMET_80kHz_vec,
                          makeConfigLabels("GEP Total MET"),
                          "GEP Total MET Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_TotalMET_multi.pdf",
                          thr_TotalMET_80kHz_vec, "Rate = 80 kHz", truthOverlay);

        if (!sig_h_gMET_JwoJSim_vec.empty()) {
            drawOverlayMulti(sig_h_gMET_JwoJSim_vec, back_h_gMET_JwoJSim_vec, labels, "gFEX MET (JwoJ Sim)",    "MET [GeV]", mDir + "gFEX_MET_JwoJSim.pdf",   signalName);
            drawOverlayMulti(sig_h_gMET_NCSim_vec,   back_h_gMET_NCSim_vec,   labels, "gFEX MET (NoiseCut Sim)","MET [GeV]", mDir + "gFEX_MET_NCSim.pdf",     signalName);
            drawOverlayMulti(sig_h_gMET_RmsSim_vec,  back_h_gMET_RmsSim_vec,  labels, "gFEX MET (Rms Sim)",     "MET [GeV]", mDir + "gFEX_MET_RmsSim.pdf",    signalName);
            drawRateVsThresholdMulti(back_hw_gMET_JwoJSim_vec, labels, "Rate vs gFEX MET (JwoJ Sim)",    "MET threshold [GeV]", mDir + "Rate_gFEX_MET_JwoJSim.pdf", signalName);
            drawRateVsThresholdMulti(back_hw_gMET_NCSim_vec,   labels, "Rate vs gFEX MET (NoiseCut Sim)","MET threshold [GeV]", mDir + "Rate_gFEX_MET_NCSim.pdf",   signalName);
            drawRateVsThresholdMulti(back_hw_gMET_RmsSim_vec,  labels, "Rate vs gFEX MET (Rms Sim)",     "MET threshold [GeV]", mDir + "Rate_gFEX_MET_RmsSim.pdf",  signalName);
            drawTurnOnOverlay(eff_gMET_JwoJSim_80kHz_vec, makeConfigLabels("gFEX JwoJ Sim"),
                              "gFEX JwoJ Sim Turn-on at 80 kHz",
                              mDir + "TurnOn_80kHz_gFEX_JwoJSim_multi.pdf",
                              thr_gMET_JwoJSim_80kHz_vec, "Rate = 80 kHz", truthOverlay);
            drawTurnOnOverlay(eff_gMET_NCSim_80kHz_vec, makeConfigLabels("gFEX NoiseCut Sim"),
                              "gFEX NoiseCut Sim Turn-on at 80 kHz",
                              mDir + "TurnOn_80kHz_gFEX_NCSim_multi.pdf",
                              thr_gMET_NCSim_80kHz_vec, "Rate = 80 kHz", truthOverlay);
            drawTurnOnOverlay(eff_gMET_RmsSim_80kHz_vec, makeConfigLabels("gFEX Rms Sim"),
                              "gFEX Rms Sim Turn-on at 80 kHz",
                              mDir + "TurnOn_80kHz_gFEX_RmsSim_multi.pdf",
                              thr_gMET_RmsSim_80kHz_vec, "Rate = 80 kHz", truthOverlay);
        }
    }

    std::cout << "Done. Plots saved to " << outputDir << "\n";
}

// -----------------------------------------------------------------------
void metAnalysisAndRates() {
    gErrorIgnoreLevel = kError;
    SetPlotStyle();

    // Each entry: { HERNTupler input ntuple, MET emulator output }
    const std::string sigInput  = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/ZvvHbb_v3/mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_resim_DAOD_NTUPLE_GEP.root";
    const std::string backInput = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/mc21_14TeV_jj_JZ_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    const std::string emuDir    = "/data/larsonma/GEPMET/outputNTuplesDev_METv2/";

    // The trailing _twrSF{x}_jetSF{y} tag encodes the (tower, jet) scale factors
    // applied in the totalMET sum (see makeOutputMETFileName). Compare 1,1 vs 1,0p4
    // to see the effect of down-weighting the jet contribution. Future: per-eta calibrated SFs.
    std::vector<std::pair<std::string, std::string>> signalFiles = {
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt10_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt10_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
    };

    std::vector<std::pair<std::string, std::string>> backgroundFiles = {
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt10_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt10_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
    };

    // Label for each file pair — shown in legend for multi-file overlays - NEEDS TO BE UPDATED AND THE SAME SIZE AS THE FILES BEING PROCESSED!
    std::vector<std::string> labels = {
        "J10_T2_EtaSK_OR_TC1_JC1",
        "J10_T2_EtaSK_OR_TC0.4_JC1",
        "J15_T2_EtaSK_OR_TC1_JC1",
        "J15_T2_EtaSK_OR_TC0.4_JC1",
        "J20_T2_EtaSK_OR_TC1_JC1",
        "J20_T2_EtaSK_OR_TC0.4_JC1",
    };

    std::string signalName = "Z #rightarrow #nu#bar{#nu}, H #rightarrow b#bar{b}";
    std::string outputDir  = "metAnalysisPlots/";
    gSystem->mkdir(outputDir.c_str(), true);

    analyze_files(signalFiles, backgroundFiles, labels, outputDir, signalName);
}
