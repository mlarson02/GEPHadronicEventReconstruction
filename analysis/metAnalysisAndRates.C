// To execute: root -b -l -q 'metAnalysisAndRates.C+'

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
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
#include "TLatex.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "analysisHelperFunctions.h"
#include "chainSource.h"

#include "/home/larsonma/atlasrootstyle/AtlasStyle.C"

// Colorblind-friendly palettes (Petroff 6/8/10).
// 6-color palette
const int kP6Blue   = TColor::GetColor("#5790fc");
const int kP6Yellow = TColor::GetColor("#f89c20");
const int kP6Red    = TColor::GetColor("#e42536");
const int kP6Grape  = TColor::GetColor("#964a8b");
const int kP6Gray   = TColor::GetColor("#9c9ca1");
const int kP6Violet = TColor::GetColor("#7a21dd");

// 8-color palette
const int kP8Blue   = TColor::GetColor("#1845fb");
const int kP8Orange = TColor::GetColor("#ff5e02");
const int kP8Red    = TColor::GetColor("#c91f16");
const int kP8Pink   = TColor::GetColor("#c849a9");
const int kP8Green  = TColor::GetColor("#adad7d");
const int kP8Cyan   = TColor::GetColor("#86c8dd");
const int kP8Azure  = TColor::GetColor("#578dff");
const int kP8Gray   = TColor::GetColor("#656364");

// 10-color palette
const int kP10Blue   = TColor::GetColor("#3f90da");
const int kP10Yellow = TColor::GetColor("#ffa90e");
const int kP10Red    = TColor::GetColor("#bd1f01");
const int kP10Gray   = TColor::GetColor("#94a4a2");
const int kP10Violet = TColor::GetColor("#832db6");
const int kP10Brown  = TColor::GetColor("#a96b59");
const int kP10Orange = TColor::GetColor("#e76300");
const int kP10Green  = TColor::GetColor("#b9ac70");
const int kP10Ash    = TColor::GetColor("#717581");
const int kP10Cyan   = TColor::GetColor("#92dadd");

// Process label (e.g. signal name) drawn at the top-right of the ATLAS label on every plot.
// Set per-file before the per-file plots and cleared for multi-file overlays. Empty = nothing drawn.
std::string gProcLabel = "";

// Average pileup quoted on the info line. Set from the sample being processed (see
// SetPileupFromPath) so PU140 samples are not labelled as PU200.
int gPileup = 200;

// Is this a PU140 path? Two conventions are recognised:
//   * the reconstruction tag in the file name — r16129 = PU140, r16130 = PU200.
//     This is the primary differentiator: HERNTupler stamps it into the ntuple
//     names and the emulation carries it through to its own outputs, which land
//     in the same directory for both pileup scenarios.
//   * a PU140 / pu140 / _140 marker in the path, e.g. the ntuples_PU140 directory.
bool IsPU140Path(const std::string& path) {
    return path.find("r16129") != std::string::npos ||
           path.find("PU140")  != std::string::npos ||
           path.find("pu140")  != std::string::npos ||
           path.find("_140")   != std::string::npos;
}

// PU140 inputs carry the pileup in their filename; PU200 inputs are unlabelled (see the
// sum-of-weights convention in HERNTupler.C). Detect from any input path.
void SetPileupFromPath(const std::string& path) {
    gPileup = IsPU140Path(path) ? 140 : 200;
}

// Draw the ATLAS "Work in progress" label (plus beam-energy / pileup info) in a
// white strip ABOVE the plot frame on the currently active canvas.
// Call after cd()'ing to the canvas and before SaveAs/Print.
// The top margin is enlarged so the frame shrinks down and leaves room above it.
void DrawATLASLabel(double x = 0.20, double /*y*/ = 0.88, const char* status = "Work in progress") {
    if (gPad) {
        gPad->SetTopMargin(0.14);   // frame top now ~0.86, leaving a white strip above
        gPad->Modified();
        gPad->Update();
    }
    const double yAtlas = 0.945;    // "ATLAS <status>" line, in the strip above the frame
    const double yInfo  = 0.895;    // beam-energy / pileup line, just below it
    TLatex l; l.SetNDC(); l.SetTextFont(72); l.SetTextColor(kBlack); l.SetTextSize(0.04);
    l.DrawLatex(x, yAtlas, "ATLAS");
    TLatex p; p.SetNDC(); p.SetTextFont(42); p.SetTextColor(kBlack); p.SetTextSize(0.04);
    p.DrawLatex(x + 0.13, yAtlas, status);
    TLatex e; e.SetNDC(); e.SetTextFont(42); e.SetTextColor(kBlack); e.SetTextSize(0.035);
    // #LT / #GT are TLatex's angle brackets; plain "<PU>" renders less-than/greater-than glyphs.
    e.DrawLatex(x, yInfo, Form("#sqrt{s} = 14 TeV, #LTPU#GT = %d", gPileup));
    // Process label at the top-right of the strip (right-aligned), to the right of "ATLAS <status>".
    if (!gProcLabel.empty()) {
        TLatex s; s.SetNDC(); s.SetTextFont(42); s.SetTextColor(kBlack); s.SetTextSize(0.038);
        s.SetTextAlign(31);
        s.DrawLatex(0.95, yAtlas, gProcLabel.c_str());
    }
}

// Background-only plots must not carry the signal process label — they are the same dijet
// sample whatever signal is being processed alongside. DrawATLASLabel reads gProcLabel from
// inside every drawing helper, so swap the global for the duration of a draw and restore it
// afterwards.
const std::string kBkgProcLabel = "QCD dijet";
struct BkgProcLabel {
    std::string saved;
    BkgProcLabel() : saved(gProcLabel) { gProcLabel = kBkgProcLabel; }
    ~BkgProcLabel() { gProcLabel = saved; }
};

const int   nColors = 7;
int cols[nColors] = { kP10Red, kP10Blue, kP10Green, kP10Violet, kP10Orange, kP10Cyan};

// Skip background events with passHSTP == false (HSTP filter removes high-energy
// pileup transients that are not modelled correctly in dijet MC).
const bool applyHSTPFilter = true;

// Progress printouts through the long event loops and the tree-opening stages, so a slow run
// can be told apart from a stuck one. Every line is flushed, because a redirected stdout is
// block-buffered and would otherwise hold the output back and look like a hang by itself.
const bool printIOProgress   = true;
const Long64_t progressEvery = 10000;   // events between progress lines

// --- Rate normalization after the HSTP filter ----------------------------------
// HERNTupler normalizes the per-event weights so that the sum over the UNFILTERED sample is
// its targetRate (30 MHz). The HSTP filter above then removes ~99% of that weight — almost
// all of JZ0, which alone is 95.4% of the target, plus 79% of JZ1 — and nothing puts it back,
// so the weighted background sums to ~323 kHz instead of 30 MHz (measured on the v4 PU200
// ntuples with checkHSTPNormalization.C, 2026-08-12).
//
// This is corrected once, at the histogram level, immediately after the background loop has
// filled them: every weighted background histogram is scaled so that its rate at threshold 0 —
// its integral, which is bin 1 of the cumulative — is targetTotalRateHz. Everything downstream
// is a cumulative sum of those same histograms, so the rate-vs-threshold curves, the
// 20/40/60/80 kHz thresholds, the rate-vs-efficiency curves and the turn-on curves evaluated at
// those thresholds all come out normalized, rather than only the plots that rescale explicitly.
// Set normalizeRateToTarget = false to go back to the weights exactly as filled.
//
// CAVEAT on the normalized curves: the surviving sample is enriched in slices that pass the
// filter (JZ2 keeps 96.8% of its weight, JZ3 99.9%, against JZ0's 0.064%), so scaling it hands
// the full collision rate to a jet-enriched population, and every threshold moves by the same
// factor — including the high-threshold region fed by JZ2/JZ3 where the normalization already
// closed. Check JZSlices/*_JZSlices_Rate.pdf for which slices populate a given working point.
// Scale target: the colliding-bunch crossing rate, 30.9 MHz — the rate an L1 trigger sees at
// zero threshold. Taken from kCrossingRateHz in analysisHelperFunctions.h so this rescale and
// the jet binomial conversion use the same number. Deliberately NOT HERNTupler's targetRate
// (30.0 MHz), which is a normalization choice for the event weights, not a crossing rate.
const double targetTotalRateHz = kCrossingRateHz;
// false = leave every histogram with the normalization the weights carry as filled.
const bool normalizeRateToTarget = true;
const unsigned int analysisPileup = 200;   // selects L_inst for the sigma x L overlay curve

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

// Threshold axis cut-off for the rate-vs-threshold plots [GeV]. The histograms run to 600 GeV,
// but the last 200 GeV carry no rate worth reading, so the curves are cropped for display. This
// only sets the axis range — the cumulative rate at every threshold still integrates the full
// spectrum above it, including the part beyond the cut-off. Callers that want a different range
// (e.g. the 200 GeV zoomed algorithm comparison) pass their own xMax.
const double kRateVsThrXMax = 400.0;

// --- MET types shared by the Z->mumu turn-ons and the MET-vs-jet-multiplicity profiles ------
// One entry per MET flavour, in a fixed order that both the per-file arrays and the multi-file
// vectors index into. Everything else (rate histograms, per-event values, labels) is looked up
// through this order, so a new MET type only has to be added here and in the two arrays that
// bind it to its histogram / branch variable inside analyze_files.
const int   nMETTypes = 7;
const char* metTypeShort[nMETTypes] = {
    "gFEX_JwoJ", "gFEX_NoiseCut", "gFEX_Rms", "jFEX", "JetMET", "TowerMET", "TotalMET"
};
const char* metTypeLabel[nMETTypes] = {
    "gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Tower MET", "GEP Total MET"
};
// The GEP entries of the table above, for the GEP-only overlays.
const int nGEPMETTypes = 3;
const int gepMETTypeIdx[nGEPMETTypes] = { 4, 5, 6 };   // Jet, Tower, Total
// Tower MET is meaningless once Overlap Removal is on, so this entry is dropped from the
// overlays for OR configs (same convention as the GEP algorithm-comparison plots below).
const int towerMETTypeIdx = 5;

// Rate points the Z->mumu dimuon-pT turn-on curves are matched to.
const int    nMuRates = 3;
const double muRatesHz[nMuRates]   = { 40e3, 60e3, 80e3 };
const char*  muRateNames[nMuRates] = { "40kHz", "60kHz", "80kHz" };

// --- Background MET vs jet multiplicity ------------------------------------------------------
// Jet multiplicity counts the truth AntiKt4 WZ-dressed jets (hard scatter) and the in-time
// AntiKt4 truth jets (pileup overlay) together, each above kNJetMinEt: an event's jet activity
// as MET sees it does not care which collision a jet came from.
const int    nNJetBins  = 20;
const double nJetAxisMax = 20.0;
const double kNJetMinEt  = 15.0;   // [GeV]

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

    sig->SetLineColor(kP10Red);   sig->SetLineWidth(2);
    back->SetLineColor(kP10Blue); back->SetLineWidth(2);
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

    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
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

    sig->SetLineColor(kP10Red);   sig->SetLineWidth(2);
    back->SetLineColor(kP10Blue); back->SetLineWidth(2);
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

    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
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
    style(sig1,  kP10Red,  1); style(back1, kP10Red,  2);
    style(sig2,  kP10Blue, 1); style(back2, kP10Blue, 2);

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

    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
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
    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
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

    int mcols[] = { kBlack, kP10Red, kP10Blue, kP10Green, kP10Violet, kP10Orange, kP10Cyan, kP10Brown };
    const int nMcols = 8;
    std::vector<TGraph*> graphs;
    for (unsigned int i = 0; i < sigs.size(); i++) {
        TGraph* g = makeEffVsThresholdGraph(sigs[i], mcols[i % nMcols]);
        g->SetTitle((title + ";" + xLabel + ";Signal Efficiency").c_str());
        graphs.push_back(g);
        g->Draw(i == 0 ? "AP" : "P SAME");
        leg.AddEntry(g, labels[i].c_str(), "l");
    }
    // Fix y-axis range after drawing
    if (!graphs.empty()) {
        graphs[0]->GetYaxis()->SetRangeUser(0.0, 1.05);
    }
    leg.Draw();
    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
    for (auto* g : graphs) delete g;
}

// -----------------------------------------------------------------------
// Build a TGraphErrors of rate vs threshold from a weighted background histogram.
// Rate error = sqrt(sum of squared bin errors) from threshold bin upward.
TGraphErrors* makeRateGraph(TH1F* h, Color_t col, Style_t markerStyle = 20) {
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
    g->SetMarkerColor(col); g->SetMarkerStyle(markerStyle); g->SetMarkerSize(0.7);
    return g;
}

// -----------------------------------------------------------------------
// Rate vs threshold from a weighted background histogram (rate in Hz)
void drawRateVsThreshold(TH1F* back_weighted, const std::string& title,
                         const std::string& xLabel, const std::string& outputPath,
                         double xMax = kRateVsThrXMax) {
    BkgProcLabel bkgProc;   // background-only plot: label as QCD dijet, not the signal
    TGraphErrors* g = makeRateGraph(back_weighted, kP10Blue);
    g->SetTitle((title + ";" + xLabel + ";Rate [Hz]").c_str());

    TCanvas c("c", title.c_str(), 700, 600);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
    gPad->SetLogy();
    g->Draw("AP");
    if (xMax > 0.0)
        g->GetXaxis()->SetLimits(back_weighted->GetXaxis()->GetXmin(), xMax);
    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
    delete g;
}

// -----------------------------------------------------------------------
// Rate vs threshold overlay for multiple algorithm configs
void drawRateVsThresholdMulti(const std::vector<TH1F*>& backs_weighted,
                              const std::vector<std::string>& labels,
                              const std::string& title, const std::string& xLabel,
                              const std::string& outputPath, const std::string& signalName = "",
                              const std::string& yLabel = "Rate [Hz]",
                              double xMax = kRateVsThrXMax, double yScale = 1.0, double yMin = -1.0,
                              double yMax = -1.0) {
    if (backs_weighted.empty()) return;
    BkgProcLabel bkgProc;   // background-only plot: label as QCD dijet, not the signal
    TCanvas c("c", title.c_str(), 700, 600);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
    gPad->SetLogy();

    int nConfigs = (int)backs_weighted.size();
    double legTop = 0.92, legH = 0.04 * (!signalName.empty() + nConfigs);
    TLegend leg(0.54, legTop - legH, 0.88, legTop);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);
    if (!signalName.empty())
        leg.AddEntry((TObject*)nullptr, signalName.c_str(), "");

    int mcols[] = { kBlack, kP10Red, kP10Blue, kP10Green, kP10Violet, kP10Orange, kP10Cyan, kP10Brown };
    const int nMcols = 8;
    const Style_t mkstyles[] = { 20, 21, 22, 23, 29, 33, 34, 47 };
    const int nStyles = 8;
    std::vector<TGraphErrors*> graphs;
    for (unsigned int i = 0; i < backs_weighted.size(); i++) {
        TGraphErrors* g = makeRateGraph(backs_weighted[i], mcols[i % nMcols], mkstyles[i % nStyles]);
        if (yScale != 1.0)
            for (int p = 0; p < g->GetN(); ++p) {
                g->SetPoint(p, g->GetX()[p], g->GetY()[p] * yScale);
                g->SetPointError(p, g->GetEX()[p], g->GetEY()[p] * yScale);
            }
        g->SetTitle((title + ";" + xLabel + ";" + yLabel).c_str());
        graphs.push_back(g);
        g->Draw(i == 0 ? "AP" : "P SAME");
        if (i == 0 && xMax > 0.0)
            g->GetXaxis()->SetLimits(backs_weighted[0]->GetXaxis()->GetXmin(), xMax);
        if (i == 0 && yMin > 0.0)
            g->SetMinimum(yMin);
        // The frame is owned by the first graph, so a later curve that starts higher would be
        // clipped off the top; callers overlaying curves of very different scale pass yMax.
        if (i == 0 && yMax > 0.0)
            g->SetMaximum(yMax);
        leg.AddEntry(g, labels[i].c_str(), "lp");
    }
    leg.Draw();
    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
    for (auto* g : graphs) delete g;
}

// -----------------------------------------------------------------------
// Per-JZ-slice rate-vs-threshold overlay drawn as points (markers, not HIST),
// mirroring OverlayAndSave's 10-colour JZ palette + right-side JZ0..JZ9 legend.
// Rate is the cumulative weighted integral above each threshold; yScale converts
// Hz -> kHz (1e-3). xMax caps the threshold axis (e.g. 200 GeV).
void OverlayRateAndSave(TH1F* hists[], int n, const char* canvasName,
                        TString pdfOut, const char* legHeader,
                        double xMax = 200.0, double yScale = 1e-3,
                        const char* yLabel = "Estimated Background Rate [kHz]",
                        double yMin = 1e-8) {
    if (n <= 0 || !hists[0]) return;
    BkgProcLabel bkgProc;   // per-JZ-slice rates are background only

    // Same palette/order as OverlayAndSave so the JZ colours stay consistent.
    const Color_t colors[10] = {
        kRed+1, kGreen+1, kBlue, kYellow+1, kMagenta+1,
        kCyan+1, kTeal+2, kViolet+1, kGray+1, kGray+3
    };
    const Style_t mkstyles[10] = { 20, 21, 22, 23, 29, 33, 34, 47, 43, 45 };

    TCanvas* c = new TCanvas(canvasName, canvasName, 900, 700);
    c->SetMargin(0.12, 0.22, 0.16, 0.06); // room for right-side legend
    c->SetTicks(1, 1);
    c->SetLogy();

    TLegend* leg = new TLegend(0.80, 0.15, 0.97, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.030);
    if (legHeader && legHeader[0] != '\0')
        leg->AddEntry((TObject*)nullptr, legHeader, "");

    std::vector<TGraphErrors*> graphs;
    double yMaxSeen = yMin;
    for (int i = 0; i < n; ++i) {
        if (!hists[i]) continue;
        TGraphErrors* g = makeRateGraph(hists[i], colors[i % 10], mkstyles[i % 10]);
        if (yScale != 1.0)
            for (int p = 0; p < g->GetN(); ++p) {
                g->SetPoint(p, g->GetX()[p], g->GetY()[p] * yScale);
                g->SetPointError(p, g->GetEX()[p], g->GetEY()[p] * yScale);
            }
        for (int p = 0; p < g->GetN(); ++p)
            if (g->GetY()[p] > yMaxSeen) yMaxSeen = g->GetY()[p];
        g->SetTitle(TString::Format(";%s;%s", hists[i]->GetXaxis()->GetTitle(), yLabel));
        graphs.push_back(g);
        g->Draw(graphs.size() == 1 ? "AP" : "P SAME");
        if (graphs.size() == 1 && xMax > 0.0)
            g->GetXaxis()->SetLimits(hists[i]->GetXaxis()->GetXmin(), xMax);
        leg->AddEntry(g, TString::Format("JZ%d", i), "lp");
    }
    // Pin the log-y frame: floor at yMin (default 1e-8), headroom above the tallest slice.
    if (!graphs.empty()) {
        graphs[0]->SetMinimum(yMin);
        graphs[0]->SetMaximum(yMaxSeen * 5.0);
    }
    leg->Draw();
    gPad->Modified(); gPad->Update(); gPad->RedrawAxis();
    c->cd(); DrawATLASLabel(); c->SaveAs(pdfOut);
    for (auto* g : graphs) delete g;
    delete c;
}

// -----------------------------------------------------------------------
// Overlay N normalized distributions on one canvas (no sig/bkg distinction)
void drawMultiDist(std::vector<TH1F*> hists, const std::vector<std::string>& labels,
                   const std::string& title, const std::string& xLabel,
                   const std::string& outputPath, bool logy = true,
                   const std::string& units = "GeV") {
    if (hists.empty()) return;
    int mcols[] = { kBlack, kP10Red, kP10Blue, kP10Green, kP10Violet,
                               kP10Orange, kP10Cyan, kP10Brown, kGray+2, kP10Yellow };
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
    TLegend leg(0.58, legTop - legH, 0.88, legTop);
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
    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
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
// xLabel defaults to the truth MET axis these curves were written for; the Z->mumu turn-ons
// pass the dimuon p_{T} axis instead.
void drawTurnOnOverlay(std::vector<TH1F*> effs, const std::vector<std::string>& labels,
                       const std::string& title, const std::string& outputPath,
                       const std::vector<double>& thresholds = {},
                       const std::string& rateLabel = "",
                       TH1F* truthDist = nullptr,
                       double legX1 = 0.45, double legY1 = 0.15,
                       const std::string& xLabel = "Truth MET_{NonInt} [GeV]") {
    if (effs.empty()) return;
    int  mcols[]    = { kBlack, kP10Red, kP10Blue, kP10Green, kP10Violet, kP10Orange, kP10Cyan};
    // Seven entries, matching mcols: both are indexed with i % 7 below, and the widest overlays
    // (all seven MET types) reach index 6.
    const Style_t  mkstyles[] = { 20, 21, 22, 23, 29, 33, 34 };

    TCanvas c("c", title.c_str(), 700, 600);
    gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);

    int nLegRows = (int)effs.size() + (!rateLabel.empty() ? 1 : 0);
    double legH = 0.058 * nLegRows;
    TLegend leg(legX1, legY1, legX1 + 0.33, legY1 + legH);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032);
    if (!rateLabel.empty())
        leg.AddEntry((TObject*)nullptr, rateLabel.c_str(), "");

    for (unsigned int i = 0; i < effs.size(); i++) {
        effs[i]->SetLineColor(mcols[i % 7]);
        effs[i]->SetMarkerColor(mcols[i % 7]);
        effs[i]->SetMarkerStyle(mkstyles[i % 7]);
        effs[i]->SetMarkerSize(0.8);
        effs[i]->SetLineWidth(2);
        effs[i]->SetTitle(title.c_str());
        effs[i]->GetXaxis()->SetTitle(xLabel.c_str());
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
    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
}

// -----------------------------------------------------------------------
// Overlay N TProfile curves (mean of y in each x bin) on one canvas, drawn as points with
// error bars. Used for the background <MET> vs jet-multiplicity plots, where every profile is
// filled with the rate weights so the mean is the weighted mean.
void drawProfileOverlay(const std::vector<TProfile*>& profs,
                        const std::vector<std::string>& labels,
                        const std::string& xLabel, const std::string& yLabel,
                        const std::string& outputPath, const std::string& legHeader = "",
                        double legX1 = 0.20, double legY1 = 0.60) {
    if (profs.empty()) return;
    int mcols[] = { kBlack, kP10Red, kP10Blue, kP10Green, kP10Violet,
                    kP10Orange, kP10Cyan, kP10Brown, kGray+2, kP10Yellow };
    const Style_t mkstyles[] = { 20, 21, 22, 23, 29, 33, 34, 47, 43, 45 };
    const int nStyles = 10;

    TCanvas c("c", "", 700, 600);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);

    int nLegRows = (int)profs.size() + (!legHeader.empty() ? 1 : 0);
    TLegend leg(legX1, legY1, legX1 + 0.36, legY1 + 0.05 * nLegRows);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);
    if (!legHeader.empty())
        leg.AddEntry((TObject*)nullptr, legHeader.c_str(), "");

    // The frame belongs to the first profile drawn, so a later curve that runs higher would be
    // clipped: take the maximum over all of them (mean + error) before drawing anything.
    double ymax = 0.0;
    for (auto* p : profs) {
        if (!p) continue;
        for (int ib = 1; ib <= p->GetNbinsX(); ++ib)
            if (p->GetBinEntries(ib) > 0)
                ymax = std::max(ymax, p->GetBinContent(ib) + p->GetBinError(ib));
    }
    if (ymax <= 0.0) ymax = 1.0;

    bool first = true;
    for (unsigned int i = 0; i < profs.size(); i++) {
        if (!profs[i]) continue;
        profs[i]->SetLineColor(mcols[i % nStyles]);
        profs[i]->SetMarkerColor(mcols[i % nStyles]);
        profs[i]->SetMarkerStyle(mkstyles[i % nStyles]);
        profs[i]->SetMarkerSize(0.9);
        profs[i]->SetLineWidth(2);
        profs[i]->SetTitle("");
        profs[i]->GetXaxis()->SetTitle(xLabel.c_str());
        profs[i]->GetYaxis()->SetTitle(yLabel.c_str());
        profs[i]->SetMinimum(0.0);
        profs[i]->SetMaximum(ymax * 1.35);
        profs[i]->Draw(first ? "E1" : "E1 SAME");
        first = false;
        leg.AddEntry(profs[i], labels[i].c_str(), "lp");
    }
    leg.Draw();
    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
}

// -----------------------------------------------------------------------
// Overlay multiple Rate-vs-Efficiency TGraph* on one canvas
void drawRateVsEffOverlay(std::vector<TGraph*> graphs,
                          const std::vector<std::string>& labels,
                          const std::string& outputPath,
                          const std::string& signalName = "") {
    if (graphs.empty()) return;
    int  mcols[]    = { kBlack, kP10Red, kP10Blue, kP10Green, kP10Violet, kP10Orange, kP10Cyan, kP10Brown };
    const Style_t  mkstyles[] = { 20, 21, 22, 23, 29, 33, 20, 21 };
    const int nStyles = 8;

    TCanvas c("c", "", 700, 600);
    gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
    gPad->SetLogy(); gPad->SetLogx();

    int nLeg = (int)(!signalName.empty()) + (int)graphs.size();
    TLegend leg(0.2, 0.45, 0.43, 0.49 + 0.05 * nLeg);
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.025);
    if (!signalName.empty())
        leg.AddEntry((TObject*)nullptr, signalName.c_str(), "");

    // Clamp the y-axis to a 10 Hz floor (the curves dip far below physical interest and
    // would otherwise push the frame down to ~1e-3 Hz, crowding the legend).
    /*double ymax = 0.0;
    for (auto* g : graphs)
        for (int p = 0; p < g->GetN(); ++p) {
            double x, y; g->GetPoint(p, x, y);
            if (y > ymax) ymax = y;
        }
    if (ymax <= 0) ymax = 1e8;
    graphs[0]->SetMinimum(10.0);
    graphs[0]->SetMaximum(ymax * 3.0);*/

    // Enforce axis floors on this log-log plot: signal efficiency (x) >= 1e-5, rate (y) >= 1e-3.
    /*double xmax = 0.0;
    for (auto* g : graphs)
        for (int p = 0; p < g->GetN(); ++p) {
            double x, y; g->GetPoint(p, x, y);
            if (x > xmax) xmax = x;
        }
    if (xmax <= 1e-5) xmax = 1.0;
    graphs[0]->SetMinimum(1e-3);   // y-axis (rate) floor*/

    for (unsigned int i = 0; i < graphs.size(); i++) {
        graphs[i]->SetLineColor(mcols[i % nStyles]);
        graphs[i]->SetMarkerColor(mcols[i % nStyles]);
        graphs[i]->SetMarkerStyle(mkstyles[i % nStyles]);
        graphs[i]->SetMarkerSize(0.8);
        graphs[i]->SetLineWidth(2);
        graphs[i]->GetXaxis()->SetTitle("Signal Efficiency");
        graphs[i]->GetYaxis()->SetTitle("Estimated Background Rate [Hz]");
        graphs[i]->Draw(i == 0 ? "AP" : "P SAME");
        //if (i == 0) graphs[0]->GetXaxis()->SetLimits(1e-5, xmax * 1.05); // x-axis (efficiency) floor
        leg.AddEntry(graphs[i], labels[i].c_str(), "lp");
    }
    gPad->Modified(); gPad->Update();
    leg.Draw();
    c.cd(); DrawATLASLabel(); c.SaveAs(outputPath.c_str());
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
// Per-event effective SoftKiller threshold = smallest non-zero Et among
// surviving towers (post-SK). gepCellsTowers{SK,EtaSK}Tree stores all towers
// with Et=0 for killed ones, so this is min(Et) over the kept subset.
// Reads only from the HERNTupler input ntuples (same input across all MET-emu
// configs), so we run this once per analyze_files call.
void plotSKThresholds(const std::string& sigPath, const std::string& backPath,
                      const std::string& outputDir) {
    if (printIOProgress)
        std::cout << "plotSKThresholds: opening signal " << sigPath << "\n" << std::flush;
    ChainSource* sigF  = ChainSource::Open(sigPath.c_str());
    if (printIOProgress)
        std::cout << "plotSKThresholds: opening background " << backPath << "\n" << std::flush;
    ChainSource* backF = ChainSource::Open(backPath.c_str());
    if (!sigF || sigF->IsZombie() || !backF || backF->IsZombie()) {
        std::cerr << "plotSKThresholds: cannot open input ntuples\n";
        return;
    }
    TTree* skSig    = (TTree*)sigF->Get("gepCellsTowersSKTree");
    TTree* skBack   = (TTree*)backF->Get("gepCellsTowersSKTree");
    TTree* etaSig   = (TTree*)sigF->Get("gepCellsTowersEtaSKTree");
    TTree* etaBack  = (TTree*)backF->Get("gepCellsTowersEtaSKTree");
    TTree* evtBack  = (TTree*)backF->Get("eventInfoTree");
    if (!skSig || !skBack || !etaSig || !etaBack || !evtBack) {
        std::cerr << "plotSKThresholds: missing SK/EtaSK/eventInfo tree\n";
        sigF->Close(); backF->Close();
        return;
    }

    std::vector<double>* skSigEt   = nullptr;
    std::vector<double>* skBackEt  = nullptr;
    std::vector<double>* etaSigEt  = nullptr;
    std::vector<double>* etaBackEt = nullptr;
    skSig->SetBranchAddress("Et",   &skSigEt);
    skBack->SetBranchAddress("Et",  &skBackEt);
    etaSig->SetBranchAddress("Et",  &etaSigEt);
    etaBack->SetBranchAddress("Et", &etaBackEt);

    std::vector<double>* eventWeightsValuesBack = nullptr;
    bool passHSTPValuesBack = true;
    evtBack->SetBranchAddress("eventWeights", &eventWeightsValuesBack);
    evtBack->SetBranchAddress("passHSTP",     &passHSTPValuesBack);

    const int nBins = 40;
    const double xlo = 0.0, xhi = 4.0; // GeV
    TH1F* sig_h_SKThresh    = new TH1F("sig_h_SKThresh",    "", nBins, xlo, xhi);
    TH1F* back_h_SKThresh   = new TH1F("back_h_SKThresh",   "", nBins, xlo, xhi);
    TH1F* sig_h_EtaSKThresh = new TH1F("sig_h_EtaSKThresh", "", nBins, xlo, xhi);
    TH1F* back_h_EtaSKThresh= new TH1F("back_h_EtaSKThresh","", nBins, xlo, xhi);
    sig_h_SKThresh->SetDirectory(0);   back_h_SKThresh->SetDirectory(0);
    sig_h_EtaSKThresh->SetDirectory(0);back_h_EtaSKThresh->SetDirectory(0);

    auto minPositive = [](const std::vector<double>* v) {
        double mn = 1e18;
        bool found = false;
        for (double et : *v) {
            if (et > 0.0 && et < mn) { mn = et; found = true; }
        }
        return found ? mn : -1.0;
    };

    Long64_t nSig = std::min(skSig->GetEntries(), etaSig->GetEntries());
    std::cout << "plotSKThresholds: signal events = " << nSig << "\n" << std::flush;
    TStopwatch swSK; swSK.Start();
    for (Long64_t i = 0; i < nSig; ++i) {
        if (printIOProgress && i > 0 && i % progressEvery == 0) {
            std::cout << "  plotSKThresholds: signal " << i << "/" << nSig
                      << " (" << swSK.RealTime() << " s)\n" << std::flush;
            swSK.Continue();   // RealTime() stops the watch
        }
        skSig->GetEntry(i);
        etaSig->GetEntry(i);
        double tSK    = minPositive(skSigEt);
        double tEtaSK = minPositive(etaSigEt);
        if (tSK    >= 0.0) sig_h_SKThresh->Fill(std::min(tSK,    xhi - 1e-9));
        if (tEtaSK >= 0.0) sig_h_EtaSKThresh->Fill(std::min(tEtaSK, xhi - 1e-9));
    }

    Long64_t nBack = std::min({skBack->GetEntries(), etaBack->GetEntries(), evtBack->GetEntries()});
    std::cout << "plotSKThresholds: background events = " << nBack << "\n" << std::flush;
    swSK.Start();
    for (Long64_t i = 0; i < nBack; ++i) {
        if (printIOProgress && i > 0 && i % progressEvery == 0) {
            std::cout << "  plotSKThresholds: background " << i << "/" << nBack
                      << " (" << swSK.RealTime() << " s)\n" << std::flush;
            swSK.Continue();
        }
        evtBack->GetEntry(i);
        if (applyHSTPFilter && !passHSTPValuesBack) continue;
        double w = eventWeightsValuesBack->at(0);
        skBack->GetEntry(i);
        etaBack->GetEntry(i);
        double tSK    = minPositive(skBackEt);
        double tEtaSK = minPositive(etaBackEt);
        if (tSK    >= 0.0) back_h_SKThresh->Fill(std::min(tSK,    xhi - 1e-9), w);
        if (tEtaSK >= 0.0) back_h_EtaSKThresh->Fill(std::min(tEtaSK, xhi - 1e-9), w);
    }

    std::string thrDir = outputDir + "SKThresholds/";
    gSystem->mkdir(thrDir.c_str(), true);

    drawOverlay(sig_h_SKThresh,    back_h_SKThresh,
                "Effective SK threshold (gepCellsTowersSK)",
                "Effective SK threshold per event [GeV]",
                thrDir + "SK_EffectiveThreshold.pdf");
    drawOverlay(sig_h_EtaSKThresh, back_h_EtaSKThresh,
                "Effective EtaSK threshold (gepCellsTowersEtaSK)",
                "Effective EtaSK threshold per event [GeV]",
                thrDir + "EtaSK_EffectiveThreshold.pdf");
    drawMultiDist({sig_h_SKThresh,  sig_h_EtaSKThresh,
                   back_h_SKThresh, back_h_EtaSKThresh},
                  {"SK signal", "EtaSK signal", "SK bkg", "EtaSK bkg"},
                  "Effective SoftKiller threshold per event",
                  "Effective threshold [GeV]",
                  thrDir + "SK_vs_EtaSK_EffectiveThreshold.pdf", true, "GeV");

    sigF->Close(); backF->Close();
}

// -----------------------------------------------------------------------
// Each signal/background entry is a pair: .first = HERNTupler input ntuple,
// .second = MET emulator output (contains metTree + emulEventInfoTree only).
void analyze_files(std::vector<std::pair<std::string, std::string>> signalFiles,
                   std::vector<std::pair<std::string, std::string>> backgroundFiles,
                   std::vector<std::string> labels,
                   std::string outputDir,
                   std::string signalName = "", std::string overlayDir = "multiFileOverlay_MET/",
                   std::vector<std::string> signalNames = {}) {
    // signalName    : legend header used for multi-file overlays (all configs).
    // signalNames   : optional per-file legend headers (parallel to signalFiles/labels);
    //                 lets one analyze_files call mix signal processes (e.g. ZvvHbb,
    //                 ttbar semilep, ttbar dilep). When empty or short, falls back to signalName.

    gSystem->mkdir(outputDir.c_str(), true);

    // Effective SoftKiller / EtaSoftKiller threshold distributions (input
    // ntuples don't vary across emu configs, so do this once).
    if (!signalFiles.empty() && !backgroundFiles.empty())
        plotSKThresholds(signalFiles[0].first, backgroundFiles[0].first, outputDir);

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
    std::vector<TH1F*> back_hw_SumJetET_vec;   // H_T
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

    // 60 kHz efficiency histograms and thresholds per file for multi-file turn-on comparison
    std::vector<TH1F*> eff_gMET_60kHz_vec;
    std::vector<TH1F*> eff_gMET_NC_60kHz_vec;
    std::vector<TH1F*> eff_gMET_Rms_60kHz_vec;
    std::vector<TH1F*> eff_jMET_60kHz_vec;
    std::vector<TH1F*> eff_JetMET_60kHz_vec;
    std::vector<TH1F*> eff_TowerMET_60kHz_vec;
    std::vector<TH1F*> eff_TotalMET_60kHz_vec;
    std::vector<double> thr_gMET_60kHz_vec;
    std::vector<double> thr_gMET_NC_60kHz_vec;
    std::vector<double> thr_gMET_Rms_60kHz_vec;
    std::vector<double> thr_jMET_60kHz_vec;
    std::vector<double> thr_JetMET_60kHz_vec;
    std::vector<double> thr_TowerMET_60kHz_vec;
    std::vector<double> thr_TotalMET_60kHz_vec;
    std::vector<TH1F*> sig_h_metTruthNonInt_unscaled_vec; // unscaled clone for truth overlay

    // Resimulated gFEX MET histogram vectors (populated only when hasGFexSimMET)
    std::vector<TH1F*> sig_h_gMET_JwoJAOD_vec, back_h_gMET_JwoJAOD_vec;
    std::vector<TH1F*> sig_h_gMET_NCAOD_vec,   back_h_gMET_NCAOD_vec;
    std::vector<TH1F*> sig_h_gMET_RmsAOD_vec,  back_h_gMET_RmsAOD_vec;
    std::vector<TH1F*> back_hw_gMET_JwoJAOD_vec;
    std::vector<TH1F*> back_hw_gMET_NCAOD_vec;
    std::vector<TH1F*> back_hw_gMET_RmsAOD_vec;
    std::vector<TH1F*> eff_gMET_JwoJAOD_80kHz_vec;
    std::vector<TH1F*> eff_gMET_NCAOD_80kHz_vec;
    std::vector<TH1F*> eff_gMET_RmsAOD_80kHz_vec;
    std::vector<double> thr_gMET_JwoJAOD_80kHz_vec;
    std::vector<double> thr_gMET_NCAOD_80kHz_vec;
    std::vector<double> thr_gMET_RmsAOD_80kHz_vec;
    std::vector<TH1F*> eff_gMET_JwoJAOD_60kHz_vec;
    std::vector<TH1F*> eff_gMET_NCAOD_60kHz_vec;
    std::vector<TH1F*> eff_gMET_RmsAOD_60kHz_vec;
    std::vector<double> thr_gMET_JwoJAOD_60kHz_vec;
    std::vector<double> thr_gMET_NCAOD_60kHz_vec;
    std::vector<double> thr_gMET_RmsAOD_60kHz_vec;

    // Z->mumu dimuon-p_{T} turn-ons, [MET type][rate point], one entry per Z->mumu input file.
    // Only Z->mumu files contribute (see isZmumuSample below), so these run parallel to
    // zmumuLabels rather than to labels — a run mixing Z->mumu with other signal processes
    // still gets a consistent multi-file overlay.
    std::vector<TH1F*>  effMu_vec[nMETTypes][nMuRates];
    std::vector<double> thrMu_vec[nMETTypes][nMuRates];
    std::vector<std::string> zmumuLabels;
    std::vector<TH1F*>  sig_h_dimuonPt_coarse_vec;   // unnormalized clone for the turn-on overlay

    // Background <MET> vs jet multiplicity, one profile per MET type per file.
    std::vector<TProfile*> back_prof_METvsNJets_vec[nMETTypes];
    std::vector<std::string> nJetProfLabels;

    bool hasSumJetET   = false;
    bool hasSumTowerET = false;

    for (unsigned int fileIt = 0; fileIt < signalFiles.size(); fileIt++) {
        std::cout << "Processing file " << fileIt << ": " << labels[fileIt] << "\n";
        // HERNTupler input ntuples: provide all upstream trees (gFEX MET, truth, gepCellsTowers jets, etc.)
        ChainSource* sigF  = ChainSource::Open(signalFiles[fileIt].first.c_str());
        ChainSource* backF = ChainSource::Open(backgroundFiles[fileIt].first.c_str());
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
        // TTree* coreEMTopoTreeSig       = (TTree*)sigF->Get("metCoreAntiKt4EMTopoTree");
        // TTree* coreEMTopoTreeBack      = (TTree*)backF->Get("metCoreAntiKt4EMTopoTree");
        // TTree* coreEMPFlowTreeSig      = (TTree*)sigF->Get("metCoreAntiKt4EMPFlowTree");
        // TTree* coreEMPFlowTreeBack     = (TTree*)backF->Get("metCoreAntiKt4EMPFlowTree");
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
        // Leading truth jet of each event, entry 0 of its own collection (HERNTupler writes it
        // separately). Feeds the truth-jet binomial correction below; if the tree is missing the
        // correction falls back to the maximum of the full Et vector.
        TTree* leadingTruthAntiKt4WZDressedJetsTreeBack = (TTree*)backF->Get("leadingTruthAntiKt4TruthDressedWZJets");

        // --- In-time pileup truth jet tree (optional) ---
        TTree* inTimeAntiKt4TruthJetsTreeSig  = (TTree*)sigF->Get("inTimeAntiKt4TruthJetsTree");
        TTree* inTimeAntiKt4TruthJetsTreeBack = (TTree*)backF->Get("inTimeAntiKt4TruthJetsTree");
        bool hasInTimeAntiKt4TruthJets = (inTimeAntiKt4TruthJetsTreeSig  != nullptr) &&
                                         (inTimeAntiKt4TruthJetsTreeBack != nullptr);

        // --- Z->mumu dimuon system (signal eventInfoTree, optional) ---
        // HERNTupler builds the dimuon system from the two leading truth muons and only
        // retrieves TruthMuons for the Zmumu sample; every other process leaves the branches at
        // their defaults (dimuonPt = -1). So a dimuon-p_{T} turn-on is only meaningful on Zmumu,
        // and both conditions are required here: the sample has to BE Zmumu and the ntuple has
        // to carry the branch. Without the name check a non-Zmumu ntuple would happily produce a
        // plot built entirely from default values.
        TTree* eventInfoTreeSig = (TTree*)sigF->Get("eventInfoTree");
        const bool isZmumuName = (signalFiles[fileIt].first.find("Zmumu")  != std::string::npos) ||
                                 (signalFiles[fileIt].second.find("Zmumu") != std::string::npos);
        const bool hasDimuonBranch = (eventInfoTreeSig != nullptr) &&
                                     (eventInfoTreeSig->GetBranch("dimuonPt") != nullptr);
        const bool isZmumuSample = isZmumuName && hasDimuonBranch;
        if (isZmumuName && !hasDimuonBranch)
            std::cout << "  Z->mumu sample but no dimuonPt branch — dimuon p_{T} turn-ons skipped\n";
        else if (isZmumuSample)
            std::cout << "  Z->mumu sample — dimuon p_{T} turn-ons enabled\n";

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
        // Nominal gFEX := the RESIMULATED gFEX (read from gFexMET*SimTree) for ALL plots; when
        // resim trees are absent we fall back to the AOD gFEX. The AOD values are kept separately
        // in the *_*AOD variables below, used only by the AOD-vs-Sim comparison block.
        TTree* gNomJwoJSig  = hasGFexSimMET ? gFexMETJwoJSimTreeSig      : gFexMETTreeSig;
        TTree* gNomJwoJBack = hasGFexSimMET ? gFexMETJwoJSimTreeBack     : gFexMETTreeBack;
        TTree* gNomNCSig    = hasGFexSimMET ? gFexMETNoiseCutSimTreeSig  : gFexMETNoiseCutTreeSig;
        TTree* gNomNCBack   = hasGFexSimMET ? gFexMETNoiseCutSimTreeBack : gFexMETNoiseCutTreeBack;
        TTree* gNomRmsSig   = hasGFexSimMET ? gFexMETRmsSimTreeSig       : gFexMETRmsTreeSig;
        TTree* gNomRmsBack  = hasGFexSimMET ? gFexMETRmsSimTreeBack      : gFexMETRmsTreeBack;

        double sig_gMET = 0.0, sig_gSumET = 0.0;
        double back_gMET = 0.0, back_gSumET = 0.0;
        double sig_gMET_NC = 0.0, sig_gSumET_NC = 0.0;
        double back_gMET_NC = 0.0, back_gSumET_NC = 0.0;
        double sig_gMET_Rms = 0.0, sig_gSumET_Rms = 0.0;
        double back_gMET_Rms = 0.0, back_gSumET_Rms = 0.0;

        gNomJwoJSig->SetBranchAddress("gMET",    &sig_gMET);
        gNomJwoJSig->SetBranchAddress("gSumET",  &sig_gSumET);
        gNomJwoJBack->SetBranchAddress("gMET",   &back_gMET);
        gNomJwoJBack->SetBranchAddress("gSumET", &back_gSumET);
        gNomNCSig->SetBranchAddress("gMET",    &sig_gMET_NC);
        gNomNCSig->SetBranchAddress("gSumET",  &sig_gSumET_NC);
        gNomNCBack->SetBranchAddress("gMET",   &back_gMET_NC);
        gNomNCBack->SetBranchAddress("gSumET", &back_gSumET_NC);
        gNomRmsSig->SetBranchAddress("gMET",    &sig_gMET_Rms);
        gNomRmsSig->SetBranchAddress("gSumET",  &sig_gSumET_Rms);
        gNomRmsBack->SetBranchAddress("gMET",   &back_gMET_Rms);
        gNomRmsBack->SetBranchAddress("gSumET", &back_gSumET_Rms);

        // --- Branch variables: jFEX MET (magnitude only) ---
        double sig_jMET = 0.0, back_jMET = 0.0;
        jFexMETTreeSig->SetBranchAddress("jMET",  &sig_jMET);
        jFexMETTreeBack->SetBranchAddress("jMET", &back_jMET);

        // --- Branch variables: AOD gFEX MET (read from gFexMETTree etc.; no X/Y components) ---
        // These hold the AOD gFEX and are used ONLY by the AOD-vs-Sim comparison block.
        double sig_gMET_JwoJAOD = 0.0, sig_gSumET_JwoJAOD = 0.0;
        double back_gMET_JwoJAOD = 0.0, back_gSumET_JwoJAOD = 0.0;
        double sig_gMET_NCAOD = 0.0, sig_gSumET_NCAOD = 0.0;
        double back_gMET_NCAOD = 0.0, back_gSumET_NCAOD = 0.0;
        double sig_gMET_RmsAOD = 0.0, sig_gSumET_RmsAOD = 0.0;
        double back_gMET_RmsAOD = 0.0, back_gSumET_RmsAOD = 0.0;
        if (hasGFexSimMET) {
            gFexMETTreeSig->SetBranchAddress("gMET",   &sig_gMET_JwoJAOD);
            gFexMETTreeSig->SetBranchAddress("gSumET", &sig_gSumET_JwoJAOD);
            gFexMETTreeBack->SetBranchAddress("gMET",   &back_gMET_JwoJAOD);
            gFexMETTreeBack->SetBranchAddress("gSumET", &back_gSumET_JwoJAOD);
            gFexMETNoiseCutTreeSig->SetBranchAddress("gMET",   &sig_gMET_NCAOD);
            gFexMETNoiseCutTreeSig->SetBranchAddress("gSumET", &sig_gSumET_NCAOD);
            gFexMETNoiseCutTreeBack->SetBranchAddress("gMET",   &back_gMET_NCAOD);
            gFexMETNoiseCutTreeBack->SetBranchAddress("gSumET", &back_gSumET_NCAOD);
            gFexMETRmsTreeSig->SetBranchAddress("gMET",   &sig_gMET_RmsAOD);
            gFexMETRmsTreeSig->SetBranchAddress("gSumET", &sig_gSumET_RmsAOD);
            gFexMETRmsTreeBack->SetBranchAddress("gMET",   &back_gMET_RmsAOD);
            gFexMETRmsTreeBack->SetBranchAddress("gSumET", &back_gSumET_RmsAOD);
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
        // double sig_coreEMTopo_SoftClus_MET   = 0.0, sig_coreEMTopo_PVSoftTrk_MET  = 0.0, sig_coreEMTopo_SoftClusEM_MET = 0.0;
        // double back_coreEMTopo_SoftClus_MET  = 0.0, back_coreEMTopo_PVSoftTrk_MET = 0.0, back_coreEMTopo_SoftClusEM_MET = 0.0;

        // coreEMTopoTreeSig->SetBranchAddress("SoftClus_MET",   &sig_coreEMTopo_SoftClus_MET);
        // coreEMTopoTreeSig->SetBranchAddress("PVSoftTrk_MET",  &sig_coreEMTopo_PVSoftTrk_MET);
        // coreEMTopoTreeSig->SetBranchAddress("SoftClusEM_MET", &sig_coreEMTopo_SoftClusEM_MET);
        // coreEMTopoTreeBack->SetBranchAddress("SoftClus_MET",   &back_coreEMTopo_SoftClus_MET);
        // coreEMTopoTreeBack->SetBranchAddress("PVSoftTrk_MET",  &back_coreEMTopo_PVSoftTrk_MET);
        // coreEMTopoTreeBack->SetBranchAddress("SoftClusEM_MET", &back_coreEMTopo_SoftClusEM_MET);

        // --- Branch variables: core MET EMPFlow ---
        // double sig_coreEMPFlow_SoftClus_MET  = 0.0, sig_coreEMPFlow_PVSoftTrk_MET  = 0.0;
        // double back_coreEMPFlow_SoftClus_MET = 0.0, back_coreEMPFlow_PVSoftTrk_MET = 0.0;

        // coreEMPFlowTreeSig->SetBranchAddress("SoftClus_MET",  &sig_coreEMPFlow_SoftClus_MET);
        // coreEMPFlowTreeSig->SetBranchAddress("PVSoftTrk_MET", &sig_coreEMPFlow_PVSoftTrk_MET);
        // coreEMPFlowTreeBack->SetBranchAddress("SoftClus_MET",  &back_coreEMPFlow_SoftClus_MET);
        // coreEMPFlowTreeBack->SetBranchAddress("PVSoftTrk_MET", &back_coreEMPFlow_PVSoftTrk_MET);

        // --- Background event weights + JZ slice index ---
        std::vector<double>* eventWeightsValuesBack = nullptr;
        int sampleJZSliceBack = -1;
        bool passHSTPValuesBack = true;
        eventInfoTreeBack->SetBranchAddress("eventWeights",  &eventWeightsValuesBack);
        eventInfoTreeBack->SetBranchAddress("sampleJZSlice", &sampleJZSliceBack);
        eventInfoTreeBack->SetBranchAddress("passHSTP",      &passHSTPValuesBack);
        // Reconstructed primary vertices with >= 2 tracks, written by HERNTupler. Guarded so
        // older ntuples without the branch still run.
        int nPrimaryVerticesBack = -1;
        const bool hasNPrimaryVertices = (eventInfoTreeBack->GetBranch("nPrimaryVertices") != nullptr);
        if (hasNPrimaryVertices)
            eventInfoTreeBack->SetBranchAddress("nPrimaryVertices", &nPrimaryVerticesBack);
        else
            std::cout << "  nPrimaryVertices branch not found — primary-vertex plot skipped\n";

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
        std::vector<double>* leadingTruthAntiKt4WZDressedJetsEtValuesBack = nullptr;
        if (hasTruthAntiKt4WZDressed) {
            truthAntiKt4WZDressedJetsTreeSig->SetBranchAddress("Et",  &truthAntiKt4WZDressedJetsEtValuesSig);
            truthAntiKt4WZDressedJetsTreeBack->SetBranchAddress("Et", &truthAntiKt4WZDressedJetsEtValuesBack);
            if (leadingTruthAntiKt4WZDressedJetsTreeBack)
                leadingTruthAntiKt4WZDressedJetsTreeBack->SetBranchAddress("Et", &leadingTruthAntiKt4WZDressedJetsEtValuesBack);
        }

        // --- Branch variables: in-time pileup truth jets (optional) ---
        std::vector<double>* inTimeAntiKt4TruthJetsEtValuesSig  = nullptr;
        std::vector<double>* inTimeAntiKt4TruthJetsEtValuesBack = nullptr;
        if (hasInTimeAntiKt4TruthJets) {
            inTimeAntiKt4TruthJetsTreeSig->SetBranchAddress("Et",  &inTimeAntiKt4TruthJetsEtValuesSig);
            inTimeAntiKt4TruthJetsTreeBack->SetBranchAddress("Et", &inTimeAntiKt4TruthJetsEtValuesBack);
        }

        // --- Branch variables: Z->mumu dimuon system (signal only, Zmumu only) ---
        double sig_dimuonPt = -1.0, sig_dimuonMass = -1.0;
        int    sig_nTruthMuons = 0;
        if (isZmumuSample) {
            eventInfoTreeSig->SetBranchAddress("dimuonPt",    &sig_dimuonPt);
            eventInfoTreeSig->SetBranchAddress("dimuonMass",  &sig_dimuonMass);
            eventInfoTreeSig->SetBranchAddress("nTruthMuons", &sig_nTruthMuons);
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
        // TH1F* sig_h_coreEMTopo_SoftClus_MET    = new TH1F(("sig_coreEMTopo_SoftClus_"  +tag).c_str(), "", 70, 0, 350);
        // TH1F* back_h_coreEMTopo_SoftClus_MET   = new TH1F(("back_coreEMTopo_SoftClus_" +tag).c_str(), "", 70, 0, 350);
        // TH1F* sig_h_coreEMTopo_PVSoftTrk_MET   = new TH1F(("sig_coreEMTopo_PVSoftTrk_" +tag).c_str(), "", 20, 0, 100);
        // TH1F* back_h_coreEMTopo_PVSoftTrk_MET  = new TH1F(("back_coreEMTopo_PVSoftTrk_"+tag).c_str(), "", 20, 0, 100);
        // TH1F* sig_h_coreEMTopo_SoftClusEM_MET  = new TH1F(("sig_coreEMTopo_SoftClusEM_" +tag).c_str(), "", 40, 0, 200);
        // TH1F* back_h_coreEMTopo_SoftClusEM_MET = new TH1F(("back_coreEMTopo_SoftClusEM_"+tag).c_str(), "", 40, 0, 200);
        // TH1F* sig_h_coreEMPFlow_SoftClus_MET   = new TH1F(("sig_coreEMPFlow_SoftClus_"  +tag).c_str(), "", 20, 0, 100);
        // TH1F* back_h_coreEMPFlow_SoftClus_MET  = new TH1F(("back_coreEMPFlow_SoftClus_" +tag).c_str(), "", 20, 0, 100);
        // TH1F* sig_h_coreEMPFlow_PVSoftTrk_MET  = new TH1F(("sig_coreEMPFlow_PVSoftTrk_" +tag).c_str(), "", 10, 0, 50);
        // TH1F* back_h_coreEMPFlow_PVSoftTrk_MET = new TH1F(("back_coreEMPFlow_PVSoftTrk_"+tag).c_str(), "", 10, 0, 50);

        // Weighted background histograms for rate plots
        TH1F* back_hw_TotalMET = new TH1F(("back_hw_TotalMET_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET     = new TH1F(("back_hw_gMET_"    +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_NC  = new TH1F(("back_hw_gMET_NC_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_Rms = new TH1F(("back_hw_gMET_Rms_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_jMET     = new TH1F(("back_hw_jMET_"    +tag).c_str(), "", nMETBins, metBinEdges);
        back_hw_jMET->Sumw2();
        TH1F* back_hw_JetMET   = new TH1F(("back_hw_JetMET_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_TowerMET = new TH1F(("back_hw_TowerMET_"+tag).c_str(), "", nMETBins, metBinEdges);
        // H_T = sum of jet E_T, for the H_T rate curve. Its own weighted histogram rather than
        // back_h_SumJetET, which the shape overlays normalize to unit area in place.
        TH1F* back_hw_SumJetET = new TH1F(("back_hw_SumJetET_"+tag).c_str(), "", 52, 0, 1040);
        back_hw_SumJetET->Sumw2();

        // Resimulated gFEX MET distribution histograms and weighted background histograms
        TH1F* sig_h_gMET_JwoJAOD  = new TH1F(("sig_h_gMET_JwoJAOD_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_gMET_JwoJAOD = new TH1F(("back_h_gMET_JwoJAOD_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_gMET_NCAOD    = new TH1F(("sig_h_gMET_NCAOD_"   +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_gMET_NCAOD   = new TH1F(("back_h_gMET_NCAOD_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_gMET_RmsAOD   = new TH1F(("sig_h_gMET_RmsAOD_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_h_gMET_RmsAOD  = new TH1F(("back_h_gMET_RmsAOD_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_JwoJAOD = new TH1F(("back_hw_gMET_JwoJAOD_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_NCAOD   = new TH1F(("back_hw_gMET_NCAOD_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_RmsAOD  = new TH1F(("back_hw_gMET_RmsAOD_" +tag).c_str(), "", nMETBins, metBinEdges);
        back_hw_gMET_JwoJAOD->Sumw2(); back_hw_gMET_NCAOD->Sumw2(); back_hw_gMET_RmsAOD->Sumw2();

        // Binomial per-crossing counterparts of the weighted background histograms above: same
        // events, same fills, but each weight multiplied by the truth-jet binomial correction
        // c(E_T^lead truth jet) built before the loop. Their cumulative is a per-crossing rate
        // rather than a collision rate, and they are deliberately left out of the rate
        // normalization further down (see the note there).
        TH1F* back_hwBin_TotalMET = new TH1F(("back_hwBin_TotalMET_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hwBin_gMET     = new TH1F(("back_hwBin_gMET_"    +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hwBin_gMET_NC  = new TH1F(("back_hwBin_gMET_NC_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hwBin_gMET_Rms = new TH1F(("back_hwBin_gMET_Rms_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hwBin_jMET     = new TH1F(("back_hwBin_jMET_"    +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hwBin_JetMET   = new TH1F(("back_hwBin_JetMET_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hwBin_TowerMET = new TH1F(("back_hwBin_TowerMET_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hwBin_gMET_JwoJAOD = new TH1F(("back_hwBin_gMET_JwoJAOD_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hwBin_gMET_NCAOD   = new TH1F(("back_hwBin_gMET_NCAOD_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hwBin_gMET_RmsAOD  = new TH1F(("back_hwBin_gMET_RmsAOD_" +tag).c_str(), "", nMETBins, metBinEdges);
        back_hwBin_TotalMET->Sumw2(); back_hwBin_gMET->Sumw2();
        back_hwBin_gMET_NC->Sumw2();  back_hwBin_gMET_Rms->Sumw2();
        back_hwBin_jMET->Sumw2();     back_hwBin_JetMET->Sumw2();
        back_hwBin_TowerMET->Sumw2();
        back_hwBin_gMET_JwoJAOD->Sumw2(); back_hwBin_gMET_NCAOD->Sumw2(); back_hwBin_gMET_RmsAOD->Sumw2();

        // Weighted background histograms, JZ0 only with no HSTP requirement (rate comparison
        // against the all-JZ + HSTP-filtered curves). Filled before the HSTP filter below.
        TH1F* back_hw_TotalMET_JZ0 = new TH1F(("back_hw_TotalMET_JZ0_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_JZ0     = new TH1F(("back_hw_gMET_JZ0_"    +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_NC_JZ0  = new TH1F(("back_hw_gMET_NC_JZ0_" +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_Rms_JZ0 = new TH1F(("back_hw_gMET_Rms_JZ0_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_jMET_JZ0     = new TH1F(("back_hw_jMET_JZ0_"    +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_JetMET_JZ0   = new TH1F(("back_hw_JetMET_JZ0_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_TowerMET_JZ0 = new TH1F(("back_hw_TowerMET_JZ0_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_JwoJAOD_JZ0 = new TH1F(("back_hw_gMET_JwoJAOD_JZ0_"+tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_NCAOD_JZ0   = new TH1F(("back_hw_gMET_NCAOD_JZ0_"  +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* back_hw_gMET_RmsAOD_JZ0  = new TH1F(("back_hw_gMET_RmsAOD_JZ0_" +tag).c_str(), "", nMETBins, metBinEdges);
        back_hw_TotalMET_JZ0->Sumw2(); back_hw_gMET_JZ0->Sumw2();
        back_hw_gMET_NC_JZ0->Sumw2();  back_hw_gMET_Rms_JZ0->Sumw2();
        back_hw_jMET_JZ0->Sumw2();     back_hw_JetMET_JZ0->Sumw2();
        back_hw_TowerMET_JZ0->Sumw2();
        back_hw_gMET_JwoJAOD_JZ0->Sumw2(); back_hw_gMET_NCAOD_JZ0->Sumw2(); back_hw_gMET_RmsAOD_JZ0->Sumw2();

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
        TH1F* h_turnOn_num_gMET_60kHz    = new TH1F(("h_turnOn_num_gMET_60kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_JetMET_20kHz  = new TH1F(("h_turnOn_num_JetMET_20kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_JetMET_40kHz  = new TH1F(("h_turnOn_num_JetMET_40kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_JetMET_80kHz  = new TH1F(("h_turnOn_num_JetMET_80kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_JetMET_60kHz  = new TH1F(("h_turnOn_num_JetMET_60kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TowerMET_20kHz = new TH1F(("h_turnOn_num_TowerMET_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TowerMET_40kHz = new TH1F(("h_turnOn_num_TowerMET_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TowerMET_80kHz = new TH1F(("h_turnOn_num_TowerMET_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TowerMET_60kHz = new TH1F(("h_turnOn_num_TowerMET_60kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TotalMET_20kHz = new TH1F(("h_turnOn_num_TotalMET_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TotalMET_40kHz = new TH1F(("h_turnOn_num_TotalMET_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TotalMET_80kHz = new TH1F(("h_turnOn_num_TotalMET_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_TotalMET_60kHz = new TH1F(("h_turnOn_num_TotalMET_60kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        // gFEX NoiseCut and Rms individual turn-on numerators
        TH1F* h_turnOn_num_gMET_NC_20kHz  = new TH1F(("h_turnOn_num_gMET_NC_20kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NC_40kHz  = new TH1F(("h_turnOn_num_gMET_NC_40kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NC_80kHz  = new TH1F(("h_turnOn_num_gMET_NC_80kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NC_60kHz  = new TH1F(("h_turnOn_num_gMET_NC_60kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_Rms_20kHz = new TH1F(("h_turnOn_num_gMET_Rms_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_Rms_40kHz = new TH1F(("h_turnOn_num_gMET_Rms_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_Rms_80kHz = new TH1F(("h_turnOn_num_gMET_Rms_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_Rms_60kHz = new TH1F(("h_turnOn_num_gMET_Rms_60kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        // jFEX MET turn-on numerators
        TH1F* h_turnOn_num_jMET_20kHz    = new TH1F(("h_turnOn_num_jMET_20kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_jMET_40kHz    = new TH1F(("h_turnOn_num_jMET_40kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_jMET_80kHz    = new TH1F(("h_turnOn_num_jMET_80kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_jMET_60kHz    = new TH1F(("h_turnOn_num_jMET_60kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        // Resimulated gFEX MET turn-on numerators
        TH1F* h_turnOn_num_gMET_JwoJAOD_20kHz = new TH1F(("h_turnOn_num_gMET_JwoJAOD_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_JwoJAOD_40kHz = new TH1F(("h_turnOn_num_gMET_JwoJAOD_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_JwoJAOD_80kHz = new TH1F(("h_turnOn_num_gMET_JwoJAOD_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_JwoJAOD_60kHz = new TH1F(("h_turnOn_num_gMET_JwoJAOD_60kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NCAOD_20kHz   = new TH1F(("h_turnOn_num_gMET_NCAOD_20kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NCAOD_40kHz   = new TH1F(("h_turnOn_num_gMET_NCAOD_40kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NCAOD_80kHz   = new TH1F(("h_turnOn_num_gMET_NCAOD_80kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_NCAOD_60kHz   = new TH1F(("h_turnOn_num_gMET_NCAOD_60kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_RmsAOD_20kHz  = new TH1F(("h_turnOn_num_gMET_RmsAOD_20kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_RmsAOD_40kHz  = new TH1F(("h_turnOn_num_gMET_RmsAOD_40kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_RmsAOD_80kHz  = new TH1F(("h_turnOn_num_gMET_RmsAOD_80kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_gMET_RmsAOD_60kHz  = new TH1F(("h_turnOn_num_gMET_RmsAOD_60kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        // Combined gFEX+GEP turn-on numerators (6 combos × 3 rate targets = 18)
        // Declared here as nullptr; created after the 2D scan finds the best thresholds
        TH1F* h_turnOn_num_combo_JwoJ_Jet_20kHz   = new TH1F(("h_turnOn_num_combo_JwoJ_Jet_20kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Jet_40kHz   = new TH1F(("h_turnOn_num_combo_JwoJ_Jet_40kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Jet_80kHz   = new TH1F(("h_turnOn_num_combo_JwoJ_Jet_80kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Jet_60kHz   = new TH1F(("h_turnOn_num_combo_JwoJ_Jet_60kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Tower_20kHz = new TH1F(("h_turnOn_num_combo_JwoJ_Tower_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Tower_40kHz = new TH1F(("h_turnOn_num_combo_JwoJ_Tower_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Tower_80kHz = new TH1F(("h_turnOn_num_combo_JwoJ_Tower_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_JwoJ_Tower_60kHz = new TH1F(("h_turnOn_num_combo_JwoJ_Tower_60kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Jet_20kHz     = new TH1F(("h_turnOn_num_combo_NC_Jet_20kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Jet_40kHz     = new TH1F(("h_turnOn_num_combo_NC_Jet_40kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Jet_80kHz     = new TH1F(("h_turnOn_num_combo_NC_Jet_80kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Jet_60kHz     = new TH1F(("h_turnOn_num_combo_NC_Jet_60kHz_"   +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Tower_20kHz   = new TH1F(("h_turnOn_num_combo_NC_Tower_20kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Tower_40kHz   = new TH1F(("h_turnOn_num_combo_NC_Tower_40kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Tower_80kHz   = new TH1F(("h_turnOn_num_combo_NC_Tower_80kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_NC_Tower_60kHz   = new TH1F(("h_turnOn_num_combo_NC_Tower_60kHz_" +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Jet_20kHz    = new TH1F(("h_turnOn_num_combo_Rms_Jet_20kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Jet_40kHz    = new TH1F(("h_turnOn_num_combo_Rms_Jet_40kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Jet_80kHz    = new TH1F(("h_turnOn_num_combo_Rms_Jet_80kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Jet_60kHz    = new TH1F(("h_turnOn_num_combo_Rms_Jet_60kHz_"  +tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Tower_20kHz  = new TH1F(("h_turnOn_num_combo_Rms_Tower_20kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Tower_40kHz  = new TH1F(("h_turnOn_num_combo_Rms_Tower_40kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Tower_80kHz  = new TH1F(("h_turnOn_num_combo_Rms_Tower_80kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOn_num_combo_Rms_Tower_60kHz  = new TH1F(("h_turnOn_num_combo_Rms_Tower_60kHz_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);

        // Z->mumu turn-on histograms vs dimuon p_{T}: one denominator plus a numerator for every
        // (MET type, rate point). The dimuon system is the well-measured proxy for the invisible
        // momentum here — muons are near-invisible to the calorimeter, so a Z->mumu event looks
        // like an event with MET = p_{T}^{mumu} to every MET algorithm in the table. The binning
        // is the same as the truth-MET turn-ons so the two can be read side by side.
        TH1F* h_turnOnMu_denom = new TH1F(("h_turnOnMu_denom_"+tag).c_str(), "", nTurnOnBins, turnOnBinEdges);
        TH1F* h_turnOnMu_num[nMETTypes][nMuRates];
        for (int iA = 0; iA < nMETTypes; ++iA)
            for (int iR = 0; iR < nMuRates; ++iR)
                h_turnOnMu_num[iA][iR] = new TH1F(
                    (std::string("h_turnOnMu_num_") + metTypeShort[iA] + "_" + muRateNames[iR] + "_" + tag).c_str(),
                    "", nTurnOnBins, turnOnBinEdges);
        // Coarse dimuon p_{T} spectrum, drawn as the shaded distribution under the turn-on curves
        // (same role sig_h_metTruthNonInt_coarse plays for the truth-MET turn-ons).
        TH1F* sig_h_dimuonPt_coarse = new TH1F(("sig_h_dimuonPt_coarse_"+tag).c_str(), "", 30, 0, 600);
        TH1F* sig_h_dimuonPt        = new TH1F(("sig_h_dimuonPt_"       +tag).c_str(), "", nMETBins, metBinEdges);
        TH1F* sig_h_dimuonMass      = new TH1F(("sig_h_dimuonMass_"     +tag).c_str(), "", 60, 60, 120);

        // Background <MET> vs jet multiplicity, one profile per MET type. Weighted with the rate
        // weights, so each point is the rate-weighted mean MET of the crossings with that many
        // jets — the same population the rate curves are built from.
        TProfile* back_prof_METvsNJets[nMETTypes];
        for (int iA = 0; iA < nMETTypes; ++iA) {
            back_prof_METvsNJets[iA] = new TProfile(
                (std::string("back_prof_METvsNJets_") + metTypeShort[iA] + "_" + tag).c_str(),
                "", nNJetBins, 0.0, nJetAxisMax);
            back_prof_METvsNJets[iA]->SetDirectory(0);
        }

        // Reconstructed primary vertices per background event. Filled both ways: raw counts
        // give the sample's vertex multiplicity, the weighted version gives the vertex
        // multiplicity of the background rate (dominated by whichever slices carry the weight).
        TH1F* back_h_nPrimaryVertices    = new TH1F(("back_h_nPrimaryVertices_"   +tag).c_str(),
                                                    "Reconstructed primary vertices;N_{primary vertices};Events",
                                                    35, 20, 160);
        TH1F* back_hw_nPrimaryVertices   = new TH1F(("back_hw_nPrimaryVertices_"  +tag).c_str(),
                                                    "Reconstructed primary vertices;N_{primary vertices};Weighted events",
                                                    35, 20, 160);
        back_hw_nPrimaryVertices->Sumw2();

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
        //TH1F* sig_sel_dPhi_GEP_gFEX[nSel80],   *back_sel_dPhi_GEP_gFEX[nSel80];
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
            //sig_sel_dPhi_GEP_gFEX[iSel]   = new TH1F(("sig_sel_dPhi_GEP_gFEX_" +st).c_str(), "", 32, 0, M_PI);
            //back_sel_dPhi_GEP_gFEX[iSel]  = new TH1F(("back_sel_dPhi_GEP_gFEX_"+st).c_str(), "", 32, 0, M_PI);
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
            back_sel_gMETsig[iSel]->Sumw2();        //back_sel_dPhi_GEP_gFEX[iSel]->Sumw2();
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
            gNomJwoJSig->GetEntry(iEvt);   // nominal gFEX (resim when available, else AOD)
            gNomNCSig->GetEntry(iEvt);
            gNomRmsSig->GetEntry(iEvt);
            metTruthTreeSig->GetEntry(iEvt);
            // metTruthNonIntX/Y already points in the standard MET direction (sum of NonInt particle
            // momenta = direction of neutrinos). The DAOD note's constraint is consistent with
            // momentum conservation, with the Int/IntOut/IntMuon terms (not NonInt) carrying the sign flip.
            // coreEMTopoTreeSig->GetEntry(iEvt);
            // coreEMPFlowTreeSig->GetEntry(iEvt);

            sig_h_TotalMET->Fill(clampVal(sig_h_TotalMET, sig_TotalMET));
            sig_h_TotalMETX->Fill(sig_TotalMETX);
            sig_h_TotalMETY->Fill(sig_TotalMETY);
            sig_h_TowerMet->Fill(clampVal(sig_h_TowerMet, sig_TowerMet));
            sig_h_JetMet->Fill(clampVal(sig_h_JetMet, sig_JetMet));
            sig_h_SumET->Fill(sig_SumET);
            if (hasSumJetET)   sig_h_SumJetET->Fill(sig_SumJetET);
            if (hasSumTowerET) sig_h_SumTowerET->Fill(sig_SumTowerET);
            //std::cout << "sig_gMET: " << sig_gMET << "\n";
            sig_h_gMET->Fill(clampVal(sig_h_gMET, sig_gMET));
            sig_h_gMET_NC->Fill(clampVal(sig_h_gMET_NC, sig_gMET_NC));
            sig_h_gMET_Rms->Fill(clampVal(sig_h_gMET_Rms, sig_gMET_Rms));
            jFexMETTreeSig->GetEntry(iEvt);
            sig_h_jMET->Fill(clampVal(sig_h_jMET, sig_jMET));
            //sig_h_gMETX->Fill(sig_gMETX); sig_h_gMETY->Fill(sig_gMETY);
            //sig_h_gMETX_NC->Fill(sig_gMETX_NC); sig_h_gMETY_NC->Fill(sig_gMETY_NC);
           // sig_h_gMETX_Rms->Fill(sig_gMETX_Rms); sig_h_gMETY_Rms->Fill(sig_gMETY_Rms);
            if (hasGFexSimMET) {
                gFexMETTreeSig->GetEntry(iEvt);          // AOD gFEX copies for the AOD-vs-Sim block
                gFexMETNoiseCutTreeSig->GetEntry(iEvt);
                gFexMETRmsTreeSig->GetEntry(iEvt);
                sig_h_gMET_JwoJAOD->Fill(clampVal(sig_h_gMET_JwoJAOD, sig_gMET_JwoJAOD));
                sig_h_gMET_NCAOD->Fill(clampVal(sig_h_gMET_NCAOD,     sig_gMET_NCAOD));
                sig_h_gMET_RmsAOD->Fill(clampVal(sig_h_gMET_RmsAOD,   sig_gMET_RmsAOD));
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
            // sig_h_coreEMTopo_SoftClus_MET->Fill(clampVal(sig_h_coreEMTopo_SoftClus_MET, sig_coreEMTopo_SoftClus_MET));
            // sig_h_coreEMTopo_PVSoftTrk_MET->Fill(clampVal(sig_h_coreEMTopo_PVSoftTrk_MET, sig_coreEMTopo_PVSoftTrk_MET));
            // sig_h_coreEMTopo_SoftClusEM_MET->Fill(clampVal(sig_h_coreEMTopo_SoftClusEM_MET, sig_coreEMTopo_SoftClusEM_MET));
            // sig_h_coreEMPFlow_SoftClus_MET->Fill(clampVal(sig_h_coreEMPFlow_SoftClus_MET, sig_coreEMPFlow_SoftClus_MET));
            // sig_h_coreEMPFlow_PVSoftTrk_MET->Fill(clampVal(sig_h_coreEMPFlow_PVSoftTrk_MET, sig_coreEMPFlow_PVSoftTrk_MET));
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

        // --- Per-JZ-slice rate bookkeeping (printed after the background loop) ---
        // Diagnoses where the background rate normalization ends up: how much weighted rate
        // each slice carries, how much of it the HSTP filter removes, and how large the
        // single-event weights are (a JZ0 event can be worth O(kHz) on its own).
        const bool printJZRateDiagnostics = true;
        long   dbgNEvents_jz[nJZSlices_]     = {0};
        long   dbgNEventsHSTP_jz[nJZSlices_] = {0};
        double dbgSumW_jz[nJZSlices_]        = {0.0};
        double dbgSumWHSTP_jz[nJZSlices_]    = {0.0};
        double dbgMaxW_jz[nJZSlices_]        = {0.0};
        long   dbgNoSliceEvents = 0;   // events whose sampleJZSlice is outside [0, nJZSlices_)
        long   dbgNoWeightEvents = 0;  // events with an empty eventWeights vector (see below)

        // --- Background event loop ---
        unsigned int nBack = metTreeBack->GetEntries();
        std::cout << "  Background events: " << nBack << "\n";

        // --- Truth-jet binomial correction, built before the main loop ---
        // The binomial pileup factorization splits "the crossing fired" into mu independent
        // collisions. That is only meaningful for an object belonging to a single collision — a
        // jet — not for MET, which is a crossing-level vector sum. So the correction is derived
        // in truth-jet space, where it is valid, and carried over to the MET spectra as a
        // per-event weight (largeRJetAnalysisAndRates.C does the jet-space part of this).
        //
        //   R_coll(x)  cumulative weighted rate of the leading AntiKt4TruthDressedWZ jet above x,
        //              all JZ slices, before the HSTP filter. HERNTupler runs in luminosity mode
        //              (useRateNormalization = false), so the weights already carry the physical
        //              collision rate sigma x L_inst and no RateModeToCollisionRateScale applies.
        //   p(x)       = R_coll(x) / (f_BX x mu), the per-collision pass probability.
        //   R_BX(x)    = f_BX (1 - (1-p)^mu), the per-crossing rate — MakeBinomialCrossingRateHist.
        //   c(x)       = R_BX(x) / R_coll(x), in (0, 1] and rising to 1 where p << 1.
        //
        // Each background event is then filled into back_hwBin_* with weight w x c(x_event),
        // x_event being that event's leading truth jet E_T. This treats an event as if it sat at
        // a threshold equal to its own leading truth jet E_T: events with little jet activity —
        // which dominate the low-MET, high-rate end where crossings saturate — are suppressed
        // most, while rare high-E_T events keep c ~ 1. That is a re-weighting of the event
        // population by saturation, and is deliberately NOT the same thing as applying the
        // binomial formula to the MET rate curve directly.
        // mu is read from THIS file pair's paths rather than from gPileup: SetPileupFromPath does
        // not run until the plotting section further down, so gPileup still holds the previous
        // file's value at this point in the loop (and the initial 200 on the first iteration).
        const double binomialPileup =
            IsPU140Path(backgroundFiles[fileIt].first + " " + backgroundFiles[fileIt].second) ? 140.0 : 200.0;
        TH1F* back_h_leadTruthWZ_Et = nullptr;   // as filled: collision rate vs truth jet E_T
        TH1* h_binomialCorr         = nullptr;   // c(x), looked up per event below
        if (hasTruthAntiKt4WZDressed) {
            back_h_leadTruthWZ_Et = new TH1F(("back_h_leadTruthWZ_Et_"+tag).c_str(),
                                             "Leading truth HS jet E_{T} (JZ0-9, no HSTP);"
                                             "E_{T} [GeV];Rate [Hz]", nMETBins, metBinEdges);
            back_h_leadTruthWZ_Et->Sumw2();
            TStopwatch swPre; swPre.Start();
            if (printIOProgress)
                std::cout << "  [binomial] truth-jet pre-pass over " << nBack << " events...\n" << std::flush;
            for (unsigned int iEvt = 0; iEvt < nBack; iEvt++) {
                if (printIOProgress && iEvt > 0 && iEvt % progressEvery == 0) {
                    std::cout << "  [binomial] pre-pass " << iEvt << "/" << nBack
                              << " (" << swPre.RealTime() << " s)\n" << std::flush;
                    swPre.Continue();
                }
                eventInfoTreeBack->GetEntry(iEvt);
                if (!eventWeightsValuesBack || eventWeightsValuesBack->empty()) continue;
                const double wTruth = eventWeightsValuesBack->at(0);
                // Events with no truth jet fill 0, so bin 1 of the cumulative stays the total
                // event rate — the same convention largeRJetAnalysisAndRates.C uses.
                double etLeadTruthWZ = 0.0;
                if (leadingTruthAntiKt4WZDressedJetsTreeBack) {
                    leadingTruthAntiKt4WZDressedJetsTreeBack->GetEntry(iEvt);
                    if (leadingTruthAntiKt4WZDressedJetsEtValuesBack &&
                        !leadingTruthAntiKt4WZDressedJetsEtValuesBack->empty())
                        etLeadTruthWZ = leadingTruthAntiKt4WZDressedJetsEtValuesBack->at(0);
                } else {
                    truthAntiKt4WZDressedJetsTreeBack->GetEntry(iEvt);
                    if (truthAntiKt4WZDressedJetsEtValuesBack &&
                        !truthAntiKt4WZDressedJetsEtValuesBack->empty())
                        etLeadTruthWZ = *std::max_element(truthAntiKt4WZDressedJetsEtValuesBack->begin(),
                                                         truthAntiKt4WZDressedJetsEtValuesBack->end());
                }
                back_h_leadTruthWZ_Et->Fill(clampVal(back_h_leadTruthWZ_Et, etLeadTruthWZ), wTruth);
            }

            TH1* hTruthColl = MakeCumulativeRateHist(back_h_leadTruthWZ_Et,
                                                     ("truthColl_"+tag).c_str());
            TH1* hTruthBX   = MakeBinomialCrossingRateHist(hTruthColl, ("truthBX_"+tag).c_str(),
                                                           /*collisionRateScale=*/1.0, binomialPileup);
            h_binomialCorr = (TH1*)hTruthBX->Clone(("binomialCorr_"+tag).c_str());
            h_binomialCorr->SetDirectory(nullptr);
            for (int ib = 1; ib <= h_binomialCorr->GetNbinsX(); ++ib) {
                const double rColl = hTruthColl->GetBinContent(ib);
                // Empty bins carry no events to weight; 1.0 leaves any that appear untouched.
                h_binomialCorr->SetBinContent(ib, rColl > 0.0 ? hTruthBX->GetBinContent(ib) / rColl : 1.0);
                h_binomialCorr->SetBinError(ib, 0.0);
            }
            printf("  [binomial] mu = %.0f, f_BX = %.1f MHz, truth-jet collision rate at E_{T} > 0"
                   " = %.4g GHz\n", binomialPileup, kCrossingRateHz / 1e6,
                   hTruthColl->GetBinContent(1) / 1e9);
            printf("  [binomial] c(E_T) at 0 / 50 / 100 / 200 GeV = %.4g / %.4g / %.4g / %.4g\n",
                   h_binomialCorr->GetBinContent(h_binomialCorr->FindBin(0.0)),
                   h_binomialCorr->GetBinContent(h_binomialCorr->FindBin(50.0)),
                   h_binomialCorr->GetBinContent(h_binomialCorr->FindBin(100.0)),
                   h_binomialCorr->GetBinContent(h_binomialCorr->FindBin(200.0)));
            delete hTruthColl; delete hTruthBX;
        } else {
            std::cout << "  Truth AntiKt4 WZ-dressed jets not found — binomial per-crossing"
                      << " curves skipped\n";
        }

        // Per-event lookup of c(x): bin of this event's leading truth jet E_T, overflow clamped
        // to the last bin. Returns 1.0 when the correction could not be built, which leaves the
        // back_hwBin_* histograms identical to the as-filled ones.
        auto binomialWeight = [&](double etLeadTruthWZ) -> double {
            if (!h_binomialCorr) return 1.0;
            int ib = h_binomialCorr->FindBin(etLeadTruthWZ);
            if (ib < 1) ib = 1;
            if (ib > h_binomialCorr->GetNbinsX()) ib = h_binomialCorr->GetNbinsX();
            return h_binomialCorr->GetBinContent(ib);
        };
        TStopwatch swBack; swBack.Start();
        for (unsigned int iEvt = 0; iEvt < nBack; iEvt++) {
            if (printIOProgress && iEvt > 0 && iEvt % progressEvery == 0) {
                std::cout << "  background loop " << iEvt << "/" << nBack
                          << " (" << swBack.RealTime() << " s)\n" << std::flush;
                swBack.Continue();
            }
            metTreeBack->GetEntry(iEvt);
            gNomJwoJBack->GetEntry(iEvt);   // nominal gFEX (resim when available, else AOD)
            gNomNCBack->GetEntry(iEvt);
            gNomRmsBack->GetEntry(iEvt);
            metTruthTreeBack->GetEntry(iEvt);
            // coreEMTopoTreeBack->GetEntry(iEvt);
            // coreEMPFlowTreeBack->GetEntry(iEvt);
            eventInfoTreeBack->GetEntry(iEvt);

            // The two blocks below run BEFORE the HSTP filter, so unlike the rest of the loop
            // they also see events the filter rejects. Some of those carry an empty
            // eventWeights vector, so read it defensively rather than with at(0).
            const bool haveWeight = (eventWeightsValuesBack && !eventWeightsValuesBack->empty());
            const double wPre = haveWeight ? eventWeightsValuesBack->at(0) : 0.0;
            if (!haveWeight) dbgNoWeightEvents++;

            // Rate bookkeeping — counted before the HSTP filter so the filter's effect is visible.
            if (haveWeight && sampleJZSliceBack >= 0 && sampleJZSliceBack < (int)nJZSlices_) {
                double wDbg = wPre;
                dbgNEvents_jz[sampleJZSliceBack]++;
                dbgSumW_jz[sampleJZSliceBack] += wDbg;
                if (wDbg > dbgMaxW_jz[sampleJZSliceBack]) dbgMaxW_jz[sampleJZSliceBack] = wDbg;
                if (!applyHSTPFilter || passHSTPValuesBack) {
                    dbgNEventsHSTP_jz[sampleJZSliceBack]++;
                    dbgSumWHSTP_jz[sampleJZSliceBack] += wDbg;
                }
            } else if (haveWeight) {
                dbgNoSliceEvents++;
            }

            // JZ0-only fills (no HSTP requirement) for rate comparison — must run before the
            // HSTP filter below. jFEX / AOD gFEX entries are loaded here since their GetEntry
            // calls sit after the filter; repeating them below is harmless.
            if (haveWeight && sampleJZSliceBack == 0) {
                double wJZ0 = wPre;
                jFexMETTreeBack->GetEntry(iEvt);
                back_hw_TotalMET_JZ0->Fill(clampVal(back_hw_TotalMET_JZ0, back_TotalMET),  wJZ0);
                back_hw_gMET_JZ0->Fill(clampVal(back_hw_gMET_JZ0,         back_gMET),      wJZ0);
                back_hw_gMET_NC_JZ0->Fill(clampVal(back_hw_gMET_NC_JZ0,   back_gMET_NC),   wJZ0);
                back_hw_gMET_Rms_JZ0->Fill(clampVal(back_hw_gMET_Rms_JZ0, back_gMET_Rms),  wJZ0);
                back_hw_jMET_JZ0->Fill(clampVal(back_hw_jMET_JZ0,         back_jMET),      wJZ0);
                back_hw_JetMET_JZ0->Fill(clampVal(back_hw_JetMET_JZ0,     back_JetMet),    wJZ0);
                back_hw_TowerMET_JZ0->Fill(clampVal(back_hw_TowerMET_JZ0, back_TowerMet),  wJZ0);
                if (hasGFexSimMET) {
                    gFexMETTreeBack->GetEntry(iEvt);
                    gFexMETNoiseCutTreeBack->GetEntry(iEvt);
                    gFexMETRmsTreeBack->GetEntry(iEvt);
                    back_hw_gMET_JwoJAOD_JZ0->Fill(clampVal(back_hw_gMET_JwoJAOD_JZ0, back_gMET_JwoJAOD), wJZ0);
                    back_hw_gMET_NCAOD_JZ0->Fill(clampVal(back_hw_gMET_NCAOD_JZ0,     back_gMET_NCAOD),   wJZ0);
                    back_hw_gMET_RmsAOD_JZ0->Fill(clampVal(back_hw_gMET_RmsAOD_JZ0,   back_gMET_RmsAOD),  wJZ0);
                }
            }

            if (applyHSTPFilter && !passHSTPValuesBack) continue;
            double w = eventWeightsValuesBack->at(0);
            if (hasNPrimaryVertices && nPrimaryVerticesBack >= 0) {
                back_h_nPrimaryVertices->Fill(clampVal(back_h_nPrimaryVertices, nPrimaryVerticesBack));
                back_hw_nPrimaryVertices->Fill(clampVal(back_hw_nPrimaryVertices, nPrimaryVerticesBack), w);
            }
            back_h_TotalMET->Fill(clampVal(back_h_TotalMET, back_TotalMET), w);
            back_h_TotalMETX->Fill(back_TotalMETX, w);
            back_h_TotalMETY->Fill(back_TotalMETY, w);
            back_h_TowerMet->Fill(clampVal(back_h_TowerMet, back_TowerMet), w);
            back_h_JetMet->Fill(clampVal(back_h_JetMet, back_JetMet), w);
            back_h_SumET->Fill(back_SumET, w);
            if (hasSumJetET)   back_h_SumJetET->Fill(back_SumJetET, w);
            if (hasSumTowerET) back_h_SumTowerET->Fill(back_SumTowerET, w);
            if (hasSumJetET)
                back_hw_SumJetET->Fill(clampVal(back_hw_SumJetET, back_SumJetET), w);
            back_h_gMET->Fill(clampVal(back_h_gMET, back_gMET), w);
            back_h_gMET_NC->Fill(clampVal(back_h_gMET_NC, back_gMET_NC), w);
            back_h_gMET_Rms->Fill(clampVal(back_h_gMET_Rms, back_gMET_Rms), w);
            jFexMETTreeBack->GetEntry(iEvt);
            back_h_jMET->Fill(clampVal(back_h_jMET, back_jMET), w);
            //back_h_gMETX->Fill(back_gMETX, w); back_h_gMETY->Fill(back_gMETY, w);
            //back_h_gMETX_NC->Fill(back_gMETX_NC, w); back_h_gMETY_NC->Fill(back_gMETY_NC, w);
            //back_h_gMETX_Rms->Fill(back_gMETX_Rms, w); back_h_gMETY_Rms->Fill(back_gMETY_Rms, w);
            back_h_metTruthNonIntX->Fill(back_metTruthNonIntX, w);
            back_h_metTruthNonIntY->Fill(back_metTruthNonIntY, w);
            back_h_JetMetX->Fill(back_JetMetX, w); back_h_JetMetY->Fill(back_JetMetY, w);
            back_h_TowerMetX->Fill(back_TowerMetX, w); back_h_TowerMetY->Fill(back_TowerMetY, w);
            back_h_metTruthNonInt->Fill(clampVal(back_h_metTruthNonInt, back_metTruthNonInt), w);
            back_h_metTruthInt->Fill(clampVal(back_h_metTruthInt, back_metTruthInt), w);
            back_h_metTruthIntOut->Fill(clampVal(back_h_metTruthIntOut, back_metTruthIntOut), w);
            // back_h_coreEMTopo_SoftClus_MET->Fill(clampVal(back_h_coreEMTopo_SoftClus_MET, back_coreEMTopo_SoftClus_MET), w);
            // back_h_coreEMTopo_PVSoftTrk_MET->Fill(clampVal(back_h_coreEMTopo_PVSoftTrk_MET, back_coreEMTopo_PVSoftTrk_MET), w);
            // back_h_coreEMTopo_SoftClusEM_MET->Fill(clampVal(back_h_coreEMTopo_SoftClusEM_MET, back_coreEMTopo_SoftClusEM_MET), w);
            // back_h_coreEMPFlow_SoftClus_MET->Fill(clampVal(back_h_coreEMPFlow_SoftClus_MET, back_coreEMPFlow_SoftClus_MET), w);
            // back_h_coreEMPFlow_PVSoftTrk_MET->Fill(clampVal(back_h_coreEMPFlow_PVSoftTrk_MET, back_coreEMPFlow_PVSoftTrk_MET), w);

            back_hw_TotalMET->Fill(clampVal(back_hw_TotalMET, back_TotalMET), w);
            back_hw_gMET->Fill(clampVal(back_hw_gMET, back_gMET), w);
            back_hw_gMET_NC->Fill(clampVal(back_hw_gMET_NC, back_gMET_NC), w);
            back_hw_gMET_Rms->Fill(clampVal(back_hw_gMET_Rms, back_gMET_Rms), w);
            back_hw_jMET->Fill(clampVal(back_hw_jMET, back_jMET), w);
            back_hw_JetMET->Fill(clampVal(back_hw_JetMET, back_JetMet), w);
            back_hw_TowerMET->Fill(clampVal(back_hw_TowerMET, back_TowerMet), w);

            // Same fills weighted by the truth-jet binomial correction at this event's leading
            // truth jet E_T, giving a per-crossing rather than a collision rate. wBin == w when
            // the correction could not be built, so these stay usable either way.
            double etLeadTruthWZEvt = 0.0;
            if (hasTruthAntiKt4WZDressed) {
                if (leadingTruthAntiKt4WZDressedJetsTreeBack) {
                    leadingTruthAntiKt4WZDressedJetsTreeBack->GetEntry(iEvt);
                    if (leadingTruthAntiKt4WZDressedJetsEtValuesBack &&
                        !leadingTruthAntiKt4WZDressedJetsEtValuesBack->empty())
                        etLeadTruthWZEvt = leadingTruthAntiKt4WZDressedJetsEtValuesBack->at(0);
                } else {
                    truthAntiKt4WZDressedJetsTreeBack->GetEntry(iEvt);
                    if (truthAntiKt4WZDressedJetsEtValuesBack &&
                        !truthAntiKt4WZDressedJetsEtValuesBack->empty())
                        etLeadTruthWZEvt = *std::max_element(truthAntiKt4WZDressedJetsEtValuesBack->begin(),
                                                             truthAntiKt4WZDressedJetsEtValuesBack->end());
                }
            }
            const double wBin = w * binomialWeight(etLeadTruthWZEvt);
            back_hwBin_TotalMET->Fill(clampVal(back_hwBin_TotalMET, back_TotalMET), wBin);
            back_hwBin_gMET->Fill(clampVal(back_hwBin_gMET, back_gMET), wBin);
            back_hwBin_gMET_NC->Fill(clampVal(back_hwBin_gMET_NC, back_gMET_NC), wBin);
            back_hwBin_gMET_Rms->Fill(clampVal(back_hwBin_gMET_Rms, back_gMET_Rms), wBin);
            back_hwBin_jMET->Fill(clampVal(back_hwBin_jMET, back_jMET), wBin);
            back_hwBin_JetMET->Fill(clampVal(back_hwBin_JetMET, back_JetMet), wBin);
            back_hwBin_TowerMET->Fill(clampVal(back_hwBin_TowerMET, back_TowerMet), wBin);
            if (hasGFexSimMET) {
                gFexMETTreeBack->GetEntry(iEvt);          // AOD gFEX copies for the AOD-vs-Sim block
                gFexMETNoiseCutTreeBack->GetEntry(iEvt);
                gFexMETRmsTreeBack->GetEntry(iEvt);
                back_h_gMET_JwoJAOD->Fill(clampVal(back_h_gMET_JwoJAOD, back_gMET_JwoJAOD), w);
                back_h_gMET_NCAOD->Fill(clampVal(back_h_gMET_NCAOD,     back_gMET_NCAOD),   w);
                back_h_gMET_RmsAOD->Fill(clampVal(back_h_gMET_RmsAOD,   back_gMET_RmsAOD),  w);
                back_hw_gMET_JwoJAOD->Fill(clampVal(back_hw_gMET_JwoJAOD, back_gMET_JwoJAOD), w);
                back_hw_gMET_NCAOD->Fill(clampVal(back_hw_gMET_NCAOD,     back_gMET_NCAOD),   w);
                back_hw_gMET_RmsAOD->Fill(clampVal(back_hw_gMET_RmsAOD,   back_gMET_RmsAOD),  w);
                back_hwBin_gMET_JwoJAOD->Fill(clampVal(back_hwBin_gMET_JwoJAOD, back_gMET_JwoJAOD), wBin);
                back_hwBin_gMET_NCAOD->Fill(clampVal(back_hwBin_gMET_NCAOD,     back_gMET_NCAOD),   wBin);
                back_hwBin_gMET_RmsAOD->Fill(clampVal(back_hwBin_gMET_RmsAOD,   back_gMET_RmsAOD),  wBin);
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

        // --- Rate normalization (see the note next to targetTotalRateHz) ---
        // Applied here, once, to the weighted background histograms themselves, so that every
        // rate number derived from them below starts from the same normalization: the rate
        // curves, the 20/40/60/80 kHz thresholds, the rate-vs-efficiency curves and — through
        // those thresholds — the turn-on curves. Each histogram is filled once per surviving
        // event, so its integral is the rate at threshold 0 and also bin 1 of the cumulative.
        //
        // Each histogram is scaled by its OWN integral, so every curve starts at exactly
        // targetTotalRateHz. The factors are identical for all of them except the AOD ones,
        // which are filled only on events that carry resimulated gFEX MET. Any constant already
        // in the weights cancels, so this is unaffected by HERNTupler's useRateNormalization
        // setting — only the as-filled scale moves with that, not the normalized one.
        //
        // NOT scaled here: the JZ0-only histograms (back_hw_*_JZ0), which drawJZ0Comparison
        // normalizes itself and compares against the as-filled all-JZ curve;
        // back_hw_nPrimaryVertices, which is a distribution rather than a rate; and the
        // binomial per-crossing histograms (back_hwBin_*), which are absolute per-crossing rates
        // — rescaling those to targetTotalRateHz at threshold 0 would erase the saturation they
        // exist to show. They are the per-crossing counterpart of the as-filled curve.
        std::map<TH1*, double> appliedRateScale;   // histogram -> factor applied, for undoing it
        auto normalizeRateHist = [&](TH1* h) -> double {
            if (!h) return 1.0;
            const double atZero = h->Integral();
            const double scale  = (normalizeRateToTarget && atZero > 0.0)
                                  ? targetTotalRateHz / atZero : 1.0;
            if (scale != 1.0) h->Scale(scale);
            appliedRateScale[h] = scale;
            return scale;
        };

        const double rateAtZero        = back_hw_gMET->Integral();   // as filled, before scaling
        const double rateToTargetScale = normalizeRateHist(back_hw_gMET);
        const double sc_gMET_NC        = normalizeRateHist(back_hw_gMET_NC);
        const double sc_gMET_Rms       = normalizeRateHist(back_hw_gMET_Rms);
        const double sc_jMET           = normalizeRateHist(back_hw_jMET);
        const double sc_JetMET         = normalizeRateHist(back_hw_JetMET);
        const double sc_TowerMET       = normalizeRateHist(back_hw_TowerMET);
        normalizeRateHist(back_hw_TotalMET);
        if (hasSumJetET) normalizeRateHist(back_hw_SumJetET);   // H_T rate curve
        normalizeRateHist(back_hw_gMET_JwoJAOD);
        normalizeRateHist(back_hw_gMET_NCAOD);
        normalizeRateHist(back_hw_gMET_RmsAOD);
        // 2D histograms behind the combined (gFEX MET, GEP MET) rate-vs-efficiency scans.
        normalizeRateHist(back_hw2_JwoJ_Jet);   normalizeRateHist(back_hw2_JwoJ_Tower);
        normalizeRateHist(back_hw2_NC_Jet);     normalizeRateHist(back_hw2_NC_Tower);
        normalizeRateHist(back_hw2_Rms_Jet);    normalizeRateHist(back_hw2_Rms_Tower);
        // Per-JZ-slice rate curves get their parent MET type's factor rather than one of their
        // own: the slices are meant to add up to the total, and normalizing each slice to
        // targetTotalRateHz would throw away exactly the composition those plots exist to show.
        for (unsigned int jz = 0; jz < nJZSlices_; ++jz) {
            if (back_h_gMET_jz[jz])     back_h_gMET_jz[jz]->Scale(rateToTargetScale);
            if (back_h_gMET_NC_jz[jz])  back_h_gMET_NC_jz[jz]->Scale(sc_gMET_NC);
            if (back_h_gMET_Rms_jz[jz]) back_h_gMET_Rms_jz[jz]->Scale(sc_gMET_Rms);
            if (back_h_jMET_jz[jz])     back_h_jMET_jz[jz]->Scale(sc_jMET);
            if (back_h_JetMET_jz[jz])   back_h_JetMET_jz[jz]->Scale(sc_JetMET);
            if (back_h_TowerMET_jz[jz]) back_h_TowerMET_jz[jz]->Scale(sc_TowerMET);
        }

        if (rateAtZero <= 0.0)
            std::cout << "  [rate] WARNING: weighted background integral is zero — rate plots will be empty\n";
        else if (normalizeRateToTarget)
            std::cout << "  [rate] rate at threshold 0 = " << rateAtZero / 1e3 << " kHz as filled;"
                      << " normalized to " << targetTotalRateHz / 1e6 << " MHz (x" << rateToTargetScale << ")\n";
        else
            std::cout << "  [rate] rate at threshold 0 = " << rateAtZero / 1e3 << " kHz;"
                      << " normalizeRateToTarget = false, using the weights as filled\n";

        if (h_binomialCorr) {
            const double binRateAtZero = back_hwBin_gMET->Integral();
            printf("  [binomial] per-crossing rate at threshold 0 = %.4g MHz (gFEX JwoJ)\n",
                   binRateAtZero / 1e6);
            if (binRateAtZero > kCrossingRateHz)
                printf("  [binomial] WARNING: that is above the %.1f MHz crossing rate — the"
                       " correction is not doing what it should\n", kCrossingRateHz / 1e6);
        }

        // --- Per-JZ-slice rate diagnostics ---
        // The per-event weight is fixed in HERNTupler (xsec x filterEff x targetRate/sigma_ref
        // / sumOfWeights), i.e. the 30 MHz target is normalized over the UNFILTERED sample. The
        // HSTP filter is applied here, afterwards, so the surviving Sum(w) is whatever is left —
        // it is not renormalized back to 30 MHz. This table shows how much rate each slice
        // carries before and after the filter.
        if (printJZRateDiagnostics) {
            double totW = 0.0, totWHSTP = 0.0;
            long   totN = 0,   totNHSTP = 0;
            for (unsigned int jz = 0; jz < nJZSlices_; ++jz) {
                totW += dbgSumW_jz[jz]; totWHSTP += dbgSumWHSTP_jz[jz];
                totN += dbgNEvents_jz[jz]; totNHSTP += dbgNEventsHSTP_jz[jz];
            }
            std::cout << "\n  ===== JZ-slice rate diagnostics (" << labels[fileIt] << ") =====\n"
                      << "  applyHSTPFilter = " << (applyHSTPFilter ? "true" : "false") << "\n"
                      << "   JZ   NEvents   Sum(w) [kHz]   NEvt HSTP   Sum(w) HSTP [kHz]"
                      << "   HSTP keep(w)   mean w [Hz]   max w [Hz]\n";
            for (unsigned int jz = 0; jz < nJZSlices_; ++jz) {
                double keep = dbgSumW_jz[jz] > 0.0 ? dbgSumWHSTP_jz[jz] / dbgSumW_jz[jz] : 0.0;
                double meanW = dbgNEvents_jz[jz] > 0 ? dbgSumW_jz[jz] / dbgNEvents_jz[jz] : 0.0;
                printf("   %2u %9ld %14.4g %11ld %19.4g %14.4f %13.4g %12.4g\n",
                       jz, dbgNEvents_jz[jz], dbgSumW_jz[jz] * 1e-3, dbgNEventsHSTP_jz[jz],
                       dbgSumWHSTP_jz[jz] * 1e-3, keep, meanW, dbgMaxW_jz[jz]);
            }
            printf("   ALL %8ld %14.4g %11ld %19.4g %14.4f\n",
                   totN, totW * 1e-3, totNHSTP, totWHSTP * 1e-3,
                   totW > 0.0 ? totWHSTP / totW : 0.0);
            if (dbgNoSliceEvents > 0)
                std::cout << "  WARNING: " << dbgNoSliceEvents
                          << " background events had sampleJZSlice outside [0," << nJZSlices_
                          << ") and were left out of the table\n";
            if (dbgNoWeightEvents > 0)
                std::cout << "  WARNING: " << dbgNoWeightEvents
                          << " background events had an empty eventWeights vector and were left out"
                          << " of the table and the JZ0 fills. These are only visible here because"
                          << " this block runs before the HSTP filter.\n";
            std::cout << "  Total rate before HSTP = " << totW * 1e-6 << " MHz"
                      << " (HERNTupler targetRate is 30 MHz over the unfiltered sample;"
                      << " a value well below that means only part of each slice was ntupled)\n"
                      << "  Total rate after  HSTP = " << totWHSTP * 1e-3 << " kHz\n";
            // Cross-check against the histograms feeding the rate curves (bin 1 = threshold 0).
            printf("  Rate at threshold 0: all-JZ gFEX JwoJ = %.4g kHz (after the x%.4g rescale),"
                   " JZ0-only no-HSTP = %.4g kHz (never rescaled)\n",
                   back_hw_gMET->Integral() * 1e-3, rateToTargetScale,
                   back_hw_gMET_JZ0->Integral() * 1e-3);
            std::cout << "  =========================================================\n\n";
        }

        // --- Threshold computation at 20 / 40 / 80 kHz ---
        double thr_gMET_20kHz      = findThreshold(back_hw_gMET,     20e3);
        double thr_gMET_40kHz      = findThreshold(back_hw_gMET,     40e3);
        double thr_gMET_80kHz      = findThreshold(back_hw_gMET,     80e3);
        double thr_gMET_60kHz      = findThreshold(back_hw_gMET,     60e3);
        double thr_gMET_NC_20kHz   = findThreshold(back_hw_gMET_NC,  20e3);
        double thr_gMET_NC_40kHz   = findThreshold(back_hw_gMET_NC,  40e3);
        double thr_gMET_NC_80kHz   = findThreshold(back_hw_gMET_NC,  80e3);
        double thr_gMET_NC_60kHz   = findThreshold(back_hw_gMET_NC,  60e3);
        double thr_gMET_Rms_20kHz  = findThreshold(back_hw_gMET_Rms, 20e3);
        double thr_gMET_Rms_40kHz  = findThreshold(back_hw_gMET_Rms, 40e3);
        double thr_gMET_Rms_80kHz  = findThreshold(back_hw_gMET_Rms, 80e3);
        double thr_gMET_Rms_60kHz  = findThreshold(back_hw_gMET_Rms, 60e3);
        double thr_jMET_20kHz      = findThreshold(back_hw_jMET,     20e3);
        double thr_jMET_40kHz      = findThreshold(back_hw_jMET,     40e3);
        double thr_jMET_80kHz      = findThreshold(back_hw_jMET,     80e3);
        double thr_jMET_60kHz      = findThreshold(back_hw_jMET,     60e3);
        double thr_JetMET_20kHz    = findThreshold(back_hw_JetMET,   20e3);
        double thr_JetMET_40kHz    = findThreshold(back_hw_JetMET,   40e3);
        double thr_JetMET_80kHz    = findThreshold(back_hw_JetMET,   80e3);
        double thr_JetMET_60kHz    = findThreshold(back_hw_JetMET,   60e3);
        double thr_TowerMET_20kHz  = findThreshold(back_hw_TowerMET,  20e3);
        double thr_TowerMET_40kHz  = findThreshold(back_hw_TowerMET,  40e3);
        double thr_TowerMET_80kHz  = findThreshold(back_hw_TowerMET,  80e3);
        double thr_TowerMET_60kHz  = findThreshold(back_hw_TowerMET,  60e3);
        double thr_TotalMET_20kHz  = findThreshold(back_hw_TotalMET,  20e3);
        double thr_TotalMET_40kHz  = findThreshold(back_hw_TotalMET,  40e3);
        double thr_TotalMET_80kHz  = findThreshold(back_hw_TotalMET,  80e3);
        double thr_TotalMET_60kHz  = findThreshold(back_hw_TotalMET,  60e3);
        double thr_gMET_JwoJAOD_20kHz = hasGFexSimMET ? findThreshold(back_hw_gMET_JwoJAOD, 20e3) : 0.0;
        double thr_gMET_JwoJAOD_40kHz = hasGFexSimMET ? findThreshold(back_hw_gMET_JwoJAOD, 40e3) : 0.0;
        double thr_gMET_JwoJAOD_80kHz = hasGFexSimMET ? findThreshold(back_hw_gMET_JwoJAOD, 80e3) : 0.0;
        double thr_gMET_JwoJAOD_60kHz = hasGFexSimMET ? findThreshold(back_hw_gMET_JwoJAOD, 60e3) : 0.0;
        double thr_gMET_NCAOD_20kHz   = hasGFexSimMET ? findThreshold(back_hw_gMET_NCAOD,   20e3) : 0.0;
        double thr_gMET_NCAOD_40kHz   = hasGFexSimMET ? findThreshold(back_hw_gMET_NCAOD,   40e3) : 0.0;
        double thr_gMET_NCAOD_80kHz   = hasGFexSimMET ? findThreshold(back_hw_gMET_NCAOD,   80e3) : 0.0;
        double thr_gMET_NCAOD_60kHz   = hasGFexSimMET ? findThreshold(back_hw_gMET_NCAOD,   60e3) : 0.0;
        double thr_gMET_RmsAOD_20kHz  = hasGFexSimMET ? findThreshold(back_hw_gMET_RmsAOD,  20e3) : 0.0;
        double thr_gMET_RmsAOD_40kHz  = hasGFexSimMET ? findThreshold(back_hw_gMET_RmsAOD,  40e3) : 0.0;
        double thr_gMET_RmsAOD_80kHz  = hasGFexSimMET ? findThreshold(back_hw_gMET_RmsAOD,  80e3) : 0.0;
        double thr_gMET_RmsAOD_60kHz  = hasGFexSimMET ? findThreshold(back_hw_gMET_RmsAOD,  60e3) : 0.0;

        // --- Rate-matched thresholds for the Z->mumu dimuon-p_{T} turn-ons ---
        // Same weighted background histograms and same findThreshold as the block above, indexed
        // through the MET-type table so the turn-on loop can stay generic. These are the
        // normalized histograms (the rescale further up has already run), so the 40/60/80 kHz
        // points mean the same thing as everywhere else in this macro.
        TH1F* back_hw_byMETType[nMETTypes] = {
            back_hw_gMET, back_hw_gMET_NC, back_hw_gMET_Rms, back_hw_jMET,
            back_hw_JetMET, back_hw_TowerMET, back_hw_TotalMET
        };
        double thrMu[nMETTypes][nMuRates];
        for (int iA = 0; iA < nMETTypes; ++iA)
            for (int iR = 0; iR < nMuRates; ++iR)
                thrMu[iA][iR] = findThreshold(back_hw_byMETType[iA], muRatesHz[iR]);

        if (hasGFexSimMET)
            std::cout << "  Thresholds [GeV] — gMET JwoJSim 20/40/80 kHz: "
                      << thr_gMET_JwoJAOD_20kHz << " / " << thr_gMET_JwoJAOD_40kHz << " / " << thr_gMET_JwoJAOD_80kHz << "\n"
                      << "                    gMET NCSim   20/40/80 kHz: "
                      << thr_gMET_NCAOD_20kHz   << " / " << thr_gMET_NCAOD_40kHz   << " / " << thr_gMET_NCAOD_80kHz   << "\n"
                      << "                    gMET RmsSim  20/40/80 kHz: "
                      << thr_gMET_RmsAOD_20kHz  << " / " << thr_gMET_RmsAOD_40kHz  << " / " << thr_gMET_RmsAOD_80kHz  << "\n";
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

        // --- Estimated background rate for fixed-threshold trigger items (kHz) ---
        // gXEJWOJ = gFEX XE (JwoJ), jXE = jFEX XE; number = fixed MET threshold [GeV].
        {
            auto rateAtThr_kHz = [](TH1F* hw, double thrGeV) {
                const int binLo = hw->FindFixBin(thrGeV);
                return hw->Integral(binLo, hw->GetNbinsX() + 1) / 1e3; // Hz -> kHz
            };
            std::cout << "  Trigger item rates [kHz]:\n"
                      << "                    gXEJWOJ110: " << rateAtThr_kHz(back_hw_gMET, 110.0) << "\n"
                      << "                    gXEJWOJ100: " << rateAtThr_kHz(back_hw_gMET, 100.0) << "\n"
                      << "                    jXE110:     " << rateAtThr_kHz(back_hw_jMET, 110.0) << "\n"
                      << "                    jXE100:     " << rateAtThr_kHz(back_hw_jMET, 100.0) << "\n";
        }

        // --- Per-JZ-slice contributions to jFEX MET above a fixed threshold ---
        {
            const double jzMETThr = 140.0; // GeV
            double jzWeightedTotal = 0.0;
            double jzWeightedPass[nJZSlices_];
            for (unsigned int jz = 0; jz < nJZSlices_; ++jz) {
                const int binLo = back_h_jMET_jz[jz]->FindFixBin(jzMETThr);
                // include the overflow bin (nbins+1) in the pass integral
                jzWeightedPass[jz] = back_h_jMET_jz[jz]->Integral(binLo, back_h_jMET_jz[jz]->GetNbinsX() + 1);
                jzWeightedTotal += jzWeightedPass[jz];
            }
            std::cout << "  --- h_rates_jfexmet slice contributions at MET > " << jzMETThr << " GeV ---\n";
            for (unsigned int jz = 0; jz < nJZSlices_; ++jz) {
                const double frac = (jzWeightedTotal > 0.0) ? (100.0 * jzWeightedPass[jz] / jzWeightedTotal) : 0.0;
                std::cout << "     JZ" << jz << ": "
                          << std::fixed << std::setw(6) << std::setprecision(2) << frac << "%   (n_events="
                          << (long long)back_h_jMET_jz[jz]->GetEntries() << ")\n"
                          << std::defaultfloat;
            }

            double jzWeightedTotalNC = 0.0;
            double jzWeightedPassNC[nJZSlices_];
            for (unsigned int jz = 0; jz < nJZSlices_; ++jz) {
                const int binLo = back_h_gMET_NC_jz[jz]->FindFixBin(jzMETThr);
                // include the overflow bin (nbins+1) in the pass integral
                jzWeightedPassNC[jz] = back_h_gMET_NC_jz[jz]->Integral(binLo, back_h_gMET_NC_jz[jz]->GetNbinsX() + 1);
                jzWeightedTotalNC += jzWeightedPassNC[jz];
            }
            std::cout << "  --- h_rates_gfexncmet slice contributions at MET > " << jzMETThr << " GeV ---\n";
            for (unsigned int jz = 0; jz < nJZSlices_; ++jz) {
                const double frac = (jzWeightedTotalNC > 0.0) ? (100.0 * jzWeightedPassNC[jz] / jzWeightedTotalNC) : 0.0;
                std::cout << "     JZ" << jz << ": "
                          << std::fixed << std::setw(6) << std::setprecision(2) << frac << "%   (n_events="
                          << (long long)back_h_gMET_NC_jz[jz]->GetEntries() << ")\n"
                          << std::defaultfloat;
            }
        }

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
        auto best_JwoJ_Jet_60   = findBestThresholds2D(out2D_JwoJ_Jet,   60e3);
        auto best_JwoJ_Tower_20 = findBestThresholds2D(out2D_JwoJ_Tower, 20e3);
        auto best_JwoJ_Tower_40 = findBestThresholds2D(out2D_JwoJ_Tower, 40e3);
        auto best_JwoJ_Tower_80 = findBestThresholds2D(out2D_JwoJ_Tower, 80e3);
        auto best_JwoJ_Tower_60 = findBestThresholds2D(out2D_JwoJ_Tower, 60e3);
        auto best_NC_Jet_20     = findBestThresholds2D(out2D_NC_Jet,     20e3);
        auto best_NC_Jet_40     = findBestThresholds2D(out2D_NC_Jet,     40e3);
        auto best_NC_Jet_80     = findBestThresholds2D(out2D_NC_Jet,     80e3);
        auto best_NC_Jet_60     = findBestThresholds2D(out2D_NC_Jet,     60e3);
        auto best_NC_Tower_20   = findBestThresholds2D(out2D_NC_Tower,   20e3);
        auto best_NC_Tower_40   = findBestThresholds2D(out2D_NC_Tower,   40e3);
        auto best_NC_Tower_80   = findBestThresholds2D(out2D_NC_Tower,   80e3);
        auto best_NC_Tower_60   = findBestThresholds2D(out2D_NC_Tower,   60e3);
        auto best_Rms_Jet_20    = findBestThresholds2D(out2D_Rms_Jet,    20e3);
        auto best_Rms_Jet_40    = findBestThresholds2D(out2D_Rms_Jet,    40e3);
        auto best_Rms_Jet_80    = findBestThresholds2D(out2D_Rms_Jet,    80e3);
        auto best_Rms_Jet_60    = findBestThresholds2D(out2D_Rms_Jet,    60e3);
        auto best_Rms_Tower_20  = findBestThresholds2D(out2D_Rms_Tower,  20e3);
        auto best_Rms_Tower_40  = findBestThresholds2D(out2D_Rms_Tower,  40e3);
        auto best_Rms_Tower_80  = findBestThresholds2D(out2D_Rms_Tower,  80e3);
        auto best_Rms_Tower_60  = findBestThresholds2D(out2D_Rms_Tower,  60e3);
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
            gNomJwoJSig->GetEntry(iEvt);   // nominal gFEX (resim when available, else AOD)
            gNomNCSig->GetEntry(iEvt);
            gNomRmsSig->GetEntry(iEvt);
            metTruthTreeSig->GetEntry(iEvt);
            double truthMET = clampVal(h_turnOn_denom, sig_metTruthNonInt);
            // Individual gFEX and GEP turn-ons
            if (sig_gMET     > thr_gMET_20kHz)     h_turnOn_num_gMET_20kHz->Fill(truthMET);
            if (sig_gMET     > thr_gMET_40kHz)     h_turnOn_num_gMET_40kHz->Fill(truthMET);
            if (sig_gMET     > thr_gMET_80kHz)     h_turnOn_num_gMET_80kHz->Fill(truthMET);
            if (sig_gMET     > thr_gMET_60kHz)     h_turnOn_num_gMET_60kHz->Fill(truthMET);
            if (sig_gMET_NC  > thr_gMET_NC_20kHz)  h_turnOn_num_gMET_NC_20kHz->Fill(truthMET);
            if (sig_gMET_NC  > thr_gMET_NC_40kHz)  h_turnOn_num_gMET_NC_40kHz->Fill(truthMET);
            if (sig_gMET_NC  > thr_gMET_NC_80kHz)  h_turnOn_num_gMET_NC_80kHz->Fill(truthMET);
            if (sig_gMET_NC  > thr_gMET_NC_60kHz)  h_turnOn_num_gMET_NC_60kHz->Fill(truthMET);
            if (sig_gMET_Rms > thr_gMET_Rms_20kHz) h_turnOn_num_gMET_Rms_20kHz->Fill(truthMET);
            if (sig_gMET_Rms > thr_gMET_Rms_40kHz) h_turnOn_num_gMET_Rms_40kHz->Fill(truthMET);
            if (sig_gMET_Rms > thr_gMET_Rms_80kHz) h_turnOn_num_gMET_Rms_80kHz->Fill(truthMET);
            if (sig_gMET_Rms > thr_gMET_Rms_60kHz) h_turnOn_num_gMET_Rms_60kHz->Fill(truthMET);
            jFexMETTreeSig->GetEntry(iEvt);
            if (sig_jMET     > thr_jMET_20kHz)     h_turnOn_num_jMET_20kHz->Fill(truthMET);
            if (sig_jMET     > thr_jMET_40kHz)     h_turnOn_num_jMET_40kHz->Fill(truthMET);
            if (sig_jMET     > thr_jMET_80kHz)     h_turnOn_num_jMET_80kHz->Fill(truthMET);
            if (sig_jMET     > thr_jMET_60kHz)     h_turnOn_num_jMET_60kHz->Fill(truthMET);
            if (sig_JetMet   > thr_JetMET_20kHz)   h_turnOn_num_JetMET_20kHz->Fill(truthMET);
            if (sig_JetMet   > thr_JetMET_40kHz)   h_turnOn_num_JetMET_40kHz->Fill(truthMET);
            if (sig_JetMet   > thr_JetMET_80kHz)   h_turnOn_num_JetMET_80kHz->Fill(truthMET);
            if (sig_JetMet   > thr_JetMET_60kHz)   h_turnOn_num_JetMET_60kHz->Fill(truthMET);
            if (sig_TowerMet  > thr_TowerMET_20kHz)  h_turnOn_num_TowerMET_20kHz->Fill(truthMET);
            if (sig_TowerMet  > thr_TowerMET_40kHz)  h_turnOn_num_TowerMET_40kHz->Fill(truthMET);
            if (sig_TowerMet  > thr_TowerMET_80kHz)  h_turnOn_num_TowerMET_80kHz->Fill(truthMET);
            if (sig_TowerMet  > thr_TowerMET_60kHz)  h_turnOn_num_TowerMET_60kHz->Fill(truthMET);
            if (sig_TotalMET  > thr_TotalMET_20kHz)  h_turnOn_num_TotalMET_20kHz->Fill(truthMET);
            if (sig_TotalMET  > thr_TotalMET_40kHz)  h_turnOn_num_TotalMET_40kHz->Fill(truthMET);
            if (sig_TotalMET  > thr_TotalMET_80kHz)  h_turnOn_num_TotalMET_80kHz->Fill(truthMET);
            if (sig_TotalMET  > thr_TotalMET_60kHz)  h_turnOn_num_TotalMET_60kHz->Fill(truthMET);

            // --- Z->mumu turn-on vs dimuon p_{T} ---
            // Denominator and numerators are filled here, in the same pass, so both see exactly
            // the same events. Events with fewer than two truth muons leave dimuonPt at its -1
            // default and are dropped from both — they carry no dimuon system to bin them by.
            if (isZmumuSample) {
                eventInfoTreeSig->GetEntry(iEvt);
                if (sig_dimuonPt >= 0.0) {
                    const double dimuonPtCl = clampVal(h_turnOnMu_denom, sig_dimuonPt);
                    h_turnOnMu_denom->Fill(dimuonPtCl);
                    sig_h_dimuonPt->Fill(clampVal(sig_h_dimuonPt, sig_dimuonPt));
                    sig_h_dimuonPt_coarse->Fill(std::min(sig_dimuonPt, 599.9));
                    sig_h_dimuonMass->Fill(sig_dimuonMass);
                    const double sigMETByType[nMETTypes] = {
                        sig_gMET, sig_gMET_NC, sig_gMET_Rms, sig_jMET,
                        sig_JetMet, sig_TowerMet, sig_TotalMET
                    };
                    for (int iA = 0; iA < nMETTypes; ++iA)
                        for (int iR = 0; iR < nMuRates; ++iR)
                            if (sigMETByType[iA] > thrMu[iA][iR])
                                h_turnOnMu_num[iA][iR]->Fill(dimuonPtCl);
                }
            }
            if (hasGFexSimMET) {
                gFexMETTreeSig->GetEntry(iEvt);          // AOD gFEX copies for the AOD-vs-Sim block
                gFexMETNoiseCutTreeSig->GetEntry(iEvt);
                gFexMETRmsTreeSig->GetEntry(iEvt);
                if (sig_gMET_JwoJAOD > thr_gMET_JwoJAOD_20kHz) h_turnOn_num_gMET_JwoJAOD_20kHz->Fill(truthMET);
                if (sig_gMET_JwoJAOD > thr_gMET_JwoJAOD_40kHz) h_turnOn_num_gMET_JwoJAOD_40kHz->Fill(truthMET);
                if (sig_gMET_JwoJAOD > thr_gMET_JwoJAOD_80kHz) h_turnOn_num_gMET_JwoJAOD_80kHz->Fill(truthMET);
                if (sig_gMET_JwoJAOD > thr_gMET_JwoJAOD_60kHz) h_turnOn_num_gMET_JwoJAOD_60kHz->Fill(truthMET);
                if (sig_gMET_NCAOD   > thr_gMET_NCAOD_20kHz)   h_turnOn_num_gMET_NCAOD_20kHz->Fill(truthMET);
                if (sig_gMET_NCAOD   > thr_gMET_NCAOD_40kHz)   h_turnOn_num_gMET_NCAOD_40kHz->Fill(truthMET);
                if (sig_gMET_NCAOD   > thr_gMET_NCAOD_80kHz)   h_turnOn_num_gMET_NCAOD_80kHz->Fill(truthMET);
                if (sig_gMET_NCAOD   > thr_gMET_NCAOD_60kHz)   h_turnOn_num_gMET_NCAOD_60kHz->Fill(truthMET);
                if (sig_gMET_RmsAOD  > thr_gMET_RmsAOD_20kHz)  h_turnOn_num_gMET_RmsAOD_20kHz->Fill(truthMET);
                if (sig_gMET_RmsAOD  > thr_gMET_RmsAOD_40kHz)  h_turnOn_num_gMET_RmsAOD_40kHz->Fill(truthMET);
                if (sig_gMET_RmsAOD  > thr_gMET_RmsAOD_80kHz)  h_turnOn_num_gMET_RmsAOD_80kHz->Fill(truthMET);
                if (sig_gMET_RmsAOD  > thr_gMET_RmsAOD_60kHz)  h_turnOn_num_gMET_RmsAOD_60kHz->Fill(truthMET);
            }
            // Combined gFEX + GEP turn-ons at best thresholds
            if (sig_gMET     > best_JwoJ_Jet_20.t1   && sig_JetMet   > best_JwoJ_Jet_20.t2)   h_turnOn_num_combo_JwoJ_Jet_20kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Jet_40.t1   && sig_JetMet   > best_JwoJ_Jet_40.t2)   h_turnOn_num_combo_JwoJ_Jet_40kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Jet_80.t1   && sig_JetMet   > best_JwoJ_Jet_80.t2)   h_turnOn_num_combo_JwoJ_Jet_80kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Jet_60.t1   && sig_JetMet   > best_JwoJ_Jet_60.t2)   h_turnOn_num_combo_JwoJ_Jet_60kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Tower_20.t1 && sig_TowerMet > best_JwoJ_Tower_20.t2) h_turnOn_num_combo_JwoJ_Tower_20kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Tower_40.t1 && sig_TowerMet > best_JwoJ_Tower_40.t2) h_turnOn_num_combo_JwoJ_Tower_40kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Tower_80.t1 && sig_TowerMet > best_JwoJ_Tower_80.t2) h_turnOn_num_combo_JwoJ_Tower_80kHz->Fill(truthMET);
            if (sig_gMET     > best_JwoJ_Tower_60.t1 && sig_TowerMet > best_JwoJ_Tower_60.t2) h_turnOn_num_combo_JwoJ_Tower_60kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Jet_20.t1     && sig_JetMet   > best_NC_Jet_20.t2)     h_turnOn_num_combo_NC_Jet_20kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Jet_40.t1     && sig_JetMet   > best_NC_Jet_40.t2)     h_turnOn_num_combo_NC_Jet_40kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Jet_80.t1     && sig_JetMet   > best_NC_Jet_80.t2)     h_turnOn_num_combo_NC_Jet_80kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Jet_60.t1     && sig_JetMet   > best_NC_Jet_60.t2)     h_turnOn_num_combo_NC_Jet_60kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Tower_20.t1   && sig_TowerMet > best_NC_Tower_20.t2)   h_turnOn_num_combo_NC_Tower_20kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Tower_40.t1   && sig_TowerMet > best_NC_Tower_40.t2)   h_turnOn_num_combo_NC_Tower_40kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Tower_80.t1   && sig_TowerMet > best_NC_Tower_80.t2)   h_turnOn_num_combo_NC_Tower_80kHz->Fill(truthMET);
            if (sig_gMET_NC  > best_NC_Tower_60.t1   && sig_TowerMet > best_NC_Tower_60.t2)   h_turnOn_num_combo_NC_Tower_60kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Jet_20.t1    && sig_JetMet   > best_Rms_Jet_20.t2)    h_turnOn_num_combo_Rms_Jet_20kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Jet_40.t1    && sig_JetMet   > best_Rms_Jet_40.t2)    h_turnOn_num_combo_Rms_Jet_40kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Jet_80.t1    && sig_JetMet   > best_Rms_Jet_80.t2)    h_turnOn_num_combo_Rms_Jet_80kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Jet_60.t1    && sig_JetMet   > best_Rms_Jet_60.t2)    h_turnOn_num_combo_Rms_Jet_60kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Tower_20.t1  && sig_TowerMet > best_Rms_Tower_20.t2)  h_turnOn_num_combo_Rms_Tower_20kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Tower_40.t1  && sig_TowerMet > best_Rms_Tower_40.t2)  h_turnOn_num_combo_Rms_Tower_40kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Tower_80.t1  && sig_TowerMet > best_Rms_Tower_80.t2)  h_turnOn_num_combo_Rms_Tower_80kHz->Fill(truthMET);
            if (sig_gMET_Rms > best_Rms_Tower_60.t1  && sig_TowerMet > best_Rms_Tower_60.t2)  h_turnOn_num_combo_Rms_Tower_60kHz->Fill(truthMET);

            // --- Event properties at 80 kHz selections ---
            if (hasWTAConeJets) gepWTAConeCellsTowersJetsTreeSig->GetEntry(iEvt);
            if (hasTruthAntiKt4WZDressed) truthAntiKt4WZDressedJetsTreeSig->GetEntry(iEvt);
            if (hasInTimeAntiKt4TruthJets) inTimeAntiKt4TruthJetsTreeSig->GetEntry(iEvt);
            {
                double phi_GEP  = std::atan2(sig_TotalMETY, sig_TotalMETX);
                //double phi_gFEX = std::atan2(sig_gMETY, sig_gMETX);
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
                /*double phi_TOBsig[nSel80] = {
                    phi_gFEX,
                    phi_gFEX,
                   // std::atan2(sig_gMETY_NC,  sig_gMETX_NC),
                    //std::atan2(sig_gMETY_Rms, sig_gMETX_Rms),
                    hasJetMetXY   ? std::atan2(sig_JetMetY,   sig_JetMetX)   : 0.0,
                    hasTowerMetXY ? std::atan2(sig_TowerMetY, sig_TowerMetX) : 0.0
                };*/
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
                    //sig_sel_dPhi_GEP_gFEX[iSel]->Fill(absDeltaPhi(phi_GEP, phi_gFEX));
                    if (hasTruthNonIntXY && (sig_metTruthNonIntX != 0.0 || sig_metTruthNonIntY != 0.0)) {
                        double phi_truth = std::atan2(sig_metTruthNonIntY, sig_metTruthNonIntX);
                        sig_sel_dPhi_GEP_truth[iSel]->Fill(absDeltaPhi(phi_GEP, phi_truth));
                        bool hasTOBphi = (iSel < 4) || (iSel == 4 && hasJetMetXY) || (iSel == 5 && hasTowerMetXY);
                        /*if (hasTOBphi) {
                            sig_sel_dPhi_truth_TOB[iSel]->Fill(absDeltaPhi(phi_truth, phi_TOBsig[iSel]));
                            sig_sel_phi2D_TOB_truth[iSel]->Fill(phi_truth, phi_TOBsig[iSel]);
                        }*/
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
            gNomJwoJBack->GetEntry(iEvt);   // nominal gFEX (resim when available, else AOD)
            gNomNCBack->GetEntry(iEvt);
            gNomRmsBack->GetEntry(iEvt);
            metTruthTreeBack->GetEntry(iEvt);
            eventInfoTreeBack->GetEntry(iEvt);
            if (applyHSTPFilter && !passHSTPValuesBack) continue;
            if (hasWTAConeJets) gepWTAConeCellsTowersJetsTreeBack->GetEntry(iEvt);
            if (hasTruthAntiKt4WZDressed) truthAntiKt4WZDressedJetsTreeBack->GetEntry(iEvt);
            if (hasInTimeAntiKt4TruthJets) inTimeAntiKt4TruthJetsTreeBack->GetEntry(iEvt);
            double w = eventWeightsValuesBack->at(0);

            double phi_GEP  = std::atan2(back_TotalMETY, back_TotalMETX);
            //double phi_gFEX = std::atan2(back_gMETY, back_gMETX);
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

            // --- <MET> vs jet multiplicity ---
            // Jets from both truth collections count towards the same multiplicity: the hard
            // scatter's AntiKt4 WZ-dressed jets and the in-time pileup overlay's AntiKt4 truth
            // jets, each above kNJetMinEt. Counted here rather than reusing nTruthJets /
            // nPileupJets above so the cut follows kNJetMinEt rather than those plots' own 15 GeV.
            if (hasTruthAntiKt4WZDressed || hasInTimeAntiKt4TruthJets) {
                int nJetsForProf = 0;
                if (hasTruthAntiKt4WZDressed && truthAntiKt4WZDressedJetsEtValuesBack)
                    for (double et : *truthAntiKt4WZDressedJetsEtValuesBack)
                        if (et > kNJetMinEt) nJetsForProf++;
                if (hasInTimeAntiKt4TruthJets && inTimeAntiKt4TruthJetsEtValuesBack)
                    for (double et : *inTimeAntiKt4TruthJetsEtValuesBack)
                        if (et > kNJetMinEt) nJetsForProf++;
                // jFEX MET is not loaded by this pass's GetEntry block above, unlike the other
                // six MET types, so pull its entry before reading back_jMET.
                jFexMETTreeBack->GetEntry(iEvt);
                const double backMETByType[nMETTypes] = {
                    back_gMET, back_gMET_NC, back_gMET_Rms, back_jMET,
                    back_JetMet, back_TowerMet, back_TotalMET
                };
                const double nJetsCl = std::min((double)nJetsForProf, nJetAxisMax - 1e-9);
                for (int iA = 0; iA < nMETTypes; ++iA)
                    back_prof_METvsNJets[iA]->Fill(nJetsCl, backMETByType[iA], w);
            }
            // Per-selection TOB MET phi: Incl/JwoJ=gFEX JwoJ, NC, Rms, Jet, Tower
            /*double phi_TOBback[nSel80] = {
                phi_gFEX,
                phi_gFEX,
                std::atan2(back_gMETY_NC,  back_gMETX_NC),
                std::atan2(back_gMETY_Rms, back_gMETX_Rms),
                hasJetMetXY   ? std::atan2(back_JetMetY,   back_JetMetX)   : 0.0,
                hasTowerMetXY ? std::atan2(back_TowerMetY, back_TowerMetX) : 0.0
            };*/
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
                //back_sel_dPhi_GEP_gFEX[iSel]->Fill(absDeltaPhi(phi_GEP, phi_gFEX), w);
                if (hasTruthNonIntXY && (back_metTruthNonIntX != 0.0 || back_metTruthNonIntY != 0.0)) {
                    double phi_truth = std::atan2(back_metTruthNonIntY, back_metTruthNonIntX);
                    back_sel_dPhi_GEP_truth[iSel]->Fill(absDeltaPhi(phi_GEP, phi_truth), w);
                    bool hasTOBphi = (iSel < 4) || (iSel == 4 && hasJetMetXY) || (iSel == 5 && hasTowerMetXY);
                    /*if (hasTOBphi) {
                        back_sel_dPhi_truth_TOB[iSel]->Fill(absDeltaPhi(phi_truth, phi_TOBback[iSel]), w);
                        back_sel_phi2D_TOB_truth[iSel]->Fill(phi_truth, phi_TOBback[iSel], w);
                    }*/
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
        TH1F* eff_gMET_60kHz      = makeEff(h_turnOn_num_gMET_60kHz,      "eff_gMET_60kHz_"     +tag);
        TH1F* eff_gMET_NC_20kHz   = makeEff(h_turnOn_num_gMET_NC_20kHz,   "eff_gMET_NC_20kHz_"  +tag);
        TH1F* eff_gMET_NC_40kHz   = makeEff(h_turnOn_num_gMET_NC_40kHz,   "eff_gMET_NC_40kHz_"  +tag);
        TH1F* eff_gMET_NC_80kHz   = makeEff(h_turnOn_num_gMET_NC_80kHz,   "eff_gMET_NC_80kHz_"  +tag);
        TH1F* eff_gMET_NC_60kHz   = makeEff(h_turnOn_num_gMET_NC_60kHz,   "eff_gMET_NC_60kHz_"  +tag);
        TH1F* eff_gMET_Rms_20kHz  = makeEff(h_turnOn_num_gMET_Rms_20kHz,  "eff_gMET_Rms_20kHz_" +tag);
        TH1F* eff_gMET_Rms_40kHz  = makeEff(h_turnOn_num_gMET_Rms_40kHz,  "eff_gMET_Rms_40kHz_" +tag);
        TH1F* eff_gMET_Rms_80kHz  = makeEff(h_turnOn_num_gMET_Rms_80kHz,  "eff_gMET_Rms_80kHz_" +tag);
        TH1F* eff_gMET_Rms_60kHz  = makeEff(h_turnOn_num_gMET_Rms_60kHz,  "eff_gMET_Rms_60kHz_" +tag);
        TH1F* eff_jMET_20kHz      = makeEff(h_turnOn_num_jMET_20kHz,      "eff_jMET_20kHz_"     +tag);
        TH1F* eff_jMET_40kHz      = makeEff(h_turnOn_num_jMET_40kHz,      "eff_jMET_40kHz_"     +tag);
        TH1F* eff_jMET_80kHz      = makeEff(h_turnOn_num_jMET_80kHz,      "eff_jMET_80kHz_"     +tag);
        TH1F* eff_jMET_60kHz      = makeEff(h_turnOn_num_jMET_60kHz,      "eff_jMET_60kHz_"     +tag);
        TH1F* eff_JetMET_20kHz    = makeEff(h_turnOn_num_JetMET_20kHz,    "eff_JetMET_20kHz_"   +tag);
        TH1F* eff_JetMET_40kHz    = makeEff(h_turnOn_num_JetMET_40kHz,    "eff_JetMET_40kHz_"   +tag);
        TH1F* eff_JetMET_80kHz    = makeEff(h_turnOn_num_JetMET_80kHz,    "eff_JetMET_80kHz_"   +tag);
        TH1F* eff_JetMET_60kHz    = makeEff(h_turnOn_num_JetMET_60kHz,    "eff_JetMET_60kHz_"   +tag);
        TH1F* eff_TowerMET_20kHz  = makeEff(h_turnOn_num_TowerMET_20kHz,  "eff_TowerMET_20kHz_" +tag);
        TH1F* eff_TowerMET_40kHz  = makeEff(h_turnOn_num_TowerMET_40kHz,  "eff_TowerMET_40kHz_" +tag);
        TH1F* eff_TowerMET_80kHz  = makeEff(h_turnOn_num_TowerMET_80kHz,  "eff_TowerMET_80kHz_" +tag);
        TH1F* eff_TowerMET_60kHz  = makeEff(h_turnOn_num_TowerMET_60kHz,  "eff_TowerMET_60kHz_" +tag);
        TH1F* eff_TotalMET_20kHz  = makeEff(h_turnOn_num_TotalMET_20kHz,  "eff_TotalMET_20kHz_" +tag);
        TH1F* eff_TotalMET_40kHz  = makeEff(h_turnOn_num_TotalMET_40kHz,  "eff_TotalMET_40kHz_" +tag);
        TH1F* eff_TotalMET_80kHz  = makeEff(h_turnOn_num_TotalMET_80kHz,  "eff_TotalMET_80kHz_" +tag);
        TH1F* eff_TotalMET_60kHz  = makeEff(h_turnOn_num_TotalMET_60kHz,  "eff_TotalMET_60kHz_" +tag);
        // Combined efficiencies
        TH1F* eff_combo_JwoJ_Jet_20kHz   = makeEff(h_turnOn_num_combo_JwoJ_Jet_20kHz,   "eff_combo_JwoJ_Jet_20kHz_"  +tag);
        TH1F* eff_combo_JwoJ_Jet_40kHz   = makeEff(h_turnOn_num_combo_JwoJ_Jet_40kHz,   "eff_combo_JwoJ_Jet_40kHz_"  +tag);
        TH1F* eff_combo_JwoJ_Jet_80kHz   = makeEff(h_turnOn_num_combo_JwoJ_Jet_80kHz,   "eff_combo_JwoJ_Jet_80kHz_"  +tag);
        TH1F* eff_combo_JwoJ_Jet_60kHz   = makeEff(h_turnOn_num_combo_JwoJ_Jet_60kHz,   "eff_combo_JwoJ_Jet_60kHz_"  +tag);
        TH1F* eff_combo_JwoJ_Tower_20kHz = makeEff(h_turnOn_num_combo_JwoJ_Tower_20kHz, "eff_combo_JwoJ_Tower_20kHz_"+tag);
        TH1F* eff_combo_JwoJ_Tower_40kHz = makeEff(h_turnOn_num_combo_JwoJ_Tower_40kHz, "eff_combo_JwoJ_Tower_40kHz_"+tag);
        TH1F* eff_combo_JwoJ_Tower_80kHz = makeEff(h_turnOn_num_combo_JwoJ_Tower_80kHz, "eff_combo_JwoJ_Tower_80kHz_"+tag);
        TH1F* eff_combo_JwoJ_Tower_60kHz = makeEff(h_turnOn_num_combo_JwoJ_Tower_60kHz, "eff_combo_JwoJ_Tower_60kHz_"+tag);
        TH1F* eff_combo_NC_Jet_20kHz     = makeEff(h_turnOn_num_combo_NC_Jet_20kHz,     "eff_combo_NC_Jet_20kHz_"    +tag);
        TH1F* eff_combo_NC_Jet_40kHz     = makeEff(h_turnOn_num_combo_NC_Jet_40kHz,     "eff_combo_NC_Jet_40kHz_"    +tag);
        TH1F* eff_combo_NC_Jet_80kHz     = makeEff(h_turnOn_num_combo_NC_Jet_80kHz,     "eff_combo_NC_Jet_80kHz_"    +tag);
        TH1F* eff_combo_NC_Jet_60kHz     = makeEff(h_turnOn_num_combo_NC_Jet_60kHz,     "eff_combo_NC_Jet_60kHz_"    +tag);
        TH1F* eff_combo_NC_Tower_20kHz   = makeEff(h_turnOn_num_combo_NC_Tower_20kHz,   "eff_combo_NC_Tower_20kHz_"  +tag);
        TH1F* eff_combo_NC_Tower_40kHz   = makeEff(h_turnOn_num_combo_NC_Tower_40kHz,   "eff_combo_NC_Tower_40kHz_"  +tag);
        TH1F* eff_combo_NC_Tower_80kHz   = makeEff(h_turnOn_num_combo_NC_Tower_80kHz,   "eff_combo_NC_Tower_80kHz_"  +tag);
        TH1F* eff_combo_NC_Tower_60kHz   = makeEff(h_turnOn_num_combo_NC_Tower_60kHz,   "eff_combo_NC_Tower_60kHz_"  +tag);
        TH1F* eff_combo_Rms_Jet_20kHz    = makeEff(h_turnOn_num_combo_Rms_Jet_20kHz,    "eff_combo_Rms_Jet_20kHz_"   +tag);
        TH1F* eff_combo_Rms_Jet_40kHz    = makeEff(h_turnOn_num_combo_Rms_Jet_40kHz,    "eff_combo_Rms_Jet_40kHz_"   +tag);
        TH1F* eff_combo_Rms_Jet_80kHz    = makeEff(h_turnOn_num_combo_Rms_Jet_80kHz,    "eff_combo_Rms_Jet_80kHz_"   +tag);
        TH1F* eff_combo_Rms_Jet_60kHz    = makeEff(h_turnOn_num_combo_Rms_Jet_60kHz,    "eff_combo_Rms_Jet_60kHz_"   +tag);
        TH1F* eff_combo_Rms_Tower_20kHz  = makeEff(h_turnOn_num_combo_Rms_Tower_20kHz,  "eff_combo_Rms_Tower_20kHz_" +tag);
        TH1F* eff_combo_Rms_Tower_40kHz  = makeEff(h_turnOn_num_combo_Rms_Tower_40kHz,  "eff_combo_Rms_Tower_40kHz_" +tag);
        TH1F* eff_combo_Rms_Tower_80kHz  = makeEff(h_turnOn_num_combo_Rms_Tower_80kHz,  "eff_combo_Rms_Tower_80kHz_" +tag);
        TH1F* eff_combo_Rms_Tower_60kHz  = makeEff(h_turnOn_num_combo_Rms_Tower_60kHz,  "eff_combo_Rms_Tower_60kHz_" +tag);
        TH1F* eff_gMET_JwoJAOD_20kHz = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_JwoJAOD_20kHz, "eff_gMET_JwoJAOD_20kHz_"+tag) : nullptr;
        TH1F* eff_gMET_JwoJAOD_40kHz = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_JwoJAOD_40kHz, "eff_gMET_JwoJAOD_40kHz_"+tag) : nullptr;
        TH1F* eff_gMET_JwoJAOD_80kHz = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_JwoJAOD_80kHz, "eff_gMET_JwoJAOD_80kHz_"+tag) : nullptr;
        TH1F* eff_gMET_JwoJAOD_60kHz = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_JwoJAOD_60kHz, "eff_gMET_JwoJAOD_60kHz_"+tag) : nullptr;
        TH1F* eff_gMET_NCAOD_20kHz   = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_NCAOD_20kHz,   "eff_gMET_NCAOD_20kHz_"+tag)   : nullptr;
        TH1F* eff_gMET_NCAOD_40kHz   = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_NCAOD_40kHz,   "eff_gMET_NCAOD_40kHz_"+tag)   : nullptr;
        TH1F* eff_gMET_NCAOD_80kHz   = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_NCAOD_80kHz,   "eff_gMET_NCAOD_80kHz_"+tag)   : nullptr;
        TH1F* eff_gMET_NCAOD_60kHz   = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_NCAOD_60kHz,   "eff_gMET_NCAOD_60kHz_"+tag)   : nullptr;
        TH1F* eff_gMET_RmsAOD_20kHz  = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_RmsAOD_20kHz,  "eff_gMET_RmsAOD_20kHz_"+tag)  : nullptr;
        TH1F* eff_gMET_RmsAOD_40kHz  = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_RmsAOD_40kHz,  "eff_gMET_RmsAOD_40kHz_"+tag)  : nullptr;
        TH1F* eff_gMET_RmsAOD_80kHz  = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_RmsAOD_80kHz,  "eff_gMET_RmsAOD_80kHz_"+tag)  : nullptr;
        TH1F* eff_gMET_RmsAOD_60kHz  = hasGFexSimMET ? makeEff(h_turnOn_num_gMET_RmsAOD_60kHz,  "eff_gMET_RmsAOD_60kHz_"+tag)  : nullptr;

        // Z->mumu dimuon-p_{T} efficiencies — same binomial division, against the dimuon
        // denominator rather than the truth-MET one.
        TH1F* effMu[nMETTypes][nMuRates] = {};
        if (isZmumuSample) {
            for (int iA = 0; iA < nMETTypes; ++iA)
                for (int iR = 0; iR < nMuRates; ++iR) {
                    TH1F* e = (TH1F*)h_turnOnMu_num[iA][iR]->Clone(
                        (std::string("effMu_") + metTypeShort[iA] + "_" + muRateNames[iR] + "_" + tag).c_str());
                    e->SetDirectory(0);
                    e->Divide(h_turnOnMu_num[iA][iR], h_turnOnMu_denom, 1.0, 1.0, "B");
                    effMu[iA][iR] = e;
                }
            std::cout << "  [Z->mumu] events with a dimuon system = "
                      << (long long)h_turnOnMu_denom->GetEntries() << " of " << nSig << "\n";
        }

        // --- Per-file signal vs background overlays ---
        // Extract signal tag from the emulation-output basename (includes config tags: N_Towers, jetEt, SK/OR)
        std::string sigPath = signalFiles[fileIt].second;
        std::string sigBase = sigPath.substr(sigPath.rfind('/') + 1);
        std::string sigTag  = sigBase.substr(0, sigBase.rfind('.'));  // strip .root
        std::string fDir = outputDir + "metPlots/" + sigTag + "_" + labels[fileIt] + "/";
        gSystem->mkdir(fDir.c_str(), true);
        // Per-file legend header — use the process-specific name when provided, else the global one.
        std::string fileSignalName = (fileIt < signalNames.size() && !signalNames[fileIt].empty())
            ? signalNames[fileIt] : signalName;
        // Per-file plots show the process name at the top-right of the ATLAS label (not in legends).
        gProcLabel = fileSignalName;
        // Quote the pileup of the sample actually being processed, not a hardcoded 200.
        SetPileupFromPath(signalFiles[fileIt].first + " " + signalFiles[fileIt].second);
        std::cout << "  Labelling plots as <PU> = " << gPileup << "\n";

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
        //drawComponentOverlay(sig_h_gMETX_Rms, back_h_gMETX_Rms, "gFEX MET_{x} (Rms)",      "gFEX MET_{x} (rms) [GeV]", fDir + "gFEX_METX_Rms.pdf");
        //drawComponentOverlay(sig_h_gMETY_Rms, back_h_gMETY_Rms, "gFEX MET_{y} (Rms)",      "gFEX MET_{y} (rms) [GeV]", fDir + "gFEX_METY_Rms.pdf");
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
                           "MET [GeV]", fDir + "GEP_vs_gFEX_MET.pdf", "");

        // --- Per-file core term sig vs bkg overlays ---
        // drawOverlay(sig_h_coreEMTopo_SoftClus_MET,   back_h_coreEMTopo_SoftClus_MET,   "Core EMTopo SoftClus MET",   "Core EMTopo SoftClus MET [GeV]", fDir + "CoreEMTopo_SoftClus_MET.pdf");
        // drawOverlay(sig_h_coreEMTopo_PVSoftTrk_MET,  back_h_coreEMTopo_PVSoftTrk_MET,  "Core EMTopo PVSoftTrk MET",  "Core EMTopo PVSoftTrk MET [GeV]", fDir + "CoreEMTopo_PVSoftTrk_MET.pdf");
        // drawOverlay(sig_h_coreEMTopo_SoftClusEM_MET, back_h_coreEMTopo_SoftClusEM_MET, "Core EMTopo SoftClusEM MET", "Core EMTopo SoftClusEM MET [GeV]", fDir + "CoreEMTopo_SoftClusEM_MET.pdf");
        // drawOverlay(sig_h_coreEMPFlow_SoftClus_MET,  back_h_coreEMPFlow_SoftClus_MET,  "Core EMPFlow SoftClus MET",  "Core EMPFlow SoftClus MET [GeV]", fDir + "CoreEMPFlow_SoftClus_MET.pdf");
        // drawOverlay(sig_h_coreEMPFlow_PVSoftTrk_MET, back_h_coreEMPFlow_PVSoftTrk_MET, "Core EMPFlow PVSoftTrk MET", "Core EMPFlow PVSoftTrk MET [GeV]", fDir + "CoreEMPFlow_PVSoftTrk_MET.pdf");

        // --- Multi-MET comparison: all types on one canvas, signal ---
        drawMultiDist({sig_h_TotalMET, sig_h_gMET, sig_h_metTruthNonInt},
                    // sig_h_coreEMTopo_SoftClus_MET, sig_h_coreEMTopo_PVSoftTrk_MET, sig_h_coreEMTopo_SoftClusEM_MET,
                    // sig_h_coreEMPFlow_SoftClus_MET, sig_h_coreEMPFlow_PVSoftTrk_MET},
                      {"GEP TotalMET", "gFEX MET", "Truth NonInt"},
                    // "EMTopo SoftClus", "EMTopo PVSoftTrk", "EMTopo SoftClusEM",
                    // "EMPFlow SoftClus", "EMPFlow PVSoftTrk"},
                      "MET comparison — signal", "MET [GeV]", fDir + "METComparison_sig.pdf");

        // --- Multi-MET comparison: all types on one canvas, background ---
        drawMultiDist({back_h_TotalMET, back_h_gMET, back_h_metTruthNonInt},
                    // back_h_coreEMTopo_SoftClus_MET, back_h_coreEMTopo_PVSoftTrk_MET, back_h_coreEMTopo_SoftClusEM_MET,
                    // back_h_coreEMPFlow_SoftClus_MET, back_h_coreEMPFlow_PVSoftTrk_MET},
                      {"GEP TotalMET", "gFEX MET", "Truth NonInt"},
                    // "EMTopo SoftClus", "EMTopo PVSoftTrk", "EMTopo SoftClusEM",
                    // "EMPFlow SoftClus", "EMPFlow PVSoftTrk"},
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
            //drawMultiDist(toVec(sig_sel_dPhi_GEP_gFEX), selLegLabels, "#Delta#phi(GEP MET, gFEX MET) — signal at 80 kHz sel.",         "|#Delta#phi(GEP MET, gFEX MET)| [rad]",                              selDir + "sig_sel_dPhi_GEP_gFEX.pdf",  false, "rad");
            if (hasTruthNonIntXY) {
                drawMultiDist(toVec(sig_sel_dPhi_GEP_truth), selLegLabels, "#Delta#phi(GEP MET, truth MET) — signal at 80 kHz sel.",    "|#Delta#phi(GEP MET, truth MET)| [rad]",                             selDir + "sig_sel_dPhi_GEP_truth.pdf",  false, "rad");
                drawMultiDist(toVec(sig_sel_dPhi_truth_TOB), selLegLabels, "#Delta#phi(truth MET, TOB MET) — signal at 80 kHz sel.",    "|#Delta#phi(truth MET, TOB MET)| [rad]",                             selDir + "sig_sel_dPhi_truth_TOB.pdf",  false, "rad");
                for (int iSel = 0; iSel < nSel80; iSel++) {
                    TCanvas c2D(("c_sig_phi2D_"+std::string(selShortNames[iSel])).c_str(), "", 700, 600);
                    c2D.SetRightMargin(0.15);
                    sig_sel_phi2D_TOB_truth[iSel]->SetTitle(Form("Signal — %s;Truth MET #phi [rad];TOB MET #phi [rad]", selLegLabels[iSel].c_str()));
                    sig_sel_phi2D_TOB_truth[iSel]->Draw("COLZ");
                    gPad->SetLogz(sig_sel_phi2D_TOB_truth[iSel]->GetMaximum() > 0);
                    c2D.cd(); DrawATLASLabel(); c2D.SaveAs((selDir + "sig_sel_phi2D_TOB_truth_" + selShortNames[iSel] + ".pdf").c_str());
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
            BkgProcLabel bkgProcSel;   // everything below this point in the block is background
            drawMultiDist(toVec(back_sel_truthMET),      selLegLabels, "Truth MET_{NonInt} — bkg at 80 kHz sel.",                       "Truth MET_{NonInt} [GeV]",                                           selDir + "back_sel_truthMET.pdf");
            drawMultiDist(toVec(back_sel_SumET),         selLegLabels, "GEP TOB #Sigma E_{T} — bkg at 80 kHz sel.",                     "GEP TOB #Sigma E_{T} [GeV]",                                         selDir + "back_sel_SumET.pdf");
            drawMultiDist(toVec(back_sel_gSumET),        selLegLabels, "gFEX TOB #Sigma E_{T} — bkg at 80 kHz sel.",                    "gFEX TOB #Sigma E_{T} [GeV]",                                        selDir + "back_sel_gSumET.pdf");
            drawMultiDist(toVec(back_sel_METsig),        selLegLabels, "GEP MET significance — bkg at 80 kHz sel.",                     "GEP MET / #sqrt{GEP TOB #Sigma E_{T}} [#sqrt{GeV}]",                 selDir + "back_sel_METsig.pdf",        true,  "#sqrt{GeV}");
            drawMultiDist(toVec(back_sel_gMETsig),       selLegLabels, "gFEX MET significance — bkg at 80 kHz sel.",                    "gFEX MET / #sqrt{gFEX TOB #Sigma E_{T}} [#sqrt{GeV}]",               selDir + "back_sel_gMETsig.pdf",       true,  "#sqrt{GeV}");
            //drawMultiDist(toVec(back_sel_dPhi_GEP_gFEX),selLegLabels, "#Delta#phi(GEP MET, gFEX MET) — bkg at 80 kHz sel.",            "|#Delta#phi(GEP MET, gFEX MET)| [rad]",                              selDir + "back_sel_dPhi_GEP_gFEX.pdf", false, "rad");
            if (hasTruthNonIntXY) {
                drawMultiDist(toVec(back_sel_dPhi_GEP_truth), selLegLabels, "#Delta#phi(GEP MET, truth MET) — bkg at 80 kHz sel.",      "|#Delta#phi(GEP MET, truth MET)| [rad]",                             selDir + "back_sel_dPhi_GEP_truth.pdf", false, "rad");
                drawMultiDist(toVec(back_sel_dPhi_truth_TOB), selLegLabels, "#Delta#phi(truth MET, TOB MET) — bkg at 80 kHz sel.",      "|#Delta#phi(truth MET, TOB MET)| [rad]",                             selDir + "back_sel_dPhi_truth_TOB.pdf", false, "rad");
                for (int iSel = 0; iSel < nSel80; iSel++) {
                    TCanvas c2D(("c_back_phi2D_"+std::string(selShortNames[iSel])).c_str(), "", 700, 600);
                    c2D.SetRightMargin(0.15);
                    back_sel_phi2D_TOB_truth[iSel]->SetTitle(Form("Background — %s;Truth MET #phi [rad];TOB MET #phi [rad]", selLegLabels[iSel].c_str()));
                    back_sel_phi2D_TOB_truth[iSel]->Draw("COLZ");
                    gPad->SetLogz(back_sel_phi2D_TOB_truth[iSel]->GetMaximum() > 0);
                    c2D.cd(); DrawATLASLabel(); c2D.SaveAs((selDir + "back_sel_phi2D_TOB_truth_" + selShortNames[iSel] + ".pdf").c_str());
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

        // --- Reconstructed primary vertices ---
        if (hasNPrimaryVertices && back_h_nPrimaryVertices->GetEntries() > 0) {
            BkgProcLabel bkgProc;   // background only
            TCanvas cPV(("c_nPrimaryVertices_"+tag).c_str(), "Primary vertices", 700, 600);
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
            back_h_nPrimaryVertices->SetLineColor(kBlack);
            back_h_nPrimaryVertices->SetLineWidth(2);
            back_h_nPrimaryVertices->Draw("HIST");
            cPV.cd(); DrawATLASLabel();
            cPV.SaveAs((fDir + "back_nPrimaryVertices.pdf").c_str());

            cPV.SetLogy();
            back_hw_nPrimaryVertices->SetLineColor(kP10Blue);
            back_hw_nPrimaryVertices->SetLineWidth(2);
            back_hw_nPrimaryVertices->Draw("HIST");
            cPV.cd(); DrawATLASLabel();
            cPV.SaveAs((fDir + "back_nPrimaryVertices_weighted.pdf").c_str());
            std::cout << "  [nPV] mean primary vertices = " << back_h_nPrimaryVertices->GetMean()
                      << " (raw), " << back_hw_nPrimaryVertices->GetMean() << " (rate weighted)\n";
        }

        // --- Per-JZ-slice MET distributions ---
        TString jzDir = (fDir + "JZSlices/").c_str();
        gSystem->mkdir(jzDir);
        {
        BkgProcLabel bkgProc;   // every JZ-slice distribution is background only
        OverlayAndSave(back_h_gMET_jz,     nJZSlices_, "c_gMET_jz",     jzDir + "gFEX_MET_JZSlices.pdf",        0);
        OverlayAndSave(back_h_gMET_NC_jz,  nJZSlices_, "c_gMET_NC_jz",  jzDir + "gFEX_NoiseCut_MET_JZSlices.pdf", 0);
        OverlayAndSave(back_h_gMET_Rms_jz, nJZSlices_, "c_gMET_Rms_jz", jzDir + "gFEX_RMS_MET_JZSlices.pdf",     0);
        OverlayAndSave(back_h_jMET_jz,     nJZSlices_, "c_jMET_jz",     jzDir + "jFEX_MET_JZSlices.pdf",        0);
        OverlayAndSave(back_h_JetMET_jz,   nJZSlices_, "c_JetMET_jz",   jzDir + "GEP_JetMET_JZSlices.pdf",     0);
        OverlayAndSave(back_h_TowerMET_jz, nJZSlices_, "c_TowerMET_jz", jzDir + "GEP_TowerMET_JZSlices.pdf",   0);
        }

        // --- Per-JZ-slice rate vs threshold (points, kHz, up to 200 GeV) ---
        OverlayRateAndSave(back_h_gMET_jz,     nJZSlices_, "c_gMET_jz_rate",     jzDir + "gFEX_MET_JZSlices_Rate.pdf",       "gFEX JwoJ");
        OverlayRateAndSave(back_h_gMET_NC_jz,  nJZSlices_, "c_gMET_NC_jz_rate",  jzDir + "gFEX_NoiseCut_MET_JZSlices_Rate.pdf", "gFEX NoiseCut");
        OverlayRateAndSave(back_h_gMET_Rms_jz, nJZSlices_, "c_gMET_Rms_jz_rate", jzDir + "gFEX_RMS_MET_JZSlices_Rate.pdf",    "gFEX Rms");
        OverlayRateAndSave(back_h_jMET_jz,     nJZSlices_, "c_jMET_jz_rate",     jzDir + "jFEX_MET_JZSlices_Rate.pdf",       "jFEX");
        OverlayRateAndSave(back_h_JetMET_jz,   nJZSlices_, "c_JetMET_jz_rate",   jzDir + "GEP_JetMET_JZSlices_Rate.pdf",     "GEP Jet MET");
        OverlayRateAndSave(back_h_TowerMET_jz, nJZSlices_, "c_TowerMET_jz_rate", jzDir + "GEP_TowerMET_JZSlices_Rate.pdf",   "GEP Tower MET");

        // --- JZ0-only (no HSTP) vs all-JZ rate vs threshold, one plot per MET type ---
        // (background only — the guard inside the block below labels these QCD dijet)
        // hAllJZ is the HSTP-filtered sum over every JZ slice, and hJZ0 the JZ0 events with the
        // HSTP requirement dropped.
        {
            BkgProcLabel bkgProc;
            // One plot per MET type, four curves:
            //   * all JZ slices (HSTP filtered), normalized so its rate at threshold 0 is
            //     targetTotalRateHz,
            //   * JZ0 alone with the HSTP requirement dropped, normalized the same way,
            //   * all JZ slices (HSTP filtered) as filled, i.e. the collision rate the weights
            //     actually carry, for reference against the two normalized curves,
            //   * the same events with the truth-jet binomial correction applied per event, i.e.
            //     the per-crossing counterpart of the as-filled curve (see the block that builds
            //     h_binomialCorr). Drawn only when the truth jets were available.
            // The as-filled JZ0 curve is deliberately not drawn: it sits ~180x above the
            // normalized ones, which no single ratio panel can span.
            auto drawJZ0Comparison = [&](TH1F* hAllJZ, TH1F* hJZ0, TH1F* hAllJZBinomial,
                                         const std::string& algoName, const std::string& fileStem) {
                if (!hAllJZ || !hJZ0) return;
                const std::string hstpSuffix = applyHSTPFilter ? " (HSTP filtered)" : "";
                const std::vector<Int_t> jz0Cols = { kBlack, kP10Red, kP10Blue, kP10Green };

                TH1* cumAll = MakeCumulativeRateHist(hAllJZ, (std::string(hAllJZ->GetName()) + "_cumAll").c_str());
                TH1* cumJZ0 = MakeCumulativeRateHist(hJZ0,   (std::string(hJZ0->GetName())   + "_cumJZ0").c_str());
                if (!cumAll || !cumJZ0) return;

                // Normalized copies: each scaled so its own first bin is targetTotalRateHz.
                // hAllJZ was already normalized upstream, so normAll is a no-op scale on it and
                // the as-filled reference curve is recovered by undoing the factor that was
                // applied there. hJZ0 is never scaled upstream, so it is normalized here.
                TH1* normAll = (TH1*)cumAll->Clone((std::string(hAllJZ->GetName()) + "_cumAllNorm").c_str());
                TH1* normJZ0 = (TH1*)cumJZ0->Clone((std::string(hJZ0->GetName())   + "_cumJZ0Norm").c_str());
                normAll->SetDirectory(nullptr); normJZ0->SetDirectory(nullptr);
                const double allAtZero = cumAll->GetBinContent(1);
                const double jz0AtZero = cumJZ0->GetBinContent(1);
                if (allAtZero > 0.0) normAll->Scale(targetTotalRateHz / allAtZero);
                if (jz0AtZero > 0.0) normJZ0->Scale(targetTotalRateHz / jz0AtZero);

                const double appliedAll = appliedRateScale.count(hAllJZ) ? appliedRateScale[hAllJZ] : 1.0;
                TH1* rawAll = (TH1*)cumAll->Clone((std::string(hAllJZ->GetName()) + "_cumAllRaw").c_str());
                rawAll->SetDirectory(nullptr);
                if (appliedAll > 0.0 && appliedAll != 1.0) rawAll->Scale(1.0 / appliedAll);
                const double allAtZeroAsFilled = appliedAll > 0.0 ? allAtZero / appliedAll : allAtZero;

                std::vector<TH1*> normCurves = { normAll, normJZ0, rawAll };
                // Scale factors go to the log, not the legend — the factor on its own is
                // unreadable, and the number worth having is where each curve started.
                printf("  [norm] %-20s %-28s from %10.4g MHz at threshold 0, scaled by %.5g\n",
                       algoName.c_str(), ("All JZ" + hstpSuffix).c_str(), allAtZeroAsFilled / 1e6,
                       allAtZeroAsFilled > 0.0 ? targetTotalRateHz / allAtZeroAsFilled : 1.0);
                printf("  [norm] %-20s %-28s from %10.4g MHz at threshold 0, scaled by %.5g\n",
                       algoName.c_str(), "JZ0 only (no HSTP)", jz0AtZero / 1e6,
                       jz0AtZero > 0.0 ? targetTotalRateHz / jz0AtZero : 1.0);

                std::vector<std::string> normLabels;
                normLabels.push_back("All JZ slices" + hstpSuffix + ", normalized");
                normLabels.push_back("JZ0 only (no HSTP), normalized");
                normLabels.push_back("All JZ slices" + hstpSuffix + ", not scaled");

                // Per-crossing curve: the as-filled events re-weighted by c(E_T^lead truth jet).
                if (h_binomialCorr && hAllJZBinomial) {
                    TH1* cumBin = MakeCumulativeRateHist(hAllJZBinomial,
                                                         (std::string(hAllJZBinomial->GetName()) + "_cumBin").c_str());
                    if (cumBin) {
                        normCurves.push_back(cumBin);
                        normLabels.push_back(Form("All JZ slices%s, binomial correction"
                                                  " (#mu=%.0f, f_{BX}=%.1f MHz)",
                                                  hstpSuffix.c_str(), binomialPileup,
                                                  kCrossingRateHz / 1e6));
                        printf("  [binomial] %-20s %-28s %10.4g MHz at threshold 0"
                               " (as filled %10.4g MHz)\n",
                               algoName.c_str(), "All JZ, per crossing",
                               cumBin->GetBinContent(1) / 1e6, allAtZeroAsFilled / 1e6);
                    }
                }

                DrawRateCurvesWithRatio(normCurves, normLabels, jz0Cols,
                                        /*refIndex=*/1, /*refShortName=*/"JZ0",
                                        "MET threshold [GeV]", "Estimated Background Rate [kHz]",
                                        algoName,
                                        TString((fDir + "JZSlices/" + fileStem + "_JZ0_vs_AllJZ_Rate.pdf").c_str()),
                                        /*yScale=*/1e-3, /*xMax=*/200.0, /*yMin=*/1e-3,
                                        /*ratioMin=*/0.0, /*ratioMax=*/3.0);
            };
            drawJZ0Comparison(back_hw_gMET,     back_hw_gMET_JZ0,     back_hwBin_gMET,     "gFEX JwoJ",     "gFEX_MET");
            drawJZ0Comparison(back_hw_gMET_NC,  back_hw_gMET_NC_JZ0,  back_hwBin_gMET_NC,  "gFEX NoiseCut", "gFEX_NoiseCut_MET");
            drawJZ0Comparison(back_hw_gMET_Rms, back_hw_gMET_Rms_JZ0, back_hwBin_gMET_Rms, "gFEX Rms",      "gFEX_RMS_MET");
            drawJZ0Comparison(back_hw_jMET,     back_hw_jMET_JZ0,     back_hwBin_jMET,     "jFEX",          "jFEX_MET");
            drawJZ0Comparison(back_hw_JetMET,   back_hw_JetMET_JZ0,   back_hwBin_JetMET,   "GEP Jet MET",   "GEP_JetMET");
            if (!hasOverlapRemoval)
                drawJZ0Comparison(back_hw_TowerMET, back_hw_TowerMET_JZ0, back_hwBin_TowerMET, "GEP Tower MET", "GEP_TowerMET");
            drawJZ0Comparison(back_hw_TotalMET, back_hw_TotalMET_JZ0, back_hwBin_TotalMET, "GEP Total MET", "GEP_TotalMET");
            if (hasGFexSimMET) {
                drawJZ0Comparison(back_hw_gMET_JwoJAOD, back_hw_gMET_JwoJAOD_JZ0, back_hwBin_gMET_JwoJAOD, "gFEX JwoJ (AOD)",     "gFEX_MET_AOD");
                drawJZ0Comparison(back_hw_gMET_NCAOD,   back_hw_gMET_NCAOD_JZ0,   back_hwBin_gMET_NCAOD,   "gFEX NoiseCut (AOD)", "gFEX_NoiseCut_MET_AOD");
                drawJZ0Comparison(back_hw_gMET_RmsAOD,  back_hw_gMET_RmsAOD_JZ0,  back_hwBin_gMET_RmsAOD,  "gFEX Rms (AOD)",      "gFEX_RMS_MET_AOD");
            }
        }

        // --- Per-file rate vs threshold plots ---
        drawRateVsThreshold(back_hw_TotalMET,  "Rate vs Emulated MET threshold",          "MET threshold [GeV]", fDir + "Rate_TotalMET.pdf");
        drawRateVsThreshold(back_hw_gMET,      "Rate vs gFEX MET threshold (JwoJ)",       "MET threshold [GeV]", fDir + "Rate_gFEX_MET_JwoJ.pdf");
        drawRateVsThreshold(back_hw_gMET_NC,   "Rate vs gFEX MET threshold (NoiseCut)",   "MET threshold [GeV]", fDir + "Rate_gFEX_MET_NoiseCut.pdf");
        drawRateVsThreshold(back_hw_gMET_Rms,  "Rate vs gFEX MET threshold (Rms)",        "MET threshold [GeV]", fDir + "Rate_gFEX_MET_Rms.pdf");
        drawRateVsThreshold(back_hw_jMET,      "Rate vs jFEX MET threshold",              "MET threshold [GeV]", fDir + "Rate_jFEX_MET.pdf");
        drawRateVsThreshold(back_hw_JetMET,    "Rate vs GEP Jet MET threshold",           "MET threshold [GeV]", fDir + "Rate_JetMET.pdf");
        drawRateVsThreshold(back_hw_TowerMET,  "Rate vs GEP Tower MET threshold",         "MET threshold [GeV]", fDir + "Rate_TowerMET.pdf");
        // H_T runs to 1040 GeV, well past the MET threshold cut-off, so it keeps its full range.
        if (hasSumJetET)
            drawRateVsThreshold(back_hw_SumJetET, "Rate vs GEP H_{T} threshold",
                                "H_{T} threshold [GeV]", fDir + "Rate_SumJetET.pdf", /*xMax=*/-1.0);

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
                                    fDir + "SigEff_vs_Threshold_AlgoComparison.pdf", "");
            // gFEX-only version
            std::vector<TH1F*> gfexHists  = {sig_h_gMET, sig_h_gMET_NC, sig_h_gMET_Rms};
            std::vector<std::string> gfexL = {"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms"};
            drawEffVsThresholdMulti(gfexHists, gfexL,
                                    "Signal Efficiency vs gFEX MET Threshold", "MET threshold [GeV]",
                                    fDir + "SigEff_vs_Threshold_gFEX_AlgoComparison.pdf", "");
            // GEP-only version
            std::vector<TH1F*> gepHists = hasOverlapRemoval
                ? std::vector<TH1F*>{sig_h_JetMet, sig_h_TotalMET}
                : std::vector<TH1F*>{sig_h_JetMet, sig_h_TowerMet, sig_h_TotalMET};
            std::vector<std::string> gepL = hasOverlapRemoval
                ? std::vector<std::string>{"GEP Jet MET", "GEP Total MET"}
                : std::vector<std::string>{"GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
            drawEffVsThresholdMulti(gepHists, gepL,
                                    "Signal Efficiency vs GEP MET Threshold", "MET threshold [GeV]",
                                    fDir + "SigEff_vs_Threshold_GEP_AlgoComparison.pdf", "");
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
                fitFn->SetLineColor(kP10Blue); fitFn->SetLineWidth(2);
                fitFn->Draw("SAME"); prof->Draw("SAME");
                double rng = std::min(xmax, ymax);
                TLine* diag = new TLine(0, 0, rng, rng);
                diag->SetLineColor(kP10Red); diag->SetLineStyle(2); diag->SetLineWidth(2);
                diag->Draw("SAME");
                TLegend leg(0.16, 0.62, 0.61, 0.88);
                leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);
                leg.AddEntry(diag, "TOB MET = Truth MET", "l");
                leg.AddEntry(fitFn, Form("Fit: y = %.3f x + %.1f GeV", slope, intercept), "l");
                leg.AddEntry((TObject*)nullptr, Form("r = %.4f", r), "");
                leg.Draw();
                cTmp.cd(); DrawATLASLabel(); cTmp.SaveAs(path.c_str());
                delete prof; delete fitFn;
            };
            // Signal
            drawTOBvsTruth(sig_h2_gJwoJ_TOBMet_vs_truthMET,    calDir + "sig_gFEX_JwoJ_TOBMet_vs_truthMET.pdf",    "gFEX JwoJ");
            drawTOBvsTruth(sig_h2_gNC_TOBMet_vs_truthMET,      calDir + "sig_gFEX_NoiseCut_TOBMet_vs_truthMET.pdf","gFEX NoiseCut");
            drawTOBvsTruth(sig_h2_gRms_TOBMet_vs_truthMET,     calDir + "sig_gFEX_Rms_TOBMet_vs_truthMET.pdf",     "gFEX Rms");
            drawTOBvsTruth(sig_h2_JetMET_TOBMet_vs_truthMET,   calDir + "sig_GEP_JetMET_TOBMet_vs_truthMET.pdf",   "GEP Jet");
            drawTOBvsTruth(sig_h2_TowerMET_TOBMet_vs_truthMET, calDir + "sig_GEP_TowerMET_TOBMet_vs_truthMET.pdf", "GEP Tower");
            // Background (lower z-floor to show high-MET tails; x capped at 300 GeV)
            {
            BkgProcLabel bkgProcCal;   // background calibration plots
            drawTOBvsTruth(back_h2_gJwoJ_TOBMet_vs_truthMET,    calDir + "back_gFEX_JwoJ_TOBMet_vs_truthMET.pdf",    "gFEX JwoJ",    1e-14, 300.0);
            drawTOBvsTruth(back_h2_gNC_TOBMet_vs_truthMET,      calDir + "back_gFEX_NoiseCut_TOBMet_vs_truthMET.pdf","gFEX NoiseCut",1e-14, 300.0);
            drawTOBvsTruth(back_h2_gRms_TOBMet_vs_truthMET,     calDir + "back_gFEX_Rms_TOBMet_vs_truthMET.pdf",     "gFEX Rms",     1e-14, 300.0);
            drawTOBvsTruth(back_h2_JetMET_TOBMet_vs_truthMET,   calDir + "back_GEP_JetMET_TOBMet_vs_truthMET.pdf",   "GEP Jet",      1e-14, 300.0);
            drawTOBvsTruth(back_h2_TowerMET_TOBMet_vs_truthMET, calDir + "back_GEP_TowerMET_TOBMet_vs_truthMET.pdf", "GEP Tower",    1e-14, 300.0);
            }

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
                h->SetLineColor(kP10Blue); h->SetLineWidth(2);
                std::string yT = Form("Fraction of Events / %.4g", h->GetBinWidth(1));
                if (!yUnit.empty()) yT += " " + yUnit;
                h->SetTitle((tobLabel + ";" + xLbl + ";" + yT).c_str());
                h->SetMinimum(1e-5);
                h->Draw("HIST");
                TLine* zero = new TLine(0, h->GetMinimum(), 0, h->GetMaximum() * 1.5);
                zero->SetLineColor(kP10Red); zero->SetLineStyle(2); zero->SetLineWidth(2);
                zero->Draw("SAME");
                cTmp.cd(); DrawATLASLabel(); cTmp.SaveAs(path.c_str());
            };
            // Helper to draw the 5-algorithm overlay for a given set of residual histograms
            auto drawRelResidualOverlay = [&](std::vector<TH1F*> resH, const std::string& path,
                                              const std::string& xLbl = "(Truth - TOB) / Truth MET_{NonInt}",
                                              const std::string& yUnit = "") {
                int rcols[] = {kBlack, kP10Red, kP10Orange, kP10Blue, kP10Green};
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
                cOvl.cd(); DrawATLASLabel(); cOvl.SaveAs(path.c_str());
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
                zero->SetLineColor(kP10Red); zero->SetLineStyle(2); zero->SetLineWidth(2);
                zero->Draw("SAME");
                TLegend leg(0.20, 0.85, 0.58, 0.93);
                leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030);
                leg.AddEntry(zero, "Perfect resolution", "l");
                leg.AddEntry(prof, "Mean profile", "lp");
                leg.Draw();
                cTmp.cd(); DrawATLASLabel(); cTmp.SaveAs(path.c_str());
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
            BkgProcLabel bkgProcRes;   // everything below here in this block is background
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

        // --- gFEX algorithm comparison (JwoJ vs NoiseCut vs Rms), jFEX, and GEP types overlaid in kHz, up to 200 GeV threshold
        {
            std::vector<TH1F*> rateHists = hasOverlapRemoval
                ? std::vector<TH1F*>{back_hw_gMET, back_hw_gMET_NC, back_hw_gMET_Rms, back_hw_jMET, back_hw_JetMET, back_hw_TotalMET}
                : std::vector<TH1F*>{back_hw_gMET, back_hw_gMET_NC, back_hw_gMET_Rms, back_hw_jMET, back_hw_JetMET, back_hw_TowerMET, back_hw_TotalMET};
            std::vector<std::string> rateLabels = hasOverlapRemoval
                ? std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Total MET"}
                : std::vector<std::string>{"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms", "jFEX", "GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
            drawRateVsThresholdMulti(rateHists, rateLabels,
                                     "Rate vs MET threshold", "MET threshold [GeV]",
                                     fDir + "Rate_AlgoComparison_kHz_200GeV.pdf", "",
                                     "Estimated Background Rate [kHz]", 200.0, 1e-3, 1e-3);
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
            "Rate = 80 kHz", sig_h_metTruthNonInt_coarse, 0.52);
        drawTurnOnOverlay(
            allTOEffs(eff_gMET_60kHz, eff_gMET_NC_60kHz, eff_gMET_Rms_60kHz, eff_jMET_60kHz, eff_JetMET_60kHz, eff_TowerMET_60kHz, eff_TotalMET_60kHz),
            allTOLabels,
            "Turn-on at 60 kHz", fDir + "TurnOn_60kHz.pdf",
            allTOThrs(thr_gMET_60kHz, thr_gMET_NC_60kHz, thr_gMET_Rms_60kHz, thr_jMET_60kHz, thr_JetMET_60kHz, thr_TowerMET_60kHz, thr_TotalMET_60kHz),
            "Rate = 60 kHz", sig_h_metTruthNonInt_coarse, 0.52);
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
        drawTurnOnOverlay(
            {eff_gMET_60kHz, eff_gMET_NC_60kHz, eff_gMET_Rms_60kHz},
            {"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms"},
            "gFEX algorithm comparison — Turn-on at 60 kHz", fDir + "TurnOn_gFEX_AlgoComparison_60kHz.pdf",
            {thr_gMET_60kHz, thr_gMET_NC_60kHz, thr_gMET_Rms_60kHz}, "Rate = 60 kHz",
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
        drawTurnOnOverlay(
            gepTOEffs(eff_JetMET_60kHz, eff_TowerMET_60kHz, eff_TotalMET_60kHz),
            gepTOLabels,
            "GEP algorithm comparison — Turn-on at 60 kHz", fDir + "TurnOn_GEP_AlgoComparison_60kHz.pdf",
            gepTOThrs(thr_JetMET_60kHz, thr_TowerMET_60kHz, thr_TotalMET_60kHz), "Rate = 60 kHz",
            sig_h_metTruthNonInt_coarse);

        // --- Z->mumu turn-on curves vs dimuon p_{T}, rate matched to 40 / 60 / 80 kHz ---
        // One canvas per rate point with every MET type overlaid, plus a GEP-only canvas.
        if (isZmumuSample) {
            std::string muDir = fDir + "ZmumuTurnOn/";
            gSystem->mkdir(muDir.c_str(), true);
            // Which MET types go on the overlays: all of them, minus Tower MET for OR configs
            // where it carries no meaning (same convention as the GEP comparison plots above).
            std::vector<int> muIdx;
            for (int iA = 0; iA < nMETTypes; ++iA)
                if (!(hasOverlapRemoval && iA == towerMETTypeIdx)) muIdx.push_back(iA);
            for (int iR = 0; iR < nMuRates; ++iR) {
                std::vector<TH1F*>       effs, gepEffs;
                std::vector<std::string> lbls, gepLbls;
                std::vector<double>      thrs, gepThrs;
                for (int iA : muIdx) {
                    effs.push_back(effMu[iA][iR]);
                    lbls.push_back(metTypeLabel[iA]);
                    thrs.push_back(thrMu[iA][iR]);
                }
                for (int g = 0; g < nGEPMETTypes; ++g) {
                    const int iA = gepMETTypeIdx[g];
                    if (hasOverlapRemoval && iA == towerMETTypeIdx) continue;
                    gepEffs.push_back(effMu[iA][iR]);
                    gepLbls.push_back(metTypeLabel[iA]);
                    gepThrs.push_back(thrMu[iA][iR]);
                }
                const std::string rateLbl = std::string("Rate = ") + muRateNames[iR];
                drawTurnOnOverlay(effs, lbls,
                                  std::string("Z #rightarrow #mu#mu turn-on at ") + muRateNames[iR],
                                  muDir + "TurnOn_dimuonPt_" + muRateNames[iR] + ".pdf",
                                  thrs, rateLbl, sig_h_dimuonPt_coarse, 0.52, 0.15,
                                  "Dimuon p_{T} [GeV]");
                drawTurnOnOverlay(gepEffs, gepLbls,
                                  std::string("GEP algorithm comparison — Z #rightarrow #mu#mu turn-on at ") + muRateNames[iR],
                                  muDir + "TurnOn_dimuonPt_GEP_AlgoComparison_" + muRateNames[iR] + ".pdf",
                                  gepThrs, rateLbl, sig_h_dimuonPt_coarse, 0.45, 0.15,
                                  "Dimuon p_{T} [GeV]");
            }
            // The dimuon system itself, as a cross-check that it looks like a Z.
            drawMultiDist({sig_h_dimuonPt}, {"Dimuon system"},
                          "Dimuon p_{T}", "Dimuon p_{T} [GeV]", muDir + "dimuonPt.pdf");
            drawMultiDist({sig_h_dimuonMass}, {"Dimuon system"},
                          "Dimuon mass", "m_{#mu#mu} [GeV]", muDir + "dimuonMass.pdf", false);
        }

        // --- Background <MET> vs jet multiplicity ---
        if (hasTruthAntiKt4WZDressed || hasInTimeAntiKt4TruthJets) {
            BkgProcLabel bkgProc;   // background-only plot
            const std::string nJetXLabel =
                Form("N_{truth jets} (E_{T} > %.0f GeV, hard scatter + in-time pileup)", kNJetMinEt);
            std::vector<TProfile*>   allProfs;
            std::vector<std::string> allLbls;
            std::vector<TProfile*>   gepProfs;
            std::vector<std::string> gepLbls;
            for (int iA = 0; iA < nMETTypes; ++iA) {
                if (hasOverlapRemoval && iA == towerMETTypeIdx) continue;
                allProfs.push_back(back_prof_METvsNJets[iA]);
                allLbls.push_back(metTypeLabel[iA]);
            }
            for (int g = 0; g < nGEPMETTypes; ++g) {
                const int iA = gepMETTypeIdx[g];
                if (hasOverlapRemoval && iA == towerMETTypeIdx) continue;
                gepProfs.push_back(back_prof_METvsNJets[iA]);
                gepLbls.push_back(metTypeLabel[iA]);
            }
            drawProfileOverlay(allProfs, allLbls, nJetXLabel, "#LTMET#GT [GeV]",
                               fDir + "MET_vs_NJets_bkg.pdf");
            drawProfileOverlay(gepProfs, gepLbls, nJetXLabel, "#LTMET#GT [GeV]",
                               fDir + "MET_vs_NJets_GEP_bkg.pdf");
        }

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
            {eff_JetMET_60kHz,
             eff_combo_JwoJ_Jet_60kHz, eff_combo_NC_Jet_60kHz, eff_combo_Rms_Jet_60kHz},
            {"GEP Jet MET",
             std::string("JwoJ+Jet (")+comboThrLabel(best_JwoJ_Jet_60)+")",
             std::string("NC+Jet (")  +comboThrLabel(best_NC_Jet_60)+")",
             std::string("Rms+Jet (") +comboThrLabel(best_Rms_Jet_60)+")"},
            "GEP Jet MET + Combined Turn-on at 60 kHz", comboDir + "TurnOn_combined_Jet_60kHz.pdf",
            {}, "Rate = 60 kHz", sig_h_metTruthNonInt_coarse, 0.38, 0.15);
        drawTurnOnOverlay(
            {eff_TowerMET_80kHz,
             eff_combo_JwoJ_Tower_80kHz, eff_combo_NC_Tower_80kHz, eff_combo_Rms_Tower_80kHz},
            {"GEP Tower MET",
             std::string("JwoJ+Tower (")+comboThrLabel(best_JwoJ_Tower_80)+")",
             std::string("NC+Tower (")  +comboThrLabel(best_NC_Tower_80)+")",
             std::string("Rms+Tower (") +comboThrLabel(best_Rms_Tower_80)+")"},
            "GEP Tower MET + Combined Turn-on at 80 kHz", comboDir + "TurnOn_combined_Tower_80kHz.pdf",
            {}, "Rate = 80 kHz", sig_h_metTruthNonInt_coarse, 0.38, 0.15);
        drawTurnOnOverlay(
            {eff_TowerMET_60kHz,
             eff_combo_JwoJ_Tower_60kHz, eff_combo_NC_Tower_60kHz, eff_combo_Rms_Tower_60kHz},
            {"GEP Tower MET",
             std::string("JwoJ+Tower (")+comboThrLabel(best_JwoJ_Tower_60)+")",
             std::string("NC+Tower (")  +comboThrLabel(best_NC_Tower_60)+")",
             std::string("Rms+Tower (") +comboThrLabel(best_Rms_Tower_60)+")"},
            "GEP Tower MET + Combined Turn-on at 60 kHz", comboDir + "TurnOn_combined_Tower_60kHz.pdf",
            {}, "Rate = 60 kHz", sig_h_metTruthNonInt_coarse, 0.38, 0.15);

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
        styleRVE(out_RateVsEff_gMET_NC.gRate_vsEff,  kP10Red);
        styleRVE(out_RateVsEff_gMET_Rms.gRate_vsEff, kP10Orange);
        styleRVE(out_RateVsEff_jMET.gRate_vsEff,     kP10Cyan);
        styleRVE(out_RateVsEff_JetMET.gRate_vsEff,   kP10Blue);
        styleRVE(out_RateVsEff_TowerMET.gRate_vsEff, kP10Green);
        styleRVE(out_RateVsEff_TotalMET.gRate_vsEff, kP10Violet);

        // Individual PDFs
        auto drawRVEsingle = [&](TGraph* g, const std::string& title, const std::string& path) {
            TCanvas c("c", "", 700, 600);
            gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.14); gPad->SetTicks(1,1);
            gPad->SetLogy();
            g->SetTitle((title+";Signal Efficiency;Estimated Background Rate [Hz]").c_str());
            g->Draw("AP");
            c.cd(); DrawATLASLabel(); c.SaveAs(path.c_str());
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
            drawRateVsEffOverlay(allRVE, allRVELabels, (fDir + "RateVsEff_overlay.pdf"), "");
        }
        // gFEX-only algorithm comparison
        drawRateVsEffOverlay(
            {out_RateVsEff_gMET.gRate_vsEff, out_RateVsEff_gMET_NC.gRate_vsEff, out_RateVsEff_gMET_Rms.gRate_vsEff},
            {"gFEX JwoJ", "gFEX NoiseCut", "gFEX Rms"},
            (fDir + "RateVsEff_gFEX_AlgoComparison.pdf"), "");
        // GEP-only algorithm comparison
        {
            std::vector<TGraph*> gepRVE = hasOverlapRemoval
                ? std::vector<TGraph*>{out_RateVsEff_JetMET.gRate_vsEff, out_RateVsEff_TotalMET.gRate_vsEff}
                : std::vector<TGraph*>{out_RateVsEff_JetMET.gRate_vsEff, out_RateVsEff_TowerMET.gRate_vsEff, out_RateVsEff_TotalMET.gRate_vsEff};
            std::vector<std::string> gepRVELabels = hasOverlapRemoval
                ? std::vector<std::string>{"GEP Jet MET", "GEP Total MET"}
                : std::vector<std::string>{"GEP Jet MET", "GEP Tower MET", "GEP Total MET"};
            drawRateVsEffOverlay(gepRVE, gepRVELabels, (fDir + "RateVsEff_GEP_AlgoComparison.pdf"), "");
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
            (comboDir + "RateVsEff_combined_frontiers.pdf"), "");
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
            drawRateVsEffOverlay(cviRVE, cviLabels, (comboDir + "RateVsEff_combined_vs_individual.pdf"), "");
        }

        // --- AOD vs Sim gFEX MET comparison (only for resim ntuples) ---
        if (hasGFexSimMET) {
            // Distribution overlays: AOD vs Sim
            drawAlgoComparison(sig_h_gMET, back_h_gMET, sig_h_gMET_JwoJAOD, back_h_gMET_JwoJAOD,
                               "gFEX JwoJ (Sim)", "gFEX JwoJ (AOD)",
                               "MET [GeV]", fDir + "gFEX_JwoJ_AOD_vs_Sim.pdf", "");
            drawAlgoComparison(sig_h_gMET_NC, back_h_gMET_NC, sig_h_gMET_NCAOD, back_h_gMET_NCAOD,
                               "gFEX NoiseCut (Sim)", "gFEX NoiseCut (AOD)",
                               "MET [GeV]", fDir + "gFEX_NC_AOD_vs_Sim.pdf", "");
            drawAlgoComparison(sig_h_gMET_Rms, back_h_gMET_Rms, sig_h_gMET_RmsAOD, back_h_gMET_RmsAOD,
                               "gFEX Rms (Sim)", "gFEX Rms (AOD)",
                               "MET [GeV]", fDir + "gFEX_Rms_AOD_vs_Sim.pdf", "");
            // Rate vs threshold: AOD vs Sim
            drawRateVsThresholdMulti({back_hw_gMET,     back_hw_gMET_JwoJAOD},
                                     {"gFEX JwoJ (Sim)", "gFEX JwoJ (AOD)"},
                                     "Rate vs gFEX JwoJ MET (AOD vs Sim)", "MET threshold [GeV]",
                                     fDir + "Rate_gFEX_JwoJ_AOD_vs_Sim.pdf", "", "Estimated Background Rate [Hz]");
            drawRateVsThresholdMulti({back_hw_gMET_NC,   back_hw_gMET_NCAOD},
                                     {"gFEX NoiseCut (Sim)", "gFEX NoiseCut (AOD)"},
                                     "Rate vs gFEX NoiseCut MET (AOD vs Sim)", "MET threshold [GeV]",
                                     fDir + "Rate_gFEX_NC_AOD_vs_Sim.pdf", "", "Estimated Background Rate [Hz]");
            drawRateVsThresholdMulti({back_hw_gMET_Rms,  back_hw_gMET_RmsAOD},
                                     {"gFEX Rms (Sim)", "gFEX Rms (AOD)"},
                                     "Rate vs gFEX Rms MET (AOD vs Sim)", "MET threshold [GeV]",
                                     fDir + "Rate_gFEX_Rms_AOD_vs_Sim.pdf", "", "Estimated Background Rate [Hz]");
            // Turn-on curves: AOD vs Sim at 20/40/80 kHz
            drawTurnOnOverlay(
                {eff_gMET_80kHz, eff_gMET_JwoJAOD_80kHz},
                {"gFEX JwoJ (Sim)", "gFEX JwoJ (AOD)"},
                "gFEX JwoJ Turn-on at 80 kHz (AOD vs Sim)", fDir + "TurnOn_gFEX_JwoJ_AOD_vs_Sim_80kHz.pdf",
                {thr_gMET_80kHz, thr_gMET_JwoJAOD_80kHz}, "Rate = 80 kHz", sig_h_metTruthNonInt_coarse);
            drawTurnOnOverlay(
                {eff_gMET_60kHz, eff_gMET_JwoJAOD_60kHz},
                {"gFEX JwoJ (Sim)", "gFEX JwoJ (AOD)"},
                "gFEX JwoJ Turn-on at 60 kHz (AOD vs Sim)", fDir + "TurnOn_gFEX_JwoJ_AOD_vs_Sim_60kHz.pdf",
                {thr_gMET_60kHz, thr_gMET_JwoJAOD_60kHz}, "Rate = 60 kHz", sig_h_metTruthNonInt_coarse);
            drawTurnOnOverlay(
                {eff_gMET_NC_80kHz, eff_gMET_NCAOD_80kHz},
                {"gFEX NoiseCut (Sim)", "gFEX NoiseCut (AOD)"},
                "gFEX NoiseCut Turn-on at 80 kHz (AOD vs Sim)", fDir + "TurnOn_gFEX_NC_AOD_vs_Sim_80kHz.pdf",
                {thr_gMET_NC_80kHz, thr_gMET_NCAOD_80kHz}, "Rate = 80 kHz", sig_h_metTruthNonInt_coarse);
            drawTurnOnOverlay(
                {eff_gMET_NC_60kHz, eff_gMET_NCAOD_60kHz},
                {"gFEX NoiseCut (Sim)", "gFEX NoiseCut (AOD)"},
                "gFEX NoiseCut Turn-on at 60 kHz (AOD vs Sim)", fDir + "TurnOn_gFEX_NC_AOD_vs_Sim_60kHz.pdf",
                {thr_gMET_NC_60kHz, thr_gMET_NCAOD_60kHz}, "Rate = 60 kHz", sig_h_metTruthNonInt_coarse);
            drawTurnOnOverlay(
                {eff_gMET_Rms_80kHz, eff_gMET_RmsAOD_80kHz},
                {"gFEX Rms (Sim)", "gFEX Rms (AOD)"},
                "gFEX Rms Turn-on at 80 kHz (AOD vs Sim)", fDir + "TurnOn_gFEX_Rms_AOD_vs_Sim_80kHz.pdf",
                {thr_gMET_Rms_80kHz, thr_gMET_RmsAOD_80kHz}, "Rate = 80 kHz", sig_h_metTruthNonInt_coarse);
            drawTurnOnOverlay(
                {eff_gMET_Rms_60kHz, eff_gMET_RmsAOD_60kHz},
                {"gFEX Rms (Sim)", "gFEX Rms (AOD)"},
                "gFEX Rms Turn-on at 60 kHz (AOD vs Sim)", fDir + "TurnOn_gFEX_Rms_AOD_vs_Sim_60kHz.pdf",
                {thr_gMET_Rms_60kHz, thr_gMET_RmsAOD_60kHz}, "Rate = 60 kHz", sig_h_metTruthNonInt_coarse);
            // Rate vs efficiency: AOD vs Sim
            auto out_RVE_JwoJAOD = MakeRateVsEff(sig_h_gMET_JwoJAOD, back_hw_gMET_JwoJAOD);
            auto out_RVE_NCAOD   = MakeRateVsEff(sig_h_gMET_NCAOD,   back_hw_gMET_NCAOD);
            auto out_RVE_RmsAOD  = MakeRateVsEff(sig_h_gMET_RmsAOD,  back_hw_gMET_RmsAOD);
            styleRVE(out_RVE_JwoJAOD.gRate_vsEff, kBlack);
            styleRVE(out_RVE_NCAOD.gRate_vsEff,   kP10Red);
            styleRVE(out_RVE_RmsAOD.gRate_vsEff,  kP10Orange);
            drawRateVsEffOverlay(
                {out_RateVsEff_gMET.gRate_vsEff,    out_RVE_JwoJAOD.gRate_vsEff},
                {"gFEX JwoJ (Sim)", "gFEX JwoJ (AOD)"},
                (fDir + "RateVsEff_gFEX_JwoJ_AOD_vs_Sim.pdf"), "");
            drawRateVsEffOverlay(
                {out_RateVsEff_gMET_NC.gRate_vsEff, out_RVE_NCAOD.gRate_vsEff},
                {"gFEX NoiseCut (Sim)", "gFEX NoiseCut (AOD)"},
                (fDir + "RateVsEff_gFEX_NC_AOD_vs_Sim.pdf"), "");
            drawRateVsEffOverlay(
                {out_RateVsEff_gMET_Rms.gRate_vsEff, out_RVE_RmsAOD.gRate_vsEff},
                {"gFEX Rms (Sim)", "gFEX Rms (AOD)"},
                (fDir + "RateVsEff_gFEX_Rms_AOD_vs_Sim.pdf"), "");
            drawRateVsEffOverlay(
                {out_RateVsEff_gMET.gRate_vsEff,    out_RateVsEff_gMET_NC.gRate_vsEff,  out_RateVsEff_gMET_Rms.gRate_vsEff,
                 out_RVE_JwoJAOD.gRate_vsEff,        out_RVE_NCAOD.gRate_vsEff,           out_RVE_RmsAOD.gRate_vsEff},
                {"gFEX JwoJ (Sim)", "gFEX NoiseCut (Sim)", "gFEX Rms (Sim)",
                 "gFEX JwoJ (AOD)", "gFEX NoiseCut (AOD)", "gFEX Rms (AOD)"},
                (fDir + "RateVsEff_gFEX_AOD_vs_Sim_all.pdf"), "");
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
        if (hasSumJetET) back_hw_SumJetET_vec.push_back(cloneDetached(back_hw_SumJetET));
        eff_gMET_80kHz_vec.push_back(cloneDetached(eff_gMET_80kHz));
        eff_gMET_60kHz_vec.push_back(cloneDetached(eff_gMET_60kHz));
        eff_gMET_NC_80kHz_vec.push_back(cloneDetached(eff_gMET_NC_80kHz));
        eff_gMET_NC_60kHz_vec.push_back(cloneDetached(eff_gMET_NC_60kHz));
        eff_gMET_Rms_80kHz_vec.push_back(cloneDetached(eff_gMET_Rms_80kHz));
        eff_gMET_Rms_60kHz_vec.push_back(cloneDetached(eff_gMET_Rms_60kHz));
        eff_jMET_80kHz_vec.push_back(cloneDetached(eff_jMET_80kHz));
        eff_jMET_60kHz_vec.push_back(cloneDetached(eff_jMET_60kHz));
        eff_JetMET_80kHz_vec.push_back(cloneDetached(eff_JetMET_80kHz));
        eff_JetMET_60kHz_vec.push_back(cloneDetached(eff_JetMET_60kHz));
        eff_TowerMET_80kHz_vec.push_back(cloneDetached(eff_TowerMET_80kHz));
        eff_TowerMET_60kHz_vec.push_back(cloneDetached(eff_TowerMET_60kHz));
        eff_TotalMET_80kHz_vec.push_back(cloneDetached(eff_TotalMET_80kHz));
        eff_TotalMET_60kHz_vec.push_back(cloneDetached(eff_TotalMET_60kHz));
        thr_gMET_80kHz_vec.push_back(thr_gMET_80kHz);
        thr_gMET_60kHz_vec.push_back(thr_gMET_60kHz);
        thr_gMET_NC_80kHz_vec.push_back(thr_gMET_NC_80kHz);
        thr_gMET_NC_60kHz_vec.push_back(thr_gMET_NC_60kHz);
        thr_gMET_Rms_80kHz_vec.push_back(thr_gMET_Rms_80kHz);
        thr_gMET_Rms_60kHz_vec.push_back(thr_gMET_Rms_60kHz);
        thr_jMET_80kHz_vec.push_back(thr_jMET_80kHz);
        thr_jMET_60kHz_vec.push_back(thr_jMET_60kHz);
        thr_JetMET_80kHz_vec.push_back(thr_JetMET_80kHz);
        thr_JetMET_60kHz_vec.push_back(thr_JetMET_60kHz);
        thr_TowerMET_80kHz_vec.push_back(thr_TowerMET_80kHz);
        thr_TowerMET_60kHz_vec.push_back(thr_TowerMET_60kHz);
        thr_TotalMET_80kHz_vec.push_back(thr_TotalMET_80kHz);
        thr_TotalMET_60kHz_vec.push_back(thr_TotalMET_60kHz);
        sig_h_metTruthNonInt_unscaled_vec.push_back(cloneDetached(sig_h_metTruthNonInt_coarse));
        // Z->mumu turn-ons and the background jet-multiplicity profiles keep their own label
        // vectors: the former only exists for Zmumu inputs, the latter only when at least one
        // truth jet collection was present, so neither is guaranteed parallel to labels.
        if (isZmumuSample) {
            zmumuLabels.push_back(labels[fileIt]);
            sig_h_dimuonPt_coarse_vec.push_back(cloneDetached(sig_h_dimuonPt_coarse));
            for (int iA = 0; iA < nMETTypes; ++iA)
                for (int iR = 0; iR < nMuRates; ++iR) {
                    effMu_vec[iA][iR].push_back(cloneDetached(effMu[iA][iR]));
                    thrMu_vec[iA][iR].push_back(thrMu[iA][iR]);
                }
        }
        if (hasTruthAntiKt4WZDressed || hasInTimeAntiKt4TruthJets) {
            nJetProfLabels.push_back(labels[fileIt]);
            for (int iA = 0; iA < nMETTypes; ++iA) {
                TProfile* p = (TProfile*)back_prof_METvsNJets[iA]->Clone(
                    (std::string("back_prof_METvsNJets_multi_") + metTypeShort[iA] + "_" + tag).c_str());
                p->SetDirectory(0);
                back_prof_METvsNJets_vec[iA].push_back(p);
            }
        }
        if (hasGFexSimMET) {
            sig_h_gMET_JwoJAOD_vec.push_back(cloneDetached(sig_h_gMET_JwoJAOD));
            back_h_gMET_JwoJAOD_vec.push_back(cloneDetached(back_h_gMET_JwoJAOD));
            sig_h_gMET_NCAOD_vec.push_back(cloneDetached(sig_h_gMET_NCAOD));
            back_h_gMET_NCAOD_vec.push_back(cloneDetached(back_h_gMET_NCAOD));
            sig_h_gMET_RmsAOD_vec.push_back(cloneDetached(sig_h_gMET_RmsAOD));
            back_h_gMET_RmsAOD_vec.push_back(cloneDetached(back_h_gMET_RmsAOD));
            back_hw_gMET_JwoJAOD_vec.push_back(cloneDetached(back_hw_gMET_JwoJAOD));
            back_hw_gMET_NCAOD_vec.push_back(cloneDetached(back_hw_gMET_NCAOD));
            back_hw_gMET_RmsAOD_vec.push_back(cloneDetached(back_hw_gMET_RmsAOD));
            eff_gMET_JwoJAOD_80kHz_vec.push_back(cloneDetached(eff_gMET_JwoJAOD_80kHz));
            eff_gMET_JwoJAOD_60kHz_vec.push_back(cloneDetached(eff_gMET_JwoJAOD_60kHz));
            eff_gMET_NCAOD_80kHz_vec.push_back(cloneDetached(eff_gMET_NCAOD_80kHz));
            eff_gMET_NCAOD_60kHz_vec.push_back(cloneDetached(eff_gMET_NCAOD_60kHz));
            eff_gMET_RmsAOD_80kHz_vec.push_back(cloneDetached(eff_gMET_RmsAOD_80kHz));
            eff_gMET_RmsAOD_60kHz_vec.push_back(cloneDetached(eff_gMET_RmsAOD_60kHz));
            thr_gMET_JwoJAOD_80kHz_vec.push_back(thr_gMET_JwoJAOD_80kHz);
            thr_gMET_JwoJAOD_60kHz_vec.push_back(thr_gMET_JwoJAOD_60kHz);
            thr_gMET_NCAOD_80kHz_vec.push_back(thr_gMET_NCAOD_80kHz);
            thr_gMET_NCAOD_60kHz_vec.push_back(thr_gMET_NCAOD_60kHz);
            thr_gMET_RmsAOD_80kHz_vec.push_back(thr_gMET_RmsAOD_80kHz);
            thr_gMET_RmsAOD_60kHz_vec.push_back(thr_gMET_RmsAOD_60kHz);
        }
        sigF->Close();
        backF->Close();
    } // file loop

    // --- Multi-file overlays (sig + bkg overlaid, dashed=bkg) ---
    // Multi-file overlays keep the process name in the legend header (via signalName), so clear
    // the top-right per-file process label.
    gProcLabel = "";
    if (signalFiles.size() > 1) {
        std::string mDir = overlayDir;
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
        if (!back_hw_SumJetET_vec.empty())
            drawRateVsThresholdMulti(back_hw_SumJetET_vec, labels, "Rate vs GEP H_{T} threshold",       "H_{T} threshold [GeV]", mDir + "Rate_SumJetET.pdf",        signalName,
                                     "Rate [Hz]", /*xMax=*/-1.0);

        // --- Multi-file signal efficiency vs threshold ---
        drawEffVsThresholdMulti(sig_h_TotalMET_vec,  labels, "Signal Efficiency vs GEP Total MET Threshold",    "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_TotalMET.pdf",   signalName);
        drawEffVsThresholdMulti(sig_h_gMET_vec,      labels, "Signal Efficiency vs gFEX MET Threshold (JwoJ)",  "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_gFEX_JwoJ.pdf", signalName);
        drawEffVsThresholdMulti(sig_h_gMET_NC_vec,   labels, "Signal Efficiency vs gFEX MET Threshold (NC)",    "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_gFEX_NC.pdf",   signalName);
        drawEffVsThresholdMulti(sig_h_gMET_Rms_vec,  labels, "Signal Efficiency vs gFEX MET Threshold (Rms)",   "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_gFEX_Rms.pdf",  signalName);
        drawEffVsThresholdMulti(sig_h_jMET_vec,      labels, "Signal Efficiency vs jFEX MET Threshold",         "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_jFEX_MET.pdf",  signalName);
        drawEffVsThresholdMulti(sig_h_JetMet_vec,    labels, "Signal Efficiency vs GEP Jet MET Threshold",      "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_JetMET.pdf",    signalName);
        drawEffVsThresholdMulti(sig_h_TowerMet_vec,  labels, "Signal Efficiency vs GEP Tower MET Threshold",    "MET threshold [GeV]", mDir + "SigEff_vs_Threshold_TowerMET.pdf",  signalName);

        // --- Rate vs Efficiency multi-file overlays ---
        // One canvas per MET type, one curve per config: the GEP types and the FEX types, so
        // every algorithm that gets a Rate_*.pdf overlay above also gets a rate-vs-efficiency one.
        {
            std::vector<TGraph*> rveJet, rveTower, rveTotal;
            std::vector<TGraph*> rveGMET, rveGMET_NC, rveGMET_Rms, rveJMET;
            for (unsigned int i = 0; i < labels.size(); i++) {
                auto outJ = MakeRateVsEff(sig_h_JetMet_vec[i],   back_hw_JetMET_vec[i]);
                auto outT = MakeRateVsEff(sig_h_TowerMet_vec[i], back_hw_TowerMET_vec[i]);
                auto outTot = MakeRateVsEff(sig_h_TotalMET_vec[i], back_hw_TotalMET_vec[i]);
                auto outG   = MakeRateVsEff(sig_h_gMET_vec[i],     back_hw_gMET_vec[i]);
                auto outGNC = MakeRateVsEff(sig_h_gMET_NC_vec[i],  back_hw_gMET_NC_vec[i]);
                auto outGRms= MakeRateVsEff(sig_h_gMET_Rms_vec[i], back_hw_gMET_Rms_vec[i]);
                auto outjM  = MakeRateVsEff(sig_h_jMET_vec[i],     back_hw_jMET_vec[i]);
                rveJet.push_back((TGraph*)outJ.gRate_vsEff->Clone(Form("rve_JetMET_%u",   i)));
                rveTower.push_back((TGraph*)outT.gRate_vsEff->Clone(Form("rve_TowerMET_%u", i)));
                rveTotal.push_back((TGraph*)outTot.gRate_vsEff->Clone(Form("rve_TotalMET_%u", i)));
                rveGMET.push_back((TGraph*)outG.gRate_vsEff->Clone(Form("rve_gMET_%u",     i)));
                rveGMET_NC.push_back((TGraph*)outGNC.gRate_vsEff->Clone(Form("rve_gMET_NC_%u",  i)));
                rveGMET_Rms.push_back((TGraph*)outGRms.gRate_vsEff->Clone(Form("rve_gMET_Rms_%u", i)));
                rveJMET.push_back((TGraph*)outjM.gRate_vsEff->Clone(Form("rve_jMET_%u",     i)));
            }
            drawRateVsEffOverlay(rveJet,   labels, mDir + "RateVsEff_JetMET_multi.pdf",   signalName);
            drawRateVsEffOverlay(rveTower, labels, mDir + "RateVsEff_TowerMET_multi.pdf",  signalName);
            drawRateVsEffOverlay(rveTotal, labels, mDir + "RateVsEff_TotalMET_multi.pdf",  signalName);
            drawRateVsEffOverlay(rveGMET,     labels, mDir + "RateVsEff_gFEX_MET_JwoJ_multi.pdf",     signalName);
            drawRateVsEffOverlay(rveGMET_NC,  labels, mDir + "RateVsEff_gFEX_MET_NoiseCut_multi.pdf", signalName);
            drawRateVsEffOverlay(rveGMET_Rms, labels, mDir + "RateVsEff_gFEX_MET_Rms_multi.pdf",      signalName);
            drawRateVsEffOverlay(rveJMET,     labels, mDir + "RateVsEff_jFEX_MET_multi.pdf",          signalName);
            for (auto* g : rveJet)   delete g;
            for (auto* g : rveTower) delete g;
            for (auto* g : rveTotal) delete g;
            for (auto* g : rveGMET)     delete g;
            for (auto* g : rveGMET_NC)  delete g;
            for (auto* g : rveGMET_Rms) delete g;
            for (auto* g : rveJMET)     delete g;
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
        drawTurnOnOverlay(eff_gMET_60kHz_vec,
                          makeConfigLabels("gFEX JwoJ"),
                          "gFEX JwoJ Turn-on at 60 kHz",
                          mDir + "TurnOn_60kHz_gFEX_JwoJ_multi.pdf",
                          thr_gMET_60kHz_vec, "Rate = 60 kHz", truthOverlay);
        drawTurnOnOverlay(eff_gMET_NC_80kHz_vec,
                          makeConfigLabels("gFEX NoiseCut"),
                          "gFEX NoiseCut Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_gFEX_NoiseCut_multi.pdf",
                          thr_gMET_NC_80kHz_vec, "Rate = 80 kHz", truthOverlay);
        drawTurnOnOverlay(eff_gMET_NC_60kHz_vec,
                          makeConfigLabels("gFEX NoiseCut"),
                          "gFEX NoiseCut Turn-on at 60 kHz",
                          mDir + "TurnOn_60kHz_gFEX_NoiseCut_multi.pdf",
                          thr_gMET_NC_60kHz_vec, "Rate = 60 kHz", truthOverlay);
        drawTurnOnOverlay(eff_gMET_Rms_80kHz_vec,
                          makeConfigLabels("gFEX Rms"),
                          "gFEX Rms Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_gFEX_Rms_multi.pdf",
                          thr_gMET_Rms_80kHz_vec, "Rate = 80 kHz", truthOverlay);
        drawTurnOnOverlay(eff_gMET_Rms_60kHz_vec,
                          makeConfigLabels("gFEX Rms"),
                          "gFEX Rms Turn-on at 60 kHz",
                          mDir + "TurnOn_60kHz_gFEX_Rms_multi.pdf",
                          thr_gMET_Rms_60kHz_vec, "Rate = 60 kHz", truthOverlay);
        drawTurnOnOverlay(eff_jMET_80kHz_vec,
                          makeConfigLabels("jFEX"),
                          "jFEX MET Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_jFEX_MET_multi.pdf",
                          thr_jMET_80kHz_vec, "Rate = 80 kHz", truthOverlay);
        drawTurnOnOverlay(eff_jMET_60kHz_vec,
                          makeConfigLabels("jFEX"),
                          "jFEX MET Turn-on at 60 kHz",
                          mDir + "TurnOn_60kHz_jFEX_MET_multi.pdf",
                          thr_jMET_60kHz_vec, "Rate = 60 kHz", truthOverlay);
        drawTurnOnOverlay(eff_JetMET_80kHz_vec,
                          makeConfigLabels("GEP Jet MET"),
                          "GEP Jet MET Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_JetMET_multi.pdf",
                          thr_JetMET_80kHz_vec, "Rate = 80 kHz", truthOverlay);
        drawTurnOnOverlay(eff_JetMET_60kHz_vec,
                          makeConfigLabels("GEP Jet MET"),
                          "GEP Jet MET Turn-on at 60 kHz",
                          mDir + "TurnOn_60kHz_JetMET_multi.pdf",
                          thr_JetMET_60kHz_vec, "Rate = 60 kHz", truthOverlay);
        drawTurnOnOverlay(eff_TowerMET_80kHz_vec,
                          makeConfigLabels("GEP Tower MET"),
                          "GEP Tower MET Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_TowerMET_multi.pdf",
                          thr_TowerMET_80kHz_vec, "Rate = 80 kHz", truthOverlay,
                          0.45, 0.20);
        drawTurnOnOverlay(eff_TowerMET_60kHz_vec,
                          makeConfigLabels("GEP Tower MET"),
                          "GEP Tower MET Turn-on at 60 kHz",
                          mDir + "TurnOn_60kHz_TowerMET_multi.pdf",
                          thr_TowerMET_60kHz_vec, "Rate = 60 kHz", truthOverlay,
                          0.45, 0.20);
        drawTurnOnOverlay(eff_TotalMET_80kHz_vec,
                          makeConfigLabels("GEP Total MET"),
                          "GEP Total MET Turn-on at 80 kHz",
                          mDir + "TurnOn_80kHz_TotalMET_multi.pdf",
                          thr_TotalMET_80kHz_vec, "Rate = 80 kHz", truthOverlay, 0.38, 0.18);
        drawTurnOnOverlay(eff_TotalMET_60kHz_vec,
                          makeConfigLabels("GEP Total MET"),
                          "GEP Total MET Turn-on at 60 kHz",
                          mDir + "TurnOn_60kHz_TotalMET_multi.pdf",
                          thr_TotalMET_60kHz_vec, "Rate = 60 kHz", truthOverlay, 0.38, 0.18);

        // --- Z->mumu dimuon-p_{T} turn-on across configs, one canvas per MET type per rate ---
        if (zmumuLabels.size() > 1) {
            auto makeZmumuLabels = [&](const std::string& algoName) {
                std::vector<std::string> v;
                for (const auto& lbl : zmumuLabels) v.push_back(algoName + " (" + lbl + ")");
                return v;
            };
            TH1F* dimuonOverlay = sig_h_dimuonPt_coarse_vec.empty() ? nullptr : sig_h_dimuonPt_coarse_vec[0];
            for (int iA = 0; iA < nMETTypes; ++iA)
                for (int iR = 0; iR < nMuRates; ++iR)
                    drawTurnOnOverlay(effMu_vec[iA][iR], makeZmumuLabels(metTypeLabel[iA]),
                                      std::string(metTypeLabel[iA]) + " Z #rightarrow #mu#mu turn-on at " + muRateNames[iR],
                                      mDir + "TurnOn_dimuonPt_" + muRateNames[iR] + "_" + metTypeShort[iA] + "_multi.pdf",
                                      thrMu_vec[iA][iR],
                                      std::string("Rate = ") + muRateNames[iR], dimuonOverlay,
                                      0.45, 0.18, "Dimuon p_{T} [GeV]");
        }

        // --- Background <MET> vs jet multiplicity across configs, one canvas per MET type ---
        if (nJetProfLabels.size() > 1) {
            BkgProcLabel bkgProc;   // background-only plots: label as QCD dijet, not the signal
            const std::string nJetXLabel =
                Form("N_{truth jets} (E_{T} > %.0f GeV, hard scatter + in-time pileup)", kNJetMinEt);
            for (int iA = 0; iA < nMETTypes; ++iA)
                drawProfileOverlay(back_prof_METvsNJets_vec[iA], nJetProfLabels,
                                   nJetXLabel,
                                   Form("#LT%s#GT [GeV]", metTypeLabel[iA]),
                                   mDir + "MET_vs_NJets_" + metTypeShort[iA] + "_multi.pdf");
        }

        // Standalone multi-file overlays of the AOD gFEX (the *_*Sim histograms now hold AOD,
        // since nominal gFEX was promoted to resim above). Labelled AOD accordingly.
        if (!sig_h_gMET_JwoJAOD_vec.empty()) {
            drawOverlayMulti(sig_h_gMET_JwoJAOD_vec, back_h_gMET_JwoJAOD_vec, labels, "gFEX MET (JwoJ AOD)",    "MET [GeV]", mDir + "gFEX_MET_JwoJAOD.pdf",   signalName);
            drawOverlayMulti(sig_h_gMET_NCAOD_vec,   back_h_gMET_NCAOD_vec,   labels, "gFEX MET (NoiseCut AOD)","MET [GeV]", mDir + "gFEX_MET_NCAOD.pdf",     signalName);
            drawOverlayMulti(sig_h_gMET_RmsAOD_vec,  back_h_gMET_RmsAOD_vec,  labels, "gFEX MET (Rms AOD)",     "MET [GeV]", mDir + "gFEX_MET_RmsAOD.pdf",    signalName);
            drawRateVsThresholdMulti(back_hw_gMET_JwoJAOD_vec, labels, "Rate vs gFEX MET (JwoJ AOD)",    "MET threshold [GeV]", mDir + "Rate_gFEX_MET_JwoJAOD.pdf", signalName);
            drawRateVsThresholdMulti(back_hw_gMET_NCAOD_vec,   labels, "Rate vs gFEX MET (NoiseCut AOD)","MET threshold [GeV]", mDir + "Rate_gFEX_MET_NCAOD.pdf",   signalName);
            drawRateVsThresholdMulti(back_hw_gMET_RmsAOD_vec,  labels, "Rate vs gFEX MET (Rms AOD)",     "MET threshold [GeV]", mDir + "Rate_gFEX_MET_RmsAOD.pdf",  signalName);
            drawTurnOnOverlay(eff_gMET_JwoJAOD_80kHz_vec, makeConfigLabels("gFEX JwoJ AOD"),
                              "gFEX JwoJ AOD Turn-on at 80 kHz",
                              mDir + "TurnOn_80kHz_gFEX_JwoJAOD_multi.pdf",
                              thr_gMET_JwoJAOD_80kHz_vec, "Rate = 80 kHz", truthOverlay);
            drawTurnOnOverlay(eff_gMET_JwoJAOD_60kHz_vec, makeConfigLabels("gFEX JwoJ AOD"),
                              "gFEX JwoJ AOD Turn-on at 60 kHz",
                              mDir + "TurnOn_60kHz_gFEX_JwoJAOD_multi.pdf",
                              thr_gMET_JwoJAOD_60kHz_vec, "Rate = 60 kHz", truthOverlay);
            drawTurnOnOverlay(eff_gMET_NCAOD_80kHz_vec, makeConfigLabels("gFEX NoiseCut AOD"),
                              "gFEX NoiseCut AOD Turn-on at 80 kHz",
                              mDir + "TurnOn_80kHz_gFEX_NCAOD_multi.pdf",
                              thr_gMET_NCAOD_80kHz_vec, "Rate = 80 kHz", truthOverlay);
            drawTurnOnOverlay(eff_gMET_NCAOD_60kHz_vec, makeConfigLabels("gFEX NoiseCut AOD"),
                              "gFEX NoiseCut AOD Turn-on at 60 kHz",
                              mDir + "TurnOn_60kHz_gFEX_NCAOD_multi.pdf",
                              thr_gMET_NCAOD_60kHz_vec, "Rate = 60 kHz", truthOverlay);
            drawTurnOnOverlay(eff_gMET_RmsAOD_80kHz_vec, makeConfigLabels("gFEX Rms AOD"),
                              "gFEX Rms AOD Turn-on at 80 kHz",
                              mDir + "TurnOn_80kHz_gFEX_RmsAOD_multi.pdf",
                              thr_gMET_RmsAOD_80kHz_vec, "Rate = 80 kHz", truthOverlay);
            drawTurnOnOverlay(eff_gMET_RmsAOD_60kHz_vec, makeConfigLabels("gFEX Rms AOD"),
                              "gFEX Rms AOD Turn-on at 60 kHz",
                              mDir + "TurnOn_60kHz_gFEX_RmsAOD_multi.pdf",
                              thr_gMET_RmsAOD_60kHz_vec, "Rate = 60 kHz", truthOverlay);
        }
    }

    std::cout << "Done. Plots saved to " << outputDir << "\n";
}

// -----------------------------------------------------------------------
void metAnalysisAndRates() {
    gErrorIgnoreLevel = kError;
    SetPlotStyle();

    // Each entry: { HERNTupler input ntuple, MET emulator output }
    // The input ntuple is per-process; set the .first of each signalFiles entry to the
    // matching process input below so multiple signal processes can be analysed together.
    const std::string sigInput  = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/ZvvHbb_v4/mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    const std::string sigInputTtbarSemilep = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/ttbar_semilep_v4/mc21_14TeV_ttbar_hdamp258p75_semilep_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    const std::string sigInputTtbarDilep   = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/ttbar_dilep_v4/mc21_14TeV_ttbar_hdamp258p75_dilep_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    // Z->mumu: the sample the dimuon-p_{T} turn-ons need. HERNTupler only retrieves TruthMuons
    // (and so only fills dimuonPt) for this sample, and analyze_files keys the dimuon plots off
    // "Zmumu" appearing in the path — so this input must be used, not a copy under another name.
    const std::string sigInputZmumu = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/Zmumu_v4/mc21_14TeV_Zmumu_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    // Glob over the ten per-slice ntuples rather than a hadd'd combination: the merged
    // file would exceed ROOT's 100 GB TTree::fgMaxTreeSize and come out truncated.
    // ChainSource turns this into one TChain per tree (see chainSource.h).
    const std::string backInput = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/QCD_Dijet_JZ*_v4/mc21_14TeV_jj_JZ*_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    // PU140 counterparts of the inputs above: same samples in the ntuples_PU140 mirror tree,
    // r16129 in place of r16130 (see submit_all_ntupler.sh --pu 140). An r16129 emulator output
    // MUST be paired with the r16129 input here — pairing it with the PU200 input silently
    // compares the PU140 emulation against PU200 truth/AOD quantities.
    const std::string sigInputPU140  = "/data/larsonma/GEPHadronicEventReconstruction/ntuples_PU140/ZvvHbb_v4/mc21_14TeV_ZvvH125_bb_e8557_s4422_r16129_DAOD_NTUPLE_GEP.root";
    const std::string sigInputTtbarSemilepPU140 = "/data/larsonma/GEPHadronicEventReconstruction/ntuples_PU140/ttbar_semilep_v4/mc21_14TeV_ttbar_hdamp258p75_semilep_e8557_s4422_r16129_DAOD_NTUPLE_GEP.root";
    const std::string sigInputTtbarDilepPU140   = "/data/larsonma/GEPHadronicEventReconstruction/ntuples_PU140/ttbar_dilep_v4/mc21_14TeV_ttbar_hdamp258p75_dilep_e8557_s4422_r16129_DAOD_NTUPLE_GEP.root";
    const std::string sigInputZmumuPU140 = "/data/larsonma/GEPHadronicEventReconstruction/ntuples_PU140/Zmumu_v4/mc21_14TeV_Zmumu_e8557_s4422_r16129_DAOD_NTUPLE_GEP.root";
    const std::string backInputPU140 = "/data/larsonma/GEPHadronicEventReconstruction/ntuples_PU140/QCD_Dijet_JZ*_v4/mc21_14TeV_jj_JZ*_e8557_s4422_r16129_DAOD_NTUPLE_GEP.root";
    const std::string emuDir    = "/data/larsonma/GEPMET/outputNTuplesDev_METv2/";

    // The trailing _twrSF{x}_jetSF{y} tag encodes the (tower, jet) scale factors
    // applied in the totalMET sum (see makeOutputMETFileName). Compare 1,1 vs 1,0p4
    // to see the effect of down-weighting the jet contribution. Future: per-eta calibrated SFs.
    std::vector<std::pair<std::string, std::string>> signalFiles = {
        /*{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt10_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt10_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },*/

        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_SK_NoOR_twrSF1_jetSF1.root"   }

       // { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt10_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },

        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF0p4_jetSF1.root"   },
        //{ sigInputTtbarSemilep, emuDir + "mc21_14TeV_ttbar_hdamp258p75_semilep_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        //{ sigInputTtbarSemilep, emuDir + "mc21_14TeV_ttbar_hdamp258p75_semilep_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF0p4_jetSF1.root"   },
        //{ sigInputTtbarDilep,   emuDir + "mc21_14TeV_ttbar_hdamp258p75_dilep_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        //{ sigInputTtbarDilep,   emuDir + "mc21_14TeV_ttbar_hdamp258p75_dilep_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF0p4_jetSF1.root"   },


        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_SK_NoOR_twrSF1_jetSF1.root"   },

        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root"   },


        // For new OR studies 08122026
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root"   },
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF0p5.root"   },

        // For new PU studies 08132026
        //{ sigInputPU140, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt0_NoSK_NoOR_twrSF1_jetSF1.root"   },
        //{ sigInputPU140, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt0_EtaSK_NoOR_twrSF1_jetSF1.root"   },
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_NoSK_NoOR_twrSF1_jetSF1.root"   },
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_EtaSK_NoOR_twrSF1_jetSF1.root"   },


        // PU 200 + varying tower coefficient studies 08182026
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p2_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p3_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p5_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p6_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p7_jetSF1.root"   },
        { sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },

        // PU 200 + varying tower coefficient studies 08182026
        /*{ sigInputPU140, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p2_jetSF1.root"   },
        { sigInputPU140, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p3_jetSF1.root"   },
        { sigInputPU140, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        { sigInputPU140, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p5_jetSF1.root"   },
        { sigInputPU140, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p6_jetSF1.root"   },
        { sigInputPU140, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p7_jetSF1.root"   },
        { sigInputPU140, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },*/

        // Z->mumu, same four configs as the PU study above (PU140/PU200 x NoSK/EtaSK). This is
        // the sample the dimuon-p_{T} turn-ons come from — they are only produced for inputs
        // with "Zmumu" in the path, since HERNTupler only fills dimuonPt for this sample.
        // Uncomment these (and the matching labels / signalName / signalNames below) to run it.
        //{ sigInputZmumuPU140, emuDir + "mc21_14TeV_Zmumu_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt0_NoSK_NoOR_twrSF1_jetSF1.root"   },
        //{ sigInputZmumuPU140, emuDir + "mc21_14TeV_Zmumu_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt0_EtaSK_NoOR_twrSF1_jetSF1.root"  },
        //{ sigInputZmumu,      emuDir + "mc21_14TeV_Zmumu_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_NoSK_NoOR_twrSF1_jetSF1.root"   },
        //{ sigInputZmumu,      emuDir + "mc21_14TeV_Zmumu_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_EtaSK_NoOR_twrSF1_jetSF1.root"  },

        //{ sigInputTtbarSemilep, emuDir + "mc21_14TeV_ttbar_hdamp258p75_semilep_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root" },
        //{ sigInputTtbarDilep,   emuDir + "mc21_14TeV_ttbar_hdamp258p75_dilep_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root"   },

        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        //{ sigInputTtbarSemilep, emuDir + "mc21_14TeV_ttbar_hdamp258p75_semilep_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
        //{ sigInputTtbarDilep,   emuDir + "mc21_14TeV_ttbar_hdamp258p75_dilep_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },


        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_SK_OR_twrSF1_jetSF1.root"   },
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_NoSK_OR_twrSF1_jetSF1.root"   },
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_SK_OR_twrSF1_jetSF1.root"   },
        ///{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_SK_OR_twrSF1_jetSF1.root"   },
        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_NoSK_OR_twrSF1_jetSF1.root"   },

        //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_SK_OR_twrSF0p4_jetSF1.root"   },
       //{ sigInput, emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_SK_OR_twrSF1_jetSF1.root"   },
    };

    std::vector<std::pair<std::string, std::string>> backgroundFiles = {
        /*{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt10_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt10_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root" },*/

        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_SK_NoOR_twrSF1_jetSF1.root"   },

       // { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt10_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },

        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF0p4_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF0p4_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF0p4_jetSF1.root"   },

        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root"   },

        // For new OR studies 08122026
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF0p5.root"   },

        // For new PU studies 08132026
        //{ backInputPU140, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt0_NoSK_NoOR_twrSF1_jetSF1.root"   },
        //{ backInputPU140, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt0_EtaSK_NoOR_twrSF1_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_NoSK_NoOR_twrSF1_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_EtaSK_NoOR_twrSF1_jetSF1.root"   },
        // The Z->mumu block above pairs with these same four background entries — the emulator
        // configs and the pileup scenarios line up one for one, so nothing changes on this side.

        // PU 200 + varying tower coefficient studies 08182026
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p2_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p3_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p5_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p6_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p7_jetSF1.root"   },
        { backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },

        // PU 140 + varying tower coefficient studies 08182026
        /*{ backInputPU140, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p2_jetSF1.root"   },
        { backInputPU140, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p3_jetSF1.root"   },
        { backInputPU140, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        { backInputPU140, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p5_jetSF1.root"   },
        { backInputPU140, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p6_jetSF1.root"   },
        { backInputPU140, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p7_jetSF1.root"   },
        { backInputPU140, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16129_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF1_jetSF1.root"   },*/

        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root"   },
        // ttbar semileptonic — same dijet background (process-independent), matching emu config
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root"   },
        // ttbar dileptonic — same dijet background (process-independent), matching emu config
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root"   },


        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        // ttbar semileptonic — same dijet background (process-independent), matching emu config
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },
        // ttbar dileptonic — same dijet background (process-independent), matching emu config
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_OR_twrSF0p4_jetSF1.root"   },


        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_SK_OR_twrSF1_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_NoSK_OR_twrSF1_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_EtaSK_OR_twrSF1_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_SK_OR_twrSF1_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt0_NoSK_OR_twrSF1_jetSF1.root"   },

        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_SK_OR_twrSF0p4_jetSF1.root"   },
        //{ backInput, emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_SK_OR_twrSF1_jetSF1.root"   },
    };

    // Label for each file pair — shown in legend for multi-file overlays - NEEDS TO BE UPDATED AND THE SAME SIZE AS THE FILES BEING PROCESSED!
    std::vector<std::string> labels = {
        /*"J10_T2_EtaSK_OR_TC1_JC1",
        "J10_T2_EtaSK_OR_TC0.4_JC1",
        "J15_T2_EtaSK_OR_TC1_JC1",
        "J15_T2_EtaSK_OR_TC0.4_JC1",
        "J20_T2_EtaSK_OR_TC1_JC1",
        "J20_T2_EtaSK_OR_TC0.4_JC1",*/
        //"NoOR_TC1_JC1",
        //"OR_TC1_JC1",
        //"NoOR_TC1_JC0.5",

        //"J15_T1_NoSK_NoOR_PU140",
        //"J15_T1_EtaSK_NoOR_PU140",
        //"J15_T1_NoSK_NoOR_PU200",
        //"J15_T1_EtaSK_NoOR_PU200",

        "TC_JC",

        // Labels for the commented-out Z->mumu block above (same four configs).
        //"Zmumu_J15_T1_NoSK_NoOR_PU140",
        //"Zmumu_J15_T1_EtaSK_NoOR_PU140",
        //"Zmumu_J15_T1_NoSK_NoOR_PU200",
        //"Zmumu_J15_T1_EtaSK_NoOR_PU200",
        //"J10_EtaSK",
        //"J15_EtaSK",
        //"J20_EtaSK",

        //"J15_T2_EtaSK_OR_0.4S_1H",
        //"J15_T2_EtaSK_OR_1S_1H",

        //"ZvvHbb_J15_T2_OR_0.4S_1H",
        //"ZvvHbb_J15_T2_NoOR_0.4S_1H",
        //"tt_1lep_J15_T2_OR_0.4S_1H",
        //"tt_1lep_J15_T2_NoOR_0.4S_1H",
        //"tt_2lep_J15_T2_OR_0.4S_1H",
        //"tt_2lep_J15_T2_NoOR_0.4S_1H",

        //"ZvvHbb_J15_T2_EtaSK_OR_0.4S_1H",
        //"tt_1lep_J15_T2_EtaSK_OR_0.4S_1H",
        //"tt_2lep_J15_T2_EtaSK_OR_0.4S_1H",

        //"J15_T2_SK_OR",
        //"J15_T2_NoSK_OR",
        //"J15_T0_EtaSK_OR",
        //"J15_T0_SK_OR",
        //"J15_T0_NoSK_OR",

        //"J15_T2_SK_OR_TC0.4_JC1",
        //"J15_T2_SK_OR_TC1_JC1",
    };

    std::string signalName = "Z #rightarrow #nu#bar{#nu}, H #rightarrow b#bar{b}";
    //std::string signalName = "Z #rightarrow #mu#mu";     // Z->mumu block above
    // Per-file legend header (parallel to signalFiles/labels) so each signal process gets
    // its own label on per-file plots. Falls back to signalName for any unset entry.
    std::vector<std::string> signalNames = {
        "Z #rightarrow #nu#bar{#nu}, H #rightarrow b#bar{b}",
        //"Z #rightarrow #mu#mu",
        //"t#bar{t} semileptonic decay",
        //"t#bar{t} dileptonic decay",
    };
    std::string outputDir  = "metAnalysisPlots/";
    std::string overlayDir = "multiFileOverlay_MET_CoeffComparison_ZvvHbb_PU200/";
    //std::string overlayDir = "multiFileOverlay_MET_Zmumu_PUComparison/";   // Z->mumu block above
    gSystem->mkdir(outputDir.c_str(), true);

    analyze_files(signalFiles, backgroundFiles, labels, outputDir, signalName, overlayDir, signalNames);
}
