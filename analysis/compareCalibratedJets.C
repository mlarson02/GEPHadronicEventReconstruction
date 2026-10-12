// Compare the JES-calibrated WTACone CellsTowers jets with the uncalibrated ones, for one signal and
// the full QCD dijet JZ0-JZ9 background of one pileup scenario. Needs HERNTupler ntuples produced with
// apply_gep_jes_calibration_ = true, which hold the *CalibJetsTree trees.
//
// For each PU-suppression variant (none / SK / EtaSK), uncalibrated vs calibrated:
//   - E_T distributions of all jets, of the leading and of the subleading jet
//   - mean number of jets per event above an E_T threshold
//   - leading-jet rate (background, Hz) and signal efficiency vs threshold
//   - calibration factor E_T(calib) / E_T(raw) vs raw E_T and vs eta
//   - closure: mean E_T / E_T(truth) vs truth E_T for jets matched to AntiKt4TruthJets with the
//     selection used to derive the LUTs, applied in truth E_T. The LUTs were derived against truth
//     p_T (HERCalibHistMaker), so the calibrated curve sits slightly below 1 where the truth jet
//     mass matters.
//   - leading WTACone jet E_T vs leading truth jet E_T, and vs leading HLT AntiKt4EMTopo jet E_T
//     (TH2F, uncalibrated and calibrated), with the slope of the least-squares line through the
//     events and their correlation coefficient r. The HLT jets are the only calibrated small-R reco
//     jets in the ntuples; there is no offline AntiKt4 collection.
//   - rate-matched signal turn-ons as in largeRJetAnalysisCore.h: leading jet at 25 kHz and 4th
//     leading jet at 100 kHz of background rate, vs the leading / 4th leading AntiKt4TruthDressedWZ
//     jet E_T. Each curve has its own threshold for that rate (printed); events with too few jets
//     enter at E_T = 0, as in the core. The core's cone-jet turn-ons use the NoSK variant.
// Truth jet E_T = E / cosh(eta), from the stored p_T, eta and mass.
// Background events are weighted as in largeRJetAnalysisAndRates.C / metAnalysisAndRates.C:
// eventWeights->at(0) with passHSTP required; signal is unweighted. The E_T distributions and the
// TH2Fs are normalized to unit area and drawn down to 1e-10 to show the high-E_T tail; the jet
// multiplicity is the mean per event and the leading-jet rate stays in Hz.
//
// The printout breaks the background down by JZ slice (events, HSTP pass fractions, share of the
// weight, mean jet multiplicity) and reports, per variant, the jets below the 10 GeV lower edge of
// the LUT derivation, which the calibration scales with the response of the 10-20 GeV bin.
//
// Output: one PDF per plot, gepJetCalibComparison/PU<pileup>_<signal>/<NoSK|SK|EtaSK>/<sig|back>_<plot>.pdf
//
// Run: root -b -q 'compareCalibratedJets.C+(140, "ggF_HHbbbb")'     signal: ggF_HHbbbb or ttbar_allhad

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TColor.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TObjArray.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TSystem.h"

const std::string calib_comparison_ntuple_base_ = "/data/larsonma/GEPHadronicEventReconstruction/ntuples";
const std::string calib_comparison_output_dir_ = "gepJetCalibComparison";
const unsigned int calib_comparison_algo_version_ = 4;

const int n_calib_variants_ = 3;
const std::string calib_variant_tags_[n_calib_variants_] = {"", "SK", "EtaSK"}; // Tree name infix
const std::string calib_variant_labels_[n_calib_variants_] = {"No PU suppression", "SK", "EtaSK"};
const std::string calib_variant_dirs_[n_calib_variants_] = {"NoSK", "SK", "EtaSK"}; // Output subdirectory

// Thresholds for the jet multiplicity, at the TProfile bin centers: 0, 5, ..., 245 GeV
const int n_multiplicity_thresholds_ = 50;
const double multiplicity_threshold_step_ = 5.0;

// Lower edge of the log axis of the unit-normalized distributions
const double normalized_log_min_ = 1e-10;

// Lowest raw E_T of the jets the LUTs were derived from (HERCalibHistMaker). CalibLUTMaker gives the
// empty 0-10 GeV LUT bin the response of the nearest filled bin, so softer jets get the 10-20 GeV factor.
const double lut_derivation_min_et_ = 10.0;

// Background rates the turn-on thresholds are set at, as in largeRJetAnalysisCore.h
const double single_jet_rate_hz_ = 25000.0;  // Leading jet
const double fourth_jet_rate_hz_ = 100000.0; // 4th leading jet

// Per-JZ-slice breakdown of the background multiplicity, no-PU-suppression variant
const int n_jz_slices_ = 10;
const int n_slice_thresholds_ = 4;
const double slice_thresholds_[n_slice_thresholds_] = {0.0, 5.0, 10.0, 20.0};

struct CalibComparisonHists {
    // [0] uncalibrated, [1] calibrated
    TH1D* allEt[2];
    TH1D* leadEt[2];
    TH1D* subleadEt[2];
    TProfile* nAbove[2];
    TProfile* factorVsEt;
    TProfile* factorVsEta;
    TProfile* closure[2][2]; // [|eta(truth)| < 2.5, 2.5-4.9][uncalibrated, calibrated]
    TH2F* leadVsTruth[2];    // Leading jet E_T vs leading truth jet E_T
    TH2F* leadVsHLT[2];      // Leading jet E_T vs leading HLT AntiKt4EMTopo jet E_T
    TH1D* leadResidual[2];   // (E_T - E_T^truth) / E_T^truth of the leading jet, as in ResHistMaker
    TH1D* leadResidualUnweighted[2]; // The same, every event with weight 1 (differs from the above for the background only)
    TH1D* fourthResidual[2];           // The same for the 4th leading jet
    TH1D* fourthResidualUnweighted[2];
    TH1D* leadEtFine[2];     // Leading / 4th leading jet E_T (0 without such a jet), finely binned,
    TH1D* fourthEtFine[2];   // for the rate-matched thresholds
};

// One signal event for the rate-matched turn-ons: leading / 4th leading AntiKt4TruthDressedWZ jet E_T
// (-1 without such a jet) and WTACone jet E_T per variant, [uncalibrated, calibrated] (0 without such a jet)
struct CalibComparisonTurnOnEvent {
    double truthLead = -1.0;
    double truth4th = -1.0;
    double lead[n_calib_variants_][2];
    double fourth[n_calib_variants_][2];
};

// Background events of one JZ slice
struct CalibComparisonSlice {
    long nEvents = 0;
    long nHSTP = 0;
    double sumW = 0.0;                                 // All events
    double sumWHSTP = 0.0;                             // passHSTP events
    double sumWNJets[2][n_slice_thresholds_] = {};     // [uncalibrated, calibrated] sum of w * jets above threshold, passHSTP
};

// Jets of one variant below the LUT derivation range, as sums of event weights
struct CalibComparisonLowEt {
    double sumWJets = 0.0;         // All jets with raw E_T > 0
    double sumWLow = 0.0;          // Raw E_T < lut_derivation_min_et_
    double sumWLowPromoted = 0.0;  // Of those, calibrated to E_T >= lut_derivation_min_et_
};

struct CalibComparisonSample {
    std::string label;
    bool isBackground;
    double sumWeights = 0.0;
    long nEvents = 0;
    CalibComparisonHists hists[n_calib_variants_];
    CalibComparisonSlice slices[n_jz_slices_];
    CalibComparisonLowEt lowEt[n_calib_variants_];
    std::vector<CalibComparisonTurnOnEvent> turnOnEvents; // Signal only
};

double calibComparisonDeltaR(double eta1, double phi1, double eta2, double phi2) {
    const double deta = eta1 - eta2;
    const double dphi = std::remainder(phi1 - phi2, 2.0 * M_PI);
    return std::sqrt(deta * deta + dphi * dphi);
}

// Truth jet E_T = E / cosh(eta), with E^2 = (p_T cosh(eta))^2 + m^2
double calibComparisonTruthEt(double pt, double eta, double mass) {
    const double c = std::cosh(eta);
    return std::sqrt(pt * pt + mass * mass / (c * c));
}

void bookCalibComparisonHists(CalibComparisonHists& h, const std::string& tag) {
    const char* kind[2] = {"raw", "calib"};
    const double closureEtEdges[] = {20, 30, 40, 50, 60, 80, 100, 140, 200, 300, 500, 1000};
    const int nClosureBins = sizeof(closureEtEdges) / sizeof(closureEtEdges[0]) - 1;
    for (int k = 0; k < 2; ++k) {
        h.allEt[k] = new TH1D(Form("allEt_%s_%s", tag.c_str(), kind[k]), ";Jet E_{T} [GeV];Jets", 100, 0, 500);
        h.leadEt[k] = new TH1D(Form("leadEt_%s_%s", tag.c_str(), kind[k]), ";Leading jet E_{T} [GeV];Events", 100, 0, 1000);
        h.subleadEt[k] = new TH1D(Form("subleadEt_%s_%s", tag.c_str(), kind[k]), ";Subleading jet E_{T} [GeV];Events", 100, 0, 1000);
        h.nAbove[k] = new TProfile(Form("nAbove_%s_%s", tag.c_str(), kind[k]), ";Jet E_{T} threshold [GeV];Mean number of jets / event above threshold",
                                   n_multiplicity_thresholds_, -0.5 * multiplicity_threshold_step_,
                                   (n_multiplicity_thresholds_ - 0.5) * multiplicity_threshold_step_);
        for (int r = 0; r < 2; ++r) {
            h.closure[r][k] = new TProfile(Form("closure_%s_eta%d_%s", tag.c_str(), r, kind[k]),
                                           ";Truth jet E_{T} [GeV];Mean E_{T} / E_{T}^{truth}", nClosureBins, closureEtEdges);
        }
        const char* yTitle = (k == 0) ? "Leading jet E_{T}, uncalibrated [GeV]" : "Leading jet E_{T}, calibrated [GeV]";
        h.leadVsTruth[k] = new TH2F(Form("leadVsTruth_%s_%s", tag.c_str(), kind[k]),
                                    Form(";Leading truth jet E_{T} [GeV];%s", yTitle), 100, 0, 1000, 100, 0, 1000);
        h.leadVsHLT[k] = new TH2F(Form("leadVsHLT_%s_%s", tag.c_str(), kind[k]),
                                  Form(";Leading HLT AntiKt4EMTopo jet E_{T} [GeV];%s", yTitle), 100, 0, 1000, 100, 0, 1000);
        h.leadResidual[k] = new TH1D(Form("leadResidual_%s_%s", tag.c_str(), kind[k]),
                                     ";Leading jet (E_{T} - E_{T}^{truth}) / E_{T}^{truth};Fraction of jets", 100, -1, 1);
        h.leadResidualUnweighted[k] = (TH1D*)h.leadResidual[k]->Clone(Form("leadResidualUnweighted_%s_%s", tag.c_str(), kind[k]));
        h.fourthResidual[k] = new TH1D(Form("fourthResidual_%s_%s", tag.c_str(), kind[k]),
                                       ";4th leading jet (E_{T} - E_{T}^{truth}) / E_{T}^{truth};Fraction of jets", 100, -1, 1);
        h.fourthResidualUnweighted[k] = (TH1D*)h.fourthResidual[k]->Clone(Form("fourthResidualUnweighted_%s_%s", tag.c_str(), kind[k]));
        h.leadEtFine[k] = new TH1D(Form("leadEtFine_%s_%s", tag.c_str(), kind[k]), ";Leading jet E_{T} [GeV];Events", 2000, 0, 1000);
        h.fourthEtFine[k] = new TH1D(Form("fourthEtFine_%s_%s", tag.c_str(), kind[k]), ";4th leading jet E_{T} [GeV];Events", 2000, 0, 200);
    }
    h.factorVsEt = new TProfile(Form("factorVsEt_%s", tag.c_str()), ";Uncalibrated jet E_{T} [GeV];Mean E_{T}^{calib} / E_{T}^{raw}", 50, 0, 500);
    h.factorVsEta = new TProfile(Form("factorVsEta_%s", tag.c_str()), ";Jet #eta;Mean E_{T}^{calib} / E_{T}^{raw} (E_{T}^{raw} > 30 GeV)", 98, -4.9, 4.9);
}

// Files of one sample: the hadd'ed ntuple of each sample directory (..._DAOD_NTUPLE_GEP.root), one per
// JZ slice for the background. The per-job ..._GEP_NNNNNN.root files are not read.
std::vector<std::string> calibComparisonFilePatterns(const std::string& sample, unsigned int pileup) {
    const std::string base = calib_comparison_ntuple_base_ + (pileup == 140 ? "_PU140" : "");
    const std::string rtag = (pileup == 140) ? "r16129" : "r16130";
    std::vector<std::string> patterns;
    if (sample == "QCD_Dijet") {
        for (unsigned int jz = 0; jz < 10; ++jz)
            patterns.push_back(Form("%s/QCD_Dijet_JZ%u_v%u/mc21_14TeV_jj_JZ%u_e8557_s4422_%s_DAOD_NTUPLE_GEP.root",
                                    base.c_str(), jz, calib_comparison_algo_version_, jz, rtag.c_str()));
    } else {
        patterns.push_back(Form("%s/%s_v%u/*_%s_DAOD_NTUPLE_GEP.root", base.c_str(), sample.c_str(),
                                calib_comparison_algo_version_, rtag.c_str()));
    }
    return patterns;
}

std::unique_ptr<TChain> makeCalibComparisonChain(const std::string& treeName, const std::vector<std::string>& patterns) {
    auto chain = std::make_unique<TChain>(treeName.c_str());
    for (const auto& pattern : patterns) chain->Add(pattern.c_str());
    chain->SetBranchStatus("*", 0);
    return chain;
}

// JZ slice of each file of a chain, from the "_JZ<n>_" in its name; -1 where there is none
std::vector<int> calibComparisonFileSlices(TChain* chain) {
    std::vector<int> slices;
    TObjArray* files = chain->GetListOfFiles();
    for (int f = 0; f < files->GetEntries(); ++f) {
        const std::string name = gSystem->BaseName(files->At(f)->GetTitle());
        const size_t pos = name.find("_JZ");
        slices.push_back(pos == std::string::npos ? -1 : std::atoi(name.c_str() + pos + 3));
    }
    return slices;
}

// Fill the jet E_T spectra and the multiplicity above threshold of one E_T-sorted jet collection
void fillCalibComparisonEventHists(const std::vector<double>& et, double w, TH1D* all, TH1D* lead, TH1D* sublead, TProfile* nAbove) {
    for (double e : et) all->Fill(e, w);
    if (et.size() > 0) lead->Fill(et[0], w);
    if (et.size() > 1) sublead->Fill(et[1], w);
    int counts[n_multiplicity_thresholds_] = {0};
    for (double e : et) {
        for (int k = 0; k < n_multiplicity_thresholds_ && k * multiplicity_threshold_step_ < e; ++k) counts[k]++;
    }
    for (int k = 0; k < n_multiplicity_thresholds_; ++k) nAbove->Fill(k * multiplicity_threshold_step_, counts[k], w);
}

bool fillCalibComparisonSample(CalibComparisonSample& s, const std::string& sample, unsigned int pileup) {
    const std::vector<std::string> patterns = calibComparisonFilePatterns(sample, pileup);

    auto info = makeCalibComparisonChain("eventInfoTree", patterns);
    const Long64_t nEntries = info->GetEntries();
    std::cout << s.label << ": " << info->GetNtrees() << " file(s), " << nEntries << " events\n";
    if (nEntries <= 0) {
        std::cerr << "[error] no events for " << sample << " in " << patterns.front() << "\n";
        return false;
    }
    std::vector<double>* eventWeights = nullptr;
    bool passHSTP = true;
    info->SetBranchStatus("eventWeights", 1);
    info->SetBranchStatus("passHSTP", 1);
    info->SetBranchAddress("eventWeights", &eventWeights);
    info->SetBranchAddress("passHSTP", &passHSTP);
    const std::vector<int> fileSlices = calibComparisonFileSlices(info.get());

    auto sameEntries = [&](TChain* chain) {
        const Long64_t n = chain->GetEntries();
        if (n == nEntries) return true;
        std::cerr << "[error] " << chain->GetName() << " has " << n << " entries in the " << sample << " ntuples, eventInfoTree "
                  << nEntries << ". The Calib trees need ntuples made with apply_gep_jes_calibration_ = true.\n";
        return false;
    };

    auto truth = makeCalibComparisonChain("antiKt4TruthJetsTree", patterns);
    if (!sameEntries(truth.get())) return false;
    std::vector<double> *truthPt = nullptr, *truthEta = nullptr, *truthPhi = nullptr, *truthMass = nullptr;
    truth->SetBranchStatus("Pt", 1); truth->SetBranchAddress("Pt", &truthPt);
    truth->SetBranchStatus("Eta", 1); truth->SetBranchAddress("Eta", &truthEta);
    truth->SetBranchStatus("Phi", 1); truth->SetBranchAddress("Phi", &truthPhi);
    truth->SetBranchStatus("Mass", 1); truth->SetBranchAddress("Mass", &truthMass);

    auto hlt = makeCalibComparisonChain("leadingHltAntiKt4EMTopoJetsTree", patterns);
    if (!sameEntries(hlt.get())) return false;
    std::vector<double>* hltLeadEt = nullptr;
    hlt->SetBranchStatus("Et", 1); hlt->SetBranchAddress("Et", &hltLeadEt);

    // Turn-on reference, as in largeRJetAnalysisCore.h: AntiKt4TruthDressedWZJets, E_T-sorted
    auto truthWZ = makeCalibComparisonChain("truthAntiKt4TruthDressedWZJets", patterns);
    if (!sameEntries(truthWZ.get())) return false;
    std::vector<double>* truthWZEt = nullptr;
    truthWZ->SetBranchStatus("Et", 1); truthWZ->SetBranchAddress("Et", &truthWZEt);

    std::unique_ptr<TChain> raw[n_calib_variants_], calib[n_calib_variants_];
    std::vector<double>* rawEt[n_calib_variants_] = {nullptr};
    std::vector<double>* calEt[n_calib_variants_] = {nullptr};
    std::vector<double>* calEtRaw[n_calib_variants_] = {nullptr};
    std::vector<double>* calEta[n_calib_variants_] = {nullptr};
    std::vector<double>* calPhi[n_calib_variants_] = {nullptr};
    for (int v = 0; v < n_calib_variants_; ++v) {
        raw[v] = makeCalibComparisonChain("gepWTAConeCellsTowers" + calib_variant_tags_[v] + "JetsTree", patterns);
        calib[v] = makeCalibComparisonChain("gepWTAConeCellsTowers" + calib_variant_tags_[v] + "CalibJetsTree", patterns);
        if (!sameEntries(raw[v].get()) || !sameEntries(calib[v].get())) return false;
        raw[v]->SetBranchStatus("Et", 1); raw[v]->SetBranchAddress("Et", &rawEt[v]);
        calib[v]->SetBranchStatus("Et", 1); calib[v]->SetBranchAddress("Et", &calEt[v]);
        calib[v]->SetBranchStatus("EtRaw", 1); calib[v]->SetBranchAddress("EtRaw", &calEtRaw[v]);
        calib[v]->SetBranchStatus("Eta", 1); calib[v]->SetBranchAddress("Eta", &calEta[v]);
        calib[v]->SetBranchStatus("Phi", 1); calib[v]->SetBranchAddress("Phi", &calPhi[v]);
    }

    long nInconsistent[n_calib_variants_] = {0};
    long nNoSlice = 0;
    std::vector<double> truthEt;
    for (Long64_t i = 0; i < nEntries; ++i) {
        if ((i + 1) % 200000 == 0) std::cout << "  " << i + 1 << " / " << nEntries << "\n";
        info->GetEntry(i);
        if (s.isBackground && (!eventWeights || eventWeights->empty())) continue;
        const double w = s.isBackground ? eventWeights->at(0) : 1.0;
        CalibComparisonSlice* slice = nullptr;
        if (s.isBackground) {
            const int jz = fileSlices[info->GetTreeNumber()];
            if (jz >= 0 && jz < n_jz_slices_) slice = &s.slices[jz];
            else nNoSlice++;
        }
        if (slice) { slice->nEvents++; slice->sumW += w; }
        if (s.isBackground && !passHSTP) continue;
        if (slice) { slice->nHSTP++; slice->sumWHSTP += w; }
        s.sumWeights += w;
        s.nEvents++;

        // Truth jets as in the LUT derivation, in E_T: E_T >= 10 GeV, |eta| <= 4.9; a reference jet
        // must have E_T >= 20 GeV and no other truth jet within dR 0.6
        truth->GetEntry(i);
        hlt->GetEntry(i);
        truthEt.resize(truthPt->size());
        double leadTruthEt = -1.0; // The tree is sorted by p_T, not E_T
        for (size_t t = 0; t < truthPt->size(); ++t) {
            truthEt[t] = calibComparisonTruthEt(truthPt->at(t), truthEta->at(t), truthMass->at(t));
            leadTruthEt = std::max(leadTruthEt, truthEt[t]);
        }
        const double leadHLTEt = hltLeadEt->empty() ? -1.0 : hltLeadEt->at(0);
        CalibComparisonTurnOnEvent turnOn;
        if (!s.isBackground) {
            truthWZ->GetEntry(i);
            if (truthWZEt->size() >= 1) turnOn.truthLead = truthWZEt->at(0);
            if (truthWZEt->size() >= 4) turnOn.truth4th = truthWZEt->at(3);
        }
        std::vector<size_t> truthSel;
        for (size_t t = 0; t < truthEt.size(); ++t)
            if (truthEt[t] >= 10.0 && std::fabs(truthEta->at(t)) <= 4.9) truthSel.push_back(t);
        std::vector<bool> truthRef(truthEt.size(), false);
        for (size_t a : truthSel) {
            if (truthEt[a] < 20.0) continue;
            bool isolated = true;
            for (size_t b : truthSel)
                if (b != a && calibComparisonDeltaR(truthEta->at(a), truthPhi->at(a), truthEta->at(b), truthPhi->at(b)) <= 0.6) isolated = false;
            truthRef[a] = isolated;
        }

        for (int v = 0; v < n_calib_variants_; ++v) {
            raw[v]->GetEntry(i);
            calib[v]->GetEntry(i);
            CalibComparisonHists& h = s.hists[v];

            // The Calib tree holds the same jets as the uncalibrated tree, re-sorted by calibrated E_T
            double sumRaw = 0.0, sumCalRaw = 0.0;
            for (double e : *rawEt[v]) sumRaw += e;
            for (double e : *calEtRaw[v]) sumCalRaw += e;
            if (rawEt[v]->size() != calEtRaw[v]->size() || std::fabs(sumRaw - sumCalRaw) > 1e-6 * std::fabs(sumRaw) + 1e-9) nInconsistent[v]++;

            fillCalibComparisonEventHists(*rawEt[v], w, h.allEt[0], h.leadEt[0], h.subleadEt[0], h.nAbove[0]);
            fillCalibComparisonEventHists(*calEt[v], w, h.allEt[1], h.leadEt[1], h.subleadEt[1], h.nAbove[1]);

            const std::vector<double>* leadSource[2] = {rawEt[v], calEt[v]};
            // Rate-matched turn-ons: leading / 4th leading jet E_T, 0 without such a jet (as in the core)
            for (int k = 0; k < 2; ++k) {
                const std::vector<double>& et = *leadSource[k];
                turnOn.lead[v][k] = et.empty() ? 0.0 : et[0];
                turnOn.fourth[v][k] = (et.size() >= 4) ? et[3] : 0.0;
                h.leadEtFine[k]->Fill(turnOn.lead[v][k], w);
                h.fourthEtFine[k]->Fill(turnOn.fourth[v][k], w);
            }
            for (int k = 0; k < 2; ++k) {
                if (leadSource[k]->empty()) continue;
                if (leadTruthEt > 0) h.leadVsTruth[k]->Fill(leadTruthEt, leadSource[k]->at(0), w);
                if (leadHLTEt > 0) h.leadVsHLT[k]->Fill(leadHLTEt, leadSource[k]->at(0), w);
            }

            if (v == 0 && slice) {
                for (int k = 0; k < 2; ++k)
                    for (double e : *leadSource[k])
                        for (int t = 0; t < n_slice_thresholds_; ++t)
                            if (e > slice_thresholds_[t]) slice->sumWNJets[k][t] += w;
            }

            // Leading / 4th leading jet residual, as in ResHistMaker (its jet_0 / jet_3): the jet of that rank
            // in raw E_T (the same jet before and after calibration), |eta| < 2.0, raw E_T >= 10 GeV, no other
            // jet with raw E_T >= 7 GeV within dR 0.4, matched to the closest reference truth jet (E_T >= 20 GeV,
            // isolated) within dR 0.15
            {
                std::vector<size_t> byRawEt(calEtRaw[v]->size());
                for (size_t j = 0; j < byRawEt.size(); ++j) byRawEt[j] = j;
                std::sort(byRawEt.begin(), byRawEt.end(), [&](size_t a, size_t b) { return calEtRaw[v]->at(a) > calEtRaw[v]->at(b); });
                const size_t ranks[2] = {0, 3};
                TH1D** residuals[2][2] = {{h.leadResidual, h.leadResidualUnweighted}, {h.fourthResidual, h.fourthResidualUnweighted}};
                for (int r = 0; r < 2; ++r) {
                    if (byRawEt.size() <= ranks[r]) continue;
                    const size_t jet = byRawEt[ranks[r]];
                    bool use = std::fabs(calEta[v]->at(jet)) < 2.0 && calEtRaw[v]->at(jet) >= 10.0;
                    for (size_t b = 0; use && b < calEtRaw[v]->size(); ++b)
                        if (b != jet && calEtRaw[v]->at(b) >= 7.0 &&
                            calibComparisonDeltaR(calEta[v]->at(jet), calPhi[v]->at(jet), calEta[v]->at(b), calPhi[v]->at(b)) <= 0.4) use = false;
                    if (!use) continue;
                    int match = -1;
                    double bestDR = 0.15;
                    for (size_t t : truthSel) {
                        if (!truthRef[t]) continue;
                        const double dR = calibComparisonDeltaR(calEta[v]->at(jet), calPhi[v]->at(jet), truthEta->at(t), truthPhi->at(t));
                        if (dR < bestDR) { bestDR = dR; match = (int)t; }
                    }
                    if (match < 0) continue;
                    const double etTruth = truthEt[match];
                    const double res[2] = {(calEtRaw[v]->at(jet) - etTruth) / etTruth, (calEt[v]->at(jet) - etTruth) / etTruth};
                    for (int k = 0; k < 2; ++k) {
                        residuals[r][0][k]->Fill(res[k], w);
                        residuals[r][1][k]->Fill(res[k]);
                    }
                }
            }

            // Jets entering the LUT derivation: raw E_T >= 10 GeV, |eta| <= 4.9, no other such jet within dR 0.4
            std::vector<size_t> jetSel;
            for (size_t j = 0; j < calEt[v]->size(); ++j) {
                const double etRaw = calEtRaw[v]->at(j);
                if (etRaw > 0) {
                    h.factorVsEt->Fill(etRaw, calEt[v]->at(j) / etRaw, w);
                    if (etRaw > 30.0) h.factorVsEta->Fill(calEta[v]->at(j), calEt[v]->at(j) / etRaw, w);
                    s.lowEt[v].sumWJets += w;
                    if (etRaw < lut_derivation_min_et_) {
                        s.lowEt[v].sumWLow += w;
                        if (calEt[v]->at(j) >= lut_derivation_min_et_) s.lowEt[v].sumWLowPromoted += w;
                    }
                }
                if (etRaw >= 10.0 && std::fabs(calEta[v]->at(j)) <= 4.9) jetSel.push_back(j);
            }
            for (size_t a : jetSel) {
                bool isolated = true;
                for (size_t b : jetSel)
                    if (b != a && calibComparisonDeltaR(calEta[v]->at(a), calPhi[v]->at(a), calEta[v]->at(b), calPhi[v]->at(b)) <= 0.4) isolated = false;
                if (!isolated) continue;
                int match = -1;
                double bestDR = 0.3;
                for (size_t t : truthSel) {
                    if (!truthRef[t]) continue;
                    const double dR = calibComparisonDeltaR(calEta[v]->at(a), calPhi[v]->at(a), truthEta->at(t), truthPhi->at(t));
                    if (dR <= bestDR) { bestDR = dR; match = (int)t; }
                }
                if (match < 0) continue;
                const double etTruth = truthEt[match];
                const int region = (std::fabs(truthEta->at(match)) < 2.5) ? 0 : 1;
                h.closure[region][0]->Fill(etTruth, calEtRaw[v]->at(a) / etTruth, w);
                h.closure[region][1]->Fill(etTruth, calEt[v]->at(a) / etTruth, w);
            }
        }
        if (!s.isBackground) s.turnOnEvents.push_back(turnOn);
    }

    std::cout << s.label << ": " << s.nEvents << " events used" << (s.isBackground ? " (passHSTP)" : "")
              << ", sum of weights " << s.sumWeights << "\n";
    for (int v = 0; v < n_calib_variants_; ++v) {
        std::cout << "  " << calib_variant_labels_[v] << ": " << nInconsistent[v]
                  << " event(s) whose Calib tree jets differ from the uncalibrated tree"
                  << (nInconsistent[v] ? "  <-- should be 0" : "") << "\n";
    }

    // Jets below the LUT derivation range, which take the response of the 10-20 GeV bin
    if (s.sumWeights > 0) {
        printf("  Jets / event (E_T > 0), and those with raw E_T < %.0f GeV, below the LUT derivation range:\n", lut_derivation_min_et_);
        for (int v = 0; v < n_calib_variants_; ++v) {
            const CalibComparisonLowEt& l = s.lowEt[v];
            printf("    %-17s  all %.3f   raw < %.0f GeV %.3f   of which calibrated to >= %.0f GeV %.3f\n",
                   calib_variant_labels_[v].c_str(), l.sumWJets / s.sumWeights, lut_derivation_min_et_,
                   l.sumWLow / s.sumWeights, lut_derivation_min_et_, l.sumWLowPromoted / s.sumWeights);
        }
    }

    if (s.isBackground) {
        if (nNoSlice > 0) std::cout << "  [warning] " << nNoSlice << " event(s) in files without a JZ slice in their name\n";
        double totWHSTP = 0.0;
        for (const CalibComparisonSlice& sl : s.slices) totWHSTP += sl.sumWHSTP;
        printf("\n  Background by JZ slice (%s). keep = passHSTP fraction, share = fraction of the passHSTP weight.\n",
               calib_variant_labels_[0].c_str());
        printf("  Mean jets / event above E_T threshold, uncalibrated | calibrated:\n");
        printf("  JZ    events  passHSTP  keep(N)   keep(w)  share(w) |     >0     >5    >10    >20 |     >0     >5    >10    >20\n");
        double totWN[2][n_slice_thresholds_] = {};
        for (int jz = 0; jz < n_jz_slices_; ++jz) {
            const CalibComparisonSlice& sl = s.slices[jz];
            if (sl.nEvents == 0) continue;
            printf("  %2d %9ld %9ld %8.4f %9.3g %9.4f |", jz, sl.nEvents, sl.nHSTP, (double)sl.nHSTP / sl.nEvents,
                   sl.sumW > 0 ? sl.sumWHSTP / sl.sumW : 0.0, totWHSTP > 0 ? sl.sumWHSTP / totWHSTP : 0.0);
            for (int k = 0; k < 2; ++k) {
                for (int t = 0; t < n_slice_thresholds_; ++t) {
                    printf(" %6.3f", sl.sumWHSTP > 0 ? sl.sumWNJets[k][t] / sl.sumWHSTP : 0.0);
                    totWN[k][t] += sl.sumWNJets[k][t];
                }
                printf(k == 0 ? " |" : "\n");
            }
        }
        // Equals the 0 / 5 / 10 / 20 GeV points of the multiplicity plot
        printf("  all%49s|", "");
        for (int k = 0; k < 2; ++k) {
            for (int t = 0; t < n_slice_thresholds_; ++t) printf(" %6.3f", totWHSTP > 0 ? totWN[k][t] / totWHSTP : 0.0);
            printf(k == 0 ? " |" : "\n\n");
        }
    }
    return true;
}

// Uncalibrated and calibrated curves in the current pad; raw may be null (single curve).
// refLineY > 0 draws a dashed horizontal line there; logMin > 0 fixes the lower edge of a log y axis.
void drawCalibComparisonPad(TH1* raw, TH1* cal, bool logy, bool logx, const std::string& line1, const std::string& line2,
                            const char* calLegend = "Calibrated", double refLineY = -1.0, double logMin = -1.0,
                            const std::string& line3 = "") {
    const int rawColor = TColor::GetColor("#5790fc"); // Petroff 6 blue
    const int calColor = TColor::GetColor("#e42536"); // Petroff 6 red
    gPad->SetLogy(logy);
    gPad->SetLogx(logx);
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);

    std::vector<TH1*> curves;
    if (raw) curves.push_back(raw);
    curves.push_back(cal);
    double yMax = 0.0, yMinPositive = 1e300;
    for (TH1* c : curves) {
        for (int b = 1; b <= c->GetNbinsX(); ++b) {
            const double y = c->GetBinContent(b);
            yMax = std::max(yMax, y);
            if (y > 0) yMinPositive = std::min(yMinPositive, y);
        }
    }
    const bool isProfile = cal->InheritsFrom(TProfile::Class());
    for (size_t i = 0; i < curves.size(); ++i) {
        TH1* c = curves[i];
        const int color = (raw && i == 0) ? rawColor : calColor;
        c->SetLineColor(color);
        c->SetMarkerColor(color);
        c->SetMarkerStyle(20);
        c->SetMarkerSize(0.6);
        c->SetLineWidth(2);
        // Headroom for the text and legend above the curves
        if (logy) {
            c->SetMaximum(yMax * 30.0);
            c->SetMinimum(logMin > 0 ? logMin : (yMinPositive > 1e299 ? 1e-3 : yMinPositive * 0.5));
        }
        else { c->SetMaximum(std::max(yMax, refLineY) * 1.4); c->SetMinimum(0.0); }
        c->Draw(i == 0 ? (isProfile ? "E" : "HIST") : (isProfile ? "E SAME" : "HIST SAME"));
    }

    if (refLineY > 0) {
        TLine* line = new TLine(cal->GetXaxis()->GetXmin(), refLineY, cal->GetXaxis()->GetXmax(), refLineY);
        line->SetLineStyle(2);
        line->SetLineColor(kGray + 2);
        line->Draw();
    }

    TLegend* leg = new TLegend(0.65, 0.79, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    if (raw) leg->AddEntry(raw, "Uncalibrated", isProfile ? "lp" : "l");
    leg->AddEntry(cal, calLegend, isProfile ? "lp" : "l");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.19, 0.85, line1.c_str());
    latex.DrawLatex(0.19, 0.80, line2.c_str());
    if (!line3.empty()) latex.DrawLatex(0.19, 0.75, line3.c_str());
}

// Copy of h normalized to unit area, with the given y (1D) or z (2D) axis title
TH1* normalizedCalibComparisonCopy(const TH1* h, const char* title) {
    TH1* copy = (TH1*)h->Clone(Form("%s_normalized", h->GetName()));
    const double integral = copy->Integral();
    if (integral > 0) copy->Scale(1.0 / integral);
    if (copy->GetDimension() == 2) copy->GetZaxis()->SetTitle(title);
    else copy->GetYaxis()->SetTitle(title);
    return copy;
}

// 2D correlation in the current pad, normalized to unit area, with the y = x diagonal and the
// least-squares line through the events. Its slope Cov(x, y) / Var(x) and the correlation
// coefficient r come from the weighted fill moments inside the axis ranges.
void drawCalibComparison2D(const TH2F* hFilled, const std::string& line1, const std::string& line2) {
    gPad->SetLogz(true);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.18); // Room for the palette labels and the z title
    gPad->SetBottomMargin(0.13);
    TH1* h = normalizedCalibComparisonCopy(hFilled, "Fraction of events");
    h->SetMinimum(normalized_log_min_);
    h->GetZaxis()->SetTitleOffset(1.4);
    h->Draw("COLZ");
    TLine* diagonal = new TLine(0, 0, h->GetXaxis()->GetXmax(), h->GetYaxis()->GetXmax());
    diagonal->SetLineStyle(2);
    diagonal->SetLineColor(kGray + 2);
    diagonal->Draw();

    const double varX = hFilled->GetCovariance(1, 1);
    const double slope = varX > 0 ? hFilled->GetCovariance(1, 2) / varX : 0.0;
    const double intercept = hFilled->GetMean(2) - slope * hFilled->GetMean(1);
    const double r = hFilled->GetCorrelationFactor(1, 2);
    TF1* fitLine = new TF1(Form("%s_fitLine", hFilled->GetName()), "pol1", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
    fitLine->SetParameters(intercept, slope);
    fitLine->SetLineColor(TColor::GetColor("#e42536")); // Petroff 6 red
    fitLine->SetLineWidth(2);
    fitLine->Draw("SAME");

    TLegend* leg = new TLegend(0.55, 0.16, 0.79, 0.26);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(diagonal, "y = x", "l");
    leg->AddEntry(fitLine, "Linear fit", "l");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.19, 0.85, line1.c_str());
    latex.DrawLatex(0.19, 0.80, line2.c_str());
    latex.DrawLatex(0.19, 0.75, Form("Linear fit: slope = %.3f, r = %.3f", slope, r));
}

// Leading-jet residual in the current pad, uncalibrated vs calibrated, unit area, each with a Gaussian
// fit over the full range (as ResHistMaker) whose mean and width go in the legend
void drawCalibComparisonResidual(const TH1D* raw, const TH1D* cal, const std::string& line1, const std::string& line2) {
    const int colors[2] = {TColor::GetColor("#5790fc"), TColor::GetColor("#e42536")}; // Petroff 6 blue, red
    const char* names[2] = {"Uncalibrated", "Calibrated"};
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    TH1* h[2] = {normalizedCalibComparisonCopy(raw, "Fraction of jets"), normalizedCalibComparisonCopy(cal, "Fraction of jets")};
    const double yMax = std::max(h[0]->GetMaximum(), h[1]->GetMaximum());

    TLegend* leg = new TLegend(0.52, 0.70, 0.89, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.032);
    for (int k = 0; k < 2; ++k) {
        h[k]->SetLineColor(colors[k]);
        h[k]->SetLineWidth(2);
        h[k]->SetMinimum(0.0);
        h[k]->SetMaximum(yMax * 1.5); // Headroom for the legend
        h[k]->Draw(k == 0 ? "HIST" : "HIST SAME");
        TF1* fit = new TF1(Form("%s_gaus", h[k]->GetName()), "gaus", h[k]->GetXaxis()->GetXmin(), h[k]->GetXaxis()->GetXmax());
        const bool ok = h[k]->Integral() > 0 && h[k]->Fit(fit, "Q0R") == 0;
        fit->SetLineColor(colors[k]);
        fit->SetLineStyle(2);
        fit->SetLineWidth(2);
        if (ok) fit->Draw("SAME");
        leg->AddEntry(h[k], names[k], "l");
        leg->AddEntry(fit, ok ? Form("Fit: mu = %.3f, sigma = %.3f", fit->GetParameter(1), fit->GetParameter(2)) : "Fit failed", "l");
    }
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.19, 0.85, line1.c_str());
    latex.DrawLatex(0.19, 0.80, line2.c_str());
}

// E_T threshold at which the weighted background rate above it equals targetHz. The rate above each
// bin's lower edge is the summed weight from that bin up, overflow included; the threshold is
// interpolated linearly between the two edges that bracket targetHz. NaN if the rate never gets there.
double calibComparisonThresholdForRate(const TH1D* hWeighted, double targetHz) {
    const int n = hWeighted->GetNbinsX();
    std::vector<double> above(n + 2, 0.0); // above[b]: rate above the lower edge of bin b
    above[n + 1] = hWeighted->GetBinContent(n + 1);
    for (int b = n; b >= 1; --b) above[b] = hWeighted->GetBinContent(b) + above[b + 1];
    for (int b = 2; b <= n; ++b) {
        if (above[b] > targetHz) continue;
        const double x0 = hWeighted->GetXaxis()->GetBinLowEdge(b - 1), x1 = hWeighted->GetXaxis()->GetBinLowEdge(b);
        const double y0 = above[b - 1], y1 = above[b];
        return (y0 == y1) ? x1 : x0 + (targetHz - y0) / (y1 - y0) * (x1 - x0);
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// Rate-matched turn-on in the current pad, uncalibrated vs calibrated, each with its threshold in the
// legend. spectrum, the truth distribution the efficiency is binned in, is drawn behind the curves in
// grey as in DrawTurnOnWithRatio (analysisHelperFunctions.h): normalized so its bins sum to one.
void drawCalibComparisonTurnOn(TH1* effRaw, TH1* effCal, double thrRaw, double thrCal, TH1* spectrum,
                               const char* spectrumLabel, const std::string& line1, const std::string& line2,
                               const std::string& line3) {
    const int colors[2] = {TColor::GetColor("#5790fc"), TColor::GetColor("#e42536")}; // Petroff 6 blue, red
    gPad->SetLeftMargin(0.15);
    gPad->SetBottomMargin(0.13);
    TH1* curves[2] = {effRaw, effCal};
    for (int k = 0; k < 2; ++k) {
        curves[k]->SetLineColor(colors[k]);
        curves[k]->SetMarkerColor(colors[k]);
        curves[k]->SetMarkerStyle(k == 0 ? 20 : 21);
        curves[k]->SetMarkerSize(0.8);
        curves[k]->SetLineWidth(2);
        curves[k]->SetMinimum(0.0);
        curves[k]->SetMaximum(1.4); // Headroom for the text above the plateau
    }
    effRaw->Draw("E1"); // Owns the frame
    spectrum->SetFillColorAlpha(kGray + 1, 0.35);
    spectrum->SetLineColor(kGray + 2);
    spectrum->SetLineWidth(1);
    spectrum->Draw("HIST SAME");
    for (int k = 0; k < 2; ++k) curves[k]->Draw("E1 SAME"); // Markers on top of the fill
    TLine* line = new TLine(effRaw->GetXaxis()->GetXmin(), 1.0, effRaw->GetXaxis()->GetXmax(), 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray + 2);
    line->Draw();

    TLegend* leg = new TLegend(0.50, 0.16, 0.89, 0.33);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(effRaw, Form("Uncalibrated, E_{T} > %.1f GeV", thrRaw), "lp");
    leg->AddEntry(effCal, Form("Calibrated, E_{T} > %.1f GeV", thrCal), "lp");
    leg->AddEntry(spectrum, spectrumLabel, "f");
    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.19, 0.85, line1.c_str());
    latex.DrawLatex(0.19, 0.80, line2.c_str());
    latex.DrawLatex(0.19, 0.75, line3.c_str());
}

void compareCalibratedJets(unsigned int pileup = 200, std::string signal = "ggF_HHbbbb") {
    gStyle->SetOptStat(0);
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2(true); // Weighted background fills

    std::string signalLabel;
    if (signal == "ggF_HHbbbb") signalLabel = "hh#rightarrow4b [ggF]";
    else if (signal == "ttbar_allhad") signalLabel = "t#bar{t} [had]";
    else {
        std::cerr << "[error] signal must be ggF_HHbbbb or ttbar_allhad (got " << signal << ")\n";
        return;
    }
    if (pileup != 140 && pileup != 200) {
        std::cerr << "[error] pileup must be 140 or 200 (got " << pileup << ")\n";
        return;
    }

    CalibComparisonSample samples[2];
    samples[0].label = signalLabel;
    samples[0].isBackground = false;
    samples[1].label = "QCD dijet JZ0-JZ9";
    samples[1].isBackground = true;
    const std::string sampleNames[2] = {signal, "QCD_Dijet"};
    for (int s = 0; s < 2; ++s) {
        for (int v = 0; v < n_calib_variants_; ++v)
            bookCalibComparisonHists(samples[s].hists[v], sampleNames[s] + "_" + std::to_string(v));
        if (!fillCalibComparisonSample(samples[s], sampleNames[s], pileup)) return;
    }

    // Rate-matched turn-ons: each curve's threshold from the weighted background, efficiency in the
    // signal vs the leading / 4th leading truth jet E_T. [variant][leading, 4th leading][uncalibrated, calibrated]
    const double turnOnRates[2] = {single_jet_rate_hz_, fourth_jet_rate_hz_};
    const char* turnOnJetNames[2] = {"leading jet", "4th leading jet"};
    TH1D* turnOnEff[n_calib_variants_][2][2];
    double turnOnThr[n_calib_variants_][2][2];
    TH1D* turnOnSpectrum[2] = {nullptr, nullptr}; // Signal truth distribution, the same for every variant
    const char* turnOnSpectrumLabels[2] = {"Leading truth jet E_{T}", "4th leading truth jet E_{T}"};
    std::cout << "\nRate-matched thresholds (" << samples[1].label << "), uncalibrated | calibrated:\n";
    for (int v = 0; v < n_calib_variants_; ++v) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
                const TH1D* hRate = (j == 0) ? samples[1].hists[v].leadEtFine[k] : samples[1].hists[v].fourthEtFine[k];
                turnOnThr[v][j][k] = calibComparisonThresholdForRate(hRate, turnOnRates[j]);
                const std::string tag = Form("turnOn_%d_%d_%d", v, j, k);
                TH1D num((tag + "_num").c_str(), "", 40, 0, (j == 0) ? 400 : 200);
                TH1D den((tag + "_den").c_str(), "", 40, 0, (j == 0) ? 400 : 200);
                for (const CalibComparisonTurnOnEvent& e : samples[0].turnOnEvents) {
                    const double truth = (j == 0) ? e.truthLead : e.truth4th;
                    if (truth < 0) continue;
                    den.Fill(truth);
                    const double reco = (j == 0) ? e.lead[v][k] : e.fourth[v][k];
                    if (reco >= turnOnThr[v][j][k]) num.Fill(truth);
                }
                if (!turnOnSpectrum[j]) {
                    turnOnSpectrum[j] = (TH1D*)den.Clone(Form("turnOnSpectrum_%d", j));
                    const double area = turnOnSpectrum[j]->Integral();
                    if (area > 0) turnOnSpectrum[j]->Scale(1.0 / area); // Bins sum to one, as in DrawTurnOnWithRatio
                }
                TH1D* eff = (TH1D*)num.Clone((tag + "_eff").c_str());
                eff->Divide(&num, &den, 1.0, 1.0, "B");
                eff->GetXaxis()->SetTitle((j == 0) ? "Leading truth jet E_{T} [GeV]" : "4th leading truth jet E_{T} [GeV]");
                eff->GetYaxis()->SetTitle("Efficiency");
                turnOnEff[v][j][k] = eff;
            }
            printf("  %-17s  %-15s @ %3.0f kHz: %7.1f | %7.1f GeV\n", calib_variant_labels_[v].c_str(), turnOnJetNames[j],
                   turnOnRates[j] / 1e3, turnOnThr[v][j][0], turnOnThr[v][j][1]);
        }
    }

    // One PDF per plot: <outputDir>/<variant>/<sig|back>_<plot>.pdf
    const std::string outputDir = calib_comparison_output_dir_ + "/PU" + std::to_string(pileup) + "_" + signal;
    TCanvas canvas("calibComparison", "", 700, 600);
    int nWritten = 0;
    for (int v = 0; v < n_calib_variants_; ++v) {
        const std::string variantDir = outputDir + "/" + calib_variant_dirs_[v];
        gSystem->mkdir(variantDir.c_str(), true);
        for (int s = 0; s < 2; ++s) {
            CalibComparisonHists& h = samples[s].hists[v];
            const std::string prefix = variantDir + "/" + (samples[s].isBackground ? "back_" : "sig_");
            const std::string& text = samples[s].label;
            const std::string text2 = calib_variant_labels_[v] + ", PU = " + std::to_string(pileup);
            // Clear() keeps the pad's log scales and margins, so reset them or the previous plot's carry over
            auto newPlot = [&]() {
                canvas.Clear();
                canvas.cd();
                canvas.SetLogx(false);
                canvas.SetLogy(false);
                canvas.SetLogz(false);
                canvas.SetRightMargin(gStyle->GetPadRightMargin());
            };
            auto save = [&](const std::string& name) { canvas.SaveAs((prefix + name + ".pdf").c_str()); nWritten++; };

            // E_T distributions, unit area; the y title carries the bin width
            struct { TH1D** hists; const char* yTitle; const char* name; } dists[3] = {
                {h.allEt, "Fraction of jets / %g GeV", "jetEt_allJets"},
                {h.leadEt, "Fraction of events / %g GeV", "jetEt_leading"},
                {h.subleadEt, "Fraction of events / %g GeV", "jetEt_subleading"}};
            for (const auto& d : dists) {
                newPlot();
                const std::string yTitle = Form(d.yTitle, d.hists[0]->GetXaxis()->GetBinWidth(1));
                drawCalibComparisonPad(normalizedCalibComparisonCopy(d.hists[0], yTitle.c_str()),
                                       normalizedCalibComparisonCopy(d.hists[1], yTitle.c_str()),
                                       true, false, text, text2, "Calibrated", -1.0, normalized_log_min_);
                save(d.name);
            }

            newPlot();
            drawCalibComparisonPad(h.nAbove[0], h.nAbove[1], true, false, text, text2);
            save("nJetsAboveThreshold");

            // Integral above threshold: rate in Hz for the background, efficiency for the signal
            newPlot();
            TH1* rate[2];
            for (int k = 0; k < 2; ++k) {
                rate[k] = h.leadEt[k]->GetCumulative(false, "_cumulative");
                if (!samples[s].isBackground && samples[s].sumWeights > 0) rate[k]->Scale(1.0 / samples[s].sumWeights);
                rate[k]->GetXaxis()->SetTitle("Leading jet E_{T} threshold [GeV]");
                rate[k]->GetYaxis()->SetTitle(samples[s].isBackground ? "Rate [Hz]" : "Signal Efficiency");
            }
            drawCalibComparisonPad(rate[0], rate[1], true, false, text, text2);
            save(samples[s].isBackground ? "leadingJetRate_vs_threshold" : "leadingJetEfficiency_vs_threshold");

            newPlot();
            drawCalibComparisonPad(nullptr, h.factorVsEt, false, false, text, text2, "E_{T}^{calib} / E_{T}^{raw}", 1.0);
            save("calibFactor_vs_rawEt");
            newPlot();
            drawCalibComparisonPad(nullptr, h.factorVsEta, false, false, text, text2, "E_{T}^{calib} / E_{T}^{raw}", 1.0);
            save("calibFactor_vs_eta");

            newPlot();
            drawCalibComparisonPad(h.closure[0][0], h.closure[0][1], false, true, text, text2, "Calibrated", 1.0, -1.0, "|#eta^{truth}| < 2.5");
            save("closure_vs_truthEt_central");
            newPlot();
            drawCalibComparisonPad(h.closure[1][0], h.closure[1][1], false, true, text, text2, "Calibrated", 1.0, -1.0, "2.5 < |#eta^{truth}| < 4.9");
            save("closure_vs_truthEt_forward");

            newPlot();
            drawCalibComparisonResidual(h.leadResidual[0], h.leadResidual[1], text, text2);
            save("leadingJetResidual");
            if (samples[s].isBackground) {
                newPlot();
                drawCalibComparisonResidual(h.leadResidualUnweighted[0], h.leadResidualUnweighted[1], text, text2 + ", unweighted");
                save("leadingJetResidual_unweighted");
            }
            newPlot();
            drawCalibComparisonResidual(h.fourthResidual[0], h.fourthResidual[1], text, text2);
            save("4thLeadingJetResidual");
            if (samples[s].isBackground) {
                newPlot();
                drawCalibComparisonResidual(h.fourthResidualUnweighted[0], h.fourthResidualUnweighted[1], text, text2 + ", unweighted");
                save("4thLeadingJetResidual_unweighted");
            }

            const char* kindName[2] = {"uncalibrated", "calibrated"};
            for (int k = 0; k < 2; ++k) {
                newPlot();
                drawCalibComparison2D(h.leadVsTruth[k], text, text2);
                save(std::string("leadingJetEt_vs_leadingTruthJetEt_") + kindName[k]);
                newPlot();
                drawCalibComparison2D(h.leadVsHLT[k], text, text2);
                save(std::string("leadingJetEt_vs_leadingHLTJetEt_") + kindName[k]);
            }
        }

        // Rate-matched signal turn-ons
        for (int j = 0; j < 2; ++j) {
            canvas.Clear();
            canvas.cd();
            canvas.SetLogx(false);
            canvas.SetLogy(false);
            canvas.SetLogz(false);
            canvas.SetRightMargin(gStyle->GetPadRightMargin());
            drawCalibComparisonTurnOn(turnOnEff[v][j][0], turnOnEff[v][j][1], turnOnThr[v][j][0], turnOnThr[v][j][1],
                                      turnOnSpectrum[j], turnOnSpectrumLabels[j],
                                      samples[0].label, calib_variant_labels_[v] + ", PU = " + std::to_string(pileup),
                                      Form("Fixed to %.0f kHz background rate", turnOnRates[j] / 1e3));
            const std::string name = Form("sig_turnOn_%s_%.0fkHz", (j == 0) ? "leadingJet" : "4thLeadingJet", turnOnRates[j] / 1e3);
            canvas.SaveAs((variantDir + "/" + name + ".pdf").c_str());
            nWritten++;
        }
    }
    std::cout << "Wrote " << nWritten << " PDFs to " << outputDir << "/{NoSK,SK,EtaSK}/\n";
}
