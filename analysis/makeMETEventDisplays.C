#include "analysisHelperFunctions.h"

#include <regex>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TEllipse.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TSystem.h"
#include "TROOT.h"

// -----------------------------------------------------------------------
// Parse "_jetEtN_" or "_towerEtN_" (N integer or decimal) from the emu file name.
// Returns fallback when no match.
static double parseEtThreshold(const std::string& fname, const std::string& key, double fallback) {
    std::regex re("_" + key + "([0-9]+(?:\\.[0-9]+)?)_");
    std::smatch m;
    if (std::regex_search(fname, m, re)) return std::stod(m[1]);
    return fallback;
}

// -----------------------------------------------------------------------
// Per-MET-algorithm display info: name, color, line style, and the
// METx/METy/MET/SumET values used both for the phi line on page 1 and
// the vector arrow on page 2 and the text dump on page 3.
struct METEntry {
    const char* name;
    Color_t     color;
    Style_t     lineStyle;   // for the eta-phi dashed line
    double      MET;
    double      METx;
    double      METy;        // for jFEX: only MET is valid; xy/sumET marked invalid
    double      SumET;
    bool        hasXY;       // false for jFEX (no phi available)
    bool        hasSumET;    // false for jFEX and truth (truth has no SumET)
};

// -----------------------------------------------------------------------
void callMakeMETEventDisplays(std::string herInputFile,
                              std::string metEmuFile,
                              int desiredJZSlice,
                              bool signalBool,
                              std::string signalString,
                              double backMinJetMETOverTruth  = 50.0,  // background-only: GEP Jet MET - truth NonInt MET > this [GeV] (0 = any over-reco)
                              double backMinGFexMETOverTruth = 150.0  // background-only: OR with above — gFEX JwoJ MET - truth NonInt MET > this [GeV]
                              ) {
    SetPlotStyle();

    // --- Hardcoded display thresholds (jet E_T cuts for eta-phi circles / transverse arrows) ---
    const double truthJetMinEt    = 20.0;   // truth WZ-dressed AntiKt4 jets
    const double inTimePUJetMinEt = 20.0;   // in-time pileup AntiKt4 truth jets
    const double ootPUJetMinEt    = 100.0;  // out-of-time pileup AntiKt4 truth jets

    // --- Parse jet / tower Et thresholds and the emu config tag from the emu filename ---
    double jetEtThreshold   = parseEtThreshold(metEmuFile, "jetEt",   10.0);
    double towerEtThreshold = parseEtThreshold(metEmuFile, "towerEt",  2.0);
    std::cout << "jetEtThreshold (from emu file): "   << jetEtThreshold   << " GeV\n";
    std::cout << "towerEtThreshold (from emu file): " << towerEtThreshold << " GeV\n";

    // Use the emu file basename (sans .root) as the config tag in the output dir
    std::string emuBase = metEmuFile.substr(metEmuFile.rfind('/') + 1);
    std::string emuTag  = emuBase.substr(0, emuBase.rfind(".root"));

    TString outDir = "metEventDisplays/" + signalString + "_" + emuTag + "/";
    gSystem->mkdir(outDir, true);

    // --- Files ---
    TFile* herFile = TFile::Open(herInputFile.c_str(), "READ");
    TFile* emuFile = TFile::Open(metEmuFile.c_str(),  "READ");
    if (!herFile || herFile->IsZombie()) { std::cerr << "Cannot open " << herInputFile << "\n"; return; }
    if (!emuFile || emuFile->IsZombie()) { std::cerr << "Cannot open " << metEmuFile  << "\n"; return; }

    // --- Trees from HER input ntuple ---
    TTree* eventInfoTree              = (TTree*)herFile->Get("eventInfoTree");
    TTree* gepCellsTowersEtaSKTree    = (TTree*)herFile->Get("gepCellsTowersEtaSKTree");
    TTree* gepWTAConeEtaSKJetsTree    = (TTree*)herFile->Get("gepWTAConeCellsTowersEtaSKJetsTree");
    TTree* gFexMETTree                = (TTree*)herFile->Get("gFexMETTree");
    TTree* gFexMETNoiseCutTree        = (TTree*)herFile->Get("gFexMETNoiseCutTree");
    TTree* gFexMETRmsTree             = (TTree*)herFile->Get("gFexMETRmsTree");
    TTree* jFexMETTree                = (TTree*)herFile->Get("jFexMETTree");
    TTree* metTruthTree               = (TTree*)herFile->Get("metTruthTree");
    TTree* inTimeAntiKt4TruthJetsTree = (TTree*)herFile->Get("inTimeAntiKt4TruthJetsTree");
    TTree* outOfTimeAntiKt4TruthJetsTree = (TTree*)herFile->Get("outOfTimeAntiKt4TruthJetsTree");
    TTree* truthAntiKt4WZDressedJetsTree = (TTree*)herFile->Get("truthAntiKt4TruthDressedWZJets");

    // --- Tree from emulator output ---
    TTree* metTree = (TTree*)emuFile->Get("metTree");

    // ---- Event info ----
    std::vector<double>* eventWeightsValues = nullptr;
    int  sampleJZSliceValues = -1;
    bool passHSTPValues = true;
    eventInfoTree->SetBranchAddress("eventWeights",  &eventWeightsValues);
    eventInfoTree->SetBranchAddress("sampleJZSlice", &sampleJZSliceValues);
    eventInfoTree->SetBranchAddress("passHSTP",      &passHSTPValues);

    // ---- EtaSK towers (input to GEP MET) ----
    std::vector<double>* towerEtValues  = nullptr;
    std::vector<double>* towerEtaValues = nullptr;
    std::vector<double>* towerPhiValues = nullptr;
    gepCellsTowersEtaSKTree->SetBranchAddress("Et",  &towerEtValues);
    gepCellsTowersEtaSKTree->SetBranchAddress("Eta", &towerEtaValues);
    gepCellsTowersEtaSKTree->SetBranchAddress("Phi", &towerPhiValues);

    // ---- EtaSK WTA-cone jets (input to GEP MET) ----
    std::vector<double>* jetEtValues  = nullptr;
    std::vector<double>* jetEtaValues = nullptr;
    std::vector<double>* jetPhiValues = nullptr;
    gepWTAConeEtaSKJetsTree->SetBranchAddress("Et",  &jetEtValues);
    gepWTAConeEtaSKJetsTree->SetBranchAddress("Eta", &jetEtaValues);
    gepWTAConeEtaSKJetsTree->SetBranchAddress("Phi", &jetPhiValues);

    // ---- Truth / in-time-PU / out-of-time-PU AntiKt4 jets ----
    std::vector<double>* truthJetEtValues  = nullptr;
    std::vector<double>* truthJetEtaValues = nullptr;
    std::vector<double>* truthJetPhiValues = nullptr;
    truthAntiKt4WZDressedJetsTree->SetBranchAddress("Et",  &truthJetEtValues);
    truthAntiKt4WZDressedJetsTree->SetBranchAddress("Eta", &truthJetEtaValues);
    truthAntiKt4WZDressedJetsTree->SetBranchAddress("Phi", &truthJetPhiValues);

    std::vector<double>* inTimePUJetEtValues  = nullptr;
    std::vector<double>* inTimePUJetEtaValues = nullptr;
    std::vector<double>* inTimePUJetPhiValues = nullptr;
    inTimeAntiKt4TruthJetsTree->SetBranchAddress("Et",  &inTimePUJetEtValues);
    inTimeAntiKt4TruthJetsTree->SetBranchAddress("Eta", &inTimePUJetEtaValues);
    inTimeAntiKt4TruthJetsTree->SetBranchAddress("Phi", &inTimePUJetPhiValues);

    std::vector<double>* ootPUJetEtValues  = nullptr;
    std::vector<double>* ootPUJetEtaValues = nullptr;
    std::vector<double>* ootPUJetPhiValues = nullptr;
    outOfTimeAntiKt4TruthJetsTree->SetBranchAddress("Et",  &ootPUJetEtValues);
    outOfTimeAntiKt4TruthJetsTree->SetBranchAddress("Eta", &ootPUJetEtaValues);
    outOfTimeAntiKt4TruthJetsTree->SetBranchAddress("Phi", &ootPUJetPhiValues);

    // ---- gFEX JwoJ MET (from HER) ----
    double sig_gMET = 0, sig_gMETX = 0, sig_gMETY = 0, sig_gSumET = 0;
    gFexMETTree->SetBranchAddress("gMET",   &sig_gMET);
    gFexMETTree->SetBranchAddress("gMETX",  &sig_gMETX);
    gFexMETTree->SetBranchAddress("gMETY",  &sig_gMETY);
    gFexMETTree->SetBranchAddress("gSumET", &sig_gSumET);

    // ---- gFEX NoiseCut MET ----
    double sig_gMET_NC = 0, sig_gMETX_NC = 0, sig_gMETY_NC = 0, sig_gSumET_NC = 0;
    gFexMETNoiseCutTree->SetBranchAddress("gMET",   &sig_gMET_NC);
    gFexMETNoiseCutTree->SetBranchAddress("gMETX",  &sig_gMETX_NC);
    gFexMETNoiseCutTree->SetBranchAddress("gMETY",  &sig_gMETY_NC);
    gFexMETNoiseCutTree->SetBranchAddress("gSumET", &sig_gSumET_NC);

    // ---- gFEX Rms MET ----
    double sig_gMET_Rms = 0, sig_gMETX_Rms = 0, sig_gMETY_Rms = 0, sig_gSumET_Rms = 0;
    gFexMETRmsTree->SetBranchAddress("gMET",   &sig_gMET_Rms);
    gFexMETRmsTree->SetBranchAddress("gMETX",  &sig_gMETX_Rms);
    gFexMETRmsTree->SetBranchAddress("gMETY",  &sig_gMETY_Rms);
    gFexMETRmsTree->SetBranchAddress("gSumET", &sig_gSumET_Rms);

    // ---- jFEX MET (magnitude only) ----
    double sig_jMET = 0;
    jFexMETTree->SetBranchAddress("jMET", &sig_jMET);

    // ---- Truth MET (NonInt + optional XY) ----
    double sig_metTruthNonInt = 0;
    double sig_metTruthNonIntX = 0, sig_metTruthNonIntY = 0;
    metTruthTree->SetBranchAddress("metTruthNonInt", &sig_metTruthNonInt);
    bool hasTruthNonIntXY = metTruthTree->FindBranch("metTruthNonIntX") != nullptr;
    if (hasTruthNonIntXY) {
        metTruthTree->SetBranchAddress("metTruthNonIntX", &sig_metTruthNonIntX);
        metTruthTree->SetBranchAddress("metTruthNonIntY", &sig_metTruthNonIntY);
    }

    // ---- GEP MET (from emu output) ----
    double sig_TotalMET = 0, sig_TotalMETX = 0, sig_TotalMETY = 0, sig_SumET = 0;
    double sig_TowerMet = 0, sig_JetMet   = 0;
    metTree->SetBranchAddress("TotalMET",  &sig_TotalMET);
    metTree->SetBranchAddress("TotalMETX", &sig_TotalMETX);
    metTree->SetBranchAddress("TotalMETY", &sig_TotalMETY);
    metTree->SetBranchAddress("TowerMet",  &sig_TowerMet);
    metTree->SetBranchAddress("JetMet",    &sig_JetMet);
    metTree->SetBranchAddress("SumET",     &sig_SumET);

    // GEP Jet/Tower MET XY components (optional in older outputs)
    double sig_JetMetX  = 0, sig_JetMetY  = 0;
    double sig_TowerMetX = 0, sig_TowerMetY = 0;
    bool hasJetMetXY   = metTree->FindBranch("JetMetX")   != nullptr;
    bool hasTowerMetXY = metTree->FindBranch("TowerMetX") != nullptr;
    if (hasJetMetXY) {
        metTree->SetBranchAddress("JetMetX", &sig_JetMetX);
        metTree->SetBranchAddress("JetMetY", &sig_JetMetY);
    }
    if (hasTowerMetXY) {
        metTree->SetBranchAddress("TowerMetX", &sig_TowerMetX);
        metTree->SetBranchAddress("TowerMetY", &sig_TowerMetY);
    }

    // GEP scalar sums (optional)
    double sig_SumJetET = 0, sig_SumTowerET = 0;
    bool hasSumJetET   = metTree->FindBranch("SumJetET")   != nullptr;
    bool hasSumTowerET = metTree->FindBranch("SumTowerET") != nullptr;
    if (hasSumJetET)   metTree->SetBranchAddress("SumJetET",   &sig_SumJetET);
    if (hasSumTowerET) metTree->SetBranchAddress("SumTowerET", &sig_SumTowerET);

    // --- Multi-page PDF setup ---
    TString pdf_ED = outDir + "METEventDisplays.pdf";
    TCanvas cEventDisplay("cED", "cED", 800, 700);
    cEventDisplay.Print(pdf_ED + "(");
    TString pdf_XY = outDir + "METTransverseDisplays.pdf";
    TCanvas cXY("cXY", "cXY", 700, 700);
    cXY.Print(pdf_XY + "(");

    const unsigned int maxDisplays = 200;
    unsigned int nAccepted = 0;
    const int nEvents = eventInfoTree->GetEntries();

    for (int iEvt = 0; iEvt < nEvents && nAccepted < maxDisplays; ++iEvt) {
        eventInfoTree->GetEntry(iEvt);

        // JZ filter for background, HSTP filter for background (same convention as metAnalysisAndRates.C)
        if (!signalBool && desiredJZSlice >= 0 && sampleJZSliceValues != desiredJZSlice) continue;
        if (!signalBool && !passHSTPValues) continue;

        // Read MET trees before the background MET/rate filter so we can decide whether to skip
        gFexMETTree->GetEntry(iEvt);
        gFexMETNoiseCutTree->GetEntry(iEvt);
        gFexMETRmsTree->GetEntry(iEvt);
        jFexMETTree->GetEntry(iEvt);
        metTruthTree->GetEntry(iEvt);
        // metTruthNonIntX/Y stores the raw sum of non-interacting particle momenta
        // (i.e. Σ p_x of neutrinos / LSPs), which already points in the standard MET
        // direction. The DAOD constraint mpx(NonInt) - mpx(Int) - mpx(IntOut) - mpx(IntMuon) = 0
        // matches global momentum conservation Σ p_x(all) = 0 with Int/IntOut/IntMuon stored
        // sign-flipped relative to raw px sums (NOT NonInt). No flip is needed here.
        metTree->GetEntry(iEvt);

        // --- Build the list of MET algorithms (per-event values) ---
        // jFEX has only a magnitude; truth has no SumET; everything else has full info.
        std::vector<METEntry> mets;
        mets.push_back({"Truth NonInt", kMagenta+1,  2, sig_metTruthNonInt, sig_metTruthNonIntX, sig_metTruthNonIntY, 0.0, hasTruthNonIntXY, false});
        mets.push_back({"gFEX JwoJ",    kRed,        2, sig_gMET,           sig_gMETX,           sig_gMETY,           sig_gSumET,     true,  true});
        mets.push_back({"gFEX NoiseCut",kPink+7,     2, sig_gMET_NC,        sig_gMETX_NC,        sig_gMETY_NC,        sig_gSumET_NC,  true,  true});
        mets.push_back({"gFEX Rms",     kOrange+1,   2, sig_gMET_Rms,       sig_gMETX_Rms,       sig_gMETY_Rms,       sig_gSumET_Rms, true,  true});
        mets.push_back({"jFEX",         kCyan+2,     2, sig_jMET,           0.0,                 0.0,                 0.0,            false, false});
        mets.push_back({"GEP Total",    kBlue,       1, sig_TotalMET,       sig_TotalMETX,       sig_TotalMETY,       sig_SumET,      true,  true});
        mets.push_back({"GEP Jet",      kBlack,      2, sig_JetMet,         sig_JetMetX,         sig_JetMetY,         sig_SumJetET,   hasJetMetXY,   hasSumJetET});
        mets.push_back({"GEP Tower",    kGreen+2,    2, sig_TowerMet,       sig_TowerMetX,       sig_TowerMetY,       sig_SumTowerET, hasTowerMetXY, hasSumTowerET});

        // Background-only: accept event if EITHER GEP Jet MET or gFEX JwoJ MET is over-reconstructed
        // relative to truth NonInt MET by at least its respective threshold.
        if (!signalBool) {
            bool passJet  = (sig_JetMet - sig_metTruthNonInt) >= backMinJetMETOverTruth;
            bool passGFex = (sig_gMET   - sig_metTruthNonInt) >= backMinGFexMETOverTruth;
            if (!passJet && !passGFex) continue;
        }

        ++nAccepted;
        if (iEvt % 10 == 0) std::cout << "iEvt: " << iEvt << "\n";

        // Heavy trees (tower / jet collections) only needed for events that pass the filter
        gepCellsTowersEtaSKTree->GetEntry(iEvt);
        gepWTAConeEtaSKJetsTree->GetEntry(iEvt);
        truthAntiKt4WZDressedJetsTree->GetEntry(iEvt);
        inTimeAntiKt4TruthJetsTree->GetEntry(iEvt);
        outOfTimeAntiKt4TruthJetsTree->GetEntry(iEvt);

        // ===================================================================
        // Page 1: eta-phi event display
        // ===================================================================
        cEventDisplay.cd();
        cEventDisplay.Clear();

        TH2F* hEvent = new TH2F(Form("hEvent_%d", iEvt),
                                "GEP MET inputs;#eta;#varphi",
                                100, -5, 5,
                                64, -3.2, 3.2);

        // Fill the 2D hist with EtaSK towers above the threshold parsed from the emu filename
        for (unsigned int iT = 0; iT < towerEtValues->size(); ++iT) {
            if (towerEtValues->at(iT) > towerEtThreshold) {
                hEvent->Fill(towerEtaValues->at(iT),
                             towerPhiValues->at(iT),
                             towerEtValues->at(iT));
            }
        }
        hEvent->GetZaxis()->SetTitle("E_{T} [GeV]");
        hEvent->Draw("COLZ");

        // WTA-cone jets above jetEtThreshold (input to GEP MET) — dashed black circles
        std::vector<std::pair<double,double>> wtaConeJetsAboveThr;
        for (unsigned int iJ = 0; iJ < jetEtValues->size(); ++iJ) {
            if (jetEtValues->at(iJ) > jetEtThreshold) {
                wtaConeJetsAboveThr.emplace_back(jetEtaValues->at(iJ), jetPhiValues->at(iJ));
                TEllipse* c = new TEllipse(jetEtaValues->at(iJ), jetPhiValues->at(iJ), 0.4, 0.4);
                c->SetLineColor(kBlack); c->SetLineWidth(2);
                c->SetFillStyle(0); c->SetLineStyle(2);
                c->Draw("same");
            }
        }

        // Truth WZ-dressed AntiKt4 jets above 15 GeV — cyan solid
        for (unsigned int iTJ = 0; iTJ < truthJetEtValues->size(); ++iTJ) {
            if (truthJetEtValues->at(iTJ) > truthJetMinEt) {
                TEllipse* c = new TEllipse(truthJetEtaValues->at(iTJ), truthJetPhiValues->at(iTJ), 0.4, 0.4);
                c->SetLineColor(kCyan); c->SetLineWidth(2);
                c->SetFillStyle(0); c->SetLineStyle(1);
                c->Draw("same");
            }
        }

        // In-time pileup truth jets above 15 GeV — green dotted
        for (unsigned int iPU = 0; iPU < inTimePUJetEtValues->size(); ++iPU) {
            if (inTimePUJetEtValues->at(iPU) > inTimePUJetMinEt) {
                TEllipse* c = new TEllipse(inTimePUJetEtaValues->at(iPU), inTimePUJetPhiValues->at(iPU), 0.4, 0.4);
                c->SetLineColor(kGreen); c->SetLineWidth(2);
                c->SetFillStyle(0); c->SetLineStyle(3);
                c->Draw("same");
            }
        }

        // Out-of-time pileup truth jets — violet dotted (matches makeJetTaggerEventDisplays convention)
        for (unsigned int iOOT = 0; iOOT < ootPUJetEtValues->size(); ++iOOT) {
            if (ootPUJetEtValues->at(iOOT) > ootPUJetMinEt) {
                TEllipse* c = new TEllipse(ootPUJetEtaValues->at(iOOT), ootPUJetPhiValues->at(iOOT), 0.4, 0.4);
                c->SetLineColor(kViolet); c->SetLineWidth(2);
                c->SetFillStyle(0); c->SetLineStyle(3);
                c->Draw("same");
            }
        }

        // Horizontal dashed line at phi_MET for each algorithm that has x/y info
        TLegend legMET(0.12, 0.74, 0.42, 0.89);
        legMET.SetBorderSize(0); legMET.SetFillStyle(0); legMET.SetTextSize(0.022);
        for (const auto& m : mets) {
            if (!m.hasXY) continue;
            if (m.METx == 0.0 && m.METy == 0.0) continue;
            double phi = std::atan2(m.METy, m.METx);
            TLine* l = new TLine(-5.0, phi, 5.0, phi);
            l->SetLineColor(m.color); l->SetLineWidth(2); l->SetLineStyle(m.lineStyle);
            l->Draw("same");
            legMET.AddEntry(l, Form("%s (#varphi=%.2f)", m.name, phi), "l");
        }
        legMET.Draw();

        // Top-of-canvas info line
        TLatex topLat;
        topLat.SetTextSize(0.028);
        topLat.DrawLatexNDC(0.10, 0.95,
            Form("Event %d   (towers > %.1f GeV, WTA-cone jets > %.1f GeV)",
                 iEvt, towerEtThreshold, jetEtThreshold));

        cEventDisplay.Print(pdf_ED);

        // ===================================================================
        // Page 2a: x-y transverse view — jets as unit vectors at (cos #varphi, sin #varphi)
        // with their E_T labeled at the edge of the unit circle. MET algorithms drawn
        // as arrows scaled to MET / MET_max. The numeric legend is on the next page.
        // ===================================================================
        cXY.cd();
        cXY.Clear();

        // Symmetric frame around the unit circle
        const double xyHalf = 1.3;
        TH2F* hFrame = new TH2F(Form("hXYFrame_%d", iEvt),
                                Form("Transverse view (Event %d);E_{x}/E_{T};E_{y}/E_{T}", iEvt),
                                10, -xyHalf, xyHalf, 10, -xyHalf, xyHalf);
        hFrame->SetStats(0);
        hFrame->Draw();

        // Unit circle for reference
        TEllipse refCircle(0.0, 0.0, 1.0, 1.0);
        refCircle.SetFillStyle(0); refCircle.SetLineColor(kGray+1); refCircle.SetLineStyle(3);
        refCircle.Draw("same");

        // Reference axes through origin
        TLine xAxis(-xyHalf, 0, xyHalf, 0); xAxis.SetLineColor(kGray+1); xAxis.SetLineStyle(2); xAxis.Draw("same");
        TLine yAxis(0, -xyHalf, 0, xyHalf); yAxis.SetLineColor(kGray+1); yAxis.SetLineStyle(2); yAxis.Draw("same");

        // Single shared scale: longest jet OR longest MET vector hits the unit circle.
        // Jets and MET arrows are both scaled by maxScale so relative magnitudes are visible.
        double maxScale = 1e-9;
        for (unsigned int iJ = 0; iJ < jetEtValues->size(); ++iJ)
            if (jetEtValues->at(iJ) > jetEtThreshold && jetEtValues->at(iJ) > maxScale)
                maxScale = jetEtValues->at(iJ);
        for (unsigned int iTJ = 0; iTJ < truthJetEtValues->size(); ++iTJ)
            if (truthJetEtValues->at(iTJ) > truthJetMinEt && truthJetEtValues->at(iTJ) > maxScale)
                maxScale = truthJetEtValues->at(iTJ);
        for (unsigned int iPU = 0; iPU < inTimePUJetEtValues->size(); ++iPU)
            if (inTimePUJetEtValues->at(iPU) > inTimePUJetMinEt && inTimePUJetEtValues->at(iPU) > maxScale)
                maxScale = inTimePUJetEtValues->at(iPU);
        for (unsigned int iOOT = 0; iOOT < ootPUJetEtValues->size(); ++iOOT)
            if (ootPUJetEtValues->at(iOOT) > ootPUJetMinEt && ootPUJetEtValues->at(iOOT) > maxScale)
                maxScale = ootPUJetEtValues->at(iOOT);
        for (const auto& m : mets)
            if (m.hasXY && m.MET > maxScale) maxScale = m.MET;

        // Jet arrows scaled by E_T/maxScale, with a small E_T label just beyond the tip (same color).
        auto drawJetArrow = [&](double phi, double Et, Color_t col, Style_t ls) {
            double r = Et / maxScale;
            TArrow* a = new TArrow(0.0, 0.0, r*std::cos(phi), r*std::sin(phi), 0.015, "|>");
            a->SetLineColor(col); a->SetLineWidth(2); a->SetLineStyle(ls); a->SetFillColor(col);
            a->Draw();
            double labelR = r + 0.10;
            TLatex* t = new TLatex(labelR*std::cos(phi), labelR*std::sin(phi),
                                   Form("%.0f GeV", std::round(Et)));
            t->SetTextSize(0.020); t->SetTextColor(col); t->SetTextAlign(22);
            t->Draw();
            return a;
        };
        for (unsigned int iJ = 0; iJ < jetEtValues->size(); ++iJ) {
            if (jetEtValues->at(iJ) > jetEtThreshold)
                drawJetArrow(jetPhiValues->at(iJ), jetEtValues->at(iJ), kBlack, 2);
        }
        for (unsigned int iTJ = 0; iTJ < truthJetEtValues->size(); ++iTJ) {
            if (truthJetEtValues->at(iTJ) > truthJetMinEt)
                drawJetArrow(truthJetPhiValues->at(iTJ), truthJetEtValues->at(iTJ), kCyan, 2);
        }
        for (unsigned int iPU = 0; iPU < inTimePUJetEtValues->size(); ++iPU) {
            if (inTimePUJetEtValues->at(iPU) > inTimePUJetMinEt)
                drawJetArrow(inTimePUJetPhiValues->at(iPU), inTimePUJetEtValues->at(iPU), kGreen, 3);
        }
        for (unsigned int iOOT = 0; iOOT < ootPUJetEtValues->size(); ++iOOT) {
            if (ootPUJetEtValues->at(iOOT) > ootPUJetMinEt)
                drawJetArrow(ootPUJetPhiValues->at(iOOT), ootPUJetEtValues->at(iOOT), kViolet, 3);
        }

        // MET arrows — color-coded per algorithm, scaled by the same maxScale as jets.
        for (const auto& m : mets) {
            if (!m.hasXY) continue;
            if (m.METx == 0.0 && m.METy == 0.0) continue;
            double phi = std::atan2(m.METy, m.METx);
            double r   = m.MET / maxScale;
            TArrow* a = new TArrow(0.0, 0.0, r * std::cos(phi), r * std::sin(phi), 0.025, "|>");
            a->SetLineColor(m.color); a->SetLineWidth(3); a->SetFillColor(m.color);
            a->Draw();
        }

        cXY.Print(pdf_XY);

        // ===================================================================
        // Page 2b: TLatex legend for the transverse view — color-coded per algorithm,
        // magnitudes and phi (in degrees). jFEX magnitude only.
        // ===================================================================
        cXY.Clear();

        TLatex xyLeg;

        // Title: centered at top, slightly larger
        xyLeg.SetTextAlign(23);   // hcenter, top
        xyLeg.SetTextSize(0.040);
        xyLeg.SetTextColor(kBlack);
        xyLeg.DrawLatexNDC(0.50, 0.96,
            Form("Event %d: jets and MET vectors (scaled / %.0f GeV)", iEvt, maxScale));

        // Body: left-aligned, smaller, color-coded per algorithm
        xyLeg.SetTextAlign(13);
        double y_1 = 0.85;
        xyLeg.SetTextSize(0.030);
        for (const auto& m : mets) {
            xyLeg.SetTextColor(m.color);
            if (!m.hasXY) {
                xyLeg.DrawLatexNDC(0.10, y_1,
                    Form("%s: |MET| = %.1f GeV   (no #varphi)", m.name, m.MET));
            } else if (m.METx == 0.0 && m.METy == 0.0) {
                // e.g. GEP Jet MET = 0 when no WTA-cone jets pass the E_T threshold
                xyLeg.DrawLatexNDC(0.10, y_1,
                    Form("%s: |MET| = %.1f GeV", m.name, m.MET));
            } else {
                double phiDeg = std::atan2(m.METy, m.METx) * 180.0 / M_PI;
                xyLeg.DrawLatexNDC(0.10, y_1,
                    Form("%s: |MET| = %.1f GeV,   #varphi = %.1f#circ", m.name, m.MET, phiDeg));
            }
            y_1 -= 0.07;
        }
        xyLeg.SetTextColor(kBlack);

        // Background-only extras: event weight (rate contribution), JZ slice,
        // and leading/subleading truth WZ-dressed AntiKt4 jet pT.
        if (!signalBool) {
            double leadTruthEt = 0.0, subleadTruthEt = 0.0;
            if (truthJetEtValues) {
                for (double et : *truthJetEtValues) {
                    if      (et > leadTruthEt)    { subleadTruthEt = leadTruthEt; leadTruthEt = et; }
                    else if (et > subleadTruthEt) { subleadTruthEt = et; }
                }
            }
            double leadInTimePUEt = 0.0;
            if (inTimePUJetEtValues)
                for (double et : *inTimePUJetEtValues)
                    if (et > leadInTimePUEt) leadInTimePUEt = et;
            double leadOOTPUEt = 0.0;
            if (ootPUJetEtValues)
                for (double et : *ootPUJetEtValues)
                    if (et > leadOOTPUEt) leadOOTPUEt = et;
            double w = (eventWeightsValues && !eventWeightsValues->empty())
                       ? eventWeightsValues->at(0) : 0.0;
            y_1 -= 0.03;
            xyLeg.SetTextColor(kGray+2);
            xyLeg.SetTextSize(0.026);
            xyLeg.DrawLatexNDC(0.10, y_1, Form("Event weight (rate contrib.): %.4g", w));
            y_1 -= 0.05;
            xyLeg.DrawLatexNDC(0.10, y_1, Form("JZ slice: %d", sampleJZSliceValues));
            y_1 -= 0.05;
            xyLeg.DrawLatexNDC(0.10, y_1,
                Form("Leading E_{T}: truth = %.1f GeV, IT PU = %.1f GeV, OOT PU = %.1f GeV",
                     leadTruthEt, leadInTimePUEt, leadOOTPUEt));
            y_1 -= 0.05;
            xyLeg.DrawLatexNDC(0.10, y_1,
                Form("Subleading truth jet E_{T}: %.1f GeV", subleadTruthEt));
            xyLeg.SetTextColor(kBlack);
        }

        cXY.Print(pdf_XY);

        // ===================================================================
        // Page 3 (on cEventDisplay): text dump of MET values + rate contribution.
        // ===================================================================
        cEventDisplay.cd();
        cEventDisplay.Clear();

        TLatex lat;
        lat.SetTextAlign(13);

        lat.SetTextSize(0.030);
        lat.DrawLatexNDC(0.05, 0.95, Form("Event %d   (signal=%s)", iEvt, signalBool ? "true" : "false"));

        // Header row
        lat.SetTextSize(0.026);
        double y_2 = 0.86;
        lat.DrawLatexNDC(0.05, y_2, "MET algorithm          MET [GeV]    METx [GeV]    METy [GeV]    SumET [GeV]");
        y_2 -= 0.025;
        lat.DrawLatexNDC(0.05, y_2, "----------------------------------------------------------------------------------");
        y_2 -= 0.04;

        for (const auto& m : mets) {
            std::string sMETx  = m.hasXY    ? Form("%8.1f", m.METx)  : std::string("    ---");
            std::string sMETy  = m.hasXY    ? Form("%8.1f", m.METy)  : std::string("    ---");
            std::string sSumET = m.hasSumET ? Form("%8.1f", m.SumET) : std::string("    ---");
            lat.SetTextColor(m.color);
            lat.DrawLatexNDC(0.05, y_2,
                Form("%-18s   %8.1f       %s       %s       %s",
                     m.name, m.MET, sMETx.c_str(), sMETy.c_str(), sSumET.c_str()));
            y_2 -= 0.045;
        }
        lat.SetTextColor(kBlack);

        // Background extras: rate contribution + JZ slice + event number
        if (!signalBool) {
            y_2 -= 0.03;
            lat.SetTextSize(0.028);
            lat.DrawLatexNDC(0.05, y_2,
                Form("Event Weight (Rate Contribution): %.4g",
                     eventWeightsValues && !eventWeightsValues->empty() ? eventWeightsValues->at(0) : 0.0));
            y_2 -= 0.05;
            lat.DrawLatexNDC(0.05, y_2,
                Form("JZ Slice: %d, Event Number: %d", sampleJZSliceValues, iEvt));
        }

        // Emu config tag at bottom
        lat.SetTextSize(0.020);
        lat.SetTextColor(kGray+2);
        lat.DrawLatexNDC(0.05, 0.05, Form("emu config: %s", emuTag.c_str()));
        lat.SetTextColor(kBlack);

        cEventDisplay.Print(pdf_ED);

        delete hEvent;
        delete hFrame;
    }

    // close multi-page PDFs
    cEventDisplay.Print(pdf_ED + ")");
    cXY.Print(pdf_XY + ")");
}

// -----------------------------------------------------------------------
void makeMETEventDisplays() {
    gErrorIgnoreLevel = kWarning;

    // Same input paths used by metAnalysisAndRates.C
    const std::string sigHER  = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/ZvvHbb_v3/mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_resim_DAOD_NTUPLE_GEP.root";
    const std::string backHER = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/mc21_14TeV_jj_JZ_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    const std::string emuDir  = "/data/larsonma/GEPMET/outputNTuplesDev_METv2/";

    // Default signal call: ZvvHbb with the (jetEt=20, towerEt=2, EtaSK_OR, twrSF=1, jetSF=1) config
    callMakeMETEventDisplays(
        sigHER,
        emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF1_jetSF1.root",
        -1, true, "ZvvHbb");

    // Background: JZ1 only. Default filter accepts events where EITHER GEP Jet MET or gFEX JwoJ MET
    // is over-reconstructed vs truth NonInt MET by its respective threshold (50 GeV / 150 GeV).
    // To override, append e.g. /*backMinJetMETOverTruth=*/30.0, /*backMinGFexMETOverTruth=*/100.0
    callMakeMETEventDisplays(
        backHER,
        emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt20_towerEt2_EtaSK_OR_twrSF1_jetSF1.root",
        1, false, "jj_1");

    gSystem->Exit(0);
}
