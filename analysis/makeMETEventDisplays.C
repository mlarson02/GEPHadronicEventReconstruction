#include "analysisHelperFunctions.h"
#include "chainSource.h"

#include <regex>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <cmath>
#include <cstdio>
#include <iterator>
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
// Tiling of one FEX tower collection (gFEX / jFEX emulated towers) in eta.
// The granularity is NOT constant: gTowers are 0.2 x 0.2 centrally but get
// coarser in the endcap / FCal, and the same is true for jTowers. Rather than
// hardcoding the forward widths, reconstruct the tiling from the tower centres:
//
//   * the smallest gap between the distinct phi centres at a given eta is d(phi)
//   * eta rings tile contiguously and have an edge at eta = 0, so walking
//     outward from 0 with  highEdge_i = 2*centre_i - lowEdge_i  (lowEdge_i =
//     highEdge_{i-1}) gives the exact eta width of every ring.
struct FexTowerRing {
    double etaCentre = 0.0;
    double etaLow    = 0.0;
    double etaHigh   = 0.0;
    double etaWidth  = 0.0;
    double phiWidth  = 0.0;
    int    nPhi      = 0;
    int    nTowers   = 0;
};

// Round to 1e-3 so float-precision jitter in the stored centres does not split rings.
static double fexTowerKey(double v) { return std::round(v * 1000.0) / 1000.0; }

// Build (and optionally dump) the eta-ring geometry of a tower collection.
// Keyed by the rounded eta centre so per-tower lookup is a single map find.
static std::map<double, FexTowerRing> buildFexTowerGeometry(const std::vector<double>& towerEta,
                                                            const std::vector<double>& towerPhi,
                                                            const std::string& label,
                                                            bool debugPrint) {
    std::map<double, std::set<double>> phisAtEta;
    std::map<double, int>              nTowersAtEta;
    for (unsigned int i = 0; i < towerEta.size(); ++i) {
        double kEta = fexTowerKey(towerEta.at(i));
        phisAtEta[kEta].insert(fexTowerKey(towerPhi.at(i)));
        ++nTowersAtEta[kEta];
    }

    std::vector<double> etaKeys;
    etaKeys.reserve(phisAtEta.size());
    for (const auto& kv : phisAtEta) etaKeys.push_back(kv.first);

    std::map<double, FexTowerRing> geom;
    if (etaKeys.empty()) return geom;

    // Index of the first ring with centre >= 0; walk the two sides outward from eta = 0.
    unsigned int iFirstPos = 0;
    while (iFirstPos < etaKeys.size() && etaKeys.at(iFirstPos) < 0.0) ++iFirstPos;

    // Fall back on the centre-to-centre spacing if the contiguous-tiling walk gives
    // something unphysical (e.g. a collection that does not tile from eta = 0).
    auto neighbourGap = [&](unsigned int i) {
        double gap = 0.0;
        if (i + 1 < etaKeys.size()) gap = std::abs(etaKeys.at(i+1) - etaKeys.at(i));
        if (i > 0) {
            double prevGap = std::abs(etaKeys.at(i) - etaKeys.at(i-1));
            if (gap <= 0.0 || prevGap < gap) gap = prevGap;
        }
        return gap > 0.0 ? gap : 0.1;
    };

    auto addRing = [&](unsigned int i, double width, double low, double high) {
        double c = etaKeys.at(i);
        FexTowerRing r;
        r.etaCentre = c;
        r.etaLow    = low;
        r.etaHigh   = high;
        r.etaWidth  = width;
        const std::set<double>& phis = phisAtEta.at(c);
        r.nPhi    = (int)phis.size();
        r.nTowers = nTowersAtEta.at(c);
        // Smallest gap between neighbouring phi centres = d(phi) of the ring; this stays
        // correct even if a tower happens to be missing from the collection (2 pi / N would not).
        double minPhiGap = 0.0;
        for (auto itPhi = phis.begin(); itPhi != phis.end(); ++itPhi) {
            auto itNext = std::next(itPhi);
            if (itNext == phis.end()) break;
            double gap = *itNext - *itPhi;
            if (gap > 1e-6 && (minPhiGap <= 0.0 || gap < minPhiGap)) minPhiGap = gap;
        }
        if (minPhiGap <= 0.0) minPhiGap = r.nPhi > 0 ? 2.0 * M_PI / r.nPhi : 0.1;
        r.phiWidth = minPhiGap;
        geom[c] = r;
    };

    double edge = 0.0;
    for (unsigned int i = iFirstPos; i < etaKeys.size(); ++i) {
        double w = 2.0 * (etaKeys.at(i) - edge);
        if (w <= 1e-6 || w > 2.0) w = neighbourGap(i);
        addRing(i, w, etaKeys.at(i) - 0.5*w, etaKeys.at(i) + 0.5*w);
        edge = etaKeys.at(i) + 0.5*w;
    }
    edge = 0.0;
    for (int i = (int)iFirstPos - 1; i >= 0; --i) {
        double w = 2.0 * (edge - etaKeys.at(i));
        if (w <= 1e-6 || w > 2.0) w = neighbourGap((unsigned int)i);
        addRing((unsigned int)i, w, etaKeys.at(i) - 0.5*w, etaKeys.at(i) + 0.5*w);
        edge = etaKeys.at(i) - 0.5*w;
    }

    if (debugPrint) {
        std::cout << "\n=== " << label << " tower geometry (derived from tower centres) ===\n";
        std::cout << "  " << towerEta.size() << " towers, " << geom.size() << " distinct eta rings\n";
        std::cout << "   eta centre   eta low   eta high     d(eta)   N(phi)     d(phi)   N(towers)\n";
        for (const auto& kv : geom) {
            const FexTowerRing& r = kv.second;
            printf("  %10.3f %9.3f %10.3f %10.3f %8d %10.3f %11d\n",
                   r.etaCentre, r.etaLow, r.etaHigh, r.etaWidth, r.nPhi, r.phiWidth, r.nTowers);
        }
        // Group the rings by (d(eta), d(phi)) so the granularity regions are easy to read off.
        std::map<std::pair<double,double>, std::pair<double,double>> regions;  // (dEta,dPhi) -> (min|eta|, max|eta|)
        for (const auto& kv : geom) {
            const FexTowerRing& r = kv.second;
            auto key = std::make_pair(fexTowerKey(r.etaWidth), fexTowerKey(r.phiWidth));
            double lo = std::min(std::abs(r.etaLow), std::abs(r.etaHigh));
            double hi = std::max(std::abs(r.etaLow), std::abs(r.etaHigh));
            auto it = regions.find(key);
            if (it == regions.end()) regions[key] = std::make_pair(lo, hi);
            else {
                it->second.first  = std::min(it->second.first,  lo);
                it->second.second = std::max(it->second.second, hi);
            }
        }
        std::cout << "  granularity regions (|eta| coverage):\n";
        for (const auto& kv : regions)
            printf("    d(eta) x d(phi) = %.3f x %.3f  for  %.3f < |eta| < %.3f\n",
                   kv.first.first, kv.first.second, kv.second.first, kv.second.second);
        std::cout << std::endl;
    }

    return geom;
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
                              double backMinGFexMETOverTruth = 150.0, // background-only: OR with above — gFEX JwoJ MET - truth NonInt MET > this [GeV]
                              bool   drawGFexTowers          = false, // separate PDF of gFEX emulated tower eta-phi displays (every event)
                              bool   drawJFexTowers          = false, // separate PDF of jFEX emulated tower eta-phi displays (every event)
                              double fexTowerZMin            = -10.0 // z-axis floor [GeV] for the g/jFEX tower displays; every tower is drawn, no E_T threshold
                              ) {
    SetPlotStyle();

    // --- Hardcoded display thresholds (jet E_T cuts for eta-phi circles / transverse arrows) ---
    const double truthJetMinEt    = 15.0;   // truth WZ-dressed AntiKt4 jets
    const double inTimePUJetMinEt = 15.0;   // in-time pileup AntiKt4 truth jets
    const double ootPUJetMinEt    = 50.0;  // out-of-time pileup AntiKt4 truth jets

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
    // HER ntuple: ChainSource so a JZ-slice glob is read as one chained set of trees,
    // matching how the analysis macros read the same input.
    ChainSource* herFile = ChainSource::Open(herInputFile.c_str());
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

    // ---- Bunch-train position (optional: only in ntuples made after it was added to HERNTupler) ----
    int distFrontBunchTrainValues = -1, distTailBunchTrainValues = -1;
    bool hasBunchTrain = eventInfoTree->FindBranch("distFrontBunchTrain") != nullptr
                      && eventInfoTree->FindBranch("distTailBunchTrain")  != nullptr;
    if (hasBunchTrain) {
        eventInfoTree->SetBranchAddress("distFrontBunchTrain", &distFrontBunchTrainValues);
        eventInfoTree->SetBranchAddress("distTailBunchTrain",  &distTailBunchTrainValues);
    } else {
        std::cout << "distFrontBunchTrain / distTailBunchTrain not present in " << herInputFile
                  << " — bunch-train position will be printed as n/a\n";
    }

    // ---- g/jFEX emulated towers (optional standalone eta-phi displays) ----
    TTree* gFexTowersTree = drawGFexTowers ? (TTree*)herFile->Get("gFexEmulatedTowersTree") : nullptr;
    TTree* jFexTowersTree = drawJFexTowers ? (TTree*)herFile->Get("jFexEmulatedTowersTree") : nullptr;
    if (drawGFexTowers && !gFexTowersTree) {
        std::cout << "gFexEmulatedTowersTree not found — skipping gFEX tower displays\n";
        drawGFexTowers = false;
    }
    if (drawJFexTowers && !jFexTowersTree) {
        std::cout << "jFexEmulatedTowersTree not found — skipping jFEX tower displays\n";
        drawJFexTowers = false;
    }

    std::vector<double>* gFexTowerEtValues  = nullptr;
    std::vector<double>* gFexTowerEtaValues = nullptr;
    std::vector<double>* gFexTowerPhiValues = nullptr;
    if (drawGFexTowers) {
        gFexTowersTree->SetBranchAddress("Et",  &gFexTowerEtValues);
        gFexTowersTree->SetBranchAddress("Eta", &gFexTowerEtaValues);
        gFexTowersTree->SetBranchAddress("Phi", &gFexTowerPhiValues);
    }

    std::vector<double>* jFexTowerEtValues  = nullptr;
    std::vector<double>* jFexTowerEtaValues = nullptr;
    std::vector<double>* jFexTowerPhiValues = nullptr;
    if (drawJFexTowers) {
        jFexTowersTree->SetBranchAddress("Et",  &jFexTowerEtValues);
        jFexTowersTree->SetBranchAddress("Eta", &jFexTowerEtaValues);
        jFexTowersTree->SetBranchAddress("Phi", &jFexTowerPhiValues);
    }

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

    // g/jFEX tower displays live in their own PDFs and are made for every event that
    // passes the JZ / HSTP selection (i.e. they ignore the background MET filter).
    TString pdf_gTow = outDir + "gFexTowerEventDisplays.pdf";
    TString pdf_jTow = outDir + "jFexTowerEventDisplays.pdf";
    TCanvas cFexTower("cFexTower", "cFexTower", 800, 700);
    if (drawGFexTowers) cFexTower.Print(pdf_gTow + "(");
    if (drawJFexTowers) cFexTower.Print(pdf_jTow + "(");

    // Set to true to dump the derived tower geometry (eta ring edges, d(eta), d(phi),
    // granularity regions) for the first event with towers in each collection.
    const bool debugFexTowerGeometry = true;

    // Tower geometry is static, so derive it once per collection and cache it.
    std::map<double, FexTowerRing> gFexTowerGeom, jFexTowerGeom;

    const unsigned int maxDisplays         = 200;
    const unsigned int maxFexTowerDisplays = 200;
    unsigned int nAccepted = 0, nFexTowerDisplays = 0;
    const int nEvents = eventInfoTree->GetEntries();

    for (int iEvt = 0;
         iEvt < nEvents && (nAccepted < maxDisplays
                            || ((drawGFexTowers || drawJFexTowers) && nFexTowerDisplays < maxFexTowerDisplays));
         ++iEvt) {
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

        // --- Overlays shared by every eta-phi display (GEP towers, gFEX towers, jFEX towers):
        //     jet circles plus a horizontal dashed line at phi_MET for each algorithm. ---
        auto drawEtaPhiOverlays = [&]() {
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
            TLegend* legMET = new TLegend(0.12, 0.74, 0.42, 0.89);
            legMET->SetBorderSize(0); legMET->SetFillStyle(0); legMET->SetTextSize(0.022);
            for (const auto& m : mets) {
                if (!m.hasXY) continue;
                if (m.METx == 0.0 && m.METy == 0.0) continue;
                double phi = std::atan2(m.METy, m.METx);
                TLine* l = new TLine(-5.0, phi, 5.0, phi);
                l->SetLineColor(m.color); l->SetLineWidth(2); l->SetLineStyle(m.lineStyle);
                l->Draw("same");
                legMET->AddEntry(l, Form("%s (#varphi=%.2f)", m.name, phi), "l");
            }
            legMET->Draw();
        };

        // --- One eta-phi display of a g/jFEX tower collection, printed to its own PDF ---
        auto drawFexTowerDisplay = [&](const std::string& label,
                                       const std::vector<double>* tEt,
                                       const std::vector<double>* tEta,
                                       const std::vector<double>* tPhi,
                                       std::map<double, FexTowerRing>& geom,
                                       const char* tag,
                                       const TString& pdfName) {
            if (!tEt || tEt->empty()) return;
            if (geom.empty()) geom = buildFexTowerGeometry(*tEta, *tPhi, label, debugFexTowerGeometry);

            cFexTower.cd();
            cFexTower.Clear();

            // Fine uniform grid; each tower is painted over its own footprint (from the
            // derived geometry) so the coarser endcap / FCal towers show at their true size.
            const int    nEtaFine = 200, nPhiFine = 128;
            const double etaRange = 5.0, phiRange = 3.2;
            TH2F* hTow = new TH2F(Form("h%sTowers_%d", tag, iEvt),
                                  Form("%s emulated towers;#eta;#varphi", label.c_str()),
                                  nEtaFine, -etaRange, etaRange, nPhiFine, -phiRange, phiRange);
            hTow->SetStats(0);
            hTow->GetZaxis()->SetTitle("E_{T} [GeV]");

            // Every tower is drawn, with no E_T threshold: gFEX towers legitimately carry
            // negative E_T (pedestal / noise subtraction can push a tower below zero) and that
            // structure is the point of the display. The overlap rule works on |E_T| — "keep
            // the largest E_T" would never let a negative tower win a bin against the
            // histogram's default content of 0.
            double maxEt = 0.0, minEt = 0.0;
            for (unsigned int iT = 0; iT < tEt->size(); ++iT) {
                double et = tEt->at(iT);
                double eta = tEta->at(iT), phi = tPhi->at(iT);
                auto itGeom = geom.find(fexTowerKey(eta));
                double dEta = (itGeom != geom.end()) ? itGeom->second.etaWidth : 0.1;
                double dPhi = (itGeom != geom.end()) ? itGeom->second.phiWidth : 0.1;
                int bxLo = hTow->GetXaxis()->FindBin(eta - 0.5*dEta + 1e-6);
                int bxHi = hTow->GetXaxis()->FindBin(eta + 0.5*dEta - 1e-6);
                int byLo = hTow->GetYaxis()->FindBin(phi - 0.5*dPhi + 1e-6);
                int byHi = hTow->GetYaxis()->FindBin(phi + 0.5*dPhi - 1e-6);
                // Overlapping entries (duplicated towers) keep the largest |E_T|, sign included,
                // rather than summing.
                for (int bx = std::max(1, bxLo); bx <= std::min(nEtaFine, bxHi); ++bx)
                    for (int by = std::max(1, byLo); by <= std::min(nPhiFine, byHi); ++by)
                        if (std::fabs(et) > std::fabs(hTow->GetBinContent(bx, by)))
                            hTow->SetBinContent(bx, by, et);
                if (et > maxEt) maxEt = et;
                if (et < minEt) minEt = et;
            }

            // Fixed negative floor so the negative towers occupy a visible part of the scale
            // rather than being squashed against the bottom of an auto-ranged axis. Anything
            // below the floor is clamped onto it, so a single very negative tower cannot
            // stretch the scale and flatten everything else.
            hTow->SetMinimum(fexTowerZMin);
            if (maxEt > 0.0) hTow->SetMaximum(maxEt);
            for (int bx = 1; bx <= nEtaFine; ++bx)
                for (int by = 1; by <= nPhiFine; ++by)
                    if (hTow->GetBinContent(bx, by) < fexTowerZMin)
                        hTow->SetBinContent(bx, by, fexTowerZMin);
            if (minEt < fexTowerZMin)
                std::cout << "  " << label << " event " << iEvt << ": most negative tower "
                          << minEt << " GeV clamped to the z-axis floor " << fexTowerZMin << " GeV\n";
            hTow->Draw("COLZ");
            drawEtaPhiOverlays();

            TLatex towLat;
            towLat.SetTextSize(0.028);
            towLat.DrawLatexNDC(0.10, 0.95,
                Form("Event %d   (all %s towers, WTA-cone jets > %.1f GeV)",
                     iEvt, label.c_str(), jetEtThreshold));

            cFexTower.Print(pdfName);
            delete hTow;
        };

        // g/jFEX tower displays are made for every event passing the JZ / HSTP selection,
        // i.e. before the background MET over-reconstruction filter below.
        if ((drawGFexTowers || drawJFexTowers) && nFexTowerDisplays < maxFexTowerDisplays) {
            // Jet collections are needed here for the overlays (cheap compared with the towers)
            gepWTAConeEtaSKJetsTree->GetEntry(iEvt);
            truthAntiKt4WZDressedJetsTree->GetEntry(iEvt);
            inTimeAntiKt4TruthJetsTree->GetEntry(iEvt);
            outOfTimeAntiKt4TruthJetsTree->GetEntry(iEvt);
            if (drawGFexTowers) {
                gFexTowersTree->GetEntry(iEvt);
                drawFexTowerDisplay("gFEX", gFexTowerEtValues, gFexTowerEtaValues, gFexTowerPhiValues,
                                    gFexTowerGeom, "gFex", pdf_gTow);
            }
            if (drawJFexTowers) {
                jFexTowersTree->GetEntry(iEvt);
                drawFexTowerDisplay("jFEX", jFexTowerEtValues, jFexTowerEtaValues, jFexTowerPhiValues,
                                    jFexTowerGeom, "jFex", pdf_jTow);
            }
            ++nFexTowerDisplays;
        }

        // Background-only: accept event if EITHER GEP Jet MET or gFEX JwoJ MET is over-reconstructed
        // relative to truth NonInt MET by at least its respective threshold.
        if (!signalBool) {
            bool passJet  = (sig_JetMet - sig_metTruthNonInt) >= backMinJetMETOverTruth;
            bool passGFex = (sig_gMET   - sig_metTruthNonInt) >= backMinGFexMETOverTruth;
            if (!passJet && !passGFex) continue;
        }

        // MET display cap reached — keep looping only to fill the g/jFEX tower PDFs
        if (nAccepted >= maxDisplays) continue;

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

        // Jet circles (GEP WTA-cone, truth, in-time PU, out-of-time PU) and the
        // per-algorithm phi_MET lines — shared with the g/jFEX tower displays
        drawEtaPhiOverlays();

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

        // Bunch-train position: distance in BCID from the tail and the head (front) of the
        // train, i.e. how much out-of-time pileup this event sees from either side.
        y_2 -= 0.03;
        lat.SetTextSize(0.028);
        if (hasBunchTrain) {
            std::string sTail = distTailBunchTrainValues  >= 0 ? Form("%d", distTailBunchTrainValues)  : std::string("n/a");
            std::string sHead = distFrontBunchTrainValues >= 0 ? Form("%d", distFrontBunchTrainValues) : std::string("n/a");
            lat.DrawLatexNDC(0.05, y_2,
                Form("Position in bunch train [BCID] (distance from tail, head): (%s, %s)",
                     sTail.c_str(), sHead.c_str()));
        } else {
            lat.DrawLatexNDC(0.05, y_2,
                "Position in bunch train [BCID] (distance from tail, head): (n/a, n/a)");
        }
        y_2 -= 0.02;

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
    if (drawGFexTowers) cFexTower.Print(pdf_gTow + ")");
    if (drawJFexTowers) cFexTower.Print(pdf_jTow + ")");

    std::cout << "MET displays: " << nAccepted << " events, g/jFEX tower displays: "
              << nFexTowerDisplays << " events\n";
}

// -----------------------------------------------------------------------
void makeMETEventDisplays() {
    gErrorIgnoreLevel = kWarning;

    // Same input paths used by metAnalysisAndRates.C
    const std::string sigHER  = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/ZvvHbb_v4/mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    const std::string backHER = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/QCD_Dijet_JZ*_v4/mc21_14TeV_jj_JZ*_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    const std::string emuDir  = "/data/larsonma/GEPMET/outputNTuplesDev_METv2/";

    // Default signal call: ZvvHbb with the (jetEt=20, towerEt=2, EtaSK_OR, twrSF=1, jetSF=1) config.
    // The last three arguments turn on the standalone gFEX / jFEX tower display PDFs
    // (made for every event, not just the MET-filtered ones) with a 1 GeV display threshold.
    callMakeMETEventDisplays(
        sigHER,
        emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root",
        -1, true, "ZvvHbb",
        /*backMinJetMETOverTruth=*/50.0, /*backMinGFexMETOverTruth=*/150.0,
        /*drawGFexTowers=*/true, /*drawJFexTowers=*/true, /*fexTowerZMin=*/-10.0);

    // Background: JZ1 only. Default filter accepts events where EITHER GEP Jet MET or gFEX JwoJ MET
    // is over-reconstructed vs truth NonInt MET by its respective threshold (50 GeV / 150 GeV).
    // To override, append e.g. /*backMinJetMETOverTruth=*/30.0, /*backMinGFexMETOverTruth=*/100.0
    callMakeMETEventDisplays(
        backHER,
        emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root",
        1, false, "jj_1",
        /*backMinJetMETOverTruth=*/50.0, /*backMinGFexMETOverTruth=*/150.0,
        /*drawGFexTowers=*/true, /*drawJFexTowers=*/true, /*fexTowerZMin=*/-10.0);

    gSystem->Exit(0);
}
