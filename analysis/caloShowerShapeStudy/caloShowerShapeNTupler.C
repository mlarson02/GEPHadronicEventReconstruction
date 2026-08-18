// caloShowerShapeNTupler.C
// ---------------------------------------------------------------------------
// Focused ntupler for a *displaced jet trigger based on shower shape* study.
//
// Adapted from HERNTupler.C, but trimmed down to only the inputs needed here.
// It reads one DAOD_JETM42 (xAOD) file paired with its matching GEPOutputReader
// ntuple ("ntuple" TTree) and writes a small set of output trees:
//
//   * jetTaggerLRJEtaSKTree       : JetTaggerLRJ EtaSK jets, taken straight from
//                                   GEPOutputReader (the same objects used for
//                                   TrigGepPerf validation). OFF by default, see
//                                   kWriteJetTaggerLRJ below.
//   * wtaConeCellsTowersEtaSKTree : WTACone (EtaSK) cells-tower jets, from
//                                   GEPOutputReader, Et-sorted.
//   * gepCellsTowersEtaSKTree     : GEPCellsTowerEtaSK towers incl. per-layer Et
//                                   (Et_l0..Et_l6), from GEPOutputReader.
//   * gepCellsTowersTree          : the SAME towers with NO soft killer
//                                   (GEPCellsTower), and
//   * gepCellsTowersSKTree        : with the plain (non-eta-dependent) SK.
//                                   EtaSK's dynamic O(1-2) GeV per-tower threshold
//                                   is right for jet finding but starves the shower
//                                   measurement -- a surviving tower spreads that Et
//                                   over up to 7 layers, leaving a jet with 2-3
//                                   towers and 1-2 lit layers, which is not enough
//                                   for the per-layer-centroid fit. Measure the
//                                   shower on these; keep the jets EtaSK.
//   * truthBSMTree                : truth BSM particles from the DAOD (TruthBSM),
//                                   4-vector + pdgId + production/decay vertices
//                                   (loop the whole collection; filter by pdgId
//                                   offline for displaced dark photon / emerging
//                                   jets).
//   * eventInfoTree               : DAOD event/run/mcChannel + weight, plus the
//                                   GEP event/run numbers for cross-checking that
//                                   the DAOD and GEP files are event-aligned.
//
// Run interpreted, exactly like HERNTupler.C (xAOD types resolve via
// autoloading in an AnalysisBase environment):
//   root -b -q 'caloShowerShapeNTupler.C(signalBool,"signalString",
//               "outputDir/","daodFile","gepFile","fileSuffix",jzSlice,pileup)'
//
// QCD dijet background is produced one JZ slice at a time (jzSlice 0..9), one
// merged ntuple per slice, and the downstream macros chain the ten of them. So
// each background event is stamped with its slice index, its cross-section
// weight and its HSTP filter decision -- eventInfoTree branches sampleJZSlice /
// eventWeights / passHSTP, named as in HERNTupler.C. Signal keeps
// sampleJZSlice = -1, eventWeights = {mcEventWeight}, passHSTP = true.
//
// We'll start with just these blocks and add more later.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <numeric>   // std::iota
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <set>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TString.h"
#include "jzSliceWeights.h"

// ---------------------------------------------------------------------------
// Shower-parent selection + debug.
// The shower-initiating LLPs (dark photon Zd -> e+ e-, diagonal dark pion 4900111
// -> q qbar) decay to Standard-Model products at a displaced vertex; that decay
// vertex is the displacement label for the regression model. Those decay vertices
// are NOT in TruthBSM but ARE in TruthBSMWithDecayParticles (BSM particles + their
// decay products) -- which is what showerParentTree reads.
//   * kShowerParentPdgIds empty -> auto-select: a BSM parent (|pdgId|>=4900000 or
//     Zd=32) with a decay vertex that decays to >=1 visible SM particle.
//   * kShowerParentPdgIds set   -> select exactly those |pdgId| (with a decay vtx).
// Read pdgIds/children off the "[BSMwithDecay dump]" printouts to verify.
// ---------------------------------------------------------------------------
static const std::set<int> kShowerParentPdgIds = {
    // optional explicit override, e.g. 32 (Zd), 4900111 (diagonal dark pion)
};
static const int kDebugTruthBSMEvents = 5;   // dump decay chains for the first N events (0 = off)

// The study is WTACone-only for now: the JetTaggerLRJ collection is empty in the
// JZ1 (v22 PU200) production and the large-R side is not what we are measuring, so
// its tree is neither bound nor written. Flip this on to bring it back.
static const bool kWriteJetTaggerLRJ = false;

void caloShowerShapeNTupler(bool signalBool,
                            std::string signalString,
                            std::string outputDir,
                            std::string daodFile,
                            std::string gepFile,
                            std::string fileSuffix = "",
                            int jzSlice = -1,
                            unsigned int pileup = 200) {

    // ---------------------------------------------------------------
    // Output file
    // ---------------------------------------------------------------
    // Background gets a per-slice tag so the ten JZ ntuples sit side by side in
    // the output directory and can be chained with a caloShowerShape_dijet_JZ[0-9]
    // glob downstream.
    std::string tag = signalString.empty() ? "background" : signalString;
    if (!signalBool && signalString.empty() && jzSlice >= 0)
        tag = "dijet_JZ" + std::to_string(jzSlice);
    if (!signalBool && jzSlice < 0)
        std::cerr << "[caloShowerShapeNTupler] WARNING: background run without a JZ slice "
                     "(jzSlice=-1) -- eventWeights fall back to mcEventWeight and the ntuple "
                     "cannot be mixed with other slices.\n";
    std::string outDir = outputDir;
    if (!outDir.empty() && outDir.back() != '/') outDir += "/";
    std::string outName = outDir + "caloShowerShape_" + tag + fileSuffix + ".root";
    TString outputFileName = outName.c_str();
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    std::cout << "[caloShowerShapeNTupler] output: " << outputFileName << "\n";
    if (!signalBool && jzSlice >= 0 && jzSlice < (int)JZSliceWeights::nJZSlices)
        std::cout << "[norm] JZ" << jzSlice << " (PU" << pileup << "): xsec x filterEff = "
                  << JZSliceWeights::crossSectionsByJZSlice[jzSlice]
                     * JZSliceWeights::filterEffienciesByJZSlice[jzSlice] << " b, sumOfWeights = "
                  << JZSliceWeights::sumOfEventWeightsForPU(jzSlice, pileup) << ", L = "
                  << JZSliceWeights::reweightLuminosityForPU(pileup) << " b^-1\n";

    // ===============================================================
    // Output trees, output vectors, and (for GEP inputs) the pointers
    // bound to the GEPOutputReader "ntuple" tree.
    // ===============================================================

    // ---- event info -----------------------------------------------
    TTree* eventInfoTree = new TTree("eventInfoTree", "event info + GEP sync");
    ULong64_t     eventNumber     = 0;
    unsigned int  runNumber       = 0;
    unsigned int  mcChannelNumber = 0;
    double        mcEventWeight   = 0.;
    int           gepEventNumber  = 0;   // from GEP ntuple (int branch)
    int           gepRunNumber    = 0;
    int           signalFlag      = signalBool ? 1 : 0;
    // JZ bookkeeping (background): slice index, the cross-section weight that makes
    // a JZ0-9 chain a physical mixture, and the HSTP filter decision. Signal writes
    // sampleJZSlice = -1, eventWeights = {mcEventWeight}, passHSTP = true.
    int                 sampleJZSlice = signalBool ? -1 : jzSlice;
    std::vector<double> eventWeights;
    bool                passHSTP      = true;
    bool                filterHSTP    = true;   // GEP ntuple branch (background only)
    eventInfoTree->Branch("eventNumber",     &eventNumber);
    eventInfoTree->Branch("runNumber",       &runNumber);
    eventInfoTree->Branch("mcChannelNumber", &mcChannelNumber);
    eventInfoTree->Branch("mcEventWeight",   &mcEventWeight);
    eventInfoTree->Branch("gepEventNumber",  &gepEventNumber);
    eventInfoTree->Branch("gepRunNumber",    &gepRunNumber);
    eventInfoTree->Branch("signal",          &signalFlag);
    eventInfoTree->Branch("sampleJZSlice",   &sampleJZSlice);
    eventInfoTree->Branch("eventWeights",    &eventWeights);
    eventInfoTree->Branch("passHSTP",        &passHSTP);

    // ---- (1) JetTaggerLRJ EtaSK jets (GEPOutputReader) -------------
    // Pass-through: each pointer is bound BOTH to the GEP input branch and to
    // the output branch, so gt->GetEntry() fills it and tree->Fill() writes it.
    TTree* jetTaggerLRJEtaSKTree = new TTree("jetTaggerLRJEtaSKTree", "JetTaggerLRJ EtaSK jets (GEPOutputReader)");
    std::vector<float>* lrj_Et         = new std::vector<float>();
    std::vector<float>* lrj_eta        = new std::vector<float>();
    std::vector<float>* lrj_phi        = new std::vector<float>();
    std::vector<float>* lrj_m          = new std::vector<float>();
    std::vector<int>*   lrj_nSubjets   = new std::vector<int>();
    std::vector<float>* lrj_psi_R      = new std::vector<float>();
    std::vector<float>* lrj_tau_1      = new std::vector<float>();
    std::vector<float>* lrj_tau_2      = new std::vector<float>();
    std::vector<float>* lrj_tau_21     = new std::vector<float>();
    std::vector<float>* lrj_massApprox = new std::vector<float>();
    std::vector<float>* lrj_subjet_Et  = new std::vector<float>();
    std::vector<float>* lrj_subjet_eta = new std::vector<float>();
    std::vector<float>* lrj_subjet_phi = new std::vector<float>();
    jetTaggerLRJEtaSKTree->Branch("Et",         &lrj_Et);
    jetTaggerLRJEtaSKTree->Branch("Eta",        &lrj_eta);
    jetTaggerLRJEtaSKTree->Branch("Phi",        &lrj_phi);
    jetTaggerLRJEtaSKTree->Branch("M",          &lrj_m);
    jetTaggerLRJEtaSKTree->Branch("NSubjets",   &lrj_nSubjets);
    jetTaggerLRJEtaSKTree->Branch("Psi_R",      &lrj_psi_R);
    jetTaggerLRJEtaSKTree->Branch("Tau_1",      &lrj_tau_1);
    jetTaggerLRJEtaSKTree->Branch("Tau_2",      &lrj_tau_2);
    jetTaggerLRJEtaSKTree->Branch("Tau_21",     &lrj_tau_21);
    jetTaggerLRJEtaSKTree->Branch("MassApprox", &lrj_massApprox);
    jetTaggerLRJEtaSKTree->Branch("SubjetEt",   &lrj_subjet_Et);
    jetTaggerLRJEtaSKTree->Branch("SubjetEta",  &lrj_subjet_eta);
    jetTaggerLRJEtaSKTree->Branch("SubjetPhi",  &lrj_subjet_phi);

    // ---- (2) WTACone EtaSK cells-tower jets (GEPOutputReader) ------
    TTree* wtaConeCellsTowersEtaSKTree = new TTree("wtaConeCellsTowersEtaSKTree", "WTACone EtaSK cells-tower jets (GEPOutputReader)");
    // output vectors (Et-sorted); pt in GeV, ring Et kept in GEP units (as HERNTupler)
    std::vector<double> wta_pt, wta_eta, wta_phi;
    std::vector<int>    wta_nConstituents;
    std::vector<double> wta_ring0Et, wta_ring1Et, wta_ring2Et, wta_ring3Et, wta_ring4Et;
    std::vector<int>    wta_totalTobN, wta_ring0TobN, wta_ring1TobN, wta_ring2TobN, wta_ring3TobN, wta_ring4TobN;
    wtaConeCellsTowersEtaSKTree->Branch("Pt",            &wta_pt);
    wtaConeCellsTowersEtaSKTree->Branch("Eta",           &wta_eta);
    wtaConeCellsTowersEtaSKTree->Branch("Phi",           &wta_phi);
    wtaConeCellsTowersEtaSKTree->Branch("NConstituents", &wta_nConstituents);
    wtaConeCellsTowersEtaSKTree->Branch("Ring0Et",       &wta_ring0Et);
    wtaConeCellsTowersEtaSKTree->Branch("Ring1Et",       &wta_ring1Et);
    wtaConeCellsTowersEtaSKTree->Branch("Ring2Et",       &wta_ring2Et);
    wtaConeCellsTowersEtaSKTree->Branch("Ring3Et",       &wta_ring3Et);
    wtaConeCellsTowersEtaSKTree->Branch("Ring4Et",       &wta_ring4Et);
    wtaConeCellsTowersEtaSKTree->Branch("TotalTobN",     &wta_totalTobN);
    wtaConeCellsTowersEtaSKTree->Branch("Ring0TobN",     &wta_ring0TobN);
    wtaConeCellsTowersEtaSKTree->Branch("Ring1TobN",     &wta_ring1TobN);
    wtaConeCellsTowersEtaSKTree->Branch("Ring2TobN",     &wta_ring2TobN);
    wtaConeCellsTowersEtaSKTree->Branch("Ring3TobN",     &wta_ring3TobN);
    wtaConeCellsTowersEtaSKTree->Branch("Ring4TobN",     &wta_ring4TobN);
    // input pointers bound to the GEP tree
    std::vector<float>* in_wta_pt = nullptr, *in_wta_eta = nullptr, *in_wta_phi = nullptr, *in_wta_nc = nullptr;
    std::vector<float>* in_wta_r0 = nullptr, *in_wta_r1 = nullptr, *in_wta_r2 = nullptr, *in_wta_r3 = nullptr, *in_wta_r4 = nullptr;
    std::vector<int>*   in_wta_totN = nullptr, *in_wta_r0N = nullptr, *in_wta_r1N = nullptr, *in_wta_r2N = nullptr, *in_wta_r3N = nullptr, *in_wta_r4N = nullptr;

    // ---- (3) GEPCellsTowerEtaSK (incl. per-layer Et) --------------
    TTree* gepCellsTowersEtaSKTree = new TTree("gepCellsTowersEtaSKTree", "GEPCellsTowerEtaSK towers incl. per-layer Et (GEPOutputReader)");
    std::vector<double> tow_Et, tow_Eta, tow_Phi;
    std::vector<double> tow_Et_l[7];
    gepCellsTowersEtaSKTree->Branch("Et",  &tow_Et);
    gepCellsTowersEtaSKTree->Branch("Eta", &tow_Eta);
    gepCellsTowersEtaSKTree->Branch("Phi", &tow_Phi);
    for (int l = 0; l < 7; ++l)
        gepCellsTowersEtaSKTree->Branch(Form("Et_l%d", l), &tow_Et_l[l]);
    // input pointers bound to the GEP tree
    std::vector<float>* in_tow_et = nullptr, *in_tow_eta = nullptr, *in_tow_phi = nullptr;
    std::vector<float>* in_tow_et_l[7] = { nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr };

    // ---- (3b) the same towers WITHOUT the eta-dependent soft killer ----------
    // EtaSK applies a dynamic O(1-2) GeV per-tower threshold from the event's energy
    // density. That is right for jet finding and rates, but it starves the shower
    // measurement: the per-layer centroid fit needs many towers per LAYER, and a
    // surviving 1-2 GeV tower spreads its Et over up to 7 layers. Writing the
    // unsuppressed (GEPCellsTower) and plain-SK (GEPCellsTowerSK) collections lets
    // the shower shape / pointing be measured on all towers while the trigger jets
    // stay EtaSK. Both are pass-throughs of the same shape as the EtaSK tree above.
    struct TowerColl {
        const char* treeName;
        const char* branchPrefix;   // GEP branch prefix, e.g. "GEPCellsTower"
        TTree*      tree;
        std::vector<double> Et, Eta, Phi;
        std::vector<double> Et_l[7];
        std::vector<float>* in_et;
        std::vector<float>* in_eta;
        std::vector<float>* in_phi;
        std::vector<float>* in_et_l[7];
    };
    // Value-initialized (new TowerColl()) so every input pointer starts null; the
    // interpreter is happier with this than with aggregate-initializing a struct
    // that mixes vectors, arrays and pointers.
    std::vector<TowerColl*> extraTowerColls;
    {
        const char* names[2]    = { "gepCellsTowersTree", "gepCellsTowersSKTree" };
        const char* prefixes[2] = { "GEPCellsTower",      "GEPCellsTowerSK"      };
        for (int i = 0; i < 2; ++i) {
            TowerColl* tc = new TowerColl();
            tc->treeName     = names[i];
            tc->branchPrefix = prefixes[i];
            tc->tree = new TTree(tc->treeName,
                                 Form("%s towers incl. per-layer Et (GEPOutputReader)", tc->branchPrefix));
            tc->tree->Branch("Et",  &tc->Et);
            tc->tree->Branch("Eta", &tc->Eta);
            tc->tree->Branch("Phi", &tc->Phi);
            for (int l = 0; l < 7; ++l) tc->tree->Branch(Form("Et_l%d", l), &tc->Et_l[l]);
            extraTowerColls.push_back(tc);
        }
    }

    // ---- (4) truth BSM particles (DAOD TruthBSM) ------------------
    TTree* truthBSMTree = new TTree("truthBSMTree", "Truth BSM particles (DAOD TruthBSM)");
    std::vector<int>    bsm_pdgId, bsm_status;
    std::vector<double> bsm_pt, bsm_eta, bsm_phi, bsm_m, bsm_e, bsm_px, bsm_py, bsm_pz;
    std::vector<double> bsm_prodVtx_x, bsm_prodVtx_y, bsm_prodVtx_z, bsm_prodVtx_t;
    std::vector<double> bsm_decayVtx_x, bsm_decayVtx_y, bsm_decayVtx_z, bsm_decayVtx_t;
    std::vector<int>    bsm_hasProdVtx, bsm_hasDecayVtx;
    truthBSMTree->Branch("pdgId",       &bsm_pdgId);
    truthBSMTree->Branch("status",      &bsm_status);
    truthBSMTree->Branch("pt",          &bsm_pt);
    truthBSMTree->Branch("eta",         &bsm_eta);
    truthBSMTree->Branch("phi",         &bsm_phi);
    truthBSMTree->Branch("m",           &bsm_m);
    truthBSMTree->Branch("e",           &bsm_e);
    truthBSMTree->Branch("px",          &bsm_px);
    truthBSMTree->Branch("py",          &bsm_py);
    truthBSMTree->Branch("pz",          &bsm_pz);
    truthBSMTree->Branch("hasProdVtx",  &bsm_hasProdVtx);
    truthBSMTree->Branch("prodVtx_x",   &bsm_prodVtx_x);
    truthBSMTree->Branch("prodVtx_y",   &bsm_prodVtx_y);
    truthBSMTree->Branch("prodVtx_z",   &bsm_prodVtx_z);
    truthBSMTree->Branch("prodVtx_t",   &bsm_prodVtx_t);
    truthBSMTree->Branch("hasDecayVtx", &bsm_hasDecayVtx);
    truthBSMTree->Branch("decayVtx_x",  &bsm_decayVtx_x);
    truthBSMTree->Branch("decayVtx_y",  &bsm_decayVtx_y);
    truthBSMTree->Branch("decayVtx_z",  &bsm_decayVtx_z);
    truthBSMTree->Branch("decayVtx_t",  &bsm_decayVtx_t);

    // ---- (5) shower-parent LLPs (selected subset of TruthBSM) --------------
    // Only the BSM particles that seed calorimeter showers (kShowerParentPdgIds);
    // their decay vertex is the displacement label for the regression model.
    TTree* showerParentTree = new TTree("showerParentTree", "Shower-initiating BSM parents (decay vertex = label)");
    std::vector<int>    sp_pdgId, sp_status, sp_hasDecayVtx;
    std::vector<double> sp_pt, sp_eta, sp_phi, sp_m;
    std::vector<double> sp_decayVtx_x, sp_decayVtx_y, sp_decayVtx_z;
    std::vector<double> sp_decayVtx_Lxy, sp_decayVtx_r3d;   // convenience (mm)
    showerParentTree->Branch("pdgId",        &sp_pdgId);
    showerParentTree->Branch("status",       &sp_status);
    showerParentTree->Branch("pt",           &sp_pt);
    showerParentTree->Branch("eta",          &sp_eta);
    showerParentTree->Branch("phi",          &sp_phi);
    showerParentTree->Branch("m",            &sp_m);
    showerParentTree->Branch("hasDecayVtx",  &sp_hasDecayVtx);
    showerParentTree->Branch("decayVtx_x",   &sp_decayVtx_x);
    showerParentTree->Branch("decayVtx_y",   &sp_decayVtx_y);
    showerParentTree->Branch("decayVtx_z",   &sp_decayVtx_z);
    showerParentTree->Branch("decayVtx_Lxy", &sp_decayVtx_Lxy);
    showerParentTree->Branch("decayVtx_r3d", &sp_decayVtx_r3d);

    // ===============================================================
    // Open inputs: DAOD (xAOD) + GEP ntuple
    // ===============================================================
    TFile* f = TFile::Open(daodFile.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "[caloShowerShapeNTupler] could not open DAOD: " << daodFile << "\n";
        return;
    }
    xAOD::TEvent event(xAOD::TEvent::kClassAccess);
    if (!event.readFrom(f).isSuccess()) {
        std::cerr << "[caloShowerShapeNTupler] cannot read xAOD from: " << daodFile << "\n";
        return;
    }

    TFile* gf = TFile::Open(gepFile.c_str());
    if (!gf || gf->IsZombie()) {
        std::cerr << "[caloShowerShapeNTupler] could not open GEP file: " << gepFile << "\n";
        return;
    }
    TTree* gt = nullptr;
    gf->GetObject("ntuple", gt);
    if (!gt) {
        std::cerr << "[caloShowerShapeNTupler] no 'ntuple' TTree in GEP file: " << gepFile << "\n";
        return;
    }

    // Read only the branches we need.
    gt->SetBranchStatus("*", 0);
    auto bind = [&](const char* name, auto addr) {
        if (!gt->GetBranch(name)) {
            std::cerr << "[caloShowerShapeNTupler] WARNING: GEP branch missing: " << name << "\n";
            return;
        }
        gt->SetBranchStatus(name, 1);
        gt->SetBranchAddress(name, addr);
    };

    // event sync
    bind("eventNumber", &gepEventNumber);
    bind("runNumber",   &gepRunNumber);

    // HSTP filter (background only): removes the JZ events whose hard-scatter truth
    // particle is out of time, the same cut metAnalysisAndRates.C / largeRJet-
    // AnalysisAndRates.C apply to the background chain.
    const bool haveFilterHSTP = !signalBool && gt->GetBranch("filterHSTP") != nullptr;
    if (haveFilterHSTP) bind("filterHSTP", &filterHSTP);
    else if (!signalBool)
        std::cerr << "[caloShowerShapeNTupler] WARNING: no filterHSTP branch in the GEP file -- "
                     "passHSTP will be true for every event.\n";

    // (1) JetTaggerLRJ EtaSK
    if (kWriteJetTaggerLRJ) {
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_Et",         &lrj_Et);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_eta",        &lrj_eta);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_phi",        &lrj_phi);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_m",          &lrj_m);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_nSubjets",   &lrj_nSubjets);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_psi_R",      &lrj_psi_R);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_tau_1",      &lrj_tau_1);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_tau_2",      &lrj_tau_2);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_tau_21",     &lrj_tau_21);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_massApprox", &lrj_massApprox);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_subjet_Et",  &lrj_subjet_Et);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_subjet_eta", &lrj_subjet_eta);
        bind("JetTaggerLRJGEPCellsTowerEtaSKJets_subjet_phi", &lrj_subjet_phi);
    }

    // (2) WTACone EtaSK cells-tower jets
    bind("WTAConeGEPCellsTowerEtaSKJets_pt",            &in_wta_pt);
    bind("WTAConeGEPCellsTowerEtaSKJets_eta",           &in_wta_eta);
    bind("WTAConeGEPCellsTowerEtaSKJets_phi",           &in_wta_phi);
    bind("WTAConeGEPCellsTowerEtaSKJets_nConstituents", &in_wta_nc);
    bind("WTAConeGEPCellsTowerEtaSKJets_ring0_Et",      &in_wta_r0);
    bind("WTAConeGEPCellsTowerEtaSKJets_ring1_Et",      &in_wta_r1);
    bind("WTAConeGEPCellsTowerEtaSKJets_ring2_Et",      &in_wta_r2);
    bind("WTAConeGEPCellsTowerEtaSKJets_ring3_Et",      &in_wta_r3);
    bind("WTAConeGEPCellsTowerEtaSKJets_ring4_Et",      &in_wta_r4);
    bind("WTAConeGEPCellsTowerEtaSKJets_total_TobN",    &in_wta_totN);
    bind("WTAConeGEPCellsTowerEtaSKJets_ring0_TobN",    &in_wta_r0N);
    bind("WTAConeGEPCellsTowerEtaSKJets_ring1_TobN",    &in_wta_r1N);
    bind("WTAConeGEPCellsTowerEtaSKJets_ring2_TobN",    &in_wta_r2N);
    bind("WTAConeGEPCellsTowerEtaSKJets_ring3_TobN",    &in_wta_r3N);
    bind("WTAConeGEPCellsTowerEtaSKJets_ring4_TobN",    &in_wta_r4N);

    // (3) GEPCellsTowerEtaSK (incl. per-layer Et)
    bind("GEPCellsTowerEtaSK_et",  &in_tow_et);
    bind("GEPCellsTowerEtaSK_eta", &in_tow_eta);
    bind("GEPCellsTowerEtaSK_phi", &in_tow_phi);
    for (int l = 0; l < 7; ++l)
        bind(Form("GEPCellsTowerEtaSK_et_l%d", l), &in_tow_et_l[l]);

    // (3b) the unsuppressed + plain-SK tower collections
    for (TowerColl* tc : extraTowerColls) {
        bind(Form("%s_et",  tc->branchPrefix), &tc->in_et);
        bind(Form("%s_eta", tc->branchPrefix), &tc->in_eta);
        bind(Form("%s_phi", tc->branchPrefix), &tc->in_phi);
        for (int l = 0; l < 7; ++l)
            bind(Form("%s_et_l%d", tc->branchPrefix, l), &tc->in_et_l[l]);
    }

    // ===============================================================
    // Event loop
    // ===============================================================
    Long64_t nDaod = event.getEntries();
    Long64_t nGep  = gt->GetEntries();
    if (nDaod != nGep) {
        std::cerr << "[caloShowerShapeNTupler] WARNING: DAOD has " << nDaod
                  << " events but GEP has " << nGep << " -- processing min().\n";
    }
    Long64_t nEvents = std::min(nDaod, nGep);
    std::cout << "[caloShowerShapeNTupler] processing " << nEvents << " events\n";

    for (Long64_t iEvt = 0; iEvt < nEvents; ++iEvt) {
        event.getEntry(iEvt);   // DAOD event
        gt->GetEntry(iEvt);     // GEP  event (fills lrj_*, in_wta_*, in_tow_* )

        // ---- clear per-event output vectors ----
        wta_pt.clear(); wta_eta.clear(); wta_phi.clear(); wta_nConstituents.clear();
        wta_ring0Et.clear(); wta_ring1Et.clear(); wta_ring2Et.clear(); wta_ring3Et.clear(); wta_ring4Et.clear();
        wta_totalTobN.clear(); wta_ring0TobN.clear(); wta_ring1TobN.clear(); wta_ring2TobN.clear(); wta_ring3TobN.clear(); wta_ring4TobN.clear();
        tow_Et.clear(); tow_Eta.clear(); tow_Phi.clear();
        for (int l = 0; l < 7; ++l) tow_Et_l[l].clear();
        for (TowerColl* tc : extraTowerColls) {
            tc->Et.clear(); tc->Eta.clear(); tc->Phi.clear();
            for (int l = 0; l < 7; ++l) tc->Et_l[l].clear();
        }
        bsm_pdgId.clear(); bsm_status.clear();
        bsm_pt.clear(); bsm_eta.clear(); bsm_phi.clear(); bsm_m.clear(); bsm_e.clear();
        bsm_px.clear(); bsm_py.clear(); bsm_pz.clear();
        bsm_hasProdVtx.clear(); bsm_prodVtx_x.clear(); bsm_prodVtx_y.clear(); bsm_prodVtx_z.clear(); bsm_prodVtx_t.clear();
        bsm_hasDecayVtx.clear(); bsm_decayVtx_x.clear(); bsm_decayVtx_y.clear(); bsm_decayVtx_z.clear(); bsm_decayVtx_t.clear();
        sp_pdgId.clear(); sp_status.clear(); sp_hasDecayVtx.clear();
        sp_pt.clear(); sp_eta.clear(); sp_phi.clear(); sp_m.clear();
        sp_decayVtx_x.clear(); sp_decayVtx_y.clear(); sp_decayVtx_z.clear();
        sp_decayVtx_Lxy.clear(); sp_decayVtx_r3d.clear();

        // ---- event info ----
        eventWeights.clear();
        const xAOD::EventInfo* evtInfo = nullptr;
        if (event.retrieve(evtInfo, "EventInfo").isSuccess()) {
            eventNumber     = evtInfo->eventNumber();
            runNumber       = evtInfo->runNumber();
            mcChannelNumber = evtInfo->mcChannelNumber();
            mcEventWeight   = evtInfo->mcEventWeight();
        } else {
            eventNumber = 0; runNumber = 0; mcChannelNumber = 0; mcEventWeight = 0.;
        }
        // Slice weight for background (mcEventWeight passed through for signal, whose
        // sampleJZSlice is -1), and the HSTP decision read from the GEP ntuple.
        eventWeights.push_back(JZSliceWeights::eventWeight(mcEventWeight, sampleJZSlice, pileup));
        passHSTP = haveFilterHSTP ? filterHSTP : true;

        // ---- (2) WTACone EtaSK jets: Et-sorted (descending) ----
        if (in_wta_pt) {
            std::vector<unsigned int> order(in_wta_pt->size());
            std::iota(order.begin(), order.end(), 0u);
            std::sort(order.begin(), order.end(),
                      [&](unsigned int a, unsigned int b){ return (*in_wta_pt)[a] > (*in_wta_pt)[b]; });
            for (unsigned int idx : order) {
                wta_pt.push_back((*in_wta_pt)[idx] / 1000.0);   // GeV
                wta_eta.push_back((*in_wta_eta)[idx]);
                wta_phi.push_back((*in_wta_phi)[idx]);
                wta_nConstituents.push_back(static_cast<int>((*in_wta_nc)[idx]));
                // ring Et kept in the GEP ntuple's native units, matching HERNTupler
                wta_ring0Et.push_back((*in_wta_r0)[idx]);
                wta_ring1Et.push_back((*in_wta_r1)[idx]);
                wta_ring2Et.push_back((*in_wta_r2)[idx]);
                wta_ring3Et.push_back((*in_wta_r3)[idx]);
                wta_ring4Et.push_back((*in_wta_r4)[idx]);
                wta_totalTobN.push_back((*in_wta_totN)[idx]);
                wta_ring0TobN.push_back((*in_wta_r0N)[idx]);
                wta_ring1TobN.push_back((*in_wta_r1N)[idx]);
                wta_ring2TobN.push_back((*in_wta_r2N)[idx]);
                wta_ring3TobN.push_back((*in_wta_r3N)[idx]);
                wta_ring4TobN.push_back((*in_wta_r4N)[idx]);
            }
        }

        // ---- (3b) unsuppressed + plain-SK towers: Et-sorted (descending) ----
        for (TowerColl* tc : extraTowerColls) {
            if (!tc->in_et) continue;
            std::vector<unsigned int> order(tc->in_et->size());
            std::iota(order.begin(), order.end(), 0u);
            std::sort(order.begin(), order.end(),
                      [&](unsigned int a, unsigned int b){ return (*tc->in_et)[a] > (*tc->in_et)[b]; });
            for (unsigned int idx : order) {
                tc->Et.push_back((*tc->in_et)[idx] / 1000.0);   // GeV
                tc->Eta.push_back((*tc->in_eta)[idx]);
                tc->Phi.push_back((*tc->in_phi)[idx]);
                for (int l = 0; l < 7; ++l) {
                    if (tc->in_et_l[l]) tc->Et_l[l].push_back((*tc->in_et_l[l])[idx] / 1000.0);
                    else                tc->Et_l[l].push_back(0.0);
                }
            }
        }

        // ---- (3) GEPCellsTowerEtaSK: Et-sorted (descending) ----
        if (in_tow_et) {
            std::vector<unsigned int> order(in_tow_et->size());
            std::iota(order.begin(), order.end(), 0u);
            std::sort(order.begin(), order.end(),
                      [&](unsigned int a, unsigned int b){ return (*in_tow_et)[a] > (*in_tow_et)[b]; });
            for (unsigned int idx : order) {
                tow_Et.push_back((*in_tow_et)[idx] / 1000.0);   // GeV
                tow_Eta.push_back((*in_tow_eta)[idx]);
                tow_Phi.push_back((*in_tow_phi)[idx]);
                for (int l = 0; l < 7; ++l) {
                    if (in_tow_et_l[l]) tow_Et_l[l].push_back((*in_tow_et_l[l])[idx] / 1000.0);
                    else                tow_Et_l[l].push_back(0.0);
                }
            }
        }

        // ---- DEBUG (event 0 only): TruthBSM carries NO decay vertices here
        //      (hasDecayVtx() is false for every particle -> decayLxy prints -1).
        //      Probe which truth containers exist and which hold the LLP decay
        //      vertices we need for the displacement label. ----
        if (kDebugTruthBSMEvents > 0 && iEvt == 0) {
            std::cout << "\n[truth containers present in this DAOD]\n";
            for (const char* name : { "TruthParticles", "TruthBSM", "TruthElectrons",
                                      "TruthMuons", "TruthPhotons", "TruthTaus",
                                      "TruthBSMWithDecayParticles", "TruthBSMWithDecayVertices",
                                      "TruthLLP", "HardScatterParticles", "TruthHFWithDecayParticles" }) {
                const xAOD::TruthParticleContainer* c = nullptr;
                if (event.retrieve(c, name).isSuccess() && c) {
                    int nDec = 0;
                    for (const auto* p : *c) if (p && p->hasDecayVtx()) ++nDec;
                    std::cout << "  particles '" << name << "': " << c->size()
                              << "  (with decayVtx: " << nDec << ")\n";
                }
            }
            for (const char* name : { "TruthVertices", "TruthPrimaryVertices",
                                      "TruthDisplacedVertices", "TruthBSMVertices" }) {
                const xAOD::TruthVertexContainer* c = nullptr;
                if (event.retrieve(c, name).isSuccess() && c)
                    std::cout << "  vertices  '" << name << "': " << c->size() << "\n";
            }
            // If the full TruthParticles record is present, show the LLP decay vertices.
            const xAOD::TruthParticleContainer* tp = nullptr;
            if (event.retrieve(tp, "TruthParticles").isSuccess() && tp) {
                std::cout << "  [TruthParticles LLP decay vertices]\n";
                int shown = 0;
                for (const auto* p : *tp) {
                    if (!p) continue;
                    int a = p->pdgId(); if (a < 0) a = -a;
                    bool llp = (a == 32 || a == 4900111 || a == 4900113 || a == 4900211 || a == 4900213);
                    if (llp && p->hasDecayVtx() && p->decayVtx() && shown < 20) {
                        const auto* v = p->decayVtx();
                        std::cout << "    pdgId=" << p->pdgId() << " status=" << p->status()
                                  << " decayLxy=" << std::sqrt(v->x() * v->x() + v->y() * v->y())
                                  << "mm decayZ=" << v->z() << "mm nChildren=" << v->nOutgoingParticles() << "\n";
                        ++shown;
                    }
                }
            }
        }

        // ---- (4) truth BSM particles (loop whole collection) ----
        const xAOD::TruthParticleContainer* truthBSM = nullptr;
        if (event.retrieve(truthBSM, "TruthBSM").isSuccess()) {
            for (const auto* p : *truthBSM) {
                if (!p) continue;
                bsm_pdgId.push_back(p->pdgId());
                bsm_status.push_back(p->status());
                bsm_pt.push_back(p->pt() / 1000.0);   // GeV
                bsm_eta.push_back(p->eta());
                bsm_phi.push_back(p->phi());
                bsm_m.push_back(p->m() / 1000.0);
                bsm_e.push_back(p->e() / 1000.0);
                bsm_px.push_back(p->px() / 1000.0);
                bsm_py.push_back(p->py() / 1000.0);
                bsm_pz.push_back(p->pz() / 1000.0);

                if (p->hasProdVtx() && p->prodVtx()) {
                    const auto* v = p->prodVtx();
                    bsm_hasProdVtx.push_back(1);
                    bsm_prodVtx_x.push_back(v->x());   // mm
                    bsm_prodVtx_y.push_back(v->y());
                    bsm_prodVtx_z.push_back(v->z());
                    bsm_prodVtx_t.push_back(v->t());
                } else {
                    bsm_hasProdVtx.push_back(0);
                    bsm_prodVtx_x.push_back(-999.); bsm_prodVtx_y.push_back(-999.);
                    bsm_prodVtx_z.push_back(-999.); bsm_prodVtx_t.push_back(-999.);
                }

                if (p->hasDecayVtx() && p->decayVtx()) {
                    const auto* v = p->decayVtx();
                    bsm_hasDecayVtx.push_back(1);
                    bsm_decayVtx_x.push_back(v->x());   // mm
                    bsm_decayVtx_y.push_back(v->y());
                    bsm_decayVtx_z.push_back(v->z());
                    bsm_decayVtx_t.push_back(v->t());
                } else {
                    bsm_hasDecayVtx.push_back(0);
                    bsm_decayVtx_x.push_back(-999.); bsm_decayVtx_y.push_back(-999.);
                    bsm_decayVtx_z.push_back(-999.); bsm_decayVtx_t.push_back(-999.);
                }
            }
        } else if (iEvt == 0) {
            std::cerr << "[caloShowerShapeNTupler] WARNING: TruthBSM not available -- truthBSMTree will be empty.\n";
        }

        // ---- (5) shower-parent LLPs from TruthBSMWithDecayParticles ----
        // TruthBSM has no decay vertices; TruthBSMWithDecayParticles keeps the BSM
        // particles AND their decay products, so the LLP decay vertex (= our label)
        // and the SM children resolve here.
        const xAOD::TruthParticleContainer* bsmDecay = nullptr;
        bool haveBsmDecay = event.retrieve(bsmDecay, "TruthBSMWithDecayParticles").isSuccess() && bsmDecay;

        // parent is BSM (dark sector 4900xxx, or Zd=32); child is a visible SM particle.
        auto isBSM = [](int pid){ int a = pid < 0 ? -pid : pid; return a >= 4900000 || a == 32; };
        auto isVisibleSM = [](int pid){ int a = pid < 0 ? -pid : pid;
            return a > 0 && a < 4900000 && a != 12 && a != 14 && a != 16; };  // exclude neutrinos
        auto seedsShower = [&](const xAOD::TruthParticle* p) -> bool {
            if (!p->hasDecayVtx() || !p->decayVtx()) return false;
            int a = p->pdgId(); if (a < 0) a = -a;
            if (!kShowerParentPdgIds.empty()) return kShowerParentPdgIds.count(a) > 0;
            if (!isBSM(p->pdgId())) return false;   // parent must be a BSM particle
            const auto* dv = p->decayVtx();
            for (size_t c = 0; c < dv->nOutgoingParticles(); ++c)
                if (dv->outgoingParticle(c) && isVisibleSM(dv->outgoingParticle(c)->pdgId())) return true;
            return false;                            // decays only to dark sector / invisible
        };

        if (haveBsmDecay) {
            for (const auto* p : *bsmDecay) {
                if (!p || !seedsShower(p)) continue;
                const auto* v = p->decayVtx();
                double dx = v->x(), dy = v->y(), dz = v->z();   // mm
                sp_pdgId.push_back(p->pdgId());
                sp_status.push_back(p->status());
                sp_pt.push_back(p->pt() / 1000.0);              // GeV
                sp_eta.push_back(p->eta());
                sp_phi.push_back(p->phi());
                sp_m.push_back(p->m() / 1000.0);
                sp_hasDecayVtx.push_back(1);
                sp_decayVtx_x.push_back(dx);
                sp_decayVtx_y.push_back(dy);
                sp_decayVtx_z.push_back(dz);
                sp_decayVtx_Lxy.push_back(std::sqrt(dx * dx + dy * dy));
                sp_decayVtx_r3d.push_back(std::sqrt(dx * dx + dy * dy + dz * dz));
            }
        } else if (iEvt == 0) {
            std::cerr << "[caloShowerShapeNTupler] WARNING: TruthBSMWithDecayParticles not found -- "
                         "showerParentTree will be empty.\n";
        }

        // DEBUG: dump the decaying particles (LLP candidates) with decay vertex +
        // children pdgIds; [SEED] marks the selected shower parents.
        if (kDebugTruthBSMEvents > 0 && iEvt < kDebugTruthBSMEvents && haveBsmDecay) {
            std::cout << "\n[BSMwithDecay dump] evt " << iEvt << "  (" << bsmDecay->size()
                      << " particles; showing those with a decay vertex)\n";
            for (const auto* p : *bsmDecay) {
                if (!p || !p->hasDecayVtx() || !p->decayVtx()) continue;
                const auto* v = p->decayVtx();
                std::cout << "  pdgId=" << p->pdgId() << " status=" << p->status()
                          << " pt=" << p->pt() / 1000.0 << "GeV"
                          << " decayLxy=" << std::sqrt(v->x() * v->x() + v->y() * v->y())
                          << "mm decayZ=" << v->z() << "mm" << (seedsShower(p) ? " [SEED]" : "")
                          << "  children:";
                for (size_t c = 0; c < v->nOutgoingParticles(); ++c)
                    if (v->outgoingParticle(c)) std::cout << " " << v->outgoingParticle(c)->pdgId();
                std::cout << "\n";
            }
        }

        // ---- fill ----
        eventInfoTree->Fill();
        if (kWriteJetTaggerLRJ) jetTaggerLRJEtaSKTree->Fill();   // pass-through (lrj_* filled by GetEntry)
        wtaConeCellsTowersEtaSKTree->Fill();
        gepCellsTowersEtaSKTree->Fill();
        for (TowerColl* tc : extraTowerColls) tc->tree->Fill();
        truthBSMTree->Fill();
        showerParentTree->Fill();
    } // event loop

    // ===============================================================
    // Write + close
    // ===============================================================
    outputFile->cd();
    eventInfoTree->Write("", TObject::kOverwrite);
    if (kWriteJetTaggerLRJ) jetTaggerLRJEtaSKTree->Write("", TObject::kOverwrite);
    wtaConeCellsTowersEtaSKTree->Write("", TObject::kOverwrite);
    gepCellsTowersEtaSKTree->Write("", TObject::kOverwrite);
    for (TowerColl* tc : extraTowerColls) tc->tree->Write("", TObject::kOverwrite);
    truthBSMTree->Write("", TObject::kOverwrite);
    showerParentTree->Write("", TObject::kOverwrite);
    outputFile->Close();

    gf->Close();
    f->Close();
    std::cout << "[caloShowerShapeNTupler] done: " << outputFileName << "\n";
}
