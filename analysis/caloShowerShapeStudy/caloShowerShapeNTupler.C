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
//                                   TrigGepPerf validation).
//   * wtaConeCellsTowersEtaSKTree : WTACone (EtaSK) cells-tower jets, from
//                                   GEPOutputReader, Et-sorted.
//   * gepCellsTowersEtaSKTree     : GEPCellsTowerEtaSK towers incl. per-layer Et
//                                   (Et_l0..Et_l6), from GEPOutputReader.
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
//               "outputDir/","daodFile","gepFile","fileSuffix")'
//
// We'll start with just these blocks and add more later.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <numeric>   // std::iota
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TString.h"

void caloShowerShapeNTupler(bool signalBool,
                            std::string signalString,
                            std::string outputDir,
                            std::string daodFile,
                            std::string gepFile,
                            std::string fileSuffix = "") {

    // ---------------------------------------------------------------
    // Output file
    // ---------------------------------------------------------------
    std::string tag = signalString.empty() ? "background" : signalString;
    std::string outDir = outputDir;
    if (!outDir.empty() && outDir.back() != '/') outDir += "/";
    std::string outName = outDir + "caloShowerShape_" + tag + fileSuffix + ".root";
    TString outputFileName = outName.c_str();
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    std::cout << "[caloShowerShapeNTupler] output: " << outputFileName << "\n";

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
    eventInfoTree->Branch("eventNumber",     &eventNumber);
    eventInfoTree->Branch("runNumber",       &runNumber);
    eventInfoTree->Branch("mcChannelNumber", &mcChannelNumber);
    eventInfoTree->Branch("mcEventWeight",   &mcEventWeight);
    eventInfoTree->Branch("gepEventNumber",  &gepEventNumber);
    eventInfoTree->Branch("gepRunNumber",    &gepRunNumber);
    eventInfoTree->Branch("signal",          &signalFlag);

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

    // (1) JetTaggerLRJ EtaSK
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
        bsm_pdgId.clear(); bsm_status.clear();
        bsm_pt.clear(); bsm_eta.clear(); bsm_phi.clear(); bsm_m.clear(); bsm_e.clear();
        bsm_px.clear(); bsm_py.clear(); bsm_pz.clear();
        bsm_hasProdVtx.clear(); bsm_prodVtx_x.clear(); bsm_prodVtx_y.clear(); bsm_prodVtx_z.clear(); bsm_prodVtx_t.clear();
        bsm_hasDecayVtx.clear(); bsm_decayVtx_x.clear(); bsm_decayVtx_y.clear(); bsm_decayVtx_z.clear(); bsm_decayVtx_t.clear();

        // ---- event info ----
        const xAOD::EventInfo* evtInfo = nullptr;
        if (event.retrieve(evtInfo, "EventInfo").isSuccess()) {
            eventNumber     = evtInfo->eventNumber();
            runNumber       = evtInfo->runNumber();
            mcChannelNumber = evtInfo->mcChannelNumber();
            mcEventWeight   = evtInfo->mcEventWeight();
        } else {
            eventNumber = 0; runNumber = 0; mcChannelNumber = 0; mcEventWeight = 0.;
        }

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

        // ---- fill ----
        eventInfoTree->Fill();
        jetTaggerLRJEtaSKTree->Fill();      // pass-through (lrj_* filled by GetEntry)
        wtaConeCellsTowersEtaSKTree->Fill();
        gepCellsTowersEtaSKTree->Fill();
        truthBSMTree->Fill();
    } // event loop

    // ===============================================================
    // Write + close
    // ===============================================================
    outputFile->cd();
    eventInfoTree->Write("", TObject::kOverwrite);
    jetTaggerLRJEtaSKTree->Write("", TObject::kOverwrite);
    wtaConeCellsTowersEtaSKTree->Write("", TObject::kOverwrite);
    gepCellsTowersEtaSKTree->Write("", TObject::kOverwrite);
    truthBSMTree->Write("", TObject::kOverwrite);
    outputFile->Close();

    gf->Close();
    f->Close();
    std::cout << "[caloShowerShapeNTupler] done: " << outputFileName << "\n";
}
