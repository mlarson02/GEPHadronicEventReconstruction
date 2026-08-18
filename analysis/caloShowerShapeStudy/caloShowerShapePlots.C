// caloShowerShapePlots.C
// ---------------------------------------------------------------------------
// Aggregate per-layer shower-shape distributions for the displaced-jet study.
// This is the *quantitative* companion to caloShowerEventDisplays.C: for each
// jet (ΔR-matched to a displaced truth-BSM particle in signal, or the two
// leading jets in dijet background) it sums the constituent-tower Et in each
// calorimeter layer l0..l6 and forms the longitudinal energy profile. Signal
// and dijet are overlaid so the layer that discriminates a displaced/late
// shower stands out.
//
//   # both signals, once each vs the QCD JZ chain, into plots/<signal>/:
//   root -b -l -q 'caloShowerShapePlots.C'
//   # or one sample explicitly (optionally vs a dijet file or JZ glob):
//   root -b -q 'caloShowerShapePlots.C("signal.root","caloShowerShape_dijet_JZ[0-9].root","plots/")'
//
// The dijet background is one ntuple per JZ slice, so any input path may be a
// glob: ChainSource (../chainSource.h) turns it into one TChain per tree, the
// same way metAnalysisAndRates.C / largeRJetAnalysisAndRates.C read the ten v4
// slices. Background jets are filled with their per-event cross-section weight
// (eventInfoTree/eventWeights, written by caloShowerShapeNTupler.C) so the ten
// slices combine into a physical mixture rather than a raw event pile.
//
// One PDF per jet collection:
//   <outDir>/caloShowerShapePlots_WTACone.pdf
//   <outDir>/caloShowerShapePlots_LRJ.pdf
//
// Distributions (signal vs dijet, area-normalized):
//   * Et fraction in each layer l0..l6
//   * longitudinal shower "depth" = Et-weighted mean layer index (0..6)
//   * EM fraction  = (l0+l1+l2+l3)/total   [displaced/hadronic -> lower]
//   * n constituent towers per jet
//   * shower-line DCA3D (mm) of the leading and subleading jet, one pad each:
//     the closest approach to the IP of the Et-weighted per-layer-centroid line
//     fit (caloShowerPointing.h, the same line the event displays draw). A prompt
//     shower's line passes through the IP; one from a displaced decay misses it.
//   * <DCA3D> profiled against |eta|, the number of layers the fit used, and jet Et
//     (page 8) -- the diagnostics that separate a geometry artefact from physics.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TSystem.h"
#include "../chainSource.h"
#include "caloShowerPointing.h"

// Default QCD dijet background: the ten per-slice ntuples, chained.
// The [0-9].root ending is load-bearing: the per-job outputs
// (caloShowerShape_dijet_JZ9_000510.root) live in the same flat directory and
// survive the hadd, so a bare JZ*.root would chain all ~1700 of them -- every
// event twice, once per-job and once merged.
static const char* kDijetJZGlob =
    "/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_dijet_JZ[0-9].root";

// Weight background jets by their JZ cross-section weight. Off = every event
// counts once, which is only meaningful within a single slice.
static const bool kApplyJZWeights = true;
// Drop background events failing the HSTP filter, as the rate macros do.
static const bool kApplyHSTPFilter = true;

static const int kNJZSlices = 10;

// Tower collection feeding the layer sums and the pointing fit:
//   "gepCellsTowersTree"       - no soft killer (DEFAULT; the fit needs the towers)
//   "gepCellsTowersSKTree"     - plain SK
//   "gepCellsTowersEtaSKTree"  - eta-dependent SK (the trigger jets' own towers)
// EtaSK's dynamic O(1-2) GeV per-tower threshold is right for jet finding but leaves
// a QCD jet with 2-3 towers over 1-2 layers -- not enough for a per-layer-centroid
// fit, which is what produced the large-DCA3D background tail. The jets stay EtaSK
// either way; only the shower measurement moves. Falls back to the EtaSK tree with a
// warning for ntuples produced before the extra collections were added.
static const char* kTowerTree         = "gepCellsTowersTree";
static const char* kTowerTreeFallback = "gepCellsTowersEtaSKTree";

// Jet collections to process. WTACone only by default: each collection costs a
// full pass over the ~850k-event JZ chain, and the whole macro is ~20 min per pass
// interpreted (run it with ACLiC, caloShowerShapePlots.C+, to cut that hugely).
static const bool kRunLRJ = false;

// Cap on events read per sample; -1 = all. The chained JZ background is far larger
// than either signal, so this is the knob to turn for a quick look. Note that a cap
// truncates the chain in file order, i.e. it drops the highest slices entirely.
static const Long64_t kMaxEventsPerSample = -1;

// Shower-pointing dca3D axis (mm): the closest approach to the IP of the Et-weighted
// per-layer-centroid line fit (caloShowerPointing.h). Prompt jets do not sit at
// zero -- tower granularity (0.1 in eta is ~15 cm at the EM barrel radius) sets a
// resolution floor -- so what matters is the shift between signal and background.
static const int    kNDcaBins = 50;
static const double kDcaMaxMM = 2000.0;

// Illustrative displaced working point used for the rate view: the point of the
// study is to lower the Et threshold for displaced objects without buying
// background rate, so the background rate curve is drawn both inclusively and
// after requiring the fitted DCA3D above this value. Scan it once real samples exist.
static const double kDcaCutMM = 300.0;

// Barrel-only cross-check: skip jets with |eta| above this (-1 = keep all). The
// nominal geometry switches from fixed-r (barrel) to fixed-z (endcap) at |eta|=1.5,
// and in the endcap r = z/sinh(eta) makes the placed radius hypersensitive to the
// eta centroid, so both effects fake a large DCA3D. Setting 1.2 removes them and
// tells you how much of the background tail is geometry rather than shower physics.
static const double kMaxJetAbsEta = 1.2;

// Leading-jet Et axis for the rate/efficiency curves.
static const int    kNEtBins  = 40;
static const double kEtMaxGeV = 800.0;

static inline double dRp(double e1,double p1,double e2,double p2){
    double dp=std::fabs(p1-p2); while(dp>M_PI) dp=std::fabs(dp-2*M_PI);
    double de=e1-e2; return std::sqrt(de*de+dp*dp);
}

// Keep overflow visible in the last bin, as the rate macros' clampVal does.
static inline double clampToAxis(const TH1F* h, double v){
    const double lo = h->GetXaxis()->GetXmin();
    const double hi = h->GetXaxis()->GetXmax();
    const double eps = 0.5 * h->GetXaxis()->GetBinWidth(h->GetNbinsX());
    return std::min(std::max(v, lo), hi - eps);
}

// Truth + rate observables, grouped so the fill signature stays readable.
// Truth comes from showerParentTree (the LLP that seeded the shower); QCD has no
// shower parents, so its truth histograms stay empty by construction — a prompt
// jet's truth DCA3D is zero.
struct ExtraHists {
    TH1F* truthDca3D[2];  // truth DCA3D (mm), leading / subleading jet
    TH1F* truthLxy;       // decay vertex Lxy (mm)
    TH1F* truthR3d;       // decay vertex |r| (mm)
    TH1F* truthPt;        // matched shower-parent kinematics
    TH1F* truthEta;
    TH1F* truthPhi;
    TH1F* truthMass;
    TH1F* dca3DRes;       // fitted - truth DCA3D (mm), leading jets
    TH2F* dca3DCorr;      // fitted vs truth DCA3D (mm), leading jets
    TH1F* jetEtAll;       // leading-jet Et [GeV], weighted for background
    TH1F* jetEtDca;       // ... after the displaced (dca3D > kDcaCutMM) requirement
    // --- DCA3D diagnostics: what is the fake tail actually correlated with? ---
    TProfile* dcaVsAbsEta;    // <DCA3D> vs |eta|  (expect a step at 1.5 + growth)
    TProfile* dcaVsNLayers;   // <DCA3D> vs number of layers the fit used
    TProfile* dcaVsJetEt;     // <DCA3D> vs jet Et (soft jets fit worst)
    TH2F*     dcaVsAbsEta2D;  // the same correlation, unprofiled
    long  nTruthMatched;
    long  nBeyondEtaCut;      // jets skipped by kMaxJetAbsEta
};

// Fill one sample's shower-shape histograms. `file` may be a glob over the JZ
// slices; every tree is then read as a TChain spanning them.
// hdca3D[0]/hdca3D[1] receive the shower-line DCA3D (mm) of the leading / subleading
// selected jet by Et; ex holds the truth and rate observables.
static long fillShowerHists(const std::string& file, bool isDijet, bool isLRJ, double Rassoc,
                            double ptMinBSM, double etMinTower, bool useTruth,
                            TH1F* hfrac[7], TH1F* hdepth, TH1F* hEMfrac, TH1F* hntow,
                            TH1F* hdca3D[2], ExtraHists& ex) {
    ChainSource* fin = ChainSource::Open(file.c_str());
    if (!fin || fin->IsZombie()) { std::cerr << "cannot open " << file << "\n"; return 0; }

    TTree* towTree = fin->Get(kTowerTree);
    const char* towTreeUsed = kTowerTree;
    if (!towTree) {
        std::cerr << "    [warn] " << kTowerTree << " not in the input -- falling back to "
                  << kTowerTreeFallback << " (re-run the ntupler for the unsuppressed towers)\n";
        towTree = fin->Get(kTowerTreeFallback);
        towTreeUsed = kTowerTreeFallback;
    }
    TTree* wtaTree = fin->Get("wtaConeCellsTowersEtaSKTree");
    TTree* lrjTree = fin->Get("jetTaggerLRJEtaSKTree");
    TTree* bsmTree = fin->Get("truthBSMTree");
    TTree* evtTree = fin->Get("eventInfoTree");
    TTree* spTree  = fin->Get("showerParentTree");   // LLPs that seeded the showers
    TTree* jetTree = isLRJ ? lrjTree : wtaTree;
    // ---- diagnostics: tree presence + entries ----
    std::cout << "[fillShowerHists] " << file << (isLRJ ? "  [LRJ]" : "  [WTACone]")
              << (useTruth ? "  useTruth" : "  leadingJets") << "\n";
    for (auto pr : {std::make_pair(towTreeUsed,               towTree),
                    std::make_pair("jetTree(selected)",       jetTree),
                    std::make_pair("truthBSMTree",            bsmTree),
                    std::make_pair("showerParentTree",        spTree),
                    std::make_pair("eventInfoTree",           evtTree)}) {
        Long64_t ne = pr.second ? pr.second->GetEntries() : -1;
        std::cout << "    " << pr.first << " entries=" << ne
                  << (!pr.second ? "  [MISSING]" : (ne == 0 ? "  [EMPTY]" : "")) << "\n";
    }
    if (!towTree || !jetTree) { std::cerr << "missing trees in " << file << "\n"; fin->Close(); return 0; }

    // ---- per-event JZ bookkeeping (background) ----
    // Older ntuples predate these branches, so each is optional: without them the
    // sample is filled unweighted and unfiltered, exactly as before.
    std::vector<double>* eventWeights = nullptr;
    int  sampleJZSlice = -1;
    bool passHSTP      = true;
    const bool haveWeights = evtTree && evtTree->GetBranch("eventWeights")  != nullptr;
    const bool haveSlice   = evtTree && evtTree->GetBranch("sampleJZSlice") != nullptr;
    const bool haveHSTP    = evtTree && evtTree->GetBranch("passHSTP")      != nullptr;
    if (haveWeights) evtTree->SetBranchAddress("eventWeights",  &eventWeights);
    if (haveSlice)   evtTree->SetBranchAddress("sampleJZSlice", &sampleJZSlice);
    if (haveHSTP)    evtTree->SetBranchAddress("passHSTP",      &passHSTP);
    const bool useWeights = isDijet && kApplyJZWeights && haveWeights;
    const bool useHSTP    = isDijet && kApplyHSTPFilter && haveHSTP;
    if (isDijet && kApplyJZWeights && !haveWeights)
        std::cerr << "    [warn] no eventWeights branch -- the JZ slices are being mixed "
                     "unweighted (re-run caloShowerShapeNTupler.C with a --jz slice)\n";
    std::cout << "    weights=" << (useWeights ? "JZ cross-section" : "1 per event")
              << "  HSTP filter=" << (useHSTP ? "on" : "off") << "\n";

    // per-slice tallies, printed after the loop
    long   nJets_jz[kNJZSlices]  = {0};
    double sumW_jz[kNJZSlices]   = {0.0};
    long   nEvtHSTPRejected = 0;
    long   nNoPointingFit   = 0;   // jets with <2 lit layers -> no dca3D

    std::vector<double> *tow_Eta=nullptr,*tow_Phi=nullptr,*tow_Et_l[7]={nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
    towTree->SetBranchAddress("Eta",&tow_Eta); towTree->SetBranchAddress("Phi",&tow_Phi);
    for(int l=0;l<7;++l) towTree->SetBranchAddress(Form("Et_l%d",l),&tow_Et_l[l]);

    std::vector<double> *wta_pt=nullptr,*wta_eta=nullptr,*wta_phi=nullptr;
    std::vector<float>  *lrj_Et=nullptr,*lrj_eta=nullptr,*lrj_phi=nullptr;
    if (isLRJ){ lrjTree->SetBranchAddress("Et",&lrj_Et); lrjTree->SetBranchAddress("Eta",&lrj_eta); lrjTree->SetBranchAddress("Phi",&lrj_phi); }
    else      { wtaTree->SetBranchAddress("Pt",&wta_pt); wtaTree->SetBranchAddress("Eta",&wta_eta); wtaTree->SetBranchAddress("Phi",&wta_phi); }

    // ---- shower-parent LLPs: the truth labels (decay vertex + kinematics) ----
    std::vector<double> *sp_pt=nullptr,*sp_eta=nullptr,*sp_phi=nullptr,*sp_m=nullptr;
    std::vector<double> *sp_dvx=nullptr,*sp_dvy=nullptr,*sp_dvz=nullptr;
    std::vector<double> *sp_lxy=nullptr,*sp_r3d=nullptr;
    const bool haveSP = spTree && spTree->GetEntries() > 0;
    if (haveSP){
        spTree->SetBranchAddress("pt",           &sp_pt);
        spTree->SetBranchAddress("eta",          &sp_eta);
        spTree->SetBranchAddress("phi",          &sp_phi);
        spTree->SetBranchAddress("m",            &sp_m);
        spTree->SetBranchAddress("decayVtx_x",   &sp_dvx);
        spTree->SetBranchAddress("decayVtx_y",   &sp_dvy);
        spTree->SetBranchAddress("decayVtx_z",   &sp_dvz);
        spTree->SetBranchAddress("decayVtx_Lxy", &sp_lxy);
        spTree->SetBranchAddress("decayVtx_r3d", &sp_r3d);
    }

    std::vector<int>    *bsm_hasDecayVtx=nullptr;
    std::vector<double> *bsm_pt=nullptr,*bsm_eta=nullptr,*bsm_phi=nullptr;
    if (bsmTree && !isDijet && useTruth){
        bsmTree->SetBranchAddress("pt",&bsm_pt); bsmTree->SetBranchAddress("eta",&bsm_eta);
        bsmTree->SetBranchAddress("phi",&bsm_phi); bsmTree->SetBranchAddress("hasDecayVtx",&bsm_hasDecayVtx);
    }

    long nJets=0;
    Long64_t nEvt=towTree->GetEntries();
    if (kMaxEventsPerSample >= 0 && nEvt > kMaxEventsPerSample){
        std::cout << "    capping at " << kMaxEventsPerSample << " of " << nEvt
                  << " events (kMaxEventsPerSample); later JZ slices are dropped\n";
        nEvt = kMaxEventsPerSample;
    }
    for (Long64_t iEvt=0; iEvt<nEvt; ++iEvt){
        towTree->GetEntry(iEvt); jetTree->GetEntry(iEvt);
        if (bsmTree && !isDijet) bsmTree->GetEntry(iEvt);
        if (haveSP) spTree->GetEntry(iEvt);
        if (evtTree) evtTree->GetEntry(iEvt);
        if (useHSTP && !passHSTP) { ++nEvtHSTPRejected; continue; }
        const double w = (useWeights && eventWeights && !eventWeights->empty())
                       ? eventWeights->front() : 1.0;

        std::vector<std::pair<double,double>> jetEP; // eta,phi
        std::vector<double> jetPt;
        if (isLRJ){ if(lrj_Et) for(size_t j=0;j<lrj_Et->size();++j){ jetEP.push_back({(double)(*lrj_eta)[j],(double)(*lrj_phi)[j]}); jetPt.push_back((*lrj_Et)[j]); } }
        else      { if(wta_pt) for(size_t j=0;j<wta_pt->size();++j){ jetEP.push_back({(*wta_eta)[j],(*wta_phi)[j]}); jetPt.push_back((*wta_pt)[j]); } }
        if (jetEP.empty()) continue;

        std::vector<int> selJets;
        if (!isDijet && useTruth && bsm_pt){
            std::vector<int> boi;
            for(size_t b=0;b<bsm_pt->size();++b)
                if (bsm_hasDecayVtx && (*bsm_hasDecayVtx)[b]==1 && (*bsm_pt)[b]>ptMinBSM) boi.push_back((int)b);
            for(size_t j=0;j<jetEP.size();++j){
                for(int b:boi) if (dRp(jetEP[j].first,jetEP[j].second,(*bsm_eta)[b],(*bsm_phi)[b])<Rassoc){ selJets.push_back((int)j); break; }
            }
        } else {
            std::vector<int> ord(jetEP.size()); for(size_t j=0;j<jetEP.size();++j) ord[j]=(int)j;
            std::sort(ord.begin(),ord.end(),[&](int a,int b){return jetPt[a]>jetPt[b];});
            for(int k=0;k<(int)ord.size() && k<2;++k) selJets.push_back(ord[k]);
        }

        // Et-rank the selected jets so the dca3D histograms mean "leading" and
        // "subleading" for every sample: the truth-matched selection above keeps
        // jet-container order, the leading-jet branch is already Et-sorted.
        std::sort(selJets.begin(), selJets.end(),
                  [&](int a,int b){ return jetPt[a] > jetPt[b]; });

        for (size_t rank=0; rank<selJets.size(); ++rank){
            const int sj = selJets[rank];
            double je=jetEP[sj].first, jp=jetEP[sj].second;
            // Barrel-only option. Deliberately a SOFT cut: forward jets are kept out
            // of every physics distribution, but they still feed the |eta|
            // diagnostics on page 8 -- otherwise a barrel-only run could never show
            // the barrel/endcap step at |eta| = 1.5 that motivates the cut.
            const bool passEta = (kMaxJetAbsEta <= 0.0 || std::fabs(je) <= kMaxJetAbsEta);
            if (!passEta) ++ex.nBeyondEtaCut;
            double sum[7]={0,0,0,0,0,0,0}; double total=0; int ntow=0;
            for (size_t t=0;t<tow_Eta->size();++t){
                if (dRp((*tow_Eta)[t],(*tow_Phi)[t],je,jp) >= Rassoc) continue;
                bool any=false;
                for(int l=0;l<7;++l){ double et = tow_Et_l[l]?(*tow_Et_l[l])[t]:0.0; if(et>etMinTower){ sum[l]+=et; total+=et; any=true; } }
                if (any) ++ntow;
            }
            if (total<=0) continue;
            double depth=0, emf=0;
            for(int l=0;l<7;++l){ double f=sum[l]/total; if(passEta) hfrac[l]->Fill(f,w); depth+=l*f; if(l<=3) emf+=f; }
            if (passEta){ hdepth->Fill(depth,w); hEMfrac->Fill(emf,w); hntow->Fill(ntow,w); }

            // Shower-pointing dca3D: the same per-layer-centroid line fit the event
            // displays draw, histogrammed in mm for the two leading jets. Jets whose
            // shower lights up fewer than two layers have no line and are skipped.
            double dca3DFitMM = -1.0;
            if (rank < 2 && hdca3D[rank]){
                double cc[3], dd[3], dcaM = 0.0;
                int nLayersUsed = 0;
                if (jetShowerPointing(je, jp, Rassoc, etMinTower, tow_Eta, tow_Phi, tow_Et_l,
                                      cc, dd, dcaM, &nLayersUsed)){
                    dca3DFitMM = dcaM * 1000.0;
                    if (passEta) hdca3D[rank]->Fill(clampToAxis(hdca3D[rank], dca3DFitMM), w);
                    // Diagnostics: is the tail geometry (|eta|), sampling (n layers)
                    // or softness (Et)?
                    const double jetEtGeV = isLRJ ? jetPt[sj]/1000.0 : jetPt[sj];
                    if (ex.dcaVsAbsEta)   ex.dcaVsAbsEta->Fill(std::fabs(je), dca3DFitMM, w);
                    if (ex.dcaVsNLayers)  ex.dcaVsNLayers->Fill(nLayersUsed, dca3DFitMM, w);
                    if (ex.dcaVsJetEt)    ex.dcaVsJetEt->Fill(clampToAxis(ex.jetEtAll, jetEtGeV), dca3DFitMM, w);
                    if (ex.dcaVsAbsEta2D) ex.dcaVsAbsEta2D->Fill(std::fabs(je),
                                                                 clampToAxis(hdca3D[rank], dca3DFitMM), w);
                } else {
                    ++nNoPointingFit;
                }
            }

            // ---- truth: the shower-initiating LLP matched to this jet ----
            // Highest-pt shower parent within the association cone; its decay vertex
            // is the displacement the fit above is trying to recover.
            double dca3DTruthMM = -1.0;
            if (passEta && rank < 2 && haveSP && sp_pt){
                int best=-1; double bestPt=-1.0;
                for (size_t s=0;s<sp_pt->size();++s){
                    if (dRp(je,jp,(*sp_eta)[s],(*sp_phi)[s]) >= Rassoc) continue;
                    if ((*sp_pt)[s] > bestPt){ bestPt=(*sp_pt)[s]; best=(int)s; }
                }
                if (best>=0){
                    dca3DTruthMM = truthShowerDca3D((*sp_dvx)[best],(*sp_dvy)[best],(*sp_dvz)[best],
                                              (*sp_eta)[best],(*sp_phi)[best]);
                    ++ex.nTruthMatched;
                    if (ex.truthDca3D[rank]) ex.truthDca3D[rank]->Fill(clampToAxis(ex.truthDca3D[rank], dca3DTruthMM), w);
                    if (ex.truthLxy)  ex.truthLxy->Fill(clampToAxis(ex.truthLxy, (*sp_lxy)[best]), w);
                    if (ex.truthR3d)  ex.truthR3d->Fill(clampToAxis(ex.truthR3d, (*sp_r3d)[best]), w);
                    if (ex.truthPt)   ex.truthPt->Fill(clampToAxis(ex.truthPt, (*sp_pt)[best]), w);
                    if (ex.truthEta)  ex.truthEta->Fill((*sp_eta)[best], w);
                    if (ex.truthPhi)  ex.truthPhi->Fill((*sp_phi)[best], w);
                    if (ex.truthMass && sp_m) ex.truthMass->Fill(clampToAxis(ex.truthMass, (*sp_m)[best]), w);
                    // fit vs truth, leading jets only (best-measured shower)
                    if (rank==0 && dca3DFitMM >= 0.0){
                        if (ex.dca3DRes)  ex.dca3DRes->Fill(clampToAxis(ex.dca3DRes, dca3DFitMM - dca3DTruthMM), w);
                        if (ex.dca3DCorr) ex.dca3DCorr->Fill(dca3DTruthMM, dca3DFitMM, w);
                    }
                }
            }

            // ---- rate view: leading-jet Et, with and without the displaced cut ----
            // Background weights are per-second, so the running sum of these bins from
            // a threshold upwards is the trigger rate in Hz at that threshold.
            if (passEta && rank==0){
                // WTACone Pt is written in GeV by the ntupler; the LRJ Et branch is a
                // GEP pass-through in MeV (same convention as make_training_parquet.py).
                const double jetEt = isLRJ ? jetPt[sj]/1000.0 : jetPt[sj];
                if (ex.jetEtAll) ex.jetEtAll->Fill(clampToAxis(ex.jetEtAll, jetEt), w);
                if (ex.jetEtDca && dca3DFitMM >= kDcaCutMM)
                    ex.jetEtDca->Fill(clampToAxis(ex.jetEtDca, jetEt), w);
            }
            if (!passEta) continue;   // counted above; keep it out of the jet tallies
            ++nJets;
            if (sampleJZSlice>=0 && sampleJZSlice<kNJZSlices){
                ++nJets_jz[sampleJZSlice];
                sumW_jz[sampleJZSlice]+=w;
            }
        }
    }

    // ---- per-slice diagnostics: which JZ slices actually populate the plots ----
    if (haveSlice){
        long   nJetsSliced=0; double sumWTotal=0;
        for(int jz=0;jz<kNJZSlices;++jz){ nJetsSliced+=nJets_jz[jz]; sumWTotal+=sumW_jz[jz]; }
        if (nJetsSliced>0){
            std::cout << "    per-JZ-slice jets (sum of weights):";
            for(int jz=0;jz<kNJZSlices;++jz)
                if (nJets_jz[jz]>0)
                    std::cout << "  JZ" << jz << "=" << nJets_jz[jz] << " (" << sumW_jz[jz] << ")";
            std::cout << "\n    total sliced jets=" << nJetsSliced << "  sum of weights=" << sumWTotal
                      << (nJets>nJetsSliced ? "  [some jets had no slice index]" : "") << "\n";
        }
    }
    if (nEvtHSTPRejected)
        std::cout << "    HSTP filter rejected " << nEvtHSTPRejected << " of " << nEvt << " events\n";
    if (nNoPointingFit)
        std::cout << "    no shower-pointing fit (<2 lit layers) for " << nNoPointingFit
                  << " of the two leading jets\n";
    std::cout << "    truth shower parents matched to a leading/subleading jet: "
              << ex.nTruthMatched << (haveSP ? "" : "  [no showerParentTree]") << "\n";

    fin->Close();
    return nJets;
}

static void styleOverlay(TH1F* hs, TH1F* hd){
    if(hs){ hs->SetLineColor(kRed+1);  hs->SetLineWidth(2); if(hs->Integral()>0) hs->Scale(1.0/hs->Integral()); }
    if(hd){ hd->SetLineColor(kBlue+1); hd->SetLineWidth(2); hd->SetLineStyle(2); if(hd->Integral()>0) hd->Scale(1.0/hd->Integral()); }
}

// ===========================================================================
void caloShowerShapePlots(std::string signalFile = "",
                          std::string dijetFile = kDijetJZGlob,
                          std::string outDir    = ".",
                          double      ptMinBSM   = 5.0,
                          // Per-LAYER zero suppression on top of whatever the tower
                          // collection already applies. 0.5 GeV was emptying the cone:
                          // ~2500 towers/event survive EtaSK and ~50 tower slots fall
                          // inside DeltaR<0.4, yet only 2-3 had any layer above 0.5 GeV,
                          // because a 1-2 GeV tower spreads that Et over up to 7 layers.
                          // Starving the fit of layer centroids is what produced the
                          // huge-DCA3D tail, so the default is now 0 (keep everything).
                          double      etMinTower = 0.0,
                          bool        useTruth   = false) {
    // No-arg default (root -b -l -q 'caloShowerShapePlots.C'): run both signals
    // once each against the chained QCD JZ0-9 background, each written into its
    // own plots/<label>/ sub-directory. Pass a signalFile explicitly to process one
    // sample (and "" as dijetFile to drop the background overlay).
    if (signalFile.empty()) {
        const std::vector<std::pair<std::string,std::string>> samples = {
            {"/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_displaced_dark_photon.root", "displaced_dark_photon"},
            {"/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_emerging_jets.root",         "emerging_jets"},
        };
        for (const auto& s : samples) {
            std::string od = "plots/" + s.second;
            gSystem->mkdir(od.c_str(), kTRUE);
            std::cout << "[caloShowerShapePlots] === " << s.second << " -> " << od << "/ ===\n";
            caloShowerShapePlots(s.first, dijetFile, od, ptMinBSM, etMinTower, useTruth);
        }
        return;
    }

    gStyle->SetOptStat(0);

    struct JetColl { const char* tag; bool isLRJ; double Rassoc; };
    std::vector<JetColl> colls = { { "WTACone", false, 0.4 } };
    // The LRJ pass doubles the run time (a second full read of the JZ chain) for a
    // collection whose low slices are nearly empty — JZ1 contributed zero LRJ jets.
    // Flip kRunLRJ back on when the large-R side is what you are after.
    if (kRunLRJ) colls.push_back({ "LRJ", true, 1.0 });

    for (const auto& coll : colls) {
        TH1F* sfrac[7]; TH1F* dfrac[7];
        for(int l=0;l<7;++l){
            sfrac[l]=new TH1F(Form("s_frac%d_%s",l,coll.tag),Form("layer %d E_{T} fraction;E_{T,l%d}/E_{T};jets (norm)",l,l),40,0,1);
            dfrac[l]=new TH1F(Form("d_frac%d_%s",l,coll.tag),"",40,0,1);
            sfrac[l]->SetDirectory(nullptr); dfrac[l]->SetDirectory(nullptr);
        }
        TH1F* sdepth=new TH1F(Form("s_depth_%s",coll.tag),"shower depth;E_{T}-weighted mean layer;jets (norm)",42,0,6.3);
        TH1F* ddepth=new TH1F(Form("d_depth_%s",coll.tag),"",42,0,6.3);
        TH1F* sEM   =new TH1F(Form("s_emf_%s",coll.tag),"EM fraction;(l0+l1+l2+l3)/E_{T};jets (norm)",40,0,1);
        TH1F* dEM   =new TH1F(Form("d_emf_%s",coll.tag),"",40,0,1);
        TH1F* sntow =new TH1F(Form("s_ntow_%s",coll.tag),"n constituent towers;n towers;jets (norm)",40,0,40);
        TH1F* dntow =new TH1F(Form("d_ntow_%s",coll.tag),"",40,0,40);
        // Shower-pointing dca3D, leading (index 0) and subleading (index 1) jet.
        TH1F* sdca3D[2]; TH1F* ddca3D[2];
        const char* dcaRank[2] = { "leading", "subleading" };
        for(int r=0;r<2;++r){
            sdca3D[r]=new TH1F(Form("s_dca3D_%d_%s",r,coll.tag),
                            Form("%s jet shower-line DCA_{3D};DCA_{3D} [mm];jets (norm)",dcaRank[r]),
                            kNDcaBins,0,kDcaMaxMM);
            ddca3D[r]=new TH1F(Form("d_dca3D_%d_%s",r,coll.tag),"",kNDcaBins,0,kDcaMaxMM);
            sdca3D[r]->SetDirectory(nullptr); ddca3D[r]->SetDirectory(nullptr);
            sdca3D[r]->Sumw2(); ddca3D[r]->Sumw2();
        }
        for (TH1F* h : {sdepth,ddepth,sEM,dEM,sntow,dntow}) h->SetDirectory(nullptr);
        // The JZ weights span orders of magnitude, so the area normalization below
        // needs proper sum-of-squares errors.
        for(int l=0;l<7;++l){ sfrac[l]->Sumw2(); dfrac[l]->Sumw2(); }
        for (TH1F* h : {sdepth,ddepth,sEM,dEM,sntow,dntow}) h->Sumw2();

        // Truth + rate observables, booked identically for signal and background so
        // they can be overlaid (QCD simply has no truth shower parents).
        ExtraHists sx, dx;
        auto bookExtra = [&](ExtraHists& e, const char* pfx){
            e = ExtraHists();
            for(int r=0;r<2;++r){
                e.truthDca3D[r]=new TH1F(Form("%s_tdca3D_%d_%s",pfx,r,coll.tag),
                                      Form("%s jet truth DCA_{3D};truth DCA_{3D} [mm];jets (norm)",dcaRank[r]),
                                      kNDcaBins,0,kDcaMaxMM);
            }
            e.truthLxy =new TH1F(Form("%s_tlxy_%s",pfx,coll.tag),"truth decay L_{xy};L_{xy} [mm];jets (norm)",50,0,4000);
            e.truthR3d =new TH1F(Form("%s_tr3d_%s",pfx,coll.tag),"truth decay |r|;|r_{decay}| [mm];jets (norm)",50,0,6000);
            e.truthPt  =new TH1F(Form("%s_tpt_%s", pfx,coll.tag),"truth parent p_{T};p_{T} [GeV];jets (norm)",50,0,500);
            e.truthEta =new TH1F(Form("%s_teta_%s",pfx,coll.tag),"truth parent #eta;#eta;jets (norm)",50,-5,5);
            e.truthPhi =new TH1F(Form("%s_tphi_%s",pfx,coll.tag),"truth parent #phi;#phi;jets (norm)",50,-M_PI,M_PI);
            e.truthMass=new TH1F(Form("%s_tm_%s",  pfx,coll.tag),"truth parent mass;m [GeV];jets (norm)",50,0,200);
            e.dca3DRes    =new TH1F(Form("%s_dca3Dres_%s",pfx,coll.tag),
                                 "fitted - truth DCA_{3D} (leading);DCA_{3D}^{fit} - DCA_{3D}^{truth} [mm];jets (norm)",
                                 60,-kDcaMaxMM,kDcaMaxMM);
            e.dca3DCorr   =new TH2F(Form("%s_dca3Dcorr_%s",pfx,coll.tag),
                                 "fitted vs truth DCA_{3D} (leading);truth DCA_{3D} [mm];fitted DCA_{3D} [mm]",
                                 40,0,kDcaMaxMM,40,0,kDcaMaxMM);
            e.jetEtAll =new TH1F(Form("%s_etall_%s",pfx,coll.tag),
                                 "leading jet E_{T};E_{T} [GeV];rate [Hz] or jets",kNEtBins,0,kEtMaxGeV);
            e.jetEtDca  =new TH1F(Form("%s_etdca_%s", pfx,coll.tag),"",kNEtBins,0,kEtMaxGeV);
            e.dcaVsAbsEta =new TProfile(Form("%s_dcaeta_%s",pfx,coll.tag),
                                        "#LTDCA_{3D}#GT vs |#eta|;|#eta_{jet}|;#LTDCA_{3D}#GT [mm]",25,0,5);
            e.dcaVsNLayers=new TProfile(Form("%s_dcanlay_%s",pfx,coll.tag),
                                        "#LTDCA_{3D}#GT vs layers used;n layers in fit;#LTDCA_{3D}#GT [mm]",7,0.5,7.5);
            e.dcaVsJetEt  =new TProfile(Form("%s_dcaet_%s",pfx,coll.tag),
                                        "#LTDCA_{3D}#GT vs jet E_{T};jet E_{T} [GeV];#LTDCA_{3D}#GT [mm]",kNEtBins,0,kEtMaxGeV);
            e.dcaVsAbsEta2D=new TH2F(Form("%s_dcaeta2d_%s",pfx,coll.tag),
                                     "DCA_{3D} vs |#eta|;|#eta_{jet}|;DCA_{3D} [mm]",25,0,5,kNDcaBins,0,kDcaMaxMM);
            for (TProfile* pr : {e.dcaVsAbsEta,e.dcaVsNLayers,e.dcaVsJetEt}) pr->SetDirectory(nullptr);
            e.dcaVsAbsEta2D->SetDirectory(nullptr); e.dcaVsAbsEta2D->Sumw2();
            for (TH1F* h : {e.truthDca3D[0],e.truthDca3D[1],e.truthLxy,e.truthR3d,e.truthPt,e.truthEta,
                            e.truthPhi,e.truthMass,e.dca3DRes,e.jetEtAll,e.jetEtDca}){
                h->SetDirectory(nullptr); h->Sumw2();
            }
            e.dca3DCorr->SetDirectory(nullptr); e.dca3DCorr->Sumw2();
        };
        bookExtra(sx,"s"); bookExtra(dx,"d");

        long nS=0,nD=0;
        if (!signalFile.empty()) nS=fillShowerHists(signalFile,false,coll.isLRJ,coll.Rassoc,ptMinBSM,etMinTower,useTruth,sfrac,sdepth,sEM,sntow,sdca3D,sx);
        if (!dijetFile.empty())  nD=fillShowerHists(dijetFile, true, coll.isLRJ,coll.Rassoc,ptMinBSM,etMinTower,useTruth,dfrac,ddepth,dEM,dntow,ddca3D,dx);
        std::cout << "[caloShowerShapePlots] " << coll.tag << ": signal jets=" << nS << "  dijet jets=" << nD << "\n";

        for(int l=0;l<7;++l) styleOverlay(sfrac[l],dfrac[l]);
        styleOverlay(sdepth,ddepth); styleOverlay(sEM,dEM); styleOverlay(sntow,dntow);
        for(int r=0;r<2;++r) styleOverlay(sdca3D[r],ddca3D[r]);

        TString pdf = TString(outDir) + "/caloShowerShapePlots_" + coll.tag + ".pdf";
        TCanvas* c = new TCanvas("cShower","shower shapes",1500,1000);
        c->Print(pdf + "[");

        auto drawPair=[&](TH1F* hs, TH1F* hd, const char* title){
            double mx=0; if(hs) mx=std::max(mx,hs->GetMaximum()); if(hd) mx=std::max(mx,hd->GetMaximum());
            if(hs){ hs->SetTitle(title); hs->SetMaximum(1.3*mx); hs->Draw("hist"); if(hd) hd->Draw("hist same"); }
            else if(hd){ hd->SetTitle(title); hd->SetMaximum(1.3*mx); hd->Draw("hist"); }
            // The counts are the number of jets ENTERING each histogram (not the
            // weighted yield), so say so — a bare number next to "dijet JZ-wgt" reads
            // like a rate otherwise.
            TLegend* lg=new TLegend(0.48,0.70,0.93,0.90);
            lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045);
            if(hs) lg->AddEntry(hs,Form("signal: %ld jets",nS),"l");
            if(hd) lg->AddEntry(hd,Form("dijet%s: %ld jets", kApplyJZWeights ? " (JZ-weighted)" : "", nD),"l");
            lg->Draw();
        };

        // page 1: the 7 layer fractions (4x2)
        c->Clear(); c->Divide(4,2);
        for(int l=0;l<7;++l){ c->cd(l+1); drawPair(sfrac[l],dfrac[l],Form("layer %d E_{T} fraction;E_{T,l%d}/E_{T};jets (norm)",l,l)); }
        c->cd(8); { TLatex tl; tl.SetNDC(); tl.SetTextSize(0.06);
            tl.DrawLatex(0.1,0.6,Form("%s jets",coll.tag)); tl.DrawLatex(0.1,0.45,"red = signal (solid)"); tl.DrawLatex(0.1,0.33,"blue = dijet (dashed)");
            if (kApplyJZWeights) { tl.SetTextSize(0.045); tl.DrawLatex(0.1,0.20,"dijet: JZ0-9 cross-section weighted"); } }
        c->Print(pdf);

        // page 2: summary observables (2x2)
        c->Clear(); c->Divide(2,2);
        c->cd(1); drawPair(sdepth,ddepth,"shower depth;E_{T}-weighted mean layer;jets (norm)");
        c->cd(2); drawPair(sEM,dEM,"EM fraction;(l0+l1+l2+l3)/E_{T};jets (norm)");
        c->cd(3); drawPair(sntow,dntow,"n constituent towers;n towers;jets (norm)");
        c->Print(pdf);

        // page 3: shower-line DCA3D of the two leading jets, one pad each. This is
        // the displacement-sensitive observable: the fitted per-layer-centroid line
        // is not forced through the IP, so a shower that started at a displaced
        // vertex misses the origin. Last bin holds the overflow.
        c->Clear(); c->Divide(2,1);
        for(int r=0;r<2;++r){
            c->cd(r+1);
            drawPair(sdca3D[r],ddca3D[r],
                     Form("%s jet shower-line DCA_{3D};DCA_{3D} [mm];jets (norm)",dcaRank[r]));
            TLatex tl; tl.SetNDC(); tl.SetTextSize(0.028);
            tl.DrawLatex(0.14,0.86,Form("%s jets, #DeltaR < %.1f",coll.tag,coll.Rassoc));
            tl.DrawLatex(0.14,0.82,"last bin = overflow");
        }
        c->Print(pdf);

        // page 4: truth displacement of the matched shower parent (signal only —
        // prompt QCD has no shower parents, so the dijet curves are empty here).
        for(int r=0;r<2;++r) styleOverlay(sx.truthDca3D[r],dx.truthDca3D[r]);
        styleOverlay(sx.truthLxy,dx.truthLxy); styleOverlay(sx.truthR3d,dx.truthR3d);
        styleOverlay(sx.dca3DRes,dx.dca3DRes);
        c->Clear(); c->Divide(2,2);
        c->cd(1); drawPair(sx.truthDca3D[0],dx.truthDca3D[0],"leading jet truth DCA_{3D};truth DCA_{3D} [mm];jets (norm)");
        c->cd(2); drawPair(sx.truthDca3D[1],dx.truthDca3D[1],"subleading jet truth DCA_{3D};truth DCA_{3D} [mm];jets (norm)");
        c->cd(3); drawPair(sx.truthLxy,dx.truthLxy,"truth decay L_{xy};L_{xy} [mm];jets (norm)");
        c->cd(4); drawPair(sx.truthR3d,dx.truthR3d,"truth decay |r|;|r_{decay}| [mm];jets (norm)");
        c->Print(pdf);

        // page 5: truth kinematics of the matched shower parent
        styleOverlay(sx.truthPt,dx.truthPt);   styleOverlay(sx.truthEta,dx.truthEta);
        styleOverlay(sx.truthPhi,dx.truthPhi); styleOverlay(sx.truthMass,dx.truthMass);
        c->Clear(); c->Divide(2,2);
        c->cd(1); drawPair(sx.truthPt,dx.truthPt,"truth parent p_{T};p_{T} [GeV];jets (norm)");
        c->cd(2); drawPair(sx.truthEta,dx.truthEta,"truth parent #eta;#eta;jets (norm)");
        c->cd(3); drawPair(sx.truthPhi,dx.truthPhi,"truth parent #phi;#phi;jets (norm)");
        c->cd(4); drawPair(sx.truthMass,dx.truthMass,"truth parent mass;m [GeV];jets (norm)");
        c->Print(pdf);

        // page 6: how well the fit recovers the truth displacement
        c->Clear(); c->Divide(2,2);
        c->cd(1); drawPair(sx.dca3DRes,dx.dca3DRes,"fitted - truth DCA_{3D} (leading);DCA_{3D}^{fit} - DCA_{3D}^{truth} [mm];jets (norm)");
        c->cd(2); if (sx.dca3DCorr->GetEntries()>0) sx.dca3DCorr->Draw("colz");
        c->cd(3); { TLatex tl; tl.SetNDC(); tl.SetTextSize(0.05);
            tl.DrawLatex(0.05,0.80,"fitted vs truth DCA_{3D}: signal, leading jets");
            tl.SetTextSize(0.038);
            tl.DrawLatex(0.05,0.68,Form("signal jets with a truth parent: %ld", sx.nTruthMatched));
            tl.DrawLatex(0.05,0.60,Form("dijet jets with a truth parent: %ld  (expect 0)", dx.nTruthMatched));
            tl.DrawLatex(0.05,0.48,"truth DCA_{3D} = |v - (v#upoint#hat{u})#hat{u}| from the");
            tl.DrawLatex(0.05,0.41,"LLP decay vertex + its direction;");
            tl.DrawLatex(0.05,0.34,"the fit uses only calorimeter towers."); }
        c->Print(pdf);

        // page 7: the point of the study — can a displacement requirement buy back
        // an E_T threshold? Left: background rate (Hz) above threshold, with and
        // without the dca3D requirement (background weights are per-second, so the
        // running sum of bins IS the rate). Right: the signal efficiency paid for it.
        {
            std::vector<double> thr, rateAll, rateDca, effAll, effDca;
            const double sTotal = sx.jetEtAll->Integral(1, kNEtBins);
            for (int b=1; b<=kNEtBins; ++b){
                thr.push_back(dx.jetEtAll->GetXaxis()->GetBinLowEdge(b));
                rateAll.push_back(dx.jetEtAll->Integral(b, kNEtBins));
                rateDca .push_back(dx.jetEtDca ->Integral(b, kNEtBins));
                effAll .push_back(sTotal>0 ? sx.jetEtAll->Integral(b,kNEtBins)/sTotal : 0.0);
                effDca  .push_back(sTotal>0 ? sx.jetEtDca ->Integral(b,kNEtBins)/sTotal : 0.0);
            }
            auto mkGraph=[&](std::vector<double>& y,int col,int style){
                TGraph* g=new TGraph((int)thr.size(), thr.data(), y.data());
                g->SetLineColor(col); g->SetLineWidth(2); g->SetLineStyle(style);
                g->SetMarkerColor(col); g->SetMarkerStyle(20); g->SetMarkerSize(0.6);
                return g;
            };
            c->Clear(); c->Divide(2,1);
            c->cd(1); gPad->SetLogy();
            TGraph* gr1=mkGraph(rateAll,kBlue+1,1); TGraph* gr2=mkGraph(rateDca,kRed+1,2);
            gr1->SetTitle("QCD rate above threshold;leading jet E_{T} threshold [GeV];rate [Hz]");
            gr1->Draw("ALP"); gr2->Draw("LP same");
            { TLegend* lg=new TLegend(0.45,0.72,0.9,0.9); lg->SetBorderSize(0);
              lg->AddEntry(gr1,"all jets","lp");
              lg->AddEntry(gr2,Form("DCA_{3D}^{fit} > %.0f mm",kDcaCutMM),"lp"); lg->Draw(); }
            c->cd(2);
            TGraph* ge1=mkGraph(effAll,kBlue+1,1); TGraph* ge2=mkGraph(effDca,kRed+1,2);
            ge1->SetTitle("signal efficiency above threshold;leading jet E_{T} threshold [GeV];fraction of signal jets");
            ge1->SetMinimum(0.0); ge1->SetMaximum(1.1);
            ge1->Draw("ALP"); ge2->Draw("LP same");
            { TLegend* lg=new TLegend(0.45,0.72,0.9,0.9); lg->SetBorderSize(0);
              lg->AddEntry(ge1,"all jets","lp");
              lg->AddEntry(ge2,Form("DCA_{3D}^{fit} > %.0f mm",kDcaCutMM),"lp"); lg->Draw(); }
            c->Print(pdf);
            std::cout << "[caloShowerShapePlots] " << coll.tag
                      << ": QCD rate at 0 GeV = " << (rateAll.empty()?0:rateAll[0]) << " Hz, with dca3D cut = "
                      << (rateDca.empty()?0:rateDca[0]) << " Hz\n";
        }

        // page 8: what is the DCA3D tail correlated with? The miss distance is
        // amplified by ~2.4 m per radian of per-layer angular drift, so one 0.1
        // tower cell of drift is already ~240 mm. These three profiles separate the
        // candidate causes: |eta| (nominal-geometry effects: the barrel/endcap
        // switch at 1.5 and the endcap r = z/sinh(eta) sensitivity), the number of
        // layers the fit had to work with, and jet softness.
        {
            auto drawProfiles=[&](TProfile* ps, TProfile* pd, const char* title){
                double mx=0;
                if(ps) mx=std::max(mx,ps->GetMaximum());
                if(pd) mx=std::max(mx,pd->GetMaximum());
                if(ps){ ps->SetLineColor(kRed+1);  ps->SetLineWidth(2); ps->SetMarkerColor(kRed+1);
                        ps->SetMarkerStyle(20); ps->SetMarkerSize(0.7); }
                if(pd){ pd->SetLineColor(kBlue+1); pd->SetLineWidth(2); pd->SetLineStyle(2);
                        pd->SetMarkerColor(kBlue+1); pd->SetMarkerStyle(24); pd->SetMarkerSize(0.7); }
                if(ps){ ps->SetTitle(title); ps->SetMinimum(0); ps->SetMaximum(1.3*mx);
                        ps->Draw("E1"); if(pd) pd->Draw("E1 same"); }
                else if(pd){ pd->SetTitle(title); pd->SetMinimum(0); pd->SetMaximum(1.3*mx); pd->Draw("E1"); }
                TLegend* lg=new TLegend(0.48,0.74,0.93,0.90);
                lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.042);
                if(ps) lg->AddEntry(ps,"signal","lp");
                if(pd) lg->AddEntry(pd,kApplyJZWeights?"dijet (JZ-weighted)":"dijet","lp");
                lg->Draw();
            };
            c->Clear(); c->Divide(2,2);
            c->cd(1); drawProfiles(sx.dcaVsAbsEta, dx.dcaVsAbsEta,
                                   "#LTDCA_{3D}#GT vs |#eta|;|#eta_{jet}|;#LTDCA_{3D}#GT [mm]");
            c->cd(2); drawProfiles(sx.dcaVsNLayers, dx.dcaVsNLayers,
                                   "#LTDCA_{3D}#GT vs layers used;n layers in fit;#LTDCA_{3D}#GT [mm]");
            c->cd(3); drawProfiles(sx.dcaVsJetEt, dx.dcaVsJetEt,
                                   "#LTDCA_{3D}#GT vs jet E_{T};jet E_{T} [GeV];#LTDCA_{3D}#GT [mm]");
            c->cd(4); if (dx.dcaVsAbsEta2D->GetEntries()>0) dx.dcaVsAbsEta2D->Draw("colz");
                      else if (sx.dcaVsAbsEta2D->GetEntries()>0) sx.dcaVsAbsEta2D->Draw("colz");
            c->Print(pdf);
            std::cout << "[caloShowerShapePlots] " << coll.tag << ": jets dropped by |eta| > "
                      << kMaxJetAbsEta << ": signal " << sx.nBeyondEtaCut
                      << ", dijet " << dx.nBeyondEtaCut
                      << (kMaxJetAbsEta > 0.0 ? "" : "  (cut disabled)") << "\n";
        }

        c->Print(pdf + "]");
        std::cout << "[caloShowerShapePlots] wrote -> " << pdf << "\n";
        delete c;
    }
}
