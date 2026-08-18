// caloShowerEventDisplays.C
// ---------------------------------------------------------------------------
// Per-event calorimeter event displays for the displaced-jet / shower-shape
// study. Reads the caloShowerShapeNTupler.C output and, for a number of events,
// draws the per-layer calorimeter energy deposits of the jets that are ΔR-matched
// to displaced truth-BSM particles (signal), or of the leading jets (dijet
// background). Three views per event: r-z, transverse (x-y), and 3D.
//
//   # both signals + the QCD dijet background, into plots/<label>/ (no args):
//   root -b -l -q 'caloShowerEventDisplays.C'
//   # or one sample explicitly:
//   root -b -q 'caloShowerEventDisplays.C("in.root","displaced_dark_photon",false,20,"plots/")'
//   root -b -q 'caloShowerEventDisplays.C("caloShowerShape_dijet_JZ4.root","dijet_JZ4",true,20,"plots/")'
//   # the whole QCD background, chained over the ten JZ slices:
//   root -b -q 'caloShowerEventDisplays.C("caloShowerShape_dijet_JZ[0-9].root","dijet_JZ0to9",true,20,"plots/")'
//
// The input path may be a glob: ChainSource (../chainSource.h) reads every tree
// as a TChain spanning the matched files, so the ten per-slice dijet ntuples need
// no hadd. For background the drawn events are then taken as a fixed quota per JZ
// slice (kEventsPerSlice from each of kDisplaySlices, JZ1-4 by default, JZ0
// excluded), found by jumping straight to each slice's entry range via the chain's
// tree offsets. Each page is labeled with the event's JZ slice.
// Note the [0-9] before .root: the per-job outputs (..._JZ9_000510.root) sit in
// the same directory after the hadd, so a bare JZ*.root would chain both the
// merged files and their inputs.
//
// Two PDFs are written per call (one per jet collection):
//   <outDir>/caloShowerEventDisplays_WTACone_<sampleLabel>.pdf
//   <outDir>/caloShowerEventDisplays_LRJ_<sampleLabel>.pdf
//
// IMPORTANT physics note -----------------------------------------------------
// GEP towers are *projective*: every layer of a tower shares one IP-pointing
// (eta,phi). We therefore place each layer's Et at the *nominal* ATLAS radius
// (barrel) or z (endcap) for that sampling -- an illustrative geometry, not a
// measured position. The calo deposits will always appear to point back to the
// origin; the displacement is shown by the TRUTH overlay (production/decay
// vertices + flight path drawn off-origin) and by the per-layer longitudinal
// profile (color = layer). For the quantitative separation see
// caloShowerShapePlots.C.
//
// Dijet background pages carry no truth overlay -- those jets are prompt by
// construction, so the fitted line SHOULD pass through the IP and whatever DCA3D
// it returns is the resolution floor the signal has to beat. Instead each page
// states the JZ slice and the event's rate contribution in Hz (its JZ weight,
// normalized to 1 s of HL-LHC luminosity): JZ0/JZ1 events are worth O(kHz) each
// while JZ9 events are worth O(1e-4) Hz, so a page is not "one event's worth" of
// background. HSTP-failed events carry no rate and are skipped by default.
//
// Truth overlay (magenta, signal only): for each drawn jet, the highest-pt LLP
// from showerParentTree inside the cone gives the decay-vertex marker, the IP ->
// decay-vertex flight path. The truth Lxy / |r| / DCA3D / parent pT are printed on
// the info pad next to the fitted DCA3D.
//
// NOTE on truth DCA3D: because the LLP starts at the IP, the line along its own
// direction through its decay vertex passes back through the IP, so this truth DCA3D
// is ~0 by construction (the few tens of mm seen is the production-vertex offset).
// It is therefore NOT the quantity the fit should be validated against -- the decay
// RADIUS (Lxy, |r|) is what is physically displaced. See the README.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMarker.h"
#include "TLine.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TObject.h"
#include "TSystem.h"
#include "TPolyLine.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include "TEllipse.h"
#include "../chainSource.h"
#include "caloShowerPointing.h"

// ---------------------------------------------------------------------------
// The nominal layer geometry (kRbarrel / kZendcap / kEtaBarrelLimit), layerXYZ()
// and the shower-pointing line fit live in caloShowerPointing.h, shared with
// caloShowerShapePlots.C so the histogrammed dca3D is the same quantity this macro
// draws. Only the drawing-specific styling stays here.
// ---------------------------------------------------------------------------
// QCD dijet background for the no-arg default: the ten per-slice ntuples chained.
// [0-9].root excludes the per-job outputs that share the directory (see above).
static const char* kDijetJZGlobDisplays =
    "/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_dijet_JZ[0-9].root";

// Background events failing the HSTP filter do not contribute to the rate, so by
// default they are not drawn either — same convention as caloShowerShapePlots.C.
static const bool kSkipHSTPFailedBackground = true;

// Background event displays: how many events to draw from each JZ slice, and which
// slices. JZ0 is excluded on purpose -- the HSTP filter removes essentially all of
// it, so its events carry no rate and a page spent on one says nothing about the
// background we actually trigger on. JZ5-9 are omitted as well: they are the hard
// tail, tiny in rate, and the interesting failure modes live in the soft slices.
static const std::vector<int> kDisplaySlices  = { 1, 2, 3, 4 };
static const int              kEventsPerSlice = 10;

// Tower collection feeding the deposits and the pointing fit:
//   "gepCellsTowersTree"       - no soft killer (DEFAULT; the fit needs the towers)
//   "gepCellsTowersSKTree"     - plain SK
//   "gepCellsTowersEtaSKTree"  - eta-dependent SK (the trigger jets' own towers)
// EtaSK leaves a QCD jet with 2-3 towers over 1-2 layers, which is not enough for a
// per-layer-centroid fit. Falls back to the EtaSK tree, with a warning, when the
// requested tree is absent -- ntuples made before the extra collections were added
// only contain the EtaSK one.
static const char* kTowerTree         = "gepCellsTowersTree";
static const char* kTowerTreeFallback = "gepCellsTowersEtaSKTree";

// Draw the LRJ collection as well as WTACone. Off: the JetTaggerLRJ collection is
// empty in the JZ1 (v22 PU200) production, and the large-R side is not what this
// study measures. Mirrors kRunLRJ in caloShowerShapePlots.C.
static const bool kRunLRJ = false;

// Depth gradient EM(inner, blue) -> HAD(outer, red) so color encodes shower depth.
static       int    kLayerColor[7] = { 616, 600, 860, 433, 418, 807, 632 };
//                                     kMag  kBlu kAz  kCy+ kGr+ kOr+ kRed
static const char*  kLayerName[7]  = { "l0 PreSamp", "l1 EM1", "l2 EM2", "l3 EM3",
                                       "l4 Tile0", "l5 Tile1", "l6 Tile2" };

static inline double dR(double e1, double p1, double e2, double p2) {
    double dp = std::fabs(p1 - p2);
    while (dp > M_PI) dp = std::fabs(dp - 2*M_PI);
    double de = e1 - e2;
    return std::sqrt(de*de + dp*dp);
}

static inline double signedR(double x, double y) {
    return std::sqrt(x*x + y*y) * (y >= 0 ? 1.0 : -1.0);
}
static inline double etMarkerSize(double et) {
    double s = 0.3 + 0.28 * std::sqrt(std::max(0.0, et));
    return std::min(3.0, s);
}

// Draw the nominal calorimeter layer geometry (black, dotted) as a background
// reference into the 3 pads: barrel radii + endcap faces in r-z, concentric rings
// in x-y, and inner/outer barrel wireframe rings in 3D. Call after the frames and
// before the deposits so the hits render on top.
static void drawCaloGeometry(TCanvas* c, std::vector<TObject*>& garbage) {
    const int col = kBlack, sty = 3;   // dotted

    // --- pad 1: r-z --- barrel = horizontal lines at +/-R; endcap = vertical faces at +/-Z
    c->cd(1);
    for (int l=0;l<7;++l) {
        double R = kRbarrel[l], zb = R * std::sinh(kEtaBarrelLimit);
        for (double s : {+1.0,-1.0}) {
            TLine* h = new TLine(-zb, s*R, zb, s*R);
            h->SetLineColor(col); h->SetLineStyle(sty); h->Draw(); garbage.push_back(h);
        }
        double Z = kZendcap[l], rEnd = Z / std::sinh(kEtaBarrelLimit);
        for (double s : {+1.0,-1.0}) {
            TLine* v = new TLine(s*Z, -rEnd, s*Z, rEnd);
            v->SetLineColor(col); v->SetLineStyle(sty); v->Draw(); garbage.push_back(v);
        }
    }

    // --- pad 2: x-y --- concentric barrel rings
    c->cd(2);
    for (int l=0;l<7;++l) {
        TEllipse* e = new TEllipse(0,0,kRbarrel[l],kRbarrel[l]);
        e->SetFillStyle(0); e->SetLineColor(col); e->SetLineStyle(sty); e->Draw(); garbage.push_back(e);
    }

    // --- pad 3: 3D --- inner + outer barrel envelope as rings at z = -zb, 0, +zb
    c->cd(3);
    const int NC=40;
    for (double R : {kRbarrel[0], kRbarrel[6]}) {
        double zb = R * std::sinh(kEtaBarrelLimit);
        for (double zc : {-zb, 0.0, zb}) {
            TPolyLine3D* ring = new TPolyLine3D(NC+1);
            for (int k=0;k<=NC;++k){ double a=2*M_PI*k/NC; ring->SetPoint(k, R*std::cos(a), R*std::sin(a), zc); }
            ring->SetLineColor(col); ring->SetLineStyle(sty); ring->Draw(); garbage.push_back(ring);
        }
    }
}

// Draw the jet association cone (DeltaR = Rassoc in eta-phi). In 3D it is the full
// boundary ring at the outer calo layer plus a few spokes back to the IP. In the 2D
// panels we draw ONLY the outermost silhouette as a single wedge from the IP -- the
// eta extent (je +/- Rassoc, at fixed phi) in r-z and the phi extent (jp +/- Rassoc,
// at fixed eta) in x-y. Drawing the full projected ring there made one cone look like
// several: the signed-r fold mirrors it about r=0, and the barrel/endcap seam makes
// the ring stride across the pad. Uses the same projective placement as the deposits.
static void drawJetCone(TCanvas* c, double je, double jp, double Rassoc, int col,
                        std::vector<TObject*>& garbage) {
    const int N = 60;

    // --- 3D: full cone = boundary ring at the outer calo layer + spokes to the IP ---
    std::vector<double> bx(N+1), by(N+1), bz(N+1);
    for (int k=0;k<=N;++k) {
        double t  = 2*M_PI*k/N;
        double eb = je + Rassoc*std::cos(t);
        double pb = jp + Rassoc*std::sin(t);
        layerXYZ(6, eb, pb, bx[k], by[k], bz[k]);
    }
    c->cd(3);
    { TPolyLine3D* r=new TPolyLine3D(N+1); for(int k=0;k<=N;++k) r->SetPoint(k,bx[k],by[k],bz[k]);
      r->SetLineColor(col); r->SetLineWidth(2); r->Draw(); garbage.push_back(r); }
    for (int k=0;k<N;k+=10) {
        TPolyLine3D* g=new TPolyLine3D(2); g->SetPoint(0,0,0,0); g->SetPoint(1,bx[k],by[k],bz[k]);
        g->SetLineColor(col); g->Draw(); garbage.push_back(g);
    }

    // --- 2D: outermost silhouette only, a wedge (two edges + closing base) from IP ---
    double ax,ay,az, bx2,by2,bz2;
    // r-z: eta extent at fixed phi = jp (single phi -> no signed-r fold)
    c->cd(1);
    layerXYZ(6, je-Rassoc, jp, ax,ay,az);
    layerXYZ(6, je+Rassoc, jp, bx2,by2,bz2);
    { TLine* e1=new TLine(0,0,az,signedR(ax,ay));         e1->SetLineColor(col); e1->SetLineWidth(2); e1->Draw(); garbage.push_back(e1);
      TLine* e2=new TLine(0,0,bz2,signedR(bx2,by2));      e2->SetLineColor(col); e2->SetLineWidth(2); e2->Draw(); garbage.push_back(e2);
      TLine* bs=new TLine(az,signedR(ax,ay),bz2,signedR(bx2,by2)); bs->SetLineColor(col); bs->SetLineWidth(2); bs->Draw(); garbage.push_back(bs); }
    // x-y: phi extent at fixed eta = je
    c->cd(2);
    layerXYZ(6, je, jp-Rassoc, ax,ay,az);
    layerXYZ(6, je, jp+Rassoc, bx2,by2,bz2);
    { TLine* e1=new TLine(0,0,ax,ay);       e1->SetLineColor(col); e1->SetLineWidth(2); e1->Draw(); garbage.push_back(e1);
      TLine* e2=new TLine(0,0,bx2,by2);     e2->SetLineColor(col); e2->SetLineWidth(2); e2->Draw(); garbage.push_back(e2);
      TLine* bs=new TLine(ax,ay,bx2,by2);   bs->SetLineColor(col); bs->SetLineWidth(2); bs->Draw(); garbage.push_back(bs); }
}

// ===========================================================================
void caloShowerEventDisplays(std::string inputFile = "",
                             std::string sampleLabel = "signal",
                             bool        isDijet     = false,
                             int         nEventsToDraw = 20,
                             std::string outDir      = ".",
                             double      ptMinBSM    = 5.0,   // GeV, BSM select
                             // Per-LAYER threshold, matching caloShowerShapePlots.C: 0
                             // keeps every layer. At 0.5 GeV a QCD jet was left with
                             // 1-2 lit layers and the pointing fit could not run.
                             double      etMinTower  = 0.0,
                             bool        useTruth    = false)  // false = leading jets, skip truth overlay
{
    // No-arg default (root -b -l -q 'caloShowerEventDisplays.C'): both signals plus
    // the QCD dijet background, each written into its own plots/<label>/ dir.
    // Pass an inputFile explicitly to process a single sample.
    if (inputFile.empty()) {
        struct Sample { const char* path; const char* label; bool isDijet; };
        const std::vector<Sample> samples = {
            {"/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_displaced_dark_photon.root", "displaced_dark_photon", false},
            {"/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_emerging_jets.root",         "emerging_jets",         false},
            // QCD dijet: the ten slices chained, so the drawn events span JZ0-9
            // (see the [0-9] note above). No truth to overlay -- these jets are
            // prompt by construction, which is exactly the point of looking at
            // them: the fitted line should pass through the IP, and whatever
            // DCA3D they show is the resolution floor the signal must beat.
            {kDijetJZGlobDisplays, "dijet_JZ0to9", true},
        };
        // Truth-BSM matching selects no jets for these samples right now, so the
        // no-arg default draws leading jets (jets + towers only, no truth overlay).
        // Flip to true to require displaced-BSM matching once truth is sorted out.
        const bool useTruthDefault = false;
        for (const auto& s : samples) {
            std::string od = std::string("plots/") + s.label;
            gSystem->mkdir(od.c_str(), kTRUE);
            std::cout << "[caloShowerEventDisplays] === " << s.label << " -> " << od
                      << "/  (isDijet=" << s.isDijet << ", useTruth=" << useTruthDefault << ") ===\n";
            caloShowerEventDisplays(s.path, s.label, s.isDijet, nEventsToDraw, od, ptMinBSM, etMinTower, useTruthDefault);
        }
        return;
    }

    gStyle->SetOptStat(0);

    ChainSource* fin = ChainSource::Open(inputFile.c_str());
    if (!fin || fin->IsZombie()) { std::cerr << "cannot open " << inputFile << "\n"; return; }

    TTree* towTree = fin->Get(kTowerTree);
    if (!towTree) {
        std::cerr << "[caloShowerEventDisplays] " << kTowerTree << " not in the input -- "
                     "falling back to " << kTowerTreeFallback << " (re-run the ntupler to get "
                     "the unsuppressed towers)\n";
        towTree = fin->Get(kTowerTreeFallback);
    } else {
        std::cout << "[caloShowerEventDisplays] tower collection: " << kTowerTree << "\n";
    }
    TTree* wtaTree = fin->Get("wtaConeCellsTowersEtaSKTree");
    TTree* lrjTree = fin->Get("jetTaggerLRJEtaSKTree");
    TTree* bsmTree = fin->Get("truthBSMTree");
    TTree* evtTree = fin->Get("eventInfoTree");
    TTree* spTree  = fin->Get("showerParentTree");   // LLPs that seeded the showers
    const size_t nInputFiles = fin->Files().size();
    // ---- diagnostics: tree presence + entries ----
    std::cout << "[caloShowerEventDisplays] input: " << inputFile
              << "  (" << nInputFiles << " file(s))"
              << (useTruth ? "  useTruth" : "  leadingJets") << "\n";
    for (auto pr : {std::make_pair("towerTree(selected)",         towTree),
                    std::make_pair("wtaConeCellsTowersEtaSKTree", wtaTree),
                    std::make_pair("jetTaggerLRJEtaSKTree",       lrjTree),
                    std::make_pair("truthBSMTree",                bsmTree),
                    std::make_pair("showerParentTree",            spTree),
                    std::make_pair("eventInfoTree",               evtTree)}) {
        Long64_t ne = pr.second ? pr.second->GetEntries() : -1;
        std::cout << "    " << pr.first << " entries=" << ne
                  << (!pr.second ? "  [MISSING]" : (ne == 0 ? "  [EMPTY]" : "")) << "\n";
    }
    if (!towTree) { std::cerr << "no tower tree in " << inputFile << "\n"; return; }

    // ---- tower branches ----
    std::vector<double> *tow_Et=nullptr, *tow_Eta=nullptr, *tow_Phi=nullptr;
    std::vector<double> *tow_Et_l[7] = {nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr};
    towTree->SetBranchAddress("Et",  &tow_Et);
    towTree->SetBranchAddress("Eta", &tow_Eta);
    towTree->SetBranchAddress("Phi", &tow_Phi);
    for (int l=0;l<7;++l) towTree->SetBranchAddress(Form("Et_l%d",l), &tow_Et_l[l]);

    // ---- WTACone jets (double) ----
    std::vector<double> *wta_pt=nullptr, *wta_eta=nullptr, *wta_phi=nullptr;
    if (wtaTree) {
        wtaTree->SetBranchAddress("Pt",  &wta_pt);
        wtaTree->SetBranchAddress("Eta", &wta_eta);
        wtaTree->SetBranchAddress("Phi", &wta_phi);
    }
    // ---- LRJ jets (float) ----
    std::vector<float> *lrj_Et=nullptr, *lrj_eta=nullptr, *lrj_phi=nullptr;
    if (lrjTree) {
        lrjTree->SetBranchAddress("Et",  &lrj_Et);
        lrjTree->SetBranchAddress("Eta", &lrj_eta);
        lrjTree->SetBranchAddress("Phi", &lrj_phi);
    }
    // ---- truth BSM ----
    std::vector<int>    *bsm_pdgId=nullptr, *bsm_hasDecayVtx=nullptr, *bsm_hasProdVtx=nullptr;
    std::vector<double> *bsm_pt=nullptr, *bsm_eta=nullptr, *bsm_phi=nullptr;
    std::vector<double> *bsm_pvx=nullptr,*bsm_pvy=nullptr,*bsm_pvz=nullptr;
    std::vector<double> *bsm_dvx=nullptr,*bsm_dvy=nullptr,*bsm_dvz=nullptr;
    if (bsmTree && !isDijet && useTruth) {
        bsmTree->SetBranchAddress("pdgId",       &bsm_pdgId);
        bsmTree->SetBranchAddress("pt",          &bsm_pt);
        bsmTree->SetBranchAddress("eta",         &bsm_eta);
        bsmTree->SetBranchAddress("phi",         &bsm_phi);
        bsmTree->SetBranchAddress("hasProdVtx",  &bsm_hasProdVtx);
        bsmTree->SetBranchAddress("prodVtx_x",   &bsm_pvx);
        bsmTree->SetBranchAddress("prodVtx_y",   &bsm_pvy);
        bsmTree->SetBranchAddress("prodVtx_z",   &bsm_pvz);
        bsmTree->SetBranchAddress("hasDecayVtx", &bsm_hasDecayVtx);
        bsmTree->SetBranchAddress("decayVtx_x",  &bsm_dvx);
        bsmTree->SetBranchAddress("decayVtx_y",  &bsm_dvy);
        bsmTree->SetBranchAddress("decayVtx_z",  &bsm_dvz);
    }

    // ---- shower-parent LLPs (showerParentTree) ----
    // TruthBSM carries no decay vertices, which is why the truthBSM overlay above
    // draws nothing and useTruth defaults to false. showerParentTree is the curated
    // collection that DOES have them (BSM parent + its visible decay products), so
    // it is what actually labels the displacement: decay vertex, flight path from
    // the IP.
    std::vector<double> *sp_pt=nullptr,*sp_eta=nullptr,*sp_phi=nullptr,*sp_m=nullptr;
    std::vector<double> *sp_dvx=nullptr,*sp_dvy=nullptr,*sp_dvz=nullptr;
    std::vector<double> *sp_lxy=nullptr,*sp_r3d=nullptr;
    const bool haveSP = spTree && spTree->GetEntries() > 0 && !isDijet;
    if (haveSP) {
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

    // ---- per-event bookkeeping, used to label the pages ----
    // For background: which JZ slice the event came from, its rate contribution
    // (eventWeights[0], in Hz because the ntupler normalizes to 1 s of HL-LHC
    // luminosity), and whether it passed the HSTP filter.
    int                  sampleJZSlice = -1;
    std::vector<double>* eventWeights  = nullptr;
    bool                 passHSTP      = true;
    const bool haveSlice   = evtTree && evtTree->GetBranch("sampleJZSlice") != nullptr;
    const bool haveWeights = evtTree && evtTree->GetBranch("eventWeights")  != nullptr;
    const bool haveHSTP    = evtTree && evtTree->GetBranch("passHSTP")      != nullptr;
    if (haveSlice)   evtTree->SetBranchAddress("sampleJZSlice", &sampleJZSlice);
    if (haveWeights) evtTree->SetBranchAddress("eventWeights",  &eventWeights);
    if (haveHSTP)    evtTree->SetBranchAddress("passHSTP",      &passHSTP);
    const bool haveEvtInfo = haveSlice || haveWeights || haveHSTP;
    const bool useHSTP     = isDijet && kSkipHSTPFailedBackground && haveHSTP;

    Long64_t nEvt = towTree->GetEntries();

    // ---- which entries to consider, in order ----------------------------------
    // Background: a fixed quota per JZ slice (kDisplaySlices x kEventsPerSlice)
    // rather than a blind stride over the chain. JZ0 is deliberately excluded --
    // the HSTP filter removes essentially all of it, so its events carry no rate.
    // Each slice lives in its own file, so the chain's tree offsets give the entry
    // range per slice directly and we never have to walk the slices we don't want.
    // Signal (one file): the old stride behaviour.
    std::vector<Long64_t> candidates;
    std::vector<int>      sliceOfCandidate;   // parallel; -1 when unknown
    const bool perSliceMode = isDijet && haveSlice && nInputFiles > 1;

    if (perSliceMode) {
        TChain* towChain = dynamic_cast<TChain*>(towTree);
        Long64_t* offs = towChain ? towChain->GetTreeOffset() : nullptr;
        const int nTrees = towChain ? towChain->GetNtrees() : 0;
        if (!offs || nTrees <= 0) {
            std::cerr << "[caloShowerEventDisplays] cannot read chain tree offsets -- "
                         "falling back to a stride over the whole chain\n";
        } else {
            // Identify each file's slice from its first entry rather than assuming
            // file order maps to JZ number.
            std::vector<int> sliceOfFile(nTrees, -1);
            for (int i = 0; i < nTrees; ++i) {
                evtTree->GetEntry(offs[i]);
                sliceOfFile[i] = sampleJZSlice;
            }
            for (int wanted : kDisplaySlices) {
                int fileIdx = -1;
                for (int i = 0; i < nTrees; ++i) if (sliceOfFile[i] == wanted) { fileIdx = i; break; }
                if (fileIdx < 0) {
                    std::cerr << "[caloShowerEventDisplays] JZ" << wanted << " not in the chain\n";
                    continue;
                }
                const Long64_t start = offs[fileIdx];
                const Long64_t end   = (fileIdx + 1 < nTrees) ? offs[fileIdx + 1] : nEvt;
                // Spread the quota across the slice, with headroom for events that
                // get skipped (no jets, HSTP-failed).
                const Long64_t span   = std::max<Long64_t>(1, end - start);
                const Long64_t stride = std::max<Long64_t>(1, span / (4 * (Long64_t)kEventsPerSlice));
                for (Long64_t e = start; e < end; e += stride) {
                    candidates.push_back(e);
                    sliceOfCandidate.push_back(wanted);
                }
                std::cout << "[caloShowerEventDisplays] JZ" << wanted << ": entries [" << start
                          << ", " << end << ")  stride " << stride << "  quota "
                          << kEventsPerSlice << "\n";
            }
        }
    }
    if (candidates.empty()) {
        // Stride over everything (signal, or background without slice info).
        Long64_t evtStride = 1;
        if (nInputFiles > 1 && nEventsToDraw > 0)
            evtStride = std::max<Long64_t>(1, nEvt / (2 * (Long64_t)nEventsToDraw));
        if (evtStride > 1)
            std::cout << "[caloShowerEventDisplays] chained input: stepping " << nEvt
                      << " events in strides of " << evtStride << " to span all files\n";
        for (Long64_t e = 0; e < nEvt; e += evtStride) {
            candidates.push_back(e);
            sliceOfCandidate.push_back(-1);
        }
    }
    const bool useSliceQuota = perSliceMode && !candidates.empty()
                               && sliceOfCandidate.front() >= 0;
    const int maxDrawTotal = useSliceQuota
                           ? (int)kDisplaySlices.size() * kEventsPerSlice
                           : nEventsToDraw;

    // Two jet collections, each -> its own PDF.
    struct JetColl { const char* tag; bool isLRJ; double Rassoc; };
    std::vector<JetColl> colls = {
        { "WTACone", false, 0.4 },   // WTACone EtaSK; Rassoc = association cone (tune)
    };
    // LRJ off: the collection is empty in the JZ1 (v22 PU200) production and the
    // large-R side is not what this study measures. Matches kRunLRJ in
    // caloShowerShapePlots.C and kWriteJetTaggerLRJ in the ntupler.
    if (kRunLRJ) colls.push_back({ "LRJ", true, 1.0 });   // JetTaggerLRJ EtaSK (large-R)

    for (const auto& coll : colls) {
        if (coll.isLRJ && !lrjTree) { std::cerr << "no LRJ tree, skipping\n"; continue; }
        if (!coll.isLRJ && !wtaTree) { std::cerr << "no WTACone tree, skipping\n"; continue; }

        std::string pdfName = outDir + "/caloShowerEventDisplays_" + coll.tag + "_" + sampleLabel + ".pdf";
        TString pdf = pdfName.c_str();
        TCanvas* c = new TCanvas("cEvt", "event display", 1500, 1100);
        c->Print(pdf + "[");   // open multipage

        int drawn = 0;
        std::map<int,int> drawnPerSlice;
        for (size_t iCand = 0; iCand < candidates.size() && drawn < maxDrawTotal; ++iCand) {
            const Long64_t iEvt = candidates[iCand];
            // Per-slice quota: stop spending pages on a slice once it has its 10.
            if (useSliceQuota && drawnPerSlice[sliceOfCandidate[iCand]] >= kEventsPerSlice)
                continue;
            towTree->GetEntry(iEvt);
            if (coll.isLRJ) lrjTree->GetEntry(iEvt); else wtaTree->GetEntry(iEvt);
            if (bsmTree && !isDijet) bsmTree->GetEntry(iEvt);
            if (haveSP) spTree->GetEntry(iEvt);
            if (haveEvtInfo) evtTree->GetEntry(iEvt);
            // An HSTP-failed background event carries no rate, so skip it rather
            // than spend a page on it.
            if (useHSTP && !passHSTP) continue;

            // ---- jet list (eta,phi,pt) ----
            std::vector<std::array<double,3>> jets; // eta,phi,pt
            if (coll.isLRJ) {
                if (lrj_Et) for (size_t j=0;j<lrj_Et->size();++j)
                    jets.push_back({ (double)(*lrj_eta)[j], (double)(*lrj_phi)[j], (double)(*lrj_Et)[j] });
            } else {
                if (wta_pt) for (size_t j=0;j<wta_pt->size();++j)
                    jets.push_back({ (*wta_eta)[j], (*wta_phi)[j], (*wta_pt)[j] });
            }
            if (jets.empty()) continue;

            // ---- pick the jets to draw ----
            // signal : jets ΔR-matched to a displaced BSM particle (pt>ptMinBSM, has decay vtx)
            // dijet  : the two leading jets
            std::vector<int> selJets;        // index into jets
            std::vector<int> matchedBSM;     // parallel: BSM index matched (or -1)
            if (!isDijet && useTruth && bsm_pt) {
                std::vector<int> bsmOfInterest;
                for (size_t b=0;b<bsm_pt->size();++b) {
                    bool disp = bsm_hasDecayVtx && (*bsm_hasDecayVtx)[b] == 1;
                    if (disp && (*bsm_pt)[b] > ptMinBSM) bsmOfInterest.push_back((int)b);
                }
                for (size_t j=0;j<jets.size();++j) {
                    int best=-1; double bestdr=coll.Rassoc;
                    for (int b : bsmOfInterest) {
                        double d = dR(jets[j][0], jets[j][1], (*bsm_eta)[b], (*bsm_phi)[b]);
                        if (d < bestdr) { bestdr=d; best=b; }
                    }
                    if (best>=0) { selJets.push_back((int)j); matchedBSM.push_back(best); }
                }
            } else {
                std::vector<int> order(jets.size());
                for (size_t j=0;j<jets.size();++j) order[j]=(int)j;
                std::sort(order.begin(),order.end(),[&](int a,int b){return jets[a][2]>jets[b][2];});
                for (int k=0;k<(int)order.size() && k<2;++k) { selJets.push_back(order[k]); matchedBSM.push_back(-1); }
            }
            if (selJets.empty()) continue;

            // =======================================================
            // Draw one page: 2x2 -> rz, xy, 3d, legend/info
            // =======================================================
            std::vector<TObject*> garbage;    // page-scoped, deleted after Print
            std::vector<double>   pointingDca3D; // per selected jet: shower-line closest approach to IP (m)
            // truth counterparts for the matched shower parent, per selected jet (mm / GeV)
            std::vector<double>   truthDca3DMM, truthLxyMM, truthR3dMM, truthPtGeV, truthEta, truthMassGeV;

            // ---- towers in each selected jet's cone, and how many layers are lit ----
            // This is the number the DCA3D fit actually has to work with: EtaSK plus
            // the per-layer etMinTower cut can leave a jet with 2-3 towers and only
            // 1-2 lit layers, at which point the fitted line is meaningless.
            std::vector<int> nTowInJet(selJets.size(), 0), nLayersInJet(selJets.size(), 0);
            for (size_t ji=0; ji<selJets.size(); ++ji) {
                const int sj = selJets[ji];
                bool layerLit[7] = {false,false,false,false,false,false,false};
                for (size_t t=0; t<tow_Et->size(); ++t) {
                    if (dR((*tow_Eta)[t],(*tow_Phi)[t],jets[sj][0],jets[sj][1]) >= coll.Rassoc) continue;
                    bool any=false;
                    for (int l=0;l<7;++l) {
                        double et = tow_Et_l[l] ? (*tow_Et_l[l])[t] : 0.0;
                        if (et >= etMinTower) { layerLit[l]=true; any=true; }
                    }
                    if (any) ++nTowInJet[ji];
                }
                for (int l=0;l<7;++l) if (layerLit[l]) ++nLayersInJet[ji];
            }
            c->Clear();
            c->Divide(2,2);

            // ---- frames (unique names, detached from file dir to avoid warnings) ----
            c->cd(1);
            TH2F* frz = new TH2F(Form("frz_%s_%lld",coll.tag,iEvt), Form("r-z  (event %lld);z [m];signed r [m]", iEvt), 10,-6.2,6.2, 10,-4.2,4.2);
            frz->SetDirectory(nullptr); frz->Draw(); garbage.push_back(frz);
            c->cd(2);
            TH2F* fxy = new TH2F(Form("fxy_%s_%lld",coll.tag,iEvt), "transverse x-y;x [m];y [m]", 10,-4.2,4.2, 10,-4.2,4.2);
            fxy->SetDirectory(nullptr); fxy->Draw(); garbage.push_back(fxy);
            c->cd(3);
            TH3F* f3 = new TH3F(Form("f3_%s_%lld",coll.tag,iEvt), ";x [m];y [m];z [m]", 1,-4.2,4.2, 1,-4.2,4.2, 1,-6.2,6.2);
            f3->SetDirectory(nullptr); f3->Draw(); garbage.push_back(f3);

            // ---- nominal calo geometry (black dotted), drawn behind the deposits ----
            drawCaloGeometry(c, garbage);

            // one TPolyMarker3D per layer (3D can't size per-point; color=layer)
            TPolyMarker3D* pm3[7];
            for (int l=0;l<7;++l){ pm3[l]=new TPolyMarker3D(); pm3[l]->SetMarkerColor(kLayerColor[l]);
                                   pm3[l]->SetMarkerStyle(20); pm3[l]->SetMarkerSize(0.6); garbage.push_back(pm3[l]); }

            // ---- deposits: towers within Rassoc of any selected jet ----
            for (size_t t=0; t<tow_Et->size(); ++t) {
                double te=(*tow_Eta)[t], tp=(*tow_Phi)[t];
                bool inJet=false;
                for (int sj : selJets) if (dR(te,tp,jets[sj][0],jets[sj][1]) < coll.Rassoc) { inJet=true; break; }
                if (!inJet) continue;
                for (int l=0;l<7;++l) {
                    double et = tow_Et_l[l] ? (*tow_Et_l[l])[t] : 0.0;
                    if (et < etMinTower) continue;
                    double x,y,z; layerXYZ(l, te, tp, x,y,z);
                    double ms = etMarkerSize(et);
                    // rz
                    c->cd(1);
                    TMarker* mrz=new TMarker(z, signedR(x,y), 20); mrz->SetMarkerColor(kLayerColor[l]);
                    mrz->SetMarkerSize(ms); mrz->Draw(); garbage.push_back(mrz);
                    // xy
                    c->cd(2);
                    TMarker* mxy=new TMarker(x, y, 20); mxy->SetMarkerColor(kLayerColor[l]);
                    mxy->SetMarkerSize(ms); mxy->Draw(); garbage.push_back(mxy);
                    // 3d
                    pm3[l]->SetNextPoint(x,y,z);
                }
            }
            c->cd(3); for (int l=0;l<7;++l) if (pm3[l]->Size()>0) pm3[l]->Draw();

            // ---- jet axes (IP -> calo at r=2m) ----
            // Leading vs subleading jet get a slightly different shade so their cone,
            // axis and shower line can be told apart (deposits stay layer-coloured).
            for (size_t ji=0; ji<selJets.size(); ++ji) {
                int sj = selJets[ji];
                int axisCol = (ji==0 ? kGray+2   : kGray+1);
                int coneCol = (ji==0 ? kViolet+1 : kViolet-4);
                double je=jets[sj][0], jp=jets[sj][1];
                double rr=2.0, xx=rr*std::cos(jp), yy=rr*std::sin(jp), zz=rr*std::sinh(je);
                c->cd(1); TLine* lrz=new TLine(0,0,zz,signedR(xx,yy)); lrz->SetLineColor(axisCol);
                lrz->SetLineStyle(2); lrz->Draw(); garbage.push_back(lrz);
                c->cd(2); TLine* lxy=new TLine(0,0,xx,yy); lxy->SetLineColor(axisCol);
                lxy->SetLineStyle(2); lxy->Draw(); garbage.push_back(lxy);
                TPolyLine3D* l3=new TPolyLine3D(2); l3->SetPoint(0,0,0,0); l3->SetPoint(1,xx,yy,zz);
                l3->SetLineColor(axisCol); l3->SetLineStyle(2); c->cd(3); l3->Draw(); garbage.push_back(l3);

                // ---- jet association cone (DeltaR = Rassoc) ----
                drawJetCone(c, je, jp, coll.Rassoc, coneCol, garbage);
            }

            // ---- shower pointing: Et-weighted per-layer-centroid line fit ----
            // Straight line through the 7 per-layer energy centroids; NOT forced
            // through the IP, so a displaced shower's line misses the origin (dca3D>0).
            for (size_t ji=0; ji<selJets.size(); ++ji) {
                int sj = selJets[ji];
                int showerCol = (ji==0 ? kGreen+2 : kGreen-6);
                double je=jets[sj][0], jp=jets[sj][1];
                double cc[3], dd[3], dca3D=0;
                if (!jetShowerPointing(je,jp,coll.Rassoc,etMinTower,tow_Eta,tow_Phi,tow_Et_l,cc,dd,dca3D)) continue;
                pointingDca3D.push_back(dca3D);
                // Draw the fit segment only where it is meaningful: from just past
                // the IP out to the outer calo layer, on the shower side. A fixed
                // +/-8 m span let a small fit tilt throw the ends across the canvas
                // (and a second jet's line into a spurious "X").
                double cd2  = cc[0]*dd[0]+cc[1]*dd[1]+cc[2]*dd[2];
                double t0   = -cd2;                       // closest-approach (dca3D) point
                double Rout = kRbarrel[6] + 0.2;
                double reach= std::sqrt(std::max(0.0, Rout*Rout - dca3D*dca3D));
                double sgn  = (cd2 >= 0 ? 1.0 : -1.0);    // extend toward the deposits
                double tlo  = t0 - 0.5*sgn, thi = t0 + sgn*(reach + 0.5);
                double ax=cc[0]+tlo*dd[0], ay=cc[1]+tlo*dd[1], az=cc[2]+tlo*dd[2];
                double bx=cc[0]+thi*dd[0], by=cc[1]+thi*dd[1], bz=cc[2]+thi*dd[2];
                // r-z: sample so the signed-r sign flip near y=0 is drawn faithfully
                const int NS=24; TPolyLine* prz=new TPolyLine(NS);
                for (int k=0;k<NS;++k){ double tt=tlo + (thi-tlo)*k/(NS-1);
                    double x=cc[0]+tt*dd[0], y=cc[1]+tt*dd[1], z=cc[2]+tt*dd[2];
                    prz->SetPoint(k, z, signedR(x,y)); }
                prz->SetLineColor(showerCol); prz->SetLineWidth(2);
                c->cd(1); prz->Draw(); garbage.push_back(prz);
                c->cd(2); TLine* pxy=new TLine(ax,ay,bx,by); pxy->SetLineColor(showerCol); pxy->SetLineWidth(2);
                pxy->Draw(); garbage.push_back(pxy);
                TPolyLine3D* p3=new TPolyLine3D(2); p3->SetPoint(0,ax,ay,az); p3->SetPoint(1,bx,by,bz);
                p3->SetLineColor(showerCol); p3->SetLineWidth(2); c->cd(3); p3->Draw(); garbage.push_back(p3);
            }

            // ---- truth overlay: IP + BSM prod/decay vertices + flight path ----
            // IP marker
            c->cd(1); { TMarker* ip=new TMarker(0,0,29); ip->SetMarkerColor(kBlack); ip->SetMarkerSize(1.4); ip->Draw(); garbage.push_back(ip);}
            c->cd(2); { TMarker* ip=new TMarker(0,0,29); ip->SetMarkerColor(kBlack); ip->SetMarkerSize(1.4); ip->Draw(); garbage.push_back(ip);}

            // ---- truth shower parents (showerParentTree), the displacement label ----
            // For each selected jet, the highest-pt LLP inside the association cone:
            // IP -> decay vertex (flight path, magenta solid) then decay vertex ->
            // magenta solid). No ray is drawn past the decay vertex: it would lie on
            // that same line (see below).
            for (size_t ji=0; ji<selJets.size() && haveSP && sp_pt; ++ji) {
                int sj = selJets[ji];
                double je=jets[sj][0], jp=jets[sj][1];
                int best=-1; double bestPt=-1;
                for (size_t s=0;s<sp_pt->size();++s){
                    if (dR(je,jp,(*sp_eta)[s],(*sp_phi)[s]) >= coll.Rassoc) continue;
                    if ((*sp_pt)[s] > bestPt){ bestPt=(*sp_pt)[s]; best=(int)s; }
                }
                if (best<0) continue;
                const int truthCol = (ji==0 ? kMagenta+2 : kMagenta-7);
                const double vx=(*sp_dvx)[best]/1000.0, vy=(*sp_dvy)[best]/1000.0, vz=(*sp_dvz)[best]/1000.0; // mm -> m
                const double se=(*sp_eta)[best], sp_=(*sp_phi)[best];
                truthDca3DMM.push_back(truthShowerDca3D((*sp_dvx)[best],(*sp_dvy)[best],(*sp_dvz)[best], se, sp_));
                truthLxyMM.push_back((*sp_lxy)[best]);
                truthR3dMM.push_back((*sp_r3d)[best]);
                truthPtGeV.push_back((*sp_pt)[best]);
                truthEta.push_back(se);
                truthMassGeV.push_back(sp_m ? (*sp_m)[best] : -1.0);
                // decay vertex marker
                c->cd(1); { TMarker* m=new TMarker(vz, signedR(vx,vy), 34); m->SetMarkerColor(truthCol);
                            m->SetMarkerSize(1.7); m->Draw(); garbage.push_back(m); }
                c->cd(2); { TMarker* m=new TMarker(vx, vy, 34); m->SetMarkerColor(truthCol);
                            m->SetMarkerSize(1.7); m->Draw(); garbage.push_back(m); }
                { TPolyMarker3D* m3=new TPolyMarker3D(1); m3->SetPoint(0,vx,vy,vz);
                  m3->SetMarkerColor(truthCol); m3->SetMarkerStyle(34); m3->SetMarkerSize(1.4);
                  c->cd(3); m3->Draw(); garbage.push_back(m3); }
                // IP -> decay vertex (the LLP flight path)
                c->cd(1); { TLine* l=new TLine(0,0,vz,signedR(vx,vy)); l->SetLineColor(truthCol);
                            l->SetLineWidth(2); l->Draw(); garbage.push_back(l); }
                c->cd(2); { TLine* l=new TLine(0,0,vx,vy); l->SetLineColor(truthCol);
                            l->SetLineWidth(2); l->Draw(); garbage.push_back(l); }
                { TPolyLine3D* l3=new TPolyLine3D(2); l3->SetPoint(0,0,0,0); l3->SetPoint(1,vx,vy,vz);
                  l3->SetLineColor(truthCol); l3->SetLineWidth(2); c->cd(3); l3->Draw(); garbage.push_back(l3); }
                // No truth "shower ray": the LLP is produced at the IP and travels
                // along u, so its decay vertex satisfies v = t*u and the ray v + s*u
                // is the SAME line as the IP -> decay-vertex flight path already
                // drawn above. It added a second copy of that line, extended well
                // past the calorimeter and off the pad, and nothing else.
            }
            if (!isDijet && useTruth && bsm_pt) {
                for (int b : matchedBSM) {
                    if (b < 0) continue;
                    // vertices in mm -> m
                    double pvx=0,pvy=0,pvz=0, dvx=0,dvy=0,dvz=0;
                    bool hasP = bsm_hasProdVtx && (*bsm_hasProdVtx)[b]==1;
                    bool hasD = bsm_hasDecayVtx && (*bsm_hasDecayVtx)[b]==1;
                    if (hasP){ pvx=(*bsm_pvx)[b]/1000.; pvy=(*bsm_pvy)[b]/1000.; pvz=(*bsm_pvz)[b]/1000.; }
                    if (hasD){ dvx=(*bsm_dvx)[b]/1000.; dvy=(*bsm_dvy)[b]/1000.; dvz=(*bsm_dvz)[b]/1000.; }
                    // decay vertex marker (star)
                    if (hasD) {
                        c->cd(1); TMarker* dm=new TMarker(dvz, signedR(dvx,dvy), 29); dm->SetMarkerColor(kBlack); dm->SetMarkerSize(1.8); dm->Draw(); garbage.push_back(dm);
                        c->cd(2); TMarker* dm2=new TMarker(dvx, dvy, 29); dm2->SetMarkerColor(kBlack); dm2->SetMarkerSize(1.8); dm2->Draw(); garbage.push_back(dm2);
                    }
                    // flight path prod -> decay (solid black)
                    if (hasP && hasD) {
                        c->cd(1); TLine* fp=new TLine(pvz,signedR(pvx,pvy),dvz,signedR(dvx,dvy)); fp->SetLineColor(kBlack); fp->SetLineWidth(2); fp->Draw(); garbage.push_back(fp);
                        c->cd(2); TLine* fp2=new TLine(pvx,pvy,dvx,dvy); fp2->SetLineColor(kBlack); fp2->SetLineWidth(2); fp2->Draw(); garbage.push_back(fp2);
                        TPolyLine3D* fp3=new TPolyLine3D(2); fp3->SetPoint(0,pvx,pvy,pvz); fp3->SetPoint(1,dvx,dvy,dvz); fp3->SetLineColor(kBlack); fp3->SetLineWidth(2); c->cd(3); fp3->Draw(); garbage.push_back(fp3);
                    }
                    // decay -> calo along parent direction (dashed) from decay vtx
                    if (hasD) {
                        double be=(*bsm_eta)[b], bp=(*bsm_phi)[b];
                        double rr=2.0, cx=rr*std::cos(bp), cy=rr*std::sin(bp), cz=rr*std::sinh(be);
                        c->cd(1); TLine* dc=new TLine(dvz,signedR(dvx,dvy),dvz+cz,signedR(dvx+cx,dvy+cy)); dc->SetLineColor(kBlack); dc->SetLineStyle(3); dc->Draw(); garbage.push_back(dc);
                        c->cd(2); TLine* dc2=new TLine(dvx,dvy,dvx+cx,dvy+cy); dc2->SetLineColor(kBlack); dc2->SetLineStyle(3); dc2->Draw(); garbage.push_back(dc2);
                    }
                }
            }

            // ---- info / legend pad ----
            c->cd(4);
            TLegend* leg = new TLegend(0.05,0.30,0.55,0.95);
            leg->SetHeader(sampleJZSlice >= 0
                           ? Form("%s  |  %s  |  JZ%d  |  event %lld", coll.tag, sampleLabel.c_str(), sampleJZSlice, iEvt)
                           : Form("%s  |  %s  |  event %lld", coll.tag, sampleLabel.c_str(), iEvt));
            for (int l=0;l<7;++l){ TMarker* m=new TMarker(0,0,20); m->SetMarkerColor(kLayerColor[l]); leg->AddEntry(m,kLayerName[l],"p"); garbage.push_back(m);}
            {   TMarker* mip=new TMarker(0,0,29); mip->SetMarkerColor(kBlack); leg->AddEntry(mip,"IP / BSM decay vtx","p"); garbage.push_back(mip);
                TLine* lf=new TLine(); lf->SetLineColor(kBlack); lf->SetLineWidth(2); leg->AddEntry(lf,"BSM flight path","l"); garbage.push_back(lf);
                TLine* ld=new TLine(); ld->SetLineColor(kBlack); ld->SetLineStyle(3); leg->AddEntry(ld,"decay #rightarrow calo dir","l"); garbage.push_back(ld);
                TLine* lj=new TLine(); lj->SetLineColor(kGray+2); lj->SetLineStyle(2); leg->AddEntry(lj,"jet axis (from IP)","l"); garbage.push_back(lj);
                TLine* lp=new TLine(); lp->SetLineColor(kGreen+2); lp->SetLineWidth(2); leg->AddEntry(lp,"shower pointing (layer fit)","l"); garbage.push_back(lp);
                TLine* lg=new TLine(); lg->SetLineColor(kBlack); lg->SetLineStyle(3); leg->AddEntry(lg,"calo layers (nominal)","l"); garbage.push_back(lg);
                TLine* lc=new TLine(); lc->SetLineColor(kViolet+1); lc->SetLineWidth(2); leg->AddEntry(lc,Form("jet cone (#DeltaR=%.1f)",coll.Rassoc),"l"); garbage.push_back(lc);
                TLine* l1=new TLine(); l1->SetLineColor(kViolet+1); l1->SetLineWidth(2); leg->AddEntry(l1,"leading jet (cone/axis/line)","l"); garbage.push_back(l1);
                TLine* l2=new TLine(); l2->SetLineColor(kViolet-4); l2->SetLineWidth(2); leg->AddEntry(l2,"subleading jet (lighter shade)","l"); garbage.push_back(l2);
                if (haveSP) {
                    TMarker* mt=new TMarker(0,0,34); mt->SetMarkerColor(kMagenta+2); leg->AddEntry(mt,"truth LLP decay vertex","p"); garbage.push_back(mt);
                    TLine* lt=new TLine(); lt->SetLineColor(kMagenta+2); lt->SetLineWidth(2); leg->AddEntry(lt,"truth LLP flight path","l"); garbage.push_back(lt);
                }
            }
            leg->SetBorderSize(0); leg->Draw(); garbage.push_back(leg);
            TLatex* note=new TLatex(); note->SetNDC(); note->SetTextSize(0.030);
            note->DrawLatex(0.05,0.26,Form("#splitline{marker size #propto layer E_{T}}{%d matched jet(s)}",(int)selJets.size()));
            // Towers / lit layers per jet, in leading-then-subleading order. Printed
            // next to the fitted DCA3D because the two go together: a jet with a
            // handful of towers spread over one or two layers has no lever arm, and
            // its DCA3D is resolution rather than displacement.
            {
                TString sn="n towers (lead, sublead):";
                for (size_t i=0;i<nTowInJet.size();++i) sn += Form("%s%d", i?", ":" ", nTowInJet[i]);
                note->DrawLatex(0.05,0.21,sn);
                TString sl="n lit layers:";
                for (size_t i=0;i<nLayersInJet.size();++i) sl += Form("%s%d/7", i?", ":" ", nLayersInJet[i]);
                note->DrawLatex(0.05,0.165,sl);
            }
            if (!pointingDca3D.empty()) {
                TString s="fitted DCA_{3D}#approx";
                for (size_t i=0;i<pointingDca3D.size();++i) s += Form("%s%.0f", i?", ":" ", pointingDca3D[i]*1000.);
                s += " mm";
                note->DrawLatex(0.05,0.135,s);
            }
            // Background: no truth to overlay, so give the event's weight in the
            // JZ mixture instead. The weight is per-second, so it is
            // literally this event's contribution to the trigger rate in Hz — a
            // JZ0/JZ1 event is worth O(kHz) on its own while a JZ9 event is worth
            // O(1e-4) Hz, which is what makes an unweighted display misleading.
            if (isDijet) {
                if (haveWeights && eventWeights && !eventWeights->empty()) {
                    TString sr = Form("JZ%d  rate contribution = %.3g Hz",
                                      sampleJZSlice, eventWeights->front());
                    note->DrawLatex(0.05,0.10,sr);
                }
                if (haveHSTP)
                    note->DrawLatex(0.05,0.055,Form("passHSTP = %s", passHSTP ? "true" : "false"));
            }
            // Truth numbers for the matched shower parents, in the same jet order, so
            // the fitted DCA3D above can be read against the truth it should recover.
            if (!truthDca3DMM.empty()) {
                TString st="truth DCA_{3D}#approx";
                for (size_t i=0;i<truthDca3DMM.size();++i) st += Form("%s%.0f", i?", ":" ", truthDca3DMM[i]);
                st += " mm";
                note->DrawLatex(0.05,0.11,st);
                TString sv="truth L_{xy}/|r|#approx";
                for (size_t i=0;i<truthLxyMM.size();++i)
                    sv += Form("%s%.0f/%.0f", i?", ":" ", truthLxyMM[i], truthR3dMM[i]);
                sv += " mm";
                note->DrawLatex(0.05,0.06,sv);
                TString sp2="truth parent p_{T}/#eta/m:";
                for (size_t i=0;i<truthPtGeV.size();++i)
                    sp2 += Form("%s%.0f GeV/%.2f/%.1f GeV", i?", ":" ",
                                truthPtGeV[i], truthEta[i], truthMassGeV[i]);
                note->DrawLatex(0.05,0.01,sp2);
            }
            garbage.push_back(note);

            c->Print(pdf);
            ++drawn;
            if (useSliceQuota) ++drawnPerSlice[sliceOfCandidate[iCand]];
            for (auto* o : garbage) delete o;
        }

        c->Print(pdf + "]");   // close multipage
        std::cout << "[caloShowerEventDisplays] wrote " << drawn << " pages -> " << pdf << "\n";
        delete c;
    }

    fin->Close();
}
