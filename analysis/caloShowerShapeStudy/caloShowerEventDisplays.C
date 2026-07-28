// caloShowerEventDisplays.C
// ---------------------------------------------------------------------------
// Per-event calorimeter event displays for the displaced-jet / shower-shape
// study. Reads the caloShowerShapeNTupler.C output and, for a number of events,
// draws the per-layer calorimeter energy deposits of the jets that are ΔR-matched
// to displaced truth-BSM particles (signal), or of the leading jets (dijet
// background). Three views per event: r-z, transverse (x-y), and 3D.
//
//   # both signals, once each, into plots/<signal>/ (no arguments needed):
//   root -b -l -q 'caloShowerEventDisplays.C'
//   # or one sample explicitly:
//   root -b -q 'caloShowerEventDisplays.C("in.root","displaced_dark_photon",false,20,"plots/")'
//   root -b -q 'caloShowerEventDisplays.C("dijet.root","dijet_JZ4",true,20,"plots/")'
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
// ---------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "TFile.h"
#include "TTree.h"
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

// ---------------------------------------------------------------------------
// Nominal calorimeter layer geometry (meters). Illustrative only.
// Layer grouping matches the ntupler (l0..l6):
//   l0 PreSampler | l1 EM1 | l2 EM2 | l3 EM3 | l4 Tile0 | l5 Tile1 | l6 Tile2
// ---------------------------------------------------------------------------
static const double kRbarrel[7] = { 1.45, 1.53, 1.63, 1.74, 2.45, 3.02, 3.63 };
static const double kZendcap[7] = { 3.60, 3.75, 3.87, 3.97, 4.40, 5.00, 5.60 };
static const double kEtaBarrelLimit = 1.5;
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

// (layer, eta, phi) -> nominal (x,y,z) in meters
static void layerXYZ(int l, double eta, double phi, double& x, double& y, double& z) {
    double r, zz;
    if (std::fabs(eta) < kEtaBarrelLimit) {
        r  = kRbarrel[l];
        zz = r * std::sinh(eta);
    } else {
        zz = (eta >= 0 ? 1.0 : -1.0) * kZendcap[l];
        double s = std::sinh(eta);
        r = (std::fabs(s) > 1e-6) ? std::fabs(kZendcap[l] / s) : kRbarrel[l];
    }
    x = r * std::cos(phi);
    y = r * std::sin(phi);
    z = zz;
}

static inline double signedR(double x, double y) {
    return std::sqrt(x*x + y*y) * (y >= 0 ? 1.0 : -1.0);
}
static inline double etMarkerSize(double et) {
    double s = 0.3 + 0.28 * std::sqrt(std::max(0.0, et));
    return std::min(3.0, s);
}

// Et-weighted per-layer-centroid straight-line fit of a jet's calorimeter shower.
// For each layer l0..l6 we take the Et-weighted mean tower position (within Rassoc
// of the jet), then fit a 3D line (weighted principal axis) through those up-to-7
// layer centroids. Unlike the jet axis this line is NOT forced through the IP, so
// for a displaced shower it points back to an offset origin. Returns false if
// fewer than two layers have energy. Fills line anchor c[3], unit direction d[3],
// and d0 = 3D closest approach of the fitted line to the origin (m).
static bool jetShowerPointing(double je, double jp, double Rassoc, double etMinTower,
                              const std::vector<double>* towEta,
                              const std::vector<double>* towPhi,
                              std::vector<double>* const towEtL[7],
                              double c[3], double d[3], double& d0)
{
    if (!towEta || !towPhi) return false;
    double Lx[7]={0,0,0,0,0,0,0}, Ly[7]={0,0,0,0,0,0,0}, Lz[7]={0,0,0,0,0,0,0}, Lw[7]={0,0,0,0,0,0,0};
    for (size_t t=0; t<towEta->size(); ++t) {
        double te=(*towEta)[t], tp=(*towPhi)[t];
        if (dR(te,tp,je,jp) >= Rassoc) continue;
        for (int l=0;l<7;++l) {
            double et = towEtL[l] ? (*towEtL[l])[t] : 0.0;
            if (et < etMinTower) continue;
            double x,y,z; layerXYZ(l,te,tp,x,y,z);
            Lx[l]+=et*x; Ly[l]+=et*y; Lz[l]+=et*z; Lw[l]+=et;
        }
    }
    std::vector<std::array<double,3>> P; std::vector<double> W;
    for (int l=0;l<7;++l) if (Lw[l]>0) { P.push_back({Lx[l]/Lw[l],Ly[l]/Lw[l],Lz[l]/Lw[l]}); W.push_back(Lw[l]); }
    if (P.size()<2) return false;

    double sw=0; c[0]=c[1]=c[2]=0;
    for (size_t i=0;i<P.size();++i){ sw+=W[i]; for(int k=0;k<3;++k) c[k]+=W[i]*P[i][k]; }
    for(int k=0;k<3;++k) c[k]/=sw;

    TMatrixDSym M(3);
    for(int a=0;a<3;++a) for(int b=0;b<3;++b) M(a,b)=0.0;
    for (size_t i=0;i<P.size();++i){
        double dp[3]={P[i][0]-c[0],P[i][1]-c[1],P[i][2]-c[2]};
        for(int a=0;a<3;++a) for(int b=0;b<3;++b) M(a,b)+=W[i]*dp[a]*dp[b];
    }
    TMatrixDSymEigen eig(M);
    TVectorD val=eig.GetEigenValues();
    TMatrixD vec=eig.GetEigenVectors();
    int imax=0; for(int i=1;i<3;++i) if(val[i]>val[imax]) imax=i;
    for(int k=0;k<3;++k) d[k]=vec(k,imax);
    double n=std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]); if(n>0) for(int k=0;k<3;++k) d[k]/=n;

    double cd=c[0]*d[0]+c[1]*d[1]+c[2]*d[2];
    double px=c[0]-cd*d[0], py=c[1]-cd*d[1], pz=c[2]-cd*d[2];
    d0=std::sqrt(px*px+py*py+pz*pz);
    return true;
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

// Draw the jet association cone (DeltaR = Rassoc in eta-phi) as a wireframe: a
// boundary ring placed at the outer calo radius plus a few generator lines back to
// the IP, in all 3 pads. Uses the same projective placement as the deposits.
static void drawJetCone(TCanvas* c, double je, double jp, double Rassoc,
                        std::vector<TObject*>& garbage) {
    const int col = kViolet+1;
    const double Rdraw = kRbarrel[6];   // boundary ring at the outer calo layer
    const int N = 60;
    std::vector<double> bx(N+1), by(N+1), bz(N+1);
    for (int k=0;k<=N;++k) {
        double t  = 2*M_PI*k/N;
        double eb = je + Rassoc*std::cos(t);
        double pb = jp + Rassoc*std::sin(t);
        bx[k]=Rdraw*std::cos(pb); by[k]=Rdraw*std::sin(pb); bz[k]=Rdraw*std::sinh(eb);
    }
    // boundary ring in each view
    c->cd(1); { TPolyLine* r=new TPolyLine(N+1); for(int k=0;k<=N;++k) r->SetPoint(k,bz[k],signedR(bx[k],by[k]));
                r->SetLineColor(col); r->SetLineWidth(2); r->Draw(); garbage.push_back(r); }
    c->cd(2); { TPolyLine* r=new TPolyLine(N+1); for(int k=0;k<=N;++k) r->SetPoint(k,bx[k],by[k]);
                r->SetLineColor(col); r->SetLineWidth(2); r->Draw(); garbage.push_back(r); }
    { TPolyLine3D* r=new TPolyLine3D(N+1); for(int k=0;k<=N;++k) r->SetPoint(k,bx[k],by[k],bz[k]);
      r->SetLineColor(col); r->SetLineWidth(2); c->cd(3); r->Draw(); garbage.push_back(r); }
    // generator lines IP -> ring (3D cone surface only, to keep the 2D views clean)
    for (int k=0;k<N;k+=10) {
        TPolyLine3D* g=new TPolyLine3D(2); g->SetPoint(0,0,0,0); g->SetPoint(1,bx[k],by[k],bz[k]);
        g->SetLineColor(col); c->cd(3); g->Draw(); garbage.push_back(g);
    }
}

// ===========================================================================
void caloShowerEventDisplays(std::string inputFile = "",
                             std::string sampleLabel = "signal",
                             bool        isDijet     = false,
                             int         nEventsToDraw = 20,
                             std::string outDir      = ".",
                             double      ptMinBSM    = 5.0,   // GeV, BSM select
                             double      etMinTower  = 0.5,   // GeV, layer deposit threshold
                             bool        useTruth    = false)  // false = leading jets, skip truth overlay
{
    // No-arg default (root -b -l -q 'caloShowerEventDisplays.C'): run both signals
    // once each (isDijet=false), each written into its own plots/<label>/ dir.
    // Pass an inputFile explicitly to process a single sample.
    if (inputFile.empty()) {
        const std::vector<std::pair<std::string,std::string>> samples = {
            {"/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_displaced_dark_photon.root", "displaced_dark_photon"},
            {"/data/larsonma/CaloShowerShapeTriggers/ntuples/caloShowerShape_emerging_jets.root",         "emerging_jets"},
        };
        // Truth-BSM matching selects no jets for these samples right now, so the
        // no-arg default draws leading jets (jets + towers only, no truth overlay).
        // Flip to true to require displaced-BSM matching once truth is sorted out.
        const bool useTruthDefault = false;
        for (const auto& s : samples) {
            std::string od = "plots/" + s.second;
            gSystem->mkdir(od.c_str(), kTRUE);
            std::cout << "[caloShowerEventDisplays] === " << s.second << " -> " << od
                      << "/  (useTruth=" << useTruthDefault << ") ===\n";
            caloShowerEventDisplays(s.first, s.second, false, nEventsToDraw, od, ptMinBSM, etMinTower, useTruthDefault);
        }
        return;
    }

    gStyle->SetOptStat(0);

    TFile* fin = TFile::Open(inputFile.c_str());
    if (!fin || fin->IsZombie()) { std::cerr << "cannot open " << inputFile << "\n"; return; }

    TTree* towTree = nullptr; fin->GetObject("gepCellsTowersEtaSKTree", towTree);
    TTree* wtaTree = nullptr; fin->GetObject("wtaConeCellsTowersEtaSKTree", wtaTree);
    TTree* lrjTree = nullptr; fin->GetObject("jetTaggerLRJEtaSKTree", lrjTree);
    TTree* bsmTree = nullptr; fin->GetObject("truthBSMTree", bsmTree);
    // ---- diagnostics: tree presence + entries ----
    std::cout << "[caloShowerEventDisplays] input: " << inputFile
              << (useTruth ? "  useTruth" : "  leadingJets") << "\n";
    for (auto pr : {std::make_pair("gepCellsTowersEtaSKTree",    towTree),
                    std::make_pair("wtaConeCellsTowersEtaSKTree", wtaTree),
                    std::make_pair("jetTaggerLRJEtaSKTree",       lrjTree),
                    std::make_pair("truthBSMTree",                bsmTree)}) {
        Long64_t ne = pr.second ? pr.second->GetEntries() : -1;
        std::cout << "    " << pr.first << " entries=" << ne
                  << (!pr.second ? "  [MISSING]" : (ne == 0 ? "  [EMPTY]" : "")) << "\n";
    }
    if (!towTree) { std::cerr << "no gepCellsTowersEtaSKTree in " << inputFile << "\n"; return; }

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

    Long64_t nEvt = towTree->GetEntries();

    // Two jet collections, each -> its own PDF.
    struct JetColl { const char* tag; bool isLRJ; double Rassoc; };
    std::vector<JetColl> colls = {
        { "WTACone", false, 0.4 },   // WTACone EtaSK; Rassoc = association cone (tune)
        { "LRJ",     true,  1.0 },   // JetTaggerLRJ EtaSK (large-R)
    };

    for (const auto& coll : colls) {
        if (coll.isLRJ && !lrjTree) { std::cerr << "no LRJ tree, skipping\n"; continue; }
        if (!coll.isLRJ && !wtaTree) { std::cerr << "no WTACone tree, skipping\n"; continue; }

        std::string pdfName = outDir + "/caloShowerEventDisplays_" + coll.tag + "_" + sampleLabel + ".pdf";
        TString pdf = pdfName.c_str();
        TCanvas* c = new TCanvas("cEvt", "event display", 1500, 1100);
        c->Print(pdf + "[");   // open multipage

        int drawn = 0;
        for (Long64_t iEvt = 0; iEvt < nEvt && drawn < nEventsToDraw; ++iEvt) {
            towTree->GetEntry(iEvt);
            if (coll.isLRJ) lrjTree->GetEntry(iEvt); else wtaTree->GetEntry(iEvt);
            if (bsmTree && !isDijet) bsmTree->GetEntry(iEvt);

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
            std::vector<double>   pointingD0; // per selected jet: shower-line closest approach to IP (m)
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
            for (int sj : selJets) {
                double je=jets[sj][0], jp=jets[sj][1];
                double rr=2.0, xx=rr*std::cos(jp), yy=rr*std::sin(jp), zz=rr*std::sinh(je);
                c->cd(1); TLine* lrz=new TLine(0,0,zz,signedR(xx,yy)); lrz->SetLineColor(kGray+2);
                lrz->SetLineStyle(2); lrz->Draw(); garbage.push_back(lrz);
                c->cd(2); TLine* lxy=new TLine(0,0,xx,yy); lxy->SetLineColor(kGray+2);
                lxy->SetLineStyle(2); lxy->Draw(); garbage.push_back(lxy);
                TPolyLine3D* l3=new TPolyLine3D(2); l3->SetPoint(0,0,0,0); l3->SetPoint(1,xx,yy,zz);
                l3->SetLineColor(kGray+2); l3->SetLineStyle(2); c->cd(3); l3->Draw(); garbage.push_back(l3);

                // ---- jet association cone (DeltaR = Rassoc) ----
                drawJetCone(c, je, jp, coll.Rassoc, garbage);
            }

            // ---- shower pointing: Et-weighted per-layer-centroid line fit ----
            // Straight line through the 7 per-layer energy centroids; NOT forced
            // through the IP, so a displaced shower's line misses the origin (d0>0).
            for (int sj : selJets) {
                double je=jets[sj][0], jp=jets[sj][1];
                double cc[3], dd[3], d0=0;
                if (!jetShowerPointing(je,jp,coll.Rassoc,etMinTower,tow_Eta,tow_Phi,tow_Et_l,cc,dd,d0)) continue;
                pointingD0.push_back(d0);
                const double L=8.0;   // half-length; the pad clips the drawn segment
                double ax=cc[0]-L*dd[0], ay=cc[1]-L*dd[1], az=cc[2]-L*dd[2];
                double bx=cc[0]+L*dd[0], by=cc[1]+L*dd[1], bz=cc[2]+L*dd[2];
                // r-z: sample so the signed-r sign flip near y=0 is drawn faithfully
                const int NS=24; TPolyLine* prz=new TPolyLine(NS);
                for (int k=0;k<NS;++k){ double tt=-L + 2.0*L*k/(NS-1);
                    double x=cc[0]+tt*dd[0], y=cc[1]+tt*dd[1], z=cc[2]+tt*dd[2];
                    prz->SetPoint(k, z, signedR(x,y)); }
                prz->SetLineColor(kGreen+2); prz->SetLineWidth(2);
                c->cd(1); prz->Draw(); garbage.push_back(prz);
                c->cd(2); TLine* pxy=new TLine(ax,ay,bx,by); pxy->SetLineColor(kGreen+2); pxy->SetLineWidth(2);
                pxy->Draw(); garbage.push_back(pxy);
                TPolyLine3D* p3=new TPolyLine3D(2); p3->SetPoint(0,ax,ay,az); p3->SetPoint(1,bx,by,bz);
                p3->SetLineColor(kGreen+2); p3->SetLineWidth(2); c->cd(3); p3->Draw(); garbage.push_back(p3);
            }

            // ---- truth overlay: IP + BSM prod/decay vertices + flight path ----
            // IP marker
            c->cd(1); { TMarker* ip=new TMarker(0,0,29); ip->SetMarkerColor(kBlack); ip->SetMarkerSize(1.4); ip->Draw(); garbage.push_back(ip);}
            c->cd(2); { TMarker* ip=new TMarker(0,0,29); ip->SetMarkerColor(kBlack); ip->SetMarkerSize(1.4); ip->Draw(); garbage.push_back(ip);}
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
            leg->SetHeader(Form("%s  |  %s  |  event %lld", coll.tag, sampleLabel.c_str(), iEvt));
            for (int l=0;l<7;++l){ TMarker* m=new TMarker(0,0,20); m->SetMarkerColor(kLayerColor[l]); leg->AddEntry(m,kLayerName[l],"p"); garbage.push_back(m);}
            {   TMarker* mip=new TMarker(0,0,29); mip->SetMarkerColor(kBlack); leg->AddEntry(mip,"IP / BSM decay vtx","p"); garbage.push_back(mip);
                TLine* lf=new TLine(); lf->SetLineColor(kBlack); lf->SetLineWidth(2); leg->AddEntry(lf,"BSM flight path","l"); garbage.push_back(lf);
                TLine* ld=new TLine(); ld->SetLineColor(kBlack); ld->SetLineStyle(3); leg->AddEntry(ld,"decay #rightarrow calo dir","l"); garbage.push_back(ld);
                TLine* lj=new TLine(); lj->SetLineColor(kGray+2); lj->SetLineStyle(2); leg->AddEntry(lj,"jet axis (from IP)","l"); garbage.push_back(lj);
                TLine* lp=new TLine(); lp->SetLineColor(kGreen+2); lp->SetLineWidth(2); leg->AddEntry(lp,"shower pointing (layer fit)","l"); garbage.push_back(lp);
                TLine* lg=new TLine(); lg->SetLineColor(kBlack); lg->SetLineStyle(3); leg->AddEntry(lg,"calo layers (nominal)","l"); garbage.push_back(lg);
                TLine* lc=new TLine(); lc->SetLineColor(kViolet+1); lc->SetLineWidth(2); leg->AddEntry(lc,Form("jet cone (#DeltaR=%.1f)",coll.Rassoc),"l"); garbage.push_back(lc);
            }
            leg->SetBorderSize(0); leg->Draw(); garbage.push_back(leg);
            TLatex* note=new TLatex(); note->SetNDC(); note->SetTextSize(0.030);
            note->DrawLatex(0.05,0.22,Form("#splitline{marker size #propto layer E_{T}}{#splitline{%d matched jet(s)}{nominal calo geometry}}",(int)selJets.size()));
            if (!pointingD0.empty()) {
                TString s="shower line d_{0}#approx";
                for (size_t i=0;i<pointingD0.size();++i) s += Form("%s%.2f", i?", ":" ", pointingD0[i]);
                s += " m";
                note->DrawLatex(0.05,0.10,s);
            }
            garbage.push_back(note);

            c->Print(pdf);
            ++drawn;
            for (auto* o : garbage) delete o;
        }

        c->Print(pdf + "]");   // close multipage
        std::cout << "[caloShowerEventDisplays] wrote " << drawn << " pages -> " << pdf << "\n";
        delete c;
    }

    fin->Close();
}
