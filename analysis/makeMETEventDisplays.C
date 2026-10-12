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
// analysisHelperFunctions.h declares DrawATLASLabel and calls it from several of its plotting
// helpers, but deliberately leaves the definition to each analysis macro (see the note at its
// declaration). This macro draws its own page labels and never calls it, but the symbol still
// has to exist or the translation unit will not link.
//
// Kept free of the gPileup / gProcLabel globals that metAnalysisAndRates.C's copy reads: those
// live in that macro, not in the header, and nothing here sets them.
void DrawATLASLabel(double x = 0.20, double /*y*/ = 0.88, const char* status = "Work in progress") {
    if (gPad) {
        gPad->SetTopMargin(0.14);   // frame top now ~0.86, leaving a white strip above
        gPad->Modified();
        gPad->Update();
    }
    const double yAtlas = 0.945;    // "ATLAS <status>" line, in the strip above the frame
    const double yInfo  = 0.895;    // beam-energy line, just below it
    TLatex l; l.SetNDC(); l.SetTextFont(72); l.SetTextColor(kBlack); l.SetTextSize(0.04);
    l.DrawLatex(x, yAtlas, "ATLAS");
    TLatex p; p.SetNDC(); p.SetTextFont(42); p.SetTextColor(kBlack); p.SetTextSize(0.04);
    p.DrawLatex(x + 0.13, yAtlas, status);
    TLatex e; e.SetNDC(); e.SetTextFont(42); e.SetTextColor(kBlack); e.SetTextSize(0.035);
    e.DrawLatex(x, yInfo, "#sqrt{s} = 14 TeV");
}

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
// One eta ring. The phi footprints are stored per tower rather than as a single pitch,
// because a ring is not always regular — the forward rings mix pitches — and a per-tower
// edge is the only thing that lets the ring tile 2 pi exactly.
struct FexTowerRing {
    double etaCentre = 0.0;
    double etaLow    = 0.0;
    double etaHigh   = 0.0;
    double etaWidth  = 0.0;
    double phiWidth  = 0.0;   // median phi pitch, for reporting only
    int    nPhi      = 0;
    int    nTowers   = 0;
    std::vector<double> phiCentres;  // ascending
    std::vector<double> phiLow;      // per-tower phi edges; may run outside +-pi (they wrap)
    std::vector<double> phiHigh;
};

// Full geometry of one tower collection: the per-eta-ring info plus the variable bin
// edges of the eta-phi histogram used to display it.
//
// The bin edges matter as much as the ring widths. Painting the towers onto a *uniform*
// fine grid aliases: the tower pitch is not a multiple of the fine bin width (gFEX is
// 2 pi / 32 = 0.19635 in phi, not 0.2, with edges at -pi + k * 0.19635), so
// FindBin(lowEdge) .. FindBin(highEdge) covers 4 fine bins on some rows and only 3 on
// others and the display comes out striped in phi — and likewise in eta wherever the
// forward granularity stops being a multiple of the fine bin width.
// Binning the histogram on the tower edges themselves removes the aliasing by
// construction: every tower covers an exact whole number of bins.
struct FexTowerGeometry {
    std::map<double, FexTowerRing> rings;     // keyed by the rounded eta centre
    std::vector<double>            etaEdges;  // tower boundaries in eta, ascending
    std::vector<double>            phiEdges;  // tower boundaries in phi, ascending
    bool                           built = false;  // derivation has been run (not "it succeeded")
    bool usable() const { return etaEdges.size() > 1 && phiEdges.size() > 1; }
};

// Round to 1e-3 so float-precision jitter in the stored centres does not split rings.
static double fexTowerKey(double v) { return std::round(v * 1000.0) / 1000.0; }

// Fold into [-pi, pi], leaving the two endpoints themselves untouched.
static double fexWrapPhi(double p) {
    while (p >  M_PI) p -= 2.0 * M_PI;
    while (p < -M_PI) p += 2.0 * M_PI;
    return p;
}

// Sort a raw edge list and merge duplicates. The tolerance is far below the finest real
// pitch (jFEX phi is 2 pi / 64 = 0.098) but well above the 1e-3 rounding of fexTowerKey,
// so it collapses jitter without merging genuinely distinct tower boundaries.
static void fexFinaliseEdges(std::vector<double>& edges, double tol = 2e-3) {
    std::sort(edges.begin(), edges.end());
    std::vector<double> out;
    out.reserve(edges.size());
    for (double e : edges)
        if (out.empty() || e - out.back() > tol) out.push_back(e);
    edges.swap(out);
}

// Collapse centres that sit closer together than tol onto their cluster mean. No real tower
// is narrower than ~0.05 in either eta or phi, so a 0.01 threshold absorbs float jitter and
// duplicated entries without merging genuinely distinct towers. Without this, two copies of
// one tower stored 1e-3 apart look like two rings (or like a 1e-3 phi pitch), which is enough
// to wreck every width derived below.
static std::vector<double> fexClusterCentres(const std::set<double>& raw, double tol = 0.01) {
    std::vector<double> centres;
    std::vector<double> run;
    for (double v : raw) {
        if (!run.empty() && v - run.back() > tol) {
            double s = 0.0;
            for (double x : run) s += x;
            centres.push_back(s / run.size());
            run.clear();
        }
        run.push_back(v);
    }
    if (!run.empty()) {
        double s = 0.0;
        for (double x : run) s += x;
        centres.push_back(s / run.size());
    }
    return centres;
}

// Snap a value onto the nearest cluster centre, or return it unchanged if nothing is within tol.
static double fexSnap(const std::vector<double>& centres, double v, double tol = 0.01) {
    if (centres.empty()) return v;
    auto it = std::lower_bound(centres.begin(), centres.end(), v);
    double best = (it == centres.end()) ? centres.back() : *it;
    if (it != centres.begin()) {
        double prev = *std::prev(it);
        if (std::abs(prev - v) < std::abs(best - v)) best = prev;
    }
    return (std::abs(best - v) <= tol) ? best : v;
}

// Median of the neighbour gaps in a sorted set of centres. The median is used rather than the
// minimum because it survives both failure modes a real collection shows: duplicated towers
// (which make a spuriously tiny gap) and missing towers (which make a doubled one). Taking the
// minimum, as this used to, lets a single duplicate collapse the whole ring's pitch.
static double fexMedianGap(const std::set<double>& centres) {
    std::vector<double> gaps;
    for (auto it = centres.begin(); it != centres.end(); ++it) {
        auto next = std::next(it);
        if (next == centres.end()) break;
        double g = *next - *it;
        if (g > 1e-6) gaps.push_back(g);
    }
    if (gaps.empty()) return 0.0;
    std::sort(gaps.begin(), gaps.end());
    return gaps.at(gaps.size() / 2);
}

// Index of the phi cell of this ring that contains phi, else the nearest centre.
// Containment first so it does not matter whether the stored phi is a cell centre or an edge.
// A linear scan is fine: a ring has at most a few tens of towers.
static int fexPhiIndex(const FexTowerRing& r, double phi) {
    const double eps = 1e-6;
    for (int i = 0; i < (int)r.phiLow.size(); ++i) {
        // The stored edges can run outside +-pi, so compare in the cell's own frame.
        double d = fexWrapPhi(phi - 0.5 * (r.phiLow.at(i) + r.phiHigh.at(i)));
        if (std::abs(d) <= 0.5 * (r.phiHigh.at(i) - r.phiLow.at(i)) + eps) return i;
    }
    int best = -1;
    double bestD = 0.0;
    for (int i = 0; i < (int)r.phiCentres.size(); ++i) {
        double d = std::abs(fexWrapPhi(phi - r.phiCentres.at(i)));
        if (best < 0 || d < bestD) { best = i; bestD = d; }
    }
    return best;
}

// Ring whose eta footprint contains eta, else the nearest one within half a width.
// Containment first, for the same reason as above, and because with an explicit grid the
// ring centres are fixed while the stored tower eta may be quoted either way.
static std::map<double, FexTowerRing>::const_iterator
fexFindRing(const std::map<double, FexTowerRing>& rings, double eta) {
    if (rings.empty()) return rings.end();
    const double eps = 1e-6;
    for (auto it = rings.begin(); it != rings.end(); ++it)
        if (eta >= it->second.etaLow - eps && eta < it->second.etaHigh + eps) return it;

    auto it = rings.lower_bound(eta);
    auto best = it;
    if (it == rings.end()) best = std::prev(rings.end());
    else if (it != rings.begin()) {
        auto prev = std::prev(it);
        if (std::abs(prev->first - eta) < std::abs(it->first - eta)) best = prev;
    }
    double tol = 0.5 * best->second.etaWidth + 0.01;
    return (std::abs(best->first - eta) <= tol) ? best : rings.end();
}

// -----------------------------------------------------------------------
// Ring-by-ring table of a tower geometry, plus the granularity regions and the resulting
// display binning. Works for both the explicit grid and a derived one.
static void printFexTowerGeometry(const FexTowerGeometry& geom, const std::string& label) {
    std::cout << "\n=== " << label << " tower geometry ===\n";
    std::cout << "  " << geom.rings.size() << " eta rings\n";
    std::cout << "   eta centre   eta low   eta high     d(eta)   N(phi)     d(phi)   N(towers)\n";
    for (const auto& kv : geom.rings) {
        const FexTowerRing& r = kv.second;
        printf("  %10.3f %9.3f %10.3f %10.3f %8d %10.3f %11d\n",
               r.etaCentre, r.etaLow, r.etaHigh, r.etaWidth, r.nPhi, r.phiWidth, r.nTowers);
    }
    // Group the rings by (d(eta), d(phi)) so the granularity regions are easy to read off.
    std::map<std::pair<double,double>, std::pair<double,double>> regions;  // (dEta,dPhi) -> (min|eta|, max|eta|)
    for (const auto& kv : geom.rings) {
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
    int nCells = 0;
    for (const auto& kv : geom.rings) nCells += kv.second.nPhi;
    printf("  %d grid cells total\n", nCells);
    if (geom.usable())
        printf("  display binning: %d eta bins over [%.3f, %.3f], %d phi bins over [%.3f, %.3f]\n",
               (int)geom.etaEdges.size() - 1, geom.etaEdges.front(), geom.etaEdges.back(),
               (int)geom.phiEdges.size() - 1, geom.phiEdges.front(), geom.phiEdges.back());
    std::cout << std::endl;
}

// -----------------------------------------------------------------------
// The gFEX tower grid, transcribed from LVL1::gTower::getEtaPhi — the function that assigns
// every gTower its eta and phi, and so the only authoritative source for it.
//
//   s_EtaCenter[40]  the eta centre of each global iEta index, reproduced verbatim below
//   phi              16 bins for iEta 0-3 and 36-39, 32 bins otherwise, centred at
//                    (iPhi + 0.5) * (2 pi / nPhi) and folded into [-pi, pi)
//
// 8 forward rings x 16 + 32 central rings x 32 = 128 + 1024 = 1152, matching reserve(1152) in
// gTowerBuilder and the "Expected towers 1152" in gTowerMakerFromGfexTowers.
//
// Everything before this was guesswork off gTowerBuilder's region comments, which are stale
// (BuildFCALgTowers is labelled "pure EMEC region") and gave bins where no tower lives, leaving
// whole eta rings empty in every event. Two things are only visible from the centre list:
//
//   * the forward/central split is by INDEX, not by an |eta| threshold. iEta 3 sits at 3.18 with
//     16 phi bins while iEta 4 sits at 3.15 with 32 — so "nPhi = 16 when |eta| > 3.2" is wrong,
//     it would put 3.18 in the central class.
//   * those two centres are 0.03 apart, i.e. the FCAL and central towers genuinely OVERLAP in
//     eta. There is no contiguous tiling to find, which is why every hand-built grid failed here.
//
// Edges are midpoints between adjacent centres. That keeps the axis contiguous and renders the
// overlap as a pair of narrow bins rather than a gap.
static FexTowerGeometry buildGFexTowerGeometrySpec() {
    static const double etaCentre[] = {
        -4.5, -3.8, -3.38, -3.18, -3.15, -3.0,
        -2.8, -2.6, -2.35, -2.1, -1.9, -1.7, -1.5, -1.3, -1.1, -0.9,
        -0.7, -0.5, -0.3, -0.1,  0.1,  0.3,  0.5,  0.7,  0.9,  1.1,
         1.3,  1.5,  1.7,  1.9,  2.1,  2.35, 2.6,  2.8,  3.0,
         3.15, 3.18, 3.38, 3.8,  4.5};
    const int nRings = (int)(sizeof(etaCentre) / sizeof(double));

    FexTowerGeometry geom;
    for (int i = 0; i < nRings; ++i) {
        double c = etaCentre[i];
        double lo = (i == 0)          ? c - 0.5 * (etaCentre[1] - etaCentre[0])
                                      : 0.5 * (etaCentre[i-1] + c);
        double hi = (i == nRings - 1) ? c + 0.5 * (etaCentre[nRings-1] - etaCentre[nRings-2])
                                      : 0.5 * (c + etaCentre[i+1]);

        FexTowerRing r;
        r.etaCentre = c;
        r.etaLow    = lo;
        r.etaHigh   = hi;
        r.etaWidth  = hi - lo;

        // Forward classes are iEta 0-3 and 36-39, exactly as gTower::getEtaPhi tests them.
        int nPhi = (i <= 3 || i >= nRings - 4) ? 16 : 32;
        double dPhi = 2.0 * M_PI / nPhi;
        r.nPhi     = nPhi;
        r.nTowers  = nPhi;
        r.phiWidth = dPhi;

        // (iPhi + 0.5) * dPhi folded into [-pi, pi), then sorted: the fold leaves the raw order
        // running 0..pi then -pi..0, and everything downstream expects ascending centres.
        std::vector<double> ctrs;
        for (int k = 0; k < nPhi; ++k) ctrs.push_back(fexWrapPhi((k + 0.5) * dPhi));
        std::sort(ctrs.begin(), ctrs.end());
        for (double ctr : ctrs) {
            r.phiCentres.push_back(ctr);
            r.phiLow.push_back(ctr - 0.5 * dPhi);
            r.phiHigh.push_back(ctr + 0.5 * dPhi);
        }

        geom.rings[c] = r;
        geom.etaEdges.push_back(lo);
        geom.etaEdges.push_back(hi);
    }
    // The 16-bin edges are every other 32-bin edge, so the union is just the 32-bin grid and a
    // forward tower spans exactly two display bins.
    double dPhiFine = 2.0 * M_PI / 32;
    for (int k = 0; k <= 32; ++k) geom.phiEdges.push_back(-M_PI + k * dPhiFine);

    fexFinaliseEdges(geom.etaEdges);
    fexFinaliseEdges(geom.phiEdges);
    geom.built = true;
    return geom;
}

// Build (and optionally dump) the geometry of a tower collection.
// Rings are keyed by the rounded eta centre so per-tower lookup is a single map find.
static FexTowerGeometry buildFexTowerGeometry(const std::vector<double>& towerEta,
                                              const std::vector<double>& towerPhi,
                                              const std::string& label,
                                              bool debugPrint) {
    // Cluster the raw centres first, so jitter / duplicate entries cannot masquerade as extra
    // rings or as an absurdly fine phi pitch.
    std::set<double> rawEta, rawPhi;
    for (unsigned int i = 0; i < towerEta.size(); ++i) {
        rawEta.insert(fexTowerKey(towerEta.at(i)));
        rawPhi.insert(fexTowerKey(towerPhi.at(i)));
    }
    std::vector<double> etaCentres = fexClusterCentres(rawEta);
    std::vector<double> phiCentres = fexClusterCentres(rawPhi);

    std::map<double, std::set<double>> phisAtEta;
    std::map<double, int>              nTowersAtEta;
    for (unsigned int i = 0; i < towerEta.size(); ++i) {
        double kEta = fexSnap(etaCentres, fexTowerKey(towerEta.at(i)));
        phisAtEta[kEta].insert(fexSnap(phiCentres, fexTowerKey(towerPhi.at(i))));
        ++nTowersAtEta[kEta];
    }

    std::vector<double> etaKeys;
    etaKeys.reserve(phisAtEta.size());
    for (const auto& kv : phisAtEta) etaKeys.push_back(kv.first);

    FexTowerGeometry geom;
    if (etaKeys.empty()) { geom.built = true; return geom; }

    // ---- Ring eta edges: midpoints between adjacent ring centres ----
    // This used to walk outward from eta = 0 taking highEdge = 2*centre - lowEdge, which
    // recovers the exact width of every ring *provided* the rings really do tile contiguously
    // from an edge at zero. The forward rings do not, and because each width is derived from the
    // running edge, one bad ring corrupts every ring outside it — which is exactly why the
    // central region looked right while the FCal came out full of holes and one-bin slivers.
    //
    // The midpoint rule gives up exactness at a granularity transition (where a narrow ring abuts
    // a wide one the boundary lands at the average of the two centres, off by a fraction of a
    // tower) in exchange for two properties worth much more in a display: it is contiguous by
    // construction, so no gap or sliver can ever appear, and each ring's width depends only on
    // its immediate neighbours, so a local oddity stays local.
    const int nRings = (int)etaKeys.size();
    for (int i = 0; i < nRings; ++i) {
        double c = etaKeys.at(i);
        double low, high;
        if (nRings == 1) { low = c - 0.1; high = c + 0.1; }
        else {
            low  = (i == 0)          ? c - 0.5 * (etaKeys.at(1) - etaKeys.at(0))
                                     : 0.5 * (etaKeys.at(i-1) + c);
            high = (i == nRings - 1) ? c + 0.5 * (etaKeys.at(nRings-1) - etaKeys.at(nRings-2))
                                     : 0.5 * (c + etaKeys.at(i+1));
        }

        FexTowerRing r;
        r.etaCentre = c;
        r.etaLow    = low;
        r.etaHigh   = high;
        r.etaWidth  = high - low;

        const std::set<double>& phis = phisAtEta.at(c);
        r.nPhi    = (int)phis.size();
        r.nTowers = nTowersAtEta.at(c);
        r.phiCentres.assign(phis.begin(), phis.end());

        // Median pitch, reported in the geometry table so an irregular ring is visible.
        double phiGap = fexMedianGap(phis);
        if (phiGap <= 0.0) phiGap = r.nPhi > 0 ? 2.0 * M_PI / r.nPhi : 0.1;
        r.phiWidth = phiGap;

        // Per-tower phi edges, again by midpoints, wrapping circularly so the ring tiles 2 pi
        // exactly. A ring with a tower genuinely missing shows up as one wide neighbour rather
        // than a hole; the geometry dump flags those so they are not mistaken for real structure.
        int n = (int)r.phiCentres.size();
        r.phiLow.resize(n);
        r.phiHigh.resize(n);
        for (int k = 0; k < n; ++k) {
            if (n == 1) { r.phiLow.at(k) = r.phiCentres.at(k) - M_PI;
                          r.phiHigh.at(k) = r.phiCentres.at(k) + M_PI; continue; }
            double prev = (k == 0)     ? r.phiCentres.at(n-1) - 2.0 * M_PI : r.phiCentres.at(k-1);
            double next = (k == n - 1) ? r.phiCentres.at(0)   + 2.0 * M_PI : r.phiCentres.at(k+1);
            r.phiLow.at(k)  = 0.5 * (prev + r.phiCentres.at(k));
            r.phiHigh.at(k) = 0.5 * (r.phiCentres.at(k) + next);
        }
        geom.rings[c] = r;
    }

    // ---- Display bin edges: the union of every tower's own boundaries ----
    // Collecting the edges of each individual tower (rather than assuming one global pitch)
    // is what makes the display exact when the phi granularity changes with eta: whatever
    // the mix of pitches, every tower's low and high edge is itself a bin edge, so each
    // tower covers a whole number of bins and no gap can open up inside its footprint.
    // Where a coarse forward ring and a finer central ring interleave, the union simply
    // carries a few extra narrow bins — harmless, since the towers still tile them exactly.
    for (const auto& kv : geom.rings) {
        geom.etaEdges.push_back(kv.second.etaLow);
        geom.etaEdges.push_back(kv.second.etaHigh);
        for (unsigned int k = 0; k < kv.second.phiLow.size(); ++k) {
            geom.phiEdges.push_back(fexWrapPhi(kv.second.phiLow.at(k)));
            geom.phiEdges.push_back(fexWrapPhi(kv.second.phiHigh.at(k)));
        }
    }
    // A tower straddling +-pi contributes wrapped edges, so pin the axis to the full range.
    geom.phiEdges.push_back(-M_PI);
    geom.phiEdges.push_back( M_PI);
    fexFinaliseEdges(geom.etaEdges);
    fexFinaliseEdges(geom.phiEdges);
    geom.built = true;

    if (debugPrint) printFexTowerGeometry(geom, label + " (derived from tower centres)");

    return geom;
}

// -----------------------------------------------------------------------
// One tower of the display grid: the summed E_T of a cell and where it sits.
struct FexCell {
    double etaCentre = 0.0;
    int    phiIndex  = 0;
    double eta       = 0.0;
    double phi       = 0.0;
    double et        = 0.0;   // sum over the entries that landed in this cell
    int    entries   = 0;
};

// Map the collection's entries onto the grid cells, keeping the largest |E_T| where several
// land in the same cell.
//
// It does NOT sum them. Summing was tried and is wrong: the stored towers are finer in eta
// than this grid, so the number of entries a cell absorbs scales with the cell's width. The
// 0.3-wide bins either side of |eta| = 2.35 and the 0.45-0.5 wide FCal bins then come out
// systematically hotter than the 0.2-wide barrel, painting vertical bands that track bin width
// rather than energy. Keeping the largest at least shows a real tower's E_T; any cell holding
// more than one entry is reported below, because that is a grid mismatch worth fixing, not
// something to paper over.
static std::vector<FexCell> fexAccumulateCells(const FexTowerGeometry& geom,
                                               const std::vector<double>& tEt,
                                               const std::vector<double>& tEta,
                                               const std::vector<double>& tPhi,
                                               unsigned int* nUnplaced = nullptr) {
    std::map<std::pair<double,int>, FexCell> cells;
    unsigned int unplaced = 0;
    for (unsigned int i = 0; i < tEt.size(); ++i) {
        double eta = tEta.at(i), phi = fexWrapPhi(tPhi.at(i));
        auto itRing = fexFindRing(geom.rings, eta);
        if (itRing == geom.rings.end()) { ++unplaced; continue; }
        const FexTowerRing& r = itRing->second;
        int k = fexPhiIndex(r, phi);
        if (k < 0) { ++unplaced; continue; }
        auto key = std::make_pair(r.etaCentre, k);
        FexCell& c = cells[key];
        if (c.entries == 0) {
            c.etaCentre = r.etaCentre;
            c.phiIndex  = k;
            c.eta       = r.etaCentre;
            c.phi       = r.phiCentres.at(k);
            c.et        = tEt.at(i);
        } else if (std::fabs(tEt.at(i)) > std::fabs(c.et)) {
            c.et = tEt.at(i);
        }
        ++c.entries;
    }
    if (nUnplaced) *nUnplaced = unplaced;

    std::vector<FexCell> out;
    out.reserve(cells.size());
    for (const auto& kv : cells) out.push_back(kv.second);
    return out;
}

// The n highest-E_T cells of an event, most energetic first.
static std::vector<FexCell> fexTopCells(const std::vector<FexCell>& cells, unsigned int n) {
    std::vector<FexCell> sorted = cells;
    std::sort(sorted.begin(), sorted.end(),
              [](const FexCell& a, const FexCell& b) { return a.et > b.et; });
    if (sorted.size() > n) sorted.resize(n);
    return sorted;
}

// How many entries landed in each grid cell, as a count per multiplicity.
//
// This is the measurement that decides whether summing is right at all. "1 x N" for every cell
// means the collection is already one object per geometric tower and the accumulation below is
// a no-op. Anything else means several objects share a tower, and only then does the question
// of summing versus picking one arise.
static void printFexCellMultiplicity(const std::string& label, int iEvt,
                                     const std::vector<FexCell>& cells, size_t nEntries) {
    std::map<int, int> byMult;
    for (const FexCell& c : cells) ++byMult[c.entries];
    printf("  %s event %d: %u entries -> %u cells; entries/cell:",
           label.c_str(), iEvt, (unsigned int)nEntries, (unsigned int)cells.size());
    for (const auto& kv : byMult) printf("  %d x %d", kv.first, kv.second);
    printf("\n");
}

// Compare the eta values the collection actually stores against the display grid.
//
// A cell holding more than one entry means the grid is coarser than the data, and that shows up
// as banding that follows bin width. This prints the container's own distinct eta values (and
// how many towers sit at each) next to the grid's bins, which is what identifies the real
// segmentation instead of inferring it from tower counts.
static void printFexEtaGranularity(const std::string& label,
                                   const std::vector<double>& tEta,
                                   const FexTowerGeometry& geom) {
    std::set<double> raw;
    for (double e : tEta) raw.insert(fexTowerKey(e));
    std::vector<double> centres = fexClusterCentres(raw);

    std::map<double, int> countAt;
    for (double e : tEta) countAt[fexSnap(centres, fexTowerKey(e))]++;

    std::cout << "\n=== " << label << " stored eta values vs display grid ===\n";
    printf("  container holds %u distinct eta values; grid has %u bins\n",
           (unsigned int)centres.size(), (unsigned int)geom.rings.size());
    printf("     eta      towers    gap to next   grid bin [low, high]   d(eta)\n");
    for (unsigned int i = 0; i < centres.size(); ++i) {
        double gap = (i + 1 < centres.size()) ? centres.at(i+1) - centres.at(i) : 0.0;
        auto it = fexFindRing(geom.rings, centres.at(i));
        if (it != geom.rings.end())
            printf("  %8.3f %10d %14.3f   [%7.3f, %7.3f] %8.3f\n",
                   centres.at(i), countAt[centres.at(i)], gap,
                   it->second.etaLow, it->second.etaHigh, it->second.etaWidth);
        else
            printf("  %8.3f %10d %14.3f   (outside the grid)\n",
                   centres.at(i), countAt[centres.at(i)], gap);
    }
    std::cout << std::endl;
}

static void printFexTopCells(const std::string& label, int iEvt, const std::vector<FexCell>& hot) {
    printf("  %s event %d: %u highest-E_T towers\n",
           label.c_str(), iEvt, (unsigned int)hot.size());
    printf("       rank     E_T [GeV]        eta        phi   entries\n");
    for (unsigned int i = 0; i < hot.size(); ++i)
        printf("     %6u %13.4f %10.3f %10.3f %9d\n",
               i + 1, hot.at(i).et, hot.at(i).eta, hot.at(i).phi, hot.at(i).entries);
}

// -----------------------------------------------------------------------
// Dump every tower of one event: position, footprint and E_T. Also reports the E_T
// range, which is the quickest way to spot a unit problem (counts vs MeV vs GeV) in
// the input branch.
static void dumpFexTowerPositions(const std::string& label, int iEvt,
                                  const std::vector<double>& tEt,
                                  const std::vector<double>& tEta,
                                  const std::vector<double>& tPhi,
                                  const FexTowerGeometry& geom) {
    std::cout << "\n=== " << label << " tower positions, event " << iEvt
              << " (" << tEt.size() << " towers) ===\n";
    std::cout << "       i        eta        phi     E_T [GeV]    d(eta)    d(phi)"
                 "        eta range             phi range\n";
    double minEt = 0.0, maxEt = 0.0, sumEt = 0.0;
    for (unsigned int i = 0; i < tEt.size(); ++i) {
        double eta = tEta.at(i), phi = tPhi.at(i), et = tEt.at(i);
        auto it = fexFindRing(geom.rings, eta);
        double dEta = (it != geom.rings.end()) ? it->second.etaWidth : 0.0;
        double dPhi = (it != geom.rings.end()) ? it->second.phiWidth : 0.0;
        double eLo  = (it != geom.rings.end()) ? it->second.etaLow  : eta;
        double eHi  = (it != geom.rings.end()) ? it->second.etaHigh : eta;
        printf("  %6u %10.4f %10.4f %12.4f %9.3f %9.3f  [%7.3f, %7.3f]  [%7.3f, %7.3f]\n",
               i, eta, phi, et, dEta, dPhi, eLo, eHi, phi - 0.5*dPhi, phi + 0.5*dPhi);
        if (i == 0 || et < minEt) minEt = et;
        if (i == 0 || et > maxEt) maxEt = et;
        sumEt += et;
    }
    printf("  %u towers, E_T range [%.4f, %.4f] GeV, scalar sum = %.2f GeV\n",
           (unsigned int)tEt.size(), minEt, maxEt, sumEt);
    // The smallest non-zero |E_T| is the collection's LSB: if it comes out at 0.2 (or some other
    // fixed step) rather than a continuous spread, the branch is in counts and the /1000 applied
    // in HERNTupler is converting the wrong unit.
    double minAbsNonZero = 0.0;
    for (double e : tEt)
        if (e != 0.0 && (minAbsNonZero == 0.0 || std::abs(e) < minAbsNonZero)) minAbsNonZero = std::abs(e);
    printf("  smallest non-zero |E_T| (apparent LSB): %.6f GeV\n", minAbsNonZero);
    std::cout << std::endl;
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
                              bool   drawGFexTowers          = false, // separate PDF of gFEX emulated tower eta-phi displays (same events as the MET displays)
                              bool   drawJFexTowers          = false, // separate PDF of jFEX emulated tower eta-phi displays (same events as the MET displays)
                              double fexTowerZMin            = 0.0,  // z-axis floor [GeV] for the g/jFEX tower displays; >= 0 means auto (most negative tower, or 0). Every tower is drawn, no E_T threshold
                              double fexTowerZMax            = 0.0,  // z-axis ceiling [GeV]; <= 0 means auto (largest tower E_T in the event)
                              int    fexTowerDumpEvent       = -1    // dump every tower position/E_T: -1 = the first displayed event, >= 0 = that iEvt, -2 = never
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

    // g/jFEX tower displays live in their own PDFs and draw the same events as
    // METEventDisplays.pdf: they are made inside the accepted-event block, after the
    // background MET over-reconstruction filter, so page N of each PDF is the same event.
    // They stop early (maxFexTowerDisplays below) to keep the files openable, so they hold
    // the first N accepted events rather than all of them.
    TString pdf_gTow = outDir + "gFexTowerEventDisplays.pdf";
    TString pdf_jTow = outDir + "jFexTowerEventDisplays.pdf";
    TCanvas cFexTower("cFexTower", "cFexTower", 800, 700);
    if (drawGFexTowers) cFexTower.Print(pdf_gTow + "(");
    if (drawJFexTowers) cFexTower.Print(pdf_jTow + "(");

    // Set to true to dump the derived tower geometry (eta ring edges, d(eta), d(phi),
    // granularity regions) for the first event with towers in each collection.
    const bool debugFexTowerGeometry = true;

    // gFEX uses the exact grid from LVL1::gTower::getEtaPhi. jFEX has no transcription here, so
    // its geometry is still derived from the first event with towers and cached.
    FexTowerGeometry gFexTowerGeom = buildGFexTowerGeometrySpec();
    FexTowerGeometry jFexTowerGeom;
    if (drawGFexTowers && debugFexTowerGeometry)
        printFexTowerGeometry(gFexTowerGeom, "gFEX (LVL1::gTower::getEtaPhi)");

    // One-shot per-tower dump (see fexTowerDumpEvent), tracked per collection.
    bool gFexDumpDone = (fexTowerDumpEvent == -2);
    bool jFexDumpDone = (fexTowerDumpEvent == -2);

    // Leading gFEX tower per event, over every accepted event rather than just the displayed
    // ones. Written to its own directory so it is not mixed in with the event displays.
    std::vector<double> gFexLeadingTowerEt, gFexLeadingTowerEta;

    // One-shot comparison of the stored eta values against the display grid.
    bool gFexEtaGranularityDone = false;

    const unsigned int maxDisplays = 200;
    // Far fewer tower pages than MET pages: a tower display paints every tower in the
    // calorimeter, so the PDF gets unmanageably large well before 200 events.
    const unsigned int maxFexTowerDisplays = 10;
    unsigned int nAccepted = 0, nFexTowerDisplays = 0;
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

        // Filled by the tower block further down, and read back by the text page so the
        // hottest towers appear under the MET table. Declared here because the text-page
        // lambda captures it by reference and is defined before that block runs.
        std::vector<FexCell> gFexTopTen;

        // --- Overlays shared by every eta-phi display (GEP towers, gFEX towers, jFEX towers):
        //     jet circles plus a horizontal dashed line at phi_MET for each algorithm. ---
        // etaLo/etaHi set how far the phi_MET lines run: the GEP display is a fixed +-5 frame,
        // while the g/jFEX tower displays are binned on the real tower edges and so span
        // whatever the collection actually covers.
        // drawMETLines is off for the tower displays: eight coloured lines plus their legend sit
        // on top of the towers and hide the structure the plot exists to show. Those numbers are
        // on the MET table page that follows each tower display instead.
        auto drawEtaPhiOverlays = [&](double etaLo = -5.0, double etaHi = 5.0,
                                      bool drawMETLines = true) {
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
            if (!drawMETLines) return;
            TLegend* legMET = new TLegend(0.12, 0.74, 0.42, 0.89);
            legMET->SetBorderSize(0); legMET->SetFillStyle(0); legMET->SetTextSize(0.022);
            for (const auto& m : mets) {
                if (!m.hasXY) continue;
                if (m.METx == 0.0 && m.METy == 0.0) continue;
                double phi = std::atan2(m.METy, m.METx);
                TLine* l = new TLine(etaLo, phi, etaHi, phi);
                l->SetLineColor(m.color); l->SetLineWidth(2); l->SetLineStyle(m.lineStyle);
                l->Draw("same");
                legMET->AddEntry(l, Form("%s (#varphi=%.2f)", m.name, phi), "l");
            }
            legMET->Draw();
        };

        // --- Text page: per-algorithm MET table plus the event-level info. Shared by
        //     METEventDisplays.pdf and the g/jFEX tower PDFs, so every display is followed
        //     by the numbers for that event. ---
        auto drawMETTextPage = [&](TCanvas& can, const TString& pdfName) {
            can.cd();
            can.Clear();

            TLatex lat;
            lat.SetTextAlign(13);

            lat.SetTextSize(0.030);
            lat.DrawLatexNDC(0.05, 0.95, Form("Event %d   (signal=%s)", iEvt, signalBool ? "true" : "false"));

            // The table is drawn one cell at a time at fixed NDC positions, with the numbers
            // right-aligned on their column. Formatting a whole row into one padded string does
            // not line anything up: TLatex renders in a proportional font, where a space is
            // narrower than a digit and every row ends up ragged.
            const double xName  = 0.04;   // left-aligned
            const double xMET   = 0.37;   // these are the RIGHT edge of each numeric column
            const double xMETx  = 0.51;
            const double xMETy  = 0.65;
            const double xSumET = 0.81;
            const double xPhi   = 0.95;

            auto drawTableRow = [&](double y, const char* name, const char* met, const char* metx,
                                    const char* mety, const char* sumet, const char* phi) {
                lat.SetTextAlign(13);                 // left, top
                lat.DrawLatexNDC(xName, y, name);
                lat.SetTextAlign(33);                 // right, top
                lat.DrawLatexNDC(xMET,   y, met);
                lat.DrawLatexNDC(xMETx,  y, metx);
                lat.DrawLatexNDC(xMETy,  y, mety);
                lat.DrawLatexNDC(xSumET, y, sumet);
                lat.DrawLatexNDC(xPhi,   y, phi);
            };

            lat.SetTextSize(0.026);
            double yText = 0.86;
            drawTableRow(yText, "MET algorithm", "MET [GeV]", "METx [GeV]",
                         "METy [GeV]", "SumET [GeV]", "#varphi [rad]");
            yText -= 0.035;
            TLine rule;
            rule.SetLineColor(kGray+1);
            rule.DrawLineNDC(xName, yText, xPhi, yText);
            yText -= 0.03;

            // phi is listed here because the tower displays no longer carry the phi_MET lines.
            for (const auto& m : mets) {
                bool hasPhi = m.hasXY && !(m.METx == 0.0 && m.METy == 0.0);
                std::string sMETx  = m.hasXY    ? Form("%.1f", m.METx)  : std::string("---");
                std::string sMETy  = m.hasXY    ? Form("%.1f", m.METy)  : std::string("---");
                std::string sSumET = m.hasSumET ? Form("%.1f", m.SumET) : std::string("---");
                std::string sPhi   = hasPhi ? Form("%.2f", std::atan2(m.METy, m.METx))
                                            : std::string("---");
                lat.SetTextColor(m.color);
                drawTableRow(yText, m.name, Form("%.1f", m.MET), sMETx.c_str(),
                             sMETy.c_str(), sSumET.c_str(), sPhi.c_str());
                yText -= 0.045;
            }
            lat.SetTextColor(kBlack);
            lat.SetTextAlign(13);

            // Ten hottest gFEX towers of this event, E_T only — their eta/phi are in the stdout
            // dump. Sits under the table so the tower pages carry the numbers behind the display.
            if (!gFexTopTen.empty()) {
                yText -= 0.015;
                lat.SetTextSize(0.024);
                lat.DrawLatexNDC(xName, yText, "10 highest gFEX tower E_{T} [GeV]:");
                yText -= 0.042;
                TString rowA, rowB;
                for (unsigned int i = 0; i < gFexTopTen.size(); ++i) {
                    TString v = Form("%.2f", gFexTopTen.at(i).et);
                    TString& row = (i < 5) ? rowA : rowB;
                    if (row.Length()) row += ",   ";
                    row += v;
                }
                lat.DrawLatexNDC(xName + 0.04, yText, rowA);
                if (rowB.Length()) {
                    yText -= 0.038;
                    lat.DrawLatexNDC(xName + 0.04, yText, rowB);
                }
                lat.SetTextSize(0.026);
            }

            // Bunch-train position: distance in BCID from the tail and the head (front) of the
            // train, i.e. how much out-of-time pileup this event sees from either side.
            yText -= 0.03;
            lat.SetTextSize(0.028);
            if (hasBunchTrain) {
                std::string sTail = distTailBunchTrainValues  >= 0 ? Form("%d", distTailBunchTrainValues)  : std::string("n/a");
                std::string sHead = distFrontBunchTrainValues >= 0 ? Form("%d", distFrontBunchTrainValues) : std::string("n/a");
                lat.DrawLatexNDC(0.05, yText,
                    Form("Position in bunch train [BCID] (distance from tail, head): (%s, %s)",
                         sTail.c_str(), sHead.c_str()));
            } else {
                lat.DrawLatexNDC(0.05, yText,
                    "Position in bunch train [BCID] (distance from tail, head): (n/a, n/a)");
            }
            yText -= 0.02;

            // Background extras: rate contribution + JZ slice + event number
            if (!signalBool) {
                yText -= 0.03;
                lat.SetTextSize(0.028);
                lat.DrawLatexNDC(0.05, yText,
                    Form("Event Weight (Rate Contribution): %.4g",
                         eventWeightsValues && !eventWeightsValues->empty() ? eventWeightsValues->at(0) : 0.0));
                yText -= 0.05;
                lat.DrawLatexNDC(0.05, yText,
                    Form("JZ Slice: %d, Event Number: %d", sampleJZSliceValues, iEvt));
            }

            // Emu config tag at bottom
            lat.SetTextSize(0.020);
            lat.SetTextColor(kGray+2);
            lat.DrawLatexNDC(0.05, 0.05, Form("emu config: %s", emuTag.c_str()));
            lat.SetTextColor(kBlack);

            can.Print(pdfName);
        };

        // --- Per-event tower E_T distributions: gFEX, jFEX and the EtaSK global towers on one
        //     plot. Puts the three collections on a common axis, which is the quickest way to
        //     see whether they share an energy scale or one of them is in different units. ---
        auto drawTowerEtDistributions = [&](TCanvas& can, const TString& pdfName) {
            // Empty towers are dropped from all three collections. Each is a fixed-size grid
            // whose unused cells sit at exactly zero; counting them buries the real distribution
            // under a single spike at the origin and inflates the tower count with cells that
            // hold nothing.
            // setsRange marks the collections the x axis is sized on: gFEX and jFEX, whichever
            // reaches further. The global EtaSK towers are still drawn but do not stretch the
            // axis, so their tail cannot squash the two FEX distributions being compared.
            struct EtColl { const char* name; const std::vector<double>* et; Color_t col; bool setsRange; };
            std::vector<EtColl> colls;
            if (gFexTowerEtValues && !gFexTowerEtValues->empty())
                colls.push_back({"gFEX towers",           gFexTowerEtValues, kRed+1,   true});
            if (jFexTowerEtValues && !jFexTowerEtValues->empty())
                colls.push_back({"jFEX towers",           jFexTowerEtValues, kAzure+2, true});
            if (towerEtValues && !towerEtValues->empty())
                colls.push_back({"Global towers (EtaSK)", towerEtValues,     kGreen+2, false});
            if (colls.empty()) return;

            double lo = 0.0, hi = 0.0;
            bool first = true;
            for (const auto& c : colls) {
                if (!c.setsRange) continue;
                for (double e : *c.et) {
                    if (e == 0.0) continue;
                    if (first) { lo = hi = e; first = false; }
                    else { lo = std::min(lo, e); hi = std::max(hi, e); }
                }
            }
            // Nothing in the FEX collections (or neither is present): fall back to everything.
            if (first) {
                for (const auto& c : colls) {
                    for (double e : *c.et) {
                        if (e == 0.0) continue;
                        if (first) { lo = hi = e; first = false; }
                        else { lo = std::min(lo, e); hi = std::max(hi, e); }
                    }
                }
            }
            if (first) return;
            double span = hi - lo;
            if (span <= 0.0) span = 1.0;
            lo -= 0.05 * span;
            hi += 0.05 * span;

            can.cd();
            can.Clear();
            can.SetLogy(true);

            // Left-hand side: the entries carry tower counts and sums, so the box is wide and
            // runs off the right edge of the pad if it is anchored there.
            TLegend* leg = new TLegend(0.18, 0.70, 0.63, 0.90);
            leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.023);

            std::vector<TH1F*> hEts;
            double yMax = 1.0;
            for (unsigned int iC = 0; iC < colls.size(); ++iC) {
                TH1F* h = new TH1F(Form("hTowEt_%d_%u", iEvt, iC),
                                   Form("Tower E_{T} distributions, event %d;E_{T} [GeV];Towers", iEvt),
                                   100, lo, hi);
                h->SetStats(0);
                h->SetLineColor(colls.at(iC).col);
                h->SetLineWidth(2);
                double sumEt = 0.0;
                int nUsed = 0;
                for (double e : *colls.at(iC).et) {
                    if (e == 0.0) continue;
                    h->Fill(e);
                    sumEt += e;
                    ++nUsed;
                }
                if (h->GetMaximum() > yMax) yMax = h->GetMaximum();
                leg->AddEntry(h, Form("%s: %d with E_{T} #neq 0, #Sigma E_{T} = %.1f GeV",
                                      colls.at(iC).name, nUsed, sumEt), "l");
                hEts.push_back(h);
            }

            hEts.front()->SetMinimum(0.5);
            hEts.front()->SetMaximum(yMax * 5.0);
            for (unsigned int iC = 0; iC < hEts.size(); ++iC)
                hEts.at(iC)->Draw(iC == 0 ? "HIST" : "HIST SAME");
            leg->Draw();

            // Units cross-check for the gFEX towers.
            //
            // towerEt() is in MeV (confirmed), so the /1000 HERNTupler applies to it is correct,
            // matching the jtowerEtMeV() used for jFEX. Kept as a standing check rather than a
            // one-off: gSumET comes from the xAOD gFexGlobalRoI SumEt(), converted on a separate
            // path, so a ratio that drifts from ~1 means the two have diverged.
            // (For reference, LVL1::gTower::setET divides by 200 to reach the gFEX digit scale —
            // that is the internal firmware representation, not what towerEt() hands back.)
            if (gFexTowerEtValues && !gFexTowerEtValues->empty() && sig_gSumET > 0.0) {
                double sumTow = 0.0;
                for (double e : *gFexTowerEtValues) sumTow += e;
                double ratio = sumTow / sig_gSumET;
                TLatex chk;
                chk.SetTextSize(0.022);
                chk.SetTextColor(kGray+3);
                chk.DrawLatexNDC(0.18, 0.66,
                    Form("#Sigma E_{T}(gFEX towers) = %.1f GeV,  gSumET (HER) = %.1f GeV,  ratio = %.3f",
                         sumTow, sig_gSumET, ratio));
                printf("  gFEX scale check, event %d: sum tower E_T = %.1f GeV, gSumET = %.1f GeV, "
                       "ratio = %.4f\n", iEvt, sumTow, sig_gSumET, ratio);
            }

            can.Print(pdfName);

            // Clear the pad before the histograms go away, then restore the linear scale for
            // the eta-phi pages that follow.
            can.Clear();
            can.SetLogy(false);
            for (TH1F* h : hEts) delete h;
        };

        // --- One eta-phi display of a g/jFEX tower collection, printed to its own PDF ---
        auto drawFexTowerDisplay = [&](const std::string& label,
                                       const std::vector<double>* tEt,
                                       const std::vector<double>* tEta,
                                       const std::vector<double>* tPhi,
                                       const std::vector<FexCell>& cells,
                                       FexTowerGeometry& geom,
                                       bool& dumpDone,
                                       const char* tag,
                                       const TString& pdfName) {
            if (!tEt || tEt->empty() || cells.empty()) return;
            if (!geom.usable()) return;

            // One-shot dump of every tower position / footprint / E_T.
            if (!dumpDone && (fexTowerDumpEvent < 0 || fexTowerDumpEvent == iEvt)) {
                dumpFexTowerPositions(label, iEvt, *tEt, *tEta, *tPhi, geom);
                dumpDone = true;
            }

            cFexTower.cd();
            cFexTower.Clear();

            // Bin on the real tower boundaries (see FexTowerGeometry): one bin per tower in the
            // uniform regions, plus the extra narrow bins where rings of different phi pitch
            // interleave. Nothing is resampled, so no tower can be clipped or split by the grid.
            const int nEtaBins = (int)geom.etaEdges.size() - 1;
            const int nPhiBins = (int)geom.phiEdges.size() - 1;
            TH2F* hTow = new TH2F(Form("h%sTowers_%d", tag, iEvt),
                                  Form("%s emulated towers;#eta;#varphi", label.c_str()),
                                  nEtaBins, geom.etaEdges.data(),
                                  nPhiBins, geom.phiEdges.data());
            hTow->SetStats(0);
            hTow->GetZaxis()->SetTitle("E_{T} [GeV]");

            // Track which bins a tower actually covers, separately from their content: E_T is
            // quantised, so a genuine zero-E_T tower is common and must not be confused with a
            // bin no tower reached.
            std::vector<char> covered((size_t)nEtaBins * nPhiBins, 0);
            const double eps = 1e-6;
            auto paintFootprint = [&](double eLo, double eHi, double pLo, double pHi, double et) {
                int bxLo = hTow->GetXaxis()->FindBin(eLo + eps);
                int bxHi = hTow->GetXaxis()->FindBin(eHi - eps);
                int byLo = hTow->GetYaxis()->FindBin(pLo + eps);
                int byHi = hTow->GetYaxis()->FindBin(pHi - eps);
                for (int bx = std::max(1, bxLo); bx <= std::min(nEtaBins, bxHi); ++bx) {
                    for (int by = std::max(1, byLo); by <= std::min(nPhiBins, byHi); ++by) {
                        size_t k = (size_t)(bx - 1) * nPhiBins + (by - 1);
                        hTow->SetBinContent(bx, by, et);
                        covered[k] = 1;
                    }
                }
            };

            // Cells are summed by the caller (fexAccumulateCells) so the same numbers feed the
            // display, the top-10 printout and the leading-tower distribution.
            // Every cell is drawn, with no E_T threshold: these towers legitimately carry negative
            // E_T (pedestal / noise subtraction can push one below zero) and that structure is the
            // point of the display.
            double maxEt = 0.0, minEt = 0.0;
            bool firstCell = true;
            for (const FexCell& c : cells) {
                auto itGeom = geom.rings.find(c.etaCentre);
                if (itGeom == geom.rings.end()) continue;
                const FexTowerRing& r = itGeom->second;
                int k = c.phiIndex;
                double et = c.et;
                if (firstCell) { minEt = maxEt = et; firstCell = false; }
                else { minEt = std::min(minEt, et); maxEt = std::max(maxEt, et); }

                // The cell's own edges, not centre +- half a nominal pitch: with the explicit
                // grid these are exact, and a forward cell is twice as wide as a central one.
                double pLo = r.phiLow.at(k), pHi = r.phiHigh.at(k);
                // A cell straddling +-pi is painted as its two pieces.
                if (pHi > M_PI) {
                    paintFootprint(r.etaLow, r.etaHigh, pLo, M_PI, et);
                    paintFootprint(r.etaLow, r.etaHigh, -M_PI, pHi - 2.0*M_PI, et);
                } else if (pLo < -M_PI) {
                    paintFootprint(r.etaLow, r.etaHigh, pLo + 2.0*M_PI, M_PI, et);
                    paintFootprint(r.etaLow, r.etaHigh, -M_PI, pHi, et);
                } else {
                    paintFootprint(r.etaLow, r.etaHigh, pLo, pHi, et);
                }
            }

            // Auto means the event's own full range, floor included: clamping the floor to 0 when
            // every tower sits above it (E_T in [1.05, 5.17], say) throws away a fifth of the
            // palette on values that never occur and flattens the contrast across the ones that do.
            double zMin = (fexTowerZMin <  0.0) ? fexTowerZMin : minEt;
            double zMax = (fexTowerZMax >  0.0) ? fexTowerZMax : maxEt;
            if (zMax <= zMin) zMax = zMin + 1.0;
            hTow->SetMinimum(zMin);
            hTow->SetMaximum(zMax);

            // COLZ skips a bin whose content is exactly zero, and draws nothing for one below the
            // axis minimum. Clamp the covered bins into range and nudge the exact zeros by an
            // amount far below one colour step, so every tower is painted and only the regions no
            // tower covers are left blank.
            const double zeroNudge = 1e-6 * std::max(1.0, zMax - zMin);
            size_t nCovered = 0;
            double coveredArea = 0.0, totalArea = 0.0;
            for (int bx = 1; bx <= nEtaBins; ++bx) {
                double wEta = hTow->GetXaxis()->GetBinWidth(bx);
                for (int by = 1; by <= nPhiBins; ++by) {
                    double area = wEta * hTow->GetYaxis()->GetBinWidth(by);
                    totalArea += area;
                    if (!covered[(size_t)(bx - 1) * nPhiBins + (by - 1)]) continue;
                    ++nCovered;
                    coveredArea += area;
                    double c = hTow->GetBinContent(bx, by);
                    if      (c < zMin) c = zMin;
                    else if (c > zMax) c = zMax;
                    if (c == 0.0) c = std::max(zMin, 0.0) + zeroNudge;
                    hTow->SetBinContent(bx, by, c);
                }
            }
            if (minEt < zMin)
                std::cout << "  " << label << " event " << iEvt << ": most negative tower "
                          << minEt << " GeV clamped to the z-axis floor " << zMin << " GeV\n";

            // Coverage is the decisive check on the derived footprints: the towers should tile
            // the eta-phi plane, so anything well below 100% means the widths are wrong and the
            // blank areas are holes rather than real structure.
            double coverFrac = (totalArea > 0.0) ? 100.0 * coveredArea / totalArea : 0.0;
            printf("  %s event %d: %d x %d display bins, %.1f%% of the eta-phi area covered\n",
                   label.c_str(), iEvt, nEtaBins, nPhiBins, coverFrac);

            hTow->Draw("COLZ");
            drawEtaPhiOverlays(geom.etaEdges.front(), geom.etaEdges.back(), /*drawMETLines=*/false);

            TLatex towLat;
            towLat.SetTextSize(0.026);
            towLat.DrawLatexNDC(0.10, 0.95,
                Form("Event %d   (%d %s towers, E_{T} #in [%.2f, %.2f] GeV, %.0f%% area covered)",
                     iEvt, (int)cells.size(), label.c_str(), minEt, maxEt, coverFrac));

            cFexTower.Print(pdfName);

            // Page 2: the same per-algorithm MET table that METEventDisplays.pdf carries, so each
            // tower display is followed by the numbers for that event. It clears the canvas first,
            // which is also what takes hTow off the pad before it is deleted.
            drawMETTextPage(cFexTower, pdfName);
            delete hTow;
        };

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
        // g/jFEX emulated towers. The per-cell sums, the top-10 printout and the leading-tower
        // distribution are done for EVERY accepted event; only the first maxFexTowerDisplays get
        // an eta-phi page, since those are what make the PDF large.
        // ===================================================================
        if (drawGFexTowers || drawJFexTowers) {
            // Load both collections up front so the E_T distribution page can overlay all three.
            if (drawGFexTowers) gFexTowersTree->GetEntry(iEvt);
            if (drawJFexTowers) jFexTowersTree->GetEntry(iEvt);

            int nGlobalNonZero = 0;
            if (towerEtValues)
                for (double e : *towerEtValues) if (e != 0.0) ++nGlobalNonZero;
            printf("Event %d tower counts: gFEX = %d, jFEX = %d, global EtaSK = %d (%d with E_T != 0)\n",
                   iEvt,
                   gFexTowerEtValues ? (int)gFexTowerEtValues->size() : -1,
                   jFexTowerEtValues ? (int)jFexTowerEtValues->size() : -1,
                   towerEtValues     ? (int)towerEtValues->size()     : -1,
                   nGlobalNonZero);

            // LVL1::gTowerBuilder reserves 1152 towers and gTowerMakerFromGfexTowers logs
            // "Expected towers 1152" for an xAOD::gFexTowerContainer, which is what
            // L1_gFexEmulatedTowers is. Anything else means the branch is not one entry per gFEX
            // tower, and no display grid can be made to fit it — 1152 is all the cells there are.
            const int kExpectedGFexTowers = 1152;
            if (gFexTowerEtValues && (int)gFexTowerEtValues->size() != kExpectedGFexTowers)
                printf("  gFEX event %d: container holds %d entries, but a gFexTowerContainer is "
                       "expected to hold %d (%.2f x) — the branch is not one entry per tower\n",
                       iEvt, (int)gFexTowerEtValues->size(), kExpectedGFexTowers,
                       (double)gFexTowerEtValues->size() / kExpectedGFexTowers);

            bool drawThisEvent = (nFexTowerDisplays < maxFexTowerDisplays);

            std::vector<FexCell> gCells, jCells;
            if (drawGFexTowers && gFexTowerEtValues && !gFexTowerEtValues->empty()) {
                unsigned int nUnplaced = 0;
                gCells = fexAccumulateCells(gFexTowerGeom, *gFexTowerEtValues,
                                            *gFexTowerEtaValues, *gFexTowerPhiValues, &nUnplaced);
                if (nUnplaced > 0)
                    printf("  gFEX event %d: %u entr(ies) at an eta/phi the grid does not cover\n",
                           iEvt, nUnplaced);
                printFexCellMultiplicity("gFEX", iEvt, gCells, gFexTowerEtValues->size());
                if (!gFexEtaGranularityDone) {
                    printFexEtaGranularity("gFEX", *gFexTowerEtaValues, gFexTowerGeom);
                    gFexEtaGranularityDone = true;
                }
                std::vector<FexCell> hot = fexTopCells(gCells, 10);
                gFexTopTen = hot;          // shown under the MET table on the text page
                printFexTopCells("gFEX", iEvt, hot);
                // Kept as raw values rather than filled straight into a histogram: the E_T scale
                // is one of the things under investigation, so the binning is chosen from the
                // data once every event has been seen.
                if (!hot.empty()) {
                    gFexLeadingTowerEt.push_back(hot.front().et);
                    gFexLeadingTowerEta.push_back(hot.front().eta);
                }
            }
            if (drawJFexTowers && jFexTowerEtValues && !jFexTowerEtValues->empty()) {
                // No explicit grid for jFEX, so derive one from the first event that has towers.
                if (!jFexTowerGeom.built) {
                    jFexTowerGeom = buildFexTowerGeometry(*jFexTowerEtaValues, *jFexTowerPhiValues,
                                                          "jFEX", debugFexTowerGeometry);
                    if (!jFexTowerGeom.usable())
                        std::cout << "  jFEX: could not derive a tower geometry — displays skipped\n";
                }
                unsigned int nUnplaced = 0;
                jCells = fexAccumulateCells(jFexTowerGeom, *jFexTowerEtValues,
                                            *jFexTowerEtaValues, *jFexTowerPhiValues, &nUnplaced);
                printFexCellMultiplicity("jFEX", iEvt, jCells, jFexTowerEtValues->size());
                printFexTopCells("jFEX", iEvt, fexTopCells(jCells, 10));
            }

            if (drawThisEvent) {
                if (drawGFexTowers) {
                    drawFexTowerDisplay("gFEX", gFexTowerEtValues, gFexTowerEtaValues, gFexTowerPhiValues,
                                        gCells, gFexTowerGeom, gFexDumpDone, "gFex", pdf_gTow);
                    drawTowerEtDistributions(cFexTower, pdf_gTow);
                }
                if (drawJFexTowers) {
                    drawFexTowerDisplay("jFEX", jFexTowerEtValues, jFexTowerEtaValues, jFexTowerPhiValues,
                                        jCells, jFexTowerGeom, jFexDumpDone, "jFex", pdf_jTow);
                    drawTowerEtDistributions(cFexTower, pdf_jTow);
                }
                ++nFexTowerDisplays;
            }
        }

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
        drawMETTextPage(cEventDisplay, pdf_ED);

        delete hEvent;
        delete hFrame;
    }

    // close multi-page PDFs
    cEventDisplay.Print(pdf_ED + ")");
    cXY.Print(pdf_XY + ")");
    if (drawGFexTowers) cFexTower.Print(pdf_gTow + ")");
    if (drawJFexTowers) cFexTower.Print(pdf_jTow + ")");

    // ===================================================================
    // Leading gFEX tower E_T over every accepted event, in its own directory.
    // ===================================================================
    if (drawGFexTowers && !gFexLeadingTowerEt.empty()) {
        TString gTowDir = "gTowerInvestigation/" + signalString + "_" + emuTag + "/";
        gSystem->mkdir(gTowDir, true);

        double loEt = gFexLeadingTowerEt.front(), hiEt = gFexLeadingTowerEt.front();
        for (double e : gFexLeadingTowerEt) { loEt = std::min(loEt, e); hiEt = std::max(hiEt, e); }
        double spanEt = hiEt - loEt;
        if (spanEt <= 0.0) spanEt = 1.0;
        loEt = std::max(0.0, loEt - 0.05 * spanEt);
        hiEt += 0.05 * spanEt;

        TH1F hLeadEt("hGFexLeadingTowerEt",
                     "Leading gFEX tower E_{T};Leading tower E_{T} [GeV];Events",
                     100, loEt, hiEt);
        TH1F hLeadEta("hGFexLeadingTowerEta",
                      "#eta of the leading gFEX tower;#eta;Events", 50, -5.0, 5.0);
        hLeadEt.SetStats(0);  hLeadEt.SetLineColor(kRed+1);  hLeadEt.SetLineWidth(2);
        hLeadEta.SetStats(0); hLeadEta.SetLineColor(kRed+1); hLeadEta.SetLineWidth(2);
        for (double e : gFexLeadingTowerEt)  hLeadEt.Fill(e);
        for (double e : gFexLeadingTowerEta) hLeadEta.Fill(e);

        TCanvas cLead("cLead", "cLead", 800, 700);
        TString pdf_lead = gTowDir + "gFexLeadingTower.pdf";
        cLead.Print(pdf_lead + "(");

        cLead.cd();
        cLead.SetLogy(true);
        hLeadEt.Draw("HIST");
        TLatex leadLat;
        leadLat.SetTextSize(0.026);
        leadLat.DrawLatexNDC(0.16, 0.96,
            Form("%s: %d events, leading gFEX tower E_{T} mean %.2f GeV, max %.2f GeV",
                 signalString.c_str(), (int)gFexLeadingTowerEt.size(), hLeadEt.GetMean(), hiEt));
        cLead.Print(pdf_lead);

        cLead.SetLogy(false);
        hLeadEta.Draw("HIST");
        cLead.Print(pdf_lead);
        cLead.Print(pdf_lead + ")");

        TFile fLead(gTowDir + "gFexLeadingTower.root", "RECREATE");
        hLeadEt.Write();
        hLeadEta.Write();
        fLead.Close();

        printf("Leading gFEX tower: %d events, E_T in [%.3f, %.3f] GeV, mean %.3f — written to %s\n",
               (int)gFexLeadingTowerEt.size(), loEt, hiEt, hLeadEt.GetMean(), gTowDir.Data());
    }

    std::cout << "MET displays: " << nAccepted << " events, g/jFEX tower displays: "
              << nFexTowerDisplays << " events (the first " << nFexTowerDisplays
              << " of the same event set)\n";
}

// -----------------------------------------------------------------------
void makeMETEventDisplays() {
    gErrorIgnoreLevel = kWarning;

    // Same input paths used by metAnalysisAndRates.C
    const std::string sigHER  = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/ZvvHbb_v4/mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    const std::string backHER = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/QCD_Dijet_JZ*_v4/mc21_14TeV_jj_JZ*_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    const std::string emuDir  = "/data/larsonma/GEPMET/outputNTuplesDev_METv3/";

    // Default signal call: ZvvHbb with the (jetEt=20, towerEt=2, EtaSK_OR, twrSF=1, jetSF=1) config.
    // The trailing arguments turn on the gFEX / jFEX tower display PDFs, which draw the first
    // maxFexTowerDisplays events of the same set, in the same order, as METEventDisplays.pdf.
    // fexTowerZMax = 0 auto-ranges the colour scale to the largest tower in each event;
    // fexTowerDumpEvent = -1 prints every tower position / footprint / E_T for the first
    // displayed event of each collection.
    callMakeMETEventDisplays(
        sigHER,
        emuDir + "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root",
        -1, true, "ZvvHbb",
        /*backMinJetMETOverTruth=*/50.0, /*backMinGFexMETOverTruth=*/150.0,
        /*drawGFexTowers=*/true, /*drawJFexTowers=*/true,
        /*fexTowerZMin=*/0.0, /*fexTowerZMax=*/0.0, /*fexTowerDumpEvent=*/-1);

    // Background: JZ1 only. Default filter accepts events where EITHER GEP Jet MET or gFEX JwoJ MET
    // is over-reconstructed vs truth NonInt MET by its respective threshold (50 GeV / 150 GeV).
    // To override, append e.g. /*backMinJetMETOverTruth=*/30.0, /*backMinGFexMETOverTruth=*/100.0
    callMakeMETEventDisplays(
        backHER,
        emuDir + "mc21_14TeV_jj_JZ_e8557_s4422_r16130_N_Towers_4096_jetEt15_towerEt2_EtaSK_NoOR_twrSF1_jetSF1.root",
        1, false, "jj_1",
        /*backMinJetMETOverTruth=*/50.0, /*backMinGFexMETOverTruth=*/150.0,
        /*drawGFexTowers=*/true, /*drawJFexTowers=*/true,
        /*fexTowerZMin=*/0.0, /*fexTowerZMax=*/0.0, /*fexTowerDumpEvent=*/-1);

    gSystem->Exit(0);
}
