// -----------------------------------------------------------------------------
// arctanPhiLUTPrototype.cc
//
// Standalone prototype for a fixed-point arctan LUT that computes the MET phi
// word (phi_MET = atan2(E_y, E_x)) in the same digitization the firmware uses.
//
// Idea:
//   * E_x, E_y arrive as signed digitized ET components (sign-magnitude, see
//     pack_signed_et in metEmulation.cc: 1 sign bit + (signed_et_bit_length_-1)
//     magnitude bits).
//   * atan2 only depends on the RATIO |E_y| / |E_x| plus the quadrant (signs).
//   * So we keep only the 5 most significant bits of |E_x| and |E_y|, aligned
//     to a common leading-bit position so the ratio is preserved, and use them
//     to index a small LUT: index = ax5 + (ay5 << REDUCE_BITS)   (i.e.
//     "E_x + E_y * 2^bitwidth", the flattened-2D layout).
//   * The LUT stores the first-quadrant reference angle as a phi code in
//     [0, g_nPhi/4]. The signs of E_x, E_y then fold it into the correct
//     quadrant of the full 0..2pi range, producing a 6-bit phi word.
//
// The prototype compares this LUT implementation against:
//   (a) the exact truth atan2, digitized, and
//   (b) the simple linear approximation phi ~= (pi/4) * (E_y / E_x).
//
// Run (interpreted or compiled):
//   root -b -l -q 'arctanPhiLUTPrototype.cc+()'
//
// This file is strictly additive: it does not modify the emulator or constants.
// -----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TSystem.h"

// Pull in the digitization constants (signed_et_bit_length_ for the E_T
// component magnitude width). The phi_missing word width is independent (below).
#include "metConstants/constants.h"

// -----------------------------------------------------------------------------
// Prototype configuration
// -----------------------------------------------------------------------------

// Number of most-significant bits of |E_x| and |E_y| kept to index the LUT.
// Configurable at runtime via the arctanPhiLUTPrototype() argument; these
// globals are set from that argument before the LUT is built. Setting it equal
// to ET_MAG_BITS (12) means no reduction, i.e. full E_T precision (and a
// LUT of 2^(2*bits) entries).
static unsigned int g_reduceBits = 5;
static unsigned int g_reduceMask = (1u << 5) - 1;

// Output phi_missing word width for this prototype. This is intentionally an
// independent quantity from phi_bit_length_ in constants.h (which digitizes
// tower/jet phi over [-3.2, 3.2]); the MET phi_missing word spans the full
// 0..2pi and is kept fully configurable here for prototyping.
constexpr unsigned int DEFAULT_PHI_MISSING_BITS = 10;

// phi-domain constants, all derived from the phi_missing word width and
// refreshed by setPhiBits() before the LUT is built.
static unsigned int g_phiBits    = DEFAULT_PHI_MISSING_BITS;
static unsigned int g_nPhi       = 1u << DEFAULT_PHI_MISSING_BITS;         // codes over 2pi
static unsigned int g_nPhiMask   = (1u << DEFAULT_PHI_MISSING_BITS) - 1;
static unsigned int g_phiQuarter = (1u << DEFAULT_PHI_MISSING_BITS) / 4;   // pi/2 in codes
static unsigned int g_phiHalf    = (1u << DEFAULT_PHI_MISSING_BITS) / 2;   // pi   in codes
static unsigned int g_phiEighth  = (1u << DEFAULT_PHI_MISSING_BITS) / 8;   // pi/4 in codes
static double       g_phiLsbRad  = (2.0 * M_PI) / double(1u << DEFAULT_PHI_MISSING_BITS);

// Refresh the derived phi-domain constants for a given phi word width.
static void setPhiBits(unsigned int phiBits) {
    g_phiBits    = phiBits;
    g_nPhi       = 1u << phiBits;
    g_nPhiMask   = (1u << phiBits) - 1;
    g_phiQuarter = (1u << phiBits) / 4;
    g_phiHalf    = (1u << phiBits) / 2;
    g_phiEighth  = (1u << phiBits) / 8;
    g_phiLsbRad  = (2.0 * M_PI) / double(1u << phiBits);
}

// Magnitude bit width of a signed digitized ET component (drop the sign bit).
constexpr unsigned int ET_MAG_BITS = signed_et_bit_length_ - 1;        // 12
constexpr uint32_t     ET_MAG_MAX  = (1u << ET_MAG_BITS) - 1;          // 4095

// -----------------------------------------------------------------------------
// Bit helpers
// -----------------------------------------------------------------------------

// Position (0-indexed) of the most significant set bit; -1 for value 0.
static inline int msbPosition(uint32_t v) {
    int pos = -1;
    for (int b = 31; b >= 0; --b) {
        if (v & (1u << b)) { pos = b; break; }
    }
    return pos;
}

// Reduce (ax, ay) to their REDUCE_BITS MSBs, aligned to the leading bit of the
// larger magnitude so the ratio ay/ax is preserved. In firmware this is a
// priority-encoder (leading-zero count) followed by a barrel shift.
static inline void reduceToMSBs(uint32_t ax, uint32_t ay,
                                uint32_t& ax_r, uint32_t& ay_r) {
    uint32_t m   = std::max(ax, ay);
    int      msb = msbPosition(m);

    // shift so the larger value's MSB lands in bit (g_reduceBits - 1).
    int shift = msb - int(g_reduceBits - 1);
    if (shift < 0) shift = 0;                 // small magnitudes: keep as-is

    ax_r = ax >> shift;
    ay_r = ay >> shift;

    if (ax_r > g_reduceMask) ax_r = g_reduceMask;  // safety clamp
    if (ay_r > g_reduceMask) ay_r = g_reduceMask;
}

// Circular difference a - b in phi codes, wrapped to (-g_nPhi/2, g_nPhi/2].
static inline int circDiffCodes(int a, int b) {
    int d = (a - b) % int(g_nPhi);
    if (d < 0)               d += g_nPhi;
    if (d > int(g_phiHalf))   d -= g_nPhi;
    return d;
}

// -----------------------------------------------------------------------------
// LUT: first-quadrant reference angle (phi code in [0, g_phiQuarter]).
// -----------------------------------------------------------------------------

static std::vector<uint8_t> g_arctanLUT;

static void buildArctanLUT() {
    g_arctanLUT.assign(size_t(1) << (2 * g_reduceBits), 0);
    for (uint32_t ay5 = 0; ay5 <= g_reduceMask; ++ay5) {
        for (uint32_t ax5 = 0; ax5 <= g_reduceMask; ++ax5) {
            double phi;                                   // reference angle [0, pi/2]
            if (ax5 == 0 && ay5 == 0) {
                phi = 0.0;                                // undefined; define as 0
            } else {
                phi = std::atan2(double(ay5), double(ax5));
            }
            // Digitize into phi codes: fraction of the full circle * g_nPhi.
            int code = int(std::lround(phi / (2.0 * M_PI) * double(g_nPhi)));
            if (code < 0)                  code = 0;
            if (code > int(g_phiQuarter))   code = g_phiQuarter;   // clamp to a quadrant
            g_arctanLUT[ax5 + (ay5 << g_reduceBits)] = uint8_t(code);
        }
    }
}

// Fold a first-quadrant reference-angle code into the full 0..2pi phi word
// using the signs of E_x, E_y.
static inline int foldQuadrant(int alphaCode, bool negX, bool negY) {
    int phi;
    if (!negX && !negY)      phi = alphaCode;                 // Q1: [0, pi/2)
    else if (negX && !negY)  phi = int(g_phiHalf)   - alphaCode; // Q2: (pi/2, pi]
    else if (negX &&  negY)  phi = int(g_phiHalf)   + alphaCode; // Q3: (pi, 3pi/2]
    else                     phi = int(g_nPhi)      - alphaCode; // Q4: (3pi/2, 2pi)
    return phi & int(g_nPhiMask);                             // wrap mod g_nPhi
}

// -----------------------------------------------------------------------------
// The three phi estimators (all return a 6-bit phi code in [0, g_nPhi)).
// -----------------------------------------------------------------------------

// (1) LUT-based fixed-point implementation.
static int phiFromLUT(int Ex, int Ey) {
    uint32_t ax = uint32_t(std::abs(Ex));
    uint32_t ay = uint32_t(std::abs(Ey));

    uint32_t ax5, ay5;
    reduceToMSBs(ax, ay, ax5, ay5);

    int alphaCode = g_arctanLUT[ax5 + (ay5 << g_reduceBits)];
    return foldQuadrant(alphaCode, Ex < 0, Ey < 0);
}

// (2) Simple linear approximation phi ~= (pi/4) * (E_y / E_x).
//     In phi codes, (pi/4) maps to g_phiEighth (= 8), so the reference angle is
//     g_phiEighth * (|E_y| / |E_x|). Since this only stays sane for
//     |E_y| <= |E_x|, we fold the octant (use the complement when |E_y| > |E_x|)
//     so the comparison is fair across the whole quadrant. The underlying
//     approximation is still the requested pi/4 * ratio.
static int phiLinear(int Ex, int Ey) {
    double ax = std::abs(double(Ex));
    double ay = std::abs(double(Ey));

    double alpha;                                            // reference angle in codes
    if (ax == 0.0 && ay == 0.0) {
        alpha = 0.0;
    } else if (ay <= ax) {
        alpha = double(g_phiEighth) * (ay / ax);             // pi/4 * (Ey/Ex)
    } else {
        alpha = double(g_phiQuarter) - double(g_phiEighth) * (ax / ay); // complement
    }

    int alphaCode = int(std::lround(alpha));
    if (alphaCode < 0)                 alphaCode = 0;
    if (alphaCode > int(g_phiQuarter))  alphaCode = g_phiQuarter;
    return foldQuadrant(alphaCode, Ex < 0, Ey < 0);
}

// (3) Truth reference: exact atan2, digitized to a 6-bit phi word.
static int phiTruth(int Ex, int Ey) {
    double a = std::atan2(double(Ey), double(Ex));          // [-pi, pi]
    if (a < 0.0) a += 2.0 * M_PI;                           // [0, 2pi)
    int code = int(std::lround(a / (2.0 * M_PI) * double(g_nPhi)));
    return code & int(g_nPhiMask);
}

// -----------------------------------------------------------------------------
// Optional: dump the LUT as a C array, mirroring the sinLUT_ style in
// constants.h, so it can be pasted into the firmware/emulation constants.
// -----------------------------------------------------------------------------

static void writeLUTHeader(const std::string& path) {
    const size_t lutSize = g_arctanLUT.size();
    std::ofstream out(path);
    out << "// Auto-generated by arctanPhiLUTPrototype.cc\n";
    out << "// First-quadrant reference-angle phi codes, indexed by\n";
    out << "//   index = ax5 + (ay5 << " << g_reduceBits << ")\n";
    out << "// where ax5, ay5 are the top " << g_reduceBits
        << " bits of |E_x|, |E_y| (ratio-aligned).\n";
    out << "static const unsigned char arctanPhiLUT_[" << lutSize << "] = {\n    ";
    for (size_t i = 0; i < lutSize; ++i) {
        out << std::setw(2) << unsigned(g_arctanLUT[i]);
        if (i + 1 < lutSize) out << ((i + 1) % g_reduceMask == 0 ? ",\n    " : ", ");
    }
    out << "\n};\n";
    out.close();
    std::cout << "Wrote LUT header: " << path << "\n";
}

// -----------------------------------------------------------------------------
// Main driver: scan random (E_x, E_y), compare estimators, make ROOT plots.
// -----------------------------------------------------------------------------

void arctanPhiLUTPrototype(unsigned int reduceBits = 5,
                           unsigned int nSamples   = 500000,
                           bool         dumpLUT     = true,
                            unsigned int phiBits    = DEFAULT_PHI_MISSING_BITS) {
    // Configure the LUT reduction bit-width. reduceBits == ET_MAG_BITS (12)
    // means no reduction, i.e. full E_T precision.
    if (reduceBits < 1)           reduceBits = 1;
    if (reduceBits > ET_MAG_BITS) reduceBits = ET_MAG_BITS;   // cannot exceed input precision
    g_reduceBits = reduceBits;
    g_reduceMask = (1u << reduceBits) - 1;

    // Configure the output phi ("phi_missing") word width and refresh the
    // phi-domain constants. Clamp to a sane range (>=3 keeps pi/4 representable).
    if (phiBits < 3)  phiBits = 3;
    if (phiBits > 16) phiBits = 16;
    setPhiBits(phiBits);

    // Output directory keyed by the E_T LUT-index bits and the phi_missing
    // word width.
    const std::string outDir =
        "arctanPhiLUTPrototype_plots_et" + std::to_string(reduceBits) +
        "b_phimiss" + std::to_string(phiBits) + "b";
    gSystem->mkdir(outDir.c_str(), true);

    buildArctanLUT();
    if (dumpLUT) writeLUTHeader(outDir + "/arctanPhiLUT.h");

    std::cout << "Configuration:\n";
    std::cout << "  signed_et_bit_length_ = " << signed_et_bit_length_ << "  (magnitude bits = " << ET_MAG_BITS << ")\n";
    std::cout << "  reduceBits (E_T LUT)  = " << g_reduceBits << "  (LUT size = " << g_arctanLUT.size() << ")\n";
    std::cout << "  phiBits (phi word)    = " << g_phiBits << "  (N_PHI = " << g_nPhi << " digitized values over 2pi)\n";
    std::cout << "  phi LSB               = " << g_phiLsbRad << " rad\n\n";

    gStyle->SetOptStat(0);

    // Residual histograms, in phi digitized values (1 digitized value = g_phiLsbRad).
    TH1F* hResLUT = new TH1F("hResLUT",
        "phi residual (estimate - truth);#Delta#phi [digitized value];entries",
        2 * g_nPhi + 1, -double(g_nPhi) - 0.5, double(g_nPhi) + 0.5);
    TH1F* hResLin = new TH1F("hResLin", "linear residual",
        2 * g_nPhi + 1, -double(g_nPhi) - 0.5, double(g_nPhi) + 0.5);
    hResLUT->SetLineColor(kAzure + 1);
    hResLUT->SetLineWidth(2);
    hResLin->SetLineColor(kOrange + 7);
    hResLin->SetLineWidth(2);

    // |error| vs true phi, to expose quadrant/octant structure.
    TProfile* pErrLUT = new TProfile("pErrLUT",
        "mean |#Delta#phi| vs true #phi;#phi_{truth} [digitized value];#LT|#Delta#phi|#GT [digitized value]",
        g_nPhi, -0.5, double(g_nPhi) - 0.5);
    TProfile* pErrLin = new TProfile("pErrLin", "linear", g_nPhi, -0.5, double(g_nPhi) - 0.5);
    pErrLUT->SetLineColor(kAzure + 1);
    pErrLUT->SetMarkerColor(kAzure + 1);
    pErrLUT->SetLineWidth(2);
    pErrLin->SetLineColor(kOrange + 7);
    pErrLin->SetMarkerColor(kOrange + 7);
    pErrLin->SetLineWidth(2);

    // LUT phi vs truth phi correlation.
    TH2F* hCorrLUT = new TH2F("hCorrLUT",
        "LUT #phi vs truth #phi;#phi_{truth} [digitized value];#phi_{LUT} [digitized value]",
        g_nPhi, -0.5, double(g_nPhi) - 0.5, g_nPhi, -0.5, double(g_nPhi) - 0.5);

    // Linear phi vs truth phi correlation.
    TH2F* hCorrLin = new TH2F("hCorrLin",
        "linear #phi vs truth #phi;#phi_{truth} [digitized value];#phi_{linear} [digitized value]",
        g_nPhi, -0.5, double(g_nPhi) - 0.5, g_nPhi, -0.5, double(g_nPhi) - 0.5);

    // Random signed components across the full magnitude range.
    std::mt19937 rng(12345);
    std::uniform_int_distribution<int> mag(0, int(ET_MAG_MAX));
    std::uniform_int_distribution<int> sgn(0, 1);

    long   nCounted   = 0;
    double sumSqLUT   = 0.0, sumSqLin = 0.0;
    int    maxAbsLUT  = 0,   maxAbsLin = 0;

    for (unsigned int i = 0; i < nSamples; ++i) {
        int Ex = mag(rng) * (sgn(rng) ? 1 : -1);
        int Ey = mag(rng) * (sgn(rng) ? 1 : -1);
        if (Ex == 0 && Ey == 0) continue;                   // phi undefined

        int pT = phiTruth(Ex, Ey);
        int pL = phiFromLUT(Ex, Ey);
        int pN = phiLinear(Ex, Ey);

        int dL = circDiffCodes(pL, pT);
        int dN = circDiffCodes(pN, pT);

        hResLUT->Fill(dL);
        hResLin->Fill(dN);
        pErrLUT->Fill(pT, std::abs(dL));
        pErrLin->Fill(pT, std::abs(dN));
        hCorrLUT->Fill(pT, pL);
        hCorrLin->Fill(pT, pN);

        sumSqLUT  += double(dL) * dL;
        sumSqLin  += double(dN) * dN;
        maxAbsLUT  = std::max(maxAbsLUT, std::abs(dL));
        maxAbsLin  = std::max(maxAbsLin, std::abs(dN));
        ++nCounted;
    }

    double rmsLUT = std::sqrt(sumSqLUT / double(nCounted));
    double rmsLin = std::sqrt(sumSqLin / double(nCounted));

    std::cout << "Samples used: " << nCounted << "\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "LUT    : RMS = " << rmsLUT << " digitized value ("
              << rmsLUT * g_phiLsbRad << " rad), max|err| = "
              << maxAbsLUT << " digitized value\n";
    std::cout << "Linear : RMS = " << rmsLin << " digitized value ("
              << rmsLin * g_phiLsbRad << " rad), max|err| = "
              << maxAbsLin << " digitized value\n";

    // ---- Plots: one PDF per panel in the per-bit-width output directory ----
    // (outDir was created at the top of the function.)

    // (1) phi residual histogram (LUT vs linear) with RMS + max-error legend.
    {
        TCanvas* c1 = new TCanvas("cPhiResidual", "phi residual", 800, 600);
        c1->SetLogy();
        double maxY = std::max(hResLUT->GetMaximum(), hResLin->GetMaximum());
        hResLUT->SetMaximum(maxY * 2.0);
        hResLUT->GetXaxis()->SetRangeUser(-double(g_phiQuarter), double(g_phiQuarter));
        hResLUT->Draw("HIST");
        hResLin->Draw("HIST SAME");
        TLegend* leg = new TLegend(0.36, 0.72, 0.88, 0.88);
        leg->SetHeader("residual = estimate - truth [digitized value]");
        leg->AddEntry(hResLUT, Form("LUT: RMS %.2f, max %d (%.3f rad)",
                                    rmsLUT, maxAbsLUT, maxAbsLUT * g_phiLsbRad), "l");
        leg->AddEntry(hResLin, Form("Linear: RMS %.2f, max %d (%.3f rad)",
                                    rmsLin, maxAbsLin, maxAbsLin * g_phiLsbRad), "l");
        leg->Draw();
        c1->SaveAs((outDir + "/phiResidual.pdf").c_str());
    }

    // (2) mean |error| vs true phi.
    {
        TCanvas* c2 = new TCanvas("cMeanAbsErr", "mean |dphi| vs phi", 800, 600);
        double maxP = std::max(pErrLUT->GetMaximum(), pErrLin->GetMaximum());
        pErrLUT->SetMaximum(maxP * 1.3);
        pErrLUT->SetMinimum(0.0);
        pErrLUT->Draw("HIST L");
        pErrLin->Draw("HIST L SAME");
        TLegend* leg = new TLegend(0.60, 0.75, 0.88, 0.88);
        leg->AddEntry(pErrLUT, "LUT", "l");
        leg->AddEntry(pErrLin, "Linear", "l");
        leg->Draw();
        c2->SaveAs((outDir + "/meanAbsErr_vs_phi.pdf").c_str());
    }

    // (3) LUT phi vs truth phi correlation.
    {
        TCanvas* c3 = new TCanvas("cCorrLUT", "LUT phi vs truth phi", 800, 600);
        hCorrLUT->Draw("COLZ");
        c3->SaveAs((outDir + "/LUTphi_vs_truthphi.pdf").c_str());
    }

    // (4) Linear phi vs truth phi correlation.
    {
        TCanvas* c4 = new TCanvas("cCorrLin", "linear phi vs truth phi", 800, 600);
        hCorrLin->Draw("COLZ");
        c4->SaveAs((outDir + "/linearphi_vs_truthphi.pdf").c_str());
    }

    // (5) LUT residual expressed in radians.
    TH1F* hResLUTrad = new TH1F("hResLUTrad",
        "LUT phi residual;#Delta#phi [rad];entries",
        2 * g_nPhi + 1, (-double(g_nPhi) - 0.5) * g_phiLsbRad, (double(g_nPhi) + 0.5) * g_phiLsbRad);
    for (int b = 1; b <= hResLUT->GetNbinsX(); ++b) {
        double code = hResLUT->GetBinCenter(b);
        hResLUTrad->Fill(code * g_phiLsbRad, hResLUT->GetBinContent(b));
    }
    hResLUTrad->SetLineColor(kAzure + 1);
    hResLUTrad->SetLineWidth(2);
    hResLUTrad->GetXaxis()->SetRangeUser(-double(g_phiQuarter) * g_phiLsbRad,
                                          double(g_phiQuarter) * g_phiLsbRad);
    {
        TCanvas* c5 = new TCanvas("cResLUTrad", "LUT phi residual [rad]", 800, 600);
        c5->SetLogy();
        hResLUTrad->Draw("HIST");
        c5->SaveAs((outDir + "/LUTphiResidual_rad.pdf").c_str());
    }

    // ---- Transfer-function overlay: reference angle vs |E_y|/|E_x| ----
    // Sweep the ratio and plot the first-quadrant reference angle three ways:
    // the full-precision arctan (the function we are approximating), the LUT
    // staircase, and the pi/4 * ratio linear approximation. This is the direct
    // "how well does the LUT track arctan" comparison. |E_x| is fixed and
    // |E_y| = r * |E_x|, both positive, so the estimators return the Q1 angle.
    // The RMS of each approximation about the exact arctan (in rad) is reported.
    // Keep pointers to the full-range graphs so they can be written to the file.
    TGraph* gArctan = nullptr;
    TGraph* gLUTtf  = nullptr;
    TGraph* gLinTf  = nullptr;

    auto makeTransferFunction = [&](double rMax, const std::string& pdfPath,
                                    const std::string& tag,
                                    TGraph*& gAout, TGraph*& gLout, TGraph*& gNout) {
        const int nSweep  = 400;
        const int ExFixed = 1000;                       // fixed |E_x|; |E_y| = r * |E_x|
        TGraph* gA = new TGraph();                       // full-precision arctan
        TGraph* gL = new TGraph();                       // LUT staircase
        TGraph* gN = new TGraph();                       // linear approximation
        double sumSqL = 0.0, sumSqN = 0.0;
        int    nPts   = 0;
        for (int i = 0; i <= nSweep; ++i) {
            double r     = rMax * double(i) / double(nSweep);
            int    Ey    = int(std::lround(r * ExFixed));
            double exact = std::atan(r);
            double lut   = phiFromLUT(ExFixed, Ey) * g_phiLsbRad;
            double lin   = phiLinear(ExFixed, Ey) * g_phiLsbRad;
            gA->SetPoint(i, r, exact);
            gL->SetPoint(i, r, lut);
            gN->SetPoint(i, r, lin);
            sumSqL += (lut - exact) * (lut - exact);
            sumSqN += (lin - exact) * (lin - exact);
            ++nPts;
        }
        double rmsL = std::sqrt(sumSqL / double(nPts));
        double rmsN = std::sqrt(sumSqN / double(nPts));

        gA->SetName(("gArctan_" + tag).c_str());
        gA->SetTitle("reference angle vs |E_{y}|/|E_{x}|;|E_{y}| / |E_{x}|;#phi [rad]");
        gA->SetLineColor(kBlack);
        gA->SetLineWidth(2);
        gL->SetName(("gLUTtf_" + tag).c_str());
        gL->SetLineColor(kAzure + 1);
        gL->SetLineWidth(2);
        gN->SetName(("gLinTf_" + tag).c_str());
        gN->SetLineColor(kOrange + 7);
        gN->SetLineWidth(2);
        gN->SetLineStyle(2);

        TCanvas* ct = new TCanvas(("cTransfer_" + tag).c_str(), "reference angle", 800, 600);
        gA->Draw("AL");
        gL->Draw("L SAME");
        gN->Draw("L SAME");
        // Legend in the bottom-right so it does not block the rising curves.
        TLegend* leg = new TLegend(0.48, 0.15, 0.88, 0.44);
        leg->SetHeader(Form("E_{x,y}: %u bit,  #phi_{miss}: %u bit",
                            g_reduceBits, g_phiBits));
        leg->AddEntry(gA, "arctan (full precision)", "l");
        leg->AddEntry(gL, Form("LUT (RMS %.3f rad)", rmsL), "l");
        leg->AddEntry(gN, Form("linear (#pi/4)#upoint r (RMS %.3f rad)", rmsN), "l");
        leg->Draw();
        ct->SaveAs(pdfPath.c_str());

        gAout = gA; gLout = gL; gNout = gN;
    };

    // Full range (out to |E_y|/|E_x| = 4) and a zoomed duplicate to the octant
    // boundary (|E_y|/|E_x| = 1), where the linear approximation is defined
    // without folding.
    constexpr double rMaxZoom = 1.0;
    makeTransferFunction(4.0, outDir + "/referenceAngle_vs_ratio.pdf", "full",
                         gArctan, gLUTtf, gLinTf);
    TGraph *gAzoom = nullptr, *gLzoom = nullptr, *gNzoom = nullptr;
    makeTransferFunction(rMaxZoom, outDir + "/referenceAngle_vs_ratio_zoom.pdf", "zoom",
                         gAzoom, gLzoom, gNzoom);

    // ---- Persist objects for later inspection ----
    TFile* out = TFile::Open((outDir + "/arctanPhiLUTPrototype.root").c_str(), "RECREATE");
    hResLUT->Write();
    hResLin->Write();
    pErrLUT->Write();
    pErrLin->Write();
    hCorrLUT->Write();
    hCorrLin->Write();
    hResLUTrad->Write();
    gArctan->Write();
    gLUTtf->Write();
    gLinTf->Write();
    gAzoom->Write();
    gLzoom->Write();
    gNzoom->Write();
    out->Close();

    std::cout << "\nWrote per-panel PDFs, LUT header, and ROOT file to "
              << outDir << "/\n";
}

// Allow running as a compiled executable too (root ... .cc+ calls the function
// of the same name; a bare a.out build uses main()).
#ifndef __CINT__
int main() {
    arctanPhiLUTPrototype();
    return 0;
}
#endif
