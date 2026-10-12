#include <string>
#include <sstream>
#include <iomanip>  // for std::setprecision
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>  // std::llabs, used by the MET TOB encoders below
#include <vector>   // for the GEP JwoJ per-eta-bin parameter table
#include "metConstants/constants.h"

// Wrap unsigned phi index into [0, two_pi_digitized_in_phi_]
inline unsigned int wrapPhiUnsigned(unsigned int phi) {
    unsigned int phiWrapped = phi;
    if (phi > two_pi_digitized_in_phi_) phiWrapped = phi - (two_pi_digitized_in_phi_ + 1);
    return phiWrapped;
}

// Wrap an unsigned HALF-code phi index into [0, two_pi_digitized_in_phi_half_].
// The half-code counterpart of wrapPhiUnsigned above, for sinLUTHalf_. One
// subtraction is enough for the same reason it is there: the only index this is
// handed is a phi half-code (at most 126) plus the pi/2 offset (32), so it can
// overshoot the range by less than one full turn.
inline unsigned int wrapPhiHalfUnsigned(unsigned int phiHalf) {
    unsigned int phiWrapped = phiHalf;
    if (phiHalf > two_pi_digitized_in_phi_half_) phiWrapped = phiHalf - (two_pi_digitized_in_phi_half_ + 1);
    return phiWrapped;
}

// Compute digitized deltaR^2 between two (eta, phi) coordinates, with phi wrap.
// The fold is around pi_digitized_in_phi_ = phi_range_ / 2, so codes on opposite
// sides of the +-pi seam come out a code apart rather than most of a turn apart.
inline unsigned int digitizedDeltaR2(unsigned int eta1, unsigned int phi1, unsigned int eta2, unsigned int phi2) {
    unsigned int uDeltaEta = static_cast<unsigned int>(std::abs(int(eta1) - int(eta2)));
    unsigned int uDeltaPhi = static_cast<unsigned int>(std::abs(int(phi1) - int(phi2)));
    if (uDeltaPhi >= pi_digitized_in_phi_) uDeltaPhi = (2 * pi_digitized_in_phi_) - uDeltaPhi;
    return uDeltaEta * uDeltaEta + uDeltaPhi * uDeltaPhi;
}

// Digitize phi onto the tower grid.
//
// Kept separate from digitize() because phi is periodic while eta/Et are not:
// the generic function saturates at range - 1, which is wrong at both ends of
// the phi axis. A value in the top half tower (phi within pi/64 of +pi) belongs
// in code 0, not in a code one past the end of the range -- which would not even
// fit the phi field. The code count is phi_range_ rather than
// 1 << phi_bit_length_ so the grid stays tied to the towers, not to the field
// width. Same definition as the jet tagger chain in emulationHelperFunctions.h.
inline unsigned int digitize_phi(double phi) {
    const int nPhi = static_cast<int>(phi_range_);
    const int code = static_cast<int>(std::lround((phi - phi_min_) / phi_granularity_));
    return static_cast<unsigned int>(((code % nPhi) + nPhi) % nPhi);
}

// The reconstruction tag encodes the pileup scenario: PU200 samples are r16130,
// PU140 samples r16129 (HERNTupler stamps this into the ntuple names). Carrying it
// through to the emulation outputs is what keeps the two scenarios distinguishable
// by file name alone, since both land in the same output directory.
// For the makeInputFileName fallback the ntuple directory is switched as well,
// PU140 ntuples living in the parallel ntuples_PU140 tree.
inline std::string applyPileupTags(std::string path, unsigned int pileup) {
    if (pileup != 140) return path;
    auto replaceFirst = [&path](const std::string& from, const std::string& to) {
        size_t pos = path.find(from);
        if (pos != std::string::npos) path.replace(pos, from.size(), to);
    };
    replaceFirst("_r16130_", "_r16129_");
    replaceFirst("/ntuples/", "/ntuples_PU140/");
    return path;
}

// Returns input NTuple file name given parameters
std::string makeInputFileName(bool signalBool, std::string signalString,
                              std::string inputRootFilePath = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/",
                              unsigned int pileup = 200) {
    std::ostringstream ss;

    if (signalBool) {
        if(signalString == "VBF_hh_bbbb_cvv0") ss << inputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv0_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if(signalString == "VBF_hh_bbbb_cvv1") ss << inputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv1_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ggF_hh_bbbb") ss << inputRootFilePath << "ggF_HHbbbb_v3/mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ZvvHbb") ss << inputRootFilePath << "ZvvHbb_v3/mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_resim_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ttbar_had") ss << inputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "Zprime_ttbar") ss << inputRootFilePath << "mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ttbar_semilep") ss << inputRootFilePath << "ttbar_semilep_v4/mc21_14TeV_ttbar_hdamp258p75_semilep_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ttbar_dilep") ss << inputRootFilePath << "ttbar_dilep_v4/mc21_14TeV_ttbar_hdamp258p75_dilep_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "Zmumu") ss << inputRootFilePath << "Zmumu_v4/mc21_14TeV_Zmumu_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    } else {
        ss << inputRootFilePath << "mc21_14TeV_jj_JZ_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    }
    return applyPileupTags(ss.str(), pileup);
}

unsigned int digitize(double value, int bit_length, double min_val, double max_val, unsigned int altRange = 0) {
    unsigned int range = (altRange == 0) ? (1u << bit_length) : altRange;
    double scale = double(range) / (max_val - min_val);
    //std::cout << "max_val - scale; " << max_val - (1/scale) << "\n";
    // Check if value is in range
    if (value < min_val) {
        value = min_val;
        //std::cout << "Warning: Value " << value
        //  << " is out of range (" << min_val
        //  << ", " << max_val << ")\n";
    }
    if (value >= max_val){
        return range - 1;
    }

    return static_cast<unsigned int>(std::round((value - min_val) * scale));
}

// --- GEP JwoJ MET parameters ------------------------------------------------------------
// The GEP JwoJ MET algorithm splits the tower collection in two, in the spirit of the gFEX
// "jets without jets" algorithm:
//
//   hard term : BLOCKS of blockSize x blockSize adjacent towers whose summed E_T is strictly
//               above hardEtThreshold, each contributing once at the block's geometric centre
//   soft term : the individual towers of every block that did not pass, each at its own
//               tower centre
//
// Blocks select, towers carry the remainder -- the same division of labour as gFEX, where
// gBlocks pick the hard term and the surviving gTowers make up the soft one. It keeps the
// diffuse soft term at full 0.1 x 0.1 resolution while the hard term is judged on local
// energy density rather than on single towers.
//
// and recombines them as
//
//   MET_{x,y} = -( hardCoeff * SUM_hard E_T{x,y} + softCoeff * SUM_soft E_T{x,y} )
//
// The three numbers are held as per-eta-bin vectors rather than as scalars so an
// eta-dependent calibration can be dropped in later without the event loop changing: the
// loop already looks each tower's parameters up by its digitized eta code. What exists today
// is the single-bin table makeUniformGEPJwoJParams builds, which spans the whole eta range
// and so reproduces the scalar behaviour exactly.
//
// The binning is on the digitized eta CODE and not on the double, for the same reason the
// rest of this chain works in codes: a bin boundary has to fall between two towers of the
// grid rather than part way through one, or the firmware and the emulator will disagree
// about which side of it a tower sits on.
struct GEPJwoJParams {
    // Inclusive upper eta code of each bin, ascending. The last entry must be eta_range_ - 1
    // so that every code on the grid lands in a bin.
    std::vector<unsigned int> etaCodeUpperEdges;
    std::vector<double>       hardEtThreshold;   // [GeV], one per bin
    std::vector<double>       hardCoeff;         // one per bin
    std::vector<double>       softCoeff;         // one per bin

    // Side of the square tower block the hard-term threshold is applied to, in towers.
    // A whole-algorithm setting rather than a per-eta-bin one, so it is a scalar here.
    //
    // Only 1 and 2 are supported, and that is a property of this grid rather than a gap in
    // the arithmetic: 1 and 2 both divide it exactly (phi 64 = 32 x 2, eta 98 = 49 x 2), so
    // every block is full and the phi ring stays rotationally symmetric. 3 divides NEITHER
    // (phi 64 = 21 x 3 + 1, eta 98 = 32 x 3 + 2), which would leave a one-tower-wide block at
    // the phi seam judged against the same threshold as its nine-tower neighbours -- a phi
    // slice systematically pushed into the soft term, i.e. a preferred MET direction. Adding
    // 3 means first deciding how that seam is handled; see isSupportedGEPJwoJBlockSize.
    //
    // blockSize 1 is the single-tower algorithm and is bit-identical to what this code did
    // before blocks existed: a 1x1 block is its own tower, and its centre half-code 2k reads
    // the same LUT value as sinLUT_[k].
    unsigned int blockSize = 1;

    unsigned int nBins() const { return static_cast<unsigned int>(etaCodeUpperEdges.size()); }

    // Bin holding this eta code. A linear scan: the table is one entry today and would be a
    // handful at most, so anything cleverer would cost more than it saves. A code past the
    // last edge -- which digitize() should never produce, since it saturates at range - 1 --
    // falls into the last bin rather than off the end of the vectors.
    unsigned int binForEtaCode(unsigned int etaCode) const {
        for (unsigned int iBin = 0; iBin + 1 < nBins(); ++iBin)
            if (etaCode <= etaCodeUpperEdges[iBin]) return iBin;
        return nBins() - 1;
    }
};

// A block size is supported only if it tiles BOTH axes of the tower grid exactly. Checked
// rather than hard-coded to {1, 2} so that the reason is in the code: change the grid and the
// supported set follows, and a size that would leave a ragged seam is rejected instead of
// silently producing one.
inline bool isSupportedGEPJwoJBlockSize(unsigned int blockSize) {
    return blockSize >= 1 && (phi_range_ % blockSize) == 0 && (eta_range_ % blockSize) == 0;
}

// One eta bin covering the whole grid: what the scalar command-line parameters mean.
inline GEPJwoJParams makeUniformGEPJwoJParams(double hardEtThreshold,
                                              double hardCoeff,
                                              double softCoeff,
                                              unsigned int blockSize = 1) {
    GEPJwoJParams params;
    params.etaCodeUpperEdges = { eta_range_ - 1 };
    params.hardEtThreshold   = { hardEtThreshold };
    params.hardCoeff         = { hardCoeff };
    params.softCoeff         = { softCoeff };
    if (!isSupportedGEPJwoJBlockSize(blockSize)) {
        std::cerr << "Error: GEP JwoJ blockSize " << blockSize << " does not tile the tower grid ("
                  << phi_range_ << " phi x " << eta_range_ << " eta) -- falling back to 1x1"
                  << std::endl;
        blockSize = 1;
    }
    params.blockSize = blockSize;
    return params;
}

// Build an eta-binned table from physical eta upper edges, for when the calibration exists.
// The edges are digitized onto the tower grid here, once, so the event loop never sees a
// double; the last bin is forced to eta_range_ - 1 whatever was passed, so no code can fall
// outside the table.
//
// Nothing calls this yet -- it is the entry point the eta-dependent calibration is meant to
// arrive through, kept beside the uniform builder so the two cannot drift apart.
inline GEPJwoJParams makeEtaBinnedGEPJwoJParams(const std::vector<double>& etaUpperEdges,
                                                const std::vector<double>& hardEtThresholds,
                                                const std::vector<double>& hardCoeffs,
                                                const std::vector<double>& softCoeffs,
                                                unsigned int blockSize = 1) {
    const size_t nBins = etaUpperEdges.size();
    if (nBins == 0 || hardEtThresholds.size() != nBins ||
        hardCoeffs.size() != nBins || softCoeffs.size() != nBins) {
        std::cerr << "Error: GEP JwoJ eta-binned parameter vectors have mismatched sizes"
                  << " -- falling back to unit coefficients with no hard term" << std::endl;
        return makeUniformGEPJwoJParams(0.0, 1.0, 1.0, blockSize);
    }
    GEPJwoJParams params;
    if (!isSupportedGEPJwoJBlockSize(blockSize)) {
        std::cerr << "Error: GEP JwoJ blockSize " << blockSize << " does not tile the tower grid ("
                  << phi_range_ << " phi x " << eta_range_ << " eta) -- falling back to 1x1"
                  << std::endl;
        blockSize = 1;
    }
    params.blockSize = blockSize;
    params.hardEtThreshold = hardEtThresholds;
    params.hardCoeff       = hardCoeffs;
    params.softCoeff       = softCoeffs;
    params.etaCodeUpperEdges.reserve(nBins);
    for (size_t iBin = 0; iBin < nBins; ++iBin)
        params.etaCodeUpperEdges.push_back(digitize(etaUpperEdges[iBin], eta_bit_length_, eta_min_, eta_max_, eta_range_));
    params.etaCodeUpperEdges.back() = eta_range_ - 1;
    return params;
}

// Sign-magnitude: MSB = sign bit, remaining signed_et_bit_length_-1 bits = magnitude
template<size_t signed_et_bit_length_>
inline double undigitize_signed_et(const std::bitset<signed_et_bit_length_>& bits) {
    std::cout << "bits: " << bits << "\n";
    int sign = bits[signed_et_bit_length_ - 1] ? -1 : 1;
    std::cout << "sign: " << sign << "\n";
    int mag = static_cast<int>(bits.to_ulong()) & ((1u << (signed_et_bit_length_ - 1)) - 1);
    std::cout << "mag: " << mag << "\n";
    std::cout << "output: " << sign * mag * et_granularity_ << "\n";
    return sign * mag * et_granularity_;
}

inline double undigitize_et(const std::bitset<et_bit_length_>& et_bits) {
    return et_bits.to_ulong() * et_granularity_;
}

// =====================================================================================
// ===                MET TOB encoding: phi, integer root, saturation                ===
// =====================================================================================
// Ports of the MET firmware modules, bit-for-bit. Each function below names the module
// it comes from; between them they produce every field of the 64-bit MET TOB word laid
// out in Table 4.5 (see the writer in metEmulation.cc).
//
// These are deliberately integer-only and free of any reference to the tower grid's
// digitization. They take E_T COUNTS in and give TOB FIELDS out, so the same functions
// serve the ROOT event loop, the memory-print validation path, and -- once it exists --
// GlobalSim, none of which need agree about where their inputs came from.

// --- missing-energy azimuth, MET_PHI_COMPARATOR.v ------------------------------------
// Divider-free comparator ladder. The first octant is quantized by met_phi_octant_
// tangent midpoints (eight at the 6-bit field width), and the octant/quadrant is then
// reconstructed from the swap flag and the two sign bits. No LUT and no division: the
// only arithmetic is compare, add and a constant multiply.
//
// CONVENTION, and it is NOT the tower grid's. This index is an azimuth measured from +x
// over [0, 2*pi) in met_phi_range_ bins, i.e. at the 6-bit field width:
//
//     bin 0 = +x,  bin 16 = +y,  bin 32 = -x,  bin 48 = -y
//
// whereas digitize_phi() puts index 0 at phi_min_ = -pi + pi/64. A MET phi index must
// never be fed to digitize_phi's inverse, nor to sinLUT_, nor compared against a tower's
// phi index. The zero vector maps to bin 0, and undigitize_met_phi() below is the only
// inverse.
//
// OUTPUT phi is sized by its FIELD WIDTH -- the opposite of how the tower grid is sized,
// and deliberately so.
//
// An INPUT phi is an index into a physical grid of phi_range_ = 64 towers, so there the
// tower count is primary and phi_bit_length_ is only what carries it: a wider field would
// not buy a finer tower. The missing-energy azimuth is not on that grid at all. It is
// arctan(Ey, Ex), which is continuous, so nothing physical caps its resolution -- the
// only thing that does is the width of the field it ships in. Widen that field and the
// ladder resolves more finely. Hence met_phi_range_ derives FROM met_phi_bit_length_
// here, and EVERYTHING below follows from it: the bin count, the quadrant and octant
// sizes, the number of ladder thresholds, and the threshold values themselves. This width
// is expected to change, so nothing may hard-code 64 or the eight constants.
//
// At the Table 4.5 width of 6 bits this happens to be 64 bins, coinciding numerically
// with the tower count. That is a coincidence and not a relationship: the origins differ
// too, and at any other width so do the bin counts.
constexpr unsigned int met_phi_bit_length_ = 6;                         // TOB phi_miss field width
constexpr unsigned int met_phi_range_      = 1u << met_phi_bit_length_; // bins spanning [0, 2*pi)
constexpr unsigned int met_phi_quadrant_   = met_phi_range_ / 4;        // bins in pi/2
constexpr unsigned int met_phi_octant_     = met_phi_range_ / 8;        // ladder thresholds needed

// Fixed-point scale of the tangent comparison, the firmware's Q10. The comparison
// minor * scale >= major * T[k] is integer at any width.
//
// THIS MUST GROW WITH met_phi_bit_length_. The thresholds are themselves rounded at this
// scale, and finer bins mean smaller, more closely spaced thresholds -- at 6 bits the
// smallest is 50, at 10 bits it would be 3. Leaving the scale at Q10 while widening the
// field spends bits on a resolution the ladder cannot actually deliver. Max deviation
// from ideal atan2 binning, in bins, where 0.5 is perfect nearest-bin behaviour:
//
//        phi bits |     Q10     Q12     Q14     Q16
//               6 |  0.5045  0.5010  0.5001  0.5000   <- Q10 is the current firmware
//               7 |  0.5091  0.5022  0.5005  0.5000
//               8 |  0.5169  0.5037  0.5011  0.5001
//               9 |  0.5371  0.5083  0.5021  0.5005
//              10 |  0.5754  0.5187  0.5042  0.5010
//
// So roughly scale >= phi bits + 4 to stay within 0.501 bins, + 6 to be exact. Keep it at
// 10 while the firmware is at 10: matching it matters more than being closer to ideal.
constexpr unsigned int met_phi_tan_scale_bit_length_ = 10;

// T[k] = round(tan((k + 0.5) * 2*pi / met_phi_range_) * 2^10), the tangent of the
// boundary between bin k and bin k+1 within the first octant.
//
// Generated rather than transcribed, for the same reason makeSinLUT() is. At the 6-bit
// width it reproduces the firmware's {50, 152, 256, 366, 484, 614, 759, 928} exactly
// (verified entry by entry); at any other width it follows instead of silently being the
// wrong table for the field.
inline std::array<unsigned int, met_phi_octant_> makeMetPhiTanThresholds() {
    std::array<unsigned int, met_phi_octant_> thresholds{};
    for (unsigned int k = 0; k < met_phi_octant_; ++k) {
        const double angle = (k + 0.5) * (2.0 * M_PI / double(met_phi_range_));
        thresholds[k] = static_cast<unsigned int>(
            std::lround(std::tan(angle) * double(1u << met_phi_tan_scale_bit_length_)));
    }
    return thresholds;
}
static const std::array<unsigned int, met_phi_octant_> met_phi_tan_thresholds_ = makeMetPhiTanThresholds();

inline unsigned int metPhiIndex(int ex, int ey) {
    const bool signX = (ex < 0);
    const bool signY = (ey < 0);

    // 64-bit throughout: the accumulator is 27 bits wide in firmware and the ladder
    // scales the minor component by 2^met_phi_tan_scale_bit_length_, so the comparison
    // needs 27 + 10 = 37 bits at the present settings and more if either grows.
    const unsigned long long absX = static_cast<unsigned long long>(std::llabs(static_cast<long long>(ex)));
    const unsigned long long absY = static_cast<unsigned long long>(std::llabs(static_cast<long long>(ey)));

    const bool swap = (absY > absX);
    const unsigned long long majorVal = swap ? absY : absX;
    const unsigned long long minorVal = swap ? absX : absY;
    if (majorVal == 0ull) return 0u; // zero vector -> bin 0

    // Thermometer code: the thresholds ascend, so the number of them passed IS the
    // in-octant bin index (the firmware's popcount over the thermometer bits).
    const unsigned long long minorScaled = minorVal << met_phi_tan_scale_bit_length_;
    unsigned int base = 0;
    for (unsigned int k = 0; k < met_phi_tan_thresholds_.size(); ++k)
        if (minorScaled >= majorVal * met_phi_tan_thresholds_[k]) ++base;

    // Fold the octant back out, then the quadrant. The wrap at the end is what turns a
    // full turn back into bin 0; it is a wrap around the azimuth range, which the firmware
    // gets for free by truncating to the field width because that range is a power of two.
    const unsigned int phiAbs = swap ? (met_phi_quadrant_ - base) : base;
    unsigned int phi;
    if      (!signX && !signY) phi = phiAbs;                          // quadrant I
    else if ( signX && !signY) phi = 2 * met_phi_quadrant_ - phiAbs;  // quadrant II
    else if ( signX &&  signY) phi = 2 * met_phi_quadrant_ + phiAbs;  // quadrant III
    else                       phi = 4 * met_phi_quadrant_ - phiAbs;  // quadrant IV
    return phi % met_phi_range_;
}

// Bin centre in radians, folded into (-pi, pi] to match every other phi in this chain.
// Bin k is centred on k * 2*pi/met_phi_range_ measured from +x, since the bin EDGES are
// the tangent midpoints at (k +- 0.5) * 2*pi/met_phi_range_.
inline double undigitize_met_phi(unsigned int phiIndex) {
    const double phi = (phiIndex % met_phi_range_) * (2.0 * M_PI / met_phi_range_);
    return (phi > M_PI) ? (phi - 2.0 * M_PI) : phi;
}

// --- square root, MET_LUT_SQRT.v -----------------------------------------------------
// Normalized-LUT approximation to floor(sqrt(radicand)) over a 26-bit radicand. NOT the
// exact integer root: this is the firmware's distributed-ROM estimator, and it is up to
// 3 output LSBs (0.75 GeV at et_granularity_) BELOW or above the exact value. That
// deviation is a property of the algorithm, not an error to be corrected -- correcting it
// here would put the emulator out of agreement with the firmware, which is the one thing
// this file exists to prevent. Anything comparing MET against a floating-point
// sqrt(Ex^2 + Ey^2) has to allow for it.
//
// What "normalized" means here: any radicand can be written as 2^exponent * f, with f
// somewhere in [1, 2). Taking the root of f is the only hard part, since the 2^exponent
// contributes a plain shift by half the exponent. So the ROM tabulates sqrt(f), and the
// LUT INDEX is simply f's leading fractional bits -- the bits sitting immediately below
// the radicand's leading one.
//
// The exponent's PARITY is folded in as the ROM's top address bit: an odd exponent leaves
// half a shift over, worth a factor sqrt(2), so the table stores that half-shift already
// applied and no multiplier is needed.
//
// The path is therefore: find the leading one, take the index bits below it, look up a
// Q1.13 coefficient at {parity, index}, shift by half the exponent, drop the fractional
// bits.
constexpr unsigned int met_sqrt_lut_index_bit_length_ = 9;
constexpr unsigned int met_sqrt_frac_bit_length_      = 13;  // Q1.13 ROM coefficients
constexpr unsigned int met_sqrt_coeff_bit_length_     = 14;
constexpr unsigned int met_sqrt_radicand_bit_length_  = 26;
constexpr unsigned int met_sqrt_scaled_bit_length_    = 26;  // width of scaled_coeff in the RTL
constexpr unsigned int met_sqrt_lut_index_count_      = 1u << met_sqrt_lut_index_bit_length_; // 512
constexpr unsigned int met_sqrt_rom_size_             = 2u * met_sqrt_lut_index_count_;       // 1024

// ROM contents, generated rather than transcribed -- for the same reason makeSinLUT()
// generates its table: a 1024-entry literal block cannot be checked by eye and drifts
// silently. This closed form reproduces all 1024 entries of sqrt_rom in MET_LUT_SQRT.v
// exactly (verified entry by entry). The +0.5 puts the entry at the MIDPOINT of the range
// of f the index stands for, rather than at its bottom edge.
inline std::array<unsigned int, met_sqrt_rom_size_> makeMetSqrtLUT() {
    std::array<unsigned int, met_sqrt_rom_size_> lut{};
    for (unsigned int addr = 0; addr < met_sqrt_rom_size_; ++addr) {
        const unsigned int parity   = addr >> met_sqrt_lut_index_bit_length_;
        const unsigned int lutIndex = addr & (met_sqrt_lut_index_count_ - 1);
        const double normalized = (1.0 + (lutIndex + 0.5) / double(met_sqrt_lut_index_count_))
                                * (parity ? 2.0 : 1.0);
        unsigned int coeff = static_cast<unsigned int>(
            std::lround(std::sqrt(normalized) * double(1u << met_sqrt_frac_bit_length_)));
        if (coeff > maskN(met_sqrt_coeff_bit_length_)) coeff = maskN(met_sqrt_coeff_bit_length_);
        lut[addr] = coeff;
    }
    return lut;
}
static const std::array<unsigned int, met_sqrt_rom_size_> metSqrtLUT_ = makeMetSqrtLUT();

inline unsigned int metLutSqrt(unsigned long long radicand) {
    if (radicand == 0ull) return 0u;

    // Leading-one position, the RTL's leading_one_pos().
    unsigned int exponent = 0;
    while ((radicand >> (exponent + 1)) != 0ull) ++exponent;

    // The LUT index is the bits immediately below the leading one -- that is, the leading
    // fractional bits of the radicand once it has been normalized into [1, 2). They are
    // zero-padded on the right when the radicand is too small to supply them all.
    // Exponent 0 yields 0, matching the RTL's unlisted case falling through to its default.
    unsigned int lutIndex;
    if (exponent >= met_sqrt_lut_index_bit_length_)
        lutIndex = static_cast<unsigned int>((radicand >> (exponent - met_sqrt_lut_index_bit_length_))
                                             & (met_sqrt_lut_index_count_ - 1));
    else
        lutIndex = static_cast<unsigned int>((radicand << (met_sqrt_lut_index_bit_length_ - exponent))
                                             & (met_sqrt_lut_index_count_ - 1));

    const unsigned int romAddr    = ((exponent & 1u) << met_sqrt_lut_index_bit_length_) | lutIndex;
    const unsigned int scaleShift = exponent >> 1;

    // The RTL's scaled_coeff is 26 bits wide, so the shift TRUNCATES there. Masking is
    // not cosmetic: without it a large coefficient and shift would carry bits the
    // firmware silently drops.
    const unsigned long long scaledMask = (1ull << met_sqrt_scaled_bit_length_) - 1ull;
    const unsigned long long scaled =
        (static_cast<unsigned long long>(metSqrtLUT_[romAddr]) << scaleShift) & scaledMask;

    return static_cast<unsigned int>((scaled >> met_sqrt_frac_bit_length_) & maskN(et_bit_length_));
}

// --- missing-energy magnitude, MET_SUM_SQUARE.v --------------------------------------
// 13-bit magnitude of the missing-energy vector, saturating at 0x1fff.
//
// Two separate overflow conditions, both from the firmware and both needed: a single
// component of 8192 or more already guarantees the magnitude overflows, and independently
// the sum of squares reaching 2^26 does (that carry bit IS the threshold, since
// sqrt(2^26) = 8192). The firmware squares only the low 13 bits of each magnitude, which
// is harmless precisely because the component test has already fired in that case.
inline unsigned int metIntegerRoot(int ex, int ey, bool& overflow) {
    const unsigned long long absX = static_cast<unsigned long long>(std::llabs(static_cast<long long>(ex)));
    const unsigned long long absY = static_cast<unsigned long long>(std::llabs(static_cast<long long>(ey)));

    const bool componentOverflow = ((absX >> et_bit_length_) != 0ull) ||
                                   ((absY >> et_bit_length_) != 0ull);

    const unsigned long long exLow = absX & maskN(et_bit_length_);
    const unsigned long long eyLow = absY & maskN(et_bit_length_);
    const unsigned long long squareSum = exLow * exLow + eyLow * eyLow;

    overflow = componentOverflow || ((squareSum >> met_sqrt_radicand_bit_length_) != 0ull);
    if (overflow) return maskN(et_bit_length_); // 0x1fff

    return metLutSqrt(squareSum);
}

// --- saturating TOB field packers ----------------------------------------------------
// The firmware SATURATES where pack_signed_et in metEmulation.cc masks. The difference
// only shows on values that do not fit, but there it is the difference between a MET of
// 2047.75 GeV flagged as overflowing and a MET that has silently wrapped to near zero,
// so these are what the TOB writer uses.
constexpr unsigned int met_tob_signed_mag_bit_length_ = signed_et_bit_length_ - 1; // 12

// 13-bit sign-magnitude {sign, magnitude[11:0]}, saturating at 0xfff.
// The sign bit is forced low when the magnitude is zero, so there is no negative zero --
// the firmware's `sign && (|magnitude)` term. Without it a small negative component that
// saturates to zero would encode as 0x1000 and read back as -0, which is a bit pattern
// the firmware can never emit.
inline unsigned int metTobSignedEt(int value, bool& overflow) {
    const unsigned long long mag = static_cast<unsigned long long>(std::llabs(static_cast<long long>(value)));
    overflow = ((mag >> met_tob_signed_mag_bit_length_) != 0ull);
    const unsigned int magOut = overflow ? maskN(met_tob_signed_mag_bit_length_)
                                         : static_cast<unsigned int>(mag);
    const unsigned int sign = (value < 0 && magOut != 0u) ? 1u : 0u;
    return (sign << met_tob_signed_mag_bit_length_) | magOut;
}

// 13-bit two's complement, saturating SYMMETRICALLY at +-4095 (0x0fff / 0x1001).
//
// This is the encoding the MET TOB Ex/Ey fields are moving to; MET_Engine.v still builds
// sign-magnitude at the time of writing, so the validation path leads the firmware here.
//
// Built on metTobSignedEt rather than alongside it, so the overflow condition and the
// saturated magnitude are identical by construction: overflow is still |value| >= 4096,
// and the only difference is how the saturated magnitude is laid out in the 13 bits.
// The consequence is that -4096 (0x1000), which two's complement could represent, is
// never emitted -- a component of exactly -4096 saturates to -4095 and raises the flag,
// the same as it does in sign-magnitude. On the firmware side that is the existing
// saturated ex_magnitude negated, with no change to ex_overflow.
inline unsigned int metTobTwosComplementEt(int value, bool& overflow) {
    const unsigned int magOut = metTobSignedEt(value, overflow) & maskN(met_tob_signed_mag_bit_length_);
    const long long signedOut = (value < 0) ? -static_cast<long long>(magOut)
                                            :  static_cast<long long>(magOut);
    return static_cast<unsigned int>(signedOut) & maskN(signed_et_bit_length_);
}

// Inverse of metTobTwosComplementEt: sign-extend a 13-bit two's complement field.
inline int metTobTwosComplementToInt(unsigned int field) {
    field &= maskN(signed_et_bit_length_);
    return ((field >> met_tob_signed_mag_bit_length_) & 1u)
               ? static_cast<int>(field) - static_cast<int>(1u << signed_et_bit_length_)
               : static_cast<int>(field);
}

// Unsigned 13-bit scalar energy, saturating at 0x1fff.
inline unsigned int metTobSumEt(unsigned long long sumEt, bool& overflow) {
    overflow = ((sumEt >> et_bit_length_) != 0ull);
    return overflow ? maskN(et_bit_length_) : static_cast<unsigned int>(sumEt);
}

// An accepted input object sits at the top of its E_T field, i.e. its digitized E_T has
// saturated. This is what feeds TOB bit [62] (upstream overflow): the input was already
// clipped before this algorithm ever saw it, so every sum built from it is a lower bound.
inline bool metInputEtSaturated(unsigned int digitizedEt) {
    return digitizedEt >= maskN(et_bit_length_);
}

// Returns output NTuple file name for MET emulation
std::string makeOutputMETFileName(unsigned int maxTowersProcessed,
                                  bool signalBool,
                                  std::string signalString,
                                  bool useSKObjects,
                                  double jetEtThreshold,
                                  double towerEtThreshold,
                                  bool doJetTowerOverlapRemoval,
                                  std::string outputRootFilePath = "/data/larsonma/GEPMET/outputNTuplesDev_METv3/",
                                  bool useEtaSKObjects = false,
                                  double towerScaleFactor = 1.0,
                                  double jetScaleFactor   = 1.0,
                                  unsigned int pileup     = 200,
                                  bool useGEPJwoJ              = false,
                                  double jwojHardEtThreshold   = 0.0,
                                  unsigned int jwojBlockSize   = 1) {
    gSystem->mkdir(outputRootFilePath.c_str());
    std::string usePUSuppress = useEtaSKObjects ? "EtaSK" : (useSKObjects ? "SK" : "NoSK");
    std::string overlapTag    = doJetTowerOverlapRemoval ? "OR" : "NoOR";
    // Format thresholds as integer if whole number, e.g. 20.0 -> "20"
    auto formatThreshold = [](double thr) -> std::string {
        std::ostringstream ss;
        int thrInt = static_cast<int>(thr);
        if (thr == thrInt) ss << thrInt;
        else ss << std::fixed << std::setprecision(1) << thr;
        return ss.str();
    };
    // Format scale factor: 1.0 -> "1", 0.4 -> "0p4", 1.2 -> "1p2", 0.25 -> "0p25".
    // Default stream formatting (%g, 6 significant digits), not a fixed one decimal: fixed
    // precision rounded 0.25 to "0p2", so twrSF 0.2 and 0.25 wrote the same output file.
    // Must match fmt_scale_factor_tag in condor/submit_met_emulation.py.
    auto formatScaleFactor = [](double sf) -> std::string {
        std::ostringstream ss;
        int sfInt = static_cast<int>(sf);
        if (sf == sfInt) { ss << sfInt; return ss.str(); }
        ss << sf;
        std::string s = ss.str();
        for (char& c : s) if (c == '.') c = 'p';
        return s;
    };
    std::string jetThrTag   = formatThreshold(jetEtThreshold);
    std::string towerThrTag = formatThreshold(towerEtThreshold);
    std::string twrSFTag    = formatScaleFactor(towerScaleFactor);
    std::string jetSFTag    = formatScaleFactor(jetScaleFactor);

    std::ostringstream ss;
    if (signalBool) {
        if      (signalString == "VBF_hh_bbbb_cvv0") ss << outputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv0_e8557_s4422_r16130_";
        else if (signalString == "VBF_hh_bbbb_cvv1") ss << outputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv1_e8557_s4422_r16130_";
        else if (signalString == "ggF_hh_bbbb")      ss << outputRootFilePath << "mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_";
        else if (signalString == "ZvvHbb")            ss << outputRootFilePath << "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_";
        else if (signalString == "ttbar_had")         ss << outputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422_r16130_";
        else if (signalString == "Zprime_ttbar")      ss << outputRootFilePath << "mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422_r16130_";
        else if (signalString == "ttbar_semilep")     ss << outputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_semilep_e8557_s4422_r16130_";
        else if (signalString == "ttbar_dilep")       ss << outputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_dilep_e8557_s4422_r16130_";
        else if (signalString == "Zmumu")             ss << outputRootFilePath << "mc21_14TeV_Zmumu_e8557_s4422_r16130_";
    } else {
        ss << outputRootFilePath << "mc21_14TeV_jj_JZ_e8557_s4422_r16130_";
    }
    ss << "N_Towers_" << maxTowersProcessed
       << "_jetEt"    << jetThrTag
       << "_towerEt"  << towerThrTag
       << "_"         << usePUSuppress
       << "_"         << overlapTag
       << "_twrSF"    << twrSFTag
       << "_jetSF"    << jetSFTag;
    // GEP JwoJ tag, written ONLY when that algorithm is enabled, so every filename produced
    // before it existed comes out byte-identical and nothing already on disk has to move.
    // The hard/soft coefficients are not repeated here: in JwoJ mode they ARE twrSF (soft)
    // and jetSF (hard), so the two tags above already carry them. The hard-term threshold is
    // the only genuinely new number, so it is the only new tag.
    // metAnalysisAndRates.C keys its GEP JwoJ path on the "_GEPJwoJ_" substring.
    if (useGEPJwoJ) ss << "_GEPJwoJ_hardEt" << formatThreshold(jwojHardEtThreshold);
    // Block size tag, appended only when it is NOT the single-tower default, so a 1x1 output
    // keeps the name it had before blocks existed -- and, since 1x1 is bit-identical to that
    // algorithm, a file already on disk under the old name is still the right file.
    if (useGEPJwoJ && jwojBlockSize != 1) ss << "_blk" << jwojBlockSize << "x" << jwojBlockSize;
    ss << ".root";
    return applyPileupTags(ss.str(), pileup);
}