// Constants used by MET emulation

#include <array>
#include <cmath>
#include <cstdint>

static inline uint32_t maskN(unsigned n) { return (n >= 32) ? 0xFFFFFFFFu : ((1u << n) - 1u); }
constexpr unsigned int maxTowersConsidered_ = 4096;
constexpr unsigned int maxJetsConsidered_ = 10;
constexpr unsigned int signed_et_bit_length_ = 13;
constexpr unsigned int et_bit_length_ = 13;
constexpr unsigned int eta_bit_length_ = 7;
constexpr unsigned int phi_bit_length_ = 6;
constexpr unsigned int sin_bit_length_ = 13; 
constexpr unsigned int padded_zeroes_length_ = 64 - et_bit_length_ - et_bit_length_ - eta_bit_length_ - phi_bit_length_;
constexpr unsigned int total_bits_output_ = padded_zeroes_length_ + et_bit_length_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
// ---- digitization grid ----
// Same physical grid as the jet tagger chain (see emulationHelperFunctions.h):
// the GEP tower grid, 98 eta towers of 0.1 spanning |eta| < 4.9 and 64 phi
// towers of pi/32 covering the full 2*pi, matching Athena's
// CaloTowerContainer::configureGrid(98, -4.9, 4.9, 64). The dynamic range and
// granularity come from the code counts eta_range_ / phi_range_, not from the
// field widths. eta_min_ / phi_min_ are the first tower *centre*, so a tower
// digitizes to an integer code instead of landing on a round-half tie halfway
// between two, and *_max_ is one LSB past the last centre. Phi codes count from
// phi = 0: index 0 = [0, pi/32).
constexpr unsigned int eta_range_ = 98;
constexpr unsigned int phi_range_ = 64;
constexpr double eta_granularity_ = 0.1;
constexpr double eta_min_ = -4.85;
constexpr double eta_max_ = eta_min_ + eta_range_ * eta_granularity_;
constexpr double phi_granularity_ = (2 * M_PI) / double(phi_range_); // pi/32
constexpr double phi_min_ = phi_granularity_ / 2;
constexpr double phi_max_ = phi_min_ + phi_range_ * phi_granularity_;
constexpr unsigned int pi_digitized_in_phi_          =  phi_range_ / 2;      // 32
constexpr unsigned int half_pi_digitized_in_phi_     =  pi_digitized_in_phi_/2; // 16, exactly pi/2
constexpr unsigned int two_pi_digitized_in_phi_      =  phi_range_ - 1;
constexpr double rCut_ = 0.4;
constexpr double r2Cut_ = rCut_ * rCut_;
constexpr unsigned int et_min_ = 0;
constexpr unsigned int et_max_ = 2048;
constexpr double et_granularity_ = (et_max_ - et_min_) / double((1 << et_bit_length_)); // 0.25 GeV
constexpr unsigned int deltaR2_length_ = 8;
constexpr double deltaR2_granularity_ = eta_granularity_ * eta_granularity_;
constexpr unsigned int digitized_delta_R2Cut_ = static_cast<unsigned int>(r2Cut_/deltaR2_granularity_ + 0.5); // round up
constexpr unsigned int phi_low_  = 0;
constexpr unsigned int phi_high_ = phi_low_ + phi_bit_length_ - 1;
constexpr unsigned int eta_low_  = phi_high_ + 1;
constexpr unsigned int eta_high_ = eta_low_ + eta_bit_length_ - 1;
constexpr unsigned int met_low_   = eta_high_ + 1;
constexpr unsigned int met_high_  = met_low_ + et_bit_length_ - 1;
constexpr unsigned int sumET_low_   = met_high_ + 1;
constexpr unsigned int sumET_high_  = sumET_low_ + et_bit_length_ - 1;


// sin(phi) for every phi code, scaled by 1 << (sin_bit_length_ - 1) = 4096 (the
// sign takes one of the sin_bit_length_ bits, leaving that as the magnitude
// range) and rounded to nearest. Built from the phi grid rather than written out
// by hand so it cannot drift from phi_min_ / phi_granularity_ / phi_range_:
// entry k is the sine of the *centre* of tower k, which is why no entry is
// exactly 0 or exactly full scale -- the tower centres straddle 0 and +-pi/2
// rather than landing on them. pi/2 is exactly half_pi_digitized_in_phi_ = 16
// codes, so cos(phi) is this same table read 16 codes along.
//
// This table AGREED BIT-FOR-BIT with the firmware trig module (MET_B_Trig in
// trig.v) on all 64 codes while the phi origin was -pi. Since the move of the phi
// codes to the [0, 2pi) origin it is the same table rotated by pi_digitized_in_phi_
// codes, and trig.v needs the same rotation to match again. Two divergences were
// recorded here previously and are both closed; they are written out because the
// note is easy to reintroduce by misreading the firmware:
//
//  1. Index origin. trig.v takes its sign from phi[5], so codes 0..31 are the
//     negative-sine half-plane and its code 0 is phi = -pi + pi/64 -- the OLD
//     phi_min_. With phi_min_ = +pi/64 codes 0..31 are the positive-sine
//     half-plane, i.e. trig.v must invert phi[5] (a 32-code rotation).
//  2. Amplitude. The table was said to sit ~0.042% low, at an effective amplitude
//     of ~4094.3, as if from an uncompensated CORDIC gain. It does not: it is
//     round(4096 * sin) exactly. Its maximum entry is 4091 rather than 4096 only
//     because the tower centres straddle pi/2 instead of landing on it -- the same
//     reason no entry here is exactly full scale.
inline std::array<int, phi_range_> makeSinLUT() {
    std::array<int, phi_range_> lut{};
    for (unsigned int iPhi = 0; iPhi < phi_range_; ++iPhi) {
        const double phi = phi_min_ + iPhi * phi_granularity_;
        lut[iPhi] = static_cast<int>(std::lround(std::sin(phi) * (1 << (sin_bit_length_ - 1))));
    }
    return lut;
}
static const std::array<int, phi_range_> sinLUT_ = makeSinLUT();

// ---- half-granularity phi LUT, for even-sized GEP JwoJ blocks -----------------
// An NxN block of towers is positioned at its geometric centre. For ODD N that
// centre is the centre of the middle tower and sinLUT_ above already has it; for
// EVEN N it falls on the CORNER where four towers meet, half a tower off every
// entry sinLUT_ holds. This table supplies those half-tower positions.
//
// The index is a HALF-code h over phi_half_range_ = 2 * phi_range_ = 128 steps of
// phi_granularity_ / 2:
//
//    h = 2k     -> the centre of tower k          (what sinLUT_[k] holds)
//    h = 2k + 1 -> the boundary between towers k and k + 1
//
// so a block spanning tower codes [a, b] sits at h = a + b, integer for every
// block size, and the parity of h alone says whether the centre is on a tower or
// on a corner.
//
// Two properties this construction is chosen for, both checked by dumpSinLUT.cc:
//
//  1. Every EVEN entry is bit-identical to sinLUT_[h / 2]. Same expression, same
//     scale, same rounding, evaluated at the same angle -- h * (granularity / 2)
//     with h = 2k is k * granularity exactly. So an odd-sized block reads the
//     value it would have read from sinLUT_, and only even-sized blocks ever
//     touch an entry that did not exist before. 1x1 output is unchanged.
//  2. pi/2 is exactly half_pi_digitized_in_phi_half_ = 32 half-codes, so cos is
//     this same table read 32 along -- the direct analogue of sinLUT_'s 16, and
//     exact to 0 LSB rather than approximate.
//
// Interpolating sinLUT_ instead was the obvious alternative and is wrong at a
// level that matters: (sin a + sin b) / 2 = sin(midpoint) * cos(pi/64), i.e.
// 0.120% low on every even-sized block, a systematic scale error on the hard term
// rather than a rounding one.
//
// For firmware this is one extra address bit on the existing table, not a second
// table: trig.v would index sin at half-tower resolution and the present
// tower-centre entries remain its even addresses. The two divergences noted above
// (index origin, amplitude) apply here unchanged.
constexpr unsigned int phi_half_range_                = 2 * phi_range_;                // 128
constexpr unsigned int half_pi_digitized_in_phi_half_ = 2 * half_pi_digitized_in_phi_; // 32, exactly pi/2
constexpr unsigned int two_pi_digitized_in_phi_half_  = phi_half_range_ - 1;           // 127
constexpr double phi_half_granularity_                = phi_granularity_ / 2;          // pi/64

inline std::array<int, phi_half_range_> makeSinLUTHalf() {
    std::array<int, phi_half_range_> lut{};
    for (unsigned int iHalf = 0; iHalf < phi_half_range_; ++iHalf) {
        const double phi = phi_min_ + iHalf * phi_half_granularity_;
        lut[iHalf] = static_cast<int>(std::lround(std::sin(phi) * (1 << (sin_bit_length_ - 1))));
    }
    return lut;
}
static const std::array<int, phi_half_range_> sinLUTHalf_ = makeSinLUTHalf();