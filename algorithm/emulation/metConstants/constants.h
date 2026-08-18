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
// between two, and *_max_ is one LSB past the last centre.
constexpr unsigned int eta_range_ = 98;
constexpr unsigned int phi_range_ = 64;
constexpr double eta_granularity_ = 0.1;
constexpr double eta_min_ = -4.85;
constexpr double eta_max_ = eta_min_ + eta_range_ * eta_granularity_;
constexpr double phi_granularity_ = (2 * M_PI) / double(phi_range_); // pi/32
constexpr double phi_min_ = -M_PI + phi_granularity_ / 2;
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
// Two known divergences from the firmware trig module (trig.v), both to be
// closed on the firmware side -- run dumpSinLUT.cc for the numbers:
//
//  1. Index origin. trig.v reads the 6-bit phi code as *signed* (o_sin's sign bit
//     is taken straight from i_index[5]), i.e. it assumes code 0 is phi = +pi/64
//     and code 32 is phi = -pi + pi/64. This grid is unsigned from -pi, so the two
//     numberings differ by a pi rotation and trig.v needs its index flipped
//     (i_index ^ 6'd32) to line up. Left this way deliberately: the LRJ chain only
//     ever uses deltaPhi, which a global rotation leaves untouched.
//  2. Amplitude. The table currently in trig.v is round-to-nearest at an effective
//     amplitude of ~4094.3, not 4096 -- 0.042% low, the signature of an
//     uncompensated CORDIC gain. No scale reproduces it with truncation. Its
//     entries therefore sit 0 to 2 LSB below this table until it is regenerated
//     from round(4096 * sin).
inline std::array<int, phi_range_> makeSinLUT() {
    std::array<int, phi_range_> lut{};
    for (unsigned int iPhi = 0; iPhi < phi_range_; ++iPhi) {
        const double phi = phi_min_ + iPhi * phi_granularity_;
        lut[iPhi] = static_cast<int>(std::lround(std::sin(phi) * (1 << (sin_bit_length_ - 1))));
    }
    return lut;
}
static const std::array<int, phi_range_> sinLUT_ = makeSinLUT();