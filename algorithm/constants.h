#ifndef CONSTANTS_H
#define CONSTANTS_H
// Constants used by SW & FW implementation

#include <cmath>

#define UNROLLFACTOR 1
#define PIPELINEII 3

#define sortSeeds_ false
#define iterative_ false
constexpr unsigned int nTotalSeeds_ = 200;
constexpr unsigned int nSeedsInput_ = 6;
constexpr unsigned int nSeedsOutput_ = 2;
constexpr unsigned int maxObjectsConsidered_ = 10;
constexpr unsigned int inputEnergyCut_ = 1;
#define useInputEnergyCut_ false
constexpr double et_granularity_ = 0.25; // 250 MeV LSB = et_max_ / (1 << et_bit_length_)
constexpr double r2Cut_ = 1.21;
constexpr double rCut_ = 1.1;
constexpr double rMergeCut_ = 0.001;
constexpr unsigned int et_bit_length_ = 13;
constexpr unsigned int eta_bit_length_ = 10;
constexpr unsigned int eta_range_ = 98;
constexpr unsigned int phi_bit_length_ = 9;
// ---- digitization grid ----
// Identical physical grid to the advanced version (constants_adv.h): the GEP
// tower grid, 98 eta towers of 0.1 spanning |eta| < 4.9 and 64 phi towers of
// pi/32 covering the full 2*pi, matching Athena's
// CaloTowerContainer::configureGrid(98, -4.9, 4.9, 64). Only the field widths
// differ here (10b eta / 9b phi, the standard TOB format), and they are just
// field widths: the dynamic range and the granularity come from the code counts
// eta_range_ / phi_range_, and the high bits the grid never reaches are padding.
// eta_min_ / phi_min_ are the first tower centre, *_max_ one LSB past the last.
constexpr unsigned int phi_range_ = 64;
constexpr double eta_granularity_ = 0.1;
constexpr double eta_min_ = -4.85;
constexpr double eta_max_ = eta_min_ + eta_range_ * eta_granularity_; // 4.95
constexpr double phi_granularity_ = (2 * M_PI) / double(phi_range_); // pi/32
constexpr double phi_min_ = -M_PI + phi_granularity_ / 2;
constexpr double phi_max_ = phi_min_ + phi_range_ * phi_granularity_;
constexpr unsigned int pi_digitized_in_phi_ = phi_range_ / 2; // 32
constexpr double deltaR2_granularity_ = eta_granularity_ * eta_granularity_;
constexpr unsigned int et_min_ = 0;
constexpr unsigned int et_max_ = 2048;
#define useMax_ false
constexpr unsigned int max_R2lut_size_ = 45056;
constexpr unsigned int max_Rlut_size_ = 1;
constexpr double deltaR_max_ = 10.48187;
constexpr unsigned int deltaR_bits_ = 8;
constexpr unsigned int max_R_8b_lut_size_ = 45056;



const unsigned int lut_size_ = (eta_range_ * (phi_range_ / 2)); // rows of |deltaPhi| codes below pi
#if !WRITE_LUT
constexpr unsigned int padded_zeroes_length_ = 64 - et_bit_length_ - eta_bit_length_ - phi_bit_length_;
constexpr unsigned int padded_zeroes_length_32b_ = 128 - et_bit_length_ - eta_bit_length_ - phi_bit_length_;
constexpr unsigned int total_bits_input_ = padded_zeroes_length_32b_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
constexpr unsigned int total_bits_output_ = padded_zeroes_length_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
typedef ap_uint<total_bits_input_> input; // need 32b input, 64b output!
typedef ap_uint<total_bits_output_> output;

// MSB -> LSB word order is phi | eta | et, i.e. et occupies the LSBs and phi the MSBs
constexpr unsigned int et_low_   = 0;
constexpr unsigned int et_high_  = et_low_ + et_bit_length_ - 1;

constexpr unsigned int eta_low_  = et_high_ + 1;
constexpr unsigned int eta_high_ = eta_low_ + eta_bit_length_ - 1;

constexpr unsigned int phi_low_  = eta_high_ + 1;
constexpr unsigned int phi_high_ = phi_low_ + phi_bit_length_ - 1;

constexpr unsigned int padded_zeroes_low_  = phi_high_ + 1;
constexpr unsigned int padded_zeroes_high_ = padded_zeroes_low_ + padded_zeroes_length_ - 1;


constexpr unsigned int nSeedsDeltaR_ = nSeedsInput_ - nSeedsOutput_;

constexpr unsigned int digitized_delta_R2_ = static_cast<unsigned int>(r2Cut_/deltaR2_granularity_ + 0.5);

// NOTE: the basic algorithm computes deltaR^2 directly with DSP multipliers (calcDeltaR2 in
// helperFunctions.h), so the lut_/lutR_/lutR_8b_ arrays are never read. They were dead code (and
// caused a size mismatch vs the generated LUT files), so they have been removed.

#endif
constexpr unsigned int deltaR_levels_ = (1 << deltaR_bits_); // 256
constexpr float deltaR_step_ = deltaR_max_ / (deltaR_levels_ - 1); // ~0.041
constexpr unsigned int rMergeConsiderCutDigitized_ = (rMergeCut_) / deltaR_step_;
        
#endif // CONSTANTS_H
