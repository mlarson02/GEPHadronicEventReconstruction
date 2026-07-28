#ifndef CONSTANTS_ADV_H
#define CONSTANTS_ADV_H
// Constants used by SW & FW implementation

#define UNROLLFACTOR 16
#define PIPELINEII 3

constexpr unsigned int nTotalSeeds_ = 10;
constexpr unsigned int nSeedsInput_ = 6;
constexpr unsigned int nSeedsOutput_ = 2;
constexpr unsigned int maxObjectsConsidered_ = 128;
constexpr double et_granularity_ = 0.125;
constexpr unsigned int subjet_et_threshold_ = 200;
constexpr double r2Cut_ = 1.21;
constexpr double rCut_ = 1.1;
constexpr double rMergeCut_ = 2.0;
constexpr unsigned int et_bit_length_ = 13;
constexpr unsigned int eta_bit_length_ = 7;
constexpr unsigned int phi_bit_length_ = 6;
constexpr unsigned int eta_range_ = 98;
constexpr unsigned int num_subjets_length_ = 2;
constexpr unsigned int deltaRBits_ = 8;
constexpr double phi_min_ = -3.2;
constexpr double phi_max_ = 3.2;
constexpr unsigned int pi_digitized_in_phi_ = 31;
constexpr double eta_min_ = -4.85;
constexpr double eta_max_ = 4.95;
constexpr double eta_granularity_ = 0.1;
constexpr double phi_granularity_ = 0.1;
constexpr unsigned int et_min_ = 0;
constexpr unsigned int et_max_ = 1024;
constexpr double phi_range_ = 6.4;


constexpr unsigned int padded_zeroes_length_ = 64 - et_bit_length_ - eta_bit_length_ - phi_bit_length_ - num_subjets_length_ - num_subjets_length_;
constexpr unsigned int padded_zeroes_length_32b_ = 32 - et_bit_length_ - eta_bit_length_ - phi_bit_length_;
constexpr unsigned int total_bits_input_ = padded_zeroes_length_32b_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
constexpr unsigned int total_bits_output_ = padded_zeroes_length_ + num_subjets_length_ + num_subjets_length_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
typedef ap_uint<total_bits_input_> input; // need 32b input, 64b output!
typedef ap_uint<total_bits_output_> output;

// MSB -> LSB word order is phi | eta | et, i.e. et occupies the LSBs and phi the MSBs
constexpr unsigned int et_low_   = 0;
constexpr unsigned int et_high_  = et_low_ + et_bit_length_ - 1;

constexpr unsigned int eta_low_  = et_high_ + 1;
constexpr unsigned int eta_high_ = eta_low_ + eta_bit_length_ - 1;

constexpr unsigned int phi_low_  = eta_high_ + 1;
constexpr unsigned int phi_high_ = phi_low_ + phi_bit_length_ - 1;

constexpr unsigned int num_subjets_low_  = phi_high_ + 1;
constexpr unsigned int num_subjets_high_ = num_subjets_low_ + num_subjets_length_ - 1;

constexpr unsigned int padded_zeroes_low_  = num_subjets_high_ + 1;
constexpr unsigned int padded_zeroes_high_ = padded_zeroes_low_ + padded_zeroes_length_ - 1;

constexpr unsigned int nSeedsDeltaR_ = nSeedsInput_ - nSeedsOutput_;

constexpr unsigned int deltaR2_bits_  = 10; // bit width of the saturated deltaR2 result (values > 1023 clamp to 1023)
constexpr double deltaR2_granularity_ = eta_granularity_ * eta_granularity_; // FIXME THIS SHOULD BE EQUIVALENT TO SQUARING PHI_GRANULARITY_ - maybe add an exception if they are not the same

constexpr unsigned int digitized_delta_R2Cut_ = static_cast<unsigned int>(r2Cut_/deltaR2_granularity_ + 0.5); //+ 0.5 for correct rounding

constexpr unsigned int digitized_d_search_squared_ = static_cast<unsigned int>(((rMergeCut_) * (rMergeCut_)/(deltaR2_granularity_)) + 0.5);

// Controls whether calcDeltaR2 (helperFunctions_adv.h) computes deltaR^2 via DSP multipliers (baseline)
// or via a precomputed LUT indexed by (deltaEta, deltaPhi). Set to 0 to use the LUT path instead.
#define USE_DSPS_ 1

#if !USE_DSPS_
// Boolean pass/fail against r2Cut_, indexed by (deltaEta, deltaPhi). Used by passesRCut().
constexpr unsigned int max_R2lut_size_ = 704;
static const bool lut_[max_R2lut_size_] =
#include "/home/mlarson/GEPHadronicEventReconstruction/algorithm/emulation/LUT_Constants_Generation/LUTs/v3/LUT_deltaR2Cut_rMerge_2_R2_1.21.h"
;

// Digitized deltaR (8-bit, full grid coverage, cut-independent) indexed by (deltaEta, deltaPhi).
// Used by passesTwoRCut() and passesSearchRadius() against the digitized thresholds below.
constexpr unsigned int max_Rlut_size_ = 1280;
static const ap_uint<deltaRBits_> lutR_[max_Rlut_size_] =
#include "/home/mlarson/GEPHadronicEventReconstruction/algorithm/emulation/LUT_Constants_Generation/LUTs/v3/LUT_deltaR_rMerge_2_R2_1.21.h"
;

constexpr double deltaR_granularity_ = 0.03993458703717951;
constexpr unsigned int digitized_two_rCut_ = 55;
constexpr unsigned int digitized_rMergeCut_ = 50;
#endif

        
#endif // CONSTANTS_ADV_H
