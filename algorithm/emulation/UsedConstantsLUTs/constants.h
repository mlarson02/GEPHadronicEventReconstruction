// Constants used by JetTagger Emulation

#include <cmath>
static inline uint32_t maskN(unsigned n) { return (n >= 32) ? 0xFFFFFFFFu : ((1u << n) - 1u); }
constexpr unsigned int algoVersion_ = 3;
constexpr unsigned int nSeedsInput_ = 10;
constexpr unsigned int nProtoSeeds_ = 6;
constexpr unsigned int nSeedsOutput_ = 2;
constexpr unsigned int max_R_8b_lut_size_ = 1408;
constexpr unsigned int maxObjectsConsidered_ = 512;
constexpr double r2Cut_ = 1.2100;
constexpr double rCut_ = 1.1000;
constexpr double rMergeCut_ = 2.0000;
constexpr unsigned int et_bit_length_ = 13;
constexpr unsigned int eta_bit_length_ = 7;
constexpr unsigned int phi_bit_length_ = 6;
constexpr unsigned int num_subjets_length_ = 2;
constexpr unsigned int mass_approx_bit_length_ = 8;
constexpr unsigned int N_subjetiness_bit_length_ = 8;
constexpr unsigned int psi_R_bit_length_ = 8;
constexpr unsigned int deltaR_lut_length_ = 8;
constexpr unsigned int padded_zeroes_length_ = 64 - et_bit_length_ - eta_bit_length_ - phi_bit_length_ - num_subjets_length_;
constexpr unsigned int total_bits_output_ = padded_zeroes_length_ + num_subjets_length_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
constexpr double phi_min_ = -3.2;
constexpr double phi_max_ = 3.2;
constexpr double eta_min_ = -4.85;
constexpr double eta_max_ = 4.95;
constexpr unsigned int et_min_ = 0;
constexpr unsigned int et_max_ = 1024;
constexpr double phi_granularity_ = (phi_max_ - phi_min_) / (1 << (phi_bit_length_));
constexpr unsigned int pi_digitized_in_phi_ = (M_PI) / (phi_granularity_);
const int PI_D     =  int(pi_digitized_in_phi_);
const int TWO_PI_D =  2 * pi_digitized_in_phi_ + 1;
constexpr unsigned int eta_range_ = (eta_max_ - eta_min_) / (phi_granularity_);
constexpr double eta_granularity_ = (eta_max_ - eta_min_) / double(eta_range_);
constexpr double et_granularity_ = (et_max_ - et_min_) / double((1 << et_bit_length_));
constexpr double deltaR_granularity_ = (2 * rCut_) / double((1 << deltaR_lut_length_));
constexpr double psi_R_granularity_ = (2 * rCut_) / double((1 << deltaR_lut_length_)) ;
constexpr double massApproxRawLSB_GeV_  = et_granularity_ * deltaR_granularity_ ;
constexpr double massApprox_max_  = 512;
constexpr double massApproxNewLSB_GeV_ = massApprox_max_ / (1 << mass_approx_bit_length_);
constexpr double massApproxScaleFactor_ = massApproxNewLSB_GeV_ / massApproxRawLSB_GeV_;
constexpr unsigned int massApproxDivisor_ = static_cast<unsigned int>(massApproxScaleFactor_ + 0.5);
constexpr double mass_approx_granularity_ = massApprox_max_ / (1 << mass_approx_bit_length_);
constexpr double N_subjetiness_granularity_ = (1.0) / (1 << N_subjetiness_bit_length_) ;
constexpr double phi_range_ = phi_max_ - phi_min_;
constexpr double deltaR2_granularity_ = eta_granularity_ * eta_granularity_;
constexpr unsigned int digitized_delta_R2Cut_ = static_cast<unsigned int>(r2Cut_/deltaR2_granularity_ + 0.5);
constexpr unsigned int digitized_d_search_squared_ = static_cast<unsigned int>(((rMergeCut_) * (rMergeCut_))/deltaR2_granularity_ + 0.5);
constexpr unsigned int total_bits_ = padded_zeroes_length_ + num_subjets_length_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
// MSB -> LSB word order is substruct_0 | phi | eta | et, i.e. et occupies the LSBs
constexpr unsigned int et_low_   = 0;
constexpr unsigned int et_high_  = et_low_ + et_bit_length_ - 1;
constexpr unsigned int eta_low_  = et_high_ + 1;
constexpr unsigned int eta_high_ = eta_low_ + eta_bit_length_ - 1;
constexpr unsigned int phi_low_  = eta_high_ + 1;
constexpr unsigned int phi_high_ = phi_low_ + phi_bit_length_ - 1;
constexpr unsigned int substruct_0_bit_low_  = phi_high_ + 1;
constexpr unsigned int substruct_0_bit_high_ = substruct_0_bit_low_ + num_subjets_length_ - 1;
constexpr unsigned int padded_zeroes_low_  = substruct_0_bit_high_ + 1;
constexpr unsigned int padded_zeroes_high_ = padded_zeroes_low_ + padded_zeroes_length_ - 1;
constexpr unsigned int nPreSeeds_ = nProtoSeeds_ - nSeedsOutput_;
static const unsigned int lutR_8b_[max_R_8b_lut_size_] = 
#include "deltaRLUT_8b.h"
;

