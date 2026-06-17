// Constants used by MET emulation

#include <cmath>

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
constexpr double phi_min_ = -3.2;
constexpr double phi_max_ = 3.2;
constexpr double eta_min_ = -4.9;
constexpr double eta_max_ = 4.9;
constexpr double rCut_ = 0.4;
constexpr double r2Cut_ = rCut_ * rCut_; 
constexpr unsigned int et_min_ = 0;
constexpr unsigned int et_max_ = 1024;
constexpr double phi_granularity_ = (phi_max_ - phi_min_) / (1 << (phi_bit_length_));
constexpr unsigned int pi_digitized_in_phi_ = (M_PI) / (phi_granularity_);
constexpr unsigned int half_pi_digitized_in_phi_     =  pi_digitized_in_phi_/2;
constexpr unsigned int two_pi_digitized_in_phi_      =  (1 << phi_bit_length_) - 1; 
constexpr unsigned int eta_range_ = (eta_max_ - eta_min_) / (phi_granularity_); // ensures phi and eta have the same granularity
constexpr double eta_granularity_ = (eta_max_ - eta_min_) / double(eta_range_);
constexpr double et_granularity_ = (et_max_ - et_min_) / double((1 << et_bit_length_));
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


static const int sinLUT_[1 << phi_bit_length_] = {
    0, -401, -799, -1189, -1567, -1930, -2275, -2598,
    -2896, -3165, -3405, -3611, -3783, -3919, -4016, -4075,
    -4095, -4075, -4016, -3919, -3783, -3611, -3405, -3165,
    -2896, -2598, -2275, -1930, -1567, -1189, -799, -401,
    0, 401, 799, 1189, 1567, 1930, 2275, 2598,
    2896, 3165, 3405, 3611, 3783, 3919, 4016, 4075,
    4095, 4075, 4016, 3919, 3783, 3611, 3405, 3165,
    2896, 2598, 2275, 1930, 1567, 1189, 799, 401
};