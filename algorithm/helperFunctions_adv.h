#ifndef HELPER_FUNCTIONS_ADV_H
#define HELPER_FUNCTIONS_ADV_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <array>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "ap_int.h"
#include "constants_adv.h"

// Shortest-arc absolute (deltaEta, deltaPhi) between two digitized (eta, phi) points.
inline void calcDeltaEtaPhi(
    ap_uint<eta_bit_length_> eta1, ap_uint<phi_bit_length_> phi1,
    ap_uint<eta_bit_length_> eta2, ap_uint<phi_bit_length_> phi2,
    ap_uint<eta_bit_length_>& uDEta, ap_uint<phi_bit_length_ - 1>& corrDPhi
) {
    #pragma HLS INLINE
    ap_int<eta_bit_length_ + 1> dEta = eta1 - eta2;
    ap_int<phi_bit_length_ + 1> dPhi = phi1 - phi2;
    uDEta = dEta[eta_bit_length_] ? static_cast<ap_uint<eta_bit_length_>>(-dEta) : static_cast<ap_uint<eta_bit_length_>>(dEta);
    ap_uint<phi_bit_length_> uDPhi = dPhi[phi_bit_length_] ? static_cast<ap_uint<phi_bit_length_>>(-dPhi) : static_cast<ap_uint<phi_bit_length_>>(dPhi);
    if (uDPhi >= pi_digitized_in_phi_) uDPhi = 2 * pi_digitized_in_phi_ - uDPhi;
    corrDPhi = uDPhi;
}

#if USE_DSPS_

// Compute saturated digitized deltaR^2 between two (eta, phi) points via DSP multipliers.
inline ap_uint<deltaR2_bits_> calcDeltaR2(
    ap_uint<eta_bit_length_> eta1, ap_uint<phi_bit_length_> phi1,
    ap_uint<eta_bit_length_> eta2, ap_uint<phi_bit_length_> phi2
) {
    #pragma HLS INLINE
    ap_uint<eta_bit_length_> uDEta;
    ap_uint<phi_bit_length_ - 1> corrDPhi;
    calcDeltaEtaPhi(eta1, phi1, eta2, phi2, uDEta, corrDPhi);
    ap_uint<2*eta_bit_length_>     etaSq    = uDEta * uDEta;
    #pragma HLS bind_op variable=etaSq op=mul impl=dsp
    // The phi square is only a 5-bit x 5-bit multiply -- a whole DSP48 per call is wasteful, and there's
    // ample LUT headroom, so map it to fabric. Halves DSP usage (one DSP per calcDeltaR2 instead of two)
    // while keeping the eta square on a genuine DSP. Bit-exact -- same arithmetic, different resource.
    ap_uint<2*(phi_bit_length_-1)> phiSq    = corrDPhi * corrDPhi;
    #pragma HLS bind_op variable=phiSq op=mul impl=fabric
    ap_uint<2*eta_bit_length_+1>   rawSum   = ap_uint<2*eta_bit_length_+1>(etaSq) + phiSq;
    return rawSum > ap_uint<deltaR2_bits_>((1 << deltaR2_bits_) - 1)
           ? ap_uint<deltaR2_bits_>((1 << deltaR2_bits_) - 1)
           : ap_uint<deltaR2_bits_>(rawSum);
}

// True if deltaR between the two points is <= rCut_ (jet building / subjet counting).
inline bool passesRCut(
    ap_uint<eta_bit_length_> eta1, ap_uint<phi_bit_length_> phi1,
    ap_uint<eta_bit_length_> eta2, ap_uint<phi_bit_length_> phi2
) {
    #pragma HLS INLINE
    return calcDeltaR2(eta1, phi1, eta2, phi2) <= digitized_delta_R2Cut_;
}

// True if deltaR between the two points is <= 2*rCut_ (overlap removal).
inline bool passesTwoRCut(
    ap_uint<eta_bit_length_> eta1, ap_uint<phi_bit_length_> phi1,
    ap_uint<eta_bit_length_> eta2, ap_uint<phi_bit_length_> phi2
) {
    #pragma HLS INLINE
    return calcDeltaR2(eta1, phi1, eta2, phi2) <= 2 * 2 * digitized_delta_R2Cut_;
}

// True if deltaR between the two points is <= rMergeCut_ (seed position search).
inline bool passesSearchRadius(
    ap_uint<eta_bit_length_> eta1, ap_uint<phi_bit_length_> phi1,
    ap_uint<eta_bit_length_> eta2, ap_uint<phi_bit_length_> phi2
) {
    #pragma HLS INLINE
    return calcDeltaR2(eta1, phi1, eta2, phi2) <= digitized_d_search_squared_;
}

#else

// Baseline (LUT) path. lut_/lutR_/max_R2lut_size_/max_Rlut_size_/digitized_two_rCut_/digitized_rMergeCut_
// are all defined in constants_adv.h, generated from LUT_deltaR2Cut_*.h / LUT_deltaR_*.h.

// True if deltaR between the two points is <= rCut_ (jet building / subjet counting).
inline bool passesRCut(
    ap_uint<eta_bit_length_> eta1, ap_uint<phi_bit_length_> phi1,
    ap_uint<eta_bit_length_> eta2, ap_uint<phi_bit_length_> phi2
) {
    #pragma HLS INLINE
    ap_uint<eta_bit_length_> uDEta;
    ap_uint<phi_bit_length_ - 1> corrDPhi;
    calcDeltaEtaPhi(eta1, phi1, eta2, phi2, uDEta, corrDPhi);
    ap_uint<eta_bit_length_ + phi_bit_length_ - 1> lutIndex = uDEta * (1 << (phi_bit_length_ - 1)) + corrDPhi;
    if (lutIndex >= max_R2lut_size_) return false; // out of table -> definitely beyond rCut_
    return lut_[lutIndex];
}

// True if deltaR between the two points is <= 2*rCut_ (overlap removal).
inline bool passesTwoRCut(
    ap_uint<eta_bit_length_> eta1, ap_uint<phi_bit_length_> phi1,
    ap_uint<eta_bit_length_> eta2, ap_uint<phi_bit_length_> phi2
) {
    #pragma HLS INLINE
    ap_uint<eta_bit_length_> uDEta;
    ap_uint<phi_bit_length_ - 1> corrDPhi;
    calcDeltaEtaPhi(eta1, phi1, eta2, phi2, uDEta, corrDPhi);
    ap_uint<eta_bit_length_ + phi_bit_length_ - 1> lutIndex = uDEta * (1 << (phi_bit_length_ - 1)) + corrDPhi;
    if (lutIndex >= max_Rlut_size_) return false; // out of table -> definitely beyond 2*rCut_
    return lutR_[lutIndex] <= digitized_two_rCut_;
}

// True if deltaR between the two points is <= rMergeCut_ (seed position search).
inline bool passesSearchRadius(
    ap_uint<eta_bit_length_> eta1, ap_uint<phi_bit_length_> phi1,
    ap_uint<eta_bit_length_> eta2, ap_uint<phi_bit_length_> phi2
) {
    #pragma HLS INLINE
    ap_uint<eta_bit_length_> uDEta;
    ap_uint<phi_bit_length_ - 1> corrDPhi;
    calcDeltaEtaPhi(eta1, phi1, eta2, phi2, uDEta, corrDPhi);
    ap_uint<eta_bit_length_ + phi_bit_length_ - 1> lutIndex = uDEta * (1 << (phi_bit_length_ - 1)) + corrDPhi;
    if (lutIndex >= max_Rlut_size_) return false; // out of table -> definitely beyond rMergeCut_
    return lutR_[lutIndex] <= digitized_rMergeCut_;
}

#endif

#endif
