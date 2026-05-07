#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H
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
#include "constants.h"

// Compute digitized deltaR^2 between two (eta, phi) points.
// Applies shortest-arc phi wrapping. Return type is wide enough to never overflow.
inline ap_uint<2*(eta_bit_length_ + phi_bit_length_)> calcDeltaR2(
    ap_uint<eta_bit_length_> eta1, ap_uint<phi_bit_length_> phi1,
    ap_uint<eta_bit_length_> eta2, ap_uint<phi_bit_length_> phi2
) {
    #pragma HLS INLINE
    ap_int<eta_bit_length_ + 1> dEta = eta1 - eta2;
    ap_int<phi_bit_length_ + 1> dPhi = phi1 - phi2;
    ap_uint<eta_bit_length_> uDEta = dEta[eta_bit_length_] ? static_cast<ap_uint<eta_bit_length_>>(-dEta) : static_cast<ap_uint<eta_bit_length_>>(dEta);
    ap_uint<phi_bit_length_> uDPhi = dPhi[phi_bit_length_] ? static_cast<ap_uint<phi_bit_length_>>(-dPhi) : static_cast<ap_uint<phi_bit_length_>>(dPhi);
    if (uDPhi >= pi_digitized_in_phi_) uDPhi = 2 * pi_digitized_in_phi_ - uDPhi;
    ap_uint<phi_bit_length_ - 1>   corrDPhi = uDPhi;
    ap_uint<2*eta_bit_length_>     etaSq    = uDEta * uDEta;
    #pragma HLS bind_op variable=etaSq op=mul impl=dsp
    ap_uint<2*(phi_bit_length_-1)> phiSq    = corrDPhi * corrDPhi;
    #pragma HLS bind_op variable=phiSq op=mul impl=dsp
    return ap_uint<2*(eta_bit_length_ + phi_bit_length_)>(etaSq) + phiSq;
}

#endif
