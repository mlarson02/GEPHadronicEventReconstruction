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
// Applies shortest-arc phi wrapping. Return type is wide enough to never overflow:
// max is (eta_range_-1)^2 + pi_digitized_in_phi_^2 = 97^2 + 32^2 = 10433, which needs 14 bits.
inline ap_uint<2*eta_code_bits_ + 1> calcDeltaR2(
    ap_uint<eta_code_bits_> eta1, ap_uint<phi_code_bits_> phi1,
    ap_uint<eta_code_bits_> eta2, ap_uint<phi_code_bits_> phi2
) {
    #pragma HLS INLINE
    ap_int<eta_code_bits_ + 1> dEta = eta1 - eta2;
    ap_int<phi_code_bits_ + 1> dPhi = phi1 - phi2;
    ap_uint<eta_code_bits_> uDEta = dEta[eta_code_bits_] ? static_cast<ap_uint<eta_code_bits_>>(-dEta) : static_cast<ap_uint<eta_code_bits_>>(dEta);
    ap_uint<phi_code_bits_> uDPhi = dPhi[phi_code_bits_] ? static_cast<ap_uint<phi_code_bits_>>(-dPhi) : static_cast<ap_uint<phi_code_bits_>>(dPhi);
    if (uDPhi >= pi_digitized_in_phi_) uDPhi = 2 * pi_digitized_in_phi_ - uDPhi;
    // phi_code_bits_ wide, NOT phi_code_bits_ - 1: the wrapped value reaches pi_digitized_in_phi_
    // (32) when two objects are exactly pi apart, which needs 6 bits. Narrowing this by one is the
    // bug that was in jet_tag_adv.cc -- 32 truncated to 0, so back-to-back objects read as coincident.
    ap_uint<phi_code_bits_>   corrDPhi = uDPhi;
    ap_uint<2*eta_code_bits_> etaSq    = uDEta * uDEta;
    #pragma HLS bind_op variable=etaSq op=mul impl=dsp
    ap_uint<2*phi_code_bits_> phiSq    = corrDPhi * corrDPhi;
    #pragma HLS bind_op variable=phiSq op=mul impl=dsp
    return ap_uint<2*eta_code_bits_ + 1>(etaSq) + phiSq;
}

#endif
