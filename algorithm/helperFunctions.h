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

// Jet eta code (10b two's complement, 0.0125 LSB) -> tower-grid index floor(eta / 0.1) + eta_range_ / 2, clamped to the grid
inline ap_uint<eta_code_bits_> jetEtaToTowerIndex(ap_uint<jet_eta_bit_length_> jetEta) {
    #pragma HLS INLINE
    ap_int<jet_eta_bit_length_> signedEta;
    signedEta.range(jet_eta_bit_length_ - 1, 0) = jetEta; // same bits, read as two's complement
    ap_int<jet_eta_bit_length_ - jet_to_tower_eta_shift_ + 1> index =
        (signedEta >> jet_to_tower_eta_shift_) + static_cast<int>(eta_range_ / 2);
    if (index < 0) index = 0;
    if (index > static_cast<int>(eta_range_) - 1) index = static_cast<int>(eta_range_) - 1;
    return ap_uint<eta_code_bits_>(index);
}

// Tower-grid index -> jet eta code of that tower's centre, 8 * (index - eta_range_ / 2) + 4
inline ap_uint<jet_eta_bit_length_> towerIndexToJetEta(ap_uint<eta_code_bits_> index) {
    #pragma HLS INLINE
    ap_int<jet_eta_bit_length_> signedEta =
        (ap_int<jet_eta_bit_length_ + 1>(index) - static_cast<int>(eta_range_ / 2)) * (1 << jet_to_tower_eta_shift_)
        + (1 << (jet_to_tower_eta_shift_ - 1));
    return signedEta.range(jet_eta_bit_length_ - 1, 0);
}

// Jet phi code (9b, LSB 2pi/512 from phi = 0) -> tower-grid index
inline ap_uint<phi_code_bits_> jetPhiToTowerIndex(ap_uint<jet_phi_bit_length_> jetPhi) {
    #pragma HLS INLINE
    return jetPhi.range(jet_phi_bit_length_ - 1, jet_to_tower_phi_shift_);
}

// Tower-grid index -> jet phi code, the index with jet_to_tower_phi_shift_ zero LSBs
inline ap_uint<jet_phi_bit_length_> towerIndexToJetPhi(ap_uint<phi_code_bits_> index) {
    #pragma HLS INLINE
    return ap_uint<jet_phi_bit_length_>(index) << jet_to_tower_phi_shift_;
}

// Tower-grid indices of an input jet's eta / phi field, per input_jets_use_jet_tob_encoding_
inline ap_uint<eta_code_bits_> inputJetEtaIndex(ap_uint<eta_bit_length_> field) {
    #pragma HLS INLINE
    return input_jets_use_jet_tob_encoding_ ? jetEtaToTowerIndex(field) : ap_uint<eta_code_bits_>(field.range(eta_code_bits_ - 1, 0));
}
inline ap_uint<phi_code_bits_> inputJetPhiIndex(ap_uint<phi_bit_length_> field) {
    #pragma HLS INLINE
    return input_jets_use_jet_tob_encoding_ ? jetPhiToTowerIndex(field) : ap_uint<phi_code_bits_>(field.range(phi_code_bits_ - 1, 0));
}

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
