#define USEFIFO 0
#if !USEFIFO
#include "jet_tag.h"
#define WRITE_LUT 0

// Main function
void jet_tag(input inputObjectValues[maxObjectsConsidered_], output (&outputJetValues)[nSeedsOutput_]){ // FIXME potentially use templated / overloaded func to deal with whether write out files while running synth or c-sim
    // Pragma for partitioning (allowing simultaneous access to) LUT array
    //#pragma HLS INTERFACE mode=ap_fifo port=inputObjectValues
    //#pragma HLS INTERFACE mode=ap_fifo port=outputJetValues
    //#pragma HLS ARRAY_PARTITION variable=lut_ complete //dim=1
    #pragma HLS ARRAY_PARTITION variable=inputObjectValues complete //dim=1
    #pragma HLS ARRAY_PARTITION variable=outputJetValues complete //dim=1
    //#pragma HLS ARRAY_PARTITION variable=lutR_8b_ cyclic factor=4 dim=1
    #pragma HLS bind_storage variable=inputObjectValues type=RAM_1P impl=lutram
    #pragma HLS bind_storage variable=outputJetValues type=RAM_1P impl=lutram
    //#pragma HLS bind_storage variable=lut_ type=RAM_1P impl=lutram // FIXME how to get these to use URAM rather than BRAM!? 
    //#pragma HLS ARRAY_PARTITION variable=inputObjectValues cyclic factor=4 dim=1 
    // PRAGMAS FOR WRITING DATA TO FPGA BRAMS (FOR STANDALONE HARDWARE IMPLEMENTATION TEST ONLY)
    // AXI4-Master interfaces for input arrays
    //#pragma HLS INTERFACE m_axi port=inputObjectValues        bundle=gmem0 offset=slave depth=maxObjectsConsidered_
    //#pragma HLS INTERFACE m_axi port=inputObjectValues bundle=gmem1 offset=slave depth=maxObjectsConsidered_
    // AXI4-Master interfaces for output arrays
    //#pragma HLS INTERFACE m_axi port=outputJetValues   bundle=gmem2 offset=slave depth=nSeedsOutput_
    // AXI4-Lite interface only for control signals (function arguments, etc.)
    //#pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    
    for (unsigned int i = 0; i < nSeedsOutput_; ++i) // NOTE: no effect on total latency in clock cycles
        #pragma HLS unroll
        outputJetValues[i] = 0;
    
    for (unsigned int iSeed = 0; iSeed < nSeedsOutput_; ++iSeed){ // Loop through seeds (which access up to Nth element of input object values now)
        #pragma HLS unroll
        ap_uint<eta_bit_length_ > seedEta = inputObjectValues[iSeed].range(eta_high_, eta_low_);
        ap_uint<phi_bit_length_ > seedPhi = inputObjectValues[iSeed].range(phi_high_, phi_low_);
        ap_uint<et_bit_length_ > outputJetEt = inputObjectValues[iSeed].range(et_high_, et_low_);
        //std::cout << "inputObjectValues[iSeed]: " << inputObjectValues[iSeed] << "\n";
        //std::cout << "et_high_: " << et_high_ << " , et_low_: " << et_low_ << "\n";
        //std::cout << "outputJetet: " << outputJetEt << "\n";
        //#pragma HLS ARRAY_PARTITION variable=seedEta complete //dim=0
        //#pragma HLS ARRAY_PARTITION variable=seedPhi complete //dim=0
        //#pragma HLS ARRAY_PARTITION variable=outputJetEt complete //dim=0
        //std::cout << "----------------------------- " << "\n";
        //std::cout << "iSeed : "<< iSeed << "\n";
        //std::cout << "seed eta: " << seedEta << "\n";
        //std::cout << "seed phi: " << seedPhi << "\n";
        //std::cout << "seed et: " << inputObjectValues[iSeed].range(et_high_, et_low_) << "\n";

        // EXPLICIT balanced adder tree (same rationale as jet_tag_adv.cc): the serial "outputJetEt +=
        // inputEt" reduction, even fully unrolled, is a loop-carried dependency that HLS wires as a
        // depth-N serial adder chain -- the latency bottleneck. Rewriting it as an explicit strided tree
        // makes it log2(N) deep. Unsigned addition is associative, so this is bit-identical to the serial
        // saturating sum: both compute min(seedEt + sum of passing inputEt, maxEt).
        // Stage 1 (all parallel): mask each input's Et by its deltaR^2 cut, promoted to the wide
        // accumulator. Indices [0, nSeedsOutput_) are the seeds themselves (not merge candidates) -> 0,
        // matching the original loop that started at iInput = nSeedsOutput_.
        ap_uint<et_bit_length_ + 11> etTree[maxObjectsConsidered_];
        #pragma HLS ARRAY_PARTITION variable=etTree complete
        for (unsigned int iInput = 0; iInput < maxObjectsConsidered_; ++iInput){
            #pragma HLS unroll
            if (iInput < nSeedsOutput_){
                etTree[iInput] = 0;
            } else {
                ap_uint<eta_bit_length_ > inputEta = inputObjectValues[iInput].range(eta_high_, eta_low_);
                ap_uint<phi_bit_length_ > inputPhi = inputObjectValues[iInput].range(phi_high_, phi_low_);
                ap_uint<et_bit_length_ > inputEt = inputObjectValues[iInput].range(et_high_, et_low_);
                ap_uint<2*(eta_bit_length_ + phi_bit_length_)> deltaR2 = calcDeltaR2(seedEta, seedPhi, inputEta, inputPhi);
                etTree[iInput] = (deltaR2 <= digitized_delta_R2_)
                               ? ap_uint<et_bit_length_ + 11>(inputEt)
                               : ap_uint<et_bit_length_ + 11>(0);
            }
        }
        // Stage 2 (log2(N) levels): strided pairwise reduction; the total lands in etTree[0].
        for (unsigned int stride = 1; stride < maxObjectsConsidered_; stride <<= 1){
            #pragma HLS unroll
            for (unsigned int i = 0; i + stride < maxObjectsConsidered_; i += (stride << 1)){
                #pragma HLS unroll
                etTree[i] += etTree[i + stride];
            }
        }
        // Add the seed's own Et (outputJetEt currently holds it) and clamp once to max Et.
        ap_uint<et_bit_length_ + 11> jetEtWide = ap_uint<et_bit_length_ + 11>(outputJetEt) + etTree[0];
        outputJetEt = (jetEtWide >= ((1 << (et_bit_length_)) - 1))
                    ? ap_uint<et_bit_length_>((1 << (et_bit_length_)) - 1) // clamp to max Et
                    : ap_uint<et_bit_length_>(jetEtWide);
        //std::cout << "final outputJetEt: " << outputJetEt << "\n";
        outputJetValues[iSeed].range(padded_zeroes_high_, padded_zeroes_low_) = 0; 
        outputJetValues[iSeed].range(et_high_, et_low_) = outputJetEt;
        outputJetValues[iSeed].range(eta_high_, eta_low_) = seedEta;
        outputJetValues[iSeed].range(phi_high_, phi_low_) = seedPhi;
    }
}

#endif

// FIFO EQUIVALENT VERSION - needs to be validated
#if USEFIFO
#include "jet_tag.h"
#include <hls_stream.h>

typedef ap_uint<total_bits_input_>  input_t;
typedef ap_uint<total_bits_output_> output_t;

void jet_tag(hls::stream<input_t>  &in,
             hls::stream<output_t> &out) {
#pragma HLS INTERFACE mode=ap_fifo port=in
#pragma HLS INTERFACE mode=ap_fifo port=out
#pragma HLS INTERFACE ap_ctrl_none port=return

    // Local event buffer so we can do random access exactly like the non-FIFO version
    input_t inputObjectValues[maxObjectsConsidered_];
#pragma HLS ARRAY_PARTITION variable=inputObjectValues complete
#pragma HLS bind_storage variable=inputObjectValues type=RAM_1P impl=lutram

    // Read one event worth of inputs
    READ_IN:
    for (unsigned int i = 0; i < maxObjectsConsidered_; ++i) {
#pragma HLS PIPELINE II=1
        inputObjectValues[i] = in.read();
    }

    output_t outputJetValues[nSeedsOutput_];
#pragma HLS ARRAY_PARTITION variable=outputJetValues complete
#pragma HLS bind_storage variable=outputJetValues type=RAM_1P impl=lutram

    INIT_OUT:
    for (unsigned int i = 0; i < nSeedsOutput_; ++i) {
#pragma HLS UNROLL
        outputJetValues[i] = 0;
    }

    SEEDS:
    for (unsigned int iSeed = 0; iSeed < nSeedsOutput_; ++iSeed) {
#pragma HLS UNROLL

        ap_uint<eta_bit_length_> seedEta =
            inputObjectValues[iSeed].range(eta_high_, eta_low_);
        ap_uint<phi_bit_length_> seedPhi =
            inputObjectValues[iSeed].range(phi_high_, phi_low_);
        ap_uint<et_bit_length_> outputJetEt =
            inputObjectValues[iSeed].range(et_high_, et_low_);

        INPUTS:
        for (unsigned int iInput = nSeedsOutput_; iInput < maxObjectsConsidered_; ++iInput) {
#pragma HLS UNROLL

            ap_uint<eta_bit_length_> inputEta =
                inputObjectValues[iInput].range(eta_high_, eta_low_);
            ap_uint<phi_bit_length_> inputPhi =
                inputObjectValues[iInput].range(phi_high_, phi_low_);
            ap_uint<et_bit_length_> inputEt =
                inputObjectValues[iInput].range(et_high_, et_low_);

            // Match non-FIFO version: signed differences
            ap_int<eta_bit_length_> deltaEta = seedEta - inputEta;
            ap_int<phi_bit_length_> deltaPhi = seedPhi - inputPhi;

            // Absolute values
            ap_uint<eta_bit_length_> uDeltaEta =
                deltaEta[eta_bit_length_ - 1]
                    ? static_cast<ap_uint<eta_bit_length_> >(-deltaEta)
                    : static_cast<ap_uint<eta_bit_length_> >( deltaEta);

            ap_uint<phi_bit_length_> uDeltaPhi =
                deltaPhi[phi_bit_length_ - 1]
                    ? static_cast<ap_uint<phi_bit_length_> >(-deltaPhi)
                    : static_cast<ap_uint<phi_bit_length_> >( deltaPhi);

            // Phi wrapping exactly as in non-FIFO version
            if (uDeltaPhi >= pi_digitized_in_phi_) {
                uDeltaPhi = 2 * pi_digitized_in_phi_ - uDeltaPhi;
            }

            ap_uint<phi_bit_length_ - 1> corrDeltaPhi = uDeltaPhi;

            ap_uint<phi_bit_length_ + eta_bit_length_> deltaR2 =
                uDeltaEta * uDeltaEta + corrDeltaPhi * corrDeltaPhi;
#pragma HLS bind_op variable=deltaR2 op=mul impl=dsp

            // Same deltaR^2 cut as non-FIFO version
            if (deltaR2 < digitized_delta_R2_) {
                if (outputJetEt + inputEt >= ((1 << et_bit_length_) - 1)) {
                    outputJetEt = ((1 << et_bit_length_) - 1);
                } else {
                    outputJetEt += inputEt;
                }
            }
        }

        outputJetValues[iSeed].range(padded_zeroes_high_, padded_zeroes_low_) = 0;
        outputJetValues[iSeed].range(et_high_, et_low_) = outputJetEt;
        outputJetValues[iSeed].range(eta_high_, eta_low_) = seedEta;
        outputJetValues[iSeed].range(phi_high_, phi_low_) = seedPhi;
    }

    WRITE_OUT:
    for (unsigned int i = 0; i < nSeedsOutput_; ++i) {
#pragma HLS PIPELINE II=1
        out.write(outputJetValues[i]);
    }
}
#endif