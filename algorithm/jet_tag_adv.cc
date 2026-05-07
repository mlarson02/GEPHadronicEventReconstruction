#include "jet_tag_adv.h"
#define WRITE_LUT 0

// Main function
void jet_tag_adv(input seedValues[nTotalSeeds_], input inputObjectValues[maxObjectsConsidered_], output (&outputJetValues)[nSeedsOutput_]){ // FIXME potentially use templated / overloaded func to deal with whether write out files while running synth or c-sim
    // Pragma for partitioning (allowing simultaneous access to) LUT array
    #pragma HLS ARRAY_PARTITION variable=seedValues complete
    #pragma HLS ARRAY_PARTITION variable=inputObjectValues complete
    #pragma HLS ARRAY_PARTITION variable=outputJetValues complete
    // PRAGMAS FOR WRITING DATA TO FPGA BRAMS (TESTING IMPLEMENTATION ONLY)
    // AXI4-Master interfaces for input arrays
    //#pragma HLS INTERFACE m_axi port=seedValues        bundle=gmem0 offset=slave depth=nSeedsInput_
    //#pragma HLS INTERFACE m_axi port=inputObjectValues bundle=gmem1 offset=slave depth=maxObjectsConsidered_
    // AXI4-Master interfaces for output arrays
    //#pragma HLS INTERFACE m_axi port=outputJetValues   bundle=gmem2 offset=slave depth=nSeedsOutput_
    // AXI4-Lite interface only for control signals (function arguments, etc.)
    //#pragma HLS INTERFACE s_axilite port=return bundle=CTRL
    /*std::cout << "--------------------" << "\n";
    std::cout << "constants: "<< "\n";
    std::cout << "deltaR2_granularity_: " << deltaR2_granularity_ << "\n";
    std::cout << "digitized_delta_R2Cut_: " << digitized_delta_R2Cut_ << "\n";
    std::cout << "r2 cut: " << r2Cut_ << "\n";*/
    for (unsigned int i = 0; i < nSeedsOutput_; ++i)
        #pragma HLS unroll
        outputJetValues[i] = 0;

    // Save original seed positions before any modifications; used for subjet finding later
    input seedValuesForSubjets[nTotalSeeds_];
    #pragma HLS ARRAY_PARTITION variable=seedValuesForSubjets complete
    for (unsigned int i = 0; i < nTotalSeeds_; ++i){
        #pragma HLS unroll
        seedValuesForSubjets[i] = seedValues[i];
    }

    // --- Overlap Removal (OR) ---
    // Ensure leading and subleading seeds don't overlap within deltaR < 2 * jet radius.
    // If they do, search seeds 3–N for the highest-Et seed that is separated from the leading
    // seed by > 2R and swap it in as the subleading seed.
    ap_uint<deltaR2_bits_> deltaR2LeadingSubleading = calcDeltaR2(
        seedValues[0].range(eta_high_, eta_low_), seedValues[0].range(phi_high_, phi_low_),
        seedValues[1].range(eta_high_, eta_low_), seedValues[1].range(phi_high_, phi_low_)
    );

    //std::cout << "seed 1 OR eta: " << seedValues[0].range(eta_high_, eta_low_) << " , phi: " << seedValues[0].range(phi_high_, phi_low_) << "\n";
    //std::cout << "seed 2 OR eta: " << seedValues[1].range(eta_high_, eta_low_) << " , phi: " << seedValues[1].range(phi_high_, phi_low_) << "\n";
    //std::cout << "deltaR2leadingsubleading: " << deltaR2LeadingSubleading << "\n";
    //std::cout << "deltaR^2 cut: " << 2 * 2 * digitized_delta_R2Cut_ << "\n";
    if(deltaR2LeadingSubleading <= 2 * 2 * digitized_delta_R2Cut_){ // (2 * R_Cut) ^ 2

        for(unsigned int iSeedOR = nSeedsOutput_; iSeedOR < nTotalSeeds_; iSeedOR++){
            #pragma HLS unroll
            //std::cout << "iSeed OR : "<< iSeedOR << "\n";
            if(seedValues[iSeedOR].range(et_high_, et_low_) == 0 && seedValues[iSeedOR].range(eta_high_, eta_low_) == 0 && seedValues[iSeedOR].range(phi_high_, phi_low_) == 0) continue; // don't consider if et, eta, phi all = 0 (no jet)
            ap_uint<deltaR2_bits_> deltaR2LeadingSubleadingNthLeading = calcDeltaR2(
                seedValues[0].range(eta_high_, eta_low_), seedValues[0].range(phi_high_, phi_low_),
                seedValues[iSeedOR].range(eta_high_, eta_low_), seedValues[iSeedOR].range(phi_high_, phi_low_)
            );
            //std::cout << "deltaR2 Nth leading: " << deltaR2LeadingSubleadingNthLeading << "\n";
            if(deltaR2LeadingSubleadingNthLeading > 2 * 2 * digitized_delta_R2Cut_){
                //std::cout << "triggering OR for iSeed OR : "<< iSeedOR << "\n";
                ap_uint<et_bit_length_ > swappedEt = seedValues[iSeedOR].range(et_high_, et_low_);
                ap_uint<eta_bit_length_ > swappedEta = seedValues[iSeedOR].range(eta_high_, eta_low_);
                ap_uint<phi_bit_length_ > swappedPhi = seedValues[iSeedOR].range(phi_high_, phi_low_);

                ap_uint<et_bit_length_ > originalSubleadingEt = seedValues[1].range(et_high_, et_low_);
                ap_uint<eta_bit_length_ > originalSubleadingEta = seedValues[1].range(eta_high_, eta_low_);
                ap_uint<phi_bit_length_ > originalSubleadingPhi = seedValues[1].range(phi_high_, phi_low_);

                // Swap the entire (Et, eta, phi) triplet for original subleading, new subleading seed
                seedValues[1].range(et_high_, et_low_) = swappedEt;
                seedValues[1].range(eta_high_, eta_low_) = swappedEta;
                seedValues[1].range(phi_high_, phi_low_) = swappedPhi;

                seedValues[iSeedOR].range(et_high_, et_low_) = originalSubleadingEt;
                seedValues[iSeedOR].range(eta_high_, eta_low_) = originalSubleadingEta;
                seedValues[iSeedOR].range(phi_high_, phi_low_) = originalSubleadingPhi;

                break; // break out of loop as we've found something to swap, thus preventing overlapping large-R jets
            } // If deltaR^2 between original leading, other proto-seed farther than 2 times jet radius
        } // Loop through potential additional seeds to find
    } // If deltaR^2 between original leading, subleading closer than 2 times jet radius


    // FIXME make this entire process more dynamic to account for nSeedsOutput_ != 2 (progressively do this for highest Et seeds rather than for 1st 2 seeds immediately)
    // --- Seed Position Optimization ---
    // For each output seed, check seeds 3–6 (nSeedsDeltaR_ candidates) to find any within
    // a search radius. If found, shift the output seed position to the midpoint between itself
    // and the closest/highest-Et candidate.
    ap_uint<4> protoSeedCounter[nSeedsOutput_] = {0, 0};

    bool indicesofProtoSeeds[nSeedsOutput_][nSeedsDeltaR_] = {false};
    ap_int<4> indices[2] = {-1, -1};
    #pragma HLS bind_storage variable=indicesofProtoSeeds type=RAM_1P impl=lutram
    //std::cout << "Seed 1 Eta, Phi: " << seedValues[0].range(eta_high_, eta_low_) << " , " << seedValues[0].range(phi_high_, phi_low_) << "\n";
    //std::cout << "Seed 2 Eta, Phi: " << seedValues[1].range(eta_high_, eta_low_) << " , " << seedValues[1].range(phi_high_, phi_low_) << "\n";
    //std::cout << " ----------------- SEED POS OPT -----------------" << "\n";
    for (unsigned int iSeed = 0; iSeed < nSeedsOutput_; iSeed++){
        #pragma HLS unroll

        //std::cout << "iSeed: " << iSeed << "\n";
        for (unsigned int iPreSeed = 0; iPreSeed < nSeedsDeltaR_; iPreSeed++){
            #pragma HLS unroll
            //std::cout << "iPreSeed: " << iPreSeed << "\n";
            //std::cout << "seed eta: " << seedValues[iSeed].range(eta_high_, eta_low_) << " , phi: " << seedValues[iSeed].range(phi_high_, phi_low_) << " , et: " << seedValues[iSeed].range(et_high_, et_low_) << "\n";
            //std::cout << "pre seed eta: " << seedValues[iPreSeed + nSeedsOutput_].range(eta_high_, eta_low_) << " , phi: " << seedValues[iPreSeed + nSeedsOutput_].range(phi_high_, phi_low_) << " , et: " << seedValues[iPreSeed + nSeedsOutput_].range(et_high_, et_low_) << "\n";
            if(seedValues[iPreSeed + nSeedsOutput_].range(et_high_, et_low_) == 0) continue;
            ap_uint<deltaR2_bits_> deltaR2 = calcDeltaR2(
                seedValues[iSeed].range(eta_high_, eta_low_), seedValues[iSeed].range(phi_high_, phi_low_),
                seedValues[iPreSeed + nSeedsOutput_].range(eta_high_, eta_low_), seedValues[iPreSeed + nSeedsOutput_].range(phi_high_, phi_low_)
            );
            //std::cout << "deltaR2: " << deltaR2 << " , digitized_d_search_squared_: " << digitized_d_search_squared_ << "\n";
            if (deltaR2 <= digitized_d_search_squared_){
                protoSeedCounter[iSeed]++;
                indicesofProtoSeeds[iSeed][iPreSeed] = true;
                indices[iSeed] = iPreSeed;
                //std::cout << "setting indices to: " << iPreSeed << " , -> indices[iSeed]: " << indices[iSeed] << "\n";
            }
        }
        //std::cout << "proto seed coutner[iseed]: " <<  protoSeedCounter[iSeed] << "\n";
    }

    // When multiple candidates are within the search radius, use the highest-Et one as the target midpoint
    // For seed 0
    if (protoSeedCounter[0] > 1) {
        //std::cout << "mnultiple proto seeds [seed 0]:" << "\n";
        ap_uint<et_bit_length_> maxEt = 0;
        for (unsigned int iPreSeed = 0; iPreSeed < nSeedsDeltaR_; iPreSeed++) {
            #pragma HLS unroll
            if (indicesofProtoSeeds[0][iPreSeed]) {
                ap_uint<et_bit_length_> et = seedValues[iPreSeed + nSeedsOutput_].range(et_high_, et_low_);
                if (et > maxEt) {
                    maxEt = et;
                    indices[0] = iPreSeed;
                }
            }
        }
    }

    // For seed 1
    if (protoSeedCounter[1] > 1) {
        //std::cout << "mnultiple proto seeds [seed 1]:" << "\n";
        ap_uint<et_bit_length_> maxEt = 0;
        for (unsigned int iPreSeed = 0; iPreSeed < nSeedsDeltaR_; iPreSeed++) {
            #pragma HLS unroll
            if (indicesofProtoSeeds[1][iPreSeed]) {
                ap_uint<et_bit_length_> et = seedValues[iPreSeed + nSeedsOutput_].range(et_high_, et_low_);
                if (et > maxEt) {
                    maxEt = et;
                    indices[1] = iPreSeed;
                }
            }
        }
    }

    // If both seeds would shift toward the same candidate, prevent the lower-Et seed from shifting
    // to avoid unnecessary overlap. FIXME implement this cleanly in HLS
    bool skipSecondSeed = false;
    if(indices[0] == indices[1] && indices[0] != -1){
        skipSecondSeed = true;
    }
    // Shift seed positions to the midpoint between the seed and its closest candidate
    for (unsigned int iSeed = 0; iSeed < nSeedsOutput_; iSeed++){
        #pragma HLS unroll
        //std::cout << "iSeed: " << iSeed << " , skipSecondSeed: " << skipSecondSeed << "\n";
        if(skipSecondSeed == true && iSeed == 1) continue;
        if(indices[iSeed] == -1) continue;

        //std::cout << "-------------- calcing mid point -----------------" << "\n";
        //std::cout << "iSeed: " << iSeed << "\n";
        //fflush(stdout);

        const ap_int<phi_bit_length_ + 2> PI_D     =  ap_int<phi_bit_length_ + 2>(pi_digitized_in_phi_);
        const ap_int<phi_bit_length_ + 2> TWO_PI_D =  ap_int<phi_bit_length_ + 2>((1 << phi_bit_length_) - 1);

        //std::cout << "pi_digitized_in_phi_: " << pi_digitized_in_phi_ << "\n";
        //std::cout << "PI_D: " << PI_D  << " , TWO_PI_D: " << TWO_PI_D << "\n";

        // Wraps a signed phi value into [-PI, PI) using digitized constants
        auto wrapSym = [&](ap_int<phi_bit_length_ + 2> x) -> ap_int<phi_bit_length_ + 2> {
            #pragma HLS INLINE
            if (x >  PI_D)   x -= TWO_PI_D;
            if (x < -PI_D)  x += TWO_PI_D;
            return x;
        };

        // Convert eta/phi from unsigned digitized format to signed centered at 0
        ap_int<eta_bit_length_ + 1> eta1 = seedValues[iSeed].range(eta_high_, eta_low_) - (1 << (eta_bit_length_ - 1));
        ap_int<eta_bit_length_ + 1> eta2 = seedValues[indices[iSeed] + nSeedsOutput_].range(eta_high_, eta_low_) - (1 << (eta_bit_length_ - 1));
        //std::cout << "eta 1 : " << eta1 << " and eta2 : " << eta2 << "\n";

        //std::cout << "phi1: " << seedValues[iSeed].range(phi_high_, phi_low_) << " and phi2: " << seedValues[indices[iSeed] + nSeedsOutput_].range(phi_high_, phi_low_) << "\n";
        ap_int<phi_bit_length_ + 2> phi1s = seedValues[iSeed].range(phi_high_, phi_low_) - (1 << (phi_bit_length_ - 1));
        ap_int<phi_bit_length_ + 2> phi2s = seedValues[indices[iSeed] + nSeedsOutput_].range(phi_high_, phi_low_) - (1 << (phi_bit_length_ - 1));
        //std::cout << "phi1s: " << phi1s << " and phi2s: " << phi2s << "\n";

        // --- Shortest-arc phi midpoint ---
        ap_int<phi_bit_length_ + 2> dphi = phi2s - phi1s;
        //std::cout << "dphi before wrap :" << dphi << "\n";
        dphi = wrapSym(dphi);                          // now in [-PI, PI)
        //std::cout << "dphi after wrap :" << dphi << "\n";
        ap_int<phi_bit_length_ + 2> phi_mid = phi1s + (dphi >> 1);  // arithmetic shift (divide by 2)
        //std::cout << "phi_mid before wrap :" << phi_mid << "\n";
        phi_mid = wrapSym(phi_mid);
        //std::cout << "phi_mid after wrap :" << phi_mid << "\n";

        // --- Unweighted eta midpoint ---
        ap_int<eta_bit_length_ + 1> eta_mid = (eta1 + eta2) >> 1;

        //std::cout << "eta mid : " << eta_mid << " phi_mid: " << phi_mid << "\n";

        // Convert midpoints back to digitized unsigned format
        ap_uint<eta_bit_length_> eta_mid_digitized = eta_mid + (1 << (eta_bit_length_ - 1));
        ap_uint<phi_bit_length_> phi_mid_digitized = phi_mid + (1 << (phi_bit_length_ - 1));

        //std::cout << "eta_mid_digitized: " << eta_mid_digitized << " phi_mid_digitized: " << phi_mid_digitized << "\n";
        seedValues[iSeed].range(eta_high_, eta_low_) = eta_mid_digitized;
        seedValues[iSeed].range(phi_high_, phi_low_) = phi_mid_digitized;
    }

    // must also account for when deltaR2ValuesSeed not set (no deltaR < 2.5 --> leave seed position as is)

    // --- Jet building: accumulate Et from input objects within deltaR < R_cut ---
    // Note: seedValues eta/phi may differ from seedValuesForSubjets after OR and seed optimization.
    // Subjet counting uses original seed positions (seedValuesForSubjets) for subjet identification.
    for (unsigned int iSeed = 0; iSeed < nSeedsOutput_; ++iSeed){
        #pragma HLS unroll
        ap_uint<et_bit_length_ > outputJetEt = 0;
        ap_uint<num_subjets_length_ > numSubjets = 0;

        //std::cout << "----------" << "\n";
        //std::cout << "iSeed: " << iSeed << "\n";
        //std::cout << "ORIGINAL seed eta: " << seedValuesForSubjets[iSeed].range(eta_high_, eta_low_) << " , phi: " << seedValuesForSubjets[iSeed].range(phi_high_, phi_low_) << " , et: " << seedValuesForSubjets[iSeed].range(et_high_, et_low_) << "\n";
        //std::cout << "seed eta: " << seedValues[iSeed].range(eta_high_, eta_low_) << " , phi: " << seedValues[iSeed].range(phi_high_, phi_low_) << " , et: " << seedValues[iSeed].range(et_high_, et_low_) << "\n";

        if(seedValues[iSeed].range(et_high_, et_low_) != 0){
            for (unsigned int iInput = 0; iInput < maxObjectsConsidered_; ++iInput){
                #pragma HLS unroll
                //std::cout << "iInput: " << iInput << "\n";
                //std::cout << "input eta: " << inputObjectValues[iInput].range(eta_high_, eta_low_) << " , phi: " << inputObjectValues[iInput].range(phi_high_, phi_low_) << " , et: " << inputObjectValues[iInput].range(et_high_, et_low_) << "\n";
                ap_uint<deltaR2_bits_> deltaR2 = calcDeltaR2(
                    seedValues[iSeed].range(eta_high_, eta_low_), seedValues[iSeed].range(phi_high_, phi_low_),
                    inputObjectValues[iInput].range(eta_high_, eta_low_), inputObjectValues[iInput].range(phi_high_, phi_low_)
                );
                //std::cout << "deltaR2 : " << deltaR2 << " , digitized deltaR2 cut: " << digitized_delta_R2Cut_ << "\n";
                if (deltaR2 <= digitized_delta_R2Cut_){ // FIXME will throw an error for deltaR!=1.0
                    if(outputJetEt + inputObjectValues[iInput][0] >= ((1 << et_bit_length_) - 1)){
                        outputJetEt = ((1 << et_bit_length_) - 1); // clamp to max Et
                    }
                    else{
                        outputJetEt += inputObjectValues[iInput].range(et_high_, et_low_);
                    }
                    //std::cout << "outputJetEt after sum: " << outputJetEt << "\n";
                }
            }

            // Count subjets: original seeds within deltaR < R_cut of the (possibly shifted) output seed center
            for(unsigned int iSubjet = 0; iSubjet < nTotalSeeds_; iSubjet++){
                #pragma HLS unroll
                //std::cout << "iSubjet: " << iSubjet << "\n";
                //std::cout << "seedValuesForSubjets[iSubjet].range(et_high_, et_low_): " << seedValuesForSubjets[iSubjet].range(et_high_, et_low_) << " , subjet_et_threshold_: " << subjet_et_threshold_ << "\n";
                if(seedValuesForSubjets[iSubjet].range(et_high_, et_low_) > subjet_et_threshold_){ // > 25 GeV nominally
                    //std::cout << " passes et threshold" << "\n";
                    ap_uint<deltaR2_bits_> deltaR2 = calcDeltaR2(
                        seedValues[iSeed].range(eta_high_, eta_low_), seedValues[iSeed].range(phi_high_, phi_low_),
                        seedValuesForSubjets[iSubjet].range(eta_high_, eta_low_), seedValuesForSubjets[iSubjet].range(phi_high_, phi_low_)
                    );
                    //std::cout << "deltaR2: " << deltaR2 << " , delta_R2Cut_: " << digitized_delta_R2Cut_ << "\n";
                    if (deltaR2 <= digitized_delta_R2Cut_){
                        //std::cout << "FOUND SUBJET" << "\n";
                        numSubjets = (numSubjets < ap_uint<num_subjets_length_>((1 << num_subjets_length_) - 1))
                                    ? ap_uint<num_subjets_length_>(numSubjets + 1)
                                    : ap_uint<num_subjets_length_>((1 << num_subjets_length_) - 1);
                    }
                }
            }
        }

        //std::cout << "numSubjets: " << numSubjets << "\n";
        outputJetValues[iSeed].range(padded_zeroes_high_, padded_zeroes_low_) = 0;
        outputJetValues[iSeed].range(num_subjets_high_, num_subjets_low_) = numSubjets;
        outputJetValues[iSeed].range(et_high_, et_low_) = outputJetEt;
        outputJetValues[iSeed].range(eta_high_, eta_low_) = seedValues[iSeed].range(eta_high_, eta_low_);
        outputJetValues[iSeed].range(phi_high_, phi_low_) = seedValues[iSeed].range(phi_high_, phi_low_);
    }
}
