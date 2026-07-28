#include "jet_tag_adv.h"
#define WRITE_LUT 0

// Main function
void jet_tag_adv(input seedValues[nTotalSeeds_], input inputObjectValues[maxObjectsConsidered_], output (&outputJetValues)[nSeedsOutput_]){ // FIXME potentially use templated / overloaded func to deal with whether write out files while running synth or c-sim
    // Pragma for partitioning (allowing simultaneous access to) LUT array
    #pragma HLS ARRAY_PARTITION variable=seedValues complete
    #pragma HLS ARRAY_PARTITION variable=inputObjectValues complete
    #pragma HLS ARRAY_PARTITION variable=outputJetValues complete
#if !USE_DSPS_
    // Fully partition the deltaR LUTs so every unrolled call site gets its own concurrent read,
    // decomposing them into individual registers rather than a single memory block (so there's no
    // BRAM to bind away in the first place).
    #pragma HLS ARRAY_PARTITION variable=lut_ complete dim=1
    #pragma HLS ARRAY_PARTITION variable=lutR_ complete dim=1
#endif
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
    //std::cout << "seed 1 OR eta: " << seedValues[0].range(eta_high_, eta_low_) << " , phi: " << seedValues[0].range(phi_high_, phi_low_) << "\n";
    //std::cout << "seed 2 OR eta: " << seedValues[1].range(eta_high_, eta_low_) << " , phi: " << seedValues[1].range(phi_high_, phi_low_) << "\n";
    // Restructured as compute-always-then-mask for deterministic (data-independent) latency: every
    // candidate check runs unconditionally, a priority encoder selects the first qualifying seed, and
    // the swap is applied via a masked write. Bit-identical to the original guarded-loop-with-break:
    // same "first non-empty seed separated from the leading seed by > 2R" selection, same single swap.
    bool orOverlap = passesTwoRCut(
        seedValues[0].range(eta_high_, eta_low_), seedValues[0].range(phi_high_, phi_low_),
        seedValues[1].range(eta_high_, eta_low_), seedValues[1].range(phi_high_, phi_low_)
    ); // leading & subleading overlap within 2 * rCut_

    // For each candidate seed (3..N): a valid swap target iff it is non-empty AND separated from the
    // leading seed by > 2R.
    bool orCand[nTotalSeeds_] = {false};
    #pragma HLS ARRAY_PARTITION variable=orCand complete
    for(unsigned int iSeedOR = nSeedsOutput_; iSeedOR < nTotalSeeds_; iSeedOR++){
        #pragma HLS unroll
        bool nonEmpty = !(seedValues[iSeedOR].range(et_high_, et_low_) == 0 && seedValues[iSeedOR].range(eta_high_, eta_low_) == 0 && seedValues[iSeedOR].range(phi_high_, phi_low_) == 0);
        orCand[iSeedOR] = nonEmpty && !passesTwoRCut(
            seedValues[0].range(eta_high_, eta_low_), seedValues[0].range(phi_high_, phi_low_),
            seedValues[iSeedOR].range(eta_high_, eta_low_), seedValues[iSeedOR].range(phi_high_, phi_low_)
        );
    }

    // Priority-select the lowest-index qualifying candidate (matches the original's first-match + break).
    bool orFound = false;
    ap_uint<4> orSwapIdx = 0;
    for(unsigned int iSeedOR = nSeedsOutput_; iSeedOR < nTotalSeeds_; iSeedOR++){
        #pragma HLS unroll
        if(orCand[iSeedOR] && !orFound){
            orFound = true;
            orSwapIdx = iSeedOR;
        }
    }
    bool doOrSwap = orOverlap && orFound;

    // Read both triplets first, then apply the swap under a single mask (a dynamic-index read/write on
    // the fully-partitioned seedValues array is a mux/demux -- fixed latency, not a control branch).
    ap_uint<et_bit_length_ > candEt  = seedValues[orSwapIdx].range(et_high_, et_low_);
    ap_uint<eta_bit_length_> candEta = seedValues[orSwapIdx].range(eta_high_, eta_low_);
    ap_uint<phi_bit_length_> candPhi = seedValues[orSwapIdx].range(phi_high_, phi_low_);
    ap_uint<et_bit_length_ > subEt  = seedValues[1].range(et_high_, et_low_);
    ap_uint<eta_bit_length_> subEta = seedValues[1].range(eta_high_, eta_low_);
    ap_uint<phi_bit_length_> subPhi = seedValues[1].range(phi_high_, phi_low_);
    if(doOrSwap){
        // Swap the (Et, eta, phi) triplet between the original subleading seed (slot 1) and the candidate
        seedValues[1].range(et_high_, et_low_)   = candEt;
        seedValues[1].range(eta_high_, eta_low_) = candEta;
        seedValues[1].range(phi_high_, phi_low_) = candPhi;
        seedValues[orSwapIdx].range(et_high_, et_low_)   = subEt;
        seedValues[orSwapIdx].range(eta_high_, eta_low_) = subEta;
        seedValues[orSwapIdx].range(phi_high_, phi_low_) = subPhi;
    }


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
            #pragma HLS unroll // fully unrolled for deterministic (feed-forward) latency
            //std::cout << "iPreSeed: " << iPreSeed << "\n";
            //std::cout << "seed eta: " << seedValues[iSeed].range(eta_high_, eta_low_) << " , phi: " << seedValues[iSeed].range(phi_high_, phi_low_) << " , et: " << seedValues[iSeed].range(et_high_, et_low_) << "\n";
            //std::cout << "pre seed eta: " << seedValues[iPreSeed + nSeedsOutput_].range(eta_high_, eta_low_) << " , phi: " << seedValues[iPreSeed + nSeedsOutput_].range(phi_high_, phi_low_) << " , et: " << seedValues[iPreSeed + nSeedsOutput_].range(et_high_, et_low_) << "\n";
            if(seedValues[iPreSeed + nSeedsOutput_].range(et_high_, et_low_) == 0) continue;
            if (passesSearchRadius(
                seedValues[iSeed].range(eta_high_, eta_low_), seedValues[iSeed].range(phi_high_, phi_low_),
                seedValues[iPreSeed + nSeedsOutput_].range(eta_high_, eta_low_), seedValues[iPreSeed + nSeedsOutput_].range(phi_high_, phi_low_)
            )){
                protoSeedCounter[iSeed]++;
                indicesofProtoSeeds[iSeed][iPreSeed] = true;
                indices[iSeed] = iPreSeed;
                //std::cout << "setting indices to: " << iPreSeed << " , -> indices[iSeed]: " << indices[iSeed] << "\n";
            }
        }
        //std::cout << "proto seed coutner[iseed]: " <<  protoSeedCounter[iSeed] << "\n";
    }

    // When multiple candidates are within the search radius, use the highest-Et one as the target
    // midpoint. Compute-always-then-mask for deterministic latency: always find the highest-Et flagged
    // candidate into a local, then override indices[iSeed] only when there were multiple candidates
    // (protoSeedCounter > 1). maxIdx defaults to the current indices[iSeed] so that if no flagged
    // candidate has Et > 0 the value is unchanged -- bit-identical to the original guarded loops.
    // For seed 0
    {
        ap_uint<et_bit_length_> maxEt = 0;
        ap_int<4> maxIdx = indices[0];
        for (unsigned int iPreSeed = 0; iPreSeed < nSeedsDeltaR_; iPreSeed++) {
            #pragma HLS unroll
            if (indicesofProtoSeeds[0][iPreSeed]) {
                ap_uint<et_bit_length_> et = seedValues[iPreSeed + nSeedsOutput_].range(et_high_, et_low_);
                if (et > maxEt) {
                    maxEt = et;
                    maxIdx = iPreSeed;
                }
            }
        }
        if (protoSeedCounter[0] > 1) indices[0] = maxIdx;
    }

    // For seed 1
    {
        ap_uint<et_bit_length_> maxEt = 0;
        ap_int<4> maxIdx = indices[1];
        for (unsigned int iPreSeed = 0; iPreSeed < nSeedsDeltaR_; iPreSeed++) {
            #pragma HLS unroll
            if (indicesofProtoSeeds[1][iPreSeed]) {
                ap_uint<et_bit_length_> et = seedValues[iPreSeed + nSeedsOutput_].range(et_high_, et_low_);
                if (et > maxEt) {
                    maxEt = et;
                    maxIdx = iPreSeed;
                }
            }
        }
        if (protoSeedCounter[1] > 1) indices[1] = maxIdx;
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
        // compute-always-then-mask: the midpoint is always computed (fixed latency) and written back
        // only when this seed actually has a candidate to shift toward. Bit-identical to the two
        // continues below (indices[iSeed] + nSeedsOutput_ is a valid index even when indices==-1 ->
        // reads slot 1, but the result is discarded since doShift is false).
        bool doShift = !(skipSecondSeed == true && iSeed == 1) && (indices[iSeed] != -1);

        // Safe candidate index for the always-on read below: indices[iSeed] is -1 (ap_int) when no
        // candidate was found, and indices[iSeed] + nSeedsOutput_ could then mix signed/unsigned into an
        // out-of-range index. Substitute a dummy in-range slot in that case -- the midpoint is discarded
        // via doShift, so which valid slot we read is irrelevant.
        ap_uint<4> candSeedIdx = (indices[iSeed] == -1)
                               ? ap_uint<4>(nSeedsOutput_)
                               : ap_uint<4>((unsigned int)indices[iSeed] + nSeedsOutput_);

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
        ap_int<eta_bit_length_ + 1> eta2 = seedValues[candSeedIdx].range(eta_high_, eta_low_) - (1 << (eta_bit_length_ - 1));
        //std::cout << "eta 1 : " << eta1 << " and eta2 : " << eta2 << "\n";

        //std::cout << "phi1: " << seedValues[iSeed].range(phi_high_, phi_low_) << " and phi2: " << seedValues[indices[iSeed] + nSeedsOutput_].range(phi_high_, phi_low_) << "\n";
        ap_int<phi_bit_length_ + 2> phi1s = seedValues[iSeed].range(phi_high_, phi_low_) - (1 << (phi_bit_length_ - 1));
        ap_int<phi_bit_length_ + 2> phi2s = seedValues[candSeedIdx].range(phi_high_, phi_low_) - (1 << (phi_bit_length_ - 1));
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
        if(doShift){
            seedValues[iSeed].range(eta_high_, eta_low_) = eta_mid_digitized;
            seedValues[iSeed].range(phi_high_, phi_low_) = phi_mid_digitized;
        }
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

        // compute-always-then-mask for deterministic latency: the Et accumulation and subjet count run
        // unconditionally, then both outputs are zeroed for an empty seed. Bit-identical to the previous
        // if(seedEt != 0) guard, which produced 0/0 for an empty seed by skipping the loops (the vars
        // init to 0). Running the loops for an empty seed only affects these locals, which are masked.
        bool seedActive = (seedValues[iSeed].range(et_high_, et_low_) != 0);

        // Accumulate Et from all passing input objects. Fully parallel: every passesRCut is evaluated
        // concurrently, the masked Et values are summed via a balanced adder tree (HLS infers the tree
        // from the unrolled reduction, no mid-loop clamp), and the result is clamped once at the end.
        // This also fixes the prior clamp, which tested inputObjectValues[iInput][0] (bit 0 only) instead
        // of the full Et; it is now a proper saturating sum. outputJetEtWide holds the sum of
        // maxObjectsConsidered_ full-width Et values without overflow before the clamp (et_bit_length_ +
        // 11 covers up to 2048 input objects).
        // EXPLICIT balanced adder tree. The plain "outputJetEtWide += et" reduction (even fully unrolled)
        // is a loop-carried dependency, so HLS wired the 128 adders in a serial chain (~depth 128), which
        // serialized the whole otherwise-parallel input loop and dominated latency (the schedule showed
        // ~56 of ~72 phase-3 cycles here). Writing the reduction as an explicit strided tree forces
        // log2(N) depth (~7 levels) instead. Unsigned addition is associative, so this is bit-identical
        // to the serial sum. FIXME (deltaR!=1.0) carried over from the passesRCut mask below.
        // Stage 1 (all parallel): mask each input's Et by its passesRCut, promoted to the wide accumulator.
        ap_uint<et_bit_length_ + 11> etTree[maxObjectsConsidered_];
        #pragma HLS ARRAY_PARTITION variable=etTree complete
        for (unsigned int iInput = 0; iInput < maxObjectsConsidered_; ++iInput){
            #pragma HLS unroll
            etTree[iInput] = passesRCut(
                seedValues[iSeed].range(eta_high_, eta_low_), seedValues[iSeed].range(phi_high_, phi_low_),
                inputObjectValues[iInput].range(eta_high_, eta_low_), inputObjectValues[iInput].range(phi_high_, phi_low_)
            ) ? ap_uint<et_bit_length_ + 11>(inputObjectValues[iInput].range(et_high_, et_low_))
              : ap_uint<et_bit_length_ + 11>(0);
        }
        // Stage 2 (log2(N) levels): strided pairwise reduction; the total lands in etTree[0].
        for (unsigned int stride = 1; stride < maxObjectsConsidered_; stride <<= 1){
            #pragma HLS unroll
            for (unsigned int i = 0; i + stride < maxObjectsConsidered_; i += (stride << 1)){
                #pragma HLS unroll
                etTree[i] += etTree[i + stride];
            }
        }
        ap_uint<et_bit_length_ + 11> outputJetEtWide = etTree[0];

        // Count subjets: original seeds within deltaR < R_cut of the (possibly shifted) output seed center
        for(unsigned int iSubjet = 0; iSubjet < nTotalSeeds_; iSubjet++){
            #pragma HLS unroll // fully unrolled for deterministic (feed-forward) latency
            //std::cout << "iSubjet: " << iSubjet << "\n";
            //std::cout << "seedValuesForSubjets[iSubjet].range(et_high_, et_low_): " << seedValuesForSubjets[iSubjet].range(et_high_, et_low_) << " , subjet_et_threshold_: " << subjet_et_threshold_ << "\n";
            if(seedValuesForSubjets[iSubjet].range(et_high_, et_low_) > subjet_et_threshold_){ // > 25 GeV nominally
                //std::cout << " passes et threshold" << "\n";
                if (passesRCut(
                    seedValues[iSeed].range(eta_high_, eta_low_), seedValues[iSeed].range(phi_high_, phi_low_),
                    seedValuesForSubjets[iSubjet].range(eta_high_, eta_low_), seedValuesForSubjets[iSubjet].range(phi_high_, phi_low_)
                )){
                    //std::cout << "FOUND SUBJET" << "\n";
                    numSubjets = (numSubjets < ap_uint<num_subjets_length_>((1 << num_subjets_length_) - 1))
                                ? ap_uint<num_subjets_length_>(numSubjets + 1)
                                : ap_uint<num_subjets_length_>((1 << num_subjets_length_) - 1);
                }
            }
        }

        // Mask both outputs to 0 for an empty seed (bit-identical to the old guarded form).
        outputJetEt = seedActive
                    ? ((outputJetEtWide >= ((1 << et_bit_length_) - 1))
                        ? ap_uint<et_bit_length_>((1 << et_bit_length_) - 1) // clamp to max Et
                        : ap_uint<et_bit_length_>(outputJetEtWide))
                    : ap_uint<et_bit_length_>(0);
        if(!seedActive) numSubjets = 0;

        //std::cout << "numSubjets: " << numSubjets << "\n";
        outputJetValues[iSeed].range(padded_zeroes_high_, padded_zeroes_low_) = 0;
        outputJetValues[iSeed].range(num_subjets_high_, num_subjets_low_) = numSubjets;
        outputJetValues[iSeed].range(et_high_, et_low_) = outputJetEt;
        outputJetValues[iSeed].range(eta_high_, eta_low_) = seedValues[iSeed].range(eta_high_, eta_low_);
        outputJetValues[iSeed].range(phi_high_, phi_low_) = seedValues[iSeed].range(phi_high_, phi_low_);
    }
}
