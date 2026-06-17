#include <string>
#include <sstream>
#include <iomanip>  // for std::setprecision
#include <fstream>
#include <iostream>
#include <cmath>
#include "UsedConstantsLUTs/constants.h"

inline double undigitize_phi(const std::bitset<phi_bit_length_>& phi_bits) { // FIXME these functions should be templated instead of having 4 copies..
    return phi_min_ + phi_bits.to_ulong() * phi_granularity_;
}

inline double undigitize_eta(const std::bitset<eta_bit_length_>& eta_bits) {
    //std::cout << "eta_min_: " << eta_min_ << " , eta_bits.to_ulong(): " << eta_bits.to_ulong() << " ,  eta_granularity_: " <<  eta_granularity_ << "\n";
    //std::cout << "returns: " << eta_min_ + eta_bits.to_ulong() * eta_granularity_ << "\n";
    return eta_min_ + eta_bits.to_ulong() * eta_granularity_;
}

inline double undigitize_et(const std::bitset<et_bit_length_>& et_bits) {
    return et_bits.to_ulong() * et_granularity_;
}

inline double undigitize_psi_R(const std::bitset<psi_R_bit_length_>& psi_R_bits) {
    return psi_R_bits.to_ulong() * psi_R_granularity_; // FIXME shouldn't be hard-coded
}

inline double undigitize_mass_approx(const std::bitset<mass_approx_bit_length_>& mass_approx_bits) {
    return mass_approx_bits.to_ulong() * mass_approx_granularity_; // FIXME shouldn't be hard-coded
}

inline double undigitize_N_subjetiness(const std::bitset<N_subjetiness_bit_length_>& N_subjetiness_bits) {
    return N_subjetiness_bits.to_ulong() * N_subjetiness_granularity_; // FIXME shouldn't be hard-coded
}


int wrapSym(int phi){
    int phiWrapped = phi;
    if(phi > PI_D) phiWrapped = phi - TWO_PI_D;
    if(phi < - PI_D) phiWrapped = phi + TWO_PI_D;
    return phiWrapped;
}

double wrapSym_dbl(double deltaPhi){
    double deltaPhiWrapped = deltaPhi; 
    if(deltaPhi > M_PI) deltaPhiWrapped -= 2 * M_PI;
    else if(deltaPhi < -M_PI) deltaPhiWrapped += 2 * M_PI;
    return deltaPhiWrapped;
}

// Compute digitized deltaR^2 between two (eta, phi) coordinates, with phi wrap
inline unsigned int digitizedDeltaR2(unsigned int eta1, unsigned int phi1, unsigned int eta2, unsigned int phi2) {
    unsigned int uDeltaEta = static_cast<unsigned int>(std::abs(int(eta1) - int(eta2)));
    unsigned int uDeltaPhi = static_cast<unsigned int>(std::abs(int(phi1) - int(phi2)));
    if (uDeltaPhi >= pi_digitized_in_phi_) uDeltaPhi = (2 * pi_digitized_in_phi_) - uDeltaPhi;
    return uDeltaEta * uDeltaEta + uDeltaPhi * uDeltaPhi;
}

// Compute LUT index from wrapped absolute (deltaEta, deltaPhi) components
inline unsigned int calcLutIndex(unsigned int uDeltaEta, unsigned int uDeltaPhi) {
    return uDeltaEta * (1u << (phi_bit_length_ - 1)) + uDeltaPhi;
}

unsigned int digitize(double value, int bit_length, double min_val, double max_val, unsigned int altRange = 0) {
    unsigned int range = (altRange == 0) ? (1u << bit_length) : altRange;
    double scale = double(range) / (max_val - min_val);
    //std::cout << "max_val - scale; " << max_val - (1/scale) << "\n";
    // Check if value is in range
    if (value < min_val) {
        value = min_val;
        //std::cout << "Warning: Value " << value
        //  << " is out of range (" << min_val
        //  << ", " << max_val << ")\n";
    }
    if (value >= max_val){
        return range - 1;
    }

    return static_cast<unsigned int>(std::round((value - min_val) * scale));
}

double calcDeltaR2(double eta1, double phi1, double eta2, double phi2) {
    //std::cout << "eta1: " << eta1 << " eta2: " << eta2 << " phi1: " << phi1 << " phi2: " << phi2 << "\n";
    double dEta = eta1 - eta2;
    double dPhi = std::abs(phi1 - phi2);
    
    // Ensure Δφ is within [-π, π] range
    if (dPhi > M_PI) {
        dPhi = 2 * M_PI - dPhi;
    }

    return dEta * dEta + dPhi * dPhi;
}

unsigned int index_of_min(unsigned int (&in)[nPreSeeds_]) { // FIXME can't use ap_uint anymore...

    unsigned int min_val = in[0];
    unsigned int min_idx = 0;

    if (in[1] < min_val) {
        min_val = in[1];
        min_idx = 1;
    }
    return min_idx;
}

unsigned int calculate_lut_max_size(double cut,
                                    unsigned int eta_bit_length,
                                    unsigned int phi_bit_length,
                                    double eta_granularity,
                                    double phi_granularity,
                                    bool deltaR2orDeltaRBool // true == computing for DeltaR2, false for DeltaR 
                                    )
{
    unsigned int last_one_index = 0;
    unsigned int idx = 0;

    // Changed to now map pairs of deltaEta, deltaPhi --> deltaR 
    // So now looping through deltaEta up to 20.0
    // deltaPhi up to 6.4 (which is wrapped in this computation, to avoid doing in FW)
    for (unsigned int etaIt = 0; etaIt < (1u << (eta_bit_length)); ++etaIt) {
        for (unsigned int phiIt = 0; phiIt < (1u << (phi_bit_length)); ++phiIt) {
            //std::cout << "etaIt: " << etaIt << " , phiIt: " << phiIt << "\n";
            double deltaPhiWrapped = wrapSym_dbl(phiIt * phi_granularity); 
            double etaSquared = std::pow(etaIt * eta_granularity, 2);
            double phiSquared = std::pow(deltaPhiWrapped, 2);
            //std::cout << "etaSquared: " << etaSquared << " , phiSquared: " << phiSquared << "\n";
            double deltaR2orR = 0; 
            if(deltaR2orDeltaRBool) deltaR2orR = etaSquared + phiSquared;
            else deltaR2orR = std::sqrt(etaSquared + phiSquared);
            //std::cout << "deltaR2orR:  "<< deltaR2orR << "\n";
            //std::cout << "last_one_index: " << last_one_index << "\n";
            if (deltaR2orR < cut) {
                //std::cout << "less than cut" << "\n";
                last_one_index = idx;
                if(cut == 0.001){
                    std::cout << "idx: " << idx << " and last_one_index: " << last_one_index << "\n";
                }
            }
            ++idx;
        }
    }
    
    unsigned int lut_max_size = last_one_index + 1;
    if(cut == 0.001){
        std::cout << "lut_max_size: " << lut_max_size << "\n";
    }
    
    return lut_max_size;
}

// Returns input NTuple file name given parameters
std::string makeInputFileName(bool signalBool, std::string signalString,
                              std::string inputRootFilePath = "/home/larsonma/GEPHadronicEventReconstruction/data/inputNTuples/") {
    std::ostringstream ss;
    
    if (signalBool) {
        if(signalString == "VBF_hh_bbbb_cvv0") ss << inputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv0_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if(signalString == "VBF_hh_bbbb_cvv1") ss << inputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv1_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ggF_hh_bbbb") ss << inputRootFilePath << "mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ZvvHbb") ss << inputRootFilePath << "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ttbar_had") ss << inputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "Zprime_ttbar") ss << inputRootFilePath << "mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";       
    } else {
        ss << inputRootFilePath << "mc21_14TeV_jj_JZ_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    }
    return ss.str();
}

// Returns output NTuple file name given parameters
std::string makeOutputFileName(double rMergeCut,
                               unsigned int NIOs,
                               unsigned int nSeeds,
                               double RSquaredCut,
                               bool signalBool,
                               std::string signalString,
                               std::string inputObjectType,
                               std::string seedObjectType,
                               bool useSKObjects,
                               unsigned int algoVersion,
                               double subjetEtThreshold = 25.0,
                               bool enableEtWeightedMidpoint = false,
                               bool minEtSeedPosOptimization = true,
                               double minEtSeedPosOptimizationCut = 20.0,
                               std::string outputRootFilePath = "/data/larsonma/LargeRadiusJets/outputNTuplesDev_CondorSubmission_NewSamples/",
                               bool useEtaSKObjects = false) {
    gSystem->mkdir(outputRootFilePath.c_str());
    std::string usePUSuppress;
    if(useEtaSKObjects){
        usePUSuppress = "EtaSK";
    }
    else if(useSKObjects){
        usePUSuppress = "SK";
    }
    else{
        usePUSuppress = "NoSK";
    }
    std::ostringstream ss;
    if (signalBool) {
        if(signalString == "VBF_hh_bbbb_cvv0") ss << outputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv0_e8557_s4422_r16130_";
        else if(signalString == "VBF_hh_bbbb_cvv1") ss << outputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv1_e8557_s4422_r16130_";
        else if (signalString == "ggF_hh_bbbb") ss << outputRootFilePath << "mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_";
        else if (signalString == "ZvvHbb") ss << outputRootFilePath << "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_";
        else if (signalString == "ttbar_had") ss << outputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422_r16130_";
        else if (signalString == "Zprime_ttbar") ss << outputRootFilePath << "mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422_r16130_";   
    } else {
        ss << outputRootFilePath << "mc21_14TeV_jj_JZ_e8557_s4422_r16130_";
    }
    
    ss << "rMerge_" << std::setprecision(3) << rMergeCut << "_"
       << "IOs_" << NIOs << "_"
       << "Seeds_" << nSeeds << "_"
       << "R2_" << std::setprecision(3) << RSquaredCut << "_IO_" << inputObjectType << "_Seed_" << seedObjectType << "_" << usePUSuppress
       << "_subjetEt" << static_cast<int>(subjetEtThreshold) << "GeV"
       << "_ewm" << static_cast<int>(enableEtWeightedMidpoint)
       << "_mep" << static_cast<int>(minEtSeedPosOptimization)
       << "_mec" << static_cast<int>(minEtSeedPosOptimizationCut) << "GeV"
       << "_v" << algoVersion << ".root";

    return ss.str();
}

std::string makeOutputTextFileName(double rMergeCut,
                               unsigned int NIOs,
                               unsigned int nSeeds,
                               double RSquaredCut,
                               bool signalBool,
                               std::string signalString,
                               std::string inputObjectType,
                               std::string seedObjectType,
                               bool useSKObjects,
                               unsigned int algoVersion,
                               std::string outputTextFilePath = "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrintsEmulation/",
                               bool useEtaSKObjects = false) {
    std::string usePUSuppress;
    if(useEtaSKObjects){
        usePUSuppress = "EtaSK";
    }
    else if(useSKObjects){
        usePUSuppress = "SK";
    }
    else{
        usePUSuppress = "NoSK";
    }                            
    std::ostringstream ss;
    if (signalBool) {
        if(signalString == "VBF_hh_bbbb_cvv0") ss << outputTextFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv0_e8557_s4422_r16130_";
        else if(signalString == "VBF_hh_bbbb_cvv1") ss << outputTextFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv1_e8557_s4422_r16130_";
        else if (signalString == "ggF_hh_bbbb") ss << outputTextFilePath << "mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_";
        else if (signalString == "ZvvHbb") ss << outputTextFilePath << "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_";
        else if (signalString == "ttbar_had") ss << outputTextFilePath << "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422_r16130_";
        else if (signalString == "Zprime_ttbar") ss << outputTextFilePath << "mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422_r16130_"; 
    } else {
        ss << outputTextFilePath << "mc21_14TeV_jj_JZ_e8557_s4422_r16130_";
    }
    ss << "rMerge_" << std::setprecision(3) << rMergeCut << "_"
       << "IOs_" << NIOs << "_"
       << "Seeds_" << nSeeds << "_"
       << "R2_" << std::setprecision(3) << RSquaredCut << "_IO_" << inputObjectType << "_Seed_" << seedObjectType << "_" << usePUSuppress << "_v" << algoVersion << ".dat";

    return ss.str();
}

// Returns LUT file name given parameters
std::string makeInputLUTFileName(double rMergeCut,
                                 double RSquaredCut,
                                 std::string LUTType,
                                 std::string inputLUTFilePath = "/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation/LUT_Constants_Generation/LUTs/") {
    std::ostringstream ss;
    ss << inputLUTFilePath << "LUT_" << LUTType << "_rMerge_" << std::setprecision(4) << rMergeCut << "_"
       << "R2_" << std::setprecision(3) << RSquaredCut << ".h";

    return ss.str();
}

// Returns Constants.h file name given parameters
std::string makeInputConstantsFileName(double rMergeCut,
                                       unsigned int NIOs,
                                       unsigned int nSeeds,
                                       double RSquaredCut,
                                       unsigned int algoVersion,
                                       std::string constantsFilePath = "/home/larsonma/GEPHadronicEventReconstruction/algorithm/emulation/LUT_Constants_Generation/constants/") {
    std::ostringstream ss;
    ss << constantsFilePath << "constants_rMerge_" << std::setprecision(4) << rMergeCut << "_"
       << "R2_" << std::setprecision(3) << RSquaredCut << "_IOs_" << NIOs << "_"
       << "Seeds_" << nSeeds << "_v" << algoVersion << ".h";

    return ss.str();
}


void write_constants_header(const std::string& header_path,
                            double r2Cut, // nominally 1.21 (R_jet = 1.1) - needs to be re-optimized
                            double rMergeCut, // 0.001 for basic (disables this), between 1.25 - 2.0 for adv 
                            unsigned int maxObjectsConsidered, // nominally 8 for basic, 8 (cone jets) or 128 (towers) for adv 
                            unsigned int nSeedsOutput, // assumed to be always 2 - 2 large-R jets should capture all relevant information in the event based on existing studies, changing it from 2 will likely break things! 
                            unsigned int nSeedsInput, // up to 10, used for subjet finding
                            unsigned int nProtoSeeds, // up to 6, used for seed position optimization
                            unsigned int total_bits_, // 64 for basic, adv for now.. 
                            unsigned int et_bit_length_, // 13 for basic & adv
                            unsigned int eta_bit_length_, // 10 for basic (standard TOB format), 7 for adv (allows 0.1 x 0.1 tower granularity)
                            unsigned int eta_range_, // 784 for basic (4.9 - -4.9 / 0.0125) to match phi granularity of 0.0125, 98 for adv (4.9 - -4.9 / 0.1) to match phi granularity of 0.1
                            unsigned int phi_bit_length_, // 9 for basic (standard TOB format), 6 for adv (allows 0.1 x 0.1 tower granularity)
                            unsigned int max_R_8b_lut_size_, // 8, only needed for advanced algorithm, provides 2 * R_cut / 256 = ~0.09 granularity in deltaR
                            unsigned int substruct_0_bit_length_, // number of subjets for large-R jet - ALREADY IMPLEMENTED IN HLS
                            unsigned int substruct_1_bit_length_, // TBD - for now LRJ 1 subjetiness (tau_1)
                            unsigned int substruct_2_bit_length_, // TBD - for now LRJ 2 subjetiness (tau_2)
                            unsigned int substruct_3_bit_length_, // TBD - for now LRJ Mass approx (deltaR subjets * E_T subjets)
                            unsigned int substruct_4_bit_length_, // TBD - for now  LRJ Psi_R
                            unsigned int algoVersion) // 2 (basic) or 3 (advanced) following naming convention of US ATLAS timeline
{
    std::ofstream out(header_path);
    if (!out) {
        std::cerr << "Failed to open " << header_path << " for writing.\n";
        return;
    }

    // choose a stable formatting for floating-point values in the header
    out.setf(std::ios::fixed);
    out << std::setprecision(4);

    // ---- file begins ----
    out << "// Constants used by JetTagger Emulation\n\n";
    out << "#include <cmath>\n";

    // Fixed (from your template)
    out << "static inline uint32_t maskN(unsigned n) { return (n >= 32) ? 0xFFFFFFFFu : ((1u << n) - 1u); }\n";
    out << "constexpr unsigned int algoVersion_ = " << algoVersion << ";\n";
    out << "constexpr unsigned int nSeedsInput_ = " << nSeedsInput << ";\n";
    out << "constexpr unsigned int nProtoSeeds_ = " << nProtoSeeds << ";\n";
    out << "constexpr unsigned int nSeedsOutput_ = " << nSeedsOutput << ";\n";
    out << "constexpr unsigned int max_R_8b_lut_size_ = " << max_R_8b_lut_size_ << ";\n";
    out << "constexpr unsigned int maxObjectsConsidered_ = " << maxObjectsConsidered << ";\n";
    out << "constexpr double r2Cut_ = " << r2Cut << ";\n";
    out << "constexpr double rCut_ = " << sqrt(r2Cut) << ";\n";
    out << "constexpr double rMergeCut_ = " << rMergeCut << ";\n";
    out << "constexpr unsigned int et_bit_length_ = " << et_bit_length_ << ";\n";
    out << "constexpr unsigned int eta_bit_length_ = " << eta_bit_length_ << ";\n";
    out << "constexpr unsigned int phi_bit_length_ = " << phi_bit_length_ << ";\n";
    out << "constexpr unsigned int num_subjets_length_ = " << substruct_0_bit_length_ << ";\n";
    out << "constexpr unsigned int mass_approx_bit_length_ = " << substruct_3_bit_length_ << ";\n";
    out << "constexpr unsigned int N_subjetiness_bit_length_ = " << substruct_1_bit_length_ << ";\n";
    out << "constexpr unsigned int psi_R_bit_length_ = " << substruct_4_bit_length_ << ";\n";
    out << "constexpr unsigned int deltaR_lut_length_ = 8;\n"; // maybe this shouldn't be hard-coded? 
    out << "constexpr unsigned int padded_zeroes_length_ = 64 - et_bit_length_ - eta_bit_length_ - phi_bit_length_ - num_subjets_length_;\n"; // For now only include substruct 0 - only this is included in HLS currently!
    out << "constexpr unsigned int total_bits_output_ = padded_zeroes_length_ + num_subjets_length_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;\n";
    out << "constexpr double phi_min_ = -3.2;\n";
    out << "constexpr double phi_max_ = 3.2;\n";
    out << "constexpr double eta_min_ = -4.9;\n";
    out << "constexpr double eta_max_ = 4.9;\n";
    out << "constexpr unsigned int et_min_ = 0;\n";
    out << "constexpr unsigned int et_max_ = 1024;\n";
    out << "constexpr double phi_granularity_ = (phi_max_ - phi_min_) / (1 << (phi_bit_length_));\n";
    out << "constexpr unsigned int pi_digitized_in_phi_ = (M_PI) / (phi_granularity_);\n";
    out << "const int PI_D     =  int(pi_digitized_in_phi_);\n"; // 31
    out << "const int TWO_PI_D =  2 * pi_digitized_in_phi_ + 1;\n"; // wrap period (63 for 6 bits)
    out << "constexpr unsigned int eta_range_ = (eta_max_ - eta_min_) / (phi_granularity_);\n"; // basically enforcing that eta granularity matches phi granularity 
    out << "constexpr double eta_granularity_ = (eta_max_ - eta_min_) / double(eta_range_);\n";
    out << "constexpr double et_granularity_ = (et_max_ - et_min_) / double((1 << et_bit_length_));\n";
    out << "constexpr double deltaR_granularity_ = (2 * rCut_) / double((1 << deltaR_lut_length_));\n"; 
    out << "constexpr double psi_R_granularity_ = (2 * rCut_) / double((1 << deltaR_lut_length_)) ;\n";
    out << "constexpr double massApproxRawLSB_GeV_  = et_granularity_ * deltaR_granularity_ ;\n";
    out << "constexpr double massApprox_max_  = 512;\n"; // in GeV
    out << "constexpr double massApproxNewLSB_GeV_ = massApprox_max_ / (1 << mass_approx_bit_length_);\n";
    out << "constexpr double massApproxScaleFactor_ = massApproxNewLSB_GeV_ / massApproxRawLSB_GeV_;\n"; 
    out << "constexpr unsigned int massApproxDivisor_ = static_cast<unsigned int>(massApproxScaleFactor_ + 0.5);\n"; // rounded
    out << "constexpr double mass_approx_granularity_ = massApprox_max_ / (1 << mass_approx_bit_length_);\n"; // 2 GeV
    out << "constexpr double N_subjetiness_granularity_ = (1.0) / (1 << N_subjetiness_bit_length_) ;\n"; // FIXME replace these with correct values! 
    out << "constexpr double phi_range_ = phi_max_ - phi_min_;\n";
    out << "constexpr double deltaR2_granularity_ = eta_granularity_ * eta_granularity_;\n"; // should be the same as phi_granularity_ squared too
    out << "constexpr unsigned int digitized_delta_R2Cut_ = static_cast<unsigned int>(r2Cut_/deltaR2_granularity_ + 0.5);\n";
    out << "constexpr unsigned int digitized_d_search_squared_ = static_cast<unsigned int>(((rMergeCut_) * (rMergeCut_))/deltaR2_granularity_ + 0.5);\n";
    out << "constexpr unsigned int total_bits_ = padded_zeroes_length_ + num_subjets_length_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;\n";
    out << "constexpr unsigned int phi_low_  = 0;\n";
    out << "constexpr unsigned int phi_high_ = phi_low_ + phi_bit_length_ - 1;\n";

    out << "constexpr unsigned int eta_low_  = phi_high_ + 1;\n";
    out << "constexpr unsigned int eta_high_ = eta_low_ + eta_bit_length_ - 1;\n";

    out << "constexpr unsigned int et_low_   = eta_high_ + 1;\n";
    out << "constexpr unsigned int et_high_  = et_low_ + et_bit_length_ - 1;\n";

    out << "constexpr unsigned int substruct_0_bit_low_  = et_high_ + 1;\n";
    out << "constexpr unsigned int substruct_0_bit_high_ = substruct_0_bit_low_ + num_subjets_length_ - 1;\n";

    //out << "constexpr unsigned int substruct_1_bit_low_  = substruct_0_bit_high_ + 1;\n";
    //out << "constexpr unsigned int substruct_1_bit_high_ = substruct_1_bit_low_ + substruct_1_bit_length_ - 1;\n";

    out << "constexpr unsigned int padded_zeroes_low_  = substruct_0_bit_high_ + 1;\n";
    out << "constexpr unsigned int padded_zeroes_high_ = padded_zeroes_low_ + padded_zeroes_length_ - 1;\n";

    out << "constexpr unsigned int nPreSeeds_ = nProtoSeeds_ - nSeedsOutput_;\n"; // FIXME unhard code this to account for when using cone jets! 

    out << "static const unsigned int lutR_8b_[max_R_8b_lut_size_] = \n";
    out << "#include \"deltaRLUT_8b.h\"\n";
    out << ";\n\n";

    // ---- file ends ----
    std::cout << "writing constants file to: " << header_path << "\n";
    out.close();
}