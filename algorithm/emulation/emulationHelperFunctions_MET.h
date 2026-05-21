#include <string>
#include <sstream>
#include <iomanip>  // for std::setprecision
#include <fstream>
#include <iostream>
#include <cmath>
#include "metConstants/constants.h"

// Wrap unsigned phi index into [0, two_pi_digitized_in_phi_]
inline unsigned int wrapPhiUnsigned(unsigned int phi) {
    unsigned int phiWrapped = phi;
    if (phi > two_pi_digitized_in_phi_) phiWrapped = phi - (two_pi_digitized_in_phi_ + 1);
    return phiWrapped;
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

// Sign-magnitude: MSB = sign bit, remaining signed_et_bit_length_-1 bits = magnitude
template<size_t signed_et_bit_length_>
inline double undigitize_signed_et(const std::bitset<signed_et_bit_length_>& bits) {
    std::cout << "bits: " << bits << "\n";
    int sign = bits[signed_et_bit_length_ - 1] ? -1 : 1;
    std::cout << "sign: " << sign << "\n";
    int mag = static_cast<int>(bits.to_ulong()) & ((1u << (signed_et_bit_length_ - 1)) - 1);
    std::cout << "mag: " << mag << "\n";
    std::cout << "output: " << sign * mag * et_granularity_ << "\n";
    return sign * mag * et_granularity_;
}

inline double undigitize_et(const std::bitset<et_bit_length_>& et_bits) {
    return et_bits.to_ulong() * et_granularity_;
}

// Returns output NTuple file name for MET emulation
std::string makeOutputMETFileName(unsigned int maxTowersProcessed,
                                  bool signalBool,
                                  std::string signalString,
                                  bool useSKObjects,
                                  double jetEtThreshold,
                                  std::string outputRootFilePath = "/data/larsonma/GEPMET/outputNTuplesDev_METv1/") {
    gSystem->mkdir(outputRootFilePath.c_str());
    std::string usePUSuppress = useSKObjects ? "SK" : "NoSK";
    // Format threshold as integer if it is a whole number, e.g. 20.0 -> "20"
    std::ostringstream thrSS;
    int thrInt = static_cast<int>(jetEtThreshold);
    if (jetEtThreshold == thrInt) thrSS << thrInt;
    else thrSS << std::fixed << std::setprecision(1) << jetEtThreshold;
    std::string thrTag = thrSS.str();

    std::ostringstream ss;
    if (signalBool) {
        if      (signalString == "VBF_hh_bbbb_cvv0") ss << outputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv0_e8557_s4422_r16130_";
        else if (signalString == "VBF_hh_bbbb_cvv1") ss << outputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv1_e8557_s4422_r16130_";
        else if (signalString == "ggF_hh_bbbb")      ss << outputRootFilePath << "mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_";
        else if (signalString == "ZvvHbb")            ss << outputRootFilePath << "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_";
        else if (signalString == "ttbar_had")         ss << outputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422_r16130_";
        else if (signalString == "Zprime_ttbar")      ss << outputRootFilePath << "mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422_r16130_";
    } else {
        ss << outputRootFilePath << "mc21_14TeV_jj_JZ_e8557_s4422_r16130_";
    }
    ss << "N_Towers_" << maxTowersProcessed << "_jetEt" << thrTag << "_" << usePUSuppress << ".root";
    return ss.str();
}