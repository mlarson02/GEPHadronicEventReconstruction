#include <string>
#include <sstream>
#include <iomanip>  // for std::setprecision
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>   // for the GEP JwoJ per-eta-bin parameter table
#include "metConstants/constants.h"

// Wrap unsigned phi index into [0, two_pi_digitized_in_phi_]
inline unsigned int wrapPhiUnsigned(unsigned int phi) {
    unsigned int phiWrapped = phi;
    if (phi > two_pi_digitized_in_phi_) phiWrapped = phi - (two_pi_digitized_in_phi_ + 1);
    return phiWrapped;
}

// Compute digitized deltaR^2 between two (eta, phi) coordinates, with phi wrap.
// The fold is around pi_digitized_in_phi_ = phi_range_ / 2, so codes on opposite
// sides of the +-pi seam come out a code apart rather than most of a turn apart.
inline unsigned int digitizedDeltaR2(unsigned int eta1, unsigned int phi1, unsigned int eta2, unsigned int phi2) {
    unsigned int uDeltaEta = static_cast<unsigned int>(std::abs(int(eta1) - int(eta2)));
    unsigned int uDeltaPhi = static_cast<unsigned int>(std::abs(int(phi1) - int(phi2)));
    if (uDeltaPhi >= pi_digitized_in_phi_) uDeltaPhi = (2 * pi_digitized_in_phi_) - uDeltaPhi;
    return uDeltaEta * uDeltaEta + uDeltaPhi * uDeltaPhi;
}

// Digitize phi onto the tower grid.
//
// Kept separate from digitize() because phi is periodic while eta/Et are not:
// the generic function saturates at range - 1, which is wrong at both ends of
// the phi axis. A value in the top half tower (phi within pi/64 of +pi) belongs
// in code 0, not in a code one past the end of the range -- which would not even
// fit the phi field. The code count is phi_range_ rather than
// 1 << phi_bit_length_ so the grid stays tied to the towers, not to the field
// width. Same definition as the jet tagger chain in emulationHelperFunctions.h.
inline unsigned int digitize_phi(double phi) {
    const int nPhi = static_cast<int>(phi_range_);
    const int code = static_cast<int>(std::lround((phi - phi_min_) / phi_granularity_));
    return static_cast<unsigned int>(((code % nPhi) + nPhi) % nPhi);
}

// The reconstruction tag encodes the pileup scenario: PU200 samples are r16130,
// PU140 samples r16129 (HERNTupler stamps this into the ntuple names). Carrying it
// through to the emulation outputs is what keeps the two scenarios distinguishable
// by file name alone, since both land in the same output directory.
// For the makeInputFileName fallback the ntuple directory is switched as well,
// PU140 ntuples living in the parallel ntuples_PU140 tree.
inline std::string applyPileupTags(std::string path, unsigned int pileup) {
    if (pileup != 140) return path;
    auto replaceFirst = [&path](const std::string& from, const std::string& to) {
        size_t pos = path.find(from);
        if (pos != std::string::npos) path.replace(pos, from.size(), to);
    };
    replaceFirst("_r16130_", "_r16129_");
    replaceFirst("/ntuples/", "/ntuples_PU140/");
    return path;
}

// Returns input NTuple file name given parameters
std::string makeInputFileName(bool signalBool, std::string signalString,
                              std::string inputRootFilePath = "/data/larsonma/GEPHadronicEventReconstruction/ntuples/",
                              unsigned int pileup = 200) {
    std::ostringstream ss;

    if (signalBool) {
        if(signalString == "VBF_hh_bbbb_cvv0") ss << inputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv0_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if(signalString == "VBF_hh_bbbb_cvv1") ss << inputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv1_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ggF_hh_bbbb") ss << inputRootFilePath << "ggF_HHbbbb_v3/mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ZvvHbb") ss << inputRootFilePath << "ZvvHbb_v3/mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_resim_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ttbar_had") ss << inputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "Zprime_ttbar") ss << inputRootFilePath << "mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ttbar_semilep") ss << inputRootFilePath << "ttbar_semilep_v4/mc21_14TeV_ttbar_hdamp258p75_semilep_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "ttbar_dilep") ss << inputRootFilePath << "ttbar_dilep_v4/mc21_14TeV_ttbar_hdamp258p75_dilep_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
        else if (signalString == "Zmumu") ss << inputRootFilePath << "Zmumu_v4/mc21_14TeV_Zmumu_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    } else {
        ss << inputRootFilePath << "mc21_14TeV_jj_JZ_e8557_s4422_r16130_DAOD_NTUPLE_GEP.root";
    }
    return applyPileupTags(ss.str(), pileup);
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

// --- GEP JwoJ MET parameters ------------------------------------------------------------
// The GEP JwoJ MET algorithm splits the tower collection in two, in the spirit of the gFEX
// "jets without jets" algorithm but on single towers rather than on gBlocks:
//
//   hard term : towers with E_T strictly above hardEtThreshold
//   soft term : every other tower that survived the tower E_T threshold
//
// and recombines them as
//
//   MET_{x,y} = -( hardCoeff * SUM_hard E_T{x,y} + softCoeff * SUM_soft E_T{x,y} )
//
// The three numbers are held as per-eta-bin vectors rather than as scalars so an
// eta-dependent calibration can be dropped in later without the event loop changing: the
// loop already looks each tower's parameters up by its digitized eta code. What exists today
// is the single-bin table makeUniformGEPJwoJParams builds, which spans the whole eta range
// and so reproduces the scalar behaviour exactly.
//
// The binning is on the digitized eta CODE and not on the double, for the same reason the
// rest of this chain works in codes: a bin boundary has to fall between two towers of the
// grid rather than part way through one, or the firmware and the emulator will disagree
// about which side of it a tower sits on.
struct GEPJwoJParams {
    // Inclusive upper eta code of each bin, ascending. The last entry must be eta_range_ - 1
    // so that every code on the grid lands in a bin.
    std::vector<unsigned int> etaCodeUpperEdges;
    std::vector<double>       hardEtThreshold;   // [GeV], one per bin
    std::vector<double>       hardCoeff;         // one per bin
    std::vector<double>       softCoeff;         // one per bin

    unsigned int nBins() const { return static_cast<unsigned int>(etaCodeUpperEdges.size()); }

    // Bin holding this eta code. A linear scan: the table is one entry today and would be a
    // handful at most, so anything cleverer would cost more than it saves. A code past the
    // last edge -- which digitize() should never produce, since it saturates at range - 1 --
    // falls into the last bin rather than off the end of the vectors.
    unsigned int binForEtaCode(unsigned int etaCode) const {
        for (unsigned int iBin = 0; iBin + 1 < nBins(); ++iBin)
            if (etaCode <= etaCodeUpperEdges[iBin]) return iBin;
        return nBins() - 1;
    }
};

// One eta bin covering the whole grid: what the scalar command-line parameters mean.
inline GEPJwoJParams makeUniformGEPJwoJParams(double hardEtThreshold,
                                              double hardCoeff,
                                              double softCoeff) {
    GEPJwoJParams params;
    params.etaCodeUpperEdges = { eta_range_ - 1 };
    params.hardEtThreshold   = { hardEtThreshold };
    params.hardCoeff         = { hardCoeff };
    params.softCoeff         = { softCoeff };
    return params;
}

// Build an eta-binned table from physical eta upper edges, for when the calibration exists.
// The edges are digitized onto the tower grid here, once, so the event loop never sees a
// double; the last bin is forced to eta_range_ - 1 whatever was passed, so no code can fall
// outside the table.
//
// Nothing calls this yet -- it is the entry point the eta-dependent calibration is meant to
// arrive through, kept beside the uniform builder so the two cannot drift apart.
inline GEPJwoJParams makeEtaBinnedGEPJwoJParams(const std::vector<double>& etaUpperEdges,
                                                const std::vector<double>& hardEtThresholds,
                                                const std::vector<double>& hardCoeffs,
                                                const std::vector<double>& softCoeffs) {
    const size_t nBins = etaUpperEdges.size();
    if (nBins == 0 || hardEtThresholds.size() != nBins ||
        hardCoeffs.size() != nBins || softCoeffs.size() != nBins) {
        std::cerr << "Error: GEP JwoJ eta-binned parameter vectors have mismatched sizes"
                  << " -- falling back to unit coefficients with no hard term" << std::endl;
        return makeUniformGEPJwoJParams(0.0, 1.0, 1.0);
    }
    GEPJwoJParams params;
    params.hardEtThreshold = hardEtThresholds;
    params.hardCoeff       = hardCoeffs;
    params.softCoeff       = softCoeffs;
    params.etaCodeUpperEdges.reserve(nBins);
    for (size_t iBin = 0; iBin < nBins; ++iBin)
        params.etaCodeUpperEdges.push_back(digitize(etaUpperEdges[iBin], eta_bit_length_, eta_min_, eta_max_));
    params.etaCodeUpperEdges.back() = eta_range_ - 1;
    return params;
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
                                  double towerEtThreshold,
                                  bool doJetTowerOverlapRemoval,
                                  std::string outputRootFilePath = "/data/larsonma/GEPMET/outputNTuplesDev_METv3/",
                                  bool useEtaSKObjects = false,
                                  double towerScaleFactor = 1.0,
                                  double jetScaleFactor   = 1.0,
                                  unsigned int pileup     = 200,
                                  bool useGEPJwoJ              = false,
                                  double jwojHardEtThreshold   = 0.0) {
    gSystem->mkdir(outputRootFilePath.c_str());
    std::string usePUSuppress = useEtaSKObjects ? "EtaSK" : (useSKObjects ? "SK" : "NoSK");
    std::string overlapTag    = doJetTowerOverlapRemoval ? "OR" : "NoOR";
    // Format thresholds as integer if whole number, e.g. 20.0 -> "20"
    auto formatThreshold = [](double thr) -> std::string {
        std::ostringstream ss;
        int thrInt = static_cast<int>(thr);
        if (thr == thrInt) ss << thrInt;
        else ss << std::fixed << std::setprecision(1) << thr;
        return ss.str();
    };
    // Format scale factor (assumed granularity 0.1): 1.0 -> "1", 0.4 -> "0p4", 1.2 -> "1p2"
    auto formatScaleFactor = [](double sf) -> std::string {
        std::ostringstream ss;
        int sfInt = static_cast<int>(sf);
        if (sf == sfInt) { ss << sfInt; return ss.str(); }
        ss << std::fixed << std::setprecision(1) << sf;
        std::string s = ss.str();
        for (char& c : s) if (c == '.') c = 'p';
        return s;
    };
    std::string jetThrTag   = formatThreshold(jetEtThreshold);
    std::string towerThrTag = formatThreshold(towerEtThreshold);
    std::string twrSFTag    = formatScaleFactor(towerScaleFactor);
    std::string jetSFTag    = formatScaleFactor(jetScaleFactor);

    std::ostringstream ss;
    if (signalBool) {
        if      (signalString == "VBF_hh_bbbb_cvv0") ss << outputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv0_e8557_s4422_r16130_";
        else if (signalString == "VBF_hh_bbbb_cvv1") ss << outputRootFilePath << "mc21_14TeV_hh_bbbb_vbf_novhh_cvv1_e8557_s4422_r16130_";
        else if (signalString == "ggF_hh_bbbb")      ss << outputRootFilePath << "mc21_14TeV_HHbbbb_HLLHC_e8564_s4422_r16130_";
        else if (signalString == "ZvvHbb")            ss << outputRootFilePath << "mc21_14TeV_ZvvH125_bb_e8557_s4422_r16130_";
        else if (signalString == "ttbar_had")         ss << outputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_allhad_e8557_s4422_r16130_";
        else if (signalString == "Zprime_ttbar")      ss << outputRootFilePath << "mc21_14TeV_flatpT_Zprime_tthad_e8557_s4422_r16130_";
        else if (signalString == "ttbar_semilep")     ss << outputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_semilep_e8557_s4422_r16130_";
        else if (signalString == "ttbar_dilep")       ss << outputRootFilePath << "mc21_14TeV_ttbar_hdamp258p75_dilep_e8557_s4422_r16130_";
        else if (signalString == "Zmumu")             ss << outputRootFilePath << "mc21_14TeV_Zmumu_e8557_s4422_r16130_";
    } else {
        ss << outputRootFilePath << "mc21_14TeV_jj_JZ_e8557_s4422_r16130_";
    }
    ss << "N_Towers_" << maxTowersProcessed
       << "_jetEt"    << jetThrTag
       << "_towerEt"  << towerThrTag
       << "_"         << usePUSuppress
       << "_"         << overlapTag
       << "_twrSF"    << twrSFTag
       << "_jetSF"    << jetSFTag;
    // GEP JwoJ tag, written ONLY when that algorithm is enabled, so every filename produced
    // before it existed comes out byte-identical and nothing already on disk has to move.
    // The hard/soft coefficients are not repeated here: in JwoJ mode they ARE twrSF (soft)
    // and jetSF (hard), so the two tags above already carry them. The hard-term threshold is
    // the only genuinely new number, so it is the only new tag.
    // metAnalysisAndRates.C keys its GEP JwoJ path on the "_GEPJwoJ_" substring.
    if (useGEPJwoJ) ss << "_GEPJwoJ_hardEt" << formatThreshold(jwojHardEtThreshold);
    ss << ".root";
    return applyPileupTags(ss.str(), pileup);
}