from pathlib import Path
#import vitis
import os
import re
import shutil
import math
import xml.etree.ElementTree as ET
import subprocess

REPO_ROOT = Path(__file__).resolve().parent.parent
EOS_DATA_ROOT = Path("/eos/user/m/mlarson/TransferMemPrintsLUTs/data")

# Base/default constants
base_constants = {
    "nTotalSeeds_": 10,
    "nSeedsInput_": 6,
    "nSeedsOutput_": 2,
    "maxObjectsConsidered_": 128,
    "et_granularity_": 0.125,
    "subjet_et_threshold_": 200, # == 25 GeV 
    "r2Cut_": 1.21,
    "rCut_": 1.1,
    "rMergeCut_": 2.0,
    "et_bit_length_": 13,
    "eta_bit_length_": 7,
    "phi_bit_length_": 6,
    "eta_range_": 98, # the range of allowed eta values - for giving the correct tower granularity (9.8 / 98 = 0.1)
    "num_subjets_length_": 2, # note: with latest format, have 36 free bits! 
    "deltaRBits_": 8,
    "phi_min_": -3.2,
    "phi_max_": 3.2,
    "pi_digitized_in_phi_": 31,
    "eta_min_": -4.85,
    "eta_max_": 4.95,
    "eta_granularity_": 0.1,
    "phi_granularity_": 0.1,
    "et_min_": 0,
    "et_max_": 1024
}

def compute_derived(constants):
    constants["phi_range_"] = constants["phi_max_"] - constants["phi_min_"]
    return constants

import os

def write_file_read_header(file_path, file_suffix, signal_bool, jzSlice):
    # Remove existing constants.h if it exists
    if os.path.exists(file_path):
        os.remove(file_path)
    # Base content up to lutPath_
    header_content = f"""#ifndef FILE_READ_H  // Check if the macro is defined
#define FILE_READ_H  // Define the macro
#include <iostream>
#include <fstream>
#include <sstream>
#include <bitset>
#include <array>
#include <ap_int.h>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

// Define constants used by testbench
const std::string memPrintsPath_ = "{EOS_DATA_ROOT}/MemPrints_v3/";
static inline uint32_t maskN(unsigned n) {{ return (n >= 32) ? 0xFFFFFFFFu : ((1u << n) - 1u); }}
"""

    # Add your file suffix line *after* lutPath_
    suffix_line = f'const std::string kFileSuffix = "{file_suffix}";\n'
    # Now create the signalBool_ line dynamically
    signal_bool_line = f'constexpr bool signalBool_ = {"true" if signal_bool else "false"};\n'

    jzSlice_line = f'constexpr unsigned int jzSlice_ = {jzSlice};\n'

    # Now continue with the rest of the file content
    #mc21_14TeV_hh_bbbb_vbf_novhh
    #mc21_14TeV_flatpT_Zprime_tthad
    remaining_content = """
const unsigned int maxEvent_ = signalBool_ ? 100 : 100;
const std::string fileName_ =
    signalBool_        ? "mc21_14TeV_HHbbbb_HLLHC" :
    (jzSlice_ == 2)    ? "mc21_14TeV_jj_JZ2" :
    (jzSlice_ == 3)    ? "mc21_14TeV_jj_JZ3" :
    (jzSlice_ == 4)    ? "mc21_14TeV_jj_JZ4" :
                         "unknown_sample";

// read values from .dat files for a provided event
template <unsigned int arraySize >
inline void extract_values_from_file(const std::string& fileName, input (&values)[arraySize], unsigned int eventToProcess) {

    std::ifstream inFile(fileName);
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return;
    }

    input zero = 0;
    std::fill(std::begin(values), std::end(values), zero);
    int iEvt = -1;
    std::string line;
    input valuesForEvent[arraySize];
    int lastEvent = 0;
    int objectIt = -1;
    int eventNumber = -1;
    int prevEventNumber = -1;
    while (std::getline(inFile, line)) {
        if (line.find("Event") != std::string::npos) {
            std::stringstream ss0(line);
            std::string temp;
            prevEventNumber = eventNumber;
            ss0 >> temp >> temp >> eventNumber;
            if (iEvt >= 0 && prevEventNumber == (int)eventToProcess){
                for (unsigned int i = 0; i < arraySize; ++i) {
                    values[i] = valuesForEvent[i];
                }
                break;
            }
            iEvt++;
            input zero = 0;
            //std::cout << "ievt in fileread: " << iEvt << std::endl;
            std::fill(std::begin(valuesForEvent), std::end(valuesForEvent), zero);
            continue;
        }
        if (iEvt > eventToProcess) break;
        if (eventNumber == eventToProcess){
            if (lastEvent != eventToProcess) objectIt = -1;
            lastEvent = eventToProcess;
            objectIt++; 
            
            std::stringstream ss(line);
            std::string index, bin, hex_word;
            ss >> index >> bin >> hex_word;

            //std::cout << "hex_word: " << hex_word << std::endl;

            size_t first_pipe = bin.find('|');
            size_t second_pipe = bin.rfind('|');
            if (first_pipe == std::string::npos || second_pipe == std::string::npos || first_pipe == second_pipe) {
                std::cerr << "Error: Malformed line -> " << line << std::endl;
                continue;
            }

            // Binary column is written MSB -> LSB as phi | eta | et, i.e. et occupies the LSBs
            std::string phi_bin = bin.substr(0, first_pipe);
            std::string eta_bin = bin.substr(first_pipe + 1, second_pipe - first_pipe - 1);
            std::string et_bin  = bin.substr(second_pipe + 1);
            //std::cout << "phi_bin : " << phi_bin << std::endl;
            //std::cout << "eta_bin : " << eta_bin << std::endl;
            //std::cout << "et_bin : " << et_bin << std::endl;

            // Prepend 5 zero bits (as MSB) to represent num_io = 0
            std::string num_io_bin = "00000";

            // Final bitstring in MSB to LSB order: num_io | phi | eta | et
            std::string full_bin = num_io_bin + phi_bin + eta_bin + et_bin;

            //std::cout << "full_bin: " << full_bin << std::endl;

            input fullInput = ap_uint<input::width>(std::stoull(full_bin, nullptr, 2));


            if (objectIt >= arraySize) continue;
            //std::cout << "full input: " << fullInput << std::endl;
            valuesForEvent[objectIt]  = fullInput;
        }
    }
    if (eventNumber == (int)eventToProcess){
        for (unsigned int i = 0; i < arraySize; ++i) {
            values[i] = valuesForEvent[i];
        }
    }

    inFile.close();
    return;
}

#endif
"""

    # Now write everything to file
    with open(file_path, 'w') as f:
        f.write(header_content)
        f.write(suffix_line)
        f.write(signal_bool_line)
        f.write(jzSlice_line)
        f.write(remaining_content)

    print(f"File {file_path} written successfully with the file suffix.")




LUT_V3_DIR = REPO_ROOT / "algorithm" / "emulation" / "LUT_Constants_Generation" / "LUTs" / "v3"


def find_lut_file(lut_type: str, r2cut: float, rmergecut: float, tol: float = 1e-6) -> Path:
    # Match by parsed numeric value rather than reconstructing the filename, since C++'s default
    # ostream float formatting (e.g. "2" not "2.0") is finicky to replicate exactly in Python.
    pattern = re.compile(rf"^LUT_{re.escape(lut_type)}_rMerge_(?P<rmerge>[0-9.]+)_R2_(?P<r2>[0-9.]+)\.h$")
    for path in LUT_V3_DIR.glob(f"LUT_{lut_type}_rMerge_*_R2_*.h"):
        m = pattern.match(path.name)
        if not m:
            continue
        if abs(float(m.group("r2")) - r2cut) < tol and abs(float(m.group("rmerge")) - rmergecut) < tol:
            return path
    raise FileNotFoundError(
        f"No {lut_type} LUT found in {LUT_V3_DIR} for r2Cut_={r2cut}, rMergeCut_={rmergecut}. "
        f"Regenerate via writeLUTsConstantsEmulation.cc first."
    )


def count_lut_entries(path: Path) -> int:
    body = path.read_text().strip()
    body = body[body.index("{") + 1 : body.rindex("}")].strip().rstrip(",")
    return len([v for v in body.split(",") if v.strip() != ""])


def write_constants_h(constants: dict, output_file: str, unroll: int, ii: int, use_dsps: bool = True):
    # Remove existing constants.h if it exists
    if os.path.exists(output_file):
        os.remove(output_file)

    deltaR2_bits = 10
    r2cut = constants["r2Cut_"]
    rmergecut = constants["rMergeCut_"]

    deltaR2Cut_lut_path = find_lut_file("deltaR2Cut", r2cut, rmergecut)
    deltaR_lut_path = find_lut_file("deltaR", r2cut, rmergecut)
    max_R2lut_size = count_lut_entries(deltaR2Cut_lut_path)
    max_Rlut_size = count_lut_entries(deltaR_lut_path)

    # Mirrors the generator's own deltaR_max/deltaR_granularity computation exactly, so digitized
    # linear thresholds compared against lutR_ entries line up with how that table was digitized.
    eta_range = constants["eta_range_"]
    eta_granularity = constants["eta_granularity_"]
    phi_bit_length = constants["phi_bit_length_"]
    phi_granularity = constants["phi_granularity_"]
    deltaR_max = math.sqrt(
        ((eta_range - 1) * eta_granularity) ** 2
        + (((1 << (phi_bit_length - 1)) - 1) * phi_granularity) ** 2
    )
    deltaR_granularity = deltaR_max / 255.0
    digitized_two_rCut = round(2 * constants["rCut_"] / deltaR_granularity)
    digitized_rMergeCut = round(rmergecut / deltaR_granularity)

    with open(output_file, "w") as f:
        f.write("#ifndef CONSTANTS_ADV_H\n")
        f.write("#define CONSTANTS_ADV_H\n")
        f.write("// Constants used by SW & FW implementation\n\n")

        # Write UNROLLFACTOR and PIPELINEII
        f.write(f"#define UNROLLFACTOR {unroll}\n")
        f.write(f"#define PIPELINEII {ii}\n\n")


        # Write all constants
        for key, value in constants.items():
            if isinstance(value, bool) or value in ["true", "false"]:
                f.write(f"#define {key} {'true' if str(value).lower() == 'true' else 'false'}\n")
            elif isinstance(value, int):
                f.write(f"constexpr unsigned int {key} = {value};\n")
            elif isinstance(value, float):
                f.write(f"constexpr double {key} = {value};\n")
            else:
                f.write(f"constexpr auto {key} = {value};\n")

        # Write the extra constants, LUT include, and struct
        f.write(f'''

constexpr unsigned int padded_zeroes_length_ = 64 - et_bit_length_ - eta_bit_length_ - phi_bit_length_ - num_subjets_length_ - num_subjets_length_;
constexpr unsigned int padded_zeroes_length_32b_ = 32 - et_bit_length_ - eta_bit_length_ - phi_bit_length_;
constexpr unsigned int total_bits_input_ = padded_zeroes_length_32b_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
constexpr unsigned int total_bits_output_ = padded_zeroes_length_ + num_subjets_length_ + num_subjets_length_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
typedef ap_uint<total_bits_input_> input; // need 32b input, 64b output!
typedef ap_uint<total_bits_output_> output;

// MSB -> LSB word order is phi | eta | et, i.e. et occupies the LSBs and phi the MSBs
constexpr unsigned int et_low_   = 0;
constexpr unsigned int et_high_  = et_low_ + et_bit_length_ - 1;

constexpr unsigned int eta_low_  = et_high_ + 1;
constexpr unsigned int eta_high_ = eta_low_ + eta_bit_length_ - 1;

constexpr unsigned int phi_low_  = eta_high_ + 1;
constexpr unsigned int phi_high_ = phi_low_ + phi_bit_length_ - 1;

constexpr unsigned int num_subjets_low_  = phi_high_ + 1;
constexpr unsigned int num_subjets_high_ = num_subjets_low_ + num_subjets_length_ - 1;

constexpr unsigned int padded_zeroes_low_  = num_subjets_high_ + 1;
constexpr unsigned int padded_zeroes_high_ = padded_zeroes_low_ + padded_zeroes_length_ - 1;

constexpr unsigned int nSeedsDeltaR_ = nSeedsInput_ - nSeedsOutput_;

constexpr unsigned int deltaR2_bits_  = {deltaR2_bits}; // bit width of the saturated deltaR2 result (values > {(1 << deltaR2_bits) - 1} clamp to {(1 << deltaR2_bits) - 1})
constexpr double deltaR2_granularity_ = eta_granularity_ * eta_granularity_; // FIXME THIS SHOULD BE EQUIVALENT TO SQUARING PHI_GRANULARITY_ - maybe add an exception if they are not the same

constexpr unsigned int digitized_delta_R2Cut_ = static_cast<unsigned int>(r2Cut_/deltaR2_granularity_ + 0.5); //+ 0.5 for correct rounding

constexpr unsigned int digitized_d_search_squared_ = static_cast<unsigned int>(((rMergeCut_) * (rMergeCut_)/(deltaR2_granularity_)) + 0.5);

// Controls whether calcDeltaR2 (helperFunctions_adv.h) computes deltaR^2 via DSP multipliers (baseline)
// or via a precomputed LUT indexed by (deltaEta, deltaPhi). Set to 0 to use the LUT path instead.
#define USE_DSPS_ {1 if use_dsps else 0}

#if !USE_DSPS_
// Boolean pass/fail against r2Cut_, indexed by (deltaEta, deltaPhi). Used by passesRCut().
constexpr unsigned int max_R2lut_size_ = {max_R2lut_size};
static const bool lut_[max_R2lut_size_] =
#include "{deltaR2Cut_lut_path}"
;

// Digitized deltaR (8-bit, full grid coverage, cut-independent) indexed by (deltaEta, deltaPhi).
// Used by passesTwoRCut() and passesSearchRadius() against the digitized thresholds below.
constexpr unsigned int max_Rlut_size_ = {max_Rlut_size};
static const ap_uint<deltaRBits_> lutR_[max_Rlut_size_] =
#include "{deltaR_lut_path}"
;

constexpr double deltaR_granularity_ = {deltaR_granularity!r};
constexpr unsigned int digitized_two_rCut_ = {digitized_two_rCut};
constexpr unsigned int digitized_rMergeCut_ = {digitized_rMergeCut};
#endif

        ''')

        f.write("\n#endif // CONSTANTS_ADV_H\n")


if __name__ == "__main__":
    # Sweep over combinations
    #nSeeds_options = [1, 2, 3, 4, 5, 6, 7, 8]
    nSeeds_options = [2]
    #r2Cut_options = [0.8, 0.9, 1.0, 1.1, 1.2]
    #r2Cut_options = [1.0, 1.21, 1.44, 1.69, 1.96]
    #r2Cut_options = [1.69]
    r2Cut_options = [1.21]
    #r2Cut_options = [1.44]
    #maxObjectsConsidered_options = [128, 256, 512, 1024]
    maxObjectsConsidered_options = [128]
    #maxObjectsConsidered_options = [128]
    rMergeCut_options = [2.0]
    #rMergeCut_options = [3.5]
    signalBool_options = [True]

    use_dsps = True # baseline is DSP-based calcDeltaR2; set False to use the LUT path instead

    jzSlices = [3]

    for jzSlice in jzSlices:
        print("jzSlice processed:", jzSlice)
        for maxObjectsConsidered in maxObjectsConsidered_options:
            print("maxObjectsConsidered processed:", maxObjectsConsidered)
            for rMergeCut in rMergeCut_options:
                print("rMergeCut processed:", rMergeCut)
                for nSeeds in nSeeds_options:
                    print("nSeeds processed:", nSeeds)
                    for r2Cut in r2Cut_options:
                        print("r2Cut processed:", r2Cut)
                        for signalBool in signalBool_options:
                            print("signalBool processed:", signalBool)
                            # Copy base constants
                            constants = base_constants.copy()
                            constants["nSeedsOutput_"] = nSeeds
                            constants["r2Cut_"] = r2Cut
                            constants["rMergeCut_"] = rMergeCut
                            constants["maxObjectsConsidered_"] = maxObjectsConsidered
                            print("nSeeds:", nSeeds)
                            print("r2Cut:", r2Cut)
                            print("rMergeCut:", rMergeCut)
                            print("maxobjectsconsidered:", maxObjectsConsidered)

                            # Calculate phi and eta granularities from base constants
                            phi_range = base_constants["phi_max_"] - base_constants["phi_min_"]
                            phi_granularity = phi_range / (1 << base_constants["phi_bit_length_"])
                            constants["phi_granularity_"] = phi_granularity

                            eta_undigi_range = base_constants["eta_max_"] - base_constants["eta_min_"]
                            eta_granularity = eta_undigi_range / (base_constants["eta_range_"])
                            #eta_granularity = eta_range / (1 << base_constants["eta_bit_length_"])
                            constants["eta_granularity_"] = eta_granularity

                            # Update derived values
                            constants = compute_derived(constants)

                            # Clean up float formatting for r2Cut
                            r2Cut_str = str(r2Cut).replace('.', 'p')
                            rMergeCut_str = str(rMergeCut).replace('.', 'p')

                            # Bool formatting for sortSeeds and inputEnergyCut
                            if signalBool:
                                signal_str = 'sig'
                            else:
                                if (jzSlice == 2):
                                    signal_str = 'back_JZ2'
                                elif (jzSlice == 3):
                                    signal_str = 'back_JZ3'
                                elif (jzSlice == 4):
                                    signal_str = 'back_JZ4'

                            constsFilename = f"../algorithm/constants_adv.h"

                            # File name based on parameter values
                            file_suffix = (
                                f"nSeeds{nSeeds}_"
                                f"r2Cut{r2Cut_str}_"
                                f"maxObj{maxObjectsConsidered}_"
                                f"rMerge{rMergeCut_str}_"
                                f"{signal_str}"
                                "_WTAConeJetsCellsTowers_Adv_ValidateEmulation_FINAL_MODIFICATIONS_DETERMINISTIC_TREEFIX"
                            )
                            print(f"Launching HLS with project name: {file_suffix}")

                            unroll = int(maxObjectsConsidered/8)

                            ii = 3

                            # Write to file
                            write_constants_h(constants, constsFilename, unroll, ii, use_dsps)
                            fileReadPath = str(REPO_ROOT / "algorithm" / "fileRead.h")  # Path to save the file
                            write_file_read_header(fileReadPath, file_suffix, signalBool, jzSlice)

                            print(f" Wrote {constsFilename}")
                            #run_lut_generator_via_root(str(REPO_ROOT / "algorithm" / "writeDeltaR2LUT_adv.cc"))

                            subprocess.run(["vitis", "-s", "jet_tagger_hls_adv.py", file_suffix, "1"], check=True)
                            xml_report_path = os.path.join('w', file_suffix, file_suffix, 'syn', 'report', 'jet_tagger_top_csynth.xml')
                            print("xml_report_path,", xml_report_path)
                                    #resources, latency = extract_hls_report(xml_report_path)
                                #print(f"Run {file_suffix}: Resources {resources}, Latency {latency}")

