from pathlib import Path
#import vitis
import os
import shutil
import math
import ROOT
import xml.etree.ElementTree as ET
import subprocess

def wrapSym(phi):
    if(phi > math.pi):
        phi -= 2 * math.pi
    elif(phi < -math.pi):
        phi += 2 * math.pi
    return phi

# Base/default constants
base_constants = {
    "sortSeeds_": "true",
    "iterative_": "false",
    "nTotalSeeds_": 200,
    "nSeedsInput_": 6,
    "nSeedsOutput_": 2,
    "maxObjectsConsidered_": 1024,
    "inputEnergyCut_": 1,
    "useInputEnergyCut_": "false", 
    "et_granularity_": 0.125,
    "r2Cut_": 1.0,
    "rCut_": 1.0,
    "rMergeCut_": 5.0,
    "et_bit_length_": 13,
    "eta_bit_length_": 10,
    "eta_range_": 784, 
    "phi_bit_length_": 9,
    "phi_min_": -3.2,
    "phi_max_": 3.2,
    "pi_digitized_in_phi_": 251, # rounding 3.1415/0.0125
    "eta_min_": -4.9,
    "eta_max_": 4.9,
    "eta_granularity_": 0.0125,
    "phi_granularity_": 0.0125, # fixme these should be calculated from range and bit length not hard-coded!
    "deltaR2_granularity_": 0.00015625,
    "et_min_": 0,
    "et_max_": 1024,
    "useMax_": "false",
    "max_R2lut_size_": 2048,
    "max_Rlut_size_": 2048,
    "deltaR_max_": 10.48187,
    "deltaR_bits_": 8
}

def calculate_lutR2_max_size(r2Cut, eta_range, phi_bit_length, eta_granularity, phi_granularity):
    last_one_index = 0
    idx = 0
    #print("r2cut used in calculating max lut size: ", r2Cut)
    for etaIt in range(eta_range): # now looping through for deltaEta, deltaPhi (before wrapped)
        for phiIt in range(1 << (phi_bit_length)):
            #print("etaIt:", etaIt)
            #print("phiIt:", phiIt)
            etaSquared = (etaIt * eta_granularity) ** 2
            phiSquared = (wrapSym(phiIt * phi_granularity)) ** 2
            #print("etaSquared:", etaSquared)
            #print("phiSquared:", phiSquared)
            
            deltaR2 = etaSquared + phiSquared
            #print("deltaR2:", deltaR2)
            #print("deltaR2:", deltaR2)
            #print("r2cut:", r2Cut)
            if deltaR2 < r2Cut:
                last_one_index = idx
                #print("last_one_index:", last_one_index)
            idx += 1
            #print("last_one_index:", last_one_index)
    lut_max_size = last_one_index + 1
    #print("lut_max_size:", lut_max_size)
    #print("lut_max_size:", lut_max_size)
    return lut_max_size

def calculate_lutR_max_size(rCut, eta_range, phi_bit_length, eta_granularity, phi_granularity):
    last_one_index = 0
    idx = 0
    #print("r2cut used in calculating max lut size: ", r2Cut)
    for etaIt in range(eta_range):
        for phiIt in range(1 << (phi_bit_length)):
            etaSquared = (etaIt * eta_granularity) ** 2
            phiSquared = (wrapSym(phiIt * phi_granularity)) ** 2
            #print("etaSquared:", etaSquared)
            #print("phiSquared:", phiSquared)
            
            deltaR = math.sqrt(etaSquared + phiSquared)
            #print("deltaR2:", deltaR2)

            if deltaR < rCut:
                last_one_index = idx
            idx += 1
            #print("last_one_index:", last_one_index)
    
    lut_max_size = last_one_index + 1
    #print("lut_max_size:", lut_max_size)
    return lut_max_size

def compute_derived(constants):
    constants["phi_range_"] = constants["phi_max_"] - constants["phi_min_"]
    return constants

def run_lut_generator_via_root(source_file: str):
    ROOT.gSystem.Load("libCore")  # Make sure libCore is loaded
    ROOT.gROOT.ProcessLine(f'.L {source_file}')
    ROOT.gROOT.ProcessLine('main();')  # Assuming your main is called "main"

def extract_hls_report(xml_path):
    if not os.path.exists(xml_path):
        raise FileNotFoundError(f"Report file not found: {xml_path}")

    tree = ET.parse(xml_path)
    root = tree.getroot()

    resources = {}
    latency = {}

    # Extract Resource Utilization
    for summary in root.findall(".//SummaryOfResources"):
        resources['LUT'] = int(summary.findtext('LUT', default='0'))
        resources['FF'] = int(summary.findtext('FF', default='0'))
        resources['BRAM_18K'] = int(summary.findtext('BRAM_18K', default='0'))
        resources['DSP48E'] = int(summary.findtext('DSP48E', default='0'))

    # Extract Latency
    for timing in root.findall(".//SummaryOfTiming/Latency"):
        latency['BestCaseCycles'] = int(timing.findtext('Best-case', default='0'))
        latency['WorstCaseCycles'] = int(timing.findtext('Worst-case', default='0'))

    return resources, latency


import os

def write_file_read_header(file_path, file_suffix, signal_bool, jzSlice):
    # Remove existing constants.h if it exists
    if os.path.exists(file_path):
        os.remove(file_path)
    # Base content up to lutPath_
    header_content = """#ifndef FILE_READ_H  // Check if the macro is defined
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
const std::string memPrintsPath_ = "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints_v2/";
const std::string lutPath_ = "/home/larsonma/GEPHadronicEventReconstruction/data/LUTs/deltaR2Cut.dat";
static inline uint32_t maskN(unsigned n) { return (n >= 32) ? 0xFFFFFFFFu : ((1u << n) - 1u); }
"""

    # Add your file suffix line *after* lutPath_
    suffix_line = f'const std::string kFileSuffix = "{file_suffix}";\n'
    # Now create the signalBool_ line dynamically
    signal_bool_line = f'constexpr bool signalBool_ = {"true" if signal_bool else "false"};\n'

    jzSlice_line = f'constexpr unsigned int jzSlice_ = {jzSlice};\n'
    #mc21_14TeV_HHbbbb_HLLHC
    #mc21_14TeV_flatpT_Zprime_tthad
    # Now continue with the rest of the file content
    remaining_content = """
const unsigned int maxEvent_ = signalBool_ ? 10000 : 10000;
const std::string fileName_ =
    signalBool_        ? "mc21_14TeV_flatpT_Zprime_tthad" :
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

            std::string et_bin  = bin.substr(0, first_pipe);
            std::string eta_bin = bin.substr(first_pipe + 1, second_pipe - first_pipe - 1);
            std::string phi_bin = bin.substr(second_pipe + 1);
            //std::cout << "et_bin : " << et_bin << std::endl;
            //std::cout << "eta_bin : " << eta_bin << std::endl;
            //std::cout << "phi_bin : " << phi_bin << std::endl;

            // Prepend 5 zero bits (as MSB) to represent num_io = 0
            std::string num_io_bin = "00000";

            // Final bitstring in MSB to LSB order: num_io | et | eta | phi
            std::string full_bin = num_io_bin + et_bin + eta_bin + phi_bin;

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




def write_constants_h(constants: dict, output_file: str, unroll: int, ii: int):
    # Remove existing constants.h if it exists
    if os.path.exists(output_file):
        os.remove(output_file)

    with open(output_file, "w") as f:
        f.write("#ifndef CONSTANTS_H\n")
        f.write("#define CONSTANTS_H\n")
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
        f.write('''


const unsigned int lut_size_ = (eta_range_ * (1 << (phi_bit_length_)));
#if !WRITE_LUT
constexpr unsigned int padded_zeroes_length_ = 64 - et_bit_length_ - eta_bit_length_ - phi_bit_length_;
constexpr unsigned int padded_zeroes_length_32b_ = 128 - et_bit_length_ - eta_bit_length_ - phi_bit_length_;
constexpr unsigned int total_bits_input_ = padded_zeroes_length_32b_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
constexpr unsigned int total_bits_output_ = padded_zeroes_length_ + et_bit_length_ + eta_bit_length_ + phi_bit_length_;
typedef ap_uint<total_bits_input_> input; // need 32b input, 64b output!
typedef ap_uint<total_bits_output_> output;

constexpr unsigned int phi_low_  = 0;
constexpr unsigned int phi_high_ = phi_low_ + phi_bit_length_ - 1;

constexpr unsigned int eta_low_  = phi_high_ + 1;
constexpr unsigned int eta_high_ = eta_low_ + eta_bit_length_ - 1;

constexpr unsigned int et_low_   = eta_high_ + 1;
constexpr unsigned int et_high_  = et_low_ + et_bit_length_ - 1;

constexpr unsigned int padded_zeroes_low_  = et_high_ + 1;
constexpr unsigned int padded_zeroes_high_ = padded_zeroes_low_ + padded_zeroes_length_ - 1;


constexpr unsigned int nSeedsDeltaR_ = nSeedsInput_ - nSeedsOutput_;

constexpr unsigned int digitized_delta_R2_ = static_cast<unsigned int>(r2Cut_/deltaR2_granularity_ + 0.5);

static const bool lut_[max_R2lut_size_] =
#include "../data/LUTs/deltaR2LUT.h"
;

static const ap_uint<deltaR_bits_ > lutR_[max_Rlut_size_] = 
#include "../data/LUTs/deltaRLUT.h"
;

static const ap_uint<deltaR_bits_ > lutR_8b_[max_R_8b_lut_size_] = 
#include "../data/LUTs/deltaRLUT_8b.h"
;

#endif
constexpr unsigned int deltaR_levels_ = (1 << deltaR_bits_); // 256
constexpr float deltaR_step_ = deltaR_max_ / (deltaR_levels_ - 1); // ~0.041
constexpr unsigned int rMergeConsiderCutDigitized_ = (rMergeCut_) / deltaR_step_;
        ''')

        f.write("\n#endif // CONSTANTS_H\n")


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
    maxObjectsConsidered_options = [10]
    #maxObjectsConsidered_options = [128]
    #sortSeeds_options = [False, True]
    rMergeCut_options = [0.001]
    #rMergeCut_options = [3.5]
    sortSeeds_options = [False] # NOW ALWAYS TRUE UNDER ASSUMPTION THAT GFEX IS DOING SORTING
    signalBool_options = [True]
    #inputEnergyCut_options = [False]
    inputEnergyCut_options = [False]
    jzSlices = [3]
    #inputEnergyCutValues_options = [0.5, 1.0, 1.5, 2.0] # only use these if inputEnergyCutOption is true
    inputEnergyCutValues_options = [4] # note 4 == 1 GeV (0.25 GeV granularity for inputobjectvalues)
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
                        for sortSeeds in sortSeeds_options:
                            print("sortSeeds processed:", sortSeeds)
                            for inputEnergyCutBool in inputEnergyCut_options:
                                print("inputEnergyCutBool processed:", inputEnergyCutBool)
                                for inputEnergyCutValues in inputEnergyCutValues_options:
                                    print("inputEnergyCutValues processed:", inputEnergyCutValues)
                                    for signalBool in signalBool_options:
                                        print("signalBool processed:", signalBool)
                                        # Copy base constants
                                        constants = base_constants.copy()
                                        constants["nSeedsOutput_"] = nSeeds
                                        constants["r2Cut_"] = r2Cut
                                        constants["rCut_"] = math.sqrt(r2Cut)
                                        constants["rMergeCut_"] = rMergeCut
                                        constants["maxObjectsConsidered_"] = maxObjectsConsidered
                                        constants["sortSeeds_"] = sortSeeds
                                        constants["useInputEnergyCut_"] = inputEnergyCutBool
                                        print("nSeeds:", nSeeds)
                                        print("r2Cut:", r2Cut)
                                        print("rMergeCut:", rMergeCut)
                                        print("maxobjectsconsidered:", maxObjectsConsidered)
                                        if inputEnergyCutBool:
                                            constants["inputEnergyCut_"] = inputEnergyCutValues

                                        # Calculate phi and eta granularities from base constants
                                        phi_range = base_constants["phi_max_"] - base_constants["phi_min_"]
                                        phi_granularity = phi_range / (1 << base_constants["phi_bit_length_"])
                                        constants["phi_granularity_"] = phi_granularity

                                        eta_range = base_constants["eta_max_"] - base_constants["eta_min_"]
                                        eta_granularity = eta_range / (base_constants["eta_range_"])
                                        constants["eta_granularity_"] = eta_granularity

                                        lut_max_size = calculate_lutR2_max_size(
                                            r2Cut=constants["r2Cut_"],
                                            eta_range=base_constants["eta_range_"],
                                            phi_bit_length=constants["phi_bit_length_"],
                                            eta_granularity=eta_granularity,
                                            phi_granularity=phi_granularity
                                        )

                                        lutR_max_size = calculate_lutR_max_size(
                                            rCut=constants["rMergeCut_"],
                                            eta_range=base_constants["eta_range_"],
                                            phi_bit_length=constants["phi_bit_length_"],
                                            eta_granularity=eta_granularity,
                                            phi_granularity=phi_granularity
                                        )

                                        lutR_8b_max_size = calculate_lutR_max_size(
                                            rCut=constants["rCut_"],
                                            eta_range=base_constants["eta_range_"],
                                            phi_bit_length=constants["phi_bit_length_"],
                                            eta_granularity=eta_granularity,
                                            phi_granularity=phi_granularity
                                        )

                                        constants["max_R2lut_size_"] = lut_max_size  # Add it to constants.h!
                                        constants["max_Rlut_size_"] = lutR_max_size  # Add it to constants.h!
                                        constants["max_R_8b_lut_size_"] = lutR_8b_max_size

                                        # Update derived values
                                        constants = compute_derived(constants)

                                        # Clean up float formatting for r2Cut and inputEnergyCutValues
                                        r2Cut_str = str(r2Cut).replace('.', 'p')
                                        rMergeCut_str = str(rMergeCut).replace('.', 'p')
                                        energyCut_str = str(inputEnergyCutValues).replace('.', 'p')

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

                                        energyCutBool_str = 'ecut' if inputEnergyCutBool else 'noecut'
                                        constsFilename = f"../algorithm/constants.h"

                                        # File name based on parameter values
                                        file_suffix = (
                                            f"nSeeds{nSeeds}_"
                                            f"r2Cut{r2Cut_str}_"
                                            f"maxObj{maxObjectsConsidered}_"
                                            f"{rMergeCut_str}"
                                            f"{signal_str}"
                                            "_WTAConeJetsCellsTowers_ValidateEmulation_FINAL"
                                            #f"{energyCutBool_str}_"
                                            #f"ecutVal{energyCut_str}"
                                        )
                                        print(f"Launching HLS with project name: {file_suffix}")

                                        unroll = int(maxObjectsConsidered/8)

                                        ii = 3

                                        # Write to file
                                        write_constants_h(constants, constsFilename, unroll, ii)
                                        fileReadPath = "/home/larsonma/GEPHadronicEventReconstruction/algorithm/fileRead.h"  # Path to save the file
                                        write_file_read_header(fileReadPath, file_suffix, signalBool, jzSlice)

                                        print(f" Wrote {constsFilename}")
                                        run_lut_generator_via_root("/home/larsonma/GEPHadronicEventReconstruction/algorithm/writeDeltaR2LUT.cc")

                                        subprocess.run(["vitis", "-s", "jet_tagger_hls.py", file_suffix, "1"], check=True)
                                        xml_report_path = os.path.join('w', file_suffix, file_suffix, 'syn', 'report', 'jet_tagger_top_csynth.xml')
                                        print("xml_report_path,", xml_report_path)
                                        #resources, latency = extract_hls_report(xml_report_path)
                                        #print(f"Run {file_suffix}: Resources {resources}, Latency {latency}")
                            
