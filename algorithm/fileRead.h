#ifndef FILE_READ_H  // Check if the macro is defined
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
const std::string memPrintsPath_ = "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints_v3/";
static inline uint32_t maskN(unsigned n) { return (n >= 32) ? 0xFFFFFFFFu : ((1u << n) - 1u); }
const std::string kFileSuffix = "nSeeds2_r2Cut1p21_maxObj128_rMerge2p0_sig_WTAConeJetsCellsTowers_Adv_ValidateEmulation_FINAL";
constexpr bool signalBool_ = true;
constexpr unsigned int jzSlice_ = 3;

const unsigned int maxEvent_ = signalBool_ ? 10000 : 10000;
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
