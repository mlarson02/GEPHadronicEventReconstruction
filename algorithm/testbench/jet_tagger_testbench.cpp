#include "/home/larsonma/LargeRadiusJets/algorithm/jet_tag.h"
#include "/home/larsonma/LargeRadiusJets/algorithm/fileRead.h"
//#include "/home/larsonma/LargeRadiusJets/algorithm/helperFunctions.h"
#include <string>
#include <iostream>
#include <array>
#include <vector>

int main() {

    input inputObjectValues[maxObjectsConsidered_];

    const std::string inputObjectFile = memPrintsPath_ + "GEPConeJetsCellsTowersSK/" + fileName_ + "_gepconejetscellstowerssk.dat";
    std::cout << "inputObjectFile: " << inputObjectFile << "\n";
    // Call the function under test
    std::string outputJetsFile = memPrintsPath_ + "largeRJets_SK_ValidateEmulation/" + fileName_ + "_largeR";
    
    outputJetsFile += kFileSuffix; 
    outputJetsFile += ".dat";
    std::cout << "output file : "<< outputJetsFile << "\n";
    std::ofstream outFile(outputJetsFile);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << outputJetsFile << std::endl;
        return 0;
    }
    std::cout << "rMergeConsiderCutDigitized_: " << rMergeConsiderCutDigitized_ << "\n";
    for (unsigned int iEvt = 0; iEvt < maxEvent_; iEvt++){
        std::cout << "----------------------------------" << "\n";
        std::cout << "iEvt: " << iEvt << "\n";
        outFile << "Event : " << std::dec << iEvt << std::endl;

        extract_values_from_file<maxObjectsConsidered_ >(inputObjectFile, inputObjectValues, iEvt);

        for (unsigned j = 0; j < maxObjectsConsidered_; ++j) {
            uint32_t w = static_cast<uint32_t>(inputObjectValues[j]);
            std::cout << "input object value: 0x"
                    << std::setw(8) << std::setfill('0') << std::uppercase
                    << std::hex << w
                    << std::dec << std::setfill(' ') << "\n";
        }
        output outputJetValues[nSeedsOutput_];

        jet_tag(inputObjectValues, outputJetValues); // top function 
        for (unsigned int iOutput = 0; iOutput < nSeedsOutput_; iOutput++){
            unsigned long long et_value = outputJetValues[iOutput].range(et_high_, et_low_).to_uint64();  // Convert to unsigned long long
            unsigned long long eta_value = outputJetValues[iOutput].range(eta_high_, eta_low_).to_uint64();
            unsigned long long phi_value = outputJetValues[iOutput].range(phi_high_, phi_low_).to_uint64();
            /*
            std::cout << "OUTPUT DATA: " << "\n";
            std::cout << "Et: " << outputJetValues[iOutput].range(et_high_, et_low_) << "\n";
            std::cout << "Et 64: " << outputJetValues[iOutput].range(et_high_, et_low_).to_uint64() << "\n";
            std::cout << "eta: " << outputJetValues[iOutput].range(eta_high_, eta_low_) << "\n";
            std::cout << "eta 64: " << eta_value << "\n";*/
            //std::cout << "et_value : " << et_value << "\n";
     
            // Convert to std::bitset
            std::bitset<et_bit_length_ > et_bitset(et_value);
            std::bitset<eta_bit_length_ > eta_bitset(eta_value);
            std::bitset<phi_bit_length_ > phi_bitset(phi_value);

            uint64_t combined_value =
                ((et_value   & maskN(et_bit_length_  )) << (eta_bit_length_ + phi_bit_length_)) |
                ((eta_value  & maskN(eta_bit_length_ )) <<  phi_bit_length_) |
                ( phi_value  & maskN(phi_bit_length_ ));

            // Convert to hexadecimal (for the last field)
            std::stringstream hex_stream;
            hex_stream << std::hex << std::nouppercase << std::setfill('0') << std::setw(total_bits_output_ / 4) << combined_value;
            std::string hexValue = hex_stream.str();

            //std::cout << "hex_stream: " << std::hex << std::setfill('0') << std::setw(total_bits_output_ / 4) << hex_stream.str() << "\n";

            //std::cout << "et_bitset.to_string() : " << et_bitset.to_string() << "\n";
            // Output in the required format
            outFile << "0x" << std::setw(2) << std::setfill('0') << std::hex << iOutput
            << " " << et_bitset.to_string() << "|" << eta_bitset.to_string() << "|" << phi_bitset.to_string()
            << " 0x" << hexValue << std::endl; 

            //totalOutputJetMergedIO[iOutput] += numio_value; 
            //if (iOutput == 0) allOutputJetMergedIOLeading.push_back(numio_value);
            //if (iOutput == 1) allOutputJetMergedIOSubLeading.push_back(numio_value);
        }
    }   

    std::cout << "Test bench completed." << std::endl;
    return 0;
}
