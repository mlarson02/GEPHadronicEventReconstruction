#include "/home/larsonma/LargeRadiusJets/algorithm/jet_tag_adv.h"
#include "/home/larsonma/LargeRadiusJets/algorithm/fileRead.h"
//#include "/home/larsonma/LargeRadiusJets/algorithm/helperFunctions.h"
#include <string>
#include <iostream>
#include <array>
#include <vector>

int main() {

    input seedValues[nTotalSeeds_]; 
    input inputObjectValues[maxObjectsConsidered_];
    //std::cout << "nseeds output test bench: " << nSeedsOutput_ << "\n";

    unsigned int totalOutputJetMergedIO[nSeedsOutput_]; 
    std::vector<unsigned int > allOutputJetMergedIOLeading; 
    std::vector<unsigned int > allOutputJetMergedIOSubLeading; 

    // FIXME add further configurability for which seeds, input objects are used!
    const std::string seedFile = memPrintsPath_ + "GEPConeJetsCellsTowersSK/" + fileName_ + "_gepconejetscellstowerssk.dat";
    //const std::string inputObjectFile = memPrintsPath_ + "CaloTopo_422/" + fileName_ + "_topo422.dat";
    const std::string inputObjectFile = memPrintsPath_ + "GEPCellsTowersSK/" + fileName_ + "_gepcellstowerssk.dat";
    std::cout << "inputObjectFile: " << inputObjectFile << "\n";
    std::cout << " seed file: " << seedFile << "\n";

    // Call the function under test
    std::string outputJetsFile = memPrintsPath_ + "largeRJets_SK_ValidateEmulation/" + fileName_ + "_largeR";
    outputJetsFile += kFileSuffix; 
    outputJetsFile += ".dat";

    std::ofstream outFile(outputJetsFile);
    if (!outFile) {
        std::cerr << "Error: Could not open file " << outputJetsFile << std::endl;
        return 0;
    }
    std::cout << "digitized_d_search_squared_: " << digitized_d_search_squared_ << "\n";
    for (unsigned int iEvt = 0; iEvt < maxEvent_; iEvt++){
        outFile << "Event : " << std::dec << iEvt << std::endl;
        //std::cout << " ---------------------------------- " << "\n";
        //std::cout << "processing event: " << std::dec << iEvt << "\n";
        extract_values_from_file<nTotalSeeds_ >(seedFile, seedValues, iEvt);

        extract_values_from_file<maxObjectsConsidered_ >(inputObjectFile, inputObjectValues, iEvt);


        output outputJetValues[nSeedsOutput_];

        jet_tag_adv(seedValues, inputObjectValues, outputJetValues); // top function 
        for (unsigned int iOutput = 0; iOutput < nSeedsOutput_; iOutput++){
            unsigned long long numsubjets_value = outputJetValues[iOutput].range(num_subjets_high_, num_subjets_low_).to_uint64();
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
            std::bitset<num_subjets_length_> io_bitset(numsubjets_value); 
            std::bitset<et_bit_length_ > et_bitset(et_value);
            std::bitset<eta_bit_length_ > eta_bitset(eta_value);
            std::bitset<phi_bit_length_ > phi_bitset(phi_value);

            uint32_t combined_value =
                ((numsubjets_value & maskN(num_subjets_length_)) << (et_bit_length_ + eta_bit_length_ + phi_bit_length_)) |
                ((et_value   & maskN(et_bit_length_  )) << (eta_bit_length_ + phi_bit_length_)) |
                ((eta_value  & maskN(eta_bit_length_ )) <<  phi_bit_length_) |
                ( phi_value  & maskN(phi_bit_length_ ));

            //std::cout << "total_bits_output_: " << std::dec << total_bits_output_ << "\n";
            //std::cout << "combined_value : " << std::hex << combined_value << "\n";
            // Convert to hexadecimal (for the last field)
            std::stringstream hex_stream;
            hex_stream << std::hex << std::nouppercase << std::setfill('0') << std::setw(total_bits_output_ / 4) << combined_value;
            std::string hexValue = hex_stream.str();

            //std::cout << "hex_stream: " << std::hex << std::setfill('0') << std::setw(total_bits_output_ / 4) << hex_stream.str() << "\n";

            //std::cout << "et_bitset.to_string() : " << et_bitset.to_string() << "\n";
            // Output in the required format
            outFile << "0x" << std::setw(2) << std::setfill('0') << std::hex << iOutput
            << " " << io_bitset.to_string() << "|" << et_bitset.to_string() << "|" << eta_bitset.to_string() << "|" << phi_bitset.to_string()
            << " 0x" << hexValue << std::endl; 

            //totalOutputJetMergedIO[iOutput] += numio_value; 
            //if (iOutput == 0) allOutputJetMergedIOLeading.push_back(numio_value);
            //if (iOutput == 1) allOutputJetMergedIOSubLeading.push_back(numio_value);
        }
    }   

    std::cout << "Test bench completed." << std::endl;
    return 0;
}
