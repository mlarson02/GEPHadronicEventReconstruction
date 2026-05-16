#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdint>
#include <vector>
#include <cassert>
#include "../emulationHelperFunctions.h"


int main() {
    gSystem->RedirectOutput("debuglog.log", "w");

    std::vector<double > rMergeCuts;
    unsigned int etBitLength;
    unsigned int etaBitLength;
    unsigned int etaRange;
    unsigned int phiBitLength;
    unsigned int substruct0BitLength;
    unsigned int substruct1BitLength;
    unsigned int substruct2BitLength;
    unsigned int substruct3BitLength;
    unsigned int substruct4BitLength;
    unsigned int totalBits;

    

    assert(etBitLength + etaBitLength + phiBitLength + substruct0BitLength + substruct1BitLength + substruct2BitLength + substruct3BitLength + substruct4BitLength <= totalBits);
    std::vector<unsigned int > algoVersions = {2, 3};
    std::vector<double > rSquaredCuts = {1.0, 1.21, 1.44};
    //std::vector<double > rSquaredCuts = {1.21};
    std::vector<unsigned int > nIOs = {8, 32, 64, 128, 256};
    std::vector<unsigned int > nSeeds = {2};
    std::vector<unsigned int > nSeedsInput = {10};
    std::vector<unsigned int > nProtoSeeds = {6};

    
    for(unsigned int algoVersion : algoVersions){

        if(algoVersion == 2){ // BASIC ALGORITHM
            etBitLength = 13;
            etaBitLength = 10;
            etaRange = 784; // FIXME this should be computed on the fly
            phiBitLength = 9;
            substruct0BitLength = 0;
            substruct1BitLength = 0;
            substruct2BitLength = 0;
            substruct3BitLength = 0;
            substruct4BitLength = 0;
            totalBits = 64;
            rMergeCuts = {0.001}; // Only allow this to be disabled
        }

        else if(algoVersion == 3){
            etBitLength = 13;
            etaBitLength = 7;
            etaRange = 98;
            phiBitLength = 6;
            substruct0BitLength = 2;
            substruct1BitLength = 8;
            substruct2BitLength = 8;
            substruct3BitLength = 8;
            substruct4BitLength = 8;
            totalBits = 64;
            rMergeCuts = {0.001, 1.25, 1.5, 1.75, 1.875, 2.0};
            //rMergeCuts = { 2.0};
        }

        double etaGranularity = 9.8 / double(etaRange);
        double phiGranularity = 6.4 / double((1 << phiBitLength));

        for(double rMergeCut : rMergeCuts){ // Loop through algorithm configurations affecting constants.h
            for(double rSquaredCut : rSquaredCuts){ // And generate different constants files
                //unsigned int max_R2lut_size = calculate_lut_max_size(rSquaredCut, etaBitLength, phiBitLength, etaGranularity, phiGranularity, true); // fixme remove unnecessary _ as no longer global variable
                //unsigned int max_Rlut_size = calculate_lut_max_size(rMergeCut, etaBitLength, phiBitLength, etaGranularity, phiGranularity, false);
                double psi_R_granularity = 2 * sqrt(rSquaredCut) / (1 << substruct4BitLength); 
                unsigned int max_R_8b_lut_size = calculate_lut_max_size(2 * std::sqrt(rSquaredCut), etaBitLength, phiBitLength, etaGranularity, phiGranularity, false); // NOTE use this now to compute distances not just from jet axis, but between subjets up to the jet diameter distance apart
                for(unsigned int nIOOption : nIOs){
                    for(unsigned int nSeedOption : nSeeds){
                        for (unsigned int nSeedInputOption : nSeedsInput){
                            for (unsigned int nProtoSeedsOption : nProtoSeeds){
                                std::cout << "---------------------" << "\n";
                                std::cout << "writing constants file with: " << "\n";
                                std::cout << "rMergeCut: " << rMergeCut << "\n";
                                std::cout << "rSquaredCut: " << rSquaredCut << "\n";
                                std::cout << "nIOs: " << nIOOption << "\n";
                                std::cout << "nSeeds: " << nSeedOption << "\n";
                                std::cout << "---------------------" << "\n";
                                std::string constantsFilePath = makeInputConstantsFileName(rMergeCut, nIOOption, nSeedOption, rSquaredCut, algoVersion);
                                write_constants_header(constantsFilePath, rSquaredCut, rMergeCut, nIOOption, nSeedOption, nSeedInputOption, nProtoSeedsOption,
                                    totalBits, etBitLength, etaBitLength, etaRange, phiBitLength, max_R_8b_lut_size, 
                                    substruct0BitLength, substruct1BitLength, substruct2BitLength, substruct3BitLength, substruct4BitLength, algoVersion);
                            }
                        }
                    } 
                }
                
                std::string lutR_8b_output_path = makeInputLUTFileName(rMergeCut, rSquaredCut, "deltaR_8b"); // Lower precision deltaR calculation
                std::ofstream outfileR8b(lutR_8b_output_path);
                unsigned int iR = 0;
                outfileR8b << "{\n    ";

                //std::cout << "r2Cut_ for LUT: " << r2Cut_ << "\n";

                for (unsigned int etaIt = 0; etaIt < etaRange; ++etaIt) { // loop through
                    for (unsigned int phiIt = 0; phiIt < (1 << (phiBitLength - 1)); ++phiIt) { // phi bit length - 1 as max deltaPhi = pi, not 2pi
                        if (iR >= max_R_8b_lut_size) break;

                        float deltaEta = etaIt * etaGranularity;
                        float deltaPhi = phiIt * phiGranularity;
                        float deltaR2 = deltaEta * deltaEta + deltaPhi * deltaPhi;
                        float deltaR = std::sqrt(deltaR2);
                        //std::cout << "deltaR: " << deltaR << "\n";
                        // Digitize deltaR to 8-bit value
                        uint8_t digitizedDeltaR = static_cast<uint8_t>(deltaR / psi_R_granularity + 0.5f); // psi_R granularity = 2 * R_jet / (1 << 8) 8 = deltaR bit width & psi_R bit width
                        if (digitizedDeltaR > ((1 << substruct4BitLength) - 1)) digitizedDeltaR = ((1 << substruct4BitLength) - 1); // Clamp (optional safety)
                        
                        outfileR8b << static_cast<unsigned int>(digitizedDeltaR); // cast for readability

                        iR++;
                        //std::cout << "digitizedDeltaR: " << std::dec << static_cast<unsigned int>(digitizedDeltaR) << "\n";
                        if (iR < max_R_8b_lut_size) outfileR8b << ", " << ((iR % 255 == 0) ? "\n    " : "");
                    } // Loop through digitized phi values 
                    if (iR >= max_R_8b_lut_size) break;
                } // loop through digitized eta values

                outfileR8b << "\n};\n";
                outfileR8b.close();
            } // Loop through jet radii squared
        } // Loop through search distances for seed optimization 
    } // Loop through algoVersions
    

    gSystem->Exit(0);
    return 0;
} // main() function definition