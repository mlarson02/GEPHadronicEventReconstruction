// File: generate_dat.C
// Usage in ROOT:
//   root [0] .L generate_dat.C+
//   root [1] main()

#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <bitset>

// ---- Configurable parameters ----
static const int N_EVENTS      = 100;   // how many events to generate
static const int N_SEEDS       = 6;    // packets per event for seeds.dat
static const int N_INPUTS      = 128;  // packets per event for inputs.dat
// ---------------------------------

// Helper: format 27-bit value into "bbbbbbbbbbbbb|bbbbbbbb|bbbbbb"
std::string format_bits27(uint32_t val) {
    // Keep only low 27 bits
    val &= 0x07FFFFFF;
    std::bitset<27> bits(val);
    std::string s = bits.to_string(); // 27-char string, MSB..LSB

    // Slice: 13|8|6
    return s.substr(0,13) + "|" + s.substr(13,8) + "|" + s.substr(21,6);
}

void write_file(const char* path, int nEvents, int nPackets) {
    std::ofstream out(path);
    if (!out) {
        std::cerr << "[ERROR] cannot open " << path << "\n";
        return;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dist(0, (1u<<27) - 1); // 27-bit

    for (int e = 0; e < nEvents; ++e) {
        out << "Event : " << e << "\n";
        for (int i = 0; i < nPackets; ++i) {
            uint32_t val = dist(gen);
            out << "0x" << std::setw(2) << std::setfill('0') << std::hex << i
                << " " << format_bits27(val)
                << " 0x" << std::setw(8) << std::setfill('0') << std::hex << std::uppercase << val
                << "\n";
        }
    }

    std::cout << "[OK] wrote " << nEvents << " events x " << nPackets
              << " packets to " << path << "\n";
}

int main() {
    write_file("/home/larsonma/GEPHadronicEventReconstruction/data/RandomMemPrints/seeds/random_seeds.dat",  N_EVENTS, N_SEEDS);
    write_file("/home/larsonma/GEPHadronicEventReconstruction/data/RandomMemPrints/inputObjects/random_inputObjects.dat", N_EVENTS, N_INPUTS);
    return 0;
}
