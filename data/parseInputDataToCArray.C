// File: parseInputDataToCArray.C
// Usage in ROOT:
//   root [0] .L parseInputDataToCArray.C+   // or .L parseInputDataToCArray.C
//   root [1] main();             // no args

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdint>
#include <cstdlib>

// ---- Edit these defaults if you want different paths/limits ----
static const char* OUTPUT_PATH = "CArrayData/input_data.h";
static const int   N_INPUT_OBJECTS     = 128;    // take first N seeds per event
static const int   MAX_EVENTS  = 100;  // cap number of events to process
// ----------------------------------------------------------------

// Trim whitespace
static inline std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, b - a + 1);
}

// Find substring starting at the last "0x" or "0X" occurrence
static inline std::string find_last_hex_token(const std::string& line) {
    for (int i = static_cast<int>(line.size()) - 2; i >= 0; --i) {
        if (line[i] == '0' && (line[i+1] == 'x' || line[i+1] == 'X')) {
            return line.substr(i);
        }
    }
    return "";
}

int main() {
    const bool physicsOrRandom = false;
    const char* INPUT_PATH  = (physicsOrRandom) ? "/home/larsonma/GEPHadronicEventReconstruction/data/MemPrints/GEPBasicClusters/mc21_14TeV_hh_bbbb_vbf_novhh_gepbasicclusters.dat" : "/home/larsonma/GEPHadronicEventReconstruction/data/RandomMemPrints/inputObjects/random_inputObjects.dat";
    std::ifstream in(INPUT_PATH);
    if (!in) {
        std::cerr << "[ERROR] Cannot open input: " << INPUT_PATH << "\n";
        return 1;
    }

    std::vector<std::vector<uint32_t>> events;
    std::vector<uint32_t> current;
    std::string line;
    bool inEvent = false;

    while (std::getline(in, line)) {
        std::string t = trim(line);
        if (t.empty()) continue;

        if (t.rfind("Event", 0) == 0) {
            if (static_cast<int>(events.size()) >= MAX_EVENTS) break; // cap
            if (inEvent) {
                events.push_back(current);
                current.clear();
            }
            inEvent = true;
            continue;
        }

        if (!inEvent) continue; // ignore preamble

        std::string lastHex = find_last_hex_token(t);
        if (!lastHex.empty()) {
            std::istringstream iss(lastHex);
            std::string tok;
            iss >> tok; // e.g. "0x00ffe554"
            try {
                uint32_t val = static_cast<uint32_t>(std::stoul(tok, nullptr, 16));
                current.push_back(val);
            } catch (...) {
                std::cerr << "[WARN] Failed to parse token: " << tok << "\n";
            }
        }
    }
    if (inEvent && static_cast<int>(events.size()) < MAX_EVENTS) {
        events.push_back(current);
    }

    if (events.empty()) {
        std::cerr << "[ERROR] No events parsed. Check format/path.\n";
        return 1;
    }

    // Validate each event has at least N_INPUT_OBJECTS items
    for (size_t e = 0; e < events.size(); ++e) {
        if (static_cast<int>(events[e].size()) < N_INPUT_OBJECTS) {
            std::cerr << "[ERROR] Event " << e << " has only "
                      << events[e].size() << " items; need at least "
                      << N_INPUT_OBJECTS << ".\n";
            return 1;
        }
    }

    std::ofstream out(OUTPUT_PATH);
    if (!out) {
        std::cerr << "[ERROR] Cannot open output: " << OUTPUT_PATH << "\n";
        return 1;
    }

    out << "#ifndef INPUT_DATA_H\n#define INPUT_DATA_H\n\n";
    out << "#define INPUT_N_EVENTS " << events.size() << "\n";
    out << "#define N_INPUT_OBJECTS  " << N_INPUT_OBJECTS << "\n\n";
    out << "static const unsigned int inputs[INPUT_N_EVENTS][N_INPUT_OBJECTS] = {\n";

    out << std::uppercase << std::hex << std::setfill('0');
    for (size_t e = 0; e < events.size(); ++e) {
        out << "  {";
        for (int i = 0; i < N_INPUT_OBJECTS; ++i) {
            out << "0x" << std::setw(8) << events[e][i];
            if (i + 1 < N_INPUT_OBJECTS) out << ", ";
        }
        out << "}";
        if (e + 1 < events.size()) out << ",";
        out << "\n";
    }
    out << "};\n\n#endif // SEEDS_DATA_H\n";

    out.close();

    std::cout << "[OK] Parsed " << events.size() << " events (limit " << MAX_EVENTS
              << "). Wrote first " << N_INPUT_OBJECTS << " seeds/event to "
              << OUTPUT_PATH << ".\n";
    std::cout << "    INPUT=  " << INPUT_PATH << "\n";

    return 0;
}
