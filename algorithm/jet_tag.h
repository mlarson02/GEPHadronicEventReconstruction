#ifndef _JET_TAG_H_
#define _JET_TAG_H_
#include "/home/larsonma/GEPHadronicEventReconstruction/algorithm/helperFunctions.h"
#include "/home/larsonma/GEPHadronicEventReconstruction/algorithm/constants.h"
#include <iostream>
#include <array>
#include <bitset>
#include <fstream>
#include <sstream>
#include <bitset>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <cmath>
//#include "helperFunctions.h"

// Main function
//void process_event(input seedValues[nSeedsInput_], input inputObjectValues[maxObjectsConsidered_], input (&outputJetValues)[nSeedsOutput_]); // FIXME potentially use templated / overloaded func to deal with whether write out files while running synth or c-sim
void jet_tag(input inputObjectValues[maxObjectsConsidered_], output (&outputJetValues)[nSeedsOutput_]); // FIXME potentially use templated / overloaded func to deal with whether write out files while running synth or c-sim
    
    
#endif