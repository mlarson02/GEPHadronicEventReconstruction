#ifndef _JET_TAG_ADV_H_
#define _JET_TAG_ADV_H_
#include "/home/larsonma/LargeRadiusJets/algorithm/helperFunctions_adv.h"
#include "/home/larsonma/LargeRadiusJets/algorithm/constants_adv.h"
#include <iostream>
#include <array>
#include <bitset>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <cmath>

// Main function
//void process_event(input seedValues[nSeedsInput_], input inputObjectValues[maxObjectsConsidered_], input (&outputJetValues)[nSeedsOutput_]); // FIXME potentially use templated / overloaded func to deal with whether write out files while running synth or c-sim
void jet_tag_adv(input seedValues[nTotalSeeds_], input inputObjectValues[maxObjectsConsidered_], output (&outputJetValues)[nSeedsOutput_]); // FIXME potentially use templated / overloaded func to deal with whether write out files while running synth or c-sim


#endif
