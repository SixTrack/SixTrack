/*
*  SixTrack Wrapper for Pythia8
*  V.K. Berglyd Olsen, BE-ABP-HSS
*  Last modified: 2018-07-30
*/

#include <iostream>
#include <string>

#include "Pythia8/Pythia.h"

extern "C" void pythiaWrapper_init(int partType);
extern "C" void pythiaWrapper_initElastic();
extern "C" void pythiaWrapper_initDiffractive();
