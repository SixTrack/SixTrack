/*
*  SixTrack Wrapper for Pythia8
*  V.K. Berglyd Olsen, BE-ABP-HSS
*  Last modified: 2018-07-30
*/

#include <iostream>
#include <string>

#include "Pythia8/Pythia.h"

extern "C" void pythiaWrapper_init();
extern "C" void pythiaWrapper_defaults();
extern "C" void pythiaWrapper_setSeed(int rndSeed);
extern "C" void pythiaWrapper_setBeam(int frameType, int idA, int idB, double eA, double eB);
extern "C" void pythiaWrapper_readFile(char* fileName);
