/*
*  SixTrack Wrapper for Pythia8
*  V.K. Berglyd Olsen, BE-ABP-HSS
*  Last modified: 2018-07-30
*/

#include <iostream>
#include <string>

#include "Pythia8/Pythia.h"

extern "C" void pythiaWrapper_init(int partType, int rndSeed);
extern "C" void pythiaWrapper_setBeamType(int frameType, int idA, int idB);
extern "C" void pythiaWrapper_setBeamCM(double eCM);
extern "C" void pythiaWrapper_setBeamEnergy(double eA, double eB);
extern "C" void pythiaWrapper_setBeamMomenta(double pxA, double pyA, double pzA, double pxB, double pyB, double pzB);
