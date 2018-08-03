/*
*  SixTrack Wrapper for Pythia8
*  V.K. Berglyd Olsen, BE-ABP-HSS
*  Last modified: 2018-07-30
*/

#include <iostream>
#include <fstream>
#include <string>

#include "Pythia8/Pythia.h"

std::ofstream pythia_log("pythia_log.dat");
Pythia8::Pythia pythia("../share/Pythia8/xmldoc",false);

extern "C" bool pythiaWrapper_init();
extern "C" bool pythiaWrapper_defaults();
extern "C" void pythiaWrapper_setProcess(bool sEL, bool sSD, bool sDD, bool sCD, bool sND);
extern "C" void pythiaWrapper_setCoulomb(bool sCMB, double tAbsMin);
extern "C" void pythiaWrapper_setSeed(int rndSeed);
extern "C" void pythiaWrapper_setBeam(int frameType, int idA, int idB, double eA, double eB);
extern "C" void pythiaWrapper_readFile(char* fileName);
extern "C" void pythiaWrapper_getCrossSection(double& sigTot, double& sigEl);
extern "C" void pythiaWrapper_getEvent(bool& status, int& code, double& t, double& xi);
