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
Pythia8::Pythia pythia("../share/Pythia8/xmldoc", false);

extern "C" {
  bool pythiaWrapper_init();
  bool pythiaWrapper_defaults();
  void pythiaWrapper_setProcess(bool sEL, bool sSD, bool sDD, bool sCD, bool sND);
  void pythiaWrapper_setModes(int sigTotMode, int sigDiffMode);
  void pythiaWrapper_setCrossSection(double csTot, double csEL, double csSD, double csDD, double csCD);
  void pythiaWrapper_setCoulomb(bool sCMB, double tAbsMin);
  void pythiaWrapper_setSeed(int rndSeed);
  void pythiaWrapper_setBeam(int frameType, int idA, int idB, double eA, double eB);
  void pythiaWrapper_readFile(char* fileName);
  void pythiaWrapper_getInitial(double& sigTot, double& sigEl, double& m0_1, double& m0_2);
  void pythiaWrapper_getEvent(bool& status, int& code, double& t, double& theta, double& dEE, double& dPP);
}
