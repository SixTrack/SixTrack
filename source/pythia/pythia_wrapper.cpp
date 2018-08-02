/*
*  SixTrack Wrapper for Pythia8
*  V.K. Berglyd Olsen, BE-ABP-HSS
*  Last modified: 2018-07-30
*/

#include "pythia_wrapper.h"

using namespace Pythia8;

extern "C" bool pythiaWrapper_init() {
  if(!pythia.init()) return false;
  pythia.settings.writeFile("pythia_settings.dat", true);
  std::cout << "PYTHIA> Done" << std::endl;
  return true;
}

extern "C" bool pythiaWrapper_defaults() {
  std::cout << "PYTHIA> Setting defaults" << std::endl;
  pythia.settings.flag("Init:showChangedSettings", true);
  pythia.settings.flag("SigmaTotal:mode", 3);
  pythia.settings.flag("SigmaDiffractive:mode", 3);
  return true;
}

extern "C" void pythiaWrapper_setProcess(bool sEL, bool sSD, bool sDD, bool sCD, bool sND) {
  std::cout << "PYTHIA> Setting processes" << std::endl;
  pythia.settings.flag("SoftQCD:elastic", sEL);
  pythia.settings.flag("SoftQCD:singleDiffractive", sSD);
  pythia.settings.flag("SoftQCD:doubleDiffractive", sDD);
  pythia.settings.flag("SoftQCD:centralDiffractive", sCD);
  pythia.settings.flag("SoftQCD:nonDiffractive", sND);
}

extern "C" void pythiaWrapper_setSeed(int rndSeed) {
  std::cout << "PYTHIA> Setting random seed" << std::endl;
  pythia.settings.mode("Random:seed", rndSeed);
}

extern "C" void pythiaWrapper_setBeam(int frameType, int idA, int idB, double eA, double eB) {
  std::cout << "PYTHIA> Setting beam parameters" << std::endl;
  pythia.settings.mode("Beams:frameType", frameType);
  pythia.settings.mode("Beams:idA", idA);
  pythia.settings.mode("Beams:idB", idB);
  pythia.settings.parm("Beams:eA", eA);
  pythia.settings.parm("Beams:eB", eB);
}

extern "C" void pythiaWrapper_readFile(char* fileName) {
  std::cout << "PYTHIA> Loading settings from external file" << std::endl;
  pythia.readFile(std::string(fileName));
}
