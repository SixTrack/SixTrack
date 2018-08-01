/*
*  SixTrack Wrapper for Pythia8
*  V.K. Berglyd Olsen, BE-ABP-HSS
*  Last modified: 2018-07-30
*/

#include "pythia_wrapper.hpp"

using namespace Pythia8;

Pythia pythia;

extern "C" void pythiaWrapper_init() {
  std::cout << "PYTHIA> Hello Kitty!" << std::endl;
  pythia.init();
}

extern "C" void pythiaWrapper_defaults() {
  pythia.settings.writeFile("pythia_log.dat");
  pythia.settings.flag("Init:showChangedSettings", true);
  pythia.settings.flag("SoftQCD:all", false);
  pythia.settings.flag("HardQCD:all", false);
}

extern "C" void pythiaWrapper_setSeed(int rndSeed) {
  pythia.settings.mode("Random:seed", rndSeed);
}

extern "C" void pythiaWrapper_setBeam(int frameType, int idA, int idB, double eA, double eB) {
  pythia.settings.mode("Beams:frameType", frameType);
  pythia.settings.mode("Beams:idA", idA);
  pythia.settings.mode("Beams:idB", idB);
  pythia.settings.parm("Beams:eA", eA);
  pythia.settings.parm("Beams:eB", eB);
}

extern "C" void pythiaWrapper_readFile(char* fileName) {
  pythia.readFile(std::string(fileName));
}
