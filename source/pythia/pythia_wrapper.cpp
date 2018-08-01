/*
*  SixTrack Wrapper for Pythia8
*  V.K. Berglyd Olsen, BE-ABP-HSS
*  Last modified: 2018-07-30
*/

#include "pythia_wrapper.hpp"

using namespace Pythia8;

Pythia* pythia;

extern "C" void pythiaWrapper_init(int partType, int rndSeed) {

  pythia->settings.parm("Random:setSeed", true);
  pythia->settings.parm("Random:seed",    rndSeed);

  pythia->init();

}

extern "C" void pythiaWrapper_setBeamType(int frameType, int idA, int idB) {
  pythia->settings.parm("Beams:frameType", frameType);
  pythia->settings.parm("Beams:idA",       idA);
  pythia->settings.parm("Beams:idB",       idB);
}

extern "C" void pythiaWrapper_setBeamCM(double eCM) {
  pythia->settings.parm("Beams:eCM",       eCM);
}

extern "C" void pythiaWrapper_setBeamEnergy(double eA, double eB) {
  pythia->settings.parm("Beams:eA",        eA);
  pythia->settings.parm("Beams:eB",        eB);
}

extern "C" void pythiaWrapper_setBeamMomenta(double pxA, double pyA, double pzA, double pxB, double pyB, double pzB) {
  pythia->settings.parm("Beams:pxA",       pxA);
  pythia->settings.parm("Beams:pyA",       pyA);
  pythia->settings.parm("Beams:pzA",       pzA);
  pythia->settings.parm("Beams:pxB",       pxB);
  pythia->settings.parm("Beams:pyB",       pyB);
  pythia->settings.parm("Beams:pzB",       pzB);
}
