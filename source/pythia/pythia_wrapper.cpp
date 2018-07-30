/*
*  SixTrack Wrapper for Pythia8
*  V.K. Berglyd Olsen, BE-ABP-HSS
*  Last modified: 2018-07-30
*/

#include "pythia_wrapper.hpp"

using namespace Pythia8;

Pythia* pythia;

extern "C" void pythiaWrapper_init(int partType) {
  pythia->init();
}
