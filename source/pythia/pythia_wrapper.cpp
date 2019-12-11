/*
*  SixTrack Wrapper for Pythia8
*  V.K. Berglyd Olsen, BE-ABP-HSS
*  Created: 2018-07-30
*  Updated: 2019-07-17
*/

#include "pythia_wrapper.h"

using namespace Pythia8;

extern "C" {

bool pythiaWrapper_init() {
  if(!pythia.init()) return false;
  pythia.settings.writeFile("pythia_settings.dat", true);
  std::cout << "PYTHIA> Done" << std::endl;
  return true;
}

bool pythiaWrapper_defaults() {
  std::cout << "PYTHIA> Setting defaults" << std::endl;
  pythia.settings.flag("Init:showChangedSettings", true);
  pythia.settings.flag("Init:showChangedParticleData", false);
  pythia.settings.flag("SigmaTotal:mode", 3);
  pythia.settings.flag("SigmaDiffractive:mode", 3);
  return true;
}

void pythiaWrapper_setProcess(bool sEL, bool sSD, bool sDD, bool sCD, bool sND) {
  std::cout << "PYTHIA> Setting processes" << std::endl;
  pythia.settings.flag("SoftQCD:elastic", sEL);
  pythia.settings.flag("SoftQCD:singleDiffractive", sSD);
  pythia.settings.flag("SoftQCD:doubleDiffractive", sDD);
  pythia.settings.flag("SoftQCD:centralDiffractive", sCD);
  pythia.settings.flag("SoftQCD:nonDiffractive", sND);
}

void pythiaWrapper_setModes(int sigTotMode, int sigDiffMode) {
  std::cout << "PYTHIA> Setting cross section modes" << std::endl;
  pythia.settings.mode("SigmaTotal:mode",       sigTotMode);
  pythia.settings.mode("SigmaDiffractive:mode", sigDiffMode);
}

void pythiaWrapper_setCrossSection(double csTot, double csEL, double csSD, double csDD, double csCD) {
  std::cout << "PYTHIA> Setting cross sections" << std::endl;
  if(csTot > 0.0) pythia.settings.parm("SigmaTotal:sigmaTot", csTot);
  if(csEL  > 0.0) pythia.settings.parm("SigmaTotal:sigmaEl",  csEL);
  if(csSD  > 0.0) pythia.settings.parm("SigmaTotal:sigmaXB",  csSD);
  if(csSD  > 0.0) pythia.settings.parm("SigmaTotal:sigmaAX",  csSD);
  if(csDD  > 0.0) pythia.settings.parm("SigmaTotal:sigmaXX",  csDD);
  if(csCD  > 0.0) pythia.settings.parm("SigmaTotal:sigmaAXB", csCD);
}

void pythiaWrapper_setCoulomb(bool sCMB, double tAbsMin) {
  pythia.settings.flag("SigmaElastic:Coulomb", sCMB);
  pythia.settings.parm("SigmaElastic:tAbsMin", tAbsMin);
}

void pythiaWrapper_setSeed(int rndSeed) {
  std::cout << "PYTHIA> Setting random seed" << std::endl;
  pythia.settings.mode("Random:seed", rndSeed);
}

void pythiaWrapper_setBeam(int frameType, int idA, int idB, double eA, double eB) {
  std::cout << "PYTHIA> Setting beam parameters" << std::endl;
  pythia.settings.mode("Beams:frameType", frameType);
  pythia.settings.mode("Beams:idA", idA);
  pythia.settings.mode("Beams:idB", idB);
  pythia.settings.parm("Beams:eA", eA);
  pythia.settings.parm("Beams:eB", eB);
  if(frameType == 3) {
    pythia.settings.flag("Beams:allowVariableEnergy", true);
  }
}

void pythiaWrapper_readFile(char* fileName) {
  std::cout << "PYTHIA> Loading settings from external file" << std::endl;
  pythia.readFile(std::string(fileName));
}

void pythiaWrapper_getInitial(double& sigTot, double& sigEl, double& m0_1, double& m0_2) {
  // sigTot = pythia.info.sigmaGen(0);
  sigTot = pythia.parm("SigmaTotal:sigmaTot");
  sigEl  = pythia.parm("SigmaTotal:sigmaEl");

  // Generate one test event, and extract the particle mass
  pythia.next();
  m0_1 = pythia.event[1].p().mCalc();
  m0_2 = pythia.event[2].p().mCalc();
}

/**
 *  SoftQCD Events
 * ================
 *  101 : Non-Diffrcative
 *  102 : Elastic             AB -> AB
 *  103 : Single Diffractive  AB -> XB
 *  104 : Single Diffractive  AB -> AX
 *  105 : Double Diffractive  AB -> XX
 *  106 : Central Diffractive AB -> AXB
 */

void pythiaWrapper_getEvent(bool& status, int& code, double& t, double& theta, double& dEE, double& dPP) {
  status = pythia.next();
  code   = pythia.info.code();
  if(!status) {
    code  = 0;
    t     = 0.0;
    theta = 0.0;
    dEE   = 0.0;
    dPP   = 0.0;
    return;
  }
  if(code == 101) {
    t     = 0.0;
    theta = 0.0;
    dEE   = 0.0;
    dPP   = 0.0;
  }
  else if(code == 102) { // Elastic
    t     =  pythia.info.tHat();
    theta =  pythia.event[3].theta();
    dEE   = (pythia.event[3].e()    - pythia.event[1].e())    / pythia.event[1].e();
    dPP   = (pythia.event[3].pAbs() - pythia.event[1].pAbs()) / pythia.event[1].pAbs();
  }
  else if(code == 104) { // Single Diffractive AB->AX
    t     =  pythia.info.tHat();
    theta =  pythia.event[3].theta();
    dEE   = (pythia.event[3].e()    - pythia.event[1].e())    / pythia.event[1].e();
    dPP   = (pythia.event[3].pAbs() - pythia.event[1].pAbs()) / pythia.event[1].pAbs();
  }
  else if(code == 106) { // Central Diffractive AB->AXB
    t     = (pythia.event[3].p()    - pythia.event[1].p()).m2Calc();
    theta =  pythia.event[3].theta();
    dEE   = (pythia.event[3].e()    - pythia.event[1].e())    / pythia.event[1].e();
    dPP   = (pythia.event[3].pAbs() - pythia.event[1].pAbs()) / pythia.event[1].pAbs();
  }
  else {
    t     = 0.0;
    theta = 0.0;
    dEE   = 0.0;
    dPP   = 0.0;
  }
}

void pythiaWrapper_getEventPVector(bool& status, int& code, double& t, double& theta, double& dEE, double& dPP, double* vecPin, double* vecPout) {

  status = pythia.next(vecPin[0],vecPin[1],vecPin[2],vecPin[3],vecPin[4],vecPin[5]);
  code   = pythia.info.code();

  t      = 0.0;
  theta  = 0.0;
  dEE    = 0.0;
  dPP    = 0.0;

  if(code == 102 || code == 104 || code == 106) {
    theta =  pythia.event[3].theta();
    dEE   = (pythia.event[3].e()    - pythia.event[1].e())    / pythia.event[1].e();
    dPP   = (pythia.event[3].pAbs() - pythia.event[1].pAbs()) / pythia.event[1].pAbs();
  }
  if(code == 102 || code == 104) {
    t = pythia.info.tHat();
  }
  if(code == 106) {
    t = (pythia.event[3].p() - pythia.event[1].p()).m2Calc();
  }

  vecPin[0]  = pythia.event[1].px();
  vecPin[1]  = pythia.event[1].py();
  vecPin[2]  = pythia.event[1].pz();
  vecPin[3]  = pythia.event[2].px();
  vecPin[4]  = pythia.event[2].py();
  vecPin[5]  = pythia.event[2].pz();
  vecPout[0] = pythia.event[3].px();
  vecPout[1] = pythia.event[3].py();
  vecPout[2] = pythia.event[3].pz();
  vecPout[3] = pythia.event[4].px();
  vecPout[4] = pythia.event[4].py();
  vecPout[5] = pythia.event[4].pz();

}

}
