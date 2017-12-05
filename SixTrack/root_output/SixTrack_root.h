#include "TROOT.h"

//General Functions
extern "C" void DoSixTrackRootInit(int eos, int run_number, char* eos_server, char* root_path, char* root_prefix, int ApertureCheck, int Collimation);
extern "C" void SixTrackRootWrite();
extern "C" void SixTrackRootExit();
//Collimation functions
extern "C" void CollimationRootInit();
//Dump functions
extern "C" void DumpRootInit();
