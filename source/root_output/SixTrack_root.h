#ifndef SixTrack_root_h
#define SixTrack_root_h

#include "TROOT.h"

//General Functions
extern "C" void DoSixTrackRootInit(int eos, int run_number, char* eos_server, char* root_path, char* root_prefix, int Accelerator, int Optics, int ApertureCheck, int Collimation, int CollimationDB, int FLUKA, int ApertureDump);
extern "C" void SixTrackRootWrite();
extern "C" void SixTrackRootExit();

//Dump functions
extern "C" void DumpRootInit();

#endif

