#include "TROOT.h"

extern "C" void CollimationRootInit();
extern "C" void CollimatorLossRootWrite(int, char*, int, int, int, double, double, double);
extern "C" void SurvivalRootWrite(int, int);
extern "C" void CollimatorDatabaseRootWrite(int, char*, int, char*, int, double, double, double, double);
