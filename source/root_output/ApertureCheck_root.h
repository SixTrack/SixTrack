#ifndef ApertureCheck_root_h
#define ApertureCheck_root_h

#include "TTree.h"
#include "TROOT.h"

extern "C" void ApertureCheckRootInit();
extern "C" void ApertureCheckWriteLossParticle(int, int, int, char*, int, double, int, double, double, double, double, double, double, double, int, int);

class ApertureCheckRootOutput
{
public:

ApertureCheckRootOutput();

void WriteLossParticle(int, int, int, char*, int, double, int, double, double, double, double, double, double, double, int, int);


private:

//Aperture check functions
Int_t turn;
Int_t i;
Int_t ix;
char name[49];
Double_t slos;
Int_t ipart;
Double_t x;
Double_t xp;
Double_t y;
Double_t yp;
Double_t ct;
Double_t p;
Double_t dp;

int na;
int nz;

TTree *ApertureLossTree;

};

#endif

