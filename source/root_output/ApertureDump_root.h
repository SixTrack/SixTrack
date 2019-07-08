#ifndef ApertureDump_root_h
#define ApertureDump_root_h

#include "TTree.h"
#include "TROOT.h"

extern "C" void root_ApertureDumpInit();
extern "C" void root_DumpAperture(char* apname_in, int apname_len, char* aptype_in, int aptype_len, double s_in, double ap1_in, double ap2_in, double ap3_in, double ap4_in, double ap5_in, double ap6_in, double ap7_in, double ap8_in, double ap9_in, double ap10_in, double ap11_in );

class ApertureDumpRootOutput
{
public:

ApertureDumpRootOutput();

void DumpAperture(char* apname_in, int apname_len, char* aptype_in, int aptype_len, double s_in, double ap1_in, double ap2_in, double ap3_in, double ap4_in, double ap5_in, double ap6_in, double ap7_in, double ap8_in, double ap9_in, double ap10_in, double ap11_in );


private:

char name[49];
char type[3];

Double_t s;
Double_t ap1;
Double_t ap2;
Double_t ap3;
Double_t ap4;
Double_t ap5;
Double_t ap6;
Double_t ap7;
Double_t ap8;
Double_t ap9;
Double_t ap10;
Double_t ap11;

TTree *ApertureDumpTree;

};

#endif

