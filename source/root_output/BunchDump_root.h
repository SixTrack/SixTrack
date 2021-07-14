#ifndef BunchDump_root_h
#define BunchDump_root_h

#include "TTree.h"
#include "TROOT.h"

extern "C" void root_BunchDumpInit();
extern "C" void root_DumpBunch(char* name_in, int name_len, int i_in, int ix_in, int turn_in, int particleID_in, int parentID_in, int pdgID_in, int16_t q_in, double weight_in, double s_in, double x_in, double xp_in, double y_in, double yp_in, double z_in, double dp_in, double sx_in, double sy_in, double sz_in, double mass_in);

class BunchDumpRootOutput
{
public:

BunchDumpRootOutput();

void DumpBunch(char* name_in, int name_len, int i_in, int ix_in, int turn_in, int particleID_in, int parentID_in, int pdgID_in, int16_t q_in, double weight_in, double s_in, double x_in, double xp_in, double y_in, double yp_in, double z_in, double dp_in, double sx_in, double sy_in, double sz_in, double mass_in);


private:

char name[49];
Int_t i;		//Current structure element (0 for StartDUMP)
Int_t ix;		//Corresponding single element (<0 for BLOC, only for ALL; 0 for StartDUMP
Int_t turn;
Int_t particleID;	//int32
Int_t parentID;		//int32
Int_t pdgID;		//int32
Short_t q;
Double_t weight;
Double_t s;
Double_t x;
Double_t xp;
Double_t y;
Double_t yp;
Double_t z;
Double_t dp;
Double_t sx;
Double_t sy;
Double_t sz;
Double_t mass;

TTree *BunchDumpTree;

};

#endif

