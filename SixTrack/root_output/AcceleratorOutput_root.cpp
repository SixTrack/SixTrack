#include "AcceleratorOutput_root.h"
#include "TTree.h"

Int_t ktrack;
Double_t position;


TTree *AcceleratorTree;

/**
*/
extern "C" void AcceleratorOutputRootInit()
{

//Tree stuff
AcceleratorTree = new TTree("Accelerator","AcceleratorTree");
AcceleratorTree->Branch("ktrack",&ktrack,"ktrack/I");
//CollimatorLossTree->Branch("name",&name,"name/C");
AcceleratorTree->Branch("position",&position,"position/D");
}

extern "C" void AcceleratorRootWrite(int ktrack_in, char* name_in, double position_in)
{
ktrack = ktrack_in;
position = position_in;
AcceleratorTree->Fill();

}

