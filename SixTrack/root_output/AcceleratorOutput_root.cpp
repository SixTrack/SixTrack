#include "AcceleratorOutput_root.h"
#include "TTree.h"

Int_t ktrack;
Double_t position;
Double_t s_int;
TTree *AcceleratorTree;

/**
* Writes the Accelerator lattice to the root file
*/
extern "C" void AcceleratorOutputRootInit()
{
    //Tree stuff
    AcceleratorTree = new TTree("Accelerator","AcceleratorTree");
    AcceleratorTree->Branch("ktrack",&ktrack,"ktrack/I");
    //CollimatorLossTree->Branch("name",&name,"name/C");
    AcceleratorTree->Branch("position",&position,"position/D");
    s_int = 0;
}

extern "C" void AcceleratorRootWrite(int ktrack_in, char* name_in, double position_in)
{
    ktrack = ktrack_in;
    position = position_in;
    AcceleratorTree->Fill();
}

