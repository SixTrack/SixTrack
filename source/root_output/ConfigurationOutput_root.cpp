#include "ConfigurationOutput_root.h"
#include "TTree.h"

Int_t npart_total;
Int_t nturns_total;
//Double_t position;
TTree *ConfigurationTree;

/**
* Writes the Accelerator lattice to the root file
*/
extern "C" void ConfigurationOutputRootInit()
{
    //Tree stuff
    ConfigurationTree = new TTree("Configuration","ConfigurationTree");
    ConfigurationTree->Branch("npart",&npart_total,"npart/I");
    ConfigurationTree->Branch("nturns",&nturns_total,"nturns/I");
}

extern "C" void ConfigurationOutputRootSet_npart(int napx_in)
{
    npart_total = napx_in;
}

extern "C" void ConfigurationOutputRootSet_nturns(int nturns_in)
{
    nturns_total = nturns_in;
}

extern "C" void ConfigurationRootWrite()
{
    ConfigurationTree->Fill();
}
