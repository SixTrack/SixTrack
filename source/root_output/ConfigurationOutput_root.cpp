#include "ConfigurationOutput_root.h"
#include "TTree.h"

Int_t npart_total;
Int_t nturns_total;
Double_t bin_size;
Double_t eref;
Double_t mref;
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
    ConfigurationTree->Branch("binsize",&bin_size,"binsize/D");
    ConfigurationTree->Branch("e0",&eref,"e0/D");
    ConfigurationTree->Branch("m0",&mref,"m0/D");
}

extern "C" void ConfigurationOutputRootSet_npart(int napx_in)
{
    npart_total = napx_in;
}

extern "C" void ConfigurationOutputRootSet_nturns(int nturns_in)
{
    nturns_total = nturns_in;
}

extern "C" void ConfigurationOutputRootSet_aperture_binsize(double bin_size_in)
{
    bin_size = bin_size_in;
}

extern "C" void ConfigurationOutputRootSet_reference_energy(double e0_in)
{
    eref = e0_in;
}

extern "C" void ConfigurationOutputRootSet_reference_mass(double nucm0_in)
{
    mref = nucm0_in;
}

extern "C" void ConfigurationRootWrite()
{
    ConfigurationTree->Fill();
}
