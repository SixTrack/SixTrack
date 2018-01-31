#include "AcceleratorOutput_root.h"
#include "TTree.h"
#include <algorithm>

char ElementName[17];
Int_t type;
Double_t value;
Double_t extra;
Double_t ElementLength;
Double_t position;

TTree *AcceleratorTree;

/**
* Writes the Accelerator lattice to the root file
* read fort.2 (or fort.3), idat -> bez = single element name,
* kz = type of element,
* ed,ek,el = strength, random error on strenght, length (can be anything)
* bbbx,bbby,bbbs = beam-beam, beam-beam parameters (will be removed soon)
* read(ch1,*) idat,kz(i),ed(i),ek(i),el(i),bbbx(i),bbby(i),bbbs(i)
*/
extern "C" void AcceleratorOutputRootInit()
{
    //Tree stuff
    AcceleratorTree = new TTree("Accelerator","AcceleratorTree");
    AcceleratorTree->Branch("name",&ElementName,"name[17]/C");
    AcceleratorTree->Branch("type",&type,"type/I");
    AcceleratorTree->Branch("value",&value,"value/D");
    AcceleratorTree->Branch("extra",&extra,"extra/D");
    AcceleratorTree->Branch("length",&ElementLength,"length/D");
    AcceleratorTree->Branch("position",&position,"position/D");
    position = 0;
}

extern "C" void AcceleratorRootWrite(char* name_in, int name_len, int ktrack_in, double value_in, double extra_in, double length_in)
{
    std::string tmpname(name_in);
    std::transform(tmpname.begin(), tmpname.end(), tmpname.begin(), ::toupper);
    strncpy(ElementName,tmpname.substr(0,name_len).c_str(),17);

    type = ktrack_in;
    value = value_in;
    extra = extra_in;

    if(std::abs(type) == 12)
    {
        ElementLength = 0.0;
    }
    else
    {
        ElementLength = length_in;
    }

    position += ElementLength;

    AcceleratorTree->Fill();
}

