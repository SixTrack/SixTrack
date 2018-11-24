#include <iostream>
#include <string>
#include <random>

#include "TFile.h"
//#include "TTree.h"

#include "SixTrack_root.h"
#include "ApertureCheck_root.h"
#include "ApertureDump_root.h"
#include "Collimation_root.h"
#include "AcceleratorOutput_root.h"
#include "ConfigurationOutput_root.h"
#include "Optics_root.h"
#include "RunTime_root.h"

TFile *RootFile;

/**
* General configuration and start up
* i.e. opening files on EOS
* Catch errors?
*/
extern "C" void DoSixTrackRootInit(int eos, int run_number, char* eos_server, char* root_path, char* root_prefix, int Accelerator, int Optics, int ApertureCheck, int Collimation, int CollimationDB, int FLUKA, int ApertureDump)
{
    std::cout << "Root output initialization" << std::endl;
    std::cout << "Accelerator is enabled: " << Accelerator << std::endl;
    std::cout << "Optics are enabled: " << Optics << std::endl;
    std::cout << "ApertureCheck is enabled: " << ApertureCheck << std::endl;
    std::cout << "Collimation is enabled: " << Collimation << std::endl;
    std::cout << "FLUKA Collimation is enabled: " << FLUKA << std::endl;
    std::cout << "Aperture Dump is enabled: " << ApertureDump << std::endl;
    std::string fname;
    if(eos)
    {
        fname = "root://" + std::string(eos_server) + "/" + std::string(root_path) + "/" + std::string(root_prefix) + std::to_string(run_number) + ".root";
    }
    else
    {
        fname = std::string(root_path) + "/" + std::string(root_prefix) + std::to_string(run_number) + ".root";
    }

    std::cout << "Opening " << fname << std::endl;

    RootFile = TFile::Open(fname.c_str(),"RECREATE");

    //Enable maximum compression
    RootFile->SetCompressionLevel(99);

    //Dumps the configuration of this simulation to the root tree
    ConfigurationOutputRootInit();

    //Stats on the run time
    RunTimeRootInit();

    //Dumps the accelerator to the root tree
    if(Accelerator)
    {
        AcceleratorOutputRootInit();
    }

    //Dumps the optics to the root tree
    if(Optics)
    {
        OpticsRootInit();
    }

    //Dumps particles lost in aperture
    if(ApertureCheck)
    {
        ApertureCheckRootInit();
    }

    //Aperture Layout

    //Dumps particles lost in collimators and other collimation info (collimator settings)
    if(Collimation)
    {
        CollimationRootInit();
		CollimationEnergyRootInit();
    }

    if(CollimationDB)
    {
        CollimationDBRootInit();
    }

    if(FLUKA)
    {
        CollimationFLUKARootInit();
    }

    if(ApertureDump)
    {
        root_ApertureDumpInit();
    }
    //Dump

    //Tracking

    //FMA

    //Write the root file to flush the headers to storage.
    //SixTrackRootWrite();
}

/**
* Just writes the root file.
* For Checkpointing etc.
*/
extern "C" void SixTrackRootWrite()
{
    RootFile->Write();
}

/**
* Writes and closes the root file
* To be called at exit
*/
extern "C" void SixTrackRootExit()
{
    RootFile->Write();
    RootFile->Close();
}

/**
* Configuration for the Dump block
*/
extern "C" void DumpRootInit()
{

}
