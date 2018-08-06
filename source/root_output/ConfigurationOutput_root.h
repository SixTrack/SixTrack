#ifndef ConfigurationOutput_root_h
#define ConfigurationOutput_root_h

extern "C" void ConfigurationOutputRootInit();

extern "C" void ConfigurationOutputRootSet_npart(int);
extern "C" void ConfigurationOutputRootSet_nturns(int);
extern "C" void ConfigurationOutputRootSet_aperture_binsize(double);
extern "C" void ConfigurationOutputRootSet_reference_energy(double);
extern "C" void ConfigurationOutputRootSet_reference_mass(double);

extern "C" void ConfigurationRootWrite();

#endif

