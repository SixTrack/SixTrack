#ifndef AcceleratorOutput_root_h
#define AcceleratorOutput_root_h

extern "C" void AcceleratorOutputRootInit();
extern "C" void AcceleratorRootWrite(char* name_in, int name_len, int ktrack_in, double value_in, double extra_in, double length_in);

#endif

