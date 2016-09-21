#include <iostream>
//#include "testclass.h"

using namespace std;

extern "C" void merlinscatter_setup_(double* Plab){
  cout << "merlinscatter_setup, plab= "<< Plab << endl;
}

extern "C" void merlinscatter_setdata_(double* pptot, double* ppel, double* ppsd){
  cout << "merlinscatter_setdata" << endl;
  *pptot = 10.0;
  *ppel  = 2.0;
  *ppsd  = 3.0;
}
