#include <cmath>
#include "testclass.h"

TestClass::TestClass(int i) {
  myInt = i;
}

void TestClass::manipulateMyInt(int j){
  myInt *= j;
}

double TestClass::getPi(){
  return 4*atan2(1.0,1.0);
}
