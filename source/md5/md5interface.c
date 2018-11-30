/**
 *  Fortran Interface for MD5
 */

#include <stdlib.h>
#include <stdio.h>
#include "md5.h"

void md5wrapper_digestFloatArray(double* inArr, int arrLen) {
  int   iSign;
  int   dPoint;
  char* dStr;
  printf("Size = %d\n", arrLen);
  for(int i=0; i<arrLen; i++) {
    printf("Element   %d = %24.17e\n", i, inArr[i]);
    printf("RoundCTL0 %d = %s\n", i, dtoa_cr(inArr[i],0,17,&iSign,&dPoint,&dStr));
    printf("RoundCTL1 %d = %s\n", i, dtoa_cr(inArr[i],1,17,&iSign,&dPoint,&dStr));
    printf("RoundCTL2 %d = %s\n", i, dtoa_cr(inArr[i],2,17,&iSign,&dPoint,&dStr));
    printf("RoundCTL3 %d = %s\n", i, dtoa_cr(inArr[i],3,17,&iSign,&dPoint,&dStr));
    printf("RoundCTL4 %d = %s\n", i, dtoa_cr(inArr[i],4,17,&iSign,&dPoint,&dStr));
    printf("RoundCTL5 %d = %s\n", i, dtoa_cr(inArr[i],5,17,&iSign,&dPoint,&dStr));
  }
}

MD5_CTX* mdCtx;

void md5wrapper_md5Init(int nInst) {
  mdCtx = malloc(sizeof(MD5_CTX)*nInst);
  for(int i=0; i<nInst; i++) {
    MD5Init(&mdCtx[i]);
  }
}

void md5wrapper_md5Update(int ctxID, unsigned char* inStr, unsigned int strLen) {
  MD5Update(&mdCtx[ctxID], inStr, strLen);
}

void md5wrapper_md5Final(int ctxID, int* md5Vals, int md5Size) {
  // md5Vals = malloc(md5Size*sizeof(int));
  MD5Final(&mdCtx[ctxID]);
  printf("C One>   '");
  for(int i = 0; i < 16; i++) {
    printf("%02x", mdCtx[ctxID].digest[i]);
    md5Vals[i] = (int)mdCtx[ctxID].digest[i];
  }
  printf("'\n");
}
