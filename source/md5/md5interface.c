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
  MD5Final(&mdCtx[ctxID]);
  for(int i = 0; i < md5Size; i++) {
    if(i < 16) {
      md5Vals[i] = (int)mdCtx[ctxID].digest[i];
    } else {
      md5Vals[i] = 0;
    }
  }
}

void md5wrapper_digestString(char* inStr, int strLen, int* md5Vals, int md5Size) {
  MD5_CTX mdCtx;
  MD5Init(&mdCtx);
  MD5Update(&mdCtx, inStr, strLen);
  MD5Final(&mdCtx);
  for(int i = 0; i < md5Size; i++) {
    if(i < 16) {
      md5Vals[i] = (int)mdCtx.digest[i];
    } else {
      md5Vals[i] = 0;
    }
  }
}

void md5wrapper_digestFile(char* fileName, int* md5Vals, int md5Size) {
  MD5_CTX mdCtx;
  int bytes;
  unsigned char data[1024];

  FILE* inFile = fopen(fileName, "rb");
  if(inFile == NULL) {
    printf("MD5> ERROR Cannot open file '%s'.\n", fileName);
    md5Vals[0] = -1;
    return;
  }

  MD5Init(&mdCtx);
  while((bytes = fread(data, 1, 1024, inFile)) != 0) {
    MD5Update(&mdCtx, data, bytes);
  }
  MD5Final(&mdCtx);
  fclose (inFile);
  for(int i = 0; i < md5Size; i++) {
    if(i < 16) {
      md5Vals[i] = (int)mdCtx.digest[i];
    } else {
      md5Vals[i] = 0;
    }
  }
}
