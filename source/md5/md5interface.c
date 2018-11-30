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

// void md5wrapper_md5Final(int ctxID, unsigned char* md5Digest) {
void md5wrapper_md5Final(int ctxID) {
  // md5Digest = malloc(32);
  unsigned char* md5Digest[32];
  MD5Final(&mdCtx[ctxID]);
  for(int i = 0; i < 16; i++) {
    printf("%02x", mdCtx[ctxID].digest[i]);
  }
  printf("\n");
  sprintf(md5Digest,"%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x",
    mdCtx[ctxID].digest[0],  mdCtx[ctxID].digest[1],
    mdCtx[ctxID].digest[2],  mdCtx[ctxID].digest[3],
    mdCtx[ctxID].digest[4],  mdCtx[ctxID].digest[5],
    mdCtx[ctxID].digest[6],  mdCtx[ctxID].digest[7],
    mdCtx[ctxID].digest[8],  mdCtx[ctxID].digest[9],
    mdCtx[ctxID].digest[10], mdCtx[ctxID].digest[11],
    mdCtx[ctxID].digest[12], mdCtx[ctxID].digest[13],
    mdCtx[ctxID].digest[14], mdCtx[ctxID].digest[15]
  );
  printf("%s\n",md5Digest);
}
