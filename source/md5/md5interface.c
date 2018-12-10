/**
 *  Fortran Interface for MD5
 */

#include <stdlib.h>
#include <stdio.h>
#include "md5.h"

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

void md5wrapper_digestFile(char* fileName, int strLen, int* md5Vals, int md5Size, int isAscii) {
  MD5_CTX mdCtx;
  int bytes;
  unsigned char data[1024];

  char* fileNameMod;
  fileNameMod = (char*) malloc(strLen+4);

#ifdef WIN32
  unsigned char dataTmp[1024];

  if (isAscii == 1) { //STRIP THE CARRIAGE RETURNS FROM THE FILE.
      snprintf(fileNameMod,strLen+4,"%s.tmp",fileName);

      FILE* inFile = fopen(fileName, "rb");
      FILE* outFile = fopen(fileNameMod, "wb");
      while((bytes = fread(data, 1, 1024, inFile)) != 0) {
          int j = 0;
          for (int i = 0; i < bytes; i++) {
              if (data[i] != '\r') {
                  dataTmp[j]=data[i];
                  j++;
              }
          }
          int err = fwrite(data, sizeof(char), j, outFile);
          if (err) {
              printf("MD5> ERROR in stripping file '%s'->'%s'.\n", fileName,fileNameMod);
          }
      }
      fclose (inFile);
      fclose (outFile);
  }
  else {
      snprintf(fileNameMod,strLen,"%s",fileName);
  }
#else
  snprintf(fileNameMod,strLen,"%s",fileName);
#endif

  FILE* inFile = fopen(fileNameMod, "rb");
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

  free (fileNameMod);
}
