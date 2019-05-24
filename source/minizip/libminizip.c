/*
  Wrapper for minizip
 ~~~~~~~~~~~~~~~~~~~~~
  V.K. Berglyd Olsen, BE-ABP-HSS
  Created: 2019-05-24

  Code extracted from minizip.c and miniunz.c
*/

#if (!defined(_WIN32)) && (!defined(WIN32)) && (!defined(__APPLE__))
#ifndef __USE_FILE_OFFSET64
#define __USE_FILE_OFFSET64
#endif
#ifndef __USE_LARGEFILE64
#define __USE_LARGEFILE64
#endif
#ifndef _LARGEFILE64_SOURCE
#define _LARGEFILE64_SOURCE
#endif
#ifndef _FILE_OFFSET_BIT
#define _FILE_OFFSET_BIT 64
#endif
#endif

#ifdef __APPLE__
// In darwin and perhaps other BSD variants off_t is a 64 bit value, hence no need for specific 64 bit functions
#define FOPEN_FUNC(filename, mode) fopen(filename, mode)
#define FTELLO_FUNC(stream) ftello(stream)
#define FSEEKO_FUNC(stream, offset, origin) fseeko(stream, offset, origin)
#else
#define FOPEN_FUNC(filename, mode) fopen64(filename, mode)
#define FTELLO_FUNC(stream) ftello64(stream)
#define FSEEKO_FUNC(stream, offset, origin) fseeko64(stream, offset, origin)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <fcntl.h>

#ifdef _WIN32
#include <direct.h>
#include <io.h>
#else
#include <unistd.h>
#include <utime.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

#include "zip.h"

#ifdef _WIN32
#define USEWIN32IOAPI
#include "iowin32.h"
#endif

#define WRITEBUFFERSIZE (16384)
#define MAXFILENAME (256)

#ifdef _WIN32
uLong filetime(f, tmzip, dt)
    char *f;        /* name of file to get info on */
    tm_zip *tmzip;  /* return value: access, modific. and creation times */
    uLong *dt;      /* dostime */
{
  int ret = 0;
  {
    FILETIME ftLocal;
    HANDLE hFind;
    WIN32_FIND_DATAA ff32;

    hFind = FindFirstFileA(f,&ff32);
    if (hFind != INVALID_HANDLE_VALUE)
    {
      FileTimeToLocalFileTime(&(ff32.ftLastWriteTime),&ftLocal);
      FileTimeToDosDateTime(&ftLocal,((LPWORD)dt)+1,((LPWORD)dt)+0);
      FindClose(hFind);
      ret = 1;
    }
  }
  return ret;
}
#else
#if defined(unix) || defined(__APPLE__)
uLong filetime(f, tmzip, dt)
    char *f;       /* name of file to get info on */
    tm_zip *tmzip; /* return value: access, modific. and creation times */
    uLong *dt;     /* dostime */
{
  int ret=0;
  struct stat s;   /* results of stat() */
  struct tm* filedate;
  time_t tm_t=0;

  if (strcmp(f,"-")!=0)
  {
    char name[MAXFILENAME+1];
    int len = strlen(f);
    if (len > MAXFILENAME)
      len = MAXFILENAME;

    strncpy(name, f,MAXFILENAME-1);
    /* strncpy doesn't append the trailing NULL if the string is too long. */
    name[ MAXFILENAME ] = '\0';

    if (name[len - 1] == '/')
      name[len - 1] = '\0';
    /* not all systems allow stat'ing a file with / appended */
    if (stat(name,&s)==0)
    {
      tm_t = s.st_mtime;
      ret = 1;
    }
  }
  filedate = localtime(&tm_t);

  tmzip->tm_sec  = filedate->tm_sec;
  tmzip->tm_min  = filedate->tm_min;
  tmzip->tm_hour = filedate->tm_hour;
  tmzip->tm_mday = filedate->tm_mday;
  tmzip->tm_mon  = filedate->tm_mon ;
  tmzip->tm_year = filedate->tm_year;

  return ret;
}
#else
uLong filetime(f, tmzip, dt)
  char *f;        /* name of file to get info on */
  tm_zip *tmzip;  /* return value: access, modific. and creation times */
  uLong *dt;      /* dostime */
{
  return 0;
}
#endif
#endif

int isLargeFile(const char* filename) {
  int i, largeFile = 0;
  ZPOS64_T pos = 0;
  FILE* pFile = FOPEN_FUNC(filename, "rb");
  if(pFile != NULL)
  {
    int n = FSEEKO_FUNC(pFile, 0, SEEK_END);
    pos = FTELLO_FUNC(pFile);
    printf("MINIZIP> File: %s ", filename);
    for(size_t i=strlen(filename); i<31; i++) printf(".");
    if(pos >= 1073741824)
      printf("  %7.2f Gb\n", filename, pos/1073741824.0);
    else if(pos >= 1048576)
      printf("  %7.2f Mb\n", filename, pos/1048576.0);
    else
      printf("  %7.2f kb\n", filename, pos/1024.0);
    if(pos >= 0xffffffff) largeFile = 1;
    fclose(pFile);
  }
  return largeFile;
}

/*
  Short version of minizip.c main()
*/

void minizip_zip(char* zipName, char* zipFiles, int nFiles, int compLevel, int* ret_err, int lenOne, int lenTwo) {

  int i, j;
  int err = 0;
  int size_buf = WRITEBUFFERSIZE;
  void* buf = NULL;

  char* archiveFile;
  archiveFile = (char*)malloc(lenOne+1);
  strncpy(archiveFile, zipName, lenOne);
  archiveFile[lenOne] = 0;

  buf = (void*)malloc(size_buf);
  if(buf == NULL) {
    printf("MINIZIP> ERROR Could not allocate buffer memory\n");
    return ZIP_INTERNALERROR;
  }

  zipFile zf;
  int errclose;
#ifdef USEWIN32IOAPI
  zlib_filefunc64_def ffunc;
  fill_win32_filefunc64A(&ffunc);
  zf = zipOpen2_64(archiveFile,0,NULL,&ffunc);
#else
  zf = zipOpen64(archiveFile,0);
#endif

  if(zf == NULL) {
    printf("MINIZIP> ERROR Could not open '%s'\n",archiveFile);
    err = ZIP_ERRNO;
  } else {
    printf("MINIZIP> Creating zip file '%s'\n",archiveFile);
  }

  for(int i=0; i<nFiles; i++) {
    char* fileNameInZip = malloc((lenTwo+1)*sizeof(char));
    strncpy(fileNameInZip,&zipFiles[lenTwo*i],lenTwo);
    fileNameInZip[lenTwo] = 0;

    // Trim trailing whitespace
    for(int j=lenTwo-1; j>=0; j--) {
      if(fileNameInZip[j] == ' ') {
        fileNameInZip[j] = 0;
      } else {
        break;
      }
    }

    FILE* fin;
    int size_read;
    char* saveFileNameInZip;
    zip_fileinfo zi;
    int zip64 = 0;

    zi.tmz_date.tm_sec  = zi.tmz_date.tm_min = zi.tmz_date.tm_hour =
    zi.tmz_date.tm_mday = zi.tmz_date.tm_mon = zi.tmz_date.tm_year = 0;
    zi.dosDate = 0;
    zi.internal_fa = 0;
    zi.external_fa = 0;
    filetime(fileNameInZip,&zi.tmz_date,&zi.dosDate);

    zip64 = isLargeFile(fileNameInZip);

    /* The path name saved, should not include a leading slash. */
    /* If it does, windows/xp and dynazip can't read the zip file. */
    saveFileNameInZip = fileNameInZip;
    while(saveFileNameInZip[0] == '\\' || saveFileNameInZip[0] == '/') {
      saveFileNameInZip++;
    }

    err = zipOpenNewFileInZip3_64(
      zf, saveFileNameInZip, &zi, NULL, 0, NULL, 0, NULL,
      (compLevel != 0) ? Z_DEFLATED : 0, compLevel,
      0, -MAX_WBITS, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY, NULL, 0, zip64
    );

    if(err != ZIP_OK) {
      printf("MINIZIP> ERROR Could not open '%s' in zip file\n",fileNameInZip);
    } else {
      fin = FOPEN_FUNC(fileNameInZip,"rb");
      if(fin == NULL) {
        err = ZIP_ERRNO;
        printf("MINIZIP> ERROR Could not open '%s' for reading\n",fileNameInZip);
      }
    }

    if(err == ZIP_OK) {
      do {
        err = ZIP_OK;
        size_read = (int)fread(buf,1,size_buf,fin);
        if(size_read < size_buf)
          if(feof(fin) == 0) {
            printf("MINIZIP> ERROR Could not read '%s'\n",fileNameInZip);
            err = ZIP_ERRNO;
          }

          if(size_read > 0) {
            err = zipWriteInFileInZip(zf,buf,size_read);
            if(err<0) {
              printf("MINIZIP> ERROR Could not write '%s' to zip file\n",fileNameInZip);
            }
          }
      } while ((err == ZIP_OK) && (size_read>0));
    }

    if(fin) {
      fclose(fin);
    }

    if(err < 0) {
      err = ZIP_ERRNO;
    } else {
      err = zipCloseFileInZip(zf);
      if(err != ZIP_OK) {
        printf("MINIZIP> ERROR Could not close '%s' in zip file\n",fileNameInZip);
      }
    }

    free(fileNameInZip);
  }

  errclose = zipClose(zf,NULL);
  if (errclose != ZIP_OK) printf("MINIZIP> ERROR Could not close '%s'\n",archiveFile);

  free(buf);
  free(archiveFile);

  *ret_err = err;

  return;
}