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

int check_exist_file(filename)
  const char* filename;
{
  FILE* ftestexist;
  int ret = 1;
  ftestexist = FOPEN_FUNC(filename,"rb");
  if (ftestexist==NULL)
    ret = 0;
  else
    fclose(ftestexist);
  return ret;
}

/* calculate the CRC32 of a file,
   because to encrypt a file, we need known the CRC32 of the file before */
int getFileCrc(const char* filenameinzip,void*buf,unsigned long size_buf,unsigned long* result_crc)
{
  unsigned long calculate_crc=0;
  int err=ZIP_OK;
  FILE * fin = FOPEN_FUNC(filenameinzip,"rb");

  unsigned long size_read = 0;
  unsigned long total_read = 0;
  if (fin==NULL)
  {
    err = ZIP_ERRNO;
  }

  if (err == ZIP_OK)
    do
    {
      err = ZIP_OK;
      size_read = (int)fread(buf,1,size_buf,fin);
      if (size_read < size_buf)
        if (feof(fin)==0)
        {
          printf("error in reading %s\n",filenameinzip);
          err = ZIP_ERRNO;
        }

        if (size_read>0)
          calculate_crc = crc32(calculate_crc,buf,size_read);
        total_read += size_read;

    } while ((err == ZIP_OK) && (size_read>0));

    if (fin)
      fclose(fin);

  *result_crc=calculate_crc;
  printf("file %s crc %lx\n", filenameinzip, calculate_crc);
  return err;
}

int isLargeFile(const char* filename)
{
  int largeFile = 0;
  ZPOS64_T pos = 0;
  FILE* pFile = FOPEN_FUNC(filename, "rb");

  if(pFile != NULL)
  {
    int n = FSEEKO_FUNC(pFile, 0, SEEK_END);
    pos = FTELLO_FUNC(pFile);

    printf("File : %s is %lld bytes\n", filename, pos);

    if(pos >= 0xffffffff)
      largeFile = 1;

    fclose(pFile);
  }

  return largeFile;
}

/*
  Short version of minizip.c main()
*/

int minizip_dozip(
  char* zipfile,
  char* filetozip,
  int opt_overwrite,
  int opt_compress_level,
  int len_zipfile,
  int len_filetozip
) {
  int i;
  char filename_try[MAXFILENAME+16];
  int zipok;
  int err = 0;
  int opt_exclude_path = 0;
  int size_buf = 0;
  void* buf = NULL;

  size_buf = WRITEBUFFERSIZE;
  buf = (void*)malloc(size_buf);
  if(buf == NULL) {
    printf("MINIZIP> ERROR Could not allocate buffer memory\n");
    return ZIP_INTERNALERROR;
  }

  zipok = 1 ;
  strncpy(filename_try, zipfile, MAXFILENAME-1);
  /* strncpy doesnt append the trailing NULL if the string is too long. */
  filename_try[ MAXFILENAME ] = '\0';

  if(opt_overwrite==2) {
    /* if the file don't exist, we not append file */
    if (check_exist_file(filename_try)==0)
      opt_overwrite=1;
  } else if (opt_overwrite==0) {
    if (check_exist_file(filename_try)!=0) opt_overwrite = 2;
  }

  zipFile zf;
  int errclose;
#ifdef USEWIN32IOAPI
  zlib_filefunc64_def ffunc;
  fill_win32_filefunc64A(&ffunc);
  zf = zipOpen2_64(filename_try,(opt_overwrite==2) ? 2 : 0,NULL,&ffunc);
#else
  zf = zipOpen64(filename_try,(opt_overwrite==2) ? 2 : 0);
#endif

  if(zf == NULL) {
    printf("error opening %s\n",filename_try);
    err = ZIP_ERRNO;
  } else {
    printf("creating %s\n",filename_try);
  }

  // for(i=zipfilenamearg+1;(i<argc) && (err==ZIP_OK);i++) {

    FILE* fin;
    int size_read;
    const char* filenameinzip = filetozip;
    const char* savefilenameinzip;
    zip_fileinfo zi;
    int zip64 = 0;

    zi.tmz_date.tm_sec = zi.tmz_date.tm_min = zi.tmz_date.tm_hour =
    zi.tmz_date.tm_mday = zi.tmz_date.tm_mon = zi.tmz_date.tm_year = 0;
    zi.dosDate = 0;
    zi.internal_fa = 0;
    zi.external_fa = 0;
    filetime(filenameinzip,&zi.tmz_date,&zi.dosDate);

    zip64 = isLargeFile(filenameinzip);

    //  /* The path name saved, should not include a leading slash. */
    //  /*if it did, windows/xp and dynazip couldn't read the zip file. */
    //    savefilenameinzip = filenameinzip;
    //    while( savefilenameinzip[0] == '\\' || savefilenameinzip[0] == '/' )
    //    {
    //        savefilenameinzip++;
    //    }

    /*should the zip file contain any path at all?*/
    // if(opt_exclude_path) {
    //   const char* tmpptr;
    //   const char* lastslash = 0;
    //   for(tmpptr = savefilenameinzip; *tmpptr; tmpptr++) {
    //     if(*tmpptr == '\\' || *tmpptr == '/') {
    //       lastslash = tmpptr;
    //     }
    //   }
    //   if(lastslash != NULL) {
    //     savefilenameinzip = lastslash+1; // base filename follows last slash.
    //   }
    // }

    err = zipOpenNewFileInZip3_64(
      zf,filenameinzip,&zi,
      NULL,0,NULL,0,NULL,
      (opt_compress_level != 0) ? Z_DEFLATED : 0,
      opt_compress_level,0,
      -MAX_WBITS, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY,
      NULL, 0, zip64
    );

    if(err != ZIP_OK) {
      printf("error in opening %s in zipfile\n",filenameinzip);
    } else {
      fin = FOPEN_FUNC(filenameinzip,"rb");
      if(fin==NULL) {
        err=ZIP_ERRNO;
        printf("error in opening %s for reading\n",filenameinzip);
      }
    }

    if(err == ZIP_OK) {
      do {
        err = ZIP_OK;
        size_read = (int)fread(buf,1,size_buf,fin);
        if(size_read < size_buf)
          if(feof(fin)==0) {
            printf("error in reading %s\n",filenameinzip);
            err = ZIP_ERRNO;
          }

          if(size_read>0) {
            err = zipWriteInFileInZip (zf,buf,size_read);
            if(err<0) {
              printf("error in writing %s in the zipfile\n",filenameinzip);
            }
          }
      } while ((err == ZIP_OK) && (size_read>0));
    }

    if(fin) fclose(fin);

    if(err<0) {
      err = ZIP_ERRNO;
    } else {
      err = zipCloseFileInZip(zf);
      if(err!=ZIP_OK) {
        printf("error in closing %s in the zipfile\n",filenameinzip);
      }
    }
  // }

  errclose = zipClose(zf,NULL);
  if (errclose != ZIP_OK) printf("error in closing %s\n",filename_try);

  free(buf);
  return 0;

}