/*
  Wrapper for minizip
 ~~~~~~~~~~~~~~~~~~~~~
  V.K. Berglyd Olsen, BE-ABP-HSS
  Created: 2019-05-24

  Code extracted from minizip.c and miniunz.c
*/

#include "utils.h"

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
    *ret_err = ZIP_INTERNALERROR;
    return;
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
    printf("MINIZIP> Compression level is %d\n",compLevel);
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

/*
  Short version of miniunz.c main()
*/

void minizip_unzip(char* zipName, char* toDir, int* ret_err, int lenOne, int lenTwo) {

  int opt_do_extract_withoutpath = 1;
  int opt_overwrite = 1;
  const char *dirname = NULL;
  unzFile uf = NULL;

  *ret_err = 0;

#ifdef USEWIN32IOAPI
  zlib_filefunc64_def ffunc;
#endif

  char* archiveFile;
  archiveFile = (char*)malloc(lenOne+1);
  strncpy(archiveFile, zipName, lenOne);
  archiveFile[lenOne] = 0;

  char* extractDir;
  extractDir = (char*)malloc(lenTwo+1);
  strncpy(extractDir, toDir, lenTwo);
  extractDir[lenTwo] = 0;

#ifdef USEWIN32IOAPI
  fill_win32_filefunc64A(&ffunc);
  uf = unzOpen2_64(archiveFile, &ffunc);
#else
  uf = unzOpen64(archiveFile);
#endif

  if(uf == NULL) {
    printf("MINIZIP> ERROR Cannot open '%s'\n", archiveFile);
    *ret_err = 1;
    return;
  }
  printf("MINIZIP> Opened file '%s'\n",archiveFile);

#ifdef _WIN32
  if(_chdir(extractDir)) {
#else
  if(chdir(extractDir)) {
#endif
    printf("MINIZIP> ERROR Changing into '%s', aborting.\n", extractDir);
    *ret_err = 1;
    return;
  }

  *ret_err = do_extract(uf, opt_do_extract_withoutpath, opt_overwrite, NULL);

  unzClose(uf);

  return;
}
