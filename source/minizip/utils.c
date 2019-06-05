/*
  Wrapper for minizip
 ~~~~~~~~~~~~~~~~~~~~~
  V.K. Berglyd Olsen, BE-ABP-HSS
  Created: 2019-05-24

  Code extracted from minizip.c and miniunz.c
*/

#include "utils.h"

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

/* change_file_date : change the date/time of a file
    filename : the filename of the file where date/time must be modified
    dosdate : the new date at the MSDos format (4 bytes)
    tmu_date : the SAME new date at the tm_unz format */
void change_file_date(filename,dosdate,tmu_date)
  const char *filename;
  uLong dosdate;
  tm_unz tmu_date;
{
#ifdef _WIN32
  HANDLE hFile;
  FILETIME ftm,ftLocal,ftCreate,ftLastAcc,ftLastWrite;

  hFile = CreateFileA(filename,GENERIC_READ | GENERIC_WRITE, 0,NULL,OPEN_EXISTING,0,NULL);
  GetFileTime(hFile,&ftCreate,&ftLastAcc,&ftLastWrite);
  DosDateTimeToFileTime((WORD)(dosdate>>16),(WORD)dosdate,&ftLocal);
  LocalFileTimeToFileTime(&ftLocal,&ftm);
  SetFileTime(hFile,&ftm,&ftLastAcc,&ftm);
  CloseHandle(hFile);
#else
#if defined(unix) || defined(__APPLE__)
  struct utimbuf ut;
  struct tm newdate;
  newdate.tm_sec = tmu_date.tm_sec;
  newdate.tm_min=tmu_date.tm_min;
  newdate.tm_hour=tmu_date.tm_hour;
  newdate.tm_mday=tmu_date.tm_mday;
  newdate.tm_mon=tmu_date.tm_mon;
  if (tmu_date.tm_year > 1900)
    newdate.tm_year=tmu_date.tm_year - 1900;
  else
    newdate.tm_year=tmu_date.tm_year ;
  newdate.tm_isdst=-1;

  ut.actime=ut.modtime=mktime(&newdate);
  utime(filename,&ut);
#endif
#endif
}

int isLargeFile(const char* filename) {
  int i, largeFile = 0;
  ZPOS64_T pos = 0;
  FILE* pFile = FOPEN_FUNC(filename, "rb");
  if(pFile != NULL)
  {
    int n = FSEEKO_FUNC(pFile, 0, SEEK_END);
    pos = FTELLO_FUNC(pFile);
    double dpos = (double)pos;
    printf("MINIZIP> File: %s ", filename);
    for(size_t i=strlen(filename); i<31; i++) printf(".");
    if(pos >= 1073741824)
      printf("  %7.2f Gb\n", dpos/1073741824.0);
    else if(pos >= 1048576)
      printf("  %7.2f Mb\n", dpos/1048576.0);
    else
      printf("  %7.2f kb\n", dpos/1024.0);
    if(pos >= 0xffffffff) largeFile = 1;
    fclose(pFile);
  }
  return largeFile;
}

int do_extract_currentfile(uf,popt_extract_without_path,popt_overwrite,password)
  unzFile uf;
  const int* popt_extract_without_path;
  int* popt_overwrite;
  const char* password;
{
  char filename_inzip[256];
  char* filename_withoutpath;
  char* p;
  int err=UNZ_OK;
  FILE *fout=NULL;
  void* buf;
  uInt size_buf;

  unz_file_info64 file_info;
  uLong ratio=0;
  err = unzGetCurrentFileInfo64(uf,&file_info,filename_inzip,sizeof(filename_inzip),NULL,0,NULL,0);

  if(err!=UNZ_OK) {
    printf("MINIZIP> ERROR unzGetCurrentFileInfo reported error %d\n",err);
    return err;
  }

  size_buf = WRITEBUFFERSIZE;
  buf = (void*)malloc(size_buf);
  if(buf==NULL) {
    printf("MINIZIP> ERROR Allocating memory\n");
    return UNZ_INTERNALERROR;
  }

  p = filename_withoutpath = filename_inzip;
  while((*p) != '\0') {
    if(((*p)=='/') || ((*p)=='\\'))
      filename_withoutpath = p+1;
    p++;
  }

  if((*filename_withoutpath) == '\0') {

    if((*popt_extract_without_path) == 0) {
      printf("MINIZIP> Creating directory '%s'\n",filename_inzip);
      mz_mkdir(filename_inzip);
    }

  } else {

    const char* write_filename;
    int skip=0;

    if((*popt_extract_without_path)==0)
      write_filename = filename_inzip;
    else
      write_filename = filename_withoutpath;

    err = unzOpenCurrentFilePassword(uf,password);
    if(err!=UNZ_OK) {
      printf("MINIZIP> ERROR unzOpenCurrentFilePassword reported error %d\n",err);
    }

    if(((*popt_overwrite)==0) && (err==UNZ_OK)) {
      char rep=0;
      FILE* ftestexist;
      ftestexist = FOPEN_FUNC(write_filename,"rb");
      if(ftestexist != NULL) {
        fclose(ftestexist);
        if(*popt_overwrite == 0) {
          printf("MINIZIP> ERROR File already exists '%s'\n",write_filename);
          return 1;
        }
      }
    }

    if((skip==0) && (err==UNZ_OK)) {
      fout = FOPEN_FUNC(write_filename,"wb");
      /* some zipfile don't contain directory alone before file */
      if((fout==NULL) && ((*popt_extract_without_path)==0) &&
                          (filename_withoutpath!=(char*)filename_inzip)) {
        char c=*(filename_withoutpath-1);
        *(filename_withoutpath-1)='\0';
        makedir(write_filename);
        *(filename_withoutpath-1)=c;
        fout=FOPEN_FUNC(write_filename,"wb");
      }

      if(fout == NULL) {
        printf("MINIZIP> ERROR Opening '%s'\n",write_filename);
      }
    }

    if(fout!=NULL) {
      printf("MINIZIP> Extracting: '%s'\n",write_filename);
      do {
        err = unzReadCurrentFile(uf,buf,size_buf);
        if(err < 0) {
          printf("MINIZIP> ERROR unzReadCurrentFile reported error %d\n",err);
          break;
        }
        if(err > 0) {
          if(fwrite(buf,err,1,fout)!=1) {
            printf("MINIZIP> ERROR Writing extracted file.\n");
            err = UNZ_ERRNO;
            break;
          }
        }
      } while (err>0);
      if(fout)
        fclose(fout);

      if(err == 0)
        change_file_date(write_filename,file_info.dosDate,file_info.tmu_date);
    }

    if(err == UNZ_OK) {
      err = unzCloseCurrentFile (uf);
      if(err != UNZ_OK) {
        printf("MINIZIP> unzCloseCurrentFile reported error %d\n",err);
      }
    } else {
      unzCloseCurrentFile(uf); /* don't lose the error */
    }
  }

  free(buf);
  return err;
}

int do_extract(uf,opt_extract_without_path,opt_overwrite,password)
  unzFile uf;
  int opt_extract_without_path;
  int opt_overwrite;
  const char* password;
{
  uLong i;
  unz_global_info64 gi;
  int err;
  FILE* fout=NULL;

  err = unzGetGlobalInfo64(uf,&gi);
  if(err!=UNZ_OK)
    printf("MINIZIP> ERROR unzGetGlobalInfo reported error %d\n",err);

  for(i=0;i<gi.number_entry;i++)
  {
    if(do_extract_currentfile(uf,&opt_extract_without_path,&opt_overwrite,password) != UNZ_OK)
      break;

    if((i+1)<gi.number_entry) {
      err = unzGoToNextFile(uf);
      if(err!=UNZ_OK) {
        printf("MINIZIP> ERROR unzGoToNextFile reported error %d\n",err);
        break;
      }
    }
  }

  return 0;
}

int mz_mkdir(dirname)
  const char* dirname;
{
  int ret=0;
#ifdef _WIN32
  ret = _mkdir(dirname);
#elif unix
  ret = mkdir (dirname,0775);
#elif __APPLE__
  ret = mkdir (dirname,0775);
#endif
  return ret;
}

int makedir(newdir)
  const char *newdir;
{
  char *buffer ;
  char *p;
  int  len = (int)strlen(newdir);

  if(len <= 0)
    return 0;

  buffer = (char*)malloc(len+1);
  if(buffer==NULL) {
    printf("MINIZIP> ERROR When allocating memory\n");
    return UNZ_INTERNALERROR;
  }
  strcpy(buffer,newdir);

  if (buffer[len-1] == '/') {
    buffer[len-1] = '\0';
  }
  if(mz_mkdir(buffer) == 0) {
    free(buffer);
    return 1;
  }

  p = buffer+1;
  while(1) {
    char hold;

    while(*p && *p != '\\' && *p != '/')
      p++;
    hold = *p;
    *p = 0;
    if((mz_mkdir(buffer) == -1) && (errno == ENOENT)) {
      printf("MINIZIP> ERROR Couldn't create directory '%s'\n",buffer);
      free(buffer);
      return 0;
    }
    if(hold == 0)
      break;
    *p++ = hold;
  }
  free(buffer);
  return 1;
}
