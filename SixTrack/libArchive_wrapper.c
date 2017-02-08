#include <stdlib.h>

#include <limits.h>
#ifndef PATH_MAX
#define PATH_MAX 4096 //Not really good... But C99 does not have this constant.
#endif

#include <fcntl.h>

#include <archive.h>
#include <archive_entry.h>

#include "libArchive_wrapper.h"

//********************************************************************************************************
void write_archive(const char* const outname, char** filename, int nFiles) {
  // Adapted from
  // https://github.com/libarchive/libarchive/wiki/Examples#A_Basic_Write_Example
  struct archive *a;
  struct archive_entry *entry;
  struct stat st;
  int err;
  char buff[8192];
  int len;
  FILE* fd;
  
  
  //printf("Writing outname='%s'\n",outname);
  
  a = archive_write_new(); //Constructs the archive in memory? If so, may be problematic for large archives.
  //For tar.gz:
  // archive_write_add_filter_gzip(a);
  // archive_write_set_format_pax_restricted(a); // Note 1
  //For .zip:
  archive_write_set_format_zip(a);

  archive_write_open_filename(a, outname);
  entry = archive_entry_new();
  for (int i=0;i<nFiles;i++){
    //printf("Compressing filename='%s'... ",filename[i]);

    //Write the header
#ifdef __WIN32
    LARGE_INTEGER filesize_union;
    HANDLE hFile = CreateFile(filename[i], GENERIC_READ,
			      FILE_SHARE_READ | FILE_SHARE_WRITE,
			      NULL, OPEN_EXISTING,
			      FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE) continue;
    if(!GetFileSizeEx(hFile, &filesize_union)) {
      printf("CRITICAL ERROR in write_archive(): GetFileSizeEx failed.");
      exit(1);
    }
    CloseHandle(hFile);
    LONGLONG filesize = filesize_union.QuadPart;
#else
    off_t filesize = 0;
    err = stat(filename[i], &st);
    //printf("stat reported err=%i\n", err);
    if(err != 0) continue;
    filesize = st.st_size;
#endif
    
    archive_entry_set_pathname(entry, filename[i]);
    archive_entry_set_size(entry, filesize);
    archive_entry_set_filetype(entry, AE_IFREG);
    archive_entry_set_perm(entry, 0644);
    archive_write_header(a, entry);
    
    //Write the data
    int len_total = 0;
    fd = fopen(filename[i], "rb");
    len = fread(buff, 1,sizeof(buff), fd);
    len_total = len;
    while ( len > 0 ) {
      err=archive_write_data(a, buff, len);
      if (err < 0){
	printf("CRITICAL ERROR in write_archive(): When writing file, got err=%i\n",err);
	printf("CRITICAL ERROR in write_archive(): %s\n",archive_error_string(a));
	exit(1);
      }
      len = fread(buff, 1,sizeof(buff), fd);
      len_total += len;
    }
    fclose(fd);
    archive_entry_clear(entry);
    
    //printf("Wrote file '%s', len_total = %i, filesize = %i\n", filename[i], len_total, filesize);
    if (len_total != filesize) {
      printf("CRITICAL ERROR in write_archive(): When writing file '%s', got len_total = %i but filesize = %i\n", filename[i], len_total, filesize);
    }
  }
  archive_entry_free(entry);
  //printf("Complete!\n");
  archive_write_close(a); // Note 4
  archive_write_free(a); // called archive_write_finish() in old versions of libarchive
}

//********************************************************************************************************
void list_archive(const char* const infile) {
  // Adapted from:
  // https://github.com/libarchive/libarchive/wiki/Examples#List_contents_of_Archive_stored_in_File
  // Note - this function will always be "chatty"; make sure to flush the stdout before using it...
  // TODO: Write a function which *returns* the list of files
  //       (by writing it into a char** supplied by the caller)
  
  printf("Opening archive '%s' for listing...\n",infile);
  
  struct archive* a = archive_read_new();
  archive_read_support_format_zip(a);
  int err = archive_read_open_filename(a, infile,10240);//Note: Blocksize isn't neccessarilly adhered to
  if (err != ARCHIVE_OK) {
    printf("CRITICAL ERROR in list_archive(): When opening archive '%s', err=%i\n",infile,err);
    printf("CRITICAL ERROR in list_archive(): %s\n",archive_error_string(a));
    exit(1);
  }
  
  struct archive_entry* entry;
  while (archive_read_next_header(a,&entry)==ARCHIVE_OK){
    printf("Found file: '%s'\n",archive_entry_pathname(entry));
    archive_read_data_skip(a);
  }
  
  archive_read_close(a);
  err = archive_read_free(a);
  if (err != ARCHIVE_OK){
    printf("CRITICAL ERROR in list_archive(): Error when calling archive_read_free(), '%s', err=%i\n",infile,err);
    printf("CRITICAL ERROR in list_archive(): %s\n",archive_error_string(a));
    exit(1);
  }
}

void list_archive_get(const char* const infile, char** filenames, int* nfiles, const int buffsize) {
  //printf("In list_archive_get\n");
  //fflush(stdout);
  
  struct archive* a = archive_read_new();
  archive_read_support_format_zip(a);
  int err = archive_read_open_filename(a, infile,10240);//Note: Blocksize isn't neccessarilly adhered to
  if (err != ARCHIVE_OK) {
    printf("CRITICAL ERROR in list_archive_get(): When opening archive '%s', err=%i\n",infile,err);
    printf("CRITICAL ERROR in list_archive_get(): %s\n",archive_error_string(a));
    exit(1);
  }
  
  const int nfiles_max = *nfiles;
  *nfiles = 0;
  
  struct archive_entry* entry;
  while (archive_read_next_header(a,&entry)==ARCHIVE_OK){
    //printf("Found file: '%s'\n",archive_entry_pathname(entry));
    //fflush(stdout);
    int buff_used = 0;
    buff_used = snprintf(filenames[*nfiles],buffsize,"%s",archive_entry_pathname(entry));
    if (buff_used >= buffsize) {
      printf("CRITICAL ERROR in list_archive_get(): When reading file '%s' from archive '%s':\n",filenames[*nfiles],infile);
      printf("CRITICAL ERROR in list_archive_get(): Buffer too small by %i characters\n",buff_used-buffsize+1);
      exit(1);
    }
    else if (buff_used < 0) {
      printf("CRITICAL ERROR in list_archive_get(): When reading file '%s' from archive '%s':\n",filenames[*nfiles],infile);
      printf("CRITICAL ERROR in list_archive_get(): Error in snprintf.\n");
      exit(1);
    }
    
    archive_read_data_skip(a);
    
    if(++(*nfiles) >= nfiles_max) {
      printf("CRITICAL ERROR in list_archive_get(): Number of files greater than nfiles_max=%i",nfiles_max);
      exit(1);
    }
  }
  
  archive_read_close(a);
  err = archive_read_free(a);
  if (err != ARCHIVE_OK){
    printf("CRITICAL ERROR in list_archive_get(): Error when calling archive_read_free(), '%s', err=%i\n",infile,err);
    printf("CRITICAL ERROR in list_archive_get(): %s\n",archive_error_string(a));
    exit(1);
  }
}

//********************************************************************************************************
void read_archive(const char* const infile, const char* const extractFolder){
  // Strongly inspired by
  // https://github.com/libarchive/libarchive/wiki/Examples#A_Complete_Extractor
  
  //printf("Opening archive '%s' for extracting to folder '%s'...\n",infile,extractFolder);

  //Check that the archive exists

  //Check that the folder exists, if not then create it
  
  struct archive* a = archive_read_new();
  archive_read_support_format_zip(a);
  struct archive* ext = archive_write_disk_new();
  archive_write_disk_set_options(ext,ARCHIVE_EXTRACT_TIME|ARCHIVE_EXTRACT_PERM|ARCHIVE_EXTRACT_ACL|ARCHIVE_EXTRACT_FFLAGS);
  archive_write_disk_set_standard_lookup(ext);
  
  int err;
  err = archive_read_open_filename(a, infile, 10240);
  if (err != ARCHIVE_OK) {
    printf("CRITICAL ERROR in read_archive(): When opening archive '%s', err=%i\n",infile,err);
    printf("CRITICAL ERROR in read_archive(): %s\n",archive_error_string(a));
    exit(1);
  }

  struct archive_entry *entry;
  
  const int fcount_max = 1000;
  char fcompleted=0; //C-Boolean
  for(int fcount=0; fcount<fcount_max;fcount++){
    err = archive_read_next_header(a,&entry);
    if (err == ARCHIVE_EOF){
      fcompleted=1;
      break;
    }
    else if (err != ARCHIVE_OK){
      printf("CRITICAL ERROR in read_archive(): When reading archive, err=%i\n",err);
      printf("CRITICAL ERROR in read_archive(): %s\n",archive_error_string(a));
      exit(1);
    }
    //printf("Found file: '%s'\n",archive_entry_pathname(entry));

    //Avoid clobbering files in current directory - solution from
    // http://stackoverflow.com/questions/4496001/libarchive-to-extract-to-a-specified-folder
    char newPath[PATH_MAX];
    if (snprintf(newPath, PATH_MAX, "%s/%s",extractFolder,archive_entry_pathname(entry)) >= PATH_MAX){
      printf("CRITICAL ERROR in read_archive(): Buffer overflow when creating the path.");
      exit(1);
    }
    archive_entry_set_pathname(entry,newPath);
    
    err = archive_write_header(ext, entry);
    if (err != ARCHIVE_OK){
      printf("CRITICAL ERROR in read_archive(): when extracting archive (creating new file), err=%i\n",err);
      printf("CRITICAL ERROR in read_archive(): %s\n",archive_error_string(ext));
      exit(1);
    }

    //Write the data!
    const void* buff;
    size_t size;
    la_int64_t offset;
    
    const int bcount_max = 100000000;
    char bcompleted = 0; //C boolean
    for (int bcount=0; bcount<bcount_max;bcount++){
      err = archive_read_data_block(a,&buff,&size, &offset);
      if ( err == ARCHIVE_EOF ) {
	bcompleted=1;
	break;
      }
      else if (err != ARCHIVE_OK){
	printf("CRITICAL ERROR in read_archive(): When extracting archive (reading data), err=%i\n",err);
	printf("CRITICAL ERROR in read_archive(): %s\n",archive_error_string(a));
	exit(1);
      }

      err = archive_write_data_block(ext,buff,size,offset);
      if (err != ARCHIVE_OK){
	printf("CRITICAL ERROR in read_archive(): When extracting archive (writing data), err=%i\n",err);
	printf("CRITICAL ERROR in read_archive(): %s\n",archive_error_string(a));
	exit(1);
      }
    }
    if (!bcompleted){
      printf("CRITICAL ERROR in read_archive(): The file writing block loop was aborted by the infinite loop guard\n");
      exit(1);
    }
    
    err=archive_write_finish_entry(ext);
    if (err != ARCHIVE_OK) {
      printf("CRITICAL ERROR in read_archive(): When extracting archive (closing new file), err=%i\n",err);
      printf("CRITICAL ERROR in read_archive(): %s\n",archive_error_string(ext));
      exit(1);
    }
  }

  archive_read_close(a);
  err=archive_read_free(a);
  if (err != ARCHIVE_OK){
    printf("CRITICAL ERROR in read_archive(): When calling archive_read_free(a), err=%i\n",err);
    printf("CRITICAL ERROR in read_archive(): %s\n",archive_error_string(a));
    exit(1);
  }
  archive_write_close(ext);
  err = archive_write_free(ext);
  if (err != ARCHIVE_OK){
    printf("CRITICAL ERROR in read_archive(): When calling archive_read_free(ext), err=%i\n",err);
    printf("CRITICAL ERROR in read_archive(): %s\n",archive_error_string(a));
    exit(1);
  }
  
  if (!fcompleted) {
    printf("CRITICAL ERROR in read_archive(): The file header loop was aborted by the infinite loop guard\n");
    exit(1);
  }
}

//********************************************************************************************************
