#include <stdlib.h>
#include <limits.h>
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
  int fd;
  
  
  printf("Writing outname='%s'\n",outname);
  
  a = archive_write_new(); //Constructs the archive in memory? If so, may be problematic for large archives.
  //For tar.gz:
  // archive_write_add_filter_gzip(a);
  // archive_write_set_format_pax_restricted(a); // Note 1
  //For .zip:
  archive_write_set_format_zip(a);

  archive_write_open_filename(a, outname);
  entry = archive_entry_new();
  for (int i=0;i<nFiles;i++){
    printf("Compressing filename='%s'... ",filename[i]);

    //Write the header
    err = stat(filename[i], &st); // POSIX only, use GetFileSizeEx on Windows.
    printf("stat reported err=%i\n", err);
    if(err != 0) continue;
    
    archive_entry_set_pathname(entry, filename[i]);
    archive_entry_set_size(entry, st.st_size);
    archive_entry_set_filetype(entry, AE_IFREG);
    archive_entry_set_perm(entry, 0644);
    archive_write_header(a, entry);
    
    //Write the data
    fd = open(filename[i], O_RDONLY);
    len = read(fd, buff, sizeof(buff));
    while ( len > 0 ) {
      err=archive_write_data(a, buff, len);
      if (err < 0){
	printf("Error when writing file, got err=%i\n",err);
	exit(1);
      }
      len = read(fd, buff, sizeof(buff));
    }
    close(fd);
    
    archive_entry_clear(entry);
  }
  archive_entry_free(entry);
  printf("Complete!\n");
  archive_write_close(a); // Note 4
  archive_write_free(a); // called archive_write_finish() in old versions of libarchive
}

//********************************************************************************************************
void list_archive(const char* const infile) {
  // Adapted from:
  // https://github.com/libarchive/libarchive/wiki/Examples#List_contents_of_Archive_stored_in_File
  
  printf("Opening archive '%s' for listing...\n",infile);
  
  struct archive* a = archive_read_new();
  archive_read_support_format_zip(a);
  int err = archive_read_open_filename(a, infile,10240);//Note: Blocksize isn't neccessarilly adhered to
  if (err != ARCHIVE_OK) {
    printf("Error when opening archive '%s', err=%i\n",infile,err);
  }
  
  struct archive_entry* entry;
  while (archive_read_next_header(a,&entry)==ARCHIVE_OK){
    printf("Found file: '%s'\n",archive_entry_pathname(entry));
    archive_read_data_skip(a);
  }
  
  archive_read_close(a);
  err = archive_read_free(a);
  if (err != ARCHIVE_OK){
    printf("Error when calling archive_read_free(), '%s', err=%i\n",infile,err);
  }
}

//********************************************************************************************************
void read_archive(const char* const infile, const char* extractFolder){
  // Strongly inspired by
  // https://github.com/libarchive/libarchive/wiki/Examples#A_Complete_Extractor
  
  printf("Opening archive '%s' for extracting to folder '%s'...\n",infile,extractFolder);

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
    printf("Error opening archive '%s', err=%i\n",infile,err);
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
      printf("Error when reading archive, err=%i\n",err);
      printf("%s\n",archive_error_string(a));
      exit(1);
    }
    printf("Found file: '%s'\n",archive_entry_pathname(entry));

    //Avoid clobbering files in current directory - solution from
    // http://stackoverflow.com/questions/4496001/libarchive-to-extract-to-a-specified-folder
    char newPath[PATH_MAX];
    snprintf(newPath, PATH_MAX, "%s/%s",extractFolder,archive_entry_pathname(entry));
    archive_entry_set_pathname(entry,newPath);
    
    err = archive_write_header(ext, entry);
    if (err != ARCHIVE_OK){
      printf("Error when extracting archive (creating new file), err=%i\n",err);
      printf("%s\n",archive_error_string(ext));
      exit(1);
    }

    //Write the data!
    const void* buff;
    size_t size;
    off_t offset;
    
    const int bcount_max = 100000000;
    char bcompleted = 0; //C boolean
    for (int bcount=0; bcount<bcount_max;bcount++){
      err = archive_read_data_block(a,&buff,&size, &offset);
      if ( err == ARCHIVE_EOF ) {
	bcompleted=1;
	break;
      }
      else if (err != ARCHIVE_OK){
	printf("Error when extracting archive (reading data), err=%i\n",err);
	printf("%s\n",archive_error_string(a));
	exit(1);
      }

      err = archive_write_data_block(ext,buff,size,offset);
      if (err != ARCHIVE_OK){
	printf("Error when extracting archive (writing data), err=%i\n",err);
	printf("%s\n",archive_error_string(a));
	exit(1);
      }
    }
    if (!bcompleted){
      printf("Error: The file writing block loop was aborted by the infinite loop guard\n");
      exit(1);
    }
    
    err=archive_write_finish_entry(ext);
    if (err != ARCHIVE_OK) {
      printf("Error when extracting archive (closing new file), err=%i\n",err);
      printf("%s\n",archive_error_string(ext));
      exit(1);
    }
  }

  archive_read_close(a);
  err=archive_read_free(a);
  if (err != ARCHIVE_OK){
    printf("Error when calling archive_read_free(a), err=%i\n",err);
  }
  archive_write_close(ext);
  err = archive_write_free(ext);
  if (err != ARCHIVE_OK){
    printf("Error when calling archive_read_free(ext), err=%i\n",err);
  }
  
  if (!fcompleted) {
    printf("Error: The file header loop was aborted by the infinite loop guard\n");
    exit(1);
  }
}

//********************************************************************************************************
