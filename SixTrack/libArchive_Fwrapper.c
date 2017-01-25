#include "libArchive_wrapper.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

//Note underscore at the end of function definitions

void f_list_archive_(const char* infile, int infile_len){
  //Make sure to zero-terminate the string!
  char* infile_c = malloc((infile_len+1)*sizeof(char));
  strncpy(infile_c,infile,infile_len);
  infile_c[infile_len]=0;
  
  list_archive(infile_c);

  free(infile_c);
  fflush(stdout);
}


void f_read_archive_(const char* const infile, const char* extractFolder,int infile_len,int extractFolder_len) {
  printf("infile_len=%i, extactFolder_len=%i\n", infile_len, extractFolder_len);
  
  //Make sure to zero-terminate the string!
  char* infile_c = malloc((infile_len+1)*sizeof(char));
  strncpy(infile_c,infile,infile_len);
  infile_c[infile_len]=0;
  
  printf("infile_c = '%s'\n",infile_c);
  
  char* extractFolder_c = malloc((extractFolder_len+1)*sizeof(char));
  strncpy(extractFolder_c,extractFolder,extractFolder_len);
  extractFolder_c[extractFolder_len]=0;
  
  printf("extractFolder_c = '%s'\n",extractFolder_c);

  read_archive(infile_c, extractFolder_c);

  free(infile_c);
  free(extractFolder_c);
  fflush(stdout);
}

void f_write_archive_(const char* const outname, const char* const filenames, int* nFiles,
		      int outname_len, int filenames_len ){
  printf("nFiles=%i\n", *nFiles);
  printf("outname_len=%i\n", outname_len);
  printf("filenames_len=%i\n", filenames_len);

  //Manually print things
  printf("outname ='");
  for (int i = 0; i<outname_len;i++){
    printf("%c",outname[i]);
  }
  printf("'\n");

  printf("filenames=\n");
  for(int j=0;j<*nFiles;j++){
    printf("'");
    for (int i = 0; i<filenames_len;i++){
      printf("%c",filenames[i+filenames_len*j]);
    }
    printf("'\n");
  }
  printf("\n");

  //Build c-style string arrays
  char* outname_c = malloc(outname_len*sizeof(char));
  strncpy(outname_c,outname,outname_len);
  outname_c[outname_len]=0;
  printf("outname_c = '%s'\n",outname_c);
  
  char** filenames_c = malloc((*nFiles)*sizeof(char*));
  for (int i=0;i<(*nFiles);i++){
    char* str2 = malloc((filenames_len+1)*sizeof(char));
    strncpy(str2,&filenames[filenames_len*i],filenames_len);
    str2[filenames_len]=0;
    filenames_c[i]=str2;
    printf("filenames_c[%i] = '%s'\n",i,filenames_c[i]);

    //Trim trailing whitespace
    for (int j=filenames_len-1;j>=0;j--){
      if (filenames_c[i][j]==' ') {
	filenames_c[i][j] = '\0';
      }
      else{
	break;
      }
    }
    printf("filenames_c[%i] = '%s'\n",i,filenames_c[i]);
  }
  
  write_archive(outname_c,filenames_c,*nFiles);

  free(outname_c);
  for (int i=0;i<*nFiles;i++){
    free(filenames_c[i]);
  }
  free(filenames_c);
  fflush(stdout);
}
