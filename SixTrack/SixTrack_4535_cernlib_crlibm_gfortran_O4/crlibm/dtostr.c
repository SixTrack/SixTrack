/*
  Compilation:
  gcc -m32 -std=c99 -W -Wall -pedantic -c dtostr.c.c
 */

// dtostr:  Fortran interface routine to sprintf

#include <errno.h>
#include <stdio.h>

#include <stdlib.h>
void enable_xp_(void);
void disable_xp_(void);
int dtostr_(double*anum,char*field,int flen)
{
  field[flen-1]='\0';
  char*format;
  format="%-23.16g";
//  printf ("dtostr: %e\n",*anum);
//  printf ("dtostr: %-23.15g\n",*anum); 
  errno=sprintf(field,format,*anum);
//  printf ("dtostr format:%s:\n",format);
//  printf ("dtostr field:%s:%i\n",field,flen);
  return errno;
}
