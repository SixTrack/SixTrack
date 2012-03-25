/*
  Compilation:
  gcc -m32 -std=c99 -W -Wall -pedantic -c round_ulp.c
 */

// round_near: Fortran interface routine to strtod
// flen is 256 (0 to 255 in C)

#include <errno.h>
#include <stdio.h>

#include <stdlib.h>

double strtod(const char *nptr, char **endptr);

double round_near_(int*ferrno,char*field,int flen)
{
  char* endptr;
  double good;
  field[flen-1]='\0';
  printf("strtod Field: %s\n",field);
  errno=0;
  good = strtod(field,&endptr);
  *ferrno=errno;
  printf("strtod Errno: %i\n",errno);
  printf("strtod returning: %21.14G\n", good);
  if (*endptr != '\0') {       /* Not necessarily an error... */
    printf("strtod Extra characters: index %i:%s\n", *endptr,endptr);
  }
  *ferrno=errno;
  return good;
}
