/*
  Compilation:
  gcc -m32 -std=c99 -W -Wall -pedantic -c dtoaf.c
 */

/*
#include <stdio.h>
*/
char *dtoa_cr(double d, int mode, int ndigits,
                        int *decpt, int *sign, char **rve);

int dtoaf_(double *d,int *mode, int *ndigits, int *decpt,
           int *sign, char *str)
{
  int len;
  char *rve;
  char *res;
  int i;

/*
  printf("%e %d %d \n",*d,*mode,*ndigits);
*/
  res=dtoa_cr(*d,*mode,*ndigits,decpt,sign,&rve);
/*
  printf("%s\n",res);
*/
  len=rve-res;
  for (i=0;i<=len-1;i++) {
    str[i]=res[i];
  }
  return len;
}
