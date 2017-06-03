/**
 * The interface routines for Fortran calls
 * An underscore appended to exp_rn and log_rn
 * The FPU flags are saved, set fro crlibm and restored
 * The argument(s) are treted as pointers

 * Author : Eric McIntosh (Eric.McIntosh@cern.ch)
 * Date of creation : 25/11/2004
 * Last Modified    : 25/11/2004
 */
#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"

/* An interface for exp_rn */
double exp_rn_(double *pointx) {
  double x,resultx;
  x=*pointx;
  resultx=exp_rn(x);
  return resultx;
}

/* An interface for log_rn */
double log_rn_(double *pointx) {
	  double x,resultx;
  x=*pointx;
  resultx=log_rn(x);
  return resultx;
}
/* An interface for log10_rn */
double log10_rn_(double *pointx) {
	  double x,resultx;
  x=*pointx;
  resultx=log10_rn(x);
  return resultx;
}
/* An interface for atan_rn */
double atan_rn_(double *pointx) {
          double x,resultx;
  x=*pointx;
  resultx=atan_rn(x);
  return resultx;
}
/* An interface for tan_rn */
double tan_rn_(double *pointx) {
          double x,resultx;
  x=*pointx;
  resultx=tan_rn(x);
  return resultx;
}
/* An interface for sin_rn */
double sin_rn_(double *pointx) {
          double x,resultx;
  x=*pointx;
  resultx=sin_rn(x);
  return resultx;
}
/* An interface for cos_rn */
double cos_rn_(double *pointx) {
          double x,resultx;
  x=*pointx;
  resultx=cos_rn(x);
  return resultx;
}
/* An interface for sinh_rn */
double sinh_rn_(double *pointx) {
          double x,resultx;
  x=*pointx;
  resultx=sinh_rn(x);
  return resultx;
}
/* An interface for cosh_rn */
double cosh_rn_(double *pointx) {
          double x,resultx;
  x=*pointx;
  resultx=cosh_rn(x);
  return resultx;
}
