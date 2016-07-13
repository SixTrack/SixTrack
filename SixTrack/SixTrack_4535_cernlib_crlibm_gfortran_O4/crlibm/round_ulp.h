#ifndef ROUND_ULP_H
#define ROUND_ULP_H

// For Fortran compatibility
#define round2ulp(a)   round2ulp_(a)
#define round4ulp(a)   round4ulp_(a)
#define round8ulp(a)   round8ulp_(a)
#define roundnulp(a,n) roundnulp_(a,n)

// Portable rounding functions to ±2, ±4 and ±8 ulps (discading 2, 3 and 4 lsbs).
// thread safe after first call.
void round2ulp(double *x);
void round4ulp(double *x);
void round8ulp(double *x);

// Portable rounding functions to ±n ulps
// not thread safe (unless n never changes).
// return 0 on success, -1 on failure (bad n)
int roundnulp(double *x, int *n); // n must be a positive power of two (checked)

#endif // ROUND_ULP_H
