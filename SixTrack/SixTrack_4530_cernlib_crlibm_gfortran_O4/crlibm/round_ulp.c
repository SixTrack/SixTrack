/*
  Compilation:
  gcc -std=c99 -O3 -W -Wall -pedantic -c round_ulp.c

  Author:
  laurent.deniau@cern.ch, 2012
 */

#include <math.h>
#include <assert.h>
#include "round_ulp.h"

// union type
typedef union {
  double             fval;
  unsigned long long ival;
} value_t;

// compile-time assertion
enum {
  invalid_size_of_ulonglong_vs_double = 1/(sizeof(double) == sizeof(unsigned long long))
};

#if defined(TEST) || defined(SEARCH)
// forward declaration
static void print(value_t x);
#endif

// ---- round_ulp -------------------------------------------------------------


// round2ulp: portable rounding function to ±2 ulps (discading the 2 lsbs).
// thread safe after first call.

void
round2ulp(double *x_)
{
  enum { ulps = 2 };
  static unsigned long long mask = 0;
  static int done = 0;

  assert(x_);
  value_t x = { *x_ };

  // do not process NaN and Inf
  if (isfinite(x.fval)) {

    // setup the mask to truncate log2(ulps)+1 lsbs
    if (!done) {
      value_t fmask = { 0.0 };

      for (int i=1; i < 2*ulps; i++)
        fmask.fval = nextafter(fmask.fval, 1);

      mask = ~fmask.ival;
      done = 0;
    }

    // add ulps-1 times 1 ulp
    double y = 2*ulps*x.fval;

    for (int i=1; i < ulps; i++)
      x.fval = nextafter(x.fval, y);

#if defined(TEST) && !defined(SEARCH)
    print(x);
#endif

    // truncate lsbs
    x.ival &= mask;
    *x_ = x.fval;
  }
}

// round4ulp: portable rounding function to ±4 ulps (discading the 3 lsbs).
// thread safe after first call.

void
round4ulp(double *x_)
{
  enum { ulps = 4 };
  static unsigned long long mask = 0;
  static int done = 0;

  assert(x_);
  value_t x = { *x_ };

  // do not process NaN and Inf
  if (isfinite(x.fval)) {

    // setup the mask to truncate log2(ulps)+1 lsbs
    if (!done) {
      value_t fmask = { 0.0 };

      for (int i=1; i < 2*ulps; i++)
        fmask.fval = nextafter(fmask.fval, 1);

      mask = ~fmask.ival;
      done = 0;
    }

    // add ulps-1 times 1 ulp
    double y = 2*ulps*x.fval;

    for (int i=1; i < ulps; i++)
      x.fval = nextafter(x.fval, y);

#if defined(TEST) && !defined(SEARCH)
    print(x);
#endif

    // truncate lsbs
    x.ival &= mask;
    *x_ = x.fval;
  }
}

// round8ulp: portable rounding function to ±8 ulps (discading the 4 lsbs).
// thread safe after first call.

void
round8ulp(double *x_)
{
  enum { ulps = 8 };
  static unsigned long long mask = 0;
  static int done = 0;

  assert(x_);
  value_t x = { *x_ };

  // do not process NaN and Inf
  if (isfinite(x.fval)) {

    // setup the mask to truncate log2(ulps)+1 lsbs
    if (!done) {
      value_t fmask = { 0.0 };

      for (int i=1; i < 2*ulps; i++)
        fmask.fval = nextafter(fmask.fval, 1);

      mask = ~fmask.ival;
      done = 0;
    }

    // add ulps-1 times 1 ulp
    double y = 2*ulps*x.fval;

    for (int i=1; i < ulps; i++)
      x.fval = nextafter(x.fval, y);

#if defined(TEST) && !defined(SEARCH)
    print(x);
#endif

    // truncate lsbs
    x.ival &= mask;
    *x_ = x.fval;
  }
}

// roundnulp: portable rounding function to ±n ulps
// not thread safe

int
roundnulp(double *x_, int *n_)
{
  static unsigned long long mask = 0;
  static int ulps = 0;

  assert(x_ && n_);
  value_t x = { *x_ };

  // do not process NaN and Inf
  if (isfinite(x.fval)) {
    int n = *n_;

    // setup the mask to truncate log2(ulps)+1 lsbs
    if (n != ulps) {
      value_t fmask = { 0.0 };

      if (n < 2 || (n & (n-1))) return -1; // must be a positive power of two

      for (int i=1; i < 2*n; i++)
        fmask.fval = nextafter(fmask.fval, 1);

      mask = ~fmask.ival;
      ulps = n;
    }

    // add ulps-1 times 1 ulp
    double y = 2*ulps*x.fval;

    for (int i=1; i < ulps; i++)
      x.fval = nextafter(x.fval, y);

#if defined(TEST) && !defined(SEARCH)
    print(x);
#endif

    // truncate lsbs
    x.ival &= mask;
    *x_ = x.fval;
  }

  return 0;
}

// ---- testsuite -------------------------------------------------------------


/*
  Testsuite compilation:
  gcc -std=c99 -O3 -W -Wall -pedantic round_ulp.c -o round_ulp -lm
  testing for special cases: add -DTEST=n where n=2,4,8
  looking for invalid cases: add -DSEARCH=n where n=2,4,8
  parallel search          : add -fopenmp
  32-bit architecture      : add -m32
*/

#if defined(TEST) || defined(SEARCH)
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>

#ifdef TEST
#define round_ulp round_ulp_name(TEST)
#else
#define round_ulp round_ulp_name(SEARCH)
#endif

#define round_ulp_name( a) round_ulp_name_(a)
#define round_ulp_name_(a) round ## a ## ulp

static void
print(value_t x)
{
  printf("%+#.13A [%016llX]", x.fval, x.ival);
}

static void
check_rounding(double a)
{
  for (int i = 0; i <= 20; i++) {
    printf("round%02d: ", i);
    value_t x = { a+(nextafter(a,2*a)-a)*i };
    value_t y = x;
    round_ulp(&y.fval);
    printf(" -> "); print(x); printf(" -> "); print(y); printf("\n");
  }
}

#ifdef SEARCH
static void
invalid_case(int i, double x, double e)
{
  printf("** invalid case detected for i=%d, x=%.17e, e=%.17e\n", i, x, e);
  check_rounding(x);
  exit(EXIT_FAILURE);
}

static void
check_case(double x, int loops)
{
  double y=0, z;

  for (int i = 0; i < loops; i++) {
    z = x; round_ulp(&z);
#ifdef XCHECK
    int n = SEARCH; 
    double z2 = x;
    assert(!roundnulp(&z2, &n));
    if (z != z2) invalid_case(i, x, x-z); 
#endif
    if (i % (SEARCH*2) == 0) y = z;
    if (y != z) invalid_case(i, x, x-y);
    x = nextafter(x,2*x);
  }
}

static double
init_case(double x)
{
  double y = x*(SEARCH*2*2);

  round_ulp(&x);
  for (int i = 0; i < SEARCH+1; i++)
    x = nextafter(x,y);

  return x;
}

static unsigned
next_seed(unsigned a)
{
  // group generator for any \frac{\setN}{2^k\setN}, k=1..32
  // must start with a=1
  return a * 2621124293u + 1;
}
#endif

int main(void)
{
  static const double pi = 3.141592653589793238462643383279502884197169399;
  static const double e  = 2.718281828459045235360287471352662497757247094;
  static const double sn = 0X1.ffffffffffff0P-1022; // close to max subnormal

#ifdef TEST
  // denormalized positive ulps
  printf("** denormalized positive ulps at +0.0\n");
  check_rounding(nextafter(+0.0,1));

  // denormalized negative ulps
  printf("** denormalized negative ulps at -0.0\n");
  check_rounding(nextafter(-0.0,-1));

  // denormalized positive ulps
  printf("** denormalized maximum positive at +sn\n");
  check_rounding(+sn);

  // denormalized negative ulps
  printf("** denormalized maximum negative at -sn\n");
  check_rounding(-sn);

  // normalized positive ulps
  printf("** normalized positive ulps at +1.0\n");
  check_rounding(+1.0);

  // normalized negative ulps
  printf("** normalized negative ulps at -1.0\n");
  check_rounding(-1.0);

  // normalized positive ulps
  printf("** normalized positive ulps at +pi\n");
  check_rounding(+pi);

  // normalized negative ulps
  printf("** normalized negative ulps at -pi\n");
  check_rounding(-pi);

  // normalized positive ulps
  printf("** normalized positive ulps at +1/pi\n");
  check_rounding(+1.0/pi);

  // normalized negative ulps
  printf("** normalized negative ulps at -1/pi\n");
  check_rounding(-1.0/pi);

  // normalized positive ulps
  printf("** normalized positive ulps at +e\n");
  check_rounding(+e);

  // normalized negative ulps
  printf("** normalized negative ulps at -e\n");
  check_rounding(-e);

  // normalized positive ulps
  printf("** normalized positive ulps at +1/e\n");
  check_rounding(+1.0/e);

  // normalized negative ulps
  printf("** normalized negative ulps at -1/e\n");
  check_rounding(-1.0/e);
#endif

#ifdef SEARCH
  // search for invalid cases (>1 million)
  enum { loops = 10000000 };

#pragma omp parallel
#pragma omp sections
{
  {
    printf("** searching for invalid cases at sn\n");
    check_case(init_case(sn), loops);
  }

#pragma omp section
  {
    printf("** searching for invalid cases at pi\n");
    check_case(init_case(pi), loops);
  }

#pragma omp section
  {
    printf("** searching for invalid cases at 1/pi\n");
    check_case(init_case(1/pi), loops);
  }

#pragma omp section
  {
    printf("** searching for invalid cases at e\n");
    check_case(init_case(e), loops);
  }

#pragma omp section
  {
    printf("** searching for invalid cases at 1/e\n");
    check_case(init_case(1/e), loops);
  }
}

  {
    enum { nperiods = 10, inloops = nperiods*SEARCH*2 }; 

    srand(time(0));
    int l = ((1<<24)-1) & rand();
    unsigned s = 1;
    while (l--) s = next_seed(s);

    printf("** searching for invalid cases at %d random positions checking %d periods\n",
           loops, nperiods);

#pragma omp parallel for shared(s)
    for (int i = 0; i < loops; i++) {
      double x = (s = next_seed(s));
      check_case(init_case(x), inloops);
      double y = 1/x;
      check_case(init_case(y), inloops);
    }
  }
#endif
}
#endif
