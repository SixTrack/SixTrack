#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct distpara
{
	struct parameters* coord;
	struct emittances* emitt;
	int dim;
	double mass;
	double momentum;
};

struct parameters
{
  double value;
  double start;
  double stop;   
  int length;
  int type;
};

struct emittances{
	double e1, e2, e3;
};


double momentum2energy(double momentum, double mass);
void  six2canonical_(double * coord, double *ref_momentum, double *mass, double *canonical);
void canonical2six_(double *canonical, double *ref_momentum, double *mass, double *coord);
void mtrx_vector_mult_(int *mp, int *np,  double mtrx_a[6][6], double mtrx_b[6], double result[6]);
void a2c_(double *e1, double *a1, double *e2, double *a2, double *e3, double *a3, double tas[6][6], double results [6]); 
void action2sixcoord_(double *e1, double *a1, double *e2, double *a2, double *e3, double *a3, double tas[6][6]
	,double *ref_momentum, double *mass, double [6]);