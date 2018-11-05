#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//extern struct distparam* dist;
//extern int dim;
struct distparam
{
	struct parameters** coord;
	struct emittances* emitt;
	double mass;
	double momentum;
	int coordtype; // This tells which type of coordinates the input is given.  // 1-Normalized 
};

struct parameters
{
  double start;
  double stop;   
  int length;
  int type; //This gives the type of distribution, constant, linear, gaussian, 
  double * values; 
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