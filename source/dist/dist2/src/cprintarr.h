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
	double **tas;
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

void a2cp(double tc[6], double cancord[6]);
void action2sixinternal_(double tc[6], double results[6]);
double momentum2energy(double momentum, double mass);
void  six2canonical_(double * coord, double *ref_momentum, double *mass, double *canonical);
void canonical2six_(double *canonical, double *ref_momentum, double *mass, double *coord);
void action2sixcoord_(double *e1, double *a1, double *e2, double *a2, double *e3, double *a3, double tas[6][6]
	,double *ref_momentum, double *mass, double [6]);