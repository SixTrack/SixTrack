#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cprintarr.h"

struct distparam* dist;

void 
six2canonical_(double * coord, double *ref_momentum, double *mass, double *canonical){
	double deltap = *(coord+5); 

	double beta0 = *ref_momentum/momentum2energy(*ref_momentum, *mass);
	double beta = (*ref_momentum+deltap)/momentum2energy((*ref_momentum)+deltap, *mass);
	double rv = beta0/beta;	

	 *(canonical+0) = *(coord+0);
	 *(canonical+1) = *(coord+1)*(1+deltap);
	 *(canonical+2) = *(coord+2);
	 *(canonical+3) = *(coord+3)*(1+deltap);
	 *(canonical+4) = *(coord+4)/rv;
	 *(canonical+5) = *(coord+5);

}

void 
canonical2six_(double *canonical, double *ref_momentum, double *mass, double *coord){
	double deltap = *(canonical+5);
	double beta0 = *ref_momentum/momentum2energy(*ref_momentum, *mass);
	double beta = (*ref_momentum+deltap)/momentum2energy((*ref_momentum)+deltap, *mass);
	double rv = beta0/beta;	
	*(coord+0) = *(canonical+0);
	*(coord+1) = *(canonical+1)/(1+deltap);
	*(coord+2) = *(canonical+2);
	*(coord+3) = *(canonical+3)/(1+deltap);
	*(coord+4) = *(canonical+4)*rv;
	*(coord+5) = *(canonical+5);	 
}

double momentum2energy(double momentum, double mass){
	return sqrt(pow(momentum,2)+pow(mass,2));
} 


/* 
    if mtrx_a is (m x n) and mtrx_b is (n x p), 
    the product is an (m x p) matrix
*/
void mtrx_vector_mult_(int *mp, int *np,  double mtrx_a[6][6], double mtrx_b[6], double result[6])
{

	int m = *mp;
	int n = *np;


    register int i=0, j=0, k=0;
    for (i = 0; i < m; i++)
    {
            for (k = 0; k < n; k++)
            {
                result [i] += mtrx_a [i][k] * mtrx_b [k];
                
            }
    }
}


void a2c_(double *e1, double *a1, double *e2, double *a2, double *e3, double *a3, double tas[6][6], double results [6]){
	double acoord[6];
	
	acoord[0]=sqrt(*e1)*cos(*a1);
	acoord[1]=-sqrt(*e1)*sin(*a1);
	acoord[2]=sqrt(*e2)*cos(*a2);
	acoord[3]=-sqrt(*e2)*sin(*a2);
	acoord[4]=sqrt(*e3)*cos(*a3);
	acoord[5]=-sqrt(*e3/1000)*sin(*a3);

	int dim, one;
	dim = 6;
	mtrx_vector_mult_(&dim,&dim, tas, acoord,results);

}
void setEmittance(double *e1, double *e2, double *e3){
		dist->emitt->e1=*e1; 
		dist->emitt->e2=*e2; 
		dist->emitt->e3=*e3; 
}
void initializedistribution_(int *numberOfDist, int *dimension){

	dist = (struct distparam*)malloc((*numberOfDist)*sizeof(struct distparam));
	int dim = *dimension;
	for(int i = 0; i <*numberOfDist; i++)
		{
		struct parameters para_tmp;
		(dist + i)->coord  =   (struct parameters**)malloc(6*sizeof(struct parameters*));
		(dist + i)->emitt =   (struct emittances*)malloc(sizeof(struct emittances));
		(dist + i)->mass  = 0;
		(dist + i)->momentum  = 0;
		(dist + i)->emitt->e1=0; 
		(dist + i)->emitt->e2=0; 
		(dist + i)->emitt->e3=0; 
		for(int j=0; j<dim; j++)
		{
			(dist +i)->coord[j] = (struct parameters*)malloc(sizeof(struct parameters));
			(dist +i)->coord[j]->start=0;
			(dist +i)->coord[j]->stop=0;
			(dist +i)->coord[j]->length=1;
			(dist +i)->coord[j]->type=-1;
		}
	}
}

void setmassmom(int ndist*, double *mass, double *momentum ){
	(dist+*ndist)-> mass = *mass;
	(dist+*ndist)-> momentum = *momentum;
}

void setParameter(int *index,  double *start, double *stop, int *length, int *type){
	struct parameters can;
	can.start = *start;
	can.stop = *stop;
	can.length = *length;
	can.type = *type;
	memcpy(dist->coord[*index], &can, sizeof(struct parameters));

}


double* createLinearSpaced(int length, double start, double stop){
	double eqspaced [length];
	double distance = (start-stop)/length;
	for(int i; i<length; i++){
		eqspaced[i] = distance*i;
	}
	return eqspaced;

}

void createTasWithNoCoupling(double betax, double alfax, double betay, double alfay, double tas[6][6]){

	    for (int i = 0; i < 6; i++)
    {
            for (int k = 0; k < 6; k++)
            {
                tas [i][k] = 0;
                
            }
    }
    tas[0][0] = sqrt(betax);
    tas[2][2] = sqrt(betay);
    tas[1][2] =-alfax/sqrt(betax);
    tas[4][3] =-alfay/sqrt(betay);

}

double RationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) / 
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double normalcdfinv_(double p)
{
	printf("printing %f", p);
    if (p <= 0.0 || p >= 1.0)
    {
    	return 0;
    }

    // See article above for explanation of this section.
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -RationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return RationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}

void action2sixcoord_(double *e1, double *a1, double *e2, double *a2, double *e3, double *a3, double tas[6][6]
	,double *ref_momentum, double *mass, double results [6]){
	double cancord[6] ;
	a2c_(e1, a1, e2, a2, e3, a3,tas, cancord);
	canonical2six_(cancord, ref_momentum, mass, results);
}
