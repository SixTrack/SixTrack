#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void
cprintarr_(float *arr, int *nx, int *ny)
{
    int i, j;
    printf("cprintarr_: arr=%p, nx=%d, ny=%d\n", arr, *nx, *ny);
    for(j = 0; j < *ny; ++j)
        for(i = 0; i < *nx; ++i)
        {
            int idx = j * *nx + i;
            printf("arr[%d] = %f\n", idx, arr[idx]);
        }
}

double momentum2energy(double momentum, double mass){
	return sqrt(pow(momentum,2)+pow(mass,2));
} 


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


/* 
    if mtrx_a is (m x n) and mtrx_b is (n x p), 
    the product is an (m x p) matrix
*/
void mtrx_vector_mult_(int *mp, int *np, int *pp, double mtrx_a[6][6], double mtrx_b[6], double result[6])
{

	int m = *mp;
	int n = *np;
	int p = *pp;


    register int i=0, j=0, k=0;
    for (i = 0; i < m; i++)
    {
            for (k = 0; k < n; k++)
            {
                result [i] += mtrx_a [i][k] * mtrx_b [k];
                
            }
        
    }

}
/*
void action2sixcoord_(double *e1, double *a1, double *e2, double *a2, double *e3, double *a3, double tas[6][6], double results [6]){
	a2c_(e1, a1, e2, a2, e3, a3,tas, results);
	canonical2six_(results, ref_momentum, mass, double canonical)
double 
}
*/

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
	one =1;
	mtrx_vector_mult_(&dim,&dim,&dim, tas, acoord,results);

} 

