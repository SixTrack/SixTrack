#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cprintarr.h"

struct distparam* dist;
int dim;
void six2canonical_(double * coord, double *ref_momentum, double *mass, double *canonical){
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

void print2file(double **coord, int length){
	FILE *fptr;

   fptr = fopen("coordinates.out", "w");
   if(fptr == NULL)
   {
      printf("Error!");
      exit(1);
   }
   for(int i=0;i < length; i++)
   {
   		fprintf(fptr,"%f  %f  %f  %f %f  %f  \n", coord[i][1], coord[i][2],coord[i][3], coord[i][4], coord[i][5], coord[i][6]);
   }
   fclose(fptr);
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
            	//printf("NOTTT matrix %f, %d, %d \n", mtrx_a [i][k], i, k );
                result [i] += mtrx_a [i][k] * mtrx_b [k];
                   //  printf("bbbbb000 %f \n", mtrx_b[k]);
                
            }
    }
      //printf(" ab %f %f %f %f %f %f \n", result[0],result[1],result[2],result[3],result[4],result[5]);
}

void mtrx_vector_mult_pointer(int mp, int np,  double **mtrx_a, double mtrx_b[6], double result[6])
{

	int m = mp;
	int n = np;


    register int i=0, j=0, k=0;
    for (i = 0; i < m; i++)
    {
            for (k = 0; k < n; k++)
            {
            	//printf("pointer matrix %f, %d, %d \n", mtrx_a [i][k], i, k );
                result [i] += mtrx_a [i][k] * mtrx_b [k];
              //  printf("bbbbb pointer %f \n", result [i]);
                
            }
    }
   // printf("%f %f %f %f %f %f \n", result[0],result[1],result[2],result[3],result[4],result[5]);
}


void a2c_(double *e1, double *a1, double *e2, double *a2, double *e3, double *a3, double tas[6][6], double results [6]){
	double acoord[6];
	
	acoord[0]=sqrt(*e1)*cos(*a1);
	acoord[1]=-sqrt(*e1)*sin(*a1);
	acoord[2]=sqrt(*e2)*cos(*a2);
	acoord[3]=-sqrt(*e2)*sin(*a2);
	acoord[4]=sqrt(*e3)*cos(*a3);
	acoord[5]=-sqrt(*e3/1000)*sin(*a3);
	printf("%f %f %f %f %f %f \n", acoord[0],acoord[1],acoord[2],acoord[3],acoord[4],acoord[5]);
	int dim, one;
	dim = 6;
	mtrx_vector_mult_(&dim,&dim, tas, acoord,results);
	//printf("%f %f %f %f %f %f \n", acoord[0],acoord[1],acoord[2],acoord[3],acoord[4],acoord[5]);

}
void a2cp(double tc[6], double cancord[6]){

	double acoord[6];
	//printf("accord %f", accord[0]);
	acoord[0]= sqrt((dist->emitt->e1)*tc[0])*cos(tc[1]);
	acoord[1]=-sqrt((dist->emitt->e1)*tc[0])*sin(tc[1]);
	acoord[2]= sqrt((dist->emitt->e2)*tc[2])*cos(tc[3]);
	acoord[3]=-sqrt((dist->emitt->e2)*tc[2])*sin(tc[3]);
	acoord[4]= sqrt((dist->emitt->e3)*tc[5])*cos(tc[4]);
	acoord[5]=-sqrt((dist->emitt->e3)*tc[5]/1000)*sin(tc[4]);
    int dima = 6;
	printf("%f %f %f %f %f %f \n", acoord[0],acoord[1],acoord[2],acoord[3],acoord[4],acoord[5]);
	mtrx_vector_mult_pointer(dima,dima, dist->tas, acoord,cancord);
	//printf("%f %f %f %f %f %f \n", cancord[0],cancord[1],cancord[2],cancord[3],cancord[4],cancord[5]);
}
void createDistribution(){

}

int getlengthdist(){
	return 100;
}
void setemittance_(double *e1, double *e2, double *e3){
		dist->emitt->e1=*e1; 
		dist->emitt->e2=*e2; 
		dist->emitt->e3=*e3; 
}

void initializedistribution_(int *numberOfDist, int *dimension){

	dist = (struct distparam*)malloc((*numberOfDist)*sizeof(struct distparam));
	dim = *dimension;
	for(int i = 0; i <*numberOfDist; i++)
		{
		struct parameters para_tmp;
		(dist + i)->coord  =   (struct parameters**)malloc(6*sizeof(struct parameters*));
		(dist + i)->emitt =   (struct emittances*)malloc(sizeof(struct emittances));
		(dist + i)->tas =(double**)malloc(dim*sizeof(double*));
		for(int k=0; k<dim;k++){
			(dist + i)->tas[k] =(double*)malloc(dim*sizeof(double));
		}
		(dist + i)->mass  = 0;
		(dist + i)->momentum  = 0;
		(dist + i)->emitt->e1=0; 
		(dist + i)->emitt->e2=0; 
		(dist + i)->emitt->e3=0;
		(dist + i)->coordtype=-1;
		for(int j=0; j<dim; j++)
		{
			(dist +i)->coord[j] = (struct parameters*)malloc(sizeof(struct parameters));
			(dist +i)->coord[j]->start=0;
			(dist +i)->coord[j]->stop=0;
			(dist +i)->coord[j]->length=1;
			(dist +i)->coord[j]->type=0;
		}
	}
}

void printdistsettings_(int *ndist){
	printf("Printing info about distribution: \n");
	printf("Distribution type (input), 1=normalized coordinates %d \n", dist->coordtype );
	printf("Mass: %f \n", dist ->mass);
	printf("Momentum: %f \n", dist ->momentum);
	printf("Emttiance, e1, e2, e3: %f, %f, %f \n", dist->emitt->e1, dist->emitt->e2, dist->emitt->e3);
	for(int j=0; j<dim; j++)
	{
		printf("coordinate index: %d \n",(j+1));
		printf("start: %f \n", (dist)->coord[j]->start);
		printf("stop: %f \n", (dist)->coord[j]->stop);
		printf("length: %d \n", (dist)->coord[j]->length);
		printf("type: %d \n", (dist)->coord[j]->type);
		printf("######################### \n");
	}

}

void settasmatrix_(double tas[6][6]){
	for(int i =0; i< dim; i++){
		for(int j =0; j< dim; j++){
			printf("howloooong,,,,%d %d %f \n", i, j, tas[i][j]);
			dist->tas[j][j] = tas[i][j];
		}
	}
}

//int arr[numRows][numCols]
void dist2sixcoord_(double **results){
	int counter = 0;
	double tc[6];
	double tmp[6];
	for(int i =0; i< dist->coord[0]->length; i++){
		for(int j =0; j< dist->coord[1]->length; j++){
			for(int k =0; k< dist->coord[2]->length; k++){
				for(int l =0; l< dist->coord[3]->length; l++){
					for(int m =0; m< dist->coord[4]->length; m++){
						for(int n =0; n< dist->coord[5]->length; n++){
							tc[0]=dist->coord[0]->values[i];
							tc[1]=dist->coord[1]->values[j];
							tc[2]=dist->coord[2]->values[k];
							tc[3]=dist->coord[3]->values[l];
							tc[4]=dist->coord[4]->values[m];
							tc[5]=dist->coord[5]->values[n];

							//sprintf("%f %f %f %f %f %f \n", tc[0],tc[1],tc[2],tc[3],tc[4],tc[5]);
							action2sixinternal_(tc, tmp);




						}
					}
				}
			}
		}	
	}
}

void setmassmom_(double *mass, double *momentum ){
	(dist)-> mass = *mass;
	(dist)-> momentum = *momentum;
}

void createLinearSpaced(int length, double start, double stop, double *eqspaced ){

	double distance = (start-stop)/length;
	for(int i; i<length; i++){
		eqspaced[i] = distance*i;
	}
	
}

void setparameter_(int *index,  double *start, double *stop, int *length, int *type){

	dist->coord[*index-1]->start = *start;
	dist->coord[*index-1]->stop = *stop;
	dist->coord[*index-1]->length = *length;
	dist->coord[*index-1]->type = *type;
	if(*type ==0){
		dist->coord[*index-1]->values = (double*)malloc(sizeof(double));
		dist->coord[*index-1]->values = start;
	}

	if(*type > 0){
		dist->coord[*index-1]->values = (double*)malloc((*length)*sizeof(double));
		memcpy(dist->coord[*index-1]->values , start, sizeof(double));  

	}
	if(*type==1){
		createLinearSpaced(*length, *start, *stop,dist->coord[*index-1]->values);
	}
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


//void action2sixcoordv_(double *acoord, double tas[6][6],double *ref_momentum, double *mass, double results [6]){
	//double cancord[6] ;
//	a2c_(acoord[0], acoord[1], acoord[2], acoord[3], acoord[4], acoord[0],tas, cancord);
//	canonical2six_(cancord, ref_momentum, mass, results);
//}
void action2sixinternal_(double tc[6], double results[6]){
	double cancord[6];
	//printf("%f %f %f %f %f %f \n", tc[0],tc[1],tc[2],tc[3],tc[4],tc[5]);
	a2cp(tc, cancord);
	//printf("%f %f %f %f %f %f \n", cancord[0],cancord[1],cancord[2],cancord[3],cancord[4],cancord[5]);

	canonical2six_(cancord, &dist->momentum, &dist->mass, results);
	//printf("%f %f %f %f %f %f \n", results[0],results[1],results[2],results[3],results[4],results[5]);

}

void action2sixcoord_(double *e1, double *a1, double *e2, double *a2, double *e3, double *a3, double tas[6][6]
	,double *ref_momentum, double *mass, double results [6]){
	double cancord[6] ;
	a2c_(e1, a1, e2, a2, e3, a3,tas, cancord);
	canonical2six_(cancord, ref_momentum, mass, results);
}
