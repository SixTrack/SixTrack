#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distr.h"


struct distparam* dist;
int dim;

void action2canonical_(double acangl[6], double cancord[6]){
	double acoord[6];
	acoord[0]= sqrt((dist->emitt->e1)*acangl[0])*cos(acangl[1]);
	acoord[1]=-sqrt((dist->emitt->e1)*acangl[0])*sin(acangl[1]);
	acoord[2]= sqrt((dist->emitt->e2)*acangl[2])*cos(acangl[3]);
	acoord[3]=-sqrt((dist->emitt->e2)*acangl[2])*sin(acangl[3]);
	acoord[4]= sqrt((dist->emitt->e3)*acangl[5])*cos(acangl[4]);
	acoord[5]=-sqrt((dist->emitt->e3)*acangl[5]/1000)*sin(acangl[4]);

	//This is the multiplication with the tas matrix 
	mtrx_vector_mult_pointer(dim,dim, dist->tas, acoord,cancord);

}
/*
Returns the total number of particles that will be created for the distribution.  
*/
int getnumberdist_(){
	double length = 1;
	for(int j =0; j<dim; j++){
		length = length*dist->coord[j]->length;
	}
	return length;
}

void setemittance12_(double *e1, double *e2){
	dist->emitt->e1=*e1; 
	dist->emitt->e2=*e2; 

}
void setemittance3_(double *e3){
	dist->emitt->e3=*e3;  	
}
void calcualteinverse_(){
	double invtas[25][25];
	double result[dim][dim];
	for(int i =0; i< dim; i++){
		for(int j =0; j< dim; j++){
			invtas[i][j] = dist->tas[i][j];
		}
	}
	cofactor(invtas, 6);
	for(int i =0; i< dim; i++){
		for(int j =0; j< dim; j++){
			dist->invtas[i][j]= invtas[i][j];
		}
	}


	printf("the inversa %f", invtas[0][0]);

}
void deltap(double *dp){
	dist->emitt->e3 = convertdp2emittance(*dp);
}

double convertdp2emittance(double dp){
	double emit=0;	
}

/*
Initalizes the distributinos 
*/
void initializedistribution_(int *numberOfDist, int *dimension){
	dist = (struct distparam*)malloc((*numberOfDist)*sizeof(struct distparam));
	dim  = *dimension;
	for(int i = 0; i <*numberOfDist; i++)
		{
		struct parameters para_tmp;
		(dist + i)->coord = (struct parameters**)malloc(dim*sizeof(struct parameters*));
		(dist + i)->emitt = (struct emittances*)malloc(sizeof(struct emittances));
		(dist + i)->tas   = (double**)malloc(dim*sizeof(double*));
		(dist + i)->invtas   = (double**)malloc(dim*sizeof(double*));
		
		for(int k=0; k<dim;k++){
			(dist + i)->tas[k] =(double*)malloc(dim*sizeof(double));
			(dist + i)->invtas[k] =(double*)malloc(dim*sizeof(double));
		}
		
		(dist + i)->mass      = 0;
		(dist + i)->momentum  = 0;
		(dist + i)->emitt->e1 = 0; 
		(dist + i)->emitt->e2 = 0; 
		(dist + i)->emitt->e3 = 0;
		(dist + i)->coordtype =-1;
		
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
			dist->tas[i][j] = tas[i][j];
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
							action2sixinternal_(tc, tmp);

						}
					}
				}
			}
		}	
	}
}

void setmassmom_(double *mass, double *momentum){
	(dist)-> mass = *mass;
	(dist)-> momentum = *momentum;
}


void setparameter_(int *index,  double *start, double *stop, int *length, int *type){

	dist->coord[*index-1]->start = *start;
	dist->coord[*index-1]->stop = *stop;
	dist->coord[*index-1]->length = *length;
	dist->coord[*index-1]->type = *type;
	if(*type ==0){ //Constant value 
		dist->coord[*index-1]->values = (double*)malloc(sizeof(double));
		dist->coord[*index-1]->values = start;
	}

	if(*type > 0){ //Allocate space for the array
		dist->coord[*index-1]->values = (double*)malloc((*length)*sizeof(double));
		memcpy(dist->coord[*index-1]->values , start, sizeof(double));  

	}
	if(*type==1){ //Linearly spaced intervalls
		createLinearSpaced(*length, *start, *stop,dist->coord[*index-1]->values);
	}
	if(*type==2){ //Exponentially spaced
		createLinearSpaced(*length, *start, *stop,dist->coord[*index-1]->values);
		for(int i=0;i <*length; i++){
			   
			dist->coord[*index-1]->values[i] = exp(dist->coord[*index-1]->values[i]);
		}
	}
	if(*type==3){ //Spaced with  ^2
		createLinearSpaced(*length, *start, *stop,dist->coord[*index-1]->values);
		for(int i=0;i <*length; i++){
			dist->coord[*index-1]->values[i] = pow(dist->coord[*index-1]->values[i],2);
		
		}

	}

}

void createtas0coupling_(double betax, double alfax, double betay, double alfay){

	    for (int i = 0; i < 6; i++)
    {
            for (int k = 0; k < 6; k++)
            {
                dist->tas [i][k] = 0;
                
            }
    }

    dist->tas[0][0] = sqrt(betax);
    dist->tas[2][2] = sqrt(betay);
    dist->tas[1][2] =-alfax/sqrt(betax);
    dist->tas[4][3] =-alfay/sqrt(betay);

}

void action2sixinternal_(double tc[6], double *results){
	double cancord[6];
	if(checkdist)
	{
		action2canonical_(tc, cancord);
		canonical2six_(cancord, &dist->momentum, &dist->mass, results);
		printf("%f %f %f %f %f %f \n", results[0],results[1],results[2],results[3],results[4],results[5]);
	}
}

int checkdist(){
	double eps =1e-16;
	if(dist->momentum < eps && dist->momentum < eps ){
		printf("Momentum and mass needs to be set");
		return 0;
	}
	if(dist->tas[0][0] < eps){
		printf("Tas matrix need to be set");
		return 0;
	}
	return 1;

}
