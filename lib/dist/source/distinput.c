#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distinput.h"
#include "distgeneration.h"
#include "outputdist.h"

/*
This allocates the the memory for the distributions
*/

void initializedistribution_(int *numberOfDist){
    dist = (struct distparam*)malloc((*numberOfDist)*sizeof(struct distparam));
    dim  = 6;
    
        for(int i = 0; i <*numberOfDist; i++)
        {
        struct parameters para_tmp;
        (dist + i)->coord = (struct parameters**)malloc(dim*sizeof(struct parameters*));
        (dist + i)->emitt = (struct emittances*)malloc(sizeof(struct emittances));
        (dist + i)->cuts2apply = (struct appliedcut*)malloc(sizeof(struct appliedcut));
        (dist + i)->cuts2apply->physical = (struct cut**)malloc(dim*sizeof(struct cut*));
        (dist + i)->cuts2apply->normalized = (struct cut**)malloc(dim*sizeof(struct cut*));
        (dist + i)->tas   = (double**)malloc(dim*sizeof(double*));
        (dist + i)->invtas   = (double**)malloc(dim*sizeof(double*));
        (dist + i)->closedorbit   = (double*)malloc(dim*sizeof(double));
        (dist + i)->isDistrcalculated = 0;
        (dist + i)->longitunalemittance = 0;

        for(int k=0; k<dim;k++){
            (dist + i)->tas[k] =(double*)malloc(dim*sizeof(double));
            (dist + i)->invtas[k] =(double*)malloc(dim*sizeof(double));
        }
        
        (dist + i)->mass      = 0;
        (dist + i)->momentum  = 0;
        (dist + i)->charge    = 1;
        (dist + i)->massnum   = 1;
        (dist + i)->atomnum   = 1;
        (dist + i)->emitt->e1 = 0; 
        (dist + i)->emitt->e2 = 0; 
        (dist + i)->emitt->e3 = 0;
        (dist + i)->emitt->dp = 0;
        (dist + i)->emitt->deltas = 0;
        (dist + i)->coordtype   =-1;
        (dist + i)->totallength = -1;
        (dist + i)->disttype = 0;
        for(int j=0; j<dim; j++)
        {
            (dist + i)->cuts2apply->physical[j] = (struct cut*)malloc(sizeof(struct cut));
            (dist + i)->cuts2apply->normalized[j] = (struct cut*)malloc(sizeof(struct cut));
            (dist +i)->coord[j] = (struct parameters*)malloc(sizeof(struct parameters));
            (dist +i)->coord[j]->start=0;
            (dist +i)->coord[j]->stop=0;
            (dist +i)->coord[j]->length=1;
            (dist +i)->coord[j]->type=0;
            (dist +i)->closedorbit[j]=0;
        }
    }
    diststart=dist;
}

void settotallength(int totallength){
	dist->totallength=totallength;
}

void setdisttype(int disttype){
	dist->disttype=disttype;
}


void setdistribution_(int *ndist){
		dist = diststart + *ndist;
}

void setemittance12_(double *e1, double *e2){
	dist->emitt->e1=*e1; 
	dist->emitt->e2=*e2; 

}

void setemittance3_(double *e3){
	dist->emitt->e3=*e3;
	dist->longitunalemittance=0;  	
}


void addclosedorbit_(double *clo){
	for(int i=0; i<dim;i++){
		dist->closedorbit[i] = clo[i];
	}
}

void setdeltap_(double *dp){
	dist->emitt->dp = *dp;
	dist->longitunalemittance=2;
	//dist->emitt->e3 = 1000*pow(*dp/dist->tas[5][5],2);
}

void settasmatrix_(double tas[6][6]){
	for(int i =0; i< dim; i++){
		for(int j =0; j< dim; j++){
			dist->tas[i][j] = tas[i][j];
		}
	}
	printmatrix(6,6, dist->tas);
}

void settasmatrix_element(double element, int i, int j){
	dist->tas[i][j] = element;

	printmatrix(6,6, dist->tas);
}

void settasmatrixpython(double **tas){
		for(int i =0; i< dim; i++){
		for(int j =0; j< dim; j++){
			dist->tas[i][j] = tas[i][j];
		}
	}
	printmatrix(6,6, dist->tas);
}
void cutnormalized(int variable, double min, double max){

}

void setmassmom_(double *mass, double *momentum){
	(dist)-> mass = *mass;
	(dist)-> momentum = *momentum;
}

/*This function sets the parameter that is used to generate the distribution later
0 -  Constant value 
1 - Linearly spaced intervalls
2 - Exponentially spaced
3 - Spaced with  ^2
4 - Uniform random */ 
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
		if(dist->disttype==1){
			if(dist->totallength>0){
				int tmp;
				tmp=dist->totallength;
				length=&tmp;
				}
			else
				printf("You need to set a totallength for disttype 1!");
		}
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
	if(*type==4){ // uniform random 
		createLinearSpaced(*length, *start, *stop,dist->coord[*index-1]->values);
		for(int i=0;i <*length; i++){
			dist->coord[*index-1]->values[i] = rand_uni(*start, *stop);
			//printf("%f \n", dist->coord[*index-1]->values[i] );
		}
	}

	if(*type==5){ // Gaussian random (Here start is mean and stop is the standard deviation)
		createLinearSpaced(*length, *start, *stop,dist->coord[*index-1]->values);
		for(int i=0;i <*length; i++){
			dist->coord[*index-1]->values[i] = randn(*start, *stop);
		}
	}

	if(*type==6){ // Rayleigh distribution
		createLinearSpaced(*length, *start, *stop,dist->coord[*index-1]->values);
		for(int i=0;i <*length; i++){
			dist->coord[*index-1]->values[i] = randray(*start, *stop);

		}
	}
}

/*Create a TAS matrix using E-T, assuming uncoupled system */
void createtas0coupling_(double betax, double alfax, double betay, double alfay, double dx, double dpx, double dy, double dpy){

	    for (int i = 0; i < 6; i++)
    {
            for (int k = 0; k < 6; k++)
            {
                dist->tas [i][k] = 0;
                
            }
    }
    dist->tas[0][0] = sqrt(betax);
    dist->tas[1][0] =-(alfax)/sqrt(betax);
    dist->tas[1][1] =-1/sqrt(betax);
    dist->tas[2][2] = sqrt(betay);
    dist->tas[3][2] =-alfay/sqrt(betay);
    dist->tas[3][3] =-(1)/sqrt(betay);

    dist->tas[0][5] = dx;
    dist->tas[1][5] = dpx;

	dist->tas[2][5] = dy;
    dist->tas[3][5] = dpy;
    //This is not really physical so it should be changed.. 
    dist->tas[5][5] =1;
    dist->tas[4][4] =1;

    printmatrix(dim,dim, dist->tas);

}


void setphysicalcut(int variable, double min, double max){
	dist->cuts2apply->isset_p=1;
	dist->cuts2apply->physical[variable-1]->min=min;
	dist->cuts2apply->physical[variable-1]->max=max;
	dist->cuts2apply->physical[variable-1]->isset=1;

}

void setnormalizedcut(int variable, double min, double max){
	dist->cuts2apply->isset_n=1;
	dist->cuts2apply->normalized[variable-1]->min=min;
	dist->cuts2apply->normalized[variable-1]->max=max;
	dist->cuts2apply->normalized[variable-1]->isset=1;

}
