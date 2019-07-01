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
        (dist + i)->ref = (struct refparam*)malloc(sizeof(struct refparam));
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
		(dist + i)->ref->e0=0;
		(dist + i)->ref->pc0=0;
		(dist + i)->ref->a0=1;
		(dist + i)->ref->z0=1;
		(dist + i)->ref->mass0=0;
		(dist + i)->ref->charge0=1;
		(dist + i)->ref->en_like=-1;



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



void setdisttype(int *disttype){
	dist->disttype=*disttype;
}


void setdistribution_(int *ndist){
		dist = diststart + *ndist;
}

void addclosedorbit_(double *clo){
	for(int i=0; i<dim;i++){
		dist->closedorbit[i] = clo[i];
	}
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