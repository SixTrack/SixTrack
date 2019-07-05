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
void initializedistribution(int numberOfDist){
    dist = (struct distparam*)malloc((numberOfDist)*sizeof(struct distparam));
    dim  = 6;
    
        for(int i = 0; i <numberOfDist; i++)
        {
        
        (dist + i)->ref = (struct refparam*)malloc(sizeof(struct refparam));
        (dist + i)->coord = (struct parameters**)malloc(dim*sizeof(struct parameters*));
        (dist + i)->cuts2apply = (struct appliedcut*)malloc(sizeof(struct appliedcut));
        (dist + i)->cuts2apply->physical = (struct cut**)malloc(dim*sizeof(struct cut*));
        (dist + i)->cuts2apply->normalized = (struct cut**)malloc(dim*sizeof(struct cut*));
        (dist + i)->tas   = (double**)malloc(dim*sizeof(double*));
        (dist + i)->invtas   = (double**)malloc(dim*sizeof(double*));
        (dist + i)->closedorbit   = (double*)malloc(dim*sizeof(double));
        (dist + i)->isDistrcalculated = 0;
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
        (dist + i)->coordtype   =-1;
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

void setdistribution(int ndist){
		dist = diststart + ndist;
}

void addclosedorbit(double *clo){
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

void get6trackcoord(double *x, double *xp, double *y, double *yp, double *sigma, double *deltap, int *totparticles){
    double tmp[6];
    if(dist->isDistrcalculated ==0){
        gen2sixcoord();
    }

    for(int i=0; i<dist->totoutcoord; i++){
        canonical2six(dist->incoord[i]->physical, dist->ref->beta0, dist->ref->pc0, dist->ref->mass0, dist->incoord[i]->mass, tmp);
        x[i]  = tmp[0];
        xp[i] = tmp[1];
        y[i]  = tmp[2];
        yp[i] = tmp[3];
        sigma[i]  = tmp[4];
        deltap[i] = tmp[5];
   }
    *totparticles=dist->totoutcoord;

}