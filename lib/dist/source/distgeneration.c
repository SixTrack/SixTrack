#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distinput.h"
#include "distgeneration.h"







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
        (dist + i)->emitt->e1 = 0; 
        (dist + i)->emitt->e2 = 0; 
        (dist + i)->emitt->e3 = 0;
        (dist + i)->emitt->dp = 0;
        (dist + i)->emitt->deltas = 0;
        (dist + i)->coordtype =-1;
        
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
/*Not fully sure what this is usefull for */
void dist2sixcoord_(){
    int counter = 0;
    double tc[6];
    double tmp[6];
    dist->distout = (double**)malloc(getnumberdist_()*sizeof(double*));
    for(int i =0; i< dist->coord[0]->length; i++){
        for(int j =0; j< dist->coord[1]->length; j++){
            for(int k =0; k< dist->coord[2]->length; k++){
                for(int l =0; l< dist->coord[3]->length; l++){
                    for(int m =0; m< dist->coord[4]->length; m++){
                        for(int n =0; n< dist->coord[5]->length; n++){
                            dist->distout[counter] = (double*)malloc(dim*sizeof(double));

                            action2sixinternal_(tc, tmp);
                            
                            if(particle_within_limits_physical(tmp)==1){
                                for(int p=0; p<dim; p++){
                                    dist->distout[counter][p] = tmp[p]+dist->closedorbit[p];
                                }
                                counter++;
                            }
                                                        
                        }
                    }
                }
            }
        }   
    }
    dist->isDistrcalculated=1;
}

void action2sixinternal_(double tc[6], double *results){
    double cancord[6];
    if(checkdist)
    {
        action2canonical_(tc, cancord);
        canonical2six_(cancord, &dist->momentum, &dist->mass, results);
    }
}

/*Checks that the necessary parameters have been set*/
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

/* This is the inverse of what is normally done*/
void canonical2emittance_(double cancord[6], double emittance[3]){
    double accord[6];
    calcualteinverse();
    mtrx_vector_mult_pointer(dim,dim, dist->invtas, cancord, accord);
    emittance[0] = pow(accord[0],2)+pow(accord[1],2);
    emittance[1] = pow(accord[2],2)+pow(accord[3],2);
    emittance[2] = pow(accord[4],2)+pow(accord[5],2);

}



