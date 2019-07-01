#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distgeneration.h"
#include "distinput.h"
#include "outputdist.h"



void gen2sixcoord(){
    int counter = 0;
    double tc[6];
    double tmp[6];
    double tmp_n[6];
    dist->distout = (double**)malloc(dist->totincoord*sizeof(double*));
    dist->distout_normalized = (double**)malloc(dist->totincoord*sizeof(double*));
    for(int i =0; i< dist->totincoord; i++){

        dist->distout[counter] = (double*)malloc(dim*sizeof(double));
        dist->distout_normalized[counter] = (double*)malloc(dim*sizeof(double));
        
        tc[0]=dist->incoord[i]->physical[0];
        tc[1]=dist->incoord[i]->physical[1];
        tc[2]=dist->incoord[i]->physical[2];
        tc[3]=dist->incoord[i]->physical[3];
        tc[4]=dist->incoord[i]->physical[4];
        tc[5]=dist->incoord[i]->physical[5];
     
        //action2sixinternal_(tc, tmp, tmp_n);

        canonical2six(tc, dist->ref->beta0, dist->ref->pc0, dist->ref->mass0, dist->incoord[i]->mass, tmp);
        if(particle_within_limits_physical(tmp)==1){
            for(int p=0; p<dim; p++){
                dist->distout[counter][p] = tmp[p]+dist->closedorbit[p];
                //dist->distout_normalized[counter][p] = tmp_n[p];
            }
            counter++;
        }
                                                                   
         
    }
    dist->totallength=counter;
    dist->isDistrcalculated=1;
}



int particle_within_limits_physical(double *physical){
    
    if(dist->cuts2apply->isset_p==0) return 1;
    for(int i=0; i<dim; i++){
        if(dist->cuts2apply->physical[i]->isset==1){
            if(physical[i] > dist->cuts2apply->physical[i]->min && physical[i] < dist->cuts2apply->physical[i]->max) return 0;
        }   
    }
    
    return 1;

}

int particle_within_limits_normalized(double *normalized){
    
    if(dist->cuts2apply->isset_n==0) return 1;
    for(int i=0; i<dim; i++){
        if(dist->cuts2apply->normalized[i]->isset==1){
            if(normalized[i] > dist->cuts2apply->normalized[i]->min && normalized[i] < dist->cuts2apply->normalized[i]->max) return 0;
        }
    }
    return 1;
}