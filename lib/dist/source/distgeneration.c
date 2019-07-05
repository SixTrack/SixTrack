#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distgeneration.h"
#include "distinput.h"
#include "outputdist.h"


/*At the moment this basically just copies (place holder for next version)*/
void gen2sixcoord(){
    int counter = 0;
    double tc[dim];
    double tmp[dim];
    double tmp_n[dim];

    for(int i =0; i< dist->totincoord; i++){
       
        tc[0]=dist->incoord[i]->physical[0];
        tc[1]=dist->incoord[i]->physical[1];
        tc[2]=dist->incoord[i]->physical[2];
        tc[3]=dist->incoord[i]->physical[3];
        tc[4]=dist->incoord[i]->physical[4];
        tc[5]=dist->incoord[i]->physical[5];
        
        if(particle_within_limits_physical(tmp)==1){
            for(int p=0; p<dim; p++){
                dist->outcoord[counter]->physical[p] = tc[p];
            }
            counter++;
        }
    }
    dist->totoutcoord=counter;
    dist->isDistrcalculated=1;
}
/*Checks if the particle is within the physical limit set by the user*/
int particle_within_limits_physical(double *physical){
    
    if(dist->cuts2apply->isset_p==0) return 1;
    for(int i=0; i<dim; i++){
        if(dist->cuts2apply->physical[i]->isset==1){
            if(physical[i] > dist->cuts2apply->physical[i]->min && physical[i] < dist->cuts2apply->physical[i]->max) return 0;
        }   
    }
    
    return 1;

}
/*Checks if the particle is within the normalized limit set by the user*/
int particle_within_limits_normalized(double *normalized){
    
    if(dist->cuts2apply->isset_n==0) return 1;
    for(int i=0; i<dim; i++){
        if(dist->cuts2apply->normalized[i]->isset==1){
            if(normalized[i] > dist->cuts2apply->normalized[i]->min && normalized[i] < dist->cuts2apply->normalized[i]->max) return 0;
        }
    }
    return 1;
}