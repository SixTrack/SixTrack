#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distgeneration.h"
#include "distinput.h"
#include "outputdist.h"

/*This is to create a grid*/
void dist2sixcoord_(){
    int counter = 0;
    double tc[6];
    double tmp[6];
    double tmp_n[6];
    dist->distout = (double**)malloc(getnumberdist_()*sizeof(double*));
    dist->distout_normalized = (double**)malloc(getnumberdist_()*sizeof(double*));
    for(int i =0; i< dist->coord[0]->length; i++){
        for(int j =0; j< dist->coord[1]->length; j++){
            for(int k =0; k< dist->coord[2]->length; k++){
                for(int l =0; l< dist->coord[3]->length; l++){
                    for(int m =0; m< dist->coord[4]->length; m++){
                        for(int n =0; n< dist->coord[5]->length; n++){
                            dist->distout[counter] = (double*)malloc(dim*sizeof(double));
                            dist->distout_normalized[counter] = (double*)malloc(dim*sizeof(double));
                            tc[0]=dist->coord[0]->values[i];
                            tc[1]=dist->coord[1]->values[j];
                            tc[2]=dist->coord[2]->values[k];
                            tc[3]=dist->coord[3]->values[l];
                            tc[4]=dist->coord[4]->values[m];
                            tc[5]=dist->coord[5]->values[n];
                           //  printf("aaaa %f \n", tc[0]);
                           // printvector("tcaa", 6, tc);
                         
                            action2sixinternal_(tc, tmp, tmp_n);
                            
                            if(particle_within_limits_physical(tmp)==1 && particle_within_limits_normalized(tmp_n)==1){
                                for(int p=0; p<dim; p++){
                                    dist->distout[counter][p] = tmp[p]+dist->closedorbit[p];
                                    dist->distout_normalized[counter][p] = tmp_n[p];
                                }
                                counter++;
                            }
                                                        
                        }
                    }
                }
            }
        }   
    }
      printf("end \n");
    dist->totallength=counter;
    dist->isDistrcalculated=1;
}
/*This is to create a random distribution*/

void createrandomdist_(){
    int counter = 0;
    double tc[6];
    double tmp[6];
    double tmp_normalized[6];
    int tempAlloc = dist->totallength;
    dist->distout = (double**)malloc(tempAlloc*sizeof(double*));
    dist->distout_normalized = (double**)malloc(tempAlloc*sizeof(double*));
    for(int i =0; i< tempAlloc; i++){
        for(int j =0; j< 6; j++){
            if(dist->coord[j]->type == 0){
                tc[j] = dist->coord[j]->values[0];
            }
            else if(dist->coord[j]->type==4 || dist->coord[j]->type==6){
                tc[j] = dist->coord[j]->values[i];
            }
            
        }
         dist->distout[i] = (double*)malloc(dim*sizeof(double));
         dist->distout_normalized[i] = (double*)malloc(dim*sizeof(double));
         action2sixinternal_(tc, tmp, tmp_normalized);
         for(int p=0; p<dim; p++){

                dist->distout[i][p] = tmp[p]+dist->closedorbit[p];
                dist->distout_normalized[i][p] = tmp_normalized[p];
        }
    }
    dist->isDistrcalculated=1;
}

void action2sixinternal_(double tc[6], double *results, double *tmp_normalized){
    double cancord[6];
    if(checkdist)
    {

        action2canonical_(tc, cancord, tmp_normalized);
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

void solve2by2eq(double a1, double b1, double c1, double a2, double b2, double c2, double *x){
    double det = a1*b2-b1*a2;
    double dx = c1*b2-b1*c2;
    double dy = a1*c2-c1*a2;
    x[0] = (dx/det);
    x[1] = (dy/det);

//    printf("this is the solution, x ,y %f %f", x[0], x[1] );

}


/*If emittance is defined it converts to canonical coordinates */
void action2canonical_(double acangl[6], double cancord[6], double acoord[6]){
    
    double dp_setting;

    acoord[0]= sqrt((dist->emitt->e1)*acangl[0]/2)*cos(acangl[1]);
    acoord[1]=-sqrt((dist->emitt->e1)*acangl[0]/2)*sin(acangl[1]);
    acoord[2]= sqrt((dist->emitt->e2)*acangl[2]/2)*cos(acangl[3]);
    acoord[3]=-sqrt((dist->emitt->e2)*acangl[2]/2)*sin(acangl[3]);
    acoord[4]= sqrt((dist->emitt->e3)*acangl[4]/2)*cos(acangl[5]);
    acoord[5]=-sqrt((dist->emitt->e3)*acangl[4]/2)/1000*sin(acangl[5]);



    if(dist->longitunalemittance==2) {
        //printf("dist22a %f \n", acangl[5]);
        //acoord[5] = acangl[5];
        double lindp = 0;
        double lindeltas=0;
        double deltap =acangl[4];
        double deltas = acangl[5];
        double *xap;
        double det = (dist->tas[4][4]*dist->tas[5][5] - dist->tas[4][5]*dist->tas[5][4]);
        for(int i=0; i<4;i++){
            lindeltas = lindeltas+dist->tas[4][i]*acoord[i];
            lindp=lindp+dist->tas[5][i]*acoord[i];
        } 

        printf("suuumm %f %f", lindeltas, lindp);
        lindp = deltap - lindp;
        lindeltas = deltas - lindeltas;

        xap = (double*)malloc(2*sizeof(double));
        solve2by2eq(dist->tas[4][4], dist->tas[4][5], lindeltas, dist->tas[5][4], dist->tas[5][5], lindp, xap );
        acoord[4] = xap[0];
        acoord[5] = xap[1];

        printf("uuua %f %f  \n",acoord[4], acoord[5] );

        mtrx_vector_mult_pointer(dim,dim, dist->tas, acoord,cancord);
        //printf("dist22b %f \n", cancord[5]);
        
        //if(fabs(cancord[5]-acoord[5])> 1e-10){
        //change_e3_to_dp_easy(cancord,acoord, acangl);

        //}
    }

    //This is the multiplication with the tas matrix 
    mtrx_vector_mult_pointer(dim,dim, dist->tas, acoord,cancord);
     
}

void change_e3_to_dp_easy(double cancord[6], double acoord[6], double acangl[6]){
 
    
    int itr = 0, maxmitr;
    double x, a, b, allerr, x1;
    double angle, emitt;
    a = -1.6;
    b=1.6;
    itr = 0;
    maxmitr = 100;
    allerr = 1e-12;
    bisection (&x, a, b, &itr);
    do
    {
        if (optideltas(cancord,acoord, acangl, a)*optideltas(cancord,acoord, acangl, x) < 0)
            b=x;
        else
            a=x;
        bisection (&x1, a, b, &itr);
        if (fabs(x1-x) < allerr)
        {
            printf("After %d iterations, root = %6.4f\n", itr, x1);
            break;
        }
        x=x1;
    }
    while (itr < maxmitr);


    toactioncord_(cancord,acoord,acangl);
        //printf("oooutt %E, %E, %E \n", cancord[4], cancord[5], acangl[5] );

    emitt =dist->emitt->e3;
    a = 0;
    b= 10000;
    itr = 0;
    

        bisection (&x, a, b, &itr);
    do
    {
        if (opemitt(cancord,acoord, acangl, a)*opemitt(cancord,acoord, acangl, x) < 0)
            b=x;
        else
            a=x;
        bisection (&x1, a, b, &itr);
        if (fabs(x1-x) < allerr)
        {
            printf("After %d iterations, root = %6.4f\n", itr, x1);
            break;
        }
        x=x1;
    }
    while (itr < maxmitr);

    toactioncord_(cancord,acoord,acangl);
}

double opemitt(double cancord[6], double acoord[6], double acangl[6], double x)
{
        dist->emitt->e3 = x;
        toactioncord_(cancord,acoord,acangl);
        return cancord[5]-dist->emitt->dp*acangl[4];
}

double createrandom(double insigma[6], double cancord[6], double normalized[6]){
    double acangl[6];
    for(int i=0; i<6; i++){
        acangl[i] = randn(0, insigma[i]);
    }
    action2canonical_(acangl, cancord, normalized);
}


double optideltas(double cancord[6], double acoord[6], double acangl[6], double x)
{
        acangl[5] = x;
        toactioncord_(cancord,acoord,acangl);
        return cancord[4];
}

double toactioncord_(double cancord[6], double acoord[6], double acangl[6])
{

        acoord[0]= sqrt((dist->emitt->e1)*acangl[0]/2)*cos(acangl[1]);
        acoord[1]=-sqrt((dist->emitt->e1)*acangl[0]/2)*sin(acangl[1]);
        acoord[2]= sqrt((dist->emitt->e2)*acangl[2]/2)*cos(acangl[3]);
        acoord[3]=-sqrt((dist->emitt->e2)*acangl[2]/2)*sin(acangl[3]);
        acoord[4]= sqrt((dist->emitt->e3)*acangl[4]/2)*cos(acangl[5]);
        acoord[5]=-sqrt((dist->emitt->e3)*acangl[4]/1000/2)*sin(acangl[5]);
        mtrx_vector_mult_pointer(dim,dim, dist->tas, acoord,cancord);

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