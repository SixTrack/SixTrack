#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distinput.h"
#include "distgeneration.h"





/*If emittance is defined it converts to canonical coordinates */
void action2canonical_(double acangl[6], double cancord[6]){
	double acoord[6];
	double dp_setting;
	acoord[0]= sqrt((dist->emitt->e1)*acangl[0])*cos(acangl[1]);
	acoord[1]=-sqrt((dist->emitt->e1)*acangl[0])*sin(acangl[1]);
	acoord[2]= sqrt((dist->emitt->e2)*acangl[2])*cos(acangl[3]);
	acoord[3]=-sqrt((dist->emitt->e2)*acangl[2])*sin(acangl[3]);
	acoord[4]= sqrt((dist->emitt->e3)*acangl[4])*cos(acangl[5]);
	acoord[5]=-sqrt((dist->emitt->e3)*acangl[4]/1000)*sin(acangl[5]);



    if(dist->longitunalemittance==2) {
    	dp_setting = acangl[4];
    	change_e3_to_dp_easy(cancord,acoord, acangl);
    	change_e3_to_dp_easy(cancord,acoord, acangl);
    	change_e3_to_dp_easy(cancord,acoord, acangl);
    }

	//This is the multiplication with the tas matrix 
	mtrx_vector_mult_pointer(dim,dim, dist->tas, acoord,cancord);
}



double createrandom(double insigma[6], double cancord[6]){
	double acangl[6];
	for(int i=0; i<6; i++){
		acangl[i] = randn(0, insigma[i]);
	}
	action2canonical_(acangl, cancord);
}


double optideltas(double cancord[6], double acoord[6], double acangl[6], double x)
{
	 	acangl[5] = x;
		toactioncord_(cancord,acoord,acangl);
		return cancord[4];
}

double toactioncord_(double cancord[6], double acoord[6], double acangl[6])
{

    	acoord[0]= sqrt((dist->emitt->e1)*acangl[0])*cos(acangl[1]);
		acoord[1]=-sqrt((dist->emitt->e1)*acangl[0])*sin(acangl[1]);
		acoord[2]= sqrt((dist->emitt->e2)*acangl[2])*cos(acangl[3]);
		acoord[3]=-sqrt((dist->emitt->e2)*acangl[2])*sin(acangl[3]);
		acoord[4]= sqrt((dist->emitt->e3)*acangl[4])*cos(acangl[5]);
		acoord[5]=-sqrt((dist->emitt->e3)*acangl[4]/1000)*sin(acangl[5]);
		mtrx_vector_mult_pointer(dim,dim, dist->tas, acoord,cancord);

}


double opemitt(double cancord[6], double acoord[6], double acangl[6], double x)
{
	 	dist->emitt->e3 = x;
		toactioncord_(cancord,acoord,acangl);
		return cancord[5]-dist->emitt->dp*acangl[4];
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

void setdistribution_(int *ndist){
		dist = diststart + *ndist;
}

void setemittance12_(double *e1, double *e2){
	dist->emitt->e1=*e1; 
	dist->emitt->e2=*e2; 

}

void addclosedorbit_(double *clo){
	for(int i=0; i<dim;i++){
		dist->closedorbit[i] = clo[i];
	}

}


void setemittance3_(double *e3){
	dist->emitt->e3=*e3;  	
}


void setdeltap_(double *dp){
	dist->emitt->dp = *dp;
	dist->longitunalemittance=2;
	dist->emitt->e3 = 1000*pow(*dp/dist->tas[5][5],2);
}

//This emittance is oversimplified but gives a good approximation. 
//void convertdp2emittance(double dp){
//	calcualteinverse();
	//dist->emitt->e3 = pow(1000*dp*dist->invtas[5][5],2);
//	dist->emitt->e3 = 1000*pow(dp/dist->tas[5][5],2);
//
//}



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

void getcoordvectors_(double *x, double *xp, double *y, double *yp, double *sigma, double *delta){
	int distlen = getnumberdist_();
	if(dist->isDistrcalculated==0){
		dist2sixcoord_();
	}
	for(int i=0; i<distlen;i++){
		x[i] = dist->distout[i][0];
		xp[i] = dist->distout[i][1];
		y[i] = dist->distout[i][2];
		yp[i] = dist->distout[i][3];
		sigma[i] = dist->distout[i][4];
		delta[i] = dist->distout[i][5];
	}
}



void getcoord_(double *coordinate, int *initial ){

	if(dist->isDistrcalculated==0){
		dist2sixcoord_();

	}
	if(*initial >= getnumberdist_()){
		printf("Not generated, total inital coordinates generated is %f:",getnumberdist_() );
	for(int i=0; i<dim; i++){
		coordinate[i] = NAN;
	}
		return;
	}
	for(int i=0; i<dim; i++){
		coordinate[i] = dist->distout[*initial][i];
	}
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

void setphysicalcut(int variable, double min, double max){
	dist->cuts2apply->isset_p=1;
	dist->cuts2apply->physical[variable-1]->min=min;
	dist->cuts2apply->physical[variable-1]->max=max;
	dist->cuts2apply->physical[variable-1]->isset=1;

}

void setnormalizedcut(int variable, double min, double max){
	dist->cuts2apply->isset_p=1;
	dist->cuts2apply->normalized[variable-1]->min=min;
	dist->cuts2apply->normalized[variable-1]->max=max;
	dist->cuts2apply->normalized[variable-1]->isset=1;

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

    printmatrix(dim,dim, dist->tas);

}


