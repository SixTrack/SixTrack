#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distr.h"


struct distparam* dist;
struct distparam* diststart;
int dim;
int distn=0;

void canonical2emittance_(double cancord[6], double emittance[3]){
	double accord[6];
	calcualteinverse();
	mtrx_vector_mult_pointer(dim,dim, dist->invtas, cancord, accord);

	emittance[0] = pow(accord[0],2)+pow(accord[1],2);
	emittance[1] = pow(accord[2],2)+pow(accord[3],2);
	emittance[2] = pow(accord[4],2)+pow(accord[5],2);

}

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

void multi_w_emittance(double acangl[6], double cancord[6]){
	double acoord[6];
	acoord[0]= sqrt((dist->emitt->e1))*(acangl[0]);
	acoord[1]= sqrt((dist->emitt->e1))*(acangl[1]);
	acoord[2]= sqrt((dist->emitt->e2))*(acangl[2]);
	acoord[3]= sqrt((dist->emitt->e2))*(acangl[3]);
	acoord[4]= sqrt((dist->emitt->e3))*(acangl[4]);
	acoord[5]= sqrt((dist->emitt->e3)/1000)*(acangl[5]);
	mtrx_vector_mult_pointer(dim,dim, dist->tas, acoord,cancord);

}

double createrandom(double insigma[6], double cancord[6]){
	double acangl[6];
	for(int i=0; i<6; i++){
		acangl[i] = randn(0, insigma[i]);
	}
	multi_w_emittance(acangl, cancord);
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



void setdistribution_(int *ndist){
		dist = diststart + *ndist;
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

void calcualteinverse(){
	double invtas[6][6];
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

}

void addclosedorbit_(double *clo){
	for(int i=0; i<dim;i++){
		dist->closedorbit[i] = clo[i];
	}

}

void setdeltap_(double *dp){
	dist->emitt->dp = *dp;
	dist->longitunalemittance=2;
	dist->emitt->e3 = 1000*pow(*dp/dist->tas[5][5],2);
}

//This emittance is oversimplified but gives a good approximation. 
void convertdp2emittance(double dp){
	calcualteinverse();
	//dist->emitt->e3 = pow(1000*dp*dist->invtas[5][5],2);
	dist->emitt->e3 = 1000*pow(dp/dist->tas[5][5],2);

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

