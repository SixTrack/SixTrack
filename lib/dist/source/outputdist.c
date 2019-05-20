#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distinput.h"
#include "distgeneration.h"
#include "outputdist.h"

void printdistsettings_(int *ndist){
	printf("Printing info about distribution: \n");
	printf("Coordianate type (input), 1=normalized coordinates %d \n", dist->coordtype );
	printf("Distribution generation type 0-Square like generation, 1-loop over every step once,  %d \n", dist->disttype );
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

/*
Returns the total number of particles that will be created for the distribution.  
*/
int getnumberdist_(){
	double length = 1;
	for(int j =0; j<dim; j++){
		length = length*dist->coord[j]->length;
	}
	dist->totallength = length;
	return length;
}

void getcoordvectors_(double *x, double *xp, double *y, double *yp, double *sigma, double *delta){
	int distlen = getnumberdist_();
	if(dist->isDistrcalculated==0){
		if(dist->disttype==0)
			dist2sixcoord_();
		else if(dist->disttype==1)
			createrandomdist_();
		else
			printf("Not a valid distribution type");
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

void getcoord_(double coordinate[6], int initial ){

	if(dist->isDistrcalculated==0){
		if(dist->disttype==0)
			dist2sixcoord_();
		else if(dist->disttype==1)
			createrandomdist_();
		else
			printf("Not a valid distribution type");
	}

	if(initial >= dist->totallength){
		printf("Not generated, total inital coordinates generated is %d :",getnumberdist_() );
	for(int i=0; i<dim; i++){
		coordinate[i] = NAN;
	}
		return;
	}

	for(int i=0; i<dim; i++){

		coordinate[i] = dist->distout[initial][i];
	}
}