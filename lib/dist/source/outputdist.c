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
	if(dist->disttype==0){
		for(int j =0; j<dim; j++){
		length = length*dist->coord[j]->length;
		}
	dist->totallength = length;
	}
	return dist->totallength;
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
void getcoord_normalized(double coordinate_normalized[6], int initial ){
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
		coordinate_normalized[i] = NAN;
	}
		return;
	}

	for(int i=0; i<dim; i++){

		coordinate_normalized[i] = dist->distout_normalized[initial][i];
	}

}
void print2file(){
   double coord[6];
   FILE *fptr;
   double x, y, xp, yp, momentum, deltas, unit;
   unit = 10e-3;

   fptr = fopen("coordinates.out", "w");
   if(fptr == NULL)
   {
      printf("Error!");
      exit(1);
   }
   for(int i=0;i < getnumberdist_(); i++)
   {
   		getcoord_(coord, i);
   		x = coord[0]*unit;
   		y = coord[2]*unit;
   		xp = coord[1]*unit;
   		yp = coord[3]*unit;
   		momentum = dist->momentum+(coord[5]*dist->momentum);
   		deltas = coord[4]*unit;
        fprintf(fptr,"%d  %d  %d  %f  %f  %d  %f  %f  %d  %d %d  %f  %f  %f \n", 
                      i,   i,  i,  x,  y,  i, xp, yp,  i, dist->massnum, dist->atomnum, dist->mass, momentum, deltas   );
   }
   fclose(fptr);
}

/*
1	particle id (not actually used)
2	parent particle id (not actually used)
3	statistical weight (not actually used)
4	x[m]
5	y [m]
6	(unused)
7	x[]
8	y′[]
9	(unused)
10	mass number
11	atomic number
12	mass [GeV/c]
13	linear momentum [GeV/c]
14	time lag [s]
*/