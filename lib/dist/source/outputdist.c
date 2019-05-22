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
	if(dist->totallength == -1){
		double length = 1;
		if(dist->disttype==0){
			for(int j =0; j<dim; j++){
			length = length*dist->coord[j]->length;
			}
		dist->totallength = length;
		}
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

void print2fort13(){

   double coord[6];
   double outputv[15];
   FILE *fptr;
   double x1, y1, xp1, yp1, energy1, deltas1, unit;
   double x2, y2, xp2, yp2, energy2, deltas2, energyref;
   unit = 1;

   fptr = fopen("fort.13", "w");
   if(fptr == NULL)
   {
      printf("Error!");
      exit(1);
   }
   for(int i=0;i < getnumberdist_()-1; i=i+2)
   {
   		getcoord_(coord, i);
   		outputv[0] = coord[0]*unit;
   		outputv[1] = coord[1]*unit;
   		outputv[2] = coord[2]*unit;
   		outputv[3] = coord[3]*unit;
   		outputv[4] = coord[4]*unit;	
   		outputv[5] = coord[5]*unit;
   	   	outputv[12] = momentum2energy(dist->momentum, dist->mass);
   	   	outputv[13] = momentum2energy((dist->momentum+(coord[5]*dist->momentum)), dist->mass);
   		getcoord_(coord, i+1);
   	   	outputv[0+6] = coord[0]*unit;
   		outputv[1+6] = coord[1]*unit;
   		outputv[2+6] = coord[2]*unit;
   		outputv[3+6] = coord[3]*unit;
   		outputv[4+6] = coord[4]*unit;
   		outputv[5+6] = coord[5]*unit;
   		outputv[14] = momentum2energy(dist->momentum+(coord[5]*dist->momentum), dist->mass);
   		
   		for(int j=0; j< 15; j++){
   			fprintf(fptr, "%f \n", outputv[j]);
   		}
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