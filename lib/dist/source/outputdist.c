#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distinput.h"
#include "distgeneration.h"
#include "outputdist.h"


void print2file(const char* nameoffile){

   FILE * fp;
   /* open the file for writing*/
   fp = fopen ("myoutfile.txt","w");
   fprintf (fp, "@ mass0 %f \n",dist->ref->mass0);
   fprintf (fp, "@ charge0 %d \n",dist->ref->charge0);
   fprintf (fp, "@ z0 %d \n",dist->ref->z0);
   fprintf (fp, "@ a0 %d \n",dist->ref->a0);
   fprintf (fp, "@ pc0 %f \n",dist->ref->pc0);

   fprintf(fp, "x px y py zeta deltap \n");
   gen2sixcoord();
   for(int i=0; i<dist->totoutcoord; i++){
   	fprintf(fp, "%.9e %.9e %.9e %.9e %.9e %.9e \n", dist->outcoord[i]->physical[0],dist->outcoord[i]->physical[1],dist->outcoord[i]->physical[2],
   		dist->outcoord[i]->physical[3],dist->outcoord[i]->physical[4],dist->outcoord[i]->physical[5]);
   }

   /* close the file*/  
   fclose (fp);
}