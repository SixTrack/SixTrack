#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distinput.h"
#include "outputdist.h"

void 
canonical2six(double *canonical, double beta0, double pc0, double mass0, double mass, double *coord){

    double deltap = *(canonical+5);
    double beta = (pc0+deltap)/momentum2energy((pc0)+deltap, mass);
    double rv = beta0/beta;
    *(coord+0) = *(canonical+0);
    *(coord+1) = *(canonical+1)/(1+deltap);
    *(coord+2) = *(canonical+2);
    *(coord+3) = *(canonical+3)/(1+deltap);
    *(coord+4) = *(canonical+4)*rv;
    *(coord+5) = *(canonical+5);  
}


int splitline(char* line, char columns[100][100], char units[100][100]  ){

  char* word;
  int count = 0;
  word = strtok(line, " ");
  strcpy(columns[count], word);
  count++;
  char *unit = NULL;
  char *columnname = NULL;
  /* the following loop gets the rest of the words until the
   * end of the message */
  while ((word = strtok(NULL, " ")) != NULL){
    unit = strstr(word, "[");
    if(unit != NULL){
      
      strcpy(units[count], unit);
      columnname=strtok(word,"[");
      strcpy(columns[count], columnname);
    }
    else{
      strcpy(units[count],  "nounit");
      strcpy(columns[count], word);
    }
   
    count ++;
  }

  return count;
}

void add2table(double table[100][100], char* line, int linenum){
  char* word;
  int count = 0;
  word = strtok(line, " ");
  table[linenum][0] = atof(word);
  count++;

  /* the following loop gets the rest of the words until the
   * end of the message */
  while ((word = strtok(NULL, " ")) != NULL){
    table[linenum][count] = atof(word);
    count ++;
  }
  for(int i=0 ;i<count; i++){
    printf("nummmss %f %d \n ", table[linenum][i], i);
  }
}

  int readfile(){
   static const char filename[] = "file.txt";
   FILE *file = fopen ( filename, "r" );

    double table [100][100];
    int linecount = -1;
    char columns[100][100];
    char units[100][100];
    int numcolum ;
// NEED to use strcp... 
   if ( file != NULL ){
      char line [ 1000 ]; /* or other suitable maximum line size */
      char tosplit [ 1000 ];
      while ( fgets ( line, sizeof line, file ) != NULL ){ /* read a line */
        if(strncmp(line, "$", 1)==0){
          char value_s[20],  shorty[20];
          char * parameter =  (char*)malloc(20*sizeof(char)); 
          char * parameter_tmp =  (char*)malloc(20*sizeof(char)); 
          char * unit =(char*)malloc(20*sizeof(char));
          char * unit_tmp =(char*)malloc(20*sizeof(char));

          double value;
          double multifactor = 1;
          sscanf( line, "%s  %s",  parameter_tmp, value_s);

          
          unit_tmp =strstr(parameter_tmp, "[");         

          if(unit_tmp != NULL){
            strcpy(unit,unit_tmp);
            strcpy(parameter,strtok(parameter_tmp,"["));
  
          }
          else{

            strcpy(unit, "nounit");
            strcpy(parameter, parameter_tmp);
          }       
          
          value = atof(value_s);

          if(strcmp(parameter, "$energy0")==0){
            multifactor = getEnergyUnit(unit);
            dist->ref->e0=value*multifactor;
          }
          else if(strcmp(parameter, "$pc0")==0){
            multifactor = getEnergyUnit(unit);
            dist->ref->pc0=value*multifactor;
          }
          else if(strcmp(parameter, "$a0")==0){
            dist->ref->a0=(int)value;
          }
          else if(strcmp(parameter, "$z0")==0){
            dist->ref->z0=(int)value;
          }
          else if(strcmp(parameter, "$mass0")==0){
            multifactor = getEnergyUnit(unit);
            dist->ref->mass0=value*multifactor;
          }
          else if(strcmp(parameter, "$charge0")==0){
            dist->ref->charge0=(int)value;
          }
          else{
            printf("Not recogniced parameter %s \n", parameter);
          }
        }
        else if(strncmp(line, "#", 1)!=0 && linecount==-1){
          strcpy(tosplit, line);
          numcolum = splitline(tosplit, columns, units);
          linecount++; 
        }
        else if(linecount>=0){
          add2table(table, line, linecount);
          linecount++;   
        }

      }

      printf("line count %d \n",linecount);
      fclose ( file );  

 
   }
  
  else
  {
     perror ( filename ); /* why didn't the file open? */
  }

  allocateincoord(linecount);
  double multifactor;
  for(int i=0; i< numcolum; i++){
     if(strcmp(columns[i], "x")==0){
      multifactor = getMetricUnit(units[i]);
      setphysical(0, i, table);
    }
    else if(strcmp(columns[i], "px")==0)
      setphysical(1, i, table);
    else if(strcmp(columns[i], "y")==0){
      multifactor = getMetricUnit(units[i]);
      setphysical(2, i, table);
    }
    else if(strcmp(columns[i], "py")==0)
      setphysical(3, i, table);
    else if(strcmp(columns[i], "zeta")==0){
      multifactor = getMetricUnit(units[i]);
      setphysical(4, i, table);
    }
    else if(strcmp(columns[i], "deltap")==0){
      checkifenergyset(5);
      setphysical(5, i, table);
    }
    else if(strcmp(columns[i], "energy")==0){
      multifactor = getEnergyUnit(units[i]);
      checkifenergyset(0);
      setnonstandard(0, i, table);
    }
     else if(strcmp(columns[i], "pc")==0){
      multifactor = getEnergyUnit(units[i]);
      checkifenergyset(1);
      setnonstandard(1, i, table);
    }
    else if(strcmp(columns[i], "psigma")==0){
      checkifenergyset(2);
      setnonstandard(2, i, table);
    }
    else if(strcmp(columns[i], "ptau")==0){
      checkifenergyset(3);
      setnonstandard(3, i, table);
    }
    else if(strcmp(columns[i], "tau")==0){
      multifactor = getMetricUnit(units[i]);
      setnonstandard(4, i, table);
    }
    else if(strcmp(columns[i], "sigma")==0){
      multifactor = getMetricUnit(units[i]);
      setnonstandard(5, i, table);
    }
    else if(strcmp(columns[i], "xp")==0)
      setnonstandard(6, i, table);
    else if(strcmp(columns[i], "yp")==0)
      setnonstandard(7, i, table);
    else if(strcmp(columns[i], "jx")==0)
      setaction(0, i, table);
    else if(strcmp(columns[i], "phix")==0)
      setaction(1, i, table);
    else if(strcmp(columns[i], "jy")==0)
      setaction(2, i, table);
    else if(strcmp(columns[i], "phiy")==0)
      setaction(3, i, table);
    else if(strcmp(columns[i], "jz")==0)
      setaction(4, i, table);
    else if(strcmp(columns[i], "phiz")==0)
      setaction(5, i, table);
    else if(strcmp(columns[i], "xn")==0)
      setnormalized(0, i, table);
    else if(strcmp(columns[i], "pxn")==0)
      setnormalized(1, i, table);
    else if(strcmp(columns[i], "yn")==0)
      setnormalized(2, i, table);
    else if(strcmp(columns[i], "pyn")==0)
      setnormalized(3, i, table);
    else if(strcmp(columns[i], "zn")==0)
      setnormalized(4, i, table);
    else if(strcmp(columns[i], "pzn")==0)
      setnormalized(5, i, table);
    else if(strcmp(columns[i], "mass")==0){
      multifactor = getEnergyUnit(units[i]);
      for(int j=0;j < dist->totincoord; j++){
        dist->incoord[j]->mass = multifactor*table[j][i];
      }
    }
    else if(strcmp(columns[i], "a")==0){
      for(int j=0;j < dist->totincoord; j++){
        dist->incoord[j]->a = table[j][i];
      }
    }
    else if(strcmp(columns[i], "z")==0){
      for(int j=0;j < dist->totincoord; j++){
        dist->incoord[j]->z = table[j][i];
      }
    }

  }

  calculaterefparam();
  convert2standard();

   return 0;
}
void
issue_error(const char* t1){
  printf("+=+=+= fatal: %s\n",t1);
  exit(1);
}

void issue_warning(const char* t1){
  printf( "+=+=+= warning:  %s\n",t1);
}

void issue_info(const char* t1){
  printf( "Info: %s\n",t1);
}

double getEnergyUnit(char *unit){
  
  double multif;
  if(strcmp(unit, "[ev]")){
    multif = 1e9;
  }
  else if(strcmp(unit, "[kev]")){
    multif = 1e6;
  } 
  else if(strcmp(unit, "[mev]")){
    multif = 1e3;
  }
  else if(strcmp(unit, "[gev]") || strcmp(unit, "nounit") ){
    multif = 1;
  }
  else if(strcmp(unit, "[tev]")){
    multif = 1e-3;
  }
  else{
    printf("Unknow unit. Must be of the form [Mev]");
  } 
}


double getMetricUnit(char *unit){
  double multif;
  if(strcmp(unit, "[mm]")){
    multif = 1e3;
  }
  else if(strcmp(unit, "[m]") || strcmp(unit, "nounit") ){
    multif = 1;
  }
  else{
    issue_warning("Unknow unit. Must be of the form [Mev]");
  } 
}

void calculaterefparam(){
  if(dist->ref->mass0 == 0){
    issue_error("A mass is needed! \n");
  }
  if(dist->ref->pc0 == 0 && dist->ref->e0==0){
    issue_error("A energy is needed! \n");    
  }
  if(dist->ref->pc0 > 0 && dist->ref->e0 > 0){
    issue_error("Can't set both energy and momentum! \n");   
  }
  if(dist->ref->pc0 > 0){
   dist->ref->e0 = momentum2energy(dist->ref->pc0,dist->ref->mass0);
  }
  if(dist->ref->e0 > 0){
   dist->ref->pc0 = energy2momentum(dist->ref->e0,dist->ref->mass0);
  }

  dist->ref->beta0 = (dist->ref->pc0)/(dist->ref->e0);
  
  if(dist->incoord[0]->mass < 1e-16){ //enough to check the first 
    for(int i=0; i< dist->totincoord; i++){
      dist->incoord[0]->mass = dist->ref->mass0; //Gives all the same mass. 
    }
  }
}

void convert2standard(){

  if(dist->ref->en_like==-1){
    issue_info("No energy variable is set. Assume 0 deviation from reference energy \n");
  }
  else if(dist->ref->en_like==0){ //energy
    for(int i=0; i< dist->totincoord; i++){
      dist->incoord[i]->physical[5] = (momentum2energy(dist->incoord[i]->nonstandard[0], dist->incoord[i]->mass)-(dist->ref->pc0))/(dist->ref->pc0);
    }  
  }
  else if(dist->ref->en_like==1){ //momentum
    for(int i=0; i< dist->totincoord; i++){
      dist->incoord[i]->physical[5] = ((dist->incoord[i]->nonstandard[1], dist->incoord[i]->mass)-(dist->ref->pc0))/(dist->ref->pc0);
    }  
  }
  else if(dist->ref->en_like==2){ //psigma
    for(int i=0; i< dist->totincoord; i++){
      dist->incoord[i]->physical[5] = psigma2deltap(dist->incoord[i]->nonstandard[2], dist->ref->beta0);
    }  
  }
  else if(dist->ref->en_like==3){ //pt
    for(int i=0; i< dist->totincoord; i++){
      dist->incoord[i]->physical[5] = pt2deltap(dist->incoord[i]->nonstandard[3], dist->ref->beta0);
    }
  }
  else{
    issue_error("Something went wrong with setting the energy! ");
  }
}

double psigma2deltap(double psigma, double beta0 ){
  return (sqrt(pow(psigma*beta0,2) +2*psigma +1)-1);
}

double pt2deltap(double pt, double beta0 ){
  return (sqrt(pow(pt,2) +2*pt*beta0 +1)-1);
}


void checkifenergyset(int entype){
  if(dist->ref->en_like==-1)
    dist->ref->en_like=entype;
  else {
    issue_error("Only allowed 1 type of energy variable!");
    
  }
}

void print2filenew(){
  if(dist->ref->en_like ==-1){

  }
   FILE * fp;
   /* open the file for writing*/
   fp = fopen ("myoutfile.txt","w");
   fprintf (fp, "mass0 %f \n",dist->ref->mass0);
   fprintf (fp, "charge0 %d \n",dist->ref->charge0);
   fprintf (fp, "z0 %d \n",dist->ref->z0);
   fprintf (fp, "a0 %d \n",dist->ref->a0);
   fprintf (fp, "pc0 %f \n",dist->ref->pc0);

   fprintf(fp, "x px y py zeta deltap");

   /* close the file*/  
   fclose (fp);


}

void allocateincoord(int linecount){
  dist->incoord = (struct coordinates**)malloc(linecount*sizeof(struct coordinates*));
  dist->totincoord = linecount;
  for(int i=0; i<linecount; i++){
    dist->incoord[i] = (struct coordinates*)malloc(sizeof(struct coordinates));
    dist->incoord[i]->physical = (double*)malloc(dim*sizeof(double));
    dist->incoord[i]->normalized = (double*)malloc(dim*sizeof(double));
    dist->incoord[i]->action = (double*)malloc(dim*sizeof(double));
    dist->incoord[i]->nonstandard = (double*)malloc(9*sizeof(double));
  }
}


void setphysical(int coordorder, int column, double table[100][100]){

  for(int i=0;i < dist->totincoord; i++){
    dist->incoord[i]->physical[coordorder] = table[i][column];
  }

}
void setnormalized(int coordorder, int column, double table[100][100]){

  for(int i=0;i < dist->totincoord; i++){
    dist->incoord[i]->normalized[coordorder] = table[i][column];
  }
}
void setaction(int coordorder, int column, double table[100][100]){

  for(int i=0;i < dist->totincoord; i++){
    dist->incoord[i]->action[coordorder] = table[i][column];
  }
}

void setnonstandard(int coordorder, int column, double table[100][100]){
  for(int i=0;i < dist->totincoord; i++){
    dist->incoord[i]->nonstandard[coordorder] = table[i][column];
  }
}

double momentum2energy(double momentum, double mass){
    return sqrt(pow(momentum,2)+pow(mass,2));
}
double energy2momentum(double energy, double mass){
    return sqrt(pow(energy,2)-pow(mass,2));
}