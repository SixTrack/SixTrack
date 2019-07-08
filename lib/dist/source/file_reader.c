#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "helper.h"
#include "distinterface.h"
#include "outputdist.h"
#include "file_reader.h"
#include "distgeneration.h"
#include "round_near.c"

int readfile(const char*  filename){

  FILE *file = fopen ( filename, "r" );
  double** table = malloc(MAX_ROWS * sizeof(double*));    // allocate the rows
  for (int r = 0; r < MAX_ROWS; ++r){
     table[r] = malloc(MAX_COLUMNS * sizeof(double));    // allocate the columns
  }

  int linecount = -1;
  int numcolum, index;
  char columns[MAX_COLUMNS][MAX_LENGTH];
  char units[MAX_COLUMNS][MAX_LENGTH];
  char line [ MAX_LENGTH*50 ]; /* or other suitable maximum line size */
  char tosplit [ MAX_LENGTH*10  ];
  char value_s[MAX_LENGTH],  shorty[MAX_LENGTH];
  char * parameter =  (char*)malloc(MAX_LENGTH*sizeof(char)); 
  char * parameter_tmp =  (char*)malloc(MAX_LENGTH*sizeof(char)); 
  char * unit =(char*)malloc(MAX_LENGTH*sizeof(char));
  char * unit_tmp =(char*)malloc(MAX_LENGTH*sizeof(char));
  char * at =(char*)malloc(MAX_LENGTH*sizeof(char));
  double value;
  char *e ;

  if ( file != NULL ){
    while ( fgets ( line, sizeof line, file ) != NULL ){ /* read a line */
      double multifactor = 1;
      int k =0;
      if(strncmp(line, "@", 1)==0){
        while( line[k] ) {
          line[k] = (tolower(line[k]));
          k++;
        }
      sscanf( line, "%s %s %s",  at, parameter_tmp, value_s);
      if(strncmp(value_s, "%", 1)!=0){
        unit_tmp =strstr(parameter_tmp, "[");         

        if(unit_tmp != NULL){
          strcpy(unit,unit_tmp);
          e = strchr(parameter_tmp, '[');
          index = (int)(e - parameter_tmp);
          strncpy(parameter, parameter_tmp, index);
        }
        else{
          strcpy(unit, "nounit");
          strcpy(parameter, parameter_tmp);
        }       
        
        value = strtod_cr(value_s,NULL);

        if(strcmp(parameter, "energy0")==0){
          multifactor = getEnergyUnit(unit);
          dist->ref->e0=value*multifactor;
        }
        else if(strcmp(parameter, "pc0")==0){
          multifactor = getEnergyUnit(unit);
          dist->ref->pc0=value*multifactor;
        }
        else if(strcmp(parameter, "a0")==0){
          dist->ref->a0=round(value);
        }
        else if(strcmp(parameter, "z0")==0){
          dist->ref->z0=round(value);
        }
        else if(strcmp(parameter, "mass0")==0){
          multifactor = getEnergyUnit(unit);
          dist->ref->mass0=value*multifactor;
        }
        else if(strcmp(parameter, "charge0")==0){
          dist->ref->charge0=round(value);
        }
        else{
          issue_error2("Not recogniced parameter!", parameter);
        }
      }
    }
    else if(strncmp(line, "!", 1)==0){
          printf("Comment line %s is read but ignored by the dist lib. \n", line );
        }
        else if(strncmp(line, "!", 1)!=0 && linecount==-1){
          
          strcpy(tosplit, line);
          numcolum = splitline(tosplit, columns, units);
          linecount++; 
        }
        else if(linecount>=0 && strncmp(line, "$", 1)!=0){
          add2table(table, line, linecount);
          linecount++;   
        }
        else{
          issue_error2("Something went wrong with reading line", line);
        }

      }
      fclose ( file );
    }
    else
    {
       perror ( filename ); 
    }
    allocateincoord(linecount); //alocates the memory
    fillcoordstructure(numcolum,columns, units, table); //fils the coordinate structures
    calculaterefparam(); //Sets up the necessary ref parameters and checks for consistensy. 
    convert2standard(); // converts to the units used internal , x, px, y, py, zeta, deltap 
    free_readin_memory(table);
    return 0;
}

// This checks that reference energy is only set once and sets up the mass in case it is not a column
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
      dist->incoord[i]->mass = dist->ref->mass0; //Gives all the same mass. 
    }
  }
  if(dist->incoord[0]->a ==0){ //enough to check the first 
    for(int i=0; i< dist->totincoord; i++){
      dist->incoord[i]->a = dist->ref->a0; //Gives all the same mass. 
    }
  }
  if(abs(dist->incoord[0]->z) ==0){ //enough to check the first 
    for(int i=0; i< dist->totincoord; i++){
      dist->incoord[i]->z = dist->ref->z0; //Gives all the same mass. 
    }
  }
}

//Converts to the internal used canonical 
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

  if(dist->ref->time_like==-1){
    issue_info("No time variable is set. Assume 0 deviation \n");
  }
  else if(dist->ref->time_like==2){ //sigma
    double beta, pc;
    for(int i=0; i< dist->totincoord; i++){
      pc = dist->ref->pc0+dist->incoord[i]->physical[5];
      beta = (pc)/momentum2energy(pc, dist->incoord[i]->mass);
      dist->incoord[i]->physical[4] = sigma2zeta(dist->incoord[i]->nonstandard[4], dist->ref->beta0, beta);
    }  
  }
  else if(dist->ref->time_like==3){ //t or tau
    for(int i=0; i< dist->totincoord; i++){
      dist->incoord[i]->physical[4] = tau2zeta(dist->incoord[i]->nonstandard[5], dist->ref->beta0);
    }
  }
  if(dist->ref->ang_like==-1){
    issue_error("No angle variable is set. \n");
  }
  else if(dist->ref->ang_like==1){ //t or tau
    for(int i=0; i< dist->totincoord; i++){
      dist->incoord[i]->physical[1] = dist->incoord[i]->physical[6]/(1+dist->incoord[i]->physical[5]);
      dist->incoord[i]->physical[3] = dist->incoord[i]->physical[7]/(1+dist->incoord[i]->physical[5]);
    }
  }
}

void allocateincoord(int linecount){
  dist->incoord = (struct coordinates**)malloc(linecount*sizeof(struct coordinates*));
  dist->outcoord = (struct coordinates**)malloc(linecount*sizeof(struct coordinates*));
  dist->totincoord = linecount;
  for(int i=0; i<linecount; i++){
    dist->incoord[i] = (struct coordinates*)malloc(sizeof(struct coordinates));
    dist->incoord[i]->physical = (double*)malloc(dim*sizeof(double));
    dist->incoord[i]->normalized = (double*)malloc(dim*sizeof(double));
    dist->incoord[i]->action = (double*)malloc(dim*sizeof(double));
    dist->incoord[i]->nonstandard = (double*)malloc(9*sizeof(double));
    
    dist->incoord[i]->mass=0;
    dist->incoord[i]->a=0;
    dist->incoord[i]->z=0;

    dist->outcoord[i] = (struct coordinates*)malloc(sizeof(struct coordinates));
    dist->outcoord[i]->physical = (double*)malloc(dim*sizeof(double));
    dist->outcoord[i]->normalized = (double*)malloc(dim*sizeof(double));
    dist->outcoord[i]->action = (double*)malloc(dim*sizeof(double));
    dist->outcoord[i]->nonstandard = (double*)malloc(9*sizeof(double));
  }
}

void setphysical(int coordorder, int column, double ** table, double multifactor){

  for(int i=0;i < dist->totincoord; i++){
    dist->incoord[i]->physical[coordorder] = multifactor*table[i][column];
  }
}
void setnormalized(int coordorder, int column, double ** table){

  for(int i=0;i < dist->totincoord; i++){
    dist->incoord[i]->normalized[coordorder] = table[i][column];
  }
}
void setaction(int coordorder, int column, double ** table){

  for(int i=0;i < dist->totincoord; i++){
    dist->incoord[i]->action[coordorder] = table[i][column];
  }
}

void setnonstandard(int coordorder, int column, double ** table){
  for(int i=0;i < dist->totincoord; i++){
    dist->incoord[i]->nonstandard[coordorder] = table[i][column];
  }
}


int splitline(char* line, char columns[MAX_COLUMNS][MAX_LENGTH], char units[MAX_COLUMNS][MAX_LENGTH]  ){
  char* word;
  int count = 0;
  int k=0;

  while( line[k] ) {
    line[k] = (tolower(line[k]));
    k++;
  }

  word = strtok(line, " ");
  if(strncmp("*", word, 1)!=0){
  add2internaltab(word, columns, units, count);
  count++;
  }
  while ((word = strtok(NULL, " ")) != NULL){
    add2internaltab(word, columns, units, count);   
    count ++;
  }

  return count;
}
void add2internaltab(char* word, char columns[MAX_COLUMNS][MAX_LENGTH], char units[MAX_COLUMNS][MAX_LENGTH] , int count){
  

  char * unit =(char*)malloc(MAX_LENGTH*sizeof(char));
  char * word_tmp =(char*)malloc(MAX_LENGTH*sizeof(char));
  char * columnname =(char*)malloc(MAX_LENGTH*sizeof(char));
  char *e;
  int index;
  strcpy(word_tmp, word);
  unit =strstr(word_tmp, "[");
      if(unit != NULL){
      
      strcpy(units[count], unit);
      e = strchr(word_tmp, '[');
      index = (int)(e - word_tmp);
      strncpy(columns[count], word_tmp, index);
      columns[count][index]='\0';
    }
    else{
      strcpy(units[count],  "nounit");
      strcpy(columns[count], word_tmp);
    }
}

void add2table(double ** table, char* line, int linenum){
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
}

void fillcoordstructure(int numcolum, char columns[MAX_COLUMNS][MAX_LENGTH], char units[MAX_COLUMNS][MAX_LENGTH], double **table){

  double multifactor;
  for(int i=0; i< numcolum; i++){
     if(strcmp(columns[i], "x")==0){
      multifactor = getMetricUnit(units[i]);
      setphysical(0, i, table, multifactor);
    }
    else if(strcmp(columns[i], "px")==0){
      checkifangset(0);
      setphysical(1, i, table,1);
    }
    else if(strcmp(columns[i], "y")==0){
      multifactor = getMetricUnit(units[i]);
      setphysical(2, i, table, multifactor);
    }
    else if(strcmp(columns[i], "py")==0){
      checkifangset(0);
      setphysical(3, i, table,1);
    }
    else if(strcmp(columns[i], "zeta")==0){
      checkiftimeset(5);
      multifactor = getMetricUnit(units[i]);
      setphysical(4, i, table, multifactor);
    }
    else if(strcmpnl(columns[i], "deltap")==0){
      checkifenergyset(5);
      setphysical(5, i, table, 1);
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
    else if(strcmp(columns[i], "sigma")==0){
      checkiftimeset(2);
      multifactor = getMetricUnit(units[i]);
      setnonstandard(4, i, table);
    }
    else if(strcmp(columns[i], "tau")==0){
      checkiftimeset(3);
      multifactor = getMetricUnit(units[i]);
      setnonstandard(5, i, table);
    }
    else if(strcmp(columns[i], "xp")==0){
      checkifangset(1);
      setnonstandard(6, i, table);
    }
    else if(strcmp(columns[i], "yp")==0){
      checkifangset(1);
      setnonstandard(7, i, table);
    }
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
}
double getEnergyUnit(char *unit){
  double multif;
  if(strcmp(unit, "[ev]")==0){
    multif = 1e-9;
  }
  else if(strcmp(unit, "[kev]")==0){
    multif = 1e-6;
  }
  else if(strcmp(unit, "[mev]")==0){
    multif = 1e-3;
  }
  else if(strcmp(unit, "[gev]")==0 || strcmp(unit, "nounit")==0){
    multif = 1;
  }
  else if(strcmp(unit, "[tev]")==0){
    multif = 1e3;
  }
  else{
    issue_error("Unknow unit. Must be of the form [Mev]");
  }
  return multif;
}

double getMetricUnit(char *unit){
  double multif;

  if(strcmp(unit, "[mm]")==0){
    multif = 1e-3;
  }
  else if(strcmp(unit, "[m]")==0 || strcmp(unit, "nounit")==0 ){
    multif = 1;
  }
  else{
    issue_warning("Unknow unit. Must be of the form [Mev]");
  }
  return multif;
}

void free_readin_memory(double **table){

  for (int r = 0; r < MAX_ROWS; ++r){
    free(table[r]);    
  }
  free(table);    
}