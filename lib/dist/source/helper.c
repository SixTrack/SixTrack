#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distinput.h"
#include "outputdist.h"




int splitline(char* line, char columns[100][100] ){

  char* word;
  int count = 0;
  word = strtok(line, " ");
  strcpy(columns[count], word);
  count++;
  /* the following loop gets the rest of the words until the
   * end of the message */
  while ((word = strtok(NULL, " ")) != NULL){
    strcpy(columns[count], word);
    count ++;
  }
  for(int i=0 ;i<count; i++){
    printf("columns %s \n ", columns[i]);
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
    int numcolum ;
// NEED to use strcp... 
   if ( file != NULL ){
      char line [ 1000 ]; /* or other suitable maximum line size */
      char tosplit [ 1000 ];
      while ( fgets ( line, sizeof line, file ) != NULL ){ /* read a line */
        if(strncmp(line, "%", 1)==0){
          char value_s[20], parameter[20], shorty[20];
          double value;
          sscanf( line, "%s %s  %s", shorty, parameter, value_s);
          value = atof(value_s);

          if(strcmp(parameter, "energy0")==0){
            dist->ref->e0=value;
          }
          else if(strcmp(parameter, "pc0")==0){
            dist->ref->pc0=value;
          }
          else if(strcmp(parameter, "a0")==0){
            dist->ref->a0=(int)value;
          }
          else if(strcmp(parameter, "z0")==0){
            dist->ref->z0=(int)value;
          }
          else if(strcmp(parameter, "mass0")==0){
            dist->ref->mass0=value;

          }
          else if(strcmp(parameter, "charge0")==0){
            dist->ref->charge0=(int)value;
          }
          else{
            printf("Not recogniced parameter %s \n", parameter);
          }
        }
        else if(strncmp(line, "#", 1)!=0 && linecount==-1){
          strcpy(tosplit, line);
          numcolum = splitline(tosplit, columns);
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
  int ab = 1;


  allocateincoord(linecount);

  for(int i=0; i< numcolum; i++){
     if(strcmp(columns[i], "x")==0)
      setphysical(0, i, table);
    else if(strcmp(columns[i], "px")==0)
      setphysical(1, i, table);
    else if(strcmp(columns[i], "y")==0)
      setphysical(2, i, table);
    else if(strcmp(columns[i], "py")==0)
      setphysical(3, i, table);
    else if(strcmp(columns[i], "zeta")==0)
      setphysical(4, i, table);
    else if(strcmp(columns[i], "deltap")==0){
      checkifenergyset(5);
      setphysical(5, i, table);
    }
    else if(strcmp(columns[i], "energy")==0){
      checkifenergyset(0);
      setnonstandard(0, i, table);
    }
     else if(strcmp(columns[i], "pc")==0){
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
    else if(strcmp(columns[i], "tau")==0)
      setnonstandard(4, i, table);
    else if(strcmp(columns[i], "sigma")==0)
      setnonstandard(5, i, table);
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
      for(int j=0;j < dist->totincoord; j++){
        dist->mass = table[j][i];
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

void calculaterefparam(){
  if(dist->ref->mass0 == 0){
    printf("A mass is needed! \n");
    exit(1);
  }
  if(dist->ref->pc0 == 0 && dist->ref->e0==0){
    printf("A energy is needed! \n");
    exit(1);    
  }
  if(dist->ref->pc0 > 0 && dist->ref->e0 > 0){
    printf("Can't set both energy and momentum! \n");
    exit(1);    
  }
  if(dist->ref->pc0 > 0){
   dist->ref->e0 = momentum2energy(dist->ref->pc0,dist->ref->mass0);
  }
  if(dist->ref->e0 > 0){
   dist->ref->pc0 = energy2momentum(dist->ref->e0,dist->ref->mass0);
  }
  dist->ref->beta0 = (dist->ref->pc0)/(dist->ref->e0);

}

void convert2standard(){


  if(dist->ref->en_like==-1){
    printf("No energy variable is set. Assume 0 deviation from reference energy \n");
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
    printf("Only allowed 1 type of energy variable! \n");
    exit(1);
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
  dist->incoord = (struct incoordinates**)malloc(linecount*sizeof(struct incoordinates*));
  dist->totincoord = linecount;
  for(int i=0; i<linecount; i++){
    dist->incoord[i] = (struct incoordinates*)malloc(sizeof(struct incoordinates));
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

void six2canonical_(double * coord, double *ref_momentum, double *mass, double *canonical){
    
    double deltap = *(coord+5); 
    double beta0 = *ref_momentum/momentum2energy(*ref_momentum, *mass);
    double beta = (*ref_momentum+deltap)/momentum2energy((*ref_momentum)+deltap, *mass);
    double rv = beta0/beta; 

     *(canonical+0) = *(coord+0);
     *(canonical+1) = *(coord+1)*(1+deltap);
     *(canonical+2) = *(coord+2);
     *(canonical+3) = *(coord+3)*(1+deltap);
     *(canonical+4) = *(coord+4)/rv;
     *(canonical+5) = *(coord+5);

}

void 
canonical2six_(double *canonical, double *ref_momentum, double *mass, double *coord){

    double deltap = *(canonical+5);
    double beta0 = *ref_momentum/momentum2energy(*ref_momentum, *mass);
    double beta = (*ref_momentum+deltap)/momentum2energy((*ref_momentum)+deltap, *mass);
    double rv = beta0/beta; 
    *(coord+0) = *(canonical+0);
    *(coord+1) = *(canonical+1)/(1+deltap);
    *(coord+2) = *(canonical+2);
    *(coord+3) = *(canonical+3)/(1+deltap);
    *(coord+4) = *(canonical+4)*rv;
    *(coord+5) = *(canonical+5);     
}

void createLinearSpaced(int length, double start, double stop, double *eqspaced ){
    
    double distance = (stop-start)/length;
    for(int i=0; i<length; i++){
        eqspaced[i] = start+distance*i;
    }
    
}

double momentum2energy(double momentum, double mass){
    return sqrt(pow(momentum,2)+pow(mass,2));
}
double energy2momentum(double energy, double mass){
    return sqrt(pow(energy,2)-pow(mass,2));
}

/* 
    if mtrx_a is (m x n) and mtrx_b is (n x p), 
    the product is an (m x p) matrix
*/
void mtrx_vector_mult_(int *mp, int *np,  double mtrx_a[6][6], double mtrx_b[6], double result[6])
{

    int m = *mp;
    int n = *np;


    register int i=0, j=0, k=0;
    for (i = 0; i < m; i++)
    {
            for (k = 0; k < n; k++)
            {
            
                result [i] += mtrx_a [i][k] * mtrx_b [k];
          
                
            }
    }
     

}

void mtrx_vector_mult_pointer(int mp, int np,  double **mtrx_a, double mtrx_b[6], double result[6])
{

    int m = mp;
    int n = np;
    result[0]=0;
    result[1]=0;
    result[2]=0;
    result[3]=0;
    result[4]=0;
    result[5]=0;

    register int i=0, j=0, k=0;
    for (i = 0; i < mp; i++)
    {
            for (k = 0; k < np; k++)
            {
                
                
                result [i] += mtrx_a [i][k] * mtrx_b [k];
                 
            }
    }

}

void solve2by2eq(double a1, double b1, double c1, double a2, double b2, double c2, double *x){
    double det = a1*b2-b1*a2;
    double dx = c1*b2-b1*c2;
    double dy = a1*c2-c1*a2;
    x[0] = (dx/det);
    x[1] = (dy/det);

//    printf("this is the solution, x ,y %f %f", x[0], x[1] );

}

/*For calculating Determinant of the Matrix */
double determinant(double a[6][6],double k)
{
  double s=1,det=0,b[6][6];
  int i,j,m,n,c;
  if (k==1)
    {
     return (a[0][0]);
    }
  else
    {
     det=0;
     for (c=0;c<k;c++)
       {
        m=0;
        n=0;
        for (i=0;i<k;i++)
          {
            for (j=0;j<k;j++)
              {
                b[i][j]=0;
                if (i != 0 && j != c)
                 {
                   b[m][n]=a[i][j];
                   if (n<(k-2))
                    n++;
                   else
                    {
                     n=0;
                     m++;
                     }
                   }
               }
             }
          det=det + s * (a[0][c] * determinant(b,k-1));
          s=-1 * s;
          }
    }
 
    return (det);
}

void cofactor(double num[6][6],double f)
{
 double b[6][6],fac[6][6];
 int p,q,m,n,i,j;
 for (q=0;q<f;q++)
 {
   for (p=0;p<f;p++)
    {
     m=0;
     n=0;
     for (i=0;i<f;i++)
     {
       for (j=0;j<f;j++)
        {
          if (i != q && j != p)
          {
            b[m][n]=num[i][j];
            if (n<(f-2))
             n++;
            else
             {
               n=0;
               m++;
               }
            }
        }
      }
      fac[q][p]=pow(-1,q + p) * determinant(b,f-1);
    }
  }
  transpose(num,fac,f);
}

/*Finding transpose of matrix*/ 
void transpose(double num[6][6],double fac[6][6],double r)
{
  int i,j;
  double b[6][6],inverse[6][6],d;
 
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
         b[i][j]=fac[j][i];
        }
    }
  d=determinant(num,r);
  for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
        inverse[i][j]=b[i][j] / d;
        }
    }
   //printf("\n\n\nThe inverse of matrix is : \n");
 
   for (i=0;i<r;i++)
    {
     for (j=0;j<r;j++)
       {
        num[i][j] = inverse[i][j];
         //printf("\t%f",inverse[i][j]);
        }
    //printf("\n");
    }
}
  void printmatrix(int m, int n, double** matrix ){
        for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%E \t", matrix[i][j]);
        }
        printf("\n");
    }
}

void printvector(const char* name, int dim, double* vector){
          printf("%s \n", name);
          for (int i = 0; i < dim; i++)
    {
            printf("%d %E \n",i, vector[i]);

    }
}


double randn(double mu, double sigma)
{
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1)
    {
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do
    {
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}

double rand_uni(double low, double high)
{

  return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

double randray(double mu, double sigma){
  double low = 0;
  double high =1;
  return pow((mu+sigma*sqrt((-2*log(rand_uni(low, high))))),2);
}