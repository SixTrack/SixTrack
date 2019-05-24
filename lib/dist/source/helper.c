#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
#include "distinput.h"
#include "outputdist.h"

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