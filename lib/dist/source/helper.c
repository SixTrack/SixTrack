#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "helper.h"
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
double rationalApproximation(double t)
{
    // Abramowitz and Stegun formula 26.2.23.
    // The absolute value of the error should be less than 4.5 e-4.
    double c[] = {2.515517, 0.802853, 0.010328};
    double d[] = {1.432788, 0.189269, 0.001308};
    return t - ((c[2]*t + c[1])*t + c[0]) / 
               (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}

double normalcdfinv_(double p)
{
    printf("printing %f", p);
    if (p <= 0.0 || p >= 1.0)
    {
        return 0;
    }

    // See article above for explanation of this section.
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        return -rationalApproximation( sqrt(-2.0*log(p)) );
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        return rationalApproximation( sqrt(-2.0*log(1-p)) );
    }
}

void createLinearSpaced(int length, double start, double stop, double *eqspaced ){

    double distance = (stop-start)/length;
    for(int i; i<length; i++){
        eqspaced[i] = distance*i;
    }
    
}

double momentum2energy(double momentum, double mass){
    return sqrt(pow(momentum,2)+pow(mass,2));
}

void print2file(double **coord, int length){
    FILE *fptr;

   fptr = fopen("coordinates.out", "w");
   if(fptr == NULL)
   {
      printf("Error!");
      exit(1);
   }
   for(int i=0;i < length; i++)
   {
        fprintf(fptr,"%f  %f  %f  %f %f  %f  \n", coord[i][1], coord[i][2],coord[i][3], coord[i][4], coord[i][5], coord[i][6]);
   }
   fclose(fptr);

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
       //  printf("\t%f",inverse[i][j]);
        }
    //printf("\n");

     }
}
