#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "helper.h"
#include "distinput.h"
#include "outputdist.h"

/*This function converts from canoncial to sixtrack tracking variables*/
void 
canonical2six(double *canonical, double beta0, double pc0, double mass0, double mass, double *coord){
    double deltap = canonical[5];
    double beta = (pc0+deltap)/momentum2energy((pc0+deltap), mass);
    double rv = beta0/beta;
    double factor = 1e3;
    coord[0] = canonical[0]*factor;
    coord[1] = canonical[1]*(factor/(1+deltap));
    coord[2] = canonical[2]*factor;
    coord[3] = canonical[3]*(factor/(1+deltap));
    coord[4] = canonical[4]*(factor*rv);
    coord[5] = canonical[5];  
}
/* Compare two strings s1 and s2, assuming either is
 * terminated by \n or a NULL, A match
 * returns 0, a non-match returns 1.
 */
int
strcmpnl(const char *s1, const char *s2)
{
  char s1c;
  char s2c;
  do
    {
      s1c = *(s1++);
      s2c = *(s2++);
      if (s1c == '\n')
          s1c = 0;
      if (s2c == '\n')
          s2c = 0;
      if (s1c != s2c)
          return 1;
    }
  while (s1c);          /* already checked *s2 is equal */
  return 0;
}
void issue_error(const char* t1){
  printf("+=+=+= fatal: %s\n",t1);
  exit(1);
}
void issue_error2(const char* t1,const char* t2){
  printf("+=+=+= fatal: %s %s\n",t1, t2);
  exit(1);
}

void issue_warning(const char* t1){
  printf( "+=+=+= warning:  %s\n",t1);
}

void issue_info(const char* t1){
  printf( "Info: %s\n",t1);
}

double psigma2deltap(double psigma, double beta0 ){
  return (sqrt(pow(psigma*beta0,2) +2*psigma +1)-1);
}
double sigma2zeta(double sigma, double beta0, double beta){
  return sigma*(beta/beta0);
}
double tau2zeta(double tau, double beta){
  return tau*beta;  
}
 
double pt2deltap(double pt, double beta0 ){
  return (sqrt(pow(pt,2) +2*pt*beta0 +1)-1);
}
// This function prevents several energy from beeing set
void checkifenergyset(int entype){
  if(dist->ref->en_like==-1)
    dist->ref->en_like=entype;
  else {
    issue_error("Only allowed 1 type of energy variable!");
  }
}

// This function prevents several energy from beeing set
void checkiftimeset(int entype){
  if(dist->ref->time_like==-1)
    dist->ref->time_like=entype;
  else {
    issue_error("Only allowed 1 type of time variable!");
  }
}

void checkifangset(int entype){
  if(dist->ref->time_like==-1)
    dist->ref->time_like=entype;
  else {
    issue_error("Only allowed 1 type of angle variable!");
  }
}

double momentum2energy(double momentum, double mass){
    return sqrt(pow(momentum,2)+pow(mass,2));
}

double energy2momentum(double energy, double mass){
    double momsq = pow(energy,2)-pow(mass,2);
    if(momsq <= 0){
      issue_error("Energy must be larger than mass!");
    } 
    return sqrt(momsq);
}