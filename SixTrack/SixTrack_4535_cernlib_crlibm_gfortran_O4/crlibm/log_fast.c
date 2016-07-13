/*
 * Function to compute the logarithm with fully exact rounding
 *
 * Author : Daramy Catherine, Florent de Dinechin
 * (Catherine.Daramy,Florent.de.Dinechin@ens-lyon.fr)
 *
 * Date of creation : 26/08/2003   
 */
#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "log_fast.h"

/* The prototypes of the second step */
extern void scs_log(scs_ptr,db_number, int);



/* switches on various printfs. Default 0 */
#define DEBUG 0




/*
 *  1) First reduction: exponent extraction      
 *         E  
 *  x = 2^   .(y)    with  1 <= y < 2
 *
 *  log(x) = E.log(2) + log(y) where:
 *     - log(2)   is tabulated
 *     - log(y) need to be evaluated 
 *  
 *
 *  2) Avoiding accuracy problem when E=-1 by testing
 *   
 *    if (ny >= sqrt(2)) then 
 *        y = z/2;  E = E+1; 
 *    and,
 *        log(x) = (E+1).log(2) + log(y/2)
 *
 *    so now:    11/16 <= sqrt(2)/2 <= y < sqrt(2) <= 23/16
 *
 *
 *  3) Second reduction: tabular reduction
 *                    
 *     The interval 1/sqrt(2) .. sqrt(2) is divided in 8 intervals.                                 
 *     So, find the interval X_i where y is.
 *     And compute z = y - middle(X_i);
 *                                    
 *  4) Computation:
 *   
 *     Polynomial evaluation of:
 *        - P(z) ~ log(z+middle(X_i))
 *
 *                   -4      -5
 *   with  |z| < 2^   or 2^   depending the considered interval.
 *
 *
 *  5) Reconstruction:
 *   log(x) = E.log(2) + P(z)
 *
 */









static void log_quick(double *pres_hi, double *pres_lo, int* prndcstindex, db_number * py, int E) {
   double ln2_times_E_HI, ln2_times_E_LO, res_hi, res_lo;
   double z, res, P_hi, P_lo;
   int k, i;
   
    res=(double)E;
    if(E<0) E=-E;

    /* find the interval including y.d */
    i = ((((*py).i[HI_ENDIAN] & 0x001F0000)>>16)-6) ;
    if (i < 10)
      i = i>>1;
    else
      i = ((i-1)>>1);
    
    z = (*py).d - (middle[i]).d;  /* (exact thanks to Sterbenz Lemma) */
    

    /* Compute ln2_times_E = E*log(2)   in double-double */
    Add12( ln2_times_E_HI, ln2_times_E_LO, res*ln2hi.d, res*ln2lo.d); 

    /* Now begin the polynomial evaluation of log(1 + z)      */

    res = (Poly_h[i][DEGREE]).d;

    for(k=DEGREE-1; k>1; k--){
      res *= z;
      res += (Poly_h[i][k]).d;
    }

    if(E <=  EMIN_FASTPATH) {
      /* Slow path */
      if(E==0) {
	*prndcstindex = 0 ;
	/* In this case we start with a double-double multiplication to get enough relative accuracy */ 
	Mul12(&P_hi, &P_lo, res, z); 
	Add22(&res_hi, &res_lo, (Poly_h[i][1]).d,  (Poly_l[i][1]).d, P_hi, P_lo);
	Mul22(&P_hi, &P_lo, res_hi, res_lo, z, 0.); 
	Add22(pres_hi, pres_lo, (Poly_h[i][0]).d, (Poly_l[i][0]).d, P_hi, P_lo);
      } 
      else
	{
	  if(E >  EMIN_MEDIUMPATH)
	    *prndcstindex = 2; 
	  else 
	    *prndcstindex =1;
	  P_hi=res*z;
	  Add12(res_hi, res_lo, (Poly_h[i][1]).d,  (Poly_l[i][1]).d + P_hi);
	  Mul22(&P_hi, &P_lo, res_hi, res_lo, z, 0.); 
	  Add22(&res_hi, &res_lo, (Poly_h[i][0]).d, (Poly_l[i][0]).d, P_hi, P_lo);

	/* Add E*log(2)  */
	  Add22(pres_hi, pres_lo, ln2_times_E_HI, ln2_times_E_LO, res_hi, res_lo);
	}
    }
    else { /* Fast path */
      
      *prndcstindex = 3 ;
      res =   z*((Poly_h[i][1]).d + z*res);
#if 1
      Add12(P_hi,P_lo,  ln2_times_E_HI, (Poly_h[i][0]).d );
      Add12(*pres_hi, *pres_lo, P_hi, (res + ((Poly_l[i][0]).d + (ln2_times_E_LO + P_lo))));
#else
      Add12(*pres_hi, *pres_lo, 
	    ln2_times_E_HI,  
	    (Poly_h[i][0]).d + (res + ((Poly_l[i][0]).d + ln2_times_E_LO)));
#endif
    }
}





/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/
 double log_rn(double x){ 
   db_number y;
   double res_hi,res_lo,roundcst;
   int E,rndcstindex;

   E=0;
   y.d=x;

   /* Filter cases */
   if (y.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
     if (((y.i[HI_ENDIAN] & 0x7fffffff)|y.i[LO_ENDIAN])==0){
       /* return -1.0/0.0; */
	  return NInf.d;
     }                    		   /* log(+/-0) = -Inf */
     if (y.i[HI_ENDIAN] < 0){ 
       /* return (x-x)/0; */                      /* log(-x) = Nan    */
	  return NaN.d;
     }
     /* Subnormal number */
     E = -52; 		
     y.d *= two52.d; 	  /* make x a normal number    */ 
   }
    
   if (y.i[HI_ENDIAN] >= 0x7ff00000){
     return  x+x;				 /* Inf or Nan       */
   }
   
   /* reduce to  y.d such that sqrt(2)/2 < y.d < sqrt(2) */
   E += (y.i[HI_ENDIAN]>>20)-1023;				/* extract the exponent */
   y.i[HI_ENDIAN] =  (y.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;	/* do exponent = 0 */
   if (y.d > SQRT_2){
     y.d *= 0.5;
     E++;
   }

   /* Call the actual computation */
   log_quick(&res_hi, &res_lo, &rndcstindex, &y, E);
   roundcst = rncst[rndcstindex];


  /* Test for rounding to the nearest */
  if(res_hi == (res_hi + (res_lo * roundcst)))
    return res_hi;
  else {
    scs_t res;
#if DEBUG
    printf("Going for Accurate Phase for x=%1.50e\n",x);
#endif
    scs_log(res, y, E);
    scs_get_d(&res_hi, res);
    return res_hi;
  }
 }














/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY		     *
 *************************************************************
 *************************************************************/
 double log_rd(double x){ 
   db_number y;
   double res_hi,res_lo,roundcst;
   int E,rndcstindex;
   db_number absyh, absyl, u, u53;

   E=0;
   y.d=x;

 /* Filter cases */
   if (y.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
     if (((y.i[HI_ENDIAN] & 0x7fffffff)|y.i[LO_ENDIAN])==0){
       /* return -1.0/0.0; */
	  return NInf.d;     
     }                    		   /* log(+/-0) = -Inf */
     if (y.i[HI_ENDIAN] < 0){ 
      /* return (x-x)/0; */                      /* log(-x) = Nan    */
	 return NaN.d;
     }
     /* Subnormal number */
     E = -52; 		
     y.d *= two52.d; 	  /* make x as normal number = x's mantissa    */ 
   }
    
   if (y.i[HI_ENDIAN] >= 0x7ff00000){
     return  x+x;				    /* Inf or Nan       */
   }
   
 
   E += (y.i[HI_ENDIAN]>>20)-1023;				/* extract the exponent */
   y.i[HI_ENDIAN] =  (y.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;	/* do exponent = 0 */
   if (y.d > SQRT_2){
     y.d *= 0.5;
     E++;
   }

  log_quick(&res_hi, &res_lo, &rndcstindex, &y, E);
   roundcst = delta[rndcstindex];
   
   /* Rounding test to + infinity */
   absyh.d=res_hi;
   absyl.d=res_lo;
   
   absyh.l = absyh.l & 0x7fffffffffffffffLL;
   absyl.l = absyl.l & 0x7fffffffffffffffLL;
   u53.l     = (absyh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;
   u.l   = u53.l - 0x0350000000000000LL;
   
   if(absyl.d > roundcst*u53.d){ 
     if(res_lo<0.)
       res_hi -= u.d;
    return res_hi;
  }else {
    scs_t res;
#if DEBUG
    printf("Going for Accurate Phase");
#endif
    scs_log(res, y, E);
    scs_get_d_minf(&res_hi, res);
    return res_hi;
  }
}







/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY		     *
 *************************************************************
 *************************************************************/
double log_ru(double x){ 
   db_number y;
   double res_hi,res_lo,roundcst;
   int E,rndcstindex;
   db_number absyh, absyl, u, u53;

   E=0;
   y.d=x;

 /* Filter cases */
   if (y.i[HI_ENDIAN] < 0x00100000){        /* x < 2^(-1022)    */
     if (((y.i[HI_ENDIAN] & 0x7fffffff)|y.i[LO_ENDIAN])==0){
       /* return -1.0/0.0; */
	  return NInf.d;     
     }                    		   /* log(+/-0) = -Inf */
     if (y.i[HI_ENDIAN] < 0){ 
      /* return (x-x)/0; */                      /* log(-x) = Nan    */
	 return NaN.d;
     }
     /* Subnormal number */
     E = -52; 		
     y.d *= two52.d; 	  /* make x as normal number = x's mantissa    */ 
   }
    
   if (y.i[HI_ENDIAN] >= 0x7ff00000){
     return  x+x;				    /* Inf or Nan       */
   }
   
    E += (y.i[HI_ENDIAN]>>20)-1023;				/* extract the exponent */
   y.i[HI_ENDIAN] =  (y.i[HI_ENDIAN] & 0x000fffff) | 0x3ff00000;	/* do exponent = 0 */
   if (y.d > SQRT_2){
     y.d *= 0.5;
     E++;
   }

   log_quick(&res_hi, &res_lo, &rndcstindex, &y, E);
   roundcst = delta[rndcstindex];


   /* Rounding test to + infinity */
   absyh.d=res_hi;
   absyl.d=res_lo;
   
   absyh.l = absyh.l & 0x7fffffffffffffffLL;
   absyl.l = absyl.l & 0x7fffffffffffffffLL;
   u53.l     = (absyh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;
   u.l   = u53.l - 0x0350000000000000LL;
   
   if(absyl.d > roundcst*u53.d){ 
     if(res_lo>0.)    res_hi += u.d;
     return res_hi;
   }else {
     scs_t res;
#if DEBUG
     printf("Going for Accurate Phase");
#endif
     scs_log(res, y, E);
     scs_get_d_pinf(&res_hi, res);
     return res_hi;
   }
}





/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  ZERO		     *
 *************************************************************
 *************************************************************/
double log_rz(double x){ 
  if(x>1)
    return log_rd(x);
  else
    return log_ru(x);
}
