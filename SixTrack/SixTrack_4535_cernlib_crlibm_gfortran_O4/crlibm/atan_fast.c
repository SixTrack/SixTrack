/*
 * Function to compute atan with fully exact rounding
 *
 * Author : Gast Nicolas
 * nicolas.gast@ens.fr
 *
 * Date : 24/06/2004
 */

#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "atan_fast.h"

/* the second step :   A VOIR PLUS TARD */ 
double scs_atan_rn(double); 
double scs_atan_rd(double); 
double scs_atan_ru(double); 

static void atan_quick(double *atanhi,double *atanlo, int *index_of_e, double x) {

  double tmphi,tmplo, x0hi,x0lo;
  double q,Xred2,x2;
  double Xredhi,Xredlo;
  double atanlolo;
  double tmphi2, tmplo2;

  
  int i;
  
  if (x > MIN_REDUCTION_NEEDED) /* test if reduction is necessary : */
    {
      /*
       * 1) Argument reduction : 
       * 
       *  tan(x) = tan( b(i) ) + tan ( (x-b(i)) / (1+x*b(i)))
       *
       *                                                     6.3
       * we choose 62 b(i) so that (x-b(i)) / (1+x*b(i)) < 2^
       */ 
      
      double xmBihi, xmBilo;
      
      if (x > arctan_table[61][B].d) {
        i=61;
        Add12( xmBihi , xmBilo , x , -arctan_table[61][B].d);
      }
      else 
        {
          /* compute i so that a[i] < x < a[i+1] */
          i=31;
          if (x < arctan_table[i][A].d) i-= 16;
          else i+=16;
          if (x < arctan_table[i][A].d) i-= 8;
          else i+= 8;
          if (x < arctan_table[i][A].d) i-= 4;
          else i+= 4;
          if (x < arctan_table[i][A].d) i-= 2;
          else i+= 2;
          if (x < arctan_table[i][A].d) i-= 1;
          else i+= 1;
          if (x < arctan_table[i][A].d) i-= 1;     
          xmBihi = x-arctan_table[i][B].d;
          xmBilo = 0.0;
        }
        
      /* we now compute Xred = ( x-b[i] ) / ( 1 + x*b[i] )
       * 
       * def : x0 := 1+x*b[i]
       *
       * 1st we compute an approximation of y = 1/x0
       * then we compute a better approx x' = y*(2-x0*y)
       * we can proove that :
       * if y = 1/x0*(1+e) 
       *     then x' = 1/x0 * (1-e^2)
       *                   
       */
      
      Mul12(&tmphi,&tmplo, x, arctan_table[i][B].d);

      if (x > 1)
        Add22(&x0hi,&x0lo,tmphi,tmplo, 1.0,0.0);
      else {Add22( &x0hi , &x0lo , 1.0,0.0,tmphi,tmplo);}

      Div22( &Xredhi, &Xredlo, xmBihi , xmBilo , x0hi,x0lo);

      /* Polynomial evaluation : 
       *  
       *  1rt compute Q(x^2) = (1 - x^2/3 + ...)
       *      then P(x) = x * Q(x^2)
       *
       */

      Xred2 = Xredhi*Xredhi;
      
      q = Xred2*(coef_poly[3]+Xred2*
                 (coef_poly[2]+Xred2*
                  (coef_poly[1]+Xred2*
                   coef_poly[0]))) ;

      /* reconstruction : atan(x) = atan(b[i]) + atan(x) */
      atanlolo = Xredlo+ arctan_table[i][ATAN_BLO].d + Xredhi*q;
      Add12( tmphi2, tmplo2, arctan_table[i][ATAN_BHI].d, Xredhi);
      Add12( *atanhi, *atanlo, tmphi2, (tmplo2+atanlolo));
      
      if (i<10)
        *index_of_e = 0;
      else
        *index_of_e = 1;
    }
  else 
    // no reduction needed
    {
      /* Polynomial evaluation : 
       *  
       *  1rt compute Q(x^2) = (1 - x^2/3 + ...)
       *      then P(x) = x * Q(x^2)
       *
       */
      
      x2 = x*x;
      q = x2*(coef_poly[3]+x2*
              (coef_poly[2]+x2*
               (coef_poly[1]+x2*
                coef_poly[0]))) ;
      Add12(*atanhi,*atanlo, x , x*q);
      
      if (x > 0x3F500000)
        *index_of_e = 2;
      else
        *index_of_e = 3;
    }
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST			     *
 *************************************************************
 *************************************************************/

/*
  static double atan_rn2(double x);
  extern double atan_rn (double x) {  
  if(x>0)
  return atan_rn2(x);
  else
  return -atan_rn2(-x);
  }
*/
extern double atan_rn(double x) {
  unsigned int hx;
  double atanhi,atanlo;
  int index_of_e;
  double sign;
  db_number x_db;
  x_db.d = x;
  hx = x_db.i[HI_ENDIAN] & 0x7fffffff; 
  if(x>=0)
    {sign=1;}
  else
    {x=-x;
    sign =-1;}
  
  /* Filter cases */
  if ( hx >= 0x43500000)           /* x >= 2^54 */
    {
      if ( (hx & 0x000fffff) == 0x00080000 && x_db.i[LO_ENDIAN] == 0)
        return x;                /* NaN */
      else
        return sign*HALFPI.d;           /* atan(x) = Pi/2 */
    }
  else
    if ( hx < 0x3E400000 )
      {return sign*x;}                   /* x<2^-27 then atan(x) =~ x */

  atan_quick(&atanhi, &atanlo,&index_of_e ,x);
  
  if (atanhi == (atanhi + (atanlo*rncst[index_of_e]))) 
    return sign*atanhi;
  else
    {
      /* more accuracy is needed , lauch accurate phase */ 
      return sign*scs_atan_rn(x);
    }
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY		     *
 *************************************************************
 *************************************************************/
extern double atan_rd(double x) {
  unsigned int hx;
  int sign;
  double atanhi,atanlo;
  db_number absyh, absyl, u, u53;
  int index_of_e;
  double roundcst;
  db_number x_db;
  x_db.d = x;
  hx = x_db.i[HI_ENDIAN] & 0x7FFFFFFF; 

  if (x<0)
    {sign = -1;
    x = -x;}
  else sign = 1;
  
  /* Filter cases */
  if ( hx >= 0x43500000)           /* x >= 2^54 */
    {
      if (  (hx & 0x000fffff) == 0x00080000 && x_db.i[LO_ENDIAN] == 0 )
        return x;                /* NaN */
      else
        if (sign>0)
          return HALFPI.d;
        else
          return -HALFPI_TO_PLUS_INFINITY.d;           /* atan(x) = Pi/2 */
    }
  else
    if ( hx < 0x3E400000 )
      {if (sign>0)
        {if(x==0)
          {x_db.i[HI_ENDIAN]  = 0x80000000;
          x_db.i[LO_ENDIAN] = 0;}
        else
          x_db.l--;
        return x_db.d;
        }
      else
        return -x;
      }
  
  atan_quick(&atanhi, &atanlo,&index_of_e, x);
  roundcst = epsilon[index_of_e];
  atanhi = sign*atanhi;
  atanlo = sign*atanlo;
  
  /* Rounding test to - infinity */ 
  
  absyh.d=atanhi;
  absyl.d=atanlo;
  
  absyh.l = absyh.l & 0x7fffffffffffffffLL;
  absyl.l = absyl.l & 0x7fffffffffffffffLL;
  u53.l     = (absyh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;
  u.l   = u53.l - 0x0350000000000000LL;
  
  if(absyl.d > roundcst*u53.d){
    if(atanlo<0.)
      {atanhi -= u.d;}
    return atanhi;
  }
  else {
    return scs_atan_rd(sign*x);
  }
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY		     *
 *************************************************************
 *************************************************************/

extern double atan_ru(double x) {
  unsigned int hx;
  db_number x_db;
  double atanhi,atanlo;
  db_number absyh, absyl, u, u53;
  int index_of_e;
  int sign;
  double roundcst;

  if (x<0)
    {sign = -1;
    x = -x;}
  else 
    sign = 1;
  
  x_db.d = x;
  hx = x_db.i[HI_ENDIAN] & 0x7FFFFFFF; 
  
  /* Filter cases */
  if ( hx >= 0x43500000)           /* x >= 2^54 */
    {
      if (  (hx & 0x000fffff) == 0x00080000 && x_db.i[LO_ENDIAN] == 0 )
        return x;                /* NaN */
      else
        {
          if (sign>0)
            return HALFPI_TO_PLUS_INFINITY.d;
        else
          return -HALFPI.d;           /* atan(x) = Pi/2 */
        }
    }
  else
    
    if ( hx < 0x3E400000 )
      {if (sign<0)
        {x_db.l--;
        return -x_db.d;
        }
      else
        if(x==0)
          return 0;
      return x;
      }                   /* x<2^-27 then atan(x) =~ x */
  
  atan_quick(&atanhi, &atanlo, &index_of_e, x);
  roundcst = epsilon[index_of_e];
  atanhi = sign*atanhi;
  atanlo = sign*atanlo;
  
  /* Rounding test to + infinity */ 
  
  absyh.d=atanhi;
  absyl.d=atanlo;
  
  absyh.l = absyh.l & 0x7fffffffffffffffLL;
  absyl.l = absyl.l & 0x7fffffffffffffffLL;
  u53.l     = (absyh.l & 0x7ff0000000000000LL) +  0x0010000000000000LL;
  u.l   = u53.l - 0x0350000000000000LL;
  
  if(absyl.d > roundcst*u53.d){
    if(atanlo>0.)
      {atanhi += u.d;}
    return atanhi;
  }
  else {
    return scs_atan_ru(sign*x);
  }
}

/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  ZERO		     *
 *************************************************************
 *************************************************************/

extern double atan_rz(double x) {
  if (x>0)
    return atan_rd(x);
  else
    return atan_ru(x);
}
