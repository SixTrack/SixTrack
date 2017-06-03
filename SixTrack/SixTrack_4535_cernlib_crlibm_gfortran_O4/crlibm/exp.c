/**
 * Function to compute the exponential on [-745, +709]
 *
 * Author : Defour David  (David.Defour@ens-lyon.fr)
 *	    Daramy Catherine (Catherine.Daramy@ens-lyon.fr)
 *
 * Date of creation : 14/03/2002   
 * Last Modified    : 28/02/2003
 */
#include <stdio.h>
#include <stdlib.h>
#include "crlibm.h"
#include "crlibm_private.h"
#include "exp.h"


/***********************************************************/
/*                      Second step                         */
/***********************************************************/



/** 
 * We use the Brent formula to evaluate exp(x)
 *
 *  1) Step 1 : range reduction
 *
 *        x  -  k * ln(2)         -ln(2)        +ln(2) 
 *  r = ---------------------     ------ <  r < ------ 
 *             512                 1024          1024
 *
 *                            (512)      (k)
 *   so that exp(x) = exp(r)^       * 2^
 *
 *
 *  2) Step 2 : polynomial evaluation
 *
 *                      2        3        4        5
 *  exp(r) ~ (1 + r + r^ /2! + r^ /3! + r^ /4! + r^ /5!
 *                      6        7        8        9  
 *                  + r^ /6! + r^ /7! + r^ /8! + r^ /9! 
 *                      10          11 
 *                  + r^  /10!  + r^  /11! )*(1+err)
 *
 *                  (-164)
 *  where |err| < 2^ 
 *
 *
 *  3) Step 3 : compute exp(r) to the power of 512
 *
 *          (512)                 2   2   2   2   2   2   2   2   2
 *   exp(r)^     = ((((((((exp(r)^  )^  )^  )^  )^  )^  )^  )^  )^
 *
 *
 *
 *  4) Step 4 : reconstruction
 *
 *                            (512)      (k)
 *            exp(x) = exp(r)^       * 2^
 *   
 *
 * Take an scs number as input and give the exponential of x
 * exactly rounded to nearest.
 *
 * Notes:
 *  - Flags need to be set to round to nearest
 */

void exp_SC(scs_ptr res_scs, double x){
  scs_t sc1, red, red_low;
  db_number db;
  int i, k;


  /* db.d = x/512  (= 2^9)  */
  
  db.d = x;                                  
  db.i[HI_ENDIAN] -= (9 << 20);              
  scs_set_d(sc1, db.d);                      
  
  
  DOUBLE2INT(k, (db.d * iln2_o512.d));       
 
  /* 1) Reduction d'argument */

  scs_set(red,     sc_ln2_o512_1_ptr);             
  scs_set(red_low, sc_ln2_o512_2_ptr);             
  if (k>0){
    scs_mul_ui(red,      (unsigned int) k);
    scs_mul_ui(red_low,  (unsigned int) k);
  }else {
    scs_mul_ui(red,      (unsigned int)(-k));
    scs_mul_ui(red_low,  (unsigned int)(-k));
    red->sign *= -1;
    red_low->sign *=-1;
  }                

  scs_sub(red, sc1, red);                    
  scs_sub(red, red, red_low);                    

  
  /* 2) Evaluation polynomiale */            
   
  scs_mul(res_scs, constant_poly_ptr[0], red);
  for(i=1; i< 10; i++){                       
    scs_add(res_scs, constant_poly_ptr[i], res_scs);
    scs_mul(res_scs, red, res_scs);            
  }
  
  scs_add(res_scs, SCS_ONE, res_scs);  
  scs_mul(res_scs, red, res_scs);
  scs_add(res_scs, SCS_ONE, res_scs);        

  /* 3) Mise a la puissance exp(r)^512  */
  for(i=0; i<9; i++)                        
    scs_square(res_scs, res_scs);
  

  res_scs->index += (int)(k/30);
  if ((k%30) > 0) 
    scs_mul_ui(res_scs, (unsigned int) (1<<((k%30))));
  else if ((k%30) < 0){
    res_scs->index --;
    scs_mul_ui(res_scs, (unsigned int) (1<<((30+(k%30)))));
  }
}


/*************************************************************
 *************************************************************
 *               ROUNDED  TO NEAREST		             *
 *************************************************************
 *************************************************************/
double scs_exp_rn(double x){ 
  scs_t res_scs;
  db_number res;
  
#if EVAL_PERF==1
  crlibm_second_step_taken++;
/*   printf("second step taken"); */
#endif


  /* 4) Reconstruction */

  exp_SC(res_scs, x);
  scs_get_d(&res.d, res_scs);

  return res.d;
}


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  -INFINITY		     *
 *************************************************************
 *************************************************************/
double scs_exp_rd(double x){ 
  scs_t res_scs;
  db_number res;
  
  exp_SC(res_scs, x);
  scs_get_d_minf(&res.d, res_scs);

  return res.d;
}


/*************************************************************
 *************************************************************
 *               ROUNDED  TOWARD  +INFINITY		     *
 *************************************************************
 *************************************************************/
double scs_exp_ru(double x){ 
  scs_t res_scs;
  db_number res;
  
  exp_SC(res_scs, x);
  scs_get_d_pinf(&res.d, res_scs);
  
  return res.d;
}

