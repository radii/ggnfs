/**************************************************************/
/* smintfact.c                                                */
/* Copyright 2004, Chris Monico.                              */
/**************************************************************/
/*  This file is part of GGNFS.
*
*   GGNFS is free software; you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation; either version 2 of the License, or
*   (at your option) any later version.
*
*   GGNFS is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with GGNFS; if not, write to the Free Software
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "ggnfs.h"


/********************************************************************/
int prho(mpz_t factor1, mpz_t factor2, mpz_t n, s32 c, s32 maxIts)
/* It is assumed that the input is composite. */
{ static mpz_t a, b, oldA, oldB, tmp, tmp2;
  static int initialized=0;
  s32   i, its;

  if (!(initialized)) {
    mpz_init(a); mpz_init(b); mpz_init(tmp); mpz_init(tmp2);
    mpz_init(oldA); mpz_init(oldB);
    initialized=1;
  }
  if (c<1) c=1;
  if (c>0x000000FF) c &= 0x000000FF;
  mpz_set_ui(a, 1);
  mpz_set_ui(b, 1);
  its=0;
  do {
    mpz_set(oldA, a); mpz_set(oldB, b);
    mpz_set_ui(tmp2, 1);
    for (i=25; i>0; i--) {
      mpz_mul(tmp, a, a); mpz_add_ui(tmp, tmp, c); mpz_mod(a, tmp, n);
      mpz_mul(tmp, b, b); mpz_add_ui(tmp, tmp, c); mpz_mod(b, tmp, n);
      mpz_mul(tmp, b, b); mpz_add_ui(tmp, tmp, c); mpz_mod(b, tmp, n);

      mpz_sub(tmp, a, b); mpz_mul(tmp2, tmp2, tmp); mpz_mod(tmp2, tmp2, n);
    }
    its += 25;
    mpz_gcd(tmp, tmp2, n);
  } while ((mpz_cmp_ui(tmp, 1)==0) && (its<maxIts));
  if (its >= maxIts) return -1;
  if (mpz_cmp(tmp, n) == 0) {
    mpz_set(a, oldA); mpz_set(b, oldB);
    do {
      mpz_mul(tmp, a, a); mpz_add_ui(tmp, tmp, c); mpz_mod(a, tmp, n);
      mpz_mul(tmp, b, b); mpz_add_ui(tmp, tmp, c); mpz_mod(b, tmp, n);
      mpz_mul(tmp, b, b); mpz_add_ui(tmp, tmp, c); mpz_mod(b, tmp, n);
      mpz_sub(tmp, a, b);
      mpz_gcd(tmp, tmp, n);
    } while (mpz_cmp_ui(tmp, 1)==0);
  }
  if (mpz_cmp(tmp, n) < 0) {
    mpz_set(factor1, tmp);
    mpz_div(factor2, n, factor1);
    return 0;
  }
  return -1; /* Unknown failure. */
}
    

/******************************************************************/    
const s32 td_primes[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53};
const s32 num_td_primes=16;
/******************************************************************/    
int doTrialDiv(s32 *factors, mpz_t n)
/******************************************************************/    
{ s32 i;
  int  numFacts=0;
  static mpz_t q, r;
  static int initialized=0;

  if (!(initialized)) {
    mpz_init(q); mpz_init(r);
    initialized=1;
  }

  for (i=0; i<num_td_primes; i++) {
    while (mpz_tdiv_qr_ui(q, r, n, td_primes[i])==0) {
      factors[numFacts++] = td_primes[i];
      mpz_set(n, q);
    }
  }
  return numFacts;
}
  
    

/* This function is used to factor relations only. It is thus
   more important that it be fast than completely correct. So
   we use only one iteration of Miler-Rabin. In fact, we could
   even consider using a faster homebrewed test at some point.
*/
/********************************************************************/
int factor(u32 *factors, mpz_t n, int useTrialDivision)
/********************************************************************/
{ int    numFactors=0, res, sorted=1, i, retVal;
  static mpz_t div1, div2, remain;
  static __mpz_struct stack[32];
  static int initialized=0, stackSize=0;
  s32 c;

  if (!(initialized)) {
    mpz_init(div1); mpz_init(div2); mpz_init(remain);
    for (i=0; i<32; i++) mpz_init(&stack[i]);
    initialized=1;
  }
  if (mpz_cmp_ui(n, 1)==0) {
    factors[0]=1;  return 0;
  }
  mpz_abs(remain, n);
  while (mpz_even_p(remain)) {
    factors[numFactors++]=2;
    mpz_tdiv_q_2exp(remain, remain, 1);
  }

  if (mpz_probab_prime_p(remain, 1)) {
    if (mpz_fits_ulong_p(remain)) {
      factors[numFactors++] = mpz_get_ui(remain);
      return numFactors;
    } else return -1;
  }
  if (useTrialDivision) {
    numFactors += doTrialDiv(&factors[numFactors], remain);
    if (mpz_cmp_ui(remain, 1)==0) return numFactors;
    if (mpz_probab_prime_p(remain, 1)) {
      if (mpz_fits_ulong_p(remain)) {
        factors[numFactors++] = mpz_get_ui(remain);
        return numFactors;
      } else return -1;
    }
  } 

  /* 50000 is sufficient for all primes below 100,000,000. */
  c = 1;
  do {
    if (c <= 2)
      res = prho(div1, div2, remain, c, 50000);
    else
      res = prho(div1, div2, remain, c, c*50000);
    c++;
  } while (res && (c <4));

  if (res) {
#define _QUIET
#ifndef _QUIET
    printf("Gave up on factoring "); mpz_out_str(stdout, 10, remain); 
    printf(" (useTrialDivision = %d)\n", useTrialDivision);
#endif
    return -1;
  }

  if (mpz_cmp(div1, div2) > 0) mpz_swap(div1, div2);
  if (mpz_probab_prime_p(div1, 1)) {
    if (mpz_fits_ulong_p(div1)) {
      factors[numFactors++] = mpz_get_ui(div1);
    } else return -1; /* Prime factor that doesn't fit in a s32. */
  } else {
    mpz_set(&stack[stackSize++], div2);
    retVal = factor(&factors[numFactors], div1, 1);
    mpz_set(div2, &stack[--stackSize]);
    if (retVal >=0) numFactors += retVal;
    else return retVal;
  }
  if (mpz_probab_prime_p(div2, 1)) {
    if (mpz_fits_ulong_p(div2)) {
      factors[numFactors++] = mpz_get_ui(div2);
    } else return -1; /* Prime factor that doesn't fit in a s32. */
  } else {
    retVal = factor(&factors[numFactors], div2, 1);
    if (retVal >=0) numFactors += retVal;
    else return retVal;
  }
    
  /* Here we should sort the factors. */
  for (i=0; i<(numFactors-1); i++) 
    if (factors[i] > factors[i+1])
      sorted=0;
  if (!(sorted)) {
    /***********************************************************/
    /* This will only happen very very rarely, so it's ok to   */
    /* use quick sort, even though it's only a small number of */
    /* integers, and we could sort it much faster in an ad-hoc */
    /* way.                                                    */
    /***********************************************************/
    qsort(factors, numFactors, sizeof(s32), cmpU32s);
  }

  return numFactors;
}  




  
  

