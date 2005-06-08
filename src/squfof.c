/**************************************************************/
/* clsieve.c                                                  */
/* Jason Papadopoulos, 2/6/05                                 */
/* Factor up to 62-bit integers using SQUFOF                  */
/**************************************************************/
/* Copyright 2004, Jason Papadopoulos                         */
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
***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "ggnfs.h"

/* This SQUFOF implementation will race as many multipliers
   of the input as will fit into 62 bits */

const u32 multipliers[] = {1, 3, 5, 7, 
                      11, 3*5, 3*7, 3*11, 
                      5*7, 5*11, 7*11, 
                      3*5*7, 3*5*11, 3*7*11, 
                      5*7*11, 3*5*7*11 };

#define MAX_MULTIPLIERS (sizeof(multipliers)/sizeof(u32))
#define QSIZE 20

/* the maximum number of inner loop iterations for all 
   multipliers combined */

#define MAX_CYCLES 100000

/* the number of iterations to do before switching to the
   next multiplier */

#define ONE_CYCLE_ITER 300

typedef struct {
  u64 n;
  u32 sqrtn[MAX_MULTIPLIERS]; 
  u32 cutoff[MAX_MULTIPLIERS];
  u32 q0[MAX_MULTIPLIERS];
  u32 p1[MAX_MULTIPLIERS];
  u32 q1[MAX_MULTIPLIERS];
  u32 num_saved[MAX_MULTIPLIERS];
  u16 saved[MAX_MULTIPLIERS][QSIZE];
  u8 failed[MAX_MULTIPLIERS];
} squfof_data_t;


/*************************************************/
/* A version of gcd() that takes unsigned inputs */
/*************************************************/
u32 gcd_ui(u32 x, u32 y) 
{ u32 tmp;

  if (y<x) {
    tmp= x; x = y; y = tmp;
  }

  while (y>0) {
    x = x%y;
    tmp = x; x = y; y = tmp;
  }
  return x;
}

/***********************************************************/
/* Perform one unit of work for SQUFOF with one multiplier */
/***********************************************************/
u32 squfof_one_cycle(squfof_data_t *data, u32 mult_idx, 
                     u32 num_iter, u32 *factor) 
{ u32 sqrtn = data->sqrtn[mult_idx];
  u32 cutoff = data->cutoff[mult_idx];
  u32 num_saved = data->num_saved[mult_idx];
  u32 multiplier = 2*multipliers[mult_idx];
  u32 coarse_cutoff = cutoff*multiplier;
  u16 *saved = data->saved[mult_idx];
  u32 q0 = data->q0[mult_idx];
  u32 p1 = data->p1[mult_idx];
  u32 q1 = data->q1[mult_idx];

  u32 i, j;
  u32 p0 = 0;
  u32 sqrtq = 0;
  u64 scaledn;

  for (i=0; i<num_iter; i++) {    
    u32 q, bits, tmp;

    /* compute (sqrtn+p1)/q1; since this is unity
       more than half the time, special case the
       division to save some effort */

    tmp = sqrtn+p1-q1;
    q = 1;
    if (tmp>=q1)
      q += tmp/q1;

    /* compute the next numerator and denominator */
    
    p0 = q*q1-p1;
    q0 = q0+(p1-p0)*q;

    /* In order to avoid trivial factorizations,
       save q1 values that are small in a list. Only
       values that do not contain factors of the
       multiplier must be saved; however, the GCD
       is 5x more expensive than all the rest of the
       work performed in the loop. To avoid most GCDs,
       only attempt one if q1 is less than the cutoff
       times the full multiplier 
       
       If the queue overflows (very rare), signal that
       this multiplier failed and move on to another one */

    if (q1<coarse_cutoff) {
      tmp = q1/gcd_ui(q1, multiplier);

      if (tmp<cutoff) {
        if (num_saved>=QSIZE) {
          data->failed[mult_idx] = 1;
          return i;
        }
        saved[num_saved++] = (u16)tmp;
      }
    }

    /* if q0 is a perfect square, then the factorization
       has probably succeeded. Most of the squareness
       tests out there require multiple divisions and 
       complicated loops. We can approximate these tests
       by doing two things: testing that the number of
       trailing zeros in q0 is even, and then testing
       if q0 shifted right this many places is 1 mod 8. */

    bits = 0; 
    tmp = q0;
    while (!(tmp&1)) {
      bits++;
      tmp >>= 1;
    }
    if (!(bits&1) && ((tmp&7) == 1)) {

      /* q0 is probably a perfect square. Take the
         square root by cheating */

      sqrtq = (u32)(sqrt((double)q0));

      if (sqrtq*sqrtq == q0) {

        /* it *is* a perfect square. If it has
           not appeared previously in the list
           for this multiplier, then we're almost
           finished */

        for (j=0; j<num_saved; j++) {
          if (saved[j]==sqrtq)
            break;
        }

        if (j==num_saved)
          break;
      }
    }

    /* perform the odd half of the SQUFOF cycle */

    tmp = sqrtn+p0-q0;
    q = 1;
    if (tmp>=q0)
      q += tmp/q0;

    p1 = q*q0-p0;
    q1 = q1+(p0-p1)*q;

    if (q0<coarse_cutoff) {
      tmp = q0/gcd_ui(q0, multiplier);

      if (tmp<cutoff) {
        if (num_saved>=QSIZE) {
          data->failed[mult_idx] = 1;
          return i;
        }
        saved[num_saved++] = (u16)tmp;
      }
    }
  }

  if (sqrtq==1) {

    /* the above found a trivial factor, so this
       multiplier has failed */

    data->failed[mult_idx] = 1;
    return i;
  }
  else if (i == num_iter) {
  
    /* no square root found; save the parameters
       and go on to the next multiplier */

    data->q0[mult_idx] = q0;
    data->p1[mult_idx] = p1;
    data->q1[mult_idx] = q1;
    data->num_saved[mult_idx] = num_saved;
    return i;
  }

  /* square root found; the algorithm cannot fail now.
     Compute the inverse quadratic form and iterate */

  q0 = sqrtq;
  p1 = p0+sqrtq*((sqrtn-p0)/sqrtq);
  scaledn = data->n * (u64)multipliers[mult_idx];
  q1 = (scaledn - (u64)p1*(u64)p1) / (u64)q0;

  while (1) {
    u32 q, tmp;

    tmp = sqrtn+p1-q1;
    q = 1;
    if (tmp>=q1)
      q += tmp/q1;

    p0 = q*q1-p1;
    q0 = q0+(p1-p0)*q;

    if (p0==p1) {
      q0 = q1;
      break;
    }

    tmp = sqrtn+p0-q0;
    q = 1;
    if (tmp >= q0)
      q += tmp/q0;

    p1 = q*q0-p0;
    q1 = q1+(p0-p1)*q;

    if (p0==p1)
      break;
  }

  /* q0 is the factor of n. Remove factors that exist
     in the multiplier and save whatever's left */

  q0 = q0/gcd_ui(q0, multiplier);
  *factor = q0;
  return i;
}

/***************************************************************/
/* Factor a number up to 62 bits in size using SQUFOF.         */
/*   For n the product of two primes, this routine will        */
/*   succeed with very high probability, although the          */
/*   likelihood of failure goes up as n increases in size.     */
/*   Empirically, 62-bit factorizations fail about 5% of the   */
/*   time; for smaller n the failure rate is nearly zero.      */
/*                                                             */
/*   If a factor is found, it is returned. If it is a trivial  */
/*   factor, the return value is 1. If SQUFOF failed, the      */
/*   return value is zero                                      */
/***************************************************************/
u32 squfof(mpz_t n) 
{ u32 num_mult;
  u32 factor_found = 0;
  u32 i, num_iter, num_failed;
  squfof_data_t data;
  static mpz_t tmp, sqrt_tmp;
  static s32 initialized = 0;
  size_t t;

  if (!initialized) { 
    initialized = 1;
    mpz_init(tmp);
    mpz_init(sqrt_tmp);
  }

  /* turn n into a u64 */
#ifdef LONG64
  data.n = mpz_get_ui(n);
#else
  mpz_export(&data.n, &t, -1, sizeof(u64), 0, 0, n);
#endif

  /* for each multiplier */

  for (i=0; i<MAX_MULTIPLIERS; i++) {

    /* use the multiplier if the multiplier times n
       will fit in 62 bits. Because multipliers are
       initialized in order of increasing size, when
       one is too big then all the rest are also too big */

    mpz_mul_ui(tmp, n, multipliers[i]);
    if (mpz_sizeinbase(tmp,2)>62)
      break;

    mpz_sqrtrem(sqrt_tmp, tmp, tmp);

    /* initialize the rest of the fields for this
       multiplier */

    data.sqrtn[i] = (u32)mpz_get_ui(sqrt_tmp);
    data.cutoff[i] = (u32)(sqrt(2.0*(double)data.sqrtn[i]));
    data.num_saved[i] = 0;
    data.failed[i] = 0;

    data.q0[i] = 1;
    data.p1[i] = data.sqrtn[i];
    data.q1[i] = (u32)mpz_get_ui(tmp);

    /* if n is a perfect square, don't run the algorithm;
       the factorization has already taken place */

    if (data.q1[i] == 0)
      return data.p1[i];
  }
  if (i==0)
    return 0;

  /* perform a block of work using each multiplier in
     turn, until our budget of work for factoring n is
     exhausted */

  num_mult = i;
  num_iter = 0;
  num_failed = 0;
  while (num_iter<MAX_CYCLES) {
    
    /* for each cycle of multipliers, begin with the
       multiplier that is largest. These have a higher
       probability of factoring n quickly */

    for (i=num_mult-1; (int)i>=0 ; i--) {
      if (data.failed[i])
        continue;

      /* work on this multiplier for a little while */
      num_iter += squfof_one_cycle(&data, i, ONE_CYCLE_ITER, &factor_found);

      /* if all multipliers have failed, then SQUFOF has failed */
      if (data.failed[i]) {
        if (++num_failed==num_mult)
          return 0;
      }
      if (factor_found)
        return factor_found;
    }
  }

  return 0;
}
