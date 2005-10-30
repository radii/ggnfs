/**************************************************************/
/* misc.c                                                     */
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

#ifdef _MSC_VER
#pragma warning (disable: 4996) /* warning C4996: 'function' was declared deprecated */
#endif

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#if !defined(_MSC_VER)
#include <sys/time.h>
#endif

#include "ggnfs.h"

#ifdef MALLOC_REPORTING
#if defined(__MINGW32__) || defined(MINGW32)
#include "malloc.h"
#else
#include <malloc.h>
#endif
#endif

#define MAX_MSG_SIZE 256


/*******************************************************/
INLINE double sTime()
/*******************************************************/
#if !defined(_MSC_VER) && !defined(__MINGW32__) && !defined(MINGW32)
{ static struct  timeval  this_tv;
  static struct  timezone dumbTZ;
  double t;
                                                                                                                                                         
  gettimeofday(&this_tv, &dumbTZ);
  t = this_tv.tv_sec + 0.000001*this_tv.tv_usec;
  return t;
}
#else
{
        return clock() / (double)CLOCKS_PER_SEC;
}
#endif

/**************************************************/
int readBinField(char *str, int size, FILE *fp)
/**************************************************/
{ int i;

  for (i=0; (i<(size-1)) && (!(feof(fp))); i++) {
    str[i] = fgetc(fp);
    if (str[i]==BF_DELIMITER) {
      str[i]=0;
      return i;
    } else if (str[i]=='\n') {
      str[i]=0;
      return i;
    }
  }
  str[i]=0;
  return i;
}

/**************************************************/
double _mpz_log(mpz_t k)
/**************************************************/
/* Approximate log(|k|).                          */
/**************************************************/
{ double res;
  long e;
  static mpz_t _absk;
  static int initialized=0;

  if (!(initialized)) {
    mpz_init(_absk); initialized=1;
  }
  mpz_abs(_absk, k);
  res = mpz_get_d_2exp(&e, _absk);
  if (res <= 0 ) {
    printf("_mpz_log() invalid argument! (");
    mpz_out_str(stdout, 10, k); printf(")\n");
    return 0.0;
  }
  return (double)e*M_LN2 + log(res);
}


/**************************************************/
int writeBinField(FILE *fp, char *str)
/**************************************************/
{
  fprintf(fp, "%s%c", str, BF_DELIMITER);
  return 0;
}

/***********************************************************************/
s32 powMod(s32 op, s32 n, s32 p)
/***********************************************************************/
/* Return op^n mod p.                                                  */
/***********************************************************************/
{ s32 remain, tmp, res;

  remain = n;
  tmp = op;
  res = 1;
  while (remain) {
    if (remain&0x01)
      res = mulmod32(res, tmp, p);
    tmp = mulmod32(tmp, tmp, p);
    remain = remain/2;
  }
  return res;
}
  
/*********************************************************************/
INLINE s32 inverseModP(s32 n, s32 p)
/*********************************************************************/
{ s32 a=n, b=p, q, r;
  s32 x, x1=0, x2=1;

  while (b>0) {
    q = a/b;
    r = a-q*b;
    x = x2 - q*x1;
    a=b;
    b=r;
    x2=x1;
    x1 = x;
  }
  x = (x2 < 0) ? (x2 + p) : x2;
  x = x%p;
  return x;
}
  
/*********************************************************************/
INLINE s32 gcd(s32 x, s32 y)
/*********************************************************************/
{ s32 a, b, r;

  a = (x>0) ? x : -x;
  b = (y>0) ? y : -y;
  
  while (b) {
    r = a%b;
    a = b;
    b = r;
  }
  return a;
}

        
/***************************************************************************/   
int sqrtModP(mpz_t res, mpz_t x2, mpz_t p)
/***************************************************************************/   
/* Given a quadratic residue $x2$ mod $p$, compute a squareroot.           */
/***************************************************************************/   
{ s32   i;
  s32   e, r, m;
  static mpz_t mz, mn, mp, ma, mx, mb, my, mt, mtmp, q;
  static int initialized=0;

  /* Use Tonelli and Shanks, ala [Cohen, "A Course in ...", Alg. 1.5.1] */
  if (!(initialized)) {
    mpz_init(mz); mpz_init(mn); mpz_init(mp); mpz_init(ma);
    mpz_init(mx); mpz_init(mb); mpz_init(my); mpz_init(mt);
    mpz_init(mtmp); mpz_init(q);
    initialized=1;
  }
  
  /* Step 1: */
  e=0;
  mpz_sub_ui(q, p, 1);
  while ((mpz_sgn(q)) && (mpz_mod_ui(mtmp, q, 2)==0)) {
    e++;
    mpz_div_ui(q, q, 2);
  }
  
  mpz_set(mp, p); mpz_set_ui(mn, 2); mpz_set(ma, x2);
  do {
    mpz_add_ui(mn, mn, 1);
  } while (mpz_jacobi(mn, mp) != -1);
  mpz_powm(mz, mn, q, mp);
  
  
  /* Step 2: */
  mpz_set(my, mz);
  r = e;
  mpz_sub_ui(mtmp, q, 1);
  mpz_div_ui(mtmp, mtmp, 2);
  mpz_powm(mx, ma, mtmp, mp);
  mpz_set(mb, ma); mpz_mul(mb, mb, mx); mpz_mul(mb, mb, mx);
  mpz_mod(mb, mb, mp);
  
  mpz_mul(mx, ma, mx); mpz_mod(mx, mx, mp);
  
  /* Step 3: */
C1_5_1_S3_mp:
  mpz_mod(mtmp, mb, p);
  if (mpz_cmp_ui(mtmp, 1)==0) {
    mpz_set(res, mx);
    goto C1_5_1_DONE_mp;
  }
  m=1;
  mpz_mul(mtmp, mb, mb); mpz_mod(mtmp, mtmp, mp);
  while (mpz_cmp_ui(mtmp, 1)) {
    m++;
    mpz_mul(mtmp, mtmp, mtmp);  mpz_mod(mtmp, mtmp, mp);
  }
  if (m==r) {
    mpz_set_ui(res, -1);
    goto C1_5_1_DONE_mp;
  }

  /* Step 4: */
  mpz_set(mt, my);
  for (i=0; i< (r-m-1); i++) {
    mpz_mul(mt, mt, mt); mpz_mod(mt, mt, mp);
  }
  mpz_mul(my, mt, mt);  mpz_mod(my, my, mp);
  r=m;
  mpz_mul(mx, mx, mt); mpz_mod(mx, mx, mp);
  mpz_mul(mb, mb, my); mpz_mod(mb, mb, mp);
  goto C1_5_1_S3_mp;

C1_5_1_DONE_mp:

  return 1;
  
}

#ifndef LONG64
void mpz_set_si64( mpz_t rop, s64 a)
{
  if( a < 0 )
  {
    a = -a;
    mpz_import( rop, 1, -1, sizeof(s64), 0, 0, &a );
    mpz_neg( rop, rop );
  }
  else
  {
    mpz_import( rop, 1, -1, sizeof(s64), 0, 0, &a );
  }
}

/****************************************************/
void mpz_mul_si64( mpz_t rop, mpz_t op1, s64 a)
{
  static int initialized = 0;
  static mpz_t temp;

  if (initialized == 0) {
    mpz_init(temp);
    initialized = 1;
  }
  mpz_set_si64( temp, a);
  mpz_mul( rop, op1, temp );
}
#endif /* !LONG64 */

/****************************************************/
int mpz_evalF(mpz_t res, s64 a, s32 b, mpz_poly f)
/****************************************************/
/* Input: integers a, b with b>0, and the poly f.   */
/* Output: res = F(a, b).                           */
/****************************************************/
{ static int initialized=0;
  static mpz_t apow, tmp;
  int    i;

  if (!(initialized)) {
    mpz_init(apow);
    mpz_init(tmp);
    initialized=1;
  }

  mpz_set(res, &f->coef[0]);
  mpz_set_si64(apow, a);
  for (i=1; i<= f->degree; i++) {
    mpz_mul_si(res, res, b);
    mpz_mul(tmp, apow, &f->coef[i]);
    mpz_add(res, res, tmp);
    mpz_mul_si64(apow, apow, a);
  }
  return 0;
}

/****************************************************/
double mpz_evalF_d(double x, double y, mpz_poly f)
/****************************************************/
/* Evaluate the homogenous polynomial               */
/*    F(x,y) = y^d * f(x/y)                         */
/* at the given point.                              */
/****************************************************/
{ static int initialized=0;
  static mpf_t apow, tmp, res, mpX, mpY;
  int    i;

  if (!(initialized)) {
    mpf_init2(apow, 256); /* Should be more than enough precision. */
    mpf_init2(tmp, 256); mpf_init2(res, 256);
    mpf_init2(mpX, 256); mpf_init2(mpY, 256);
    initialized=1;
  }

  mpf_set_z(res, &f->coef[0]);
  mpf_set_d(apow, x);
  mpf_set_d(mpX, x);
  mpf_set_d(mpY, y);
  for (i=1; i<= f->degree; i++) {
    mpf_mul(res, res, mpY);
    mpf_set_z(tmp, &f->coef[i]);
    mpf_mul(tmp, apow, tmp);
    mpf_add(res, res, tmp);
    mpf_mul(apow, apow, mpX);
  }
  return mpf_get_d(res);
}

/******************************************************/
INLINE int fplog_evalF(s32 a, s32 b, nfs_fb_t *FB)
/******************************************************/
/* Input: integers a, b, with b>0, and the FB.        */
/* Return value: floor[ log|F(a,b)|/FB->log_alb ].    */
/******************************************************/
{ static int initialized=0;
  static mpz_t norm;

  if (!(initialized)) {
    mpz_init(norm);
    initialized=1;
  }

  mpz_evalF(norm, a, b, FB->f);
  assert(mpz_sizeinbase(norm, 2)/(M_LOG2E*FB->log_alb) <= (double)INT_MAX);
  return (int)(mpz_sizeinbase(norm, 2)/(M_LOG2E*FB->log_alb));
}

        
/******************************************************/
INLINE int fplog_mpz(mpz_t k, double log_of_base)
/******************************************************/
/* Input: An mpz_t integer k.                         */
/* Return value: round[ log_2(k)]                     */
/******************************************************/
{
  assert(mpz_sizeinbase(k, 2)/(M_LOG2E*log_of_base) <= (double)INT_MAX);
  return (int)(mpz_sizeinbase(k, 2)/(M_LOG2E*log_of_base));
}

/******************************************************/
INLINE int fplog(s32 k, double log_of_base)
/******************************************************/
{ int t;
        
  t = (int)(0.5 + (log((double)k)/log_of_base));
  return t;
}

/***************************************************************/
void mpz_nearest_int(mpz_t q, mpz_t num, mpz_t den)
/* Return the nearest integer to num/den.                      */
{ static mpz_t tmpn, tmpd;
  static int initialized=0;

  if (!(initialized)) {
    mpz_init(tmpn);
    mpz_init(tmpd);
    initialized = 1;
  }
  mpz_set(tmpn, num);
  mpz_set(tmpd, den);
  if (mpz_sgn(tmpn) < 0)
    mpz_neg(tmpn, tmpn);
  if (mpz_sgn(tmpd) < 0)
    mpz_neg(tmpd, tmpd);
  mpz_mul_2exp(tmpn, tmpn, 1);
  mpz_add(tmpn, tmpn, tmpd);
  mpz_mul_2exp(tmpd, tmpd, 1);

  mpz_tdiv_q(q, tmpn, tmpd);
  if (mpz_sgn(num) != mpz_sgn(den))
    mpz_neg(q, q);
}

/********************************************/
void msgLog(char *fName, char *fmt, ...)
{ char tmp[4*MAX_MSG_SIZE], timestr[64], *f;
  va_list ap;
  FILE   *fp;
  time_t  now;
  
  /* Grab the message that was passed: It is VERY IMPORTANT */
  /* to have the va_end(ap) after!                          */
  va_start(ap, fmt);
  vsnprintf(tmp, 4*MAX_MSG_SIZE, fmt, ap);
  va_end(ap);

  time(&now); 
  strftime(timestr, 32, "[%m/%d %X]", localtime(&now));
  if ((fName!=NULL)&&fName[0]) f=fName;
  else f=GGNFS_LOG_NAME;
  if ((fp = fopen(f, "a"))) {
    fprintf(fp, "%s %s\n", timestr, tmp);
    fclose(fp);
  }
}

/********************************************/
void printTmp(char *fmt, ...)
{ char    tmp[4*MAX_MSG_SIZE];
  static  size_t lastSize=0;
  va_list ap;
  size_t     i;  

  /* Grab the message that was passed: It is VERY IMPORTANT */
  /* to have the va_end(ap) after!                          */
  va_start(ap, fmt);
  vsnprintf(tmp, 4*MAX_MSG_SIZE, fmt, ap);
  va_end(ap);

  printf("%s", tmp);
  for (i=strlen(tmp); i<=lastSize; i++)
    printf(" ");
  printf("\r");
  lastSize = strlen(tmp);
  fflush(stdout);
}


/********************************************/
void mpz_fact_init(mpz_fact_t *F)
{
  mpz_init(F->N);
  mpz_init(F->unfactored);
  F->factors = NULL;
  F->exponents = NULL;
  F->size = 0;
}
/********************************************/
void mpz_fact_clear(mpz_fact_t *F)
{
  mpz_clear(F->N);
  mpz_clear(F->unfactored);
  if (F->factors != NULL) {
    free(F->factors); F->factors = NULL;
  }
  if (F->exponents != NULL) {
    free(F->exponents); F->exponents = NULL;
  }
  F->size = 0;
}

/************************************************************/
void mpz_fact_add_factor(mpz_fact_t *F, mpz_t factor, int exponent)
/************************************************************/
/* Given a sorted (possibly empty) list of factors, add the */
/* new factor in sorted order.                              */
/************************************************************/
{ int i, loc;
  int c;
  __mpz_struct *tmpFact;
  int          *tmpExp;
	
  for (i=0, loc=F->size; (unsigned int)i<F->size; i++) {
    c = mpz_cmp(factor, &F->factors[i]);
    if (c==0) {
      F->exponents[i] += exponent;
      return;
    } else if (c < 0)
      loc = i;
  }

  tmpFact = (__mpz_struct *)malloc((F->size+1)*sizeof(__mpz_struct));
  tmpExp = (int *)malloc((F->size+1)*sizeof(int));

  
  for (i=F->size; i>=0; i--)
    mpz_init(&tmpFact[i]);

  for (i=0; i<loc; i++) {
    mpz_set(&tmpFact[i], &F->factors[i]);
    tmpExp[i] = F->exponents[i];
  }
  mpz_set(&tmpFact[loc], factor);
  tmpExp[loc] = exponent;
  for (i=loc+1; (unsigned int)i < (F->size+1); i++) {  
    mpz_set(&tmpFact[i], &F->factors[i-1]);
    tmpExp[i] = F->exponents[i-1];
  }
  for (i=0; (unsigned int)i<F->size; i++)
    mpz_clear(&F->factors[i]);
#ifdef _OLD  
  if (F->factors != NULL)  
    free(F->factors);
  if (F->exponents != NULL) 
    free(F->exponents);
#else
  if (F->size > 0) {
    free(F->factors); free(F->exponents);
  }
#endif

  F->factors = tmpFact;
  F->exponents = tmpExp;
  F->size += 1;
  return;
}

/************************************************************/
#ifdef _MAXIMAL_ORDER_ONLY
const int  _numLevels=7;
#else
const int  _numLevels=1;
#endif
const s32 _B1Sizes[]={2000, 5000, 10000, 50000, 250000, 1000000, 3000000};
const u32 _numCurves[]={50, 150, 200, 400, 600, 1000, 1000};

static int mpz_fact_factorRealWork_rec(mpz_fact_t *F, int doRealWork,
                                       u32 numSmallP, u32 *smallP,
                                       mpz_t tmp1, mpz_t tmp2, char str[256])
/****************************************************/
/* Recursive routine meant to be called from        */
/* mpz_fact_factorEasy().  Does a power check and,  */
/* if doRealWork is positive, ECM.  Assumes that    */
/* trial division and DISC_P_DIV factors have       */
/* already been divided out (so no need to repeat)  */
/****************************************************/
{
  int        top_e, e, giveUp, retVal, iter, level;
  mpz_t      s;
  s32        B1;
  u32        i;
  double     B2;

  if (doRealWork < 0) return -1;

  mpz_init(s);

  giveUp = 0;
  iter=10;
  level=0;
  B1 = _B1Sizes[0];
  B2 = 0;
  mpz_set_ui(s, 0);
  top_e = 1;
  while ((!giveUp) && mpz_cmp_ui(F->unfactored, 1)) {
    if (mpz_probab_prime_p(F->unfactored, 20)) {
      msgLog(GGNFS_LOG_NAME, "DISC_P_DIV=%s",
             mpz_get_str(str, 10, F->unfactored));
      mpz_fact_add_factor(F, F->unfactored, 1);
      mpz_set_ui(F->unfactored, 1);
    } else {
      /** Do a power check; avoid mpz_perfect_power_p() because **/
      /** it redoes trial division (as of GMP 4.1.4)            **/
      e = 1;
      while (mpz_perfect_square_p(F->unfactored)) {
        e *= 2;
        mpz_sqrt(F->unfactored, F->unfactored);
      }
      for (i=1; i<numSmallP; i++) {
        while (mpz_root(tmp1, F->unfactored, smallP[i])) {
          e *= smallP[i];
          mpz_set (F->unfactored, tmp1);
        }
        if (mpz_cmp_ui(tmp1, smallP[numSmallP-1]) < 0) break;
      } 
      if (e > 1) {
        top_e *= e;
        if (mpz_probab_prime_p(F->unfactored, 20)) {
//        printf("found factor with multiplicity %ld\n", e);
          msgLog(GGNFS_LOG_NAME, "DISC_P_DIV=%s",
                 mpz_get_str(str, 10, F->unfactored));
          mpz_fact_add_factor(F, F->unfactored, top_e);
          mpz_set_ui(F->unfactored, 1);
        }
      }
      /** Now do ECM on the remaining part **/
      if (mpz_cmp_ui(F->unfactored, 1) && (doRealWork > 0)) {
        mpz_set(tmp2, F->unfactored);
        i=0;
        do {
          retVal = ecmFactor(tmp1, tmp2, B1, B2, iter, s);
          i++;
          if ((i%10)== 0) {
            printf("Attempt %" PRId32 " / %" PRId32 " for: ", i, _numCurves[level]); 
            mpz_out_str(stdout, 10, tmp2); 
            printf("\n");
          }
          /* if ECM found input number, try again */
          if ((retVal > 0) && (mpz_cmp(tmp1, F->unfactored) == 0)) retVal = 0;
        } while ((retVal==0) && (i<_numCurves[level]));
        if (retVal <= 0) {
          if (++level < _numLevels)
            B1 = _B1Sizes[level];
	  else
            giveUp = 1;
        } else {
          mpz_divexact(F->unfactored, F->unfactored, tmp1);
          if (mpz_probab_prime_p(tmp1, 20)) {
            e = 1;
            mpz_mod(tmp2, F->unfactored, tmp1);
            while (mpz_sgn(tmp2)==0) {
              e++;
              mpz_divexact(F->unfactored, F->unfactored, tmp1);
              mpz_mod(tmp2, F->unfactored, tmp1);
            }
            msgLog(GGNFS_LOG_NAME, "DISC_P_DIV=%s", mpz_get_str(str, 10, tmp1));
            mpz_fact_add_factor(F, tmp1, top_e*e);
	  } else {
            mpz_fact_t TmpF;

            mpz_fact_init(&TmpF);
            mpz_set(TmpF.N, tmp1);
            mpz_set(TmpF.unfactored, tmp1);
            TmpF.sign = 1;
	    retVal = mpz_fact_factorRealWork_rec(&TmpF, doRealWork, numSmallP,
                                                 smallP, tmp1, tmp2, str);
            if (retVal==0) {
              for (i=0; i<TmpF.size; i++)
                mpz_fact_add_factor(F, &TmpF.factors[i],
                                    top_e*TmpF.exponents[i]);
            }
            else
              doRealWork = 0; /* give up on ECM, but do one more power check */
            mpz_fact_clear(&TmpF);
          }
        }  
      } /* end ECM block */
      else
        giveUp = 1; /* no more power checks, no more ECM */ 
    }
  }

  mpz_clear(s);
  if (mpz_cmp_ui(F->unfactored, 1)==0)
    return 0;
  return -1;

}

int mpz_fact_factorEasy(mpz_fact_t *F, mpz_t N, int doRealWork)
/************************************************************/
/* Attempt to factor N using trial division and ECM.        */
/* 'F' should already be initialized.                       */
/* prandseed() should be used first, if you want to set it. */
/* Return value: 0 <==> factored completely.                */
/* If doRealWork=0, we will not try too hard. If it is -1   */
/* we will do only trial division and divisors in ggnfs.log */
/************************************************************/
{ int        e;
  u32        numSmallP, *smallP;
  mpz_t      tmp1, tmp2;
  u32        i;
  FILE       *fp;
  char       str[256], *loc;
  
  mpz_init(tmp1);
  mpz_init(tmp2);

  mpz_set(F->N, N);
  mpz_set(F->unfactored, N);
  if (F->factors != NULL) {
    for (i=0; i<F->size; i++)
      mpz_clear(&F->factors[i]); 
    free(F->factors);
  }
  if (F->exponents != NULL) 
    free(F->exponents);
  F->factors = NULL;
  F->exponents = NULL;
  F->size = 0;

  F->sign = 1;
  if (mpz_sgn(F->unfactored) < 0) {
    F->sign = -1;
    mpz_abs(F->unfactored, F->unfactored);
  } 

  
  /** Do trial division. **/
  numSmallP = 50000;
  smallP = getPList(&numSmallP);
//  printf("Doing trial division with %ld small primes.\n", numSmallP);
  i=0;
  while ((i<numSmallP) && mpz_cmp_ui(F->unfactored, 1)) {
    if (mpz_tdiv_q_ui(tmp1, F->unfactored, smallP[i])==0) {
      e = 1;
      mpz_set(F->unfactored, tmp1);
      while (mpz_tdiv_q_ui(tmp1, F->unfactored, smallP[i])==0) {
        e++;
        mpz_set(F->unfactored, tmp1);
      }
      mpz_set_ui(tmp1, smallP[i]);
      mpz_fact_add_factor(F, tmp1, e);
      if (mpz_probab_prime_p(F->unfactored, 20)) {
        mpz_fact_add_factor(F, F->unfactored, 1);
	mpz_set_ui(F->unfactored, 1);
      }
    }
    i++;
  }
  

  /***************************************************/
  /* First, scan the log file to see if we'd already */
  /* factored this number. i.e., just scan through   */
  /* for DISC_P_DIV=<n> lines, and see if any of them*/
  /* are factors.                                    */
  /***************************************************/
  if (mpz_cmp_ui(F->unfactored, 1) && (fp = fopen(GGNFS_LOG_NAME, "r"))) {
    while (!(feof(fp)) && mpz_cmp_ui(F->unfactored, 1)) {
      fgets(str, 256, fp);
      loc = strstr(str, "DISC_P_DIV=");
      if (loc != NULL) {
        mpz_set_str(tmp1, &loc[11], 10);
        if (mpz_probab_prime_p(tmp1, 20)) {
          e = 0;
          mpz_mod(tmp2, F->unfactored, tmp1);
          while (mpz_sgn(tmp2)==0) {
            e++;
            mpz_divexact(F->unfactored, F->unfactored, tmp1);
            mpz_mod(tmp2, F->unfactored, tmp1);
          }
          if (e)
	  {
	    mpz_fact_add_factor(F, tmp1, e);
            if (mpz_cmp_ui(F->unfactored, 1)==0) {
              msgLog(GGNFS_LOG_NAME, "Discriminant factorization re-read from log file.");
            }
	  }
        }
      }
    }
    fclose(fp);
  }

  if (mpz_cmp_ui(F->unfactored, 1) && (doRealWork >= 0)) {
    mpz_fact_factorRealWork_rec(F, doRealWork, numSmallP, smallP, tmp1, tmp2, str);
  }
  free(smallP);  
  mpz_clear(tmp1);
  mpz_clear(tmp2);
  if (mpz_cmp_ui(F->unfactored, 1)==0)
    return 0;

  return -1;

}
	

/*****************************************************************/
int mpz_fact_check(mpz_fact_t *D, int doRealWork)
/*****************************************************************/
/* Check that the factorization is correct and complete. If not, */
/* we will attempt to finish it.                                 */
/* Return value:  0, D was/is ok.                                */
/*             <> 0, D is not complete.                          */
/*****************************************************************/
{ unsigned int i; 
  int        e, retVal;
  mpz_t      tmp1, tmp2;
  mpz_fact_t Tmp;
	
  mpz_init(tmp1);
  mpz_init(tmp2);
  
  mpz_abs(D->unfactored, D->N);
  for (i=0; i<D->size; i++) {
    e = 0;
    mpz_tdiv_qr(tmp1, tmp2, D->unfactored, &D->factors[i]);
    while (mpz_sgn(tmp2)==0) {
      e++;
      mpz_set(D->unfactored, tmp1);
      mpz_tdiv_qr(tmp1, tmp2, D->unfactored, &D->factors[i]);
    }
    D->exponents[i] = e;
  }
  if (mpz_cmp_ui(D->unfactored, 1)==0) {
    mpz_clear(tmp1); mpz_clear(tmp2);
    return 0;
  }
  mpz_fact_init(&Tmp);
  retVal = mpz_fact_factorEasy(&Tmp, D->unfactored, doRealWork);
  if (retVal == 0) {
    for (i=0; i<Tmp.size; i++)
      mpz_fact_add_factor(D, &Tmp.factors[i], Tmp.exponents[i]);
    mpz_fact_clear(&Tmp);
    mpz_clear(tmp1); mpz_clear(tmp2);
    return mpz_fact_check(D, doRealWork);
  } 
  /* Even if a failure occurred, still add the factors we did find. */
  for (i=0; i<Tmp.size; i++)
    mpz_fact_add_factor(D, &Tmp.factors[i], Tmp.exponents[i]);

  mpz_clear(tmp1); mpz_clear(tmp2);
  mpz_fact_clear(&Tmp);
  return -1;
}

/*****************************************************************/
int mpz_fact_removeSF(mpz_fact_t *S, mpz_fact_t *F)
/*****************************************************************/
/* Input: F is the discriminant of a poly.                       */
/* Output: S will be F with the fundamental discriminant part    */
/*         removed.                                              */
/* Note: Right now, we do only S <-- F/X, where X is the         */
/*       maximal squarefree divisor of S.                        */
/*****************************************************************/
{ unsigned int   i;
  int j, e, twoLoc=-1;

  mpz_fact_clear(S);
  mpz_fact_init(S);

  for (i=0; i<F->size; i++) {
    e = F->exponents[i];
    if (e > 1) {
      mpz_fact_add_factor(S, &F->factors[i], e - e%2);
    } else if (mpz_cmp_ui(&F->factors[i], 2) && (e >= 1)) {
      twoLoc = i;
    }
  }
  /* Later, we should do something with the two's!!!! */


  mpz_set_ui(S->N, 1);
  for (i=0; i<S->size; i++) {
    e = S->exponents[i];
    for (j=0; j<e; j++)
      mpz_mul(S->N, S->N, &S->factors[i]);
  }
  mpz_set_ui(S->unfactored, 1);
  
  return 0;
}
  
  
/*****************************************************************/
void initNF(nf_t *N)
{ int i;

  mpz_poly_init(N->f);
  mpz_poly_init(N->T);
  N->degree = 0;
  N->W = (mpz_mat_t *)malloc(sizeof(mpz_mat_t));
  N->W_inv = (mpz_mat_t *)malloc(sizeof(mpz_mat_t));
  mpz_mat_init2(N->W, MAXPOLYDEGREE, MAXPOLYDEGREE);
  mpz_mat_init2(N->W_inv, MAXPOLYDEGREE, MAXPOLYDEGREE);
  mpz_init(N->W_d);
  mpz_init(N->Kdisc);
  mpz_init(N->Tdisc);
  mpz_init(N->index);
  for (i=0; i<MAXPOLYDEGREE; i++) {
    mpf_init2(N->fZeros[i].mpr, 128);  
    mpf_init2(N->fZeros[i].mpi, 128);  
    mpf_init2(N->TZeros[i].mpr, 128);  
    mpf_init2(N->TZeros[i].mpi, 128);  
  }
  N->sPrimes = NULL;
  N->v_cd_sPrimes = NULL;
  N->numSPrimes = 0;
  N->Mt = NULL;
  mpz_poly_init(N->Sk_ib);
}

/*****************************************************************/
void clearNF(nf_t *N)
{ int i;

  mpz_poly_clear(N->f);
  mpz_poly_clear(N->T);
  N->degree = 0;
  mpz_mat_clear(N->W);
  free(N->W);
  mpz_mat_clear(N->W_inv);
  free(N->W_inv);
  mpz_clear(N->W_d);
  mpz_clear(N->Kdisc);
  mpz_clear(N->Tdisc);
  mpz_clear(N->index);
  for (i=0; i<MAXPOLYDEGREE; i++) {
    mpf_clear(N->fZeros[i].mpr);  
    mpf_clear(N->fZeros[i].mpi);  
    mpf_clear(N->TZeros[i].mpr);  
    mpf_clear(N->TZeros[i].mpi);  
  }
  mpz_mat_clear(N->W);
  mpz_poly_clear(N->Sk_ib);
}

/*****************************************************************/
void initIdeal(prime_id_t *I)
{
  mpz_init(I->p);
  mpz_poly_init(I->alpha);
  mpz_poly_init(I->beta);
  mpz_mat_init(&I->betaMat);
  I->e = I->f = 0;
}

/*****************************************************************/
void clearIdeal(prime_id_t *I)
{
  mpz_clear(I->p);
  mpz_poly_clear(I->alpha);
  mpz_poly_clear(I->beta);
  mpz_mat_clear(&I->betaMat);
  I->e = I->f = 0;
}


/*********************************************************/
int cmpCplx(const void *a, const void *b)
/*********************************************************/
/* Compare so that:
   (1) A real number is always < a non-real number.
   (2) If 'a' and 'b' are both real (non-real), then
       use a lex order on (real part, complex part).
*/
{ nfs_complex_t *A = (nfs_complex_t *)a, *B=(nfs_complex_t *)b;

  if ((fabs(A->i) < 0.00000001) && (fabs(B->i) > 0.00000001))
    return -1;
  if ((fabs(A->i) > 0.00000001) && (fabs(B->i) < 0.00000001))
    return 1;

  if (fabs(A->r - B->r) > 0.00000001) {
    if (A->r > B->r)
      return 1;
    return -1;
  }
  if (fabs(A->i - B->i) > 0.00000001) {
    if (A->i > B->i)
      return 1;
    return -1;
  }
  return 0;
}


/*********************************************************/
void reorderRoots(nfs_complex_t *z, int n)
/*********************************************************/
{
  qsort(z, n, sizeof(nfs_complex_t), cmpCplx);
}

/*********************************************************/
void reorderRoots2(nfs_complex_t *z, int n)
/*********************************************************/
/* Reorder an array of complex numbers so that:          */
/*  (1) Real ones come first, followed by complex.       */
/*  (2) If there are exactly 'r' real and 's' complex    */
/*      values, then z[r+i] is the conjugate of z[r+i+s].*/
//********************************************************/
{ int    i, j, r, s;
  double t;
  mpf_t  mt;
                                                                                                                             
  mpf_init2(mt, mpf_get_prec(z[0].mpr));
  r = 0;
  for (i=0; i<n; i++) {
    if (z[i].i < 0.000000000001) {
      mpf_set(mt, z[r].mpr); mpf_set(z[r].mpr, z[i].mpr); mpf_set(z[i].mpr, mt);
      mpf_set(mt, z[r].mpi); mpf_set(z[r].mpi, z[i].mpi); mpf_set(z[i].mpi, mt);
      t = z[r].r; z[r].r = z[i].r; z[i].r = t;
      t = z[r].i; z[r].i = z[i].i; z[i].i = t;
      r++;
    }
  }
  s = (n-r)/2;
  if (s==0) {
    mpf_clear(mt);
    return ;
  }
  for (i=0; i<s; i++) {
    /* Find the zero that is conjugate to z[r+i]. */
    j=1;
    while (((r+i+j)<n) && (fabs(z[r+i].r - z[r+i+j].r) > 0.00000000001))
      j++;
    if (j==s) {
      fprintf(stderr, "reorderZeros(): Error - could not find matching conjugate root!\n");
      mpf_clear(mt);
      return ;
    }
    /* Now swap z[r+i+j] and z[r+i+s]. */
    mpf_set(mt, z[r+i+j].mpr); mpf_set(z[r+i+j].mpr, z[r+i+s].mpr); mpf_set(z[r+i+s].mpr, mt);
    mpf_set(mt, z[r+i+j].mpi); mpf_set(z[r+i+j].mpi, z[r+i+s].mpi); mpf_set(z[r+i+s].mpi, mt);
    t = z[r+i+j].r; z[r+i+j].r = z[r+i+s].r; z[r+i+s].r = t;
    t = z[r+i+j].i; z[r+i+j].i = z[r+i+s].i; z[r+i+s].i = t;
  }
  mpf_clear(mt);
}                                          

/**************************************************************/
int cmpS32s(const void *a, const void *b)
/**************************************************************/
{ s32 *A = (s32 *)a, *B=(s32 *)b;

  if (*A < *B) return -1;
  if (*A > *B) return 1;
  return 0;
}

/**************************************************************/
int cmpU32s(const void *a, const void *b)
/**************************************************************/
{ u32 *A = (u32 *)a, *B=(u32 *)b;

  if (*A < *B) return -1;
  if (*A > *B) return 1;
  return 0;
}

/****************************************************/
double L_n(double n, double a, double c)
/****************************************************/
/* Computes exp(c * (log n)^a * (loglog n)^(1-a) ). */
/****************************************************/
{ double t1, t2;

  t1 = log(n);
  t2 = log(t1);
  t1 = pow(t1, a);
  t2 = pow(t2, 1.0-a);
  return exp(c*t1*t2);
}

/*************************************************************/
s32 removeS32Pairs(s32 *L, s32 size)
/*************************************************************/
/* Given a list of integers, remove duplicates in pairs, so  */
/* that if there are an even number of occurrences of some   */
/* integer, none will be left. If there are an odd number of */
/* occurrences, exactly one will be left.                    */
/*************************************************************/
{ s32 i, j, newEndLoc;
                                                                                                        
  if (size==0) return 0;
  qsort(L, size, sizeof(s32), cmpS32s);
  /* Now, remove pairs of duplicates: */
  newEndLoc=0;
  i=0;
  while (i<size) {
    /* How many copies of L[i] are there? */
    j=i+1;
    while ((j<size) && (L[j]==L[i]))
      j++;
    if ((j-i)%2 == 0) {
      /* An even number of them, so omit them altogether. */
      i=j;
    } else {
      /* Copy one of them over. */
      L[newEndLoc++] = L[i];
      i=j;
    }
  }
  return newEndLoc;
}

/*************************************************************/
int cmp_lpair_t(const void *a, const void *b)
/*************************************************************/
{ lpair_t *A=(lpair_t *)a, *B=(lpair_t *)b;
                                                                                                  
  if (A->x < B->x) return -1;
  if (A->x > B->x) return 1;
  if (A->y < B->y) return -1;
  if (A->y > B->y) return 1;
  return 0;
}

/************************************************************/
int mallocReport() 
/************************************************************/
/* Returns the amount of RAM currently malloc'd by this process in megs. */
{ 
#ifdef MALLOC_REPORTING
  struct mallinfo mallInf;
  
  mallInf = mallinfo();
  /* Strange - as far as I can tell, uordblks is supposed to
     be the total. But at least in my version of glibc (2.3.3-27,
     Fedora core 2), I need to add these like so:
  */
  return (mallInf.uordblks + mallInf.hblkhd)/1048576;
#else
  return 0;
#endif
}


static int mi_maxHeapUseage=0;
static int mi_errs=0;
static int mi_allocs=0, mi_reallocs=0;
/*****************************************************************/
void *lxmalloc(size_t n, int fatal)
/*****************************************************************/
{ void *p;
  int   hu;

  mi_allocs++;
  p = malloc(n);
#ifdef MALLOC_DEBUG
  FILE* f;
  f = fopen ("m_info.txt","a+");
  fprintf(f,"malloc: %ld\n",n);
  fclose(f);
#endif
  if (p==NULL) {
    mi_errs++;
    msgLog("", "Memory allocation error (%" PRIu32 " bytes requested).", n);
    printf("Memory allocation error (%" PRIu32 " bytes requested).", (u32)n);
    if (!(fatal)) return NULL;
    printf("Fatal error. Terminating...\n");
    exit(-1);
  }
  hu = mallocReport();
  mi_maxHeapUseage = MAX(mi_maxHeapUseage, hu);
  return p;
}
/*****************************************************************/
void *lxcalloc(size_t n, int fatal)
/*****************************************************************/
{ void *p;
  int   hu;

  mi_allocs++;
  p = malloc(n);
#ifdef MALLOC_DEBUG
  FILE* f;
  f = fopen ("m_info.txt","a+");
  fprintf(f,"calloc: %ld\n",n);
  fclose(f);
#endif
  if (p==NULL) {
    mi_errs++;
    msgLog("", "Memory allocation error (%" PRIu32 " bytes requested).", n);
    printf("Memory allocation error (%" PRIu32 " bytes requested).", (u32)n);
    if (!(fatal)) return NULL;
    printf("Fatal error. Terminating...\n");
    exit(-1);
  }
  hu = mallocReport();
  mi_maxHeapUseage = MAX(mi_maxHeapUseage, hu);
  memset(p, 0x00, n);
  return p;
}

/*****************************************************************/
void *lxrealloc(void *x, size_t n, int fatal)
/*****************************************************************/
{ void *p;
  int   hu;

  mi_reallocs++;
#ifdef MALLOC_DEBUG
  FILE* f;
  f = fopen ("m_info.txt","a+");
  fprintf(f,"realloc: %ld\n",n);
  fclose(f);
#endif
  p = realloc(x, n);
  if (p==NULL) {
    mi_errs++;
    msgLog("", "Memory allocation error (%" PRIu32 " bytes requested).", n);
    printf("Memory allocation error (%" PRIu32 " bytes requested).", (u32)n);
    if (!(fatal)) return NULL;
    printf("Fatal error. Terminating...\n");
    exit(-1);
  }
  hu = mallocReport();
  mi_maxHeapUseage = MAX(mi_maxHeapUseage, hu);
  return p;
}

/*****************************************************************/
int getHeapStats(int *maxUseage, int *errs, int *allocs, int *reallocs)
{
  *maxUseage = mi_maxHeapUseage;
  *errs = mi_errs;
  *allocs = mi_allocs;
  *reallocs = mi_reallocs;
  return mi_maxHeapUseage;
}

/*****************************************************************/
void logHeapStats()
{
  msgLog("", "Max heap usage: %d MB", mi_maxHeapUseage);
  msgLog("", "malloc/realloc errors: %d ", mi_errs);
  msgLog("", "total malloc's : %d ", mi_allocs);
  msgLog("", "total realloc's: %d ", mi_reallocs);

}
