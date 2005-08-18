/*********************************************************/
/* polyselect.c                                          */
/* Copyright 2004, Chris Monico.                         */
/*********************************************************/
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
//#include <math.h>
#include <sys/stat.h>
#include "ggnfs.h"
#include "prand.h"

#define START_MSG \
"\n"\
" __________________________________________________________ \n"\
"|      This is the polyselect program for GGNFS.           |\n"\
"| Version: %-25s                            |\n"\
"| This program is copyright 2004, Chris Monico, and subject|\n"\
"| to the terms of the GNU General Public License version 2.|\n"\
"|__________________________________________________________|\n"\
"* This program is neither complete nor necessarily stable! *\n"\
"*         It is still very much under development          *\n"\
"************************************************************\n"

/*#define _VERY_LOUD */


/* This is the number of pieces into which we divide the half circle [0,\pi],
   to estimate F(cos(\theta), \sin(theta)) on the unit circle. Murphy denotes
   it by K, and takes it to be 1000 which seems reasonable. */
#define HALF_CIRCLE_PIECES 1000

/* An upper bound for primes to be used in computing the rating of
   the polynomials: If you change this, be aware that it changes
   the rating scheme so that an E(*,*) score with a CONT_PRIME_BOUND
   of 2000, say, cannot be (fairly) compared to score obtained with
   3000.
*/
#define CONT_PRIME_BOUND 2000
#define MAX_CONTPRIMES 1000  /* This can be larger than pi(2000), so be safe! */

/* Primes upto SAMPLE_PRIME_BOUND will have their content estimated by
   actual computation, for more accuracy. */
#define SAMPLE_PRIME_BOUND 100

/* CONT_SAMPLE_SIZE is how many samples to gather for the purpose of estimating
   cont_p from computed values of F(x,y). Even with 1000 samples, there is
   a good deal of fluctuation!
*/
#define DEFAULT_CONT_SAMPLE_SIZE 2000


#define MIN_A (-0x7FFFFFFF)
#define MAX_A (0x7FFFFFFF)
#define MAX_B (0x7FFFFFFF)

#define MAX_REASONABLE_SKEW 2000


#define MAX_NFS_POLY_DEGREE 6
#define MAX_MV_POLY_SIZE  (8+MAX_NFS_POLY_DEGREE+1)

/* This is nowhere near done yet! It is still very much a work in
   progress, but hopefully it will become a good poly selection
   program. But for now, it is just a collection of useful functions
   along with some code to test them. */

/* This is the single term of a multivariate polynomial.
   The variables are fixed. In the notation of Murphy's thesis, they are:
   (The last DEGREE of them are not really variables - they're more like parameters
   that are convenient to treat the same way some times; they correspond to coefficients
   of the starting polynomial.)
     exps[0] <--> c0
     exps[1] <--> c1
     exps[2] <--> c2
     exps[3] <-->  t
     exps[4] <-->  s
     exps[5] <-->  x
     exps[6] <-->  y
     exps[7] <-->  m
     exps[8] <--> a0 
      ...
     exps[7+k] <--> ak

   Where f_m(x) = a0 + a1*x + ... + a_d*x^d. 
   A `s32' coefficient will suffice for our purposes for 
   number field polynomials upto degree 6 (verified via Maple), 
   but will almost certainly overflow 32 bits for degree 7 (I didn't bother checking,
   for obvious reasons.). To keep integration simple, though, I'm going to 
   use double coefficients which should also suffice. */
#define C0_LOC 0
#define C1_LOC 1
#define C2_LOC 2
#define T_LOC  3
#define S_LOC  4
#define X_LOC  5
#define Y_LOC  6
#define M_LOC  7
#define A_START 8

#define DEFAULT_J0 100
#define DEFAULT_J1 10

/* If we are enumerating polynomials with a specified divisor of c_d,
   this is the default number of them to go through.
*/
#define DEFAULT_ENUMSIZE 1000000

typedef struct {
  double coef; 
  float  exps[MAX_MV_POLY_SIZE];
} mv_term_t;

typedef struct {
  mv_term_t *terms;
  int        numTerms;
} mv_poly_t;

typedef struct {
  mpf_t    c0, c1, c2, t, s;
  mpz_t    m, n;
  mpz_poly f;
  double   logSize;
  double   score;
} Iparam_t;

#define  DEFAULT_ALLNAME   "all.poly"
#define  DEFAULT_BESTNAME  "best.poly"
#define USAGE \
"[OPTIONS]\n"\
"--help              : Show this help and exit.\n"\
"-if <file>          : <file> has the number and poly. degree\n"\
"-score <file>       : Just score the polynomial in this file; do nothing else.\n"\
"-css   <int>        : Sample size for scoring.\n"\
"     The following are accepted options from an input file:\n"\
"n: <mp integer>     : The number being factored.\n"\
"name: <string>      : A descriptive name string for the number.\n"\
"deg: <int>          : Desired polynomial degree.\n"\
"af: <file>          : Append the current best polynomial to <file>\n"\
"bf: <file>          : Keep the best poly so far in <file>\n"\
"seed: <int>         : seed the PRNG with <int>\n"\
"j0: <int>           : Size of Murphy's J0 for sieving-for-root-properties.\n"\
"j1: <int>           : Size of Murphy's J1 for sieving-for-root-properties.\n"\
"                      ( J1 <<  J0 is a good idea).\n"\
"maxs1: <float>      : Throw out polynomials whose stage 1 'I' score is\n"\
"                      larger than this.\n"\
"maxskew: <float>    : Max acceptable skew.\n"\
"lc1: <float>        : Max relative size of leading coefficient (in (0,1)).\n"\
"lcp: <int>          : Number of primes to choose from for l.c. divisors.\n"\
"leave: <int>        : Leave this many bits of the l.c. to chance.\n"\
"lcd: <long>         : Consider only l.c. divisible by this integer.\n"\
"enum: <lcd>         : Test polynomials w/ l.c. e0*lcd, (e0+1)*lcd,..., \n"\
"                      e1*lcd.\n"\
"e0: <big int>       : See `enum'.\n"\
"e1: <big int>       : See `enum'.\n"\
"css: <int>          : Contribution sample size (higher for more accurate scoring).\n"\
"                      If you are scoring a single polynomial only and want a more\n"\
"                      accurate score, 100000 is an acceptable choice. Otherwise,\n"\
"                      leave it alone - the default is ok within a significant digit.\n"\
"cutoff: <double>    : Keep any polynomial scoring better than this.\n"\
"examinefrac: <double> : A number in the range (0,1]. The lower this number is, the\n"\
"                        stricter the test used in looking at the 3rd coefficient. Thus,\n"\
"                        the smaller the number, the fewer intermediate polynomials which\n"\
"                        will be expensively examined.\n"


/*********** Globals *************/
int       d; /* degree of poly we're looking for. */
mv_poly_t I; /* The double integral of F^2(x,y) over the rectangle S. */
mv_poly_t dIdC0, dIdC1, dIdC2, dIdT, dIdS; /* Partial derivatives of I. */
int       maxPowC0, maxPowC1, maxPowC2, maxPowT, maxPowS, 
          maxPowM, maxPowA[MAX_NFS_POLY_DEGREE+1];
int       sieveJ0=DEFAULT_J0, sieveJ1=DEFAULT_J1;
int       lcp=7, leave=9;

double    minStage1Size=10000, userSetMin=0.0;
double    maxSkew = MAX_REASONABLE_SKEW;
s32      iteration=0, lcd=1;
s32      contSampleSize=DEFAULT_CONT_SAMPLE_SIZE;
mpz_t     lastLC, N, enumLCD, e0, e1, enumCt;
char      ifname[256], allName[256], bestName[256];
char      name[MAX_NAMESIZE], scoreName[512];
s32      seed;
double    lc0=0.05, lc1=0.5;
double    cutoff=1, examineFrac=0.5;
/*********************************/


void mv_poly_cp(mv_poly_t *dest, mv_poly_t *src);
void Iparam_init(Iparam_t *par);
void Iparam_cp(Iparam_t *dest, Iparam_t *src);

/******************************************************/
double randDouble(double min, double max)
{ double x = prand_01();

  return min + (max-min)*x;
}


/******************************************************/
void mv_term_cp(mv_term_t *dest, mv_term_t *src)
{ int i;

  dest->coef = src->coef;
  for (i=0; i<MAX_MV_POLY_SIZE; i++)
    dest->exps[i] = src->exps[i];
}

/******************************************************/
void mv_term_mult(mv_term_t *res, mv_term_t *op1, mv_term_t *op2)
{ int i;

  res->coef = op1->coef*op2->coef;
  for (i=0; i<MAX_MV_POLY_SIZE; i++)
    res->exps[i] = op1->exps[i] + op2->exps[i];
}

/******************************************************/
void mv_term_setone(mv_term_t *res)
{ int i;

  for (i=0; i<MAX_MV_POLY_SIZE; i++)
    res->exps[i] = 0;
  res->coef=1;
}

/******************************************************/
void mv_poly_simplify(mv_poly_t *res, mv_poly_t *op1)
{ int     i, j, k, l, found;
  mv_poly_t tmp;

  tmp.terms = (mv_term_t *)malloc(op1->numTerms*sizeof(mv_term_t));
  tmp.numTerms=0;
  mv_term_cp(&tmp.terms[0], &op1->terms[0]);
  i=1; j=0;
  if (fabs(tmp.terms[0].coef) > 0.0000001)
    j=1;
  while (i<op1->numTerms) {
    found=-1;
    for (k=0; (k<j)&&(found==-1); k++) {
      found = k;
      for (l=0; l<MAX_MV_POLY_SIZE; l++) 
        if (tmp.terms[k].exps[l] != op1->terms[i].exps[l])
          found=-1;
    }
    if (found > -1) {
      /* Combine these terms. */
      tmp.terms[found].coef += op1->terms[i].coef;
    } else if (fabs(op1->terms[i].coef) > 0.0000001) { /* Sloppy, but whatever. */
      /* Add this term to res. */
      mv_term_cp(&tmp.terms[j], &op1->terms[i]);
      j++;
    }
    i++;
  }
  tmp.numTerms = j;
  mv_poly_cp(res, &tmp);
  free(tmp.terms);
}

/******************************************************/
void mv_poly_cp(mv_poly_t *dest, mv_poly_t *src)
{ int i;

  if (dest->terms != NULL) free(dest->terms);
  dest->terms = (mv_term_t *)malloc(src->numTerms*sizeof(mv_term_t));
  dest->numTerms = src->numTerms;
  for (i=0; i<src->numTerms; i++)
    mv_term_cp(&dest->terms[i], &src->terms[i]);
}

/**************************************************************/
void mv_poly_add(mv_poly_t *res, mv_poly_t *op1, mv_poly_t *op2)
{ mv_poly_t tmp;
  int i;

  tmp.numTerms=op1->numTerms + op2->numTerms;
  tmp.terms = (mv_term_t *)malloc(tmp.numTerms*sizeof(mv_term_t));
  for (i=0; i<op1->numTerms; i++)
    mv_term_cp(&tmp.terms[i], &op1->terms[i]);
  for (i=0; i<op2->numTerms; i++)
    mv_term_cp(&tmp.terms[op1->numTerms+ i], &op2->terms[i]);
  mv_poly_simplify(res, &tmp);
  free(tmp.terms);
}

/**************************************************************/
void mv_poly_mul(mv_poly_t *res, mv_poly_t *op1, mv_poly_t *op2)
{ mv_poly_t tmp;
  int i, j, k;

  tmp.numTerms=op1->numTerms*op2->numTerms;
  tmp.terms = (mv_term_t *)malloc(tmp.numTerms*sizeof(mv_term_t));

  k=0;
  for (i=0; i<op1->numTerms; i++) {
    for (j=0; j<op2->numTerms; j++) {
      mv_term_mult(&tmp.terms[k], &op1->terms[i], &op2->terms[j]);
      k++;
    }
  }
  mv_poly_simplify(res, &tmp);
  free(tmp.terms);
}

//#define _LATEX_OUTPUT
/* This is actually Maple-esque output now. */
/********************************************************************************/
void mv_poly_print(FILE *fp, mv_poly_t *F)
{ int i, j;
  float e;
  double c;

#ifdef _LATEX_OUTPUT
  fprintf(fp, "\\begin{eqnarray*}\n  &=&");
#endif
  for (i=0; i<F->numTerms; i++) {
    c = F->terms[i].coef;
//    if (c != 1.0)
      fprintf(fp, "%1.3lf", c);
    for (j=0; j<MAX_MV_POLY_SIZE; j++) {
      e = F->terms[i].exps[j];
      if (e && (e!=1)) {
        if (j==C0_LOC) fprintf(fp, "*c0^(%1.1f)", e);
        else if (j==C1_LOC) fprintf(fp, "*c1^(%1.1f)", e);
        else if (j==C2_LOC) fprintf(fp, "*c2^(%1.1f)", e);
        else if (j==T_LOC) fprintf(fp, "*t^(%1.1f)", e);
        else if (j==S_LOC) fprintf(fp, "*s^(%1.1f)", e);
        else if (j==M_LOC) fprintf(fp, "*m^(%1.1f)", e);
        else if (j==X_LOC) fprintf(fp, "*x^(%1.1f)", e);
        else if (j==Y_LOC) fprintf(fp, "*y^(%1.1f)", e);
        else fprintf(fp, "*a%d^(%1.1f)", j-A_START, e);
      } else if (e) {
        if (j==C0_LOC) fprintf(fp, "*c0");
        else if (j==C1_LOC) fprintf(fp, "*c1");
        else if (j==C2_LOC) fprintf(fp, "*c2");
        else if (j==T_LOC) fprintf(fp, "*t");
        else if (j==S_LOC) fprintf(fp, "*s");
        else if (j==M_LOC) fprintf(fp, "*m");
        else if (j==X_LOC) fprintf(fp, "*x");
        else if (j==Y_LOC) fprintf(fp, "*y");
        else fprintf(fp, "*a%d", j-A_START);

      }
    }
    if (i < (F->numTerms-1)) 
      if (F->terms[i+1].coef > 0 ) fprintf(fp, " + ");
    if (i%4==0)  {
#ifdef _LATEX_OUTPUT
      fprintf(fp, "\\\\\n & & ");
#else
      fprintf(fp, "\n");
#endif
    }
  }
#ifdef _LATEX_OUTPUT
  fprintf(fp, "\\end{eqnarray*}");
#endif
  fprintf(fp, "\n");
}

/********************************************************************************/
void computeF(mv_poly_t *F, int d)
/* Compute the homogeneous polynomial F, with parameters c0, c1, c2, s, t, as in
   Murphy. */
{ mv_poly_t tmp1, tmp2, tmp3, tmp4;
  int i, j;

  tmp1.terms = (mv_term_t *)malloc(15*sizeof(mv_term_t));
  tmp2.terms = (mv_term_t *)malloc(15*sizeof(mv_term_t));
  tmp3.terms = (mv_term_t *)malloc(15*sizeof(mv_term_t));
  tmp4.terms = (mv_term_t *)malloc(15*sizeof(mv_term_t));
  if (F->terms != NULL) free(F->terms);
  F->terms = (mv_term_t *)malloc(10*sizeof(mv_term_t));

  /* F <-- a_0y^d. */
  F->numTerms = 1;
  mv_term_setone(&F->terms[0]);
  F->terms[0].exps[Y_LOC] = d;
  F->terms[0].exps[A_START + 0] = 1;
 
  /* tmp1 <-- x-yt. */
  tmp1.numTerms = 2;
  mv_term_setone(&tmp1.terms[0]); mv_term_setone(&tmp1.terms[1]);
  tmp1.terms[0].exps[X_LOC]=1;
  tmp1.terms[1].exps[Y_LOC]=1; tmp1.terms[1].exps[T_LOC]=1; tmp1.terms[1].coef = -1;
  mv_poly_cp(&tmp3, &tmp1);

  for (i=1; i<=d; i++) {
    /* We need to do F <-- F + a_i*y^{d-i}*tmp1. */
    mv_poly_cp(&tmp2, &tmp1);
    for (j=0; j<tmp2.numTerms; j++) {
      tmp2.terms[j].exps[Y_LOC] += d-i;
      tmp2.terms[j].exps[A_START + i] += 1;
    }
    mv_poly_add(F, F, &tmp2);

    /* Now, tmp1 <-- tmp1*tmp3. */
    mv_poly_mul(&tmp1, &tmp1, &tmp3);
  }
  /* Finally, the rest of the junk. I believe there is a typo at this
     point in Murphy's thesis, on p. 84. I am understanding the notation
     to mean: x_t = x-t, and m_t = m+t.
     In that case, f(x) should surely be:
       f(x) = f_m(x_t) + (c_2x_t^2 + c_1x_t + c_0)*(x - m_t)
     whereas it says               ...  (x_t - m_t).
     But in order that f(m_t) == 0 ( mod N), we clearly need the first of these. */

  tmp4.numTerms = 15;
  for (i=0; i<15; i++)
    mv_term_setone(&tmp4.terms[i]);
  tmp4.terms[0].exps[C1_LOC] = 1; tmp4.terms[0].exps[X_LOC]=2;

  tmp4.terms[1].coef = -1;
  tmp4.terms[1].exps[C1_LOC] = 1; tmp4.terms[1].exps[X_LOC]=1;
  tmp4.terms[1].exps[M_LOC] = 1;  

  tmp4.terms[2].coef = -2;
  tmp4.terms[2].exps[C1_LOC] = 1; tmp4.terms[2].exps[X_LOC]=1;
  tmp4.terms[2].exps[T_LOC] = 1; 

  tmp4.terms[3].exps[C1_LOC] = 1;  tmp4.terms[3].exps[T_LOC] = 1;  
  tmp4.terms[3].exps[M_LOC] = 1; 

  tmp4.terms[4].exps[C1_LOC] = 1; tmp4.terms[4].exps[T_LOC] = 2;  

  tmp4.terms[5].exps[C0_LOC] = 1; tmp4.terms[5].exps[X_LOC] = 1;  

  tmp4.terms[6].coef = -1;       tmp4.terms[6].exps[C0_LOC] = 1; 
  tmp4.terms[6].exps[M_LOC] = 1; 

  tmp4.terms[7].coef = -1;       tmp4.terms[7].exps[C0_LOC] = 1; 
  tmp4.terms[7].exps[T_LOC] = 1; 

  tmp4.terms[8].coef = -1;       tmp4.terms[8].exps[C2_LOC] = 1; 
  tmp4.terms[8].exps[T_LOC] = 3; 

  tmp4.terms[9].coef = -1;       tmp4.terms[9].exps[C2_LOC] = 1; 
  tmp4.terms[9].exps[T_LOC] = 2; tmp4.terms[9].exps[M_LOC] = 1;

  tmp4.terms[10].coef = 3;       tmp4.terms[10].exps[C2_LOC] = 1; 
  tmp4.terms[10].exps[T_LOC] = 2; tmp4.terms[10].exps[X_LOC] = 1;

  tmp4.terms[11].coef = 2;        tmp4.terms[11].exps[C2_LOC] = 1; 
  tmp4.terms[11].exps[T_LOC] = 1; tmp4.terms[11].exps[X_LOC] = 1;
  tmp4.terms[11].exps[M_LOC] = 1; 

  tmp4.terms[12].coef = -3;        tmp4.terms[12].exps[C2_LOC] = 1; 
  tmp4.terms[12].exps[T_LOC] = 1; tmp4.terms[12].exps[X_LOC] = 2;

  tmp4.terms[13].coef = -1;        tmp4.terms[13].exps[C2_LOC] = 1; 
  tmp4.terms[13].exps[X_LOC] = 2; tmp4.terms[13].exps[M_LOC] = 1;

  tmp4.terms[14].coef = 1;        tmp4.terms[14].exps[C2_LOC] = 1; 
  tmp4.terms[14].exps[X_LOC]=3;


  /* And homogenize: */
  for (i=0; i<15; i++)
    tmp4.terms[i].exps[Y_LOC] = d - tmp4.terms[i].exps[X_LOC];

  mv_poly_add(&tmp2, F, &tmp4);
  mv_poly_simplify(F, &tmp2);

  free(tmp1.terms); free(tmp2.terms); free(tmp3.terms); free(tmp4.terms);
}
/********************************************************************************/
void integrateDx(mv_poly_t *F)
/* Integrate the given polynomial wrt x, then evaluate from -s^{0.5} to s^{0.5}. */
{ int t, tmpE;
  mv_poly_t tmp1, tmp2;
  

  tmp1.terms = (mv_term_t *)malloc(F->numTerms*sizeof(mv_term_t));
  tmp2.terms = (mv_term_t *)malloc(F->numTerms*sizeof(mv_term_t));
  tmp1.numTerms = tmp2.numTerms = F->numTerms;

  for (t=0; t<F->numTerms; t++) {
    mv_term_cp(&tmp1.terms[t], &F->terms[t]);
    tmp1.terms[t].exps[X_LOC] += 1;
    tmp1.terms[t].coef /= tmp1.terms[t].exps[X_LOC];
    mv_term_cp(&tmp2.terms[t], &tmp1.terms[t]);
  }

  /* Now, evaluate tmp1 at x=s^{0.5}. */
  for (t=0; t<tmp1.numTerms; t++) {
    tmp1.terms[t].exps[S_LOC] += tmp1.terms[t].exps[X_LOC]*0.5;
    tmp1.terms[t].exps[X_LOC] = 0;
  }

  /* Evaluate tmp2 at x = -s^{0.5} and negate leading coefficients. */
  for (t=0; t<tmp1.numTerms; t++) {
    tmp2.terms[t].exps[S_LOC] += tmp2.terms[t].exps[X_LOC]*0.5;
    tmpE = (int) tmp2.terms[t].exps[X_LOC];
    if (tmpE %2 == 0)
      tmp2.terms[t].coef *= -1;
    tmp2.terms[t].exps[X_LOC] = 0;
  }

  mv_poly_add(F, &tmp1, &tmp2);
  free(tmp1.terms); free(tmp2.terms);
}

/********************************************************************************/
void integrateDy(mv_poly_t *F)
/* Integrate the given polynomial wrt y, then evaluate from -s^{-0.5} to s^{-0.5}. */
{ int t, tmpE;
  mv_poly_t tmp1, tmp2;
  

  tmp1.terms = (mv_term_t *)malloc(F->numTerms*sizeof(mv_term_t));
  tmp2.terms = (mv_term_t *)malloc(F->numTerms*sizeof(mv_term_t));
  tmp1.numTerms = tmp2.numTerms = F->numTerms;

  for (t=0; t<F->numTerms; t++) {
    mv_term_cp(&tmp1.terms[t], &F->terms[t]);
    tmp1.terms[t].exps[Y_LOC] += 1;
    tmp1.terms[t].coef /= tmp1.terms[t].exps[Y_LOC];
    mv_term_cp(&tmp2.terms[t], &tmp1.terms[t]);
  }

  /* Now, evaluate tmp1 at y=s^{-0.5}. */
  for (t=0; t<tmp1.numTerms; t++) {
    tmp1.terms[t].exps[S_LOC] += tmp1.terms[t].exps[Y_LOC]*(-0.5);
    tmp1.terms[t].exps[Y_LOC] = 0;
  }

  /* Evaluate tmp2 at y = -s^{-0.5} and negate leading coefficients. */
  for (t=0; t<tmp1.numTerms; t++) {
    tmp2.terms[t].exps[S_LOC] += tmp2.terms[t].exps[Y_LOC]*(-0.5);
    tmpE = (int) tmp2.terms[t].exps[Y_LOC];
    if (tmpE %2 == 0)
      tmp2.terms[t].coef *= -1;
    tmp2.terms[t].exps[Y_LOC] = 0;
  }

  mv_poly_add(F, &tmp1, &tmp2);
  free(tmp1.terms); free(tmp2.terms);
}




/********************************************************************************/
void estimate_contp_sampling(double *contp, s32 *p, int numP, mpz_poly f,
                                s32 A0, s32 A1, s32 B1, s32 sampleSize)
/* As in Murphy, estimate cont_p(v) by sampling in the range A0<=a<=A1, 0<b<=B1. */
{ s32         a, b, aSize, bSize;
  s32         i, j, g;
  static mpz_t eval, tmp;
  static int   initialized=0;

  if (!initialized) {
    mpz_init(eval); mpz_init(tmp);
    initialized=1;
  }

  aSize = A1-A0+1; bSize = B1+1;
  for (j=0; j<numP; j++) contp[j]=0.0;
  i=0;
  do {
    a = A0 + prand()%aSize; b = 1 + prand()%bSize;
    g = gcd(a,b);
    if ((g==1)||(g==-1)) {
      mpz_evalF(eval, a, b, f);
      mpz_abs(eval, eval);
      for (j=0; j<numP; j++) {
        mpz_set_ui(tmp, p[j]);
        contp[j] += mpz_remove(eval, eval, tmp);
      }
      i++;
    }
  } while (i<sampleSize);
  for (j=0; j<numP; j++)
    contp[j] /= (double)sampleSize;
}

/*************************************************************/
double estimate_contp_approx(mpz_poly f, s32 p)
{ s32   zeros[MAXPOLYDEGREE];
  int    numZeros;
  poly_t f_;

  mpz_poly_modp(f_, f, p);
  numZeros = poly_getZeros(zeros, f_, p);
  return  (double)numZeros/((double)p-1.0);
}

/**************************************************************/
double estimateAlpha(mpz_poly f, int sampleBound, int B)
/**************************************************************/
/* Estimate alpha, by sample for primes p<=sampleBound, and   */
/* by approximation for sampleBound < p <= B.                 */
/**************************************************************/
{ double         contp[MAX_CONTPRIMES], alpha;
  int            i;
  static s32   *p=NULL, numP=0, numP1=0, lastB=0;
  static int     initialized=0;

  if (B != lastB) {
    free(p);
    initialized=0;
  }
  if (!initialized) {
    numP = getMaxP(1, B);
    p = getPList(&numP);
    lastB = B;
    initialized=1;
  }
  for (i=0; i<numP; i++) {
    if (p[i] < sampleBound) numP1 = i;
  }
  estimate_contp_sampling(contp, p, numP1, f, MIN_A, MAX_A, MAX_B, contSampleSize);

  for (i=numP1; i<numP; i++)
    contp[i] = estimate_contp_approx(f, p[i]);
  alpha=0.0;
  for (i=0; i<numP; i++)
    alpha += ( 1/((double)p[i]-1.0) - contp[i])*log((double)p[i]);

  return alpha;
}

/***************************************************************/
double est_rating_non_skewed(mpz_poly f)
/* Compute the estimated rating $\mathbb{E}(F_1)$, as in Murphy,
   Eq. (5.6).
*/
{ int    i, K=HALF_CIRCLE_PIECES;
  s32   B=CONT_PRIME_BOUND;
  double theta, alpha, logB, est, u;

  alpha = estimateAlpha(f, SAMPLE_PRIME_BOUND, B);

  theta = M_PI/(2*K);
  logB = log((double)B);
  est = 0.0;
  for (i=0; i<K; i++) {
    u = (log(fabs(mpz_evalF_d(cos(theta), sin(theta), f))) + alpha)/logB;
    est += dickman(u);
    theta += M_PI/((double)K);
  }
  return est;
}

/***************************************************************/
double est_rating_skewed(mpz_poly f, mpz_t m, double s)
/* Compute the estimated rating $\mathbb{E}(F_1, F_2)$, as in Murphy,
   Eq. (5.6).
*/
{ int    i, K=HALF_CIRCLE_PIECES;
  s32   B=CONT_PRIME_BOUND;
  double theta, alpha1, alpha2, logB, u1, u2;
  long   double est;
  double s1=sqrt(s), s2=1.0/sqrt(s), dM = mpz_get_d(m);
  static mpz_poly f2;
  static int initialized=0;

  if (!initialized) {
    mpz_poly_init(f2);
    initialized=1;
  }
  f2->degree=1;
  mpz_set_ui(&f2->coef[1], 1);
  mpz_neg(&f2->coef[0], m);

  alpha1 = estimateAlpha(f, SAMPLE_PRIME_BOUND, B);
  alpha2 = estimateAlpha(f2, SAMPLE_PRIME_BOUND, B);

  theta = M_PI/(2*K);
  logB = log((double)B);
  est = 0.0;
  for (i=0; i<K; i++) {
    u1 = (log(fabs(mpz_evalF_d(s1*cos(theta), s2*sin(theta), f))) + alpha1)/logB;
    /* F_2(x,y) = x - y*m, so:                        */
    /*   F_2(s1*cos(theta), s2*sin(theta)) = s1*cos(theta) - m*s2*sin(theta). */
    u2 = (log(fabs(s1*cos(theta) - dM*s2*sin(theta))) + alpha2)/logB;
    est += (long double)dickman(u1)*(long double)dickman(u2);
    /* est += (long double)dickman_old(u1,55)*(long double)dickman_old(u2,55); */
    theta += M_PI/((double)K);
  }
  return (double)est;
}
  
  
/***************************************************************/
void setI(int degree)
{ mv_poly_t F;
  int i, j;
  static int initialized=0;

  if (!initialized) {
    I.terms = NULL;
    initialized=1;
  } else {
    free(I.terms);
    I.terms = NULL;
  }
  d = degree;
  F.terms = NULL;
  computeF(&F, d);

#ifdef _VERY_LOUD
  printf("Degree %d, F = \n", degree);
  mv_poly_print(stdout, &F);
#endif

  mv_poly_mul(&I, &F, &F);

  integrateDx(&I);
  integrateDy(&I);

#ifdef _VERY_LOUD
  printf("Degree %d, I = \n", degree);
  mv_poly_print(stdout, &I);
#endif

  /* These things help speed up evaluation: */
  maxPowC0=maxPowC1=maxPowT=maxPowS=maxPowM=0;
  for (i=0; i<=d; i++)
    maxPowA[i]=0;

  for (i=0; i<I.numTerms; i++) {
    maxPowC0 = MAX(maxPowC0, I.terms[i].exps[C0_LOC]);
    maxPowC1 = MAX(maxPowC1, I.terms[i].exps[C1_LOC]);
    maxPowC2 = MAX(maxPowC2, I.terms[i].exps[C2_LOC]);
    maxPowT =  MAX(maxPowT,  I.terms[i].exps[T_LOC]);
    maxPowS =  MAX(maxPowS,  I.terms[i].exps[S_LOC]);
    maxPowM =  MAX(maxPowM,  I.terms[i].exps[M_LOC]);
    for (j=0; j<=d; j++) 
      maxPowA[j] = MAX(maxPowA[j], I.terms[i].exps[A_START+j]);
  }
}

/***************************************************************/
void mv_poly_eval(mpf_t res, mv_poly_t *F, __mpf_struct *X, int numVars)
{ static mpf_t thisTerm;
  static __mpf_struct X_pows[MAX_MV_POLY_SIZE*12];
  static int initialized=0;
  int    e, i, t, eX[MAX_MV_POLY_SIZE], div;

  if (!initialized) {
    mpf_init2(thisTerm, 256);
    for (i=0;i<MAX_MV_POLY_SIZE*12; i++)
      mpf_init2(&X_pows[i], 256);
    initialized=1;
  }

  for (i=0; i<MAX_MV_POLY_SIZE; i++) {
    mpf_set_ui(&X_pows[12*i], 1);
    eX[i]=0;
  }
  mpf_set_ui(res, 0);
  for (t=0; t<F->numTerms; t++) {
    mpf_set_d(thisTerm, F->terms[t].coef);
    for (i=0; i<numVars; i++) {
      e = (int)F->terms[t].exps[i];
      div=0;
      if (e<0) { e=-e;div=1; }
      while (e > eX[i]) {
        mpf_mul(&X_pows[12*i + (eX[i]+1)], &X_pows[12*i + (eX[i])], &X[i]);
        eX[i] += 1;
      }
      if (e) {
        if (div) {
          if (mpf_sgn(&X_pows[12*i+e]) != 0) {
            mpf_div(thisTerm, thisTerm, &X_pows[12*i + e]);
          } else {
          /* There really should be something smarter here! */
            mpf_mul_ui(thisTerm, thisTerm, 1000000);
          }
        } else
          mpf_mul(thisTerm, thisTerm, &X_pows[12*i + e]);
      }
    }
    mpf_add(res, res, thisTerm);
  }
}  

/***************************************************************/
void Iparam_to_point(__mpf_struct *X, Iparam_t *p)
{ int i;

  mpf_set(&X[C0_LOC], p->c0); mpf_set(&X[C1_LOC], p->c1);
  mpf_set(&X[C2_LOC], p->c2);
  mpf_set(&X[T_LOC], p->t); mpf_set(&X[S_LOC], p->s);
  mpf_set_z(&X[M_LOC], p->m);

  mpf_set_ui(&X[X_LOC], 0); mpf_set_ui(&X[Y_LOC], 0);

  for (i=0; i<=p->f->degree; i++)
    mpf_set_z(&X[A_START+i], &p->f->coef[i]);
}

/***************************************************************/
void point_to_Iparam(__mpf_struct *X, Iparam_t *p)
{ int i;
  static mpf_t tmp1;
  static int initialized=0;

  if (!initialized) {
    mpf_init2(tmp1, 256); initialized=1;
  }
  
  mpf_set(p->c0, &X[C0_LOC]); mpf_set(p->c1, &X[C1_LOC]);
  mpf_set(p->c2, &X[C2_LOC]);
  mpf_set(p->t, &X[T_LOC]); mpf_set(p->s, &X[S_LOC]);
  mpf_set_d(tmp1, 0.5); mpf_add(tmp1, tmp1, &X[M_LOC]); mpz_set_f(p->m, tmp1);

  for (i=0; i<=p->f->degree; i++) {
    mpf_set_d(tmp1, 0.5); mpf_add(tmp1, tmp1, &X[A_START+i]);
    mpz_set_f(&p->f->coef[i], tmp1);
  }
}

/***************************************************************/
void Iparam_cp(Iparam_t *dest, Iparam_t *src)
/***************************************************************/
{
  mpf_set(dest->c0, src->c0);
  mpf_set(dest->c1, src->c1);
  mpf_set(dest->c2, src->c2);
  mpz_set(dest->m, src->m);
  mpz_set(dest->n, src->n);
  mpz_poly_cp(dest->f, src->f);
  dest->logSize = src->logSize;
  dest->score = src->score;
}  

/***************************************************************/
void diff(mv_poly_t *dF, mv_poly_t *F, int var)
/***************************************************************/
/* Compute the partial derivative of F  wrt  x_{var}.          */
/***************************************************************/
{ int    t;
  double e;
  mv_poly_t tmp;

  tmp.terms = NULL;
  mv_poly_cp(&tmp, F);
  for (t=0; t<tmp.numTerms; t++) {
    e = tmp.terms[t].exps[var];
    if (fabs(e) < 0.00001) {
      tmp.terms[t].coef = 0;
    } else {
      tmp.terms[t].coef *= e;
      tmp.terms[t].exps[var] -= 1;
    }
  }
  mv_poly_simplify(dF, &tmp);
  if (tmp.terms != NULL) free(tmp.terms);
}

/***************************************************************/
void getGradient()
/* Compute the gradient of I in (c0, c1, c2, s, t) space.*/
{ static int initialized=0;

  if (!initialized) {
    dIdC0.terms = NULL; dIdC1.terms = NULL; dIdC2.terms=NULL;
    dIdT.terms = NULL;  dIdS.terms = NULL;
    initialized=1;
  }
  diff(&dIdC0, &I, C0_LOC);
  diff(&dIdC1, &I, C1_LOC);
  diff(&dIdC2, &I, C2_LOC);
  diff(&dIdT, &I, T_LOC);
  diff(&dIdS, &I, S_LOC);
}

/******************************************************/
void evalGradient(__mpf_struct *grad, __mpf_struct *X)
/* The x&y entries of X should exist, but be zero.    */
{
  mv_poly_eval(&grad[0], &dIdC0, X, 8+d+1);
  mv_poly_eval(&grad[1], &dIdC1, X, 8+d+1);
  mv_poly_eval(&grad[2], &dIdC2, X, 8+d+1);
  mv_poly_eval(&grad[3], &dIdT, X, 8+d+1);
  mv_poly_eval(&grad[4], &dIdS, X, 8+d+1);
}

#define STEP_RATIO 2.0
/******************************************************/
double steepestDescentMinimize(__mpf_struct *Res, __mpf_struct *X, int level)
/* The higher the level, the harder we'll try to minimize.
   level=1 is the lowest. It will give something ressembling a min,
           but not too good. But it will go quickly.
   level=10 is a reasonable upper bound. Higher than that may take
           very long and not get you too much better. It will often
           take some time at that level.
*/
{ double R=32.0, RStop, r, Int, ir;
  int    i, consecGood, cont, it=0, smallImprove=0;
  static __mpf_struct grad[5], P[MAX_MV_POLY_SIZE], P1[MAX_MV_POLY_SIZE];
  static mpf_t tmp1, tmp2, evalP, evalP1;
  static int initialized=0;

//  RStop = R/(20*(MAX(level, 1)));
  RStop = 1.0/(10*(MAX(level, 1)));
  if (!initialized) {
    for (i=0; i<5; i++)
      mpf_init2(&grad[i], 256);
    for (i=0; i<MAX_MV_POLY_SIZE; i++) {
      mpf_init2(&P[i], 256);
      mpf_init2(&P1[i], 256);
    }
    mpf_init2(tmp1, 256); mpf_init2(tmp2, 256); mpf_init2(evalP, 256); mpf_init2(evalP1, 256);
    initialized=1;
  }

  for (i=0; i<MAX_MV_POLY_SIZE; i++) {
    mpf_set(&P[i], &X[i]);
    mpf_set(&P1[i], &X[i]); /* Remember : after the first 5 coordinates, the rest never change! */
  }

  consecGood=0;
  do {
    evalGradient(grad, P);
    /* Normalize the gradient. */
    mpf_set_ui(tmp1, 0);
    for (i=0; i<5; i++) {
      mpf_mul(tmp2, &grad[i], &grad[i]);
      mpf_add(tmp1, tmp1, tmp2);
    }
    mpf_sqrt(tmp1, tmp1);
    for (i=0; i<5; i++)
      mpf_div(&grad[i], &grad[i], tmp1);

    mv_poly_eval(evalP, &I, P, 8+d+1);
    do {
      cont=0;
      for (i=0; i<5; i++) {
        mpf_set_d(&P1[i], -R);
        mpf_mul(&P1[i], &P1[i], &grad[i]);
        mpf_add(&P1[i], &P[i], &P1[i]);
      }
      mpf_abs(&P1[S_LOC], &P1[S_LOC]); /* With large leaps, 's' could have become negative. */
      mv_poly_eval(evalP1, &I, P1, 8+d+1);

      if (mpf_cmp(evalP1, evalP)>0) {
        /* The function increased, so we surely jumped too far. */
//printf("At R=%1.5lf, function increased.\n", R);
        R = R/STEP_RATIO;
        consecGood=0;
        cont=1;
      } else {
//printf("At R=%1.5lf, function decreased.\n", R);
        consecGood++;
        r = log(sqrt(mpf_get_d(evalP))) - log(sqrt(mpf_get_d(evalP1)));
        ir = fabs(r/(0.5*log(mpf_get_d(evalP))));
        if (ir < 0.0002) {
          smallImprove++;
          if (ir < 0.000001)
            smallImprove += 10; /* Much larger penalty for this small of an improvement. */
        } else smallImprove=0;
      }
      if (consecGood >= 4) {
        R *= STEP_RATIO;
        consecGood = 0;
      }
      if (fabs(R) < 0.01) cont=0;  /* 0.1 */
      if (smallImprove > 20*level) cont=0; /* 10 */
    } while (cont);
    
    /* Check tolerence, and do P <-- P1, if ok. */
    if (mpf_cmp(evalP1, evalP)<0) {
      for (i=0; i<MAX_MV_POLY_SIZE; i++) 
        mpf_set(&P[i], &P1[i]);
    } else R=0;
    if (++it%500==0) { 
      Int = log(fabs(mpf_get_d(evalP)))/(2.0*M_LN2);
    }
  } while ((R > RStop)&&(smallImprove < 100*level));
  Int = log(fabs(mpf_get_d(evalP)))/(2.0*M_LN2);
  printTmp("k: %ld, log_2(I) = %1.2lf", iteration, Int);
  for (i=0; i<MAX_MV_POLY_SIZE; i++) 
    mpf_set(&Res[i], &P[i]);
//  mv_poly_eval(evalP, &I, P, 7+d+1);
  return Int;

}

/******************************************************************/  
void fixParams(Iparam_t *Opt, Iparam_t *Orig)
/* Round the 5 optimized fields and compute the polynomial. The original
   polynomial and it's 'm' value are in Orig.
*/
{ static mpf_t tmp1, tmp2, tmp3, tmp4;
  static mpz_poly tpol1, tpol2, tpol3, tpol4;
  static int initialized=0;
  int    i, j;

  if (!initialized) {
    mpf_init2(tmp1, 256); mpf_init2(tmp2, 256);
    mpf_init2(tmp3, 256); mpf_init2(tmp4, 256);
    mpz_poly_init(tpol1); mpz_poly_init(tpol2);
    mpz_poly_init(tpol3); mpz_poly_init(tpol3);
    initialized=1;
  }
  mpf_set_d(tmp1, 0.5);
  mpf_add(Opt->c0, Opt->c0, tmp1); mpf_floor(Opt->c0, Opt->c0);
  mpf_add(Opt->c1, Opt->c1, tmp1); mpf_floor(Opt->c1, Opt->c1);
  mpf_add(Opt->c2, Opt->c2, tmp1); mpf_floor(Opt->c2, Opt->c2);
  mpf_add(Opt->t, Opt->t, tmp1); mpf_floor(Opt->t, Opt->t);


  tpol1->degree = 2;

  mpz_set_f(&tpol1->coef[1], Opt->c1);
  mpf_mul(tmp1, Opt->c1, Opt->t);
  mpf_sub(tmp1, Opt->c0, tmp1);
  mpz_set_f(&tpol1->coef[0], tmp1);
  /* And the c2 stuff: */
  mpz_set_f(&tpol1->coef[2], Opt->c2);
  mpf_set_si(tmp1, -2); mpf_mul(tmp1, tmp1, Opt->c2);
  mpf_mul(tmp1, tmp1, Opt->t); mpf_set_z(tmp2, &tpol1->coef[1]);
  mpf_add(tmp2, tmp2, tmp1);  mpz_set_f(&tpol1->coef[1], tmp2);

  mpf_mul(tmp1, Opt->t, Opt->t); mpf_mul(tmp1, tmp1, Opt->c2);
  mpf_set_z(tmp2, &tpol1->coef[0]); mpf_add(tmp2, tmp2, tmp1);
  mpz_set_f(&tpol1->coef[0], tmp2);

  


  tpol2->degree = 1;
  mpz_set_ui(&tpol2->coef[1], 1);
  mpz_set_f(&tpol2->coef[0], Opt->t);
  mpz_add(Opt->m, &tpol2->coef[0], Orig->m); 
  mpz_neg(&tpol2->coef[0], Opt->m);
  mpz_poly_mul(tpol1, tpol1, tpol2);

  /* Now, compute tpol1 += f(x-t). 
     First, set tpol2 <-- x-t. We will keep tpol3 = (x-t)^i. */

  tpol2->degree = 1;
  mpz_set_ui(&tpol2->coef[1], 1);
  mpz_set_f(&tpol2->coef[0], Opt->t); mpz_neg(&tpol2->coef[0], &tpol2->coef[0]);
  tpol3->degree = 0;
  mpz_set_ui(&tpol3->coef[0], 1);
  for (i=0; i<=Orig->f->degree; i++) {
    tpol4->degree = i;
    for (j=0; j<=i; j++)
      mpz_mul(&tpol4->coef[j], &tpol3->coef[j], &Orig->f->coef[i]);
    mpz_poly_add(tpol1, tpol1, tpol4);
    /* Now tpol3 *= tpol2. */
    mpz_poly_mul(tpol3, tpol3, tpol2);
  }
  mpz_poly_cp(Opt->f, tpol1);

/* CJM 5/28/04 */
  mpf_set_ui(Opt->c0, 0); mpf_set_ui(Opt->c1, 0); mpf_set_ui(Opt->c2, 0);
  mpf_set_ui(Opt->t, 0); mpf_set(Opt->s, Orig->s);

}

/******************************************************************/
double minimize_wrt_s(Iparam_t *Param)
/******************************************************************/
/* This function is used for recomputing the skew when the other  */
/* parameters are fixed (for example, after they've been rounded  */
/* to integers).                                                  */
/* It uses only Param->f and Param->m, so these should already    */
/* have been adjusted to account for any nonzero c0 and c1.       */
/******************************************************************/
/* CJM, 5/28/04                                                   */
/* This function has also been verified to work to at least a few */
/* significant digits.                                            */
/******************************************************************/
{ static __mpf_struct X[MAX_MV_POLY_SIZE];
  static mpf_t        tmp1, R, s0, s1, evalS0, evalS1;
  static int          initialized=0;
  int                 i, cont, consecGood=0;
  

  if (!initialized) {
    mpf_init2(tmp1, 256); mpf_init2(s0, 256); mpf_init2(s1, 256);
    mpf_init2(R, 256); mpf_init2(evalS0, 256); mpf_init2(evalS1, 256);
    for (i=0; i<MAX_MV_POLY_SIZE; i++) 
      mpf_init2(&X[i], 256);
    /* These spots will always be zero. */
    mpf_set_ui(&X[C2_LOC], 0);
    mpf_set_ui(&X[C1_LOC], 0);
    mpf_set_ui(&X[C0_LOC], 0);
    mpf_set_ui(&X[T_LOC], 0);
    initialized=1;
  }

  /* Everything but 's' is fixed, so set them now: */
  mpf_set_z(&X[M_LOC], Param->m);
  for (i=0; i<=d; i++) 
    mpf_set_z(&X[A_START+i], &Param->f->coef[i]);
  /* This is the initial estimate for the optimal 's': */
  mpf_set(&X[S_LOC], Param->s);
  mpf_set(s0, Param->s);

  /* Ad-hoc method to find the local min: */
  cont=1;
  if (mpf_cmp_ui(s0, 1)<=0) {
    mpf_set_d(R, 0.25);
    mpf_set_ui(s0, 100); /* To be sure: could be a bad local min? */
  }
  else mpf_set_ui(R, 1.0);
  mv_poly_eval(evalS0, &I, X, 8+1+d);
  do {
    /* Compute s1 <-- s0 +/- R, depending on the sign of dI/dS.  */
    mv_poly_eval(tmp1, &dIdS, X, 8+1+d);
    if (mpf_sgn(tmp1)<0)
      mpf_add(s1, s0, R);
    else mpf_sub(s1, s0, R); 

    if (mpf_sgn(s1) < 0) {
      while (mpf_cmp(R, s0) > 0)
        mpf_div_ui(R, R, 2);
      mpf_div_ui(s1, s0, 2); 
    }

    if (mpf_cmp_ui(s1, 1) <= 0)
      mpf_set_d(s1, 1.01);

//printf("s1 = "); mpf_out_str(stdout, 10, 10, s1); printf("\n");
    if (mpf_get_d(s1) > maxSkew + 50) {
      cont=0;
    } else if (mpf_sgn(s1) < 0) cont=0;
    else {
      /* Now check: */
      mpf_set(&X[S_LOC], s1); 
      mv_poly_eval(evalS1, &I, X, 8+1+d);
      if (mpf_cmp(evalS1, evalS0) < 0) {
        consecGood++;
        mpf_set(s0, s1);
        mpf_set(evalS0, evalS1);
      } else {
        consecGood=0;
        mpf_div_ui(R, R, 2);
      }
      if (consecGood >= 3) {
        consecGood=0;
        mpf_mul_ui(R, R, 2);
      }
      if (mpf_get_d(R) < 0.001)
        cont=0;
    }
  } while (cont); 

  mpf_set(Param->s, s0);
  return mpf_get_d(s0);
}

/******************************************************************/
int minimize_wrt_t(Iparam_t *Param)
/******************************************************************/
/* This function is used for recomputing the t when the other     */
/* parameters are fixed (for example, after we've done the        */
/* sieving for root properties which will have changed the size   */
/* properties a little, use the function to try to get back some  */
/* of the size without altering the roots).                       */
/*   On start, the polynomial should already have been adjusted   */
/* according to any c0, c1 and earlier t. On completion from here,*/
/* these fields will all be zero and the polynomial will have     */
/* already been shifted using the 't' found.                      */
/******************************************************************/
{ static __mpf_struct X[MAX_MV_POLY_SIZE];
  static mpf_t        tmp1, R, t0, t1, evalT0, evalT1;
  static int          initialized=0;
  int                 i, cont, consecGood=0;
  

  if (!initialized) {
    mpf_init2(tmp1, 256); mpf_init2(t0, 256); mpf_init2(t1, 256);
    mpf_init2(R, 256); mpf_init2(evalT0, 256); mpf_init2(evalT1, 256);
    for (i=0; i<MAX_MV_POLY_SIZE; i++) 
      mpf_init2(&X[i], 256);
    initialized=1;
  }
  mpf_set_ui(&X[C2_LOC], 0);
  mpf_set_ui(&X[C1_LOC], 0);
  mpf_set_ui(&X[C0_LOC], 0);

  /* All but 't' are fixed, so set them now: */
  mpf_set(&X[S_LOC], Param->s);
  mpf_set_z(&X[M_LOC], Param->m);
  for (i=0; i<=d; i++) 
    mpf_set_z(&X[A_START+i], &Param->f->coef[i]);
  /* This is the initial estimate for the optimal 't': */
  mpf_set_ui(&X[T_LOC], 0);
  mpf_set_ui(t0, 0);

  /* Ad-hoc method to find the local min: */
  cont=1;
  mpf_set_ui(R, 1);
  mv_poly_eval(evalT0, &I, X, 8+1+d);
  do {
    /* Compute t1 <-- t0 +/- R, depending on the sign of dI/dT.  */
    mv_poly_eval(tmp1, &dIdT, X, 8+1+d);
    if (mpf_sgn(tmp1)<0)
      mpf_add(t1, t0, R);
    else mpf_sub(t1, t0, R); 

    /* Now check: */
    mpf_set(&X[T_LOC], t1); mv_poly_eval(evalT1, &I, X, 8+1+d);
    if (mpf_cmp(evalT1, evalT0) < 0) {
      consecGood++;
      mpf_set(t0, t1);
      mpf_set(evalT0, evalT1);
    } else {
      consecGood=0;
      mpf_div_ui(R, R, 2);
    }
    if (consecGood >= 4) {
      consecGood=0;
      mpf_mul_ui(R, R, 2);
    }
    if (mpf_get_d(R) < 0.1)
      cont=0;
  } while (cont); 

#if 0
  mpf_set(Param->c0, &X[C0_LOC]); 
  mpf_set(Param->c1, &X[C1_LOC]); 
  mpf_set(Param->c2, &X[C2_LOC]);
#else
  mpf_set_ui(Param->c0, 0); mpf_set_ui(Param->c1, 0); mpf_set_ui(Param->c2,0);
#endif

  mpf_set(Param->t, t0);
  fixParams(Param, Param);
  
  return 0;
}

/******************************************************************/
int optimizeParameters(Iparam_t *Opt, Iparam_t *Orig)
/******************************************************************/
/* The only fields used from Orig are `m' and `f'. But all fields */
/* of `Opt' will be filled in.                                    */
/******************************************************************/
{ static __mpf_struct X[MAX_MV_POLY_SIZE];
  static mpf_t        tmp1;
  static mpz_poly     tpol1, tpol2, tpol3, tpol4;
  static Iparam_t     ThisOpt;
  static int          initialized=0;
  int                 i;
  double              thisMin;
  

  if (!initialized) {
    mpf_init2(tmp1, 256);
    mpz_poly_init(tpol1); mpz_poly_init(tpol2); 
    mpz_poly_init(tpol3); mpz_poly_init(tpol4);
    for (i=0; i<MAX_MV_POLY_SIZE; i++) 
      mpf_init2(&X[i], 256);
    Iparam_init(&ThisOpt);
    initialized=1;
  }
  mpf_set_z(&X[M_LOC], Orig->m);
  for (i=0; i<=d; i++) 
    mpf_set_z(&X[A_START+i], &Orig->f->coef[i]);
  mpz_poly_cp(ThisOpt.f, Orig->f);
    

  mpf_set_ui(&X[C0_LOC], 0); mpf_set_ui(&X[C1_LOC], 0); 
  mpf_set_ui(&X[C2_LOC],0);
  mpf_set_ui(&X[T_LOC], 0);   mpf_set_ui(&X[S_LOC], 500);
 
  thisMin = steepestDescentMinimize(X, X, 2);
//printf("c0 = "); mpf_out_str(stdout, 10, 0, &X[C0_LOC]); printf("\n");
//printf("c1 = "); mpf_out_str(stdout, 10, 0, &X[C1_LOC]); printf("\n");
//printf("c2 = "); mpf_out_str(stdout, 10, 0, &X[C2_LOC]); printf("\n\n");
  if ((thisMin < minStage1Size)&& (mpf_get_d(&X[S_LOC]) < maxSkew)) {
    thisMin = steepestDescentMinimize(X, X, 30);
  } else {
    Opt->logSize = thisMin;
    return thisMin;
  }
  mpf_set(Opt->c0, &X[C0_LOC]);  mpf_set(Opt->c1, &X[C1_LOC]);
  mpf_set(Opt->c2, &X[C2_LOC]);
  mpf_set(Opt->t, &X[T_LOC]);    mpf_set(Opt->s, &X[S_LOC]);    

  /* Recompute the integral. */
  mpf_set(Opt->c0, &X[C0_LOC]);  mpf_set(Opt->c1, &X[C1_LOC]);
  mpf_set(Opt->c2, &X[C2_LOC]);
  mpf_set(Opt->t, &X[T_LOC]);  mpf_set(Opt->s, &X[S_LOC]);
  fixParams(Opt, Orig);
  /* fixParams will have rounded the arguments and recomputed the polynomial. */
  /* So now we should recompute the skew & do the evaluation again.  */
  minimize_wrt_s(Opt);

//printf("Orig->f: "); mpz_poly_print(stdout, "", Orig->f);
//printf(" Opt->f: "); mpz_poly_print(stdout, "", Opt->f);

  mpf_set(&X[C0_LOC], Opt->c0);  mpf_set(&X[C1_LOC], Opt->c1);
  mpf_set(&X[C2_LOC], Opt->c2);
  mpf_set(&X[T_LOC], Opt->t);  mpf_set(&X[S_LOC], Opt->s);
  for (i=0; i<=d; i++) 
    mpf_set_z(&X[A_START+i], &Opt->f->coef[i]);
  mv_poly_eval(tmp1, &I, X, 8+d+1);
  if (mpz_sgn(tmp1) < 0) {
    printf("Some error occurred: evaluation of I is negative!\n");
    for (i=0; i<8+d+1; i++) {
      printf("X[%d] = ", i); 
      mpf_out_str(stdout, 10, 0, &X[i]);
      printf("\n");
    }
    return 0;
  }
  mpf_sqrt(tmp1, tmp1);
  Opt->logSize = log(mpf_get_d(tmp1))/M_LN2;
  /* CJM, 5/28/04:
     The evaluation of the double integral seems to be right -
     I verified a typical case with Maple and it matched to about
     8 significant digits or so (rounding error - I only used 5
     correct digits of s, or it would have probably matched better).
   */
  return 0;
}

/******************************************************************/
int optimizeParameters2(Iparam_t *Opt, Iparam_t *Orig)
/******************************************************************/
/* Similar to above, but we will try to perturb some values to see*/
/* if we can find a better local min.                             */
/* 'Opt' should be the output of optimizeParameters(Opt, Orig).   */
/* So, this function is not intended to replace the above one -   */
/* just to suplement it.                                          */
/******************************************************************/
{ static __mpf_struct X[MAX_MV_POLY_SIZE];
  static mpf_t        tmp1;
  static mpz_poly     tpol1, tpol2, tpol3, tpol4;
  static Iparam_t     ThisOpt;
  static int          initialized=0;
  int                 i, cont=1, try=0;
  s32                dc0=100, dc1=4, dt=20, ds=300;
  double              thisMin;
  

  if (!initialized) {
    mpf_init2(tmp1, 256);
    mpz_poly_init(tpol1); mpz_poly_init(tpol2); 
    mpz_poly_init(tpol3); mpz_poly_init(tpol4);
    for (i=0; i<MAX_MV_POLY_SIZE; i++) 
      mpf_init2(&X[i], 256);
    Iparam_init(&ThisOpt);
    initialized=1;
  }
  mpf_set_z(&X[M_LOC], Orig->m);
  for (i=0; i<=d; i++) 
    mpf_set_z(&X[A_START+i], &Orig->f->coef[i]);
  Iparam_cp(&ThisOpt, Orig);
//  Iparam_cp(Opt, Orig);
  do {
    Iparam_cp(&ThisOpt, Orig);
    /* Right here, try perturbing these a little: */
    if (try==0) {
      mpf_set_ui(&X[C0_LOC], 0); 
      mpf_set_ui(&X[C1_LOC], 0); 
      mpf_set_ui(&X[C2_LOC], 0);
      mpf_set_ui(&X[T_LOC], 0);   
      mpf_set_ui(&X[S_LOC], 0);
    } else {
      mpf_set_ui(&X[C0_LOC], prand()%dc0 - dc0/2); 
      mpf_set_ui(&X[C1_LOC], prand()%dc1 - dc1/2); 
      mpf_set_ui(&X[C2_LOC], 0);
      mpf_set_ui(&X[T_LOC], prand()%dt - dt/2);   
      mpf_set_ui(&X[S_LOC], prand()%ds);
    }

    thisMin = steepestDescentMinimize(X, X, 2);
    if ((thisMin < minStage1Size)&& (mpf_get_d(&X[S_LOC]) < maxSkew)) {
      thisMin = steepestDescentMinimize(X, X, 30);
      mpf_set(ThisOpt.c0, &X[C0_LOC]);  mpf_set(ThisOpt.c1, &X[C1_LOC]);
      mpf_set(ThisOpt.c2, &X[C2_LOC]);
      mpf_set(ThisOpt.t, &X[T_LOC]);    mpf_set(ThisOpt.s, &X[S_LOC]);    

      /* Recompute the integral. */
      fixParams(&ThisOpt, Orig);
      /* fixParams will have rounded the arguments and recomputed the polynomial. */
      /* So now we should recompute the skew & do the evaluation again.  */
      minimize_wrt_s(&ThisOpt);


      mpf_set(&X[C0_LOC], ThisOpt.c0);  mpf_set(&X[C1_LOC], ThisOpt.c1);
      mpf_set(&X[C2_LOC], ThisOpt.c2);
      mpf_set(&X[T_LOC], ThisOpt.t);  mpf_set(&X[S_LOC], ThisOpt.s);
      for (i=0; i<=d; i++) 
        mpf_set_z(&X[A_START+i], &ThisOpt.f->coef[i]);
      mv_poly_eval(tmp1, &I, X, 8+d+1);
      if (mpz_sgn(tmp1) < 0) {
        printf("Some error occurred: evaluation of I is negative!\n");
        for (i=0; i<8+d+1; i++) {
          printf("X[%d] = ", i); 
          mpf_out_str(stdout, 10, 0, &X[i]);
          printf("\n");
        }
      }
      mpf_sqrt(tmp1, tmp1);
      ThisOpt.logSize = log(mpf_get_d(tmp1))/M_LN2;
      if ((ThisOpt.logSize < Opt->logSize) || (try==0)) {
        Iparam_cp(Opt, &ThisOpt);
        Opt->logSize = ThisOpt.logSize;
      }
    }
    if (++try >= 5) cont=0;
  } while (cont);
  /* CJM, 5/28/04:
     The evaluation of the double integral seems to be right -
     I verified a typical case with Maple and it matched to about
     8 significant digits or so (rounding error - I only used 5
     correct digits of s, or it would have probably matched better).
   */
  return 0;
}

/**********************************************************/
int adjustPolForParams(Iparam_t *res, Iparam_t *par, s32 J0, s32 J1)
/**********************************************************/
/* To par->f, do:                                         */
/* f <-- f + (J1*x - J0)(x-m).                            */
/**********************************************************/
{ static int initialized=0;
  static mpz_t mp1, mp2;

  if (!(initialized)) {
    mpz_init(mp1); mpz_init(mp2);
    initialized=1;
  }
  mpz_poly_cp(res->f, par->f);
  mpz_set(res->m, par->m);
  mpz_set(res->n, par->n);
  mpf_set(res->s, par->s);
  mpf_set(res->t, par->t);
  mpf_set(res->c0, par->c0);
  mpf_set(res->c1, par->c1);
  mpf_set(res->c2, par->c2);
  res->logSize = par->logSize;
  res->score = par->score;

  mpz_set_si(mp1, J1);
  mpz_add(&res->f->coef[2], &res->f->coef[2], mp1);
  mpz_set_si(mp1, J1);
  mpz_mul(mp1, mp1, res->m);
  mpz_set_si(mp2, J0);
  mpz_add(mp1, mp1, mp2);
  mpz_sub(&res->f->coef[1], &res->f->coef[1], mp1);
  mpz_set_si(mp1, J0);
  mpz_mul(mp1, mp1, res->m);
  mpz_add(&res->f->coef[0], &res->f->coef[0], mp1);
  return 0;
}


#define MAX_PRIME_POWER 1024
const s32 sievePrimes[]={
  2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
  101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,
  211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
  293,307,311,313,317,331,337,347,
  349,353,359,367,373,379,383,389,
  397,401,409,419,421,431,433,439,
  443,449,457,461,463,467,479,487,
  491,499};
const int numSievePrimes=82; 
#define MAX_KEEP 5
/* This function is screwy : I'm pretty darned sure that there's something wrong. */
/***************************************************************/
int sieveOverSimilarSize(Iparam_t *par, s32 *candJ0, s32 *candJ1, double *logAdj, int numKeep)
/***************************************************************/
/* Sieve over some polynomials with size similar to the given  */
/* one to find the one with best root properties.              */
/***************************************************************/
{ static  mpz_t mp1, mp2;
  static  int initialized=0;
  static  double *Jarray;
  static  s32    bestJ0[MAX_KEEP], bestJ1[MAX_KEEP];
  static  double  bestAlpha[MAX_KEEP];
  static  int     numBest;
  s32    J0=sieveJ0, J1=sieveJ1, j0, j1;
  double  tmpArray[MAX_PRIME_POWER], tmpd1, tmpd2;
  s32    p, l, mRes, fl, tmp1;
  int     i, index, k;

  if (!initialized) {
    if (!(Jarray = (double *)malloc((2*J0+1)*sizeof(double)))) {
      fprintf(stderr, "sieveOverSimilarSize() Fatal memory allocation error!\n");
      exit(-1);
    }
    mpz_init(mp1);
    mpz_init(mp2);
    initialized=1;
  }
  numBest = 0;
  bestAlpha[0]=0.0;

  for (j1=-J1; j1<=J1; j1++) {
    for (j0=0; j0<2*J0+1; j0++)
      Jarray[j0] = 0.0;
    for (i=0; i<numSievePrimes; i++) {
      p = sievePrimes[i]; 
      mRes = mpz_tdiv_ui(par->m, p);

      for (l=0; l<p; l++)
        tmpArray[l] = 0.0;

      for (l=0; l<p; l++) {
        fl = mpz_poly_evalModp(par->f, p, l);
        fl = (fl + mulmod32(j1, l*l, p))%p; /* l^2 << 2^32. */
        fl = (fl - mulmod32(j1, mulmod32(mRes, l, p), p))%p;
        if (fl < 0) fl += p;
        tmp1 = (p+ l - mRes )%p;
        if (tmp1) { /* CJM, 5/28/04: Aha! Take that! */
          tmp1 = inverseModP(tmp1, p);
          j0 = mulmod32(fl, tmp1, p);
        
          if (j0 < 0) j0 += p;
        /* So f_{j1,j0}(l) == 0 (mod p).  */
          tmpArray[j0] += 1.0;
        }
	else{ 
	  /* It's possible that p divides fl as well. */
	  /* In this case p divides the polynomial for all j.*/ 
	  /* Added by EJL 02/14/05. */ 
	  if(fl == 0){
	    for(j0 = 0; j0 < p; j0++)
	      tmpArray[j0] += 1.0;
	  }
	}
      }
      /* Now tmpArray[j0] is the number of roots of f_{j1,j0} mod p.    */
      /* Replace each entry with it's contribution to alpha(F_{j1,j0}). */
      tmpd1 = (double)p/(p+1.0);
      tmpd2 = log((double)p)/(double)(p-1);
      for (j0=0; j0<p; j0++)
        tmpArray[j0] = (1.0 - tmpArray[j0]*tmpd1)*tmpd2;

      /* Now add these value to the J0 array in the appropriate places. */
/* Here: check that index is right! */
      index = (p - (J0%p))%p;
      for (j0=0; j0<2*J0+1; j0++) {
        Jarray[j0] += tmpArray[index];
        index = (index+1)%p;
      }
    }
    /* Make a pass through the J0 array and see if we can beat the best so far. */
/* Below here seems to be okay, believe it or not. It's the actual
   sieving above that I'm skeptical of. It just seems to give (0,0)
   way too often. This was the case, though, even before we were keeping
   the several best candidates around. */
    if (numBest==0) {
      bestAlpha[0] = Jarray[0];
      bestJ0[0] = -J0;
      bestJ1[0] = j1;
      numBest=1;
    }
    for (i=0; i<2*J0+1; i++) {
      if (Jarray[i] < bestAlpha[numBest-1]) {
        /* Find where it should go: */
        l=numBest-1;
        while ((l>0) && (Jarray[i] < bestAlpha[l-1]))
          l--;
        numBest = MIN(numBest+1, MAX_KEEP);
        for (k=numBest-1; k>l; k--) {
          bestAlpha[k] = bestAlpha[k-1];
          bestJ0[k] = bestJ0[k-1];
          bestJ1[k] = bestJ1[k-1];
        }
        bestAlpha[l] = Jarray[i];
        bestJ0[l] = -J0+i;
        bestJ1[l] = j1;
      }
    }
  }
  for (i=0; i<MIN(numKeep, MAX_KEEP); i++) {
    candJ0[i]=bestJ0[i];
    candJ1[i]=bestJ1[i];
    logAdj[i] = bestAlpha[i]/M_LN2;
  }

  return 0;
}



/***************************************************************/
int initialize(int degree)
{
  if (degree > MAX_NFS_POLY_DEGREE) {
    printf("initialize() Error: degree=%d, max degree is %d.\n", degree, MAX_NFS_POLY_DEGREE);
    return -1;
  }
  d=degree; /* CJM, 7/30/04. */
  setI(degree);
  getGradient();
  printf("Initialization done:\n");
#ifdef _VERY_LOUD
  printf("I = "); mv_poly_print(stdout, &I);
  printf("\ndIdC0 = "); mv_poly_print(stdout, &dIdC0);
  printf("\ndIdC1 = "); mv_poly_print(stdout, &dIdC1);
  printf("\ndIdC2 = "); mv_poly_print(stdout, &dIdC2);
  printf("\ndIdT = "); mv_poly_print(stdout, &dIdT);
  printf("\ndIdS = "); mv_poly_print(stdout, &dIdS);
#endif
  return 0;
}

/******************************************************/
void Iparam_init(Iparam_t *par)
{
  mpz_poly_init(par->f);
  mpz_init(par->m); mpz_init(par->n);
  mpf_init2(par->c0, 256); mpf_init2(par->c1, 256);
  mpf_init2(par->c2, 256);
  mpf_init2(par->s, 256); mpf_init2(par->t, 256);
}

/*******************************************************************/
const int smallPrimes[]={
  2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
  73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179};
void choose_lc(mpz_t lc, double approx_lc_log, int degree)
{ int k;
  int e[41];
  s32 remain;
  double l;


  if (mpz_cmp_ui(enumLCD,0) > 0) {
    mpz_add(lastLC, lastLC, enumLCD);
    mpz_set(lc, lastLC);
    mpz_add_ui(enumCt, enumCt, 1);
    if (mpz_cmp(enumCt, e1) > 0)
      exit(0); /* Not so clean, but whatever. */
    return;
  }

  lcp = MAX(3, MIN(lcp, 41));
  leave = MAX(leave, 0);
  mpz_set_ui(lc, lcd);
  l = log((double)lcd)/M_LN2;
  memset(e, 0x00, 41*sizeof(int));
  while (l < approx_lc_log - leave) {
      k = (rand()^prand())%lcp;
    e[k]+=1;
    mpz_mul_ui(lc, lc, smallPrimes[k]);
    l += log((double)smallPrimes[k])/M_LN2;
  }
  remain = (s32)pow(2.0, approx_lc_log - l);
  if (remain > 1)
    mpz_mul_ui(lc, lc, 1 + prand()%remain);

}

/*******************************************************************/
void choose_m(mpz_t m, mpz_t N, double approx_lc_log, int degree)
{ static mpz_t lc, tmp;
  static mpz_poly  tmpPol;
  static mpf_t     tmp1;
  static int initialized=0;

  if (!initialized) {
    mpz_init(lc); mpz_init(tmp);
    mpz_poly_init(tmpPol); mpf_init2(tmp1, 1024);
    initialized=1;
  }
  choose_lc(lc, approx_lc_log, degree);
  mpz_set(tmp, lc);
  mpz_div(tmp, N, tmp);
  mpz_root(m, tmp, degree);
  mpz_sub_ui(m, m, 1);
  
} 


/******************************************************/
double doTest(double minLClog, double maxLClog, Iparam_t *param, int degree)
{ double approx_lc_log, E, oldE, oldLogSize;
  static Iparam_t startParam, thisParam, testParam;
  static mpz_poly oldF;
  static mpz_t    oldM;
  static double   log2N, examineC3Lim;
  static int initialized=0;
  int    cont, try, tryNextM, mTry=0;
  s32   J0[MAX_KEEP], J1[MAX_KEEP], i;
  double logAdj[MAX_KEEP], maxThirdLog2, lcLog2;

  if (!initialized) {
    Iparam_init(&startParam);
    Iparam_init(&thisParam);
    Iparam_init(&testParam);
    mpz_poly_init(oldF); mpz_init(oldM);
    initialize(degree);
    log2N = mpz_sizeinbase(param->n, 2);
    examineC3Lim = 2.5 + (degree - 2.5)*examineFrac;
    initialized=1;
  }
  tryNextM = 0;
  
  mpz_poly_cp(startParam.f, param->f);
  mpz_set(startParam.n, param->n);
  mpz_poly_cp(thisParam.f, param->f);
  mpz_set(thisParam.n, param->n);
  param->score = 0.0;
  approx_lc_log = minLClog + prand_01()*(maxLClog - minLClog);

  do {
      choose_m(startParam.m, startParam.n, approx_lc_log, degree);
      mTry=0;
    /* Is the third coefficient too large to even consider? */
    mpz_poly_getBaseM(startParam.f, startParam.n, startParam.m, 0);
    lcLog2 = mpz_sizeinbase(&startParam.f->coef[degree], 2);
    maxThirdLog2 = mpz_sizeinbase(startParam.m,2) - lcLog2;
//    maxThirdLog2 = lcLog2 + (maxThirdLog2/degree)*3.15;
    maxThirdLog2 = lcLog2 + (maxThirdLog2/degree)*examineC3Lim;
    if (mpz_sizeinbase(&startParam.f->coef[degree-2],2) < maxThirdLog2) {
      optimizeParameters(&thisParam, &startParam);
      if (thisParam.logSize < minStage1Size)  {
        optimizeParameters2(&thisParam, &startParam);

        oldLogSize = thisParam.logSize;
        printTmp("k: %ld, log_2(I) = %1.2lf, computing E... (mTry=%d)",
                    iteration, thisParam.logSize, mTry);
        oldE = est_rating_skewed(thisParam.f, thisParam.m, mpf_get_d(thisParam.s));
        if (oldE > param->score) {
          mpz_poly_cp(param->f, thisParam.f);
          mpz_set(param->m, thisParam.m);
          param->logSize = thisParam.logSize;
          param->score = oldE;
        }
        mpz_poly_cp(oldF, thisParam.f);
        mpz_set(oldM, thisParam.m);
        printTmp("k: %ld, log_2(I) = %1.2lf, sieving... (mTry=%d)",
                    iteration, thisParam.logSize, mTry);
   
        /* This is the current best. */
        mpz_poly_cp(param->f, thisParam.f);
        mpz_set(param->m, thisParam.m);
        param->score = oldE;
        param->logSize = thisParam.logSize;
        mpf_set(param->s, thisParam.s);
  
        cont=1; try=0;
        sieveOverSimilarSize(&thisParam, J0, J1, logAdj, MAX_KEEP);
        for (i=0; i<MAX_KEEP; i++) {
          if (thisParam.logSize + logAdj[i] < minStage1Size) {
            adjustPolForParams(&testParam, &thisParam, J0[i], J1[i]);
            testParam.logSize += logAdj[i];
            minimize_wrt_t(&testParam);  /* Shift, if it will help. */
            minimize_wrt_s(&testParam);  /* Recompute the skew.     */
            if (mpf_get_d(testParam.s) < maxSkew) {
              printTmp("k: %ld, log_2(I) = %1.2lf, recomputing E... (mTry=%d)",
                       iteration, testParam.logSize, mTry);
              E = est_rating_skewed(testParam.f, testParam.m, mpf_get_d(testParam.s));
              if (E > param->score)  {
                mpz_poly_cp(param->f, testParam.f);
                mpz_set(param->m, testParam.m);
                param->logSize = testParam.logSize;
                param->score = E;
              }
            }
          }
        } 
      }
    } 
    if (mTry > 0) /* Some wierd stall condition sometimes.
                       I think it would be best to allow it
                       to continue, but it doesn't ever seem
                       to really help, so whatever. */
      tryNextM = 0;
  } while (tryNextM);

//  printf("Score = %e    \n\n", param->score);
  return param->score;
}

/******************************************************/
int parseInput(char *input)
/******************************************************/
{ char token[128], value[512], thisLine[1024];
  s32 size=strlen(input), loc;
  int i;

  loc=0;
  while (loc<size) {
    for (i=0; (i+loc)<size && (input[i+loc]!='\n'); i++)
      thisLine[i] = input[i+loc];
    thisLine[i]=0;
    loc += i+1;
    if (thisLine[0] != '#') {
      token[0] = value[0] = 0;
      sscanf(thisLine, "%128s %512s", token, value);
      if (strncmp(token, "n:", 2)==0) {
        mpz_set_str(N, value, 10);
      } else if (strncmp(token, "deg:", 4)==0)  {
        d = atoi(value);
      } else if (strncmp(token, "name:", 5)==0) {
        strncpy(name, value, MAX_NAMESIZE);
      } else if (strncmp(token, "if:", 3)==0) {
        strncpy(ifname, value, MAXFNAMESIZE);
      } else if (strncmp(token, "af:",3)==0) {
        strncpy(allName, value, MAXFNAMESIZE);
      } else if (strncmp(token, "bf:",3)==0) {
        strcpy(bestName, value);
      } else if (strncmp(token, "seed:",5)==0) {
        seed = atol(value);
      } else if (strncmp(token, "maxs1:",6)==0) {
        minStage1Size = userSetMin = atof(value);
      } else if (strncmp(token, "maxskew:",8)==0) {
        maxSkew = atof(value);
      } else if (strncmp(token, "lc1:",4)==0) {
        lc1 = atof(value);
      } else if (strncmp(token, "lcd:",4)==0) {
        lcd = MAX(atol(value), 1);
      } else if ((strncmp(token, "enum:",5)==0)) {
        mpz_set_str(enumLCD, value, 10);
      } else if (strncmp(token, "e0:",3)==0) {
        mpz_set_str(e0, value, 10);
      } else if (strncmp(token, "e1:",3)==0) {
        mpz_set_str(e1, value, 10);
      } else if (strncmp(token, "leave:",6)==0) {
        leave = atoi(value);
      } else if (strncmp(token, "lcp:",4)==0) {
        lcp = atoi(value);
      } else if (strncmp(token, "j0:",3)==0) {
        sieveJ0 = atoi(value);
      } else if (strncmp(token, "j1:",3)==0) {
        sieveJ1 = atoi(value);
      } else if (strncmp(token, "css:",4)==0) {
        contSampleSize = atol(value);
      } else if (strncmp(token, "cutoff:",4)==0) {
        cutoff = atof(value);
      } else if (strncmp(token, "examinefrac:",12)==0) {
        examineFrac = atof(value);
      }
    } 
  }
  if (mpz_cmp_ui(N, 1)==0) return -1;
  if (d<=1) {
    printf("Invalid degree (d=%d)!\n", d);
    return -1;
  }
  return 0;
}

/******************************************************/
char *getInputFromFile(char *fName)
/******************************************************/
/* Read the (ASCII formatted) file into a string and  */
/* return the result (with newlines and all).         */
/******************************************************/
{ struct stat fileInfo;
  FILE *fp;
  s32  size, index;
  char  *input, tmpStr[256];

  if (stat(fName, &fileInfo)) {
    fprintf(stderr, "Error: Could not stat input file %s.\n", fName);
    return "";
  }
  size = 512 + fileInfo.st_size;
  if (!(input = (char *)malloc(size*sizeof(char)))) {
    fprintf(stderr, "Memory allocation error (%" PRId32 " bytes for file %s).\n",
            size, fName);
    return "";
  }
  index=0;
  input[0]=0;
  fp = fopen(fName, "r");
  while (!(feof(fp))) {
    tmpStr[0]=0;
    fgets(tmpStr, 256, fp);
    if (strlen(tmpStr)) {
      strcat(input, tmpStr);
    }
  }
  fclose(fp);
  return input;
}

/******************************************************/
int scorePol(char *fName)
{ FILE *fp;
  static nfs_fb_t FB;
  static int initialized=0;
  int    res;
  double E, s;
  Iparam_t P;


  if (!initialized) {
    initFB(&FB); 
    initialized=1;
    Iparam_init(&P);
  }
  if (!(fp = fopen(fName, "r"))) {
    fprintf(stderr, "Error opening %s for read!\n", fName);
    return -1;
  }
  res = readPoly(fp, &FB);
  if (res) {
    fprintf(stderr, "Could not read a valid polynomial from %s! (res=%d)\n", 
             fName, res);
    return -1;
  }
  fclose(fp);
  maxSkew = 10e10;
  initialize(FB.f->degree);
  
  mpz_poly_cp(P.f, FB.f);
  mpz_set(P.m, FB.m);
  mpz_set(P.n, FB.n);
  mpf_set_d(P.s, FB.skew);
  s = minimize_wrt_s(&P);
  
  printf("If this is an SNFS number/poly, take these figures with a grain of salt!\n");
  printf("The computed skew will usually be close to optimal for the siever, but\n");
  printf("may differ! You should determine the best choice by experimentation.\n");
  E = est_rating_skewed(FB.f, FB.m, FB.skew);
  printf("Input polynomial has specified skew: %1.3lf\nE = %e\n", FB.skew, E);
  printf("Skew is computed to be: %1.3lf.\n", s);
  E = est_rating_skewed(FB.f, FB.m, s);
  printf("With the computed skew of %1.3lf:\nE = %e\n", s, E);
  
  return 0;
}

/******************************************************/
int isSNFS(mpz_t N)
/******************************************************/
/* Attempt to decide if this is an SNFS number, by    */
/* looking at its m-adic expansion for some small m.  */
/******************************************************/
/* Write me! */
{
  return 0;
}

/******************************************************/
void writeFactPars(FILE *fp, Iparam_t *P)
/* At some point, this function may make some default 
   choices for these paramaters. But for now, they will
   still have to be manually set by the user.
*/
{
  fprintf(fp, "# These parameters should be manually set:\n");
  fprintf(fp, "rlim: \n");
  fprintf(fp, "alim: \n");
  fprintf(fp, "lpbr: \n");
  fprintf(fp, "lpba: \n");
  fprintf(fp, "mfbr: \n");
  fprintf(fp, "mfba: \n");
  fprintf(fp, "rlambda: \n");
  fprintf(fp, "alambda: \n");
  fprintf(fp, "qintsize: \n");
}


/******************************************************/
int main(int argC, char *args[])
/******************************************************/
{ Iparam_t bestParam;
  char    *input, str[512];
  int      i;
  FILE    *fp;
  time_t   t;
  double   minLClog, maxLClog, bestScore=0, thisScore;
  printf(START_MSG, GGNFS_VERSION);

  seed=time(&t);
  mpz_init_set_ui(lastLC, 0);
  mpz_init_set_ui(N, 1);
  mpz_init_set_ui(enumLCD, 0);
  Iparam_init(&bestParam);
  ifname[0]=0;
  strcpy(allName, DEFAULT_ALLNAME);
  strcpy(bestName, DEFAULT_BESTNAME);
  sprintf(name, "unknown");
  scoreName[0]=0;
  mpz_init_set_ui(e0, 1);
  mpz_init_set_ui(e1, DEFAULT_ENUMSIZE);
  mpz_init_set(enumCt, e0);


  for (i=1; i<argC; i++) {
    if (strcmp(args[i], "-if")==0) {
      if ((++i) < argC)
        strcpy(ifname, args[i]);
    } else if (strcmp(args[i], "-score")==0) {
      if ((++i) < argC)
        strncpy(scoreName, args[i], 512);
    } else if (strcmp(args[i], "-css")==0) {
      if ((++i) < argC)
        contSampleSize = atol(args[i]);
    } else if (strcmp(args[i], "--help")==0) {
      printf("Usage: %s %s\n", args[0], USAGE);
      exit(0);
    }
  }
  /* Defaults which can be overridden: */
  srand(seed);
  prandseed(seed, 7313*seed+5, 2*seed+1);
  lc1 = MAX(lc1, 0.2);
  lc1 = MIN(lc1, 1.0); 
  /* Are we just scoring a single polynomial and nothing more? */
  if (scoreName[0]) {
    return scorePol(scoreName);
  }

  /* Get the input for this job. */
  /* We will add some network code here to allow the input to
     be obtained from a server in the same format.
  */
  input = getInputFromFile(ifname);
  parseInput(input);

  mpz_set(bestParam.n, N);
  if (mpz_cmp_ui(enumLCD,0) > 0) {
    printf("Enumerating leading coefficients c_d = k*");
    mpz_out_str(stdout, 10, enumLCD); printf(", ");
    mpz_out_str(stdout, 10, e0); printf(" <= k <= ");
    mpz_out_str(stdout, 10, e1); printf("\n");
    mpz_mul(lastLC, enumLCD, e0); mpz_sub(lastLC, lastLC, enumLCD);
    mpz_set(enumCt, e0);
  }
  if (ifname[0]==0) {
    printf("Usage: %s %s\n", args[0], USAGE);
    exit(0);
  }

  maxLClog = (1/(double)(d+1))*(double)mpz_sizeinbase(bestParam.n,2);
  minLClog = lc0*maxLClog;
  maxLClog = lc1*maxLClog;


  while (1) {
    iteration++;
    thisScore = doTest(minLClog, maxLClog, &bestParam, d);
    if (thisScore > bestScore) {
      /* Get a more accurate score to be sure. */
      contSampleSize *= 10;
      thisScore = bestParam.score = est_rating_skewed(bestParam.f, bestParam.m, mpf_get_d(bestParam.s));
      contSampleSize /= 10;
    }
    if (thisScore > bestScore) {
      printTmp(" ");
      printf("Score: %e (adj. I=%1.4lf, iteration %" PRId32 ", minStage1=%1.2lf)\n", 
             thisScore, bestParam.logSize, iteration, minStage1Size);
      bestScore = thisScore;

      if (!(fp = fopen(bestName, "w"))) {
        fprintf(stderr, "Error opening %s for write!\n", bestName);
        exit(-1);
      }
      fprintf(fp, "name: %s\n", name);
      fprintf(fp, "n:  "); mpz_out_str(fp, 10, bestParam.n); fprintf(fp, "\n");
      fprintf(fp, "m:  "); mpz_out_str(fp, 10, bestParam.m); 
      fprintf(fp, "\ndeg: %d\n", bestParam.f->degree);
      for (i=bestParam.f->degree; i>=0; i--) {
        fprintf(fp, "c%d: ", i);
        mpz_out_str(fp, 10, &bestParam.f->coef[i]);
        fprintf(fp, "\n");
      }
      fprintf(fp, "skew: %1.3lf\n", mpf_get_d(bestParam.s));
      fprintf(fp, "type: gnfs\n");
      fprintf(fp, "# adj. I(F,S) = %1.3lf\n", bestParam.logSize);
      fprintf(fp, "# E(F1,F2) = %e\n", bestParam.score);
      fprintf(fp, "# GGNFS version %s polyselect.\n", GGNFS_VERSION);
      fprintf(fp, "# Options were: \n");
      fprintf(fp, "# lcd=%" PRId32 ", enumLCD=%s, maxS1=%1.8lf, seed=%" PRId32 ".\n",
                   lcd, mpz_get_str(str, 10, enumLCD), minStage1Size, seed);
      fprintf(fp, "# maxskew=%1.1lf\n", maxSkew);
      writeFactPars(fp, &bestParam);
      fprintf(fp, "\n");
      fclose(fp);
    }
    if (thisScore > cutoff) {
     /* Append it to the `all file'. */
      if (!(fp = fopen(allName, "a"))) {
        fprintf(stderr, "Error opening %s for write!\n", bestName);
        exit(-1);
      }
      fprintf(fp, "name: %s\n", name);
      fprintf(fp, "n: "); mpz_out_str(fp, 10, bestParam.n); fprintf(fp, "\n");
      fprintf(fp, "m:  "); mpz_out_str(fp, 10, bestParam.m); 
      fprintf(fp, "\ndeg: %d\n", bestParam.f->degree);
      for (i=bestParam.f->degree; i>=0; i--) {
        fprintf(fp, "c%d: ", i);
        mpz_out_str(fp, 10, &bestParam.f->coef[i]);
        fprintf(fp, "\n");
      }
      fprintf(fp, "skew: %1.3lf/\n", mpf_get_d(bestParam.s));
      fprintf(fp, "type: gnfs\n");
      fprintf(fp, "# adj. I(F,S) = %1.3lf\n", bestParam.logSize);
      fprintf(fp, "# E(F1,F2) = %e\n", bestParam.score);
      fprintf(fp, "\n");
      fclose(fp);
    } 
  }

  return 0;
}


