/**************************************************************/
/* mpz_poly.c                                                 */
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

/********************************************************/
/* Functions for performing various tasks with our      */
/* number field defining polynomial, f.                 */
/* In particular, there are functions here to do        */
/* various types of evaluation, find the complex zeros, */
/* and find the base-m representation.                  */
/********************************************************/
/* Throughout,                                          */
/* 'f' is the polynomial f(x)                           */
/* 'F' is the homogenized version: F(x,y) = (y^d)f(x/y).*/
/* 'g' is the monic version, g(x) = c_d^{d-1}f(x/c_d).  */
/********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "ggnfs.h"



/****************************************/
void mpz_poly_init(mpz_poly f)
/****************************************/
{ int i;
	
  for (i=MAXPOLYDEGREE; i>=0; i--)
    mpz_init(&(f->coef[i]));
  f->degree=0;
}

/****************************************/
void mpz_poly_clear(mpz_poly f)
/****************************************/
{ int i;
	
  for (i=MAXPOLYDEGREE; i>=0; i--)
    mpz_clear(&(f->coef[i]));
  f->degree=0;
}

/************************************************/
void mpz_poly_cp(mpz_poly dest, mpz_poly src)
/************************************************/
{ int i;

  for (i=dest->degree=src->degree; i>=0; i--)
    mpz_set(&dest->coef[i], &src->coef[i]);
}
	
/************************************************/
int mpz_poly_fixDeg(mpz_poly f)
/************************************************/
{ int d = f->degree;

  while ((d>0) && (mpz_sgn(&f->coef[d])==0))
    d--;
  f->degree = d;
  return d;
}
	
	
/************************************************/
int mpz_poly_diff(mpz_poly res, mpz_poly q)
/************************************************/
/* Compute res <-- q'.                          */
/************************************************/
{ int i;
  
#ifdef _OLD 
  d = q->degree - 1;
  if (d < 0) {
    res->degree = 0;
    mpz_set_ui(&res->coef[0], 0);
    return 0;
  }
  res->degree = d;
  mpz_init(tmp);
  for (i=0; i<=d; i++) {
    mpz_set_si(tmp, (s32)(i+1));
    j=i+1;
    mpz_mul_ui(tmp, tmp, &q->coef[j]);
    mpz_set(&res->coef[i], tmp);
  }
  mpz_clear(tmp);
#else
  if (q->degree <= 0) {
    mpz_set_ui(&res->coef[0], 0);
    res->degree = 0;
    return 0;
  } 

  mpz_poly_cp(res, q);
  for (i=0; i<res->degree; i++)
    mpz_mul_ui(&res->coef[i], &res->coef[i+1], (i+1));
  res->degree -= 1;
#endif

  return 0;
}


/************************************************/
/* This structure is only used locally. */
typedef struct {
  __mpf_struct c_r[MAXPOLYDEGREE+1];
  __mpf_struct c_i[MAXPOLYDEGREE+1];
  int degree;
} mpf_cpoly_t;
typedef mpf_cpoly_t mpf_cpoly[1];

/**** Prototypes for locally used functions. ****/
void   initCPoly(mpf_cpoly g);
void   clearCPoly(mpf_cpoly g);
void   evalCPoly(mpf_t vr, mpf_t vi, mpf_cpoly Q, mpf_t xr, mpf_t xi);
void   mulC(mpf_t yr, mpf_t yi, mpf_t ar, mpf_t ai, mpf_t br, mpf_t bi);
void   divC(mpf_t yr, mpf_t yi, mpf_t ar, mpf_t ai, mpf_t br, mpf_t bi);
double absC_d(mpf_t re, mpf_t im);
/************************************************/



/******************************************************/
int mpz_poly_getComplexZeros(nfs_complex_t *Z, mpz_poly f)
/******************************************************/
/* Get approximations to the zeros of 'F'. Of course, */
/* there will be deg(F) of them, so 'Z' should be     */
/* allocated accordingly.                             */
/******************************************************/
/* This is Algorithm 3.6.6 (Complex Roots) from       */
/* H. Cohen, ``A Course in Computational Algebraic    */
/* Number Theory''.                                   */
/******************************************************/
/* Note: The defining polynomial is probably well     */
/* behaved enough that we could do most of this with  */
/* just double's. But since we link in GMP anyway,    */
/* what the heck.                                     */
/******************************************************/
{ mpf_cpoly P, dP, Q, dQ;
  mpf_t xr, xi, vr, vi, dxr, dxi, m, m1, t1, t2;
  mpf_t yr, yi, v1r, v1i, a, b;
  int n, i, c, res=0, prec;

  prec = 512; /* Should be plenty. */
  /* Initialize the various GMP variables. */  
  initCPoly(P); initCPoly(dP); initCPoly(Q); initCPoly(dQ);
  mpf_init2(xr, prec);  mpf_init2(xi, prec);
  mpf_init2(yr, prec);  mpf_init2(yi, prec);
  mpf_init2(dxr, prec); mpf_init2(dxi, prec);
  mpf_init2(vr, prec);  mpf_init2(vi, prec);
  mpf_init2(v1r, prec); mpf_init2(v1i, prec);
  mpf_init2(t1, prec);  mpf_init2(t2, prec);
  mpf_init2(m, prec);   mpf_init2(m1, prec);
  mpf_init2(a, prec);   mpf_init2(b, prec);
  
  n = f->degree;
  P->degree = n;
  Q->degree = n;
  /* P <-- F, Q <-- f */
  for (i=0; i<=n; i++) {
    mpf_set_z(&(P->c_r[i]), &(f->coef[i]));
    mpf_set_ui(&(P->c_i[i]), 0);
    mpf_set(&(Q->c_r[i]), &(P->c_r[i]));
    mpf_set(&(Q->c_i[i]), &(P->c_i[i]));
  }
  dP->degree = n-1;
  /* dP <-- P' */
  for (i=0; i<=(n-1); i++) {
    mpf_set(&(dP->c_r[i]), &(P->c_r[i+1]));
    mpf_mul_ui(&(dP->c_r[i]), &(dP->c_r[i]), (i+1));
    mpf_set_ui(&(dP->c_i[i]), 0);
  }
  
  /* Step 2 */
ALG3_6_6_STEP2:
  /* do DQ <-- Q'. */
  dQ->degree = n-1;
  for (i=0; i<=dQ->degree; i++) {
    mpf_mul_ui(&(dQ->c_r[i]), &(Q->c_r[i+1]), i+1);
    mpf_mul_ui(&(dQ->c_i[i]), &(Q->c_i[i+1]), i+1);
  }
  
  /* Arbitrary - but this is Cohen's suggestion, since  */
  /* it is not too close to a trivial algebraic number. */
  mpf_set_d(xr, 1.3);
  mpf_set_d(xi, 0.314159);
  evalCPoly(vr, vi, Q, xr, xi);
  mpf_mul(m, vr, vr); mpf_mul(t1, vi, vi); mpf_add(m, m, t1);

  /* Step 3: */
ALG3_6_6_STEP3:
  c=0;
  evalCPoly(t1, t2, dQ, xr, xi);
  if ((fabs(mpf_get_d(t1))>0.0) || (fabs(mpf_get_d(t1))>0.0)) {
    divC(dxr, dxi, vr, vi, t1, t2);
    if (absC_d(dxr, dxi) < 0.0000000001) 
      goto ALG3_6_6_STEP5;
  } else {
    printf("ALG3_6_6_STEP3 failed: Q'(x) = 0!\n");
    return -1;
  }

  /* Step 4: */
ALG3_6_6_STEP4:
  mpf_sub(yr, xr, dxr);
  mpf_sub(yi, xi, dxi);
  evalCPoly(v1r, v1i,  Q, yr, yi);
  mpf_mul(m1, v1r, v1r); mpf_mul(t1, v1i, v1i); mpf_add(m1, m1, t1);
  if (mpf_cmp(m1, m) < 0) {
    mpf_set(xr, yr); mpf_set(xi, yi);
    mpf_set(vr, v1r); mpf_set(vi, v1i);
    mpf_set(m, m1);
    goto ALG3_6_6_STEP3;
  }
  c++;
  mpf_div_ui(dxr, dxr, 4); mpf_div_ui(dxi, dxi, 4);
  if (c<20)
    goto ALG3_6_6_STEP4;
  else {
    res = -1;
    goto ALG3_6_6_STOP;
  }

  

ALG3_6_6_STEP5:
  /* Three Newton iterations. */

  /* Note: We should really be a bit more clever here,
  *  and look at the difference between the last two
  *  Newton iterations to make a conclusion about
  *  how many bits can be presumed accurate. But
  *  I think this should give plenty of accuracy
  *  for now.
  */
  
  for (i=0; i<3; i++) {
    evalCPoly(t1, t2, P, xr, xi);
    evalCPoly(yr, yi, dP, xr, xi);
    divC(vr, vi, t1, t2, yr, yi);
    mpf_sub(xr, xr, vr); mpf_sub(xi, xi, vi);
  }
  if (fabs(mpf_get_d(xi)) < 0.000000000001) {
    mpf_set_d(xi, 0.0);
    mpf_set(Z[n-1].mpr, xr);
    mpf_set(Z[n-1].mpi, xi);
    Z[n-1].r = mpf_get_d(xr);
    Z[n-1].i = 0.0;
    /* Q <-- Q/(X - (xr,xi)). */
    /* We will use dQ for intermediate storage, since it */
    /* will be recomputed anyway.                        */
    mpf_set(&(dQ->c_r[n-1]), &(Q->c_r[n]));
    mpf_set(&(dQ->c_i[n-1]), &(Q->c_i[n]));
    for (i=n-1; i>=1; i--) {
      mulC(&(dQ->c_r[i-1]), &(dQ->c_i[i-1]), xr, xi, &(dQ->c_r[i]), &(dQ->c_i[i]));
      mpf_add(&(dQ->c_r[i-1]), &(dQ->c_r[i-1]),  &(Q->c_r[i]));
      mpf_add(&(dQ->c_i[i-1]), &(dQ->c_i[i-1]),  &(Q->c_i[i]));
    }
    n--;
    for (i=0; i<=n; i++) {
      mpf_set(&(Q->c_r[i]), &(dQ->c_r[i]));
      mpf_set(&(Q->c_i[i]), &(dQ->c_i[i]));
    }
    Q->degree = n;
  } else {
    Z[n-1].r = mpf_get_d(xr);
    Z[n-1].i = mpf_get_d(xi);
    mpf_set(Z[n-1].mpr, xr);
    mpf_set(Z[n-1].mpi, xi);
    Z[n-2].r = Z[n-1].r;
    Z[n-2].i = -Z[n-1].i;
    mpf_set(Z[n-2].mpr, Z[n-1].mpr);
    mpf_neg(Z[n-2].mpi, Z[n-1].mpi);
//    mpf_neg(Z[n-2].mpi, Z[n-2].mpi);
    /* Now, Q <-- Q/(X^2 - 2*Re(x)X + |x|^2). */
    mpf_add(a, xr, xr);
    mpf_mul(b, xr, xr); mpf_mul(t1, xi, xi); mpf_add(b, b, t1);
    
    mpf_set(&(dQ->c_r[n-2]), &(Q->c_r[n]));      mpf_set(&(dQ->c_i[n-2]), &(Q->c_i[n]));
    if (n>2) {
      mpf_mul(&(dQ->c_r[n-3]), &(dQ->c_r[n-2]), a); mpf_mul(&(dQ->c_i[n-3]), &(dQ->c_i[n-2]), a);
      mpf_add(&(dQ->c_r[n-3]), &(dQ->c_r[n-3]), &(Q->c_r[n-1]));
      mpf_add(&(dQ->c_i[n-3]), &(dQ->c_i[n-3]), &(Q->c_i[n-1]));
    }
    for (i=n-2; i>=2; i--) {
      mpf_mul(&(dQ->c_r[i-2]), a, &(dQ->c_r[i-1])); mpf_mul(&(dQ->c_i[i-2]), a, &(dQ->c_i[i-1])); 
      mpf_add(&(dQ->c_r[i-2]), &(dQ->c_r[i-2]), &(Q->c_r[i]));
      mpf_add(&(dQ->c_i[i-2]), &(dQ->c_i[i-2]), &(Q->c_i[i]));
      mpf_mul(t1, b, &(dQ->c_r[i]));
      mpf_mul(t2, b, &(dQ->c_i[i]));
      mpf_sub(&(dQ->c_r[i-2]), &(dQ->c_r[i-2]), t1);
      mpf_sub(&(dQ->c_i[i-2]), &(dQ->c_i[i-2]), t2);
    }
    n -= 2;
    for (i=0; i<=n; i++) {
      mpf_set(&(Q->c_r[i]), &(dQ->c_r[i]));
      mpf_set(&(Q->c_i[i]), &(dQ->c_i[i]));
    }
    Q->degree = n;
  }
  if (n==1) {
    divC(xr, xi, &(Q->c_r[0]), &(Q->c_i[0]), &(Q->c_r[1]), &(Q->c_i[1]));
    Z[0].r = -mpf_get_d(xr);
    Z[0].i = -mpf_get_d(xi);
//    mpf_set(Z[0].mpr, xr);
//    mpf_neg(Z[0].mpr, Z[0].mpr);
//    mpf_set(Z[0].mpi, xi);
//    mpf_neg(Z[0].mpi, Z[0].mpi);
    mpf_neg(Z[0].mpr, xr);
    mpf_neg(Z[0].mpi, xi);
  } else if (n>0)
    goto ALG3_6_6_STEP2;
  res = 0;  

ALG3_6_6_STOP:
  
  
  clearCPoly(P); clearCPoly(dP); clearCPoly(Q); clearCPoly(dQ);
  mpf_clear(xr);  mpf_clear(xi);
  mpf_clear(yr);  mpf_clear(yi);
  mpf_clear(dxr); mpf_clear(dxi);
  mpf_clear(vr);  mpf_clear(vi);
  mpf_clear(v1r); mpf_clear(v1i);
  mpf_clear(t1);  mpf_clear(t2);
  mpf_clear(m);   mpf_clear(m1);
  mpf_clear(a);   mpf_clear(b);
  return res;
  
}


/****************************************/
void mpz_poly_print(FILE *fp, char *str, mpz_poly f)
{ int i, flag, d=f->degree;

  fprintf(fp, "%s", str);
  if ((d==0) && mpz_sgn(&f->coef[0])==0) {
    fprintf(fp, "0\n");
    return ;
  }
  for (i=flag=0; i<=d; i++) {
    if (mpz_sgn(&f->coef[i])) {
      if (flag && mpz_sgn(&f->coef[i]) > 0) {
        fprintf(fp, " + ");
        mpz_out_str(fp, 10, &(f->coef[i]));
      } else {
        mpz_out_str(fp, 10, &(f->coef[i]));
	flag=1;
      }
      if (i>0) 
        fprintf(fp, "X");
      if (i>1)
        fprintf(fp, "^%d", i);
    }
  }
  fprintf(fp, "\n");
}
     

  
/***************************************************************/
int mpz_poly_getBaseM(mpz_poly f, mpz_t n, mpz_t m, int init_poly)
/***************************************************************/
/* Input: Integers n, m and init_poly.                         */
/* Output: 'f' will be the base-m representation of 'n',       */
/*         with integer coefficients in [-m/2, m/2].           */
/*         If 'init_poly' is nonzero, 'f' will be initialized. */
/*         Note: If 'f' has already been initialized, make     */
/*         init_poly zero, or you'll be leaking memory!.       */
/* Return value: 0, if we succeeded. Nonzero on error.         */
/***************************************************************/
{ int i, res=0;
  static mpz_t temp, half_m;
  static int initialized=0;

  if (mpz_sgn(m)==0) {
    fprintf(stderr, "mpz_poly_getBaseM(): Division by zero!!\n");
    return -1;
  }
  if (!(initialized)) {
    mpz_init(temp);
    mpz_init(half_m);
    initialized=1;
  }

  /*************************************/ 
  /* initialize the mpz_t coefficients */
  /*************************************/ 
  if (init_poly) 
    for (i=MAXPOLYDEGREE; i>=0; i--)
      mpz_init(&(f->coef[i]));
    
  mpz_set(temp, n);

  i=0;

  while ((i<=MAXPOLYDEGREE)&&(mpz_sgn(temp)>0)) {
    mpz_fdiv_qr(temp, &(f->coef[i]), temp, m);
    i++;
  }
  f->degree = i-1;

  if (mpz_sgn(temp))
    res = -1;
  
  /**************************************/
  /* Correct for any coefficients > m/2 */
  /**************************************/
  mpz_div_ui(half_m, m, 2);
  for (i=0; i< f->degree; i++) {
    if (mpz_cmp(&(f->coef[i]), half_m) > 0) {
      mpz_sub(&(f->coef[i]), &(f->coef[i]), m);
      mpz_add_ui(&(f->coef[i+1]), &(f->coef[i+1]), 1);
    }
  }
  return res;
}

/***************************************************************/
int mpz_poly_getBaseM2(mpz_poly f, mpz_t n, mpz_t m, mpz_t modulus, mpz_t residue)
/***************************************************************/
/* Same as above, but insist that the leading coefficient be   */
/* ==residue (mod modulus).                                    */
/***************************************************************/
{ int i, res=0;
  static mpz_t temp, half_m, tmp1, tmp2, mPow;
  static int initialized=0;

  if (mpz_sgn(m)==0) {
    fprintf(stderr, "mpz_poly_getBaseM(): Division by zero!!\n");
    return -1;
  }
  if (!(initialized)) {
    mpz_init(temp); mpz_init(tmp1); mpz_init(tmp2);
    mpz_init(half_m); mpz_init(mPow);
    initialized=1;
  }

  i=0;
  mpz_set(mPow, m); i=1;
  while (mpz_cmp(mPow, n) < 0) {
    mpz_mul(mPow, mPow, m);
    i++;
  }
  f->degree = --i;
  mpz_div(mPow, mPow, m); /* temp <-- m^d. */
  mpz_div(&f->coef[i], n, mPow);
  mpz_mod(tmp1, &f->coef[i], modulus);
  mpz_sub(tmp1, residue, tmp1);
  mpz_add(&f->coef[i], &f->coef[i], tmp1);
  mpz_mul(tmp2, mPow, &f->coef[i]);
  mpz_sub(temp, n, tmp2);
  i--;
  mpz_div(mPow, mPow, m);
  while ((i>=0) && (mpz_sgn(temp))) {
    mpz_div(&f->coef[i], temp, mPow);
    mpz_mul(tmp1, mPow, &f->coef[i]);
    mpz_sub(temp, temp, tmp1);
    mpz_div(mPow, mPow, m);
    i--;
  }

  if (mpz_sgn(temp))
    res = -1;
  
  /**************************************/
  /* Correct for any coefficients > m/2 */
  /**************************************/
  mpz_div_ui(half_m, m, 2);
  for (i=0; i< (f->degree-1); i++) {
    if (mpz_cmp(&(f->coef[i]), half_m) > 0) {
      mpz_sub(&(f->coef[i]), &(f->coef[i]), m);
      mpz_add_ui(&(f->coef[i+1]), &(f->coef[i+1]), 1);
    }
  }
  return res;
}



/******************************************************/
void initCPoly(mpf_cpoly g)
/******************************************************/
{ int i;

  g->degree = 0;
  for (i=0; i<MAXPOLYDEGREE+1; i++) {
    mpf_init2(&(g->c_r[i]), 128); /* This should be plenty. */
    mpf_init2(&(g->c_i[i]), 128); 
  }
  /* Note: We are not expecting to factor numbers  */
  /* larger than, say, 2^{600}. With a 5-th degree */
  /* poly, this means coefficients in the range of */
  /* 2^{100}, so 128 is plenty.                    */
}

/******************************************************/
void clearCPoly(mpf_cpoly g)
/******************************************************/
{ int i;

  g->degree = 0;
  for (i=0; i<MAXPOLYDEGREE+1; i++) {
    mpf_clear(&(g->c_r[i])); 
    mpf_clear(&(g->c_i[i])); 
  }
}


/******************************************************/
void evalCPoly(mpf_t vr, mpf_t vi, mpf_cpoly Q, mpf_t xr, mpf_t xi)
/******************************************************/
/* Compute (vr, vi) <-- Q(xr, xi).                    */
/******************************************************/
{ int i;
  mpf_t tmp_r, tmp_i, t1, t2, nr;

  mpf_init2(tmp_r, 128); mpf_init2(tmp_i, 128);
  mpf_init2(t1, 128);    mpf_init2(t2, 128);
  mpf_init2(nr, 128);
  
  mpf_set(tmp_r, xr); mpf_set(tmp_i, xi);
  mpf_set(vr, &(Q->c_r[0]));
  mpf_set(vi, &(Q->c_i[0]));
  for (i=1; i<=Q->degree; i++) {
    /* Compute the i-th term of Q(xr, xi). */
    mpf_mul(t1, &(Q->c_r[i]), tmp_r);
    mpf_mul(t2, &(Q->c_i[i]), tmp_i);
    mpf_sub(t1, t1, t2);
    mpf_add(vr, vr, t1);
    mpf_mul(t1, &(Q->c_i[i]), tmp_r);
    mpf_mul(t2, &(Q->c_r[i]), tmp_i);
    mpf_add(t1, t1, t2);
    mpf_add(vi, vi, t1);
    if (i<Q->degree) {
      /* Do (tmp_r, tmp_i) <-- (tmp_r, tmp_i)*(xr, xi) */
      mpf_mul(t1, tmp_r, xr);
      mpf_mul(t2, tmp_i, xi);
      mpf_sub(nr, t1, t2);
      mpf_mul(t1, tmp_r, xi);
      mpf_mul(t2, tmp_i, xr);
      mpf_add(tmp_i, t1, t2);
      mpf_set(tmp_r, nr);
    }
  }
  mpf_clear(tmp_r); mpf_clear(tmp_i);
  mpf_clear(t1); mpf_clear(t2);
  mpf_clear(nr);
}
    
    
/******************************************************/
void mulC(mpf_t yr, mpf_t yi, mpf_t ar, mpf_t ai, mpf_t br, mpf_t bi)
/******************************************************/
/* Compute (yr,yi) <-- (ar,ai)*(br,bi).               */
/******************************************************/
{ mpf_t t1;

  mpf_init2(t1, mpf_get_prec(ar));
  mpf_mul(yr, ar, br);
  mpf_mul(t1, ai, bi);
  mpf_sub(yr, yr, t1);

  mpf_mul(yi, ai, br);
  mpf_mul(t1, ar, bi);
  mpf_add(yi, yi, t1);
  mpf_clear(t1);
}
        
#ifdef _THIS_METHOD_SUCKS
/* I don't know where I read that this is supposed to be
 * more stable, but in our application, it sucks really
 * bad. It took me forever to figure out what was going
 * wrong. Perhaps I didn't implement it right. In either
 * case, the straightforward division technique works
 * much better.
*/
/******************************************************/
void divC(mpf_t yr, mpf_t yi, mpf_t ar, mpf_t ai, mpf_t br, mpf_t bi)
/******************************************************/
/* Compute (yr,yi) <-- (ar,ai)/(br,bi).               */
/******************************************************/
/* We do it a little better than the obvious way,     */
/* which is prone to error when |br| and |bi| are far */
/* apart.                                             */
/******************************************************/
{ mpf_t h, d;

  mpf_init2(d, mpf_get_prec(br)); mpf_init2(h, mpf_get_prec(br));
  if (fabs(mpf_get_d(bi)) < fabs(mpf_get_d(br))) {
    if (fabs(mpf_get_d(bi)) < 0.00000000000000001) {
      mpf_div(yr, ar, br);
      mpf_div(yi, ai, br);
    } else {
      mpf_div(h, bi, br);
      mpf_mul(d, h, bi); mpf_add(d, d, br);
      mpf_mul(yr, h, ai); mpf_add(yr, yr, ar); mpf_div(yr, yr, d);
      mpf_mul(yi, h, ar); mpf_sub(yi, ai, yi); mpf_div(yi, yi, d);
    }
//  } else if (fabs(mpf_get_d(bi)) > 0) {
  } else if (mpf_sgn(bi) > 0) {
    mpf_div(h, br, bi);
    mpf_mul(d, h, br); mpf_add(d, d, bi);
    mpf_mul(yr, ar, h); mpf_add(yr, yr, yi); mpf_div(yr, yr, d);
    mpf_mul(yi, ai, h); mpf_sub(yi, yi, ar); mpf_div(yi, yi, d);
  } else {
    printf("divC() Divide by zero!\n");
  }
  mpf_clear(d); mpf_clear(h);
}   
#else
/*************************************************************/
void divC(mpf_t rx, mpf_t ry, mpf_t ax, mpf_t ay, mpf_t bx, mpf_t by)
/*************************************************************/
{ mpf_t d, t1, t2;
  
  mpf_init2(d, mpf_get_prec(ax));
  mpf_init2(t1, mpf_get_prec(ax));
  mpf_init2(t2, mpf_get_prec(ax));
  mpf_mul(d, bx, bx);
  mpf_mul(t1, by, by);
  mpf_add(d, d, t1);

  mpf_div(t1, bx, d);
  mpf_div(t2, by, d);
  mpf_neg(t2, t2);

  /* Now just multiply: (rx, ry) <-- (ax, ay)(t1, t2). */
  mpf_mul(rx, ax, t1);  mpf_mul(d, ay, t2);  mpf_sub(rx, rx, d);
  mpf_mul(ry, ax, t2); mpf_mul(d, ay, t1); mpf_add(ry, ry, d);

  mpf_clear(d);
  mpf_clear(t1);
  mpf_clear(t2);
}
#endif




/************************************/
double absC_d(mpf_t re, mpf_t im)
/************************************/
{ double r, i;

  /* Close enough. */
  r = mpf_get_d(re);
  i = mpf_get_d(im);
  return r*r + i*i;
}


/************************************/
int mpz_poly_discrim(mpz_t disc, mpz_poly A)
/************************************/
/* Compute the discriminant of 'A'. */
/************************************/
{ mpz_poly dA;
  int i, m, sign;
  
  mpz_poly_init(dA);

  m = A->degree;
  /* Following Cohen,  compute */
  /* disc <-- (-1)^(m(m-1)/2) * R(A, A')/l(A). */
  
  mpz_poly_diff(dA, A); /* dA <-- A'. */
  i = m*(m-1)/2;
  sign = (i%2) ? -1 : 1;
  
  mpz_poly_resultant(disc, A, dA);
  mpz_div(disc, disc, &(A->coef[m]));
  mpz_mul_si(disc, disc, sign);

  mpz_poly_clear(dA);
  return 0;
}

/************************************/
void mpz_poly_resultant(mpz_t res, mpz_poly _A, mpz_poly _B)
/************************************/
/* Cohen, Alg. 3.3.7.               */
/************************************/
{ mpz_poly A, B, Q, R;
  mpz_t    a, b, g, h, s, t, tmp1;
  int      i, delta;
  
  if ((_A->degree == 0) && (mpz_cmp_ui(&(_A->coef[0]),0)==0)) {
    mpz_set_ui(res, 0);
    return;
  }
  if ((_B->degree == 0) && (mpz_cmp_ui(&(_B->coef[0]),0)==0)) {
    mpz_set_ui(res, 0);
    return;
  }
  
  mpz_poly_init(A);
  mpz_poly_init(B);
  mpz_poly_init(Q);
  mpz_poly_init(R);
  mpz_init(a); mpz_init(b);
  mpz_init(g); mpz_init(h);
  mpz_init(s); mpz_init(t);
  mpz_init(tmp1);
  
  for (i=0; i<= _A->degree; i++)
    mpz_set(&(A->coef[i]), &(_A->coef[i]));
  for (i=0; i<= _B->degree; i++)
    mpz_set(&(B->coef[i]), &(_B->coef[i]));
  A->degree = _A->degree;
  B->degree = _B->degree;

  /* a <-- cont(A). */
  mpz_set(a, &(A->coef[0]));
  for (i=1; i<=A->degree; i++)
    mpz_gcd(a, a, &(A->coef[i]));
    
  /* b <-- cont(B). */
  mpz_set(b, &(B->coef[0]));
  for (i=1; i<=B->degree; i++)
    mpz_gcd(b, b, &(B->coef[i]));

  /* A <-- A/a,  B <-- B/b. */
  for (i=0; i<=A->degree; i++)
    mpz_div(&(A->coef[i]), &(A->coef[i]), a);
  for (i=0; i<=B->degree; i++)
    mpz_div(&(B->coef[i]), &(B->coef[i]), b);

  mpz_set_ui(g, 1); mpz_set_ui(h, 1); mpz_set_ui(s, 1);
  mpz_pow_ui(t, a, B->degree);
  mpz_pow_ui(tmp1, b, A->degree);
  mpz_mul(t, t, tmp1);

  if (A->degree < B->degree) {
    for (i=0; i<=MAXPOLYDEGREE; i++) {
      mpz_set(tmp1, &(A->coef[i]));
      mpz_set(&(A->coef[i]), &(B->coef[i]));
      mpz_set(&(B->coef[i]), tmp1);
    }
    i = A->degree; A->degree = B->degree; B->degree = i;
    if ((A->degree%2)&&(B->degree%2))
      mpz_set_si(s, -1);
  }

ALG_3_3_7_STEP2:
  delta = A->degree - B->degree;
  if ((A->degree%2)&&(B->degree%2))
    mpz_neg(s, s);
  
  mpz_poly_psuedoDiv(Q, R, A, B);

  /********** Step 3 ***********/
  /* A <-- B, B <-- R/(g*h^{delta}). */
  A->degree = B->degree;
  for (i=0; i<=B->degree; i++)
    mpz_set(&(A->coef[i]), &(B->coef[i]));
  mpz_pow_ui(tmp1, h, delta);
  mpz_mul(tmp1, tmp1, g);
  B->degree = R->degree;
  for (i=0; i<=R->degree; i++) 
    mpz_div(&(B->coef[i]), &(R->coef[i]), tmp1);
  
  
  /********** Step 4 ***********/
  mpz_set(g, &(A->coef[A->degree]));
  mpz_pow_ui(tmp1, g, delta);
  for (i=delta-1; i>0; i--)
    mpz_div(tmp1, tmp1, h);
  mpz_set(h, tmp1);

  if (B->degree > 0) {
    goto ALG_3_3_7_STEP2;
  }

  mpz_pow_ui(tmp1, &(B->coef[B->degree]), A->degree);
  for (i=A->degree-1; i>0; i--)
    mpz_div(tmp1, tmp1, h);
  mpz_set(h, tmp1);

  mpz_mul(res, s, t); mpz_mul(res, res, h);
  
  
  mpz_poly_clear(A);
  mpz_poly_clear(B);
  mpz_poly_clear(Q);
  mpz_poly_clear(R);
  mpz_clear(a); mpz_clear(b);
  mpz_clear(g); mpz_clear(h);
  mpz_clear(s); mpz_clear(t);
  mpz_clear(tmp1);
}

/**************************************************************/
void mpz_poly_psuedoDiv(mpz_poly Q, mpz_poly R, mpz_poly A, mpz_poly B)
/**************************************************************/
{ int m=A->degree, n=B->degree, e, i, s_deg;
  mpz_t d, tmp, s;

  mpz_init(tmp);
  mpz_init(s);
  mpz_init_set(d, &(B->coef[n]));
  
  R->degree = A->degree;
  for (i=0; i<=A->degree; i++)
    mpz_set(&(R->coef[i]), &(A->coef[i]));
  Q->degree = 0;
  mpz_set_ui(&(Q->coef[0]), 0);
  e = m-n+1;

ALG_3_1_2_STEP2:
  if (R->degree < B->degree) {
    mpz_pow_ui(tmp, d, e);
    for (i=0; i<=Q->degree; i++)
      mpz_mul(&(Q->coef[i]), &(Q->coef[i]), tmp);
    for (i=0; i<=R->degree; i++)
      mpz_mul(&(R->coef[i]), &(R->coef[i]), tmp);
    mpz_clear(d);
    mpz_clear(s);
    mpz_clear(tmp);
    return;
  }

  /* Step 3: */
  s_deg = R->degree - B->degree;
  mpz_set(s, &(R->coef[R->degree]));

  /* Q <-- d*Q + s*x^{s_deg} */
  for (i=0; i<=Q->degree; i++)
    mpz_mul(&(Q->coef[i]), &(Q->coef[i]), d);
  for (i=Q->degree+1; i<=s_deg; i++)
    mpz_set_ui(&(Q->coef[i]), 0);
  Q->degree = MAX(s_deg, Q->degree);
  mpz_add(&(Q->coef[s_deg]), &(Q->coef[s_deg]), s);
  while ((Q->degree>0) && (mpz_cmp_ui(&(Q->coef[Q->degree]), 0)==0))
    Q->degree -= 1;
    
  /* R <-- d*R - B*s*x^{s_deg} */
  for (i=0; i<=R->degree; i++)
    mpz_mul(&(R->coef[i]), &(R->coef[i]), d);
  for (i=R->degree+1; i<=(s_deg+B->degree); i++)
    mpz_set_ui(&(R->coef[i]), 0);

  for (i=0; i<=B->degree; i++) {
    mpz_mul(tmp, &(B->coef[i]), s);
    mpz_sub(&(R->coef[s_deg+i]), &(R->coef[s_deg+i]), tmp);
  }
  R->degree = MAX(R->degree, s_deg+B->degree);
  while ((R->degree>0) && (mpz_cmp_ui(&(R->coef[R->degree]), 0)==0))
    R->degree -= 1;
  e--;
  goto ALG_3_1_2_STEP2;
  
}
	
/**************************************************************/
int mpz_poly_div(mpz_poly _q, mpz_poly f, mpz_poly g, mpz_t p)
/**************************************************************/
/* Compute q so that q*g = f (mod p). Result is undefined if  */
/* g does not divide f.                                       */
/**************************************************************/
{ int      i, rd;
  static mpz_poly r, q;
  static mpz_t    tmp1, tmp2, lgi;
  static int initialized=0;

  if (!(initialized)) {
    mpz_poly_init(r);
    mpz_poly_init(q);
    mpz_init(tmp1);
    mpz_init(tmp2);
    mpz_init(lgi);
    initialized=1;
  }

  /* Special case: */
  if (g->degree == 0) {
    if (mpz_sgn(&g->coef[0])==0) {
      fprintf(stderr, "mpz_poly_div() Division by zero!\n");
      return -1;
    }
    mpz_invert(lgi, &g->coef[0], p);
    for (i=f->degree; i>=0; i--) {
      mpz_mul(&_q->coef[i], &f->coef[i], lgi);
      mpz_mod(&_q->coef[i], &_q->coef[i], p);
    }
    _q->degree = f->degree;
    return 0;
  }

    
  
  mpz_poly_cp(r, f);
  q->degree = f->degree - g->degree;
  for (i=q->degree; i>=0; i--)
    mpz_set_ui(&q->coef[i], 0);
  
  mpz_invert(lgi, &g->coef[g->degree], p);
  
    
	  
  
  rd = r->degree;
  while (rd >= g->degree) {
    mpz_mul(tmp1, lgi, &r->coef[rd]);
    mpz_mod(tmp1, tmp1, p);
    mpz_set(&q->coef[rd - g->degree], tmp1);
    for (i=g->degree; i>=0; i--) {
      mpz_mul(tmp2, tmp1, &g->coef[i]);
      mpz_sub(&r->coef[rd - (g->degree - i)], &r->coef[rd - (g->degree - i)], tmp2);
      mpz_mod(&r->coef[rd - (g->degree - i)], &r->coef[rd - (g->degree - i)], p);
    }
    while ((rd>0) && (mpz_sgn(&r->coef[rd])==0))
      rd--;
    r->degree = rd;
  }
  mpz_poly_cp(_q, q);
  return rd;
} 
  
/******************************************************/
int mpz_poly_div_qr(mpz_poly _q, mpz_poly _r, mpz_poly x, mpz_poly y, mpz_t p)
/******************************************************/
/* Get q,r so that x = qy + r, deg(r) < deg(y).       */
/* i.e., compute x/y.                                 */
/******************************************************/
{ int    i, j, n;
  static mpz_t  c, lr, ly_inv, tmp1;
  static mpz_poly q, r;
  static int initialized=0;
          
  if (!(initialized)) {
    mpz_init(c); mpz_init(lr); mpz_init(ly_inv); mpz_init(tmp1);
    mpz_poly_init(q); mpz_poly_init(r); 
    initialized=1;
  }

  q->degree = MAX(x->degree - y->degree, 0);
  for (i=0; i<=q->degree; i++)
    mpz_set_ui(&q->coef[i],0);
  mpz_poly_cp(r, x);
                                                                                                 
  mpz_poly_fixDeg(r);
                                                                                                 
  i=q->degree;
  mpz_invert(ly_inv, &y->coef[y->degree], p);
  while ((r->degree >= y->degree) && !((r->degree==0)&&(mpz_sgn(&r->coef[0])==0))) {
    mpz_set(lr, &r->coef[r->degree]);
    mpz_mul(c, ly_inv, lr); mpz_mod(c, c, p);
    n = r->degree - y->degree;
    mpz_set(&q->coef[n], c);
    for (j=r->degree; j>=(r->degree - y->degree); j--) {
      mpz_mul(tmp1, &y->coef[j-n], c); 
      mpz_sub(&r->coef[j], &r->coef[j], tmp1);
      mpz_mod(&r->coef[j], &r->coef[j], p);
//      r->coef[j] = (p+(r->coef[j] - mulmod32(y->coef[j-n], c, p)))%p;
    }
    mpz_poly_fixDeg(r);
  }
  mpz_poly_cp(_q, q); mpz_poly_cp(_r, r);
  return 0;
}



/**************************************************************/
double log_evalPoly(double x, double y, mpz_poly f)
/**************************************************************/
/* Let z = f(x+iy). Return value is an approximation of       */
/* log( |z|).                                                 */
/**************************************************************/
{ double zr, zi, rr, ri, c, tr, ti;
  int i;

  rr = ri =0;
  zr = 1.0;
  zi = 0.0;

  /* Throughout, (zr, zi) = (z,y)^i */
  for (i=0; i<=f->degree; i++) {
    c = mpz_get_d(&(f->coef[i]));
    rr += c*zr;
    ri += c*zi;

    tr = zr*x - zi*y; ti = zr*y + zi*x;
    zr = tr; zi = ti;
  }
  return log(sqrt(rr*rr + ri*ri));
}
    
/**************************************************************/
void mpz_poly_eval(mpz_t res, mpz_poly f, mpz_t x)
{ int i, d=f->degree;
  static int initialized=0;
  static mpz_t tmp, xPow;

  if (!(initialized)) {
    mpz_init(tmp);
    mpz_init(xPow);
    initialized=1;
  }
  mpz_set_ui(xPow, 1);
  mpz_set_ui(res, 0);
  for (i=0; i<=d; i++) {
    mpz_mul(tmp, xPow, &(f->coef[i]));
    mpz_add(res, res, tmp);
    mpz_mul(xPow, xPow, x);
  }
}

/**************************************************************/
void mpz_poly_eval2(mpz_t res, mpz_poly f, mpz_t x, mpz_t y1)
{ int i, j, d=f->degree;
  static int initialized=0;
  static mpz_t tmp, xPow;

  if (!(initialized)) {
    mpz_init(tmp);
    mpz_init(xPow);
    initialized=1;
  }
  mpz_set_ui(xPow, 1);
  mpz_set_ui(res, 0);
  for (i=0; i<=d; i++) {
    mpz_mul(tmp, xPow, &(f->coef[i]));
    for (j=i; j<d; j++) mpz_mul(tmp, tmp, y1);
    mpz_add(res, res, tmp);
    mpz_mul(xPow, xPow, x);
  }
}

/**************************************************************/
void mpz_poly_evalD(mpz_t res, mpz_poly f, mpz_t x)
{ int i, d=f->degree;
  static int initialized=0;
  static mpz_t tmp, xPow;

  if (!(initialized)) {
    mpz_init(tmp);
    mpz_init(xPow);
    initialized=1;
  }
  mpz_set_ui(xPow, 1);
  mpz_set_ui(res, 0);
  for (i=1; i<=d; i++) {
    mpz_mul(tmp, xPow, &(f->coef[i]));
    mpz_mul_ui(tmp, tmp, i);
    mpz_add(res, res, tmp);
    mpz_mul(xPow, xPow, x);
  }
}
  
    
/**************************************************************/
void mpz_poly_HenselPow2(mpz_t res, mpz_poly f, s32 a, s32 p, int n)
/**************************************************************/
/* Given a simple zero of f mod p, lift it to a zero mod p^n, */
/* for n a power of 2. See Bach & Shallit.                    */
/**************************************************************/
{ mpz_t x0, x1, tmp1, tmp2, pn_2;
	
  if (n==1) {
    mpz_set_si(res, a);
    return;
  }
  mpz_init(x0); mpz_init(x1);
  mpz_init(tmp1); mpz_init(tmp2);
  mpz_init(pn_2);
  
  mpz_poly_HenselPow2(x0, f, a, p, n/2);

  mpz_poly_eval(tmp1, f, x0);
  mpz_ui_pow_ui(pn_2, p, n/2);
  
  /* tmp1 <-- f(x0)/p^{n/2}. */
  mpz_div(tmp1, tmp1, pn_2);

  /* tmp1 <-- tmp1 mod p^n.  */
  mpz_ui_pow_ui(tmp2, p, n);
  mpz_mod(tmp1, tmp1, tmp2);

  /* tmp2 <-- f'(x0)^{-1} mod p^{n/2}. */
  mpz_poly_evalD(tmp2, f, x0);
  mpz_mod(tmp2, tmp2, pn_2);
  if (mpz_sgn(tmp2)<0)
    mpz_add(tmp2, tmp2, pn_2);
		
  mpz_invert(tmp2, tmp2, pn_2);

  /* x1 <-- -(tmp1)(tmp2) mod p^{n/2}. */
  mpz_mul(tmp1, tmp1, tmp2);
  mpz_neg(tmp1, tmp1);
  mpz_mod(x1, tmp1, pn_2);
  if (mpz_sgn(x1) < 0)
    mpz_add(x1, x1, pn_2);
	
  

  /* res <-- x0 + x1*p^{n/2}. */
  mpz_mul(res, x1, pn_2);
  mpz_add(res, res, x0);

  mpz_poly_eval(tmp1, f, res);
  mpz_mul(pn_2, pn_2, pn_2);
  mpz_mod(tmp1, tmp1, pn_2);

  mpz_clear(x0); mpz_clear(x1);
  mpz_clear(tmp1); mpz_clear(tmp2);
  mpz_clear(pn_2);
}  
  
/********************************************************/    
int mpz_poly_HenselLift(mpz_t res, mpz_poly f, s32 a, s32 p, int k)
/********************************************************/    
/* Same as above, but for any exponent. We will lift it */
/* as far as it goes (i.e., if it is not a simple zero, */
/* it might only lift to p^m for some small m).         */
/* Return value: The largest m<=k so that 'res' is a    */
/* zero mod p^m.                                        */
/********************************************************/    
{ mpz_t modulus, eval;
  int  n, m;

  mpz_init(modulus);
  mpz_init(eval);
  n=1;
  while (n<k)
    n *= 2;

  mpz_poly_HenselPow2(res, f, a, p, n);
  mpz_ui_pow_ui(modulus, p, k);
  mpz_mod(res, res, modulus);

  mpz_poly_eval(eval, f, res);
  m = k;
  do {
    mpz_mod(eval, eval, modulus);
    if (mpz_sgn(eval)) {
      m--;
      mpz_div_ui(modulus, modulus, p);
    }
  } while (mpz_sgn(eval));
  mpz_mod(res, res, modulus);
  
  mpz_clear(modulus); mpz_clear(eval);
  return m;
}

        
/******************************************************************/
int mpz_poly_pow_mod(mpz_poly res, mpz_poly f, mpz_t n, mpz_poly mod)
/******************************************************************/
{ mpz_poly g, h;
  int      i;
  mpz_t remain;

  mpz_poly_init(g);
  mpz_poly_init(h);
  mpz_init(remain);

  g->degree = f->degree;
  for (i=0; i<=g->degree; i++)
    mpz_set(&(g->coef[i]), &(f->coef[i]));
  
  h->degree=0;  mpz_set_ui(&(h->coef[0]), 1);

  mpz_set(remain, n);
  while (mpz_sgn(remain)) {
    if (mpz_odd_p(remain)) {
      mpz_poly_mul(h, h, g);
      mpz_poly_mod(h, h, mod);
    }
    mpz_div_2exp(remain, remain, 1);
    mpz_poly_mul(g, g, g);
    mpz_poly_mod(g, g, mod);
  }

  for (i=h->degree; i>=0; i--)
    mpz_set(&(res->coef[i]), &(h->coef[i]));
  res->degree = h->degree;

  mpz_poly_clear(g);
  mpz_poly_clear(h);
  mpz_clear(remain);
  
  return res->degree;
}

/******************************************************************/
int mpz_poly_pow_mod_pp(mpz_poly res, mpz_poly f, mpz_t n, mpz_poly mod, mpz_t p)
/******************************************************************/
{ mpz_poly g, h;
  int      i;
  mpz_t remain;

  mpz_poly_init(g);
  mpz_poly_init(h);
  mpz_init(remain);

  g->degree = f->degree;
  for (i=0; i<=g->degree; i++)
    mpz_set(&(g->coef[i]), &(f->coef[i]));
  
  h->degree=0;  mpz_set_ui(&(h->coef[0]), 1);

  mpz_set(remain, n);
  while (mpz_sgn(remain)) {
    if (mpz_odd_p(remain)) {
      mpz_poly_mulmod_pp(h, h, g, mod, p);
    }
    mpz_div_2exp(remain, remain, 1);
    mpz_poly_mulmod_pp(g, g, g, mod, p);
  }

  for (i=h->degree; i>=0; i--)
    mpz_set(&(res->coef[i]), &(h->coef[i]));
  res->degree = h->degree;

  mpz_poly_clear(g);
  mpz_poly_clear(h);
  mpz_clear(remain);
  
  return res->degree;
}


/******************************************************/    
int mpz_poly_mod(mpz_poly r, mpz_poly a, mpz_poly b)
/******************************************************/
/* r <-- a mod b. b is assumed monic!                 */
/******************************************************/
{ mpz_poly R;
  mpz_t    c, tmp;
  int      i, d;

  /* Take care of some special cases. */
  if (b->degree == 0) {
    r->degree=0;
    mpz_set_ui(&r->coef[0], 0);
    return 0;
  }  
  mpz_poly_init(R);
  mpz_init(c);
  mpz_init(tmp);
  
  R->degree = a->degree;
  for (i=0; i<=a->degree; i++)
    mpz_set(&R->coef[i], &a->coef[i]);

  while (R->degree >= b->degree) {
    mpz_set(c, &R->coef[R->degree]);
    /* Now do R <-- R - c*b. */
    d = R->degree;
    for (i=b->degree; i>=0; i--) {
      mpz_mul(tmp, c, &b->coef[i]); 
      mpz_sub(&R->coef[d-(b->degree-i)], &R->coef[d-(b->degree-i)], tmp);
    }
    while ((R->degree >0) && (mpz_sgn(&R->coef[R->degree])==0))
      R->degree -= 1;
  }
  r->degree = R->degree;
  for (i=0; i<=r->degree; i++)
    mpz_set(&r->coef[i], &R->coef[i]);

  mpz_poly_clear(R);
  mpz_clear(c);
  mpz_clear(tmp);

  return 0;
}



/******************************************************/    
int mpz_poly_mod_pp(mpz_poly r, mpz_poly a, mpz_poly b, mpz_t p)
/******************************************************/    
{ mpz_poly R;
  mpz_t    c, lci, tmp;
  int      i, d;

  /* Take care of some special cases. */
  if (b->degree == 0) {
    r->degree=0;
    mpz_set_ui(&(r->coef[0]), 0);
    return 0;
  }  
	  
  
  mpz_poly_init(R);
  mpz_init(c);
  mpz_init(lci);
  mpz_init(tmp);

  mpz_invert(lci, &(b->coef[b->degree]), p);
  
  
  R->degree = a->degree;
  for (i=0; i<=a->degree; i++)
    mpz_set(&(R->coef[i]), &(a->coef[i]));

  while (R->degree >= b->degree) {
    mpz_mul(c, &(R->coef[R->degree]), lci);
    mpz_mod(c, c, p);
    /* Now do R <-- R - c*b. */
    d = R->degree;
    for (i=b->degree; i>=0; i--) {
      mpz_mul(tmp, c, &(b->coef[i])); mpz_mod(tmp, tmp, p);
      
      mpz_sub(&(R->coef[d-(b->degree-i)]), &(R->coef[d-(b->degree-i)]), tmp);
      mpz_mod(&(R->coef[d-(b->degree-i)]), &(R->coef[d-(b->degree-i)]), p);
    }
    while ((R->degree >0) && (mpz_sgn(&(R->coef[R->degree]))==0))
      R->degree -= 1;
  }
  r->degree = R->degree;
  for (i=0; i<=r->degree; i++)
    mpz_set(&(r->coef[i]), &(R->coef[i]));

  mpz_poly_clear(R);
  mpz_clear(c);
  mpz_clear(lci);
  mpz_clear(tmp);

  return 0;
}

	  

/******************************************************/    
int mpz_poly_gcd(mpz_poly g, mpz_poly h, mpz_t p)
/******************************************************/    
/* g <-- gcd(g,h)                                     */
/* g and h will both be modified.                     */
/******************************************************/    
{ static mpz_poly r;
  static mpz_t    tmp;
  static int      initialized=0;
  int i;

  if (!(initialized)) {
    mpz_poly_init(r);
    mpz_init(tmp);
    initialized=1;
  }

  /*************************************/
  /* Make sure the degrees are correct */
  /* or horrible things will happen.   */
  /*************************************/
  for (i=0; i<=h->degree; i++)
    mpz_mod(&(h->coef[i]), &(h->coef[i]), p);
  for (i=0; i<=g->degree; i++)
    mpz_mod(&(g->coef[i]), &(g->coef[i]), p);

  while ((h->degree >0) && (mpz_sgn(&(h->coef[h->degree]))==0))
    h->degree -= 1;
  while ((g->degree >0) && (mpz_sgn(&(g->coef[g->degree]))==0))
    g->degree -= 1;

    
  while ((h->degree>0) || (mpz_sgn(&(h->coef[h->degree])))) {
    mpz_poly_mod_pp(r, g, h, p);
    for (i=h->degree; i>=0; i--) 
      mpz_set(&(g->coef[i]), &(h->coef[i]));
    g->degree = h->degree;
    for (i=r->degree; i>=0; i--) 
      mpz_set(&(h->coef[i]), &(r->coef[i]));
    h->degree = r->degree;
    
  }
  while ((g->degree >0) && (mpz_sgn(&(g->coef[g->degree]))==0))
    g->degree -= 1;
  
  mpz_invert(tmp, &g->coef[g->degree], p);
  for (i=g->degree; i>=0; i--) {
    mpz_mul(&g->coef[i], &g->coef[i], tmp);
    mpz_mod(&g->coef[i], &g->coef[i], p);
  }
  return 1;
}


/************************************************/
INLINE int mpz_poly_modp(poly_t res, mpz_poly q, s32 p)
/************************************************/
/* Comute the reduction: res <-- q(x) mod p.    */
/************************************************/
{ int i;

  res->degree = q->degree;
  for (i=0; i<=res->degree; i++) {
    res->coef[i] = mpz_fdiv_ui(&(q->coef[i]), p);
    if (res->coef[i] < 0)
      res->coef[i] += p;
  }
  poly_fixDeg(res);
  
  return 0;
}



/**************************************************************************/
int mpz_polyCoprime(mpz_poly f, mpz_poly g, mpz_t p)
/**************************************************************************/
/* Return 1 if gcd(f,g)=1, where f and g are polys mod p.                 */
/**************************************************************************/
{ mpz_poly tmp1, tmp2;
  int i, retVal;

  mpz_poly_init(tmp1);
  mpz_poly_init(tmp2);
  
  tmp1->degree = f->degree;
  for (i=0; i<=f->degree; i++)
    mpz_set(&(tmp1->coef[i]), &(f->coef[i]));
  
  tmp2->degree = g->degree;
  for (i=0; i<=g->degree; i++)
    mpz_set(&(tmp2->coef[i]), &(g->coef[i]));
 
  mpz_poly_gcd(tmp1, tmp2, p);

  retVal = 1;
  if (tmp1->degree > 0)
    retVal = 0;

  mpz_poly_clear(tmp1);
  mpz_poly_clear(tmp2);
  return retVal;
}


/**************************************************************************/
INLINE s32 mpz_poly_evalModp(mpz_poly h, s32 p, s32 x)
{ int  i, d=h->degree;
  s32 xPow=1, res, tmp;

  res = 0;
  for (i=0; i<=d; i++) {
    tmp = mpz_fdiv_ui(&(h->coef[i]), p);
    tmp = (tmp < 0) ? (tmp+p) : tmp;
    res = (res + mulmod32(tmp, xPow, p))%p;
    xPow = mulmod32(xPow, x, p);
  }
  return res;
}

/**************************************************************************/
int mpz_poly_cmp(mpz_poly A, mpz_poly B)
/**************************************************************************/
{ int i, dA=A->degree, dB=B->degree;
  int d = MIN(dA, dB), r;

  i=0;
  while (i<=d) {
    r = mpz_cmp(&A->coef[i], &B->coef[i]);
    if (r)
      return r;
    i++;
  }
  if (dA < dB)
    return -1;
  else if (dA > dB)
    return 1;
  return 0;
}

/****************************************************************************/
void mpz_poly_mul(mpz_poly res, mpz_poly op1, mpz_poly op2)
{ int i, j, d1=op1->degree, d2=op2->degree, d;
  static mpz_poly Tmp;
  static mpz_t    tmp1;
  static int initialized=0;
  
  if (!(initialized)) {
    mpz_poly_init(Tmp);
    mpz_init(tmp1);
    initialized=1;
  }
  
  d = d1+d2;
  for (i=d; i>=0; i--)
    mpz_set_ui(&Tmp->coef[i], 0);

  for (i=0; i<=d1; i++) {
    for (j=0; j<=d2; j++) {
      mpz_mul(tmp1, &op1->coef[i], &op2->coef[j]);
      mpz_add(&Tmp->coef[i+j], &Tmp->coef[i+j], tmp1);
    }
  }
  Tmp->degree = d;
  mpz_poly_fixDeg(Tmp);
  mpz_poly_cp(res, Tmp);
}

/****************************************************************************/
void mpz_poly_add(mpz_poly res, mpz_poly op1, mpz_poly op2)
{ int i, d;
  static mpz_poly Tmp;
  static int initialized=0;
  
  if (!(initialized)) {
    mpz_poly_init(Tmp);
    initialized=1;
  }
  Tmp->degree = MAX(op1->degree, op2->degree);
  d = MIN(op1->degree, op2->degree);
  for (i=0; i<=d; i++)
    mpz_add(&Tmp->coef[i], &op1->coef[i], &op2->coef[i]);
  while (i<=op1->degree) {
    mpz_set(&Tmp->coef[i], &op1->coef[i]);
    i++;
  }
  while (i<=op2->degree) {
    mpz_set(&Tmp->coef[i], &op2->coef[i]);
    i++;
  }
  mpz_poly_fixDeg(Tmp);
  mpz_poly_cp(res, Tmp);
}
  

/****************************************************************************/
void mpz_poly_mulmod(mpz_poly res, mpz_poly op1, mpz_poly op2, mpz_poly mod)
{ mpz_poly t;

  mpz_poly_init(t);
  mpz_poly_mul(t, op1, op2);
  mpz_poly_mod(t, t, mod);
  mpz_poly_cp(res, t);
  mpz_poly_clear(t);
}
	
/****************************************************************************/
void mpz_poly_mulmod_p(mpz_poly res, mpz_poly op1, mpz_poly op2, mpz_t p)
{ int i, j, d1=op1->degree, d2=op2->degree, d;
  static mpz_poly Tmp;
  static mpz_t    tmp1;
  static int initialized=0;
  
  if (!(initialized)) {
    mpz_poly_init(Tmp);
    mpz_init(tmp1);
    initialized=1;
  }
  
  d = d1+d2;
  for (i=d; i>=0; i--)
    mpz_set_ui(&Tmp->coef[i], 0);

  for (i=0; i<=d1; i++) {
    for (j=0; j<=d2; j++) {
      mpz_mul(tmp1, &op1->coef[i], &op2->coef[j]);
      mpz_add(&Tmp->coef[i+j], &Tmp->coef[i+j], tmp1);
      mpz_mod(&Tmp->coef[i+j], &Tmp->coef[i+j], p);
    }
  }
  Tmp->degree = d;
  mpz_poly_fixDeg(Tmp);
  mpz_poly_cp(res, Tmp);
}

/****************************************************************************/
void mpz_poly_mulmod_pp(mpz_poly res, mpz_poly op1, mpz_poly op2, mpz_poly mod, mpz_t p)
/* 'mod' is assumed to be monic!. */
{ static mpz_poly xPow, tRes;
  static mpz_t    t1, t2;
  static int initialized=0;
  int    i, j, k, d=mod->degree, d1=op1->degree, d2=op2->degree;

  if (!(initialized)) {
    mpz_poly_init(xPow); mpz_poly_init(tRes);
    mpz_init(t1);        mpz_init(t2);
    initialized=1;
  }

  for (i=0; i<d; i++) {
    mpz_set_ui(&(xPow->coef[i]), 0);
    mpz_set_ui(&(tRes->coef[i]), 0);
  }
  mpz_set_ui(&(xPow->coef[0]), 1);
  
  for (i=0; i<=(d1+d2); i++) {
    /* Find the coefficient of x^i in the product op1*op2. */
    mpz_set_ui(t1, 0);
    j = MAX(0, i-d2);
    j = MIN(j, d1); /* Is this necessary?? */
    while ((j<=d1) && ((i-j) >=0)) {
      mpz_mul(t2, &(op1->coef[j]), &(op2->coef[i-j]));
      mpz_add(t1, t1, t2);
      mpz_mod(t1, t1, p);
      j++;
    }

    /* tRes <-- tRes + t1*xPow. */
    for (k=0; k<d; k++) {
      mpz_mul(t2, t1, &(xPow->coef[k]));
      mpz_add(&(tRes->coef[k]), &(tRes->coef[k]), t2);
      mpz_mod(&(tRes->coef[k]), &(tRes->coef[k]), p);
    }

    /* xPow <-- xPow*x (modulo 'mod'). */
    mpz_set(t2, &(xPow->coef[d-1]));
    for (k=d-1; k>0; k--)
      mpz_set(&(xPow->coef[k]), &(xPow->coef[k-1]));
    mpz_set_ui(&(xPow->coef[0]), 0);
    for (k=0; k<d; k++) {
      mpz_mul(t1, t2, &(mod->coef[k]));
      mpz_sub(&(xPow->coef[k]), &(xPow->coef[k]), t1);
      mpz_mod(&(xPow->coef[k]), &(xPow->coef[k]), p);
    }
  }
  for (i=0; i<d; i++)
    mpz_set(&(res->coef[i]), &(tRes->coef[i]));
  i=d-1;
  while ((i>0) && (mpz_sgn(&(res->coef[i]))==0))
    i--;
  res->degree = i;
}

/**************************************************************************/
int mpz_poly_irreduciblelike_modp(mpz_poly _f, mpz_t p)
/**************************************************************************/
/* Input: A polynomial 'f' of degree 'deg_f', and a prime integer 'p'.    */
/* Return value: 1, if 'f' has no roots mod p and x^{p^deg(f)}==x mod f.  */
/* 0, otherwise.                                                          */
/**************************************************************************/
{ mpz_poly xPow, tmp, f;
  mpz_t    c;
  int      i, retVal=1;

  mpz_poly_init(xPow); mpz_poly_init(tmp); mpz_poly_init(f);
  mpz_init(c);

  mpz_mod(c, &_f->coef[_f->degree], p);
  if (mpz_sgn(c)==0) {
    retVal = 0;
    goto MP_IR_DONE;
  }

  
  f->degree = _f->degree;
  mpz_invert(c, &_f->coef[_f->degree], p);
  
  /* Make 'f' monic. */
  for (i=0; i<=f->degree; i++) {
    mpz_mul(&f->coef[i], &_f->coef[i], c);
    mpz_mod(&f->coef[i], &f->coef[i], p);
  }
  
  /* First, compute x^p mod f. */
  for (i=0; i<= f->degree; i++)
    mpz_set_ui(&xPow->coef[i], 0);
  mpz_set_ui(&xPow->coef[1], 1);
  xPow->degree=1;
  

  mpz_poly_pow_mod_pp(xPow, xPow, p, f, p);
  
  /* Verify that (x^p - x, f) = 1. */
  
  
  tmp->degree = MAX(xPow->degree, 1);
  for (i=0; i<= xPow->degree; i++)
    mpz_set(&tmp->coef[i], &xPow->coef[i]);
  mpz_sub_ui(&tmp->coef[1], &tmp->coef[1], 1);
  if (mpz_sgn(&tmp->coef[1]) < 0)
    mpz_add(&tmp->coef[1], &tmp->coef[1], p);
  
  while ((tmp->degree >0) && (mpz_sgn(&tmp->coef[tmp->degree])==0))
    tmp->degree -= 1;
  
  if (!(mpz_polyCoprime(tmp, f, p))) {
    retVal=0;
    goto MP_IR_DONE;
  }

  /* Now, compute x^{p^{deg_f}} mod f. */
  for (i=2; i<=f->degree; i++)
    mpz_poly_pow_mod_pp(xPow, xPow, p, f, p);
  
  /* Check that x^{p^{deg_f}} == x (mod f). */
  if (xPow->degree != 1) {
    retVal = 0;
  }
  else {
    if (mpz_cmp_ui(&xPow->coef[1], 1))
      retVal = 0;
    if (mpz_sgn(&xPow->coef[0]))
      retVal = 0;
  }
  
MP_IR_DONE:
  mpz_poly_clear(xPow);
  mpz_poly_clear(tmp);
  mpz_poly_clear(f);
  mpz_clear(c);
  return retVal;
}
  
/**************************************************************************/
int mpz_poly_irreducible_modp(mpz_poly _f, mpz_t p)
/**************************************************************************/
/* Input: A polynomial 'f' of degree 'deg_f', and a prime integer 'p'.    */
/* Return value: 1, if 'f' is irreducible mod p. 0 otherwise.             */
/**************************************************************************/
{ static mpz_poly factors[MAXPOLYDEGREE];
  static int initialized=0;
  int i, n, e[MAXPOLYDEGREE];

  if (!initialized) {
    for (i=0; i<MAXPOLYDEGREE; i++)
      mpz_poly_init(factors[i]);
    initialized=1;
  }
  n=mpz_poly_fact(factors, e, _f, p);
  if ((n==1) && (e[0]==1))
    return 1;
  return 0;
}




/**************************************************************************/
int mpz_poly_sqf_fact(mpz_poly *Ai, int *exps, mpz_poly A, mpz_t p)
/**************************************************************************/
/* Compute the squarefree factorization of A mod p via Cohen alg. 3.4.2.  */
/* The Ai and i must already be allocated and initialized.                */
/* Note: We assume that any exponents are small enough to fit in an 'int'.*/
/* Return value: Number of factors.                                       */
/**************************************************************************/
{ int i, e, k, numFactors=0, p_;
  mpz_poly T, T0, V, W, Tmp;
  
  
  mpz_poly_init(T);
  mpz_poly_init(T0);
  mpz_poly_init(V);
  mpz_poly_init(W);
  mpz_poly_init(Tmp);

  /* This is only used if p is small anyway, so it's ok. */
  p_ = (int)mpz_get_ui(p);
  
  /* Step 1. */
  e=1;
  mpz_poly_cp(T0, A);
 
 
  /* Step 2. */
ALG_3_4_2_STEP2:
  if (T0->degree == 0) {
    mpz_poly_clear(T); mpz_poly_clear(T0); mpz_poly_clear(V);
    mpz_poly_clear(W); mpz_poly_clear(Tmp);
    return numFactors;
  }
  mpz_poly_diff(Tmp, T0);
  mpz_poly_cp(T, T0);
  mpz_poly_gcd(T, Tmp, p); /* T <-- (T0, T0') mod p. */


  mpz_poly_div(V, T0, T, p);
  k = 0;

  /* Step 3. */
ALG_3_4_2_STEP3:
  
  if (V->degree == 0) {
    /* Then p must divide the exponents of all nonzero terms of T. */
    /* By assumption, it must be the case that p is small.         */
    if (T->degree % p_)
      fprintf(stderr, "mpz_poly_sqf_fact() bad error at step 3!\n");
    T0->degree = T->degree/p_;
    for (i=T0->degree; i>=0; i--)
      mpz_set(&T0->coef[i], &T->coef[i*p_]);
    e *= p_;
    goto ALG_3_4_2_STEP2;
  }

  /* Step 4. */
  k++;
  if ((mpz_cmp_ui(p, k) <= 0) && (k%p_ == 0)) {
    mpz_poly_div(T, T, V, p);
    k++;
  }

  /* Step 5. */
  mpz_poly_cp(W, T);
  mpz_poly_cp(Tmp, V);
  mpz_poly_gcd(W, Tmp, p); /* W <-- (T, V) mod p. */
  
  mpz_poly_div(Tmp, V, W, p);
  mpz_poly_cp(V, W);
  mpz_poly_div(T, T, V, p);
  
  if (Tmp->degree > 0) {
    mpz_poly_cp(Ai[numFactors], Tmp);
    exps[numFactors] = e*k;
    numFactors++;
  }
  
  goto ALG_3_4_2_STEP3;
}  

/**************************************************************************/
int mpz_poly_Berlekamp(mpz_poly *E, mpz_poly A, mpz_t p)
#define SMALL_P_BOUND 13
/**************************************************************************/
/* Input: A squarefree monic polynomial A(x)\in F_p[x].                   */
/* Output: The factorization A = E[0]*...E[k-1].                          */
/* Return value: k>0 on success, negative on error.                       */
/**************************************************************************/
/* This is Cohen Alg. 3.4.10 (when p is small), 3.4.11 (when p is larger).*/
/**************************************************************************/
{ int      n=A->degree, i, j, k, l, r, smallP, a, fSize, found;
  s32     p_=1, s;
  mpz_poly xPow, x_p, T, B, D, Tmp, F[SMALL_P_BOUND];
  mpz_mat_t  Q, V;
  mpz_t    tmp1;
  
  smallP = (mpz_cmp_ui(p, SMALL_P_BOUND) <= 0) ? 1 : 0;
  if (smallP) {
    p_ = mpz_get_ui(p);
    for (i=0; i<p_; i++)
      mpz_poly_init(F[i]);
  }
  
  mpz_poly_init(xPow); mpz_poly_init(x_p); mpz_poly_init(Tmp);
  mpz_poly_init(T); mpz_poly_init(B); mpz_poly_init(D);
  mpz_mat_init(&Q); mpz_mat_init(&V);
  mpz_init(tmp1);
  
  /** Step 1. **/

  /* xPow <-- 1. */
  mpz_set_ui(&xPow->coef[0], 1); xPow->degree = 0;
    /* x_p <-- x^p mod <A, p>. */
  mpz_set_ui(&x_p->coef[0], 0); mpz_set_ui(&x_p->coef[1], 1); x_p->degree = 1;
  mpz_poly_pow_mod_pp(x_p, x_p, p, A, p);

  Q.rows = n;
  Q.cols = n;
  for (k=0; k<n; k++) {
    /* The k-th column of Q <-- xPow. */
    for (i=0; i<=xPow->degree; i++)
      mpz_set(&Q.entry[i][k], &xPow->coef[i]);
    for ( ; i<n; i++)
      mpz_set_ui(&Q.entry[i][k], 0);

    mpz_poly_mulmod_pp(xPow, xPow, x_p, A, p);
  }


  /** Step 2 **/

  /* First Q <-- Q-I. */
  for (i=0; i<n; i++) {
    mpz_sub_ui(&Q.entry[i][i], &Q.entry[i][i], 1);
    if (mpz_sgn(&Q.entry[i][i]) < 0)
      mpz_add(&Q.entry[i][i], &Q.entry[i][i], p);
  }
  r = mpz_mat_getKernel_modp(&V, &Q, p);
  mpz_poly_cp(E[0], A);
  k=1;  
  j=1;

  
  /** Step 3 **/
ALG_3_4_1X_STEP3:
  if (k==r)
    goto ALG_3_4_1X_DONE;
  if (smallP) {
    j++;
    T->degree=n-1;
    for (i=0; i<n; i++)
      mpz_set(&T->coef[i], &V.entry[i][j-1]);
    mpz_poly_fixDeg(T);
  } else {
    /* T <-- ai*Vi, for random ai. */
    T->degree = n-1;
    for (i=0; i<n; i++)
      mpz_set_ui(&T->coef[i], 0);
    for (i=0; i<r; i++) {
      a = rand();
      for (l=0; l<n; l++) {
        mpz_mul_ui(tmp1, &V.entry[l][i], a);
	mpz_add(&T->coef[l], &T->coef[l], tmp1);
	mpz_mod(&T->coef[l], &T->coef[l], p);
      }
    }
    mpz_poly_fixDeg(T);
  }
      
  /** Step 4 **/
  if (smallP) {
    for (i=0; i<k; i++) {
      mpz_poly_cp(B, E[i]);
      fSize=0;
      if (B->degree > 1) {
        for (s=0; s<p_; s++) {
          mpz_poly_cp(D, B);
  	  mpz_poly_cp(Tmp, T);
          mpz_sub_ui(&Tmp->coef[0], &Tmp->coef[0], s);
          mpz_mod(&Tmp->coef[0], &Tmp->coef[0], p);
          mpz_poly_gcd(D, Tmp, p);
          if (D->degree >= 1) {
            for (l=found=0; l<fSize; l++)
              if (mpz_poly_cmp(D, F[l])==0)
                found = 1;
            if (!(found)) {
              mpz_poly_cp(F[fSize], D);
              fSize++;
            }
          }
        }
      }
      if (fSize) {
        mpz_poly_cp(E[i], F[0]);
	for (l=1; l<fSize; l++)
          mpz_poly_cp(E[k+l-1], F[l]);
	k = k + fSize - 1;
	if (k==r)
          goto ALG_3_4_1X_DONE;
      }
    }
    goto ALG_3_4_1X_STEP3;
  } else {
    mpz_sub_ui(tmp1, p, 1); mpz_div_2exp(tmp1, tmp1, 1); /* tmp1 <-- (p-1)/2. */
    
    for (i=0; i<k; i++) {
	    
      mpz_poly_cp(B, E[i]);
      if (B->degree > 1) {
        mpz_poly_pow_mod_pp(Tmp, T, tmp1, B, p);
	mpz_sub_ui(&Tmp->coef[0], &Tmp->coef[0], 1);
	if (mpz_sgn(&Tmp->coef[0]) < 0)
          mpz_add(&Tmp->coef[0], &Tmp->coef[0], p);
        mpz_poly_cp(D, B);
        mpz_poly_gcd(D, Tmp, p);
        if ((D->degree > 0) && (D->degree < B->degree)) {
          mpz_poly_cp(E[i], D);
	  mpz_poly_div(E[k], B, D, p);
	  k++;
	  if (k==r)
            goto ALG_3_4_1X_DONE;
	}
      }
    }
    goto ALG_3_4_1X_STEP3;
  }    
	  

ALG_3_4_1X_DONE:
  mpz_poly_clear(xPow); mpz_poly_clear(x_p); mpz_poly_clear(Tmp);
  mpz_poly_clear(T); mpz_poly_clear(B); mpz_poly_clear(D);
  mpz_mat_clear(&Q); mpz_mat_clear(&V);
  mpz_init(tmp1);
  if (smallP)
    for (i=0; i<p_; i++)
      mpz_poly_clear(F[i]);
  return k;
}
  
  
  

	
/**************************************************************************/
int mpz_poly_fact(mpz_poly *Pi, int *exps, mpz_poly _A, mpz_t p)
/**************************************************************************/
/* Input: A nonzero polynomial A\in F_p[x].                               */
/* Output: Irredicuble polys Pi[0],...,Pi[k-1] and exps[0],...,exps[k-1]  */
/*     s.t. A = (Pi[0]^exps[0])...(Pi[k-1]^exps[k-1]) is the factorization*/
/*     of A mod p.                                                        */
/* Return value: k>0 on success, negative on error.                       */
/**************************************************************************/
{ int      d=_A->degree, i, j, k;
  mpz_t    lc, lci;
  mpz_poly A, *Sqf;
  int      *sqfExps, sqfSize, s;

  /* Special case: */
  if (_A->degree <= 1) {
    mpz_poly_cp(Pi[0], _A);
    return 1;
  }
  mpz_poly_init(A);
  mpz_poly_cp(A, _A);

  Sqf = malloc(d * sizeof(mpz_poly));
  sqfExps = malloc(d * sizeof(int));
  for (i=0; i<d; i++)
    mpz_poly_init(Sqf[i]);
  mpz_poly_init(A);
  mpz_init(lc);
  mpz_init(lci);

  mpz_poly_cp(A, _A);
  /* First, make A monic. */
  mpz_set(lc, &A->coef[A->degree]);
  if (mpz_cmp_ui(lc, 1)) {
    mpz_invert(lci, lc, p);
    for (i=A->degree; i>=0; i--) {
      mpz_mul(&A->coef[i], &A->coef[i], lci);
      mpz_mod(&A->coef[i], &A->coef[i], p);
    }
  }

  
  /* Get the squarefree factorization. */
  sqfSize = mpz_poly_sqf_fact(Sqf, sqfExps, A, p);
  k=0;

  /* Factor each squarefree part. */
  for (i=0; i<sqfSize; i++) {
    /* Now Sqf[i]^exps[i] divides A and Sqf[i] is squarefree. */
    /* So factor it and adjust the exponents.                 */
    s = mpz_poly_Berlekamp(&Pi[k], Sqf[i], p);
    for (j=0; j<s; j++)
      exps[k+j] = sqfExps[i];
    k += s;
  }
  if (mpz_cmp_ui(lc, 1)) {
    mpz_set(&(Pi[k]->coef[0]), lc);
    Pi[k]->degree = 0;
    exps[k] = 1;
    k++;
  }
  for (i=0; i<d; i++)
    mpz_poly_clear(Sqf[i]);
  mpz_poly_clear(A);
  mpz_clear(lc);
  mpz_clear(lci);
  free(Sqf);
  free(sqfExps);

  return k;
}


double getC_(mpz_poly f, mpz_poly a);
/****************************************************************************/
void getSqrtP(mpz_t p, mpz_poly a, mpz_poly f)
/* Choose a prime for the Zalpha_sqrt function. */
{ double C;
  int    attempt=0, e[MAXPOLYDEGREE], i, n, sqf;
  mpz_poly factors[MAXPOLYDEGREE];
	
  C = getC_(f, a);
  mpz_set_d(p, C);
  do {
    mpz_nextprime(p, p);
    attempt++;
  } while (!(mpz_poly_irreducible_modp(f, p)) && (attempt < 500));
  if (attempt == 500) {
    msgLog(NULL, "Warning: getSqrtP() - f appears to split modulo all primes!");
    msgLog(NULL, "grabbing a prime for which f remains squarefree.");
    for (i=0; i<MAXPOLYDEGREE; i++)
      mpz_poly_init(factors[i]);
    do {
      mpz_nextprime(p, p);
      n = mpz_poly_fact(factors, e, f, p);
      sqf=1;
      for (i=0; i<n; i++)
        if (e[i]>1) sqf=0;
    } while (!sqf);
  }
}

/****************************************************************************/
double getC_(mpz_poly f, mpz_poly a)
/* Compute s*||M||^-1, where M is the matrix \Theta in Cohen, section 4.2.4 */
/* and s = sqrt( sum |a(\alpha_i)|^2).                                      */
/* For now, just return some safe upper bound.                              */
{ int           i, j, d=f->degree;
  double        s, sx, sy, x, y, c, tx, ty;
  nfs_complex_t *Z;

  Z = malloc(d * sizeof(nfs_complex_t));
  for (i=0; i<d; i++) {
    mpf_init2(Z[i].mpr, 128);
    mpf_init2(Z[i].mpi, 128);
  }
  mpz_poly_getComplexZeros(Z, f);
  s = 0.0;
  sx = 0.0; sy = 0.0;
  
  for (i=0; i<d; i++) {
    x = 1.0; y = 0.0;
    for (j=0; j<=a->degree; j++) {
      c = mpz_get_d(&(a->coef[i]));
      sx += c*x;
      sy += c*y;
      tx = x*Z[i].r - y*Z[i].i;
      ty = x*Z[i].i + y*Z[i].r;
      x = tx; y = ty;
    }
    s += sx*sx + sy*sy;
  }
  //s = sqrt(s);
  free(Z);
  return MAX(10000.0*s, 10e10);
}


/******************************************************/    
int mpz_poly_inv(mpz_poly inv, mpz_poly f, mpz_poly m, mpz_t p)
/******************************************************/    
/* s <-- f^{-1} mod <p, m(x)>.                        */
/******************************************************/    
{ int    i; 
  static mpz_poly r, g, h, t, s1, s2, t1, t2, q, s;
  static mpz_t    c;
  static int initialized=0;

  if (!(initialized)) {
    mpz_poly_init(r); mpz_poly_init(g); mpz_poly_init(h);  
    mpz_poly_init(t); mpz_poly_init(s1); mpz_poly_init(s2);  
    mpz_poly_init(t1); mpz_poly_init(t2); mpz_poly_init(q);  
    mpz_poly_init(s); mpz_init(c);
    initialized=1;
  }


  mpz_poly_cp(g, f);
  mpz_poly_cp(h, m);

  /*************************************/
  /* Make sure the degrees are correct */
  /* or horrible things will happen.   */
  /*************************************/
  mpz_poly_fixDeg(h);
  mpz_poly_fixDeg(g);

  if ((h->degree==0) && (mpz_sgn(&h->coef[0])==0)) {
    mpz_set_ui(&inv->coef[0], 1);
    inv->degree = 0;
    return 0;
  }
  
  mpz_set_ui(&s2->coef[0],1); s2->degree=0;
  mpz_set_ui(&s1->coef[0],0); s1->degree=0;
  mpz_set_ui(&t2->coef[0],0); t2->degree=0;
  mpz_set_ui(&t1->coef[0],1); t1->degree=0;
  

  while ((h->degree>0) || (mpz_sgn(&h->coef[0]))) {
//printf("  ................\n");
//printf("p = "); mpz_out_str(stdout, 10, p); printf("\n");
//mpz_poly_print(stdout, "g = ", g);
//mpz_poly_print(stdout, "h = ", h);

    mpz_poly_div_qr(q, r, g, h, p);

//mpz_poly_print(stdout, "q = ", q);
//mpz_poly_print(stdout, "r = ", r);


//    mpz_poly_psuedoDiv(q, r, g, h);
//    for (i=0; i<=q->degree; i++)
//      mpz_mod(&q->coef[i], &q->coef[i], p);
//    for (i=0; i<=r->degree; i++)
//      mpz_mod(&r->coef[i], &r->coef[i], p);

    /* s <-- s2 - q*s1. */ 
//printf("\n");
//mpz_poly_print(stdout, "q = ", q);
//mpz_poly_print(stdout, "s1 = ", s1);
    mpz_poly_mulmod_p(s, q, s1, p);   
    for (i=0; i<=s->degree; i++) {
      mpz_neg(&s->coef[i], &s->coef[i]);
      mpz_mod(&s->coef[i], &s->coef[i], p);
    }
    while (i <= s2->degree) {
      mpz_set_ui(&s->coef[i], 0);
      s->degree += 1;
      i++;
    }
    for (i=0; i<=s2->degree; i++) {
      mpz_add(&s->coef[i], &s->coef[i], &s2->coef[i]);
      mpz_mod(&s->coef[i], &s->coef[i], p);
    }
    mpz_poly_fixDeg(s);

    /* t <-- t2 - q*t1 */
    mpz_poly_mulmod_p(t, q, t1, p);   
    for (i=0; i<=t->degree; i++) {
      mpz_neg(&t->coef[i], &t->coef[i]);
      mpz_mod(&t->coef[i], &t->coef[i], p);
    }
    while (i <= t2->degree) {
      mpz_set_ui(&t->coef[i], 0);
      t->degree += 1;
      i++;
    }
    for (i=0; i<=t2->degree; i++) {
      mpz_add(&t->coef[i], &t->coef[i], &t2->coef[i]);
      mpz_mod(&t->coef[i], &t->coef[i], p);
    }
    mpz_poly_fixDeg(t);
    
    /* g <-- h, h <-- r */
    mpz_poly_cp(g, h);
    mpz_poly_cp(h, r);
    /* s2 <-- s1, s1 <-- s */
    mpz_poly_cp(s2, s1);
    mpz_poly_cp(s1, s);
    /* t2 <-- t1, t1 <-- t */
    mpz_poly_cp(t2, t1);
    mpz_poly_cp(t1, t);
  }
//printf("Main loop done.\n");
  if (g->degree >0) {
    fprintf(stderr, "mpz_poly_inv(): Error - polynomial was not invertible!\n");
    exit(-1);
  }
  mpz_invert(c, &g->coef[0], p);

  inv->degree = s2->degree;
  for (i=0; i<=s2->degree; i++) {
    mpz_mul(&inv->coef[i], c, &s2->coef[i]);
    mpz_mod(&inv->coef[i], &inv->coef[i], p);
  }
  return 0;
}

/****************************************************************************/
int Zalpha_sqrt(mpz_poly res, mpz_poly _a, mpz_poly _f, mpz_t N, mpz_t _p)
/****************************************************************************/
/* Compute the square root of a\in Z[\alpha], by projecting to a suitable   */
/* finite field Z[x]/<p, f>, and applying Cipolla (Bach & Shallit, sec 7.2) */
/* Careful, though: Our 'f' defines the finite field, and so is different   */
/* than their 'f'.                                                          */
/* If 'p' is not specified (i.e., _p=0), we will choose one.                */
/****************************************************************************/
{ mpz_t     p,  c, exponent, c1;
  mpz_poly  f, tmp, b0, b1, xPow0, xPow1, y0, y1, z, a, t;
  mpz_poly  factors[MAXPOLYDEGREE], sqrts[MAXPOLYDEGREE];
  int       d=_f->degree, i, j, attempt=0, equal, e[MAXPOLYDEGREE],n, signs;
  int       flip;
  
  mpz_init(p); mpz_init(c); mpz_init(c1); mpz_init(exponent);
  mpz_poly_init(f); mpz_poly_init(t); mpz_poly_init(z);
  mpz_poly_init(tmp); mpz_poly_init(b0);  mpz_poly_init(b1);
  mpz_poly_init(y0);  mpz_poly_init(y1);  mpz_poly_init(xPow0);
  mpz_poly_init(xPow1); mpz_poly_init(a);
  for (i=0; i<MAXPOLYDEGREE; i++)
    mpz_poly_init(factors[i]);
  for (i=0; i<MAXPOLYDEGREE; i++)
    mpz_poly_init(sqrts[i]);

  if (mpz_sgn(_p)==0)
    /* This will determine an appropriate prime p. */
    getSqrtP(p, _a, _f);
  else
    mpz_set(p, _p);

  mpz_poly_mod_pp(a, _a, _f, p);
  printf("Attempting to compute a square root of:\na=");
  mpz_poly_print(stdout, "", a);
  printf("in Z[x]/<f, p>, where:\nf=");
  mpz_poly_print(stdout, "", _f);
  printf("\np="); mpz_out_str(stdout, 10, p);
  printf("\n");

  printf("p = "); mpz_out_str(stdout, 10, p); printf("\n");
  n = mpz_poly_fact(factors, e, _f, p);
  if (n>1) {
    printf("f splits mod p:\n");
    for (i=0; i<n; i++) {
      printf("f_%d = ", i);
      mpz_poly_print(stdout, "", factors[i]);
      printf("\n");
    }
    printf("Using CRT mode for Zalpha_sqrt()...\n");
    printf("------------------------------------\n");
    mpz_set_ui(&res->coef[0], 0); res->degree = 0;
    for (i=0; i<n; i++) {
      if (Zalpha_sqrt(sqrts[i], a, factors[i], p, p)) {
        printf("An intermediate squareroot computation failed. Cannot continue!\n");
        msgLog(NULL, "Zalpha_sqrt(): intermediate sqrt computation failed. Cannot continue!");
        equal=-1;
        goto ZASQRT_CLEANUP;
      }
      printf("Sqrt succeeded modulo factor %d / %d :\n", i+1, n);
      mpz_poly_print(stdout, "  sqrt_i = ", sqrts[i]);
      printf("\n");
    }
    /* This is a bit tricky now: we need to try all possible
       combinations of square roots until we get a compatible combination.
       That is, (-1)^{e1}sqrts[0],... , (-1)^{en}sqrts[n-1].
       We'll know the right one when we see it as it should be the only
       one which holds in Z[\alpha].
    */
    
    signs=0;
    while (signs < BIT(n)) {    
      printf("Trying lift combination %d of %d...\n", signs+1, BIT(n));
      if (signs > 0) {
        /* Which signs need to be flipped? */
        flip = (signs-1)^signs;
        for (i=0; i<n; i++) {
          if (flip&BIT(i)) {
            for (j=0; j<=sqrts[i]->degree; j++)
              mpz_neg(&sqrts[i]->coef[j], &sqrts[i]->coef[j]);
          }
        }
      }
      /* Do CRT: */
      mpz_set_ui(&res->coef[0], 0); res->degree=0;
      for (i=0; i<n; i++) {      
        /* Compute tmp <-- (factors[0])...(factors[n])/factors[i] (=f/factors[i]). */
        mpz_set_ui(&tmp->coef[0], 1);
        tmp->degree = 0;
        for (j=0; j<n; j++) {
          if (j!=i)
            mpz_poly_mul(tmp, tmp, factors[j]);
        }
        /* z <-- tmp^(-1) mod factors[i]. */
        mpz_poly_inv(z, tmp, factors[i], p);
        mpz_poly_mulmod_p(t, sqrts[i], tmp, p);
        mpz_poly_mulmod_p(t, t, z, p);
        mpz_poly_add(res, res, t);
        for (j=0; j<=res->degree; j++)
          mpz_mod(&res->coef[j], &res->coef[j], p);
      }
      mpz_poly_mod(res, res, _f);
      /* Check to see if this is the right one: */
      equal = 1;
      mpz_div_2exp(c, p, 1);
      for (i=0; i<=res->degree; i++) 
        if (mpz_cmp(&(res->coef[i]), c) > 0)
          mpz_sub(&(res->coef[i]), &(res->coef[i]), p);
      printf("Final square root is figured to be:\n");
      mpz_poly_print(stdout, "", res); printf("\n");
      /* Check to see if this is correct: */
      mpz_poly_mulmod(z, res, res, _f);
      for (i=0; i<=_a->degree; i++) {
        mpz_sub(c1, &_a->coef[i], &z->coef[i]);
        mpz_mod(c1, c1, N);
        if (mpz_sgn(c1))
          equal = 0;
      }
      if (equal)
        goto ZASQRT_CLEANUP;
      signs++;
    }
    equal = 0;
    goto ZASQRT_CLEANUP;
  }
  
  
  /* Compute f <-- _f mod p. */
  mpz_poly_init(f);
  f->degree = d;
  for (i=0; i<=d; i++) 
    mpz_mod(&(f->coef[i]), &(_f->coef[i]), p);
  
  /* Make f monic. */
  mpz_invert(c, &(f->coef[d]), p);
  for (i=0; i<=d; i++)  {
    mpz_mul(&(f->coef[i]), &(f->coef[i]), c);
    mpz_mod(&(f->coef[i]), &(f->coef[i]), p);
  }


  /* Choose a random 't' */  
PICK_T:
  for (i=0; i<d; i++) {
    mpz_set_ui(&(t->coef[i]), rand());
    mpz_mod(&(t->coef[i]), &(t->coef[i]), p);
  }
  while (mpz_sgn(&(t->coef[d-1]))==0) {
    mpz_set_ui(&(t->coef[d-1]), rand());
    mpz_mod(&(t->coef[d-1]), &(t->coef[d-1]), p);
  }
  t->degree = d-1;
  
  
  /* Rather than test if t^2-4a is a square or not, we will */
  /* just continue on as if it were a non-q.residue, and    */
  /* check the result at the end.                           */

  
  /* Compute X^{(q+1)/2} mod <X^2 - tX + a>                 */
  mpz_pow_ui(exponent, p, d);
  mpz_add_ui(exponent, exponent, 1);
  mpz_div_2exp(exponent, exponent, 1);

  for (i=0; i<d; i++) {
    mpz_set_ui(&(b0->coef[i]), 0);
    mpz_set_ui(&(b1->coef[i]), 0);
    mpz_set_ui(&(xPow0->coef[i]), 0);
    mpz_set_ui(&(xPow1->coef[i]), 0);
  }
  b0->degree = b1->degree = 0;
  mpz_set_ui(&(b0->coef[0]), 1);
  xPow0->degree = xPow1->degree = 0;
  mpz_set_ui(&(xPow1->coef[0]), 1);

  while (mpz_sgn(exponent)) {
    if (mpz_odd_p(exponent)) {
      /* (y0,y1) <-- (b0,b1)*(xPow0, xPow1) mod <X^2 - tX + a>.                      */
      /*            = (b0*xPow0) + (b0*xPow1 + b1*xPow0)X + (b1*xPow1)X^2            */	    
      /*            = (b0*xPow0 - b1*xPow1*a) + (b0*xPow1 + b1*xPow0 + b1*xPow1*t)X. */
      mpz_poly_mulmod_pp(y0, b0, xPow0, f, p);
      mpz_poly_mulmod_pp(y1, b0, xPow1, f, p);
      mpz_poly_mulmod_pp(tmp, b1, xPow0, f, p);
      for (i=0; i<=tmp->degree; i++) {
        mpz_add(&(y1->coef[i]), &(y1->coef[i]), &(tmp->coef[i]));
	mpz_mod(&(y1->coef[i]), &(y1->coef[i]), p);
      }
      y1->degree = MAX(y1->degree, tmp->degree);
      while ((y1->degree >0) && (mpz_sgn(&(y1->coef[y1->degree]))==0))
        y1->degree -= 1;
      
      mpz_poly_mulmod_pp(tmp, b1, xPow1, f, p);
      
      mpz_poly_mulmod_pp(z, tmp, a, f, p);
      for (i=0; i<=z->degree; i++) {
        mpz_sub(&(y0->coef[i]), &(y0->coef[i]), &(z->coef[i]));
	mpz_mod(&(y0->coef[i]), &(y0->coef[i]), p);
      }
      y0->degree = MAX(y0->degree, z->degree);
      while ((y0->degree >0) && (mpz_sgn(&(y0->coef[y0->degree]))==0))
        y0->degree -= 1;

      mpz_poly_mulmod_pp(z, tmp, t, f, p);
      for (i=0; i<=z->degree; i++) {
        mpz_add(&(y1->coef[i]), &(y1->coef[i]), &(z->coef[i]));
	mpz_mod(&(y1->coef[i]), &(y1->coef[i]), p);
      }
      y1->degree = MAX(y1->degree, z->degree);
      while ((y1->degree >0) && (mpz_sgn(&(y1->coef[y1->degree]))==0))
        y1->degree -= 1;
      
      /* Now (b0, b1) <-- (y0, y1). */
      for (i=0; i<d; i++) {
        mpz_set(&(b0->coef[i]), &(y0->coef[i]));
        mpz_set(&(b1->coef[i]), &(y1->coef[i]));
      }
      i = y0->degree;
      while ((i>0) && (mpz_sgn(&(b0->coef[i]))==0))
        i--;
      b0->degree = i;
      i = y1->degree;
      while ((i>0) && (mpz_sgn(&(b1->coef[i]))==0))
        i--;
      b1->degree = i;
    }
    /* (y0,y1) <-- (xPow0,xPow1)*(xPow0, xPow1) mod <X^2 - tX + a>.                */
    /*            = (xPow0)^2 + (2*xPow0*xPow1)X + (xPow1^2)X^2                    */	    
    /*            = (xPow0^2 - xPow1^2*a) + (2*xPow0*xPow1 + xPow1^2*t)X.          */
    mpz_poly_mulmod_pp(y0, xPow0, xPow0, f, p);
    mpz_poly_mulmod_pp(y1, xPow0, xPow1, f, p);
    for (i=0; i<=y1->degree; i++) {
      mpz_add(&(y1->coef[i]), &(y1->coef[i]), &(y1->coef[i]));
      mpz_mod(&(y1->coef[i]), &(y1->coef[i]), p);
    }
      
    mpz_poly_mulmod_pp(tmp, xPow1, xPow1, f, p);
      
    mpz_poly_mulmod_pp(z, tmp, a, f, p);
    for (i=0; i<=z->degree; i++) {
      mpz_sub(&(y0->coef[i]), &(y0->coef[i]), &(z->coef[i]));
      mpz_mod(&(y0->coef[i]), &(y0->coef[i]), p);
    }
    y0->degree = MAX(y0->degree, z->degree);
    while ((y0->degree >0) && (mpz_sgn(&(y0->coef[y0->degree]))==0))
      y0->degree -= 1;
    
    
    mpz_poly_mulmod_pp(z, tmp, t, f, p);
    for (i=0; i<=z->degree; i++) {
      mpz_add(&(y1->coef[i]), &(y1->coef[i]), &(z->coef[i]));
      mpz_mod(&(y1->coef[i]), &(y1->coef[i]), p);
    }
    y1->degree = MAX(y1->degree, z->degree);
    while ((y1->degree >0) && (mpz_sgn(&(y1->coef[y1->degree]))==0))
      y1->degree -= 1;

    /* (xPow0, xPow1) <-- (y0, y1). */
    for (i=0; i<d; i++) {
      mpz_set(&(xPow0->coef[i]), &(y0->coef[i]));
      mpz_set(&(xPow1->coef[i]), &(y1->coef[i]));
    }
    i = y0->degree;
    while ((i>0) && (mpz_sgn(&(xPow0->coef[i]))==0))
      i--;
    xPow0->degree = i;
    i = y1->degree;
    while ((i>0) && (mpz_sgn(&(xPow1->coef[i]))==0))
      i--;
    xPow1->degree = i;

    mpz_div_2exp(exponent, exponent, 1);
  }
  
  mpz_poly_mulmod_pp(z, b0, b0, f, p);
  res->degree = b0->degree;
  for (i=0; i<=res->degree; i++)
    mpz_set(&(res->coef[i]), &(b0->coef[i]));
  
  equal = 1;
  mpz_div_2exp(c, p, 1);
  for (i=0; i<=res->degree; i++) 
    if (mpz_cmp(&(res->coef[i]), c) > 0)
      mpz_sub(&(res->coef[i]), &(res->coef[i]), p);
  for (i=0; i<=a->degree; i++) {
    mpz_sub(c1, &(a->coef[i]), &(z->coef[i]));
    mpz_mod(c1, c1, p);
    if (mpz_sgn(c1))
      equal = 0;
  }
  if (!(equal) && ((++attempt) < 16)) {
    printf("Zalpha: square root mod p failed. Trying again...\n");
    goto PICK_T;
  }
  if (equal) {
    printf("Zalpha: SQRT mod p succeeded!\n");
    printf("res = "); mpz_poly_print(stdout, "", res); printf("\n");
    mpz_invert(c, &(_f->coef[d]), N);
    for (i=0; i<=_f->degree; i++) {
      mpz_mul(&(f->coef[i]), &(_f->coef[i]), c);
      mpz_mod(&(f->coef[i]), &(f->coef[i]), N);
    }
    f->degree = _f->degree;
    mpz_poly_mulmod_pp(z, res, res, f, N);
    equal = 1;
    for (i=0; i<=a->degree; i++) {
      mpz_sub(c1, &(z->coef[i]), &(a->coef[i]));
      mpz_mod(c1, c1, N);
      if (mpz_sgn(c1))
        equal = 0;
    }
    if (equal)
      printf("Computed square root holds mod N.\n");
    else
      printf("Computed square root is invalid mod N.\n");
/* Try the other square root, though it shouldn't matter. */
  }

ZASQRT_CLEANUP:  
  for (i=0; i<MAXPOLYDEGREE; i++)
    mpz_poly_init(factors[i]);
  mpz_clear(p); mpz_clear(c); mpz_clear(c1); mpz_clear(exponent);
  mpz_poly_clear(f); mpz_poly_clear(t); mpz_poly_clear(z);
  mpz_poly_clear(tmp); mpz_poly_clear(b0);  mpz_poly_clear(b1);
  mpz_poly_clear(y0);  mpz_poly_clear(y1);  mpz_poly_clear(xPow0);
  mpz_poly_clear(xPow1); mpz_poly_clear(a);
  if (!(equal))
    return -1;
  return 0;
  
}



