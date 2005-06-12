/**************************************************************/
/* nfmisc.c                                                   */
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

/* This code is still in debugging stages, and not ready to be cleaned
   up yet!
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ggnfs.h"


#define MSG_LOG GGNFS_LOG_NAME


#define _HEAVY_DEBUG

void getTraceConstants(nf_t *N);

/****************************************************************/
int smallestFactor(mpz_fact_t *F)
/****************************************************************/
{ int minLoc, i;

  minLoc = -1;
  for (i=0; i<F->size; i++)
    if ((F->exponents[i]>0) && 
        ((minLoc==-1) || (mpz_cmp(&F->factors[i], &F->factors[minLoc]) < 0)))
      minLoc = i;
  return minLoc;
}

/****************************************************************/
void fixOrder(mpz_mat_t *O, mpz_t d)
/****************************************************************/
/* Make the denominator right.                                  */
/****************************************************************/
{ int   i, j, n=O->rows;
  mpz_t g;

  mpz_init_set(g, d);
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      mpz_gcd(g, g, &O->entry[i][j]);
  if (mpz_cmp_si(g, 1)) {
    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
        mpz_div(&O->entry[i][j], &O->entry[i][j], g);
    mpz_div(d, d, g);
  }
  mpz_clear(g);
}

/****************************************************************/
mpz_poly *initMultTable(int n)
/****************************************************************/
{ int       i;
  mpz_poly *M;

  M = (mpz_poly *)malloc(sizeof(mpz_poly)*n*n);
  for (i=0; i<n*n; i++)
    mpz_poly_init(M[i]);
  return M;
}

/****************************************************************/
void clearMultTable(mpz_poly *M, int n)
/****************************************************************/
{ int i;

  for (i=0; i<n*n; i++)
    mpz_poly_clear(M[i]);
  free(M);
}

/****************************************************************/
void computeMultTable(mpz_poly *Mt, mpz_mat_t *M, mpz_t d, mpz_poly T)
/****************************************************************/
/* Given an order, specified by a Z-basis M in HNF and          */
/* denominator d, compute a multiplication table for the        */
/* generators. That is, if M_j is the j-th column of M, w_j the */
/* corresponding polynomial in Z[alpha] and omega_j = w_j/d,    */
/* Then Mt[i*n +j] = omega_i * omega_j, expressed as a poly in  */
/* the omega_k.                                                 */
/****************************************************************/
{ int      i, j, k, l, n=T->degree;
  mpz_poly prod, w1, w2, wprod;
  mpz_t    tmp, tmp2;

  mpz_poly_init(prod); mpz_poly_init(w1); mpz_poly_init(w2);
  mpz_poly_init(wprod); 
  mpz_init(tmp); mpz_init(tmp2);


  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      getCol(w1, M, i);
      getCol(w2, M, j);
      mpz_poly_mulmod(prod, w1, w2, T);
      /**************************************************************/
      /* A little easy algebra shows that if, as above,             */
      /* omega_j = w_j/d, where w_j is a poly with integer entries, */
      /* and (omega_i)(omega_j) = a_1*omega_1 + ... + a_n*omega_n,  */
      /* then (w_i)(w_j)/d = a_1*w_1 + ... + a_n*w_n.               */
      /* In particular, the product w_i*w_j had better be divisible */
      /* by d!                                                      */
      /**************************************************************/
      for (k=0; k<=prod->degree; k++) {
        mpz_mod(tmp, &prod->coef[k], d);
        if (mpz_sgn(tmp)) {
          printf("computeMultTable(): Sever error - product not divisible by d!\n");
          exit(-1);
        }
        mpz_div(&prod->coef[k], &prod->coef[k], d);
      }
      /************************************************/ 
      /* Now express prod as a polynomial in the w_k. */
      /* Since M is in HNF, this is easy.             */
      /************************************************/ 
      wprod->degree = prod->degree;
      for (k=prod->degree; k>=0; k--) {
        mpz_mod(tmp, &prod->coef[k], &M->entry[k][k]);
        if (mpz_sgn(tmp)) {
          printf("Some error occurred: prod does not appear to be in O!\n");
          exit(-1);
        }
        mpz_div(tmp, &prod->coef[k], &M->entry[k][k]);
        for (l=0; l<n; l++) {
          mpz_mul(tmp2, tmp, &M->entry[l][k]);
          mpz_sub(&prod->coef[l], &prod->coef[l], tmp2);
        }
        mpz_set(&wprod->coef[k], tmp);
      }
      while ((wprod->degree > 0) && (mpz_sgn(&wprod->coef[wprod->degree])==0))
        wprod->degree -= 1;
      mpz_poly_cp(Mt[i*n+j], wprod);
    }
  }
  mpz_poly_clear(prod); mpz_poly_clear(w1); mpz_poly_clear(w2);
  mpz_poly_clear(wprod);
  mpz_clear(tmp); mpz_clear(tmp2);
}

/**************************************************************/
int precomputeSpPowers(nf_t *N)
/**************************************************************/
/* Precompute some data for powers of special ideals, to make */
/* computation of the valuations of <a-b\alpha> at these      */
/* ideals faster.                                             */
/**************************************************************/
/* We assume N->sPowMats has not been allocated.              */
/**************************************************************/
{ int       i, e;
  mpz_mat_t T;

  mpz_mat_init(&T);
  
  for (i=0; i<N->numSPrimes; i++) {
    for (e=1; e<=MAX_SP_POWERS; e++) {
      getIdealHNF_ib(&T, N->sPrimes[i].p, N->sPrimes[i].alpha, e, N);
      /* Now if P is this special ideal, we have that an
         <a-b\alpha> ideal is divisible by P^e iff
         (c_d*a, -b,0,...,0) is in the column space of this matrix T.
         Luckily, T is upper triangular, so we need remember only
         the upper-left 2x2 submatrix to be able to test this condition
         later.
      */
      mpz_mat_init2(&N->sPowMats[i][e-1], 2, 2);
      mpz_set(&N->sPowMats[i][e-1].entry[0][0], &T.entry[0][0]);
      mpz_set(&N->sPowMats[i][e-1].entry[0][1], &T.entry[0][1]);
      mpz_set(&N->sPowMats[i][e-1].entry[1][0], &T.entry[1][0]); /* should always be zero. */
      mpz_set(&N->sPowMats[i][e-1].entry[1][1], &T.entry[1][1]);
      N->sPowMats[i][e-1].rows = N->sPowMats[i][e-1].cols = 2;
    }
  }
  mpz_mat_clear(&T);
  return 0;
} 


/****************************************************************/
int stdtoib(mpz_poly res, mpz_poly op, nf_t *N)
/****************************************************************/
/* Input: A polynomial op in alpha (where the number field is   */
/*        Q(alpha). The number field, N, with all good fields   */
/*        already filled in.                                    */
/* Output: res expressed as a polynomial in the omega_i (the    */
/*        integral basis). It is assumed that op is an          */
/*        algebraic integer so that this is possible.           */
/****************************************************************/
{ int      i, j;
  static   mpz_poly tpol1, tpol2;
  static   mpz_t    c, tmp;
  static   int initialized=0;
 
  if (!(initialized)) {
    mpz_poly_init(tpol1); mpz_poly_init(tpol2);
    mpz_init(c); mpz_init(tmp);
    initialized=1;
  }

  /* Let omega_i = w_i/d. Then the easiest way to proceed is */
  /* to write d*op = a_1*w_1 + ... + a_n*w_n.                */
  mpz_poly_cp(tpol1, op);
  for (i=0; i<=tpol1->degree; i++) 
    mpz_mul(&tpol1->coef[i], &tpol1->coef[i], N->W_d);

  tpol2->degree = op->degree;
  for (i=0; i<=tpol2->degree; i++)
    mpz_set_ui(&tpol2->coef[i], 0);

  for (i=tpol1->degree; i>=0; i--) {
#ifdef TEST_DIVISIONS
    mpz_mod(c, &tpol1->coef[i], &N->W->entry[i][i]);
    if (mpz_sgn(c)) {
      printf("stdtoib(): Error - op is not an algbraic integer!\n");
      return -1;
    }
#endif
    mpz_div(c, &tpol1->coef[i], &N->W->entry[i][i]);
    mpz_set(&tpol2->coef[i], c);
    for (j=0; j<=i; j++) {
      mpz_mul(tmp, c, &N->W->entry[j][i]);
      mpz_sub(&tpol1->coef[j], &tpol1->coef[j], tmp);
    }
  }
  mpz_poly_fixDeg(tpol1);
  if ((tpol1->degree > 0) || mpz_sgn(&tpol1->coef[0]))  {
    printf("stdtoib(): Error - op is not an algbraic integer!\n");
    return -1;
  }

  mpz_poly_cp(res, tpol2);
  return 0;
}  



/****************************************************************/
void basisPolyMult(mpz_poly res, mpz_poly op1, mpz_poly op2, mpz_poly *Mt, int n)
/****************************************************************/
/* Given polynomials op1 and op2 in the omega_i, compute the    */
/* product, res <-- op1*op2, expressed as a polynomial in the   */
/* omega_i. Mt is the multiplication table of the omega_i.      */
/****************************************************************/
{ int      i, j, k, s;
  static   mpz_poly tRes, tpol1;
  static   mpz_t    c, tmp;
  static   int initialized=0;

  if (!(initialized)) {
    mpz_poly_init(tRes);
    mpz_poly_init(tpol1);
    mpz_init(c); mpz_init(tmp);
    initialized=1;
  }


  tRes->degree = n-1;
  for (i=0; i<= tRes->degree; i++)
    mpz_set_ui(&tRes->coef[i], 0);

  for (i=0; i<=op1->degree; i++) {
    for (j=0; j<=op2->degree; j++) {
      mpz_mul(c, &op1->coef[i], &op2->coef[j]);
      /* It's worth checking, since we often use this */
      /* function with at least one sparse polynomial. */
      if (mpz_sgn(c)) {
        s = Mt[i*n+j]->degree;
        for (k=0; k<=s; k++) {
          mpz_mul(tmp, c, &Mt[i*n+j]->coef[k]);
          mpz_add(&tRes->coef[k], &tRes->coef[k], tmp);
        }
      }
    }
  }
  mpz_poly_cp(res, tRes);

}


/****************************************************************/
void basisPolyPowM(mpz_poly res, mpz_poly base, mpz_t e, mpz_poly *Mt, mpz_t p, int n)
/****************************************************************/
{ mpz_t    remain;
  mpz_poly tRes, tpol1;
  int      i;

  mpz_poly_init(tRes);
  mpz_poly_init(tpol1);
  mpz_init(remain);

  mpz_set_ui(&tRes->coef[0], 1);
  tRes->degree = 0;
  mpz_set(remain, e);
  mpz_poly_cp(tpol1, base);

  while (mpz_sgn(remain)) {
    if (mpz_odd_p(remain)) {
      basisPolyMult(tRes, tRes, tpol1, Mt, n);
      for (i=tRes->degree; i>=0; i--)
        mpz_mod(&tRes->coef[i], &tRes->coef[i], p);
      mpz_poly_fixDeg(tRes);
      mpz_sub_ui(remain, remain, 1);
    }
    mpz_div_2exp(remain, remain, 1);
    basisPolyMult(tpol1, tpol1, tpol1, Mt, n);
    for (i=tpol1->degree; i>=0; i--)
      mpz_mod(&tpol1->coef[i], &tpol1->coef[i], p);
    mpz_poly_fixDeg(tpol1);
  }
  mpz_poly_cp(res, tRes);

  mpz_poly_clear(tRes);
  mpz_poly_clear(tpol1);
  mpz_clear(remain);
}

/****************************************************************/
void idealHNF_ib(mpz_mat_t *H, mpz_poly a, nf_t *N)
/****************************************************************/
/* Compute the HNF representation, wrt the integral basis, of   */
/* the principal ideal <a>O, where 'a' is given wrt the I.B.    */
/****************************************************************/
{ int      i, n=N->degree;
  static mpz_poly tpol1, tpol2;
  static int initialized=0;

  if (!(initialized)) {
    mpz_poly_init(tpol1); mpz_poly_init(tpol2);
    initialized=1;
  }
  for (i=0; i<n; i++) {
    /* tpol1 <-- omega_i, wrt the omega_i basis. */
    if (i>0)
      mpz_set_ui(&tpol1->coef[i-1], 0);
    mpz_set_ui(&tpol1->coef[i], 1);
    tpol1->degree = i;
    basisPolyMult(tpol2, tpol1, a, N->Mt, n);
    setCol(H, i, tpol2);
  }
  H->rows = H->cols = n;
  mpz_mat_getHNF(H, H);
}

/****************************************************************/
void idealHNF_ib_ab(mpz_mat_t *H, s64 a, s32 b, nf_t *N)
/****************************************************************/
/* Compute the HNF representation, wrt the integral basis, of   */
/* the principal ideal <a*c_d-b\hat{\alpha}>O.                  */
/****************************************************************/
{ int      i, j, n=N->degree, s1, s2, c;
  static mpz_poly tpol1, tpol2;
  static mpz_t tmp, tmp2, bmultiplier;
  static mpz_mat_t tmpH, t2;
  static int initialized=0;
  __mpz_struct *q;

  if (!(initialized)) {
    mpz_poly_init(tpol1); mpz_poly_init(tpol2);
    mpz_init(tmp); mpz_init(tmp2); mpz_mat_init(&tmpH);
    mpz_mat_init(&t2);
    mpz_init_set(bmultiplier, N->W_d);
    mpz_mod(tmp, N->W_d, &N->W->entry[1][1]);
    if (mpz_sgn(tmp)) {
      msgLog(GGNFS_LOG_NAME,"Hmm. Something strange is afoot with bmultiplier!");
      mpz_set_ui(bmultiplier,1);
    } else
      mpz_div(bmultiplier, bmultiplier, &N->W->entry[1][1]);
    initialized=1;
  }
/* CHANGE: right-rotate the columns by one here.
   It should probably make the HNF reduction
   more efficient.
*/

/* CJM, 12/1. */
mpz_mul_si(tmp2, bmultiplier, -b);

  mpz_mul_si64(tmp, &N->f->coef[n], a);
  for (i=0; i<n; i++) {
    /* Set the i-th column of H equal to             */
    /* omega_i*(a - b*bmultiplier\omega_1) =         */
    /* a*omega_i - b*bmultiplier(omega_i*omega_1).   */
    c = i;
    s1 = i*n+1;
    s2 = N->Mt[s1]->degree;
    q = N->Mt[s1]->coef;
    for (j=0; j<=s2; j++) {
/* CJM, 12/1. */
#if 0
      mpz_mul_si(&tmpH.entry[j][c], &q[j], -b);
      mpz_mul(&tmpH.entry[j][c], &tmpH.entry[j][c], bmultiplier);
#else
      mpz_mul(&tmpH.entry[j][c], &q[j], tmp2);
#endif
    }
    for (; j<n; j++)
      mpz_set_ui(&tmpH.entry[j][c], 0);
    mpz_add(&tmpH.entry[i][c], &tmpH.entry[i][c], tmp);
  }
  tmpH.rows = tmpH.cols = n;
  mpz_mat_getHNF(H, &tmpH);

}


/****************************************************************/
int computePRadical(mpz_mat_t *Ip, mpz_t Ip_d, mpz_mat_t *O, mpz_t d, 
                    mpz_poly *Mt, mpz_poly T, mpz_t p)
/****************************************************************/
/* Compute the p-radical as suggested on p. 309 of Cohen.       */
/* Input: (O, d) is an order. Mt is the multiplication table    */
/*        for that order. T is the defining, monic poly, and    */
/*        p is obviously the prime for which the p-radical is   */
/*        needed.                                               */
/* Output: (Ip, Ip_d) is the p-radical of the given order.      */
/****************************************************************/
{ int       i, j, k, l, n=T->degree;
  mpz_t     q, tmp1, tmp2, modulus, d_q_1;
  mpz_poly  apol, tpol1;
  mpz_mat_t A, B;

  mpz_init(q); mpz_init(tmp1); mpz_init(tmp2);
  mpz_init(d_q_1);
  mpz_init(modulus);
  mpz_poly_init(apol);
  mpz_poly_init(tpol1);
  mpz_mat_init2(&A, n, n);
  mpz_mat_init2(&B, n, 2*n);

  mpz_set(q, p);
  while (mpz_cmp_ui(q, n) < 0)
    mpz_mul(q, q, p);

  for (j=0; j<n; j++) {
    /***************************************************/
    /* Compute a_j <-- (omega_j)^q, expressed          */
    /* as a polynomial in the omega_i. More precisely, */
    /* we only need the coefficients of a_j modulo p.  */
    /* Then, set the j-th column of A, A_j <-- a_j.    */ 
    /***************************************************/
    tpol1->degree = j;
    for (k=0; k<j; k++)
      mpz_set_ui(&tpol1->coef[k], 0);
    mpz_set_ui(&tpol1->coef[j], 1);
    basisPolyPowM(apol, tpol1, q, Mt, p, n);
    setCol(&A, j, apol);
  }
  A.rows = A.cols = n;
  mpz_mat_getKernel_modp(&B, &A, p);

  mpz_mat_mul(&B, O, &B);

  /**********************************************/
  /* Now, B is Ip/pO. Append with the p*omega_i */
  /* entries and do HNF to get Ip.              */
  /**********************************************/
  l = B.cols;
  for (j=0; j<n; j++)
    for (i=0; i<n; i++)
      mpz_mul(&B.entry[i][j+l], &O->entry[i][j], p);
  B.cols += n;

  mpz_mat_getHNF(Ip, &B);
  mpz_set(Ip_d, d);
  fixOrder(Ip, Ip_d);

  mpz_clear(q); mpz_clear(tmp1); mpz_clear(tmp2);
  mpz_clear(d_q_1); mpz_clear(modulus);
  mpz_poly_clear(apol);
  mpz_poly_clear(tpol1);
  mpz_mat_clear(&A);
  mpz_mat_clear(&B);
  return 0;
}

/****************************************************************/
void computePRadical2(mpz_mat_t *Ip, mpz_poly *Mt, mpz_poly T, mpz_t p)
/****************************************************************/
/* Like above, with 2 major differences:                        */
/* (1) Compute the p-radical modulo pO.                         */
/* (2) The basis returned is an F_p basis w.r.t O.              */
/****************************************************************/
{ int       j, k, n=T->degree;
  mpz_t     q;
  mpz_poly  apol, tpol1;
  mpz_mat_t A;

  mpz_init(q);
  mpz_poly_init(apol);
  mpz_poly_init(tpol1);
  mpz_mat_init2(&A, n, n);

  mpz_set(q, p);
  while (mpz_cmp_ui(q, n) < 0)
    mpz_mul(q, q, p);

  for (j=0; j<n; j++) {
    /***************************************************/
    /* Compute a_j <-- (omega_j)^q, expressed          */
    /* as a polynomial in the omega_i. More precisely, */
    /* we only need the coefficients of a_j modulo p.  */
    /* Then, set the j-th column of A, A_j <-- a_j.    */ 
    /***************************************************/
    tpol1->degree = j;
    for (k=0; k<j; k++)
      mpz_set_ui(&tpol1->coef[k], 0);
    mpz_set_ui(&tpol1->coef[j], 1);
    basisPolyPowM(apol, tpol1, q, Mt, p, n);
    setCol(&A, j, apol);
  }
  A.rows = A.cols = n;
  mpz_mat_getKernel_modp(Ip, &A, p);

  mpz_clear(q);
  mpz_poly_clear(apol);
  mpz_poly_clear(tpol1);
  mpz_mat_clear(&A);
}



/****************************************************************/
int thm6_1_3(mpz_mat_t *O, mpz_t d, mpz_mat_t *Op, mpz_t dp, 
             mpz_poly T, mpz_t p)
/****************************************************************/
/* Given an order Op,dp and a prime p, compute O' as in Cohen,  */
/* Thm. 6.1.3.                                                  */
/****************************************************************/
{ mpz_mat_t Ip, A, U;
  mpz_t     t1, t2, Ip_d;
  int       i, j, k, l, n = T->degree;
  mpz_poly  w_k, gam_i, pol1, pol2, *Mt;

  Mt = initMultTable(n);
  computeMultTable(Mt, Op, dp, T);

  mpz_mat_init2(&Ip, n, n);
  mpz_mat_init2(&A, n*n, n);
  mpz_mat_init2(&U, n, n*(n+1));
  mpz_init(t1); mpz_init(t2); mpz_init(Ip_d);
  mpz_poly_init(w_k);
  mpz_poly_init(gam_i);
  mpz_poly_init(pol1);
  mpz_poly_init(pol2);

  computePRadical(&Ip, Ip_d, Op, dp, Mt, T, p);

  /**********************************************************/
  /* We know now the p-radical, Ip.                         */
  /* Use this in Lemma 6.1.7 of Cohen. In particular,       */
  /* note that the n^2 x n matrix should be easy to compute */
  /* because we have a nice triangular basis (this is       */
  /* the alternative method he mentions @ the bottom        */
  /*  of p. 310.                                            */
  /**********************************************************/
  for (i=0; i<n*n; i++)
    for (j=0; j<n; j++)
      mpz_set_ui(&A.entry[i][j], 0);

  /* Compute the big matrix. */
  for (k=0; k<n; k++) {
    getCol(w_k, Op, k);
    for (i=0; i<n; i++) {
      getCol(gam_i, &Ip, i);
      mpz_poly_mul(pol1, w_k, gam_i);
      mpz_poly_mod(pol1, pol1, T);

      /********************************************************/
      /* Now, observe: If  omega_k = w_k/dw and               */
      /* gamma_j = c_j/dc, and                                */
      /* omega_k * gamma_j = a_1*gamma_1 + ... + a_n*gamma_n, */
      /* Then we have:                                        */
      /* (w_k * c_j)/(dw*dc) = (a_1*c_1 + ... + a_n*c_n)/dc.  */
      /* Hence,                                               */
      /* (w_k * c_j)/dw = a_1*c_1 + ... + a_n*c_n.            */
      /********************************************************/
      for (j=0; j<=pol1->degree; j++) {
        mpz_mod(t1, &pol1->coef[j], dp);
        if (mpz_sgn(t1)) {
          printf("thm6_1_3(): Severe error - t1 != 0.\n");
          exit(-1);
        }
        mpz_div(&pol1->coef[j], &pol1->coef[j], dp);
      }
      /*********************************************************/
      /* Now, find the representation of pol1 on the basis Ip. */
      /* Since Ip is in HNF, this is easy.                     */
      /*********************************************************/
      for (j=pol1->degree; j>=0; j--) {
      mpz_mod(t1, &pol1->coef[j], &Ip.entry[j][j]);
      if (mpz_sgn(t1)) {
        printf("Some error occurred: pol1 does not appear to be in Ip!\n");
        exit(-1);
      }
      mpz_div(t1, &pol1->coef[j], &Ip.entry[j][j]);
      for (l=0; l<n; l++) {
        mpz_mul(t2, t1, &Ip.entry[l][j]);
        mpz_sub(&pol1->coef[l], &pol1->coef[l], t2);
      }
      mpz_mod(&A.entry[i*n+j][k], t1, p);
      }
    }
  }
  A.rows = n*n; 
  A.cols = n;
  mpz_mat_getKernel_modp(&U, &A, p);         
  /***************************************************************/
  /* Now U should be the kernel of the map in Lemma 6.1.7 mod p. */
  /* To get the kernel we need, supplement with the p*w_j and    */
  /* get the HNF.                                                */
  /***************************************************************/
  mpz_mat_mul(&U, Op, &U);
  l = U.cols;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      mpz_mul(&U.entry[i][j+l], &Op->entry[i][j], p);
  U.cols += n;
  mpz_mat_getHNF(O, &U);
  mpz_mul(d, dp, p);
  fixOrder(O, d);
  
  clearMultTable(Mt, n);
  mpz_mat_clear(&Ip);
  mpz_mat_clear(&A);
  mpz_mat_clear(&U);
  mpz_clear(t1); mpz_clear(t2); mpz_clear(Ip_d);
  mpz_poly_clear(w_k);
  mpz_poly_clear(gam_i);
  mpz_poly_clear(pol1);
  mpz_poly_clear(pol2);
  return 0;
}

/****************************************************************/
int getPMaximal(mpz_mat_t *Op, mpz_t dp, mpz_poly T, mpz_t p)
/****************************************************************/
/* Input: An irreducible (not necessarily monic) polynomial, T. */
/*        and a prime p dividing disc(T).                       */
/* Output: A matrix Op and integer d containing the HNF rep.    */
/*        of a p-maximal order O_p.                             */
/****************************************************************/
/* We will compute Op by repeated application of Cohen,         */
/* Thm 6.1.3.                                                   */
/****************************************************************/
{ mpz_mat_t O;
  mpz_t     d;
  int       i, j, n=T->degree, pMaximal;

  mpz_mat_init2(&O, n, n);
  mpz_init(d);
  /* Initialize Op. */
  Op->rows = Op->cols = n;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      mpz_set_si(&Op->entry[i][j], (i==j));
  mpz_set_ui(dp, 1);
  pMaximal=0;
  do {
    thm6_1_3(&O, d, Op, dp, T, p);
    if ((mpz_cmp(d, dp)==0) && mpz_mat_equal(&O, Op))
      pMaximal = 1;
    else {
      mpz_mat_cp(Op, &O);
      mpz_set(dp, d);
    }
  } while (!(pMaximal));
  mpz_mat_clear(&O);
  mpz_clear(d);
  return 1;
}

/****************************************************************/
int getIntegralBasis(nf_t *N, mpz_fact_t *D, int tryHard)
/****************************************************************/
/* Input: N should have all it's fields initialized, and have   */
/*        an irreducible monic polynomial, T.                   */
/*        The discriminant, D, of T.                            */
/* Output: An basis for the ring of integers of K=Q[x]/<T>,     */
/*         and the discriminant of K.                           */
/* See Cohen, Alg. 6.1.8. Then take away the Dedekind test, and */
/*         be imaginative - the entire ``detailed'' description */
/*         of Alg. 6.1.8 does not mention the denominator once! */
/*         I tried like heck, but couldn't get that alg to work */
/*         out of the box. Adding the Dedekind test now would   */
/*         be a trivial matter, but I'm not too concerned with  */
/*         it at the moment. Things are working just fine as-is.*/
/****************************************************************/
/* If tryHard is zero, the returned basis might not be maximal, */
/* but it should be p-maximal with respect to all reasonably    */
/* sized primes (probably upto 15 digits or so.)                */
/****************************************************************/
{ int        retVal, n=N->T->degree, i, j, l, pLoc, m;
  mpz_fact_t F;
  mpz_t      p, d_Op;
  mpz_mat_t  Op;
  mpf_t      cdf;
  mpz_poly   tpol1, tpol2;

  mpz_fact_init(&F);
  mpz_init(p); mpz_init(d_Op);
  mpz_mat_init2(&Op, n, n);
  mpf_init2(cdf, 256);
  mpz_poly_init(tpol1); mpz_poly_init(tpol2);

  retVal = mpz_fact_check(D, tryHard);
  if (retVal && tryHard) {
    msgLog(MSG_LOG, "getIntegralBasis(): Warning - Could not factor discriminant!");
#ifdef _MAXIMAL_ORDER_ONLY
    mpz_fact_clear(&F); mpz_clear(p); mpz_clear(d_Op);
    mpz_mat_clear(&Op); mpf_clear(cdf);
    return -1;
#else
    msgLog(MSG_LOG, "Proceeding anyway, assuming the p-maximal order is maximal!");
#endif
  }
  mpz_set(N->Tdisc, D->N);
  mpz_fact_removeSF(&F, D);
  
  /* Initialize with W <-- Z[theta]. */
  mpz_mat_setID(N->W, n);
  mpz_set_ui(N->W_d, 1);


  while (mpz_cmp_ui(F.N, 1)) {
    /* Get the next prime to maximize at: */
    pLoc = smallestFactor(&F);
    if (pLoc < 0) {
      printf("Some error occurred getting next prime to maximize at!\n");
      mpz_fact_clear(&F); mpz_clear(p); mpz_clear(d_Op);
      mpz_mat_clear(&Op); mpf_clear(cdf);
      mpz_poly_clear(tpol1); mpz_poly_clear(tpol2);
      return -1;
    }
    mpz_set(p, &F.factors[pLoc]);
    /* Get a p-maximal order, Op. */
    getPMaximal(&Op, d_Op, N->T, p);
    /* Do W <-- W + Op */
    mpz_mat_moduleAdd(N->W, N->W_d, N->W, N->W_d, &Op, d_Op);
    /* Remove occurrences of p from F.N, since W is presumably p-maximal now. */
    mpz_remove(F.N, F.N, p);
    F.exponents[pLoc] = 0;
  }
    
  /* Compute the field discriminant. */
  if (mpz_sgn(N->Tdisc)>0)
    mpz_set_si(N->Kdisc, 1);
  else
    mpz_set_si(N->Kdisc, -1);
  for (i=0; i<n; i++) {
    mpz_mul(N->Kdisc, N->Kdisc, N->W_d);
    mpz_div(N->Kdisc, N->Kdisc, &N->W->entry[i][i]);
  }
  mpz_mul(N->Kdisc, N->Kdisc, N->Kdisc);
  mpz_div(N->Kdisc, N->Tdisc, N->Kdisc);


  mpz_div(N->index, N->Tdisc, N->Kdisc);
  mpz_sqrt(N->index, N->index);
  N->degree = N->T->degree;
  N->Mt = initMultTable(n);
  computeMultTable(N->Mt, N->W, N->W_d, N->T);

  /* Get the zeros of f. */
  mpz_poly_getComplexZeros(N->fZeros, N->f);

  /* Order them in the expected way: */
  reorderRoots(N->fZeros, n);
  /* And compute the zeros of T as tZero = c_d*fZero. */
  mpf_set_z(cdf, &N->f->coef[n]);
  for (i=0; i<n; i++) {
    mpf_mul(N->TZeros[i].mpr, N->fZeros[i].mpr, cdf);
    mpf_mul(N->TZeros[i].mpi, N->fZeros[i].mpi, cdf);
    N->TZeros[i].r = mpf_get_d(N->TZeros[i].mpr);
    N->TZeros[i].i = mpf_get_d(N->TZeros[i].mpi);
  }

  /* Compute the decompositions of the special primes: */
  mpz_set(F.N, N->index);
  F.size = 0;
  retVal = mpz_fact_check(&F, tryHard);
#if 0
  if ((retVal) && (tryHard)) {
    msgLog(MSG_LOG, "getIntegralBasis(): Fatal - Could not factor index!");
    mpz_fact_clear(&F); mpz_clear(p); mpz_clear(d_Op);
    mpz_mat_clear(&Op); mpf_clear(cdf);
    mpz_poly_clear(tpol1); mpz_poly_clear(tpol2);
    return -1;
  }
#endif
  /* The maximum number of special prime ideals is: */
  m = F.size * n;
  if (!(N->sPrimes = (prime_id_t *)malloc((1+m)*sizeof(prime_id_t)))) {
    fprintf(stderr, "getIntegralBasis(): Memory allocation error for sPrimes!\n");
    mpz_fact_clear(&F); mpz_clear(p); mpz_clear(d_Op);
    mpz_mat_clear(&Op); mpf_clear(cdf);
    mpz_poly_clear(tpol1); mpz_poly_clear(tpol2);
    return -1;
  }
  for (i=0; i<m; i++)
    initIdeal(&N->sPrimes[i]);
  for (i=m=0; i<F.size; i++) { 
    m += factorPrime(&N->sPrimes[m], &F.factors[i], N);
  }
  N->numSPrimes = m;
  /* Some constants used to compute valuations. */
  precomputeSpPowers(N);

/***/
  for (i=0; i<m; i++) {
    /* Compute the betaMat's. */
    N->sPrimes[i].betaMat.rows = N->sPrimes[i].betaMat.cols = n;
    for (j=0; j<n; j++) {
      /* The j-th column of betaMat is beta*alpha^j. */
      for (l=0; l<j; l++)
        mpz_set_ui(&tpol1->coef[l], 0);
      mpz_set_ui(&tpol1->coef[j], 1);
      tpol1->degree = j;
      basisPolyMult(tpol2, tpol1, N->sPrimes[i].beta, N->Mt, n);
      setCol(&N->sPrimes[i].betaMat, j, tpol2);
    }
  }


  if (!(N->v_cd_sPrimes = (int *)malloc((1+m)*sizeof(int)))) {
    fprintf(stderr, "getIntegralBasis(): Memory allocation error for v_cd_ePrimes!\n");
    return -1;
  }
  /* Compute the valuation of c_d at the i-th special prime. */
  mpz_mat_setID(&Op, n);
  /* Get the HNF of <c_d>. */
  for (i=0; i<n; i++)
    mpz_set(&Op.entry[i][i], &N->f->coef[n]);
  for (i=0; i<m; i++) {
    /* Compute the valuation of the i-th exceptional prime at <c_d>. */
    N->v_cd_sPrimes[i] = valuation(&Op, &N->sPrimes[i], N);
  }



  mpz_mat_pseudoInvert(N->W_inv, p, N->W);
  getTraceConstants(N);


  mpz_fact_clear(&F);
  mpz_clear(p); mpz_clear(d_Op);
  mpz_mat_clear(&Op);
  mpf_clear(cdf);
  mpz_poly_clear(tpol1); mpz_poly_clear(tpol2);
  return 0;
}

/********************************************************************/
void getTraceConstants(nf_t *N)
/********************************************************************/
{ int      i, k, n=N->degree;
  mpz_poly sk_ah;
  mpz_t    tmp1, tmp2;

  mpz_poly_init(sk_ah);
  mpz_init(tmp1); mpz_init(tmp2);

  /* First, compute the constants in \hat{\alpha}. */
  mpz_set_ui(&sk_ah->coef[0], n);
  for (k=1; k<n; k++) {
    mpz_mul_ui(tmp1, &N->T->coef[n-k], k);
    mpz_neg(&sk_ah->coef[k], tmp1);
    for (i=1; i<k; i++) {
      mpz_mul(tmp1, &N->T->coef[n-i], &sk_ah->coef[k-i]);
      mpz_sub(&sk_ah->coef[k], &sk_ah->coef[k], tmp1);
    }
  }

  /* Bingo. Now compute Tr(omega_i). */
  for (i=0; i<n; i++) {
    mpz_set_ui(tmp1, 0);
    for (k=0; k<=i; k++) {
      mpz_mul(tmp2, &N->W->entry[k][i], &sk_ah->coef[k]);
      mpz_add(tmp1, tmp1, tmp2);
    }
    /* Now, tmp1 better be divisible by W_d! */
    mpz_tdiv_qr(&N->Sk_ib->coef[i], tmp2, tmp1, N->W_d);
    if (mpz_sgn(tmp2)) {
      fprintf(stderr, "getTraceConstants(): Error - ?!?!\n");
      exit(-1);
    }
//    printf("Trace of omega_%d = ", i);
//    mpz_out_str(stdout, 10, &N->Sk_ib->coef[i]);
//    printf("\n");
  }
  N->Sk_ib->degree = n-1; /* Just a formality, in case we'd like to print it maybe. */
  mpz_poly_clear(sk_ah);
  mpz_clear(tmp1); mpz_clear(tmp2);
}





/********************************************************************/
void getTrace_ib(mpz_t t, mpz_poly a, nf_t *N)
/********************************************************************/
/* Input: An algebraic integer 'a' given on the integral basis.     */
/* Output: t <-- Trace(a) in the number field N over Q.             */
/********************************************************************/
{ int   i;
  mpz_t tmp;

  mpz_init(tmp);
  mpz_set_ui(t, 0);
  for (i=0; i<=a->degree; i++) {
    mpz_mul(tmp, &a->coef[i], &N->Sk_ib->coef[i]);
    mpz_add(t, t, tmp);
  }
  mpz_clear(tmp);
}


/********************************************************************/
void norm_std(mpz_t n, mpz_poly a, nf_t *N)
/********************************************************************/
/* Input: An algebraic integer, 'a' given by it's standard          */
/*        representation (a polynomial in \alpha = RootOf(T) ).     */
/*        The number field, N.                                      */
/* Output: n <-- The norm of 'a'.                                   */
/********************************************************************/
{
  mpz_poly_resultant(n, N->T, a);
}

/********************************************************************/
void norm_ib(mpz_t n, mpz_poly a, nf_t *N)
/********************************************************************/
/* Input: An algebraic integer, 'a' given in terms of the integral  */
/*        basis: a = \sum a->coef[i]*\omega_i.                      */
/*        The number field, N.                                      */
/* Output: n <-- The norm of 'a'.                                   */
/********************************************************************/
{ static mpz_poly tpol1;
  static mpz_t    tmp;
  static int      initialized=0;
  int    i, j, d=N->degree;

  if (!(initialized)) {
    mpz_poly_init(tpol1);
    mpz_init(tmp);
    initialized=1;
  }

  /* First, get the standard representation of (W_d^n)*a. */
  for (i=0; i<d; i++)
    mpz_set_ui(&tpol1->coef[i], 0);
  for (i=0; i<d; i++) {
    for (j=0; j<=i; j++) {
      mpz_mul(tmp, &a->coef[i], &N->W->entry[j][i]);
      mpz_add(&tpol1->coef[j], &tpol1->coef[j], tmp);
    }
  }
  tpol1->degree = a->degree;
  if (a->degree > 0)
    mpz_poly_resultant(n, N->T, tpol1);
  else
    mpz_pow_ui(n, &tpol1->coef[0], N->degree);
  mpz_pow_ui(tmp, N->W_d, d);
  mpz_tdiv_qr(n, tmp, n, tmp);
  if (mpz_sgn(tmp))
    printf("norm_ib(): Unexpected error: resultant not divisible by W_d^n!\n");

}

/********************************************************************/
int twoEltRep(mpz_poly a, mpz_mat_t *I, mpz_t p, int f, nf_t *N)
/********************************************************************/
/* Cohen, Alg. 4.7.10.                                              */
/* Input: A prime ideal over p with norm p^f,  given by its basis   */
/*        w.r.t. the omega_i.                                       */               
/*        The number field N with all good fields filled in.        */
/*        That is, I is the HNF of the ideal w.r.t the integral     */
/*        basis.                                                    */
/* Output: 'a' is an algebraic integer, represented as a polynomial */
/*         in the omega_i, such that (p, a) is the two-elt. rep.    */
/*         of the given prime ideal.                                */
/* We assume that I_1 = p.                                          */
/********************************************************************/
{ int       i, j, l, R, lam[MAXPOLYDEGREE], k=N->degree;
  mpz_poly  alpha, tpol2;
  mpz_t     tmp, n, tmp2, G_d;
  mpz_mat_t G;

  if (I->cols == 1) {
    getCol(a, I, 0);
    return 0;
  }

  mpz_poly_init(alpha); mpz_poly_init(tpol2);
  mpz_init(tmp); mpz_init(n); mpz_init(tmp2);
  mpz_init(G_d);
  mpz_mat_init2(&G, k, k);

  mpz_mat_mul(&G, N->W, I);
  mpz_set(G_d, N->W_d);


  fixOrder(&G, G_d);
  /* Now, G is the HNF w.r.t. the standard basis, and it */
  /* has denominator G_d.                                */

#ifdef __HEAVY_DEBUG
  printf("WRT standard basis, I=\n");
  mpz_mat_print(stdout, &G);
  printf("Denominator: "); mpz_out_str(stdout, 10, G_d); printf("\n");
  printf(" f=%d\n", f);
#endif

  R = 1;
ALG_4_7_10_S2:
  for (i=1; i<k; i++)
    lam[i] = R;
  
ALG_4_7_10_S3:
  for (i=0; i<k; i++)
    mpz_set_si(&alpha->coef[i], 0);
  alpha->degree = k-1;
  for (i=1; i<k; i++) {
    for (l=0; l<k; l++)   {
      mpz_mul_si(tmp, &G.entry[l][i], lam[i]);
      mpz_add(&alpha->coef[l], &alpha->coef[l], tmp);
    }
  }
  mpz_poly_fixDeg(alpha);
  mpz_poly_resultant(n, N->T, alpha);

  mpz_pow_ui(tmp, G_d, k);
  mpz_mod(tmp2, n, tmp);
  if (mpz_sgn(tmp2))  {
    printf("twoEltRep(): Something is wrong - subresultant not divisible by d^n.\n");
    printf(" alpha = "); mpz_poly_print(stdout, "", alpha);
    mpz_poly_resultant(n, N->T, alpha);
    printf(" resultant = "); mpz_out_str(stdout, 10, n); printf("\n");
    printf(" G_d^k = "); mpz_out_str(stdout, 10, tmp); printf("\n\n");
  } else
    mpz_div(n, n, tmp); 
#ifdef __HEAVY_DEBUG
  printf("N(alpha) = "); mpz_out_str(stdout, 10, n); printf("\n");
#endif

  mpz_pow_ui(tmp, p, f);
  mpz_div(n, n, tmp);
  mpz_mod(tmp2, n, p);
  if (mpz_sgn(tmp2)) {
    /* p doesn't divide n anymore. */
    mpz_poly_cp(a, alpha);
    goto ALG_4_7_10_DONE;
  }

  mpz_mul(tmp, p, G_d);
  mpz_add(&alpha->coef[0], &alpha->coef[0], tmp);
  mpz_poly_resultant(n, N->T, alpha);
  mpz_sub(&alpha->coef[0], &alpha->coef[0], tmp);

  /* This needs to be double checked! */
  mpz_pow_ui(tmp, G_d, k);
  mpz_mod(tmp2, n, tmp);
  if (mpz_sgn(tmp2))  {
    printf("twoEltRep(): Something is wrong - subresultant not divisible by d^n.\n");
    printf("  alpha = "); mpz_poly_print(stdout, "", alpha); 
    mpz_poly_resultant(n, N->T, alpha);
    printf(" resultant = "); mpz_out_str(stdout, 10, n); printf("\n");
    printf("  G_d^k = "); mpz_out_str(stdout, 10, tmp); printf("\n");
  } else
    mpz_div(n, n, tmp); 
#ifdef __HEAVY_DEBUG
  printf("N(alpha) = "); mpz_out_str(stdout, 10, n); printf("\n");
#endif

  mpz_pow_ui(tmp, p, f);
  mpz_div(n, n, tmp);
  mpz_mod(tmp2, n, p);
  if (mpz_sgn(tmp2)) {
    /* p doesn't divide n anymore. */
    mpz_poly_cp(a, alpha);
    goto ALG_4_7_10_DONE;
  }

  /* Step 4: */
  j=k-1;
  while ((j>0) && (lam[j]==-R))
    j--;

  lam[j] -= 1;
  for (i=j+1; i<k; i++)
    lam[i] = R;

  /* Step 5: */
  j=1;
  while ((j<k) && (lam[j]==0))
    j++;
  if (j >= k) {
    R++;
    goto ALG_4_7_10_S2;
  }
  goto ALG_4_7_10_S3;


ALG_4_7_10_DONE:
  /*****************************************************************/
  /* We're almost there. What we have is that alpha = a(theta)/G_d */
  /* does the job. Since this is an algebraic integer, we can      */
  /* just convert a(theta) to its representation w.r.t. the        */
  /* integral basis, and divide the coefficients by G_d.           */
  /*****************************************************************/
  stdtoib(a, a, N);
  for (i=0; i<=a->degree; i++) {
    mpz_mod(tmp, &a->coef[i], G_d);
    if (mpz_sgn(tmp)) {
      printf("twoEltRep() Error - inconsistent result!\n");
    }
    mpz_div(&a->coef[i], &a->coef[i], G_d);
  }
  

  mpz_poly_clear(alpha); mpz_poly_clear(tpol2);
  mpz_clear(tmp); mpz_clear(n); mpz_clear(tmp2); mpz_clear(G_d);
  mpz_mat_clear(&G);
  return 0;
}


/****************************************************************/
void getIdealHNF_ib(mpz_mat_t *H, mpz_t p, mpz_poly a, s32 e, nf_t *N)
/****************************************************************/
/* Input: A prime, p, and an algebraic integer 'a' given as a   */
/*        polynomial on the integral basis in 'N'. An exponent  */
/*        e.                                                    */
/* Output: The HNF representation of the ideal (p, a)^e w.r.t   */
/*        the integral basis.                                   */
/****************************************************************/
/* Check me!!!!
*/
{ int       i, j, n=N->degree;
  s32      tmpE;
  static mpz_mat_t M, Mpow, Me;
  static mpz_poly  tpol1, tpol2;
  static int initialized=0;

  if (!(initialized)) {
    mpz_mat_init2(&M, n, 2*n); mpz_mat_init2(&Mpow, n, n); 
    mpz_mat_init2(&Me, n, n);
    mpz_poly_init(tpol1); mpz_poly_init(tpol2);
    initialized=1;
  }

  M.rows = n; M.cols = 2*n;
  /* Compute the p*omega_i columns. */
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      if (i!=j)
        mpz_set_si(&M.entry[i][n+j], 0);
      else
        mpz_set(&M.entry[i][n+j], p);
  /* Compute the a*omega_i columns. */
  for (i=0; i<n; i++) {
    tpol1->degree=i;
    if (i>0)
      mpz_set_si(&tpol1->coef[i-1], 0);
    mpz_set_si(&tpol1->coef[i], 1);
    basisPolyMult(tpol2, tpol1, a, N->Mt, n);
    setCol(&M, i, tpol2);
  }
  if (e==1) {
    /* Whew - that was easy :) */
    mpz_mat_getHNF(H, &M);
    return;
  }
  
  mpz_mat_setID(&Me, n);
  mpz_mat_getHNF(&M, &M);
  mpz_mat_cp(&Mpow, &M);

  tmpE = e;
  /* Now we have Me = product so far,  */
  /*           Mpow = M^{2^j}.         */
  while (tmpE) {
    if (tmpE&0x01) 
      idealHNFMul_ib(&Me, &Me, &Mpow, N);
    idealHNFMul_ib(&Mpow, &Mpow, &Mpow, N);
    tmpE = (tmpE>>1);
  } 
  mpz_mat_cp(H, &Me);
 
}

/****************************************************************/
void idealHNFMul_ib(mpz_mat_t *A, mpz_mat_t *B, mpz_mat_t *C, nf_t *N)
/****************************************************************/
/* B and C are ideals, represented wrt the integral basis in N. */
/* Do A <-- B*C.                                                */
/****************************************************************/
{ int       i, j, n=N->degree;
  static mpz_poly  tpol1, bi, cj;
  static mpz_mat_t Big;
  static int initialized=0;

  if (!(initialized)) {
    mpz_mat_init2(&Big, n, n*n);
    mpz_poly_init(tpol1); mpz_poly_init(bi); mpz_poly_init(cj);
    initialized=1;
  }
  for (i=0; i<n; i++) {
    getCol(bi, B, i);
    for (j=0; j<n; j++) {
      getCol(cj, C, j);
      basisPolyMult(tpol1, bi, cj, N->Mt, n);
      setCol(&Big, (i*n+j), tpol1);
    }
  }
  Big.rows = n;
  Big.cols = n*n;
  mpz_mat_getHNF(A, &Big);
}


/****************************************************************/
void idealHNFPow_ib(mpz_mat_t *A, mpz_mat_t *B, s32 e, nf_t *N)
/****************************************************************/
/* B and C are ideals, represented wrt the integral basis in N. */
/* Do A <-- B^e.                                                */
/****************************************************************/
{ int    n=N->degree;
  s32   remain;
  static mpz_mat_t Bpow2k, Tmp;
  static int initialized=0;

  if (!(initialized)) {
    mpz_mat_init2(&Bpow2k, n, n);
    mpz_mat_init2(&Tmp, n, n);
    initialized=1;
  }

  mpz_mat_setID(&Tmp, n);
  mpz_mat_cp(&Bpow2k, B);
  remain = e;
  while (remain) {
    if (remain%2) {
      idealHNFMul_ib(&Tmp, &Tmp, &Bpow2k, N);
    }
    remain = (remain>>1);
    idealHNFMul_ib(&Bpow2k, &Bpow2k, &Bpow2k, N);
  }

  mpz_mat_cp(A, &Tmp);
}

/****************************************************************/
int cmpIdeals(const void *A, const void *B)
/****************************************************************/
{ prime_id_t *a=(prime_id_t *)A, *b=(prime_id_t *)B;
  int res, i;

  res = mpz_cmp(a->p, b->p);
  if (res) return res;
  if (a->alpha->degree < b->alpha->degree) return -1;
  if (a->alpha->degree > b->alpha->degree) return 1;
  for (i=0; i<=a->alpha->degree; i++) {
    res = mpz_cmp(&a->alpha->coef[i], &b->alpha->coef[i]);
    if (res) return res;
  }
  return 0;
}
  

/****************************************************************/
void idealSort(prime_id_t *I, int Isize)
/****************************************************************/
{ 
  qsort(I, Isize, sizeof(prime_id_t), cmpIdeals);
}


/****************************************************************/
int factorPrime(prime_id_t *I, mpz_t p, nf_t *N)
/****************************************************************/
/* Don't look at this function! Between the really involved math */
/* and the really terrible code, it's quite possible that your   */
/* head might actually implode.                                  */
/****************************************************************/
/* Input: I must be allocated already.                          */
/*        'p' is a prime integer.                               */
/*        'N' is a number field, as obtained from               */
/*            getIntegralBasis().                               */
/* Output: 'I' will be a list of prime ideals and exponents     */
/*         so that <p> = I[0]^e[0]...I[k]^e[k].                 */
/* Return value: The number of distinct prime ideals, k.        */
/****************************************************************/
/* Make sure all these mpz_??? things get cleared!!!!
   And fix me!!!!
   There is a serious problem somewhere, and I think it's likely
   in here. Without the bizarre fprintf()'s below (or some other
   suitable call to a printf or fprintf function), we get a 
   seg fault. That suggests that we might be overwriting some
   important part of the stack or something. But I haven't yet
   figured it out.
*/
{ int       n = N->degree, c=0, Isize=0;
  int       *ej, i=0, j=0, k=0, g=0, f=0, retVal=0, r=0;
  int       index1=0, index2=0, *m_e, s=0;
  mpz_poly  T, gam_i, gam_j, tpol1, tpol2, tpol3, m, alpha;
  mpz_poly  *Tj, *m_i, *Gam_ij;
  mpz_t     t1, Ip_d;
  mpz_mat_t Ip, *L, Beta, Gamma;
  mpz_mat_t T1, T2, M, M1, H;


  ej = (int *)malloc(n*sizeof(int));
  m_e = (int *)malloc(n*sizeof(int));
  mpz_poly_init(T);
  mpz_poly_init(gam_i);  mpz_poly_init(gam_j);
  mpz_poly_init(tpol1);  mpz_poly_init(tpol2); mpz_poly_init(tpol3);
  mpz_poly_init(m);
  mpz_poly_init(alpha);
  Tj = (mpz_poly *)malloc(n*sizeof(mpz_poly));
  m_i = (mpz_poly *)malloc(n*sizeof(mpz_poly));
  Gam_ij = (mpz_poly *)malloc(n*n*sizeof(mpz_poly));
  L = (mpz_mat_t *)malloc(n*sizeof(mpz_mat_t));
  for (i=0; i<n; i++) {
    mpz_poly_init(Tj[i]);
    mpz_poly_init(m_i[i]);
    mpz_mat_init2(&L[i], n, n);
  }
  for (j=0; j<n*n; j++) 
    mpz_poly_init(Gam_ij[j]);
  mpz_init(t1); mpz_init(Ip_d);
  mpz_mat_init2(&Ip, n, n);
  mpz_mat_init2(&Beta, n, n+1);
  mpz_mat_init2(&Gamma, n, n);
  mpz_mat_init2(&T1, 2*n, 2*n);
  mpz_mat_init2(&T2, 2*n, 2*n);
  mpz_mat_init2(&M, n, 2*n);
  mpz_mat_init2(&M1, n, n);
  mpz_mat_init2(&H, n, 2*n);


  mpz_mod(t1, N->index, p);
  if (mpz_sgn(t1)) {
    /* Since p does not divide the index, life is easy! */
    T->degree = 0;
    for (i=0; i<=n; i++) {
      mpz_mod(&T->coef[i], &N->T->coef[i], p);
      if (mpz_sgn(&T->coef[i]))
        T->degree = i;
    }
    g = mpz_poly_fact(Tj, ej, T, p);
    for (i=0; i<g; i++) {
      mpz_set(I[i].p, p);
      mpz_poly_cp(I[i].alpha, Tj[i]);
      computeVconst(I[i].beta, p, I[i].alpha, N);
      I[i].e = ej[i];
      I[i].f = Tj[i]->degree;
    }
    Isize = g;
    goto ALG_6_2_9_DONE;
  }  

  /* Step 2: */
  computePRadical2(&Ip, N->Mt, N->T, p);
  if (Ip.rows == 0)
    Ip.rows = n; /* Special case, where Ip = p*Op. */
  c = 1;
  mpz_mat_cp(&L[0], &Ip);

  /* Step 8: */
ALG_6_2_9_S8:
  mpz_mat_cp(&Beta, &L[c-1]); /* Always choose the last one - later we'll overwrite it! */
  mpz_mat_cp(&H, &Beta); /* for later use. */
  r = Beta.cols;
  Beta.rows = n; /* If L[0] is the zero ideal, L[0] could be 0x0 (a bug). */
  mpz_set_si(&Beta.entry[0][r], 1);
  for (index1=1; index1<n; index1++)
    mpz_set_si(&Beta.entry[index1][r], 0);
  Beta.cols = r+1;
  if (r+1 < n) /* Bug in suppBasis when Beta is nxn. */
    mpz_mat_suppBasis_modp(&Beta, &Beta, p);
  /* The F_p-basis of A is now given by Beta_{r},...,Beta_{n-1}. */

  /* Step 9: */
  f = n-r;
//  printf("f=%d\n", f);
  /* Note: This could be made much simpler, but I want
     to follow Cohen as much as possible for now.
  */
  /* We are going to compute the mult table of the gamma_i. */
  for (index1=0; index1<f; index1++) 
    for (index2=0; index2<n; index2++) 
      mpz_set(&Gamma.entry[index2][index1], &Beta.entry[index2][r+index1]);
  for (index1=f; index1<n; index1++) 
    for (index2=0; index2<n; index2++) 
      mpz_set(&Gamma.entry[index2][index1], &Beta.entry[index2][index1-f]);
  Gamma.rows = n; Gamma.cols = n;

  for (index1=0; index1<f; index1++) {
    getCol(gam_i, &Gamma, index1);
    for (index2=0; index2<f; index2++) {
      getCol(gam_j, &Gamma, index2);
      basisPolyMult(tpol1, gam_i, gam_j, N->Mt, n);
      T1.rows=n; T1.cols=1;
      setCol(&T1, 0, tpol1);
      if (mpz_mat_invIm(&T2, &T1, &Gamma, p)) {
        printf("Some error occurred computing mult. table of the gamma_i!\n");
        retVal = -1;
        goto ALG_6_2_9_DONE;
      }
      /* We now have that T1 = Gam*T2. */
      getCol(Gam_ij[index1*f + index2], &T2, 0);
      /* But ignore the Beta coefficients. */
      Gam_ij[index1*f + index2]->degree = f-1;
      mpz_poly_fixDeg(Gam_ij[index1*f + index2]);
    }
  }

  /* Step 10: */
  M.rows = f; M.cols = f;
  for (index1=0; index1<f; index1++) {
    /* tpol1 <-- gamma_{index1}. */
    for (index2=0; index2<f; index2++) 
      mpz_set_si(&tpol1->coef[index2], 0);
    mpz_set_si(&tpol1->coef[index1], 1);
    tpol1->degree = index1;
    
    /* tpol2 <-- tpol1^{p} */
    basisPolyPowM(tpol2, tpol1, p, Gam_ij, p, f);

    /* tpol2 <-- tpol2 - tpol1. */
    for (index2=tpol2->degree+1; index2<=tpol1->degree; index2++) {
      mpz_set_si(&tpol2->coef[index2], 0);
      tpol2->degree = index2;
    }
    mpz_sub_ui(&tpol2->coef[index1], &tpol2->coef[index1], 1);
    mpz_mod(&tpol2->coef[index1], &tpol2->coef[index1], p);
    mpz_poly_fixDeg(tpol2);
    setCol(&M, index1, tpol2);
  }
  mpz_mat_getKernel_modp(&M1, &M, p);
//printf("Step 10: M1 = \n"); mpz_mat_print(stdout, &M1);

  /* Step 11: */
  if (M1.cols > 1)
    goto ALG_6_2_9_S12;
  mpz_set(I[Isize].p, p);
  for (index1=H.cols; index1>0; index1--)
    for (index2=0; index2<n; index2++)
      mpz_set(&H.entry[index2][index1], &H.entry[index2][index1-1]);
  mpz_set(&H.entry[0][0], p);
  for (index2=1; index2<n; index2++)
    mpz_set_ui(&H.entry[index2][0], 0);
  H.cols += 1;
//  mpz_mat_getHNF(&H, &H);
//printf("Calling twoEltRep() on H=\n");
//mpz_mat_print(stdout, &H);

  twoEltRep(tpol1, &H, p, f, N);
  mpz_set(I[Isize].p, p);
  mpz_poly_cp(I[Isize].alpha, tpol1);
  computeVconst(I[Isize].beta, p, I[Isize].alpha, N);
  mpz_mat_setID(&T1, n);
  for (index1=0; index1<n; index1++)
    mpz_set(&T1.entry[index1][index1], p);
  I[Isize].e = valuation(&T1, &I[Isize], N);
  if (I[Isize].e<=0) {
    fprintf(stderr, "Error: factorPrime() is broken!\n");
    exit(-1);
  }
  I[Isize].f = f;
  Isize++;
  c--;
  if (c > 0)
    goto ALG_6_2_9_S8;
 goto ALG_6_2_9_DONE; 



ALG_6_2_9_S12:
  /* Step 12: */
  alpha->degree = 0;
  for (index1=0; (index1 < M1.cols)&&(alpha->degree==0); index1++)
    for (index2=1; index2 < M1.rows; index2++)
      if (mpz_sgn(&M1.entry[index2][index1])) {
        getCol(alpha, &M1, index1);
        break;
      }
  /* Now, we need the minimal polynomial m(X) of alpha. */
  mpz_set_ui(&tpol2->coef[0], 1);
  tpol2->degree = 0;
  T1.rows = f;
  T1.cols = 0;
  index1=0;
  do {
    T1.cols += 1;
    setCol(&T1, index1, tpol2);
    basisPolyMult(tpol2, tpol2, alpha, Gam_ij, f);
    mpz_mat_getKernel_modp(&T2, &T1, p);
    index1++;
  } while (T2.cols <= 0);
  getCol(m, &T2, 0);
  /* And make sure the leading coefficient is one: */
  if (mpz_invert(t1, &m->coef[m->degree], p)) {
    for (index1=0; index1<=m->degree; index1++) {
      mpz_mul(&m->coef[index1], &m->coef[index1], t1);
      mpz_mod(&m->coef[index1], &m->coef[index1], p);
    }
  }
  mpz_poly_fixDeg(m);

  /* There is either something wrong with this function (the
     stack being illegally overwritten) or a bug in gcc 3.2.2.
     Whatever it is/was, a printf/fprintf call seemed to correct
     for it (I guess, by magically re-aligning the stack pointer
     or something). Anyway, it was totally not clear what was
     going on. I tried very hard to debug it with ddd/gdb, and
     various memory-management tools, but to no avail. Anyway,
     it seems ok now. But if something happens and the code acts
     in some very unexpected way, like an infinite loop, try adding
     this statement back in.
  */
//fprintf(stderr, "111.\n");

  /* Step 13: */
  k = mpz_poly_fact(m_i, m_e, m, p);
  /* Step 14: */

  r = H.cols;
  for (s=0; s<k; s++) {
    mpz_mat_cp(&M, &H);
    M.rows = n;
    /* We need to compute m_i(alpha), expressed on the integral basis. */
    /* But remember : alpha is expressed on the gamma basis!           */
    /* And m_i(s) is linear, so this is easy.                          */
    for (index1=0; index1<n; index1++)
      mpz_set_ui(&tpol1->coef[index1], 0);

    for (index1=0; index1<=alpha->degree; index1++) {
      getCol(gam_j, &Gamma, index1);
      for (index2=0; index2<=gam_j->degree; index2++) {
        mpz_mul(t1, &gam_j->coef[index2], &alpha->coef[index1]);
        mpz_add(&tpol1->coef[index2], &tpol1->coef[index2], t1);
      }
    }
    for (i=0; i<=tpol1->degree; i++)
      mpz_mul(&tpol1->coef[i], &tpol1->coef[i], &m_i[s]->coef[1]);
    mpz_add(&tpol1->coef[0], &tpol1->coef[0], &m_i[s]->coef[0]);
    tpol1->degree = n-1;
    mpz_poly_fixDeg(tpol1);

    for (index1=0; index1<n; index1++) {
      /* tpol2 <-- omega_{index1} */
      for (index2=0; index2<index1; index2++)
        mpz_set_ui(&tpol2->coef[index2], 0);
      mpz_set_ui(&tpol2->coef[index1], 1);
      tpol2->degree = index1;
      /* tpol3 <-- tpol1 * tpol2, expressed on the omega_i. */
      basisPolyMult(tpol3, tpol1, tpol2, N->Mt, n);
      setCol(&M, M.cols, tpol3);
      M.cols += 1;
    }
    /* Step 15, integrated. */
    mpz_mat_getImage_modp(&L[c-1+s], &M, p);
  }
  c = c+k-1;
  goto ALG_6_2_9_S8;


ALG_6_2_9_DONE:
  /* Sort the special ideals so they always appear in the same order. */
  if (Isize > 1)
    idealSort(I, Isize);

  free(ej); free(m_e);
  mpz_poly_clear(T);
  mpz_poly_clear(m);
  mpz_poly_clear(gam_i); mpz_poly_clear(gam_j);
  mpz_poly_clear(tpol1); mpz_poly_clear(tpol2); mpz_poly_clear(tpol3);
  mpz_poly_clear(alpha);
  for (i=0; i<N->degree; i++)
    mpz_poly_clear(Tj[i]);
  mpz_clear(t1); mpz_clear(Ip_d);
  mpz_mat_clear(&Ip);
  mpz_mat_clear(&Beta);
  mpz_mat_clear(&Gamma);
  mpz_mat_clear(&T1);
  mpz_mat_clear(&T2);
  mpz_mat_clear(&M);
  mpz_mat_clear(&M1);
  mpz_mat_clear(&H);
  for (i=0; i<n; i++) {
    mpz_mat_clear(&L[i]);
    mpz_poly_clear(m_i[i]);
    for (j=0; j<n; j++)
      mpz_poly_clear(Gam_ij[i*n + j]);
  }
  free(Tj); free(m_i); free(Gam_ij); free(L);
  if (retVal < 0)
    return retVal;
  return Isize;
}


/****************************************************************/
int computeVconst(mpz_poly Beta, mpz_t p, mpz_poly alpha, nf_t *N)
/****************************************************************/
/* As in Cohen Alg. 4.8.17, steps 1-2, compute the value Beta   */
/* for the ideal <p, alpha> to help in computing valuations at  */
/* <p, alpha>.                                                  */
/* alpha and beta are given w.r.t. the integral basis in N.     */
/****************************************************************/
{ int       i, k, n=N->degree, retVal=0;
  mpz_mat_t A, B;
  mpz_poly  w_i, g_j, tpol1;

  mpz_mat_init2(&A, 2*n, n);
  mpz_mat_init2(&B, n, n);
  mpz_poly_init(w_i); mpz_poly_init(g_j);
  mpz_poly_init(tpol1);


  for (i=0; i<n; i++) {
    if (i>0)
      mpz_set_ui(&w_i->coef[i-1], 0);
    mpz_set_ui(&w_i->coef[i], 1); 
    w_i->degree = i;
    
    g_j->degree = 0;
    mpz_set(&g_j->coef[0], p);
    
    basisPolyMult(tpol1, w_i, g_j, N->Mt, n);
    for (k=0; k<n; k++) {
      if (k <= tpol1->degree)
        mpz_set(&A.entry[k][i], &tpol1->coef[k]);
      else
        mpz_set_ui(&A.entry[k][i], 0);
    }

    mpz_poly_cp(g_j, alpha);
 
    basisPolyMult(tpol1, w_i, g_j, N->Mt, n);
    for (k=0; k<n; k++) {
      if (k <= tpol1->degree)
        mpz_set(&A.entry[n+k][i], &tpol1->coef[k]);
      else
        mpz_set_ui(&A.entry[k][i], 0);
    }
  }
  A.cols = n; A.rows = 2*n;
  mpz_mat_getKernel_modp(&B, &A, p);
  if (B.cols > 0)
    getCol(Beta, &B, 0);
  else {
    printf("computeVconst(): Error - kernel is zero!\n");
    retVal = -1;
  }

  mpz_mat_clear(&A);
  mpz_mat_clear(&B);
  mpz_poly_clear(w_i); mpz_poly_clear(g_j);
  mpz_poly_clear(tpol1);
  return retVal;
}


/**********************************************************************/
int valuation(mpz_mat_t *_A, prime_id_t *I, nf_t *N)
/*********************************************************************/
/* Input: A prime, I->p, and an algebraic integer I->alpha given as  */
/*        a polynomial in the omega_i, so that <I->p, I->alpha> is a */
/*        prime ideal. I->Beta is the structure constant already     */
/*        computed by computeVconst() for this ideal. _A is an       */
/*        HNF form of an ideal w.r.t. the N->W basis.                */
/* Return value: The valuation of the ideal defined by _A at the     */
/*        prime ideal <p, alpha>.                                    */
/*        negative, on error (_A is assumed an integral ideal).      */
/*********************************************************************/
{ int    i, j, v=-1, n=N->degree;
  static mpz_mat_t A;
  static mpz_poly  tpol1, tpol2;
  static mpz_t     norm, tmp;
  static int initialized=0;

  if (!(initialized)) {
    mpz_mat_init2(&A, n, n);
    mpz_poly_init(tpol1); mpz_poly_init(tpol2);
    mpz_init(norm); mpz_init(tmp);
    initialized=1;
  }

  mpz_set_ui(norm, 1);
  for (i=0; i<n; i++) {
    mpz_mul(norm, norm, &_A->entry[i][i]);
    mpz_mod(norm, norm, I->p);
  }
  v = 0;
  if (mpz_sgn(norm)) 
    goto ALG_4_8_17_DONE;
  mpz_mat_cp(&A, _A);

  /* Step 4: */
ALG_4_8_17_S4:
  for (i=0; i<n; i++) {
    getCol(tpol1, &A, i);
    basisPolyMult(tpol2, tpol1, I->beta, N->Mt, n);
    setCol(&A, i, tpol2);
  }

//  mpz_mat_getHNF(&A, &A);
  /* Step 5: */
  mpz_mod(tmp, &A.entry[n-1][n-1], I->p);
  if (mpz_sgn(tmp))
    goto ALG_4_8_17_DONE;
  mpz_mod(tmp, N->index, I->p);
  if (mpz_sgn(tmp)) {
    v++;
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
        mpz_tdiv_qr(&A.entry[i][j], tmp, &A.entry[i][j], I->p);
        if (mpz_sgn(tmp)) {
          printf("valuation() Error : 'A' is supposed to be divisible by p, but is not!\n");
          v=-1; goto ALG_4_8_17_DONE;
        }
      }
    }
    goto ALG_4_8_17_S4;
  }

  /* Step 6: */
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      mpz_tdiv_qr(&A.entry[i][j], tmp, &A.entry[i][j], I->p);
      if (mpz_sgn(tmp)) {
        goto ALG_4_8_17_DONE;
      }
    }
  }
  v++;
  goto ALG_4_8_17_S4;

ALG_4_8_17_DONE:

  return v;
}
  
/**********************************************************************/
int valuation2(mpz_mat_t *_A, prime_id_t *I, nf_t *N)
/*********************************************************************/
/* Same as above, but a little faster. The caveat is that this one   */
/* apparently cannot be used from within the `factorPrime' function -*/
/* at that point, something hasn't been initialized yet or something */
/* which causes a failure. This one is intended for use in the       */
/* relation factorization code, so this is the one that should be    */
/* optimized.                                                        */
/*********************************************************************/
{ int    i, j, l, v=-1, n=N->degree;
  static mpz_mat_t A;
  static mpz_poly  tpol1, tpol2;
  static mpz_t     norm, tmp;
  static int initialized=0;

  if (!(initialized)) {
    mpz_mat_init2(&A, n, n);
    mpz_poly_init(tpol1); mpz_poly_init(tpol2);
    mpz_init(norm); mpz_init(tmp);
    initialized=1;
  }

  mpz_set_ui(norm, 1);
  for (i=0; i<n; i++) {
    mpz_mul(norm, norm, &_A->entry[i][i]);
    mpz_mod(norm, norm, I->p);
  }
  v = 0;
  if (mpz_sgn(norm)) 
    goto ALG_4_8_17_DONE;
  mpz_mat_cp(&A, _A);


  /* Step 4: */
ALG_4_8_17_S4:
  for (i=0; i<n; i++) {
    getCol(tpol1, &A, i);
    /* Do tpol2 <-- tpol1*beta, using betaMat to shortcut the computation. */
    for (j=0; j<n; j++)
      mpz_set_ui(&tpol2->coef[j], 0);
    for (j=0; j<=tpol1->degree; j++) {
      /* tpol2 <-- tpol2 + tpol1[j]*betaMat_col_j. */
      for (l=0; l<n; l++) {
        mpz_mul(tmp, &tpol1->coef[j], &I->betaMat.entry[l][j]);
        mpz_add(&tpol2->coef[l], &tpol2->coef[l], tmp);
      }
    }
    tpol2->degree = n-1;
    setCol(&A, i, tpol2);
  } /* end bottleneck. */

/* This should be necessary, but things have been okay without it. Why?
   Does it just so happen that the matrix is in HNF because of the way we are computing it?
*/
#if 0
   mpz_mat_getHNF(&A, &A);
#endif
  /* Step 5: */
  mpz_mod(tmp, &A.entry[n-1][n-1], I->p);
  if (mpz_sgn(tmp))
    goto ALG_4_8_17_DONE;
  mpz_mod(tmp, N->index, I->p);
  if (mpz_sgn(tmp)) {
    v++;
    for (i=0; i<n; i++) {
      for (j=0; j<n; j++) {
        mpz_tdiv_qr(&A.entry[i][j], tmp, &A.entry[i][j], I->p);
        if (mpz_sgn(tmp)) {
          printf("valuation() Error : 'A' is supposed to be divisible by p, but is not!\n");
          v=-1; goto ALG_4_8_17_DONE;
        }
      }
    }
    goto ALG_4_8_17_S4;
  }

  /* Step 6: */
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      mpz_tdiv_qr(&A.entry[i][j], tmp, &A.entry[i][j], I->p);
      if (mpz_sgn(tmp)) {
        goto ALG_4_8_17_DONE;
      }
    }
  }
  v++;
  goto ALG_4_8_17_S4;

ALG_4_8_17_DONE:

  return v;
}

/**********************************************************************/
int valuation_ab(mpz_t a, mpz_t b, int k, nf_t *N)
/*********************************************************************/
/* Compute the valuation of <a-b\omega_1> at special ideal #k.       */
/* Return value: the valuation ( >=0 ) on success, -255 on error.    */
/*********************************************************************/
{ static mpz_t q, r, tmp;
  static mpz_mat_t H;
  static int initialized=0;
  int    e, divisible;

  if (!(initialized)) {
    mpz_init(q); mpz_init(r); mpz_init(tmp);
    mpz_mat_init(&H);
    initialized=1;
  }

  e=0;
  do {
    divisible=0;
    /* Is (a,b)^T in the column space of N->sPowMats[k][e] ? */
    mpz_tdiv_qr(q, r, b, &N->sPowMats[k][e].entry[1][1]);
    if (mpz_sgn(r)==0) {
      mpz_neg(q, q); /* Since it's a <a MINUS b\alpha>. */
      mpz_mul(tmp, q, &N->sPowMats[k][e].entry[0][1]);
      mpz_sub(tmp, a, tmp);
      mpz_mod(tmp, tmp, &N->sPowMats[k][e].entry[0][0]);
      if (mpz_sgn(tmp)==0) {
        divisible=1;
        e++;
      }
    }
  } while (divisible && (e<MAX_SP_POWERS));
  if (e >= MAX_SP_POWERS) {
    return -255;
  }
  return e;
}
