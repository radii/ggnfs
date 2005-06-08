/**************************************************************/
/* mpz_mat.c                                                  */
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

#include <stdio.h>
#include <stdlib.h>
#include "ggnfs.h"



/* Support functions. */
int  redI(int k, int l, mpz_mat_t *B, mpz_mat_t *H);
int  swapI(int k, int k_max, mpz_mat_t *B, mpz_mat_t *H);
void mpz_nearest_int(mpz_t q, mpz_t num, mpz_t den);
int  swapK(int k, int k_max, int *f, mpz_mat_t *B, mpz_mat_t *H);


/* Globals */
static mpz_mat_t lam, d;
static mpz_t     u, q, t, tmp1, tmp2, btmp, ltmp;
static int       mpz_mat_initialized=0;


/************************************************************/
void initGlobals()
{
  if (mpz_mat_initialized)
    return;
  mpz_mat_init(&lam); mpz_mat_init(&d);
  mpz_init(u); mpz_init(q); mpz_init(t);
  mpz_init(tmp1); mpz_init(tmp2);
  mpz_init(btmp); mpz_init(ltmp);
  mpz_mat_initialized = 1;
}


/************************************************************/
void mpz_mat_init(mpz_mat_t *M)
{ int i, j;

  M->maxRows = DEFAULT_MPZ_MAT_ROWS;
  M->maxCols = DEFAULT_MPZ_MAT_COLS;
  M->entry = (__mpz_struct **)malloc(M->maxRows*sizeof(__mpz_struct *));
  for (i=0; i<M->maxRows; i++) {
    M->entry[i] = (__mpz_struct *)malloc(M->maxCols*sizeof(__mpz_struct));
    for (j=0; j<M->maxCols; j++) 
      mpz_init(&M->entry[i][j]);
  }
  M->rows = M->cols = 0;
}
/************************************************************/
void mpz_mat_init2(mpz_mat_t *M, int maxRows, int maxCols)
{ int i, j;

  M->maxRows = MAX(maxRows, 1);
  M->maxCols = MAX(maxCols, 1);
  M->entry = (__mpz_struct **)malloc(M->maxRows*sizeof(__mpz_struct *));
  for (i=0; i<M->maxRows; i++) {
    M->entry[i] = (__mpz_struct *)malloc(M->maxCols*sizeof(__mpz_struct));
    for (j=0; j<M->maxCols; j++) 
      mpz_init(&M->entry[i][j]);
  }
  M->rows = M->cols = 0;
}
/************************************************************/
void mpz_mat_clear(mpz_mat_t *M)
{ int i, j;

  for (i=0; i<M->maxRows; i++) {
    for (j=0; j<M->maxCols; j++) 
      mpz_clear(&M->entry[i][j]);
    free(M->entry[i]);
  }
  free(M->entry);
  M->rows = M->cols = 0;
  M->maxRows = M->maxCols = 0;
}

/************************************************************/
void mpz_mat_print(FILE *fp, mpz_mat_t *M)
{ int i, j;

  if (M->rows == 0 || M->cols == 0) {
    fprintf(fp, "[]\n");
    return;
  }
  for (i=0; i<M->rows; i++) {
    for (j=0; j<M->cols; j++) {
      mpz_out_str(fp, 10, &M->entry[i][j]);
      fprintf(fp, "  ");
    }
    fprintf(fp, "\n");
  }
}

/***************************************************************/
int mpz_mat_equal(mpz_mat_t *A, mpz_mat_t *B)
{ int i, j;

  if ((A->rows != B->rows) || (A->cols != B->cols))
    return 0;
  for (i=0; i< A->rows; i++)
    for (j=0; j<A->cols; j++)
      if (mpz_cmp(&A->entry[i][j], &B->entry[i][j]))
        return 0;
  return 1;
}

/***************************************************************/
int mpz_mat_iszero(mpz_mat_t *M)
{ int i, j;

  for (i=0; i<M->rows; i++)
    for (j=0; j<M->cols; j++)
      if (mpz_sgn(&M->entry[i][j]))
        return 0;
  return 1;
}

/***************************************************************/
int mpz_mat_isID(mpz_mat_t *M)
{ int i, j;

  for (i=0; i<M->rows; i++)
    for (j=0; j<M->cols; j++)
      if (mpz_cmp_ui(&M->entry[i][j], (i==j)) )
        return 0;
  return 1;
}


/***************************************************************/
void mpz_mat_setID(mpz_mat_t *M, int n)
{ int i, j;

  M->rows = M->cols = n;
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      mpz_set_ui(&M->entry[i][j], (i==j));
}

/****************************************************************/
void getCol(mpz_poly h, mpz_mat_t *M, int col)
/****************************************************************/
{ int i;
                                                                                                      
  for (i=M->rows-1; i>=0; i--)
    mpz_set(&h->coef[i], &M->entry[i][col]);
  h->degree = M->rows-1;
  mpz_poly_fixDeg(h);
}
                                                                                                      
/****************************************************************/
void setCol(mpz_mat_t *M, int col, mpz_poly h)
/****************************************************************/
{ int i;
                                                                                                      
  for (i=0; i<=h->degree; i++)
    mpz_set(&M->entry[i][col], &h->coef[i]);
  for ( ; i<M->rows; i++)
    mpz_set_ui(&M->entry[i][col], 0);
}

/***************************************************************/
void mpz_mat_dotCols(mpz_t res, mpz_mat_t *M, int col1, int col2)
{ int i;
  mpz_t tmp;
  
  mpz_init(tmp);
  mpz_set_ui(res, 0);
  for (i=0; i<M->rows; i++) {
    mpz_mul(tmp, &M->entry[i][col1], &M->entry[i][col2]);
    mpz_add(res, res, tmp);
  }
  mpz_clear(tmp);
}

/***************************************************************/
int mpz_mat_LLL(mpz_mat_t *B, mpz_mat_t *H, mpz_mat_t *A)
/***************************************************************/
/* Input: 'A' is a matrix of column vectors to be LLL-reduced. */
/* Output: 'B' is a matrix whose column vectors are an         */
/*          LLL-reduced basis for colsp(A).                    */
/*         'H' is the transition matrix so that B=A*H.         */
/***************************************************************/
{ int       i, j, k, k_max, l, res=0, n=A->cols;
  
  initGlobals();
  
  B->rows = A->rows;
  B->cols = A->cols;
  for (i=0; i<B->rows; i++)
    for (j=0; j<B->cols; j++)
      mpz_set(&B->entry[i][j], &A->entry[i][j]);

  
  /* (1) Initialize. */
  k=2; k_max=1;
  mpz_set_ui(&d.entry[0][0], 1);
  mpz_mat_dotCols(&d.entry[0][1], B, 0, 0);
		  
  if (H) {
    H->rows = H->cols = A->cols;
    for (i=0; i<H->rows; i++)
      for (j=0; j<H->cols; j++)
        if (i!=j)
          mpz_set_ui(&H->entry[i][j], 0);
        else
          mpz_set_ui(&H->entry[i][j], 1);

  }
  
ALG_2_6_7_STEP2:
  /* (2) Incremental Gram-Schmidt. */
  if (k <= k_max)
    goto ALG_2_6_7_STEP3;
  k_max = k;
  for (j=1; j <= k; j++) {
    mpz_mat_dotCols(u, B, k-1, j-1);
    for (i=1; i<=j-1; i++) {
      mpz_mul(tmp1, &d.entry[0][i], u);
      mpz_mul(tmp2, &lam.entry[k-1][i-1], &lam.entry[j-1][i-1]);
      mpz_sub(u, tmp1, tmp2);
      mpz_div(u, u, &d.entry[0][i-1]);
    }
    if (j < k) 
      mpz_set(&lam.entry[k-1][j-1], u);
    if (j == k) {
      mpz_set(&d.entry[0][k], u);
      if (mpz_cmp_ui(u, 0)==0) {
        fprintf(stderr, "mpz_mat_LLL() Error - input was not a basis!\n");
	res = -1;
	goto ALG_2_6_7_DONE;
      }
    }
  }
    
ALG_2_6_7_STEP3:
  /* (3) Test LLL condition. */
  redI(k, k-1, B, H);
  mpz_mul(tmp1, &d.entry[0][k-1], &d.entry[0][k-1]);
  mpz_mul_ui(tmp1, tmp1, 3);
  mpz_div_ui(u, tmp1, 4);
  mpz_mul(tmp2, &lam.entry[k-1][k-2], &lam.entry[k-1][k-2]);
  mpz_sub(u, u, tmp2); 
  mpz_mul(tmp1, &d.entry[0][k], &d.entry[0][k-2]);
  if (mpz_cmp(tmp1, u) < 0) {
    swapI(k, k_max, B, H);
    k = MAX(2, k-1);
    goto ALG_2_6_7_STEP3;
  }
  for (l=k-2; l>=1; l--)
    redI(k, l, B, H);
  k++;

  /* (4) Finished? */
  if (k <= n)
    goto ALG_2_6_7_STEP2;

ALG_2_6_7_DONE:
  return res;
}

/***************************************************************/
int redI(int k, int l, mpz_mat_t *B, mpz_mat_t *H)
{ int i;
	
  mpz_mul_2exp(tmp1, &lam.entry[k-1][l-1], 1);
  if (mpz_sgn(tmp1) < 0)
    mpz_neg(tmp1, tmp1);
  if (mpz_cmp(tmp1, &d.entry[0][l]) <= 0)
    return 0;
  mpz_nearest_int(q, &lam.entry[k-1][l-1], &d.entry[0][l]);
  if (H) {
    for (i=0; i<H->rows; i++) {
      mpz_mul(tmp1, q, &H->entry[i][l-1]);
      mpz_sub(&H->entry[i][k-1], &H->entry[i][k-1], tmp1);
    }
  }
  for (i=0; i<B->rows; i++) {
    mpz_mul(tmp1, q, &B->entry[i][l-1]);
    mpz_sub(&B->entry[i][k-1], &B->entry[i][k-1], tmp1);
  }
  mpz_mul(tmp1, q, &d.entry[0][l]);
  mpz_sub(&lam.entry[k-1][l-1], &lam.entry[k-1][l-1], tmp1);
  for (i=1; i <= l-1; i++) {
    mpz_mul(tmp1, q, &lam.entry[l-1][i-1]);
    mpz_sub(&lam.entry[k-1][i-1], &lam.entry[k-1][i-1], tmp1);
  }
  return 0;
}
  

/***************************************************************/
int swapI(int k, int k_max, mpz_mat_t *B, mpz_mat_t *H)
{ int i, j;
 
  if (H) {	
    for (i=0; i<H->rows; i++) {
      mpz_set(tmp1, &H->entry[i][k-1]);
      mpz_set(&H->entry[i][k-1], &H->entry[i][k-2]);
      mpz_set(&H->entry[i][k-2], tmp1);
    }
  }
  for (i=0; i<B->rows; i++) {
    mpz_set(tmp1, &B->entry[i][k-1]);
    mpz_set(&B->entry[i][k-1], &B->entry[i][k-2]);
    mpz_set(&B->entry[i][k-2], tmp1);
  }
  for (j=1; j <= k-2; j++) {
    mpz_set(tmp1, &lam.entry[k-1][j-1]);
    mpz_set(&lam.entry[k-1][j-1], &lam.entry[k-2][j-1]);
    mpz_set(&lam.entry[k-2][j-1], tmp1);
  }
  mpz_set(ltmp, &lam.entry[k-1][k-2]);
  
  mpz_mul(btmp, &d.entry[0][k-2], &d.entry[0][k]);
  mpz_mul(tmp1, ltmp, ltmp);
  mpz_add(btmp, btmp, tmp1);
  mpz_div(btmp, btmp, &d.entry[0][k-1]);

  for (i=k+1; i<= k_max; i++) {
    mpz_set(t, &lam.entry[i-1][k-1]);
    mpz_mul(tmp1, &d.entry[0][k], &lam.entry[i-1][k-2]);
    mpz_mul(tmp2, ltmp, t);
    mpz_sub(&lam.entry[i-1][k-1], tmp1, tmp2);
    mpz_div(&lam.entry[i-1][k-1], &lam.entry[i-1][k-1], &d.entry[0][k-1]);

    mpz_mul(tmp1, btmp, t);
    mpz_mul(tmp2,  ltmp, &lam.entry[i-1][k-1]);
    mpz_add(&lam.entry[i-1][k-2], tmp1, tmp2);
    mpz_div(&lam.entry[i-1][k-2], &lam.entry[i-1][k-2], &d.entry[0][k]);
  }
  mpz_set(&d.entry[0][k-1], btmp);
  return 0;
  
}

/***************************************************************/
void mpz_mat_mul(mpz_mat_t *R, mpz_mat_t *X, mpz_mat_t *Y)
{ int i, j, k;
  mpz_t tmp;
  mpz_mat_t Res;

  mpz_mat_init2(&Res, X->rows, Y->cols);	
  mpz_init(tmp);
  Res.rows = X->rows;
  Res.cols = Y->cols;
  for (i=0; i<Res.rows; i++)
    for (j=0; j<Res.cols; j++)
      mpz_set_ui(&Res.entry[i][j], 0);
  
  for (i=0; i<X->rows; i++) {
    for (j=0; j<X->cols; j++) {
      for (k=0; k<Y->cols; k++) {
        mpz_mul(tmp, &X->entry[i][j], &Y->entry[j][k]);
	mpz_add(&Res.entry[i][k], &Res.entry[i][k], tmp);
      }
    }
  }
  mpz_mat_cp(R, &Res);
  mpz_clear(tmp);
  mpz_mat_clear(&Res);
}

/***************************************************************/
void mpz_mat_cp(mpz_mat_t *dest, mpz_mat_t *src)
{ int i, j, r, c;

  r = dest->rows = src->rows;
  c = dest->cols = src->cols;
  for (i=0; i<r; i++)
    for (j=0; j<c; j++)
      mpz_set(&dest->entry[i][j], &src->entry[i][j]);
}

/***************************************************************/
void mpz_mat_cat(mpz_mat_t *dest, mpz_mat_t *src1, mpz_mat_t *src2)
{ int i, j, r, c1, c2;

  r = dest->rows = src1->rows;
  c1 = src1->cols;
  c2 = src2->cols;
  dest->cols = c1+c2;
  for (i=0; i<r; i++) {
    for (j=0; j<c1; j++)
      mpz_set(&dest->entry[i][j], &src1->entry[i][j]);
    for (j=0; j<c2; j++)
      mpz_set(&dest->entry[i][c1+j], &src2->entry[i][j]);
  }
}


/***************************************************************/
int mpz_mat_LatticeIntersection(mpz_mat_t *I, mpz_mat_t *B_i, int s)
/***************************************************************/
/* Given triangular lattice bases B_i[0],...,B_i[s-1],         */
/* compute an LLL-reduced basis for the intersection and store */
/* it in I.                                                    */
/***************************************************************/
{ mpz_mat_t B, T, U;
  int i;
  
  mpz_mat_init(&B); mpz_mat_init(&T); mpz_mat_init(&U);
  mpz_mat_cp(I, &B_i[0]);
  
  for (i=1; i<s; i++) {
    mpz_mat_cat(&B, &B_i[i], I);
    mpz_mat_getKernel(&T, &B);
    mpz_mat_mul(&U, &B_i[i], &T);
    mpz_mat_LLL(I, NULL, &U);
  }
  if (s==1) {
    mpz_mat_cp(&U, I);
    mpz_mat_LLL(I, NULL, &U);
  }

  mpz_mat_clear(&B); mpz_mat_clear(&T); mpz_mat_clear(&U);
  return 0;
}
    

/***************************************************************/
int mpz_mat_getKernel(mpz_mat_t *I, mpz_mat_t *B)
/***************************************************************/
/* Compute a basis for the Z-kernel of 'B'.                    */
/***************************************************************/
{ int       i, j, k, k_max, l, res=0, m, n, f[256], r;
  mpz_mat_t H;
  
  initGlobals();
  mpz_mat_init(&H);
  
  m = B->rows;
  n = B->cols;

  for (i=0; i<=n; i++)
    f[i] = 0;
  
  /* (1) Initialize. */
  k=2; k_max=1;
  mpz_set_ui(&d.entry[0][0], 1);
  mpz_mat_dotCols(t, B, 0, 0);
  if (mpz_cmp_ui(t,0)) {
    mpz_set(&d.entry[0][1], t);
    mpz_set_ui(&lam.entry[0][0], 1);
    f[1] = 1;
  } else {
    mpz_set_ui(&d.entry[0][1], 1);
    f[1] = 0;
  }
		  
  H.rows = H.cols = n;
  for (i=0; i<H.rows; i++)
    for (j=0; j<H.cols; j++) {
      if (i!=j)
        mpz_set_ui(&H.entry[i][j], 0);
      else
        mpz_set_ui(&H.entry[i][j], 1);
      mpz_set_ui(&lam.entry[i][j], 0);
    }
  
ALG_2_7_2_STEP2:
  /* (2) Incremental Gram-Schmidt. */
  if (k <= k_max)
    goto ALG_2_7_2_STEP3;
  k_max = k;
  for (j=1; j <= k; j++) {
    if ((f[j]==0)&&(j<k)) {
      mpz_set_ui(&lam.entry[k-1][j-1], 0);
    } else {
      mpz_mat_dotCols(u, B, k-1, j-1);
      for (i=1; i<=j-1; i++) {
        if (f[i]) {
          mpz_mul(tmp1, &d.entry[0][i], u);
          mpz_mul(tmp2, &lam.entry[k-1][i-1], &lam.entry[j-1][i-1]);
          mpz_sub(u, tmp1, tmp2);
          mpz_div(u, u, &d.entry[0][i-1]);
	}
      }
      if (j < k) 
        mpz_set(&lam.entry[k-1][j-1], u);
      else {
        if (mpz_cmp_ui(u, 0)) {
          mpz_set(&d.entry[0][k], u);
	  mpz_set_ui(&lam.entry[k-1][k-1], 1);
          f[k] = 1;
        } else {
          mpz_set(&d.entry[0][k], &d.entry[0][k-1]);
          f[k]=0;
        }
      }
    }
  }
  
ALG_2_7_2_STEP3:
  /* (3) Test f_k = 0 and f_{k-1} != 0 */
  if (f[k-1] && (f[k]==0)) {
    redI(k, k-1, B, &H);
    swapK(k, k_max, f, B, &H);
    k = MAX(2, k-1);
    goto ALG_2_7_2_STEP3;
  }

  
  for (l=k-1; l>=1; l--)
    if (f[l]) 
      redI(k, l, B, &H);
  k++;

  /* (4) Finished? */
  if (k <= n)
    goto ALG_2_7_2_STEP2;

  r=1;
  while ((r <= n) && (f[r]==0))
    r++;
  r--;
  
  H.cols = r;
  mpz_mat_cp(I, &H);
  
  mpz_mat_clear(&H);
  return res;
}

/*********************************************************/
int swapK(int k, int k_max, int *f, mpz_mat_t *B, mpz_mat_t *H)
/*********************************************************/
{ int i, j;
 
  /* Note the error in Cohen's description right here:
   * We must, in addition, swap the cols of 'B' and
   * more importantly, set lambda before swapping the
   * lambda_{i,j}!
  */
  mpz_set(ltmp, &lam.entry[k-1][k-2]);
  for (i=0; i<H->rows; i++) {
    mpz_set(tmp1, &H->entry[i][k-1]);
    mpz_set(&H->entry[i][k-1], &H->entry[i][k-2]);
    mpz_set(&H->entry[i][k-2], tmp1);
  }
  for (i=0; i<B->rows; i++) {
    mpz_set(tmp1, &B->entry[i][k-1]);
    mpz_set(&B->entry[i][k-1], &B->entry[i][k-2]);
    mpz_set(&B->entry[i][k-2], tmp1);
  }

  for (j=1; j <= k-2; j++) {
    mpz_set(tmp1, &lam.entry[k-1][j-1]);
    mpz_set(&lam.entry[k-1][j-1], &lam.entry[k-2][j-1]);
    mpz_set(&lam.entry[k-2][j-1], tmp1);
  }
  if (mpz_sgn(ltmp)==0) {
    mpz_set(&d.entry[0][k-1], &d.entry[0][k-2]);
    f[k-1] = 0; f[k] = 1;
    mpz_set_ui(&lam.entry[k-1][k-2], 0);
    for (i=k+1; i<=k_max; i++) {
      mpz_set(&lam.entry[i-1][k-1], &lam.entry[i-1][k-2]);
      mpz_set_ui(&lam.entry[i-1][k-2], 0);
    }
  } else {
    for (i=k+1; i<= k_max; i++) {
      mpz_mul(&lam.entry[i-1][k-2], ltmp, &lam.entry[i-1][k-2]);
      mpz_div(&lam.entry[i-1][k-2], &lam.entry[i-1][k-2], &d.entry[0][k-1]);
    }
    mpz_set(t, &d.entry[0][k]);
    mpz_mul(tmp1, ltmp, ltmp);
    mpz_div(&d.entry[0][k-1], tmp1, &d.entry[0][k-1]);
    mpz_set(&d.entry[0][k], &d.entry[0][k-1]);
    for (j=k+1; j < k_max; j++) {
      for (i=j+1; i <= k_max; i++) {
        mpz_mul(&lam.entry[i-1][j-1], &lam.entry[i-1][j-1], &d.entry[0][k-1]);
        mpz_div(&lam.entry[i-1][j-1], &lam.entry[i-1][j-1], t);
      }
    }
    for (j=k+1; j <= k_max; j++) {
      mpz_mul(&d.entry[0][j], &d.entry[0][j], &d.entry[0][k-1]);
      mpz_div(&d.entry[0][j], &d.entry[0][j], t);
    }
  }
  return 0;
  
}

#define MAX_HNF_MAT_SIZE 36
/**********************************************************/
int mpz_mat_getHNF(mpz_mat_t *W, mpz_mat_t *_A)
/**********************************************************/
/* Compute the Hermite Normal Form of a given matrix. See */
/* Cohen, Algorithm 2.4.4.                                */
/**********************************************************/
/* I'm not entirely sure this works if _A->rows > _A->cols*/
/* but we don't really need it for that case anyway.      */
/**********************************************************/
{ int i, j, k, l, index, minLoc;
  int m = _A->rows, n = _A->cols;
  static mpz_mat_t A;
  static mpz_t tmp, minEntry, b, q;
  static int initialized=0;  

  if (!(initialized)) {
    mpz_mat_init2(&A, MAX_HNF_MAT_SIZE, MAX_HNF_MAT_SIZE);
    mpz_init(tmp); mpz_init(minEntry); mpz_init(b); mpz_init(q);
    initialized=1;
  }
		  
  mpz_mat_cp(&A, _A);
  
  i = m; k = n;
  if (m <= n)
    l=1;
  else
    l = m-n+1;

ALG_2_4_4_STEP2:
  
  minLoc=-1;
  for (j=1; j < k; j++) {
    if (mpz_sgn(&A.entry[i-1][j-1])) {
      if (minLoc < 0) {
        minLoc = j;
	mpz_abs(minEntry, &A.entry[i-1][j-1]);
      }
      else {
        mpz_abs(tmp, &A.entry[i-1][j-1]);
        if (mpz_cmp(tmp, minEntry) < 0) {
          mpz_set(minEntry, tmp);
          minLoc = j;
        }
      }
    }
  }
  
  if (minLoc==-1) {
    if (mpz_sgn(&A.entry[i-1][k-1]) < 0) {
      for (j=1; j<=A.rows; j++)
        mpz_neg(&A.entry[j-1][k-1], &A.entry[j-1][k-1]);
    }
    goto ALG_2_4_4_STEP5;
  }
  

  /* Step 3: */
  mpz_abs(tmp, &A.entry[i-1][k-1]);
  if (mpz_sgn(tmp) && ((minLoc < 0) || (mpz_cmp(tmp, minEntry) < 0)))
    minLoc = k;

  if (minLoc < k) {
    /* Swap columns k and minLoc. */
    for (j=1; j<=A.rows; j++) 
      mpz_swap(&A.entry[j-1][k-1], &A.entry[j-1][minLoc-1]);
    if (mpz_sgn(&A.entry[i-1][k-1]) < 0)
      for (j=1; j<=A.rows; j++)
        mpz_neg(&A.entry[j-1][k-1], &A.entry[j-1][k-1]);
  }
  mpz_set(b, &A.entry[i-1][k-1]);

  /* Step 4 */
  for (j=1; j<k; j++) {
    mpz_nearest_int(q, &A.entry[i-1][j-1], b);
    for (index=1; index<=A.rows; index++) {
      mpz_mul(tmp, q, &A.entry[index-1][k-1]);
      mpz_sub(&A.entry[index-1][j-1], &A.entry[index-1][j-1], tmp);
    }
  }
  goto ALG_2_4_4_STEP2;

ALG_2_4_4_STEP5:
  mpz_set(b, &A.entry[i-1][k-1]);
  if (mpz_sgn(b)==0) {
    k++;
    goto ALG_2_4_4_STEP6;  
  }
  for (j=k+1; j<=A.cols; j++) {
    mpz_fdiv_q(q, &A.entry[i-1][j-1], b);  
    for (index=1; index<=A.rows; index++) {
      mpz_mul(tmp, q, &A.entry[index-1][k-1]);
      mpz_sub(&A.entry[index-1][j-1], &A.entry[index-1][j-1], tmp);
    }
  }

ALG_2_4_4_STEP6:
  if (i==l) {
    for (j=1; j<=(n-k+1); j++) {
      for (index=1; index<=A.rows; index++)
        mpz_set(&W->entry[index-1][j-1], &A.entry[index-1][j+k-1 - 1]);
    }
    W->rows = A.rows;
    W->cols = n-k+1;
  } else {
    i--;
    k--;
    goto ALG_2_4_4_STEP2;
  }
  return n-k;
}

/*************************************************************/
int mpz_mat_invIm(mpz_mat_t *X, mpz_mat_t *V, mpz_mat_t *_M, mpz_t p)
/*************************************************************/
/* Cohen, Alg. 2.3.5. Solve V=MX mod p.                      */
{ int       i, j, k, l, m=_M->rows, n=_M->cols, r=V->cols, retVal=0;
  mpz_mat_t M, B, C;
  mpz_t     d, tmp, tmpSum;

  mpz_mat_init2(&M, m, n); mpz_mat_cp(&M, _M);
  mpz_mat_init2(&B, m, r); 
  mpz_mat_init2(&C, m, 1); 
  mpz_init(d); mpz_init(tmp); mpz_init(tmpSum);

  /* Step 1: */
  j=0;
  mpz_mat_cp(&B, V);

  /* Step 2: */
ALG_2_3_5_S2:
  j++;
  if (j>n)
    goto ALG_2_3_5_S6;

  /* Step 3: */
  for (i=m; i>=j; i--)
    if (mpz_sgn(&M.entry[i-1][j-1]))
      break;
  if (i<j) {
    printf("mpz_mat_invIm(): Error - col's of M are not lin. ind.!\n");
    retVal = -1;
    goto ALG_2_3_5_DONE;
  }

  /* Step 4: */
  if (i>j) {
    for (l=j; l<=n; l++)
      mpz_swap(&M.entry[i-1][l-1], &M.entry[j-1][l-1]);
    for (l=1; l<=B.cols; l++)
      mpz_swap(&B.entry[i-1][l-1], &B.entry[j-1][l-1]);
  }

  /* Step 5: */
  mpz_invert(d, &M.entry[j-1][j-1], p);
  for (k=m; k>j; k--)
    mpz_mul(&C.entry[k-1][0], d, &M.entry[k-1][j-1]);
  for (k=m; k>j; k--) {
    for (l=n; l>j; l--) {
      mpz_mul(tmp, &C.entry[k-1][0], &M.entry[j-1][l-1]);
      mpz_sub(&M.entry[k-1][l-1], &M.entry[k-1][l-1], tmp);
      mpz_mod(&M.entry[k-1][l-1], &M.entry[k-1][l-1], p);
    }
  }
  for (k=m; k>j; k--) {
    for (l=1; l<=B.cols; l++) {
      mpz_mul(tmp, &C.entry[k-1][0], &B.entry[j-1][l-1]);
      mpz_sub(&B.entry[k-1][l-1], &B.entry[k-1][l-1], tmp);
      mpz_mod(&B.entry[k-1][l-1], &B.entry[k-1][l-1], p);
    }
  }
  goto ALG_2_3_5_S2;

  /* Step 6: */
ALG_2_3_5_S6:
  X->rows = n; X->cols = r;

  for (i=n; i>=1; i--) {
    /* First, X_i' <-- B_i' */
    for (l=1; l<=X->cols; l++)
      mpz_set(&X->entry[i-1][l-1], &B.entry[i-1][l-1]);

    for (j=i+1; j<=n; j++) {
      /* Now, X_i' <-- X_i' - m_{ij}X_j' */ 
      mpz_set(d, &M.entry[i-1][j-1]);
      for (l=1; l<=X->cols; l++) {
        mpz_mul(tmp, d, &X->entry[j-1][l-1]);
        mpz_sub(&X->entry[i-1][l-1], &X->entry[i-1][l-1], tmp);
        mpz_mod(&X->entry[i-1][l-1], &X->entry[i-1][l-1], p);
      }
    }
    /* Finally, X_i' <-- X_i' / m_{ii}. */
    mpz_invert(d, &M.entry[i-1][i-1], p);
    for (l=1; l<=X->cols; l++) {
      mpz_mul(&X->entry[i-1][l-1], &X->entry[i-1][l-1], d);
      mpz_mod(&X->entry[i-1][l-1], &X->entry[i-1][l-1], p);
    }
  }


  
  /* Step 7: (verify) */
  mpz_mat_mul(&B, _M, X);
  for (i=0; i<B.rows; i++) {
    for (j=0; j<B.cols; j++) {
      mpz_sub(tmp, &B.entry[i][j], &V->entry[i][j]);
      mpz_mod(tmp, tmp, p);
      if (mpz_sgn(tmp)) {
        retVal = -1;
        printf("mpz_mat_invIm() : Error: V is not in Colsp(M)!\n");
        printf("V =\n"); mpz_mat_print(stdout, V);
        printf("M =\n"); mpz_mat_print(stdout, _M);
        printf("Attempted solution of V = MX was:\nX = \n");
        mpz_mat_print(stdout, X);
        printf("But MX = \n");
        mpz_mat_print(stdout, &B);
        goto ALG_2_3_5_DONE;
      }
    }
  }

ALG_2_3_5_DONE:
  mpz_mat_clear(&M); mpz_mat_clear(&B);
  mpz_mat_clear(&C);
  mpz_clear(d); mpz_clear(tmp); mpz_clear(tmpSum);
  return retVal;
}

  
/*************************************************************/
int mpz_mat_suppSubspace_modp(mpz_mat_t *Z, mpz_mat_t *V, mpz_mat_t *M, mpz_t p)
/*************************************************************/
/* Cohen, Alg. 2.3.7.                                        */
/* Find a basis for the supplement of Colsp(V) in Colsp(M).  */
/*************************************************************/
{ int       i, j, r=V->cols, n=M->cols, retVal=0;
  mpz_mat_t X, C;

  mpz_mat_init2(&X, n, MAX(r,n));
  mpz_mat_init2(&C, n, n-r);

  if (mpz_mat_invIm(&X, V, M, p)) {
    printf("mpz_mat_suppSubspace_modp(): Error - Colsp(V) not in Colsp(M)!\n");
    retVal = -1; goto ALG_2_3_7_DONE;
  }

  mpz_mat_suppBasis_modp(&X, &X, p);

  for (j=0; j<n-r; j++)
    for (i=0; i<n; i++)
      mpz_set(&C.entry[i][j], &X.entry[i][2*r-n+j+1]);
  C.rows = n; C.cols = n-r;

  mpz_mat_mul(Z, M, &C);
  for (i=0; i<Z->rows; i++)
    for (j=0; j<Z->cols; j++)
      mpz_mod(&Z->entry[i][j], &Z->entry[i][j], p);


ALG_2_3_7_DONE:
  mpz_mat_clear(&X); mpz_mat_clear(&C);
  return retVal;
}
  
	
    
/*************************************************************/
int mpz_mat_getHNF_mod_D(mpz_mat_t *W, mpz_mat_t *_A, mpz_t D)
/*************************************************************/
/* HNF when it's known that the determinant of the generated */
/* Z-module divides D. See Cohen, Alg. 2.4.8.                */
/* A must have full rank!                                    */
/*************************************************************/
{ int       m=_A->rows, n=_A->cols;
  mpz_mat_t A, B;
  int       i, j, k, l;
  mpz_t     R, R_2, u, v, d, tmp1, tmp2, tmp3, q;

  mpz_mat_init2(&A, m, n);
  mpz_mat_init2(&B, m, 1);
  mpz_mat_cp(&A, _A);
  mpz_init(R); mpz_init(R_2); mpz_init(u); mpz_init(v); mpz_init(d);
  mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3); mpz_init(q);
  
  /* Step 1: */
  i=m; j=n; k=n; mpz_set(R, D); mpz_div_2exp(R_2, R, 1);
  /* Step 2: */
ALG_2_4_8_STEP2:
  if (j==1)
    goto ALG_2_4_8_STEP4;
  j--;
  if (mpz_sgn(&A.entry[i-1][j-1])==0)
    goto ALG_2_4_8_STEP2;
  
  /* Step 3: */
    /* Ext. gcd, with minimal u, v. */
  mpz_gcdext(d, u, v, &A.entry[i-1][k-1], &A.entry[i-1][j-1]);
  mpz_lcm(tmp1, &A.entry[i-1][k-1], &A.entry[i-1][j-1]);
  mpz_abs(tmp1, tmp1);
  mpz_mod(u, u, tmp1); mpz_mod(v, v, tmp1);
  mpz_div_2exp(tmp2, tmp1, 1);
  if (mpz_cmp(u, tmp2) > 0)
    mpz_sub(u, u, tmp1);  
  if (mpz_cmp(v, tmp2) > 0)
    mpz_sub(v, v, tmp1);  
    /* B <-- u*A_k + v*A_j. */
  for (l=0; l<m; l++) { 
    mpz_mul(&B.entry[l][0], u, &A.entry[l][k-1]);
    mpz_mul(tmp1, v, &A.entry[l][j-1]);
    mpz_add(&B.entry[l][0], &B.entry[l][0], tmp1);
  }
  
  mpz_div(tmp1, &A.entry[i-1][k-1], d);
  mpz_div(tmp2, &A.entry[i-1][j-1], d);
  mpz_neg(tmp2, tmp2);
  /* Now A_j <-- tmp1*A_j + tmp2*A_k mod R. */
  for (l=1; l<=m; l++) {
    mpz_mul(tmp3, tmp2, &A.entry[l-1][k-1]);
    mpz_mul(&A.entry[l-1][j-1], &A.entry[l-1][j-1], tmp1);
    mpz_add(&A.entry[l-1][j-1], &A.entry[l-1][j-1], tmp3);
    mpz_mod(&A.entry[l-1][j-1], &A.entry[l-1][j-1], R);
    if (mpz_cmp(&A.entry[l-1][j-1], R_2)>0)
      mpz_sub(&A.entry[l-1][j-1], &A.entry[l-1][j-1], R);
  }
  /* A_k <-- B mod R. */
  for (l=1; l<=m; l++) {
    mpz_set(&A.entry[l-1][k-1], &B.entry[l-1][0]);
    if (mpz_cmp(&A.entry[l-1][k-1], R_2) > 0)
      mpz_sub(&A.entry[l-1][k-1], &A.entry[l-1][k-1], R);
  }
  goto ALG_2_4_8_STEP2;

  
  /* Step 4: */
ALG_2_4_8_STEP4:

  mpz_gcdext(d, u, v, &A.entry[i-1][k-1], R);
  for (l=1; l<=m; l++) {
    mpz_mul(&W->entry[l-1][i-1], u, &A.entry[l-1][k-1]);
    mpz_mod(&W->entry[l-1][i-1], &W->entry[l-1][i-1], R);
  }

  if (mpz_sgn(&W->entry[i-1][i-1])==0)
    mpz_set(&W->entry[i-1][i-1], R);
  
  for (j=i+1; j<=m; j++) {
    mpz_fdiv_q(tmp1, &W->entry[i-1][j-1], &W->entry[i-1][i-1]);
    for (l=1; l<=m; l++) {
      mpz_mul(tmp2, tmp1, &W->entry[l-1][i-1]);
      mpz_sub(&W->entry[l-1][j-1], &W->entry[l-1][j-1], tmp2);
      mpz_mod(&W->entry[l-1][j-1], &W->entry[l-1][j-1], R);
      if (mpz_cmp(&W->entry[l-1][j-1], R_2) > 0)
        mpz_sub(&W->entry[l-1][j-1], &W->entry[l-1][j-1], R);
    }
  }
  if (i != 1) {
    mpz_div(R, R, d);
    mpz_div_2exp(R_2, R, 1);
    i--;
    k--;
    j=k;
    if (mpz_sgn(&A.entry[i-1][k-1])==0) {
      mpz_set(&A.entry[i-1][k-1], R);
    }
    goto ALG_2_4_8_STEP2;
  }
  
  W->rows = W->cols = m;
  
  mpz_mat_clear(&A);  mpz_mat_clear(&B);
  mpz_clear(R); mpz_clear(R_2); mpz_clear(u); mpz_clear(v); mpz_clear(d);
  mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3); mpz_clear(q);

  return 0;
}
	
/***************************************************************/
int mpz_mat_getKernel_modp(mpz_mat_t *X, mpz_mat_t *B, mpz_t p)
/***************************************************************/
/* Cohen, Algorithm 2.3.1.                                     */
/***************************************************************/
{ mpz_mat_t     M;
  mpz_t         d, tmp;
  int           i, j, k, r, s, m=B->rows, n=B->cols;
  int           *ci, *di, col;


  ci = malloc(m * sizeof(int)); 
  di = malloc(n * sizeof(int));
  mpz_mat_init2(&M, B->rows, B->cols); mpz_init(d); mpz_init(tmp);

  mpz_mat_cp(&M, B);
  for (i=0; i<M.rows; i++)
    for (j=0; j<M.cols; j++)
      mpz_mod(&M.entry[i][j], &M.entry[i][j], p);

  /* Step 1. */
  r=0; k=1;
  for (i=1; i<=m; i++)
    ci[i-1] = 0;

ALG_2_3_1_STEP2:
  /* Step 2. */
  j = 1;
  while ((j<=m) && ((mpz_sgn(&M.entry[j-1][k-1])==0) || (ci[j-1])))
    j++;
  if (j>m) {
    r++;
    di[k-1] = 0;
    goto ALG_2_3_1_STEP4;
  }

  /* Step 3. */
  mpz_invert(d, &M.entry[j-1][k-1], p);
  mpz_neg(d, d);  mpz_mod(d, d, p);
  
  mpz_set(&M.entry[j-1][k-1], p);
  mpz_sub_ui(&M.entry[j-1][k-1], &M.entry[j-1][k-1], 1);
  
  for (s=k+1; s<=n; s++) {
    mpz_mul(&M.entry[j-1][s-1], d, &M.entry[j-1][s-1]);
    mpz_mod(&M.entry[j-1][s-1], &M.entry[j-1][s-1], p);
  }
  for (i=1; i<=m; i++) {
    if (i != j) {
      mpz_set(d, &M.entry[i-1][k-1]);
      mpz_set_ui(&M.entry[i-1][k-1], 0);
      for (s=k+1; s<=n; s++) {
        mpz_mul(tmp, d, &M.entry[j-1][s-1]);
	mpz_add(&M.entry[i-1][s-1], &M.entry[i-1][s-1], tmp);
	mpz_mod(&M.entry[i-1][s-1], &M.entry[i-1][s-1], p);
      }
    }
  }
  ci[j-1] = k;
  di[k-1] = j;

ALG_2_3_1_STEP4:
  /* Step 4. */
  if (k<n) {
    k++;
    goto ALG_2_3_1_STEP2;
  }

  /* Step 5. */
  X->rows = n;
  col=0;
  for (k=1; k<=n; k++) {
    if (di[k-1]==0) {
      for (i=1; i<=n; i++) {
        if (di[i-1] > 0) {
          mpz_set(&X->entry[i-1][col], &M.entry[di[i-1]-1][k-1]);
	} else if (i==k) {
          mpz_set_ui(&X->entry[i-1][col], 1);
	} else {
          mpz_set_ui(&X->entry[i-1][col], 0);
        }
      }
      col++;
    }
  }
  X->cols = col;
  if (col != r) {
    printf("mpz_mat_getKernel_modp(): Some error occurred. col=%d, r=%d.\n", col, r);
  }
  mpz_mat_clear(&M); mpz_clear(d); mpz_clear(tmp);
  free(ci); free(di);

  return col;
}

/***************************************************************/
int mpz_mat_getImage_modp(mpz_mat_t *Y, mpz_mat_t *B, mpz_t p)
/***************************************************************/
/* Cohen, Algorithm 2.3.2.                                     */
/***************************************************************/
{ mpz_mat_t     M;
  mpz_t         d, tmp;
  int           i, j, k, r, s, m=B->rows, n=B->cols;
  int           *ci, *di, col;


  ci = malloc(m * sizeof(int)); 
  di = malloc(n * sizeof(int));
  mpz_mat_init2(&M, m, n); mpz_init(d); mpz_init(tmp);

  mpz_mat_cp(&M, B);
  for (i=0; i<M.rows; i++)
    for (j=0; j<M.cols; j++)
      mpz_mod(&M.entry[i][j], &M.entry[i][j], p);

  /* Step 1. */
  r=0; k=1;
  for (i=1; i<=m; i++)
    ci[i-1] = 0;

ALG_2_3_2_STEP2:
  /* Step 2. */
  j = 1;
  while ((j<=m) && ((mpz_sgn(&M.entry[j-1][k-1])==0) || (ci[j-1]))) {
    j++;
  }
  if (j>m) {
    r++;
    di[k-1] = 0;
    goto ALG_2_3_2_STEP4;
  }

  /* Step 3. */
  mpz_invert(d, &M.entry[j-1][k-1], p);
  mpz_neg(d, d);  mpz_mod(d, d, p);
  
  mpz_set(&M.entry[j-1][k-1], p);
  mpz_sub_ui(&M.entry[j-1][k-1], &M.entry[j-1][k-1], 1);
  
  for (s=k+1; s<=n; s++) {
    mpz_mul(&M.entry[j-1][s-1], d, &M.entry[j-1][s-1]);
    mpz_mod(&M.entry[j-1][s-1], &M.entry[j-1][s-1], p);
  }
  for (i=1; i<=m; i++) {
    if (i != j) {
      mpz_set(d, &M.entry[i-1][k-1]);
      mpz_set_ui(&M.entry[i-1][k-1], 0);
      for (s=k+1; s<=n; s++) {
        mpz_mul(tmp, d, &M.entry[j-1][s-1]);
	mpz_add(&M.entry[i-1][s-1], &M.entry[i-1][s-1], tmp);
	mpz_mod(&M.entry[i-1][s-1], &M.entry[i-1][s-1], p);
      }
    }
  }
  ci[j-1] = k;
  di[k-1] = j;

ALG_2_3_2_STEP4:
  /* Step 4. */
  if (k<n) {
    k++;
    goto ALG_2_3_2_STEP2;
  }

  /* Step 5. */
  col=0;
  for (j=1; j<=m; j++) {
    if (ci[j-1]) {
      for (i=1; i<=m; i++) {
//         printf("i=%d, col=%d, j=%d, ci[j-1]=%d\n", i,col,j,ci[j-1]);
         mpz_mod(&Y->entry[i-1][col], &B->entry[i-1][ci[j-1]-1], p);
      }
      col++;
    }
  }
  Y->rows = m;
  Y->cols = col;

  mpz_mat_clear(&M); mpz_clear(d); mpz_clear(tmp);
  free(ci); free(di);

  return col;
}
  
/***************************************************************/
int mpz_mat_moduleAdd(mpz_mat_t *H, mpz_t dH, 
                      mpz_mat_t *W, mpz_t d, mpz_mat_t *W_, mpz_t d_)
/***************************************************************/
/* Input: HNF's (W, d) and (W_, d_).                           */
/* Output: The HNF (H, dH) of the module sum (W,d) + (W_,d_).  */
/***************************************************************/
{ mpz_mat_t A;
  mpz_t     D, tmp;
  int       i, j, m=W->rows;
  
  mpz_mat_init2(&A, m, 2*m);
  mpz_init(D); mpz_init(tmp);

  mpz_lcm(D, d, d_);

  mpz_div(tmp, D, d);
  for (i=0; i<m; i++)
    for (j=0; j<m; j++)
      mpz_mul(&A.entry[i][j], &W->entry[i][j], tmp);

  mpz_div(tmp, D, d_);
  for (i=0; i<m; i++)
    for (j=0; j<m; j++)
      mpz_mul(&A.entry[i][m+j], &W_->entry[i][j], tmp);
  A.rows = m;
  A.cols = 2*m;
 
  mpz_mat_getHNF(H, &A);
  
  mpz_set(tmp, D);
  for (i=0; i<m; i++)
    for (j=0; j<m; j++)
      mpz_gcd(tmp, tmp, &H->entry[i][j]);
  mpz_div(dH, D, tmp);
  for (i=0; i<m; i++)
    for (j=0; j<m; j++)
      mpz_div(&H->entry[i][j], &H->entry[i][j], tmp);

  mpz_mat_clear(&A);
  mpz_clear(D); mpz_clear(tmp);
  return 0;
}

/***************************************************************/
void mpz_mat_pseudoInvert(mpz_mat_t *I, mpz_t d, mpz_mat_t *J)
/***************************************************************/
/* Input: A Q-invertible matrix J with integer entries.        */
/* Output: A matrix I with integer entries and an integer d    */
/*         so that I*J = d*Identity.                           */
/* You probably wouldn't like to use this function on anything */
/* more than 5x5, since it is not terribly smart and will have */
/* intermediate coefficient explosion.                         */
/***************************************************************/
{ int       i, j, k, l, n=J->rows;
  mpz_mat_t U;
  mpz_t     tmp1, tmp2, tmp3;

  mpz_mat_init2(&U, n, 2*n);
  mpz_init(tmp1); mpz_init(tmp2); mpz_init(tmp3);

  mpz_set_ui(d, 1);
  for (i=0;i<n; i++) {
    for (j=0; j<n; j++)
      mpz_set(&U.entry[i][j], &J->entry[i][j]);
    for (j=0; j<n; j++)
      mpz_set_ui(&U.entry[i][j+n], (i==j));
  }
  U.rows = n; U.cols = 2*n;
  
  /* Now, just do elementary row ops to U. */
  for (j=0; j<n; j++) {
    /* Find a nonzero entry in column j. */
    for (i=j; i<n; i++)
      if (mpz_sgn(&U.entry[i][j]))
        break;
    if (i>=n) {
      fprintf(stderr, "pseudoInvert() Error - matrix not Q-invertible!\n");
      exit(-1);
    }
    /* Swap rows i and j: */
    for (k=0; k<2*n; k++)
      mpz_swap(&U.entry[i][k], &U.entry[j][k]);

    mpz_set(tmp1, &U.entry[j][j]);
    /* ``divide'' row j by tmp1. */
    mpz_mul(d, d, tmp1);
    for (k=0; k<n; k++) {
      if (k!=j)
        for (l=0; l<2*n; l++)
          mpz_mul(&U.entry[k][l], &U.entry[k][l], tmp1);
    }

    /* Zero out the rest of the j-th column. */
    for (k=0; k<n; k++) {
      if (k!=j) {
        mpz_div(tmp2, &U.entry[k][j], &U.entry[j][j]);
        /* Row k <-- (Row k) - tmp2*(Row j). */
        for (l=0; l<2*n; l++) {
          mpz_mul(tmp3, tmp2, &U.entry[j][l]);
          mpz_sub(&U.entry[k][l], &U.entry[k][l], tmp3);
        }
      }
    }
    /* Perhaps we can reduce some coefficients now? */
    mpz_set(tmp1, d);
    for (k=0; k<n; k++) {
      for (l=0; l<2*n; l++)
        mpz_gcd(tmp1, tmp1, &U.entry[k][l]);
      if (mpz_cmp_ui(tmp1, 1)==0)
        break;
    }
    if (mpz_cmp_ui(tmp1, 1)) {
      mpz_div(d, d, tmp1);
      for (k=0; k<n; k++)
        for (l=0; l<2*n; l++)
          mpz_div(&U.entry[k][l], &U.entry[k][l], tmp1);
    }
  }
  /* Now we have a diagonalization, so there's a little
     multiplying left to do:
  */


  for (i=0; i<n; i++) 
    if (mpz_sgn(&U.entry[i][i]) < 0)
      for (j=0; j<2*n; j++)
        mpz_neg(&U.entry[i][j], &U.entry[i][j]);


  mpz_set_ui(tmp1, 1);
  for (i=0; i<n; i++)
    mpz_lcm(tmp1, tmp1, &U.entry[i][i]);

  mpz_set(d, tmp1);
  for (i=0; i<n; i++) {
    mpz_div(tmp2, tmp1, &U.entry[i][i]);
    for (j=0; j<n; j++) {
      mpz_mul(&I->entry[i][j], tmp2, &U.entry[i][j+n]);
    }
  }
  I->rows = I->cols = n;

  mpz_set(tmp1, d);
  for (k=0; k<n; k++) {
    for (l=0; l<n; l++)
      mpz_gcd(tmp1, tmp1, &I->entry[k][l]);
    if (mpz_cmp_ui(tmp1, 1)==0)
      break;
  }
  mpz_abs(tmp1, tmp1);
  if (mpz_cmp_ui(tmp1, 1)) {
    mpz_div(d, d, tmp1);
    for (k=0; k<n; k++)
      for (l=0; l<n; l++)
        mpz_div(&I->entry[k][l], &I->entry[k][l], tmp1);
  }

  mpz_mat_clear(&U);
  mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(tmp3);
}
 
    

#ifdef _TEST_MAIN		  
int main()
{ mpz_mat_t Y, X;
  mpz_t p;

  mpz_mat_init(&X);
  mpz_mat_init(&Y);
  mpz_init_set_ui(p, 7);

  Y.rows = 2;
  Y.cols = 3;
  mpz_set_ui(&Y.entry[0][0], 1);
  mpz_set_ui(&Y.entry[1][0], 3);
  mpz_set_ui(&Y.entry[0][1], 2);
  mpz_set_ui(&Y.entry[1][1], 3);
  mpz_set_ui(&Y.entry[0][2], 3);
  mpz_set_ui(&Y.entry[1][2], 5);

  printf("Computing mod 7 kernel of:\n");
  mpz_mat_print(stdout, &Y);
  mpz_mat_getKernel_modp(&X, &Y, p);
  printf("kernel:\n");
  mpz_mat_print(stdout, &X);
}
#endif 
/*************************************************************/
int mpz_mat_suppBasis_modp(mpz_mat_t *B, mpz_mat_t *_M, mpz_t p)
/*************************************************************/
/* Given an nxk matrix _M with rank k <= n, find B=[_M|U]    */
/* an nxn invertible matrix. Cohen Alg. 2.3.6.               */
/* Return value: 0 on success, negative on failure (i.e.,    */
/* _M did not have full rank).                               */
/* Note: I had to tweak this - there seems to be a typo in   */
/* Cohen, at Step 4, that causes failure in most cases.      */
/*************************************************************/
{ int       n=_M->rows, k=_M->cols;
  mpz_mat_t M;
  mpz_t     d, tmp1;
  int       r, i, c, j, t, found;

  mpz_mat_init2(&M, n, n); 
  mpz_mat_cp(&M,_M);
  mpz_init(d); mpz_init(tmp1);

  /* Disaster could occur if M is not reduced mod p! */
  for (i=0; i<n; i++)
    for (j=0; j<k; j++)
      mpz_mod(&M.entry[i][j], &M.entry[i][j], p);
  mpz_mat_cp(B, &M);
  for (i=0; i<n; i++)
    for (j=k; j<n; j++)
      mpz_set_ui(&B->entry[i][j], 0);
  B->rows = B->cols = n;



//printf("Before CREF, M=\n");
//mpz_mat_print(stdout, &M);

  /* Do column ops on M to get it in CREF. */
  i=0;
  c=0;
  while (c<k) {
    /* Find the first column after c with a nonzero entry in row i, if there is one. */
    do {
      for (found=-1, t=c; t<k; t++)
        if (mpz_sgn(&M.entry[i][t])) found = t;
      if (found==-1) i++;
    } while ((found==-1) && (i<n));
    if (i>=n) {
      fprintf(stderr, "Error: M did not have full rank!\n");
      exit(-1); /* testing. */
      return -1;
    }
    /* Swap columns `found' and c. */
    if (found != c) {
      for (t=0; t<n; t++) mpz_swap(&M.entry[t][found], &M.entry[t][c]);
    }
    /* Use column c to zero out the rest of row i. */
    mpz_invert(d, &M.entry[i][c], p);
    for (r=0; r<n; r++) {
      mpz_mul(&M.entry[r][c], &M.entry[r][c], d);
      mpz_mod(&M.entry[r][c], &M.entry[r][c], p);
    }
    for (t=0; t<k; t++) {
      if (t != c) {
        mpz_set(d, &M.entry[i][t]);
        for (r=0; r<n; r++) {
          mpz_mul(tmp1, &M.entry[r][c], d); mpz_mod(tmp1, tmp1, p);
          mpz_sub(&M.entry[r][t], &M.entry[r][t], tmp1);
          mpz_mod(&M.entry[r][t], &M.entry[r][t], p);
        }
      }
    }
    c++;
  }
//  printf("After CREF, M=\n");
//  mpz_mat_print(stdout, &M);

  r=0; c=k;
  for (j=0; j<k; j++) {
    while (mpz_sgn(&M.entry[r][j])==0) {
      mpz_set_ui(&B->entry[r][c], 1);
      r++; c++;
    }
    r++;
  }
  while (r<n) 
    mpz_set_ui(&B->entry[r++][c++], 1);
 
//  printf("After supplementing, B=\n");
//  mpz_mat_print(stdout, B);
    



  mpz_mat_clear(&M); 
  mpz_clear(d); mpz_clear(tmp1);
  return 0;
}  
    

