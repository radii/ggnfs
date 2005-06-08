/**************************************************************/
/* blanczos64-no-mmx.c                                        */
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
#include <string.h>
#include "prand.h"
#include "ggnfs.h"


#ifdef __GNUC__
#define ALIGNED16 __attribute__ ((aligned (16)))
#else
#define ALIGNED16
#endif

#define XOR64(_a, _b) \
  { *(_a) ^= *(_b); }


static const u64 bit64[] ALIGNED16 ={
0x0000000000000001ULL,0x0000000000000002ULL,0x0000000000000004ULL,0x0000000000000008ULL,
0x0000000000000010ULL,0x0000000000000020ULL,0x0000000000000040ULL,0x0000000000000080ULL,
0x0000000000000100ULL,0x0000000000000200ULL,0x0000000000000400ULL,0x0000000000000800ULL,
0x0000000000001000ULL,0x0000000000002000ULL,0x0000000000004000ULL,0x0000000000008000ULL,
0x0000000000010000ULL,0x0000000000020000ULL,0x0000000000040000ULL,0x0000000000080000ULL,
0x0000000000100000ULL,0x0000000000200000ULL,0x0000000000400000ULL,0x0000000000800000ULL,
0x0000000001000000ULL,0x0000000002000000ULL,0x0000000004000000ULL,0x0000000008000000ULL,
0x0000000010000000ULL,0x0000000020000000ULL,0x0000000040000000ULL,0x0000000080000000ULL,
0x0000000100000000ULL,0x0000000200000000ULL,0x0000000400000000ULL,0x0000000800000000ULL,
0x0000001000000000ULL,0x0000002000000000ULL,0x0000004000000000ULL,0x0000008000000000ULL,
0x0000010000000000ULL,0x0000020000000000ULL,0x0000040000000000ULL,0x0000080000000000ULL,
0x0000100000000000ULL,0x0000200000000000ULL,0x0000400000000000ULL,0x0000800000000000ULL,
0x0001000000000000ULL,0x0002000000000000ULL,0x0004000000000000ULL,0x0008000000000000ULL,
0x0010000000000000ULL,0x0020000000000000ULL,0x0040000000000000ULL,0x0080000000000000ULL,
0x0100000000000000ULL,0x0200000000000000ULL,0x0400000000000000ULL,0x0800000000000000ULL,
0x1000000000000000ULL,0x2000000000000000ULL,0x4000000000000000ULL,0x8000000000000000ULL};

#ifdef BIT
#undef BIT
#endif
#ifdef BIT64
#undef BIT64
#endif
#define BIT64(_i) bit64[(_i)]


/*************************************************************/ 
/* Implementation of the Block Lanczos algorithm for finding */
/* vectors in the kernel of a large sparse matrix.           */
/* The implementation here is as descibed in:                */
/* ``A Block Lanczos Algorithm for Finding Dependencies over */
/*   GF(2)'', Peter Montgomery.                              */
/* The paper I have does not have a publication source on it,*/
/* so it's probably a preprint.                              */
/*************************************************************/ 


/***** Prototypes for locally used functions. *****/
void MultB64(u64 *Product, u64 *x, void *P);
void MultB_T64(u64 *Product, u64 *x, void *P);
void multT(u64 *res, u64 *A, u64 *B, s32 n);
void multS(u64 *D, int *S);
void mult64x64(u64 *res, u64 *A, u64 *B);
void preMult(u64 *A, u64 *B);
void multnx64(u64 *C_n, u64 *A_n, u64 *B, s32 n);
void addmultnx64(u64 *C_n, u64 *A_n, u64 *B, s32 n);
void getW_S(u64 *Wi, int *Si, u64 *T, int *Si_1);
int  isZeroV(u64 *A, s32 size);
int  doColumnOps(u64 *A, u64 *B, s32 n);
int  doColumnOps64(u64 *A, u64 *B, s32 n);
/**************************************************/


void seedBlockLanczos(s32 seed)
{ prandseed(seed, 712*seed + 21283, seed^0xF3C91D1A);
}

void MultB64(u64 *Product, u64 *x, void *P) {
  nfs_sparse_mat_t *M = (nfs_sparse_mat_t *)P;
  memset(Product, 0, M->numCols * sizeof(u64)); 
  {
    int i;
    for (i = 0; i < M->numDenseBlocks; i++) {
      multT(Product + M->denseBlockIndex[i], M->denseBlocks[i], x, M->numCols);
    }
  }
#if defined(L2_CACHE_SIZE) && (L2_CACHE_SIZE > 0)
// L2_CACHE_SIZE has to be a power of 2.
// MULTB64_PAGESIZE is a half of L2 cache size.
#define MULTB64_PAGESIZE (L2_CACHE_SIZE * 1024 / 2 / sizeof(u64))
#define MULTB64_PAGEMASK (-MULTB64_PAGESIZE)
  {
    s32 n = M->numCols;
    u32 pagestart;
    for (pagestart = 0; pagestart < M->numCols; pagestart += MULTB64_PAGESIZE) {
      u32 *p = M->cEntry;
      s32 i;
      for (i = 0; i < n; i++) {
        u64 t = x[i];
        u32 *s = p + ((M->cIndex[i + 1] - M->cIndex[i]) & -16);
        for (; p < s; p += 16) {
          if ((p[0] & MULTB64_PAGEMASK) == pagestart) Product[p[0]] ^= t;
          if ((p[1] & MULTB64_PAGEMASK) == pagestart) Product[p[1]] ^= t;
          if ((p[2] & MULTB64_PAGEMASK) == pagestart) Product[p[2]] ^= t;
          if ((p[3] & MULTB64_PAGEMASK) == pagestart) Product[p[3]] ^= t;
          if ((p[4] & MULTB64_PAGEMASK) == pagestart) Product[p[4]] ^= t;
          if ((p[5] & MULTB64_PAGEMASK) == pagestart) Product[p[5]] ^= t;
          if ((p[6] & MULTB64_PAGEMASK) == pagestart) Product[p[6]] ^= t;
          if ((p[7] & MULTB64_PAGEMASK) == pagestart) Product[p[7]] ^= t;
          if ((p[8] & MULTB64_PAGEMASK) == pagestart) Product[p[8]] ^= t;
          if ((p[9] & MULTB64_PAGEMASK) == pagestart) Product[p[9]] ^= t;
          if ((p[10] & MULTB64_PAGEMASK) == pagestart) Product[p[10]] ^= t;
          if ((p[11] & MULTB64_PAGEMASK) == pagestart) Product[p[11]] ^= t;
          if ((p[12] & MULTB64_PAGEMASK) == pagestart) Product[p[12]] ^= t;
          if ((p[13] & MULTB64_PAGEMASK) == pagestart) Product[p[13]] ^= t;
          if ((p[14] & MULTB64_PAGEMASK) == pagestart) Product[p[14]] ^= t;
          if ((p[15] & MULTB64_PAGEMASK) == pagestart) Product[p[15]] ^= t;
        }
        s = M->cEntry + M->cIndex[i + 1];
        for (; p < s; p++) {
          if ((p[0] & MULTB64_PAGEMASK) == pagestart) Product[p[0]] ^= t;
        }
      }
    }
  }
#else
  {
    s32 n = M->numCols;
    u32 *p = M->cEntry;
    s32 i;
    for (i = 0; i < n; i++) {
      u64 t = x[i];
      u32 *s = p + ((M->cIndex[i + 1] - M->cIndex[i]) & -16);
      for (; p < s; p += 16) {
        Product[p[0]] ^= t;
        Product[p[1]] ^= t;
        Product[p[2]] ^= t;
        Product[p[3]] ^= t;
        Product[p[4]] ^= t;
        Product[p[5]] ^= t;
        Product[p[6]] ^= t;
        Product[p[7]] ^= t;
        Product[p[8]] ^= t;
        Product[p[9]] ^= t;
        Product[p[10]] ^= t;
        Product[p[11]] ^= t;
        Product[p[12]] ^= t;
        Product[p[13]] ^= t;
        Product[p[14]] ^= t;
        Product[p[15]] ^= t;
      }
      s = M->cEntry + M->cIndex[i + 1];
      for (; p < s; p++) {
        Product[p[0]] ^= t;
      }
    }
  }
#endif
}

void MultB_T64(u64 *Product, u64 *x, void *P) {
  nfs_sparse_mat_t *M = (nfs_sparse_mat_t *)P;
  memset(Product, 0, M->numCols * sizeof(u64));
  {
    int i;
    for (i = 0; i < M->numDenseBlocks; i++) {
      addmultnx64(Product, M->denseBlocks[i], x + M->denseBlockIndex[i], M->numCols);
    }
  }
#if defined(L2_CACHE_SIZE) && (L2_CACHE_SIZE > 0)
// L2_CACHE_SIZE has to be a power of 2.
// MULTB_T64_PAGESIZE is a half of L2 cache size.
#define MULTB_T64_PAGESIZE (L2_CACHE_SIZE * 1024 / 2 / sizeof(u64))
#define MULTB_T64_PAGEMASK (-MULTB_T64_PAGESIZE)
  {
    s32 n = M->numCols;
    u32 pagestart;
    for (pagestart = 0; pagestart < n; pagestart += MULTB_T64_PAGESIZE) {
      u32 *p = M->cEntry;
      s32 i;
      for (i = 0; i < n; i++) {
        u64 t = Product[i];
        u32 *s = p + ((M->cIndex[i + 1] - M->cIndex[i]) & -16);
        for (; p < s; p += 16) {
          if ((p[0] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[0]];
          if ((p[1] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[1]];
          if ((p[2] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[2]];
          if ((p[3] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[3]];
          if ((p[4] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[4]];
          if ((p[5] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[5]];
          if ((p[6] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[6]];
          if ((p[7] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[7]];
          if ((p[8] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[8]];
          if ((p[9] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[9]];
          if ((p[10] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[10]];
          if ((p[11] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[11]];
          if ((p[12] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[12]];
          if ((p[13] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[13]];
          if ((p[14] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[14]];
          if ((p[15] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[15]];
        }
        s = M->cEntry + M->cIndex[i + 1];
        for (; p < s; p++) {
          if ((p[0] & MULTB_T64_PAGEMASK) == pagestart) t ^= x[p[0]];
        }
        Product[i] = t;
      }
    }
  }
#else
  {
    s32 n = M->numCols;
    u32 *p = M->cEntry;
    s32 i;
    for (i = 0; i < n; i++) {
      u64 t = Product[i];
      u32 *s = p + ((M->cIndex[i + 1] - M->cIndex[i]) & -16);
      for (; p < s; p += 16) {
        t ^= x[p[0]];
        t ^= x[p[1]];
        t ^= x[p[2]];
        t ^= x[p[3]];
        t ^= x[p[4]];
        t ^= x[p[5]];
        t ^= x[p[6]];
        t ^= x[p[7]];
        t ^= x[p[8]];
        t ^= x[p[9]];
        t ^= x[p[10]];
        t ^= x[p[11]];
        t ^= x[p[12]];
        t ^= x[p[13]];
        t ^= x[p[14]];
        t ^= x[p[15]];
      }
      s = M->cEntry + M->cIndex[i + 1];
      for (; p < s; p++) {
        t ^= x[p[0]];
      }
      Product[i] = t;
    }
  }
#endif
}


/**********************************************************************/
int blockLanczos64(u64 *deps, MAT_MULT_FUNC_PTR64 MultB, 
                  MAT_MULT_FUNC_PTR64 MultB_T, void *P, s32 n)
/**********************************************************************/
{ u64 *Y=NULL, *X=NULL, *Vi=NULL, *Vi_1=NULL, *Vi_2=NULL, *tmp_n=NULL, *tmp2_n=NULL;
  u64 *V0=NULL, *tmp3_n=NULL, *Z=NULL, *AZ=NULL;
  u64 D[64] ALIGNED16, E[64] ALIGNED16, F[64] ALIGNED16, Wi[64] ALIGNED16;
  u64 Wi_1[64] ALIGNED16, Wi_2[64] ALIGNED16, T[64] ALIGNED16, T_1[64] ALIGNED16;
  u64 tmp[64] ALIGNED16;
  u64 U[64] ALIGNED16, U_1[64] ALIGNED16, tmp2[64] ALIGNED16;
  int  Si[64], Si_1[64];
  u64 i, j, m, mask, isZero, r1,r2;
  u32  iterations;
  int  errs=0, numDeps=-1, cont, s;

  
  /* Memory allocation: */
  if (!(Y = (u64 *)malloc(n*sizeof(u64))))    errs++;
  if (!(X = (u64 *)malloc(n*sizeof(u64))))    errs++;
  if (!(Vi = (u64 *)malloc(n*sizeof(u64))))   errs++;
  if (!(V0 = (u64 *)malloc(n*sizeof(u64))))   errs++;
  if (!(Vi_1 = (u64 *)malloc(n*sizeof(u64)))) errs++;
  if (!(Vi_2 = (u64 *)malloc(n*sizeof(u64)))) errs++;
  if (!(tmp_n = (u64 *)malloc(n*sizeof(u64)))) errs++;
  if (!(tmp2_n = (u64 *)malloc(n*sizeof(u64)))) errs++;
  if (!(tmp3_n = (u64 *)malloc(n*sizeof(u64)))) errs++;

  if (errs)
    goto SHORT_CIRC_STOP;
  
  /******************************************************************/
  /* Throughout, 'A' means the matrix A := (B^T)B. In fact, all the */
  /* notation is directly from Montgomery's paper, except that I    */
  /* don't bother indexing the D_i, E_i, F_i. They are simply D,E,F.*/
  /******************************************************************/
  /* Initialization: */
  for (j=0; j<n; j++) {
    r1 = prand(); r2 = prand(); 
    Y[j] = r1^(r2<<32);
    X[j] = Vi_1[j] = Vi_2[j] = tmp_n[j] = tmp2_n[j] = tmp3_n[j] = 0;
  }
  for (i=0; i<64; i++) {
    Wi[i] = Wi_1[i] = Wi_2[i] = 0;
    T[i] = T_1[i] = U[i] = U_1[i] = tmp[i] = 0;
    Si[i] = Si_1[i] = i; /* Si and Si_1 are both I_64. */
  }
  MultB(tmp_n, Y, P); 
  MultB_T(V0, tmp_n, P);
  memcpy(Vi, V0, n*sizeof(u64));

  /* Prime 'T' for the loop, so that we always have */
  /* T = (Vi^T)A(Vi) and T_1 = (Vi_1^T)A(Vi_1).     */
  
  MultB(tmp2_n, Vi, P);
  MultB_T(tmp_n, tmp2_n, P);  /* tmp_n <-- A*Vi */
  multT(T, Vi, tmp_n, n);      /* T <-- (Vi^T)(tmp_n) = (Vi^T)A(Vi) */
  
  cont = 1;
  i = 0;
  getW_S(Wi, Si, T, Si_1); /* Compute W0 and S0. */
  /* Initialize X <-- (V0)(W0)(V0^T)(V0). */
  multT(tmp, V0, V0, n);    /* tmp <-- (V0^T)(V0). */
  mult64x64(tmp2, Wi, tmp); /* tmp2 <-- (W0)(tmp) = (W0)(V0^T)(V0). */
  multnx64(X, V0, tmp2, n); /* X <-- V0(tmp2). */
  iterations = 0;
  do {
    /* Iteration step. */
    iterations++;

    /********** Compute D_{i+1}. **********/
    /* tmp_n = A*Vi from initialization, or the previous iteration. */

    multT(U, tmp_n, tmp_n, n);  /* U <-- (tmp_n)^T(tmp_n) = (Vi^T)(A^2)(Vi) */
    multS(U, Si);            /* U <-- (Vi^T)(A^2)(Vi)(Si)(Si^T). */
    memcpy(D, U, 64*sizeof(u64)); /* D <-- U. */
    for (j=0; j<64; j++)
      D[j] ^= T[j]; /* D <-- D + (Vi^T)A(Vi). */
    preMult(D, Wi); /* D <-- (Wi)D. */
    for (j=0; j<64; j++)
      D[j] ^= BIT64(j); /* D <-- D + I_{64}. */

    /********** Compute E_{i+1}. **********/
    mult64x64(E, Wi_1, T); /* E <-- (Wi_1)(Vi^T)A(Vi). */
    multS(E, Si);          /* E <-- E(Si)(Si^T).       */

    /**************** Compute F_{i+1}. *************/
    /* Note: We should, at this point, always have */
    /* T_1 = (Vi_1^T)A(Vi_1) and                   */
    /* U_1 = (Vi_1^T)A^2(Vi_1)(Si_1)(Si_1^T)       */
    /***********************************************/
    for (j=0; j<64; j++)
      F[j] = U_1[j] ^ T_1[j];
    multS(F, Si);

    mult64x64(tmp, T_1, Wi_1); /* tmp <-- (Vi_1^T)A(Vi_1)(Wi_1). */
    for (j=0; j<64; j++)
      tmp[j] ^= BIT64(j); /* tmp <-- tmp + I_64. */
    preMult(tmp, Wi_2);
    preMult(F, tmp);
    /*** Done computing 'F'. */

    /* Finally, do the actual iteration step, putting V_{i+1} in 'tmp_n'. */
    /* Right now, we still have tmp_n = A*Vi, so we should first use it.  */
    
    /* Do tmp_n <-- tmp_n*(Si)(Si^T). */
    mask = 0x00000000;
    for (j=0; j<64; j++) {
      s = Si[j];
      if ((s>=0) && (s<64))
        mask |= BIT64(s);
    }
    for (j=0; j<n; j++) 
      tmp_n[j] &= mask;
    
    addmultnx64(tmp_n, Vi, D, n);   /* tmp_n <-- tmp_n + (Vi)D.        */

    addmultnx64(tmp_n, Vi_1, E, n); /* tmp_n <-- tmp_n + (Vi_1)E */

    addmultnx64(tmp_n, Vi_2, F, n); /* tmp_n <-- tmp_n + (Vi_2)F */
    /*** Done. We now have tmp_n = V_{i+1}. ***/
   

    i++;
    memcpy(Vi_2, Vi_1, n*sizeof(u64));
    memcpy(Vi_1, Vi, n*sizeof(u64));
    memcpy(Vi, tmp_n, n*sizeof(u64));
    memcpy(Wi_2, Wi_1, 64*sizeof(u64));
    memcpy(Wi_1, Wi, 64*sizeof(u64));
    memcpy(T_1, T, 64*sizeof(u64));
    memcpy(U_1, U, 64*sizeof(u64));
    memcpy(Si_1, Si, 64*sizeof(int));
    
    /******** My `step 3'. ********/
    /* Is (Vi^T)(A)(Vi) == 0 ? */
    MultB(tmp2_n, Vi, P);
    MultB_T(tmp_n, tmp2_n, P);  /* tmp_n <-- A*Vi       */

    multT(T, Vi, tmp_n, n); /* T <-- (Vi^T)A(Vi). */
    getW_S(Wi, Si, T, Si_1);
    
    if (!(isZeroV(T, 64))) {
      /* X <-- X + (Vi)(Wi)(Vi^T)(V_0) */
      multT(tmp, Vi, V0, n); /* tmp <-- (Vi^T)(V0). */
      preMult(tmp, Wi);   /* tmp <-- (Wi)(tmp).  */
      
      addmultnx64(X, Vi, tmp, n);   /* X <-- X + (Vi)(tmp)    */
    } else {
      cont=0;
      m = i;
    }
    printTmp("Lanczos: Estimate %1.1lf%% complete...",
              (double)100.0*64.0*iterations/n);  
    if ((double)100.0*64.0*iterations/n > 250) {
      fprintf(stderr, "Some error has occurred: Lanczos is not converging!\n");
      fprintf(stderr, "Number of iterations is %ld.\n", iterations);
      /* Add some debugging stuff here! */
      fprintf(stderr, "Terminating...\n");
      exit(-1);
    }
  } while (cont);
  printf("\nBlock Lanczos used %ld iterations.\n", iterations);

          
  Z = (u64 *)malloc(2*n*sizeof(u64));
  AZ = (u64 *)malloc(2*n*sizeof(u64));
  if (!(Z&&AZ)) {
    fprintf(stderr, "blanczos(): Memory allocation error!\n");
    goto SHORT_CIRC_STOP;
  }
  
  if (isZeroV(Vi, n)) {
    /* Then <X+Y> < ker(A). Later, we will expect AX=0, so do X <-- X+Y. */
    printf("After Block Lanczos iteration, Vm=0.\n");

    for (i=0; i<n; i++)
      X[i] ^= Y[i];
  } else {
    printf("After Block Lanczos iteration, Vm is nonzero. Finishing...\n");
    /* We need more memory (for convenience), so free off what we don't need. */
    free(V0); free(Vi_1); free(Vi_2); free(tmp3_n);
    V0 = Vi_1 = Vi_2 = tmp3_n = NULL;

    /* Construct Z=[ X+Y | Vi] and compute AZ=[ A(X+Y) | AVi ] */
    /* X <-- X+Y, for convenience.                     */
    for (j=0; j<n; j++) {
      X[j] ^= Y[j];
      Z[2*j] = X[j];
      Z[2*j+1] = Vi[j];
    }
    /* Now, compute AZ = [ AX | AVi ]. */
    MultB(tmp2_n, X, P);  
    MultB_T(tmp_n, tmp2_n, P); /* tmp_n <-- AX       */
    for (j=0; j<n; j++) 
      AZ[2*j] = tmp_n[j];
    MultB(tmp2_n, Vi, P);  
    MultB_T(tmp_n, tmp2_n, P); /* tmp_n <-- AVi       */
    for (j=0; j<n; j++) 
      AZ[2*j+1] = tmp_n[j];
    /* Now, AZ should have small rank, and we should be  */
    /* able to do column ops to make zero vectors. Doing */
    /* the same ops to 'Z', we get vectors in ker(A).    */
    doColumnOps(Z, AZ, n);
    /* Now, look for zero columns in AZ. If a column is zero, */
    /* copy the corresponding column of Z into 'X'.           */
    for (i=0; i<n; i++)
      X[i] = 0;
    numDeps=0;
    for (i=0; (i<64) && (numDeps < 64); i++) {
      j=0;
      while ((j<n) && ((AZ[2*j + i/64]&BIT64(i%64))==0))
        j++;
      if (j==n) {
        /* Copy column i of Z into column 'numDeps' of X. */
        for (j=0; j<n; j++) {
          if (Z[2*j + i/64]&BIT64(i%64))
            X[j] ^= BIT64(numDeps);
        }
        numDeps++;
      }
    }
    printf("Found %d dependencies for A=(B^T)B.\n", numDeps);
  }
  j=0;
  while ((j<n)&&(X[j]==0))
    j++;
  if (j==n)
    printf("Probable error: The matrix X is identically zero!!!\n");


  printf("Getting dependencies for original matrix, B...\n");
  /***************************************************/
  /* At this point, we should have AX = (B^T)BX = 0. */
  /* Do simultaneous column ops to BX and X, putting */
  /* BX in RREF, to get vectors in ker(B).           */
  /***************************************************/

  for (i=0; i<n; i++)
    deps[i] = 0;
  numDeps=0;

  MultB(tmp_n, X, P);  /* tmp_n <-- BX. */
  doColumnOps64(X, tmp_n, n);
  /* We only want 32 of the dependencies. */
  for (i=0; i<32; i++) {
    for (j=0, isZero=1; j<n; j++)
      if (tmp_n[j]&BIT64(i)) {
        isZero=0;
      }
    if (isZero) {
      for (j=0, isZero=1; j<n; j++) {
        if (X[j]&BIT64(i)) {
          deps[j] ^= BIT64(numDeps);
          isZero=0;
        }
      }
      if (!(isZero))
        numDeps++;
    }
  }

  if (numDeps) {
    printf("Found %d dependencies for 'B'. Verifying...\n", numDeps);
  } else {
    printf("Some error occurred: all dependencies found seem to be trivial!\n");
    printf("Is Rank(A) = Rank((B^T)B) too small?\n");
    goto SHORT_CIRC_STOP;
  }

  MultB(tmp_n, deps, P);
  for (i=0, isZero=1; (i<n)&&(isZero); i++)
    if (tmp_n[i])
      isZero = 0;
  if (!(isZero))
    printf("Some error occurred: Final product (B)(deps) is nonzero (i=%ld)!\n", (s32)i);
  else
    printf("Verified.\n"); 



SHORT_CIRC_STOP:
  
  if (Y != NULL)      free(Y);
  if (X != NULL)      free(X);
  if (Vi != NULL)     free(Vi);
  if (V0 != NULL)     free(V0);
  if (Vi_1 != NULL)   free(Vi_1);
  if (Vi_2 != NULL)   free(Vi_2);
  if (tmp_n != NULL)  free(tmp_n);
  if (tmp2_n != NULL) free(tmp2_n);
  if (tmp3_n != NULL) free(tmp3_n);
  if (Z != NULL)      free(Z);
  if (AZ != NULL)     free(AZ);
  return numDeps;
}


u64 mult_w[2048] ALIGNED16;

void multT(u64 *c, u64 *a, u64 *b, s32 n) {
  memset(mult_w, 0, sizeof(u64) * 256 * 8);
  {
    int i;
    for (i = 0; i < n; i++) {
      u64 u = b[i];
      u32 t = (u32)a[i];
      mult_w[t & 255] ^= u;
      mult_w[256 + ((t >> 8) & 255)] ^= u;
      mult_w[512 + ((t >> 16) & 255)] ^= u;
      mult_w[768 + (t >> 24)] ^= u;
      t = (u32)(a[i] >> 32);
      mult_w[1024 + (t & 255)] ^= u;
      mult_w[1280 + ((t >> 8) & 255)] ^= u;
      mult_w[1536 + ((t >> 16) & 255)] ^= u;
      mult_w[1792 + (t >> 24)] ^= u;
    }
  }
  {
    int i;
    for (i = 0; i < 64; i += 8) {
      u64 *a = mult_w + (i << 5), *b = c + i;
      u64 r0, r1, r2, r3, r4;
	r0 = a[255];	//255->254
	r0 ^= a[254];
	a[254] = r0;
	r0 ^= a[252];	//254->252
	r0 ^= a[253];	//253->252
	a[252] = r0;
	r1 = a[251];	//251->250
	r1 ^= a[250];
	a[250] = r1;
	r0 ^= a[248];	//252->248
	r0 ^= r1;	//250->248
	r0 ^= a[249];	//249->248
	a[248] = r0;
	r1 = a[247];	//247->246
	r1 ^= a[246];
	a[246] = r1;
	r1 ^= a[244];	//246->244
	r1 ^= a[245];	//245->244
	a[244] = r1;
	r2 = a[243];	//243->242
	r2 ^= a[242];
	a[242] = r2;
	r0 ^= a[240];	//248->240
	r0 ^= r1;	//244->240
	r0 ^= r2;	//242->240
	r0 ^= a[241];	//241->240
	a[240] = r0;
	r1 = a[239];	//239->238
	r1 ^= a[238];
	a[238] = r1;
	r1 ^= a[236];	//238->236
	r1 ^= a[237];	//237->236
	a[236] = r1;
	r2 = a[235];	//235->234
	r2 ^= a[234];
	a[234] = r2;
	r1 ^= a[232];	//236->232
	r1 ^= r2;	//234->232
	r1 ^= a[233];	//233->232
	a[232] = r1;
	r2 = a[231];	//231->230
	r2 ^= a[230];
	a[230] = r2;
	r2 ^= a[228];	//230->228
	r2 ^= a[229];	//229->228
	a[228] = r2;
	r3 = a[227];	//227->226
	r3 ^= a[226];
	a[226] = r3;
	r0 ^= a[224];	//240->224
	r0 ^= r1;	//232->224
	r0 ^= r2;	//228->224
	r0 ^= r3;	//226->224
	r0 ^= a[225];	//225->224
	a[224] = r0;
	r1 = a[223];	//223->222
	r1 ^= a[222];
	a[222] = r1;
	r1 ^= a[220];	//222->220
	r1 ^= a[221];	//221->220
	a[220] = r1;
	r2 = a[219];	//219->218
	r2 ^= a[218];
	a[218] = r2;
	r1 ^= a[216];	//220->216
	r1 ^= r2;	//218->216
	r1 ^= a[217];	//217->216
	a[216] = r1;
	r2 = a[215];	//215->214
	r2 ^= a[214];
	a[214] = r2;
	r2 ^= a[212];	//214->212
	r2 ^= a[213];	//213->212
	a[212] = r2;
	r3 = a[211];	//211->210
	r3 ^= a[210];
	a[210] = r3;
	r1 ^= a[208];	//216->208
	r1 ^= r2;	//212->208
	r1 ^= r3;	//210->208
	r1 ^= a[209];	//209->208
	a[208] = r1;
	r2 = a[207];	//207->206
	r2 ^= a[206];
	a[206] = r2;
	r2 ^= a[204];	//206->204
	r2 ^= a[205];	//205->204
	a[204] = r2;
	r3 = a[203];	//203->202
	r3 ^= a[202];
	a[202] = r3;
	r2 ^= a[200];	//204->200
	r2 ^= r3;	//202->200
	r2 ^= a[201];	//201->200
	a[200] = r2;
	r3 = a[199];	//199->198
	r3 ^= a[198];
	a[198] = r3;
	r3 ^= a[196];	//198->196
	r3 ^= a[197];	//197->196
	a[196] = r3;
	r4 = a[195];	//195->194
	r4 ^= a[194];
	a[194] = r4;
	r0 ^= a[192];	//224->192
	r0 ^= r1;	//208->192
	r0 ^= r2;	//200->192
	r0 ^= r3;	//196->192
	r0 ^= r4;	//194->192
	a[192] = r0;
	r0 ^= a[193];	//193->192
	a[192] = r0;
	r0 = a[191];	//191->190
	r0 ^= a[190];
	a[190] = r0;
	r0 ^= a[188];	//190->188
	r0 ^= a[189];	//189->188
	a[188] = r0;
	r1 = a[187];	//187->186
	r1 ^= a[186];
	a[186] = r1;
	r0 ^= a[184];	//188->184
	r0 ^= r1;	//186->184
	r0 ^= a[185];	//185->184
	a[184] = r0;
	r1 = a[183];	//183->182
	r1 ^= a[182];
	a[182] = r1;
	r1 ^= a[180];	//182->180
	r1 ^= a[181];	//181->180
	a[180] = r1;
	r2 = a[179];	//179->178
	r2 ^= a[178];
	a[178] = r2;
	r0 ^= a[176];	//184->176
	r0 ^= r1;	//180->176
	r0 ^= r2;	//178->176
	r0 ^= a[177];	//177->176
	a[176] = r0;
	r1 = a[175];	//175->174
	r1 ^= a[174];
	a[174] = r1;
	r1 ^= a[172];	//174->172
	r1 ^= a[173];	//173->172
	a[172] = r1;
	r2 = a[171];	//171->170
	r2 ^= a[170];
	a[170] = r2;
	r1 ^= a[168];	//172->168
	r1 ^= r2;	//170->168
	r1 ^= a[169];	//169->168
	a[168] = r1;
	r2 = a[167];	//167->166
	r2 ^= a[166];
	a[166] = r2;
	r2 ^= a[164];	//166->164
	r2 ^= a[165];	//165->164
	a[164] = r2;
	r3 = a[163];	//163->162
	r3 ^= a[162];
	a[162] = r3;
	r0 ^= a[160];	//176->160
	r0 ^= r1;	//168->160
	r0 ^= r2;	//164->160
	r0 ^= r3;	//162->160
	a[160] = r0;
	r0 ^= a[161];	//161->160
	a[160] = r0;
	r0 = a[159];	//159->158
	r0 ^= a[158];
	a[158] = r0;
	r0 ^= a[156];	//158->156
	r0 ^= a[157];	//157->156
	a[156] = r0;
	r1 = a[155];	//155->154
	r1 ^= a[154];
	a[154] = r1;
	r0 ^= a[152];	//156->152
	r0 ^= r1;	//154->152
	r0 ^= a[153];	//153->152
	a[152] = r0;
	r1 = a[151];	//151->150
	r1 ^= a[150];
	a[150] = r1;
	r1 ^= a[148];	//150->148
	r1 ^= a[149];	//149->148
	a[148] = r1;
	r2 = a[147];	//147->146
	r2 ^= a[146];
	a[146] = r2;
	r0 ^= a[144];	//152->144
	r0 ^= r1;	//148->144
	r0 ^= r2;	//146->144
	a[144] = r0;
	r0 ^= a[145];	//145->144
	a[144] = r0;
	r0 = a[143];	//143->142
	r0 ^= a[142];
	a[142] = r0;
	r0 ^= a[140];	//142->140
	r0 ^= a[141];	//141->140
	a[140] = r0;
	r1 = a[139];	//139->138
	r1 ^= a[138];
	a[138] = r1;
	r0 ^= a[136];	//140->136
	r0 ^= r1;	//138->136
	a[136] = r0;
	r0 ^= a[137];	//137->136
	a[136] = r0;
	r0 = a[135];	//135->134
	r0 ^= a[134];
	a[134] = r0;
	r0 ^= a[132];	//134->132
	a[132] = r0;
	r0 ^= a[133];	//133->132
	a[132] = r0;
	r0 = a[131];	//131->130
	r0 ^= a[130];
	a[130] = r0;
	r0 = a[127];	//127->126
	r0 ^= a[126];
	a[126] = r0;
	r0 ^= a[124];	//126->124
	r0 ^= a[125];	//125->124
	a[124] = r0;
	r1 = a[123];	//123->122
	r1 ^= a[122];
	a[122] = r1;
	r0 ^= a[120];	//124->120
	r0 ^= r1;	//122->120
	r0 ^= a[121];	//121->120
	a[120] = r0;
	r1 = a[119];	//119->118
	r1 ^= a[118];
	a[118] = r1;
	r1 ^= a[116];	//118->116
	r1 ^= a[117];	//117->116
	a[116] = r1;
	r2 = a[115];	//115->114
	r2 ^= a[114];
	a[114] = r2;
	r0 ^= a[112];	//120->112
	r0 ^= r1;	//116->112
	r0 ^= r2;	//114->112
	r0 ^= a[113];	//113->112
	a[112] = r0;
	r1 = a[111];	//111->110
	r1 ^= a[110];
	a[110] = r1;
	r1 ^= a[108];	//110->108
	r1 ^= a[109];	//109->108
	a[108] = r1;
	r2 = a[107];	//107->106
	r2 ^= a[106];
	a[106] = r2;
	r1 ^= a[104];	//108->104
	r1 ^= r2;	//106->104
	r1 ^= a[105];	//105->104
	a[104] = r1;
	r2 = a[103];	//103->102
	r2 ^= a[102];
	a[102] = r2;
	r2 ^= a[100];	//102->100
	r2 ^= a[101];	//101->100
	a[100] = r2;
	r3 = a[99];	//99->98
	r3 ^= a[98];
	a[98] = r3;
	r0 ^= a[96];	//112->96
	r0 ^= r1;	//104->96
	r0 ^= r2;	//100->96
	r0 ^= r3;	//98->96
	a[96] = r0;
	r0 ^= a[97];	//97->96
	a[96] = r0;
	r0 = a[95];	//95->94
	r0 ^= a[94];
	a[94] = r0;
	r0 ^= a[92];	//94->92
	r0 ^= a[93];	//93->92
	a[92] = r0;
	r1 = a[91];	//91->90
	r1 ^= a[90];
	a[90] = r1;
	r0 ^= a[88];	//92->88
	r0 ^= r1;	//90->88
	r0 ^= a[89];	//89->88
	a[88] = r0;
	r1 = a[87];	//87->86
	r1 ^= a[86];
	a[86] = r1;
	r1 ^= a[84];	//86->84
	r1 ^= a[85];	//85->84
	a[84] = r1;
	r2 = a[83];	//83->82
	r2 ^= a[82];
	a[82] = r2;
	r0 ^= a[80];	//88->80
	r0 ^= r1;	//84->80
	r0 ^= r2;	//82->80
	a[80] = r0;
	r0 ^= a[81];	//81->80
	a[80] = r0;
	r0 = a[79];	//79->78
	r0 ^= a[78];
	a[78] = r0;
	r0 ^= a[76];	//78->76
	r0 ^= a[77];	//77->76
	a[76] = r0;
	r1 = a[75];	//75->74
	r1 ^= a[74];
	a[74] = r1;
	r0 ^= a[72];	//76->72
	r0 ^= r1;	//74->72
	a[72] = r0;
	r0 ^= a[73];	//73->72
	a[72] = r0;
	r0 = a[71];	//71->70
	r0 ^= a[70];
	a[70] = r0;
	r0 ^= a[68];	//70->68
	a[68] = r0;
	r0 ^= a[69];	//69->68
	a[68] = r0;
	r0 = a[67];	//67->66
	r0 ^= a[66];
	a[66] = r0;
	r0 = a[63];	//63->62
	r0 ^= a[62];
	a[62] = r0;
	r0 ^= a[60];	//62->60
	r0 ^= a[61];	//61->60
	a[60] = r0;
	r1 = a[59];	//59->58
	r1 ^= a[58];
	a[58] = r1;
	r0 ^= a[56];	//60->56
	r0 ^= r1;	//58->56
	r0 ^= a[57];	//57->56
	a[56] = r0;
	r1 = a[55];	//55->54
	r1 ^= a[54];
	a[54] = r1;
	r1 ^= a[52];	//54->52
	r1 ^= a[53];	//53->52
	a[52] = r1;
	r2 = a[51];	//51->50
	r2 ^= a[50];
	a[50] = r2;
	r0 ^= a[48];	//56->48
	r0 ^= r1;	//52->48
	r0 ^= r2;	//50->48
	a[48] = r0;
	r0 ^= a[49];	//49->48
	a[48] = r0;
	r0 = a[47];	//47->46
	r0 ^= a[46];
	a[46] = r0;
	r0 ^= a[44];	//46->44
	r0 ^= a[45];	//45->44
	a[44] = r0;
	r1 = a[43];	//43->42
	r1 ^= a[42];
	a[42] = r1;
	r0 ^= a[40];	//44->40
	r0 ^= r1;	//42->40
	a[40] = r0;
	r0 ^= a[41];	//41->40
	a[40] = r0;
	r0 = a[39];	//39->38
	r0 ^= a[38];
	a[38] = r0;
	r0 ^= a[36];	//38->36
	a[36] = r0;
	r0 ^= a[37];	//37->36
	a[36] = r0;
	r0 = a[35];	//35->34
	r0 ^= a[34];
	a[34] = r0;
	r0 = a[31];	//31->30
	r0 ^= a[30];
	a[30] = r0;
	r0 ^= a[28];	//30->28
	r0 ^= a[29];	//29->28
	a[28] = r0;
	r1 = a[27];	//27->26
	r1 ^= a[26];
	a[26] = r1;
	r0 ^= a[24];	//28->24
	r0 ^= r1;	//26->24
	a[24] = r0;
	r0 ^= a[25];	//25->24
	a[24] = r0;
	r0 = a[23];	//23->22
	r0 ^= a[22];
	a[22] = r0;
	r0 ^= a[20];	//22->20
	a[20] = r0;
	r0 ^= a[21];	//21->20
	a[20] = r0;
	r0 = a[19];	//19->18
	r0 ^= a[18];
	a[18] = r0;
	r0 = a[15];	//15->14
	r0 ^= a[14];
	a[14] = r0;
	r0 ^= a[12];	//14->12
	a[12] = r0;
	r0 ^= a[13];	//13->12
	a[12] = r0;
	r0 = a[11];	//11->10
	r0 ^= a[10];
	a[10] = r0;
	r0 = a[7];	//7->6
	r0 ^= a[6];
	a[6] = r0;
	r0 = a[1];
	r0 ^= a[3];
	r0 ^= a[5];
	r0 ^= a[7];
	r0 ^= a[9];
	r0 ^= a[11];
	r0 ^= a[13];
	r0 ^= a[15];
	r0 ^= a[17];
	r0 ^= a[19];
	r0 ^= a[21];
	r0 ^= a[23];
	r0 ^= a[25];
	r0 ^= a[27];
	r0 ^= a[29];
	r0 ^= a[31];
	r0 ^= a[33];
	r0 ^= a[35];
	r0 ^= a[37];
	r0 ^= a[39];
	r0 ^= a[41];
	r0 ^= a[43];
	r0 ^= a[45];
	r0 ^= a[47];
	r0 ^= a[49];
	r0 ^= a[51];
	r0 ^= a[53];
	r0 ^= a[55];
	r0 ^= a[57];
	r0 ^= a[59];
	r0 ^= a[61];
	r0 ^= a[63];
	r0 ^= a[65];
	r0 ^= a[67];
	r0 ^= a[69];
	r0 ^= a[71];
	r0 ^= a[73];
	r0 ^= a[75];
	r0 ^= a[77];
	r0 ^= a[79];
	r0 ^= a[81];
	r0 ^= a[83];
	r0 ^= a[85];
	r0 ^= a[87];
	r0 ^= a[89];
	r0 ^= a[91];
	r0 ^= a[93];
	r0 ^= a[95];
	r0 ^= a[97];
	r0 ^= a[99];
	r0 ^= a[101];
	r0 ^= a[103];
	r0 ^= a[105];
	r0 ^= a[107];
	r0 ^= a[109];
	r0 ^= a[111];
	r0 ^= a[113];
	r0 ^= a[115];
	r0 ^= a[117];
	r0 ^= a[119];
	r0 ^= a[121];
	r0 ^= a[123];
	r0 ^= a[125];
	r0 ^= a[127];
	r0 ^= a[129];
	r0 ^= a[131];
	r0 ^= a[133];
	r0 ^= a[135];
	r0 ^= a[137];
	r0 ^= a[139];
	r0 ^= a[141];
	r0 ^= a[143];
	r0 ^= a[145];
	r0 ^= a[147];
	r0 ^= a[149];
	r0 ^= a[151];
	r0 ^= a[153];
	r0 ^= a[155];
	r0 ^= a[157];
	r0 ^= a[159];
	r0 ^= a[161];
	r0 ^= a[163];
	r0 ^= a[165];
	r0 ^= a[167];
	r0 ^= a[169];
	r0 ^= a[171];
	r0 ^= a[173];
	r0 ^= a[175];
	r0 ^= a[177];
	r0 ^= a[179];
	r0 ^= a[181];
	r0 ^= a[183];
	r0 ^= a[185];
	r0 ^= a[187];
	r0 ^= a[189];
	r0 ^= a[191];
	r0 ^= a[193];
	r0 ^= a[195];
	r0 ^= a[197];
	r0 ^= a[199];
	r0 ^= a[201];
	r0 ^= a[203];
	r0 ^= a[205];
	r0 ^= a[207];
	r0 ^= a[209];
	r0 ^= a[211];
	r0 ^= a[213];
	r0 ^= a[215];
	r0 ^= a[217];
	r0 ^= a[219];
	r0 ^= a[221];
	r0 ^= a[223];
	r0 ^= a[225];
	r0 ^= a[227];
	r0 ^= a[229];
	r0 ^= a[231];
	r0 ^= a[233];
	r0 ^= a[235];
	r0 ^= a[237];
	r0 ^= a[239];
	r0 ^= a[241];
	r0 ^= a[243];
	r0 ^= a[245];
	r0 ^= a[247];
	r0 ^= a[249];
	r0 ^= a[251];
	r0 ^= a[253];
	r0 ^= a[255];
	b[0] = r0;
	r0 = a[2];
	r0 ^= a[3];
	r0 ^= a[6];
	r0 ^= a[10];
	r0 ^= a[14];
	r0 ^= a[18];
	r0 ^= a[22];
	r0 ^= a[26];
	r0 ^= a[30];
	r0 ^= a[34];
	r0 ^= a[38];
	r0 ^= a[42];
	r0 ^= a[46];
	r0 ^= a[50];
	r0 ^= a[54];
	r0 ^= a[58];
	r0 ^= a[62];
	r0 ^= a[66];
	r0 ^= a[70];
	r0 ^= a[74];
	r0 ^= a[78];
	r0 ^= a[82];
	r0 ^= a[86];
	r0 ^= a[90];
	r0 ^= a[94];
	r0 ^= a[98];
	r0 ^= a[102];
	r0 ^= a[106];
	r0 ^= a[110];
	r0 ^= a[114];
	r0 ^= a[118];
	r0 ^= a[122];
	r0 ^= a[126];
	r0 ^= a[130];
	r0 ^= a[134];
	r0 ^= a[138];
	r0 ^= a[142];
	r0 ^= a[146];
	r0 ^= a[150];
	r0 ^= a[154];
	r0 ^= a[158];
	r0 ^= a[162];
	r0 ^= a[166];
	r0 ^= a[170];
	r0 ^= a[174];
	r0 ^= a[178];
	r0 ^= a[182];
	r0 ^= a[186];
	r0 ^= a[190];
	r0 ^= a[194];
	r0 ^= a[198];
	r0 ^= a[202];
	r0 ^= a[206];
	r0 ^= a[210];
	r0 ^= a[214];
	r0 ^= a[218];
	r0 ^= a[222];
	r0 ^= a[226];
	r0 ^= a[230];
	r0 ^= a[234];
	r0 ^= a[238];
	r0 ^= a[242];
	r0 ^= a[246];
	r0 ^= a[250];
	r0 ^= a[254];
	b[1] = r0;
	r0 = a[4];
	r0 ^= a[5];
	r0 ^= a[6];
	r0 ^= a[12];
	r0 ^= a[20];
	r0 ^= a[28];
	r0 ^= a[36];
	r0 ^= a[44];
	r0 ^= a[52];
	r0 ^= a[60];
	r0 ^= a[68];
	r0 ^= a[76];
	r0 ^= a[84];
	r0 ^= a[92];
	r0 ^= a[100];
	r0 ^= a[108];
	r0 ^= a[116];
	r0 ^= a[124];
	r0 ^= a[132];
	r0 ^= a[140];
	r0 ^= a[148];
	r0 ^= a[156];
	r0 ^= a[164];
	r0 ^= a[172];
	r0 ^= a[180];
	r0 ^= a[188];
	r0 ^= a[196];
	r0 ^= a[204];
	r0 ^= a[212];
	r0 ^= a[220];
	r0 ^= a[228];
	r0 ^= a[236];
	r0 ^= a[244];
	r0 ^= a[252];
	b[2] = r0;
	r0 = a[8];
	r0 ^= a[9];
	r0 ^= a[10];
	r0 ^= a[12];
	r0 ^= a[24];
	r0 ^= a[40];
	r0 ^= a[56];
	r0 ^= a[72];
	r0 ^= a[88];
	r0 ^= a[104];
	r0 ^= a[120];
	r0 ^= a[136];
	r0 ^= a[152];
	r0 ^= a[168];
	r0 ^= a[184];
	r0 ^= a[200];
	r0 ^= a[216];
	r0 ^= a[232];
	r0 ^= a[248];
	b[3] = r0;
	r0 = a[16];
	r0 ^= a[17];
	r0 ^= a[18];
	r0 ^= a[20];
	r0 ^= a[24];
	r0 ^= a[48];
	r0 ^= a[80];
	r0 ^= a[112];
	r0 ^= a[144];
	r0 ^= a[176];
	r0 ^= a[208];
	r0 ^= a[240];
	b[4] = r0;
	r0 = a[32];
	r0 ^= a[33];
	r0 ^= a[34];
	r0 ^= a[36];
	r0 ^= a[40];
	r0 ^= a[48];
	r0 ^= a[96];
	r0 ^= a[160];
	r0 ^= a[224];
	b[5] = r0;
	r0 = a[64];
	r0 ^= a[65];
	r0 ^= a[66];
	r0 ^= a[68];
	r0 ^= a[72];
	r0 ^= a[80];
	r0 ^= a[96];
	r0 ^= a[192];
	b[6] = r0;
	r0 = a[128];
	r0 ^= a[129];
	r0 ^= a[130];
	r0 ^= a[132];
	r0 ^= a[136];
	r0 ^= a[144];
	r0 ^= a[160];
	r0 ^= a[192];
	b[7] = r0;
    }
  }
}

  
/*********************************************/
void multS(u64 *D, int *S)
/*********************************************/
/* Input: D is a 64x64 matrix,               */
/*        'S' is a subset of the columns,    */
/*        {0,1,2,...,31}.                    */
/* Output: Columns of 'D' not in 'S' are     */
/*         zeroed, and the others unchanged. */
/*********************************************/
{ int i, s;
  u64 mask;

  mask = 0x00000000;
  for (i=0; i<64; i++) {
    s = S[i];
    if ((s>=0) && (s<64))
      mask |= BIT64(s);
  }
  for (i=0; i<64; i++) 
    D[i] &= mask;
}

void mult64x64(u64 *c, u64 *a, u64 *b) {
  int i;
  for (i = 0; i < 64; i++) {
    u64 u = 0;
    u32 t = (u32)a[i];
    if (t & 1) u = b[0];
    if (t & (1 << 1)) u ^= b[1];
    if (t & (1 << 2)) u ^= b[2];
    if (t & (1 << 3)) u ^= b[3];
    if (t & (1 << 4)) u ^= b[4];
    if (t & (1 << 5)) u ^= b[5];
    if (t & (1 << 6)) u ^= b[6];
    if (t & (1 << 7)) u ^= b[7];
    if (t & (1 << 8)) u ^= b[8];
    if (t & (1 << 9)) u ^= b[9];
    if (t & (1 << 10)) u ^= b[10];
    if (t & (1 << 11)) u ^= b[11];
    if (t & (1 << 12)) u ^= b[12];
    if (t & (1 << 13)) u ^= b[13];
    if (t & (1 << 14)) u ^= b[14];
    if (t & (1 << 15)) u ^= b[15];
    if (t & (1 << 16)) u ^= b[16];
    if (t & (1 << 17)) u ^= b[17];
    if (t & (1 << 18)) u ^= b[18];
    if (t & (1 << 19)) u ^= b[19];
    if (t & (1 << 20)) u ^= b[20];
    if (t & (1 << 21)) u ^= b[21];
    if (t & (1 << 22)) u ^= b[22];
    if (t & (1 << 23)) u ^= b[23];
    if (t & (1 << 24)) u ^= b[24];
    if (t & (1 << 25)) u ^= b[25];
    if (t & (1 << 26)) u ^= b[26];
    if (t & (1 << 27)) u ^= b[27];
    if (t & (1 << 28)) u ^= b[28];
    if (t & (1 << 29)) u ^= b[29];
    if (t & (1 << 30)) u ^= b[30];
    if ((s32)t < 0) u ^= b[31];
    t = (u32)(a[i] >> 32);
    if (t & 1) u ^= b[32];
    if (t & (1 << 1)) u ^= b[33];
    if (t & (1 << 2)) u ^= b[34];
    if (t & (1 << 3)) u ^= b[35];
    if (t & (1 << 4)) u ^= b[36];
    if (t & (1 << 5)) u ^= b[37];
    if (t & (1 << 6)) u ^= b[38];
    if (t & (1 << 7)) u ^= b[39];
    if (t & (1 << 8)) u ^= b[40];
    if (t & (1 << 9)) u ^= b[41];
    if (t & (1 << 10)) u ^= b[42];
    if (t & (1 << 11)) u ^= b[43];
    if (t & (1 << 12)) u ^= b[44];
    if (t & (1 << 13)) u ^= b[45];
    if (t & (1 << 14)) u ^= b[46];
    if (t & (1 << 15)) u ^= b[47];
    if (t & (1 << 16)) u ^= b[48];
    if (t & (1 << 17)) u ^= b[49];
    if (t & (1 << 18)) u ^= b[50];
    if (t & (1 << 19)) u ^= b[51];
    if (t & (1 << 20)) u ^= b[52];
    if (t & (1 << 21)) u ^= b[53];
    if (t & (1 << 22)) u ^= b[54];
    if (t & (1 << 23)) u ^= b[55];
    if (t & (1 << 24)) u ^= b[56];
    if (t & (1 << 25)) u ^= b[57];
    if (t & (1 << 26)) u ^= b[58];
    if (t & (1 << 27)) u ^= b[59];
    if (t & (1 << 28)) u ^= b[60];
    if (t & (1 << 29)) u ^= b[61];
    if (t & (1 << 30)) u ^= b[62];
    if ((s32)t < 0) u ^= b[63];
    c[i] = u;
  }
}

/******************************************/
void preMult(u64 *A, u64 *B)
/******************************************/
/* Input: 'A' and 'B' are 64x64 matrices. */
/* Output: A <-- B*A.                     */
/******************************************/
{ u64 res[64] ALIGNED16;

  mult64x64(res, B, A);
  memcpy(A, res, 64*sizeof(u64));
}
  
void multnx64(u64 *c, u64 *a, u64 *b, s32 n) {
  int i, j;
  for (i = 0, j = 0; i < 64; i += 8, j += 256) {
    u64 b0 = b[i];
    u64 b1 = b[i + 1];
    u64 b2 = b[i + 2];
    u64 t;
    mult_w[j] = 0;
    mult_w[j + 1] = t = b0;
    mult_w[j + 3] = t ^= b1;
    mult_w[j + 2] = t ^= b0;
    mult_w[j + 6] = t ^= b2;
    mult_w[j + 7] = t ^= b0;
    mult_w[j + 5] = t ^= b1;
    mult_w[j + 4] = t ^= b0;
    mult_w[j + 12] = t ^= b[i + 3];
    mult_w[j + 13] = t ^= b0;
    mult_w[j + 15] = t ^= b1;
    mult_w[j + 14] = t ^= b0;
    mult_w[j + 10] = t ^= b2;
    mult_w[j + 11] = t ^= b0;
    mult_w[j + 9] = t ^= b1;
    mult_w[j + 8] = t ^= b0;
    mult_w[j + 24] = t ^= b[i + 4];
    mult_w[j + 25] = t ^= b0;
    mult_w[j + 27] = t ^= b1;
    mult_w[j + 26] = t ^= b0;
    mult_w[j + 30] = t ^= b2;
    mult_w[j + 31] = t ^= b0;
    mult_w[j + 29] = t ^= b1;
    mult_w[j + 28] = t ^= b0;
    mult_w[j + 20] = t ^= b[i + 3];
    mult_w[j + 21] = t ^= b0;
    mult_w[j + 23] = t ^= b1;
    mult_w[j + 22] = t ^= b0;
    mult_w[j + 18] = t ^= b2;
    mult_w[j + 19] = t ^= b0;
    mult_w[j + 17] = t ^= b1;
    mult_w[j + 16] = t ^= b0;
    mult_w[j + 48] = t ^= b[i + 5];
    mult_w[j + 49] = t ^= b0;
    mult_w[j + 51] = t ^= b1;
    mult_w[j + 50] = t ^= b0;
    mult_w[j + 54] = t ^= b2;
    mult_w[j + 55] = t ^= b0;
    mult_w[j + 53] = t ^= b1;
    mult_w[j + 52] = t ^= b0;
    mult_w[j + 60] = t ^= b[i + 3];
    mult_w[j + 61] = t ^= b0;
    mult_w[j + 63] = t ^= b1;
    mult_w[j + 62] = t ^= b0;
    mult_w[j + 58] = t ^= b2;
    mult_w[j + 59] = t ^= b0;
    mult_w[j + 57] = t ^= b1;
    mult_w[j + 56] = t ^= b0;
    mult_w[j + 40] = t ^= b[i + 4];
    mult_w[j + 41] = t ^= b0;
    mult_w[j + 43] = t ^= b1;
    mult_w[j + 42] = t ^= b0;
    mult_w[j + 46] = t ^= b2;
    mult_w[j + 47] = t ^= b0;
    mult_w[j + 45] = t ^= b1;
    mult_w[j + 44] = t ^= b0;
    mult_w[j + 36] = t ^= b[i + 3];
    mult_w[j + 37] = t ^= b0;
    mult_w[j + 39] = t ^= b1;
    mult_w[j + 38] = t ^= b0;
    mult_w[j + 34] = t ^= b2;
    mult_w[j + 35] = t ^= b0;
    mult_w[j + 33] = t ^= b1;
    mult_w[j + 32] = t ^= b0;
    mult_w[j + 96] = t ^= b[i + 6];
    mult_w[j + 97] = t ^= b0;
    mult_w[j + 99] = t ^= b1;
    mult_w[j + 98] = t ^= b0;
    mult_w[j + 102] = t ^= b2;
    mult_w[j + 103] = t ^= b0;
    mult_w[j + 101] = t ^= b1;
    mult_w[j + 100] = t ^= b0;
    mult_w[j + 108] = t ^= b[i + 3];
    mult_w[j + 109] = t ^= b0;
    mult_w[j + 111] = t ^= b1;
    mult_w[j + 110] = t ^= b0;
    mult_w[j + 106] = t ^= b2;
    mult_w[j + 107] = t ^= b0;
    mult_w[j + 105] = t ^= b1;
    mult_w[j + 104] = t ^= b0;
    mult_w[j + 120] = t ^= b[i + 4];
    mult_w[j + 121] = t ^= b0;
    mult_w[j + 123] = t ^= b1;
    mult_w[j + 122] = t ^= b0;
    mult_w[j + 126] = t ^= b2;
    mult_w[j + 127] = t ^= b0;
    mult_w[j + 125] = t ^= b1;
    mult_w[j + 124] = t ^= b0;
    mult_w[j + 116] = t ^= b[i + 3];
    mult_w[j + 117] = t ^= b0;
    mult_w[j + 119] = t ^= b1;
    mult_w[j + 118] = t ^= b0;
    mult_w[j + 114] = t ^= b2;
    mult_w[j + 115] = t ^= b0;
    mult_w[j + 113] = t ^= b1;
    mult_w[j + 112] = t ^= b0;
    mult_w[j + 80] = t ^= b[i + 5];
    mult_w[j + 81] = t ^= b0;
    mult_w[j + 83] = t ^= b1;
    mult_w[j + 82] = t ^= b0;
    mult_w[j + 86] = t ^= b2;
    mult_w[j + 87] = t ^= b0;
    mult_w[j + 85] = t ^= b1;
    mult_w[j + 84] = t ^= b0;
    mult_w[j + 92] = t ^= b[i + 3];
    mult_w[j + 93] = t ^= b0;
    mult_w[j + 95] = t ^= b1;
    mult_w[j + 94] = t ^= b0;
    mult_w[j + 90] = t ^= b2;
    mult_w[j + 91] = t ^= b0;
    mult_w[j + 89] = t ^= b1;
    mult_w[j + 88] = t ^= b0;
    mult_w[j + 72] = t ^= b[i + 4];
    mult_w[j + 73] = t ^= b0;
    mult_w[j + 75] = t ^= b1;
    mult_w[j + 74] = t ^= b0;
    mult_w[j + 78] = t ^= b2;
    mult_w[j + 79] = t ^= b0;
    mult_w[j + 77] = t ^= b1;
    mult_w[j + 76] = t ^= b0;
    mult_w[j + 68] = t ^= b[i + 3];
    mult_w[j + 69] = t ^= b0;
    mult_w[j + 71] = t ^= b1;
    mult_w[j + 70] = t ^= b0;
    mult_w[j + 66] = t ^= b2;
    mult_w[j + 67] = t ^= b0;
    mult_w[j + 65] = t ^= b1;
    mult_w[j + 64] = t ^= b0;
    mult_w[j + 192] = t ^= b[i + 7];
    mult_w[j + 193] = t ^= b0;
    mult_w[j + 195] = t ^= b1;
    mult_w[j + 194] = t ^= b0;
    mult_w[j + 198] = t ^= b2;
    mult_w[j + 199] = t ^= b0;
    mult_w[j + 197] = t ^= b1;
    mult_w[j + 196] = t ^= b0;
    mult_w[j + 204] = t ^= b[i + 3];
    mult_w[j + 205] = t ^= b0;
    mult_w[j + 207] = t ^= b1;
    mult_w[j + 206] = t ^= b0;
    mult_w[j + 202] = t ^= b2;
    mult_w[j + 203] = t ^= b0;
    mult_w[j + 201] = t ^= b1;
    mult_w[j + 200] = t ^= b0;
    mult_w[j + 216] = t ^= b[i + 4];
    mult_w[j + 217] = t ^= b0;
    mult_w[j + 219] = t ^= b1;
    mult_w[j + 218] = t ^= b0;
    mult_w[j + 222] = t ^= b2;
    mult_w[j + 223] = t ^= b0;
    mult_w[j + 221] = t ^= b1;
    mult_w[j + 220] = t ^= b0;
    mult_w[j + 212] = t ^= b[i + 3];
    mult_w[j + 213] = t ^= b0;
    mult_w[j + 215] = t ^= b1;
    mult_w[j + 214] = t ^= b0;
    mult_w[j + 210] = t ^= b2;
    mult_w[j + 211] = t ^= b0;
    mult_w[j + 209] = t ^= b1;
    mult_w[j + 208] = t ^= b0;
    mult_w[j + 240] = t ^= b[i + 5];
    mult_w[j + 241] = t ^= b0;
    mult_w[j + 243] = t ^= b1;
    mult_w[j + 242] = t ^= b0;
    mult_w[j + 246] = t ^= b2;
    mult_w[j + 247] = t ^= b0;
    mult_w[j + 245] = t ^= b1;
    mult_w[j + 244] = t ^= b0;
    mult_w[j + 252] = t ^= b[i + 3];
    mult_w[j + 253] = t ^= b0;
    mult_w[j + 255] = t ^= b1;
    mult_w[j + 254] = t ^= b0;
    mult_w[j + 250] = t ^= b2;
    mult_w[j + 251] = t ^= b0;
    mult_w[j + 249] = t ^= b1;
    mult_w[j + 248] = t ^= b0;
    mult_w[j + 232] = t ^= b[i + 4];
    mult_w[j + 233] = t ^= b0;
    mult_w[j + 235] = t ^= b1;
    mult_w[j + 234] = t ^= b0;
    mult_w[j + 238] = t ^= b2;
    mult_w[j + 239] = t ^= b0;
    mult_w[j + 237] = t ^= b1;
    mult_w[j + 236] = t ^= b0;
    mult_w[j + 228] = t ^= b[i + 3];
    mult_w[j + 229] = t ^= b0;
    mult_w[j + 231] = t ^= b1;
    mult_w[j + 230] = t ^= b0;
    mult_w[j + 226] = t ^= b2;
    mult_w[j + 227] = t ^= b0;
    mult_w[j + 225] = t ^= b1;
    mult_w[j + 224] = t ^= b0;
    mult_w[j + 160] = t ^= b[i + 6];
    mult_w[j + 161] = t ^= b0;
    mult_w[j + 163] = t ^= b1;
    mult_w[j + 162] = t ^= b0;
    mult_w[j + 166] = t ^= b2;
    mult_w[j + 167] = t ^= b0;
    mult_w[j + 165] = t ^= b1;
    mult_w[j + 164] = t ^= b0;
    mult_w[j + 172] = t ^= b[i + 3];
    mult_w[j + 173] = t ^= b0;
    mult_w[j + 175] = t ^= b1;
    mult_w[j + 174] = t ^= b0;
    mult_w[j + 170] = t ^= b2;
    mult_w[j + 171] = t ^= b0;
    mult_w[j + 169] = t ^= b1;
    mult_w[j + 168] = t ^= b0;
    mult_w[j + 184] = t ^= b[i + 4];
    mult_w[j + 185] = t ^= b0;
    mult_w[j + 187] = t ^= b1;
    mult_w[j + 186] = t ^= b0;
    mult_w[j + 190] = t ^= b2;
    mult_w[j + 191] = t ^= b0;
    mult_w[j + 189] = t ^= b1;
    mult_w[j + 188] = t ^= b0;
    mult_w[j + 180] = t ^= b[i + 3];
    mult_w[j + 181] = t ^= b0;
    mult_w[j + 183] = t ^= b1;
    mult_w[j + 182] = t ^= b0;
    mult_w[j + 178] = t ^= b2;
    mult_w[j + 179] = t ^= b0;
    mult_w[j + 177] = t ^= b1;
    mult_w[j + 176] = t ^= b0;
    mult_w[j + 144] = t ^= b[i + 5];
    mult_w[j + 145] = t ^= b0;
    mult_w[j + 147] = t ^= b1;
    mult_w[j + 146] = t ^= b0;
    mult_w[j + 150] = t ^= b2;
    mult_w[j + 151] = t ^= b0;
    mult_w[j + 149] = t ^= b1;
    mult_w[j + 148] = t ^= b0;
    mult_w[j + 156] = t ^= b[i + 3];
    mult_w[j + 157] = t ^= b0;
    mult_w[j + 159] = t ^= b1;
    mult_w[j + 158] = t ^= b0;
    mult_w[j + 154] = t ^= b2;
    mult_w[j + 155] = t ^= b0;
    mult_w[j + 153] = t ^= b1;
    mult_w[j + 152] = t ^= b0;
    mult_w[j + 136] = t ^= b[i + 4];
    mult_w[j + 137] = t ^= b0;
    mult_w[j + 139] = t ^= b1;
    mult_w[j + 138] = t ^= b0;
    mult_w[j + 142] = t ^= b2;
    mult_w[j + 143] = t ^= b0;
    mult_w[j + 141] = t ^= b1;
    mult_w[j + 140] = t ^= b0;
    mult_w[j + 132] = t ^= b[i + 3];
    mult_w[j + 133] = t ^= b0;
    mult_w[j + 135] = t ^= b1;
    mult_w[j + 134] = t ^= b0;
    mult_w[j + 130] = t ^= b2;
    mult_w[j + 131] = t ^= b0;
    mult_w[j + 129] = t ^= b1;
    mult_w[j + 128] = t ^ b0;
  }
  {
    s32 i;
    for (i = 0; i < n; i++) {
      u32 t = (u32)a[i], u = (u32)(a[i] >> 32);
      c[i] = mult_w[t & 255] ^
             mult_w[256 + ((t >> 8) & 255)] ^
             mult_w[512 + ((t >> 16) & 255)] ^
             mult_w[768 + (t >> 24)] ^
             mult_w[1024 + (u & 255)] ^
             mult_w[1280 + ((u >> 8) & 255)] ^
             mult_w[1536 + ((u >> 16) & 255)] ^
             mult_w[1792 + (u >> 24)];
    }
  }
}

void addmultnx64(u64 *c, u64 *a, u64 *b, s32 n) {
  int i, j;
  for (i = 0, j = 0; i < 64; i += 8, j += 256) {
    u64 b0 = b[i];
    u64 b1 = b[i + 1];
    u64 b2 = b[i + 2];
    u64 t;
    mult_w[j] = 0;
    mult_w[j + 1] = t = b0;
    mult_w[j + 3] = t ^= b1;
    mult_w[j + 2] = t ^= b0;
    mult_w[j + 6] = t ^= b2;
    mult_w[j + 7] = t ^= b0;
    mult_w[j + 5] = t ^= b1;
    mult_w[j + 4] = t ^= b0;
    mult_w[j + 12] = t ^= b[i + 3];
    mult_w[j + 13] = t ^= b0;
    mult_w[j + 15] = t ^= b1;
    mult_w[j + 14] = t ^= b0;
    mult_w[j + 10] = t ^= b2;
    mult_w[j + 11] = t ^= b0;
    mult_w[j + 9] = t ^= b1;
    mult_w[j + 8] = t ^= b0;
    mult_w[j + 24] = t ^= b[i + 4];
    mult_w[j + 25] = t ^= b0;
    mult_w[j + 27] = t ^= b1;
    mult_w[j + 26] = t ^= b0;
    mult_w[j + 30] = t ^= b2;
    mult_w[j + 31] = t ^= b0;
    mult_w[j + 29] = t ^= b1;
    mult_w[j + 28] = t ^= b0;
    mult_w[j + 20] = t ^= b[i + 3];
    mult_w[j + 21] = t ^= b0;
    mult_w[j + 23] = t ^= b1;
    mult_w[j + 22] = t ^= b0;
    mult_w[j + 18] = t ^= b2;
    mult_w[j + 19] = t ^= b0;
    mult_w[j + 17] = t ^= b1;
    mult_w[j + 16] = t ^= b0;
    mult_w[j + 48] = t ^= b[i + 5];
    mult_w[j + 49] = t ^= b0;
    mult_w[j + 51] = t ^= b1;
    mult_w[j + 50] = t ^= b0;
    mult_w[j + 54] = t ^= b2;
    mult_w[j + 55] = t ^= b0;
    mult_w[j + 53] = t ^= b1;
    mult_w[j + 52] = t ^= b0;
    mult_w[j + 60] = t ^= b[i + 3];
    mult_w[j + 61] = t ^= b0;
    mult_w[j + 63] = t ^= b1;
    mult_w[j + 62] = t ^= b0;
    mult_w[j + 58] = t ^= b2;
    mult_w[j + 59] = t ^= b0;
    mult_w[j + 57] = t ^= b1;
    mult_w[j + 56] = t ^= b0;
    mult_w[j + 40] = t ^= b[i + 4];
    mult_w[j + 41] = t ^= b0;
    mult_w[j + 43] = t ^= b1;
    mult_w[j + 42] = t ^= b0;
    mult_w[j + 46] = t ^= b2;
    mult_w[j + 47] = t ^= b0;
    mult_w[j + 45] = t ^= b1;
    mult_w[j + 44] = t ^= b0;
    mult_w[j + 36] = t ^= b[i + 3];
    mult_w[j + 37] = t ^= b0;
    mult_w[j + 39] = t ^= b1;
    mult_w[j + 38] = t ^= b0;
    mult_w[j + 34] = t ^= b2;
    mult_w[j + 35] = t ^= b0;
    mult_w[j + 33] = t ^= b1;
    mult_w[j + 32] = t ^= b0;
    mult_w[j + 96] = t ^= b[i + 6];
    mult_w[j + 97] = t ^= b0;
    mult_w[j + 99] = t ^= b1;
    mult_w[j + 98] = t ^= b0;
    mult_w[j + 102] = t ^= b2;
    mult_w[j + 103] = t ^= b0;
    mult_w[j + 101] = t ^= b1;
    mult_w[j + 100] = t ^= b0;
    mult_w[j + 108] = t ^= b[i + 3];
    mult_w[j + 109] = t ^= b0;
    mult_w[j + 111] = t ^= b1;
    mult_w[j + 110] = t ^= b0;
    mult_w[j + 106] = t ^= b2;
    mult_w[j + 107] = t ^= b0;
    mult_w[j + 105] = t ^= b1;
    mult_w[j + 104] = t ^= b0;
    mult_w[j + 120] = t ^= b[i + 4];
    mult_w[j + 121] = t ^= b0;
    mult_w[j + 123] = t ^= b1;
    mult_w[j + 122] = t ^= b0;
    mult_w[j + 126] = t ^= b2;
    mult_w[j + 127] = t ^= b0;
    mult_w[j + 125] = t ^= b1;
    mult_w[j + 124] = t ^= b0;
    mult_w[j + 116] = t ^= b[i + 3];
    mult_w[j + 117] = t ^= b0;
    mult_w[j + 119] = t ^= b1;
    mult_w[j + 118] = t ^= b0;
    mult_w[j + 114] = t ^= b2;
    mult_w[j + 115] = t ^= b0;
    mult_w[j + 113] = t ^= b1;
    mult_w[j + 112] = t ^= b0;
    mult_w[j + 80] = t ^= b[i + 5];
    mult_w[j + 81] = t ^= b0;
    mult_w[j + 83] = t ^= b1;
    mult_w[j + 82] = t ^= b0;
    mult_w[j + 86] = t ^= b2;
    mult_w[j + 87] = t ^= b0;
    mult_w[j + 85] = t ^= b1;
    mult_w[j + 84] = t ^= b0;
    mult_w[j + 92] = t ^= b[i + 3];
    mult_w[j + 93] = t ^= b0;
    mult_w[j + 95] = t ^= b1;
    mult_w[j + 94] = t ^= b0;
    mult_w[j + 90] = t ^= b2;
    mult_w[j + 91] = t ^= b0;
    mult_w[j + 89] = t ^= b1;
    mult_w[j + 88] = t ^= b0;
    mult_w[j + 72] = t ^= b[i + 4];
    mult_w[j + 73] = t ^= b0;
    mult_w[j + 75] = t ^= b1;
    mult_w[j + 74] = t ^= b0;
    mult_w[j + 78] = t ^= b2;
    mult_w[j + 79] = t ^= b0;
    mult_w[j + 77] = t ^= b1;
    mult_w[j + 76] = t ^= b0;
    mult_w[j + 68] = t ^= b[i + 3];
    mult_w[j + 69] = t ^= b0;
    mult_w[j + 71] = t ^= b1;
    mult_w[j + 70] = t ^= b0;
    mult_w[j + 66] = t ^= b2;
    mult_w[j + 67] = t ^= b0;
    mult_w[j + 65] = t ^= b1;
    mult_w[j + 64] = t ^= b0;
    mult_w[j + 192] = t ^= b[i + 7];
    mult_w[j + 193] = t ^= b0;
    mult_w[j + 195] = t ^= b1;
    mult_w[j + 194] = t ^= b0;
    mult_w[j + 198] = t ^= b2;
    mult_w[j + 199] = t ^= b0;
    mult_w[j + 197] = t ^= b1;
    mult_w[j + 196] = t ^= b0;
    mult_w[j + 204] = t ^= b[i + 3];
    mult_w[j + 205] = t ^= b0;
    mult_w[j + 207] = t ^= b1;
    mult_w[j + 206] = t ^= b0;
    mult_w[j + 202] = t ^= b2;
    mult_w[j + 203] = t ^= b0;
    mult_w[j + 201] = t ^= b1;
    mult_w[j + 200] = t ^= b0;
    mult_w[j + 216] = t ^= b[i + 4];
    mult_w[j + 217] = t ^= b0;
    mult_w[j + 219] = t ^= b1;
    mult_w[j + 218] = t ^= b0;
    mult_w[j + 222] = t ^= b2;
    mult_w[j + 223] = t ^= b0;
    mult_w[j + 221] = t ^= b1;
    mult_w[j + 220] = t ^= b0;
    mult_w[j + 212] = t ^= b[i + 3];
    mult_w[j + 213] = t ^= b0;
    mult_w[j + 215] = t ^= b1;
    mult_w[j + 214] = t ^= b0;
    mult_w[j + 210] = t ^= b2;
    mult_w[j + 211] = t ^= b0;
    mult_w[j + 209] = t ^= b1;
    mult_w[j + 208] = t ^= b0;
    mult_w[j + 240] = t ^= b[i + 5];
    mult_w[j + 241] = t ^= b0;
    mult_w[j + 243] = t ^= b1;
    mult_w[j + 242] = t ^= b0;
    mult_w[j + 246] = t ^= b2;
    mult_w[j + 247] = t ^= b0;
    mult_w[j + 245] = t ^= b1;
    mult_w[j + 244] = t ^= b0;
    mult_w[j + 252] = t ^= b[i + 3];
    mult_w[j + 253] = t ^= b0;
    mult_w[j + 255] = t ^= b1;
    mult_w[j + 254] = t ^= b0;
    mult_w[j + 250] = t ^= b2;
    mult_w[j + 251] = t ^= b0;
    mult_w[j + 249] = t ^= b1;
    mult_w[j + 248] = t ^= b0;
    mult_w[j + 232] = t ^= b[i + 4];
    mult_w[j + 233] = t ^= b0;
    mult_w[j + 235] = t ^= b1;
    mult_w[j + 234] = t ^= b0;
    mult_w[j + 238] = t ^= b2;
    mult_w[j + 239] = t ^= b0;
    mult_w[j + 237] = t ^= b1;
    mult_w[j + 236] = t ^= b0;
    mult_w[j + 228] = t ^= b[i + 3];
    mult_w[j + 229] = t ^= b0;
    mult_w[j + 231] = t ^= b1;
    mult_w[j + 230] = t ^= b0;
    mult_w[j + 226] = t ^= b2;
    mult_w[j + 227] = t ^= b0;
    mult_w[j + 225] = t ^= b1;
    mult_w[j + 224] = t ^= b0;
    mult_w[j + 160] = t ^= b[i + 6];
    mult_w[j + 161] = t ^= b0;
    mult_w[j + 163] = t ^= b1;
    mult_w[j + 162] = t ^= b0;
    mult_w[j + 166] = t ^= b2;
    mult_w[j + 167] = t ^= b0;
    mult_w[j + 165] = t ^= b1;
    mult_w[j + 164] = t ^= b0;
    mult_w[j + 172] = t ^= b[i + 3];
    mult_w[j + 173] = t ^= b0;
    mult_w[j + 175] = t ^= b1;
    mult_w[j + 174] = t ^= b0;
    mult_w[j + 170] = t ^= b2;
    mult_w[j + 171] = t ^= b0;
    mult_w[j + 169] = t ^= b1;
    mult_w[j + 168] = t ^= b0;
    mult_w[j + 184] = t ^= b[i + 4];
    mult_w[j + 185] = t ^= b0;
    mult_w[j + 187] = t ^= b1;
    mult_w[j + 186] = t ^= b0;
    mult_w[j + 190] = t ^= b2;
    mult_w[j + 191] = t ^= b0;
    mult_w[j + 189] = t ^= b1;
    mult_w[j + 188] = t ^= b0;
    mult_w[j + 180] = t ^= b[i + 3];
    mult_w[j + 181] = t ^= b0;
    mult_w[j + 183] = t ^= b1;
    mult_w[j + 182] = t ^= b0;
    mult_w[j + 178] = t ^= b2;
    mult_w[j + 179] = t ^= b0;
    mult_w[j + 177] = t ^= b1;
    mult_w[j + 176] = t ^= b0;
    mult_w[j + 144] = t ^= b[i + 5];
    mult_w[j + 145] = t ^= b0;
    mult_w[j + 147] = t ^= b1;
    mult_w[j + 146] = t ^= b0;
    mult_w[j + 150] = t ^= b2;
    mult_w[j + 151] = t ^= b0;
    mult_w[j + 149] = t ^= b1;
    mult_w[j + 148] = t ^= b0;
    mult_w[j + 156] = t ^= b[i + 3];
    mult_w[j + 157] = t ^= b0;
    mult_w[j + 159] = t ^= b1;
    mult_w[j + 158] = t ^= b0;
    mult_w[j + 154] = t ^= b2;
    mult_w[j + 155] = t ^= b0;
    mult_w[j + 153] = t ^= b1;
    mult_w[j + 152] = t ^= b0;
    mult_w[j + 136] = t ^= b[i + 4];
    mult_w[j + 137] = t ^= b0;
    mult_w[j + 139] = t ^= b1;
    mult_w[j + 138] = t ^= b0;
    mult_w[j + 142] = t ^= b2;
    mult_w[j + 143] = t ^= b0;
    mult_w[j + 141] = t ^= b1;
    mult_w[j + 140] = t ^= b0;
    mult_w[j + 132] = t ^= b[i + 3];
    mult_w[j + 133] = t ^= b0;
    mult_w[j + 135] = t ^= b1;
    mult_w[j + 134] = t ^= b0;
    mult_w[j + 130] = t ^= b2;
    mult_w[j + 131] = t ^= b0;
    mult_w[j + 129] = t ^= b1;
    mult_w[j + 128] = t ^ b0;
  }
  {
    s32 i;
    for (i = 0; i < n; i++) {
      u32 t = (u32)a[i], u = (u32)(a[i] >> 32);
      c[i] ^= mult_w[t & 255] ^
              mult_w[256 + ((t >> 8) & 255)] ^
              mult_w[512 + ((t >> 16) & 255)] ^
              mult_w[768 + (t >> 24)] ^
              mult_w[1024 + (u & 255)] ^
              mult_w[1280 + ((u >> 8) & 255)] ^
              mult_w[1536 + ((u >> 16) & 255)] ^
              mult_w[1792 + (u >> 24)];
    }
  }
}


/*********************************************/ 
int isZeroV(u64 *A, s32 size)
/*********************************************/ 
{ u64 i;

  for (i=0; i<size; i++)
    if (A[i])
      return 0;
  return 1;
}

/*********************************************/ 
void getW_S(u64 *Wi, int *Si, u64 *T, int *Si_1)
{ u64 M[128] ALIGNED16, mask, t0, t1;
  int  c[64], i, j, k, sSize, s;
  
  /* Number the columns, with those in Si_1 coming last. */  
  for (sSize=0, mask=0x00000000, i=0; i<64; i++) {
    s = Si_1[i];
    if ((s>=0) && (s<64)) {
      mask |= BIT64(s);
      sSize++;
    }
  }
  for (i=0, j=0, k=0; i<64; i++) {
    if (mask&BIT64(i)) 
      c[64 - sSize + j++] = i;
    else
      c[k++] = i;
  }
  /* S <-- empty set. */
  sSize=0;
  for (i=0; i<64; i++)
    Si[i] = -1;

  /* Set M <-- [ T | I_64] */
  for (i=0; i<64; i++) {
    M[2*i] = T[i];
    M[2*i+1] = BIT64(i);
  }
  for (j=0; j<64; j++) {
    for (k=j; (k<64)&& ( (M[2*c[j]]&BIT64(c[j]))==0); k++) {
      if (M[2*c[k]]&BIT64(c[j])) {
        /* Exchange rows c[j] and c[k] of M */
        t0 = M[2*c[j]];  t1 = M[2*c[j]+1];
        M[2*c[j]] = M[2*c[k]]; M[2*c[j]+1] = M[2*c[k]+1];
        M[2*c[k]] = t0; M[2*c[k]+1] = t1;
      }
    }
    if ((M[2*c[j]]&BIT64(c[j]))) {
      Si[sSize++] = c[j];
      /* Now, add row c[j] to other rows as necessary, to zero out */
      /* the rest of column c[j].                                  */
      for (k=0; k<64; k++) {
        if ((k != c[j]) && (M[2*k]&BIT64(c[j]))) {
          M[2*k] ^= M[2*c[j]];
          M[2*k+1] ^= M[2*c[j]+1];
        }
      }
    } else {
      for (k=j; (k<64) && ( (M[2*c[j]+1]&BIT64(c[j]))==0); k++) {
        if (M[2*c[k]+1]&BIT64(c[j])) {
          /* Exchange rows c[j] and c[k] of M */
          t0 = M[2*c[j]];  t1 = M[2*c[j]+1];
          M[2*c[j]] = M[2*c[k]]; M[2*c[j]+1] = M[2*c[k]+1];
          M[2*c[k]] = t0; M[2*c[k]+1] = t1;
        }
      }
      /* Now, add row c[j] to other rows as necessary, to zero out */
      /* the rest of column c[j]+64                                */
      for (k=0; k<64; k++) {
        if ((k != c[j]) && (M[2*k+1]&BIT64(c[j]))) {
          M[2*k] ^= M[2*c[j]];
          M[2*k+1] ^= M[2*c[j]+1];
        }
      }
      /* Zero out row c[j]. */
      M[2*c[j]] = 0; M[2*c[j]+1] = 0;
    } /* end 'else'  */
  } /* end 'for j' */
  for (j=0; j<64; j++)
    Wi[j] = M[2*j+1];
}

/*******************************************************/
int doColumnOps(u64 *A, u64 *B, s32 n)
/*******************************************************/
/* Do column operations (Gaussian elimination) on 'B', */
/* applying the same operations to 'A' as we go.       */
/* Input: 'A' and 'B' are nx64 matrices.               */
/* Output: 'A' and 'B' are changed, as above. i.e., 'B'*/
/*         is in column-reduced echelon form, and the  */
/*         same transformation was applied to 'A'.     */
/*******************************************************/
/* Remark: This is horribly inefficient. There is      */
/* surely a better way to do it. But the easiest is to */
/* just transpose the matrices before and after.       */
/* But it doesn't much matter either way - this guy is */
/* not called from inside any loops, and it's not a    */
/* bottleneck.                                         */
/*******************************************************/
{ u64 *AT, *BT;
  u64 i, j, k, rowSize, col, row, t;
  int cont;

  rowSize = n/64;
  if (n%64)
    rowSize += 1;
  
  if (!(AT = (u64 *)malloc(64*rowSize*sizeof(u64)))) {
    printf("doColumnOps(): Memory allocation error!\n");
    return -1;
  }
  if (!(BT = (u64 *)malloc(64*rowSize*sizeof(u64)))) {
    printf("doColumnOps(): Memory allocation error!\n");
    free(AT);
    return -1;
  }
  /* Compute AT = A^T and BT = B^T. */
  memset(AT, 0x00, 64*rowSize*sizeof(u64));
  memset(BT, 0x00, 64*rowSize*sizeof(u64));

  for (row=0; row<n; row++) {
    for (col=0; col<64; col++) {
      if (A[2*row]&BIT64(col&0x3F))
        AT[col*rowSize + row/64] ^= BIT64(row&0x3F);
      if (B[2*row]&BIT64(col&0x3F))
        BT[col*rowSize + row/64] ^= BIT64(row&0x3F);
    }
    for (col=64; col<64; col++) {
      if (A[2*row + 1]&BIT64(col&0x3F))
        AT[col*rowSize + row/64] ^= BIT64(row&0x3F);
      if (B[2*row + 1]&BIT64(col&0x3F))
        BT[col*rowSize + row/64] ^= BIT64(row&0x3F);
    }
  }

  /* Now, apply row ops to get 'BT' into RREF. */
  col=0;
  for (row=0; row<64; row++) {
    /* Is there a row at or after 'row' with a '1' in column 'col'? */
    /* Consider adding a short circuit, testing whole words at a time */
    /* since we are generally expecting small rank for our matrices.  */
    cont=1;
    do {
      k=row;
      while ((k<64)&&((BT[rowSize*k + col/64]&BIT64(col&0x3F))==0))
        k++;
      if (k < 64) {
        /* Swap rows 'row' and 'k'. */
        for (i=0; i<rowSize; i++) {
          t = BT[rowSize*row + i]; 
          BT[rowSize*row + i] = BT[rowSize*k + i];
          BT[rowSize*k + i] = t;
          t = AT[rowSize*row + i]; 
          AT[rowSize*row + i] = AT[rowSize*k + i];
          AT[rowSize*k + i] = t;
        }
        /* Use row 'row' to zero out the rest of this column. */
        for (i=0; i<64; i++) {
          if ((BT[i*rowSize + col/64]&BIT64(col&0x3F)) && (i != row)) {
            for (j=0; j<rowSize; j++) {
              BT[i*rowSize + j] ^= BT[row*rowSize + j];
              AT[i*rowSize + j] ^= AT[row*rowSize + j];
            }
          }
        } 
        cont=0;
      } else {
        /* Column 'col' is zero at and below row 'row'. */
        cont=1;
      }
      col++;
    } while ((cont) && (col < n));
    if (col >= n)
      break; /* CJM, 10/20/03. */
  }

  /* Compute A = AT^T and B = BT^T. */
  for (i=0; i<2*n; i++)
    A[i] = B[i] = 0;

  for (row=0; row<64; row++) {
    for (col=0; col<n; col++) {
      if (AT[row*rowSize + col/64]&BIT64(col&0x3F))
        A[2*col + row/64] ^= BIT64(row&0x3F);
      if (BT[row*rowSize + col/64]&BIT64(col&0x3F))
        B[2*col + row/64] ^= BIT64(row&0x3F);
    }
  }
  free(AT);
  free(BT);
  return 0;
}

/*******************************************************/
int doColumnOps64(u64 *A, u64 *B, s32 n)
/*******************************************************/
/* Do column operations (Gaussian elimination) on 'B', */
/* applying the same operations to 'A' as we go.       */
/* Input: 'A' and 'B' are nx64 matrices.               */
/* Output: 'A' and 'B' are changed, as above. i.e., 'B'*/
/*         is in column-reduced echelon form, and the  */
/*         same transformation was applied to 'A'.     */
/*******************************************************/
/* Remark: This is horribly inefficient. There is      */
/* surely a better way to do it. But the easiest is to */
/* just transpose the matrices before and after.       */
/* But it doesn't much matter either way - this guy is */
/* not called from inside any loops, and it's not a    */
/* bottleneck.                                         */
/*******************************************************/
{ u64 *AT, *BT;
  u64 i, j, k, rowSize, col, row, t;
  int cont;

  rowSize = n/64;
  if (n%64)
    rowSize += 1;
  
  if (!(AT = (u64 *)malloc(64*rowSize*sizeof(u64)))) {
    printf("doColumnOps(): Memory allocation error!\n");
    return -1;
  }
  if (!(BT = (u64 *)malloc(64*rowSize*sizeof(u64)))) {
    printf("doColumnOps(): Memory allocation error!\n");
    free(AT);
    return -1;
  }
  /* Compute AT = A^T and BT = B^T. */
  memset(AT, 0x00, 64*rowSize*sizeof(u64));
  memset(BT, 0x00, 64*rowSize*sizeof(u64));

  for (row=0; row<n; row++) {
    for (col=0; col<64; col++) {
      if (A[row]&BIT64(col&0x3F))
        AT[col*rowSize + row/64] ^= BIT64(row&0x3F);
      if (B[row]&BIT64(col&0x3F))
        BT[col*rowSize + row/64] ^= BIT64(row&0x3F);
    }
  }

  /* Now, apply row ops to get 'BT' into RREF. */
  col=0;
  for (row=0; row<64; row++) {
    /* Is there a row at or after 'row' with a '1' in column 'col'? */
    /* Consider adding a short circuit, testing whole words at a time */
    /* since we are generally expecting small rank for our matrices.  */
    cont=1;
    do {
      k=row;
      while ((k<64)&&((BT[rowSize*k + col/64]&BIT64(col&0x3F))==0))
        k++;
      if (k < 64) {
        /* Swap rows 'row' and 'k'. */
        for (i=0; i<rowSize; i++) {
          t = BT[rowSize*row + i]; 
          BT[rowSize*row + i] = BT[rowSize*k + i];
          BT[rowSize*k + i] = t;
          t = AT[rowSize*row + i]; 
          AT[rowSize*row + i] = AT[rowSize*k + i];
          AT[rowSize*k + i] = t;
        }
        /* Use row 'row' to zero out the rest of this column. */
        for (i=0; i<64; i++) {
          if ((BT[i*rowSize + col/64]&BIT64(col&0x3F)) && (i != row)) {
            for (j=0; j<rowSize; j++) {
              BT[i*rowSize + j] ^= BT[row*rowSize + j];
              AT[i*rowSize + j] ^= AT[row*rowSize + j];
            }
          }
        } 
        cont=0;
      } else {
        /* Column 'col' is zero at and below row 'row'. */
        cont=1;
      }
      col++;
    } while ((cont) && (col < n));
    if (col >= n)
      break; /* CJM, 10/20/03. */
  }
  /* Compute A = AT^T and B = BT^T. */
  for (i=0; i<n; i++)
    A[i] = B[i] = 0;

  for (row=0; row<64; row++) {
    for (col=0; col<n; col++) {
      if (AT[row*rowSize + col/64]&BIT64(col&0x3F))
        A[col] ^= BIT64(row);
      if (BT[row*rowSize + col/64]&BIT64(col&0x3F))
        B[col] ^= BIT64(row);
    }
  }
  free(AT);
  free(BT);
  return 0;
}

