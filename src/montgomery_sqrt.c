/**************************************************************/
/* montgomery_sqrt.c                                          */
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
#include <sys/types.h>
#include <sys/stat.h>

#include "ggnfs.h"

#define MAX_IPBSIZE  100

#define DEFAULT_LLLMAX_LOG 401.0*log(10.0)

#define MAX_CD_EXTRA_PRIMES 20

//#define _DEBUG

#ifdef GMP_BUG
#define patched_mpz_set_d(z, d) { double d_ = (d);  \
  if (-1.0 < d_ && d_ < 1.0) mpz_set_ui((z), 0); else mpz_set_d((z), d_); }
#else 
#define patched_mpz_set_d(z, d) mpz_set_d(z, d)
#endif

#ifdef _DEBUG
#define MAX_IDEAL_STR 4096
char idealSelStr[MAX_IDEAL_STR+1];
#endif

typedef struct {
  s32 p, r;
} afb_elt_t;

#define MAX_DIST_FACTS 32

long ABexponentSum=0;
typedef struct {
  mpz_t        c_d;     /* Leading coefficient of 'f'. */
  int          numFactors;
  __mpz_struct p[MAX_DIST_FACTS]; /* The distinct primes dividing c_d. */
               /* If p[j] is in the AFB, it's first occurrence is at aI[j].
                  Otherwise, aI[j] < 0. */
  s32         aI[MAX_DIST_FACTS]; 
               /* If p[j] is an exceptional prime, it's first occurrence
                  is at eI[j]. Otherwise, eI[j] < 0. */
  s32         eI[MAX_DIST_FACTS]; 
} cd_info;

typedef struct {
  nf_t       *N;  /* The number field. */
  nfs_fb_t   *FB; /* The factor base.  */
  afb_elt_t  *AFB;  /* ideals: from AFB and large leftover primes. */
  prime_id_t *sPrimes; /* Special primes. */
  cd_info     Cd;      /* Various junk needed to handle c_d. */
  int        *v_cd_sPrimes; /* Valuations of c_d at the sPrimes. */
  s32         aMax, aSize, spMax, spSize; /* allocated size / actual size. */
  s32        *aExp, aExpLast; /* Valuations of <gamma_i> at the AFB primes. 
                               Note: these can be negative since <gamma_i> is
                               a fractional ideal!  */
  s32        *spExp;          /* Valuations of <gamma_i> at special primes. */
  mpz_mat_t   Hnum, Hden;     /* HNF of stuff leftover from previous step. */
  double      logNormGNum, logNormGDen; /* Log of N(Numer(<gamma_i>)) and N(Denom(<gamma_i>)) */
  double      log_eps[MAXPOLYDEGREE]; /* log of the embeddings sizes.      */
   __mpz_struct q[MAX_IPBSIZE];
  int         ipbSize;
  mpz_poly    crtRes[MAX_IPBSIZE]; /* Residues mod q, with \hat{alpha} representation. */
  mpz_mat_t   Beta; /* The beta_j, as defined in Nguyen. */
  mpz_t       kappa; /* The index of the maximal order over the 'A' order, in Nguyen. */
  mpz_t       ratSqrt; /* Obviously, the rational square root. */
  double      LLL_max_log;
  __mpz_struct  omegaEvalM[MAXPOLYDEGREE]; /* The omega_i evaluated at 'alpha=m' modulo 'n'. */
} msqrt_t;

typedef struct {
  s32 p;
  s32 e;
} rat_p_t;


/* Prototypes for locally used stuff. */
s32    locateP(s32 p, msqrt_t *M);
void   updateEps_ab(msqrt_t *M, s64 a, s64 b, int exponent);
void   updateEps(msqrt_t *M, mpz_poly delta, int exponent);
void   updateCRT(msqrt_t *M, mpz_poly delta, mpz_t denom, int exponent);
void   updateFactorization_ab(msqrt_t *M, relation_t *R, int sl);
double logNorm(mpz_mat_t *HNF, msqrt_t *M);
int    choose_ab_exponent(msqrt_t *M, relation_t *R);
void   idealDiv(mpz_mat_t *Res, mpz_mat_t *I, mpz_mat_t *J, nf_t *N);

/*********************************************************************/
int afb_elt_cmp(const void *A, const void *B)
/*********************************************************************/
{ afb_elt_t *a=(afb_elt_t *)A, *b=(afb_elt_t *)B;

  if (a->p < b->p) return -1;
  if (a->p > b->p) return 1;
  if (a->r < b->r) return -1;
  if (a->r > b->r) return 1;
  return 0;
}

/*********************************************************************/
int cmp_rat_p_t(const void *A, const void *B)
/*********************************************************************/
{ rat_p_t *a=(rat_p_t *)A, *b=(rat_p_t *)B;

  if (a->p < b->p) return -1;
  if (a->p > b->p) return 1;
  return 0;
}

/*********************************************************************/
void fix_for_cd(msqrt_t *M)
/*********************************************************************/
/* As I recall, this function is adjusting the factor bases so that  */
/* all factors of the leading coefficient c_d are in there somewhere.*/
/*********************************************************************/
{ mpz_fact_t F;
  int        retVal, i, j, addToFB, numZeros, k, reSort=0;
  s32       p, Zeros[MAXPOLYDEGREE];
  poly_t     f_;


  mpz_fact_init(&F);
  retVal = mpz_fact_factorEasy(&F, M->Cd.c_d, 1);
  if (retVal) {
    msgLog("", "Error: Unable to completely factor c_d!\n");
    fprintf(stderr, "fix_for_cd(): Unable to completely factor c_d!\n");
    return;
  }
  printf("fix_for_cd(): c_d factored.\n");
  M->Cd.numFactors = F.size;
  if (F.size > MAX_DIST_FACTS) {
    fprintf(stderr, "fix_for_cd() Error: MAX_DIST_FACTS < %d! Increase and recompile.\n", F.size);
    exit(-1);
  }
  for (i=0; i<F.size; i++) {
    printf("Checking on factor "); mpz_out_str(stdout, 10, &F.factors[i]); 
    printf("..."); fflush(stdout);
    mpz_set(&M->Cd.p[i], &F.factors[i]);
    M->Cd.aI[i] = -1;
    M->Cd.eI[i] = -1;
    addToFB = 1;
    /* Is this factor of c_d in the AFB? */
    if (mpz_fits_slong_p(&M->Cd.p[i])) {
      p = mpz_get_si(&M->Cd.p[i]);
      j = locateP(p, M);
      /* Really, I think we should check that all primes over 'p' are in here. */
      if (M->AFB[j].p == p) {
        M->Cd.aI[i] = j;
        addToFB=0;
      }
    }
    /* Is it special? */
    for (j=0; j<M->spSize; j++) {
      if (mpz_cmp(&M->Cd.p[i], M->sPrimes[j].p)==0) {
        M->Cd.eI[i] = j;
        M->Cd.aI[i] = 0; /* It could be in the AFB also, in which case we should ignore it from there. */
        addToFB = 0;
        break;
      }
    }
    if (addToFB) {
      /* It should be added to the FB. Hopefully, it fits in a s32, */
      /* because we cannot otherwise handle it yet.                  */
      if (mpz_fits_slong_p(&M->Cd.p[i])) {
        p = mpz_get_si(&M->Cd.p[i]);
        j = locateP(p, M);
        /* Now, 'j' is the index for where it should go.        */    
        /* So, find all the zeros of 'f' mod p and add them in. */    
        mpz_poly_modp(f_, M->FB->f, p);
        numZeros = poly_getZeros(Zeros, f_, p);
        Zeros[numZeros++] = p; /* We know this prime divides c_d ;) */
        /* Shift the existing primes down so we can add these in. */
//        memmove(&M->AFB[j+numZeros], &M->AFB[j], (M->aSize-j)*sizeof(afb_elt_t));
//        for (k=0; k<numZeros; k++) {
//          M->AFB[j+k].p = p; M->AFB[j+k].r = Zeros[k];
//        }
        for (k=0; k<numZeros; k++) {
          M->AFB[M->aSize+k].p = p; M->AFB[M->aSize+k].r = Zeros[k];
        }
        M->aSize += numZeros;
        reSort = 1;
      } else {
        /* We're S.O.L. */
      }
      
    }
    if (M->Cd.aI[i] > 0)
      printf("in AFB at index %" PRId32 "...", M->Cd.aI[i]);
    if (M->Cd.eI[i] > 0)
      printf("is exceptional with index %" PRId32 "...", M->Cd.eI[i]);
    printf("\n");

  }
  if (reSort)  {
    qsort(M->AFB, M->aSize, sizeof(afb_elt_t), afb_elt_cmp);
    fix_for_cd(M); /* This is inefficient, but the indicies may have changed - */
                   /* We need to recheck them. */
  }

  mpz_fact_clear(&F);
}


/****************************************************************/
void initOmegaEvalM(msqrt_t *M)
/****************************************************************/
/* If the integral basis is omega_0,...,omega_{d-1}, this       */
/* function is computing the image of the omega_j under the     */
/* homomorphism \psi: O --> Z/NZ given by alpha |--> m.         */
/****************************************************************/
{ int    i, j, d=M->N->degree;
  mpz_t  cdm, cdmPow, tmp1;

  for (i=0; i<d; i++)
    mpz_init_set_ui(&M->omegaEvalM[i], 0);
  mpz_init(cdm); mpz_init_set_ui(cdmPow, 1);
  mpz_init(tmp1);

  mpz_invert(cdm, M->FB->y1, M->FB->n);
  mpz_mul(cdm, cdm, M->FB->y0);
  mpz_mul(cdm, cdm, &M->N->f->coef[d]);
  mpz_neg(cdm, cdm);
  mpz_mod(cdm, cdm, M->FB->n);
  for (i=0; i<d; i++) {
    for (j=0; j<d; j++) {
      mpz_mul(tmp1, cdmPow, &M->N->W->entry[i][j]);
      mpz_add(&M->omegaEvalM[j], &M->omegaEvalM[j], tmp1);
      mpz_mod(&M->omegaEvalM[j], &M->omegaEvalM[j], M->FB->n);
    }
    mpz_mul(cdmPow, cdmPow, cdm);
    mpz_mod(cdmPow, cdmPow, M->FB->n);
  }
  mpz_invert(tmp1, M->N->W_d, M->FB->n);
  for (i=0; i<d; i++) {
    mpz_mul(&M->omegaEvalM[i], &M->omegaEvalM[i], tmp1);
    mpz_mod(&M->omegaEvalM[i], &M->omegaEvalM[i], M->FB->n);
  }
  mpz_clear(cdm); mpz_clear(cdmPow); mpz_clear(tmp1);
}

/****************************************************************/
int alphatoib(mpz_poly res, mpz_poly op, msqrt_t *M)
/****************************************************************/
/* Input: A polynomial op in alpha, where alpha is a root of the*/
/*        non-monic polynomial, f.                              */
/* Output: res expressed as a polynomial in the omega_i (the    */
/*        integral basis). It is assumed that op is an          */
/*        algebraic integer so that this is possible.           */
/****************************************************************/
{ static mpz_poly tpol;
  static mpz_t    cd_pow, tmp;
  static int      initialized=0;
  int    i, d=M->N->degree, retVal;

  if (!(initialized)) {
    mpz_poly_init(tpol);
    mpz_init(cd_pow);
    mpz_init(tmp);
    initialized=1;
  }
  
  mpz_set_ui(cd_pow, 1);
  for (i=op->degree; i>=0; i--) {
    mpz_mul(&tpol->coef[i], &op->coef[i], cd_pow);
    mpz_mul(cd_pow, cd_pow, &M->FB->f->coef[d]);
  }
  tpol->degree = op->degree;
  retVal = stdtoib(res, tpol, M->N);
  if (!retVal) {
    /* We should be able to divide back out by c_d^{op->degree}. */
    mpz_pow_ui(cd_pow, &M->FB->f->coef[d], op->degree);
    for (i=0; i<=res->degree; i++) {
      mpz_tdiv_qr(&res->coef[i], tmp, &res->coef[i], cd_pow);
      if (mpz_sgn(tmp)) {
        printf("alphatoib() Error while correcting back for c_d!\n");
        exit(-1);
      }
    }
  }
  return retVal;
}

/****************************************************************/
int computeBeta(msqrt_t *M)
/****************************************************************/
/* Set the Beta field of M. The columns of the matrix Beta are  */
/* those corresponding to the Beta_i polynomials in Nguyen, but */
/* with respect to the integral basis instead of alpha.         */
/****************************************************************/
{ mpz_poly beta_a, beta_ib;
  int      i, j, d=M->N->degree;

  mpz_poly_init(beta_a); mpz_poly_init(beta_ib);
  M->Beta.rows = d;
  M->Beta.cols = d;
  for (i=0; i<d-1; i++) {
    for (j=i+1; j<=d; j++) 
      mpz_set(&beta_a->coef[j-i-1], &M->FB->f->coef[j]);
    beta_a->degree = d-i-1; 
    mpz_poly_print(stdout, "", beta_a);
    if (alphatoib(beta_ib, beta_a, M)) {
      fprintf(stderr, "computeBeta(): Severe error while computing beta_%d!\n", i);
      exit(-1);
    }
    setCol(&M->Beta, i, beta_ib);
  }  
  /* If we wanted Beta to be a representation of the corresponding
     order, all we need to do is to throw in the (1,0,...,0)^T column.
  */
  mpz_poly_clear(beta_a); mpz_poly_clear(beta_ib);

  return 0;
}
  
/*********************************************************************/
s64 findLP(s32 p, s32 r, msqrt_t *M)
/*********************************************************************/
{ afb_elt_t thisP, *res;
  s64 index;

  thisP.p = p; thisP.r = r;
  res = (afb_elt_t *)bsearch(&thisP, M->AFB, M->aSize, sizeof(afb_elt_t), afb_elt_cmp);
  if (res==NULL)
    return -1;
  index = (res - M->AFB);
  return index;
}

/*********************************************************************/
s32 locateP(s32 p, msqrt_t *M)
/*********************************************************************/
/* Find the first occurrence of a (p',r') pair in AFB with p'=p.     */
/* return value is the index where it was found (or should be).      */
/*********************************************************************/
/* This is terrible. We could probably use locateAFB instead, but    */
/* it's not a time-critical function (only used a handful of times)  */
/* and so it's not worth the overhead of initializing locateAFB or   */
/* the memory used by locateAFB to build a hash table.               */
/*********************************************************************/
{ s32 index, offset, diff, s=M->aSize;

  for (index=0; index<s; index++) 
    if (M->AFB[index].p>=p)
      return index;

  index = s/2;
  offset = index+1;

  diff = M->AFB[index].p - p;
  while (diff && (offset > 0)) {
    if (diff < 0) {
      index -= offset;
    } else if (diff > 0)
      index += offset;
    if (index < 0) index = 0;
    if (index >=s) index=s-1;
    if (offset > 1)
      offset = (1+offset)/2;
    else offset=0;
    diff = M->AFB[index].p - p;
  }
  while ((index >0) && (M->AFB[index-1].p == p))
    index--;
  return index;
}

#define AFB_TMP_SIZE 20
/*********************************************************************/
int setupPrimes(msqrt_t *M, multi_file_t *lpF)
/*********************************************************************/
/* Allocate for and initialize AFB so that it contains all the       */
/* non-special primes we will need.                                  */
/* Also, identify the special primes and store copies in sPrimes.    */
/*********************************************************************/
{ off_t      t1Size, maxSize;
  s32        p1, r1, h1, ct1, i1, i, k;
  afb_elt_t  *T1;
  mpz_t      tmp;
  struct     stat fileInfo; 
  char       fName[64];
  FILE       *fp;

  mpz_init(tmp);
  /*****************************************************************/
  /* Handle the exceptional primes first. Since we already have    */
  /* decompositions of all primes dividing the index [K:Z[alpha]], */
  /* this is fairly easy. Just find them and copy 'em over.        */
  /*****************************************************************/
  M->spMax = M->N->numSPrimes;
  if (!(M->sPrimes = (prime_id_t *)malloc(M->spMax*sizeof(prime_id_t)))) {
    fprintf(stderr, "setupPrimes(): Mem. allocation error for M->sPrimes!\n");
    return -1;
  }
  k = 0;
  for (i=0; i<M->N->numSPrimes; i++) {  
    initIdeal(&M->sPrimes[k]);
    mpz_set(M->sPrimes[k].p, M->N->sPrimes[i].p);
    mpz_poly_cp(M->sPrimes[k].alpha, M->N->sPrimes[i].alpha);
    mpz_poly_cp(M->sPrimes[k].beta, M->N->sPrimes[i].beta);
    M->sPrimes[k].e = M->N->sPrimes[i].e;
    M->sPrimes[k].f = M->N->sPrimes[i].f;
    k++;
  }
  M->spSize = k;
  printf("There are %" PRId32 " special prime ideals:\n", k);
  for (i=0; i<k; i++) {
    printf("Sp. %" PRId32 ": ", i);
    mpz_out_str(stdout, 10, M->sPrimes[i].p); printf(", ");
    mpz_poly_print(stdout, "",M->sPrimes[i].alpha); 
    printf("        Ramification index: %d\n", M->sPrimes[i].e);
    printf("        Norm: p^%d\n", M->sPrimes[i].f);
  }

  maxSize = 0;
  for (i=0; i<=lpF->numFiles; i++) {
    if (i < lpF->numFiles)
      sprintf(fName, "%s.%" PRId32, lpF->prefix, i);
    else
      sprintf(fName, "%s.L", lpF->prefix);
    if (stat(fName, &fileInfo)) {
      fprintf(stderr, "setupPrimes() Error: could not stat file %s!\n", fName);
      exit(-1);
    }
    if (i < lpF->numFiles)
      maxSize += fileInfo.st_size/(5*sizeof(s32));
    else
      maxSize += fileInfo.st_size/(4*sizeof(s32));
  }
  maxSize += 1000; /* For safety. */
  if (!(T1 = (afb_elt_t *)malloc(maxSize*sizeof(afb_elt_t)))) {
    fprintf(stderr, "setupPrimes() Error allocating %" PRIu32 " bytes for T1!\n", 
            (u32)(maxSize*sizeof(afb_elt_t)) );
    exit(-1);
  }

  /* Now, read in the large algebraic primes. */
  t1Size = 0;
  for (i=0; i<=lpF->numFiles; i++) {
    if (i < lpF->numFiles)
      sprintf(fName, "%s.%" PRId32, lpF->prefix, i);
    else
      sprintf(fName, "%s.L", lpF->prefix);
    if (!(fp = fopen(fName, "rb"))) {
      fprintf(stderr, "setupPrimes() Error: could not open file %s for read!\n", fName);
      exit(-1);
    }
    while (!(feof(fp)) && (t1Size < maxSize)) {
      r1 = -1;
      if (i < lpF->numFiles) {
        readRaw32(&h1, fp);
        readRaw32(&p1, fp);
        readRaw32(&r1, fp);
        readRaw32(&i1, fp);
        readRaw32(&ct1, fp);
      } else {
        readRaw32(&p1, fp);
        readRaw32(&r1, fp);
        readRaw32(&i1, fp);
        readRaw32(&ct1, fp);
      }
      if (r1 >= 0) {
        T1[t1Size].p = p1; T1[t1Size].r = r1;
        t1Size++;
      }
    }
    fclose(fp);
  }
  if (t1Size >= maxSize) {
    fprintf(stderr, "setupPrimes() severe error! maxSize=%" PRIu64 " exceeded!\n", (u64)maxSize);
    exit(-1);
  
  }
  /* And sort them. */
  qsort(T1, t1Size, sizeof(afb_elt_t), afb_elt_cmp);

  printf("Found %" PRIu64 " large primes total.\n", (u64)t1Size);
  /* Finally, the large primes we needed are in T1 and */
  /* there are t1Size of them.                         */
  M->aSize = t1Size + M->FB->afb_size;
  M->aMax = M->aSize + MAX_CD_EXTRA_PRIMES;
  /* Allocate, and copy everything over into M->AFB.   */
  if (!(M->AFB = (afb_elt_t *)malloc(M->aMax*sizeof(afb_elt_t)))) {
    fprintf(stderr, "setupAFB() Mem. allocation error for M->AFB!\n");
    if (T1 != NULL) free(T1); 
    return -1;
  }
  k = M->FB->afb_size;
  for (i=0; i<k; i++) {
    M->AFB[i].p = M->FB->afb[2*i];
    M->AFB[i].r = M->FB->afb[2*i+1];
  }
  for (i=0; i<t1Size; i++) {
    M->AFB[k+i].p = T1[i].p;
    M->AFB[k+i].r = T1[i].r;
  }
  if (T1 != NULL)
    free(T1);

  /* This isn't quite complete - we may still need to add
     some more things to the factor base. In particular,
     if c_d has any more prime factors which are not
     exceptional primes and are not in the AFB, we need
     to be able to reference them somehow! The following
     function doesn't do that right now - it only factors
     c_d and looks to see if the factors are in the AFB
     or the exceptional prime list. Otherwise, it does
     nothing with the remaining factors.
  */
  fix_for_cd(M);
  mpz_clear(tmp);
  return 0;
}

/*********************************************************************/
int getIPB_mpz(msqrt_t *M, double logTotal)
/*********************************************************************/
/* Get some large inert primes.                                      */
/* 'q[i]' should not be initialized. We will initialize              */
/* them here since it isn't really known in advance how many are     */
/* needed.                                                           */
/* Return value: number of primes obtained on success.               */
/*               negative, on error.                                 */
/*********************************************************************/
{ mpz_poly f;
  mpz_t    p, tmp;
  int      baseSize, i, n=M->N->degree;
  double   sumLogs;

  mpz_poly_init(f);
  mpz_init(tmp);
  mpz_init_set_str(p, "6000000000", 10); /* Arbitrary, but > 2^32. */
  sumLogs = 0;
  baseSize = 0;
  M->ipbSize = 0;
  while ((sumLogs < logTotal) && (baseSize < MAX_IPBSIZE)) {
    mpz_nextprime(p, p);
    mpz_mod(tmp, M->N->Tdisc, p);
    if (mpz_sgn(tmp)) {
      for (i=0; i<=n; i++)
        mpz_mod(&f->coef[i], &M->N->T->coef[i], p);
      f->degree = n;
      mpz_poly_fixDeg(f);
      if ((f->degree==n) && (mpz_poly_irreduciblelike_modp(f, p))) {
        mpz_init_set(&M->q[baseSize++], p);
        sumLogs += _mpz_log(p);
      }
    }
    mpz_add_ui(p, p, 1);
  }
  mpz_clear(tmp); mpz_clear(p);
  mpz_poly_clear(f);
  if ((sumLogs < logTotal) && (baseSize >= MAX_IPBSIZE)) {
    printf("getIPB_mpz() : MAX_IPBSIZE is too small - increase and recompile!\n");
    return -1;
  }
  M->ipbSize = baseSize;
  return baseSize;
}

/*
STEN: ((NR_MULTIPLIER * x )/ 10) is actually calculated. This approach
	  allows us to avoid 'double' to 'INT32' conversion warning message
	  (NR_MULTIPLIER was defined as 1.3 previously)
*/
#define NR_MULTIPLIER 13  
						   
/*******************************************************/
int ratSqrt(relation_t *R, int e, s32 depSize, msqrt_t *M)
/*******************************************************/
/* R is an (a,b) pair in the product, or NULL.         */
/* 'e' is the exponent of (a-bm) (e = +/-1).           */
/* `depSize' is the total number of dependencies. This */
/*           is needed to estimate the number of large */
/*           primes we will have to deal with.         */
/* M is self-explanatory since there is only one       */
/*           msqrt_t structure hanging around.         */
/* If R is NULL, we will compute the final square root */
/*           and un-initialize (freeing up memory).    */
/*******************************************************/
{ static int         initialized=0;
  static rat_p_t    *ratHashList, *ratLeftoverList;
  static s32        ratLeftoverSize, *RFB_exps, nr;
  s32               i, j, p, h;
  mpz_t              tmp;
  int                res=0;

  if (!initialized) {
    /**************************************************************************/
    /* Some temporary structures to help compute the rational square root.    */
    /* This works as follows:                                                 */
    /* (1) Find nr = an upper bound on the number of large rational primes.   */
    /* (2) Make two tables : ratHashList, ratLeftoverList of sizes 2*nr,      */
    /*     0.5*nr repsectively (this is actually overkill on the second one). */
    /* (3) Initialize all entries to have a prime field of 0, exponent 0.     */
    /* (4) For each large rational prime p and exponent e (+/- 1) we find in  */
    /*     the given (a,b) pair, do as follows:                               */
    /*     (4.1) Compute h <-- NFS_HASH(p, 2*nr).                             */
    /*     (4.2) If ratHashList[h].p == p, do ratHashList[h].e += e.          */
    /*           Otherwise, add (p,e) to ratLeftoverList.                     */
    /* (5) Sort ratLeftoverList, and remove duplicates by adding exponents.   */
    /* (6) Compute the rational square root modulo N by simply running through*/
    /*     the list and using 1/2 the exponent for each prime (they should    */
    /*     all be even if it's a good dependence).                            */
    /* (7) Free the tables - we won't need them no mo'.                       */
    /**************************************************************************/
    nr = (NR_MULTIPLIER*depSize*M->N->FB->maxLP) / 10;
    if (!(ratHashList = (rat_p_t *)malloc(2*nr*sizeof(rat_p_t)))) {
      fprintf(stderr, "ratSqrt() Error allocating %" PRIu32 " bytes for ratHashList!\n",
              (u32)(2*nr*sizeof(rat_p_t)) );
      exit(-1);
    }
    if (!(ratLeftoverList = (rat_p_t *)malloc((nr/2)*sizeof(rat_p_t)))) {
      fprintf(stderr, "ratSqrt() Error allocating %" PRIu32 " bytes for ratLeftoverList!\n",
              (u32)((nr/2)*sizeof(rat_p_t)) );
      free(ratHashList); exit(-1);
    }
    if (!(RFB_exps = (s32 *)malloc(M->N->FB->rfb_size*sizeof(s32)))) {
      fprintf(stderr, "ratSqrt() Error allocating %" PRIu32 " bytes for RFB_exps!\n", 
              (u32)(M->N->FB->rfb_size*sizeof(s32)) );
      free(ratHashList); free(ratLeftoverList); exit(-1);
    }
    for (i=0; i<2*nr; i++) {
      ratHashList[i].p = ratHashList[i].e = 0;
    }
    for (i=0, ratLeftoverSize=0; i<nr/2; i++) {
      ratLeftoverList[i].p = ratLeftoverList[i].e = 0;
    }
    for (i=0; i<M->N->FB->rfb_size; i++)
      RFB_exps[i] = 0;
  
    initialized=1;
  }

  if (R != NULL) {
    /* First, the easy part: count the exponents of */
    /* primes from RFB.                             */
    for (j=0; j<R->rFSize; j++)
      RFB_exps[R->rFactors[j]] += R->rExps[j]*e;
    /* Now, the large primes: */
    j=0;
    while ((j<MAX_LARGE_RAT_PRIMES) && (j < M->FB->maxLP)) {
      p = R->p[j];
      if (p>1) {
        h = NFS_HASH(p, 0, 2*nr);
        if (ratHashList[h].p==0) {
          ratHashList[h].p = p; ratHashList[h].e = e;
        } else if (ratHashList[h].p==p) {
          ratHashList[h].e += e;
        } else { /* There is something else in our table entry, so add to leftovers. */
          ratLeftoverList[ratLeftoverSize].p = p;
          ratLeftoverList[ratLeftoverSize].e = e;
          if (++ratLeftoverSize >= (nr/2)) {
            fprintf(stderr, "ratSqrt() Error: ratLeftoverList exceeded 2*nr!\n");
            fprintf(stderr, "Increase NR_MULTIPLIER, recompile, and try again.\n");
            msgLog("", "ratSqrt() Error: ratLeftoverList exceeded 2*nr!\n");
            msgLog("", "Increase NR_MULTIPLIER, recompile, and try again.\n");
            free(ratHashList); free(ratLeftoverList); initialized=0;
            return -1;
          }
        }
      }
      j++;
    }
  } else {
    mpz_init(tmp);
    /* Compute the rational square root: */
    /* First, we must sort and combine the leftover large primes: */
    qsort(ratLeftoverList, ratLeftoverSize, sizeof(rat_p_t), cmp_rat_p_t);
    for (i=1, j=0; i<ratLeftoverSize; i++) {
      if ((ratLeftoverList[i].p == ratLeftoverList[j].p)&&(i!=j)) {
        ratLeftoverList[j].e += ratLeftoverList[i].e;
        ratLeftoverList[i].e = 0;
      } else {
        /* Find the next leftover prime with a positive exponent. */
        do {
          j++;
        } while ((j < ratLeftoverSize)&&(ratLeftoverList[j].e == 0));
        if (i<j) i=j;
      }
    }
    /* Enough already - let's compute the darned thing! */
    mpz_init_set_ui(M->ratSqrt, 1);
    for (i=0; i<M->FB->rfb_size; i++) {
      if (RFB_exps[i]%2) {
        fprintf(stderr, "Error: RFB[%" PRId32 "] has odd exponent %" PRId32 "!\n",
                 i, RFB_exps[i]);
        exit(-1);
      }
      e = RFB_exps[i]/2;
  
      if (e) {
        mpz_set_ui(tmp, M->FB->rfb[2*i]);
        if (e > 1) 
          mpz_powm_ui(tmp, tmp, e, M->FB->n);
        else if (e < 0) {
          if (mpz_invert(tmp, tmp, M->FB->n))
            mpz_powm_ui(tmp, tmp, -e, M->FB->n);
          else {
            printf("Error: Inverse of %" PRId32 " does not exist mod n!", M->FB->rfb[2*i]); 
            printf("If this is an intentionally placed factor, re-run with -knowndiv.\n");
            exit(-1);
          }
        }
        mpz_mul(M->ratSqrt, M->ratSqrt, tmp);
        mpz_mod(M->ratSqrt, M->ratSqrt, M->FB->n);
      }
    }
    for (i=0; i<2*nr; i++) {
      if (ratHashList[i].e) {
        mpz_set_ui(tmp, ratHashList[i].p);
        e = ratHashList[i].e/2;
        if (ratHashList[i].e%2) {
          fprintf(stderr, "Error: Rational prime %" PRId32 " has odd exponent %" PRId32 "!\n",
                  ratHashList[i].p, ratHashList[i].e);
          msgLog(NULL, "Error: Rational prime %" PRId32 " has odd exponent %" PRId32 "!\n",
                 ratHashList[i].p, ratHashList[i].e);
          res=-1; goto RSQRT_CLEANUP;
        }
        if (e > 1) 
          mpz_powm_ui(tmp, tmp, e, M->FB->n);
        else if (e < 0) {
          mpz_invert(tmp, tmp, M->FB->n);
          mpz_powm_ui(tmp, tmp, -e, M->FB->n);
        }
        mpz_mul(M->ratSqrt, M->ratSqrt, tmp);
        mpz_mod(M->ratSqrt, M->ratSqrt, M->FB->n);
      }
    }
    for (i=0; i<ratLeftoverSize; i++) {
      if (ratLeftoverList[i].e) {
        mpz_set_ui(tmp, ratLeftoverList[i].p);
        e = ratLeftoverList[i].e/2;
        if (ratLeftoverList[i].e%2) {
          fprintf(stderr, "Error: Rational prime %" PRId32 " (leftover) has odd exponent %" PRId32 "!\n",
                  ratLeftoverList[i].p, ratLeftoverList[i].e);
          msgLog(NULL, "Error: Rational prime %" PRId32 " (leftover) has odd exponent %" PRId32 "!\n",
                  ratLeftoverList[i].p, ratLeftoverList[i].e);
          res=-1; goto RSQRT_CLEANUP;
        }
        if (e > 1)
          mpz_powm_ui(tmp, tmp, e, M->FB->n);
        else if (e < 0) {
          mpz_invert(tmp, tmp, M->FB->n);
          mpz_powm_ui(tmp, tmp, -e, M->FB->n);
        }
        mpz_mul(M->ratSqrt, M->ratSqrt, tmp);
        mpz_mod(M->ratSqrt, M->ratSqrt, M->FB->n);
      }
    }
RSQRT_CLEANUP:
    free(ratHashList); free(ratLeftoverList); free(RFB_exps);
    mpz_clear(tmp);
    initialized=0;
  }
  return res;
}

/*********************************************************************/
int initMsqrt(msqrt_t *M,  s32 *relsInDep, multi_file_t *prelF, multi_file_t *lpF)
/*********************************************************************/
/* M->N and M->FB must already be set.                               */
/*********************************************************************/
{ s32        i, j, k, depSize, R0, R1, Rindex;
  s64        a, b, rel;
  int         d=M->N->degree, e, fileNum;
  double      xr, xi, zr, zi, zpr, zpi, tr, ti, c;
  mpz_t       cd, tmp, Zsquare, tmp2, tmp3, bmultiplier;
  mpz_poly    tpol1;
  mpz_mat_t   H;
  char        fName[64], str[256];
  relation_t  R;
  FILE       *fp;
  s32        numPairs;
#ifdef _LOUD_DEBUG
  FILE *ofp;
#endif

ABexponentSum=0;
  mpz_init_set(cd, &M->FB->f->coef[d]);
  mpz_init_set(M->Cd.c_d, cd);
  mpz_init(Zsquare); mpz_init(tmp2); mpz_init(tmp3);
  for (i=0; i<MAX_DIST_FACTS; i++)
    mpz_init(&M->Cd.p[i]);
  mpz_init(tmp); mpz_init(bmultiplier);
  mpz_poly_init(tpol1);
  mpz_mat_init2(&H, d, d);
  mpz_mat_init2(&M->Beta, d, d);

  mpz_set(bmultiplier, M->N->W_d);
  mpz_mod(tmp, bmultiplier, &M->N->W->entry[1][1]);
  if (mpz_sgn(tmp)) {
    msgLog("", "Unexpected error computing bmultiplier! Sqrt computation will probably fail!");
    mpz_set_ui(bmultiplier, 1);
  } else
    mpz_div(bmultiplier, bmultiplier, &M->N->W->entry[1][1]);
  msgLog("", "bmultiplier=%s", mpz_get_str(str, 10, bmultiplier));

  /* Compute the index of the maximal order over the 'A' order. */
  mpz_init(M->kappa);
  mpz_pow_ui(tmp, cd, (d-1)*(d-2)/2);
  mpz_div(M->kappa, M->N->index, tmp); 
  
  /* Get an inert prime base for Chinese remaindering, and */
  /* initialize the CRT residues.                          */
  if (getIPB_mpz(M, 128.0) <= 0) {
    fprintf(stderr, "initMsqrt(): Some error occurred getting IPB!\n");
    return -1;
  }
  for (i=0; i<M->ipbSize; i++) {
    mpz_poly_init(M->crtRes[i]);
    mpz_set_ui(&M->crtRes[i]->coef[0], 1);
    M->crtRes[i]->degree = 0;
  }
  /* H is what was leftover from the previous step. That is,
     when we choose some delta\in I, chances are very good
     that <delta> has some prime factors not in I. Rather
     than bothering to factor what's left, we will just compute
     the ideal quotient (using the HNF's), and leave these
     extra factors in H, to be taken as part of I in the
     next step. Thus, before the first step, we should have
     H = <1>, whose HNF is the identity matrix.
     Note that there are actually 2 H's of course. We need
     to know whether the leftover stuff was in the numerator
     or denominator.
  */
  mpz_mat_init2(&M->Hnum, d, d); mpz_mat_setID(&M->Hnum, d);
  mpz_mat_init2(&M->Hden, d, d); mpz_mat_setID(&M->Hden, d);

  M->LLL_max_log = DEFAULT_LLLMAX_LOG;
  computeBeta(M);

  depSize=0;
  while (relsInDep[depSize]>=0)
    depSize++;
  /************************************************************/
  /* Setup M->AFB, so it has the regular AFB + large primes.  */
  /************************************************************/
  setupPrimes(M, lpF);
  if (!(M->aExp = (s32 *)malloc(M->aSize*sizeof(s32)))) {
    fprintf(stderr, "initMsqrt() memory allocation error for aExp!\n");
    return -1;
  }
  /* Allocate for exponents of the special primes. */
  k = M->spSize;
  if (!(M->spExp = (s32 *)malloc(k*sizeof(s32)))) {
    fprintf(stderr, "initMsqrt() memory allocation error for spExp!\n");
    return -1;
  }
  /* Compute some constants to help compute valuations at exceptional primes. */
  if (!(M->v_cd_sPrimes = (int *)malloc(k*sizeof(int)))) {
    fprintf(stderr, "initMsqrt() memory allocation error for v_cd_sPrimes!\n");
    return -1;
  }
  mpz_mat_setID(&H, d);
  for (i=0; i<k; i++) {
    /* Get the HNF of <c_d>. */
    for (j=0; j<d; j++)
      mpz_set(&H.entry[j][j], cd);
    /* Compute the valuation of the i-th special prime at <c_d>. */
    M->v_cd_sPrimes[i] = valuation(&H, &M->sPrimes[i], M->N);
  }

  printf("The zeros of f have been computed as:\n");
  for (i=0; i<d; i++)
    printf("z%" PRId32 " = %1.15lf + I*%1.15lf\n",i,M->N->fZeros[i].r, M->N->fZeros[i].i);

  /* Some precomputation to save work later. */
  for (i=0; i<d; i++)  {
    for (j=0; j<d; j++) {
      /*****************************************************/
      /* Compute \omega_i evaluated at the j-th zero of f. */
      /* This is needed to speed up other things.          */
      /*****************************************************/
      xr = xi = 0.0;
      zr = M->N->TZeros[j].r; zi = M->N->TZeros[j].i;
      zpr = 1.0; zpi = 0.0;
      for (k=0; k<=i; k++) {
        c = mpz_get_d(&M->N->W->entry[k][i]);
        xr += c*zpr; xi += c*zpi;
        tr = zpr*zr - zpi*zi; ti = zpr*zi + zpi*zr;
        zpr = tr; zpi = ti;
      }
      /* Divide by the denominator, and we're set. */
      c = mpz_get_d(M->N->W_d);
      xr /= c; xi /= c;
      M->N->WevalZ_r[i][j] = xr;
      M->N->WevalZ_i[i][j] = xi;
    }
  }

  /*****************************************************/
  /* Intialize everything. Zero out all the exponents, */
  /* the norm sizes,...                                */
  /*****************************************************/
  for (i=0; i<M->aSize; i++)
    M->aExp[i]=0;
  M->aExpLast=0;
  for (i=0; i<M->spSize; i++)
    M->spExp[i] = 0;
  M->logNormGNum = M->logNormGDen = 0.0;

  /*******************************************************/
  /* Now, compute gamma_1 <-- Prod(c_d*a_i - b_i*alpha), */
  /* where alpha is a root of T (the monic-ized version  */
  /* of the original polynomial).                        */
  /*******************************************************/
  for (j=0; j<d; j++)
    M->log_eps[j] = 0.0;
  mpz_set(cd, &M->FB->f->coef[d]);

  printf("Reading relations and computing initial <gamma> factorization...\n");
  printf("depSize = %" PRId32 ".\n", depSize);
  /* Prime the loop by opening the first relation file. */
  sprintf(fName, "%s.0", prelF->prefix);
  if (!(fp = fopen(fName, "rb"))) {
    fprintf(stderr, "initMsqrt() Fatal error: could not open %s for read!\n", fName);
    exit(-1);
  }
  printf("Reading relations from %s...\n", fName);
  R0 = 0;
  readRaw32(&R1, fp);
  Rindex = 0; fileNum = 0;
  e=-1;
  /* Throughout this loop: the current file has relations [R0, R1). */
  /* The relation fp is currently looking at is Rindex.             */
  mpz_set_ui(Zsquare, 1); numPairs = 0;
  for (i=0; i<depSize; i++) {
    rel = relsInDep[i];
    while (Rindex <= rel) {
      if (rel >= R1) {
        fclose(fp); fileNum++;
        sprintf(fName, "%s.%d", prelF->prefix, fileNum);
        if (!(fp = fopen(fName, "rb"))) {
          fprintf(stderr, "initMsqrt() Fatal error: could not open %s for read!\n", fName);
          exit(-1);
        }
        printf("Reading relations from %s...\n", fName);
        R0 = R1;
        readRaw32(&R1, fp); R1 += R0;
        Rindex = R0;
      }
      if (readRel(&R, fp)==0) {
        Rindex++;
      }
    }
    /* Finally, 'R' should have the proper relation. */
    a = R.a; b = -R.b;
    numPairs++;
#define _DEBUG_RAT_SQRT
#ifdef _DEBUG_RAT_SQRT
    mpz_set_ui(tmp2, 1);
    for (j=0; j<R.rFSize; j++) {
      mpz_ui_pow_ui(tmp, M->FB->rfb[2*R.rFactors[j]], R.rExps[j]);
      mpz_mul(tmp2, tmp2, tmp);
    }
    for (j=0; j<MAX_LARGE_RAT_PRIMES; j++) {
      mpz_mul_ui(tmp2, tmp2, R.p[j]);
    }
    mpz_set_si64(tmp, b);
    mpz_mul(tmp, tmp, M->FB->y0);
    mpz_set_si64(tmp3, a);
    mpz_mul(tmp3, tmp3, M->FB->y1);
    mpz_sub(tmp, tmp3, tmp);
    mpz_abs(tmp, tmp);
    mpz_abs(tmp2, tmp2);
    if (mpz_cmp(tmp2, tmp)) {
      s32 pFacts[128];
      int numpFacts;
      printf("Factorization of relation %" PRId32 " is wrong:\n", i);
      printf("a=%" PRId64 ", b=%" PRId64 "\n", a, -b);
      printf("a-bm = "); mpz_out_str(stdout, 10, tmp); printf("\n");
      printf("Stored large primes are: %" PRIu32 " %" PRIu32 ".\n", (u32)R.p[0], (u32)R.p[1]);
      printf("Product of factors gives:\n       ");
      mpz_out_str(stdout, 10, tmp2); printf("\n");
      numpFacts = factor(pFacts, tmp, 1);
      printf("factor() returned %d and :\n", numpFacts);
      for (j=0; j<numpFacts; j++)
        printf("%" PRId32 " ", pFacts[j]);
      printf("\n");
      exit(-1);
    }
#endif


    /***********************************************************/
    /* tpol1 must be given with respect to the integral basis: */
    /* tpol1 <-- (a+b\alpha) = (c_d*a + b*\hat{alpha})/c_d     */
    /*                       = c_d*a + b*bmultiplier* omega_1  */
    /***********************************************************/
    mpz_set_si64(&tpol1->coef[0], a); 
    mpz_mul(&tpol1->coef[0], &tpol1->coef[0], cd);
    mpz_set_si64(&tpol1->coef[1], b);
    mpz_mul(&tpol1->coef[1], &tpol1->coef[1], bmultiplier);
    tpol1->degree = 1;
//    e *= -1;
    /* Choose the exponent for this (a,b) pair, using a greedy strategy. */
    e = choose_ab_exponent(M, &R);
ABexponentSum += e;

    updateEps_ab(M, a, -b, e); 
    updateCRT(M, tpol1, cd, e);
    updateFactorization_ab(M, &R, e);
    /* Update the rational square root info. */
    if (ratSqrt(&R, e, depSize, M))
      return -1;
 
    /* Keep track of what the final square should be. */
#if 1
  /* This is the original code: */
    mpz_set_si64(tmp, b);
    mpz_mul(tmp, tmp, M->FB->y0);
    mpz_set_si64(tmp2, a);
    mpz_mul(tmp2, tmp2, M->FB->y1);
    mpz_sub(tmp, tmp2, tmp);
mpz_abs(tmp, tmp);
#else
    mpz_set_si64(tmp, b);
    mpz_mul(tmp, tmp, M->FB->y0);
    mpz_set_si64(tmp2, a);
    mpz_mul(tmp2, tmp2, M->FB->y1);
    mpz_sub(tmp, tmp2, tmp);
    
//mpz_abs(tmp, tmp);
#endif

    if (e == -1)
      mpz_invert(tmp, tmp, M->FB->n);
    mpz_mul(Zsquare, Zsquare, tmp);
    mpz_mod(Zsquare, Zsquare, M->FB->n);

  }
  printf("The final square should be: ");
  mpz_out_str(stdout, 10, Zsquare);
  printf("\nWe used %" PRId32 " (a,b) pairs.\n", numPairs);
  if (Rindex < R1) fclose(fp);
  i=M->aSize-1;
  while ((i>=0) && (M->aExp[i]==0))
    i--;
  M->aExpLast = i;

  for (i=0; i<=M->aExpLast; i++) {
    if (M->aExp[i]%2) {
      fprintf(stderr, "Error: Odd exponent found: AFB[%" PRId32 "] has exponent %" PRId32 "!\n",
              i, M->aExp[i]);
      exit(-1);
    }
  }
  for (i=0; i<M->spSize; i++) {
    if (M->spExp[i]%2) {
      fprintf(stderr, "Error: Odd exponent found: Sp[%" PRId32 "] has exponent %" PRId32 "!\n",
              i, M->spExp[i]);
      exit(-1);
    }
  }
  printf("Done. Computing rational square root...\n");
  if (ratSqrt(NULL, 0, 0, M))
    exit(-1);
  printf("Rational square root: "); mpz_out_str(stdout, 10, M->ratSqrt);
  printf("\n");
  mpz_mul(tmp, M->ratSqrt, M->ratSqrt);
  mpz_mod(tmp, tmp, M->FB->n);
  if (mpz_cmp(tmp, Zsquare)) {
    printf("Severe error: rational square root is broken!\n");
    printf("(M->ratSqrt)^2 = ");
    mpz_out_str(stdout, 10, tmp);
    printf("\nShould be :      ");
    mpz_out_str(stdout, 10, Zsquare);
    printf("\n");
    exit(-1);
  }


  M->logNormGNum = M->logNormGDen = 0.0;
  for (i=0; i<=M->aExpLast; i++) {
    e = M->aExp[i];
    if (e > 0)
      M->logNormGNum += e*log((double)M->AFB[i].p);
    else if (e < 0)
      M->logNormGDen += -e*log((double)M->AFB[i].p);
  }
  for (i=0; i<M->spSize; i++) {
    e = M->spExp[i]*M->sPrimes[i].f;
    if (e > 0)
      M->logNormGNum += e*_mpz_log(M->sPrimes[i].p); 
    else if (e < 0)
      M->logNormGDen += -e*_mpz_log(M->sPrimes[i].p); 
  }

  printf("Initialization done.\n");
  printf("logNormGNum = %1.5lf\n", M->logNormGNum);
  printf("logNormGDen = %1.5lf\n", M->logNormGDen);

  printf("CRT residues:\n");
  for (i=0; i<M->ipbSize; i++) {
    printf("q=");
    mpz_out_str(stdout, 10, &M->q[i]);
    printf(", residue: ");
    mpz_poly_print(stdout, "", M->crtRes[i]);
  }

  i=0;
  j=0;
  printf("There are %" PRId32 " exceptional prime ideals:\n", M->spSize);
  for (i=0; i<M->spSize; i++) {
    if (M->spExp[i] ) 
      printf("[Sp. %" PRId32 "]^%" PRId32 " * ", i, M->spExp[i]);
  }
  printf("\n");
  initOmegaEvalM(M); /* Initialize the evaluation constants. */

  mpz_clear(Zsquare); mpz_clear(tmp2);
  mpz_clear(cd); mpz_clear(tmp); mpz_clear(bmultiplier);
  mpz_poly_clear(tpol1);
  mpz_mat_clear(&H);

  return 0;
}

/***************************************************************************/
void updateEps_ab(msqrt_t *M, s64 a, s64 b, int exponent)
/***************************************************************************/
/* Update the embedding sizes with gamma <-- gamma*(a-b\alpha)^{exponent}. */
/***************************************************************************/
/* This is now correct. The initialization ends up with the correct        */
/* values for log_eps, so don't mess with this!
*/ 
{ int    i, d=M->N->degree;
  double xr, xi, c;

  for (i=0; i<d; i++) {
    xr = (double)a - (double)b*M->N->fZeros[i].r;
    xi = -(double)b*M->N->fZeros[i].i;
    c = 0.5*exponent*log(xr*xr + xi*xi);
    M->log_eps[i] += c;
  }
}

//#define _MP_EPS
/***********************************************************************/
void updateEps(msqrt_t *M, mpz_poly delta, int exponent)
/***********************************************************************/
/* Update the embedding sizes with gamma <-- gamma*delta^{exponent}    */
/* gamma is a polynomial in the \omega_i (the integral basis).         */
/***********************************************************************/
#ifdef _MP_EPS
{ int     i, j, d=M->N->degree;
  static mpf_t  xr, xi, c, tmp1;
  static int initialized=0;

  if (!(initialized)) {
    mpf_init2(xr, 256); mpf_init2(xi, 256); mpf_init2(c, 256);
    mpf_init2(tmp1, 256);
    initialized=1;
  }

  for (j=0; j<d; j++) {
    /* Compute delta(\sigma_j(alpha)). */
    mpf_set_ui(xr,0); mpf_set_ui(xi, 0);
    for (i=0; i<=delta->degree; i++) {
      /* x += delta->coef[i] * \omega_i(\sigma_j(alpha)). */
      mpf_set_z(c, &delta->coef[i]);
      mpf_set_d(tmp1, M->N->WevalZ_r[i][j]); mpf_mul(tmp1, tmp1, c); mpf_add(xr, xr, tmp1);
      mpf_set_d(tmp1, M->N->WevalZ_i[i][j]); mpf_mul(tmp1, tmp1, c); mpf_add(xi, xi, tmp1);
    }
    mpf_mul(c, xr, xr);
    mpf_mul(tmp1, xi, xi);
    mpf_add(c, c, tmp1);
    mpf_sqrt(c, c);
    M->log_eps[j] += log(mpf_get_d(c))*exponent;
  }
}
#else
{ int     i, j, d=M->N->degree;
  double  xr, xi, c;

  for (j=0; j<d; j++) {
    /* Compute delta(\sigma_j(alpha)). */
    xr=0.0; xi = 0.0;
    for (i=0; i<=delta->degree; i++) {
      /* x += delta->coef[i] * \omega_i(\sigma_j(alpha)). */
      c = mpz_get_d(&delta->coef[i]);
      xr += c*M->N->WevalZ_r[i][j];
      xi += c*M->N->WevalZ_i[i][j];
    }
    c = 0.5*log(xr*xr + xi*xi);
    M->log_eps[j] += c*exponent;
  }
}
#endif

/***********************************************************************/
void updateCRT(msqrt_t *M, mpz_poly delta, mpz_t denom, int exponent)
/***********************************************************************/
/* Update the CRT residues with gamma <-- gamma*delta^{exponent}       */
/* gamma is a polynomial in the \omega_i (the integral basis).         */
/***********************************************************************/
{ int    i, j;
  static mpz_poly tpol1, delta_a;
  static mpz_t    dinv, denom_a;
  static mpz_mat_t v;
  static int initialized=0;

  if (!(initialized)) {
    mpz_poly_init(tpol1); mpz_poly_init(delta_a);
    mpz_init(dinv); mpz_init(denom_a); mpz_mat_init(&v);
    initialized=1;
  }

  /* Convert delta to the \hat{\alpha} representation. */
  setCol(&v, 0, delta);
  v.rows = M->N->degree;
  v.cols = 1;
  mpz_mat_mul(&v, M->N->W, &v);
  getCol(delta_a, &v, 0);
  mpz_mul(denom_a, M->N->W_d, denom);  


  for (j=0; j<M->ipbSize; j++) {
    if (exponent < 0) {
      mpz_poly_inv(tpol1, delta_a, M->N->T, &M->q[j]);
      for (i=0; i<=tpol1->degree; i++) {
        mpz_mul(&tpol1->coef[i], &tpol1->coef[i], denom_a);
        mpz_mod(&tpol1->coef[i], &tpol1->coef[i], &M->q[j]);
      }
      for (i=0; i<-exponent; i++) 
        mpz_poly_mulmod_pp(M->crtRes[j], M->crtRes[j], tpol1, M->N->T, &M->q[j]);
    } else {
       mpz_invert(dinv, denom_a, &M->q[j]);
       for (i=0; i<=delta_a->degree; i++) {
         mpz_mul(&tpol1->coef[i], &delta_a->coef[i], dinv);
         mpz_mod(&tpol1->coef[i], &tpol1->coef[i], &M->q[j]);
       }
       tpol1->degree = delta_a->degree;
       for (i=0; i<exponent; i++)
         mpz_poly_mulmod_pp(M->crtRes[j], M->crtRes[j], tpol1, M->N->T, &M->q[j]);
    }
  }
}

/************************************************************************/
void idealDivExact(mpz_mat_t *Res, mpz_mat_t *X, mpz_mat_t *Y, nf_t *N)
/************************************************************************/
{ static mpz_t g, tmp, tmp2, normX, normY;
  static mpz_mat_t I, J, T;
  static int initialized=0;
  int    i, n=N->degree;

  if (!(initialized)) {
    mpz_init(g); mpz_init(tmp); mpz_init(tmp2);
    mpz_init(normX); mpz_init(normY);
    mpz_mat_init(&I);
    mpz_mat_init(&J);
    mpz_mat_init(&T);
    initialized=1;
  }

  mpz_set_ui(normX, 1); mpz_set_ui(normY, 1);
  for (i=0; i<n; i++) {
    mpz_mul(normX, normX, &X->entry[i][i]);
    mpz_mul(normY, normY, &Y->entry[i][i]);
  }

  mpz_gcd(g, normX, normY);
  do {
    mpz_div(tmp, normX, g);
    mpz_gcd(tmp2, tmp, g);
    mpz_div(g, g, tmp2);
  } while (mpz_cmp_ui(tmp2, 1));

  mpz_mat_setID(&T, n);
  mpz_div(tmp, normX, g);

  for (i=0; i<n; i++)
    mpz_set(&T.entry[i][i], tmp);

  mpz_set_ui(tmp, 1);
  mpz_mat_moduleAdd(&I, tmp, X, tmp, &T, tmp);
  mpz_mat_setID(&T, n);
  mpz_div(tmp, normY, g);

  for (i=0; i<n; i++)
    mpz_set(&T.entry[i][i], tmp);

  mpz_set_ui(tmp, 1);
  mpz_mat_moduleAdd(&J, tmp, Y, tmp, &T, tmp);
    
  idealDiv(Res, &I, &J, N);
}


/************************************************************************/
void idealDiv(mpz_mat_t *Res, mpz_mat_t *I, mpz_mat_t *J, nf_t *N)
/************************************************************************/
/* It is assumed that this function will be used for the same number    */
/* field every time it's called. So we'll save some cost by remembering */
/* some of what is computed the first time.                             */
/************************************************************************/
/* I think we can actually do this much faster in our case since I is   */
/* <delta>. But for now, slow will do just fine.                        */
/************************************************************************/
{ int    i, j, n=N->degree;
  static mpz_mat_t T, T_inv, H, P;
  static mpz_t     d, tmp1, T_invden, Pden;
  static mpz_poly  g_i, d_j, r;
  static int initialized=0;


  /* First, use Cohen Alg. 4.8.21 to invert J. */
  if (!(initialized)) {
    mpz_mat_init2(&T, n, n); mpz_mat_init2(&T_inv, n, n);
    mpz_mat_init2(&P, n, n);
    mpz_mat_init2(&H, n, n*n);
    mpz_init(d); mpz_init(tmp1); mpz_init(T_invden); mpz_init(Pden);
    mpz_poly_init(g_i); mpz_poly_init(d_j); mpz_poly_init(r);
    mpz_set(d, N->Kdisc);

    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
        getTrace_ib(&T.entry[i][j], N->Mt[i*n+j], N);
    T.rows = T.cols = n;
    mpz_mat_pseudoInvert(&T_inv, T_invden, &T);
    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
        mpz_mul(&T_inv.entry[i][j], &T_inv.entry[i][j], d);
    mpz_set(tmp1, T_invden);
    for (i=0; i<n; i++)
      for (j=0; j<n; j++) 
        mpz_gcd(tmp1, tmp1, &T_inv.entry[i][j]);
    mpz_div(T_invden, T_invden, tmp1);
    for (i=0; i<n; i++)
      for (j=0; j<n; j++) 
        mpz_div(&T_inv.entry[i][j], &T_inv.entry[i][j], tmp1);
    initialized=1;
  }

  for (i=0; i<n; i++) {
    getCol(g_i, J, i);
    for (j=0; j<n; j++) {
      getCol(d_j, &T_inv, j);
      basisPolyMult(r, g_i, d_j, N->Mt, n);
      setCol(&H, i*n+j, r);
    }
  }
  H.rows = n; H.cols = n*n;
  mpz_mat_getHNF(&H, &H);
  /* Now, we need the transpose. */
  for (i=0; i<n; i++)
    for (j=i+1; j<n; j++)
      mpz_swap(&H.entry[i][j], &H.entry[j][i]);

  mpz_mat_pseudoInvert(&P, Pden, &H);
  mpz_mat_mul(&P, &T_inv, &P);
  mpz_mul(Pden, Pden, T_invden);
  /* Finally, check the denominator: */
  mpz_set(tmp1, Pden);
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++)
      mpz_gcd(tmp1, tmp1, &P.entry[i][j]);
    if (mpz_cmp_ui(tmp1, 1)==0)
      break;
  }
  if (mpz_cmp_ui(tmp1, 1)) {
    mpz_div(Pden, Pden, tmp1);
    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
        mpz_div(&P.entry[i][j], &P.entry[i][j], tmp1);
  }
  mpz_mat_getHNF(&P, &P);

  /* Now all we need is the product of I by P, Pden. */
  mpz_mat_cp(&H, I);
  idealHNFMul_ib(Res, &P, &H, N);
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
      mpz_tdiv_qr(&Res->entry[i][j], tmp1, &Res->entry[i][j], Pden);
      if (mpz_sgn(tmp1)) {
        fprintf(stderr, "Error: expected divisibility is not there!\n");
        exit(-1);
      }
    }
  }

}


/***********************************************************************/
void updateFactorization(msqrt_t *M, mpz_poly delta, mpz_mat_t *I, int sl)
/***********************************************************************/
/* Update the factorization of <gamma_l> with                          */
/* gamma <-- gamma*delta^{-2*sl}, where                                */
/* gamma is a polynomial in the \omega_i (the integral basis).         */
/* It is assumed that I is the ideal from which delta was chosen and   */
/* that the factorization has already been updated with                */
/* gamma <-- gamma*I^{-2*sl}, so that all we really need to do is find */
/* the quotient <delta>/I and put this in Hnum or Hden as appropriate. */
/***********************************************************************/
{ int    d=M->N->degree;
  static mpz_mat_t Delta, NewH;
  static int initialized=0;

  if (!(initialized)) {
    mpz_mat_init2(&Delta, d, d);
    mpz_mat_init2(&NewH, d, d);
    initialized=1;
  }
  idealHNF_ib(&Delta, delta, M->N);

  /* Update the embedding sizes. */
  updateEps(M, delta, -2*sl);
  if (sl > 0) {
    /* We are simplifying the numerator, by doing
       <gamma> <-- <gamma>/I^2. We try to accomplish this with:
        gamma  <--  gamma/delta^2, so if <delta>=IJ, we get
       <gamma>/<delta>^2 = <gamma>/(I^2J^2).
       Thus, find J and put it in Hden, so it can be cancelled the
       next time we work on the denominator.
    */
    idealDivExact(&NewH, &Delta, I, M->N);
    idealHNFMul_ib(&M->Hden, &NewH, &M->Hden, M->N);
    M->logNormGDen += 2*logNorm(&NewH, M);
  } else {     
    idealDivExact(&NewH, &Delta, I, M->N);
    idealHNFMul_ib(&M->Hnum, &NewH, &M->Hnum, M->N);
    M->logNormGNum += 2*logNorm(&NewH, M);
  }
}
       
/***********************************************************************/
int choose_ab_exponent(msqrt_t *M, relation_t *R)
/***********************************************************************/
/* Choose an exponent for this (a,b) pair greedily.                    */
/***********************************************************************/
{ int  i;
  s64 pLoc;
  s32 e, p, r;
  double netNum=0.0, netDen=0.0, l;

  /* First, handle the `regular' factors. */
  for (i=0; i<R->aFSize; i++) {
    e = R->aExps[i];
    pLoc = R->aFactors[i];
    l = log((double)M->AFB[pLoc].p);
    netNum += fabs((M->aExp[pLoc] + e)*l);
    netDen += fabs((M->aExp[pLoc] - e)*l);
  }

  for (i=0; i<R->spSize; i++) {
    e = R->spExps[i];
    pLoc = R->spFactors[i];
    l = _mpz_log(M->sPrimes[i].p);
    netNum += fabs((M->spExp[i] + e)*l);
    netDen += fabs((M->spExp[i] - e)*l);
  }  

  /* And the large primes: */
  for (i=0; i<MAX_LARGE_ALG_PRIMES; i++) {
    p = R->a_p[i]; r = R->a_r[i];
    if (p > 1) {
      pLoc = findLP(p, r, M);
      if ((M->AFB[pLoc].p != p)||(M->AFB[pLoc].r != r))  {
        fprintf(stderr, "Error: Couldn't find (%" PRId32 ", %" PRId32 ") in the primes list!\n", p, r);
        exit(-1);
      }
      l = log((double)M->AFB[pLoc].p);
      netNum += fabs((M->aExp[pLoc] + 1)*l);
      netDen += fabs((M->aExp[pLoc] - 1)*l);
    }
  }
  if (netDen < netNum) return -1;
  return 1;
}


/***********************************************************************/
void updateFactorization_ab(msqrt_t *M, relation_t *R, int exponent)
/***********************************************************************/
/* Update the factorization of <gamma_l> with                          */
/* gamma <-- gamma*(a - b\alpha)^exponent, where                       */
/* alpha = RootOf(M->FB->f), and exponent=+/-1.                        */
/* We will use the factorization from R, so it must have already been  */
/* put through factRel()!                                              */
/***********************************************************************/
{ int  i;
  s64 pLoc;
  s32 e, p, r;

  /* First, handle the `regular' factors. */
  for (i=0; i<R->aFSize; i++) {
    e = R->aExps[i];
    pLoc = R->aFactors[i];
    M->aExp[pLoc] += exponent*e;
  }

  for (i=0; i<R->spSize; i++) {
    pLoc = R->spFactors[i];
    M->spExp[pLoc] += exponent*R->spExps[i];
  }  

  /* And the large primes: */
  for (i=0; i<MAX_LARGE_ALG_PRIMES; i++) {
    p = R->a_p[i]; r = R->a_r[i];
    if (p > 1) {
      pLoc = findLP(p, r, M);
      if ((M->AFB[pLoc].p != p)||(M->AFB[pLoc].r != r))  {
        fprintf(stderr, "Error: Couldn't find (%" PRId32 ", %" PRId32 ") in the primes list!\n", p, r);
        exit(-1);
      }
      M->aExp[pLoc] += exponent;
    }
  }
}


/*********************************************************************/
void getIdealHNF_ib_pre(mpz_mat_t *Ipe, s32 p, s32 r, s32 e, msqrt_t *M)
/*********************************************************************/
/* Get the HNF, w.r.t. the integral basis, of the 1st-degree prime   */
/* ideal corresponding to (p,r) with exponent 'e'.                   */
/* It is assumed that this ideal is not exceptional!                 */
/*********************************************************************/
{ static mpz_mat_t I, Tmp;
  static mpz_poly  alpha;
  static mpz_t     p_, tmp, bmultiplier;
  static int initialized=0;
  int    i, j, d=M->N->degree;

  if (!(initialized)) {
    mpz_mat_init2(&I, d, d);
    mpz_mat_init2(&Tmp, d, d);
    mpz_poly_init(alpha);
    mpz_init(p_); mpz_init(tmp);
    mpz_init_set(bmultiplier, M->N->W_d);
    mpz_div(bmultiplier, bmultiplier, &M->N->W->entry[1][1]);
    initialized=1;
  }

  if (r != p) {
    /***************************************************/
    /* So r is finite. If p does not divide the index, */
    /* this is easy.                                   */
    /***************************************************/
    if (mpz_tdiv_ui(M->N->index, p)) {
      mpz_set_ui(p_, p);
      mpz_set_si(&alpha->coef[0], -r);
      mpz_mul(&alpha->coef[0], &alpha->coef[0], &M->FB->f->coef[d]);
      mpz_mod(&alpha->coef[0], &alpha->coef[0], p_);
      mpz_set(&alpha->coef[1], bmultiplier);
      alpha->degree = 1;
      getIdealHNF_ib(Ipe, p_, alpha, e, M->N);
      return;
    }
    /************************************************************/
    /* p divides the index, so it is one of the M->n->sPrimes.  */
    /* It could be the case that it's not exceptional, but this */
    /* seems not to happen for some reason. (?? why not ??).    */
    /************************************************************/
    fprintf(stderr, "getIdealHNF_ib_pre() : Error - cannot handle ramified primes!\n");
    mpz_mat_setID(Ipe, d);
    return ;
#ifdef _TRY_RAMIFIED
    for (i=0; i<M->N->numSPrimes; i++) {
      if ((mpz_cmp_ui(M->N->sPrimes[i].p, p)==0) && (M->N->sPrimes[i].f==1)) {
        /* The norm and degree are right. Check if this is really the one. */
        /* This remains to be done. But I don't think it is even needed.   */
      }
    }
#endif

  } else {
    /*****************************************/
    /* r==p, so this is a prime (p, \infty). */
    /*****************************************/
    /* I think this is right.
    */
    mpz_set_ui(&I.entry[0][0], p);
    for (i=1; i<d; i++)
      mpz_set_ui(&I.entry[i][0], 0);
    for (j=1; j<d; j++) 
      for (i=0; i<d; i++)
        mpz_set(&I.entry[i][j], &M->Beta.entry[i][j-1]);
    I.cols = I.rows = d;
    mpz_mat_getHNF(&I, &I);
    idealHNFPow_ib(Ipe, &I, e, M->N);
  }

} 



/**************************************************/
double logNorm(mpz_mat_t *HNF, msqrt_t *M)
/**************************************************/
/* Given a matrix in HNF, return log(|det(HNF)|). */
/**************************************************/
{ double res=0;
  int i;
  static mpz_t t;
  static int initialized=0;

  if (!(initialized)) {
    mpz_init(t);
    initialized=1;
  }

  for (i=HNF->rows-1; i>=0; i--) {
    mpz_abs(t, &HNF->entry[i][i]);
    res += _mpz_log(t);
  }
  return res;
}

/*********************************************************************/
int chooseIdeal(mpz_mat_t *I, msqrt_t *M, int sl)
/*********************************************************************/
/* Choose an ideal I with I^2 dividing                               */
/*             Num(<gamma_l>) and (Hnum)^2 if sl=1, or               */
/*             Den(<gamma_l>) and (Hden)^2 if sl=-1.                 */
/* The ideal should have norm close to exp(M->LLL_max_log).          */
/* At the same time, update the factorization of <gamma_l> with      */
/* <gamma_l> <-- <gamma_l>*I^{-2*sl}.                                */
/* We prefer to do this update here because we already know the      */
/* factorization of I, so it's more convenient. Later, after delta   */
/* is actually chosen and we need to multiply by <delta>^2, all that */
/* must be done is to see how <delta> differs from I and put that    */
/* in Hnum or Hden as appropriate (i.e., H <-- <delta>/I.)           */
/* Important: 'I' will be an ideal w.r.t the integral basis, not the */
/* standard basis!                                                   */
/*********************************************************************/
{ double    lognormI, dP;
  s32      p, r, index, indexE, e=0;
  int       i, d=M->N->degree, cont, nextIndex;
  static    mpz_mat_t Ipe;
  static    int initialized=0;
#ifdef _DEBUG
  char tmpStr[1024];
#endif

  if (!(initialized)) {
    mpz_mat_init2(&Ipe, d, d);
    initialized=1;
  }
#ifdef _DEBUG
  idealSelStr[0]=0;
#endif

  if (sl == 1) {
    mpz_mat_cp(I, &M->Hnum);
    M->logNormGNum -= 2.0*logNorm(&M->Hnum, M);
    mpz_mat_setID(&M->Hnum, d);
#ifdef _DEBUG
    sprintf(tmpStr, "M->Hnum added.\n");
    strncat(idealSelStr, tmpStr, MAX_IDEAL_STR-strlen(idealSelStr));
#endif
  } else {
    mpz_mat_cp(I, &M->Hden);
    M->logNormGDen -= 2.0*logNorm(&M->Hden, M);
    mpz_mat_setID(&M->Hden, d);
#ifdef _DEBUG
    sprintf(tmpStr, "M->Hden added.\n");
    strncat(idealSelStr, tmpStr, MAX_IDEAL_STR-strlen(idealSelStr));
#endif
  }
  /******************************************************************/
  /* Hnum and Hden are in HNF with respect to the integral basis.   */
  /* So, we obtain the norm by multiplying the diagonal entries,    */
  /* and multiplying by det(omega_j basis) = M->N->index.           */
  /******************************************************************/
  lognormI = _mpz_log(M->N->index);
  for (i=0; i<d; i++)
    lognormI += _mpz_log(&I->entry[i][i]);

  index = M->aExpLast;
  cont = 1;



  while (cont) {
    if (sl==1) {
      while ((index>=0) && (M->aExp[index] <= 0))
        index--;
      if (index >= 0) {
        /* Select this ideal. */
        p = M->AFB[index].p;
        r = M->AFB[index].r;
        if (M->aExp[index]%2) {
          printf("chooseIdeal() sever error: odd exponent found: (%" PRId32 ", %" PRId32 ") e=%" PRId32 "!\n",p,r,(s32)M->aExp[index]);
          return -1;
        }
        e=0;
        while ((2*e < M->aExp[index]) && (lognormI < M->LLL_max_log)) {
          e++;
          lognormI += log((double)p);
        }
#ifdef _DEBUG
        sprintf(tmpStr, "Ideal %ld : (%ld, %ld) e=%ld/%ld\n",index,p,r,e,M->aExp[index]);
        strncat(idealSelStr, tmpStr, MAX_IDEAL_STR-strlen(idealSelStr));
#endif
        /* Now we will actually use the (p,r) ideal with exponent e. */
        if (e > 0) {
          M->aExp[index] -= 2*e;
          getIdealHNF_ib_pre(&Ipe, p, r, e, M);
          idealHNFMul_ib(I, I, &Ipe, M->N);
          /* Take it out of the norm now: updateFactorization() expects this to be done here! */
          M->logNormGNum -= 2*e*log((double)p);
          cont = (20+lognormI < M->LLL_max_log);
        } else cont=0;
      } else {
        /* There are no (p,r) ideals left. How about exceptional ideals? */
        indexE=M->spSize-1;
        while ((indexE >=0) && (M->spExp[indexE] <= 0))
          indexE--;
        if (indexE < 0)
          cont = 0;
        else {
          if (M->spExp[indexE]%2) {
            printf("chooseIdeal() severe error: odd exponent found! (Special ideal %" PRId32 ", e=%" PRId32 ")\n", indexE,M->spExp[indexE]);
            return -1;
          }
          e=0;
          dP = mpz_get_d(M->sPrimes[indexE].p); 
          while ((2*e < M->spExp[indexE]) && (lognormI < M->LLL_max_log)) {
            e++;
            lognormI += _mpz_log(M->sPrimes[indexE].p);
          }
#ifdef _DEBUG
        sprintf(tmpStr, "EIdeal %" PRId32 " : e=%" PRId32 "/%" PRId32 "\n",indexE,e,M->spExp[indexE]);
        strncat(idealSelStr, tmpStr, MAX_IDEAL_STR-strlen(idealSelStr));
#endif
          M->spExp[indexE] -= 2*e;
          getIdealHNF_ib(&Ipe, M->sPrimes[indexE].p, M->sPrimes[indexE].alpha, e, M->N);
          idealHNFMul_ib(I, I, &Ipe, M->N);
          /* Take it out of the norm now: updateFactorization() expects this to be done here! */
          M->logNormGNum -= 2*e*log(dP)*M->sPrimes[indexE].f;
          cont = (10+lognormI < M->LLL_max_log);
        } 
      }
    } else { /* sl == -1. */
      while ((index>=0) && (M->aExp[index] >= 0))
        index--;
      if (index >= 0) {
        /* Select this ideal. */
        p = M->AFB[index].p;
        r = M->AFB[index].r;
        if (M->aExp[index]%2) {
          printf("chooseIdeal() severe error: odd exponent found: (%" PRId32 ", %" PRId32 ") e=%" PRId32 "!\n",p,r,(s32)M->aExp[index]);
          return -1;
        }
        e=0;
        while ((2*e < -M->aExp[index]) && (lognormI < M->LLL_max_log)) {
          e++;
          lognormI += log((double)p);
        }
#ifdef _DEBUG
        sprintf(tmpStr, "Ideal %" PRId32 " : (%" PRId32 ", %" PRId32 ") e=%" PRId32 "/%" PRId32 "\n",index,p,r,e,M->aExp[index]);
        strncat(idealSelStr, tmpStr, MAX_IDEAL_STR-strlen(idealSelStr));
#endif
        if (e > 0 ) {
          /* Now we will actually use the (p,r) ideal with exponent e. */
          M->aExp[index] += 2*e;
          getIdealHNF_ib_pre(&Ipe, p, r, e, M);
          idealHNFMul_ib(I, I, &Ipe, M->N);
          M->logNormGDen -= 2*e*log((double)p);
          cont = (10+lognormI < M->LLL_max_log);
        } else cont=0;
      } else {
        /* There are no (p,r) ideals left. How about special ideals? */
        indexE=M->spSize-1;
        while ((indexE >=0) && (M->spExp[indexE] >= 0))
          indexE--;
        if (indexE < 0) 
          cont=0;
        else {
          if (M->spExp[indexE]%2) {
            printf("chooseIdeal() severe error: odd exponent found! (Special ideal %" PRId32 ", e=%" PRId32 ")\n", indexE,M->spExp[indexE]);
            return -1;
          }
          e=0;
          dP = mpz_get_d(M->sPrimes[indexE].p); 
          while ((2*e < -M->spExp[indexE]) && (lognormI < M->LLL_max_log)) {
            e++;
            lognormI += _mpz_log(M->sPrimes[indexE].p);
          }
#ifdef _DEBUG
        sprintf(tmpStr, "EIdeal %ld : e=%ld/%ld\n",indexE,e,M->spExp[indexE]);
        strncat(idealSelStr, tmpStr, MAX_IDEAL_STR-strlen(idealSelStr));
#endif
          M->spExp[indexE] += 2*e;
          getIdealHNF_ib(&Ipe, M->sPrimes[indexE].p, M->sPrimes[indexE].alpha, e, M->N);
          idealHNFMul_ib(I, I, &Ipe, M->N);
          /* Take it out of the norm now: updateFactorization() expects this to be done here! */
          M->logNormGDen -= 2*e*log(dP)*M->sPrimes[indexE].f;
          cont = (10+lognormI < M->LLL_max_log);
        } 
      }
    }
    if (e>0) {
      /* Don't do another ideal over the same prime! */
      nextIndex=index;
      while ((nextIndex>=0) && (M->AFB[nextIndex].p == M->AFB[index].p))
        nextIndex--;
      index=nextIndex;
    }
  }
  return 0;
}

/*********************************************************************/
void computeEmbedding_ib(double *sr, double *si, mpz_poly v, int i, msqrt_t *M)
/*********************************************************************/
/* Input: v is a polynomial in the integral basis, omega_j.          */
/*        i is an integer 0<=i<d specifying a zero of f.             */
/* Output: (sr, si) <-- sigma_i(v).                                  */
/*        This is an approximation, which could possible go way off  */
/*        under rare circumstances (i.e., if there is alot of        */
/*        additive cancellation).                                    */
/*********************************************************************/
{ long double tr, ti, c;
  int    j;

  tr=ti=0.0;
  for (j=0; j<=v->degree; j++) {
    c = mpz_get_d(&v->coef[j]);
    tr += M->N->WevalZ_r[j][i]*c;
    ti += M->N->WevalZ_i[j][i]*c;
  }
  *sr = (double)tr;
  *si = (double)ti;
}


/*********************************************************************/
int chooseDelta(mpz_poly delta, mpz_mat_t *I, int sl, msqrt_t *M)
/*********************************************************************/
/* I is the HNF of the ideal from which delta should be chosen. So,  */
/* we will do an LLL reduction, add the `embedding-bounding' entries,*/
/* do another LLL, and choose an element of the ideal.               */
/*********************************************************************/
{ int    i, j, d=M->N->degree, lowestCol;
  double c, lambda[MAXPOLYDEGREE], zi, entry, sr, si, lowestWt, thisWt, logNormI;
  static mpz_mat_t H, U;
  static mpz_poly  v;
  static mpz_t     tmp;
  static int initialized=0;

  if (!(initialized)) {
    mpz_mat_init2(&H, 2*d, d);
    mpz_mat_init2(&U, 2*d, d);
    mpz_poly_init(v); mpz_init(tmp);
    initialized=1;
  }

  logNormI = logNorm(I, M);
  /* LLL reduction on I. */
//  mpz_mat_LLL(&H, &U, I);
  mpz_mat_LLL(&H, NULL, I);

  /*******************************************/
  /* Now add the embedding-bounding entries. */
  /* First, compute the constant 'c' needed  */
  /* for all entries.                        */
  /*******************************************/
  c = sl*fabs((M->logNormGNum - M->logNormGDen));
  mpz_abs(tmp, M->N->Kdisc);
  c -= _mpz_log(tmp);
  c = 0.5*c + M->LLL_max_log - logNormI;
  c = c/(double)d; /* Leave c logarithmic for now. */
  /* Now the lambda_j : */
  for (j=0; j<d; j++) {
    lambda[j] = c - 0.5*sl*M->log_eps[j];
  }

  /* Finally, we're ready to add the needed entries: */
  for (i=0; i<d; i++) {
    /* Is the i-th root of f real or complex? */
    zi = M->N->WevalZ_i[1][i]; /* Im(omega_1 at the i-th zero) */
    if (fabs(zi) < 0.000001) {
      /* Quite probably a real root. */
      for (j=0; j<d; j++) {
        getCol(v, &H, j);
        computeEmbedding_ib(&sr, &si, v, i, M); /* (sr, si) <-- sigma_i(v). */
        entry = exp(lambda[i] + log(fabs(sr))); if (sr < 0) entry *= -1.0;
        patched_mpz_set_d(&H.entry[d+i][j], entry);
      }
    } else {
      /* Quite probably a complex root. */
      for (j=0; j<d; j++) {
        getCol(v, &H, j);
        computeEmbedding_ib(&sr, &si, v, i, M); /* (sr, si) <-- sigma_i(v). */
        entry = exp(lambda[i]) * sr * M_SQRT2;
        patched_mpz_set_d(&H.entry[d+i][j], entry);
        entry = exp(lambda[i]) * si * M_SQRT2;
        patched_mpz_set_d(&H.entry[d+i+1][j], entry);
      }
      /***************************************************/
      /* Since the roots have been ordered, the next one */
      /* should be the conjugate of this one, so we      */
      /* should skip over it.                            */
      /***************************************************/
      i++;
    }
  }
  H.rows += d;  

//  mpz_mat_LLL(&H, &U, &H);
  mpz_mat_LLL(&H, NULL, &H);
  
  /* Now, select Delta to be the lowest weight column of H. */
  lowestCol=-1; lowestWt=-1;
  for (i=0; i<H.cols; i++) {
    thisWt=0.0;
    for (j=0; j<H.rows; j++) {
      c = mpz_get_d(&H.entry[j][i]);
      thisWt += c*c;
    }
    if ((lowestCol <0) || (thisWt < lowestWt)) {
      lowestCol = i;
      lowestWt = thisWt;
    }
  }
  /* Retrieve the part of delta corresponding to the element itself. */
  H.rows = d;
  getCol(delta, &H, lowestCol);
  /******************************************************************/
  /* CJM, 2/3/04.                                                   */
  /* We might have chosen an element whose norm is negative. If so, */
  /* we'll need to multiply it by -1.                               */
  /******************************************************************/
  /* CJM: 11/29/04: This seems fishy. Is it really necessary?!?     */
  norm_ib(tmp, delta, M->N);
#if 0
  if (mpz_sgn(tmp) < 0) {
    for (i=0; i<=delta->degree; i++)
      mpz_neg(&delta->coef[i], &delta->coef[i]);
  }
#endif
  return 0;
}

/*********************************************************************/
int crtLift(mpz_poly gamma, mpz_t denom, msqrt_t *M)
/*********************************************************************/
/* Lift the CRT residues to an algebraic integer. Hopefully, things  */
/* worked well enough that the resulting integer, gamma, is actually */
/* what's left of the real gamma = \prod (a-b\alpha).                */
/*********************************************************************/
{ mpz_t m, t1, t2, c, W_d;
  int   i, j, d=M->N->degree;

  mpz_init(t1); mpz_init(t2); mpz_init(c);
  mpz_init_set_ui(m, 1); mpz_init_set(W_d, M->N->W_d);
  /* Remember: The residues are in terms of \hat{alpha}, */
  /* a root of the monic polynomial. Thus, we will first */
  /* multiply through by the common denominator of the   */
  /* integral basis, to get just a numerator.            */
 

  for (i=0; i<d; i++)
    mpz_set_ui(&gamma->coef[i], 0);
  gamma->degree = d-1; /* To be adjusted later. */

  for (i=0; i<M->ipbSize; i++)
    mpz_mul(m, m, &M->q[i]);
  
  for (j=0; j<M->ipbSize; j++) {
    mpz_div(t1, m, &M->q[j]);
    mpz_invert(c, t1, &M->q[j]);
    mpz_mul(t2, t1, c);
    for (i=0; i<=M->crtRes[j]->degree; i++) {
      mpz_mul(t1, t2, &M->crtRes[j]->coef[i]);
      mpz_mul(t1, t1, W_d);
      mpz_add(&gamma->coef[i], &gamma->coef[i], t1);
      mpz_mod(&gamma->coef[i], &gamma->coef[i], m);
    }
  }
  mpz_poly_fixDeg(gamma);

  mpz_div_2exp(t1, m, 1);
  for (i=0; i<=gamma->degree; i++) {
    if (mpz_cmp(&gamma->coef[i], t1)>=0)
      mpz_sub(&gamma->coef[i], &gamma->coef[i], m);
  }

  mpz_set(denom, W_d);
  mpz_set(t1, denom);
  for (i=0; i<=gamma->degree; i++)
    mpz_gcd(t1, t1, &gamma->coef[i]);
  mpz_div(denom, denom, t1);
  for (i=0; i<=gamma->degree; i++)
    mpz_div(&gamma->coef[i], &gamma->coef[i], t1);

  mpz_clear(t1); mpz_clear(t2); mpz_clear(c);
  mpz_clear(m);
  return 0;
}


/*********************************************************************/
int finalSquareRoot(mpz_poly beta, mpz_poly gamma, msqrt_t *M)
/*********************************************************************/
{ mpz_t tmp;
  int   res;

  mpz_init(tmp);
  if (gamma->degree == 0) {
    beta->degree = 0;
    mpz_sqrtrem(&beta->coef[0], tmp, &gamma->coef[0]);
    if (mpz_sgn(tmp)) {
      printf("Error: Something failed horribly: gamma is not a square!\n");
      mpz_clear(tmp);
      return -1;
    }
    mpz_clear(tmp);
    return 0;
  }
  mpz_set_ui(tmp, 0);
  res = Zalpha_sqrt(beta, gamma, M->N->T, M->FB->n, tmp);
  mpz_clear(tmp);
  return res;
}

/*********************************************************************/
void evalOmegaPoly(mpz_t eval, mpz_poly h, mpz_t x, mpz_t modulus, msqrt_t *M)
/*********************************************************************/
/* Given a polynomial, h, in the integral basis, evaluate it at      */
/* \hat{\alpha} = x, mod modulus.                                    */
/*********************************************************************/
{ int       i;
  mpz_t     tmp1, res;

  mpz_init(tmp1); mpz_init_set_ui(res, 0); 
  for (i=0; i<=h->degree; i++) {
    mpz_mul(tmp1, &h->coef[i], &M->omegaEvalM[i]);
    mpz_add(res, res, tmp1);
    mpz_mod(res, res, M->FB->n);
  }
  mpz_set(eval, res);
  mpz_clear(tmp1); mpz_clear(res); 
}
  

/*********************************************************************/
/*
   STEN: g_M was previously defined inside montgomerySqrt() function as
         local variable. However, msqrt_t type is too huge to be placed
		 on stack. I decided to make it a global variable, thus making
		 montgomerySqrt() function non-reentrant (it cann't be called from
		 separate threads simultaneously).
*/
static msqrt_t g_M; 

int montgomerySqrt(mpz_t rSqrt, mpz_t aSqrt, s32 *relsInDep, multi_file_t *prelF, 
                   multi_file_t *lpF, nfs_fb_t *FB, nf_t *N)
/*********************************************************************/
/* Compute the square root in Z[\theta] from the dependence relsInDep*/
/* on the relations R in prelF.                                      */
/* relsInDep is stored as a list of indices to the 'R', terminated   */
/* by a negative number.                                             */
/* Rather than return the square root, we simply apply the           */
/* homomorphism to Z/NZ and return the resulting value there.        */
/*********************************************************************/
/* N should have all its fields filled in; in particular, we need    */
/* to already have an integral basis, the field discriminant, the    */
/* index (these are all done by getIntegralBasis()) and the zeros.   */
/*********************************************************************/
/* STEN: Note that function is not reentrant as it uses some global  */
/*       variables (see static definitons above).                    */
/*********************************************************************/

{ s32          l;
  int           i, d = N->degree, retVal, sl, cont, finalTries;
  mpz_poly      delta, beta, gamma;
  mpz_mat_t     I;
  mpz_t         tmp1, tmp2, gamma_denom, cdm, res, tmp3;
  double        normA, normN, now, lastReportTime;
  double        diff, lastDiff;
  double        lastNormA, lastNormN, lastNormNum, lastNormDen;
  double        lastHNormNum, lastHNormDen, improve, newTotal;
  double        totalSize[20];
#ifdef _LOUD_DEBUG
  FILE *ofp;
#endif
#ifdef _DEBUG
  int empty;
#endif

  /* Init variables. */
  mpz_poly_init(delta); mpz_poly_init(beta); mpz_poly_init(gamma);
  mpz_mat_init2(&I,d,d);
  mpz_init_set_ui(tmp1, 1); mpz_init(tmp2); mpz_init(tmp3);
  mpz_init(cdm); mpz_init(gamma_denom);


  lastReportTime = sTime() - 10.0;
  /****** Initialize M. ******/
  g_M.N = N; g_M.FB = FB;

  mpz_mod(tmp1, g_M.FB->n, g_M.FB->knownDiv);
  if (mpz_sgn(tmp1)==0)
    mpz_div(g_M.FB->n, g_M.FB->n, g_M.FB->knownDiv);
  mpz_set_ui(tmp1, 1);

  if ((retVal = initMsqrt(&g_M, relsInDep, prelF, lpF)))
    return retVal;
  mpz_mul(cdm, &g_M.N->f->coef[g_M.N->degree], g_M.FB->m);

  mpz_init_set_ui(res, 1);
  mpz_set_ui(rSqrt, 0);
  mpz_set_ui(aSqrt, 0);

  printf("Using monic polynomial T="); mpz_poly_print(stdout, "", g_M.N->T);
  printf("Integral basis is given by:\n");
  mpz_mat_print(stdout, N->W);
  printf("with denominator: "); mpz_out_str(stdout, 10, N->W_d);
  printf("\n");

  cont=1;
  l=0;
  finalTries=0;
  lastDiff=0.0;
  normA = g_M.logNormGNum - g_M.logNormGDen;
  normN = 0.0;
  for (i=0; i<N->degree; i++) 
    normN += g_M.log_eps[i];
  lastNormA = normA; lastNormN = normN;
  lastNormNum = g_M.logNormGNum; lastNormDen = g_M.logNormGDen;

  while (cont) {

    /* Simplify numerator or denominator? */
    if (fabs(g_M.logNormGNum) > fabs(g_M.logNormGDen))
      sl = 1;
    else
      sl = -1;
    /* On the other hand, we cannot let the HNum, HDen matrices */
    /* grow out of control:                                     */
    lastHNormNum = logNorm(&g_M.Hnum, &g_M);
    lastHNormDen = logNorm(&g_M.Hden, &g_M);
    if ((lastHNormNum > 0.9*g_M.LLL_max_log) && (g_M.logNormGNum > 1.0))
      sl = 1;
    if ((lastHNormDen > 0.9*g_M.LLL_max_log) && (g_M.logNormGDen > 1.0))
      sl = -1;
    l++;
    now = sTime();

    if ((now > (lastReportTime + 5.0)) || 
         ((g_M.logNormGDen < 1000.0) && (g_M.logNormGNum < 1000.0))) {
      printf("step %" PRId32 ", sl=%2d, logGam=%1.2lf/%1.2lf, ", 
              l, sl, g_M.logNormGNum, g_M.logNormGDen);
      printf("emb: ");
      for (i=0; i<N->degree; i++) 
        printf("%1.2lf ", g_M.log_eps[i]);
      normA = g_M.logNormGNum - g_M.logNormGDen;
      normN = 0.0;
      for (i=0; i<N->degree; i++) 
        normN += g_M.log_eps[i];
      diff = normA - normN;
      printf(", diff=%1.4lf\n", diff);
      lastReportTime = now;
    }

    chooseIdeal(&I, &g_M, sl);
    chooseDelta(delta, &I, sl, &g_M);
#ifdef _LOUD_DEBUG
    ofp = fopen("maple.txt", "a");
    fprintf(ofp, "*(");
    mpz_poly_print(ofp, "", delta);
    fprintf(ofp, ")^(%d)\n", sl);
    fclose(ofp);
#endif
    updateFactorization(&g_M, delta, &I, sl);
    mpz_set_ui(tmp1, 1);
    updateCRT(&g_M, delta, tmp1, -2*sl);

    evalOmegaPoly(tmp2, delta, cdm, g_M.FB->n, &g_M);
    if (sl == -1) {
      if (mpz_invert(tmp2, tmp2, g_M.FB->n)==0) {
        evalOmegaPoly(tmp2, delta, cdm, g_M.FB->n, &g_M);
        printf("Error: "); mpz_out_str(stdout, 10, tmp2);
        printf(" is not invertible!\n");
        printf("gcd with n is : ");
        mpz_gcd(tmp2, tmp2, g_M.FB->n);
        mpz_out_str(stdout, 10, tmp2); printf("\n");
        exit(-1);
      }
    }
    mpz_mul(res, res, tmp2);
    mpz_mod(res, res, g_M.FB->n);
    newTotal = fabs(g_M.logNormGDen) + fabs(g_M.logNormGNum) + fabs(logNorm(&g_M.Hnum, &g_M))
               + fabs(logNorm(&g_M.Hden, &g_M));
    totalSize[l%20]=newTotal;
    if (l>20) {
      improve = totalSize[(l+1)%20]-totalSize[l%20];
      if (improve < g_M.LLL_max_log) {
        printf("Warning: Small improvement (%1.4lf) over 20 iterations!\n", improve);
        g_M.LLL_max_log *= 0.9;
      } else {
        g_M.LLL_max_log = DEFAULT_LLLMAX_LOG;
      }
    }

    cont= ((g_M.logNormGDen>=0.5) || (g_M.logNormGNum > 0.5) ||
          !(mpz_mat_isID(&g_M.Hnum)) || !(mpz_mat_isID(&g_M.Hden)));
    /* If we're close, or there was too little improvement last time,
       increment the counter which will trigger an update of LLL_MAX.
    */
    if (fabs(g_M.logNormGDen) + fabs(g_M.logNormGNum) < 20.0*g_M.LLL_max_log) {
      finalTries++;
    } 
    /* If we get caught in a loop of the same leftover ideals,
       try to sneak out by mixing it up a bit.
    */
    if (finalTries > 50) {
      g_M.LLL_max_log = DEFAULT_LLLMAX_LOG/2;
    }
    if (finalTries > 100) {
      g_M.LLL_max_log = DEFAULT_LLLMAX_LOG/4;
    }
    if (finalTries > 150) {
      g_M.LLL_max_log = DEFAULT_LLLMAX_LOG/10;
    }
    if (((finalTries>250)&&(g_M.logNormGDen < 0.01)) || (finalTries > 270)) 
      /* Excessive - it shouldn't need nearly this many. */
      cont=0;
#ifdef _DEBUG
    normA = g_M.logNormGNum - g_M.logNormGDen;
    normN = 0.0;
    for (i=0; i<N->degree; i++) 
      normN += g_M.log_eps[i];
    diff = normA - normN;
    if (fabs(diff - lastDiff) > 2.0) {
      printf("**** Possible error: diff=%1.8lf vs. lastDiff=%1.8lf!\n", diff,lastDiff);
      printf("  Step %ld.\n", l);
      printf("  last normA = %1.5lf,   last normN = %1.5lf\n", lastNormA, lastNormN);
      printf("   new normA = %1.5lf,    new normN = %1.5lf\n\n", normA, normN);
      printf("last normNum = %1.5lf, last normDen = %1.5lf\n", lastNormNum, lastNormDen);
      printf(" new normNum = %1.5lf,  new normDen = %1.5lf\n\n", g_M.logNormGNum, g_M.logNormGDen);
      printf("last normHNum= %1.5lf, last normHDen= %1.5lf\n", lastHNormNum, lastHNormDen);
      printf(" new normHNum= %1.5lf,  new normHDen= %1.5lf\n\n", logNorm(&g_M.Hnum, &g_M),
              logNorm(&g_M.Hden, &g_M));
      if (sl==1) printf("The numerator was being simplified.\n");
      else printf("The denominator was being simplified.\n");
      printf("delta = "); mpz_poly_print(stdout, "", delta); printf("\n");
      norm_ib(tmp1, delta, g_M.N);
      printf("norm(delta)="); mpz_out_str(stdout, 10, tmp1); 
      printf(", log=%1.5lf\n", _mpz_log(tmp1));
      printf("Ideal selection was as follows:\n%s\n", idealSelStr);

    } 
    lastDiff=diff;
    lastNormA = normA; 
    lastNormN = normN;
    lastHNormNum = logNorm(&g_M.Hnum, &g_M); 
    lastHNormDen = logNorm(&g_M.Hden, &g_M);
    lastNormNum = g_M.logNormGNum; 
    lastNormDen = g_M.logNormGDen;
#endif
  }
  printf("-------------------------------------------------\n");
  printf("Iterative portion of square root computation done.\n");
  printf("step %" PRId32 ", sl=%2d, logGam=%1.2lf/%1.2lf, ", 
          l, sl, g_M.logNormGNum, g_M.logNormGDen);
  printf("emb: ");
  for (i=0; i<N->degree; i++) 
    printf("%1.2lf ", g_M.log_eps[i]);
  normA = g_M.logNormGNum - g_M.logNormGDen;
  normN = 0.0;
  for (i=0; i<N->degree; i++) 
    normN += g_M.log_eps[i];
  printf(", diff=%1.4lf\n", normA - normN);

  printf("Final CRT residues:\n");
  for (i=0; i<g_M.ipbSize; i++) {
    mpz_out_str(stdout, 10, &g_M.q[i]);
    printf(" : ");
    mpz_poly_print(stdout, "", g_M.crtRes[i]);
  }
#ifdef _DEBUG
  printf("Remaining ideals:\n");
  empty=1;
  for (i=g_M.aSize-1; i>=0; i--) {
    if (g_M.aExp[i] != 0) {
      printf("(%" PRId32 ", %" PRId32 ")^%" PRId32 "\n", g_M.AFB[i].p, g_M.AFB[i].r, g_M.aExp[i]);
      empty=0;
    }
  }
  for (i=g_M.spSize-1; i>=0; i--) {
    if (g_M.spExp[i] != 0) {
      printf("(EIdeal %d)^%" PRId32 "\n", i, g_M.spExp[i]);
      empty=0;
    }
  }
  if (empty) 
    printf("None.\n\n");
#endif


  /* Now everything has been simplified enough that we should
     be able to lift what's left to Z_K and compute the square root.
  */
  crtLift(gamma, gamma_denom, &g_M);

  printf("Residues lifted to a remaining gamma = ");
  printf("(1/"); mpz_out_str(stdout, 10, gamma_denom); printf(") * (");
  mpz_poly_print(stdout, "", gamma); printf("\n");

  mpz_cdiv_qr(tmp2, tmp1, g_M.N->W_d, gamma_denom);
  if (mpz_sgn(tmp1)) {
    msgLog("", "Warning: W_d not divisible by gamma_denom!");
    printf("Warning: W_d not divisible by gamma_denom!\n");
    printf("This computation will most likely fail.\n");
  }
  mpz_mul(tmp2, tmp2, g_M.N->W_d);
  for (i=0; i<=gamma->degree; i++) {
    mpz_mul(&gamma->coef[i], &gamma->coef[i], tmp2);
  }
  if (finalSquareRoot(beta, gamma, &g_M)) {
    printf("*** Some serious error occurred computing remaining square root.\n");
    printf("*** This run will most likely fail!\n");
  }
  mpz_poly_print(stdout, "final beta =(1/k)*", beta);

  /*************************************************/
  /* Now, Beta is a polynomial in \hat{\alpha}, so */
  /* it needs to be evaluated at c_d*m.            */
  /*************************************************/
  /* First correct for the factor of W_d^2 we multiplied in. */
  mpz_invert(tmp2, g_M.N->W_d, g_M.FB->n);  
  for (i=0; i<=beta->degree; i++) {
    mpz_mul(&beta->coef[i], &beta->coef[i], tmp2);
    mpz_mod(&beta->coef[i], &beta->coef[i], g_M.FB->n);
  }
  mpz_poly_eval2(tmp2, beta, cdm, g_M.FB->y1);
  mpz_mul(res, res, tmp2);
  mpz_mod(res, res, g_M.FB->n);

  printf("Beta evaluated at m = "); mpz_out_str(stdout, 10, res); printf("\n");
  /* Now, we must correct for the extra factors of Y1 in the denominator
     from the rational side. This is a sloppy way to do it, but hey - it works.
  */
  printf("ABexponentSum = %ld\n", ABexponentSum);
  if (ABexponentSum < 0) {
    mpz_invert(tmp2, g_M.FB->y1, g_M.FB->n);
    mpz_powm_ui(tmp2, tmp2, -ABexponentSum/2, g_M.FB->n);
  } else {
    mpz_mod(tmp2, g_M.FB->y1, g_M.FB->n);
    mpz_powm_ui(tmp2, tmp2, ABexponentSum/2, g_M.FB->n);
  }
  mpz_mul(res, res, tmp2);
  printf("Y1 corrected Beta   = "); mpz_out_str(stdout, 10, res); printf("\n");


  mpz_set(aSqrt, res);
  mpz_mul(tmp2, res, res); mpz_mod(tmp2, tmp2, g_M.FB->n);
  printf("Beta^2 == "); mpz_out_str(stdout, 10, tmp2); printf("\n");
  printf("Rational square root: "); mpz_out_str(stdout, 10, g_M.ratSqrt); printf("\n");
  mpz_mul(tmp2, g_M.ratSqrt, g_M.ratSqrt); mpz_mod(tmp2, tmp2, g_M.FB->n);
  printf("    ^2 == "); mpz_out_str(stdout, 10, tmp2); printf("\n");
  mpz_set(rSqrt, g_M.ratSqrt);


  mpz_poly_clear(delta); mpz_poly_clear(beta); mpz_poly_clear(gamma);
  mpz_mat_clear(&I);
  mpz_clear(tmp1); mpz_clear(tmp2); mpz_clear(gamma_denom); mpz_clear(cdm);
  mpz_clear(tmp3);

  free(g_M.AFB);
  for (i=0; i<g_M.spSize; i++)
    clearIdeal(&g_M.sPrimes[i]);
  free(g_M.sPrimes);
  free(g_M.v_cd_sPrimes);
  for (i=0; i<MAX_DIST_FACTS; i++)
    mpz_clear(&g_M.Cd.p[i]);
  free(g_M.aExp); free(g_M.spExp);
  mpz_mat_clear(&g_M.Hnum); mpz_mat_clear(&g_M.Hden);
  for (i=0; i<g_M.ipbSize; i++) {
    mpz_clear(&g_M.q[i]);
    mpz_poly_clear(g_M.crtRes[i]);
  }
  mpz_mat_clear(&g_M.Beta);
  mpz_clear(g_M.kappa);
  for (i=0; i<d; i++)
    mpz_clear(&g_M.omegaEvalM[i]);
  return 0;  
}  

