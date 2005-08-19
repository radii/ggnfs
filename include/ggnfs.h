/**************************************************************/
/* ggnfs.h                                                    */
/* Copyright 2003, Chris Monico.                              */
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
*   along with Foobar; if not, write to the Free Software
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _GGNFS_H
#define _GGNFS_H


#if defined (__MINGW32__) || defined (MINGW32)
#define _MSC_VER 1300
#endif

#include <gmp.h>

#if defined (__cplusplus)
extern "C" {
#endif

typedef unsigned char uchar;

#include <stdarg.h>
#define __STDC_FORMAT_MACROS
#if defined (_MSC_VER) && !defined(__MINGW32__)
#include <stdio.h>
#include <basetsd.h>

#define int8_t  INT8
#define uint8_t UINT8

#define int16_t  INT16
#define uint16_t UINT16

#define int32_t INT32
#define uint32_t UINT32

#define int64_t INT64
#define uint64_t UINT64

#define PRId8 "I8d"
#define PRIi8 "I8i"
#define PRIo8 "I8o"
#define PRIu8 "I8u"
#define PRIx8 "I8x"
#define PRIX8 "I8X"

#define PRId16 "I16d"
#define PRIi16 "I16i"
#define PRIo16 "I16o"
#define PRIu16 "I16u"
#define PRIx16 "I16x"
#define PRIX16 "I16X"

#define PRId32 "I32d"
#define PRIi32 "I32i"
#define PRIo32 "I32o"
#define PRIu32 "I32u"
#define PRIx32 "I32x"
#define PRIX32 "I32X"

#define PRId64 "I64d"
#define PRIi64 "I64i"
#define PRIo64 "I64o"
#define PRIu64 "I64u"
#define PRIx64 "I64x"
#define PRIX64 "I64X"

#define SCNd32 "I32d"
#define SCNi32 "I32i"
#define SCNo32 "I32o"
#define SCNu32 "I32u"
#define SCNx32 "I32x"
#define SCNX32 "I32X"

#define SCNd64 "I64d"
#define SCNi64 "I64i"
#define SCNo64 "I64o"
#define SCNu64 "I64u"
#define SCNx64 "I64x"
#define SCNX64 "I64X"

#else
#include <inttypes.h>
#endif
#define _USE_MATH_DEFINES
#include <math.h>
#include "version.h"

/* This is just for a procrels hack. */
#define BAD_LP_INDEX 0xFFFFFFFE

#define GGNFS_LOG_NAME "ggnfs.log"

#define MAX_NAMESIZE 256
#define MAX_LARGE_RAT_PRIMES  3
#define MAX_LARGE_ALG_PRIMES  3
/* Overkill, but it we don't allocate many of the structures with these. */
#define MAX_RAT_FACTORS      300
#define MAX_ALG_FACTORS      300
#define MAX_SP_FACTORS       70

/* For the most part, the <a-b\alpha> shouldn't contain more than
   this many factors of any one special prime:
*/
#define MAX_SP_POWERS 50

#define MAX_ROWS_IN_COL 2048
#define MAX_RELS_IN_FF 48

/* Clients need not report the parts of
   the relation factorizations occurring in the
   first x of entries of the factor base(s).
*/
#define CLIENT_SKIP_R_PRIMES 50
#define CLIENT_SKIP_A_PRIMES 50

#define MAXFNAMESIZE  127

#define ENDRFB "#ENDRFB"
#define ENDAFB "#ENDAFB"
#define ENDQCB "#ENDQCB"

#define BF_DELIMITER 0x0a

#ifdef _MSC_VER 
#define INLINE 
#define inline __inline
#define vsnprintf _vsnprintf
#else
#define INLINE inline
#endif

/* This is the max number of pairs per call to some of
   the classical sieve and isSmooth (withInfo) functions.
*/
#define MAX_PAIRS_PER_CALL 50000

#define BIT(_i) ((0x00000001)<<(_i))
#define MAX(_a, _b) ((_a) > (_b) ? (_a) : (_b))
#define MIN(_a, _b) ((_a) < (_b) ? (_a) : (_b))

/* should be defined in inttypes.h.  If your system doesn't
   have this, you'll need to define these and the printf/scanf
   format specifiers (e.g. PRId32, SCNu64, etc) yourself */
#define s32 int32_t
#define s64 int64_t
#define u32 uint32_t
#define u64 uint64_t
#define u8  uint8_t
#define u16 uint16_t


/************************************************/
typedef struct {
  char prefix[60];
  int  numFiles;
  s32 data1, data2; /* Used for whatever. */
} multi_file_t;


#define MAXPOLYSIZE 25
typedef struct {
  s32 coef[2*MAXPOLYSIZE+1];
  int  degree;
} _poly_t;
typedef _poly_t poly_t[1];

#define MAXPOLYDEGREE 11
typedef struct {
  __mpz_struct coef[MAXPOLYDEGREE+1];
  int degree;
} mpz_poly_t;
typedef mpz_poly_t mpz_poly[1];

typedef struct {
  double r;
  double i;
  mpf_t  mpr;
  mpf_t  mpi;
} nfs_complex_t;

/*****************************************************************/
/* This structure is useful for certain purposes, like           */
/* computing the factorizations of (a-bm) and <a-b\alpha>.       */
/* Yet, when it comes time to complete a factorization using     */
/* a large number of relations, it won't do. That's what         */
/* the `rel_list' is for.                                        */
/*****************************************************************/
typedef struct {
//  s32 a;
  s64 a;
  s32 b;
  /* The fields below are filled n by getRFact() */
  s32  rFactors[MAX_RAT_FACTORS]; /* Rational factors dividing (a+bm). */
  char  rExps[MAX_RAT_FACTORS];    /* exponents.                        */
  int   rFSize;
  s32  aFactors[MAX_ALG_FACTORS]; /* Algebraic factors dividing <a+b\alpha>. */
  char  aExps[MAX_ALG_FACTORS];    /* exponents.                              */
  int   aFSize;
  int   spFactors[MAX_SP_FACTORS];   /* ``special primes'' dividing <a+b\alpha>. */
  int   spExps[MAX_SP_FACTORS];      /* exponents. */
  int   spSize;
  s32  qcbBits[2]; /* bit-packed version of the QCB columns.  */
  s32 p[MAX_LARGE_RAT_PRIMES];   /* Large rational primes */
  s32 a_p[MAX_LARGE_ALG_PRIMES]; /* Large algebraic primes: 'p' part. */
  s32 a_r[MAX_LARGE_ALG_PRIMES]; /* Large algebraic primes: 'q' part. */
} relation_t;


typedef struct {
  s32  numRels, maxRels, maxDataSize;
  s32 *relIndex; /* relIndex[i] is the index of the start of relation 'i' in
                     the following array: */
  s32 *relData;  /* This is where all the relation data is held, in a very
                     specific format: see rels.c for the format description.
                   */
} rel_list;

/* This is a more generic structure, (which will be) used
   for building full relations and matrix operations.
*/
typedef struct {
  s32 *data;
  s32 *index;
  s32  maxDataSize, maxFields;
  s32  numFields;
} llist_t;


typedef struct {
  s32 x, y;
} lpair_t;



typedef struct {
  char       name[MAX_NAMESIZE];
  mpz_poly   f;
  double     skew;
  mpz_poly   g; /* g(x) := c_d^{d-1}f(x/c_d). */
  nfs_complex_t zeros[MAXPOLYDEGREE];
  mpz_t      n;
  mpz_t      m;
  mpz_t      y0;
  mpz_t      y1;
  mpz_t      disc;     /* not set unless it's needed. */
  mpz_t      cd;
  mpz_t      knownDiv;
  s32       rfb_size;   /* rational factor base size  */
  s32       rLim;       /* rational factor base limit. */
  s32      *rfb;        /* rational factor base: rfb[2k]=p_k, rfb[2k+1] = m (mod p_k). */
  double     rfb_log_base, log_rlb;
  int       *rfb_log;    /* Approximations of logs to the above base. */
  double     rfb_lambda; /* Franke's lambda; I'll explain later. */

  double     afb_log_base, log_alb;
  s32       afb_size;   /* algebraic factor base size  */
  double     afb_lambda; /* Franke's lambda; I'll explain later. */
  s32      *afb;        /* algebraic factor base: afb[2k]=p_k, afb[2k+1]=r_k (root of 'f' mod p_k). */
  s32       aLim;       /* algebraic factor base limit. */
  int       *afb_log;    /* fixed-point log approximations */
  int        afb_log_ff; /* fudge-factor for sieving.   */
  s32       qcb_size;   /* quadratic character base size */
  s32      *qcb;        /* quadratic character base. Formatted same as afb. */
  s32       maxP_r;     /* Max large (leftover) rational prime */
  s32       maxP_a;     /* Max large (leftover) algebraic prime size */
  int        maxLP;      /* Max # of large rational primes. */
  int        maxLPA;     /* Max # of large algebraic primes. */
  int        MFB_r;      /* Bits in largest value to factor into large primes */
  int        MFB_a;      /* Bits in largest value to factor into large primes */
} nfs_fb_t;


typedef struct {
  nfs_fb_t *FB;
  relation_t *R;
  s32 numR;
  s32 **relsInFF;
  s32 numFF;
  s32 minFF;
  s32 *colEntries;
  s32 *colIndex;
  s32 n;
  s32 m;
} nfs_mat_t;

typedef struct {
  mpz_t         N;
  mpz_t         unfactored;
  __mpz_struct *factors;
  int          *exponents;
  int           size;
  int           sign;
} mpz_fact_t;


#define SIEVE_UNKNOWN   0
#define SIEVE_CLASSICAL 1
#define SIEVE_LATTICE   2

typedef struct {
  int     type;
  double  rfb_lambda;
  double  afb_lambda;
  char    outName[MAXFNAMESIZE+1];
  /* Classical specific: */
  s32 a0;
  s32 a1;
  s32 b0;
  s32 b1;
  /* Lattic specific: */
  s32 qIndex0;
  s32 qIndex1;
  double sieveArea;
  s32 cacheSize;
  nfs_fb_t FB;
  char     fbName[MAXFNAMESIZE+1];
  int      nDigits;
} nfs_sieve_job_t;

  
typedef struct {
  s32 p;
  s32 r;
  s32 e;
  /* p==r represents (p,\infty). */
} factor_t;


#define DEFAULT_MPZ_MAT_ROWS (8*MAXPOLYDEGREE+2)
#define DEFAULT_MPZ_MAT_COLS (8*MAXPOLYDEGREE+2)
typedef struct {
  __mpz_struct **entry;
  int rows, cols;
  int maxRows, maxCols;
} mpz_mat_t;

typedef struct {
  mpz_t    p;
  mpz_poly alpha; /* The prime ideal is <p, alpha> */
  mpz_poly beta;  /* structure constant for computing valuations. */
  mpz_mat_t betaMat; /* cols are beta, beta*alpha,..., beta*alpha^{d-1}. */
  int      e;     /* valuation of <p> at this ideal. */
  int      f;     /* This ideal has norm p^f. */
} prime_id_t;


typedef struct {
  mpz_poly      f;      /* NFS polynomial, f=c_0 + c_1x + ... + c_dx^d. */
  mpz_poly      T;      /* Defining monic polynomial: T(x) = c_d^{d-1}f(x/c_d). */
  int           degree; /* Degree of the field extension over Q. */
  mpz_mat_t     *W;     /* Z-Basis for the ring of integers, in HNF. */
  mpz_poly      *Mt;    /* Multiplication table for the omega_i. */
  mpz_t         W_d;    /* Denominator for the above. */
  mpz_t         Kdisc;  /* Discriminant of the field. */
  mpz_t         Tdisc;  /* Discriminant of T.         */
  mpz_t         index;  /* index = sqrt(disc(T)/Kdisc). */
  mpz_mat_t     *W_inv; /* W_inv*W/d = Identity.        */
  nfs_complex_t fZeros[MAXPOLYDEGREE]; /* Zeros of f. */
  nfs_complex_t TZeros[MAXPOLYDEGREE]; /* Zeros of T. */
  double        WevalZ_r[MAXPOLYDEGREE][MAXPOLYDEGREE], WevalZ_i[MAXPOLYDEGREE][MAXPOLYDEGREE];
                /* The elements of the integral basis evaluated at the zeros. That is,
                   WevalZ[i][j] is \omega_i evaluated at \alpha = \sigma_j(\alpha).
                */
  prime_id_t   *sPrimes;       /* Primes dividing the index [Z_k : Z[alpha]]. */
  mpz_mat_t     sPowMats[MAX_SP_FACTORS][MAX_SP_POWERS];
  int          *v_cd_sPrimes;  /* Valuations of c_d at the sPrimes. */
  int           numSPrimes;
  mpz_poly      Sk_ib;         /* Sk_ib[j] = Trace(omega_j). */
  nfs_fb_t     *FB;            /* Not always needed, but handy to have around. */

} nf_t;

#define MAX_DENSE_BLOCKS 64
typedef struct {
  s32 numRows, numCols, maxDataSize;
  s32 *cEntry;
  s32 *cIndex;
  u64 *denseBlocks[MAX_DENSE_BLOCKS];
  s32  denseBlockIndex[MAX_DENSE_BLOCKS];
  s32  numDenseBlocks;
} nfs_sparse_mat_t;




/* matstuff.c */
int getDependencies(nfs_sparse_mat_t *M, llist_t *C, s32 *deps);
int writeSparseMat(char *fname, nfs_sparse_mat_t *M);
int readSparseMat(nfs_sparse_mat_t *M, char *fname);
int cm_init(llist_t *C, nfs_sparse_mat_t *M);
int cm_removeCols(llist_t *C, s32 *cols, s32 numCols);
int removeCols(nfs_sparse_mat_t *M, llist_t *C, s32 *cols, s32 m);
int pruneMatrix(nfs_sparse_mat_t *M, s32 minExtraCols, double wtFactor, llist_t *C);


/* mpz_poly.c */
void   mpz_poly_init(mpz_poly f);
void   mpz_poly_clear(mpz_poly f);
void   mpz_poly_cp(mpz_poly dest, mpz_poly src);
int    mpz_poly_fixDeg(mpz_poly f);
int    mpz_poly_cmp(mpz_poly A, mpz_poly B);
int    mpz_poly_diff(mpz_poly res, mpz_poly q);
int    mpz_poly_getComplexZeros(nfs_complex_t *Z, mpz_poly f);
void   mpz_poly_print(FILE *fp, char *str, mpz_poly f);
int    mpz_poly_getBaseM(mpz_poly f, mpz_t n, mpz_t m, int init_poly);
int    mpz_poly_getBaseM2(mpz_poly f, mpz_t n, mpz_t m, mpz_t modulus, mpz_t residue);
int    mpz_poly_discrim(mpz_t disc, mpz_poly A);
void   mpz_poly_resultant(mpz_t res, mpz_poly A, mpz_poly B);
void   mpz_poly_psuedoDiv(mpz_poly Q, mpz_poly R, mpz_poly A, mpz_poly B);
void   mpz_poly_eval(mpz_t res, mpz_poly f, mpz_t x);
void   mpz_poly_eval2(mpz_t res, mpz_poly f, mpz_t x, mpz_t y1);
void   mpz_poly_evalD(mpz_t res, mpz_poly f, mpz_t x);
int    mpz_poly_HenselLift(mpz_t res, mpz_poly f, s32 a, s32 p, int k);
int    mpz_poly_pow_mod(mpz_poly res, mpz_poly f, mpz_t n, mpz_poly mod);
int    mpz_poly_pow_mod_pp(mpz_poly res, mpz_poly f, mpz_t n, mpz_poly mod, mpz_t p);
int    mpz_poly_mod(mpz_poly r, mpz_poly a, mpz_poly b);
int    mpz_poly_mod_pp(mpz_poly r, mpz_poly a, mpz_poly b, mpz_t p);
int    mpz_poly_gcd(mpz_poly g, mpz_poly h, mpz_t p);
int    mpz_poly_modp(poly_t res, mpz_poly q, s32 p);
int    mpz_polyCoprime(mpz_poly f, mpz_poly g, mpz_t p);
s32   mpz_poly_evalModp(mpz_poly h, s32 p, s32 x);
void   mpz_poly_mul(mpz_poly res, mpz_poly op1, mpz_poly op2);
void   mpz_poly_add(mpz_poly res, mpz_poly op1, mpz_poly op2);
void   mpz_poly_mulmod(mpz_poly res, mpz_poly op1, mpz_poly op2, mpz_poly mod);
void   mpz_poly_mulmod_p(mpz_poly res, mpz_poly op1, mpz_poly op2, mpz_t p);
void   mpz_poly_mulmod_pp(mpz_poly res, mpz_poly op1, mpz_poly op2, mpz_poly mod, mpz_t p);
int    mpz_poly_irreducible_modp(mpz_poly _f, mpz_t p);
int    mpz_poly_irreduciblelike_modp(mpz_poly _f, mpz_t p);
int    mpz_poly_div(mpz_poly q, mpz_poly x, mpz_poly y, mpz_t p);  
int    mpz_poly_fact(mpz_poly *Pi, int *exps, mpz_poly A, mpz_t p);
int    mpz_poly_inv(mpz_poly inv, mpz_poly f, mpz_poly m, mpz_t p);


int    Zalpha_sqrt(mpz_poly res, mpz_poly a, mpz_poly _f, mpz_t N, mpz_t p);

/* poly.c */
void poly_cp(poly_t dest, poly_t src);
int  poly_fixDeg(poly_t op);
int  poly_Jacobi(poly_t f, poly_t g, s32 p);
int  poly_mul_modp(poly_t res, poly_t op1, poly_t op2, s32 p);
int  poly_mulmodpp(poly_t res, poly_t op1, poly_t op2, poly_t mod, s32 p);
int  poly_mod(poly_t res, poly_t op, poly_t _mod, s32 p);
int  poly_xpow_modpp(poly_t res, s32 n, poly_t mod, s32 p);
int  poly_pow_modpp(poly_t res, poly_t f, s32 n, poly_t mod, s32 p);
int  poly_powmpz_mod(poly_t res, poly_t f, mpz_t n, poly_t mod, s32 p);
int  poly_gcd(poly_t g, poly_t h, s32 p);
int  poly_div(poly_t _q, poly_t _r, poly_t x, poly_t y, s32 p);
int  poly_inv(poly_t inv, poly_t f, poly_t m, s32 p);
int  poly_getZeros(s32 *zeros, poly_t _f, s32 p);
int  poly_isSimpleZero(poly_t f, s32 p, s32 r);
int  poly_coprime(poly_t f, poly_t g, s32 p);
int  poly_irreducible_modp(poly_t _f, s32 p);


/* mpz_mat.c */
void mpz_mat_init(mpz_mat_t *M);
void mpz_mat_init2(mpz_mat_t *M, int maxRows, int maxCols);
void mpz_mat_clear(mpz_mat_t *M);
void getCol(mpz_poly f, mpz_mat_t *M, int col);
void setCol(mpz_mat_t *M, int col, mpz_poly f);
void mpz_mat_print(FILE *fp, mpz_mat_t *M);
int  mpz_mat_equal(mpz_mat_t *A, mpz_mat_t *B);
int  mpz_mat_iszero(mpz_mat_t *M);
int  mpz_mat_isID(mpz_mat_t *M);
void mpz_mat_setID(mpz_mat_t *M, int n);
void mpz_mat_dotCols(mpz_t res, mpz_mat_t *M, int col1, int col2);
int  mpz_mat_LLL(mpz_mat_t *B, mpz_mat_t *H, mpz_mat_t *A);
int  mpz_mat_LatticeIntersection(mpz_mat_t *I, mpz_mat_t *B_i, int s);
void mpz_mat_mul(mpz_mat_t *Res, mpz_mat_t *X, mpz_mat_t *Y);
void mpz_mat_cp(mpz_mat_t *dest, mpz_mat_t *src);
void mpz_mat_cat(mpz_mat_t *dest, mpz_mat_t *src1, mpz_mat_t *src2);
int  mpz_mat_getHNF(mpz_mat_t *W, mpz_mat_t *A);
int  mpz_mat_getHNF_mod_D(mpz_mat_t *W, mpz_mat_t *A, mpz_t D);
int  mpz_mat_getKernel(mpz_mat_t *I, mpz_mat_t *B);
int  mpz_mat_getKernel_modp(mpz_mat_t *X, mpz_mat_t *B, mpz_t p);
int  mpz_mat_getImage_modp(mpz_mat_t *Y, mpz_mat_t *N, mpz_t p);
int  mpz_mat_invIm(mpz_mat_t *X, mpz_mat_t *V, mpz_mat_t *_M, mpz_t p);
int  mpz_mat_suppBasis_modp(mpz_mat_t *B, mpz_mat_t *_M, mpz_t p);
int  mpz_mat_suppSubspace_modp(mpz_mat_t *Z, mpz_mat_t *V, mpz_mat_t *M, mpz_t p);
int  mpz_mat_moduleAdd(mpz_mat_t *H, mpz_t dH,
                       mpz_mat_t *W, mpz_t d, mpz_mat_t *W_, mpz_t d_);
void mpz_mat_pseudoInvert(mpz_mat_t *I, mpz_t d, mpz_mat_t *J);

int mpz_mat_getHNF_sp(mpz_mat_t *W, mpz_mat_t *A);





/* misc.c */
double sTime();
double _mpz_log(mpz_t k);
int    readBinField(char *str, int size, FILE *fp);
int    writeBinField(FILE *fp, char *str);
s32   powMod(s32 op, s32 n, s32 p);
s32   inverseModP(s32 n, s32 p);
s32   gcd(s32 x, s32 y);
int    sqrtModP(mpz_t res, mpz_t x2, mpz_t p);

//int    mpz_evalF(mpz_t res, s32 a, s32 b, mpz_poly f);
int    mpz_evalF(mpz_t res, s64 a, s32 b, mpz_poly f);

#ifdef LONG64
#define mpz_mul_si64(mp1, mp2, x) mpz_mul_si(mp1, mp2, x)
#define mpz_set_si64(mp, x) mpz_set_si(mp, x)
#else
void   mpz_mul_si64( mpz_t rop, mpz_t op1, s64 a);
void   mpz_set_si64( mpz_t rop, s64 a);
#endif

double mpz_evalF_d(double x, double y, mpz_poly f);
int    fplog_evalF(s32 a, s32 b, nfs_fb_t *FB);
int    fplog_mpz(mpz_t k, double log_of_base);
int    fplog(s32 k, double log_of_base);
double log_evalPoly(double x, double y, mpz_poly f);
void   mpz_nearest_int(mpz_t q, mpz_t num, mpz_t den);
void   msgLog(char *fName, char *fmt, ...);
void   printTmp(char *fmt, ...);
double L_n(double n, double a, double c);
s32    removeS32Pairs(s32 *L, s32 size);
int    mallocReport();
void  *lxmalloc(size_t n, int fatal);
void  *lxcalloc(size_t n, int fatal);
void  *lxrealloc(void *x, size_t n, int fatal);
int    getHeapStats(int *maxUseage, int *errs, int *allocs, int *reallocs);
void   logHeapStats();


       /* These are used for factoring the discriminant. */
void   mpz_fact_init(mpz_fact_t *F);
void   mpz_fact_clear(mpz_fact_t *F);
void   mpz_fact_add_factor(mpz_fact_t *F, mpz_t factor, int exponent);
int    mpz_fact_check(mpz_fact_t *D, int doRealWork);
int    mpz_fact_removeSF(mpz_fact_t *S, mpz_fact_t *F);
int    mpz_fact_factorEasy(mpz_fact_t *F, mpz_t N, int doRealWork);
void   initNF(nf_t *N);
void   clearNF(nf_t *N);
void   initIdeal(prime_id_t *I);
void   clearIdeal(prime_id_t *I);
void   reorderRoots(nfs_complex_t *z, int n);
int    cmpS32s(const void *a, const void *b);
int    cmpU32s(const void *a, const void *b);
int    cmp_lpair_t(const void *a, const void *b);


/* makefb.c */
int createFB(nfs_fb_t *FB, char *ofname);

	
	
/* muldmod32.s */
#ifdef __ppc__
#include <stdint.h>
static inline s32 mulmod32(uint32_t x, uint32_t y, uint32_t m)
{
   return ((uint64_t)x*(uint64_t)y%m);
}
#elif defined( _MSC_VER )
#else
extern s32 mulmod32(s32 op1, s32 op2, s32 modulus) asm("mulmod32");
#endif



/* clsieve.c */
static int forceStop;
int clSieve(nfs_sieve_job_t *J);

/* squfof.c */
u32 squfof(mpz_t n);

/* latsieve.c */
int latSieve(s32 *candidates, int maxCands, nfs_fb_t *FB,
             s32 qIndex, double sieveArea, s32 RAMtoUse);

        
/* smintfact.c */
int   smallIntFactor(u32 *factors, int *numFactors, mpz_t _t);
int   factor_using_pollard_rho (u32 *factors, int *numFactors, mpz_t n, int a_int);
int   factor(s32 *factors, mpz_t n, int useTrialDivision);


/* fbmisc.c */
void   initFB(nfs_fb_t *FB);
int    writePoly(FILE *fp, nfs_fb_t *FB);
int    readPoly(FILE *fp, nfs_fb_t *FB);
int    saveFB(char *fName, nfs_fb_t *FB);
int    loadFB(char *fName, nfs_fb_t *FB);
int    setLogs(nfs_fb_t *FB, s32 a0, s32 a1, s32 b0, s32 b1);
int    isSmooth_rat(s32 a, s32 b, nfs_fb_t *FB);
int    isSmooth_alg(s32 a, s32 b, nfs_fb_t *FB);
int    isSmooth_alg_q(s32 a, s32 b, s32 qIndex, nfs_fb_t *FB);
int    isSmooth_rat_withInfo(relation_t *R, nfs_fb_t *FB);
int    isSmooth_rat_withInfo_par(relation_t *R, int numRels, nfs_fb_t *FB);
int    isSmooth_alg_withInfo(relation_t *R, nfs_fb_t *FB);
int    isSmooth_alg_withInfo_par(relation_t *R, int numRels, nfs_fb_t *FB);

int    generateAFB(nfs_fb_t *FB, int verbose);
int    generateQCB(nfs_fb_t *FB, int size );
double getIPrimes(s32 *res, s32 numP, s32 upperBound, nfs_fb_t *FB);
s32  *getIPB_old(s32 *size, double lambda, nfs_fb_t *FB);
int    get_g(mpz_poly g, nfs_fb_t *FB);


/* ecm4c.c */
int ecmFactor(mpz_t p, mpz_t n, s32 _B1, double B2, int iter, mpz_t s);

/* getprimes.c */
u32   pSieve(u32 *p, u32 Psize, u32 a, u32 b);
u32   *getPList(u32 *numPrimes);
u32   getNextPrime(u32 n);
u32   getPrevPrime(u32 n);
u32   getMaxP(u32 a, u32 b);
u32   approxPi_x(u32 b);

/* llist.c */
int  ll_init(llist_t *L, s32 maxDataSize, s32 maxFields);
void ll_clear(llist_t *L);
int  ll_verify(llist_t *L);
void ll_resize(llist_t *L, s32 maxDataSize);
int  ll_appendField(llist_t *L, s32 *entries, s32 numEntries);
int  ll_deleteFields(llist_t *L, s32 *fields, s32 numFields);
int  ll_catFields(llist_t *L, lpair_t *pairs, s32 numPairs, int mod2);
int  ll_getsortOnFieldSize(s32 *fields, llist_t *L);
int  ll_numCommonEntries(llist_t *L, s32 field0, s32 field1);
int  ll_write(char *fname, llist_t *C);
int  ll_read(llist_t *C, char *fname);


/* combparts.c */
s32 combParts(llist_t *R, llist_t *P, int maxRelsInFF, s32 minFF);




/* rels.c */
#define GETNUMRFB(_s) ((_s)>>24)
#define GETNUMAFB(_s) (((_s)&0x00FF0000)>>16)
#define GETNUMSPB(_s)  (((_s)&0x0000FF00)>>8)
#define GETNUMLRP(_s)  (((_s)&0x000000C0)>>6)
#define GETNUMLAP(_s)  (((_s)&0x00000030)>>4)
                                                                                                            
#define SETNUMRFB(_s,_n) (_s = (_s&0x00FFFFFF)^((((s32)(_n)&0x000000FF)<<24)))
#define SETNUMAFB(_s,_n) (_s = (_s&0xFF00FFFF)^((((s32)(_n)&0x000000FF)<<16)))
#define SETNUMSPB(_s,_n) (_s = (_s&0xFFFF00FF)^((((s32)(_n)&0x000000FF)<<8)))
#define SETNUMLRP(_s,_n) (_s = (_s&0xFFFFFF3F)^((((s32)(_n)&0x00000003)<<6)))
#define SETNUMLAP(_s,_n) (_s = (_s&0xFFFFFFCF)^((((s32)(_n)&0x00000003)<<4)))
#define S32S_IN_ENTRY(_s) (2*GETNUMRFB(_s)+2*GETNUMAFB(_s)+2*GETNUMSPB(_s)+GETNUMLRP(_s)+2*GETNUMLAP(_s) + 6)
/* In the above, the +6 is to count also the: (1 size field) + (a,b fields) + (2 qcb fields) */
int   factRel(relation_t *R, nf_t *N);
int   factRelQ(relation_t *R, nf_t *N, s32 FBIndex);
int   factRels_clsieved(relation_t *R, int numRels, nf_t *N);
int   completePartialRelFact(relation_t *R, nf_t *N, s32 rTDiv, s32 aTDiv);
int   relConvertToData(s32 *data, relation_t *R);
int   dataConvertToRel(relation_t *R, s32 *data);
void  readRawS32(s32 *w, FILE *fp);
void  writeRawS32(FILE *fp, s32 *w);
int   writeRel(FILE *fp, s32 *relData);
int   writeRelList(char *fname, rel_list *L);
int   readRel(relation_t *R, FILE *fp);
int   readRelList(rel_list *L, char *fname);
int   parseOutputLine(relation_t *R, char *str, nfs_fb_t *FB);
void  makeOutputLine(char *str, relation_t *R, nfs_fb_t *FB);
s32  lookupRFB(s32 p, nfs_fb_t *FB);
s32  lookupAFB(s32 p, s32 r, nfs_fb_t *FB);




/* blanczos.c */
typedef void (* MAT_MULT_FUNC_PTR32)(u32 *Product, u32 *x, void *P);
int    blockLanczos32(u32 *deps, MAT_MULT_FUNC_PTR32 LeftMul, 
                      MAT_MULT_FUNC_PTR32 RightMul, void *P, s32 n);
/* I will do away with these prototypes later. */
void MultB32(u32 *Product, u32 *x, void *P);
void MultB_T32(u32 *Product, u32 *x, void *P);
#define BIT64(_i) ((0x0000000000000001ULL)<<(_i))
#define NUM_QCB_BLOCKS 1


typedef void (* MAT_MULT_FUNC_PTR64)(u64 *Product, u64 *x, void *P);
int    blockLanczos64(u64 *deps, MAT_MULT_FUNC_PTR64 LeftMul, 
                      MAT_MULT_FUNC_PTR64 RightMul, void *P, s32 n);
void   seedBlockLanczos(s32 seed);
/* I will do away with these prototypes later. */
void MultB64(u64 *Product, u64 *x, void *P);
void MultB_T64(u64 *Product, u64 *x, void *P);

        
/* nfmisc.c */
/* These aren't quite named consistently yet, I think. I want functions
   that expect or return things w.r.t. the integral basis to end with "_ib".
*/
int    getIntegralBasis(nf_t *N, mpz_fact_t *D, int tryHard);
void   basisPolyMult(mpz_poly res, mpz_poly op1, mpz_poly op2, mpz_poly *Mt, int n);
int    factorPrime(prime_id_t *I, mpz_t p, nf_t *N);
int    computeVconst(mpz_poly Beta, mpz_t p, mpz_poly alpha, nf_t *N);
int    valuation(mpz_mat_t *A, prime_id_t *I, nf_t *N);
int    valuation2(mpz_mat_t *A, prime_id_t *I, nf_t *N);
int    valuation_ab(mpz_t a, mpz_t b, int k, nf_t *N);
void   idealHNFMul_ib(mpz_mat_t *A, mpz_mat_t *B, mpz_mat_t *C, nf_t *N);
void   idealHNFPow_ib(mpz_mat_t *A, mpz_mat_t *B, s32 e, nf_t *N);
void   getIdealHNF_ib(mpz_mat_t *H, mpz_t p, mpz_poly a, s32 e, nf_t *N);
void   norm_std(mpz_t n, mpz_poly a, nf_t *N);
void   norm_ib(mpz_t n, mpz_poly a, nf_t *N);
void   idealHNF_ib(mpz_mat_t *H, mpz_poly a, nf_t *N);
void   idealHNF_ib_ab(mpz_mat_t *H, s64 a, s32 b, nf_t *N);
int    stdtoib(mpz_poly res, mpz_poly op, nf_t *N);
void   getTrace_ib(mpz_t t, mpz_poly a, nf_t *N);


/* fbgen.c */
u32 root_finder(u32 * root_buf, mpz_t * A, u32 adeg, u32 p);

int    factorN(mpz_t p, mpz_t q, s32 *dep, relation_t *R, nfs_fb_t *FB, nf_t *N);
int    montgomerySqrt(mpz_t rSqrt, mpz_t aSqrt, s32 *relsInDep, multi_file_t *prelF,
                      multi_file_t *lpF, nfs_fb_t *FB, nf_t *N);
double dickman(double x);
double dickmanStrong(double x, int numTerms);

/* assess.c */
void init_assess(double b0, double b1, double area, int pb);
unsigned int invert(unsigned int a, unsigned int p);  /* 0<b<p */
void murphy_en(double *me, int deg0, double *dbl_coeff0, int deg1, double *dbl_coeff1,
               double alpha0, double alpha1, double skewness, int nsm);
void murphy_e(double *me, int deg0, double *dbl_coeff0, int deg1, double *dbl_coeff1,
              double alpha0, double alpha1, double skewness);
int compute_alpha(double *alpha, int deg, unsigned int **coeffmod, mpz_t *gmp_coeff, double alpha_targ);
void compute_alpha_exact(double *alpha, int deg, unsigned int **coeffmod, mpz_t *gmp_coeff, unsigned int pb);
  	       
/* roots.c */
int find_optima(int *deg, double **coeff, double skewness, double **optima);

/* primes.c */
void prime_table_init();
unsigned int get_next_prime();


/*******************************************************/
/* We need a macro for hashing (s32, s32) --> Z/kZ.    */
/* It doesn't have to be a very good hash - so long as */
/* it's reasonably uniform, and doesn't change during  */
/* the course of a particular factorization.           */
/* It is used to decide which file we should put a     */
/* factored relation into. i.e., we don't necessarily  */
/* want them all in one file, or we may exceed an      */
/* OS-imposed filesize limit. On the other hand, we    */
/* don't want to randomly divy them up, or it will be  */
/* harder to look for collisions.                      */
/*******************************************************/

#define NFS_HASH(_a, _b, _k) (((u32)(2718281*(_a))^(3141592*(_b)^((_a)>>2)\
                              ^((_b)>>3)))%(_k))

#define NFS_HASH_Q(_a, _b, _k) (((u32)(3*(_a/2)+(_b)))%(_k))

#if defined(_MSC_VER) || defined(__ppc__)
#define MULMOD32(_res, _op1, _op2, _mod) _res = mulmod32(_op1, _op2, _mod)
#else
#define MULMOD32(_res, _op1, _op2, _mod) \
          asm(" imull %2\n\t idivl %3\n\t movl %%edx, %0"  \
               : "=mq" (_res)  \
               : "a" (_op1) , "mb" (_op2) , "mc" (_mod) : "%edx");     
#endif

#if defined (__cplusplus)
}
#endif

#endif
