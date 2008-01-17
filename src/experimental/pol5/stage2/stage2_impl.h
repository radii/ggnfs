#ifndef _STAGE2_IMPL_H_
#define _STAGE2_IMPL_H_

#include <stdio.h>
#include <string.h>
#include "ggnfs.h"
#include "if.h"
#include "profile.h"
#include "stage2_inline.h"
#include "poly_stage2.h"

#ifdef __cplusplus
extern "C" {
#endif

#define  MAX_PRIME_PROJ       100
#define  MAX_PRIME_AFF        200
#define  NPROJ_PRIMES          25
#define  NAFF_PRIMES           46
#define  MAX_X       1000000	/* !!! */
#define  MAX_Y           100	/* !!! */
#define  SIEVELEN           8192

/*-----------------------------------------------------------------------*/
/* profiling and statistics */

typedef struct {
	profile_t profile;
} stage2_stat_t;

#define DO_PROFILE

#ifdef DO_PROFILE
	#define profile_start(i) PROFILE_START(&stats->profile,i)
	#define profile_stop(i) PROFILE_STOP(&stats->profile,i)
#else
	#define profile_start(i) /* nothing */
	#define profile_stop(i) /* nothing */
#endif

enum profile_counters {
	PROF_ALL = 0,
	PROF_SIEVE_ALL,
	PROF_INIT_SIEVE,
	PROF_SIEVE,
	PROF_ALPHA1,
	PROF_ALPHA2,
	PROF_EVAL,
	PROF_INITIAL_OPTIMIZE,
	PROF_OPTIMIZE2,
	PROF_OPTIMIZE3,
	PROF_ROOTS,
	PROF_MURPHY,
	PROF_MAX	/* must be last */
};

void stage2_stat_init(stage2_stat_t *stats);
void stage2_stat_free(stage2_stat_t *stats);

/*-----------------------------------------------------------------------*/
/* data used in the current polynomial */

typedef struct {
	mpz_t gmp_a[6];
	mpz_t gmp_b[6];
	mpz_t gmp_lina[2];
	mpz_t gmp_linb[2];
	mpz_t gmp_help1;
	mpz_t gmp_help2;
	mpz_t gmp_help3;
	mpz_t gmp_help4;
	mpz_t gmp_p;
	mpz_t gmp_d;
	mpz_t gmp_m;
	mpz_t gmp_mb;
} curr_poly_t;

void curr_poly_init(curr_poly_t *poly);
void curr_poly_free(curr_poly_t *poly);

/*-----------------------------------------------------------------------*/
/* routines for translating polynomials */

void optimize_1(curr_poly_t *c, double *skewness, 
		double *pol_norm, double *alpha_proj);
void optimize_2(curr_poly_t *c, double skewness, 
		double *new_skewness, double *norm_ptr);
void optimize_3(curr_poly_t *c, poly_stage2_t *data, 
		double skewness, double *new_skewness, 
		double *norm_ptr, double *eptr, double *alphaptr);

/*-----------------------------------------------------------------------*/
/* data for the root sieve */

typedef struct {
	int xmin;
	int xmax;
	int ymin;
	int ymax;
} rs_bound_t;

typedef struct {
	int primepower;		/* p^k */
	int rooti;		/* zero =i mod p^k */
	int J;			/* f(i)+J*(i-m)=0 mod p^k */
	int step;		/* stepwidth p^(k-l) */
	int start;		/* sievearray[start] has zero =i mod p^k */
	int lineinc;		/* increment for change y->y+1 */
	unsigned short v;	/* scaled value */
	double value;		/* log(p)/(p^(k-1)*(p+1)) */
} primelist_t;

typedef struct {
	primelist_t *pl;
	int primelistlen;
	double limit; 
	double st_alpha;
	unsigned short *sievearray;
        unsigned short default_cut;
	int sievelen; 
	int nsubsieves;

	int prep_len; 
	int prep_p_len[NAFF_PRIMES];
	unsigned int *prep_p_begin[NAFF_PRIMES];
	unsigned int *prep_p[NAFF_PRIMES];
} root_sieve_t;

void root_sieve_init(root_sieve_t *rs);
void root_sieve_free(root_sieve_t *rs);
double compute_proj_alpha(curr_poly_t *c);
void root_sieve_run(curr_poly_t *c, double log_max_norm_2, 
		poly_stage2_t *data, root_sieve_t *rs,
		stage2_stat_t *stats, double skewness, 
		double pol_norm, double alpha_proj);

#ifdef HAVE_ASM_INTEL
void asm_root_sieve8(unsigned int **p1, unsigned int *p2, int l1,
		     unsigned int *p4, int l2);
#endif

/*-----------------------------------------------------------------------*/
/* routines for rating polynomial yield */

void check(int x, int y, curr_poly_t *c, 
		poly_stage2_t *data, stage2_stat_t *stats,
		double skewness);
void murphy_e_core(double *me, int deg0, double *dbl_coeff0, 
		int deg1, double *dbl_coeff1, double alpha0, 
		double alpha1, double skewness, int nsm);
void murphy_e(double *me, int deg0, double *dbl_coeff0, 
		int deg1, double *dbl_coeff1, double alpha0, 
		double alpha1, double skewness);

#ifdef __cplusplus
}
#endif

#endif /* !_STAGE2_IMPL_H_ */
