#ifndef _STAGE1_IMPL_H_
#define _STAGE1_IMPL_H_

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "ggnfs.h"
#include "poly_stage1.h"
#include "profile.h"
#include "if.h"

#if defined(HAVE_ASM_INTEL)
#define HAVE_ASM
#endif

/* We need to write a floorl() for Cygwin. In the meantime: 
   MinGW *must* use HAVE_FLOAT64, and we *must* currently compile 
   with -O0 (at least for gcc 3.4). This needs to be investigated */

#if defined(__CYGWIN__)
#undef HAVE_FLOAT64
#endif

#define P0_MAX             46300	/* <=2^15.5 */
#define DEG 5
#define MULTIPLIER            60	/* 2*2*3*5 */
#define NPR5 50
#define HASHSHIFT  52 /*dont change this (or disable/adapt the asm-functions) */
#define NHASH   (1<<(64-HASHSHIFT))
#define HASHSHIFT32  (HASHSHIFT-32)

#include "stage1_inline.h"

#ifdef __cplusplus
extern "C" {
#endif

/*-----------------------------------------------------------------------*/

/* profiling and statistics */

typedef struct {
	u64 n_a5;
	u64 n_pr;
	u64 n_p;
	u64 n_p_p0;
	u64 n_raw;
	u64 n_polexpand;
	u64 n_poltranslate;
	u64 n_survivors;
	profile_t profile;
} stage1_stat_t;

void stage1_stat_init(stage1_stat_t *stats);
void stage1_stat_free(stage1_stat_t *stats);

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
	PROF_COMPUTE_ROOT,
	PROF_NEXT_P,
	PROF_INIT_P_PRIMES,
	PROF_INIT_AUX_PRIMES,
	PROF_KNAPSACK_ALL,
	PROF_INIT_KNAPSACK_RAW,
	PROF_KNAPSACK_RAW,
	PROF_KNAP_RAW_COMBINE,
	PROF_CRT_AUX,
	PROF_KAPPA,
	PROF_KAPPA2,
	PROF_COMPUTE_KNAP_VALS,
	PROF_RAW_HASH1,
	PROF_RAW_HASH2,
	PROF_RAW_HASH3,
	PROF_RAW_HASH4,
	PROF_INIT_KNAPSACK_EXACT,
	PROF_KNAPSACK_EXACT,
	PROF_KNAP_EXACT_CHECK,
	PROF_KNAP_RAW_CHECK,
	PROF_CHECK_POLY,
	PROF_MAX 		/* must be last */
};

/*-----------------------------------------------------------------------*/
/* search bounds */

typedef struct {
	mpz_t gmp_a5_begin;
	mpz_t gmp_a5_end;
	double norm_max; 
	double skewness_min; 
	double skewness_max;
	double a3_max;
	double p_size_max;
	unsigned int p0_limit;
} bounds_t;

void stage1_bounds_init(bounds_t *bounds, poly_stage1_t *data);
void stage1_bounds_free(bounds_t *bounds);
void stage1_bounds_update(bounds_t *bounds, double N, double a5);

/*-----------------------------------------------------------------------*/

/* data for the current polynomial */

typedef struct {
	mpz_t gmp_N; 
	mpz_t gmp_prod; 
	mpz_t gmp_prod_base; 
	mpz_t gmp_root;
	mpz_t gmp_help1;
	mpz_t gmp_help2;
	mpz_t gmp_help3;
	mpz_t gmp_help4;
	mpz_t gmp_a5;
	mpz_t gmp_d;
	mpz_t gmp_m0;
	FILE *outputfile;
} curr_poly_t;

void curr_poly_init(curr_poly_t *poly, bounds_t *bounds, 
			poly_stage1_t *data);
void curr_poly_free(curr_poly_t *poly);

/*-----------------------------------------------------------------------*/

/* data for the primes making up the leading
   rational polynomial coefficient */

typedef struct {
	u32 pr;
	u32 pr_start;
	u32 pr_step;
	double p_log;
} rat_prime_t;

typedef struct {
	rat_prime_t primes[NPR5];
	unsigned char *bitarray_5power[NPR5];
} rat_factor_t;

void rat_factor_init(rat_factor_t *r, mpz_t N, mpz_t a5_start);
void rat_factor_free(rat_factor_t *r);

/*-----------------------------------------------------------------------*/

/* data used for the current rational polynomial */

typedef struct {
	u32 p_pr;
	u32 p_fr[DEG];
	double p_log;
	u32 p_N_inv;
	u32 p_N_mod_p2;
	u32 p_a5_mod_p2;
	u32 p_minus5a5_mod_p;
	u64 ull_p_inv;
#ifdef HAVE_FLOAT64
	long double ld_p_inv;
#endif
} curr_rat_prime_t;

typedef struct {
	u32 p;
	u32 r;
} aux_factor_t;

typedef struct {
	u32 npr_in_p;
	u32 npr_total;
	u32 p_ind[NPR5];
	u32 p_inv_table[NPR5][NPR5];
	curr_rat_prime_t primes[NPR5];

	u32 p0_list_len;
	u32 p0_list_len_max;
	u32 last_p0;
	aux_factor_t *p0_list;
} curr_rat_factor_t;

void curr_rat_factor_init(curr_rat_factor_t *c, u32 npr_in_p);
void curr_rat_factor_free(curr_rat_factor_t *c);
u32 curr_rat_factor_find(rat_factor_t *rat_factors,
			  curr_rat_factor_t *curr_rat_factors,
			  double p_size_max);
void curr_rat_factor_fill(curr_rat_factor_t *c, mpz_t N, mpz_t a5);
void curr_rat_factor_fill_aux(curr_rat_factor_t *c, mpz_t N, mpz_t a5,
				double p0_max, u32 p0_limit);

/*-----------------------------------------------------------------------*/
/* data for the knapsack code */

typedef struct {
	mpz_t gmp_disp;
	mpz_t gmp_D0;
	mpz_t gmp_D[NPR5][5];
	mpz_t gmp_kappa0;
	mpz_t gmp_kappa[NPR5][5];
	mpz_t gmp_kappa_help[NPR5 + 1];
	double dbl_kappa_help[NPR5 + 1];
	
	uint64_t ull_kappa0;
	uint64_t ull_kappa[NPR5][5];
	uint64_t raw_ull_kappa[NPR5][5];
	unsigned int ul_d_help[NPR5][5];
	uint64_t *s1;
	uint64_t *s2;
	uint64_t *s1sort;
	uint64_t *s2sort;
	int s1len;
	int s2len;
	uint64_t ull_bound;
	uint64_t lambda0;
	uint64_t raw_ull_bound;
	uint64_t *stored_pairs;
	uint64_t *raw_stored_pairs;
	int store_len;
	int nstore;
	int raw_store_len;
	int nraw_store;
	int *raw_cand;
	int nraw_cand;
	
	unsigned int *raw_cand_hash;
	
	uint64_t ull_kappa_p_inv[NPR5 + 1];
	unsigned int ull_kappa_help0[NPR5 + 1];
	
	unsigned int prep_p[NPR5];
	unsigned int prep_5a5[NPR5];
	unsigned int prep_N_mod_p2[NPR5];
	unsigned int prep_other_mod_p2[NPR5];
	unsigned int prep_inv_table[NPR5][NPR5];
	
	unsigned int *hashdata;
	unsigned int *s1sortl, *s2sortl;
	
	uint64_t *shelp;
	int len1;
	int len2;
	int len11;
	int len12;
	int len21;
	int len22;
	
	unsigned int **hashptr32;
	int *sort;
	uint64_t **hashptr;
	
	unsigned int ul_d_help_1[NPR5][5];
	unsigned int p_mod_p2[NPR5];
	unsigned int p0_inv_p[NPR5];
	unsigned int p_inv_p0[NPR5];
	
	unsigned int raw_bound;
	int *raw_cand_ptr;
	unsigned int *hashdataptr;
	unsigned int *s11l;
	unsigned int *s12l;
	unsigned int *s21l;
	unsigned int *s22l;
	unsigned int s11len;
	unsigned int s12len;
	unsigned int s21len;
	unsigned int s22len;
} knapsack_data_t;

void knapsack_data_init(knapsack_data_t *knap, u32 npr_in_p);
void knapsack_data_free(knapsack_data_t *knap);

/*-----------------------------------------------------------------------*/

int check_pol(curr_poly_t *curr, bounds_t *bounds, stage1_stat_t *stats);

/*
Given p (product of small primes) the following functions search for 
suitable d such that we can find a polynomial of degree 5 with d/p as
zero and 'small' a5, a4 and a3.
With reasonable parameters for most p there does not exist a suitable
d. Therefore the search is first done in an approximate manner
(knapsack_raw) and then exact (knapsack_exact).
*/

void init_knapsack_raw(curr_poly_t *curr, knapsack_data_t *knap,
			curr_rat_factor_t *curr_factors,
			double a3_max, stage1_stat_t *stats);
void init_knapsack_raw_p0(curr_poly_t *curr, knapsack_data_t *knap,
			curr_rat_factor_t *curr_factors,
			unsigned int p0, unsigned int r0, 
			double a3_max, stage1_stat_t *stats);
int knapsack_raw(knapsack_data_t *knap, stage1_stat_t *stats, int npr_in_p);
void check_raw(curr_poly_t *curr, knapsack_data_t *knap,
		bounds_t *bounds, stage1_stat_t *stats);

void init_knapsack_exact(curr_poly_t *curr, knapsack_data_t *knap,
			curr_rat_factor_t *curr_factors,
			double a3_max);
void init_knapsack_exact_p0(curr_poly_t *curr, knapsack_data_t *knap,
			curr_rat_factor_t *curr_factors,
			unsigned int p0, unsigned int r0, double a3_max);
void knapsack_exact(curr_poly_t *curr, knapsack_data_t *knap,  
			bounds_t *bounds, stage1_stat_t *stats);

#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_IMPL_H_ */
