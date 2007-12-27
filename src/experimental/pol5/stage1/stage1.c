#include "stage1_impl.h"

/*
Given n primes p_i such that a_5*x^5=N has 5 five solution modulo p_i
we want to compute all 5^n solutions of this congruence modulo the product
of the p_i and find an expansion N=a_5*d^5+a_4*p*d^4+a_3*p^2*d^3 where
a3 is 'small'. Let prod be the product of the p_i and fr_{ij} a fifth root
modulo p_i of N/a_5 (0<=j<5). Let d_{ij}=h_{ij}*prod/p_i such that
d_{ij}=fr_{ij} modulo p_i and 0<=d_{ij}<prod. Then the fifth roots
modulo prod of N/a5 are: sum_i d_{i,j_i}. Let D_{ij}=d_{ij}-d_{i0}.
We choose D_0 near (N/a_5)^0.2 such that D_0=sum_i d_{i0} modulo prod.
Then we can use d=D_0+sum_i D_{i,j_i} for the expansion mentioned above.
The value of a_4 corresponding to (j_i) can be written as
kappa_0+sum_i kappa_{i,j_i} modulo prod with kappa_{i0}=0. We consider 
these kappa's as integers between 0 and prod. They can be computed as follows:
kappa_{i,j}=-5*a_5*(h_{ij}-h_{i0})/p_i modulo p_k, k!=i and
kappa_{i,j}=a_5/N*(p_i/prod)*
           ((N*D_{i,j}-a_5*((d_0+D_{i,j})^6-D_0^6))/p_i) modulo p_i
and applying the chinese remainder theorem.
Now we can estimate the a_3's:
Let lambda_0 be such that:
  ((a_5*D_0^5-N)/prod+kappa_0*D_0^4)/prod=lambda_0*D_0^4
(we can replace the two D_0^4 by (N/a_5)^0.8 without loosing too much
precision).
For the other d with corresponding kappa we have:
((a_5*d^5-N)/prod+kappa*d^4)=
 (lambda_0+(kappa-kappa_0)/prod-5*a_5*(d-D_0)/prod^2+
  10*a_5*(d-D_0)^2/prod^2/d+...)*d^4
where ... are very small terms (of order n^2*a_5/d^2).
Since (kappa-kappa_0) and (d-D_0) depend 'linearly' on (j_i) this
approximation of a_3/d can be written as lambda_0+sum_i lambda_{i,j_i}
and we have to solve this knapsack problem. */

/*------------------------------------------------------------------------*/
static void
search_p(curr_poly_t *curr, curr_rat_factor_t *curr_factors,
		knapsack_data_t *knap, bounds_t *bounds, 
		stage1_stat_t *stats, unsigned int p0_max)
{
	int i, j;
	int passed;
	unsigned int r0;
	unsigned int p0 = 1;
	u32 npr_in_p = curr_factors->npr_in_p;

	for (i = 0; i < curr_factors->p0_list_len; i++) {
		curr_factors->last_p0 = p0;
		p0 = curr_factors->p0_list[i].p;
		r0 = curr_factors->p0_list[i].r;
		if (p0 > p0_max)
			break;

		if (verbose > 3) {
			u32 *p_ind = curr_factors->p_ind;
			curr_rat_prime_t *primes = curr_factors->primes;
			printf("search for: ");
			for (j = 0; j < npr_in_p; j++)
				printf("%d ", primes[p_ind[j]].p_pr);
			if (p0 != 1)
				printf("   %u (%u)", p0, r0);
			printf("\n");
		}
		stats->n_p++;
		mpz_mul_ui(curr->gmp_prod, curr->gmp_prod_base, p0);

		profile_start(PROF_INIT_KNAPSACK_RAW);
		if (p0 == 1) {
			init_knapsack_raw(curr, knap, curr_factors,
					bounds->a3_max, stats);
		}
		else {
			stats->n_p_p0++;
			init_knapsack_raw_p0(curr, knap, curr_factors,
					p0, r0, bounds->a3_max, stats);
		}
		profile_stop(PROF_INIT_KNAPSACK_RAW);

		profile_start(PROF_KNAPSACK_RAW);
		passed = knapsack_raw(knap, stats, npr_in_p);
		profile_stop(PROF_KNAPSACK_RAW);

		if (passed > 0) {
			profile_start(PROF_INIT_KNAPSACK_EXACT);
			if (p0 == 1) {
				init_knapsack_exact(curr, knap, curr_factors,
						bounds->a3_max);
			}
			else {
				init_knapsack_exact_p0(curr, knap, curr_factors,
						p0, r0, bounds->a3_max);
			}
			profile_stop(PROF_INIT_KNAPSACK_EXACT);
			profile_start(PROF_KNAP_RAW_CHECK);
			check_raw(curr, knap, bounds, stats);
			profile_stop(PROF_KNAP_RAW_CHECK);

		}
		else if (passed < 0) {
			profile_start(PROF_INIT_KNAPSACK_EXACT);
			if (p0 == 1) {
				init_knapsack_exact(curr, knap, curr_factors,
						bounds->a3_max);
			}
			else {
				init_knapsack_exact_p0(curr, knap, curr_factors,
						p0, r0, bounds->a3_max);
			}
			profile_stop(PROF_INIT_KNAPSACK_EXACT);
			profile_start(PROF_KNAPSACK_EXACT);
			knapsack_exact(curr, knap, bounds, stats);
			profile_stop(PROF_KNAPSACK_EXACT);
		}
	}
}

/*------------------------------------------------------------------------*/
static void
search_a5(curr_poly_t *curr, rat_factor_t *rat_factors,
		curr_rat_factor_t *curr_factors,
		knapsack_data_t *knap,
		bounds_t *bounds, stage1_stat_t *stats)
{
	int i, j;
	double p_size, dp0_max, dbl_a5, dbl_N;
	unsigned int p0_max;
	unsigned int shift = 0;
	mpz_t gmp_tmp;
	u32 npr_in_p;
	u32 npr_excess;
	u32 *p_ind;
	curr_rat_prime_t *primes;

	profile_start(PROF_ALL);
	dbl_N = mpz_get_d(curr->gmp_N);
	dbl_a5 = mpz_get_d(curr->gmp_a5);
	stage1_bounds_update(bounds, dbl_N, dbl_a5);
	mpz_init(gmp_tmp);

	while (1) {
		stats->n_a5++;
		shift++;
		if (shift >= 1000) {
			/* update some variables, not very precise 
			   at the beginning of a5-range */
			mpz_set_ui(gmp_tmp, shift);
			mpz_mul_ui(gmp_tmp, gmp_tmp, MULTIPLIER);
			mpz_add(curr->gmp_a5, curr->gmp_a5, gmp_tmp);
			if (mpz_cmp(curr->gmp_a5, bounds->gmp_a5_end) > 0)
				break;
			dbl_a5 = mpz_get_d(curr->gmp_a5);
			stage1_bounds_update(bounds, dbl_N, dbl_a5);
			shift = 0;
		}
		if (!curr_rat_factor_find(rat_factors, 
					curr_factors, bounds->p_size_max)) {
			continue;
		}

		stats->n_pr++;

		/* it actually is much faster to let GMP compute
		   the fifth root from scratch than to use the
		   previous root as a starting point and compute
		   the Newton iteration manually */

		profile_start(PROF_COMPUTE_ROOT);
		mpz_set_ui(gmp_tmp, shift);
		mpz_mul_ui(gmp_tmp, gmp_tmp, MULTIPLIER);
		mpz_add(curr->gmp_a5, curr->gmp_a5, gmp_tmp);
		shift = 0;

		mpz_fdiv_q(curr->gmp_root, curr->gmp_N, curr->gmp_a5);
		mpz_root(curr->gmp_root, curr->gmp_root, 5);
		mpz_set(curr->gmp_m0, curr->gmp_root);
		profile_stop(PROF_COMPUTE_ROOT);

		if (mpz_cmp(curr->gmp_a5, bounds->gmp_a5_end) > 0)
			break;
		if (verbose > 3) {
			mpz_out_str(stdout, 10, curr->gmp_a5);
			printf(" \n");
		}
		profile_start(PROF_INIT_P_PRIMES);
		curr_rat_factor_fill(curr_factors, curr->gmp_N,
						curr->gmp_a5);
		profile_stop(PROF_INIT_P_PRIMES);

		p_size = 0;
		dbl_a5 = mpz_get_d(curr->gmp_a5);
		stage1_bounds_update(bounds, dbl_N, dbl_a5);

		p_size = 0;
		p_ind = curr_factors->p_ind;
		primes = curr_factors->primes;
		npr_in_p = curr_factors->npr_in_p;
		npr_excess = curr_factors->npr_total - npr_in_p;
		for (i = 0; i < npr_in_p; i++) {
			p_ind[i] = i;
			p_size += primes[i].p_log;
		}

		profile_start(PROF_INIT_AUX_PRIMES);
		curr_rat_factor_fill_aux(curr_factors, 
				curr->gmp_N, curr->gmp_a5,
				exp(bounds->p_size_max - p_size),
				 bounds->p0_limit);
		profile_stop(PROF_INIT_AUX_PRIMES);

		while (1) {
			dp0_max = exp(bounds->p_size_max - p_size);
			if (dp0_max > (double) (bounds->p0_limit))
				p0_max = bounds->p0_limit;
			else
				p0_max = (unsigned int) (rint(dp0_max) - 1);
			if (p0_max < 1)
				p0_max = 1;

			mpz_set_ui(curr->gmp_prod_base, 1);
			for (i = 0; i < npr_in_p; i++) {
				mpz_mul_ui(curr->gmp_prod_base, 
					curr->gmp_prod_base, 
					primes[p_ind[i]].p_pr);
			}

			profile_start(PROF_KNAPSACK_ALL);
			search_p(curr, curr_factors, knap,
					bounds, stats, p0_max);
			profile_stop(PROF_KNAPSACK_ALL);

			profile_start(PROF_NEXT_P);
			i = npr_in_p - 1;
			while (i >= 0) {
				if (p_size > bounds->p_size_max) {
					i--;
					while ((i >= 0) && (p_ind[i] + 1 ==
							p_ind[i + 1])) {
						i--;
					}
					if (i < 0)
						break;
					p_ind[i]++;
					for (j = i + 1; j < npr_in_p; j++)
						p_ind[j] = p_ind[j - 1] + 1;
					p_size = 0;
					for (j = 0; j < npr_in_p; j++)
						p_size += primes[p_ind[j]].p_log;
					if (p_size <= bounds->p_size_max)
						break;
					continue;
				}
				while ((i >= 0) && (p_ind[i] >= i + npr_excess))
					i--;
				if (i < 0)
					break;
				p_ind[i]++;
				for (j = i + 1; j < npr_in_p; j++)
					p_ind[j] = p_ind[j - 1] + 1;
				p_size = 0;
				for (j = 0; j < npr_in_p; j++)
					p_size += primes[p_ind[j]].p_log;
				if (p_size <= bounds->p_size_max)
					break;
			}
			profile_stop(PROF_NEXT_P);
			if (i < 0)
				break;
		}
	}

	mpz_clear(gmp_tmp);
	profile_stop(PROF_ALL);
}

/*------------------------------------------------------------------------*/
void
poly_stage1_init(poly_stage1_t *data)
{
	mpz_init_set_ui(data->gmp_N, 0);
	mpz_init_set_ui(data->gmp_a5_begin, 0);
	mpz_init_set_ui(data->gmp_a5_end, 0);
	data->npr_in_p = 0;
	data->norm_max = 0;
	data->p0_limit = 0;
}

/*------------------------------------------------------------------------*/
void
poly_stage1_free(poly_stage1_t *data)
{
	mpz_clear(data->gmp_N);
	mpz_clear(data->gmp_a5_begin);
	mpz_clear(data->gmp_a5_end);
}

/*------------------------------------------------------------------------*/
s32
poly_stage1_run(poly_stage1_t *data)
{
	u32 success = 0;
	stage1_stat_t stats;
	bounds_t bounds;
	curr_poly_t currpoly;
	rat_factor_t rat_factors;
	curr_rat_factor_t curr_rat_factors;
	knapsack_data_t knapdata;

	stage1_stat_init(&stats);
	stage1_bounds_init(&bounds, data);
	curr_poly_init(&currpoly, &bounds, data);
	rat_factor_init(&rat_factors, currpoly.gmp_N, currpoly.gmp_a5);
	curr_rat_factor_init(&curr_rat_factors, data->npr_in_p);
	knapsack_data_init(&knapdata, data->npr_in_p);

	search_a5(&currpoly, &rat_factors, &curr_rat_factors, &knapdata,
			&bounds, &stats);

	if (stats.n_survivors)
		success = 1;

	if (verbose) {
		printf("Statistics:\n");
		printf("a5-values: %" PRIu64 ", suitable for %d primes: %"
		       PRIu64 "\n", stats.n_a5, data->npr_in_p, stats.n_pr);
		printf("raw checks of (a5,p): %" PRIu64 " (%" PRIu64
		       "),  fine checks of (a5,p): %" PRIu64 "\n", stats.n_p,
		       stats.n_p_p0, stats.n_raw);
		printf("polynomials computed: %" PRIu64 "(%" PRIu64 "),  "
		       "survivors: %" PRIu64 "\n", 
		       stats.n_polexpand, stats.n_poltranslate,
		       stats.n_survivors);
		printf("Total number of checked polynomials: %" PRIu64 "\n",
		       (uint64_t) (pow(5.0, (double)data->npr_in_p)) * 
		       stats.n_p);
	}

#ifdef DO_PROFILE
	profile_done(&stats.profile);
	printf("\nTiming:");
	printf("\nTotal            "); 
	profile_print(&stats.profile, PROF_ALL);
	printf("\n  Root             "); 
	profile_print(&stats.profile, PROF_COMPUTE_ROOT);
	printf("\n  Init primes      "); 
	profile_print(&stats.profile, PROF_INIT_P_PRIMES);
	printf("\n  Init aux primes      "); 
	profile_print(&stats.profile, PROF_INIT_AUX_PRIMES);
	printf("\n  Knapsack         "); 
	profile_print(&stats.profile, PROF_KNAPSACK_ALL);
	printf("\n    Init knaps. raw  "); 
	profile_print(&stats.profile, PROF_INIT_KNAPSACK_RAW);
	printf("\n      D's              "); 
	profile_print(&stats.profile, PROF_CRT_AUX);
	printf("\n      kappa1          "); 
	profile_print(&stats.profile, PROF_KAPPA);
	printf("\n        kappa2          "); 
	profile_print(&stats.profile, PROF_KAPPA2);
	printf("\n      rest             ");
       	profile_print(&stats.profile, PROF_COMPUTE_KNAP_VALS);
	printf("\n    Knapsack raw     "); 
	profile_print(&stats.profile, PROF_KNAPSACK_RAW);
	printf("\n      combine          "); 
	profile_print(&stats.profile, PROF_KNAP_RAW_COMBINE);
	printf("\n      hash1      "); 
	profile_print(&stats.profile, PROF_RAW_HASH1);
	printf("\n      hash2      "); 
	profile_print(&stats.profile, PROF_RAW_HASH2);
	printf("\n      hash3      "); 
	profile_print(&stats.profile, PROF_RAW_HASH3);
	printf("\n      hash4      "); 
	profile_print(&stats.profile, PROF_RAW_HASH4);
	printf("\n      check      "); 
	profile_print(&stats.profile, PROF_KNAP_RAW_CHECK);
	printf("\n    Init knaps. ex.  "); 
	profile_print(&stats.profile, PROF_INIT_KNAPSACK_EXACT);
	printf("\n    Knapsack ex.       "); 
	profile_print(&stats.profile, PROF_KNAPSACK_EXACT);
	printf("\n      check            "); 
	profile_print(&stats.profile, PROF_KNAP_EXACT_CHECK);
	printf("\n    Check pol        "); 
	profile_print(&stats.profile, PROF_CHECK_POLY);
	printf("\n  Next p           "); 
	profile_print(&stats.profile, PROF_NEXT_P);
	printf("\n");
#endif
	stage1_bounds_free(&bounds);
	curr_poly_free(&currpoly);
	rat_factor_free(&rat_factors);
	curr_rat_factor_free(&curr_rat_factors);
	stage1_stat_free(&stats);
	knapsack_data_free(&knapdata);
	return success;
}
