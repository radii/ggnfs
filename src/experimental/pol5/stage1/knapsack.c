#include "stage1_impl.h"

/*
The following is the exact version of the 'raw' functions.
This is very time-critical and not optimized.
In fact this is quite old and should be revised.
*/

/*------------------------------------------------------------------------*/
static uint64_t
renorm(curr_poly_t *curr, mpz_t n)
{
	uint64_t res;

	mpz_mul_2exp(curr->gmp_help1, n, 64);
	mpz_fdiv_q(curr->gmp_help1, curr->gmp_help1, curr->gmp_prod);
	if (mpz_sizeinbase(curr->gmp_help1, 2) > 64)
		complain("renorm\n");
	mpz_get_ull_64(&res, curr->gmp_help1);
	return res;
}

/*------------------------------------------------------------------------*/
static uint64_t
renorm2(curr_poly_t *curr, mpz_t n)
{
	uint64_t res;

	mpz_mul_2exp(curr->gmp_help1, n, 64);
	mpz_mul(curr->gmp_help1, curr->gmp_help1, curr->gmp_a5);
	mpz_mul_ui(curr->gmp_help1, curr->gmp_help1, 5);
	mpz_fdiv_q(curr->gmp_help1, curr->gmp_help1, curr->gmp_prod);
	mpz_fdiv_q(curr->gmp_help1, curr->gmp_help1, curr->gmp_prod);
	mpz_get_ull_64(&res, curr->gmp_help1);
	return res;
}

/*------------------------------------------------------------------------*/
static void
check_stored_pairs(curr_poly_t *curr, knapsack_data_t *knap,
		bounds_t *bounds, stage1_stat_t *stats)
{
	int i, j, i1, i2;
	uint64_t *stored_pairs = knap->stored_pairs;
	int nstore = knap->nstore;
	uint64_t *s1 = knap->s1;
	uint64_t *s2 = knap->s2;
	int s1len = knap->s1len;
	int s2len = knap->s2len;
	int len1 = knap->len1;
	int len2 = knap->len2;

	i = 0;
	while (i < nstore - 1) {
		for (j = i + 1; j < nstore; j++)
			if ((stored_pairs[2 * i] == stored_pairs[2 * j]) &&
			    (stored_pairs[2 * i + 1] ==
			     stored_pairs[2 * j + 1]))
				break;
		if (j < knap->nstore) {
			stored_pairs[2 * j] = stored_pairs[2 * nstore - 2];
			stored_pairs[2 * j + 1] = stored_pairs[2 * nstore - 1];
			nstore--;
		}
		i++;
	}			/* slow algorithm since nstore usually is 1 */
	knap->nstore = nstore;
	for (i = 0; i < nstore; i++) {

		/* TODO: there might be more than one i1 
		   satisfying s1[i1]==stored_pairs[2*i] */

		for (i1 = 0; i1 < s1len; i1++) {
			if (stored_pairs[2 * i] == s1[i1])
				break;
		}
		if (i1 >= s1len)
			complain("check_stored_pairs.1\n");
		for (i2 = 0; i2 < s2len; i2++) {
			if (stored_pairs[2 * i + 1] == s2[i2])
				break;
		}
		if (i2 >= s2len)
			complain("check_stored_pairs.2\n");

		mpz_set_ui(curr->gmp_d, 0);
		mpz_set_ui(curr->gmp_help4, 0);
		for (j = 0; j < len1; j++) {
			mpz_add(curr->gmp_d, curr->gmp_d, knap->gmp_D[j][i1 % 5]);
			mpz_add(curr->gmp_help4, curr->gmp_help4, 
					knap->gmp_kappa[j][i1 % 5]);
			i1 -= (i1 % 5);
			i1 /= 5;
		}
		for (j = len1; j < len1 + len2; j++) {
			mpz_add(curr->gmp_d, curr->gmp_d, knap->gmp_D[j][i2 % 5]);
			mpz_add(curr->gmp_help4, curr->gmp_help4, 
					knap->gmp_kappa[j][i2 % 5]);
			i2 -= (i2 % 5);
			i2 /= 5;
		}
		mpz_add(curr->gmp_d, curr->gmp_d, knap->gmp_disp);
		mpz_add(curr->gmp_help4, curr->gmp_help4, knap->gmp_kappa0);
		/* check */
		mpz_mul(curr->gmp_help3, curr->gmp_d, curr->gmp_d);
		mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_help3);
		mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_d);
		mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_a5);
		mpz_sub(curr->gmp_help3, curr->gmp_help3, curr->gmp_N);
		mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help2, curr->gmp_help3, curr->gmp_prod);
		if (mpz_sgn(curr->gmp_help2))
			complain("check2\n");
		mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_d);
		mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_a5);
		mpz_mul(curr->gmp_help1, curr->gmp_N, curr->gmp_help4);
		mpz_add(curr->gmp_help3, curr->gmp_help3, curr->gmp_help1);
		mpz_fdiv_r(curr->gmp_help3, curr->gmp_help3, curr->gmp_prod);

		if (mpz_sgn(curr->gmp_help3)) {
			mpz_out_str(stdout, 10, curr->gmp_d);
			printf("\n");
			mpz_out_str(stdout, 10, curr->gmp_prod);
			printf("\n");
			mpz_out_str(stdout, 10, curr->gmp_help4);
			printf("\n");
			complain("check_stored_pairs check\n");
		}

		if (check_pol(curr, bounds, stats)) {

			mpz_out_str(curr->outputfile, 10, curr->gmp_a5);
			fprintf(curr->outputfile, " ");
			mpz_out_str(curr->outputfile, 10, curr->gmp_prod);
			fprintf(curr->outputfile, " ");
			mpz_out_str(curr->outputfile, 10, curr->gmp_d);
			fprintf(curr->outputfile, "\n");

			if (verbose > 1) {
				mpz_out_str(stdout, 10, curr->gmp_a5);
				printf(" ");
				mpz_out_str(stdout, 10, curr->gmp_d);
				printf(" ");
				mpz_out_str(stdout, 10, curr->gmp_prod);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_kappa0);
				printf(" ");
				mpz_out_str(stdout, 10, knap->gmp_D0);
				printf("\n");
			}
		}
	}
	if (verbose > 2)
		printf("%d pairs checked\n", nstore);
}

/*------------------------------------------------------------------------*/
static void
store(knapsack_data_t *knap, int i1, int i2)
{
	while (knap->nstore + 1 >= knap->store_len) {
		knap->store_len += 16;
		knap->stored_pairs = (uint64_t *)xrealloc(
					knap->stored_pairs,
					2 * knap->store_len * sizeof(uint64_t));
	}
	knap->stored_pairs[2 * knap->nstore] = knap->s1sort[i1];
	knap->stored_pairs[2 * knap->nstore + 1] = knap->s2sort[i2];
	knap->nstore++;
}

/*------------------------------------------------------------------------*/
static void
combine0(knapsack_data_t *knap, 
	uint64_t * targ, int len, int ind)
{
	int i, j, disp;
	uint64_t add;

	if (knap->ull_kappa[ind][0])
		complain("combine0\n");
	disp = len;
	for (j = 1; j < 5; j++) {
		add = knap->ull_kappa[ind][j];
		for (i = 0; i < len; i++)
			targ[i + disp] = add + targ[i];
		disp += len;
	}
}

/*------------------------------------------------------------------------*/
static void
combinelast0(knapsack_data_t *knap,
		uint64_t * targ, int len, int ind)
{
	int i, j, disp;
	uint64_t add, h;
	int *sort = knap->sort;

	if (knap->ull_kappa[ind][0])
		complain("combinelast0\n");
	memset(sort, 0, NHASH * sizeof(int));
	disp = 0;
	for (j = 0; j < 5; j++) {
		add = knap->ull_kappa[ind][j];
		for (i = 0; i < len; i++) {
			h = add + targ[i];
			targ[i + disp] = h;
			h >>= HASHSHIFT;
			sort[h]++;
		}
		disp += len;
	}
}

/*------------------------------------------------------------------------*/
static void
combine(knapsack_data_t *knap, uint64_t * targ, 
		uint64_t * targsort, int i0, int i1)
{
	int i, len;
	uint64_t h, hh, *ptr;
	int *sort = knap->sort;
	uint64_t **hashptr = knap->hashptr;

	if (i1 - i0 < 1)
		complain("combine\n");
	for (i = 0; i < 5; i++)
		targ[i] = knap->ull_kappa[i0][i];
	len = 5;
	for (i = i0 + 1; i < i1 - 1; i++) {
		combine0(knap, targ, len, i);
		len *= 5;
	}
	combinelast0(knap, targ, len, i1 - 1);
	len *= 5;
	hashptr[0] = targsort;
	for (i = 0; i < NHASH - 1; i++) {
		hashptr[i + 1] = hashptr[i] + sort[i];
	}
	for (i = 0; i < len; i++) {
		h = targ[i];
		hh = h >> HASHSHIFT;
		*(hashptr[hh]) = h;
		hashptr[hh]++;
	}
	for (i = 0; i < NHASH; i++) {
		if (sort[i] > 3) {
			ptr = hashptr[i] - sort[i];
			while (ptr < hashptr[i] - 1) {
				if (ptr[0] > ptr[1]) {
					h = ptr[0];
					ptr[0] = ptr[1];
					ptr[1] = h;
					if (ptr + sort[i] > hashptr[i])
						ptr--;
				}
				else {
					ptr++;
				}
			}
			continue;
		}
		if (sort[i] > 1) {
			ptr = hashptr[i] - sort[i];
			if (sort[i] == 2) {
				if (ptr[0] > ptr[1]) {
					h = ptr[0];
					ptr[0] = ptr[1];
					ptr[1] = h;
				}
			}
			else {
				if (ptr[0] > ptr[1]) {
					h = ptr[0];
					ptr[0] = ptr[1];
					ptr[1] = h;
				}
				if (ptr[1] > ptr[2]) {
					if (ptr[0] > ptr[2]) {
						h = ptr[2];
						ptr[2] = ptr[1];
						ptr[1] = ptr[0];
						ptr[0] = h;
					}
					else {
						h = ptr[1];
						ptr[1] = ptr[2];
						ptr[2] = h;
					}
				}
			}
		}
	}
}

/*------------------------------------------------------------------------*/
void
init_knapsack_exact(curr_poly_t *curr, 
			knapsack_data_t *knap,
			curr_rat_factor_t *curr_factors,
			double a3_max)
{
	int i, j, l;
	unsigned int p, h, hh, inv;
	uint64_t lambda0;
	unsigned int dp, d0p, p2, hp, N_mod_p2, an_mod_p2;
	double db, dq;
	int npr_in_p = curr_factors->npr_in_p;
	curr_rat_prime_t *primes = curr_factors->primes;
	u32 *p_ind = curr_factors->p_ind;

	/* compute D_{i,j} */
	p = 1;
	db = 0.;
	for (i = 0; i < npr_in_p; i++) {
		curr_rat_prime_t *c = primes + p_ind[i];
		p = c->p_pr;
		if (mpz_fdiv_q_ui(curr->gmp_help1, curr->gmp_prod, p))
			complain("init_knapsack_exact.prod\n");
		h = mpz_fdiv_ui(curr->gmp_help1, p);
		inv = mp_modinv_1(h, p);
		mpz_mul_ui(knap->gmp_kappa_help[i], curr->gmp_help1, inv);
		knap->p_mod_p2[i] = mpz_fdiv_ui(curr->gmp_help1, p * p);

		knap->dbl_kappa_help[i] = (double) inv / (double) p;

		h = inv * c->p_fr[0];
		h %= p;
		knap->ul_d_help[i][0] = h;
		mpz_mul_ui(knap->gmp_D[i][0], curr->gmp_help1, h);
		db += (double) h / (double) p;

		for (j = 1; j < DEG; j++) {
			h = inv * c->p_fr[j];
			h %= p;
			h += (p - knap->ul_d_help[i][0]);
			if (h >= p)
				h -= p;
			knap->ul_d_help[i][j] = h;
			h += knap->ul_d_help[i][0];
			mpz_mul_ui(knap->gmp_D[i][j], curr->gmp_help1, h);
		}
	}
	db += 0.5 * (double) npr_in_p;

	mpz_set_ui(knap->gmp_D0, 0);
	for (i = 0; i < npr_in_p; i++)
		mpz_add(knap->gmp_D0, knap->gmp_D0, knap->gmp_D[i][0]);
	mpz_fdiv_r(curr->gmp_help1, curr->gmp_m0, curr->gmp_prod);
	mpz_sub(knap->gmp_disp, curr->gmp_m0, curr->gmp_help1);
	mpz_mul_si(curr->gmp_help1, curr->gmp_prod, (int) db);
	mpz_sub(knap->gmp_disp, knap->gmp_disp, curr->gmp_help1);
	mpz_add(knap->gmp_D0, knap->gmp_D0, knap->gmp_disp);

	/* compute kappa_0 and kappa_{i,j} */
	mpz_set_ui(curr->gmp_help1, 0);
	for (i = 0; i < npr_in_p; i++) {
		curr_rat_prime_t *c = primes + p_ind[i];

		mpz_set_ui(knap->gmp_kappa[i][0], 0);
		p = c->p_pr;
		p2 = p * p;
		d0p = mpz_fdiv_ui(knap->gmp_D0, p2);
		an_mod_p2 = c->p_a5_mod_p2;
		N_mod_p2 = c->p_N_mod_p2;

		hp = mp_expo_1(d0p, DEG, p2);
		hp = mp_modmul_1(hp, an_mod_p2, p2);
		hp = p2 - hp + N_mod_p2;
		if (hp % p)
			complain("%u %u %u neu-1e\n", p, d0p, hp);
		h = hp / p;
		inv = mp_modmul_1(an_mod_p2, c->p_N_inv, p);
		for (j = 0; j < i; j++) {
			inv = mp_modmul_1(inv, 
				curr_factors->p_inv_table[p_ind[j]][p_ind[i]], p);
		}
		for (j = i + 1; j < npr_in_p; j++) {
			inv = mp_modmul_1(inv, 
				curr_factors->p_inv_table[p_ind[j]][p_ind[i]], p);
		}
		h = mp_modmul_1(h, inv, p);
		h = mp_modmul_1(h, d0p, p);
		mpz_mul_ui(curr->gmp_help2, knap->gmp_kappa_help[i], h);
		mpz_add(curr->gmp_help1, curr->gmp_help1, curr->gmp_help2);

		for (j = 1; j < DEG; j++) {
			dp = mp_modmul_1(knap->p_mod_p2[i], 
					knap->ul_d_help[i][j], p2);
			hh = d0p + dp;
			h = mp_expo_1(hh, DEG + 1, p2);
			hh = mp_expo_1(d0p, DEG + 1, p2);
			h = mp_modmul_1(hh + (p2 - h), an_mod_p2, p2);
			h += mp_modmul_1(N_mod_p2, dp, p2);

			if (h % p)
				complain("%u %u %u neu1\n", p, d0p, dp);
			h /= p;
			h = mp_modmul_1(h, inv, p);
			mpz_mul_ui(curr->gmp_help4, knap->gmp_kappa_help[i], h);
			dq = (double) h * knap->dbl_kappa_help[i];

			for (l = 0; l < npr_in_p; l++) {
				curr_rat_prime_t *cl = primes + p_ind[l];
				if (l == i)
					continue;
				h = mp_modmul_1(cl->p_minus5a5_mod_p * 
						knap->ul_d_help[i][j],
						curr_factors->p_inv_table[p_ind[i]][p_ind[l]],
						cl->p_pr);
				mpz_mul_ui(curr->gmp_help3, 
					   knap->gmp_kappa_help[l],
					   h);
				mpz_add(curr->gmp_help4, curr->gmp_help4,
					curr->gmp_help3);
				dq += (((double) h) *
				       knap->dbl_kappa_help[l]);
			}

			h = (unsigned int) (dq);
			mpz_mul_ui(curr->gmp_help3, curr->gmp_prod, h);
			mpz_sub(knap->gmp_kappa[i][j], 
					curr->gmp_help4, curr->gmp_help3);
			if (mpz_sgn(knap->gmp_kappa[i][j]) < 0 || 
					mpz_cmp(knap->gmp_kappa[i][j], 
						curr->gmp_prod) >= 0) {
				mpz_fdiv_r(knap->gmp_kappa[i][j], 
					   knap->gmp_kappa[i][j],
					   curr->gmp_prod);
			}
		}
	}
	mpz_fdiv_r(knap->gmp_kappa0, curr->gmp_help1, curr->gmp_prod);

	/* check */
	for (i = 0; i < npr_in_p; i++) {
		for (j = 1; j < DEG; j++) {
			mpz_sub(curr->gmp_help2, 
					knap->gmp_D[i][j], knap->gmp_D[i][0]);
			mpz_add(curr->gmp_help1, knap->gmp_D0, curr->gmp_help2);
			mpz_pow_ui(curr->gmp_help3, curr->gmp_help1, DEG);
			mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_a5);
			mpz_sub(curr->gmp_help3, curr->gmp_help3, curr->gmp_N);
			mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help4, 
					curr->gmp_help3, curr->gmp_prod);
			if (mpz_sgn(curr->gmp_help4))
				complain("neu2\n");
			mpz_mul(curr->gmp_help3, curr->gmp_help3, 
					curr->gmp_help1);
			mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_a5);
			mpz_add(curr->gmp_help1, knap->gmp_kappa0, 
					knap->gmp_kappa[i][j]);
			mpz_mul(curr->gmp_help2, curr->gmp_help1, curr->gmp_N);
			mpz_add(curr->gmp_help3, curr->gmp_help3, 
					curr->gmp_help2);
			mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help1, 
					curr->gmp_help3, curr->gmp_prod);
			if (mpz_sgn(curr->gmp_help1)) {
				mpz_out_str(stdout, 10, curr->gmp_prod);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_kappa0);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_kappa[i][j]);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_D0);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_D[i][0]);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_D[i][j]);
				printf("\n");
				mpz_out_str(stdout, 10, curr->gmp_N);
				printf("\n");
				mpz_out_str(stdout, 10, curr->gmp_help1);
				printf("\n");
				mpz_out_str(stdout, 10, curr->gmp_help4);
				printf("\n");
				complain("check at %d %d\n", i, j);
			}
		}
	}

	mpz_pow_ui(curr->gmp_help1, knap->gmp_D0, DEG - 1);
	mpz_mul(curr->gmp_help3, curr->gmp_help1, knap->gmp_D0);
	mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_a5);
	mpz_sub(curr->gmp_help3, curr->gmp_help3, curr->gmp_N);
	mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help4, 
			curr->gmp_help3, curr->gmp_prod);
	if (mpz_sgn(curr->gmp_help4))
		complain("neu4\n");
	mpz_mul(curr->gmp_help2, knap->gmp_kappa0, curr->gmp_help1);
	mpz_add(curr->gmp_help3, curr->gmp_help3, curr->gmp_help2);
	mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help4, 
			curr->gmp_help3, curr->gmp_prod);
	if (mpz_sgn(curr->gmp_help4))
		complain("neu5\n");
	mpz_fdiv_r(curr->gmp_help3, curr->gmp_help3, 
			curr->gmp_help1);	/* modulo D0^(DEG-1) */
	mpz_mul_2exp(curr->gmp_help3, curr->gmp_help3, 64);
	mpz_fdiv_q(curr->gmp_help3, curr->gmp_help3, curr->gmp_help1);
	if (mpz_sizeinbase(curr->gmp_help3, 2) > 64)
		complain("neu6\n");
	mpz_get_ull_64(&lambda0, curr->gmp_help3);

	/* renormalize kappa's to 2^64-range */
	knap->ull_kappa0 = renorm(curr, knap->gmp_kappa0);

	for (i = 0; i < npr_in_p; i++) {
		knap->ull_kappa[i][0] = 0;
		for (j = 1; j < DEG; j++) {
			knap->ull_kappa[i][j] =
				renorm(curr, knap->gmp_kappa[i][j]) + 
				renorm2(curr, knap->gmp_D[i][j]) -
				renorm2(curr, knap->gmp_D[i][0]);
		}
	}
	for (j = 0; j < DEG; j++)
		knap->ull_kappa[0][j] += lambda0;

	db = a3_max / mpz_get_d(curr->gmp_m0);
	if ((db >= 1.) || (db < 0.))
		complain("a_n-2-bound so high that all pols will pass: %f\n",
			 db);
	db *= 18446744073709551616.;
	knap->ull_bound = (uint64_t) db;
	if (verbose > 3)
		printf("ull_bound: %" PRIu64 "\n", knap->ull_bound);

	/*  bound=ull_bound+1;*/

	for (i = knap->len1; i < npr_in_p; i++) {
		for (j = 0; j < DEG; j++) {
			knap->ull_kappa[i][j] = -(knap->ull_kappa[i][j]);
		}
	}
}

/*------------------------------------------------------------------------*/
void 
init_knapsack_exact_p0(curr_poly_t *curr, 
		knapsack_data_t *knap,
		curr_rat_factor_t *curr_factors,
		unsigned int p0, unsigned int r0, double a3_max)
{
	int i, j, l;
	unsigned int p, h, hh, inv;
	uint64_t lambda0;
	unsigned int dp, d0p, p2, hp, N_mod_p2, an_mod_p2;
	double db, dq;
	int npr_in_p = curr_factors->npr_in_p;
	curr_rat_prime_t *primes = curr_factors->primes;
	u32 *p_ind = curr_factors->p_ind;

/* compute D_{i,j} */
	if (mpz_fdiv_q_ui(curr->gmp_help1, curr->gmp_prod, p0))
		complain("init_knapsack_exact.p0.prod\n");
	h = mpz_fdiv_ui(curr->gmp_help1, p0);
	inv = mp_modinv_1(h, p0);
	mpz_mul_ui(knap->gmp_kappa_help[npr_in_p], curr->gmp_help1, inv);
	knap->p_mod_p2[npr_in_p] = mpz_fdiv_ui(curr->gmp_help1, p0 * p0);
	knap->dbl_kappa_help[npr_in_p] = (double) inv / (double) p0;

	p = 1;
	db = (double) h / (double) p0;
	for (i = 0; i < npr_in_p; i++) {
		curr_rat_prime_t *c = primes + p_ind[i];
		p = c->p_pr;
		if (mpz_fdiv_q_ui(curr->gmp_help1, curr->gmp_prod, p))
			complain("init_knapsack_exact.prod\n");
		h = mpz_fdiv_ui(curr->gmp_help1, p);
		inv = mp_modinv_1(h, p);
		mpz_mul_ui(knap->gmp_kappa_help[i], curr->gmp_help1, inv);
		knap->p_mod_p2[i] = mpz_fdiv_ui(curr->gmp_help1, p * p);

		knap->dbl_kappa_help[i] = (double) inv / (double) p;

		h = inv * c->p_fr[0];
		h %= p;
		knap->ul_d_help[i][0] = h;
		mpz_mul_ui(knap->gmp_D[i][0], curr->gmp_help1, h);
		db += (double) h / (double) p;

		for (j = 1; j < DEG; j++) {
			h = inv * c->p_fr[j];
			h %= p;
			h += (p - knap->ul_d_help[i][0]);
			if (h >= p)
				h -= p;
			knap->ul_d_help[i][j] = h;
			h += knap->ul_d_help[i][0];
			mpz_mul_ui(knap->gmp_D[i][j], curr->gmp_help1, h);
		}
	}
	db += 0.5 * (double) npr_in_p;

	if (mpz_fdiv_q_ui(curr->gmp_help1, curr->gmp_prod, p0))
		complain("search.prod.p0\n");
	h = mpz_fdiv_ui(curr->gmp_help1, p0);
	h = mp_modinv_1(h, p0);
	h = (h * r0) % p0;
	mpz_mul_ui(knap->gmp_D0, curr->gmp_help1, h);
	mpz_set(curr->gmp_help2, knap->gmp_D0);
	for (i = 0; i < npr_in_p; i++)
		mpz_add(knap->gmp_D0, knap->gmp_D0, knap->gmp_D[i][0]);
	mpz_fdiv_r(curr->gmp_help1, curr->gmp_m0, curr->gmp_prod);
	mpz_sub(knap->gmp_disp, curr->gmp_m0, curr->gmp_help1);
	mpz_mul_si(curr->gmp_help1, curr->gmp_prod, (int) db);
	mpz_sub(knap->gmp_disp, knap->gmp_disp, curr->gmp_help1);
	mpz_add(knap->gmp_D0, knap->gmp_D0, knap->gmp_disp);
	mpz_add(knap->gmp_disp, knap->gmp_disp, curr->gmp_help2);

	/* compute kappa_0 and kappa_{i,j} */
	p = p0;
	p2 = p * p;
	d0p = mpz_fdiv_ui(knap->gmp_D0, p2);
	an_mod_p2 = mpz_fdiv_ui(curr->gmp_a5, p2);
	N_mod_p2 = mpz_fdiv_ui(curr->gmp_N, p2);

	hp = mp_expo_1(d0p, DEG, p2);
	hp = mp_modmul_1(hp, an_mod_p2, p2);
	hp = p2 - hp + N_mod_p2;
	if (hp % p)
		complain("%u %u %u neu-1re\n", p, d0p, hp);
	h = hp / p;
	inv = N_mod_p2;
	for (j = 0; j < npr_in_p; j++)
		inv = mp_modmul_1(inv, primes[p_ind[j]].p_pr, p);
	inv = mp_modinv_1(inv, p);
	h = mp_modmul_1(h, inv, p);
	h = mp_modmul_1(h, an_mod_p2, p);
	h = mp_modmul_1(h, d0p, p);
	mpz_mul_ui(curr->gmp_help1, knap->gmp_kappa_help[npr_in_p], h);

	for (i = 0; i < npr_in_p; i++) {
		curr_rat_prime_t *c = primes + p_ind[i];

		mpz_set_ui(knap->gmp_kappa[i][0], 0);
		p = c->p_pr;
		p2 = p * p;
		d0p = mpz_fdiv_ui(knap->gmp_D0, p2);
		an_mod_p2 = c->p_a5_mod_p2;
		N_mod_p2 = c->p_N_mod_p2;

		hp = mp_expo_1(d0p, DEG, p2);
		hp = mp_modmul_1(hp, an_mod_p2, p2);
		hp = p2 - hp + N_mod_p2;
		if (hp % p)
			complain("%u %u %u neu-1re\n", p, d0p, hp);
		h = hp / p;
		inv = mp_modmul_1(an_mod_p2, c->p_N_inv, p);
		for (j = 0; j < i; j++) {
			inv = mp_modmul_1(inv, curr_factors->p_inv_table[p_ind[j]][p_ind[i]], p);
		}
		for (j = i + 1; j < npr_in_p; j++) {
			inv = mp_modmul_1(inv, curr_factors->p_inv_table[p_ind[j]][p_ind[i]], p);
		}
		inv = mp_modmul_1(inv, mp_modinv_1(p0 % p, p), p);

		h = mp_modmul_1(h, inv, p);
		h = mp_modmul_1(h, d0p, p);
		mpz_mul_ui(curr->gmp_help2, knap->gmp_kappa_help[i], h);
		mpz_add(curr->gmp_help1, curr->gmp_help1, curr->gmp_help2);

		for (j = 1; j < DEG; j++) {
			dp = mp_modmul_1(knap->p_mod_p2[i], 
					knap->ul_d_help[i][j], p2);
			hh = d0p + dp;
			h = mp_expo_1(hh, DEG + 1, p2);
			hh = mp_expo_1(d0p, DEG + 1, p2);
			h = mp_modmul_1(hh + (p2 - h), an_mod_p2, p2);
			h += mp_modmul_1(N_mod_p2, dp, p2);

			if (h % p)
				complain("%u %u %u neu1\n", p, d0p, dp);
			h /= p;
			h = mp_modmul_1(h, inv, p);
			mpz_mul_ui(curr->gmp_help4, knap->gmp_kappa_help[i], h);
			dq = (double) h * knap->dbl_kappa_help[i];

			h = mpz_fdiv_ui(curr->gmp_a5, p0);
			h = mp_modmul_1(h, 5 * p0 - 5, p0);
			h = mp_modmul_1(h, knap->ul_d_help[i][j], p0);
			h = mp_modmul_1(h, mp_modinv_1(p % p0, p0), p0);
			mpz_mul_ui(curr->gmp_help3, 
					knap->gmp_kappa_help[npr_in_p], h);
			mpz_add(curr->gmp_help4, curr->gmp_help4, 
					curr->gmp_help3);
			dq += (((double) h) * knap->dbl_kappa_help[npr_in_p]);

			for (l = 0; l < npr_in_p; l++) {
				curr_rat_prime_t *cl = primes + p_ind[l];

				if (l == i)
					continue;
				h = mp_modmul_1(cl->p_minus5a5_mod_p * 
						knap->ul_d_help[i][j],
						curr_factors->p_inv_table[p_ind[i]][p_ind[l]],
						cl->p_pr);
				mpz_mul_ui(curr->gmp_help3, 
					   knap->gmp_kappa_help[l], h);
				mpz_add(curr->gmp_help4, 
					curr->gmp_help4,
					curr->gmp_help3);
				dq += (((double) h) *
				       knap->dbl_kappa_help[l]);
			}

			h = (unsigned int) (dq);
			mpz_mul_ui(curr->gmp_help3, curr->gmp_prod, h);
			mpz_sub(knap->gmp_kappa[i][j], 
				curr->gmp_help4, curr->gmp_help3);
			if (mpz_sgn(knap->gmp_kappa[i][j]) < 0 ||
					mpz_cmp(knap->gmp_kappa[i][j], 
					curr->gmp_prod) >= 0) {
				mpz_fdiv_r(knap->gmp_kappa[i][j], 
					   knap->gmp_kappa[i][j],
					   curr->gmp_prod);
			}
		}
	}
	mpz_fdiv_r(knap->gmp_kappa0, curr->gmp_help1, curr->gmp_prod);

	/* check */
	for (i = 0; i < npr_in_p; i++) {
		for (j = 1; j < DEG; j++) {
			mpz_sub(curr->gmp_help2, knap->gmp_D[i][j], 
					knap->gmp_D[i][0]);
			mpz_add(curr->gmp_help1, knap->gmp_D0, curr->gmp_help2);
			mpz_pow_ui(curr->gmp_help3, curr->gmp_help1, DEG);
			mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_a5);
			mpz_sub(curr->gmp_help3, curr->gmp_help3, curr->gmp_N);
			mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help4, 
					curr->gmp_help3, curr->gmp_prod);
			if (mpz_sgn(curr->gmp_help4))
				complain("neu2\n");
			mpz_mul(curr->gmp_help3, curr->gmp_help3, 
					curr->gmp_help1);
			mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_a5);
			mpz_add(curr->gmp_help1, knap->gmp_kappa0, 
					knap->gmp_kappa[i][j]);
			mpz_mul(curr->gmp_help2, curr->gmp_help1, curr->gmp_N);
			mpz_add(curr->gmp_help3, curr->gmp_help3, 
					curr->gmp_help2);
			mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help1, 
					curr->gmp_help3, curr->gmp_prod);
			if (mpz_sgn(curr->gmp_help1)) {
				mpz_out_str(stdout, 10, curr->gmp_a5);
				printf("\n");
				mpz_out_str(stdout, 10, curr->gmp_prod);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_kappa0);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_kappa[i][j]);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_D0);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_D[i][0]);
				printf("\n");
				mpz_out_str(stdout, 10, knap->gmp_D[i][j]);
				printf("\n");
				mpz_out_str(stdout, 10, curr->gmp_N);
				printf("\n");
				mpz_out_str(stdout, 10, curr->gmp_help1);
				printf("\n");
				mpz_out_str(stdout, 10, curr->gmp_help4);
				printf("\n");
				complain("check at %d %d   p0: %u\n", i, j, p0);
			}
		}
	}

	mpz_pow_ui(curr->gmp_help1, knap->gmp_D0, DEG - 1);
	mpz_mul(curr->gmp_help3, curr->gmp_help1, knap->gmp_D0);
	mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_a5);
	mpz_sub(curr->gmp_help3, curr->gmp_help3, curr->gmp_N);
	mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help4, 
			curr->gmp_help3, curr->gmp_prod);
	if (mpz_sgn(curr->gmp_help4))
		complain("neu4\n");
	mpz_mul(curr->gmp_help2, knap->gmp_kappa0, curr->gmp_help1);
	mpz_add(curr->gmp_help3, curr->gmp_help3, curr->gmp_help2);
	mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help4, 
			curr->gmp_help3, curr->gmp_prod);
	if (mpz_sgn(curr->gmp_help4))
		complain("neu5\n");
	mpz_fdiv_r(curr->gmp_help3, curr->gmp_help3, curr->gmp_help1);	/* modulo D0^(DEG-1) */
	mpz_mul_2exp(curr->gmp_help3, curr->gmp_help3, 64);
	mpz_fdiv_q(curr->gmp_help3, curr->gmp_help3, curr->gmp_help1);
	if (mpz_sizeinbase(curr->gmp_help3, 2) > 64)
		complain("neu6\n");
	mpz_get_ull_64(&lambda0, curr->gmp_help3);

	/* renormalize kappa's to 2^64-range */
	knap->ull_kappa0 = renorm(curr, knap->gmp_kappa0);

	for (i = 0; i < npr_in_p; i++) {
		knap->ull_kappa[i][0] = 0;
		for (j = 1; j < DEG; j++) {
			knap->ull_kappa[i][j] =
				renorm(curr, knap->gmp_kappa[i][j]) + 
				renorm2(curr, knap->gmp_D[i][j]) -
				renorm2(curr, knap->gmp_D[i][0]);
		}
	}
	for (j = 0; j < DEG; j++) {
		knap->ull_kappa[0][j] += lambda0;
	}

	db = a3_max / mpz_get_d(curr->gmp_m0);
	if ((db >= 1.) || (db < 0.)) {
		complain("a_n-2-bound so high that all pols will pass: %f\n",
			 db);
	}
	db *= 18446744073709551616.;
	knap->ull_bound = (uint64_t) db;
	if (verbose > 3)
		printf("ull_bound: %" PRIu64 "\n", knap->ull_bound);

	/*  bound=ull_bound+1;*/

	for (i = knap->len1; i < npr_in_p; i++)
		for (j = 0; j < DEG; j++)
			knap->ull_kappa[i][j] = -(knap->ull_kappa[i][j]);
}

/*------------------------------------------------------------------------*/
void
knapsack_exact(curr_poly_t *curr, knapsack_data_t *knap,
		bounds_t *bounds, stage1_stat_t *stats)
{
	int i, i1, i2;
	uint64_t bound, v1, v2;
	uint64_t *s1 = knap->s1;
	uint64_t *s2 = knap->s2;
	uint64_t *s1sort = knap->s1sort;
	uint64_t *s2sort = knap->s2sort;
	int len1 = knap->len1;
	int len2 = knap->len2;
	int s1len = knap->s1len;
	int s2len = knap->s2len;

	combine(knap, s1, s1sort, 0, len1);
	combine(knap, s2, s2sort, len1, len1 + len2);

	bound = knap->ull_bound + 1;
	knap->nstore = 0;
	if ((s1sort[0] - s2sort[s2len - 1]) < bound)
		store(knap, 0, s2len - 1);
	if ((-s1sort[s1len - 1] + s2sort[0]) < bound)
		store(knap, s1len - 1, 0);
	i1 = 0;
	i2 = 0;
	while (1) {
		v1 = s1sort[i1];
		v2 = s2sort[i2];
		if (v1 < v2) {
			if ((v2 - v1) < bound) {
				store(knap, i1, i2);
				for (i = i2 + 1; i < s2len; i++) {
					if ((s2sort[i] - v1) < bound)
						store(knap, i1, i);
					else
						break;
				}
			}
			if (i1 < s1len - 1)
				i1++;
			else
				break;
		}
		else {
			if ((v1 - v2) < bound) {
				store(knap, i1, i2);
				for (i = i1 + 1; i < s1len; i++) {
					if ((s1sort[i] - v2) < bound)
						store(knap, i, i2);
					else
						break;
				}
			}
			if (i2 < s2len - 1)
				i2++;
			else
				break;
		}
	}
	if (knap->nstore) {
		profile_start(PROF_KNAP_EXACT_CHECK);
		check_stored_pairs(curr, knap, bounds, stats);
		profile_stop(PROF_KNAP_EXACT_CHECK);
	}
}
