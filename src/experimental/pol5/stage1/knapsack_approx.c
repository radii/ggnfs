#include "stage1_impl.h"

/*------------------------------------------------------------------------*/
static void
combine_raw0(knapsack_data_t *knap, int len, int ind)
{
	int i, j, disp;
	uint64_t add;
	uint64_t *shelp = knap->shelp;

	if (knap->raw_ull_kappa[ind][0])
		complain("combine_raw0\n");
	disp = len;
	for (j = 1; j < 5; j++) {
		add = knap->raw_ull_kappa[ind][j];
		for (i = 0; i < len; i++)
			shelp[i + disp] = add + shelp[i];
		disp += len;
	}
}

/*------------------------------------------------------------------------*/
static void
combine_raw0last(knapsack_data_t *knap, 
		unsigned int *targ, int len, int ind)
{
	int i, j, disp;
	uint64_t add;
	uint64_t *shelp = knap->shelp;

	if (knap->raw_ull_kappa[ind][0])
		complain("combine_raw0last\n");
	for (i = 0; i < len; i++)
		targ[i] = (unsigned int) (shelp[i] >> 32);
	disp = len;
	for (j = 1; j < 5; j++) {
		add = knap->raw_ull_kappa[ind][j];
		for (i = 0; i < len; i++)
			targ[i + disp] =
				(unsigned int) ((add + shelp[i]) >> 32);
		disp += len;
	}
}

/*------------------------------------------------------------------------*/
static void
combine_raw(knapsack_data_t *knap, 
		unsigned int *targ, int i0, int i1)
{
	int i, len;

	if (i1 - i0 < 1)
		complain("combine_raw\n");
	if (i1 - i0 == 1) {
		for (i = 0; i < 5; i++) {
			targ[i] = (unsigned int) 
				(knap->raw_ull_kappa[i0][i] >> 32);
		}
		return;
	}
	for (i = 0; i < 5; i++)
		knap->shelp[i] = knap->raw_ull_kappa[i0][i];
	len = 5;
	for (i = i0 + 1; i < i1 - 1; i++) {
		combine_raw0(knap, len, i);
		len *= 5;
	}
	combine_raw0last(knap, targ, len, i1 - 1);
}

/*------------------------------------------------------------------------*/
static void
raw_store(knapsack_data_t *knap, int i1, int i2)
{
	while (knap->nraw_store + 1 >= knap->raw_store_len) {
		knap->raw_store_len += 16;
		knap->raw_stored_pairs =
			(uint64_t *) xrealloc(knap->raw_stored_pairs,
					      2 * knap->raw_store_len *
					      sizeof(uint64_t));
	}
	knap->raw_stored_pairs[2 * knap->nraw_store] = i1;
	knap->raw_stored_pairs[2 * knap->nraw_store + 1] = i2;
	knap->nraw_store++;
}

/*------------------------------------------------------------------------*/
#ifndef HAVE_ASM

static int
raw_hash_1(knapsack_data_t *knap)
{
	unsigned int i, j;
	unsigned int ind, h;
	unsigned int add, *hash;
	unsigned char *sort;
	unsigned int *hashdata = knap->hashdata;
	unsigned int s11len = knap->s11len;
	unsigned int *s11l = knap->s11l;
	unsigned int s12len = knap->s12len;
	unsigned int *s12l = knap->s12l;

	memset(hashdata, 0, NHASH * sizeof(unsigned char));
	sort = (unsigned char *) hashdata;
	hash = hashdata + (NHASH >> 2);
	for (i = 0; i < s11len; i++) {
		add = s11l[i];
		for (j = 0; j < s12len; j++) {
			h = add + s12l[j];
			ind = h >> HASHSHIFT32;
			if (sort[ind] >= 32)
				return 1;
			hash[sort[ind] * NHASH + ind] = h;
			sort[ind]++;
		}
	}
	return 0;
}

static void
raw_hash_2(knapsack_data_t *knap)
{
	unsigned int i, j;
	int k;
	unsigned int ind, h;
	unsigned int add, *hash;
	unsigned char *sort;
	unsigned int *hashdata = knap->hashdata;
	unsigned int s21len = knap->s21len;
	unsigned int *s21l = knap->s21l;
	unsigned int s22len = knap->s22len;
	unsigned int *s22l = knap->s22l;
	unsigned int raw_bound = knap->raw_bound;

	sort = (unsigned char *) hashdata;
	hash = hashdata + (NHASH >> 2);
	for (i = 0; i < s21len; i++) {
		add = s21l[i];
		for (j = 0; j < s22len; j++) {
			h = add + s22l[j];
			ind = h >> HASHSHIFT32;
			for (k = 0; k < sort[ind]; k++) {
				if (hash[k * NHASH + ind] - h < raw_bound) {
					/* assumes npr_in_p<25 */
					*(knap->raw_cand_ptr) = (i + (j << 16));
					knap->raw_cand_ptr++;
					k = sort[ind];	/* avoid duplicates */
				}
			}
		}
	}
}

#endif /* !HAVE_ASM */

/*------------------------------------------------------------------------*/
static int
raw_hash_3(knapsack_data_t *knap)
{
	int k, i, j;
	unsigned int h, ind, *hash;
	unsigned char *sort;
	int nraw_cand;
	unsigned int *hashdata = knap->hashdata;
	unsigned int *s21l = knap->s21l;
	unsigned int *s22l = knap->s22l;
	int *raw_cand = knap->raw_cand;
	unsigned int *raw_cand_hash = knap->raw_cand_hash;

	nraw_cand = knap->nraw_cand = knap->raw_cand_ptr - knap->raw_cand;
	memset(hashdata, 0, NHASH * sizeof(unsigned char));
	sort = (unsigned char *) hashdata;
	hash = hashdata + (NHASH >> 2);
	for (k = 0; k < nraw_cand; k++) {
		i = raw_cand[k] & 0x0000ffff;
		j = raw_cand[k] >> 16;
		h = s21l[i] + s22l[j];
		raw_cand_hash[k] = h;
		ind = h >> HASHSHIFT32;
		if (sort[ind] >= 32)
			return 1;
		hash[sort[ind] * NHASH + ind] = h;
		sort[ind]++;
	}
	return 0;
}

/*------------------------------------------------------------------------*/
static void
raw_hash_4(knapsack_data_t *knap)
{
	int i, j, k, l, m, i1, i2;
	unsigned int ind, h, h1, hh;
	unsigned int add, *hash;
	unsigned char *sort;
	uint64_t sum1, sum2;
	unsigned int *hashdata = knap->hashdata;
	unsigned int len1 = knap->len1;
	unsigned int len2 = knap->len2;
	unsigned int s21len = knap->s21len;
	unsigned int *s11l = knap->s11l;
	unsigned int s11len = knap->s11len;
	unsigned int *s12l = knap->s12l;
	unsigned int s12len = knap->s12len;
	unsigned int raw_bound = knap->raw_bound;
	unsigned int *raw_cand_hash = knap->raw_cand_hash;
	int *raw_cand = knap->raw_cand;
	int nraw_cand = knap->nraw_cand;

	sort = (unsigned char *) hashdata;
	hash = hashdata + (NHASH >> 2);
	for (i = 0; i < s11len; i++) {
		add = s11l[i];
		for (j = 0; j < s12len; j++) {
			h = add + s12l[j];
			ind = h >> HASHSHIFT32;
			for (k = 0; k < sort[ind]; k++) {
				if ((h - hash[k * NHASH + ind]) >= raw_bound)
					continue;

				h1 = hash[k * NHASH + ind];
				for (l = 0; l < nraw_cand; l++) {
					if (h1 != raw_cand_hash[l])
						continue;

					i1 = i + j * s11len;
					i2 = (raw_cand[l] & 0xffff) + 
						(raw_cand[l] >> 16) * s21len;
					sum1 = 0;
					sum2 = 0;
					for (m = 0; m < len1; m++) {
						sum1 += knap->raw_ull_kappa[m][i1 % 5];
						i1 -= (i1 % 5);
						i1 /= 5;
					}
					for (m = 0; m < len2; m++) {
						sum2 += knap->raw_ull_kappa[len1 + m][i2 % 5];
						i2 -= (i2 % 5);
						i2 /= 5;
					}
#if 1
					/* check */
					hh = (unsigned int)((sum1 +
							(8589934592ULL +
							knap->raw_ull_bound)) >> 32);
					if (hh > h) {
						if (hh - h > 3)
							complain("hash4.0\n");
					}
					else {
						if (h - hh > 3)
							complain("%u %u %u %u %" PRIu64 " hash4.1\n", h, hh, h1, raw_bound, knap->raw_ull_bound);
					}
					hh = (unsigned int)(sum2 >> 32);
					if (hh > h1) {
						if (hh - h1 > 3)
							complain("%u %u %u %u hash4.2\n", h, hh, h1, raw_bound);
					}
					else if (h1 - hh > 3) {
						complain("%u %u %u %u %u %u hash4.3\n", h, hh, h1, raw_bound, l, nraw_cand);
					}
#endif
					if (((sum1 - sum2) < knap->raw_ull_bound) ||
					    ((sum2 - sum1) < knap->raw_ull_bound)) {
						i1 = i + j * s11len;
						i2 = (raw_cand[l] & 0xffff) +
							(raw_cand[l] >> 16) * 
							s21len;
						raw_store(knap, i1, i2);
					}
				}
			}
		}
	}
}

/*------------------------------------------------------------------------*/
void
init_knapsack_raw(curr_poly_t *curr, 
		knapsack_data_t *knap,
		curr_rat_factor_t *curr_factors,
		double a3_max, stage1_stat_t *stats)
{
	int i, j, l;
	unsigned int p, h, hh, inv, h6;
	unsigned int dp, d0p, p2, hp, N_mod_p2, a5_mod_p2;
	double db;
	int npr_in_p = curr_factors->npr_in_p;
	u32 *p_ind = curr_factors->p_ind;
	curr_rat_prime_t *primes = curr_factors->primes;

#ifdef HAVE_FLOAT64
	long double dbl_prod;
	long double dbl_5a5p;
#else
	uint64_t ull_5a5p;
#endif

	/* preparation */
#ifdef HAVE_FLOAT64
	dbl_prod = mpz_get_ld(curr->gmp_prod);
	dbl_5a5p = (5.L * mpz_get_ld(curr->gmp_a5)) / dbl_prod;
#else
	mpz_mul_ui(curr->gmp_help1, curr->gmp_a5, 5);
	mpz_mul_2exp(curr->gmp_help1, curr->gmp_help1, 64);
	mpz_fdiv_q(curr->gmp_help4, curr->gmp_help1, curr->gmp_prod);
	if (mpz_sizeinbase(curr->gmp_help4, 2) > 64)
		complain("ull5a5p-comp\n");
	mpz_get_ull_64(&ull_5a5p, curr->gmp_help4);
#endif

	for (i = 0; i < npr_in_p; i++) {
		for (j = 0; j < npr_in_p; j++) {
			knap->prep_inv_table[i][j] = 
				curr_factors->p_inv_table[p_ind[i]][p_ind[j]];
		}
	}

	profile_start(PROF_CRT_AUX);
	/* computation of D_{i,j} */
	db = 0.;
	for (i = 0; i < npr_in_p; i++) {
		curr_rat_prime_t *c = primes + p_ind[i];
		p = c->p_pr;
		knap->prep_p[i] = p;
		knap->prep_5a5[i] = c->p_minus5a5_mod_p;
		knap->prep_N_mod_p2[i] = c->p_N_mod_p2;

		if (mpz_fdiv_q_ui(curr->gmp_help1, curr->gmp_prod, p))
			complain("search.prod raw\n");
		knap->p_mod_p2[i] = mpz_fdiv_ui(curr->gmp_help1, p * p);
		inv = mp_modinv_1(knap->p_mod_p2[i] % p, p);
		knap->ull_kappa_p_inv[i] = c->ull_p_inv;

		knap->ull_kappa_help0[i] = mp_modmul_1(inv, inv, p);
		knap->prep_5a5[i] = mp_modmul_1(knap->prep_5a5[i], inv, p);
		h = mp_modmul_1(inv, c->p_fr[0], p);
		knap->ul_d_help[i][0] = h;
		mpz_mul_ui(knap->gmp_D[i][0], curr->gmp_help1, h);
		db += (double) h / (double) p;

		for (j = 1; j < 5; j++) {
			knap->ul_d_help[i][j] = mp_modmul_1(inv,
					 c->p_fr[j] + (p - c->p_fr[0]), p);
		}
	}

	for (l = 0; l < npr_in_p; l++) {
		for (i = 0; i < npr_in_p; i++) {
			knap->prep_inv_table[i][l] =
				knap->prep_5a5[l] * knap->prep_inv_table[i][l];
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

	profile_stop(PROF_CRT_AUX);
	profile_start(PROF_KAPPA);
	for (i = 0; i < npr_in_p; i++) {
		for (j = 0; j < 5; j++) {
			knap->raw_ull_kappa[i][j] = 0;
		}
	}

	/* compute kappa_0 and kappa_{i,j} */
	knap->ull_kappa0 = 0;
	for (i = 0; i < npr_in_p; i++) {
		curr_rat_prime_t *c = primes + p_ind[i];
		p = c->p_pr;
		p2 = p * p;
		d0p = mpz_fdiv_ui(knap->gmp_D0, p2);
		a5_mod_p2 = c->p_a5_mod_p2;
		N_mod_p2 = c->p_N_mod_p2;	/* ca. 270-290 Tz */

		hp = mp_modmul_1(d0p, d0p, p2);
		hp = mp_modmul_1(hp, hp, p2);
		hp = mp_modmul_1(hp, d0p, p2);
		hp = mp_modmul_1(hp, a5_mod_p2, p2);
		hp = p2 - hp + N_mod_p2;
		if (hp % p)
			complain("%u %u %u neu-1\n", p, d0p, hp);
		h = hp / p;
		inv = mp_modmul_1(a5_mod_p2, c->p_N_inv * 
					knap->ull_kappa_help0[i], p);
		h = mp_modmul_1(h, inv * d0p, p);
		ulladdmul(&knap->ull_kappa0, h, knap->ull_kappa_p_inv + i);

		hh = mp_modmul_1(d0p, d0p, p2);
		hh = mp_modmul_1(hh, d0p, p2);
		h6 = mp_modmul_1(hh, hh, p2);	/* h6=d0p^6 */

		profile_start(PROF_KAPPA2);
		for (j = 1; j < 5; j++) {
			dp = mp_modmul_1(knap->p_mod_p2[i], 
					knap->ul_d_help[i][j], p2);
			hh = d0p + dp;
			h = mp_modmul_1(hh, hh, p2);
			h = mp_modmul_1(hh, h, p2);
			h = mp_modmul_1(h, h, p2);	/* h=(d0p+dp)^6 */
			h = mp_modmul_1(h6 + (p2 - h), a5_mod_p2, p2);
			h += mp_modmul_1(N_mod_p2, dp, p2);

			if (h % p)
				complain("%u %u %u neu1\n", p, d0p, dp);
			h /= p;
			h = mp_modmul_1(h, inv, p);

			ulladdmul(&(knap->raw_ull_kappa[i][j]), h,
				  knap->ull_kappa_p_inv + i);
		}
		for (l = 0; l < i; l++) {
			for (j = 1; j < 5; j++) {
				h = mp_modmul_1(knap->ul_d_help[i][j],
					     knap->prep_inv_table[i][l],
					     knap->prep_p[l]);
				ulladdmul(&(knap->raw_ull_kappa[i][j]), h,
					  knap->ull_kappa_p_inv + l);
			}
		}
		for (l = i + 1; l < npr_in_p; l++) {
			for (j = 1; j < 5; j++) {
				h = mp_modmul_1(knap->ul_d_help[i][j],
					     knap->prep_inv_table[i][l],
					     knap->prep_p[l]);
				ulladdmul(&(knap->raw_ull_kappa[i][j]), h,
					  knap->ull_kappa_p_inv + l);
			}
		}
		profile_stop(PROF_KAPPA2);
	}

	profile_stop(PROF_KAPPA);
	profile_start(PROF_COMPUTE_KNAP_VALS);

#ifdef HAVE_FLOAT64
	{
		/* compute kappa0/prod - 5*a5*(D0-root)/prod^2 + 
		  		10*a5*(D0-root)^2/prod^2/D0 */
		uint64_t ullh;
		long double ddd, dde;

		mpz_sub(curr->gmp_help2, curr->gmp_root, knap->gmp_D0);
		dde = -mpz_get_ld(curr->gmp_help2);

		ddd = dde / dbl_prod;
		ddd *= dbl_5a5p;
		ullh = (uint64_t) (18446744073709551616.L *
				   (ddd - floorl(ddd)));
		ddd *= 2. * dde;
		ddd /= mpz_get_ld(curr->gmp_root);
		ullh -= (uint64_t) (18446744073709551616.L *
				    (ddd - floorl(ddd)));
		knap->lambda0 = ullh + knap->ull_kappa0;
	}
#else
	{
		uint64_t ullh1, ullh2;

		mpz_sub(curr->gmp_help2, knap->gmp_D0, curr->gmp_root);	/* D0-root */
		mpz_mul_ui(curr->gmp_help1, curr->gmp_a5, 5);	/* 5*a5 */
		mpz_mul(curr->gmp_help1, curr->gmp_help1, curr->gmp_help2);/* 5*a5*(D0-root) */
		mpz_mul(curr->gmp_help3, curr->gmp_prod, curr->gmp_prod);	/* prod^2 */
		mpz_mul_2exp(curr->gmp_help1, curr->gmp_help1, 64);
		mpz_fdiv_q(curr->gmp_help4, curr->gmp_help1, curr->gmp_help3);
		if (mpz_sizeinbase(curr->gmp_help4, 2) > 64)
			complain("lambda-comp\n");
		mpz_get_ull_64(&ullh1, curr->gmp_help4);

		mpz_add(curr->gmp_help2, curr->gmp_help2, curr->gmp_help2);	/* 10*a5*(D0-root)^2 * 2^64 */
		mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_help2);
		mpz_fdiv_q(curr->gmp_help4, curr->gmp_help4, knap->gmp_D0);

		if (mpz_sizeinbase(curr->gmp_help4, 2) > 64)
			complain("lambda-comp\n");
		mpz_get_ull_64(&ullh2, curr->gmp_help4);
		knap->lambda0 = -ullh1 + ullh2 + knap->ull_kappa0;
	}
#endif
	profile_stop(PROF_COMPUTE_KNAP_VALS);

	{
		uint64_t ullh;

		for (i = 0; i < npr_in_p; i++) {
			curr_rat_prime_t *c = primes + p_ind[i];
#ifdef HAVE_FLOAT64
			ullh = dull(dbl_5a5p * c->ld_p_inv);
#else
			ull_mulh(&ullh, &(c->ull_p_inv), &(ull_5a5p));
#endif
			for (j = 1; j < 5; j++) {
				ulladdmul(&(knap->raw_ull_kappa[i][j]),
					  knap->ul_d_help[i][j], &ullh);
			}
		}
		for (j = 0; j < 5; j++)
			knap->raw_ull_kappa[0][j] += knap->lambda0;
	}

	for (i = knap->len1; i < npr_in_p; i++) {
		for (j = 0; j < 5; j++) {
			knap->raw_ull_kappa[i][j] = 
				-(knap->raw_ull_kappa[i][j]);
		}
	}

	db = a3_max / mpz_get_d(curr->gmp_m0);
	if ((db >= 1.) || (db < 0.))
		complain("a3-bound so high that all pols will pass: %f\n", db);
	db *= 18446744073709551616.;
	knap->raw_ull_bound = (uint64_t) db;
	knap->raw_ull_bound += (uint64_t) ((npr_in_p + 2) * npr_in_p * 
				primes[p_ind[npr_in_p - 1]].p_pr);
	if (verbose > 3)
		printf("raw_ull_bound: %" PRIu64 "\n", knap->raw_ull_bound);
	knap->raw_bound = (unsigned int) (knap->raw_ull_bound >> 31) + 1 + 4;
}

/*------------------------------------------------------------------------*/
void
init_knapsack_raw_p0(curr_poly_t *curr, 
			knapsack_data_t *knap,
			curr_rat_factor_t *curr_factors,
			unsigned int p0, unsigned int r0, 
			double a3_max, stage1_stat_t *stats)
{
	int i, j, l;
	unsigned int p, h, hh, h6, inv;
	unsigned int dp, d0p, p2, hp, N_mod_p2, a5_mod_p2;
	double db;
	int npr_in_p = curr_factors->npr_in_p;
	u32 *p_ind = curr_factors->p_ind;
	curr_rat_prime_t *primes = curr_factors->primes;

#ifdef HAVE_FLOAT64
	long double dbl_prod;
	long double dbl_5a5p;
#else
	uint64_t ull_5a5p;
#endif
	unsigned int p_minus5a5_mod_p0;

	if (curr_factors->last_p0 == 1) {
		for (i = 0; i < npr_in_p; i++) {
			for (j = 1; j < 5; j++) {
				knap->ul_d_help_1[i][j] = knap->ul_d_help[i][j];
			}
		}
	}

	if (curr_factors->last_p0 != p0) {
		for (i = 0; i < npr_in_p; i++) {
			curr_rat_prime_t *c = primes + p_ind[i];
			knap->p0_inv_p[i] = mp_modinv_1(p0 % c->p_pr, c->p_pr);
			knap->p_inv_p0[i] = mp_modinv_1(c->p_pr % p0, p0);
		}
	}
	knap->prep_p[npr_in_p] = p0;
	p_minus5a5_mod_p0 = mp_modmul_1(mpz_fdiv_ui(curr->gmp_a5, p0), 
					5 * p0 - 5, p0);
	knap->ull_kappa_p_inv[npr_in_p] = 18446744073709551615ULL / 
					((uint64_t) p0);

	for (i = 0; i < npr_in_p; i++) {
		for (j = 0; j < npr_in_p; j++) {
			knap->prep_inv_table[i][j] = 
				curr_factors->p_inv_table[p_ind[i]][p_ind[j]];
		}
	}
	for (i = 0; i < npr_in_p; i++) {
		knap->prep_inv_table[i][npr_in_p] = knap->p_inv_p0[i];
		knap->prep_inv_table[npr_in_p][i] = knap->p0_inv_p[i];
	}

	for (i = 0; i < npr_in_p; i++) {
		curr_rat_prime_t *c = primes + p_ind[i];
		knap->prep_5a5[i] = c->p_minus5a5_mod_p;
		knap->prep_N_mod_p2[i] = c->p_N_mod_p2;
	}
	profile_start(PROF_CRT_AUX);

	/* compute D_{i,j} */
	p = 1;

#ifdef HAVE_FLOAT64
	dbl_prod = mpz_get_ld(curr->gmp_prod);
	dbl_5a5p = (5.L * mpz_get_ld(curr->gmp_a5)) / dbl_prod;
#else
	mpz_mul_ui(curr->gmp_help1, curr->gmp_a5, 5);
	mpz_mul_2exp(curr->gmp_help1, curr->gmp_help1, 64);
	mpz_fdiv_q(curr->gmp_help4, curr->gmp_help1, curr->gmp_prod);
	if (mpz_sizeinbase(curr->gmp_help4, 2) > 64)
		complain("ull5a5p-comp\n");
	mpz_get_ull_64(&ull_5a5p, curr->gmp_help4);
#endif

	h = mpz_fdiv_ui(curr->gmp_prod_base, p0);
	inv = mp_modinv_1(h, p0);
	knap->p_mod_p2[npr_in_p] = mpz_fdiv_ui(curr->gmp_prod_base, p0 * p0);
	knap->ull_kappa_help0[npr_in_p] = mp_modmul_1(inv, inv, p0);
	knap->prep_5a5[npr_in_p] = mp_modmul_1(p_minus5a5_mod_p0, inv, p0);
	h = mp_modmul_1(inv, r0, p0);
	mpz_mul_ui(knap->gmp_D0, curr->gmp_prod_base, h);
	db = (double) h / (double) p0;

	for (i = 0; i < npr_in_p; i++) {
		curr_rat_prime_t *c = primes + p_ind[i];
		p = c->p_pr;
		if (mpz_fdiv_q_ui(curr->gmp_help1, curr->gmp_prod, p))
			complain("search.prod raw p0\n");
		knap->p_mod_p2[i] = mpz_fdiv_ui(curr->gmp_help1, p * p);
		inv = mp_modinv_1(knap->p_mod_p2[i] % p, p);

		knap->ull_kappa_help0[i] = mp_modmul_1(inv, inv, p);
		knap->prep_5a5[i] = mp_modmul_1(knap->prep_5a5[i], inv, p);
		h = mp_modmul_1(inv, c->p_fr[0], p);
		knap->ul_d_help[i][0] = h;
		mpz_mul_ui(knap->gmp_D[i][0], curr->gmp_help1, h);
		db += (double) h / (double) p;

		for (j = 1; j < 5; j++) {
			knap->ul_d_help[i][j] = mp_modmul_1(knap->ul_d_help_1[i][j], 
							knap->p0_inv_p[i], p);
		}
	}

	for (l = 0; l < npr_in_p; l++) {
		for (i = 0; i < npr_in_p; i++) {
			knap->prep_inv_table[i][l] =
				knap->prep_5a5[l] * knap->prep_inv_table[i][l];
		}
	}

	for (i = 0; i < npr_in_p; i++) {
		knap->prep_inv_table[i][npr_in_p] = mp_modmul_1(
						knap->prep_5a5[npr_in_p], 
						knap->p_inv_p0[i], p0);
	}

	db += 0.5 * (double) npr_in_p;
	for (i = 0; i < npr_in_p; i++)
		mpz_add(knap->gmp_D0, knap->gmp_D0, knap->gmp_D[i][0]);
	mpz_fdiv_r(curr->gmp_help1, curr->gmp_m0, curr->gmp_prod);
	mpz_sub(knap->gmp_disp, curr->gmp_m0, curr->gmp_help1);
	mpz_mul_si(curr->gmp_help1, curr->gmp_prod, (int) db);
	mpz_sub(knap->gmp_disp, knap->gmp_disp, curr->gmp_help1);
	mpz_add(knap->gmp_D0, knap->gmp_D0, knap->gmp_disp);

	profile_stop(PROF_CRT_AUX);
	profile_start(PROF_KAPPA);
	for (i = 0; i < npr_in_p; i++) {
		for (j = 0; j < 5; j++) {
			knap->raw_ull_kappa[i][j] = 0;
		}
	}

	/* compute kappa_0 and kappa_{i,j} */
	knap->ull_kappa0 = 0;

	p = p0;
	p2 = p * p;
	d0p = mpz_fdiv_ui(knap->gmp_D0, p2);
	a5_mod_p2 = mpz_fdiv_ui(curr->gmp_a5, p2);
	N_mod_p2 = mpz_fdiv_ui(curr->gmp_N, p2);

	hp = mp_expo_1(d0p, 5, p2);
	hp = mp_modmul_1(hp, a5_mod_p2, p2);
	hp = p2 - hp + N_mod_p2;
	if (hp % p)
		complain("%u %u %u neu-1r\n", p, d0p, hp);
	h = hp / p;
	inv = mp_modinv_1(N_mod_p2 % p0, p0);
	inv = mp_modmul_1(inv, knap->ull_kappa_help0[npr_in_p], p);
	h = mp_modmul_1(h, inv, p);
	h = mp_modmul_1(h, a5_mod_p2, p);
	h = mp_modmul_1(h, d0p, p);
	ulladdmul(&knap->ull_kappa0, h, knap->ull_kappa_p_inv + npr_in_p);

	for (i = 0; i < npr_in_p; i++) {
		curr_rat_prime_t *c = primes + p_ind[i];
		p = c->p_pr;
		p2 = p * p;
		d0p = mpz_fdiv_ui(knap->gmp_D0, p2);
		a5_mod_p2 = c->p_a5_mod_p2;
		N_mod_p2 = c->p_N_mod_p2;

		hp = mp_modmul_1(d0p, d0p, p2);
		hp = mp_modmul_1(hp, hp, p2);
		hp = mp_modmul_1(hp, d0p, p2);
		hp = mp_modmul_1(hp, a5_mod_p2, p2);
		hp = p2 - hp + N_mod_p2;
		if (hp % p)
			complain("%u %u %u neu-1r\n", p, d0p, hp);
		h = hp / p;
		inv = mp_modmul_1(a5_mod_p2, c->p_N_inv * 
			       knap->ull_kappa_help0[i], p);
		h = mp_modmul_1(h, inv * d0p, p);
		ulladdmul(&knap->ull_kappa0, h, knap->ull_kappa_p_inv + i);

		hh = mp_modmul_1(d0p, d0p, p2);
		hh = mp_modmul_1(hh, d0p, p2);
		h6 = mp_modmul_1(hh, hh, p2);	/* h6=d0p^6 */

		profile_start(PROF_KAPPA2);
		for (j = 1; j < 5; j++) {
			dp = mp_modmul_1(knap->p_mod_p2[i], 
					knap->ul_d_help[i][j], p2);
			hh = d0p + dp;
			h = mp_modmul_1(hh, hh, p2);
			h = mp_modmul_1(hh, h, p2);
			h = mp_modmul_1(h, h, p2);	/* h=(d0p+dp)^6 */
			h = mp_modmul_1(h6 + (p2 - h), a5_mod_p2, p2);
			h += mp_modmul_1(N_mod_p2, dp, p2);

			if (h % p)
				complain("%u %u %u neu1\n", p, d0p, dp);
			h /= p;
			h = mp_modmul_1(h, inv, p);
			ulladdmul(&(knap->raw_ull_kappa[i][j]), h,
				  knap->ull_kappa_p_inv + i);
		}
		for (l = 0; l < i; l++) {
			for (j = 1; j < 5; j++) {
				h = mp_modmul_1(knap->ul_d_help[i][j],
					     knap->prep_inv_table[i][l],
					     knap->prep_p[l]);
				ulladdmul(&(knap->raw_ull_kappa[i][j]), h,
					  knap->ull_kappa_p_inv + l);
			}
		}
		for (l = i + 1; l < npr_in_p + 1; l++) {
			for (j = 1; j < 5; j++) {
				h = mp_modmul_1(knap->ul_d_help[i][j],
					     knap->prep_inv_table[i][l],
					     knap->prep_p[l]);
				ulladdmul(&(knap->raw_ull_kappa[i][j]), h,
					  knap->ull_kappa_p_inv + l);
			}
		}
		profile_stop(PROF_KAPPA2);
	}

	profile_stop(PROF_KAPPA);
	profile_start(PROF_COMPUTE_KNAP_VALS);
#ifdef HAVE_FLOAT64
	{
		/* compute kappa0/prod - 5*a5*(D0-root)/prod^2 + 
		  		10*a5*(D0-root)^2/prod^2/D0 */
		uint64_t ullh;
		long double ddd, dde;

		mpz_sub(curr->gmp_help2, curr->gmp_root, knap->gmp_D0);
		dde = -mpz_get_ld(curr->gmp_help2);

		ddd = dde / dbl_prod;
		ddd *= dbl_5a5p;
		ullh = (uint64_t) (18446744073709551616.L *
				   (ddd - floorl(ddd)));
		ddd *= 2. * dde;
		ddd /= mpz_get_ld(curr->gmp_root);
		ullh -= (uint64_t) (18446744073709551616.L *
				    (ddd - floorl(ddd)));
		knap->lambda0 = ullh + knap->ull_kappa0;
	}
#else
	{
		uint64_t ullh1, ullh2;

		mpz_sub(curr->gmp_help2, knap->gmp_D0, curr->gmp_root);	/* D0-root */
		mpz_mul_ui(curr->gmp_help1, curr->gmp_a5, 5);	/* 5*a5 */
		mpz_mul(curr->gmp_help1, curr->gmp_help1, curr->gmp_help2);/* 5*a5*(D0-root) */
		mpz_mul(curr->gmp_help3, curr->gmp_prod, curr->gmp_prod);	/* prod^2 */
		mpz_mul_2exp(curr->gmp_help1, curr->gmp_help1, 64);
		mpz_fdiv_q(curr->gmp_help4, curr->gmp_help1, curr->gmp_help3);
		if (mpz_sizeinbase(curr->gmp_help4, 2) > 64)
			complain("lambda-comp\n");
		mpz_get_ull_64(&ullh1, curr->gmp_help4);

		mpz_add(curr->gmp_help2, curr->gmp_help2, curr->gmp_help2);	/* 10*a5*(D0-root)^2 * 2^64 */
		mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_help2);
		mpz_fdiv_q(curr->gmp_help4, curr->gmp_help4, knap->gmp_D0);

		if (mpz_sizeinbase(curr->gmp_help4, 2) > 64)
			complain("lambda-comp\n");
		mpz_get_ull_64(&ullh2, curr->gmp_help4);
		knap->lambda0 = -ullh1 + ullh2 + knap->ull_kappa0;
	}
#endif

	{
		uint64_t ullh;

		for (i = 0; i < npr_in_p; i++) {
#ifdef HAVE_FLOAT64
			ullh = dull(dbl_5a5p * primes[p_ind[i]].ld_p_inv);
#else
			ull_mulh(&ullh, &(primes[p_ind[i]].ull_p_inv), 
					&(ull_5a5p));
#endif
			for (j = 1; j < 5; j++) {
				ulladdmul(&(knap->raw_ull_kappa[i][j]),
					  knap->ul_d_help[i][j], &ullh);
			}
		}
		for (j = 0; j < 5; j++)
			knap->raw_ull_kappa[0][j] += knap->lambda0;
	}

	for (i = knap->len1; i < npr_in_p; i++) {
		for (j = 0; j < 5; j++) {
			knap->raw_ull_kappa[i][j] = -(knap->raw_ull_kappa[i][j]);
		}
	}

	db = a3_max / mpz_get_d(curr->gmp_m0);
	if ((db >= 1.) || (db < 0.))
		complain("a3-bound so high that all pols will pass: %f\n", db);
	db *= 18446744073709551616.;
	knap->raw_ull_bound = (uint64_t) db;
	knap->raw_ull_bound +=
		(uint64_t) ((npr_in_p + 2) * npr_in_p * 
				primes[p_ind[npr_in_p - 1]].p_pr);
	if (verbose > 3)
		printf("raw_ull_bound: %" PRIu64 "\n", knap->raw_ull_bound);
	knap->raw_bound = (unsigned int) (knap->raw_ull_bound >> 31) + 1 + 4;
	profile_stop(PROF_COMPUTE_KNAP_VALS);
}

/*------------------------------------------------------------------------*/
#ifdef HAVE_ASM
extern int asm_hash1();
extern int asm_hash2();

unsigned int raw_bound asm("raw_bound");
int *raw_cand_ptr asm("raw_cand_ptr");
unsigned int *hashdataptr asm("hashdataptr");
unsigned int *s11l asm("s11l");
unsigned int *s12l asm("s12l");
unsigned int *s21l asm("s21l");
unsigned int *s22l asm("s22l");
unsigned int s11len asm("s11len");
unsigned int s12len asm("s12len");
unsigned int s21len asm("s21len");
unsigned int s22len asm("s22len");
#endif

/*------------------------------------------------------------------------*/
int
knapsack_raw(knapsack_data_t *knap,
		stage1_stat_t *stats, int npr_in_p)
{
	int j;
	int res;

	profile_start(PROF_KNAP_RAW_COMBINE);
	for (j = 0; j < 5; j++)
		knap->raw_ull_kappa[0][j] += (8589934592ULL + knap->raw_ull_bound);

	combine_raw(knap, knap->s11l, 0, knap->len11);
	combine_raw(knap, knap->s12l, knap->len11, knap->len1);
	combine_raw(knap, knap->s21l, knap->len1, knap->len1 + knap->len21);
	combine_raw(knap, knap->s22l, knap->len1 + knap->len21, npr_in_p);

	for (j = 0; j < 5; j++)
		knap->raw_ull_kappa[0][j] -= (8589934592ULL + knap->raw_ull_bound);
	profile_stop(PROF_KNAP_RAW_COMBINE);


	knap->nraw_store = 0;
	profile_start(PROF_RAW_HASH1);
#ifdef HAVE_ASM
	raw_bound = knap->raw_bound;
	hashdataptr = knap->hashdata;
	s11l = knap->s11l;
	s11len = knap->s11len;
	s12l = knap->s12l;
	s12len = knap->s12len;
	s21l = knap->s21l;
	s21len = knap->s21len;
	s22l = knap->s22l;
	s22len = knap->s22len;
	res = asm_hash1();
#else
	res = raw_hash_1(knap);
#endif
	profile_stop(PROF_RAW_HASH1);
	if (res)
		return -1;	/* too many entries in one hash-line */


	knap->raw_cand_ptr = knap->raw_cand;
	profile_start(PROF_RAW_HASH2);
#ifdef HAVE_ASM
	raw_cand_ptr = knap->raw_cand_ptr;
	asm_hash2();
	knap->raw_cand_ptr = raw_cand_ptr;
#else
	raw_hash_2(knap);
#endif
	profile_stop(PROF_RAW_HASH2);
	if (knap->raw_cand_ptr == knap->raw_cand)
		return 0;


	profile_start(PROF_RAW_HASH3);
	res = raw_hash_3(knap);
	profile_stop(PROF_RAW_HASH3);
	if (res)
		return -1;	/* should happen only for large npr_in_p */


	profile_start(PROF_RAW_HASH4);
	raw_hash_4(knap);
	profile_stop(PROF_RAW_HASH4);
	return knap->nraw_store;
}

/*------------------------------------------------------------------------*/
void
check_raw(curr_poly_t *curr, knapsack_data_t *knap,
		bounds_t *bounds, stage1_stat_t *stats)
{
	int i, j, i1, i2;
	uint64_t sum1, sum2;
	int nraw_store = knap->nraw_store;
	uint64_t *raw_stored_pairs = knap->raw_stored_pairs;
	int len1 = knap->len1;
	int len2 = knap->len2;
	int s1len = knap->s1len;
	int s2len = knap->s2len;
	uint64_t ull_bound = knap->ull_bound;

	i = 0;			/* remove duplicates: */
	while (i < nraw_store - 1) {
		for (j = i + 1; j < nraw_store; j++) {
			if ((raw_stored_pairs[2 * i] == raw_stored_pairs[2 * j])
			    && (raw_stored_pairs[2 * i + 1] ==
				raw_stored_pairs[2 * j + 1]))
				break;
		}
		if (j < nraw_store) {
			raw_stored_pairs[2 * j] =
				raw_stored_pairs[2 * nraw_store - 2];
			raw_stored_pairs[2 * j + 1] =
				raw_stored_pairs[2 * nraw_store - 1];
			nraw_store--;
		}
		i++;
	}
	knap->nraw_store = nraw_store;

	for (i = 0; i < nraw_store; i++) {
		i1 = raw_stored_pairs[2 * i];
		if (i1 >= s1len)
			complain("check_raw_stored_pairs.1 %d %d \n", i1, i);
		i2 = raw_stored_pairs[2 * i + 1];
		if (i2 >= s2len)
			complain("check_raw_stored_pairs.2 %d %d \n", i2, i);

		sum1 = 0;
		sum2 = 0;
		for (j = 0; j < len1; j++) {
			sum1 += knap->ull_kappa[j][i1 % DEG];
			i1 -= (i1 % DEG);
			i1 /= DEG;
		}
		for (j = len1; j < len1 + len2; j++) {
			sum2 += knap->ull_kappa[j][i2 % DEG];
			i2 -= (i2 % DEG);
			i2 /= DEG;
		}
		if ((sum1 - sum2 > ull_bound) && (sum2 - sum1 > ull_bound)) {
#if 1
			if (sum1 > sum2)
				printf("%" PRIu64 " %" PRIu64 "! ", sum1 - sum2,
				       ull_bound);
			else
				printf("%" PRIu64 " %" PRIu64 "! ", sum2 - sum1,
				       ull_bound);
#else
			printf("!");
#endif
			continue;
		}

		i1 = raw_stored_pairs[2 * i];
		i2 = raw_stored_pairs[2 * i + 1];
		mpz_set_ui(curr->gmp_d, 0);
		mpz_set_ui(curr->gmp_help4, 0);
		for (j = 0; j < len1; j++) {
			mpz_add(curr->gmp_d, curr->gmp_d, knap->gmp_D[j][i1 % DEG]);
			mpz_add(curr->gmp_help4, curr->gmp_help4, knap->gmp_kappa[j][i1 % DEG]);
			i1 -= (i1 % DEG);
			i1 /= DEG;
		}
		for (j = len1; j < len1 + len2; j++) {
			mpz_add(curr->gmp_d, curr->gmp_d, knap->gmp_D[j][i2 % DEG]);
			mpz_add(curr->gmp_help4, curr->gmp_help4, knap->gmp_kappa[j][i2 % DEG]);
			i2 -= (i2 % DEG);
			i2 /= DEG;
		}
		mpz_add(curr->gmp_d, curr->gmp_d, knap->gmp_disp);
		mpz_add(curr->gmp_help4, curr->gmp_help4, knap->gmp_kappa0);
		/* check */
		mpz_pow_ui(curr->gmp_help3, curr->gmp_d, DEG);
		mpz_mul(curr->gmp_help3, curr->gmp_help3, curr->gmp_a5);
		mpz_sub(curr->gmp_help3, curr->gmp_help3, curr->gmp_N);
		mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help2, curr->gmp_help3, curr->gmp_prod);
		if (mpz_sgn(curr->gmp_help2)) {
			printf("\n");
			mpz_out_str(stdout, 10, curr->gmp_d);
			printf("\n");
			mpz_out_str(stdout, 10, curr->gmp_prod);
			printf("\n");
			mpz_out_str(stdout, 10, curr->gmp_a5);
			printf("\n");
			mpz_out_str(stdout, 10, curr->gmp_help3);
			printf("\n");
			complain("check2.raw\n");
		}
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
			complain("check_raw_stored_pairs check\n");
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
}
