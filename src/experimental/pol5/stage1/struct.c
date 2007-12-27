#include "stage1_impl.h"

static const u32 pr_mod5[NPR5] = {
	11, 31, 41, 61, 71, 101, 131, 151, 181, 191,
	211, 241, 251, 271, 281, 311, 331, 401, 421, 431,
	461, 491, 521, 541, 571, 601, 631, 641, 661, 691,
	701, 751, 761, 811, 821, 881, 911, 941, 971, 991,
	1021, 1031, 1051, 1061, 1091, 1151, 1171, 1181, 1201, 1231
};				/* primes =1 mod 5, < 2^(31/3) */

static const unsigned char ucmask[8] = { 
	0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 
};

/*------------------------------------------------------------------------*/
void
stage1_stat_init(stage1_stat_t *stats)
{
	memset(stats, 0, sizeof(stage1_stat_t));
	profile_init(&stats->profile, PROF_MAX);
}

/*------------------------------------------------------------------------*/
void
stage1_stat_free(stage1_stat_t *stats)
{
	profile_free(&stats->profile);
	memset(stats, 0, sizeof(stage1_stat_t));
}

/*------------------------------------------------------------------------*/
void 
stage1_bounds_init(bounds_t *bounds, poly_stage1_t *data)
{
	mpz_init_set(bounds->gmp_a5_begin, data->gmp_a5_begin);
	mpz_init_set(bounds->gmp_a5_end, data->gmp_a5_end);
	bounds->norm_max = data->norm_max;
	bounds->p0_limit = data->p0_limit;
	if (bounds->p0_limit == 0 || bounds->p0_limit > P0_MAX)
		bounds->p0_limit = P0_MAX;
}

/*------------------------------------------------------------------------*/
void 
stage1_bounds_free(bounds_t *bounds)
{
	mpz_clear(bounds->gmp_a5_begin);
	mpz_clear(bounds->gmp_a5_end);
}

/*------------------------------------------------------------------------*/
void 
stage1_bounds_update(bounds_t *bounds, double N, double a5)
{
	bounds->skewness_max = pow(bounds->norm_max / a5, 0.4);
	bounds->skewness_min = pow(N / a5 / 
				pow(bounds->norm_max, 5.), 2./15.);
	bounds->a3_max = bounds->norm_max / 
				sqrt(bounds->skewness_min);
	bounds->p_size_max = log(bounds->a3_max / 
				bounds->skewness_min);
}

/*------------------------------------------------------------------------*/
void
curr_poly_init(curr_poly_t *curr, bounds_t *bounds, poly_stage1_t *data)
{
	int i;
	double dbl_a5, dbl_N, dbl_a5_min, dbl_a5_max;

	mpz_init_set(curr->gmp_N, data->gmp_N);
	mpz_init(curr->gmp_root);
	mpz_init(curr->gmp_prod);
	mpz_init(curr->gmp_prod_base);
	mpz_init(curr->gmp_help1);
	mpz_init(curr->gmp_help2);
	mpz_init(curr->gmp_help3);
	mpz_init(curr->gmp_help4);
	mpz_init(curr->gmp_a5);
	mpz_init(curr->gmp_d);
	mpz_init(curr->gmp_m0);
	curr->outputfile = data->outfile;

	mpz_fdiv_q_ui(curr->gmp_a5, bounds->gmp_a5_begin, MULTIPLIER);
	mpz_mul_ui(curr->gmp_a5, curr->gmp_a5, MULTIPLIER);
	if (mpz_cmp(curr->gmp_a5, bounds->gmp_a5_begin) < 0)
		mpz_add_ui(curr->gmp_a5, curr->gmp_a5, MULTIPLIER);

	dbl_N = mpz_get_d(curr->gmp_N);
	dbl_a5 = mpz_get_d(curr->gmp_a5);

	dbl_a5_max = pow(bounds->norm_max, 8.) / dbl_N;
	dbl_a5_max = sqrt(dbl_a5_max);
	dbl_a5_min = 0.;
	for (i = 0; i < data->npr_in_p; i++)
		dbl_a5_min += log((double) pr_mod5[i]);
	dbl_a5_min *= 5.;
	dbl_a5_min -= 10. * log(bounds->norm_max);
	dbl_a5_min += log(dbl_N);
	dbl_a5_min = exp(dbl_a5_min);
	dbl_a5_min = floor(dbl_a5_min) + 1.;

	printf("Parameters: number: %.3e, norm_max: %.3e, a5: %.3e - %.3e\n",
	       dbl_N, bounds->norm_max, dbl_a5, mpz_get_d(bounds->gmp_a5_end));
	printf("Accepted a5-range: %f - %f \n", dbl_a5_min, dbl_a5_max);

	if (dbl_a5_min > dbl_a5)
		complain("a5_begin smaller than %f\n", dbl_a5_min);
	if (dbl_a5_max < mpz_get_d(bounds->gmp_a5_end))
		complain("a5_end bigger than %f\n", dbl_a5_max);

	printf("Searching in subinterval: ");
	mpz_out_str(stdout, 10, curr->gmp_a5);
	printf(" - ");
	mpz_out_str(stdout, 10, bounds->gmp_a5_end);
	printf("\n");

	mpz_sub_ui(curr->gmp_a5, curr->gmp_a5, MULTIPLIER);
}

/*------------------------------------------------------------------------*/
void
curr_poly_free(curr_poly_t *curr)
{
	mpz_clear(curr->gmp_N);
	mpz_clear(curr->gmp_prod);
	mpz_clear(curr->gmp_prod_base);
	mpz_clear(curr->gmp_root);
	mpz_clear(curr->gmp_help1);
	mpz_clear(curr->gmp_help2);
	mpz_clear(curr->gmp_help3);
	mpz_clear(curr->gmp_help4);
	mpz_clear(curr->gmp_a5);
	mpz_clear(curr->gmp_d);
	mpz_clear(curr->gmp_m0);
}

/*------------------------------------------------------------------------*/
static u32
is_5power(u32 a, u32 p)
{			/* p!=5, 0 is not considered as fifth power */
	if (p % 5 != 1)
		return 0;
	if (mp_expo_1(a, (p - 1) / 5, p) == 1)
		return 1;
	return 0;
}

void
rat_factor_init(rat_factor_t *rat_factors, mpz_t N, mpz_t a5_start)
{
	u32 i, j;
	u32 len;

	for (i = 0; i < NPR5; i++) {
		rat_prime_t *r = rat_factors->primes + i;
		u32 p = pr_mod5[i];
		u32 N_mod_pr;
		u32 N_inv_mod_pr;
		u32 a5_mod_pr;

		N_mod_pr = mpz_fdiv_ui(N, p);
		N_inv_mod_pr = mp_modinv_1(N_mod_pr, p);
		a5_mod_pr = mpz_fdiv_ui(a5_start, p);
		r->pr = p;
		r->pr_step = mp_modmul_1(MULTIPLIER, N_inv_mod_pr, p);
		r->pr_start = mp_modmul_1(a5_mod_pr, N_inv_mod_pr, p);
		r->p_log = log((double)p);
	}

	for (i = len = 0; i < NPR5; i++)
		len += (pr_mod5[i] / 8 + 1);
	rat_factors->bitarray_5power[0] = (unsigned char *)xcalloc(len,
						sizeof(unsigned char));
	for (i = 1; i < NPR5; i++) {
		rat_factors->bitarray_5power[i] =
			rat_factors->bitarray_5power[i - 1] + 
				(pr_mod5[i - 1] / 8 + 1);
	}
	for (i = 0; i < NPR5; i++) {
		u32 p = pr_mod5[i];
		unsigned char *bitarray = rat_factors->bitarray_5power[i];

		for (j = 1; j < p; j++) {
			if (is_5power(j, p)) {
				bitarray[j / 8] |= ucmask[j % 8];
			}
		}
	}
}

/*------------------------------------------------------------------------*/
void
rat_factor_free(rat_factor_t *rat_factors)
{
	free(rat_factors->bitarray_5power[0]);
}

/*------------------------------------------------------------------------*/
void 
curr_rat_factor_init(curr_rat_factor_t *c, u32 npr_in_p)
{
	memset(c, 0, sizeof(curr_rat_factor_t));
	c->npr_in_p = npr_in_p;
	c->p0_list_len_max = 256;
	c->p0_list = (aux_factor_t *)xmalloc(c->p0_list_len_max *
					 sizeof(aux_factor_t));
}

/*------------------------------------------------------------------------*/
void 
curr_rat_factor_free(curr_rat_factor_t *c)
{
	free(c->p0_list);
}

/*------------------------------------------------------------------------*/
static int
coprime(unsigned int a, unsigned int b)
{
	unsigned int r, bb = b, aa = a;

	while (bb) {
		r = aa % bb;
		aa = bb;
		bb = r;
	}
	return (aa == 1);
}

/*------------------------------------------------------------------------*/
static void
aux_factor_store(curr_rat_factor_t *curr,
		unsigned int n, unsigned int r)
{
	if (curr->p0_list_len >= curr->p0_list_len_max) {
		curr->p0_list_len_max *= 2;
		curr->p0_list = (aux_factor_t *) xrealloc(curr->p0_list,
						  curr->p0_list_len_max *
						  sizeof(aux_factor_t));
	}
	curr->p0_list[curr->p0_list_len].p = n;
	curr->p0_list[curr->p0_list_len].r = r;
	curr->p0_list_len++;
}

/*------------------------------------------------------------------------*/
static const unsigned int smallprimes[] = {
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
	73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
	127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
	179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
	233, 239, 241, 251
};

#define NUM_SMALL_PRIMES (sizeof(smallprimes) / sizeof(smallprimes[0]))

static unsigned int
aux_factor_phi(unsigned int n)
{				/* n<=2^16 */
	unsigned int i;
	unsigned int phi, d, p;

	if (n > 65536)
		complain("aux_factor_phi\n");
	d = n;
	phi = 1;
	for (i = 0; i < NUM_SMALL_PRIMES; i++) {
		p = smallprimes[i];
		if (d % p == 0) {
			phi *= (p - 1);
			d /= p;
			while (d % p == 0) {
				phi *= p;
				d /= p;
			}
			if (d == 1)
				break;
		}
	}
	if (d > 1)
		phi *= (d - 1);
	return phi;
}

/*------------------------------------------------------------------------*/
static void
aux_factor_roots(curr_rat_factor_t *curr_factors, 
		unsigned int n, mpz_t gmp_N, mpz_t gmp_a5)
{
	int i;
	unsigned int h, r, a5, N, p;

	a5 = mpz_fdiv_ui(gmp_a5, n);
	if (!coprime(a5, n))
		return;
	N = mpz_fdiv_ui(gmp_N, n);
	if (n < 20) {
		if (n == 11)
			return;
		for (r = 1; r < n; r++) {
			h = mp_expo_1(r, 5, n);
			h = mp_modmul_1(h, a5, n);
			if (h == N)
				aux_factor_store(curr_factors, n, r);
		}
		return;
	}

	for (i = 0; i < NPR5; i++) {
		p = pr_mod5[i];
		if (p > n)
			break;
		if (n % p == 0)
			return;
	}

	/* for deg=5 exactly one root exists */
	h = aux_factor_phi(n);
	if (h % 5 == 0)
		return;

	r = 1 + h;
	while (r % 5)
		r += h;
	r /= 5;

	h = mp_modinv_1(a5, n);
	h = mp_modmul_1(h, N, n);
	r = mp_expo_1(h, r, n);
	aux_factor_store(curr_factors, n, r);
}

/*------------------------------------------------------------------------*/
void
curr_rat_factor_fill_aux(curr_rat_factor_t *curr_factors, 
			mpz_t N, mpz_t a5,
			double dbl_p0_max, u32 p0_limit)
{
	unsigned int n, p0m;

	curr_factors->p0_list[0].p = 1;
	curr_factors->p0_list[0].r = 0;
	curr_factors->p0_list_len = 1;
	if (dbl_p0_max > (double) p0_limit)
		p0m = p0_limit;
	else
		p0m = (unsigned int) rint(dbl_p0_max);

	if (p0m < 2)
		return;
	for (n = 2; n <= p0m; n++) {
		if (coprime(n, MULTIPLIER))
			aux_factor_roots(curr_factors, n, N, a5);
	}
}

/*------------------------------------------------------------------------*/
u32
curr_rat_factor_find(rat_factor_t *rat_factors, 
			curr_rat_factor_t *curr_factors,
			double p_size_max)
{
	int i, j, ind;
	unsigned int h, p;
	double size;
	u32 npr_in_p = curr_factors->npr_in_p;

	for (i = ind = 0; i < NPR5; i++) {
		rat_prime_t *r = rat_factors->primes + i;
		unsigned char *bitarray = rat_factors->bitarray_5power[i];

		p = r->pr;
		r->pr_start += r->pr_step;
		if (r->pr_start >= p)
			r->pr_start -= p;
		h = r->pr_start;
		if (bitarray[h / 8] & ucmask[h % 8]) {
			curr_rat_prime_t *curr = curr_factors->primes + ind++;
			curr->p_pr = p;
			curr->p_log = r->p_log;
			if (ind == npr_in_p)
				break;
		}
	}
	i++;
	if (ind < npr_in_p)
		return 0;
	size = 0.;
	for (j = 0; j < npr_in_p; j++) {
		size += curr_factors->primes[j].p_log;
	}
	if (size > p_size_max) {	/* just updating */
		for (; i < NPR5; i++) {
			rat_prime_t *r = rat_factors->primes + i;
			p = r->pr;
			r->pr_start += r->pr_step;
			if (r->pr_start >= p)
				r->pr_start -= p;
		}
		return 0;
	}			/* else search for further candidates */
	for (; i < NPR5; i++) {
		rat_prime_t *r = rat_factors->primes + i;
		unsigned char *bitarray = rat_factors->bitarray_5power[i];

		p = r->pr;
		r->pr_start += r->pr_step;
		if (r->pr_start >= p)
			r->pr_start -= p;
		h = r->pr_start;
		if (bitarray[h / 8] & ucmask[h % 8]) {
			curr_rat_prime_t *curr = curr_factors->primes + ind++;
			curr->p_pr = p;
			curr->p_log = r->p_log;
		}
	}
	curr_factors->npr_total = ind;
	return 1;
}

/*------------------------------------------------------------------------*/
static void
fifth_root(u32 *ro, u32 a, u32 p)
{				/* perhaps improve this */
	u32 ind;
	u32 h, r;

	ind = 0;
	for (r = 1; r < p; r++) {
		h = r * r;
		h %= p;
		h = h * h;
		h %= p;
		h *= r;
		h %= p;
		if (h == a) {
			ro[ind++] = r;
			if (ind == 5)
				break;
		}
	}
	if (ind < 5)
		complain("cannot find enough fifth roots %d %d\n", a, p);
}

/*------------------------------------------------------------------------*/
void
curr_rat_factor_fill(curr_rat_factor_t *curr_factors,
			mpz_t N, mpz_t a5)
{
	u32 i, j;
	u32 h, p;
	u32 npr_total = curr_factors->npr_total;
	curr_rat_prime_t *primes = curr_factors->primes;

	for (i = 0; i < npr_total; i++) {
		curr_rat_prime_t *c = primes + i;
		p = c->p_pr;
		c->p_N_mod_p2 = mpz_fdiv_ui(N, p * p);
		c->p_a5_mod_p2 = mpz_fdiv_ui(a5, p * p);
		c->p_minus5a5_mod_p = mp_modmul_1(p - 5, c->p_a5_mod_p2, p);
		c->p_N_inv = mp_modinv_1(c->p_N_mod_p2 % p, p);
		h = c->p_a5_mod_p2 % p;
		h = mp_modmul_1(mp_modinv_1(h, p), c->p_N_mod_p2, p);
#ifdef HAVE_FLOAT64
		c->ld_p_inv = 1.L / ((long double) p);
#endif
		c->ull_p_inv = 18446744073709551615ULL / ((uint64_t) p);
		fifth_root(c->p_fr, h, p);
	}
	for (i = 0; i < npr_total; i++) {
		for (j = 0; j < npr_total; j++) {
			if (i != j) {
				u32 pi = primes[i].p_pr;
				u32 pj = primes[j].p_pr;
				curr_factors->p_inv_table[i][j] =
						mp_modinv_1(pi % pj, pj);
			}
		}
	}
}

/*------------------------------------------------------------------------*/
static int
five_pow(int n)
{
	int res, i;

	res = 1;
	for (i = 0; i < n; i++)
		res *= 5;
	return res;
}

void
knapsack_data_init(knapsack_data_t *knap, u32 npr_in_p)
{
	int i, j;

	memset(knap, 0, sizeof(knapsack_data_t));

	mpz_init(knap->gmp_disp);
	mpz_init(knap->gmp_D0);

	for (i = 0; i < NPR5; i++) {
		for (j = 0; j < 5; j++) {
			mpz_init(knap->gmp_D[i][j]);
		}
	}
	mpz_init(knap->gmp_kappa0);
	for (i = 0; i < NPR5; i++) {
		for (j = 0; j < 5; j++) {
			mpz_init(knap->gmp_kappa[i][j]);
		}
	}
	for (i = 0; i < NPR5 + 1; i++) {
		mpz_init(knap->gmp_kappa_help[i]);
	}

	knap->len1 = npr_in_p / 2;
	knap->len2 = npr_in_p - knap->len1;	/* len1<=len2 */
	knap->len12 = knap->len1 / 2;
	knap->len11 = knap->len1 - knap->len12;
	knap->len22 = knap->len2 / 2;
	knap->len21 = knap->len2 - knap->len22;
	knap->s11len = five_pow(knap->len11);
	knap->s12len = five_pow(knap->len12);
	knap->s21len = five_pow(knap->len21);
	knap->s22len = five_pow(knap->len22);
	knap->s1len = knap->s11len * knap->s12len;
	knap->s2len = knap->s21len * knap->s22len;
	knap->shelp = (uint64_t *) xmalloc(knap->s21len * 
				sizeof(uint64_t)); /* s21len is maximum */
	knap->s11l = (unsigned int *) xmalloc(knap->s11len * 
						sizeof(unsigned int));
	knap->s12l = (unsigned int *) xmalloc(knap->s12len * 
						sizeof(unsigned int));
	knap->s21l = (unsigned int *) xmalloc(knap->s21len * 
						sizeof(unsigned int));
	knap->s22l = (unsigned int *) xmalloc(knap->s22len * 
						sizeof(unsigned int));
	knap->s1sortl = (unsigned int *) xmalloc(knap->s1len * 
						sizeof(unsigned int));
	knap->s2sortl = (unsigned int *) xmalloc(knap->s2len * 
						sizeof(unsigned int));
	knap->s1 = (uint64_t *) xmalloc(knap->s1len * sizeof(uint64_t));
	knap->s2 = (uint64_t *) xmalloc(knap->s2len * sizeof(uint64_t));
	knap->s1sort = (uint64_t *) xmalloc(knap->s1len * sizeof(uint64_t));
	knap->s2sort = (uint64_t *) xmalloc(knap->s2len * sizeof(uint64_t));
	knap->raw_cand = (int *) xmalloc(knap->s2len * sizeof(int));
	knap->raw_cand_hash = (unsigned int *) xmalloc(knap->s2len * 
						sizeof(unsigned int));
	knap->hashdata = (unsigned int *) xmalloc( (32 * NHASH + (NHASH >> 2)) *
						sizeof(unsigned int));
	knap->hashptr32 = (unsigned int **) xmalloc(NHASH * 
						sizeof(unsigned int *));
	knap->sort = (int *) xmalloc(NHASH * sizeof(int));
	knap->hashptr = (uint64_t **) xmalloc(NHASH * sizeof(uint64_t *));
	knap->stored_pairs = NULL;
	knap->raw_stored_pairs = NULL;
}

/*------------------------------------------------------------------------*/
void
knapsack_data_free(knapsack_data_t *knap)
{
	int i, j;

	mpz_clear(knap->gmp_disp);
	mpz_clear(knap->gmp_D0);

	for (i = 0; i < NPR5; i++) {
		for (j = 0; j < 5; j++) {
			mpz_clear(knap->gmp_D[i][j]);
		}
	}
	mpz_clear(knap->gmp_kappa0);
	for (i = 0; i < NPR5; i++) {
		for (j = 0; j < 5; j++) {
			mpz_clear(knap->gmp_kappa[i][j]);
		}
	}
	for (i = 0; i < NPR5 + 1; i++) {
		mpz_clear(knap->gmp_kappa_help[i]);
	}

	free(knap->shelp);
	free(knap->s11l);
	free(knap->s12l);
	free(knap->s21l);
	free(knap->s22l);
	free(knap->s1sortl);
	free(knap->s2sortl);
	free(knap->s1);
	free(knap->s2);
	free(knap->s1sort);
	free(knap->s2sort);
	free(knap->raw_cand);
	free(knap->raw_cand_hash);
	free(knap->hashdata);
	free(knap->hashptr32);
	free(knap->sort);
	free(knap->hashptr);
	free(knap->stored_pairs);
	free(knap->raw_stored_pairs);
}

