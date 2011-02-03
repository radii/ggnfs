#include "stage2_impl.h"

#define  MAX_PRIME_PROJ       100
#define  MAX_PRIME_AFF        200
#define  NPROJ_PRIMES          25
#define  MAX_X       1000000	/* !!! */
#define  MAX_Y           100	/* !!! */
#define  SIEVELEN           8192

static const unsigned int primes[NAFF_PRIMES] = {
	2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
	73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
	127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
	179, 181, 191, 193, 197, 199
};

/*-------------------------------------------------------------------------*/
static int
divmod(int a, int b, int pk)
{				/* a/b mod pk , a, b < 2^15 */
	int inv;

	inv = invert(b, pk);
	if (!inv)
		Schlendrian("divmod\n");
	inv = inv % pk;
	if (inv < 0)
		inv += pk;
	inv *= (a % pk);
	inv = inv % pk;
	if (inv < 0)
		inv += pk;
	return inv;
}

/*-------------------------------------------------------------------------*/
static double
compute_standard_alpha(void)
{
	int p, pk, k;
	double dp, dpk, dlp, dlog, standardalpha;

	standardalpha = 0.;
	for (k = 0; k < NAFF_PRIMES; k++) {
		p = primes[k];
		dp = (double) p;
		dlp = (double) (log(dp));
		pk = p;
		while (pk < MAX_PRIME_AFF) {
			dpk = (double) pk;
			dlog = dlp / dpk * dp / (dp + 1.);
			standardalpha -= dlog;
			pk *= p;
		}
	}
	return standardalpha;
}

/*-------------------------------------------------------------------------*/
double
compute_proj_alpha(curr_poly_t *c)
{
	unsigned int i, j, k;
	unsigned int p, p2, p3;
	unsigned int w, b3, b4, b5;
	double value, dp, dl;
	double table[MAX_PRIME_PROJ];

	value = 0.;
	for (k = 0; k < NPROJ_PRIMES; k++) {
		p = primes[k];
		if (mpz_mod_ui(c->gmp_help1, c->gmp_a[5], p) == 0) {
			dp = (double) p;
			value -= log(dp) / (dp + 1.);
			if (MAX_PRIME_PROJ / p / p < p) {/* consider only p^2 */
				p2 = p * p;
				if (mpz_mod_ui(c->gmp_help1, 
						c->gmp_a[5], p2) == 0) {
					if (mpz_mod_ui(c->gmp_help1, 
							c->gmp_a[4], p) == 0)
						value -= log(dp) * 
							dp / (dp * dp - 1);
					else
						value -= log(dp) / 
							(dp * dp - 1);
				}
				else {
					if (mpz_mod_ui(c->gmp_help1, 
							c->gmp_a[4], p) == 0)
						value -= log(dp) / 
							(dp * dp - 1);
				}
			}
			else {
				p2 = p * p;
				p3 = p * p2;
				for (i = 0; i < p2; i++)
					table[i] = 0.;
				dl = log(dp) / ((dp * dp - 1.) * dp * dp * dp);
				for (i = 0; i < p; i++)
					table[i * p] = dl;
				table[0] =
					log(dp) * (2. + 1. / (dp - 1.)) / 
						((dp * dp - 1.) * dp * dp * dp);
				b5 = mpz_mod_ui(c->gmp_help1, c->gmp_a[5], p3);
				b5 /= p;
				b4 = mpz_mod_ui(c->gmp_help1, c->gmp_a[4], p2);
				b3 = mpz_mod_ui(c->gmp_help1, c->gmp_a[3], p2);
/* (a5/p)*x^2 + a4*x*(y/p) + a3*p*(y/p)^2, b5=(a5/p),b4=a4,b3=a3,j=y/p,i=x */
				for (j = 0; j < p2; j++) {
					for (i = 0; i < p2; i++) {
						if (i % p) {
							w = b5 * i * i +
								b4 * i * j +
								b3 * p * j * j;
							w = w % p2;
							value -= table[w];
						}
					}
				}
			}
		}
	}
	return value;
}

/*-------------------------------------------------------------------------*/
static void
compute_sieving_area(double max_norm_1, double *dbl_a, double dbl_d,
			rs_bound_t *rs_bound, double skewness,
			double alpha_proj)
{
	int x0, x1, xx, y0, y1, yy;
	double da1, da0, dd;
	double v, v0;

	v = max_norm_1 * exp(-alpha_proj);

	da1 = dbl_a[1];
	da0 = dbl_a[0];
	dd = -dbl_d;
	y0 = 0;
	y1 = MAX_Y;
	dbl_a[1] = da1 - ((double) y1) * dd;
	v0 = sqrt(ifs(dbl_a, find_best_skewness(dbl_a, skewness)));
	if (v0 > v) {
		while (y1 - y0 > 1) {
			yy = (y0 + y1) / 2;
			dbl_a[1] = da1 - ((double) yy) * dd;
			v0 = sqrt(ifs
				  (dbl_a, find_best_skewness(dbl_a, skewness)));
			if (v0 > v)
				y1 = yy;
			else
				y0 = yy;
		}
	}
	rs_bound->ymin = -y1;
	y0 = 0;
	y1 = MAX_Y;
	dbl_a[1] = da1 + ((double) y1) * dd;
	v0 = sqrt(ifs(dbl_a, find_best_skewness(dbl_a, skewness)));
	if (v0 > v) {
		while (y1 - y0 > 1) {
			yy = (y0 + y1) / 2;
			dbl_a[1] = da1 + ((double) yy) * dd;
			v0 = sqrt(ifs
				  (dbl_a, find_best_skewness(dbl_a, skewness)));
			if (v0 > v)
				y1 = yy;
			else
				y0 = yy;
		}
	}
	rs_bound->ymax = y1;
	dbl_a[1] = da1;

	x0 = 0;
	x1 = MAX_X;
	dbl_a[0] = da0 - ((double) x1) * dd;
	v0 = sqrt(ifs(dbl_a, find_best_skewness(dbl_a, skewness)));
	if (v0 > v) {
		while (x1 - x0 > 1) {
			xx = (x0 + x1) / 2;
			dbl_a[0] = da0 - ((double) xx) * dd;
			v0 = sqrt(ifs
				  (dbl_a, find_best_skewness(dbl_a, skewness)));
			if (v0 > v)
				x1 = xx;
			else
				x0 = xx;
		}
	}
	rs_bound->xmin = -x1;
	x0 = 0;
	x1 = MAX_X;
	dbl_a[0] = da0 + ((double) x1) * dd;
	v0 = sqrt(ifs(dbl_a, find_best_skewness(dbl_a, skewness)));
	if (v0 > v) {
		while (x1 - x0 > 1) {
			xx = (x0 + x1) / 2;
			dbl_a[0] = da0 + ((double) xx) * dd;
			v0 = sqrt(ifs(dbl_a, 
					find_best_skewness(dbl_a, skewness)));
			if (v0 > v)
				x1 = xx;
			else
				x0 = xx;
		}
	}
	rs_bound->xmax = x1;
	dbl_a[0] = da0;
}

/*-------------------------------------------------------------------------*/
static void
initsieve(curr_poly_t *c, root_sieve_t *rs, 
		double limit, rs_bound_t *rs_bounds)
{
	int exc;
	int p;
	int i, j, co[6], md, mp, len, pk, fi, mi, l, k;
	int sieveevaldiff[6], fis[6];
	double dp, dpk, dlp, dlog;
	int J, step, help;
	int nsubsieves;
	int sievelen;
	primelist_t *pl = rs->pl;

	for (i = 0; i < 6; i++) {
		mpz_out_str(stdout, 10, c->gmp_a[5 - i]);
		printf(" ");
	}
	printf("\n");
	mpz_out_str(stdout, 10, c->gmp_p);
	printf(" ");
	mpz_out_str(stdout, 10, c->gmp_d);
	printf("\n\n");

	sievelen = SIEVELEN;
	if (2 * (rs_bounds->xmax - rs_bounds->xmin) < sievelen) {
		sievelen /= 2;
		exc = sievelen - (rs_bounds->xmax - rs_bounds->xmin);
		exc /= 2;
		rs_bounds->xmin -= exc;
		rs_bounds->xmax = sievelen - 1 + rs_bounds->xmin;
		nsubsieves = 1;
	}
	else {
		nsubsieves = (rs_bounds->xmax - rs_bounds->xmin) / 
						sievelen + 1;
		if ((nsubsieves & 1) == 0)
			nsubsieves++;	/* odd number */
		exc = nsubsieves * sievelen - 
				(rs_bounds->xmax - rs_bounds->xmin);
		exc /= 2;
		rs_bounds->xmin -= exc;
		rs_bounds->xmax = nsubsieves * sievelen - 1 + 
					rs_bounds->xmin;
	}
	len = 0;
	for (k = 0; k < NAFF_PRIMES; k++) {
		p = primes[k];
		dp = (double) p;
		dlp = (double) (log(dp));
		pk = p;
		while (pk < MAX_PRIME_AFF) {
			dpk = (double) pk;
			dlog = dlp / dpk * dp / (dp + 1.);
			limit += dlog;
			for (i = 0; i <= 5; i++)
				co[i] = mpz_mod_ui(c->gmp_help1, 
						c->gmp_a[i], pk);
			md = mpz_mod_ui(c->gmp_help1, c->gmp_d, pk);
			mp = mpz_mod_ui(c->gmp_help1, c->gmp_p, pk);
			if (pk < 8) {
				for (i = 0; i < pk; i++) {
					fi = 0;
					for (j = 0; j <= 5; j++) {
						fi *= i;
						fi += co[5 - j];
						fi %= pk;
					}
					mi = md - ((i * mp) % pk);
					if (mi < 0)
						mi += pk;
					if (mi == 0) {
						if (fi == 0) {
							limit -= dlog;	/* p^k | fi+j*(i*P-d) for all j */
						}
					}
					else {
						l = 0;
						step = pk;
						while ((mi % p == 0) &&
						       (fi % p == 0)) {
							l++;
							mi /= p;
							fi /= p;
							step /= p;
							if (l > 32)
								Schlendrian
									("initsieve\n");
						}
						if (mi % p) {
							J = divmod(fi, mi,
								   step);
							pl[len].primepower = pk;
							pl[len].rooti = i;
							pl[len].J = J;
							pl[len].step = step;
							help = (J - 
							     rs_bounds->ymin * 
							     i - rs_bounds->xmin) % step;
							if (help < 0)
								help += step;
							pl[len].start = help;
							help = (-i +
								nsubsieves *
								sievelen) %
								step;
							if (help < 0)
								help += step;
							pl[len].lineinc = help;
							pl[len].value = dlog;
							len++;
						}	/* otherwise no solutions */
					}
				}
			}
			else {
				for (i = 0; i < 6; i++) {
					fi = 0;
					for (j = 0; j <= 5; j++) {
						fi *= i;
						fi += co[5 - j];
						fi %= pk;
					}
					fis[i] = fi;
				}
				for (i = 0; i < 6; i++) {
					sieveevaldiff[i] = fis[0];
					for (j = 0; j < 5 - i; j++)
						fis[j] = fis[j + 1] - fis[j];
				}
				for (i = 0; i < pk; i++) {
					mi = md - ((i * mp) % pk);
					if (mi < 0)
						mi += pk;
					fi = sieveevaldiff[0];
					for (j = 0; j < 5; j++)
						sieveevaldiff[j] +=
							sieveevaldiff[j + 1];
					for (j = 0; j < 5; j++)
						sieveevaldiff[j] %= pk;
					if (mi == 0) {
						if (fi == 0)
							limit -= dlog;	/* p^k | fi+j*(i-m) for all j */
					}
					else {
						l = 0;
						step = pk;
						while ((mi % p == 0) &&
						       (fi % p == 0)) {
							l++;
							mi /= p;
							fi /= p;
							step /= p;
							if (l > 32)
								Schlendrian
									("initsieve\n");
						}
						if (mi % p) {
							J = divmod(fi, mi,
								   step);
							pl[len].primepower = pk;
							pl[len].rooti = i;
							pl[len].J = J;
							pl[len].step = step;
							help = (J - rs_bounds->ymin * i -
								rs_bounds->xmin) % step;
							if (help < 0)
								help += step;
							pl[len].start = help;
							help = (-i +
								nsubsieves *
								sievelen) %
								step;
							if (help < 0)
								help += step;
							pl[len].lineinc = help;
							pl[len].value = dlog;
							len++;
						}	/* otherwise no solutions */
					}
				}
			}
			pk *= p;
		}
	}

	for (i = 0; i < len; i++) {
		pl[i].v = (unsigned short) (pl[i].value * 1000. + 0.5);
		if (pl[i].v == 0)
			pl[i].v = 1;
	}
	rs->default_cut = (unsigned short) (limit * 1000. + 0.5);
	if (rs->default_cut == 0)
		rs->default_cut = 1;

	rs->primelistlen = len;
	rs->limit = limit;
	rs->nsubsieves = nsubsieves;
	rs->sievelen = sievelen;
}

/*-------------------------------------------------------------------------*/
static void
prepare_sieve(root_sieve_t *rs)
{
	int l, i0, i1, i, k, p, pk, ii, st;
	unsigned short *us_ptr;
	primelist_t *pl = rs->pl;
	int primelistlen = rs->primelistlen;
	unsigned int **prep_p_begin = rs->prep_p_begin;

	l = 0;
	i0 = 0;
	for (k = 0; k < NAFF_PRIMES; k++) {
		p = primes[k];
		for (i = i0; i < primelistlen; i++)
			if (pl[i].primepower == p)
				break;
		if (i == primelistlen)
			continue;	/* p not in pl-list */
		i0 = i;
		for (i = i0 + 1; i < primelistlen; i++)
			if (pl[i].primepower % p)
				break;
		i1 = i;		/* [i0,i1[: p-part of pl-list */
		pk = p;
		for (i = i0; i < i1; i++)
			if (pl[i].primepower > pk)
				pk = pl[i].primepower;	/* maximal primepower */

		for (i = 0; i < pk; i++)
			prep_p_begin[l][i] = 0;
		us_ptr = (unsigned short *) (prep_p_begin[l]);
		for (i = i0; i < i1; i++) {
			st = pl[i].step;
			ii = pl[i].start;
			while (ii < 2 * pk) {
				us_ptr[ii] += pl[i].v;
				ii += st;
			}
		}

		for (i = 0; i < pk; i++)
			prep_p_begin[l][i + pk] = prep_p_begin[l][i];

		rs->prep_p[l] = prep_p_begin[l];
		rs->prep_p_len[l] = pk;
		l++;
	}
	rs->prep_len = l;
}

/*-------------------------------------------------------------------------*/
static void
finish_sieve(root_sieve_t *rs)
{
	int i, ind, add;
	int sievelen = rs->sievelen;
	int nsubsieves = rs->nsubsieves;

	for (i = 0; i < rs->primelistlen; i++) {
		primelist_t *pl = rs->pl + i;
		ind = pl->start;
		add = pl->step;
		ind += add;
		ind -= ((sievelen * nsubsieves) % add);
		ind %= add;
		pl->start = ind;
	}
}

/*-------------------------------------------------------------------------*/
static void
sieve_new(root_sieve_t *rs)
{
	int i, len2;
	unsigned int *ul_sv;
	int len = rs->sievelen;
	int *prep_p_len = rs->prep_p_len;
	unsigned int **prep_p_begin = rs->prep_p_begin;
	unsigned int **prep_p = rs->prep_p;
	unsigned short *sievearray = rs->sievearray;

#ifndef HAVE_ASM_INTEL
	int ind, end;
	unsigned int *ptr, *ptrbegin, *ptrend;
#endif

	len2 = len / 2;
	ul_sv = (unsigned int *)sievearray;
	memset(sievearray, 0, len * sizeof(unsigned short));
	for (i = 0; i < rs->prep_len; i++) {
#ifdef HAVE_ASM_INTEL
		asm_root_sieve8(&(prep_p[i]), prep_p_begin[i], prep_p_len[i],
				ul_sv, len / 4);
#else
		ptr = prep_p[i];
		ptrbegin = prep_p_begin[i];
		ptrend = ptrbegin + prep_p_len[i];
		ind = 0;
		end = len2 - prep_p_len[i];
		while (ptr < ptrend)
			ul_sv[ind++] += *ptr++;
		while (ind < end) {
			for (ptr = ptrbegin; ptr < ptrend; ptr++, ind++)
				ul_sv[ind] += *ptr;
		}
		for (ptr = ptrbegin; ind < len2; ptr++, ind++)
			ul_sv[ind] += *ptr;
		prep_p[i] = ptr;
#endif
	}
}

/*-------------------------------------------------------------------------*/
static void
advancesieve(root_sieve_t *rs)
{
	int i, s, add;

	for (i = 0; i < rs->primelistlen; i++) {
		primelist_t *pl = rs->pl + i;
		s = pl->start + pl->lineinc;
		add = pl->step;
		if (s >= add)
			s -= add;
		pl->start = s;
	}
}

/*-------------------------------------------------------------------------*/
void
root_sieve_run(curr_poly_t *c, double log_max_norm_2, 
		poly_stage2_t *data, root_sieve_t *rs,
		assess_t *assess, stage2_stat_t *stats,
		double skewness, double pol_norm,
		double alpha_proj)
{
	int x, y, i, j;
	double lim0;
	double sk, v, dbl_sv[6];
	double dbl_p;
	double dbl_d;
	rs_bound_t rs_bounds;
	unsigned short cut;

	profile_start(PROF_SIEVE_ALL);

	profile_start(PROF_INIT_SIEVE);

	dbl_p = mpz_get_d(c->gmp_p);
	dbl_d = mpz_get_d(c->gmp_d);
	for (i = 0; i <= 5; i++)
		dbl_sv[i] = mpz_get_d(c->gmp_a[i]);

	compute_sieving_area(data->max_norm_1, dbl_sv, 
				dbl_d, &rs_bounds, skewness,
				alpha_proj);
	lim0 = log(pol_norm);
	initsieve(c, rs, 
		lim0 - log_max_norm_2 + alpha_proj,
		&rs_bounds);
	if (verbose > 1) {
		printf("area: [%d,%d]x[%d,%d]\n", 
				rs_bounds.xmin, rs_bounds.xmax, 
				rs_bounds.ymin, rs_bounds.ymax);
		if (verbose > 2)
			printf("limit: %f ", rs->limit);
	}

	if (rs->limit < 0.) {	/* accept any (x,y) in sieving area */
		profile_stop(PROF_INIT_SIEVE);
		profile_start(PROF_EVAL);
		for (y = rs_bounds.ymin; y <= rs_bounds.ymax; y++) {
			for (x = rs_bounds.xmin; x <= rs_bounds.xmax; x++) {
				check(x, y, c, data, assess, stats, skewness);
			}
		}
		profile_stop(PROF_EVAL);
		profile_stop(PROF_SIEVE_ALL);
		return;
	}
	dbl_sv[2] += ((double) rs_bounds.ymin) * dbl_p;
	dbl_sv[1] -= ((double) rs_bounds.ymin) * dbl_d;
	sk = skewness;
	cut = rs->default_cut;

	printf("norm: %.3e  alpha_proj: %.3f  skewness: %.3f\n", 
			pol_norm, alpha_proj, sk);
	printf("limit: %f, lim0: %f   cut: %d\n", 
			rs->limit, lim0, (int) cut);

	profile_stop(PROF_INIT_SIEVE);
	for (y = rs_bounds.ymin; y <= rs_bounds.ymax; y++) {
		sk = find_best_skewness(dbl_sv, sk);
		v = ifs(dbl_sv, sk);
		dbl_sv[2] += dbl_p;
		dbl_sv[1] -= dbl_d;
		v = log(v) / 2.;
		cut = (unsigned short) ((rs->limit + v - lim0) * 1000 + 0.5);
		if (cut == 0)
			cut = 1;
		if (verbose > 2) {
			printf("y: %d, lim: %f cut: %d,  sk: %.2f\n", 
					y, v, (int) cut, sk);
		}

		x = rs_bounds.xmin;

		profile_start(PROF_INIT_SIEVE);
		prepare_sieve(rs);
		profile_stop(PROF_INIT_SIEVE);

		for (i = 0; i < rs->nsubsieves; i++) {

			profile_start(PROF_SIEVE);
			sieve_new(rs);
			profile_stop(PROF_SIEVE);

			profile_start(PROF_EVAL);
			for (j = 0; j < rs->sievelen; j++) {
				if (rs->sievearray[j] >= cut) {
					if (verbose > 2) {
						printf("sv: %d, cut: %d  "
							"at x: %d  ", 
							rs->sievearray[j], 
							(int)cut, x + j);
					}
					check(x + j, y, c, data, 
						assess, stats, skewness);
				}
			}
			x += rs->sievelen;
			profile_stop(PROF_EVAL);
		}
		finish_sieve(rs);
		if (y < rs_bounds.ymax)
			advancesieve(rs);
		if (verbose > 1)
			printf(".");
	}
	profile_stop(PROF_SIEVE_ALL);
}

/*-------------------------------------------------------------------------*/
void root_sieve_init(root_sieve_t *rs)
{
	int i;
	memset(rs, 0, sizeof(root_sieve_t));

	rs->st_alpha = compute_standard_alpha();
	rs->pl = (primelist_t *)xmalloc(MAX_PRIME_AFF * MAX_PRIME_AFF *
					sizeof(primelist_t));
	rs->sievearray = (unsigned short *)xmalloc(SIEVELEN *
						sizeof(unsigned short));
	for (i = 0; i < NAFF_PRIMES; i++) {
		rs->prep_p_begin[i] = (unsigned int *)xmalloc(
						2 * MAX_PRIME_AFF *
						sizeof(unsigned int));
	}
}

/*-------------------------------------------------------------------------*/
void root_sieve_free(root_sieve_t *rs)
{
	int i;
	memset(rs, 0, sizeof(root_sieve_t));

	free(rs->pl);
	free(rs->sievearray);
	for (i = 0; i < NAFF_PRIMES; i++)
		free(rs->prep_p_begin[i]);
}
