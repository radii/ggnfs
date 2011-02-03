#include "stage1_impl.h"

/*
The following functions check whether given (a5,p,d) is useful, i. e.
whether we can find a4,...,a0 and a skewness such that the de-skewed
sup-norm is less than normmax.
This is called very seldom and the following implementation is slow.
*/

/*------------------------------------------------------------------------*/
static double
check_a2_eval(double b5, double b4, double b3, double b2, double sk,
		bounds_t *bounds)
{
	double sk_sq = sqrt(sk), m, m0;

	m = fabs(b5 * sk * sk * sk_sq);
	m0 = fabs(b4 * sk * sk_sq);
	if (m0 > m)
		m = m0;
	m0 = fabs(b3 * sk_sq);
	if (m0 > m)
		m = m0;
	m0 = fabs(b2 / sk_sq);
	if (m0 > m)
		m = m0;
	return m / bounds->norm_max;
}

/*------------------------------------------------------------------------*/
#define TRANSLATE k=(int)dk; dh=(double)k; d5=c5; d4=c4+5.*dh*c5; \
                  d3=c3+4.*dh*c4+10.*dh*dh*c5; \
                  d2=c2+3.*dh*c3+6.*dh*dh*c4+10.*dh*dh*dh*c5

static int
check_a2(double b5, double b4, double b3, double b2, double sk,
		bounds_t *bounds)
		  /* very bad!!!! */
{
	int i, k;
	double c2 = b2, c3 = b3, c4 = b4, c5 = b5, d2, d3, d4, d5;
	double s = sk, s0, ds, dk, dh, val, val0;

	ds = s / 10.;
	dk = ds;
	val = check_a2_eval(c5, c4, c3, c2, s, bounds);
	for (i = 0; i < 1000; i++) {
/* skewness */
		if (ds > 2.) {
			s0 = s + ds;
			if (s0 > bounds->skewness_max)
				s0 = bounds->skewness_max;
			val0 = check_a2_eval(c5, c4, c3, c2, s0, bounds);
			if (val0 < val) {
				s = s0;
				ds *= 1.4;
				val = val0;
			}
			else {
				s0 = s - ds;
				if (s0 < bounds->skewness_min)
					s0 = bounds->skewness_min;
				val0 = check_a2_eval(c5, c4, c3, c2, s0, bounds);
				if (val0 < val) {
					s = s0;
					ds *= 1.4;
					val = val0;
				}
				else
					ds /= 1.4;
			}
			if (val < 1.)
				return 1;
		}
		/* translation */
		if (dk > 2.) {
			TRANSLATE;
			val0 = check_a2_eval(d5, d4, d3, d2, s, bounds);
			if (val0 < val) {
				c2 = d2;
				c3 = d3;
				c4 = d4;
				c5 = d5;
				dk *= 1.4;
				val = val0;
			}
			else {
				dk = -dk;
				TRANSLATE;
				dk = -dk;
				val0 = check_a2_eval(d5, d4, d3, d2, s, bounds);
				if (val0 < val) {
					c2 = d2;
					c3 = d3;
					c4 = d4;
					c5 = d5;
					dk *= 1.4;
					val = val0;
				}
				else
					dk /= 1.4;
			}
			if (val < 1.)
				return 1;
		}
		if ((ds < 10.) && (dk < 10.))
			break;
	}
	if (verbose > 1)
		if (val < 100.)
			printf(".");
	return 0;
}

/*------------------------------------------------------------------------*/
int
check_pol(curr_poly_t *curr, bounds_t *bounds, stage1_stat_t *stats)
{
	double sk, s;
	double dbl_a5, dbl_a4, dbl_a3, dbl_a2;
	mpz_t gmp_a0, gmp_a1, gmp_a2, gmp_a3, gmp_a4;
	int ok;

	profile_start(PROF_CHECK_POLY);
	stats->n_polexpand++;
	mpz_init(gmp_a0);
	mpz_init(gmp_a1);
	mpz_init(gmp_a2);
	mpz_init(gmp_a3);
	mpz_init(gmp_a4);

	/* compute coefficients */
	mpz_mul(curr->gmp_help4, curr->gmp_d, curr->gmp_d);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_help4);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_d);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_a5);
	mpz_sub(curr->gmp_help3, curr->gmp_N, curr->gmp_help4);
	mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help1, curr->gmp_help3, curr->gmp_prod);
	if (mpz_sgn(curr->gmp_help1))
		complain("pol-dev5\n");

	if (!mpz_invert(curr->gmp_help2, curr->gmp_d, curr->gmp_prod))
		complain("pol-dev.i\n");
	mpz_mul(curr->gmp_help4, curr->gmp_help2, curr->gmp_help2);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_help4);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_help3);
	mpz_fdiv_r(gmp_a4, curr->gmp_help4, curr->gmp_prod);
	mpz_mul(curr->gmp_help4, curr->gmp_d, curr->gmp_d);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_help4);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, gmp_a4);
	mpz_sub(curr->gmp_help3, curr->gmp_help3, curr->gmp_help4);
	mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help1, curr->gmp_help3, curr->gmp_prod);
	if (mpz_sgn(curr->gmp_help1))
		complain("pol-dev4\n");

	mpz_mul(curr->gmp_help4, curr->gmp_d, curr->gmp_d);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_d);
	mpz_fdiv_q(gmp_a3, curr->gmp_help3, curr->gmp_help4);
	mpz_fdiv_q(gmp_a3, gmp_a3, curr->gmp_prod);
	mpz_mul(gmp_a3, gmp_a3, curr->gmp_prod);
	mpz_mul(curr->gmp_help4, curr->gmp_help2, curr->gmp_help2);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_help2);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_help3);
	mpz_fdiv_r(curr->gmp_help1, curr->gmp_help4, curr->gmp_prod);
	mpz_add(gmp_a3, gmp_a3, curr->gmp_help1);
	mpz_mul(curr->gmp_help4, curr->gmp_d, curr->gmp_d);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_d);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, gmp_a3);
	mpz_sub(curr->gmp_help3, curr->gmp_help3, curr->gmp_help4);
	mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help1, curr->gmp_help3, curr->gmp_prod);
	if (mpz_sgn(curr->gmp_help1))
		complain("pol-dev3\n");

	mpz_mul(curr->gmp_help4, curr->gmp_d, curr->gmp_d);
	mpz_fdiv_q(gmp_a2, curr->gmp_help3, curr->gmp_help4);
	mpz_fdiv_q(gmp_a2, gmp_a2, curr->gmp_prod);
	mpz_mul(gmp_a2, gmp_a2, curr->gmp_prod);
	mpz_mul(curr->gmp_help4, curr->gmp_help2, curr->gmp_help2);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, curr->gmp_help3);
	mpz_fdiv_r(curr->gmp_help1, curr->gmp_help4, curr->gmp_prod);
	mpz_add(gmp_a2, gmp_a2, curr->gmp_help1);
	mpz_mul(curr->gmp_help4, curr->gmp_d, curr->gmp_d);
	mpz_mul(curr->gmp_help4, curr->gmp_help4, gmp_a2);
	mpz_sub(curr->gmp_help3, curr->gmp_help3, curr->gmp_help4);
	mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help1, curr->gmp_help3, curr->gmp_prod);
	if (mpz_sgn(curr->gmp_help1))
		complain("pol-dev2\n");

	mpz_fdiv_q(gmp_a1, curr->gmp_help3, curr->gmp_d);
	mpz_fdiv_q(gmp_a1, gmp_a1, curr->gmp_prod);
	mpz_mul(gmp_a1, gmp_a1, curr->gmp_prod);
	mpz_mul(curr->gmp_help4, curr->gmp_help3, curr->gmp_help2);
	mpz_fdiv_r(curr->gmp_help1, curr->gmp_help4, curr->gmp_prod);
	mpz_add(gmp_a1, gmp_a1, curr->gmp_help1);
	mpz_mul(curr->gmp_help4, curr->gmp_d, gmp_a1);
	mpz_sub(curr->gmp_help3, curr->gmp_help3, curr->gmp_help4);
	mpz_fdiv_qr(curr->gmp_help3, curr->gmp_help1, curr->gmp_help3, curr->gmp_prod);
	if (mpz_sgn(curr->gmp_help1))
		complain("pol-dev1\n");

	mpz_set(gmp_a0, curr->gmp_help3);

	mpz_fdiv_qr(curr->gmp_help1, gmp_a3, gmp_a3, curr->gmp_d);
	mpz_add(curr->gmp_help2, gmp_a3, gmp_a3);
	if (mpz_cmp(curr->gmp_d, curr->gmp_help2) < 0) {
		mpz_sub(gmp_a3, gmp_a3, curr->gmp_d);
		mpz_add_ui(curr->gmp_help1, curr->gmp_help1, 1);
	}
	mpz_mul(curr->gmp_help1, curr->gmp_help1, curr->gmp_prod);
	mpz_add(gmp_a4, gmp_a4, curr->gmp_help1);

	mpz_fdiv_qr(curr->gmp_help1, gmp_a2, gmp_a2, curr->gmp_d);
	mpz_add(curr->gmp_help2, gmp_a2, gmp_a2);
	if (mpz_cmp(curr->gmp_d, curr->gmp_help2) < 0) {
		mpz_sub(gmp_a2, gmp_a2, curr->gmp_d);
		mpz_add_ui(curr->gmp_help1, curr->gmp_help1, 1);
	}
	mpz_mul(curr->gmp_help1, curr->gmp_help1, curr->gmp_prod);
	mpz_add(gmp_a3, gmp_a3, curr->gmp_help1);

	if (verbose > 2) {
		printf("polynomial:\n");
		mpz_out_str(stdout, 10, curr->gmp_a5); printf("\n");
		mpz_out_str(stdout, 10, gmp_a4); printf("\n");
		mpz_out_str(stdout, 10, gmp_a3); printf("\n");
		mpz_out_str(stdout, 10, gmp_a2); printf("\n");
		mpz_out_str(stdout, 10, gmp_a1); printf("\n");
		mpz_out_str(stdout, 10, gmp_a0); printf("\n\n");
	}

	dbl_a5 = mpz_get_d(curr->gmp_a5);
	dbl_a4 = mpz_get_d(gmp_a4);
	dbl_a3 = mpz_get_d(gmp_a3);
	dbl_a2 = mpz_get_d(gmp_a2);

	sk = bounds->skewness_max;
	s = pow(bounds->norm_max / fabs(dbl_a4), 2. / 3.);
	if (s < sk)
		sk = s;
	s = pow(bounds->norm_max / fabs(dbl_a3), 2.);
	if (s < sk)
		sk = s;

	if (verbose > 2) {
		printf("quality: %f %f %f %f\n",
		       dbl_a5 / bounds->norm_max * sk * sk * sqrt(sk),
		       dbl_a4 / bounds->norm_max * sk * sqrt(sk),
		       dbl_a3 / bounds->norm_max * sqrt(sk),
		       dbl_a2 / bounds->norm_max / sqrt(sk));
	}
	if (fabs(dbl_a2 / bounds->norm_max / sqrt(sk)) <= 1.) {
		ok = 1;
	}
	else {
		stats->n_poltranslate++;
		ok = check_a2(dbl_a5, dbl_a4, dbl_a3, dbl_a2, sk, bounds);
	}

	mpz_clear(gmp_a0);
	mpz_clear(gmp_a1);
	mpz_clear(gmp_a2);
	mpz_clear(gmp_a3);
	mpz_clear(gmp_a4);
	profile_stop(PROF_CHECK_POLY);
	if (ok)
		stats->n_survivors++;
	return ok;
}
