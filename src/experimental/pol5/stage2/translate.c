#include "stage2_impl.h"

/*-------------------------------------------------------------------------*/
#ifdef ULL_NO_UL
#define BITS_PER_ULONG 32
static void mpz_set_ull(mpz_t targ, uint64_t src)
{
  mpz_set_ui(targ,(uint32_t)(src>>32));
  mpz_mul_2exp(targ, targ, BITS_PER_ULONG);
  mpz_add_ui(targ, targ, (uint32_t)(src&0xFFFFFFFF));
}

static void mpz_set_sll(mpz_t x, int64_t src)
{
  if (src < 0) {
    mpz_set_ull(x, (uint64_t) (-src));
    mpz_neg(x, x);
  } else
    mpz_set_ull(x, (uint64_t) src);
}
#endif

/*-------------------------------------------------------------------------*/
static void
translate_dbl(double *dtarg, double *dsrc, int k)
{
	int i;
	double d, dk;

	for (i = 0; i < 6; i++)
		dtarg[i] = dsrc[i];
	dk = (double) (-k);
	for (i = 0; i < 5; i++)
		dtarg[i] += (double) (i + 1) * dsrc[i + 1] * dk;
	d = dk * dk;
	dtarg[0] += d * dsrc[2];
	dtarg[1] += 3 * d * dsrc[3];
	dtarg[2] += 6 * d * dsrc[4];
	dtarg[3] += 10 * d * dsrc[5];
	d *= dk;
	dtarg[0] += d * dsrc[3];
	dtarg[1] += 4 * d * dsrc[4];
	dtarg[2] += 10 * d * dsrc[5];
	d *= dk;
	dtarg[0] += d * dsrc[4];
	dtarg[1] += 5 * d * dsrc[5];
	d *= dk;
	dtarg[0] += d * dsrc[5];
}

/*-------------------------------------------------------------------------*/
static void
translate_gmp(curr_poly_t *c, mpz_t * gmp_z, 
		mpz_t * lin, mpz_t gmp_mz, int k)
{
	if (k >= 0) {
		mpz_set_ui(c->gmp_help1, k);
		mpz_neg(c->gmp_help1, c->gmp_help1);
	}
	else {
		mpz_set_ui(c->gmp_help1, -k);	/* help1=-k */
	}
	mpz_set(c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[1], c->gmp_help3);
	mpz_add(gmp_z[0], gmp_z[0], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[2], c->gmp_help3);
	mpz_add(gmp_z[0], gmp_z[0], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[3], c->gmp_help3);
	mpz_add(gmp_z[0], gmp_z[0], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[4], c->gmp_help3);
	mpz_add(gmp_z[0], gmp_z[0], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[5], c->gmp_help3);
	mpz_add(gmp_z[0], gmp_z[0], c->gmp_help2);
/* a0<-a0-a1*k+a2*k^2-a3*k^3+a4*k^4-a5*k^5 */

	mpz_set(c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[2], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, 2);
	mpz_add(gmp_z[1], gmp_z[1], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[3], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, 3);
	mpz_add(gmp_z[1], gmp_z[1], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[4], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, 4);
	mpz_add(gmp_z[1], gmp_z[1], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[5], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, 5);
	mpz_add(gmp_z[1], gmp_z[1], c->gmp_help2);
/* a1<-a1-2*a2*k+3*a3*k^2-4*a4*k^3+5*a5*k^4 */

	mpz_set(c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[3], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, 3);
	mpz_add(gmp_z[2], gmp_z[2], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[4], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, 6);
	mpz_add(gmp_z[2], gmp_z[2], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[5], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, 10);
	mpz_add(gmp_z[2], gmp_z[2], c->gmp_help2);
/* a2<-a2-3*a3*k+6*a4*k^2-10*a5*k^3 */

	mpz_set(c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[4], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, 4);
	mpz_add(gmp_z[3], gmp_z[3], c->gmp_help2);
	mpz_mul(c->gmp_help3, c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[5], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, 10);
	mpz_add(gmp_z[3], gmp_z[3], c->gmp_help2);
/* a3<-a3-4*a4*k+10*a5*k^2 */

	mpz_set(c->gmp_help3, c->gmp_help1);
	mpz_mul(c->gmp_help2, gmp_z[5], c->gmp_help3);
	mpz_mul_ui(c->gmp_help2, c->gmp_help2, 5);
	mpz_add(gmp_z[4], gmp_z[4], c->gmp_help2);
/* a4<-a4-5*a5*k */

	mpz_sub(gmp_mz, gmp_mz, c->gmp_help1);
/* m<-m+k */

	mpz_mul(c->gmp_help2, lin[1], c->gmp_help1);
	mpz_add(lin[0], lin[0], c->gmp_help2);
/* lin0<-lin0-lin1*k */
}

/*-------------------------------------------------------------------------*/
void
optimize_1(curr_poly_t *c,
		double *skewness, double *pol_norm,
		double *alpha_proj)
{
	int dk, i, niter;
	int64_t di1, di0;
	double dbl_a[6];
	double dbl_a0[6];
	double dbl_d, dbl_p;
	double value, v0;
	double s, s0, ds, d;

	dk = 16;
	niter = 0;
	for (i = 0; i < 6; i++)
		dbl_a[i] = mpz_get_d(c->gmp_a[i]);
	dbl_d = mpz_get_d(c->gmp_d);
	s = 1.;
	if (dbl_a[4] != 0)
		s = sqrt(fabs(dbl_a[2] / dbl_a[4]));
	s0 = 1.;
	if (dbl_a[3] != 0)
		s0 = fabs(dbl_a[2] / dbl_a[3]);
	if (s0 > s)
		s = s0;
	ds = 2.;
	di1 = (int64_t) (dbl_a[1] / dbl_d);
	di0 = (int64_t) (dbl_a[0] / dbl_d);
	value = ifs(dbl_a, s);
	while (1) {
		if (niter > 10000) {
			fprintf(stderr, "too many iterations in optimize_1\n");
			break;
		}
		if ((dk < 2) && (di1 < 2) && (di0 < 2) && (ds < 0.001))
			break;
/* skewness */
		s0 = s * (1. + ds);
		v0 = ifs(dbl_a, s0);
		if (v0 < value) {
			s = s0;
			value = v0;
			ds *= 1.1;
		}
		else {
			s0 = s / (1. + ds);
			v0 = ifs(dbl_a, s0);
			if (v0 < value) {
				s = s0;
				value = v0;
				ds *= 1.1;
			}
			else
				ds /= 1.1;
		}
/* translation */
		translate_dbl(dbl_a0, dbl_a, dk);
		v0 = ifs(dbl_a0, s);
		if (v0 < value) {
			translate_gmp(c, c->gmp_a, 
					c->gmp_lina, c->gmp_m, dk);
			value = v0;
			dk = (int) (1 + 1.1 * (double) dk);
		}
		else {
			translate_dbl(dbl_a0, dbl_a, -dk);
			v0 = ifs(dbl_a0, s);
			if (v0 < value) {
				translate_gmp(c, c->gmp_a, 
						c->gmp_lina, c->gmp_m, -dk);
				value = v0;
				dk = (int) (1 + 1.1 * (double) dk);
			}
			else {
				dk = (int) (0.9 * (double) dk - 1);
				if (dk < 1)
					dk = 1;
			}
		}
		mpz_set(c->gmp_p, c->gmp_lina[1]);
		mpz_neg(c->gmp_d, c->gmp_lina[0]);
		for (i = 0; i < 6; i++)
			dbl_a[i] = mpz_get_d(c->gmp_a[i]);
		dbl_d = mpz_get_d(c->gmp_d);
		dbl_p = mpz_get_d(c->gmp_p);
/* "rotation" with i1*x*(px-m) */
		d = (double) di1;
		dbl_a[2] += d * dbl_p;
		dbl_a[1] -= d * dbl_d;
		v0 = ifs(dbl_a, s);
		if (v0 < value) {
			mpz_set_sll(c->gmp_help1, di1);
			mpz_mul(c->gmp_help2, c->gmp_help1, c->gmp_p);
			mpz_add(c->gmp_a[2], c->gmp_a[2], c->gmp_help2);
			mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_d);
			mpz_sub(c->gmp_a[1], c->gmp_a[1], c->gmp_help1);
			value = v0;
			di1 = (int64_t) (1 + 1.1 * (double) di1);
		}
		else {
			dbl_a[2] -= 2 * d * dbl_p;
			dbl_a[1] += 2 * d * dbl_d;
			v0 = ifs(dbl_a, s);
			if (v0 < value) {
				mpz_set_sll(c->gmp_help1, -di1);
				mpz_mul(c->gmp_help2, c->gmp_help1, c->gmp_p);
				mpz_add(c->gmp_a[2], c->gmp_a[2], c->gmp_help2);
				mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_d);
				mpz_sub(c->gmp_a[1], c->gmp_a[1], c->gmp_help1);
				value = v0;
				di1 = (int64_t) (1 + 1.1 * (double) di1);
			}
			else {
				dbl_a[2] += d * dbl_p;
				dbl_a[1] -= d * dbl_d;	/* set it back */
				di1 = (int64_t) (0.9 * (double) di1 - 1);
				if (di1 < 1)
					di1 = 1;
			}
		}
/* "rotation" with i0*(px-m) */
		d = (double) di0;
		dbl_a[1] += d * dbl_p;
		dbl_a[0] -= d * dbl_d;
		v0 = ifs(dbl_a, s);
		if (v0 < value) {
			mpz_set_sll(c->gmp_help1, di0);
			mpz_mul(c->gmp_help2, c->gmp_help1, c->gmp_p);
			mpz_add(c->gmp_a[1], c->gmp_a[1], c->gmp_help2);
			mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_d);
			mpz_sub(c->gmp_a[0], c->gmp_a[0], c->gmp_help1);
			value = v0;
			di0 = (int64_t) (1 + 1.1 * (double) di0);
		}
		else {
			dbl_a[1] -= 2 * d * dbl_p;
			dbl_a[0] += 2 * d * dbl_d;
			v0 = ifs(dbl_a, s);
			if (v0 < value) {
				mpz_set_sll(c->gmp_help1, -di0);
				mpz_mul(c->gmp_help2, c->gmp_help1, c->gmp_p);
				mpz_add(c->gmp_a[1], c->gmp_a[1], c->gmp_help2);
				mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_d);
				mpz_sub(c->gmp_a[0], c->gmp_a[0], c->gmp_help1);
				value = v0;
				di0 = (int64_t) (1 + 1.1 * (double) di0);
			}
			else {
				dbl_a[1] += d * dbl_p;
				dbl_a[0] -= d * dbl_d;	/* set it back */
				di0 = (int64_t) (0.9 * (double) di0 - 1);
				if (di0 < 1)
					di0 = 1;
			}
		}
		niter++;
	}

	*skewness = s;
	*pol_norm = sqrt(ifs(dbl_a, s));
	*alpha_proj = compute_proj_alpha(c);

	if (verbose > 2) {
		int i;

		printf("optimize 1\nskewness: %.2f norm: %.4e "
			"alpha %.2f\npol0: ",
		       *skewness, *pol_norm, *alpha_proj);
		for (i = 5; i >= 0; i--) {
			mpz_out_str(stdout, 10, c->gmp_a[i]);
			printf(" ");
		}
		printf("\npol1: ");
		mpz_out_str(stdout, 10, c->gmp_lina[1]);
		printf(" ");
		mpz_out_str(stdout, 10, c->gmp_lina[0]);
		printf("\n\n");
	}
}

/*-------------------------------------------------------------------------*/
void
optimize_2(curr_poly_t *c, double skewness, 
		double *new_skewness, double *norm_ptr)
{
	int dk, i, niter;
	double dbl_b[6];
	double dbl_b0[6];
	double value, v0;
	double s, s0, ds;

	dk = 16;
	ds = 2.;
	niter = 0;
	for (i = 0; i < 6; i++)
		dbl_b[i] = mpz_get_d(c->gmp_b[i]);
	s = skewness;
	value = ifs(dbl_b, s);
	while (1) {
		if (niter > 10000) {
			fprintf(stderr, "too many iterations in optimize_2\n");
			break;
		}
		if ((dk < 2) && (ds < 1.001))
			break;
/* skewness */
		s0 = s * ds;
		v0 = ifs(dbl_b, s0);
		if (v0 < value) {
			s = s0;
			value = v0;
			ds *= 1.1;
		}
		else {
			s0 = s / ds;
			v0 = ifs(dbl_b, s0);
			if (v0 < value) {
				s = s0;
				value = v0;
				ds *= 1.1;
			}
			else
				ds = 1. + (ds - 1.) / 1.1;
		}
/* translation */
		translate_dbl(dbl_b0, dbl_b, dk);
		v0 = ifs(dbl_b0, s);
		if (v0 < value) {
			translate_gmp(c, c->gmp_b, c->gmp_linb, c->gmp_mb, dk);
			value = v0;
			dk = (int) (1 + 1.1 * (double) dk);
		}
		else {
			translate_dbl(dbl_b0, dbl_b, -dk);
			v0 = ifs(dbl_b0, s);
			if (v0 < value) {
				translate_gmp(c, c->gmp_b, c->gmp_linb, 
						c->gmp_mb, -dk);
				value = v0;
				dk = (int) (1 + 1.1 * (double) dk);
			}
			else {
				dk = (int) (0.9 * (double) dk - 1);
				if (dk < 1)
					dk = 1;
			}
		}
		for (i = 0; i < 6; i++)
			dbl_b[i] = mpz_get_d(c->gmp_b[i]);
		niter++;
	}

	*new_skewness = s;
	*norm_ptr = sqrt(value);
	if (verbose > 3) {
		int i;

		printf("optimize 2\nskewness: %.2f norm: %.4e\npol0: ", 
				*new_skewness, *norm_ptr);
		for (i = 5; i >= 0; i--) {
			mpz_out_str(stdout, 10, c->gmp_b[i]);
			printf(" ");
		}
		printf("\npol1: ");
		mpz_out_str(stdout, 10, c->gmp_linb[1]);
		printf(" ");
		mpz_out_str(stdout, 10, c->gmp_linb[0]);
		printf("\n\n");
	}
}

/*-------------------------------------------------------------------------*/
void
optimize_3(curr_poly_t *c, poly_stage2_t *data, assess_t *assess,
		stage2_stat_t *stats, double skewness, double *new_skewness, 
		double *norm_ptr, double *eptr, double *alphaptr)
{
	int dk, i, niter;
	double dbl_b[6];
	double dbl_b0[6];
	double e, new_e, alpha;
	double s, s0, ds;
	double dbl_linb[2];

	alpha = *alphaptr;
	dk = 16;
	ds = 2.;
	niter = 0;
	for (i = 0; i < 6; i++)
		dbl_b[i] = mpz_get_d(c->gmp_b[i]);

	s = skewness;
	dbl_linb[1] = mpz_get_d(c->gmp_linb[1]);
	dbl_linb[0] = mpz_get_d(c->gmp_linb[0]);
	profile_start(PROF_MURPHY_E);
	murphy_e_score(&e, 5, dbl_b, 1, dbl_linb, 
			alpha, 0., s, 1000, assess);
	profile_stop(PROF_MURPHY_E);
	if (e < data->min_e) {
		*new_skewness = s;
		*eptr = e;
		*norm_ptr = sqrt(ifs(dbl_b, s));
		return;
	}
	while (1) {
		if (niter > 10000) {
			fprintf(stderr, "too many iterations in optimize_3\n");
			break;
		}
		if ((dk < 8) && (ds < 1.001))
			break;
/* skewness */
		s0 = s * ds;
		profile_start(PROF_MURPHY_E);
		murphy_e_score(&new_e, 5, dbl_b, 1, dbl_linb, 
				alpha, 0., s0, 1000, assess);
		profile_stop(PROF_MURPHY_E);
		if (new_e > e) {
			s = s0;
			e = new_e;
			ds *= 1.1;
		}
		else {
			s0 = s / ds;
			profile_start(PROF_MURPHY_E);
			murphy_e_score(&new_e, 5, dbl_b, 1, dbl_linb, 
					alpha, 0., s0, 1000, assess);
			profile_stop(PROF_MURPHY_E);
			if (new_e > e) {
				s = s0;
				e = new_e;
				ds *= 1.1;
			}
			else {
				ds = 1. + (ds - 1.) / 1.1;
			}
		}
/* translation */
		translate_dbl(dbl_b0, dbl_b, dk);
		dbl_linb[0] -= dk * dbl_linb[1];
		profile_start(PROF_MURPHY_E);
		murphy_e_score(&new_e, 5, dbl_b0, 1, dbl_linb, 
				alpha, 0., s, 1000, assess);
		profile_stop(PROF_MURPHY_E);
		dbl_linb[0] += dk * dbl_linb[1];
		if (new_e > e) {
			translate_gmp(c, c->gmp_b, 
					c->gmp_linb, c->gmp_mb, dk);
			dbl_linb[0] -= dk * dbl_linb[1];
			e = new_e;
			dk = (int) (1 + 1.1 * (double) dk);
		}
		else {
			dbl_linb[0] += dk * dbl_linb[1];
			translate_dbl(dbl_b0, dbl_b, -dk);
			profile_start(PROF_MURPHY_E);
			murphy_e_score(&new_e, 5, dbl_b0, 1, dbl_linb, 
					alpha, 0., s, 1000, assess);
			profile_stop(PROF_MURPHY_E);
			dbl_linb[0] -= dk * dbl_linb[1];
			if (new_e > e) {
				translate_gmp(c, c->gmp_b, c->gmp_linb, 
						c->gmp_mb, -dk);
				dbl_linb[0] += dk * dbl_linb[1];
				e = new_e;
				dk = (int) (1 + 1.1 * (double) dk);
			}
			else {
				dk = (int) (0.9 * (double) dk - 1);
				if (dk < 1)
					dk = 1;
			}
		}
		for (i = 0; i < 6; i++)
			dbl_b[i] = mpz_get_d(c->gmp_b[i]);
		niter++;
	}

	profile_start(PROF_MURPHY_E);
	murphy_e_score(&e, 5, dbl_b0, 1, dbl_linb, 
			alpha, 0., s, 10000, assess);
	profile_stop(PROF_MURPHY_E);
	if (e >= data->min_e) {
		if (data->p_bound < 2000) {
			if (verbose > 2)
				printf("(%f,%g) -> ", alpha, e);

			profile_start(PROF_MURPHY_ROOTS);
			murphy_alpha_exact(&alpha, 5, assess, c->gmp_b, 2000);
			profile_stop(PROF_MURPHY_ROOTS);

			profile_start(PROF_MURPHY_E);
			murphy_e_score(&e, 5, dbl_b0, 1, dbl_linb, 
					alpha, 0., s, 10000, assess);
			profile_stop(PROF_MURPHY_E);
			if (verbose > 2)
				printf("(%f,%g)\n", alpha, e);
		}
	}

	*new_skewness = s;
	*alphaptr = alpha;
	*eptr = e;
	*norm_ptr = sqrt(ifs(dbl_b, s));

	if (verbose > 2) {
		int i;

		printf("optimize 3\nskewness: %.2f norm: %.4e "
			"alpha %.2f Murphy_E: %.3e\npol0: ", 
			*new_skewness, *norm_ptr, *alphaptr, *eptr);
		for (i = 5; i >= 0; i--) {
			mpz_out_str(stdout, 10, c->gmp_b[i]);
			printf(" ");
		}
		printf("\npol1: ");
		mpz_out_str(stdout, 10, c->gmp_linb[1]);
		printf(" ");
		mpz_out_str(stdout, 10, c->gmp_linb[0]);
		printf("\n\n");
	}
}

