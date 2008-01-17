#include "stage2_impl.h"

/*-------------------------------------------------------------------------*/
static void
write_polynomial_51(FILE * fi, int deg, mpz_t * coeff1, mpz_t * coeff2, mpz_t m,
		    double skewness, double norm, double alpha, double murphy_e)
{
	int i;

	fprintf(fi, "BEGIN POLY ");
	fprintf(fi, "#skewness %.2f norm %.2e alpha %.2f Murphy_E %.2e\n",
		skewness, norm, alpha, murphy_e);
	for (i = 0; i < deg + 1; i++) {
		fprintf(fi, "X%d ", deg - i);
		mpz_out_str(fi, 10, coeff1[deg - i]);
		fprintf(fi, "\n");
	}
	if (mpz_cmp_ui(coeff2[1], 1)) {
		fprintf(fi, "Y1 ");
		mpz_out_str(fi, 10, coeff2[1]);
		fprintf(fi, "\nY0 ");
		mpz_out_str(fi, 10, coeff2[0]);
		fprintf(fi, "\n");
	}
	fprintf(fi, "M ");
	mpz_out_str(fi, 10, m);
	fprintf(fi, "\nEND POLY\n");
}

/*-------------------------------------------------------------------------*/
void
check(int x, int y, curr_poly_t *c, 
	poly_stage2_t *data, stage2_stat_t *stats,
	double skewness)
{
	int i;
	double alpha, norm, murphye, alpha_max;

	profile_start(PROF_ALPHA1);
	if (verbose > 3)
		printf("check at i=%d, j=%d\n", x, y);

	for (i = 0; i < 6; i++)
		mpz_set(c->gmp_b[i], c->gmp_a[i]);
	mpz_set_si(c->gmp_help1, y);
	mpz_mul(c->gmp_help2, c->gmp_help1, c->gmp_p);
	mpz_add(c->gmp_b[2], c->gmp_b[2], c->gmp_help2);
	mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_d);
	mpz_sub(c->gmp_b[1], c->gmp_b[1], c->gmp_help1);
	mpz_set_si(c->gmp_help1, x);
	mpz_mul(c->gmp_help2, c->gmp_help1, c->gmp_p);
	mpz_add(c->gmp_b[1], c->gmp_b[1], c->gmp_help2);
	mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_d);
	mpz_sub(c->gmp_b[0], c->gmp_b[0], c->gmp_help1);
	mpz_set(c->gmp_mb, c->gmp_m);
	for (i = 0; i < 2; i++)
		mpz_set(c->gmp_linb[i], c->gmp_lina[i]);
	profile_stop(PROF_ALPHA1);

	profile_start(PROF_OPTIMIZE2);
	optimize_2(c, skewness, &skewness, &norm);
	profile_stop(PROF_OPTIMIZE2);

	profile_start(PROF_ALPHA2);
	alpha_max = log(data->max_norm_2 / norm);
	if (compute_alpha(&alpha, 5, NULL, c->gmp_b, alpha_max)) {
		profile_stop(PROF_ALPHA2);
		if (verbose > 2)
			printf("failed\n");
		return;
	}
	if (verbose > 2)
		printf("alpha: %.3f\n", alpha);
	profile_stop(PROF_ALPHA2);

	if (norm * exp(alpha) > data->max_norm_2) {
		if (verbose > 2)
			printf("failed\n");
		return;
	}

	if (verbose > 2)
		printf("i: %d, j: %d\n", x, y);
	profile_start(PROF_OPTIMIZE3);
	optimize_3(c, data, skewness, &skewness, &norm, &murphye, &alpha);
	profile_stop(PROF_OPTIMIZE3);
	if (murphye < data->min_e) {
		if (verbose > 1)
			printf("f");
		return;
	}
	if (verbose > 1)
		printf("P");
	write_polynomial_51(data->outfile, 5, c->gmp_b, c->gmp_linb, 
				c->gmp_mb, skewness, norm, alpha, murphye);
}

