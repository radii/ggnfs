#include "stage2_impl.h"

/*-------------------------------------------------------------------------*/
static int
pol_expand(curr_poly_t *c, mpz_t gmp_N)
{
/* compute coefficients */
	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help4);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_a[5]);
	mpz_sub(c->gmp_help3, gmp_N, c->gmp_help4);
	mpz_fdiv_qr(c->gmp_help3, c->gmp_help1, c->gmp_help3, c->gmp_p);
	if (mpz_sgn(c->gmp_help1))
		return 0;

	if (mpz_cmp_ui(c->gmp_p, 1) == 0)
		mpz_set_ui(c->gmp_help2, 1);
	else {
		if (!mpz_invert(c->gmp_help2, c->gmp_d, c->gmp_p))
			return 0;
	}
	mpz_mul(c->gmp_help4, c->gmp_help2, c->gmp_help2);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help4);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help3);
	mpz_fdiv_r(c->gmp_a[4], c->gmp_help4, c->gmp_p);
	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help4);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_a[4]);
	mpz_sub(c->gmp_help3, c->gmp_help3, c->gmp_help4);
	mpz_fdiv_qr(c->gmp_help3, c->gmp_help1, c->gmp_help3, c->gmp_p);
	if (mpz_sgn(c->gmp_help1))
		return 0;

	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_d);
	mpz_fdiv_q(c->gmp_a[3], c->gmp_help3, c->gmp_help4);
	mpz_fdiv_q(c->gmp_a[3], c->gmp_a[3], c->gmp_p);
	mpz_mul(c->gmp_a[3], c->gmp_a[3], c->gmp_p);
	mpz_mul(c->gmp_help4, c->gmp_help2, c->gmp_help2);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help2);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help3);
	mpz_fdiv_r(c->gmp_help1, c->gmp_help4, c->gmp_p);
	mpz_add(c->gmp_a[3], c->gmp_a[3], c->gmp_help1);
	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_a[3]);
	mpz_sub(c->gmp_help3, c->gmp_help3, c->gmp_help4);
	mpz_fdiv_qr(c->gmp_help3, c->gmp_help1, c->gmp_help3, c->gmp_p);
	if (mpz_sgn(c->gmp_help1))
		return 0;

	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_fdiv_q(c->gmp_a[2], c->gmp_help3, c->gmp_help4);
	mpz_fdiv_q(c->gmp_a[2], c->gmp_a[2], c->gmp_p);
	mpz_mul(c->gmp_a[2], c->gmp_a[2], c->gmp_p);
	mpz_mul(c->gmp_help4, c->gmp_help2, c->gmp_help2);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_help3);
	mpz_fdiv_r(c->gmp_help1, c->gmp_help4, c->gmp_p);
	mpz_add(c->gmp_a[2], c->gmp_a[2], c->gmp_help1);
	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_d);
	mpz_mul(c->gmp_help4, c->gmp_help4, c->gmp_a[2]);
	mpz_sub(c->gmp_help3, c->gmp_help3, c->gmp_help4);
	mpz_fdiv_qr(c->gmp_help3, c->gmp_help1, c->gmp_help3, c->gmp_p);
	if (mpz_sgn(c->gmp_help1))
		return 0;

	mpz_fdiv_q(c->gmp_a[1], c->gmp_help3, c->gmp_d);
	mpz_fdiv_q(c->gmp_a[1], c->gmp_a[1], c->gmp_p);
	mpz_mul(c->gmp_a[1], c->gmp_a[1], c->gmp_p);
	mpz_mul(c->gmp_help4, c->gmp_help3, c->gmp_help2);
	mpz_fdiv_r(c->gmp_help1, c->gmp_help4, c->gmp_p);
	mpz_add(c->gmp_a[1], c->gmp_a[1], c->gmp_help1);
	mpz_mul(c->gmp_help4, c->gmp_d, c->gmp_a[1]);
	mpz_sub(c->gmp_help3, c->gmp_help3, c->gmp_help4);
	mpz_fdiv_qr(c->gmp_help3, c->gmp_help1, c->gmp_help3, c->gmp_p);
	if (mpz_sgn(c->gmp_help1))
		return 0;

	mpz_set(c->gmp_a[0], c->gmp_help3);

	mpz_fdiv_qr(c->gmp_help1, c->gmp_a[3], c->gmp_a[3], c->gmp_d);
	mpz_add(c->gmp_help2, c->gmp_a[3], c->gmp_a[3]);
	if (mpz_cmp(c->gmp_d, c->gmp_help2) < 0) {
		mpz_sub(c->gmp_a[3], c->gmp_a[3], c->gmp_d);
		mpz_add_ui(c->gmp_help1, c->gmp_help1, 1);
	}
	mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_p);
	mpz_add(c->gmp_a[4], c->gmp_a[4], c->gmp_help1);

	mpz_fdiv_qr(c->gmp_help1, c->gmp_a[2], c->gmp_a[2], c->gmp_d);
	mpz_add(c->gmp_help2, c->gmp_a[2], c->gmp_a[2]);
	if (mpz_cmp(c->gmp_d, c->gmp_help2) < 0) {
		mpz_sub(c->gmp_a[2], c->gmp_a[2], c->gmp_d);
		mpz_add_ui(c->gmp_help1, c->gmp_help1, 1);
	}
	mpz_mul(c->gmp_help1, c->gmp_help1, c->gmp_p);
	mpz_add(c->gmp_a[3], c->gmp_a[3], c->gmp_help1);

	mpz_set(c->gmp_lina[1], c->gmp_p);
	mpz_neg(c->gmp_lina[0], c->gmp_d);
	mpz_invert(c->gmp_m, c->gmp_p, gmp_N);
	mpz_mul(c->gmp_m, c->gmp_m, c->gmp_d);
	mpz_mod(c->gmp_m, c->gmp_m, gmp_N);

	if (verbose > 2) {
		int i;

		printf("pol-expand\npol0: ");
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

	return 1;
}

/*-------------------------------------------------------------------------*/
void
curr_poly_init(curr_poly_t *c)
{
	int i;

	mpz_init(c->gmp_m);
	mpz_init(c->gmp_mb);
	mpz_init(c->gmp_p);
	mpz_init(c->gmp_d);
	mpz_init(c->gmp_lina[0]);
	mpz_init(c->gmp_lina[1]);
	mpz_init(c->gmp_linb[0]);
	mpz_init(c->gmp_linb[1]);
	for (i = 0; i < 6; i++)
		mpz_init(c->gmp_a[i]);
	for (i = 0; i < 6; i++)
		mpz_init(c->gmp_b[i]);
	mpz_init(c->gmp_help1);
	mpz_init(c->gmp_help2);
	mpz_init(c->gmp_help3);
	mpz_init(c->gmp_help4);
}

/*-------------------------------------------------------------------------*/
void
curr_poly_free(curr_poly_t *c)
{
	int i;

	mpz_clear(c->gmp_m);
	mpz_clear(c->gmp_mb);
	mpz_clear(c->gmp_p);
	mpz_clear(c->gmp_d);
	mpz_clear(c->gmp_lina[0]);
	mpz_clear(c->gmp_lina[1]);
	mpz_clear(c->gmp_linb[0]);
	mpz_clear(c->gmp_linb[1]);
	for (i = 0; i < 6; i++)
		mpz_clear(c->gmp_a[i]);
	for (i = 0; i < 6; i++)
		mpz_clear(c->gmp_b[i]);
	mpz_clear(c->gmp_help1);
	mpz_clear(c->gmp_help2);
	mpz_clear(c->gmp_help3);
	mpz_clear(c->gmp_help4);
}

/*-------------------------------------------------------------------------*/
static int
read_a5pd(curr_poly_t *c, poly_stage2_t *data)
{
	char input_line[256];
	char *line, *end;

	input_line[0] = 0;
	fgets(input_line, sizeof(input_line), data->infile);
	line = input_line;
	end = strchr(line, ' ');
	if (end == NULL)
		return -1;
	*end = 0;
	if (mpz_set_str(c->gmp_a[5], line, 10))
		return -1;
	line = end + 1;
	end = strchr(line, ' ');
	if (end == NULL)
		return -1;
	*end = 0;
	if (mpz_set_str(c->gmp_p, line, 10))
		return -1;
	line = end + 1;
	if (mpz_set_str(c->gmp_d, line, 10))
		return -1;
	return 1;
}

/*-------------------------------------------------------------------------*/
static void
optimize(curr_poly_t *c, poly_stage2_t *data, 
		root_sieve_t *rs, stage2_stat_t *stats)
{
	int err;
	double skewness;
	double pol_norm;
	double alpha_proj;
	double log_max_norm_2 = log(data->max_norm_2);

	profile_start(PROF_ALL);

	while (1) {
		err = read_a5pd(c, data);
		if (err < 0) {
			if (feof(data->infile))
				break;
			continue;
		}
		if (!pol_expand(c, data->gmp_N)) {
			mpz_out_str(stdout, 10, c->gmp_a[5]);
			fprintf(stderr, "expand failed\n");
			continue;
		}
		profile_start(PROF_INITIAL_OPTIMIZE);
		optimize_1(c, &skewness, &pol_norm, &alpha_proj);
		profile_stop(PROF_INITIAL_OPTIMIZE);

		if (pol_norm * exp(alpha_proj) > data->max_norm_1)
			continue;
		root_sieve_run(c, log_max_norm_2, data, rs, stats,
				skewness, pol_norm, alpha_proj);
	}

	profile_stop(PROF_ALL);
}

/*-------------------------------------------------------------------------*/
void
stage2_stat_init(stage2_stat_t *stats)
{
	memset(stats, 0, sizeof(stage2_stat_t));
	profile_init(&stats->profile, PROF_MAX);
}

/*------------------------------------------------------------------------*/
void
stage2_stat_free(stage2_stat_t *stats)
{
	profile_free(&stats->profile);
	memset(stats, 0, sizeof(stage2_stat_t));
}

/*-------------------------------------------------------------------------*/
void
poly_stage2_init(poly_stage2_t *data)
{
	memset(data, 0, sizeof(poly_stage2_t));
	mpz_init(data->gmp_N);
	data->max_norm_1 = 1e20;
	data->max_norm_2 = 1e18;
	data->min_e = 0;
	data->p_bound = 2000;
	data->bound0 = 1e7;
	data->bound1 = 5e6;
	data->area = 1e16;
}

/*-------------------------------------------------------------------------*/
void
poly_stage2_free(poly_stage2_t *data)
{
	mpz_clear(data->gmp_N);
}

/*-------------------------------------------------------------------------*/
s32
poly_stage2_run(poly_stage2_t *data)
{
	stage2_stat_t stats;
	curr_poly_t curr_poly;
	root_sieve_t root_sieve;

	stage2_stat_init(&stats);
	curr_poly_init(&curr_poly);
	root_sieve_init(&root_sieve);
	init_assess(data->bound0, data->bound1, data->area, data->p_bound);

	optimize(&curr_poly, data, &root_sieve, &stats);

#ifdef DO_PROFILE
	profile_done(&stats.profile);
	printf("\nTiming:");
	printf("\nTotal            ");
	profile_print(&stats.profile, PROF_ALL);
	printf("\n  1. Optimization  ");
	profile_print(&stats.profile, PROF_INITIAL_OPTIMIZE);
	printf("\n  Sieve all        ");
	profile_print(&stats.profile, PROF_SIEVE_ALL);
	printf("\n    Init sieve       ");
	profile_print(&stats.profile, PROF_INIT_SIEVE);
	printf("\n    Sieve            ");
	profile_print(&stats.profile, PROF_SIEVE);
	printf("\n    eval             ");
	profile_print(&stats.profile, PROF_EVAL);
	printf("\n      gmp/alpha1        ");
	profile_print(&stats.profile, PROF_ALPHA1);
	printf("\n      gmp/alpha2        ");
	profile_print(&stats.profile, PROF_ALPHA2);
	printf("\n      2. Optimization  ");
	profile_print(&stats.profile, PROF_OPTIMIZE2);
	printf("\n      3. Optimization  ");
	profile_print(&stats.profile, PROF_OPTIMIZE3);
	printf("\n        polroots         ");
	profile_print(&stats.profile, PROF_ROOTS);
	printf("\n        murphy-e sum     ");
	profile_print(&stats.profile, PROF_MURPHY);
	printf("\n");
#endif

	stage2_stat_free(&stats);
	curr_poly_free(&curr_poly);
	root_sieve_free(&root_sieve);
	return 1;
}
