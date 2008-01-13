/*----------------------------------------------------------------------
Copyright 2007, Jason Papadopoulos

This file is part of GGNFS.

GGNFS is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GGNFS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GGNFS; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
----------------------------------------------------------------------*/

#include "sieve.h"

/*------------------------------------------------------------------*/
void print_relation(msieve_obj *obj, int64 a, uint32 b, 
			uint32 *factors_r, uint32 num_factors_r, 
			uint32 *large_prime_r,
			uint32 *factors_a, uint32 num_factors_a, 
			uint32 *large_prime_a) {
	
	uint32 i, j;
	char buf[LINE_BUF_SIZE];
	char *tmp = buf;

	tmp += sprintf(buf, "%" PRId64 ",%u", a, b);
	for (i = 0; i < num_factors_r; i++) {
		if (i == 0)
			tmp += sprintf(tmp, ":%x", factors_r[i]);
		else
			tmp += sprintf(tmp, ",%x", factors_r[i]);
	}
	for (j = 0; j < MAX_LARGE_PRIMES; j++) {
		if (large_prime_r[j] == 1)
			continue;
		if (i == 0)
			tmp += sprintf(tmp, ":%x", large_prime_r[j]);
		else
			tmp += sprintf(tmp, ",%x", large_prime_r[j]);
		i++;
	}

	for (i = 0; i < num_factors_a; i++) {
		if (i == 0)
			tmp += sprintf(tmp, ":%x", factors_a[i]);
		else
			tmp += sprintf(tmp, ",%x", factors_a[i]);
	}
	for (j = 0; j < MAX_LARGE_PRIMES; j++) {
		if (large_prime_a[j] == 1)
			continue;
		if (i == 0)
			tmp += sprintf(tmp, ":%x", large_prime_a[j]);
		else
			tmp += sprintf(tmp, ",%x", large_prime_a[j]);
		i++;
	}
	sprintf(tmp, "\n");
	savefile_write_line(&obj->savefile, buf);
}

/*------------------------------------------------------------------*/
uint32 fplog(uint32 k, double log_of_base) {

	/* express k in a different base */

	return (uint32)(log((double)k)/log_of_base + 0.5);
}

/*------------------------------------------------------------------*/
int32 fplog_eval_poly(int64 a, uint32 b, 
			mp_poly_t *f, double log_base,
			uint32 *bits) { 

	/* Compute f(a,b) and take its logarithm in the
	   current base of sieve logarithms, rounding down */

	signed_mp_t norm;
	eval_poly(&norm, a, b, f);
	*bits = mp_bits(&norm.num);
	return (uint32)((*bits) * M_LN2 / log_base);
}

/*------------------------------------------------------------------*/
/* Bases will be chosen for logarithms with the goal of
   hitting this (base-2) target value at most */

#define LOG_TARGET 220

double get_log_base(mp_poly_t *poly, 
			int64 a0, int64 a1, uint32 b) { 

	/* Decide on a base for the logs of one polynomial. 

	   The rational poly is assumed linear, its maximum value 
	   occurs at one of the endpoints of the sieve interval

	   Rigorously finding the extreme values of the algebraic poly
	   would require finding the minima and maxima, and comparing
	   the polynomial values there to the values at a0 and a1. However, 
	   experiments show that for small b the values at the endpoints
	   are much larger than those at the extreme values, and for large
	   b the values are close but the extreme values tend to be outside
	   the sieve interval. Hence we cheat and just do the same as with 
	   the rational poly */

	double t;
	signed_mp_t tmp1, tmp2;

	eval_poly(&tmp1, a0, b, poly);
	eval_poly(&tmp2, a1, b, poly);

	if (mp_cmp(&tmp2.num, &tmp1.num) > 0) 
		t = mp_bits(&tmp1.num);
	else
		t = mp_bits(&tmp2.num);

	/* the base to use is a number x such that 
	   log_x (2^t) = LOG_TARGET. After much re-
	   arranging, x amounts to: */

	return pow(2.0, t / LOG_TARGET);
}

/*------------------------------------------------------------------*/
uint32 read_last_line(msieve_obj *obj, mp_t *n) {

	uint32 last_line = 0;
	char buf[LINE_BUF_SIZE];
	FILE *linefile;
	mp_t read_n;

	sprintf(buf, "%s.line", obj->savefile.name);
	linefile = fopen(buf, "r");
	if (linefile == NULL)
		return last_line;

	fgets(buf, (int)sizeof(buf), linefile);
	mp_clear(&read_n);
	if (buf[0] == 'N')
		mp_str2mp(buf + 2, &read_n, 10);
	if (mp_cmp(n, &read_n) == 0) {
		fgets(buf, (int)sizeof(buf), linefile);
		last_line = atoi(buf);
	}

	fclose(linefile);
	return last_line;
}

/*------------------------------------------------------------------*/
void write_last_line(msieve_obj *obj, mp_t *n, uint32 b) {

	char buf[LINE_BUF_SIZE];
	FILE *linefile;

	sprintf(buf, "%s.line", obj->savefile.name);
	linefile = fopen(buf, "w");
	if (linefile == NULL) {
		printf("error: cannot open linefile '%s'\n", buf);
		exit(-1);
	}

	fprintf(linefile, "N %s\n", mp_sprintf(n, 10, obj->mp_sprintf_buf));
	fprintf(linefile, "%u\n", b);
	fclose(linefile);
}

/*------------------------------------------------------------------*/
uint32 add_free_relations(msieve_obj *obj, 
			sieve_param_t *params, mp_t *n) {

	uint32 alg_degree;
	uint32 rat_degree;
	factor_base_t fb;
	char buf[256];
	uint32 num_relations = 0;
	prime_sieve_t prime_sieve;

	/* read in the NFS polynomials */

	memset(&fb, 0, sizeof(fb));
	if (read_poly(obj, n, &fb.rfb.poly, &fb.afb.poly)) {
		printf("error: failed to read NFS polynomials\n");
		exit(-1);
	}

	alg_degree = fb.afb.poly.degree;
	rat_degree = fb.rfb.poly.degree;
	init_prime_sieve(&prime_sieve, 100, params->afb_limit);

	while (1) {
		uint32 dummy_roots[MAX_POLY_DEGREE];
		uint32 high_coeff;
		uint32 num_roots;
		uint32 p = get_next_prime(&prime_sieve);

		if (p > params->afb_limit)
			break;

		num_roots = poly_get_zeros(dummy_roots, 
					&fb.afb.poly,
					p, &high_coeff, 1);

		if (num_roots == alg_degree && high_coeff != 0) {

			num_roots = poly_get_zeros(dummy_roots, 
						&fb.rfb.poly,
						p, &high_coeff, 1);

			if (num_roots == rat_degree && high_coeff != 0) {
				sprintf(buf, "%u,0:\n", p);
				savefile_write_line(&obj->savefile, buf);
				num_relations++;
			}
		}
	}

	logprintf(obj, "added %u free relations\n", num_relations);
	free_prime_sieve(&prime_sieve);
	savefile_flush(&obj->savefile);
	return num_relations;
}
