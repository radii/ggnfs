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

#include "sqrt.h"

#define LOG2_IDEAL_HASHTABLE_SIZE 21

/* we will need to find primes q for which f(x) mod q
   is irreducible. This is the maximum number of q to try */

#define NUM_PRIME_RETRIES 100

/*--------------------------------------------------------------------*/
uint32 get_prime_for_sqrt(mp_poly_t *alg_poly,
			  uint32 min_value,
			  uint32 *q_out) {

	uint32 i;
	uint32 status = 0;
	uint32 q = 0;
	prime_sieve_t prime_sieve;

	init_prime_sieve(&prime_sieve, min_value, (uint32)(-1));
	for (i = 0; i < NUM_PRIME_RETRIES; i++) {
		uint32 tmp_q = get_next_prime(&prime_sieve);
		if (is_irreducible(alg_poly, tmp_q)) {
			q = tmp_q;
			break;
		}
		else if (q == 0) {
			/* in rare cases, there is no q for which alg_poly
			   mod q is irreducible. Life becomes much more
			   difficult in this case, since square roots mod
			   q will not be unique. The only alternative when
			   this happens is to pick a q such that alg_poly 
			   mod q has no linear roots (so that all of the
			   relations mod q are relatively prime to alg_poly), 
			   then keep trying dependencies until by luck we 
			   start with the correct initial square root for 
			   the Newton iteration to succeed.

			   Buhler et. al. show that for polynomial degree d,
			   on average one in 2^(d/2) dependencies will lead to
			   a congruence of squares (and about half of those
			   will lead to a factor). Technically we also need to
			   check that alg_poly mod q is squarefree, but that
			   would require a full polynomial factoring routine;
			   I'm gambling that being squarefree is not rare.
			  
			   Note that experience with choosing suboptimal q
			   shows that life isn't too much more difficult than
			   the ordinary case. Experimentally, the square
			   root algorithm fails to converge approximately
			   half the time */

			uint32 roots[MAX_POLY_DEGREE];
			uint32 high_coeff;
			uint32 num_roots = poly_get_zeros(roots, alg_poly,
						tmp_q, &high_coeff, 1);
			if (high_coeff != 0 && num_roots == 0)
				q = tmp_q;
		}
	}
	free_prime_sieve(&prime_sieve);

	if (i == NUM_PRIME_RETRIES)
		status = 1;

	*q_out = q;
	return status;
}

/*--------------------------------------------------------------------*/
static void eval_poly_derivative(mp_poly_t *poly, 
				signed_mp_t *m1, signed_mp_t *m0, 
				mp_t *n, mp_t *res) {
	uint32 i;
	mp_t next_coeff;
	mp_t m1_pow, m1_tmp, m0_tmp;

	mp_copy(&m1->num, &m1_tmp);
	if (m1->sign == NEGATIVE)
		mp_sub(n, &m1_tmp, &m1_tmp);
	mp_copy(&m1_tmp, &m1_pow);

	mp_copy(&m0->num, &m0_tmp);
	if (m0->sign == POSITIVE)
		mp_sub(n, &m0_tmp, &m0_tmp);

	i = poly->degree;
	mp_mul_1(&poly->coeff[i].num, i, res);

	while (--i) {
		signed_mp_t *coeff = poly->coeff + i;

		mp_modmul(res, &m0_tmp, n, res);

		mp_mul_1(&coeff->num, i, &next_coeff);
		if (coeff->sign == NEGATIVE)
			mp_sub(n, &next_coeff, &next_coeff);
		mp_modmul(&next_coeff, &m1_pow, n, &next_coeff);
		mp_add(res, &next_coeff, res);

		if (i > 1)
			mp_modmul(&m1_pow, &m1_tmp, n, &m1_pow);
	}

	if (mp_cmp(res, n) >= 0)
		mp_sub(res, n, res);
}

/*--------------------------------------------------------------------*/
static uint32 rat_square_root(relation_t *rlist, uint32 num_relations,
				mp_t *n, mp_t *sqrt_r) {
	uint32 i, j;
	hashtable_t h;
	uint32 already_seen;
	mp_t base, exponent, tmp;
	uint32 status = 0;

	/* count up the number of times each prime factor in
	   rlist occurs */

	hashtable_init(&h, LOG2_IDEAL_HASHTABLE_SIZE, 10000, 1);
	for (i = 0; i < num_relations; i++) {
		relation_t *r = rlist + i;

		if (r->refcnt % 2 == 0)
			continue;

		for (j = 0; j < r->num_factors_r; j++) {
			hash_t *curr_ideal = hashtable_find(&h,
						r->factors + j, 
						&already_seen);
			if (!already_seen)
				curr_ideal->payload[1] = 1;
			else
				curr_ideal->payload[1]++;
		}
	}

	/* verify all such counts are even, and form the 
	   rational square root */

	mp_clear(&base); base.nwords = 1;
	mp_clear(&exponent); exponent.nwords = 1;
	mp_clear(sqrt_r); sqrt_r->nwords = sqrt_r->val[0] = 1;

	for (i = 1; i < h.match_array_size; i++) {
		hash_t *curr_ideal = h.match_array + i;
		if (curr_ideal->payload[1] % 2) {
			status = 1;
			break;
		}
		if (curr_ideal->payload[0]) {
			base.val[0] = curr_ideal->payload[0];
			exponent.val[0] = curr_ideal->payload[1] / 2;
			mp_expo(&base, &exponent, n, &tmp);
			mp_modmul(sqrt_r, &tmp, n, sqrt_r);
		}
	}

	hashtable_free(&h);
	return status;
}

/*--------------------------------------------------------------------*/
static uint32 verify_alg_ideal_powers(relation_t *rlist, 
					uint32 num_relations,
					uint32 *num_free_relations) {
	uint32 i, j, k;
	hashtable_t h;
	uint32 already_seen;
	uint32 *counts = NULL;
	uint32 count_alloc = 0;
	uint32 num_odd;

	/* the return value is the total number of relations
	   counted, which will be smaller than num_relations
	   since only odd multiplicity relations appear, and
	   even then only once */

	/* count the multiplicity of each algebraic ideal (not
	   just the prime to which the ideal corresponds) in
	   rlist. This is unnecessarily messy because we have to
	   store the counts in a separate array */

	*num_free_relations = 0;

	hashtable_init(&h, LOG2_IDEAL_HASHTABLE_SIZE, 10000, 
			(uint32)(sizeof(ideal_t) / sizeof(uint32)));

	for (i = num_odd = 0; i < num_relations; i++) {
		relation_t *r = rlist + i;
		relation_lp_t rlp;

		if (r->refcnt % 2 == 0)
			continue;

		find_large_ideals(r, &rlp, 0, 0);

		for (j = 0; j < rlp.ideal_count; j++) {
			ideal_t *curr_ideal = rlp.ideal_list + j;
			hash_t *curr_ideal_hash;

			if (curr_ideal->i.rat_or_alg == RATIONAL_IDEAL)
				continue;

			curr_ideal_hash = hashtable_find(&h,
					curr_ideal, &already_seen);
			k = curr_ideal_hash - h.match_array;

			if (!already_seen && 
			    count_alloc < h.match_array_alloc) {
				counts = (uint32 *)xrealloc(counts,
						h.match_array_alloc *
						sizeof(uint32));
				memset(counts + count_alloc, 0,
						(h.match_array_alloc -
				       		count_alloc) * sizeof(uint32));
				count_alloc = h.match_array_alloc;
			}
			counts[k]++;
		}

		num_odd++;
		if (r->b == 0)
			(*num_free_relations)++;
	}

	/* verify each ideal occurs an even number of times */

	for (i = 1; i < h.match_array_size; i++) {
		if (counts[i] % 2) {
			num_odd = 0;
			break;
		}
	}

	hashtable_free(&h);
	free(counts);
	return num_odd;
}

/*--------------------------------------------------------------------*/
uint32 nfs_find_factors(msieve_obj *obj, mp_t *n, 
			factor_list_t *factor_list) {

	/* external interface for the NFS square root */

	uint32 i, j;
	uint32 check_q;
	factor_base_t fb;
	mp_poly_t monic_alg_poly;
	mp_t sqrt_r, sqrt_a;
	mp_t c, tmp1, tmp2;
	signed_mp_t *m0;
	signed_mp_t *m1;
	uint32 dep_lower = 1;
	uint32 dep_upper = 64;
	uint32 factor_found = 0;

	logprintf(obj, "\n");
	logprintf(obj, "commencing square root phase\n");

	/* read in the factor base */

	memset(&fb, 0, sizeof(fb));
	if (read_poly(obj, n, &fb.rfb.poly, &fb.afb.poly)) {
		logprintf(obj, "polynomials not found\n");
		return 0;
	}

	/* find the values needed to convert the algebraic 
	   square root back to an integer */

	if (fb.rfb.poly.degree != 1) {
		logprintf(obj, "cannot handle non-linear polynomials\n");
		return 0;
	}
	m0 = &fb.rfb.poly.coeff[0];
	m1 = &fb.rfb.poly.coeff[1];

	/* construct a monic version of the algebraic poly,
	   saving off the leading coefficient separately */

	j = fb.afb.poly.degree;
	if (fb.afb.poly.coeff[j].sign == NEGATIVE) {
		logprintf(obj, "cannot handle negative leading "
				"algebraic polynomial coefficient\n");
		return 0;
	}
	memset(&monic_alg_poly, 0, sizeof(mp_poly_t));
	mp_copy(&fb.afb.poly.coeff[j].num, &c);
	mp_copy(&c, &tmp1);
	monic_alg_poly.coeff[j-1] = fb.afb.poly.coeff[j-1];
	for (i = j - 2; (int32)i >= 0; i--) {
		mp_mul(&fb.afb.poly.coeff[i].num, &tmp1,
				&monic_alg_poly.coeff[i].num);
		monic_alg_poly.coeff[i].sign = fb.afb.poly.coeff[i].sign;

		if (i > 0) {
			mp_mul(&c, &tmp1, &tmp2);
			mp_copy(&tmp2, &tmp1);
		}
	}
	monic_alg_poly.coeff[j].num.nwords = 1;
	monic_alg_poly.coeff[j].num.val[0] = 1;
	monic_alg_poly.degree = fb.afb.poly.degree;
	get_prime_for_sqrt(&monic_alg_poly, (uint32)0x80000000, &check_q);

	/* determine the list of dependencies to compute */

	if (obj->nfs_lower && obj->nfs_upper) {
		dep_lower = MIN(obj->nfs_lower, 64);
		dep_upper = MIN(obj->nfs_upper, 64);
		dep_upper = MAX(dep_lower, dep_upper);
	}

	/* for each dependency */

	for (i = dep_lower; i <= dep_upper; i++) {
		uint32 k;
		uint32 num_relations;
		uint32 full_num_relations;
		uint32 num_free_relations;
		relation_t *rlist;
		abpair_t *abpairs;
		uint32 num_cycles;
		la_col_t *cycle_list;
		mp_t exponent;

		logprintf(obj, "reading relations for dependency %u\n", i);

		/* read in only the relations for dependency i */

		nfs_read_cycles(obj, &fb, &num_cycles, &cycle_list, 
				&num_relations, &rlist, 0, i);

		if (num_cycles == 0)
			continue;

		/* we don't care about the cycles that the dependency
		   needs, only the raw count of each relation that occurs */

		free_cycle_list(cycle_list, num_cycles);

		/* do some sanity checking, performing increasing
		   amounts of work as the dependency proves itself
		   to be valid */

		full_num_relations = verify_alg_ideal_powers(
					rlist, num_relations,
					&num_free_relations);
		if (full_num_relations == 0) {
			logprintf(obj, "algebraic side is not a square!\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (full_num_relations % 2) {
			/* the number of relations in the dependency must
			   be even, because each relation represents a
			   degree-1 polynomial, and the product of these
			   relations will not have a square root unless the
			   degree of the product is even */
			logprintf(obj, "number of relations is not even\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}
		if (num_free_relations % 2) {
			logprintf(obj, "number of free relations (%u) is "
					"not even\n", num_free_relations);
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}

		if (rat_square_root(rlist, num_relations, n, &sqrt_r) != 0) {
			logprintf(obj, "rational side is not a square!\n");
			nfs_free_relation_list(rlist, num_relations);
			continue;
		}

		/* flatten the list of relations; each occurrence of
		   a relation gets its own abpair_t */

		abpairs = (abpair_t *)xmalloc(full_num_relations *
						sizeof(abpair_t));
		for (j = k = 0; j < num_relations; j++) {
			relation_t *r = rlist + j;
			if (r->refcnt % 2) {
				abpairs[k].a = r->a;
				abpairs[k].b = r->b;
				k++;
			}
		}
		nfs_free_relation_list(rlist, num_relations);

		/* perform the major work: the algebraic square root.
		   Note that to conserve memory, abpairs is freed in
		   the following call */

		mp_clear(&sqrt_a);
		alg_square_root(obj, &monic_alg_poly, n, &c, m1, m0, abpairs, 
					full_num_relations, check_q, &sqrt_a);
		if (mp_is_zero(&sqrt_a)) {
			logprintf(obj, "algebraic square root failed\n");
			continue;
		}

		/* an algebraic square root is available; move on
		   to the final congruence of squares. The arithmetic
		   is as given in Buhler et. al. with one exception:
		   when the rational poly is nonmonic there is a 
		   correction to the final square root value but the 
		   free relations *do not* figure into it. This latter
		   point is completely ignored in the literature! */

		eval_poly_derivative(&fb.afb.poly, m1, m0, n, &tmp1);
		mp_modmul(&tmp1, &sqrt_r, n, &sqrt_r);

		mp_clear(&exponent);
		exponent.nwords = 1;
		if (!mp_is_one(&c)) {
			exponent.val[0] = full_num_relations/ 2 + 
						fb.afb.poly.degree - 2;
			mp_expo(&c, &exponent, n, &tmp2);
			mp_modmul(&tmp2, &sqrt_r, n, &sqrt_r);
		}

		if (!mp_is_one(&m1->num)) {
			exponent.val[0] = (full_num_relations -
				       		num_free_relations) / 2;
			mp_copy(&m1->num, &tmp1);
			if (m1->sign == NEGATIVE)
				mp_sub(n, &tmp1, &tmp1);
			mp_expo(&tmp1, &exponent, n, &tmp2);
			mp_modmul(&tmp2, &sqrt_a, n, &sqrt_a);
		}

		/* a final sanity check: square the rational and algebraic 
		   square roots, expecting the same value modulo n */

		mp_modmul(&sqrt_r, &sqrt_r, n, &tmp1);
		mp_modmul(&sqrt_a, &sqrt_a, n, &tmp2);
		if (mp_cmp(&tmp1, &tmp2) != 0) {
			logprintf(obj, "dependency does not form a "
					"congruence of squares!\n");
			continue;
		}

		/* look for a nontrivial factor of n */

		mp_add(&sqrt_r, &sqrt_a, &tmp1);
		mp_gcd(&tmp1, n, &tmp1);
		if (!mp_is_one(&tmp1) && mp_cmp(n, &tmp1) != 0) {
			factor_found = 1;
			if (factor_list_add(obj, factor_list, &tmp1) == 0)
				break;
		}
	}

	return factor_found;
}
