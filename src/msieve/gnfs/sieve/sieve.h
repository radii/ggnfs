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

$Id: sieve.h,v 1.1 2007-12-16 03:54:47 jasonp_sf Exp $
----------------------------------------------------------------------*/

#ifndef _SIEVE_H_
#define _SIEVE_H_

#include "gnfs.h"
#include <ap.h>

#ifdef __cplusplus
extern "C" {
#endif

/* factors smaller than the following are not printed */

#define MAX_SKIPPED_FACTOR 256

/* dump one NFS relation to the savefile */

void print_relation(msieve_obj *obj, int64 a, uint32 b, 
		uint32 *factors_r, uint32 num_factors_r, 
		uint32 *large_prime_r,
		uint32 *factors_a, uint32 num_factors_a, 
		uint32 *large_prime_a);
	
/* convert k bits to another base */

uint32 fplog(uint32 k, double log_of_base);

/* find the size of f(a,b) both in bits and in base log_base */

int32 fplog_eval_poly(int64 a, uint32 b, 
			mp_poly_t *f, double log_base,
			uint32 *bits);

/* compute a base of logarithms suitable for the current
   sieve line. The base chosen should be larger when the
   size of sieve values are smaller, in order to use up more
   of the dynamic range in one byte of the sieve array */

double get_log_base(mp_poly_t *poly, 
			int64 a0, int64 a1, uint32 b);

uint32 read_last_line(msieve_obj *obj, mp_t *n);

void write_last_line(msieve_obj *obj, mp_t *n, uint32 b);

/* prototypes for the subsystem that batch-factors the
   portions of relations that contain large primes. The
   current system is designed for cofactors containing up to
   three large primes, though relations with any number of
   large primes will be batched and factored */

#define MAX_LARGE_PRIMES 3

/* simplified representation of one NFS relation. Note that
   any of the factors may be trivial */

typedef struct {
	int64 a;
	uint32 b;
	uint8 num_factors_r;     /* doesn't include large primes */
	uint8 num_factors_a;     /* doesn't include large primes */
	uint8 lp_r_num_words;     /* one number with this many words */
	uint8 lp_a_num_words;     /* one number with this many words */
	uint32 factor_list_word;  /* offset in an array of factors where
				     all of the above above appear, in order */
} cofactor_t;

/* main structure controlling batch factoring. The main goal
   of batch factoring is to compute gcd(each_relation, 
   product_of_many_primes) much faster than running a conventional
   factoring algorithm on each relation. We hope that the gcd
   splits a given relation into two approximately equal-size halves,
   so that each part is trivial (has one factor) or is easy (has two
   factors) to decompose; failing that, we hope that the gcd provides
   enough information to tell us whether going to the trouble of
   fully factoring the relation is unlikely to yield something useful.

   The algorithm depends critically on the fact that gcd(A,B) =
   gcd(A mod B, B), so that when B is massively smaller than A the
   former is easy to compute. We don't actually compute the gcd
   until each relation B has (A mod B) available, and it is the latter
   that Bernstein's algorithm computes */

typedef struct {
	fastmult_info_t fastmult_info; /* for large-integer multiplies */

	ap_t prime_product;  /* product of primes used in the gcd */

	uint32 num_success;       /* number of surviving relations */
	uint32 target_relations;  /* number of relations to batch up */
	uint32 lp_cutoff_r;       /* maximum size of rational factors */
	mp_t lp_cutoff_r2;        /* square of lp_cutoff_r */
	uint32 lp_cutoff_a;       /* maximum size of algebraic factors */
	mp_t lp_cutoff_a2         /* square of lp_cutoff_a */;
	mp_t max_prime2;          /* the square of the largest prime that
				     occurs in prime_product */

	uint32 num_relations;     /* number relations currently batched */
	uint32 num_relations_alloc;     /* relations allocated */
	cofactor_t *relations;    /* relations currently batched */

	uint32 num_factors;       /* words for factors of batched relations */
	uint32 num_factors_alloc; /* space for batched factors */
	uint32 *factors;          /* factors of batched relations */
} relation_batch_t;

/* initialize the relation batch. Batch factoring uses all
   the primes from min_prime to max_prime. We assume min_prime
   is the smallest element possible in unfactored parts of
   relations, i.e. MIN(max factor base prime within rational
   or algebraic factor bases).
   
   In general we do not want max_prime as large as the 
   large prime cutoffs; making it smaller allows the batch 
   factoring to split most of the cofactors in relations that 
   contain large primes, or at least prove most relations to 
   be not worth the trouble to do so manually */

void relation_batch_init(msieve_obj *obj, relation_batch_t *rb,
			uint32 min_prime, uint32 max_prime,
			uint32 lp_cutoff_r, uint32 lp_cutoff_a);

void relation_batch_free(relation_batch_t *rb);

/* add one relation to the batch. Note that the relation may
   contain any number of large primes, as long as the unfactored
   parts of the relation fit into an mp_t. Modifying the code to
   handle more than three large primes within unfactored_[ra]
   just involves modifying the base case of the recursion; maybe 
   that should be made into a callback */

void relation_batch_add(int64 a, uint32 b, 
			uint32 *factors_r, uint32 num_factors_r, 
			mp_t *unfactored_r,
			uint32 *factors_a, uint32 num_factors_a, 
			mp_t *unfactored_a,
			relation_batch_t *rb);
	
/* factor all the batched relations, saving all the ones whose
   largest {rational|algebraic} factors are all less than
   lp_cutoff_[ra] */

uint32 relation_batch_run(msieve_obj *obj, relation_batch_t *rb);

#ifdef __cplusplus
}
#endif

#endif /* _SIEVE_H_ */
