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

#include "filter.h"

/*--------------------------------------------------------------------*/
static void free_relsets(merge_t *merge) {

	uint32 i;
	relation_set_t *relset_array = merge->relset_array;
	uint32 num_relsets = merge->num_relsets;

	for (i = 0; i < num_relsets; i++) {
		relation_set_t *r = relset_array + i;
		free(r->data);
	}
	free(merge->relset_array);
	merge->relset_array = NULL;
}

/*--------------------------------------------------------------------*/
static void dump_cycles(msieve_obj *obj, merge_t *merge) {

	uint32 i;
	relation_set_t *relset_array = merge->relset_array;
	uint32 num_relsets = merge->num_relsets;
	char buf[256];
	FILE *cycle_fp;

	sprintf(buf, "%s.cyc", obj->savefile.name);
	cycle_fp = fopen(buf, "wb");
	if (cycle_fp == NULL) {
		logprintf(obj, "error: can't open cycle file\n");
		exit(-1);
	}

	fwrite(&num_relsets, sizeof(uint32), (size_t)1, cycle_fp);

	for (i = 0; i < num_relsets; i++) {
		relation_set_t *r = relset_array + i;
		uint32 num = r->num_relations;

		fwrite(&num, sizeof(uint32), (size_t)1, cycle_fp);
		fwrite(r->data, sizeof(uint32),
				(size_t)num, cycle_fp);
	}
	fclose(cycle_fp);
}

/*--------------------------------------------------------------------*/
static void find_fb_size(factor_base_t *fb, 
			uint32 limit_r, uint32 limit_a,
			uint32 *entries_r_out, uint32 *entries_a_out) {

	prime_sieve_t prime_sieve;
	uint32 entries_r = 0;
	uint32 entries_a = 0;

	init_prime_sieve(&prime_sieve, 0, 
			MAX(limit_r, limit_a) + 1000);

	while (1) {
		uint32 p = get_next_prime(&prime_sieve);
		uint32 num_roots;
		uint32 high_coeff;
		uint32 roots[MAX_POLY_DEGREE + 1];

		if (p >= limit_r && p >= limit_a)
			break;

		if (p < limit_r) {
			num_roots = poly_get_zeros(roots, &fb->rfb.poly, p, 
							&high_coeff, 1);
			if (high_coeff == 0)
				num_roots++;
			entries_r += num_roots;
		}

		if (p < limit_a) {
			num_roots = poly_get_zeros(roots, &fb->afb.poly, p, 
							&high_coeff, 1);
			if (high_coeff == 0)
				num_roots++;
			entries_a += num_roots;
		}
	}

	free_prime_sieve(&prime_sieve);
	*entries_r_out = entries_r;
	*entries_a_out = entries_a;
}

/*--------------------------------------------------------------------*/
uint32 nfs_filter_relations(msieve_obj *obj, mp_t *n) {

	filter_t filter;
	merge_t merge;
	uint32 filtmin_r;
	uint32 filtmin_a;
	uint32 entries_r;
	uint32 entries_a;
	uint32 extra_needed;
	uint32 max_weight;
	factor_base_t fb;

	logprintf(obj, "\n");
	logprintf(obj, "commencing relation filtering\n");

	memset(&fb, 0, sizeof(fb));
	if (read_poly(obj, n, &fb.rfb.poly, &fb.afb.poly)) {
		printf("filtering failed to read polynomials\n");
		exit(-1);
	}

	/* delete duplicate relations, and determine the cutoff
	   size of large primes used in the rest of the filtering.
	   We do not just use the factor base, since it may be too 
	   small or too large for this particular dataset
	  
	   Also determine the number of ideals that would be ignored
	   given the cutoff from the duplicate phase; this puts a
	   lower bound on the matrix size */

	filtmin_r = filtmin_a = nfs_purge_duplicates(obj, &fb, 0);
	if (obj->nfs_lower)
		filtmin_r = obj->nfs_lower;
	if (obj->nfs_upper)
		filtmin_a = obj->nfs_upper;

	/* calculate the initial excess to look for */

	find_fb_size(&fb, filtmin_r, filtmin_a, &entries_r, &entries_a);

	/* singleton removal happens in at least two phases. 
	   Phase 1 does the singleton removal and clique 
	   processing for all the ideals larger than filtmin, 
	   regardless of ideal weight. In the first phase we 
	   demand a lot of excess relations. 
	   
	   We avoid pruning too many cliques, to leave some room
	   for phase 2+. Phase 2 always runs, even if there is 
	   not enough excess for phase 1; the first pass is just
	   to remove the majority of the excess, in case there
	   are many more relations than necesary */

	memset(&filter, 0, sizeof(filter));
	filter.filtmin_r = filtmin_r;
	filter.filtmin_a = filtmin_a;
	filter.target_excess = (uint32)(1.5 * (entries_r + entries_a));
	logprintf(obj, "ignoring smallest %u rational and %u "
			"algebraic ideals\n", entries_r, entries_a);

	extra_needed = nfs_purge_singletons(obj, &fb, &filter, 0);

	/* throw away the current relation list */

	free(filter.relation_array);
	filter.relation_array = NULL;

	/* perform the other filtering passes; these use a very small
	   filtering bound, and rely on most unneeded relations
	   having been deleted already. The set of relations remaining
	   will be forwarded to the final merge phase */

	filtmin_r = MIN(filtmin_r / 2, 750000);
	filtmin_a = MIN(filtmin_a / 2, 750000);
	find_fb_size(&fb, filtmin_r, filtmin_a, &entries_r, &entries_a);
	filter.filtmin_r = filtmin_r;
	filter.filtmin_a = filtmin_a;

	for (max_weight = 20; max_weight < 40; max_weight += 5) {

		filter.target_excess = entries_r + entries_a;
		extra_needed = nfs_purge_singletons(obj, &fb, 
						&filter, max_weight);

		/* give up if there is not enough excess 
		   to form a matrix */

		if (filter.num_relations < filter.num_ideals ||
		    filter.num_relations - filter.num_ideals <
						FINAL_EXCESS_FRACTION *
						filter.target_excess) {
			free(filter.relation_array);
			filter.relation_array = NULL;
			return extra_needed;
		}

		/* perform the merge phase */

		extra_needed = filter.target_excess;
		nfs_merge_init(obj, &filter);
		nfs_merge_2way(obj, &filter, &merge);
		nfs_merge_full(obj, &merge, extra_needed);

		/* do not accept the collection of generated cycles
		   unless the matrix they form is dense enough. It's
		   unfortunate, but different datasets could look the
		   same going into the merge phase but lead to very
		   different matrix densities. The problem is that if
		   the matrix is not dense enough experience shows
		   that the linear algebra sometimes is unable to find
		   nontrivial dependencies. This is especially common
		   when the bound on large primes is much larger than
		   is sensible. In fact, I haven't seen a matrix
		   with average cycle weight under 50.0 get solved
		   successfully */

		if (merge.avg_cycle_weight > 60.0)
			break;

		logprintf(obj, "matrix not dense enough, retrying\n");
		free_relsets(&merge);
	}

	if (max_weight == 40) {
		printf("error: too many merge attempts\n");
		exit(-1);
	}

	/* optimize the collection of cycles 
	   and save the result */

	nfs_merge_post(obj, &merge);
	dump_cycles(obj, &merge);
	free_relsets(&merge);
	return 0;
}
