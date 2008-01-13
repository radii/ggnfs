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
	uint32 filtmin_r_pass2;
	uint32 filtmin_a_pass2;
	uint32 entries_r;
	uint32 entries_a;
	uint32 entries_r_pass2;
	uint32 entries_a_pass2;
	uint32 extra_needed;
	uint32 num_passes;
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

	/* singleton removal happens in several phases. The first
	   phase does the singleton removal and clique processing
	   for all the ideals larger than filtmin. In the first
	   phase we demand a lot of excess relations. 
	   
	   We avoid pruning too many cliques, to leave some room
	   for succeeding phases */

	memset(&filter, 0, sizeof(filter));
	filter.filtmin_r = filtmin_r;
	filter.filtmin_a = filtmin_a;
	filter.target_excess = (uint32)(1.7 * (entries_r + entries_a));
	logprintf(obj, "ignoring smallest %u rational and %u "
			"algebraic ideals\n", entries_r, entries_a);

	extra_needed = nfs_purge_singletons(obj, &fb, &filter, 1);

	/* give up if there is not enough excess to form a matrix
	   (we require 8% more excess than is strictly needed to
	   construct a matrix) */

	if (filter.num_relations < filter.num_ideals ||
	    filter.num_relations - filter.num_ideals <
	    		1.08 * (entries_r + entries_a)) {
		free(filter.relation_array);
		filter.relation_array = NULL;
		return extra_needed;
	}

	/* We have excess relations to spare. This means the singleton
	   removal can be repeated with a smaller factor base bound,
	   so that more large ideals can be explicitly tracked. It also
	   means we have more room to delete cliques. Unfortunately,
	   it's not clear up front *how much* more room there is; the
	   dataset must have a minimum amount of excess to form a matrix
	   at all, but the only way to size up the dataset in a single
	   pass is to assume an unmanageably low factor base bound,
	   do the singleton removal, then sacrifice ideals until the
	   minimum amount of matrix excess is achieved. That in turn
	   probably requires too much memory.

	   Instead, if the current dataset has sufficient excess, we
	   iteratively lower the factor base bound a little (10% at 
	   a time), then repeat the singleton and clique removal with 
	   a lower (10% at a time) amount of desired excess.

	   We will keep information on ideals larger than filtmin_pass2, 
	   not just ideals larger than filtmin. This technically means 
	   the merge phase can produce an extremely small matrix, but 
	   that matrix would be hopelessly dense.  Instead, the merge 
	   code will use the extra information to skip merging ideals 
	   that appear in many relations. The ignored ideals increase 
	   the number of cycles that must be found, until the process 
	   converges to a moderate size matrix that is still tolerably
	   sparse. This multi-step process is how the pros filter :) */

	num_passes = 1;
	filtmin_r_pass2 = filtmin_r;
	filtmin_a_pass2 = filtmin_a;
	entries_r_pass2 = entries_r;
	entries_a_pass2 = entries_a;

	while (1) {
		double excess_fraction = (double)(filter.num_relations - 
				filter.num_ideals) / (entries_r + entries_a);

		logprintf(obj, "dataset has %3.1f%% excess relations\n",
				100.0 * (excess_fraction - 1.0));

		/* throw away the current relation list */

		free(filter.relation_array);
		filter.relation_array = NULL;

		/* when singleton removal completes the code 
		   remembers the maximum number of relations 
		   an ideal appears in. When this number is large 
		   (20-something or more) the dataset is sufficiently 
		   dense that the heuristics below and in the merging 
		   code work correctly. If it is smaller, then the
		   current value of filtmin is set too conservatively,
		   and we should lower it a little bit */

		if (filter.max_ideal_degree < 22) {
			filtmin_r = filtmin_r_pass2;
			filtmin_a = filtmin_a_pass2;
			entries_r = entries_r_pass2;
			entries_a = entries_a_pass2;
			num_passes = 1;
		}

		/* reduce the bound on large ideals, so that there
		   are more ideals explicitly tracked */

		filtmin_r_pass2 = filtmin_r * (1.0 - 0.1 * num_passes);
		filtmin_a_pass2 = filtmin_a * (1.0 - 0.1 * num_passes);

		/* find the new number of entries needed for a matrix */

		find_fb_size(&fb, filtmin_r_pass2, filtmin_a_pass2,
				&entries_r_pass2, &entries_a_pass2);
		logprintf(obj, "ignoring smallest %u rational and %u "
					"algebraic ideals\n", 
					entries_r_pass2, entries_a_pass2);
		excess_fraction = MAX(1.08, excess_fraction * (1.0 -
					0.1 * num_passes));

		/* repeat the in-memory singleton removal with the more
		   stringent filtering bound. The number of large ideals 
		   will increase, so if there was too much excess before
		   there is less now. Any singletons or cliques removed
		   in previous passes would still need to be removed in this
		   pass, so the temporary files generated in previous passes
		   are incrementally modified */

		filter.filtmin_r = filtmin_r_pass2;
		filter.filtmin_a = filtmin_a_pass2;
		filter.target_excess = excess_fraction * 
					(entries_r + entries_a);

		extra_needed = nfs_purge_singletons(obj, &fb, &filter, 0);

		/* preceed to the merge phase if there is only a small
		   amount of excess left */

		if (filter.num_relations < filter.num_ideals ||
		    filter.num_relations - filter.num_ideals <=
				1.10 * (entries_r + entries_a)) {

			logprintf(obj, "dataset has %3.1f%% excess relations\n",
					100.0 * (
					(double)(filter.num_relations - 
					filter.num_ideals) / (entries_r + 
					entries_a) - 1.0));

			extra_needed = entries_r_pass2 + entries_a_pass2;
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

		if (++num_passes > 8) {
			printf("error: too many filtering passes\n");
			exit(-1);
		}
	}

	nfs_merge_post(obj, &merge);
	dump_cycles(obj, &merge);
	free_relsets(&merge);
	return 0;
}
