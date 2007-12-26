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

/* Create <savefile_name>.s, a binary file containing the line
   numbers of relations that survived the singleton and clique
   removal phases.

   Removal of singletons happens in one of two modes, disk-based
   mode and in-memory mode. The disk-based mode uses a constant-size
   hashtable of bytes to store the counts of large ideals. Large
   ideals that hash to the same bin add their counts together, and
   singleton removal proceeds based on the hashtable until the number
   of singletons found decreases below a cutoff. This allows the first
   few singleton removal passes to run with constant (low) memory use, but
   some singletons may be missed.

   The in-memory mode runs either by itself or after the disk-based
   removal process completes. All the surviving large ideals get read
   in and assigned unique integers, and we use a perfect hash
   instead of an ordinary hash for the remaining singleton removal
   passes. In-memory singleton removal is rigorous, and using a perfect
   hashtable is much faster and more memory-efficient than the on-disk
   hashtable. The only drawback to in-memory mode is the memory use
   associated with the initial reading-in of large ideals */

/*--------------------------------------------------------------------*/
#define NUM_IDEAL_BINS 6

static uint8 *purge_singletons_pass1(msieve_obj *obj, factor_base_t *fb,
					filter_t *filter, 
					uint32 log2_hashtable_size,
					uint32 hash_mask, uint32 pass,
					uint32 *hash_bins_filled_out) {

	/* handle the initial phase of disk-based singleton
	   removal. This fills the hashtable initially */

	FILE *savefile_fp;
	FILE *good_relation_fp;
	uint32 i;
	char buf[256];
	uint32 next_relation;
	uint32 curr_relation;
	uint32 num_relations;
	uint8 *hashtable;
	uint32 hash_bins_filled;
	uint32 ideal_count_bins[NUM_IDEAL_BINS+2] = {0};

	uint32 tmp_factors[TEMP_FACTOR_LIST_SIZE];
	relation_t tmp_relation;

	tmp_relation.factors = tmp_factors;

	logprintf(obj, "commencing singleton removal, pass %u\n", pass);

	savefile_fp = fopen(obj->savefile_name, "r");
	if (savefile_fp == NULL) {
		logprintf(obj, "error: singleton1 can't open savefile\n");
		exit(-1);
	}

	if (pass == 1)
		sprintf(buf, "%s.d", obj->savefile_name);
	else
		sprintf(buf, "%s.s", obj->savefile_name);
	good_relation_fp = fopen(buf, "rb");
	if (good_relation_fp == NULL) {
		logprintf(obj, "error: singleton1 can't open rel file\n");
		exit(-1);
	}
	hashtable = (uint8 *)xcalloc(
			(size_t)1 << log2_hashtable_size, (size_t)1);

	/* for each relation declared good by the duplicate
	   removal phase */

	curr_relation = (uint32)(-1);
	num_relations = 0;
	hash_bins_filled = 0;
	fgets(buf, (int)sizeof(buf), savefile_fp);
	fread(&next_relation, (size_t)1, sizeof(uint32), good_relation_fp);
	while (!feof(savefile_fp) && !feof(good_relation_fp)) {
		
		uint32 num_ideals;
		relation_lp_t tmp_ideal;

		if (buf[0] != '-' && !isdigit(buf[0])) {
			fgets(buf, (int)sizeof(buf), savefile_fp);
			continue;
		}
		if (++curr_relation < next_relation) {
			fgets(buf, (int)sizeof(buf), savefile_fp);
			continue;
		}

		/* read it in */

		nfs_read_relation(buf, fb, &tmp_relation, 1);
		num_relations++;

		/* tabulate the large ideals, and skip the
		   relation if there are too many of them */

		num_ideals = find_large_ideals(&tmp_relation, 
				&tmp_ideal, filter->filtmin);

		/* increment the number of occurrences of each
		   large ideal in the hashtable. Each bin will
		   saturate at a count of 255 or more */

		for (i = 0; i < num_ideals; i++) {
			ideal_t *ideal = tmp_ideal.ideal_list + i;
			uint32 hashval = (HASH1(ideal->blob[0] ^ hash_mask) ^
					  HASH2(ideal->blob[1] ^ hash_mask)) >>
					(32 - log2_hashtable_size);

			if (hashtable[hashval] == 0)
				hash_bins_filled++;
			if (hashtable[hashval] < 0xff)
				hashtable[hashval]++;
		}

		/* update the histogram of ideal counts */

		if (num_ideals > NUM_IDEAL_BINS)
			ideal_count_bins[NUM_IDEAL_BINS+1]++;
		else
			ideal_count_bins[num_ideals]++;

		/* get next line number to look for and the
		   next line from the savefile */

		fread(&next_relation, (size_t)1, 
				sizeof(uint32), good_relation_fp);
		fgets(buf, (int)sizeof(buf), savefile_fp);
	}

	/* print some statistics */

	for (i = 0; i < NUM_IDEAL_BINS+1; i++) {
		logprintf(obj, "relations with %u large ideals: %u\n", 
					i, ideal_count_bins[i]);
	}
	logprintf(obj, "relations with %u+ large ideals: %u\n", 
				i, ideal_count_bins[i]);
	logprintf(obj, "%u relations and about %u large ideals\n", 
				num_relations, hash_bins_filled);

	fclose(savefile_fp);
	fclose(good_relation_fp);
	*hash_bins_filled_out = hash_bins_filled;
	return hashtable;
}

/*--------------------------------------------------------------------*/
static void purge_singletons_pass2(msieve_obj *obj, factor_base_t *fb,
				uint8 *hashtable, filter_t *filter,
				uint32 log2_hashtable_size, uint32 hash_mask,
				uint32 pass) {

	/* continue the disk-based singleton removal process */

	FILE *savefile_fp;
	FILE *good_relation_fp;
	FILE *out_fp;
	uint32 i;
	char buf[256];
	char buf2[256];
	uint32 next_relation;
	uint32 curr_relation;
	uint32 num_relations;
	uint32 num_singletons;
	uint32 num_ideals;

	uint32 tmp_factors[TEMP_FACTOR_LIST_SIZE];
	relation_t tmp_relation;

	tmp_relation.factors = tmp_factors;

	logprintf(obj, "commencing singleton removal, pass %u\n", pass);

	/* we start from either the file from the duplicate
	   removal or the file from the previous singleton
	   removal pass */

	savefile_fp = fopen(obj->savefile_name, "r");
	if (savefile_fp == NULL) {
		logprintf(obj, "error: singleton2 can't open savefile\n");
		exit(-1);
	}
	if (pass == 2)
		sprintf(buf, "%s.d", obj->savefile_name);
	else
		sprintf(buf, "%s.s", obj->savefile_name);
	good_relation_fp = fopen(buf, "rb");
	if (good_relation_fp == NULL) {
		logprintf(obj, "error: singleton2 can't open rel file\n");
		exit(-1);
	}
	sprintf(buf, "%s.s0", obj->savefile_name);
	out_fp = fopen(buf, "wb");
	if (good_relation_fp == NULL) {
		logprintf(obj, "error: singleton2 can't open out file\n");
		exit(-1);
	}

	/* for each relation */

	curr_relation = (uint32)(-1);
	num_relations = 0;
	num_singletons = 0;
	fgets(buf, (int)sizeof(buf), savefile_fp);
	fread(&next_relation, (size_t)1, 
			sizeof(uint32), good_relation_fp);
	while (!feof(savefile_fp) && !feof(good_relation_fp)) {
		
		uint32 is_singleton;
		uint32 num_ideals;
		relation_lp_t tmp_ideal;

		if (buf[0] != '-' && !isdigit(buf[0])) {
			fgets(buf, (int)sizeof(buf), savefile_fp);
			continue;
		}
		if (++curr_relation < next_relation) {
			fgets(buf, (int)sizeof(buf), savefile_fp);
			continue;
		}

		/* read it in */

		nfs_read_relation(buf, fb, &tmp_relation, 1);
		num_relations++;

		/* find the large ideals */

		num_ideals = find_large_ideals(&tmp_relation, 
					&tmp_ideal, 
					filter->filtmin);

		/* check the count of each large ideal */

		for (i = is_singleton = 0; i < num_ideals; i++) {
			ideal_t *ideal = tmp_ideal.ideal_list + i;
			uint32 hashval = (HASH1(ideal->blob[0] ^ hash_mask) ^
					  HASH2(ideal->blob[1] ^ hash_mask)) >>
					(32 - log2_hashtable_size);

			/* cache the hash values in case they 
			   are needed later */

			tmp_factors[i] = hashval;
			if (hashtable[hashval] <= 1)
				is_singleton = 1;
		}

		if (is_singleton) {

			/* decrement the frequency of all of the 
			   large ideals (unless a counter was 
			   saturated) */

			for (i = 0; i < num_ideals; i++) {
				uint32 hashval = tmp_factors[i];
				if (hashtable[hashval] < 0xff)
					hashtable[hashval]--;
			}
			num_singletons++;
		}
		else {
			/* relation survived this pass */

			fwrite(&curr_relation, (size_t)1, 
				sizeof(uint32), out_fp);
		}

		/* get next line number to look for and the
		   next line from the savefile */

		fread(&next_relation, (size_t)1, 
				sizeof(uint32), good_relation_fp);
		fgets(buf, (int)sizeof(buf), savefile_fp);
	}

	/* find the number of ideals remaining and
	   print some statistics */

	for (i = num_ideals = 0; 
			i < ((uint32)1 << log2_hashtable_size); i++) {
		if (hashtable[i] > 0)
			num_ideals++;
	}

	logprintf(obj, "found %u singletons\n", num_singletons);
	logprintf(obj, "current dataset: %u relations and "
			"about %u large ideals\n",
			num_relations - num_singletons, num_ideals);

	/* clean up and create the next singleton file */

	fclose(savefile_fp);
	fclose(good_relation_fp);
	fclose(out_fp);
	if (pass == 2) {
		sprintf(buf, "%s.d", obj->savefile_name);
		remove(buf);
	}
	sprintf(buf, "%s.s0", obj->savefile_name);
	sprintf(buf2, "%s.s", obj->savefile_name);
	remove(buf2);
	if (rename(buf, buf2) != 0) {
		logprintf(obj, "error: singleton2 can't rename output file\n");
		exit(-1);
	}

	filter->num_relations = num_relations - num_singletons;
}

/*--------------------------------------------------------------------*/
static void purge_singletons_final_core(msieve_obj *obj, 
					filter_t *filter) {

	/* main routine for performing in-memory singleton
	   removal. We iterate until there are no more singletons */

	uint32 i, j;
	uint32 *freqtable;
	relation_ideal_t *relation_array;
	relation_ideal_t *curr_relation;
	relation_ideal_t *old_relation;
	uint32 orig_num_ideals;
	uint32 num_passes;
	uint32 num_relations;
	uint32 num_ideals;
	uint32 new_num_relations;

	logprintf(obj, "commencing in-memory singleton removal\n");

	num_relations = filter->num_relations;
	orig_num_ideals = num_ideals = filter->num_ideals;
	relation_array = filter->relation_array;
	freqtable = (uint32 *)xcalloc((size_t)num_ideals, sizeof(uint32));

	/* count the number of times each ideal occurs. Note
	   that since we know the exact number of ideals, we
	   don't need a hashtable to store the counts, just an
	   ordinary random-access array (i.e. a perfect hashtable) */

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {
		for (j = 0; j < curr_relation->ideal_count; j++) {
			uint32 ideal = curr_relation->ideal_list[j];
			freqtable[ideal]++;
		}
		curr_relation = next_relation_ptr(curr_relation);
	}

	logprintf(obj, "begin with %u relations and %u unique ideals\n", 
					num_relations, num_ideals);

	/* while singletons were found */

	num_passes = 0;
	new_num_relations = num_relations;
	do {
		num_relations = new_num_relations;
		new_num_relations = 0;
		curr_relation = relation_array;
		old_relation = relation_array;

		for (i = 0; i < num_relations; i++) {
			uint32 curr_num_ideals = curr_relation->ideal_count;
			uint32 ideal;
			relation_ideal_t *next_relation;

			next_relation = next_relation_ptr(curr_relation);

			/* check the count of each ideal */

			for (j = 0; j < curr_num_ideals; j++) {
				ideal = curr_relation->ideal_list[j];
				if (freqtable[ideal] <= 1)
					break;
			}

			if (j < curr_num_ideals) {

				/* relation is a singleton; decrement the
				   count of each of its ideals and skip it */

				for (j = 0; j < curr_num_ideals; j++) {
					ideal = curr_relation->ideal_list[j];
					freqtable[ideal]--;
				}
			}
			else {
				/* relation survived this pass; append it to
				   the list of survivors */

				old_relation->rel_index = 
						curr_relation->rel_index;
				old_relation->gf2_factors = 
						curr_relation->gf2_factors;
				old_relation->ideal_count = curr_num_ideals;
				for (j = 0; j < curr_num_ideals; j++) {
					old_relation->ideal_list[j] =
						curr_relation->ideal_list[j];
				}
				new_num_relations++;
				old_relation = next_relation_ptr(old_relation);
			}

			curr_relation = next_relation;
		}

		num_passes++;
	} while (new_num_relations != num_relations);

	/* find the ideal that occurs in the most
	   relations, and renumber the ideals to ignore
	   any that have a count of zero */

	num_ideals = 0;
	for (i = j = 0; i < orig_num_ideals; i++) {
		if (freqtable[i]) {
			j = MAX(j, freqtable[i]);
			freqtable[i] = num_ideals++;
		}
	}

	logprintf(obj, "reduce to %u relations and %u ideals in %u passes\n", 
				num_relations, num_ideals, num_passes);
	logprintf(obj, "max relations containing the same ideal: %u\n", j);
	
	/* save the current state */

	filter->num_relations = num_relations;
	filter->num_ideals = num_ideals;
	filter->max_ideal_degree = j;
	filter->relation_array = relation_array = 
			(relation_ideal_t *)xrealloc(relation_array,
				(curr_relation - relation_array + 1) *
				sizeof(relation_ideal_t));

	curr_relation = relation_array;
	for (i = 0; i < num_relations; i++) {
		for (j = 0; j < curr_relation->ideal_count; j++) {
			uint32 ideal = curr_relation->ideal_list[j];
			curr_relation->ideal_list[j] = freqtable[ideal];
		}
		curr_relation = next_relation_ptr(curr_relation);
	}
	free(freqtable);
}

/*--------------------------------------------------------------------*/
static void purge_singletons_final(msieve_obj *obj, 
				factor_base_t *fb,
				filter_t *filter,
				uint32 log2_hashtable_size) {

	/* the last disk-based pass through the relation
	   file; its job is to form a packed array of 
	   relation_ideal_t structures */

	uint32 i;
	FILE *savefile_fp;
	FILE *good_relation_fp;
	char buf[256];
	uint32 next_relation;
	uint32 curr_relation;
	uint32 num_relations;
	hashtable_t unique_ideals;

	relation_ideal_t *relation_array;
	uint32 relation_array_alloc;
	uint32 relation_array_word;

	uint32 tmp_factors[TEMP_FACTOR_LIST_SIZE];
	relation_t tmp_relation;

	tmp_relation.factors = tmp_factors;

	logprintf(obj, "commencing singleton removal, final pass\n");

	savefile_fp = fopen(obj->savefile_name, "r");
	if (savefile_fp == NULL) {
		logprintf(obj, "error: singleton3 can't open savefile\n");
		exit(-1);
	}
	sprintf(buf, "%s.s", obj->savefile_name);
	good_relation_fp = fopen(buf, "rb");
	if (good_relation_fp == NULL) {
		logprintf(obj, "error: singleton3 can't open rel file\n");
		exit(-1);
	}
	relation_array_word = 0;
	relation_array_alloc = 5000;
	relation_array = (relation_ideal_t *)xmalloc(relation_array_alloc *
						sizeof(relation_ideal_t));

	hashtable_init(&unique_ideals, log2_hashtable_size,
			10000, (uint32)(sizeof(ideal_t)/sizeof(uint32)));

	/* for each relation that survived previous 
	   singleton filtering passes */

	curr_relation = (uint32)(-1);
	num_relations = 0;
	fgets(buf, (int)sizeof(buf), savefile_fp);
	fread(&next_relation, (size_t)1, sizeof(uint32), good_relation_fp);
	while (!feof(savefile_fp) && !feof(good_relation_fp)) {
		
		int32 status;

		if (buf[0] != '-' && !isdigit(buf[0])) {
			fgets(buf, (int)sizeof(buf), savefile_fp);
			continue;
		}
		if (++curr_relation < next_relation) {
			fgets(buf, (int)sizeof(buf), savefile_fp);
			continue;
		}

		/* read it in */

		status = nfs_read_relation(buf, fb, &tmp_relation, 1);

		if (status == 0) {
			relation_lp_t tmp_ideal;
			relation_ideal_t packed_ideal;
			num_relations++;

			/* get the large ideals */

			find_large_ideals(&tmp_relation, &tmp_ideal, 
						filter->filtmin);

			packed_ideal.rel_index = curr_relation;
			packed_ideal.gf2_factors = tmp_ideal.gf2_factors;
			packed_ideal.ideal_count = tmp_ideal.ideal_count;

			/* map each ideal to a unique integer */

			for (i = 0; i < tmp_ideal.ideal_count; i++) {
				ideal_t *ideal = tmp_ideal.ideal_list + i;
				hash_t *entry = hashtable_find(&unique_ideals,
								ideal->blob,
								NULL);
				packed_ideal.ideal_list[i] = 
						hashtable_get_offset(
							&unique_ideals, entry);
			}

			/* make sure the relation array has room for the
			   new relation. Be careful increasing the array
			   size, since this is probably the largest array
			   in the NFS code */

			if (relation_array_word * sizeof(uint32) >= 
					(relation_array_alloc-1) * 
					sizeof(relation_ideal_t)) {
				relation_array_alloc = 1.4 * 
						relation_array_alloc;
				relation_array = (relation_ideal_t *)xrealloc(
						relation_array, 
						relation_array_alloc *
						sizeof(relation_ideal_t));
			}

			/* copy exactly the number of bytes that the new
			   relation needs, then point beyond the end of
			   the copied relation */

			memcpy((uint32 *)relation_array + relation_array_word,
					&packed_ideal, sizeof(packed_ideal));
			relation_array_word += (sizeof(packed_ideal) -
					(TEMP_FACTOR_LIST_SIZE - 
					 	packed_ideal.ideal_count) * 
					sizeof(packed_ideal.ideal_list[0])) /
					sizeof(uint32);
		}

		/* get next line number to look for and the
		   next line from the savefile */

		fread(&next_relation, (size_t)1, 
				sizeof(uint32), good_relation_fp);
		fgets(buf, (int)sizeof(buf), savefile_fp);
	}

	logprintf(obj, "memory use: %.1f MB\n", (double)
			(relation_array_alloc * sizeof(relation_ideal_t) +
			 (sizeof(uint32) << log2_hashtable_size) +
			 (unique_ideals.match_array_alloc * sizeof(hash_t))) /
			 1048576);

	fclose(savefile_fp);
	fclose(good_relation_fp);

	/* save the data */

	filter->num_relations = num_relations;
	filter->num_ideals = hashtable_getall(&unique_ideals, NULL);
	hashtable_free(&unique_ideals);
	filter->relation_array = (relation_ideal_t *)xrealloc(relation_array,
						relation_array_word * 
						sizeof(uint32));
}

/*--------------------------------------------------------------------*/
static void dump_relations(msieve_obj *obj, filter_t *filter) {

	uint32 i;
	char buf[256];
	FILE *relation_fp;
	relation_ideal_t *r = filter->relation_array;

	sprintf(buf, "%s.s", obj->savefile_name);
	relation_fp = fopen(buf, "wb");
	if (relation_fp == NULL) {
		logprintf(obj, "error: reldump can't open out file\n");
		exit(-1);
	}

	for (i = 0; i < filter->num_relations; i++) {
		fwrite(&r->rel_index, (size_t)1, sizeof(uint32), relation_fp);
		r = next_relation_ptr(r);
	}
	fclose(relation_fp);
}

/*--------------------------------------------------------------------*/
uint32 nfs_purge_singletons(msieve_obj *obj, factor_base_t *fb,
				filter_t *filter, uint32 disk_based) {

	uint32 clique_heap_size;
	uint32 hash_mask = 0x0f0f0f0f;
	uint32 log2_hashtable_size = 27;

	logprintf(obj, "filtering ideals above %u\n", filter->filtmin);
	logprintf(obj, "need %u more relations than ideals\n", 
					filter->target_excess);
	
	if (disk_based) {
		uint32 hash_bins_filled = 0;
		uint32 curr_pass = 1;
		uint8 *hashtable = purge_singletons_pass1(obj, fb, filter, 
							log2_hashtable_size,
							hash_mask, curr_pass++, 
							&hash_bins_filled);
		uint32 num_relations;

		/* perform disk-based singleton removal until the
		   number of relations is small or the number of 
		   singletons found becomes small */

		filter->num_relations = 0xffffffff;
		do {
			num_relations = filter->num_relations;
			purge_singletons_pass2(obj, fb, hashtable, 
						filter, log2_hashtable_size,
						hash_mask, curr_pass++);

			/* if the hashtable was saturated going 
			   into the above singleton removal pass, 
			   rebuild it with more entries and change
			   the hash function */

			if (hash_bins_filled > 
				0.3 * (1 << log2_hashtable_size)) {
				free(hashtable);
				log2_hashtable_size = MIN(30,
						log2_hashtable_size + 1);
				hash_mask = (hash_mask << 1) |
						(hash_mask >> 31);
				hashtable = purge_singletons_pass1(
							obj, fb, filter, 
							log2_hashtable_size,
							hash_mask, curr_pass++, 
							&hash_bins_filled);
			}
		} while (filter->num_relations > 200000 &&
			num_relations - filter->num_relations > 500000);

		free(hashtable);
	}

	purge_singletons_final(obj, fb, filter, log2_hashtable_size - 6);
	purge_singletons_final_core(obj, filter);

	/* make sure there are enough excess relations
	   to proceed with the rest of the filtering.
	   If there aren't enough relations, return a
	   crude estimate of how many more would be needed */

	if (filter->num_relations < filter->num_ideals ||
	    filter->num_relations - filter->num_ideals <
	    			filter->target_excess) {
		uint32 relations_needed = 1000000;

		if (filter->num_relations > filter->num_ideals) {
			relations_needed = 3 * (filter->target_excess -
						(filter->num_relations - 
					 	 filter->num_ideals));
			relations_needed = MAX(relations_needed, 20000);
		}
		else {
			free(filter->relation_array);
			filter->relation_array = NULL;
		}

		return relations_needed;
	}

	/* iteratively delete cliques until there are
	   just enough excess relations for the merge 
	   phase to run. We choose the number of cliques
	   to delete so that there are always several
	   passes to perform, allowing us to discover new
	   cliques to prune */

	clique_heap_size = ((filter->num_relations - filter->num_ideals) -
				filter->target_excess) / 2;
	clique_heap_size = MIN(clique_heap_size, 400000);

	while (filter->num_relations - filter->num_ideals >
			filter->target_excess + 100) {
		if (nfs_purge_cliques(obj, filter, 
				clique_heap_size,
				(filter->num_relations - filter->num_ideals) - 
				filter->target_excess) == 0) {
			break;
		}

		/* the above got rid of relations; now get rid
		   of the ideals from those relations and remove
		   any additional singletons */

		purge_singletons_final_core(obj, filter);
	}

	dump_relations(obj, filter);
	return 0;
}
