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

	savefile_t *savefile = &obj->savefile;
	FILE *bad_relation_fp;
	uint32 i;
	char buf[LINE_BUF_SIZE];
	uint32 next_bad_relation;
	uint32 curr_relation;
	uint32 num_relations;
	uint8 *hashtable;
	uint32 hash_bins_filled;
	uint32 ideal_count_bins[NUM_IDEAL_BINS+2] = {0};

	uint32 tmp_factors[TEMP_FACTOR_LIST_SIZE];
	relation_t tmp_relation;

	tmp_relation.factors = tmp_factors;

	logprintf(obj, "commencing singleton removal, pass %u\n", pass);

	savefile_open(savefile, SAVEFILE_READ);

	if (pass == 1)
		sprintf(buf, "%s.d", savefile->name);
	else
		sprintf(buf, "%s.s", savefile->name);
	bad_relation_fp = fopen(buf, "rb");
	if (bad_relation_fp == NULL) {
		logprintf(obj, "error: singleton1 can't open rel file\n");
		exit(-1);
	}
	hashtable = (uint8 *)xcalloc(
			(size_t)1 << log2_hashtable_size, (size_t)1);

	/* for each relation declared good by the duplicate
	   removal phase */

	num_relations = 0;
	hash_bins_filled = 0;
	curr_relation = (uint32)(-1);
	next_bad_relation = (uint32)(-1);
	fread(&next_bad_relation, (size_t)1, 
			sizeof(uint32), bad_relation_fp);
	savefile_read_line(buf, sizeof(buf), savefile);

	while (!savefile_eof(savefile)) {
		
		uint32 num_ideals;
		relation_lp_t tmp_ideal;

		if (buf[0] != '-' && !isdigit(buf[0])) {
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}
		if (++curr_relation == next_bad_relation) {
			fread(&next_bad_relation, (size_t)1, 
					sizeof(uint32), bad_relation_fp);
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		/* read it in */

		nfs_read_relation(buf, fb, &tmp_relation, 1);
		num_relations++;

		/* tabulate the large ideals */

		num_ideals = find_large_ideals(&tmp_relation, &tmp_ideal, 
						filter->filtmin_r,
						filter->filtmin_a);

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

		/* get the next line from the savefile */

		savefile_read_line(buf, sizeof(buf), savefile);
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

	savefile_close(savefile);
	fclose(bad_relation_fp);
	*hash_bins_filled_out = hash_bins_filled;
	return hashtable;
}

/*--------------------------------------------------------------------*/
static void purge_singletons_pass2(msieve_obj *obj, factor_base_t *fb,
				uint8 *hashtable, filter_t *filter,
				uint32 log2_hashtable_size, uint32 hash_mask,
				uint32 pass) {

	/* continue the disk-based singleton removal process */

	savefile_t *savefile = &obj->savefile;
	FILE *bad_relation_fp;
	FILE *out_fp;
	uint32 i;
	char buf[LINE_BUF_SIZE];
	char buf2[256];
	uint32 next_bad_relation;
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

	savefile_open(savefile, SAVEFILE_READ);
	if (pass == 2)
		sprintf(buf, "%s.d", savefile->name);
	else
		sprintf(buf, "%s.s", savefile->name);
	bad_relation_fp = fopen(buf, "rb");
	if (bad_relation_fp == NULL) {
		logprintf(obj, "error: singleton2 can't open rel file\n");
		exit(-1);
	}
	sprintf(buf, "%s.s0", savefile->name);
	out_fp = fopen(buf, "wb");
	if (out_fp == NULL) {
		logprintf(obj, "error: singleton2 can't open out file\n");
		exit(-1);
	}

	/* for each relation */

	num_relations = 0;
	num_singletons = 0;
	curr_relation = (uint32)(-1);
	next_bad_relation = (uint32)(-1);
	fread(&next_bad_relation, (size_t)1, 
			sizeof(uint32), bad_relation_fp);
	savefile_read_line(buf, sizeof(buf), savefile);

	while (!savefile_eof(savefile)) {
		
		uint32 is_singleton;
		uint32 num_ideals;
		relation_lp_t tmp_ideal;

		if (buf[0] != '-' && !isdigit(buf[0])) {
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}
		if (++curr_relation == next_bad_relation) {
			fwrite(&curr_relation, (size_t)1, 
					sizeof(uint32), out_fp);
			fread(&next_bad_relation, (size_t)1, 
					sizeof(uint32), bad_relation_fp);
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		/* read it in */

		nfs_read_relation(buf, fb, &tmp_relation, 1);
		num_relations++;

		/* find the large ideals */

		num_ideals = find_large_ideals(&tmp_relation, 
					&tmp_ideal, 
					filter->filtmin_r,
					filter->filtmin_a);

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
			fwrite(&curr_relation, (size_t)1, 
				sizeof(uint32), out_fp);
		}

		savefile_read_line(buf, sizeof(buf), savefile);
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

	savefile_close(savefile);
	fclose(bad_relation_fp);
	fclose(out_fp);
	if (pass == 2) {
		sprintf(buf, "%s.d", savefile->name);
		remove(buf);
	}
	sprintf(buf, "%s.s0", savefile->name);
	sprintf(buf2, "%s.s", savefile->name);
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
static size_t read_lp_file(msieve_obj *obj, filter_t *filter, FILE *fp,
				size_t mem_use, 
				uint32 max_ideal_weight) {
	uint32 i, j, k;
	size_t header_size;
	relation_ideal_t tmp;
	relation_ideal_t *relation_array;
	uint32 curr_word;
	size_t num_relation_alloc;
	uint32 *counts;
	uint32 num_relations = filter->num_relations;
	uint32 num_ideals = filter->num_ideals;

	/* read in the relations from fp, as well as the large
	   ideals they contain. Do not save ideals that occur
	   more than max_ideal_weight times in the dataset */

	header_size = (sizeof(relation_ideal_t) - 
			sizeof(tmp.ideal_list)) / sizeof(uint32);
	counts = (uint32 *)xcalloc((size_t)num_ideals, sizeof(uint32));

	/* first build a frequency table for the large ideals */

	for (i = 0; i < num_relations; i++) {

		fread(&tmp, sizeof(uint32), header_size, fp);

		for (j = 0; j < tmp.ideal_count; j++) {
			uint32 curr_ideal;

			fread(&curr_ideal, sizeof(uint32), (size_t)1, fp);
			counts[curr_ideal]++;
		}
	}

	/* renumber the ideals to ignore the ones that occur
	   too often */

	for (i = j = 0; i < num_ideals; i++) {
		if (counts[i] <= max_ideal_weight)
			counts[i] = j++;
		else
			counts[i] = (uint32)(-1);
	}
	filter->target_excess += i - j;
	filter->num_ideals = j;
	logprintf(obj, "keeping %u ideals with weight <= %u, "
			"new excess is %u\n",
			j, max_ideal_weight, 
			filter->target_excess);

	/* reread the relation list, saving the sparse ideals */

	rewind(fp);
	num_relation_alloc = 10000;
	curr_word = 0;
	relation_array = (relation_ideal_t *)xmalloc(
					num_relation_alloc *
					sizeof(relation_ideal_t));
	for (i = 0; i < num_relations; i++) {

		relation_ideal_t *r;

		/* make sure the relation array has room for the
		   new relation. Be careful increasing the array
		   size, since this is probably the largest array
		   in the NFS code */

		if (curr_word * sizeof(uint32) >=
				(num_relation_alloc-1) * 
				sizeof(relation_ideal_t)) {

			num_relation_alloc = 1.4 * num_relation_alloc;
			relation_array = (relation_ideal_t *)xrealloc(
					relation_array, 
					num_relation_alloc *
					sizeof(relation_ideal_t));
		}

		r = (relation_ideal_t *)(
			(uint32 *)relation_array + curr_word);
		fread(r, sizeof(uint32), header_size, fp);

		for (j = k = 0; j < r->ideal_count; j++) {

			uint32 curr_ideal;

			fread(&curr_ideal, sizeof(uint32), (size_t)1, fp);
			curr_ideal = counts[curr_ideal];
			if (curr_ideal != (uint32)(-1))
				r->ideal_list[k++] = curr_ideal;
		}
		r->gf2_factors += j - k;
		r->ideal_count = k;
		curr_word += header_size + k;
	}

	/* finish up: trim the allocated relation array */

	filter->relation_array = (relation_ideal_t *)xrealloc(
						relation_array, 
						curr_word * 
						sizeof(uint32));
	free(counts);
	return MAX(mem_use, num_ideals * sizeof(uint32) +
			    num_relation_alloc * 
			    sizeof(relation_ideal_t));
}

/*--------------------------------------------------------------------*/
static void purge_singletons_final(msieve_obj *obj, 
				factor_base_t *fb,
				filter_t *filter,
				uint32 log2_hashtable_size,
				uint32 max_ideal_weight) {

	/* the last disk-based pass through the relation
	   file; its job is to form a packed array of 
	   relation_ideal_t structures */

	uint32 i;
	savefile_t *savefile = &obj->savefile;
	FILE *relation_fp;
	FILE *final_fp;
	char buf[LINE_BUF_SIZE];
	uint32 next_relation;
	uint32 curr_relation;
	uint32 num_relations;
	hashtable_t unique_ideals;
	size_t relation_array_words;
	size_t mem_use;
	uint32 tmp_factors[TEMP_FACTOR_LIST_SIZE];
	relation_t tmp_relation;
	uint32 have_bad_relation_list = (max_ideal_weight == 0);

	tmp_relation.factors = tmp_factors;

	logprintf(obj, "commencing singleton removal, final pass\n");

	savefile_open(savefile, SAVEFILE_READ);
	sprintf(buf, "%s.s", savefile->name);
	relation_fp = fopen(buf, "rb");
	if (relation_fp == NULL) {
		logprintf(obj, "error: singleton3 can't open rel file\n");
		exit(-1);
	}
	sprintf(buf, "%s.lp", savefile->name);
	final_fp = fopen(buf, "w+b");
	if (final_fp == NULL) {
		logprintf(obj, "error: singleton3 can't open out file\n");
		exit(-1);
	}
	relation_array_words = 0;

	hashtable_init(&unique_ideals, log2_hashtable_size,
			10000, (uint32)(sizeof(ideal_t)/sizeof(uint32)));

	/* for each relation that survived previous 
	   singleton filtering passes */

	curr_relation = (uint32)(-1);
	next_relation = (uint32)(-1);
	num_relations = 0;
	fread(&next_relation, (size_t)1, 
			sizeof(uint32), relation_fp);
	savefile_read_line(buf, sizeof(buf), savefile);

	while (!savefile_eof(savefile)) {
		
		int32 status;
		size_t curr_relation_words;

		if (buf[0] != '-' && !isdigit(buf[0])) {
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}
		if (have_bad_relation_list) {
			if (++curr_relation == next_relation) {
				fread(&next_relation, (size_t)1, 
						sizeof(uint32), relation_fp);
				savefile_read_line(buf, sizeof(buf), savefile);
				continue;
			}
		}
		else {
			if (++curr_relation < next_relation) {
				savefile_read_line(buf, sizeof(buf), savefile);
				continue;
			}
			fread(&next_relation, (size_t)1, 
					sizeof(uint32), relation_fp);
		}

		/* read it in */

		status = nfs_read_relation(buf, fb, &tmp_relation, 1);

		if (status == 0) {
			relation_lp_t tmp_ideal;
			relation_ideal_t packed_ideal;
			num_relations++;

			/* get the large ideals */

			find_large_ideals(&tmp_relation, &tmp_ideal, 
						filter->filtmin_r,
						filter->filtmin_a);

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

			/* dump the relation to disk */

			curr_relation_words = (sizeof(packed_ideal) -
					(TEMP_FACTOR_LIST_SIZE - 
					 	packed_ideal.ideal_count) * 
					sizeof(packed_ideal.ideal_list[0])) /
					sizeof(uint32);
			fwrite(&packed_ideal, curr_relation_words,
				sizeof(uint32), final_fp);
			relation_array_words += curr_relation_words;
		}

		savefile_read_line(buf, sizeof(buf), savefile);
	}

	savefile_close(savefile);
	fclose(relation_fp);

	/* finish up; the disk-based pass simply reads all
	   the large ideals from disk, but subsequent passes
	   filter out large ideals that occur too often */

	filter->num_relations = num_relations;
	filter->num_ideals = hashtable_getall(&unique_ideals, NULL);
	mem_use = hashtable_sizeof(&unique_ideals);
	hashtable_free(&unique_ideals);

	rewind(final_fp);
	if (max_ideal_weight == 0) {
		filter->relation_array = (relation_ideal_t *)xmalloc(
						relation_array_words * 
						sizeof(uint32));
		fread(filter->relation_array, relation_array_words,
				sizeof(uint32), final_fp);
		mem_use = MAX(mem_use, 
				relation_array_words * sizeof(uint32));
	}
	else {
		mem_use = read_lp_file(obj, filter, final_fp, 
					mem_use, max_ideal_weight);
	}

	logprintf(obj, "memory use: %.1f MB\n",
			(double)mem_use / 1048576);
	fclose(final_fp);
	sprintf(buf, "%s.lp", savefile->name);
	remove(buf);
}

/*--------------------------------------------------------------------*/
static void dump_relations(msieve_obj *obj, filter_t *filter) {

	uint32 i;
	char buf[256];
	FILE *relation_fp;
	relation_ideal_t *r = filter->relation_array;

	sprintf(buf, "%s.s", obj->savefile.name);
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
				filter_t *filter, uint32 max_ideal_weight) {

	uint32 clique_heap_size;
	uint32 hash_mask = 0x0f0f0f0f;
	uint32 log2_hashtable_size = 27;
	uint32 disk_based = (max_ideal_weight == 0);
	uint32 max_clique_relations;
	uint32 target_excess;

	logprintf(obj, "filtering rational ideals above %u\n", 
					filter->filtmin_r);
	logprintf(obj, "filtering algebraic ideals above %u\n", 
					filter->filtmin_a);
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

	purge_singletons_final(obj, fb, filter, 
			log2_hashtable_size - 6, max_ideal_weight);
	purge_singletons_final_core(obj, filter);

	/* save the relation list once it is free of singletons.
	   If the disk-based pass above occurred, this will 
	   convert the list of relations-to-skip into a list
	   of relations-to-keep */

	dump_relations(obj, filter);

	/* make sure there are enough excess relations
	   to proceed with the rest of the filtering.
	   If there aren't enough relations, return a
	   crude estimate of how many more would be needed */

	target_excess = filter->target_excess;
	if (!disk_based)
		target_excess = FINAL_EXCESS_FRACTION * target_excess;

	if (filter->num_relations < filter->num_ideals ||
	    filter->num_relations - filter->num_ideals < target_excess) {
		uint32 relations_needed = 1000000;

		if (filter->num_relations > filter->num_ideals) {
			relations_needed = 3 * (target_excess -
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

	clique_heap_size = ((filter->num_relations - 
				filter->num_ideals) - target_excess) / 2;
	clique_heap_size = MIN(clique_heap_size, 400000);
	max_clique_relations = 2;

	while (1) {
		/* iteratively delete relations containing ideals
		   appearing in exactly max_clique_relations relations */

		while (filter->num_relations - filter->num_ideals >
				target_excess + 100) {
			if (nfs_purge_cliques(obj, filter, clique_heap_size,
					max_clique_relations,
					(filter->num_relations - 
					 filter->num_ideals) - 
					target_excess) == 0) {
				break;
			}

			/* the above got rid of relations; now get rid
			   of the ideals from those relations and remove
			   any additional singletons */

			purge_singletons_final_core(obj, filter);
		}

		/* in the vast majority of cases, max_clique_relations = 2
		   is enough to burn up all of the excess in the dataset.
		   However, when there is a really huge amount of 
		   initial excess, we may run out of cliques before 
		   running out of excess. In that case, look for larger 
		   cliques to delete */

		if (filter->num_relations - filter->num_ideals <=
					target_excess + 100000)
			break;

		logprintf(obj, "too much excess; switching to %u-cliques\n",
					++max_clique_relations);
	}

	/* do not save the list with cliques removed; if we have
	   to rerun the singleton removal, then more excess relations
	   may be required */

	return 0;
}
