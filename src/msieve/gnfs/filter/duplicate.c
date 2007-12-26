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

/* produce <savefile_name>.d, a binary file containing the
   line numbers of relations in the savefile that are not 
   duplicates.

   This code has to touch all of the relations, and when the
   dataset is large avoiding excessive memory use is tricky.
   Cavallar's paper suggests dropping relations into a hashtable
   and ignoring relations that map to the same hash bin. This keeps
   memory use down but causes good relations to be thrown away.
   Another option is to put the (a,b) values of relations into the
   hashtable, so that hash collisions can be resolved rigorously.
   Unfortunately this means we have to budget 12 bytes for each
   unique relation, and there could be tens of millions of them.

   The implementation here is a compromise: we do duplicate removal
   in two passes. The first pass maps relations into a hashtable
   of bits, and we save (on disk) the list of hash bins where two 
   or more relations collide. The second pass refills the hashtable
   of bits with just these entries, then reads through the complete
   dataset again and saves the (a,b) values of any relation that
   maps to one of the filled-in hash bins. The memory use in the
   first pass is constant, and the memory use of the second pass
   is 12 bytes per duplicate relation. Assuming unique relations
   greatly outnumber duplicates, this solution finds all the duplicates
   with no false positives, and the memory use is low enough so
   that singleton filtering is a larger memory bottleneck 
   
   One useful optimization for really big problems would turn the
   first-pass hashtable into a Bloom filter using several hash
   functions. This would make it much more effective at avoiding
   false positives as the hashtable gets more congested */

#define LOG2_DUP_HASHTABLE1_SIZE 28   /* first pass hashtable size */
#define LOG2_DUP_HASHTABLE2_SIZE 20   /* second pass hashtable size */

static void purge_duplicates_pass2(msieve_obj *obj) {

	FILE *savefile_fp;
	FILE *good_relation_fp;
	FILE *collision_fp;
	FILE *out_fp;
	uint32 i;
	char buf[256];
	uint32 num_duplicates;
	uint32 num_relations;
	uint32 next_relation;
	uint32 curr_relation;
	uint8 *bit_table;
	hashtable_t duplicates;

	logprintf(obj, "commencing duplicate removal, pass 2\n");

	/* fill in the list of hash collisions */

	sprintf(buf, "%s.hc", obj->savefile_name);
	collision_fp = fopen(buf, "rb");
	if (collision_fp == NULL) {
		logprintf(obj, "error: dup2 can't open collision file\n");
		exit(-1);
	}
	bit_table = (uint8 *)xcalloc(
			(size_t)1 << (LOG2_DUP_HASHTABLE1_SIZE - 3), 
			(size_t)1);

	while (fread(&i, (size_t)1, sizeof(uint32), collision_fp) != 0) {
		if (i < (1 << LOG2_DUP_HASHTABLE1_SIZE)) {
			bit_table[i / 8] |= 1 << (i % 8);
		}
	}
	fclose(collision_fp);

	/* set up for reading the list of relations */

	savefile_fp = fopen(obj->savefile_name, "r");
	if (savefile_fp == NULL) {
		logprintf(obj, "error: dup2 can't open savefile\n");
		exit(-1);
	}
	sprintf(buf, "%s.gr", obj->savefile_name);
	good_relation_fp = fopen(buf, "rb");
	if (good_relation_fp == NULL) {
		logprintf(obj, "error: dup2 can't open rel file\n");
		exit(-1);
	}
	sprintf(buf, "%s.d", obj->savefile_name);
	out_fp = fopen(buf, "wb");
	if (out_fp == NULL) {
		logprintf(obj, "error: dup2 can't open output file\n");
		exit(-1);
	}
	hashtable_init(&duplicates, LOG2_DUP_HASHTABLE2_SIZE, 10000, 2);

	num_duplicates = 0;
	num_relations = 0;
	curr_relation = (uint32)(-1);
	fgets(buf, (int)sizeof(buf), savefile_fp);
	fread(&next_relation, (size_t)1, sizeof(uint32), good_relation_fp);
	while (!feof(savefile_fp) && !feof(good_relation_fp)) {
		
		uint32 hashval;
		int64 a;
		uint32 b;
		uint32 key[2];
		char *next_field;

		if (buf[0] != '-' && !isdigit(buf[0])) {

			/* no relation on this line */

			fgets(buf, (int)sizeof(buf), savefile_fp);
			continue;
		}
		if (++curr_relation < next_relation) {

			/* this relation isn't valid */

			fgets(buf, (int)sizeof(buf), savefile_fp);
			continue;
		}

		/* determine if the (a,b) coordinates of the
		   relation collide in the table of bits */

		a = (int64)strtod(buf, &next_field);
		b = strtoul(next_field + 1, NULL, 10);
		key[0] = (uint32)a;
		key[1] = ((a >> 32) & 0xff) | (b << 8);

		hashval = (HASH1(key[0]) ^ HASH2(key[1])) >>
				(32 - LOG2_DUP_HASHTABLE1_SIZE);

		if (bit_table[hashval/8] & (1 << (hashval % 8))) {

			/* relation collides in the first hashtable;
			   use the second hashtable to determine 
			   rigorously if the relation was previously seen */

			uint32 is_dup;
			hashtable_find(&duplicates, key, &is_dup);

			if (!is_dup) {

				/* relation was seen for the first time;
				   doesn't count as a duplicate */

				fwrite(&curr_relation, (size_t)1, 
						sizeof(uint32), out_fp);
				num_relations++;
			}
			else {
				/* relation was previously seen; this
				   time it's a duplicate */

				num_duplicates++;
			}
		}
		else {
			/* no collision; relation is unique */

			fwrite(&curr_relation, (size_t)1, 
					sizeof(uint32), out_fp);
			num_relations++;
		}

		/* get next line number to look for, and
		   next line of relation file */

		fread(&next_relation, (size_t)1, 
				sizeof(uint32), good_relation_fp);
		fgets(buf, (int)sizeof(buf), savefile_fp);
	}

	logprintf(obj, "found %u duplicates and %u unique relations\n", 
				num_duplicates, num_relations);
	logprintf(obj, "memory use: %.1f MB\n", 
			(double)((1 << (LOG2_DUP_HASHTABLE1_SIZE-3)) +
			(sizeof(uint32) << LOG2_DUP_HASHTABLE2_SIZE) +
			duplicates.match_array_alloc * sizeof(hash_t)) / 
			1048576);

	/* clean up and finish */

	fclose(savefile_fp);
	fclose(good_relation_fp);
	fclose(out_fp);
	sprintf(buf, "%s.hc", obj->savefile_name);
	remove(buf);
	sprintf(buf, "%s.gr", obj->savefile_name);
	remove(buf);

	free(bit_table);
	hashtable_free(&duplicates);
}

/*--------------------------------------------------------------------*/
uint32 nfs_purge_duplicates(msieve_obj *obj, factor_base_t *fb,
				uint32 max_relations) {

	uint32 i;
	FILE *savefile_fp;
	FILE *good_relation_fp;
	FILE *collision_fp;
	uint32 curr_relation;
	char buf[256];
	uint32 num_relations;
	uint32 num_collisions;
	uint8 *hashtable;
	uint32 blob[2];

	double total_primes;
	uint32 num_primes;

	uint32 tmp_factors[TEMP_FACTOR_LIST_SIZE];
	relation_t tmp_relation;

	tmp_relation.factors = tmp_factors;

	logprintf(obj, "commencing duplicate removal, pass 1\n");

	savefile_fp = fopen(obj->savefile_name, "r");
	if (savefile_fp == NULL) {
		logprintf(obj, "error: dup1 can't open savefile\n");
		exit(-1);
	}
	sprintf(buf, "%s.gr", obj->savefile_name);
	good_relation_fp = fopen(buf, "wb");
	if (good_relation_fp == NULL) {
		logprintf(obj, "error: dup1 can't open relation file\n");
		exit(-1);
	}
	sprintf(buf, "%s.hc", obj->savefile_name);
	collision_fp = fopen(buf, "wb");
	if (collision_fp == NULL) {
		logprintf(obj, "error: dup1 can't open collision file\n");
		exit(-1);
	}
	hashtable = (uint8 *)xcalloc(
			(size_t)1 << (LOG2_DUP_HASHTABLE1_SIZE - 3), 
			(size_t)1);

	total_primes = 0;
	num_primes = 0;
	curr_relation = (uint32)(-1);
	num_relations = 0;
	num_collisions = 0;
	fgets(buf, (int)sizeof(buf), savefile_fp);
	while (!feof(savefile_fp)) {

		int32 status;
		uint32 hashval;

		if (buf[0] != '-' && !isdigit(buf[0])) {

			/* no relation on this line */

			fgets(buf, (int)sizeof(buf), savefile_fp);
			continue;
		}

		/* read and verify the relation */

		curr_relation++;
		if (max_relations && curr_relation >= max_relations)
			break;

		status = nfs_read_relation(buf, fb, &tmp_relation, 1);
		if (status != 0) {
			logprintf(obj, "error %d reading relation %u\n",
					status, curr_relation);
			fgets(buf, (int)sizeof(buf), savefile_fp);
			continue;
		}

		/* relation is good; find the value to which it
		   hashes. Note that only the bottom 40 bits of 'a'
		   and the bottom 24 bits of 'b' figure into the hash,
		   so that spurious hash collisions are possible
		   (though highly unlikely) */

		num_relations++;
		blob[0] = (uint32)tmp_relation.a;
		blob[1] = ((tmp_relation.a >> 32) & 0xff) |
			  (tmp_relation.b << 8);

		hashval = (HASH1(blob[0]) ^ HASH2(blob[1])) >>
			   (32 - LOG2_DUP_HASHTABLE1_SIZE);

		/* save the hash bucket if there's a collision */

		if (hashtable[hashval / 8] & (1 << (hashval % 8))) {
			fwrite(&hashval, (size_t)1, 
					sizeof(uint32), collision_fp);
			num_collisions++;
		}
		hashtable[hashval / 8] |= 1 << (hashval % 8);

		/* add the factors of tmp_relation to the counts of primes */
		   
		for (i = 0; i < tmp_relation.num_factors_r +
				tmp_relation.num_factors_a; i++) {
			total_primes += tmp_relation.factors[i];
		}
		num_primes += tmp_relation.num_factors_r +
				tmp_relation.num_factors_a;


		/* save the relation number and get the next line */

		fwrite(&curr_relation, (size_t)1, 
				sizeof(uint32), good_relation_fp);
		fgets(buf, (int)sizeof(buf), savefile_fp);
	}

	free(hashtable);
	fclose(savefile_fp);
	fclose(good_relation_fp);
	fclose(collision_fp);
		
	logprintf(obj, "found %u hash collisions in %u relations\n", 
				num_collisions, num_relations);

	if (num_collisions == 0) {

		/* no duplicates; no second pass is necessary */

		char buf2[256];
		sprintf(buf, "%s.hc", obj->savefile_name);
		remove(buf);
		sprintf(buf, "%s.gr", obj->savefile_name);
		sprintf(buf2, "%s.d", obj->savefile_name);
		if (rename(buf, buf2) != 0) {
			logprintf(obj, "error: dup1 can't rename outfile\n");
			exit(-1);
		}
	}
	else {
		purge_duplicates_pass2(obj);
	}

	/* the large prime cutoff for the rest of the filtering
	   process is arbitrarily chosen as the 'center of mass'
	   of the distribution of prime factors of relations.
	   This seems to choose a reasonable bound, and has the
	   advantage of depending only on the relations that are
	   actually present */

	return (uint32)(total_primes / num_primes);
}
