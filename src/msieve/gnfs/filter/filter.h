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

/* implementation of number field sieve filtering */

#ifndef _FILTER_H_
#define _FILTER_H_

#include <common.h>
#include "gnfs.h"

#ifdef __cplusplus
extern "C" {
#endif

/* the singleton removal phase uses a packed representation
   of each relation. The following maps relations to large ideals;
   the ideals themselves don't matter, only the unique number
   that each is assigned */

typedef struct {
	uint32 rel_index;      /* savefile line number where relation occurs */
	uint8 ideal_count;     /* number of large ideals */
	uint8 gf2_factors;     /* number of small factors that are ignored */
	uint8 connected;       /* scratch space used in clique formation */
	uint32 ideal_list[TEMP_FACTOR_LIST_SIZE];  /* the large ideals */
} relation_ideal_t;

/* relations can have between 0 and TEMP_FACTOR_LIST_SIZE 
   large ideals, and the average number is very small (2 or 3). 
   When processing large numbers of relations we do not 
   assume a full relation_ideal_t struct is allocated for each. 
   Instead all such structures are packed together as
   tightly as possible, and we always iterate sequentially
   through the array. The following abstracts away the process
   of figuring out the number of bytes used in a current
   relation_ideal_t in the list, then pointing beyond it */

static INLINE relation_ideal_t *next_relation_ptr(relation_ideal_t *r) {
	return (relation_ideal_t *)((uint8 *)r +
			sizeof(relation_ideal_t) - 
			sizeof(r->ideal_list[0]) * 
			(TEMP_FACTOR_LIST_SIZE - r->ideal_count));
}

/* structure for the mapping between large ideals 
   and relations (used during clique removal) */

typedef struct {
	uint32 relation_array_word;  /* 32-bit word offset into relation
					array where the relation starts */
	uint32 next;		     /* next relation containing this ideal */
} ideal_relation_t;

/* structure used to map between a large ideal and a 
   linked list of relations that use that ideal */

typedef struct {
	uint32 payload : 30;	/* offset in list of ideal_relation_t
				   structures where the linked list of
				   ideal_relation_t's for this ideal starts */
	uint32 clique : 1;      /* nonzero if this ideal can participate in
				   a clique */
	uint32 connected : 1;   /* nonzero if this ideal has already been
				   added to a clique under construction */
} ideal_map_t;

/* main structure controlling relation filtering */

typedef struct {
	relation_ideal_t *relation_array;  /* relations after singleton phase */
	uint32 num_relations;       /* current number of relations */
	uint32 num_ideals;          /* current number of unique large ideals */
	uint32 filtmin_r;           /* min. value a rational ideal needs 
				       to be tracked during filtering */
	uint32 filtmin_a;           /* min. value an algebraic ideal needs 
				       to be tracked during filtering */
	uint32 target_excess;      /* how many more relations than ideals
					are required for filtering to proceed */
} filter_t;

/* the multiple of the amount of excess needed for
   merging to proceed */

#define FINAL_EXCESS_FRACTION 1.16

/* representation of a 'relation set', i.e. a group
   of relations that will all participate in the same
   column when the final matrix is formed */

typedef struct {
	uint16 num_relations;     /* number of relations in this relation set */
	uint16 num_small_ideals;  /* number of ideals that are not tracked */
	uint16 num_large_ideals;  /* number of ideals that are tracked */
	uint16 num_active_ideals; /* number of ideals eligible for merging
				     (0 means relset is a cycle) */
	uint32 *data;             /* the first num_relations elements are
				     the relation numbers that participate in
				     this relation set, while the last 
				     num_large_ideals elements give the ideals
				     that are tracked for merging. Both lists
				     are assumed sorted in ascending order */
} relation_set_t;

/* the above simulates matrix rows; the following simulates
   the matrix columns, mapping ideals to the relation sets
   containing those ideals */

typedef struct ideal_set_t {
	uint16 num_relsets;     /* the number of relation sets
				   containing this ideal */
	uint16 num_relsets_alloc; /* maximum number of relset numbers the
				     'relsets' array can hold */
	uint16 active;          /* 1 if ideal is active, 0 if inactive */
	uint16 min_relset_size; /* the number of ideals in the
				   lightest relation set that
				   contains this ideal */
	uint32 *relsets;        /* list of members in an array of 
				   relation sets that contain this
				   ideal (no ordering assumed) */
	struct ideal_set_t *next;
	struct ideal_set_t *prev;  /* used to build circular linked lists */
} ideal_set_t;

/* relation sets with more than this many relations are deleted */

#define MAX_RELSET_SIZE 28

/* structure controlling the merge process */

typedef struct {
	uint32 num_relsets;           /* current number of relation sets */
	uint32 num_ideals;            /* number of unique ideals that must
					 be merged */
	double avg_cycle_weight;      /* the avg number of nonzeros per cycle */
	relation_set_t *relset_array; /* current list of relation sets */
} merge_t;


/* create '<savefile_name>.d', a binary file containing
   the line numbers of unique relations. The return value is
   the large prime bound to use for the rest of the filtering.
   Duplicate removal only applies to the first max_relations
   relations found (or all relations if zero) */

uint32 nfs_purge_duplicates(msieve_obj *obj, factor_base_t *fb,
				uint32 max_relations); 

/* read '<savefile_name>.d' and create '<savefile_name>.s', a 
   binary file containing the line numbers of relations that
   are not singletons. All ideals larger than the bounds specified
   in 'filter' are tacked */
   
void nfs_purge_singletons_initial(msieve_obj *obj, 
			factor_base_t *fb, filter_t *filter);

/* read and modify '<savefile_name>.s', to account for ideals
   larger than the filtering bound in 'filter' and occurring in
   at most max_ideal_weight relations (or all relations if zero) */
   
void nfs_purge_singletons(msieve_obj *obj, factor_base_t *fb,
			filter_t *filter, uint32 max_ideal_weight);

/* remove the singletons in a memory-resident set of relations */

void nfs_purge_singletons_core(msieve_obj *obj, filter_t *filter);

/* perform clique removal on the current set of relations */

void nfs_purge_cliques(msieve_obj *obj, filter_t *filter);

/* initialize the merge process */

void nfs_merge_init(msieve_obj *obj, filter_t *filter); 

/* perform all 2-way merges, converting the results into
   relation-sets that the main merge routine operates on */

void nfs_merge_2way(msieve_obj *obj, filter_t *filter, merge_t *merge);

/* do the rest of the merging. min_cycles is the minimum number
   of cycles that the input collection of relation-sets must
   produce, corresponding to the smallest matrix that can be
   built (the actual matrix is expected to be much larger than 
   this). */

void nfs_merge_full(msieve_obj *obj, merge_t *merge, uint32 min_cycles);

/* perform post-processing optimizations on the collection of cycles
   found by the merge phase */

void nfs_merge_post(msieve_obj *obj, merge_t *merge);

#ifdef __cplusplus
}
#endif

#endif /* _FILTER_H_ */
