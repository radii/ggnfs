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

#ifndef _LANCZOS_H_
#define _LANCZOS_H_

#include <common.h>

#ifdef __cplusplus
extern "C" {
#endif

/* routines for cache-efficient multiplication of
   sparse matrices */

/* the smallest number of columns that will be
   converted to packed format */

#define MIN_NCOLS_TO_PACK 50000

/* the number of moderately dense rows that are
   packed less tightly */

#define NUM_MEDIUM_ROWS 3000

/* structure representing a nonzero element of
   the matrix after packing into block format. 
   The two fields are the row and column offsets
   from the top left corner of the block */

typedef struct {
	uint16 row_off;
	uint16 col_off;
} entry_idx_t;

/* struct representing one block */

typedef struct {
	uint32 start_row;
	uint32 start_col;         /* coordinates of top left corner */
	uint32 num_rows;
	uint32 num_entries;       /* number of nonzero matrix entries */
	uint32 num_entries_alloc; /* nonzero matrix entries allocated */
	entry_idx_t *entries;     /* nonzero entries */
	uint16 *med_entries;      /* nonzero entries for medium dense rows */
} packed_block_t;

enum thread_command {
	COMMAND_WAIT,
	COMMAND_RUN,
	COMMAND_RUN_TRANS,
	COMMAND_END
};

/* struct used by threads for computing partial
   matrix multiplies */

typedef struct {
	uint32 ncols;		/* number of columns used by this thread */
	uint32 num_dense_rows;  /* number of rows packed by dense_blocks */
	uint64 **dense_blocks;  /* for holding dense matrix rows; 
				   dense_blocks[i] holds the i_th batch of
				   64 matrix rows */
	uint32 num_blocks;
	uint64 *x;
	uint64 *b;
	packed_block_t *blocks; /* sparse part of matrix, in block format */

	/* fields for thread pool synchronization */

	enum thread_command command;

#if defined(WIN32) || defined(_WIN64)
	HANDLE thread_id;
	HANDLE run_event;
	HANDLE finish_event;
#else
	pthread_t thread_id;
	pthread_mutex_t run_lock;
	pthread_cond_t run_cond;
#endif

} thread_data_t;

#define MAX_THREADS 32
#define MIN_NCOLS_TO_THREAD 250000

/* struct representing a packed matrix */

typedef struct {
	uint32 nrows;
	uint32 ncols;
	uint32 num_dense_rows;
	uint32 num_threads;

	la_col_t *unpacked_cols;  /* used if no packing takes place */

	thread_data_t thread_data[MAX_THREADS];

} packed_matrix_t;

void packed_matrix_init(msieve_obj *obj, 
			packed_matrix_t *packed_matrix,
			la_col_t *A, uint32 nrows, uint32 ncols,
			uint32 num_dense_rows);

void packed_matrix_free(packed_matrix_t *packed_matrix);

size_t packed_matrix_sizeof(packed_matrix_t *packed_matrix);

void mul_MxN_Nx64(packed_matrix_t *A, uint64 *x, uint64 *b);

void mul_trans_MxN_Nx64(packed_matrix_t *A, uint64 *x, uint64 *b);

void mul_Nx64_64x64_acc(uint64 *v, uint64 *x, uint64 *y, uint32 n);

void mul_64xN_Nx64(uint64 *x, uint64 *y, uint64 *xy, uint32 n);

/* for big jobs, we use a multithreaded framework that calls
   these two routines for the heavy lifting */

void mul_packed_core(thread_data_t *t);

void mul_trans_packed_core(thread_data_t *t);

#ifdef __cplusplus
}
#endif

#endif /* !_LANCZOS_H_ */
