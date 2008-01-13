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

#ifndef _FASTMULT_H_
#define _FASTMULT_H_

#include <util.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAX_FHT_POWER 27

#define HUGE_TWIDDLE_CUTOFF 20

typedef struct {
        double *small;
        double *large;
} huge_twiddle_t;

typedef struct {
	uint32 log2_runlength;
	volatile double round_constant[2];

	double *twiddle[HUGE_TWIDDLE_CUTOFF + 1];

	huge_twiddle_t huge_twiddle[MAX_FHT_POWER + 1 - HUGE_TWIDDLE_CUTOFF];
} fastmult_info_t;

void fastmult_info_init(fastmult_info_t *info);
void fastmult_info_free(fastmult_info_t *info);

void fastmult(uint32 *a, uint32 awords, 
		uint32 *b, uint32 bwords,
		uint32 *prod, fastmult_info_t *info);

#ifdef __cplusplus
}
#endif

#endif /* _FASTMULT_H_ */
