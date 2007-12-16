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

$Id: mp_int.h,v 1.1 2007-12-16 03:54:47 jasonp_sf Exp $
----------------------------------------------------------------------*/

#ifndef _MP_INT_H_
#define _MP_INT_H_

#include <mp.h>

#ifdef __cplusplus
extern "C" {
#endif
	/* routines that are used within the multiple precision
	   code but are also useful by themselves and as building
	   blocks in other libraries */

	/* return the index of the first nonzero word in
	   x, searching backwards from max_words */

static INLINE uint32 num_nonzero_words(uint32 *x, uint32 max_words) {

	uint32 i;
	for (i = max_words; i && !x[i-1]; i--)
		;
	return i;
}

	/* Internal structure used by routines that need
	   to do division */

typedef struct {
	uint32 nwords;
	uint32 val[2 * MAX_MP_WORDS];
} big_mp_t;

	/* internal division routine; the input is twice the
	   size of an mp_t to allow for a quotient and remainder
	   each that large. Note that no check is made that
	   the quotient fits in an mp_t */

void mp_divrem_core(big_mp_t *num, mp_t *denom, mp_t *quot, mp_t *rem);

#ifdef __cplusplus
}
#endif

#endif /* !_MP_INT_H_ */
