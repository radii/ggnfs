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

$Id: sqrt.h,v 1.1 2007-12-16 03:54:47 jasonp_sf Exp $
----------------------------------------------------------------------*/

#ifndef _SQRT_H_
#define _SQRT_H_

#include "gnfs.h"
#include <ap.h>

#ifdef __cplusplus
extern "C" {
#endif

uint32 get_prime_for_sqrt(mp_poly_t *alg_poly,
			  uint32 min_value,
			  uint32 *q_out); 

void alg_square_root(msieve_obj *obj, mp_poly_t *monic_alg_poly, 
			mp_t *n, mp_t *c, signed_mp_t *m1, signed_mp_t *m0,
			abpair_t *rlist, uint32 num_relations, uint32 check_q,
			mp_t *sqrt_a);

#ifdef __cplusplus
}
#endif

#endif /* _SQRT_H_ */
