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

/*------------------------------------------------------------------------*/
#if defined(WIN32) || defined(WIN64)

#include <float.h>
typedef uint32 precision_t;
static INLINE uint32 precision_is_ieee(void) {
	precision_t prec = _control87(0, 0);
	return  ((prec & _MCW_PC) == _PC_53) ? 1 : 0;
}
static INLINE precision_t set_precision(void) {
	precision_t old_prec = _control87(0, 0);
	_control87(_PC_53, _MCW_PC);
	return old_prec;
}
static INLINE void clear_precision(precision_t old_prec) {
	_control87(old_prec, 0xffffffff);
}

#elif (defined(__GNUC__) || defined(__ICL)) && \
	(defined(__i386__) || defined(__x86_64__))

#include <float.h>
typedef uint16 precision_t;
static INLINE uint32 precision_is_ieee(void) {
	precision_t prec;
	asm volatile ("fnstcw %0":"=m"(prec));
	return ((prec & ~0x0300) == 0x0200) ? 1 : 0;
}
static INLINE precision_t set_precision(void) {
	precision_t old_prec, new_prec;
	asm volatile ("fnstcw %0":"=m"(old_prec));
	new_prec = (old_prec & ~0x0300) | 0x0200;
	asm volatile ("fldcw %0": :"m"(new_prec));
	return old_prec;
}
static INLINE void clear_precision(precision_t old_prec) {
	asm volatile ("fldcw %0": :"m"(old_prec));
}

#else

static INLINE uint32 precision_is_ieee(void) {
	return 1;
}
static INLINE precision_t set_precision(void) {
	return 0;
}
static INLINE void clear_precision(precision_t old_prec) {
}

#endif
/*------------------------------------------------------------------------*/

#define MAX_FHT_POWER 27

#define HUGE_TWIDDLE_CUTOFF 20

typedef struct {
        double *small;
        double *large;
} huge_twiddle_t;

typedef struct {
	uint32 log2_runlength;
	volatile double round_constant[2];

	uint32 precision_changed;
	precision_t old_precision;

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
