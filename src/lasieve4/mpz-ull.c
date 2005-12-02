/* mpz-ull.c
  By Jens Franke.
  6/13/04: Hacked up for use in GGNFS by Chris Monico.
                                                                                                                                                                                                             
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
                                                                                                                                                                                                             
  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include <sys/types.h>
#include <limits.h>
#include <gmp.h>

#include "lasieve.h"

#ifdef ULL_NO_UL
static unsigned long have_init = 0;
static mpz_t auxz, auxz2;

#ifndef ULLONG_MAX
	#define ULLONG_MAX (~0ULL)
#endif

/****************************************************/
void mpz_ull_init()
/****************************************************/
{
  if (have_init != 0)
    return;
  mpz_init(auxz);
  mpz_init(auxz2);
  have_init = 1;
}

#define BITS_PER_ULONG (sizeof(unsigned long)*CHAR_BIT)
/****************************************************/
void mpz_set_ull(mpz_t targ, uint64_t src)
/****************************************************/
{
  mpz_set_ui(targ,(uint32_t)(src>>32));
  mpz_mul_2exp(targ, targ, BITS_PER_ULONG);
  mpz_add_ui(targ, targ, (uint32_t)(src&0xFFFFFFFF));
}

/****************************************************/
uint64_t mpz_get_ull(mpz_t src)
/****************************************************/
{
  mpz_fdiv_q_2exp(auxz, src, BITS_PER_ULONG);
  return (((uint64_t)mpz_get_ui(auxz)) << BITS_PER_ULONG) | mpz_get_ui(src);
}

/****************************************************/
int mpz_cmp_ull(mpz_t op1, uint64_t op2)
/****************************************************/
{
  mpz_set_ull(auxz, op2);
  return mpz_cmp(op1, auxz);
}

/****************************************************/
void mpz_mul_ull(mpz_t rop, mpz_t op1, uint64_t op2)
/****************************************************/
{
  mpz_set_ull(auxz, op2);
  mpz_mul(rop, op1, auxz);
}

/****************************************************/
int64_t mpz_get_sll(mpz_t x)
/****************************************************/
{
  if (mpz_sgn(x) < 0) {
    mpz_neg(auxz2, x);
    return -((int64_t) mpz_get_ull(auxz2));
  } else
    return mpz_get_ull(x);
}

/****************************************************/
void mpz_set_sll(mpz_t x, int64_t src)
/****************************************************/
{
  if (src < 0) {
    mpz_set_ull(x, (uint64_t) (-src));
    mpz_neg(x, x);
  } else
    mpz_set_ull(x, (uint64_t) src);
}

/****************************************************/
int mpz_fits_uint64_t_p(mpz_t x)
/****************************************************/
{
  if (mpz_sgn(x) < 0)
    return 0;
  if (mpz_sizeinbase(x, 2) > CHAR_BIT * sizeof(uint64_t))
    return 0;
  return 1;
}

/****************************************************/
int mpz_fits_sllong_p(mpz_t x)
/****************************************************/
{
  if (mpz_sgn(x) > 0)
    return mpz_get_ull(x) < ULLONG_MAX / 2;
  else {
    mpz_neg(auxz2, x);
    return mpz_get_ull(auxz2) <= ULLONG_MAX / 2;
  }
}
#endif
