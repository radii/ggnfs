/* gmp-aux.c
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
#include <string.h>
#include <gmp.h>
#ifdef __ppc__
#include "ppc32/siever-config.h"
#else
#include "asm/lasieve-asm.h"
#include "lasieve.h"
#endif

/********************************************************/
void adjust_mpz_bufsize(mpz_t ** x, size_t * alloc_ptr, size_t size,
                        size_t increment)
/********************************************************/
{
  size_t old_alloc;

  old_alloc = *alloc_ptr;
  adjust_bufsize((void **)x, alloc_ptr, size, increment, sizeof(**x));
  while (old_alloc < *alloc_ptr)
    mpz_init((*x)[old_alloc++]);
}

/********************************************************/
int string2mpz(mpz_t rop, char *x, int base)
/********************************************************/
{
  size_t l;
  char *y;
  int rv;

  x += strspn(x, " \t+");
  if (strlen(x) == 0)
    mpz_set_ui(rop, 0);
  y = strdup(x);
  for (l = strlen(y) - 1; y[l] == '\n'; l--) {
    y[l] = '\0';
    if (l == 0)
      break;
  }
  rv = mpz_set_str(rop, y, base);
  free(y);
  return rv;
}

#ifdef NEED_MPZ_MUL_SI
/********************************************************/
void mpz_mul_si(mpz_t x, mpz_t y, long int z)
/********************************************************/
{
  if (z < 0) {
    mpz_mul_ui(x, y, (unsigned long) (-z));
    mpz_neg(x, x);
  } else {
    mpz_mul_ui(x, y, (unsigned long) z);
  }
}
#endif
