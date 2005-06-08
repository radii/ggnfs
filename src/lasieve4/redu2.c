/* redu2.c
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

#include <math.h>
#ifdef __ppc__
#include "ppc32/siever-config.h"
#else
#include "asm/lasieve-asm.h"
#include "lasieve.h"
#endif
i32_t n_iter = 0;

/**************************************************************/
void reduce2(i32_t * a0_ptr, i32_t * b0_ptr, i32_t * a1_ptr, i32_t * b1_ptr,
             i32_t a0, i32_t b0, i32_t a1, i32_t b1, float sigma)
/**************************************************************/
{
  float a0sq, a1sq, s;

  a0sq = ((float) a0) * a0;
  a0sq += sigma * ((float) b0) * b0;
  a1sq = ((float) a1) * a1;
  a1sq += sigma * ((float) b1) * b1;

  for (;;) {
    s = ((float) a0) * a1;
    s += sigma * ((float) b0) * b1;

    n_iter++;
    if (a0sq < a1sq) {
      i32_t k;

      k = rint(s / a0sq);
      if (k == 0)
        break;
      a1 -= k * a0;
      b1 -= k * b0;
      a1sq = ((float) a1) * a1;
      a1sq += sigma * ((float) b1) * b1;
    } else {
      i32_t k;

      k = rint(s / a1sq);
      if (k == 0)
        break;
      a0 -= k * a1;
      b0 -= k * b1;
      a0sq = ((float) a0) * a0;
      a0sq += sigma * ((float) b0) * b0;
    }
  }

  if (b0 < 0) {
    b0 = -b0;
    a0 = -a0;
  }
  if (b1 < 0) {
    b1 = -b1;
    a1 = -a1;
  }
  if (a0sq > a1sq) {
    *a0_ptr = a0;
    *b0_ptr = b0;
    *a1_ptr = a1;
    *b1_ptr = b1;
  } else {
    *a0_ptr = a1;
    *b0_ptr = b1;
    *a1_ptr = a0;
    *b1_ptr = b0;
  }
}
