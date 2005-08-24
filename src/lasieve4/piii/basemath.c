/* basemath.c
  Written by T. Kleinjung. 
  6/13/04: Hacked up for inclusion in GGNFS by Chris Monico. 

  Copyright (C) 2001 Jens Franke, T. Kleinjung.
  This file is part of gnfs4linux, distributed under the terms of the 
  GNU General Public Licence and WITHOUT ANY WARRANTY.
 
  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "siever-config.h"

extern uint32_t montgomery_inv_n;
extern uint32_t *montgomery_modulo_n;
extern uint32_t montgomery_modulo_R2[NMAX_ULONGS];
extern uint32_t montgomery_modulo_R4[NMAX_ULONGS];
extern uint32_t montgomery_ulongs;

extern void (*asm_mulmod) (uint32_t *, uint32_t *, uint32_t *);
extern void (*asm_add2) (uint32_t *, uint32_t *);
extern void (*asm_diff) (uint32_t *, uint32_t *, uint32_t *);
extern void (*asm_zero) (uint32_t *);
extern void (*asm_copy) (uint32_t *, uint32_t *);
extern void (*asm_add2_ui) (uint32_t *, uint32_t);
extern void (*asm_sub) (uint32_t *, uint32_t *, uint32_t *);
extern void (*asm_sub_n) (uint32_t *, uint32_t *);
extern void (*asm_half) (uint32_t *);

/**************************************************/
int asm_cmp64(uint32_t * a, uint32_t * b)
/**************************************************/
{
  if (a[0] != b[0])
    return 1;
  if (a[1] != b[1])
    return 1;
  return 0;
}

/**************************************************/
int asm_cmp(uint32_t * a, uint32_t * b)
/**************************************************/
{ long i;

  for (i = 0; i < montgomery_ulongs; i++)
    if (a[i] != b[i])
      return 1;
  return 0;
}

/**************************************************/
void gcd64(uint32_t * gcd, uint32_t * a, uint32_t * b)
/**************************************************/
{ uint32_t r[2], bb[2], aa[2];

  if (!(a[0] | a[1])) {
    gcd[0] = b[0];
    gcd[1] = b[1];
    return;
  }
  if (!(b[0] | b[1])) {
    gcd[0] = a[0];
    gcd[1] = a[1];
    return;
  }
  bb[0] = b[0];
  bb[1] = b[1];
  aa[0] = a[0];
  aa[1] = a[1];
  while (!(bb[0] & 1)) {
    bb[0] = (bb[0] >> 1) | (bb[1] << 31);
    bb[1] >>= 1;
  }
  while (!(aa[0] & 1)) {
    aa[0] = (aa[0] >> 1) | (aa[1] << 31);
    aa[1] >>= 1;
  }
  while (1) {
    asm_diff64(aa, bb, r);
    if (!(r[0] | r[1]))
      break;
    while (!(r[0] & 1)) {
      r[0] = (r[0] >> 1) | (r[1] << 31);
      r[1] >>= 1;
    }
    if ((aa[1] > bb[1]) || ((aa[1] == bb[1]) && (aa[0] > bb[0]))) {
      aa[1] = r[1];
      aa[0] = r[0];
    } else {
      bb[1] = r[1];
      bb[0] = r[0];
    }
  }
  gcd[0] = aa[0];
  gcd[1] = aa[1];
}

/**************************************************
void gcd(uint32_t * gcd, uint32_t * a, uint32_t * b)
**************************************************
{ uint32_t r[NMAX_ULONGS], bb[NMAX_ULONGS], aa[NMAX_ULONGS];
  long i;

  for (i = 0; i < montgomery_ulongs; i++)
    if (a[i])
      break;
  if (i >= montgomery_ulongs) {
    asm_copy(gcd, b);
    return;
  }
  for (i = 0; i < montgomery_ulongs; i++)
    if (b[i])
      break;
  if (i >= montgomery_ulongs) {
    asm_copy(gcd, a);
    return;
  }

  asm_copy(bb, b);
  asm_copy(aa, a);
  while (!(bb[0] & 1)) {
    for (i = 0; i < montgomery_ulongs - 1; i++)
      bb[i] = (bb[i] >> 1) | (bb[i + 1] << 31);
    bb[montgomery_ulongs - 1] >>= 1;
  }
  while (!(aa[0] & 1)) {
    for (i = 0; i < montgomery_ulongs - 1; i++)
      aa[i] = (aa[i] >> 1) | (aa[i + 1] << 31);
    aa[montgomery_ulongs - 1] >>= 1;
  }
  while (1) {
    asm_diff(r, aa, bb);
    for (i = 0; i < montgomery_ulongs; i++)
      if (r[i])
        break;
    if (i >= montgomery_ulongs)
      break;
    while (!(r[0] & 1)) {
      for (i = 0; i < montgomery_ulongs - 1; i++)
        r[i] = (r[i] >> 1) | (r[i + 1] << 31);
      r[montgomery_ulongs - 1] >>= 1;
    }
    for (i = montgomery_ulongs - 1; i >= 0; i--)
      if (aa[i] != bb[i])
        break;
    if ((i >= 0) && (aa[i] > bb[i]))
      asm_copy(aa, r);
    else
      asm_copy(bb, r);
  }
  asm_copy(gcd, aa);
}*/

/**************************************************/
void asm_half_old(uint32_t * a)
/**************************************************/
{ uint32_t c, n_half[NMAX_ULONGS];
  long i;

  for (i = 0; i < montgomery_ulongs - 1; i++)
    n_half[i] =
      (montgomery_modulo_n[i] >> 1) | (montgomery_modulo_n[i + 1] << 31);
  n_half[montgomery_ulongs - 1] =
    montgomery_modulo_n[montgomery_ulongs - 1] >> 1;
  asm_add2_ui(n_half, 1);       /* (N+1)/2 */

  c = a[0] & 1;
  for (i = 0; i < montgomery_ulongs - 1; i++)
    a[i] = (a[i] >> 1) | (a[i + 1] << 31);
  a[montgomery_ulongs - 1] >>= 1;
  if (c)
    asm_add2(a, n_half);
}

/**************************************************/
int asm_invert(uint32_t * res, uint32_t * b)
/**************************************************/ 
{ long i, f1, len;
  uint32_t t1[NMAX_ULONGS], t2[NMAX_ULONGS];
  uint32_t v1[NMAX_ULONGS], v2[NMAX_ULONGS];
  uint32_t n_half[NMAX_ULONGS];

  for (i = 0; i < montgomery_ulongs; i++)
    if (b[i])
      break;
  if (i >= montgomery_ulongs)
    return 0;
  if (b[0] & 1) {
    asm_copy(t1, b);
    f1 = 0;
  } else {
    asm_sub(t1, montgomery_modulo_n, b);
    f1 = 1;
  }
  asm_zero(t2);
  t2[0] = 1;
  asm_copy(v1, montgomery_modulo_n);
  asm_zero(v2);
#if 0
  for (i = 0; i < montgomery_ulongs - 1; i++)
    n_half[i] =
      (montgomery_modulo_n[i] >> 1) | (montgomery_modulo_n[i + 1] << 31);
  n_half[montgomery_ulongs - 1] =
    montgomery_modulo_n[montgomery_ulongs - 1] >> 1;
  asm_add2_ui(n_half, 1);       /* (N+1)/2 */
#endif
  len = montgomery_ulongs - 1;
  while (1) {
    if (!(t1[len] | v1[len]))
      len--;
    for (i = len; i >= 0; i--)
      if (t1[i] != v1[i])
        break;
    if (i < 0)
      break;
    if (t1[i] > v1[i]) {        /* t1>v1 */
      asm_sub_n(t1, v1);        /* t1 even */
      asm_sub(t2, t2, v2);
      do {
#if 1
        for (i = 0; i < len; i++)
          t1[i] = (t1[i] >> 1) | (t1[i + 1] << 31);
        t1[len] >>= 1;
#else
        asm_half(t1);
#endif
        asm_half(t2);
      } while (!(t1[0] & 1));
    } else {                    /* v1>t1 */
      asm_sub_n(v1, t1);        /* t1 even */
      asm_sub(v2, v2, t2);
      do {
        for (i = 0; i < len; i++)
          v1[i] = (v1[i] >> 1) | (v1[i + 1] << 31);
        v1[len] >>= 1;
        asm_half(v2);
      } while (!(v1[0] & 1));
    }
  }
  if (t1[0] != 1)
    return 0;
  for (i = 1; i < montgomery_ulongs; i++)
    if (t1[i])
      return 0;
  if (f1)
    asm_sub(res, montgomery_modulo_n, t2);
  else
    asm_copy(res, t2);
  asm_mulmod(res, montgomery_modulo_R4, res);
  return 1;
}
