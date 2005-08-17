/* psp.c
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
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>
#include "lasieve-asm.h"

typedef unsigned long long ull;

extern ulong montgomery_inv_n;
extern ulong *montgomery_modulo_n;
extern ulong montgomery_modulo_R2[NMAX_ULONGS];
extern ulong montgomery_ulongs;

extern void (*asm_mulmod) (ulong *, ulong *, ulong *);
extern void (*asm_squmod) (ulong *, ulong *);
extern void (*asm_add2) (ulong *, ulong *);
extern void (*asm_diff) (ulong *, ulong *, ulong *);

#if 0
int asm_cmp64(ulong * a, ulong * b)
{
  if (a[0] != b[0])
    return 1;
  if (a[1] != b[1])
    return 1;
  return 0;
}
#endif

/*******************************************************/
int psp64()
/*******************************************************/
{ ulong x[2], ex[2], one[2], s;
  long e, i, b=0, v;

  if (!(montgomery_modulo_n[0] & 1))
    return 0;

  ex[0] = montgomery_modulo_n[0] - 1;
  ex[1] = montgomery_modulo_n[1];
  if (!(ex[0] | ex[1]))
    return 0;
  e = 0;
  while (!(ex[0] & 1)) {
    ex[0] = (ex[0] >> 1) | (ex[1] << 31);
    ex[1] >>= 1;
    e++;
  }
  x[0] = 1;
  x[1] = 0;
  asm_mulm64(x, montgomery_modulo_R2, x);
  one[0] = x[0];
  one[1] = x[1];
  if (ex[1]) {
    for (i = 0; i < 32; i++)
      if (ex[1] & (1 << i))
        b = i;
    b += 32;
  } else {
    for (i = 0; i < 32; i++)
      if (ex[0] & (1 << i))
        b = i;
  }
  s = 1;
  for (i = b - 1; i >= 0; i--) {
    if (s >= 32)
      break;
    s <<= 1;
    if (i >= 32)
      v = ex[1] & (1 << (i - 32));
    else
      v = ex[0] & (1 << i);
    if (v)
      s++;
  }
  x[0] = 0;
  x[1] = 0;
  x[s >> 5] = 1 << (s & 31);
  asm_mulm64(x, montgomery_modulo_R2, x);
  for (; i >= 0; i--) {
    asm_sqm64(x, x);
    if (i >= 32)
      v = ex[1] & (1 << (i - 32));
    else
      v = ex[0] & (1 << i);
    if (v)
      asm_add64(x, x);
  }                             /* x=2^ex */
  if (asm_cmp64(x, one) == 0)
    return 1;
  asm_diff64(one, montgomery_modulo_n, one);
  if (asm_cmp64(x, one) == 0)
    return 1;
  for (i = 0; i < e - 1; i++) {
    asm_sqm64(x, x);
    if (asm_cmp64(x, one) == 0)
      return 1;
  }
  return 0;
}

#if 0
int asm_cmp(ulong * a, ulong * b)
{
  long i;

  for (i = 0; i < montgomery_ulongs; i++)
    if (a[i] != b[i])
      return 1;
  return 0;
}
#endif

/*******************************************************/
int psp(mpz_t n)
/*******************************************************/
{ ulong x[NMAX_ULONGS], ex[NMAX_ULONGS], one[NMAX_ULONGS], s;
  long e, i, b, v;

  if (!set_montgomery_multiplication(n))
    return -1;
  if (montgomery_ulongs == 2)
    return psp64();

  if (!(montgomery_modulo_n[0] & 1))
    return 0;                   /* number is even */
  ex[0] = montgomery_modulo_n[0] - 1;
  for (i = 1; i < montgomery_ulongs; i++)
    ex[i] = montgomery_modulo_n[i];
  for (i = 0; i < montgomery_ulongs; i++)
    if (ex[i])
      break;
  if (i >= montgomery_ulongs)
    return 0;                   /* number is 1 */

  e = 0;
  while (!(ex[0] & 1)) {
    for (i = 0; i < montgomery_ulongs - 1; i++)
      ex[i] = (ex[i] >> 1) | (ex[i + 1] << 31);
    ex[montgomery_ulongs - 1] >>= 1;
    e++;
  }
  one[0] = 1;
  for (i = 1; i < montgomery_ulongs; i++)
    one[i] = 0;
  asm_mulmod(one, montgomery_modulo_R2, one);

  for (i = montgomery_ulongs - 1; i >= 0; i--)
    if (ex[i])
      break;
  if (i < 0)
    complain("psp\n");
  b = 32 * i + 31;
  v = ex[i];
  while (v >= 0) {
    b--;
    v <<= 1;
  }

  for (i = 0; i < montgomery_ulongs; i++)
    x[i] = 0;
  s = 1;
  for (i = b - 1; i >= 0; i--) {
    if (s >= 16 * montgomery_ulongs)
      break;
    s <<= 1;
    v = ex[i >> 5] & (1 << (i & 31));
    if (v)
      s++;
  }
  x[s >> 5] = 1 << (s & 31);
  asm_mulmod(x, montgomery_modulo_R2, x);

  if (i >= 0) {
    v = (ex[i >> 5] << (31 - (i & 31)));
    while (1) {
#if 1
      asm_squmod(x, x);
#else
      asm_mulmod(x, x, x);
#endif
      if (v < 0)
        asm_add2(x, x);
      if (i & 31) {
        i--;
        v <<= 1;
      } else {
        if (!i)
          break;
        i--;
        v = ex[i >> 5];
      }
    }
  }
  if (asm_cmp(x, one) == 0)
    return 1;
  asm_diff(one, montgomery_modulo_n, one);
  if (asm_cmp(x, one) == 0)
    return 1;
  for (i = 0; i < e - 1; i++) {
#if 1
    asm_squmod(x, x);
#else
    asm_mulmod(x, x, x);
#endif
    if (asm_cmp(x, one) == 0)
      return 1;
  }
  return 0;
}
