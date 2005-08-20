/* real-poly-aux.c
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

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/types.h>

#include "lasieve.h"

static double *aux;
static size_t aux_alloc = 0;

/* Action of GL(2,Z) on homogeneous polys. */

/************************************************************/
void tpol(double *rop, double *op, i32_t deg, i32_t x0, i32_t x1,
          i32_t y0, i32_t y1)
/************************************************************/
{
  i32_t d;

  if (deg == UINT_MAX)
    complain("Degree too large\n");
  if (deg == 0) {
    rop[0] = op[0];
    return;
  }
  if (aux_alloc < deg + 1) {
    if (aux_alloc > 0)
      free(aux);
    aux_alloc = deg + 1;
    aux = xmalloc(aux_alloc * sizeof(*aux));
  }
  rop[0] = op[0] * y1;
  rop[1] = op[0] * y0;
  aux[0] = x1;
  aux[1] = x0;
  for (d = 1;; d++) {
    i32_t i;

    for (i = 0; i <= d; i++)
      rop[i] += op[d] * aux[i];
    if (d == deg)
      return;
    /* Multiply rop by (y0*T+y1), where T is the free variable. */
    rop[d + 1] = y0 * rop[d];
    for (i = d; i > 0; i--)
      rop[i] = y1 * rop[i] + y0 * rop[i - 1];
    rop[0] *= y1;
    /* Multiply aux by (x0*T+x1). */
    aux[d + 1] = x0 * aux[d];
    for (i = d; i > 0; i--)
      aux[i] = x1 * aux[i] + x0 * aux[i - 1];
    aux[0] *= x1;
  }
}

/************************************************************/
double rpol_eval(double *p, i32_t d, double x, double y)
/************************************************************/
{
  double v, z;
  i32_t i;

  for (i = 0, v = p[0], z = x; i < d;) {
    v = y * v + p[++i] * z;
    z *= x;
  }
  return v;
}

/************************************************************/
double rpol_lb(double *pol, i32_t poldeg, double a, double b)
/************************************************************/
{
  i32_t i;
  double m1, m2, s1, s2;

  if (poldeg == UINT_MAX)
    complain("Degree too large\n");
  if (aux_alloc < poldeg + 1) {
    if (aux_alloc > 0)
      free(aux);
    aux_alloc = poldeg + 1;
    aux = xmalloc(aux_alloc * sizeof(*aux));
  }

  if (a < 0) {
    double rv1, rv2;
    i32_t s;

    for (i = 0, s = 1; i <= poldeg; i++, s = -s)
      aux[i] = pol[i] * s;
    if (b < 0)
      return rpol_lb(aux, poldeg, -b, -a);
    rv1 = rpol_lb(aux, poldeg, 0, -a);
    if ((rv2 = rpol_lb(pol, poldeg, 0, b)) < rv1)
      return rv2;
    return rv1;
  }

  if (pol[poldeg] < 0) {
    s1 = -pol[poldeg];
    s2 = -pol[poldeg];
    m1 = 0;
    m2 = 0;
  } else {
    m1 = pol[poldeg];
    m2 = pol[poldeg];
    s1 = 0;
    s2 = 0;
  }

  for (i = 1; i <= poldeg; i++) {
    if (pol[poldeg - i] < 0) {
      s1 = s1 * a - pol[poldeg - i];
      s2 = s2 * b - pol[poldeg - i];
      m1 *= a;
      m2 *= b;
    } else {
      s1 *= a;
      s2 *= b;
      m1 = m1 * a + pol[poldeg - i];
      m2 = m2 * b + pol[poldeg - i];
    }
  }
  if (s2 < m1)
    return m1 - s2;
  if (s1 > m2)
    return s1 - m2;
  return -1;
}

#if 0
/* Use this to debug tpol and rpol_eval. */
void Usage()
{
  exit(1);
}

main()
{
  i32_t d = 3, i;
  double p[4], q[4], a, b;
  i32_t x0, x1, y0, y1;

#if 1
  p[0] = 43;
  p[1] = 15;
  p[2] = 18;
  p[3] = -100;
  x0 = 5;
  x1 = -3;
  y0 = 3;
  y1 = -2;
#else
  p[0] = 0;
  p[1] = 0;
  p[2] = 0;
  p[3] = 1;
  x0 = 1;
  x1 = 1;
  y0 = 1;
  y1 = 0;
#endif
  tpol(q, p, d, x0, x1, y0, y1);
  scanf("%lf %lf", &a, &b);
  printf("%f == %f ?\n", rpol_eval(p, d, x0 * a + x1 * b, y0 * a + y1 * b),
         rpol_eval(q, d, a, b));
  for (i = 0; i <= d; i++)
    printf(i == 0 ? "%f" : " %f", q[i]);
  printf("\n");
  exit(0);
}
#endif
