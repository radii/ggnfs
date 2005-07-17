/* fbgen.c 
   Originally written by F. Bahr. For use in the quality estimator for GNFS
   polynomials, adjustments (sometimes crude adjustments) have been made by
   J. Franke.
   6/13/04: Hacked up and uglied by Chris Monico for use in GGNFS

  Copyright (C) 2001 Friedrich Bahr

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/
/* CJM: Rather than bothering to try and optimize my root finding code, I'm
   going to just use this one for FB generation since it already appears
   to be pretty well optimized.
*/

#include <limits.h>
#include <sys/types.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif
#include <math.h>
#include <gmp.h>
#include "ggnfs.h"

#define HAVE_CMOV



#define U32_MAX 0xffffffff
#define I32_MAX INT_MAX
#define ULL_NO_UL



/* 32bit.h */
#ifdef __ppc__
#include "lasieve4/ppc32/siever-config.h"
#else
#define modinv32(x) asm_modinv32(x)

#if defined (_MSC_VER) && !defined(__MINGW32__)

#define inline __inline
volatile extern u32 modulo32;
u32 asm_modinv32(u32 x);

static inline u32 modsq32(u32 x)
{	u32 res;

	__asm
	{	mov		eax,[x]
		mul		[x]
		div		[modulo32]
		mov		res,edx
	}
	return res;
}

static inline u32 modmul32(u32 x, u32 y)
{	u32 res;
	
	__asm
	{	mov		eax,[x]
		mul		[y]
		div		[modulo32]
		mov		res,edx
	}
	return res;
}

#ifdef HAVE_CMOV

static inline u32 modadd32(u32 x, u32 y)
{	u32 res;
	
	__asm
	{
		xor		edx,edx
		mov		ecx,[y]
		add		ecx,[x]
		cmovc	edx,[modulo32]
	    cmp		ecx,[modulo32]
		cmovae	edx,[modulo32]
		sub		ecx,edx
		mov		[res],ecx
	}
	return res;
}

static inline u32 modsub32(u32 subtrahend, u32 minuend)
{	u32 res;
  
	__asm
	{	xor		edx,edx
		mov		ecx,[subtrahend]
		sub		ecx,[minuend]
		cmovb	edx,[modulo32]
		add		ecx,edx
		mov		[res],ecx
	}
	return res;
}

#else

static inline u32 modadd32(u32 x,u32 y)
{	u32 res;
	
	__asm
	{	
		mov		ecx,[x]
		add		ecx,[y]
		jc		l1f
		cmp		ecx,[modulo32]
		jb		l2f
l1f:	sub		ecx,[modulo32]
l2f:	mov		[res],ecx
	}
	return res;
}

static inline u32 modsub32(u32 subtrahend,u32 minuend)
{	u32 res;

	__asm
	{
		mov		ecx,[subtrahend]
		sub		ecx,[minuend]
		jae		l1f
		addl	ecx,[modulo32]
l1f:	mov		[res],ecx
	}
	return res;
}

#endif

#else

volatile extern u32 modulo32 asm("modulo32");
u32 asm_modinv32(u32 x) asm("asm_modinv32");

static inline u32 modsq32(u32 x)
{ u32 res,clobber;
  __asm__ volatile ("mull %%eax\n"
	   "divl modulo32" : "=d" (res), "=a" (clobber) : "a" (x) : "cc" );
  return res;
}

static inline u32 modmul32(u32 x,u32 y)
{ u32 res,clobber;
  __asm__ volatile ("mull %%ecx\n"
	   "divl modulo32" : "=d" (res), "=a" (clobber) : "a" (x), "c" (y) :
	   "cc");
  return res;
}

static inline u32 modadd32(u32 x,u32 y)
{ u32 res;
#ifdef HAVE_CMOV
  __asm__ volatile ("xorl %%edx,%%edx\n"
	   "addl %%eax,%%ecx\n"
	   "cmovc modulo32,%%edx\n"
	   "cmpl modulo32,%%ecx\n"
	   "cmovae modulo32,%%edx\n"
	   "subl %%edx,%%ecx\n"
	   "2:\n" : "=c" (res) : "a" (x), "c" (y) : "%edx", "cc");
#else
  __asm__ volatile ("addl %%eax,%%ecx\n"
	   "jc 1f\n"
	   "cmpl modulo32,%%ecx\n"
	   "jb 2f\n"
	   "1:\n"
	   "subl modulo32,%%ecx\n"
	   "2:\n" : "=c" (res) : "a" (x), "c" (y) : "cc");
#endif
  return res;
}

static inline u32 modsub32(u32 subtrahend,u32 minuend)
{ u32 res;
#ifdef HAVE_CMOV
  __asm__ volatile ("xorl %%edx,%%edx\n"
	   "subl %%eax,%%ecx\n"
	   "cmovbl modulo32,%%edx\n"
	   "addl %%edx,%%ecx\n"
	   "1:" : "=c" (res) : "a" (minuend), "c" (subtrahend) : "%edx", "cc");
#else
  __asm__ volatile ("subl %%eax,%%ecx\n"
	   "jae 1f\n"
	   "addl modulo32,%%ecx\n"
	   "1:" : "=c" (res) : "a" (minuend), "c" (subtrahend) : "cc" );
#endif
  return res;
}

#endif
volatile u32 modulo32;
#endif 
u32 modulo32hbit;
s32 modulo32bit[32], modulo32ebits;
u32 *polr, t, polr_alloc = 0;

/****************************************************/
void *xmalloc(size_t size)
/****************************************************/
{ char *x;

  if (size == 0)
    return NULL;
  if ((x = malloc(size)) == NULL)
    fprintf(stderr, "xmalloc() memory allocation error!\n");
  return x;
}



/*****************************************************/
inline u32 polmodsq32(u32 * P, u32 dP)
/*****************************************************/
{
  u32 i, j, a, m;

  m = P[dP];
  P[dP << 1] = modsq32(m);
  for (i = 0; i < dP; i++) {
    a = modmul32(m, P[i]);
    P[dP + i] = modadd32(a, a);
  }
  for (j = dP; j > 1;) {
    m = P[--j];
    P[j << 1] = modadd32(P[j << 1], modsq32(m));
    a = modmul32(m, P[0]);
    P[j] = modadd32(a, a);
    for (i = 1; i < j; i++) {
      a = modmul32(m, P[i]);
      P[j + i] = modadd32(P[j + i], modadd32(a, a));
    }
  }
  if (dP != 0)
    P[0] = modsq32(P[0]);
  dP = dP << 1;
  return dP;
}

/*****************************************************/
u32 poldivmod32(u32 * P, u32 dP, u32 * D, u32 dD)
/*****************************************************/
{
  u32 i, a, b, j;

  /* D!=0 ! vorausgesetzt */
  while (dP >= dD) {
    a = P[dP], b = D[dD];
    for (i = j = dP - dD; i < dP; i++)
      P[i] = modsub32(modmul32(P[i], b), modmul32(D[i - j], a));
    for (i = 0; i < j; i++)
      P[i] = modmul32(P[i], b);
    for (dP--; P[dP] == 0 && dP > 0; dP--);
  }
  return dP;
}

/*****************************************************/
u32 polcpymod32(u32 * Q, u32 * P, u32 dP)
/*****************************************************/
{
  memcpy(Q, P, (dP + 1) * sizeof(u32));
  return dP;
}

/*****************************************************/
u32 polgcdmod32(u32 * P, u32 dP, u32 * Q, u32 dQ)
/*****************************************************/
{
  u32 *S, d, *R = Q;

  /* hier wird dP > dQ angenommen */
  /* stets ist Q das Ergebnis */

  while (dQ) {
    d = poldivmod32(P, dP, Q, dQ);
    S = P, P = Q, Q = S, dP = dQ, dQ = d;
  }
  if (Q[0] == 0) {
    if (P != R)
      polcpymod32(R, P, dP);
    return dP;
  }
  R[0] = 1;
  return 0;
}

/*****************************************************/
inline u32 poldivnormmod32(u32 * P, u32 dP, u32 * D, u32 dD)
/*****************************************************/
{
  u32 i, a, j;

  /* D 
     normiert */
  while (dP >= dD) {
    a = P[dP];
    for (i = j = dP - dD; i < dP; i++)
      P[i] = modsub32(P[i], modmul32(D[i - j], a));
    for (dP--; P[dP] == 0 && dP > 0; dP--);
  }
  return dP;
}

/*****************************************************/
u32 polnormmod32(u32 * P, u32 dP)
/*****************************************************/
{
  u32 a, i;

  if (P[dP] != 1 && P[dP] != 0) {
    a = modinv32(P[dP]);
    for (i = 0; i < dP; i++)
      P[i] = modmul32(P[i], a);
    P[dP] = 1;
  }
  return dP;
}

/* CAVE Was passiert, falls die Eingabe durch |modulo32| teilbar ist? */
/*****************************************************/
u32 polredmod32(u32 * T, __mpz_struct * A, u32 dA)
/*****************************************************/
{
  u32 i;

  for (i = 0; i <= dA; i++)
    T[i] = mpz_fdiv_ui(&A[i], modulo32);
  while (T[dA] == 0)
    dA--;
  return polnormmod32(T, dA);;
}

/*****************************************************/
u32 polXpotmodPmod32(u32 * Q, u32 * P, u32 dP)
/*****************************************************/
{
  u32 dQ = modulo32bit[modulo32hbit], k = 0;

  /* Q = X hoch (modulo32-1)/2 - 1 modulo P */
  Q += modulo32ebits;
  for (; k < dQ; k++)
    Q[k] = 0;
  Q[dQ] = 1;
  if (dQ >= dP)
    dQ = poldivnormmod32(Q, dQ, P, dP);
  for (k = modulo32hbit; k > 0;) {
    k--;
    dQ = polmodsq32(Q, dQ);
    if (modulo32bit[k]) {       /* Multiplikation mit X */
      Q--;
      dQ++;
      Q[0] = 0;
    }
    if (dQ >= dP)
      dQ = poldivnormmod32(Q, dQ, P, dP);
  }
  return dQ;
}

/*****************************************************/
void polchcomod32(u32 * P, u32 dP)
/*****************************************************/
{
  u32 i = dP, j;

  for (; i > 0;) {
    for (i--, j = i; j < dP; j++)
      P[j] = modadd32(P[j], P[j + 1]);
  }
}

/*****************************************************/
u32 polrootrecmod32(u32 * T, u32 dT, u32 * r, u32 * Q)
/*****************************************************/
{
  u32 *P, dP, dQ, i = 0, b;

  /* Ersetze T(X) durch T(X+1) */
  if (dT > 1) {
    polchcomod32(T, dT);
    dQ = polXpotmodPmod32(Q, T, dT);
    Q[0] = modadd32(Q[0], 1);
    P = Q + dT + 1;
    dP = polcpymod32(P, T, dT);
    dQ = polgcdmod32(P, dP, Q, dQ);
    b = r[0];
    r[0] = modadd32(b, 1);
    if (dQ > 0 && dQ < dT) {
      dQ = polnormmod32(Q, dQ);
      dP = polcpymod32(P, T, dT);
      dP = poldivnormmod32(P, dP, Q, dQ);
      dP = dT - dQ, P += dQ;
      i = polrootrecmod32(Q, dQ, r, P + dP + 1);        /*i=dQ */
      r[i] = modadd32(b, 1);
      i += polrootrecmod32(P, dP, r + i, P + dP + 1);   /*i=dT */
      return i;
    } else
      return polrootrecmod32(T, dT, r, Q);
  } else {
    r[0] = modsub32(r[0], T[0]);
    return 1;
  }
}

/*****************************************************/
void polprintmod32(u32 * P, u32 dP, char *c)
/*****************************************************/
{
  u32 i;

  printf("\n%s %" PRIu32 "\n", c, dP);
  for (i = 0; i <= dP; i++)
    printf("%" PRIu32 " ", P[i]);
  printf("\n");
}

/*****************************************************/
void mod32multab(dT)
/*****************************************************/
{
  u32 i, c;

  modulo32ebits = 0;
  for (i = 0, c = modulo32; (c >>= 1) > dT; i++) {
    if (c & 1) {
      modulo32bit[i] = 1;
      modulo32ebits++;
    } else
      modulo32bit[i] = 0;
  }
  modulo32bit[i] = c;
  modulo32hbit = i;             /* Anzahl der notwendigen Quadratbildungen */
}

/*****************************************************/
int uintcmp(const void *x, const void *y)
/*****************************************************/
{
  const u32 ux = *(const u32 *) x;
  const u32 uy = *(const u32 *) y;

  return ux >= uy ? 1 : -1;
}

/*****************************************************/
u32 polrootmod32(__mpz_struct  * A, u32 dT, u32 * r)
/*****************************************************/
{
  u32 i = 0, j, *P, dP, *Q, dQ, *S, dS, *T = polr;

  dT = polredmod32(T, A, dT);
  if (T[0] == 0) {              /* Nullstelle 0 */
    for (i = 1; T[i] == 0 && i < dT; i++);
    T += i;
    dT -= i;
    r[0] = 0, i = 1;
  }
  if (dT == 0)
    return i;
  mod32multab(dT);
  /* Zerlegung Q = (P,X^(modulo32-1)/2+1); S = (P,X^(modulo32-1)/2+1) */
  P = T + dT + 1;
  dP = polcpymod32(P, T, dT);
  Q = P + dP + 1;
  dQ = polXpotmodPmod32(Q, P, dP);
  S = Q + dP + 1;               /* die ggT-Berechnung kann P ergeben */
  dS = polcpymod32(S, Q, dQ);
  Q[0] = modadd32(Q[0], 1);
  dQ = polgcdmod32(P, dP, Q, dQ);
  dP = polcpymod32(P, T, dT);
  S[0] = modsub32(S[0], 1);
  dS = polgcdmod32(P, dP, S, dS);
  if (dQ > 0) {
    dQ = polnormmod32(Q, dQ);
    r[i] = 0;
    j = polrootrecmod32(Q, dQ, r + i, S + dS + 1);      /*j=dQ */
    if (j != dQ)
      fprintf(stderr, "Falsche Nullstellenanzahl Q, P %" PRIu32 "\n", modulo32);
    i += j;
  }
  if (dS > 0) {
    dS = polnormmod32(S, dS);
    r[i] = 0;
    j = polrootrecmod32(S, dS, r + i, S + dS + 1);      /*j=dS */
    if (j != dS)
      fprintf(stderr, "Falsche Nullstellenanzahl S, P %" PRIu32 "\n", modulo32);
    i += j;
  }
  qsort(r, i, sizeof(u32), uintcmp);
  return i;
}

/*****************************************************/
s32 root_finder(s32 * root_buf, __mpz_struct * A, s32 adeg, s32 p)
/*****************************************************/
{
  u32 res;

  if (polr_alloc < 300 * adeg * adeg) {
    /* CAVE falls adeg==0? */
    if (polr_alloc > 0)
      free(polr);
    polr_alloc = 300 * adeg * adeg;
    polr = (u32 *)xmalloc(polr_alloc * sizeof(*polr));
  }
  if (p == 2) {
    u32 i, pv;

    res = 0;
    if (mpz_fdiv_ui(&A[0], p) == 0)
      root_buf[res++] = 0;
    for (i = 0, pv = 0; i <= adeg; i++)
      if (mpz_fdiv_ui(&A[i], p) != 0)
        pv ^= 1;
    if (pv == 0)
      root_buf[res++] = 1;
    if (mpz_fdiv_ui(&A[adeg], p) == 0)
      root_buf[res++] = 2;
    return res;
  }
  modulo32 = p;
  res = polrootmod32(A, adeg, root_buf);
  /* CAVE: Diese Division erfolgt wohl ein zweites Mal. */
  if (mpz_fdiv_ui(&A[adeg], p) == 0)
    root_buf[res++] = p;
  return res;
}
