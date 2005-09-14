/**************************************************************/
/* getprimes.c                                                */
/* Copyright 2004, Chris Monico.                              */
/**************************************************************/
/*  This file is part of GGNFS.
*
*   GGNFS is free software; you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation; either version 2 of the License, or
*   (at your option) any later version.
*
*   GGNFS is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with GGNFS; if not, write to the Free Software
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ggnfs.h"

#define BIT(_i) ((0x00000001)<<(_i))


/* Here's what we provide: */
u32 *getPList(u32 *numPrimes);
u32 getNextPrime(u32 n);
u32 getPrevPrime(u32 n);
u32 pSieve(u32 *res, u32 size, u32 a, u32 b);
u32 getMaxP(u32 a, u32 b);
/***************************/

/********************************************************/
u32 *getPList(u32 *numPrimes)
/********************************************************/
/* Input: An integer 'numPrimes'.                       */
/* Output: 'numPrimes' is the actual number of primes   */
/*          generated.                                  */
/* Return value: A pointer to the primes generated. Be  */
/*          sure to deallocate when they are no longer  */
/*          needed!                                     */
/********************************************************/
{ u32 a, b, *primes=NULL, *tmp=NULL;
  u32 numP, thisP, i, tmpSize;
    
  if (!(primes = (u32 *)malloc(*numPrimes*sizeof(u32)))) {
    fprintf(stderr, "Memory allocation error in getPList()!\n");
    return NULL;
  }
  a=1; b=100000;
  tmpSize = getMaxP(a, b);
  tmp = (u32 *)malloc(tmpSize*sizeof(u32));
  numP = 0;
  while (numP < *numPrimes) {
    thisP = pSieve(tmp, tmpSize, a, b);
    for (i=0 ; (numP < *numPrimes) && (i<thisP); i++, numP++) {
      primes[numP] = tmp[i];
    }
    a=b; /* Always, 'b' is composite, so no overlap. */
    b += 100000;
  }
  free(tmp);
  *numPrimes = numP;
  return primes;
}


/*******************************************************************/
u32 getNextPrime(u32 n)
/*******************************************************************/
/* Get the smallest prime greater than 'n'.                        */
/* Crude, but whatever.                                            */
/*******************************************************************/
{ u32 tmp[125], a, b, res;

  a = n+1;
  b = a+500;
  res = pSieve(tmp, 125, a, b);
  if (res < 1) {
    fprintf(stderr, "getNextPrime(): failure!\n");
    exit(-1);
  }
  return tmp[0];
}
  
/*******************************************************************/
u32 getPrevPrime(u32 n)
/*******************************************************************/
/* Get the largest prime less than 'n'.                            */
/* This could be much better and smarter, but it doesn't matter.   */
/*******************************************************************/
{ u32 tmp[700], a, b, res;

  if (n<= 2)
    return 0;
  a = n-1001;
  b = n-1;
  if (a < 1)
    a = 1;
  res = pSieve(tmp, 700, a, b);
  if (res < 1) {
    fprintf(stderr, "getPrevPrime(): failure!\n");
    exit(-1);
  }
  return tmp[res-1];
}

u32 getMaxP(u32 a, u32 b);
static unsigned char pDiff[6542];
static int pDiffInit=0;
static u32 numDiffs=0;
/**********************************************/
const u32 _smallP[] = {
  2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53, 
  59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131, 
  137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,
  227,229,233,239,241,251}; 
void makeDiffs()
/***********************************************************/
/* Make a list of the difference between successive primes */
/* between 2 and 2^16.                                     */
/***********************************************************/
{ unsigned char *S;
  int i, entry;
  u32 index;
  u32  p, q;

  if (pDiffInit) return;
  S = (unsigned char *)malloc((32768)*sizeof(unsigned char));
  memset(S, 0x00, 32768);
  for (i=1; i<54; i++) {
    p = _smallP[i];
    entry = (3*p-1)/2; /* Don't sieve out p. No need to sieve out 2p. */
    while (entry < 32768) {
      S[entry] = 0x01;
      entry += p;
    }
  }
  q=2; index=0;
  for (i=1; i<32768; i++) {
    if (S[i]==0) {
      p = 2*i+1;
      pDiff[index++] = p-q;
      q=p;
    }
  }
  numDiffs = index;
  if (numDiffs != 6541)
    printf("Error: makeDiffs() is broken! (numDiffs = %d)\n", numDiffs);
  free(S);
  pDiffInit=1;
}

/**********************************************/
u32 pSieve(u32 *p, u32 Psize, u32 a, u32 b)
/**********************************************/
/* Find all primes in the range [a, b], where */
/* 1 < a < b.                                 */
/**********************************************/
{ u32 k0, k1, s, q, B, a_, b_, S;
  u32 offset, numP=0, i, j;
  int  qNum;
  unsigned char *E;

  makeDiffs();
  B = (u32)sqrt((double)b);
  a_ = a +  (1-(a&0x01));
  b_ = b -  (1-(b&0x01));
  q = 2;
  if (a_ <= B) {
    i=0;
    while ((q <= b_) && (i < numDiffs)) {
      if (q >= a_)
        p[numP++] = a_ = q;
      q += pDiff[i++];
    }
    if (q&0x01) a_ += 2;
    else a_++;
  }
  if (q > b_) return numP;
  k0 = (a_ -1)>>1;
  k1 = (b_ -1)>>1;
  S = k1-k0+1;
  s = S/8 + 1;
  if (!(E = (unsigned char *)malloc(s*sizeof(char)))) {
    fprintf(stderr, "esieve32() memory allocation error!\n");
    return 0;
  }
  memset(E, 0x00, s); /* Initialize the array to zero. */
  q=3;
  qNum=1; /* Indexing from zero. */
  while (q <= B) {
    offset = (q - ((k0 + ((q+1)>>1))%q))%q;
    while (offset <S) {
      E[offset>>3] |= BIT(offset&0x07);
      offset += q;
    }
    q += pDiff[qNum++];
  }
  /* Now pull out the primes. */
  for (i=0; i<s; i++) {
    if (E[i]^0xFF) {
      for (j=0; j<8; j++) {
        if ((E[i]&BIT(j))==0) {
          p[numP++] = (2*(k0 + 8*i + j)) + 1;
	}
      }
    }
  }
  while (p[numP-1] > b)
    numP--;
  free(E);
  return numP;
}

/*********************************************/
u32 getMaxP(u32 a, u32 b)
/*********************************************/
/* Estimate an upper bound for the number of */
/* primes in the interval [a,b].             */
/*********************************************/
{ double x;

  if (a > 30000)
    x = 1.11*(double)(b-a)/log((double)a);
  else 
    x = 1.11*(double)(b)/log((double)b);
  return (u32)x + 25;
}


/*********************************************/
u32 approxPi_x(u32 b)
/*********************************************/
{
  return (u32)(1.075*b/log((double)b));
}
