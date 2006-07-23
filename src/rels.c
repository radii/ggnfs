/**************************************************************/
/* rels.c                                                     */
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

#ifdef _MSC_VER
#pragma warning (disable: 4996) /* warning C4996: 'function' was declared deprecated */
#endif

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <gmp.h>
#include <limits.h>
#include "ggnfs.h"

#define FB_TRIAL_DIV_FRAC 0.35
#define MAX_TRIAL_DIV_SIZE 50000


#define MAX_FACTORS MAX(MAX_RAT_FACTORS, MAX_ALG_FACTORS)

/* The following is the format for processed relations (in memory and on disk): */
/* [size field][a][b][RFB entries][AFB entries][Sp entries][QCB 1][QCB 2][Large rat. primes][Large alg. primes]
   where we have:
   [size field] = [# RFB entries (high 8 bits)][# AFB (next 8 bits)][# Sp (8 bits)][# LRP (2 bits)][# LAP (2 bits)]
   The 4 least significant bits are currently unused, but reservered for possible use later.
   Then, the RFB/AFB/Sp entries look like [index (32 bits) | exponent (32 bits)].
   Although we don't need 32 bits for the exponent, it does make it much easier and faster
   to pack in and out of memory and to and from disk.
   QCB 1 & 2 are each 32 bits. The large rat. primes are each 32 bits. The large alg. primes
   are each 64 bits, looking like [p | r].
   
   Note: In a file, these fields are preceeded by the number of relations
         in the file. Specifically, the first 4 bytes of the file is a s32
         holding the count.
*/
/* These are defined in ggnfs.h, but cited here for easy access. 
  #define GETNUMRFB(_s) ((_s)>>24)
  #define GETNUMAFB(_s) (((_s)&0x00FF0000)>>16)
  #define GETNUMSPB(_s)  (((_s)&0x0000FF00)>>8)
  #define GETNUMLRP(_s)  (((_s)&0x000000C0)>>6)
  #define GETNUMLAP(_s)  (((_s)&0x00000030)>>4)

  #define SETNUMRFB(_s,_n) (_s = (_s&0x00FFFFFF)^((((s32)(_n)&0x000000FF)<<24)))
  #define SETNUMAFB(_s,_n) (_s = (_s&0xFF00FFFF)^((((s32)(_n)&0x000000FF)<<16)))
  #define SETNUMSPB(_s,_n) (_s = (_s&0xFFFF00FF)^((((s32)(_n)&0x000000FF)<<8)))
  #define SETNUMLRP(_s,_n) (_s = (_s&0xFFFFFF3F)^((((s32)(_n)&0x00000003)<<6)))
  #define SETNUMLAP(_s,_n) (_s = (_s&0xFFFFFFCF)^((((s32)(_n)&0x00000003)<<4)))
  #define S32S_IN_ENTRY(_s) (2*GETNUMRFB(_s)+2*GETNUMAFB(_s)+2*GETNUMSPB(_s)+GETNUMLRP(_s)+2*GETNUMLAP(_s) + 6)
  
*/

#define FHASH_SIZE 4000000
/************************************************************/
s32 lookupRFB(s32 P, nfs_fb_t *FB)
/************************************************************/
/* Lookup the given prime in the RFB. Return it's index, if */
/* found, or -1 otherwise.                                  */
/************************************************************/
{ s32   p, h, *loc;
  u32   k;
  static s32 *rfbHash;
  static int   initialized=0;
  
  if (!initialized) {
    rfbHash = (s32 *)malloc(FHASH_SIZE*sizeof(s32));
    if (rfbHash == NULL) {
      printf("lookupRFB() fatal memory allocation error!\n");
      exit(-1);
    }
    memset(rfbHash, 0x00, FHASH_SIZE*sizeof(s32));
    for (k=1; k<FB->rfb_size; k++) {
      p = FB->rfb[2*k];
      h = NFS_HASH(p, 0, FHASH_SIZE);
      if (rfbHash[h]==0)
        rfbHash[h] = k;
    }
    initialized=1;
  }
  p=P;
  h = NFS_HASH(p, 0, FHASH_SIZE);
  k = rfbHash[h];
  if (FB->rfb[2*k] == p) return k;

  loc = (s32 *)bsearch(&p, FB->rfb, FB->rfb_size, 2*sizeof(s32), cmpS32s);
  if (loc != NULL) {
    assert((loc - FB->rfb) <= INT_MAX);
    k = (s32)(loc - FB->rfb);
    k /= 2; /* Each entry in the array is actually 2 s32's wide. */
    return k;
  }
  return -1;
}

/************************************************************/
s32 lookupAFB(s32 p, s32 r, nfs_fb_t *FB)
/************************************************************/
/* Lookup the given prime in the AFB. Return it's index, if */
/* found, or -1 otherwise.                                  */
/************************************************************/
{ s32 _p, _r, h, *loc;
  u32 k;
  static s32 *afbHash;
  static int   initialized=0;

  if (!initialized) {
    afbHash = (s32 *)malloc(FHASH_SIZE*sizeof(s32));
    if (afbHash == NULL) {
      printf("lookupAFB() fatal memory allocation error!\n");
      exit(-1);
    }
    memset(afbHash, 0x00, FHASH_SIZE*sizeof(s32));
    k=1;
    while (k<FB->afb_size) {
      _p = FB->afb[2*k];
      h = NFS_HASH(_p, 0, FHASH_SIZE);
      if (afbHash[h]==0)
        afbHash[h] = k;
      k++;
    }
    initialized=1;
  }
  _p=p; _r=r;
  h = NFS_HASH(_p, 0, FHASH_SIZE);
  k = afbHash[h];

  if (FB->afb[2*k] != _p) {
    loc = (s32 *)bsearch(&_p, FB->afb, FB->afb_size, 2*sizeof(s32), cmpS32s);
    if (loc != NULL) {
	  assert((loc - FB->afb) <= INT_MAX);
      k = (s32)(loc - FB->afb);
      k /= 2; /* Each entry in the array is actually 2 s32's wide. */
      while ((k>0) && (FB->afb[2*(k-1)]==_p))
        k--;
    }
  }

  while ((k < FB->afb_size) && (FB->afb[2*k]==_p) && (FB->afb[2*k+1] != _r)) {
    k++;
  }
  if ((FB->afb[2*k]==_p) && (FB->afb[2*k+1]==_r))
    return k;
  return -1;
}
  

/****************************************************************************/
int relConvertToData(s32 *data, relation_t *R)
/****************************************************************************/
/* Convert the given relation into the 'data' format, specified at the top  */
/* of this file.                                                            */
/* It is assumed that there is sufficient memory allocated at 'data'.       */
/* Return value: The number of s32s used.                                  */
/****************************************************************************/
{ int i, size, numLRP, numLAP;
  s32 sF;

  size=1;
  sF = 0;
  *( (s64*)&data[size]) = R->a;
  //data[size++] = R->a;
  size += 2;
  data[size++] = R->b;
  for (i=0; i<R->rFSize; i++) {
    data[size++] = R->rFactors[i];
    data[size++] = R->rExps[i];
  }
  for (i=0; i<R->aFSize; i++) {
    data[size++] = R->aFactors[i];
    data[size++] = R->aExps[i];
  }
  for (i=0; i<R->spSize; i++) {
    data[size++] = R->spFactors[i];
    data[size++] = R->spExps[i];
  }
  data[size++] = R->qcbBits[0];  
  data[size++] = R->qcbBits[1];  
  for (numLRP=0; (numLRP<MAX_LARGE_RAT_PRIMES) && (R->p[numLRP] > 1); numLRP++)
    data[size++] = R->p[numLRP];

  for (numLAP=0; (numLAP<MAX_LARGE_ALG_PRIMES) && (R->a_p[numLAP] > 1); numLAP++) {
    data[size++] = R->a_p[numLAP];
    data[size++] = R->a_r[numLAP];
  }
  SETNUMRFB(sF, R->rFSize); SETNUMAFB(sF, R->aFSize); SETNUMSPB(sF, R->spSize);
  SETNUMLRP(sF, numLRP); SETNUMLAP(sF, numLAP);
  data[0] = sF;

  return size;
}

  
/*********************************************************************/
int dataConvertToRel(relation_t *R, s32 *data)
/*********************************************************************/
/* Convert the given data (in the format specified in ggnfs.h) to a  */
/* relation_t.                                                       */
/* Return value: The number of s32s used.                           */
/*********************************************************************/
{ int  size, i, numLAP, numLRP;
  s32 sF;

  sF = data[0];
  size=1;
  R->a = *( (s64*) &data[size] );
  size += 2;

  R->b = data[size++];
  /** Read RFB entries **/
  R->rFSize = (int)GETNUMRFB(sF);
  if (R->rFSize > MAX_RAT_FACTORS) {
    fprintf(stderr, "dataConvertToRel() Error: MAX_RAT_FACTORS too small to hold %d factors!\n",
            R->rFSize);
    return -1;
  }
  for (i=0; i<R->rFSize; i++)  {
    R->rFactors[i] = data[size++];
    R->rExps[i] = data[size++];
  }
  /** Read AFB entries **/
  R->aFSize = (int)GETNUMAFB(sF);
  if (R->aFSize > MAX_ALG_FACTORS) {
    fprintf(stderr, "dataConvertToRel() Error: MAX_ALG_FACTORS too small to hold %d factors!\n",
            R->aFSize);
    fprintf(stderr, "sF field was: %8.8" PRIx32 "\n", sF);
    return -1;
  }
  for (i=0; i<R->aFSize; i++)  {
    R->aFactors[i] = data[size++];
    R->aExps[i] = data[size++];
  }
  /** Read SPB entries **/
  R->spSize = (int)GETNUMSPB(sF);
  if (R->spSize > MAX_SP_FACTORS) {
    fprintf(stderr, "dataConvertToRel() Error: MAX_SP_FACTORS too small to hold %d factors!\n",
            R->spSize);
    return -1;
  }
  for (i=0; i<R->spSize; i++)  {
    R->spFactors[i] = data[size++];
    R->spExps[i] = data[size++];
  }
  /* Read the QCB entries. */ 
  R->qcbBits[0] = data[size++];
  R->qcbBits[1] = data[size++];

  /* Read the large rational primes. */
  numLRP = (int)GETNUMLRP(sF);
  if (numLRP > MAX_LARGE_RAT_PRIMES) {
    fprintf(stderr, "dataConvertToRel() Error: MAX_LARGE_RAT_PRIMES too small to hold %d large primes!\n",
            numLRP);
    return -1;
  }
  for (i=0; i<numLRP; i++)  {
    R->p[i] = data[size++];
  }
  for ( ; i<MAX_LARGE_RAT_PRIMES; i++)
    R->p[i]=1;

  /* Read the large algebraic primes. */
  numLAP = (int)GETNUMLAP(sF);
  if (numLAP > MAX_LARGE_ALG_PRIMES) {
    fprintf(stderr, "dataConvertToRel() Error: MAX_LARGE_ALG_PRIMES too small to hold %d large primes!\n",
            numLAP);
    return -1;
  }
  for (i=0; i<numLAP; i++)  {
    R->a_p[i] = data[size++];
    R->a_r[i] = data[size++];
  }
  for ( ; i<MAX_LARGE_ALG_PRIMES; i++) {
    R->a_p[i] = 1;
    R->a_r[i] = 0; /* Maybe -1 here? */
  }

  return size;
}
 
/***********************************************/
int writeRel(FILE *fp, s32 *relData)
/***********************************************/
/* Write a single relation out to file.        */
/***********************************************/
{ s32 size, sF;

  sF = relData[0];
  size = S32S_IN_ENTRY(sF);
  fwrite(relData, sizeof(s32), size, fp);
  return 0;
}

int readRel(relation_t *R, FILE *fp)
/*****************************************************/
/* Read a relation, in binary format, from fp.       */
/* Return 0 on success, -1 on failure.               */
/* In fact, we only really care about the large prime*/
/* fields, so this could probably be sped up.        */
/*****************************************************/
{ s32 data[1024], dataSize=0, sF;

  readRaw32(&sF, fp);
  if (feof(fp))
    return -1;
  data[0] = sF;
  dataSize = S32S_IN_ENTRY(sF) - 1; /* We already read the first one! */
  if (dataSize > 1024) {
    fprintf(stderr, "readRel() Error: dataSize = %" PRId32 " > 1024! Increase and recompile!\n",
            dataSize);
    return -1;
  }
  fread(&data[1], sizeof(s32), dataSize, fp);
  if (dataConvertToRel(R, data)>0) {
    return 0;
  }
  return -1;
}

/*********************************************************************/
int writeRelList(char *fname, rel_list *L)
/*********************************************************************/
/* Write the rel_list structure out to file, in a specific format.   */
/* It is assumed that we will not exceed the OS-imposed maximum      */
/* file size.                                                        */
/* Note: The function filterRelations() in procrels.c also writes out*/
/* relations to file in this format. So if you change this, change   */
/* that one to match also (the point is that it needs to write them  */
/* out in a different order).                                        */
/*********************************************************************/
{ s32  index, i;
  FILE *fp;
  
  if (!(fp = fopen(fname, "wb"))) {
    fprintf(stderr, "Error writing to file %s!\n", fname);
    return -1;
  }
  /* The first 4 bytes give the # of relations in the file. */
  writeRaw32(fp, &L->numRels);
  index = 0;
  for (i=0; i<L->numRels; i++) {
    writeRel(fp, &L->relData[L->relIndex[i]]);
  }
  fclose(fp);
  return 0;
}

/*****************************************************************************/
int readRelList(rel_list *L, char *fname)
/*****************************************************************************/
/* Read a rel_list structure from file in the same format as writeRelList(). */
/* The caller should already have allocated room in 'L' and set the maxRels, */
/* maxDataSize fields.                                                       */
/*****************************************************************************/
{ s32  size, relsRead, sF, relSize;
  s32  *buf, bufMax;
  size_t bufSize, bufIndex;
  FILE *fp;

  L->numRels = 0; L->relIndex[0]=0;
  if (!(fp = fopen(fname, "rb"))) {
    return -1;
  }
  /* The first 4 bytes of the file always give the
     number of relations in the file.
  */
  readRaw32(&L->numRels, fp);
  
  printf("\n");
  if (L->numRels > L->maxRels) {
    fprintf(stderr, "readRelList() Error: File contains %" PRIu32 " relations vs. maxRels=%" PRIu32 "\n",
            L->numRels, L->maxRels);
    fclose(fp);
    return -1;
  }
  size = relsRead = 0;
  bufMax = 8388608; bufSize=0; bufIndex=0;
  buf = (s32 *)malloc(bufMax*sizeof(s32));
  while ((size < L->maxDataSize) && (relsRead < L->numRels)) {
    /* Do we need to load more data into the buffer? */
    if ((bufIndex + 1024 > bufSize)&&(!feof(fp))) {
      printTmp("readRelList : [ %3.3d%% done]", 100*relsRead/L->numRels);
      memmove(buf, &buf[bufIndex], (bufSize-bufIndex)*sizeof(s32));
      bufSize -= bufIndex;
      bufIndex=0;
      bufSize += fread(&buf[bufSize], sizeof(s32), bufMax-bufSize, fp);
    }

    L->relIndex[relsRead] = size;
    /* Read the relation size field. */
    sF = buf[bufIndex++];
    L->relData[size++] = sF;
    relSize = S32S_IN_ENTRY(sF) - 1; /* We already read the first field! */
    memcpy(&L->relData[size], &buf[bufIndex], relSize*sizeof(s32));
    bufIndex += relSize;
    size += relSize;
    relsRead++;
  }
  fclose(fp);
  free(buf);
  printf("readRelList : [   100%% done]\n");

  if ((relsRead < L->numRels)&&(size < L->maxDataSize)) {
    L->numRels = relsRead;
    fprintf(stderr, "readRelList() Error: File %s appears corrupted!\n", fname);
    return -1;
  }
  if (size >= L->maxDataSize) {
    fprintf(stderr, "  (read %" PRId32 " of %" PRId32 " relations).\n", relsRead, L->numRels);
    L->numRels = relsRead;
    fprintf(stderr, "readRelList() Error: L->relData is not large enough to handle %s!\n", fname);
    return -1;
  }
  L->relIndex[relsRead] = size;
  return 0;
}

/***********************************************************************/
int factRel(relation_t *R, nf_t *N)
/***********************************************************************/
/* Input: A relation (R->a, R->b), and the factor base 'FB'.           */
/* Output: The relation_t fields: (isNeg), rFactors, rExps, rFSize,    */
/*         aFactors, aExps, aFSize, p1, p2, a1_p, a1_r, a2_p, a2_r     */
/*         will be filled.                                             */
/* Return value: 0 if (R->a, R-b) is a valid relation.                 */
/*        Otherwise, some error occurred (not smooth maybe?).          */
/***********************************************************************/
{ static int  initialized=0, qSize=0;
  static s32 afbTrialSize, rfbTrialSize, afbSize, rfbSize;
  static mpz_t     temp1, temp2, norm;
  static mpz_poly  delta;
  static mpz_mat_t deltaHNF;
  s32   i, factors[MAX_FACTORS+1], b;
  s64   a;
  s32   *loc, locIndex, r;
  s32   pFacts[10*MAX_FACTORS], p;
  int    numpFacts;
  char   exponents[MAX_FACTORS+1];
  int    rSize=0, aSize=0, spSize=0;
  int    e, numLarge, numSp=N->numSPrimes;
  nfs_fb_t *FB = N->FB;

  if (!(initialized)) {
    mpz_init(temp1); mpz_init(temp2);
    mpz_init(norm);  mpz_poly_init(delta);
    mpz_mat_init(&deltaHNF);

    qSize = FB->qcb_size/32;
    if (FB->qcb_size % 32)
      qSize++;
    if (qSize > 2) {
      qSize = 2;
      /* Let's take some liberties: */
      FB->qcb_size = 62;
    }
    rfbSize = FB->rfb_size;
    afbSize = FB->afb_size;
    rfbTrialSize = MIN((s32)(FB_TRIAL_DIV_FRAC*rfbSize), MAX_TRIAL_DIV_SIZE);
    afbTrialSize = MIN((s32)(FB_TRIAL_DIV_FRAC*afbSize), MAX_TRIAL_DIV_SIZE);
    initialized=1;
  }
  /********** Get the RFB part. **********/  
  /* Do temp1 <-- a - bm  */
  a = R->a; b = R->b;
  mpz_mul_si64(temp2,FB->y1,a);
  mpz_mul_si(temp1,FB->y0,b);
  mpz_add(temp1,temp2,temp1);
  mpz_abs(temp1,temp1);
  rSize=0;
  for (i=0; i<rfbTrialSize; i++) {
    e=0;
    p = FB->rfb[2*i];
    r = FB->rfb[2*i+1];
    if (((p!=r)&&((a - mulmod32(b, r, p))%p ==0))||((p==r)&&(b%p==0)))
      do {
        mpz_tdiv_q_ui(temp1, temp1, p);
        e++;
      } while (mpz_fdiv_ui(temp1, p)==0) ;
    if (e&&(rSize < MAX_FACTORS)) {
      factors[rSize] = i;
      exponents[rSize] = e;
      rSize++;
    }
  }
  for (i=0; i<MAX_LARGE_RAT_PRIMES; i++)
    R->p[i] = 1;
  numLarge=0; numpFacts=0;
  if (mpz_cmp_ui(temp1, 1)>0) {
    /* There are still factors left. Instead of trial dividing through */
    /* the rest of the factor base, just factor what's left, and find  */
    /* the factors.                                                    */
    numpFacts = factor(pFacts, temp1, 0);
    if (numpFacts <= 0)
      return numpFacts; /* There was a factor larger than 2^32. */
    for (i=0; i<numpFacts; i++) {
      /* Find this factor in the RFB. */
      loc = (s32 *)bsearch(&pFacts[i], FB->rfb, rfbSize, 2*sizeof(s32), cmpS32s);
      if (loc != NULL) {
	    assert((loc - FB->rfb) <= INT_MAX);
        locIndex = (s32)(loc - FB->rfb);
        locIndex /= 2; /* Each entry in the array is actually 2 s32's wide. */
        e=1;
        while ((i<(numpFacts-1)) && (pFacts[i] == pFacts[i+1])) {
          e++; i++;
        }
        factors[rSize] = locIndex;
        exponents[rSize] = e;
        rSize++;
      } else if ((numLarge < MAX_LARGE_RAT_PRIMES) && ((u32)pFacts[i]< FB->maxP_r)) {
        R->p[numLarge++] = pFacts[i];
      } else {
#ifdef GGNFS_VERBOSE
        fprintf(stderr, "%" PRId32 "\n", p);
#endif
        return -177;
      }
    }
  }
  if (rSize >= MAX_RAT_FACTORS)
    return -2;
  R->rFSize = rSize;
#ifdef _OLD
  for (i=0; i<rSize; i++) {
    R->rFactors[i] = factors[i];
    R->rExps[i] = exponents[i];
  }
#else
  memcpy(R->rFactors, factors, rSize*sizeof(s32));
  memcpy(R->rExps, exponents, rSize*sizeof(char));
#endif

  /******** Do the special primes ********/
  mpz_evalF(norm, R->a, R->b, FB->f);
  mpz_abs(norm, norm);

  idealHNF_ib_ab(&deltaHNF, R->a, R->b, N);

  /* Since we are actually getting valuations of (c_d*a - b\hat{alpha})
     = c_d(a - b\alpha)
     We have to subtract off the valuations of c_d there.
  */
  for (i=0; i<numSp; i++) {
    mpz_remove(norm, norm, N->sPrimes[i].p);
    e = valuation(&deltaHNF, &N->sPrimes[i], N) - N->v_cd_sPrimes[i];
    if (e) {
      factors[spSize] = i;
      exponents[spSize++] = e;   
    }
  }
  if (spSize > MAX_SP_FACTORS)
    return -5;
  R->spSize = spSize;
  for (i=0; i<spSize; i++) {
    R->spFactors[i] = factors[i];
    R->spExps[i] = exponents[i];
  }

  
  /********** Get the AFB part. **********/  
  i=0;
  while (i<afbTrialSize) {
    e=0;
    p = FB->afb[2*i];
    r = FB->afb[2*i+1];
    if ((r!=p)&& ((a-mulmod32(b,r,p))%p == 0)) {
      /**********************************************/
      /* This is a plain ol' (p,r) prime. It divides*/
      /* a-b\alpha iff a-br == 0 (mod p).           */
      /**********************************************/
      while (mpz_fdiv_q_ui(temp1, norm, p)==0) {
        e++;
        mpz_set(norm, temp1);
      }
    } else if ((r==p) && ((b%p)==0)) {
      /**********************************************/
      /* This is a `prime @ infinity', (p, \infty). */
      /* It divides a-b\alpha iff p|b.              */
      /**********************************************/
      while (mpz_fdiv_q_ui(temp1, norm, p)==0) {
        e++;
        mpz_set(norm, temp1);
      }
    }
    if (e&&(aSize < MAX_FACTORS)) {
      factors[aSize] = i;
      exponents[aSize++] = e;
    }
    i++;
  }
  for (i=0; i<MAX_LARGE_ALG_PRIMES; i++) {
    R->a_p[i] = 1; R->a_r[i] = 0;
  }
  numLarge=0; numpFacts=0;
  if (mpz_cmp_ui(norm, 1)>0) {
    /* There are still factors left. Instead of trial dividing through */
    /* the rest of the factor base, just factor what's left, and find  */
    /* the factors.                                                    */
    numpFacts = factor(pFacts, norm, 0);
    if (numpFacts <= 0)
      return numpFacts; /* There was a factor larger than 2^32. */
    for (i=0; i<numpFacts; i++) {
      /* Find this factor in the AFB. */
      p = pFacts[i];
	  assert(p > 0);
      if ((u32)p > FB->maxP_a) return -1923;
      loc = (s32 *)bsearch(&p, FB->afb, afbSize, 2*sizeof(s32), cmpS32s);
      /* There is only one alg. prime with this norm dividing <a-b\alpha>, */
      /* so  the exponent is easy to find:                                 */
      e=1;
      while ((i<(numpFacts-1))&&(pFacts[i+1]==p)) {
        e++; i++;
      }
      /* Find the corresponding 'r': */
      if (R->b%p==0) r=p; /* prime @ infty. */
      else r = mulmod32((p+(s32)(R->a%p))%p, inverseModP(R->b, p), p);
      if (loc != NULL) {
        /* Find the first AFB element corresponding to this 'p'. */
	    assert((loc - FB->afb) <= INT_MAX);
        locIndex = (s32)(loc - FB->afb);
        locIndex /= 2;
        while ((locIndex >0) && (FB->afb[2*(locIndex-1)]==p))
          locIndex--;
        while ((FB->afb[2*locIndex]==p) && (FB->afb[2*locIndex+1]!=r))
          locIndex++;
        if ((FB->afb[2*locIndex]==p) && (FB->afb[2*locIndex+1]==r)) {
          factors[aSize]=locIndex;
          exponents[aSize]=e;
          aSize++;
        } else {
          return -191; /* Why isn't it in the AFB? */
        }
      } else if ((numLarge < MAX_LARGE_ALG_PRIMES)&&((u32)p<FB->maxP_a)) {
        /* It's a large prime. */
        R->a_p[numLarge] = p;
        R->a_r[numLarge] = r;
        numLarge++;
      } else {
        return -19932;
      }
    }
  }
  if (aSize > MAX_ALG_FACTORS)
    return -5;
  R->aFSize = aSize;
  memcpy(R->aFactors, factors, aSize*sizeof(s32));
  memcpy(R->aExps, exponents, aSize*sizeof(char));
  
  /********** Get the QCB part. **********/  
  for (i=0; i<qSize; i++)
    R->qcbBits[i] = 0x00000000;
  
  for (i=0; i<FB->qcb_size; i++) {
    mpz_set_si(temp1, R->b);
    mpz_mul_ui(temp1, temp1, FB->qcb[2*i+1]);
    mpz_set_si64(temp2, R->a);
    mpz_sub(temp2, temp2, temp1);
    mpz_set_si(temp1, FB->qcb[2*i]);
    e = mpz_legendre(temp2, temp1);
    if (e==-1) 
      R->qcbBits[i/32] ^= BIT(i&0x1F);
  }
mpz_set_si(temp1, R->b);
mpz_mul(temp1, temp1, FB->y0);
mpz_set_si64(temp2, R->a);
mpz_mul(temp2, temp2, FB->y1);
mpz_sub(temp1, temp2, temp1);
if (mpz_sgn(temp1)<0)
  R->qcbBits[1] |= 0x10000000; 
else
  R->qcbBits[1] &= 0xEFFFFFFF; 
 
#if 0  
  R->qcbBits[1] |= 0x20000000; /* Cheap workaround for the sign of (a-bm). */
  R->qcbBits[0] |= 0x00000001; /* Cheap workaround for the sign of (a-bm). */
#endif
  return 0;
}


/***********************************************************************/
int completeRelFact(relation_t *R, nf_t *N)
/***********************************************************************/
/* R has it's (a,b) fields and rFactors, aFactors fields filled in.    */
/* This function will fill in the exponents, special primes, large     */
/* primes and quadratic characters.                                    */
/* Return value: 0 if it was a good relation. nonzero otherwise.       */
/***********************************************************************/
{ static int  initialized=0, qSize=0;
  static mpz_t     temp1, temp2, norm;
  static mpz_poly  delta;
  static mpz_mat_t deltaHNF;
  static s32 maxRFBPrime, maxAFBPrime, rfbSize, afbSize;
  s32   i, factors[MAX_FACTORS+1], b, fact;
  s64   a;
  s32   *loc, locIndex, r;
  s32   pFacts[10*MAX_FACTORS], p;
  int    numpFacts;
  char   exponents[MAX_FACTORS+1];
  int    spSize=0;
  int    e, numLarge, numSp=N->numSPrimes;
  nfs_fb_t *FB = N->FB;

  if (!(initialized)) {
    mpz_init(temp1); mpz_init(temp2);
    mpz_init(norm);  mpz_poly_init(delta);
    mpz_mat_init(&deltaHNF);
    rfbSize = FB->rfb_size;
    afbSize = FB->afb_size;
    qSize = FB->qcb_size/32;
    if (FB->qcb_size % 32)
      qSize++;
    if (qSize > 2) {
      qSize = 2;
      /* Let's take some liberties: */
      FB->qcb_size = 62;
    }
    maxRFBPrime = FB->rfb[2*(rfbSize-1)];
    maxAFBPrime = FB->afb[2*(afbSize-1)];
    initialized=1;
  }

  /********** Get the RFB part. **********/  
  /* Do temp1 <-- a - bm  */
  a = R->a; b = R->b;
  mpz_mul_si64(temp2,FB->y1,a);
  mpz_mul_si(temp1,FB->y0,b);
  mpz_add(temp1,temp2,temp1);
  mpz_abs(temp1,temp1);

  for (i=0; i<R->rFSize; i++) {
    e=0;
    fact = R->rFactors[i];
    p = FB->rfb[2*fact];
    r = FB->rfb[2*fact+1];
    if (((p!=r)&&((a - mulmod32(b, r, p))%p ==0))||((p==r)&&(b%p==0)))  {
      do {
        mpz_tdiv_q_ui(temp1, temp1, p);
        e++;
      } while (mpz_fdiv_ui(temp1, p)==0) ;
    }
    R->rExps[i] = e;
  }
  for (i=0; i<MAX_LARGE_RAT_PRIMES; i++)
    R->p[i] = 1;
  numLarge=0; numpFacts=0;
  if (mpz_cmp_ui(temp1, 1)>0) {
    numpFacts = factor(pFacts, temp1, 0);
    if (numpFacts <= 0)
      return numpFacts; /* There was a factor larger than 2^32. */
    for (i=0; i<numpFacts; i++) {
      /* Find this factor in the RFB. */
      if ((u32)pFacts[i] < (u32)maxRFBPrime) {
        loc = (s32 *)bsearch(&pFacts[i], FB->rfb, rfbSize, 2*sizeof(s32), cmpS32s);
        if (loc != NULL) {
		  assert((loc - FB->rfb) <= INT_MAX);
          locIndex = (s32)(loc - FB->rfb);
          locIndex /= 2; /* Each entry in the array is actually 2 s32's wide. */
          e=1;
          while ((i<(numpFacts-1)) && (pFacts[i] == pFacts[i+1])) {
            e++; i++;
          }
          R->rFactors[R->rFSize] = locIndex;
          R->rExps[R->rFSize] = e;
          R->rFSize++;
        }
      } else if ((numLarge < MAX_LARGE_RAT_PRIMES) && ((u32)pFacts[i]< (u32)FB->maxP_r)) {
        R->p[numLarge++] = pFacts[i];
      } else {
#ifdef GGNFS_VERBOSE
        fprintf(stderr, "%" PRId32 "\n", pFacts[i]);
#endif
        return -177;
      }
    }
  }

  /******** Do the special primes ********/
  mpz_evalF(norm, R->a, R->b, FB->f);
  mpz_abs(norm, norm);

  idealHNF_ib_ab(&deltaHNF, R->a, R->b, N);

  for (i=0; i<numSp; i++) {
    mpz_remove(norm, norm, N->sPrimes[i].p);
    e = valuation2(&deltaHNF, &N->sPrimes[i], N) - N->v_cd_sPrimes[i];
    if (e) {
      factors[spSize] = i;
      exponents[spSize++] = e;   
    }
  }
  if (spSize > MAX_SP_FACTORS)
    return -5;
  R->spSize = spSize;
  for (i=0; i<spSize; i++) {
    R->spFactors[i] = factors[i];
    R->spExps[i] = exponents[i];
  }


  /********** Get the AFB part. **********/  
if (mpz_sgn(norm)==0) {
printf("Error: norm=0!!!!\n");
exit(-1);
}
  for (i=0; i<R->aFSize; i++) {
    e=0;
    fact = R->aFactors[i];
    p = FB->afb[2*fact];
    r = FB->afb[2*fact+1];
    if (((p!=r)&&((a - mulmod32(b, r, p))%p ==0))||((p==r)&&(b%p==0)))  {
      while (mpz_fdiv_q_ui(temp1, norm, p)==0) {
        e++;
        mpz_set(norm, temp1);
      }
    }
    R->aExps[i] = e;
  }

  for (i=0; i<MAX_LARGE_ALG_PRIMES; i++) {
    R->a_p[i] = 1; R->a_r[i] = 0;
  }
  numLarge=0; numpFacts=0;
  if (mpz_cmp_ui(norm, 1)>0) {
    numpFacts = factor(pFacts, norm, 0);
    if (numpFacts <= 0)
      return numpFacts; /* There was a factor larger than 2^32. */
    for (i=0; i<numpFacts; i++) {
      /* Find this factor in the AFB. */
      p = pFacts[i];
      if ((u32)p> (u32)FB->maxP_a) return -1923;
      if ((u32)p < (u32)maxAFBPrime)
        loc = (s32 *)bsearch(&p, FB->afb, afbSize, 2*sizeof(s32), cmpS32s);
      else loc=NULL;
      /* There is only one alg. prime with this norm dividing <a-b\alpha>, */
      /* so  the exponent is easy to find:                                 */
      e=1;
      while ((i<(numpFacts-1))&&(pFacts[i+1]==p)) {
        e++; i++;
      }
      /* Find the corresponding 'r': */
      if (R->b%p==0) r=p; /* prime @ infty. */
      else r = mulmod32((p+(s32)(R->a%p))%p, inverseModP(R->b, p), p);
      if (loc != NULL) {
        /* Find the first AFB element corresponding to this 'p'. */
	    assert((loc - FB->afb) <= INT_MAX);
        locIndex = (s32)(loc - FB->afb);
        locIndex /= 2;
        while ((locIndex >0) && (FB->afb[2*(locIndex-1)]==p))
          locIndex--;
        while ((FB->afb[2*locIndex]==p) && (FB->afb[2*locIndex+1]!=r))
          locIndex++;
        if ((FB->afb[2*locIndex]==p) && (FB->afb[2*locIndex+1]==r)) {
          R->aFactors[R->aFSize] = locIndex;
          R->aExps[R->aFSize] = e;
          R->aFSize++;
        } else {
          return -191; /* Why isn't it in the AFB? */
        }
      } else if ((numLarge < MAX_LARGE_ALG_PRIMES)&&((u32)p<(u32)FB->maxP_a)) {
        /* It's a large prime. */
        R->a_p[numLarge] = p;
        R->a_r[numLarge] = r;
        numLarge++;
      } else {
        return -19932;
      }
    }
  }
  
  /********** Get the QCB part. **********/  
  for (i=0; i<qSize; i++)
    R->qcbBits[i] = 0x00000000;
  
  for (i=0; i<FB->qcb_size; i++) {
    mpz_set_si(temp1, R->b);
    mpz_mul_ui(temp1, temp1, FB->qcb[2*i+1]);
    mpz_set_si64(temp2, R->a);
    mpz_sub(temp2, temp2, temp1);
    mpz_set_si(temp1, FB->qcb[2*i]);
    e = mpz_legendre(temp2, temp1);
    if (e==-1) 
      R->qcbBits[i/32] ^= BIT(i&0x1F);
  }
mpz_set_si(temp1, R->b);
mpz_mul(temp1, temp1, FB->y0);
mpz_set_si64(temp2, R->a);
mpz_mul(temp2, temp2, FB->y1);
mpz_sub(temp1, temp2, temp1);
if (mpz_sgn(temp1)<0)
  R->qcbBits[1] |= 0x10000000; 
else
  R->qcbBits[1] &= 0xEFFFFFFF; 
#if 0
  R->qcbBits[1] |= 0x20000000; /* Cheap workaround for the sign of (a-bm). */
  R->qcbBits[0] |= 0x00000001; /* Cheap workaround for the sign of (a-bm). */
#endif
  return 0;
}

/***********************************************************************/
int cmpRelsAB(const void *X, const void *Y)
{ relation_t *x=(relation_t *)X, *y=(relation_t *)Y;

  if (x->b < y->b) return -1;
  if (x->b > y->b) return 1;
  if (x->a < y->a) return -1;
  if (x->a > y->a) return 1;
  return 0;
}

/************************************************************************/
int factRels_clsieved(relation_t *R, size_t numRels, nf_t *N)
/************************************************************************/
/* R[0], ..., R[numRels-1] have (a,b) fields that should be factored.   */
/* The b's should be all in a fairly narrow range (i.e., perhaps from a */
/* few lines of classical sieving). This function will compute the      */
/* factorizations of them by a sieving process. Any pairs which turn    */
/* out to be useless (not almost-smooth) will have their (a,b) fields   */
/* set to zero.                                                         */
/************************************************************************/
{ int  sorted=1, e, res;
  size_t relNum, numFact;
  unsigned int i, j;
  s32 b, p, r;
  s64 a;
  u32 rCutoff=300, aCutoff=300;
  s32 residue;
  s64 minA, maxA, minB, maxB;
  nfs_fb_t *FB = N->FB;

  /* Sort the relations lexicographically ascending with b>a. */
  /* Since we are using quick sort and there is a good chance     */
  /* that they are already sorted (which is bad for quicksort,    */
  /* I think), we should check to see if they are already sorted. */
  minA = maxA = R[0].a;
  minB = maxB = R[0].b;
  for (i=0; i<(numRels-1); i++)  {
    if ((R[i].b > R[i+1].b) || ((R[i].b==R[i+1].b) && (R[i].a > R[i+1].a)))
      sorted=0;
    minA = MIN(R[i].a, minA);
    maxA = MAX(R[i].a, maxA);
    minB = MIN(R[i].b, minB);
    maxB = MAX(R[i].b, maxB);
  }
  if (!sorted)
    qsort(R, numRels, sizeof(relation_t), cmpRelsAB);


  /* Initialize the factorizations. */
  for (i=0; i<numRels; i++) {
    R[i].rFSize = R[i].aFSize = 0;
    for (j=0; j<MAX_LARGE_RAT_PRIMES; j++)
      R[i].p[j] = 1;
    for (j=0; j<MAX_LARGE_ALG_PRIMES; j++)
      R[i].a_p[j] = R[i].a_r[j] = 1;
  }

  /* Trial division upto the cutoff for rational factors. */
  for (i=0; i<rCutoff; i++) {
    p = FB->rfb[2*i];
    r = FB->rfb[2*i+1];
    e=0;
    b = R[0].b;
    residue = mulmod32(b, r, p)%p;
    for (j=0; j<numRels; j++) {
      if (R[j].b != b) {
        b = R[j].b;
        residue = mulmod32(b, r, p)%p;
      }
      if ( (R[j].a - residue)%p ==0) {
        R[j].rFactors[R[j].rFSize] = i;
        R[j].rExps[R[j].rFSize] = 1;
        R[j].rFSize++;
      }
    }
  }
  /* Trial division upto the cutoff for algebraic factors. */
  for (i=0; i<aCutoff; i++) {
    p = FB->afb[2*i];
    r = FB->afb[2*i+1];
    e=0;
    b = R[0].b;
    residue = mulmod32(b, r, p)%p;
    for (j=0; j<numRels; j++) {
      if (R[j].b != b) {
        b = R[j].b;
        residue = mulmod32(b, r, p)%p;
      }
      if ( (R[j].a - residue)%p ==0) {
        R[j].aFactors[R[j].aFSize] = i;
        R[j].aExps[R[j].aFSize] = 1;
        R[j].aFSize++;
      }
    }
  }
  
  /**********************************************************************/
  /* These next two parts are only temporary - we will experiment with  */
  /* hash-table sieving and possibly other techniques to get it as fast */
  /* as possible. It is the bottleneck, though, as confirmed by gprof   */
  /* on some small numbers.                                             */
  /**********************************************************************/
  
  /********** Get the RFB parts. **********/
  for (i=rCutoff; i<FB->rfb_size; i++) {
    e=0;
    p = FB->rfb[2*i];
    r = FB->rfb[2*i+1];
    relNum = 0;

    /* Find the first a >= R[0].a so that (a, R[0].b) is divisible by (p,r). */
    b = R[0].b;
    
    a = R[0].a + (mulmod32(b, r, p) - (s32)(R[0].a%p) + p)%p;
    while (relNum < numRels) {
      if ((R[relNum].a==a) && (R[relNum].b==b))  {
        numFact = R[relNum].rFSize;
        R[relNum].rFactors[numFact] = i;
        R[relNum].rExps[numFact]=1; /* Another function will adjust exponents. */
        R[relNum].rFSize += 1;
        a += p; relNum++;
      }
      while ((relNum < numRels) && (R[relNum].a != a)) {
        if (R[relNum].a < a) relNum++;
        if (relNum < numRels) {
          if (R[relNum].b != b)  { /* Re-adjust a. */
            b = R[relNum].b;
            a = R[relNum].a + (mulmod32(b, r, p) - (s32)(R[relNum].a%p) + p)%p;
          }
          if (R[relNum].a > a) {
            do {
              a += p;
            } while (a < R[relNum].a);
          }
        }
      }
    }
  }


  /********** Get the AFB parts. **********/
  for (i=aCutoff; i<FB->afb_size; i++) {
    e=0;
    p = FB->afb[2*i];
    r = FB->afb[2*i+1];
    relNum = 0;

    /* Find the first a >= R[0].a so that (a, R[0].b) is divisible by (p,r). */
    b = R[0].b;
    
    a = R[0].a + (mulmod32(b, r, p) - (s32)(R[0].a%p) + p)%p;
    while (relNum < numRels) {
      while ((relNum < numRels) && (R[relNum].a != a)) {
        if (R[relNum].a < a) relNum++;
        if (relNum < numRels) {
          if (R[relNum].b != b)  { /* Re-adjust a. */
            b = R[relNum].b;
            a = R[relNum].a + (mulmod32(b, r, p) - (s32)(R[relNum].a%p) + p)%p;
          }
          if (R[relNum].a > a)
            do {
              a += p;
            } while (a < R[relNum].a);
        }
      }
      if ((relNum < numRels) && (R[relNum].a==a) && (R[relNum].b==b))  {
        numFact = R[relNum].aFSize;
        R[relNum].aFactors[numFact] = i;
        R[relNum].aExps[numFact]=1; /* Another function will adjust exponents. */
        R[relNum].aFSize += 1;
        a += p; relNum++;
      }
    }
  }

  /*******************************************/
  /* Finish the factorizations individually. */
  /* This will take care of the exponents on */
  /* factor base primes, the special primes  */
  /* and the large primes. It will zero out  */
  /* the `b' field of any relations that are */
  /* not (almost) smooth.                    */
  /*******************************************/
  for (i=0; i<numRels; i++) {
    res = completeRelFact(&R[i], N);
    if (res)  { /* Bad relation. */
      R[i].a = R[i].b = 0;
    }
  }

  return 0;
}


/***********************************************************************/
int completePartialRelFact(relation_t *R, nf_t *N, s32 rTDiv, s32 aTDiv)
/************************************************************************/
/* R has it's (a,b) fields and rFactors, aFactors fields filled in with */
/* a partial factorization (i.e., all of the high-end factors, but      */
/* probably none of the small ones.                                     */
/* This function will fill in the exponents, special primes, large      */
/* primes, quadratic characters, and missing factors.                   */
/* Return value: 0 if it was a good relation. nonzero otherwise.        */
/* We will only look for missing factors upto RFB[rTDiv] and AFB[aTDiv].*/
/************************************************************************/
{ static int  initialized=0, qSize=0;
  static mpz_t     temp1, temp2, norm, bMult;
  static mpz_poly  delta;
  static mpz_mat_t deltaHNF;
  static s32 maxRFBPrime, maxAFBPrime, rfbSize, afbSize;
  s32   i, factors[MAX_FACTORS+1], b, fact, r;
  s64   a;
  s32   locIndex;
  s32   pFacts[10*MAX_FACTORS], p;
  int    numpFacts, numFactors, kk;
  char   exponents[MAX_FACTORS+1];
  int    spSize=0, j;
  int    e, numLarge, numSp=N->numSPrimes;
  nfs_fb_t *FB = N->FB;

  if (!(initialized)) {
    mpz_init(temp1); mpz_init(temp2);
    mpz_init(norm);  mpz_poly_init(delta);
    mpz_init(bMult);
    mpz_mat_init(&deltaHNF);
    rfbSize = FB->rfb_size;
    afbSize = FB->afb_size;
    qSize = FB->qcb_size/32;
    if (FB->qcb_size % 32)
      qSize++;
    if (qSize > 2) {
      qSize = 2;
      /* Let's take some liberties: */
      FB->qcb_size = 62;
    }
    maxRFBPrime = FB->rfb[2*(rfbSize-1)];
    maxAFBPrime = FB->afb[2*(afbSize-1)];

    mpz_set(bMult, N->W_d);
    mpz_div(bMult, bMult, &N->W->entry[1][1]);

    initialized=1;
  }

  /********** Get the RFB part. **********/  
  /* Do temp1 <-- a - bm  */
  a = R->a; b = R->b;
  mpz_mul_si64(temp2,FB->y1,a);
  //printf("a = %" PRId64 "\n", a );
  //gmp_printf("y1*a=%Zd\n", temp2 );
 
  mpz_mul_si(temp1,FB->y0,b);
  //gmp_printf("y0*b=%Zd\n", temp1 );

  mpz_add(temp1,temp2,temp1);
  //gmp_printf("sum=%Zd\n", temp1 );

  mpz_abs(temp1,temp1);
  //gmp_printf("abs=%Zd\n", temp1 );

 
  /* Trial division to find the small factors: */
  numFactors=0;
  for (i=0; i<rTDiv; i++) {
    e=0;
    p = FB->rfb[2*i];
    r = FB->rfb[2*i+1];
    if (((p!=r)&&((a - mulmod32(b, r, p))%p ==0))||((p==r)&&(b%p==0)))  {
      while (mpz_fdiv_ui(temp1, p)==0) {
        mpz_tdiv_q_ui(temp1, temp1, p);
        e++;
      }
      if (e>0) {
        factors[numFactors] = i;
        exponents[numFactors]=e;
        numFactors++;
      }
    }
  }
  /* Now, check out the exponents for the factors which were given. */
  for (i=0; i<R->rFSize; i++) {
    e=0;
    fact = R->rFactors[i];
    p = FB->rfb[2*fact];
    r = FB->rfb[2*fact+1];
    if (((p!=r)&&((a - mulmod32(b, r, p))%p ==0))||((p==r)&&(b%p==0)))  {
      while (mpz_fdiv_ui(temp1, p)==0) {
        mpz_tdiv_q_ui(temp1, temp1, p);
        e++;
      }
      if (e>0) {
        factors[numFactors] = fact;
        exponents[numFactors] = e;
        numFactors++;
      }
    }
  }

  /* And the given large primes: */
  if (mpz_cmp_ui(temp1, 1)==0) {
    for (j=0; j<MAX_LARGE_RAT_PRIMES; j++)
      R->p[j] = 1;
  } else {
    for (j=0,kk=0,numpFacts=0; j<MAX_LARGE_RAT_PRIMES; j++) {
      pFacts[j] = R->p[j]; 
      numpFacts += (pFacts[j]>1);
      kk += (pFacts[j]>1);
      if (pFacts[j]==0) pFacts[j]=1;
      if (mpz_fdiv_ui(temp1, pFacts[j]) != 0) { 
#ifdef GGNFS_VERBOSE
        gmp_printf("temp1=%Zd, pFacts[%d]=%ld\n", temp1, j, pFacts[j]); 
#endif
	return -191; 
      }
      mpz_tdiv_q_ui(temp1, temp1, pFacts[j]);
    }
    if (mpz_cmp_ui(temp1, 1)) {
      /* This generally shouldn't happen - it's a kludge to fix my screw up, where I did a -dump
         but didn't dump the 3rd large prime from each side :)
      */
      if (mpz_probab_prime_p(temp1, 4) && (mpz_sizeinbase(temp1, 2) <= 32) && (kk < MAX_LARGE_RAT_PRIMES)) {
        pFacts[kk] = mpz_get_ui(temp1);
        R->p[kk] = pFacts[kk];
        numpFacts++;
        mpz_set_ui(temp1, 1);
      } else {
        return -193;
      }
    }
    /* Ok. so we know that these large primes do divide it,
       and that it will be a complete factorization after that.
       Let's just make sure, first, that they are actually large
       primes, and not actually from the RFB.
    */
    numLarge=0;
    for (i=0; i<numpFacts; i++) {
      /* Find this factor in the RFB. */
      if ((u32)pFacts[i] < (u32)maxRFBPrime) {
        locIndex = lookupRFB(pFacts[i], FB);
        if (locIndex >= 0) {
          e=1;
          while ((i<(numpFacts-1)) && (pFacts[i] == pFacts[i+1])) {
            e++; i++;
          }
          R->rFactors[R->rFSize] = locIndex;
          R->rExps[R->rFSize] = e;
          R->rFSize++;
        }

      } else if ((numLarge < MAX_LARGE_RAT_PRIMES) && ((u32)pFacts[i]< (u32)FB->maxP_r)) {
        R->p[numLarge++] = pFacts[i];
      } else {
#ifdef GGNFS_VERBOSE
        fprintf(stderr, "%" PRId32 " (numLarge=%d)\n", pFacts[i], numLarge);
#endif
        return -177;
      }
    }
  }
  /* Copy the factors and exponents back into R: */
  R->rFSize = numFactors;
  for (i=0; i<numFactors; i++) {
    R->rFactors[i] = factors[i];
    R->rExps[i] = exponents[i];
  }


  /******** Do the special primes ********/
  mpz_evalF(norm, R->a, R->b, FB->f);
  mpz_abs(norm, norm);

  mpz_set_si64(temp1, R->a);
  mpz_mul(temp1, temp1, &FB->f->coef[FB->f->degree]);
  mpz_mul_ui(temp2, bMult, R->b);

  for (i=0; i<numSp; i++) {
    mpz_remove(norm, norm, N->sPrimes[i].p);
    e = valuation_ab(temp1, temp2, i, N);
    if (e == -255) {
      /* This signals an error. Resort to the old slow version: */
      idealHNF_ib_ab(&deltaHNF, R->a, R->b, N);
      e = valuation2(&deltaHNF, &N->sPrimes[i], N);
    } 
    e -= N->v_cd_sPrimes[i];
    if (e) {
      factors[spSize] = i;
      exponents[spSize++] = e;   
    }
  }

  if (spSize > MAX_SP_FACTORS)
    return -5;
  R->spSize = spSize;
  for (i=0; i<spSize; i++) {
    R->spFactors[i] = factors[i];
    R->spExps[i] = exponents[i];
  }

  /********** Get the AFB part. **********/  

  /* Trial division to find the small factors: */
  numFactors=0;
  for (i=0; i<aTDiv; i++) {
    e=0;
    p = FB->afb[2*i];
    r = FB->afb[2*i+1];
    if ( (a - mulmod32(b, r, p))%p ==0)  {
      while (mpz_fdiv_ui(norm, p)==0) {
        mpz_tdiv_q_ui(norm, norm, p);
        e++;
      }
      if (e>0) {
        factors[numFactors] = i;
        exponents[numFactors]=e;
        numFactors++;
      }
    }
  }
  /* Now, check out the exponents for the factors which were given. */
  for (i=0; i<R->aFSize; i++) {
    e=0;
    fact = R->aFactors[i];
    p = FB->afb[2*fact];
    r = FB->afb[2*fact+1];
    if ( (a - mulmod32(b, r, p))%p ==0)  {
      while (mpz_fdiv_ui(norm, p)==0) {
        mpz_tdiv_q_ui(norm, norm, p);
        e++;
      }
      if (e>0) {
        factors[numFactors] = fact;
        exponents[numFactors] = e;
        numFactors++;
      }
    }
  }

  /* And the given large primes: */
  if (mpz_cmp_ui(norm, 1)==0) {
    for (j=0; j<MAX_LARGE_ALG_PRIMES; j++)
      R->a_p[j] = 1;
  } else {
    for (j=0, kk=0, numpFacts=0; j<MAX_LARGE_ALG_PRIMES; j++) {
      pFacts[j] = R->a_p[j]; 
      numpFacts += (pFacts[j]>1);
      kk += (pFacts[j] > 1);
      if (pFacts[j]==0) pFacts[j]=1;
      if (mpz_fdiv_ui(norm, pFacts[j]) != 0) { 
#ifdef GGNFS_VERBOSE
        gmp_printf("norm=%Zd, pFacts[%d]=%ld\n", norm, j, pFacts[j]); 
#endif
	return -291; 
      }
      mpz_tdiv_q_ui(norm, norm, pFacts[j]);
    }
    if (mpz_cmp_ui(norm, 1)) {
      /* This generally shouldn't happen - it's a kludge to fix my screw up, where I did a -dump
         but didn't dump the 3rd large prime from each side :)
      */
      if (mpz_probab_prime_p(norm, 4) && (mpz_sizeinbase(norm, 2) <= 32) && (kk < MAX_LARGE_ALG_PRIMES)) {
        pFacts[kk] = mpz_get_ui(norm);
        R->a_p[kk] = mpz_get_ui(norm);
        numpFacts++;
        mpz_set_ui(norm, 1);
      } else {
        printf("Algebraic failure: (%" PRId64 ", %" PRId32 ")\n", R->a, R->b);
        printf("  leftover norm is: "); mpz_out_str(stdout, 10, norm); printf("\n");
        printf("The following were the siever-supplied primes (norms):\n");
        for (i=0; i<R->aFSize; i++)
          printf("%" PRIx32 " (%" PRId32 ") ", R->aFactors[i], FB->afb[2*R->aFactors[i]]);
        for (i=0; i<MAX_LARGE_ALG_PRIMES; i++)
          printf("(%" PRId32 ")", R->a_p[i]);
        printf("\n");
        return -293;
      }
    }
    /* Ok. so we know that these large primes do divide it,
       and that it will be a complete factorization after that.
       Let's just make sure, first, that they are actually large
       primes, and not actually from the AFB.
    */
    numLarge=0;
    for (i=0; i<numpFacts; i++) {
      /* Find this factor in the AFB. */
      p = pFacts[i];
      if ((u32)p> (u32)FB->maxP_a) { 
#ifdef GGNFS_VERBOSE
        fprintf(stderr, "%" PRId32 "\n", p); 
#endif
	return -1923; 
      }
      /* Find the corresponding 'r': */
      if (R->b%p==0) r=p; /* prime @ infty. */
      else r = mulmod32((p+(s32)(R->a%p))%p, inverseModP(R->b, p), p);

      if ((u32)p < (u32)maxAFBPrime)
        locIndex = lookupAFB(p, r, FB);
      else locIndex=-1;
      /* There is only one alg. prime with this norm dividing <a-b\alpha>, */
      /* so  the exponent is easy to find:                                 */
      e=1;
      while ((i<(numpFacts-1))&&(pFacts[i+1]==p)) {
        e++; i++;
      }
      if (locIndex >= 0) {
        R->aFactors[R->aFSize] = locIndex;
        R->aExps[R->aFSize] = e;
        R->aFSize++;
      } else if ((numLarge < MAX_LARGE_ALG_PRIMES)&&((u32)p<(u32)FB->maxP_a)) {
        /* It's a large prime. */
        R->a_p[numLarge] = p;
        R->a_r[numLarge] = r;
        numLarge++;
      } else {
        return -19932;
      }
    }
  }
  /* Copy the factors and exponents back into R: */
  R->aFSize = numFactors;
  for (i=0; i<numFactors; i++) {
    R->aFactors[i] = factors[i];
    R->aExps[i] = exponents[i];
  }


  /********** Get the QCB part. **********/  
  for (i=0; i<qSize; i++)
    R->qcbBits[i] = 0x00000000;
  
  for (i=0; i<FB->qcb_size; i++) {
    mpz_set_si(temp1, R->b);
    mpz_mul_ui(temp1, temp1, FB->qcb[2*i+1]);
    mpz_set_si64(temp2, R->a);
    mpz_sub(temp2, temp2, temp1);
    mpz_set_si(temp1, FB->qcb[2*i]);
    e = mpz_legendre(temp2, temp1);
    if (e==-1) 
      R->qcbBits[i/32] ^= BIT(i&0x1F);
  }
  mpz_set_si(temp1, R->b);
  mpz_mul(temp1, temp1, FB->y0);
  mpz_set_si64(temp2, R->a);
  mpz_mul(temp2, temp2, FB->y1);
  mpz_sub(temp1, temp2, temp1);
  if (mpz_sgn(temp1)<0)
    R->qcbBits[1] |= 0x10000000; 
  else
    R->qcbBits[1] &= 0xEFFFFFFF; 
  /* And make sure we wind up with an even number of relations. */
  R->qcbBits[0] |= 0x00000001; 

#if 0
  R->qcbBits[1] |= 0x20000000; /* Cheap workaround for the sign of (a-bm). */
  R->qcbBits[0] |= 0x00000001; /* Cheap workaround for the sign of (a-bm). */
#endif
  return 0;
}

/**************************************************************/
void makeOutputLine(char *str, relation_t *R, nfs_fb_t *FB, int short_form)
/**************************************************************/
/* Turn a good relation with partial factorization into a     */
/* line of output, parseable by procrels.c (using the function*/
/* below).                                                    */
/**************************************************************/
{ int i, numR=0, numA=0;
  char s[128];

  sprintf(str, "%" PRId64 ",%" PRId32, R->a, R->b);

  if (short_form)
      return;

  strcat(str, ":");

  for (i=0; i<R->rFSize; i++) {
    if (R->rFactors[i] >= CLIENT_SKIP_R_PRIMES) {
      if (numR==0)
        sprintf(s, "%" PRIx32, FB->rfb[2*R->rFactors[i]]);
      else sprintf(s, ",%" PRIx32, FB->rfb[2*R->rFactors[i]]);
      strcat(str, s);
      numR++;
    }
  }

  for (i=0; i<MAX_LARGE_RAT_PRIMES; i++) {
    if (R->p[i] > 1) {
      if (numR>0)
        sprintf(s, ",%" PRIx32, R->p[i]);
      else
        sprintf(s, "%" PRIx32, R->p[i]);
      strcat(str, s);
      numR++;
    }
  }
  strcat(str, ":");
  for (i=0; i<R->aFSize; i++) {
    if (R->aFactors[i] >= CLIENT_SKIP_A_PRIMES) {
      if (numA==0)
        sprintf(s, "%" PRIx32, FB->afb[2*R->aFactors[i]]);
      else sprintf(s, ",%" PRIx32, FB->afb[2*R->aFactors[i]]);
      strcat(str, s);
      numA++;
    }
  }

  for (i=0; i<MAX_LARGE_ALG_PRIMES; i++) {
    if (R->a_p[i] > 1) {
      if (numA>0)
        sprintf(s, ",%" PRIx32, R->a_p[i]);
      else
        sprintf(s, "%" PRIx32, R->a_p[i]);
      strcat(str, s);
      numA++;
    }
  }
} 


/*************************************************************/
int parseOutputLine(relation_t *R, char *str, nfs_fb_t *FB)
/*************************************************************/
/* Take a line of input, in the siever-output format, and    */
/* peel off the information about (a,b) and known factors.   */
/* Put all the known info in R.                              */
/* Return value: 0 if everything is fine and this is         */
/*   potentially a good relation. Nonzero if it's bad for    */
/*   some reason.                                            */
/*************************************************************/
{ 
    size_t len = strlen(str);
    size_t i;
    size_t j;
    size_t size;
  
    int largeAlg, largeRat;
    char ab[128], rfb[256], afb[256];
  
    s32 p, r, k;
    s32 maxRFB, maxAFB;

    maxRFB = FB->rfb[2*(FB->rfb_size - 1)];
    maxAFB = FB->afb[2*(FB->afb_size - 1)];
 
    /* It's a shame - we could have parsed this with something similar
       to sscanf(str, "%[^:]:%[^:]:%[^:]", s1, s2, s3) ,
       but it doesn't handle the empty case well. There is surely some
       other standard library function capable of doing this parsing,
       but oh well.
    */
  
    i = 0;
    j = 0;
  
    while ((i < len) && (j < 127) && (str[i] != ':'))
    {
        ab[j++] = str[i++];
    }

    ab[j] = 0; 
    i++;

    j = 0;
  
    while ((i < len) && (j < 255) && (str[i] != ':'))
    {
        rfb[j++] = str[i++];
    }

    rfb[j] = 0; 
    i++;

    j = 0;
  
    while ((i < len) && (j < 255) && (str[i] != ':'))
    {
      afb[j++] = str[i++];
    }
  
    afb[j] = 0; 
    i++;

    if (sscanf(ab, "%" SCNd64 ",%" SCNd32, &R->a, &R->b) != 2) 
        return -1;

    /* Rational primes: */
    largeRat = 0;
  
    R->p[0] = 1;
    R->p[1] = 1;
  
    R->rFSize = 0; 
    j = 0;
  
    size = strlen(rfb);
  
    while ((j < size) && (sscanf(rfb + j,"%" SCNd32, &p) == 1)) 
    {
        k = lookupRFB(p, FB);
        if (k >= 0) 
        {
            R->rFactors[R->rFSize] = k;
            R->rFSize+=1;
        } 
        else 
        if ((p > maxRFB) && (largeRat < FB->maxLP)) 
        {
            R->p[largeRat++] = p;
        }
        while (isxdigit(rfb[j]))
        {
            j++;
        }

        j++; /* Pass the separator. */
    }

    /* Algebraic primes: */
    largeAlg = 0;
  
    R->a_p[0] = 1;
    R->a_p[1] = 1;

    R->a_r[0] = 1;
    R->a_r[1] = 1;

    R->aFSize = 0;
    j = 0;
  
    size = strlen(afb);

    while ((j < size) && (sscanf(afb + j,"%" SCNd32, &p) == 1)) 
    {
        if (R->b % p) 
        {
            r = mulmod32(p + (s32)(R->a%p), inverseModP(R->b, p), p);
            k = lookupAFB(p, r, FB);

            //printf("\n\n(%ld, %ld) : Look up of algebraic factor (%ld, %ld) gave k=%ld\n",
            //        R->a, R->b, p, r, k);
         
            if (k>=0) 
            {
                R->aFactors[R->aFSize] = k;
                R->aFSize++;
            } 
            else 
            if ((p > maxAFB) && (largeAlg < FB->maxLPA)) 
            {
                R->a_p[largeAlg] = p; 
                R->a_r[largeAlg++] = r;
            }
        } 

        while (isxdigit(afb[j])) 
        {
            j++;
        }

        j++; /* Pass the separator. */
    }

    return 0;
}

