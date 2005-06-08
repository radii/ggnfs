/**************************************************************/
/* sqrt.c                                                     */
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
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef _MSC_VER  
#include <sys/time.h>
#endif
#include "ggnfs.h"

#define QCB_SIZE 50


#define USAGE " -fb <fname> -prel <prefix> -deps <fname> -depnum <int>\n"\
"-fb     <fname>  : File containing the factor base.\n"\
"-deps   <fname>  : Output of `procrels'. i.e., the file with the dependencies.\n"\
"-depnum <int>    : Which dependency to try.\n"\
"-knowndiv <int>  : The product of known small divisors of n\n"\
"-nodfactor       : Don't try to factor the discriminant again.\n"

#define START_MSG \
" __________________________________________________________ \n"\
"|          This is the sqrt program for GGNFS.             |\n"\
"| Version: %-15s                                 |\n"\
"| This program is copyright 2004, Chris Monico, and subject|\n"\
"| to the terms of the GNU General Public License version 2.|\n"\
"|__________________________________________________________|\n"




typedef struct {
  s32 numRels;
  s32 numPrimes;
  s32 Rels[MAX_RELS_IN_FF];
  s32 QCB[2];
  char sign;
  s32 rows[MAX_ROWS_IN_COL];
} column_t;

/******************************************************/
void readColIndex(column_t *C, FILE *fp)
/******************************************************/
/* This is the indicies of the relations contributing */
/* to the column. We don't care about anything else.  */
/******************************************************/
{
  fread(&C->numRels, sizeof(s32), 1, fp);
  if (C->numRels > 0)
    fread(&C->Rels, sizeof(s32), C->numRels, fp);
}

/*****************************************************/
s32 reduceToOdd(s32 *x, s32 size)
/*****************************************************/
/* Given a list of s32's, remove repeated ones in   */
/* pairs, so that we get only the ones that occurred */
/* an odd number of times.                           */
/* Return value: number of entries remaining.        */
/*****************************************************/
{ int j, k, unique;

  qsort(x, size, sizeof(s32), cmpS32s);
  
  j=k=unique=0;
  while (j<size) {
    k=1;
    while (((j+k)<size) && (x[j] == x[j+k]))
      k++;
    if (k&0x01) {
      /* It occurrs an odd number of times, so keep it. */
      x[unique++] = x[j];
    }
    j += k;
  }
  return unique;
}

/**********************************************/
int cmpRels(const void *a, const void *b)
/**********************************************/
/* Compare two relations w.r.t revlex on their*/
/* corresponding (a,b) pairs (for filtering   */
/* out duplicates).                           */
/**********************************************/
{ relation_t *A = (relation_t *)a, *B = (relation_t *)b;

  if (A->b < B->b) return -1;
  else if (A->b > B->b) return 1;
  else if (A->a < B->a) return -1;
  else if (A->a > B->a) return 1;
  return 0;
}
  
/****************************************************/
int main(int argC, char *args[])
/****************************************************/
{ char       fbName[64], depName[64], colIndex[64];
  char       *rid_hash, str[1024], token[512], value[512];
  mpz_t      p, q, rSqrt, aSqrt, kDiv;
  double     startTime, now;
  FILE       *fp;
  nf_t       N;
  mpz_fact_t D;
  int        depNum=-1, res=0, cont;
  s32      *colsInDep=NULL, maxCols, numCols, t, i, j, k;
  s32      *relsInDep=NULL, maxRels, numRels;
  struct     stat fileInfo;
  multi_file_t prelF, lpF;
  column_t   C;
  int        discFact=0;

  fbName[0] = prelF.prefix[0] = 0;
  mpz_init_set_ui(kDiv, 1); 

  printf(START_MSG, GGNFS_VERSION);
 
  for (i=1; i<argC; i++) {
    if (strcmp(args[i], "-fb")==0) {
      if ((++i) < argC) 
        strncpy(fbName, args[i], 64);
    } else if (strcmp(args[i], "-deps")==0) {
      if ((++i) < argC) 
        strncpy(depName, args[i], 64);
    } else if (strcmp(args[i], "-depnum")==0) {
      if ((++i) < argC) 
        depNum = atoi(args[i]);
    } else if (strcmp(args[i], "-knowndiv")==0) {
      if ((++i) < argC)
        mpz_set_str(kDiv, args[i], 10);
    } else if (strcmp(args[i], "-nodfactor")==0) {
      discFact=0;
    }
  }
 
  if ((fbName[0]==0) || (depNum < 0)) {
    printf("USAGE: %s %s\n", args[0], USAGE);
    exit(0);
  }
  if (depNum > 31) {
    printf("depNum=%d is invalid. It should be in [0,31].\n", depNum);
    exit(0);
  }
  msgLog("", "GGNFS-%s : sqrt", GGNFS_VERSION);


  mpz_fact_init(&D);
  initNF(&N);
  if (!(N.FB = (nfs_fb_t *)malloc(sizeof(nfs_fb_t)))) {
    fprintf(stderr, "Error: Could not allocate %d bytes for N.FB!\n",
            sizeof(nfs_fb_t));
    clearNF(&N); mpz_fact_clear(&D); exit(-4);
  }
  initFB(N.FB);
  if (loadFB(fbName, N.FB)) {
    printf("Could not load FB from %s!\n", fbName);
    clearNF(&N); mpz_fact_clear(&D); exit(-4);
  }
  mpz_set(N.FB->knownDiv, kDiv);

  startTime = sTime();
  printf("Setting up nf_t object...\n");
  mpz_poly_print(stdout, "F_1 = ", N.FB->f);
  printf("F_2 = ");
  mpz_out_str(stdout, 10, N.FB->y0);
  if (mpz_sgn(N.FB->y1)>0) printf(" + ");
  mpz_out_str(stdout, 10, N.FB->y1);
  printf("X\n");


  mpz_poly_cp(N.f, N.FB->f);
  get_g(N.T, N.FB);
  mpz_poly_discrim(D.N, N.T);
  mpz_fact_factorEasy(&D, D.N, discFact);
  getIntegralBasis(&N, &D, discFact);

  printf("Reading dependency %d from file %s...\n", depNum, depName);
  if (!(fp = fopen(depName, "rb"))) {
    fprintf(stderr, "Error opening %s for read!\n", depName);
    res = -1; goto SS_DONE;
  }

  cont=1;
  while (cont) {
    str[0]=token[0]=value[0]=0;
    readBinField(str, 512, fp);
    sscanf(str, "%256s %256s", token, value);
    if (strncmp(token, "NUMCOLS:", 8)==0) sscanf(value, "%lx", &maxCols);
    else if (strncmp(token, "COLNAME:", 8)==0) sscanf(value, "%s", colIndex);
    else if (strncmp(token, "MAXRELS:", 8)==0) sscanf(value, "%lx", &maxRels);
    else if (strncmp(token, "RELPREFIX:", 10)==0) sscanf(value, "%s", prelF.prefix);
    else if (strncmp(token, "RELFILES:", 9)==0) sscanf(value, "%x", &prelF.numFiles);
    else if (strncmp(token, "LPFPREFIX:", 10)==0) sscanf(value, "%s", lpF.prefix);
    else if (strncmp(token, "LPFFILES:", 8)==0) sscanf(value, "%x", &lpF.numFiles);
    else if (strncmp(token, "END_HEADER",10)==0) cont=0;
    if (feof(fp)) cont=0;
  } 

  if (!(colsInDep = (s32 *)malloc(maxCols*sizeof(s32)))) {
    fclose(fp);
    fprintf(stderr, "Error allocating %ld bytes for columns in dependency!\n", 
            maxCols*sizeof(s32));
    res = -1; goto SS_DONE;
  }
  for (i=0, numCols=0; i<maxCols; i++) {
    fread(&t, sizeof(s32), 1, fp);
    if (t&BIT(depNum)) 
      colsInDep[numCols++] = i;
  }
  fclose(fp);

  printf("NUMCOLS = %ld\n", maxCols);
  printf("COLNAME = %s\n", colIndex);
  printf("MAXRELS = %ld\n", maxRels);
  printf("RELPREFIX = %s\n", prelF.prefix);
  printf("RELFILES = %d\n", prelF.numFiles);
  printf("LPFPREFIX = %s\n", lpF.prefix);
  printf("LPFFILES = %d\n", lpF.numFiles);
  printf("There are %ld columns in this dependency. Getting corresponding (a,b) pairs...\n", numCols);

  if (stat(colIndex, &fileInfo)) {
    fprintf(stderr, "Could not stat column index file %s!\n", colIndex);
    res=-1; goto SS_DONE;
  }

  /********************************************************************/
  /* Now, open and scan the colIndex file to find out which relations */
  /* go with the columns in our dependency - do this by hash, so that */
  /* we don't need to consider an (a,b) pair more than once. That is, */
  /* it's possible that the dependency uses an (a,b) pair multiple    */
  /* times, but it always suffices to reduce this to 0 or 1 times.    */
  /********************************************************************/
  if (!(rid_hash = (char *)malloc(maxRels*sizeof(char)))) {
    fprintf(stderr, "Error allocating %ld bytes for rid_hash!\n", maxRels*sizeof(char));
    res = -1; goto SS_DONE;
  }
  for (i=0; i<maxRels; i++)
    rid_hash[i] = 0x00;
  if (!(fp = fopen(colIndex, "rb"))) {
    fprintf(stderr, "Error opening column index file %s for read!\n", colIndex);
    res = -1; goto SS_DONE;
  }
  for (i=0, j=0; i<numCols; i++) {
    /* j is the column number waiting on 'fp', and colsInDep[i] is         */
    /* the next column we need (already sorted b/c the way it was read in. */
    while (j <= colsInDep[i]) {
      readColIndex(&C, fp);
      j++;
    }
    for (k=0; k<C.numRels; k++) {
      if (C.Rels[k] < maxRels)
        rid_hash[C.Rels[k]] ^= 0x01;
      else {
        fprintf(stderr, "Error: Column claims use of relation %ld (maxRels = %ld)!\n",
                C.Rels[k], maxRels);
  
        exit(-4);
      }
    }
  }
  fclose(fp);
  free(colsInDep); colsInDep = NULL; /* No longer needed. */
  numRels = 0;
  for (i=0; i<maxRels; i++) {
    if (rid_hash[i]==0x01)
      numRels++;
  }
  printf("This dependency consists of %ld (a,b) pairs.\n", numRels);  
  
  if (!(relsInDep = (s32 *)malloc(sizeof(s32)*(numRels+1)))) {
    fprintf(stderr, "Error allocating %ld bytes for relsInDep!\n", (numRels+1)*sizeof(s32));
    res = -1; goto SS_DONE;
  }
  numRels = 0;
  for (i=0; i<maxRels; i++) {
    if (rid_hash[i]==0x01) 
      relsInDep[numRels++] = i;
  }
  relsInDep[numRels] = -1; /* Terminator. */
  free(rid_hash); /* No longer needed. */


  mpz_init(rSqrt); mpz_init(aSqrt);
  /* The montgomerySqrt() call goes here. */
  montgomerySqrt(rSqrt, aSqrt, relsInDep, &prelF, &lpF, N.FB, &N);
  mpz_init(p); mpz_init(q);
  mpz_sub(p, rSqrt, aSqrt);
  mpz_gcd(p, p, N.FB->n);
  mpz_div(q, N.FB->n, p);
  printf("Square root computations result in N=(r1)(r2) where:\n");
  printf("r1 = "); mpz_out_str(stdout, 10, p); printf("\n");
  printf("r2 = "); mpz_out_str(stdout, 10, q); printf("\n");
  now = sTime();
  msgLog("", "From dependence %d, sqrt obtained:", depNum);
/*
   res = -2 : not factored
   res =  0 : completely factored
   res =  1 : incompletely factored
*/
  res = -2;
  mpz_get_str(str, 10, p);
  if ((mpz_cmp_ui(p, 1)>0) && (mpz_cmp(p, N.FB->n)<0)) {
    res = 0;
    if (mpz_probab_prime_p(p, 10))
      sprintf(str, "%s (pp%d)", str,strlen(str));
    else {
      sprintf(str, "%s (c%d)", str,strlen(str));
      res=1;
    }
  }
  msgLog("", "  r1=%s", str,strlen(str));
  mpz_get_str(str, 10, q);
  if ((mpz_cmp_ui(q, 1)>0) && (mpz_cmp(q, N.FB->n)<0)) {
    if (mpz_probab_prime_p(q, 10)) 
      sprintf(str, "%s (pp%d)", str,strlen(str));
    else {
      sprintf(str, "%s (c%d)", str,strlen(str));
      res=1;
    }
  }
  msgLog("", "  r2=%s", str);
  if (res>=0) 
    msgLog("", "(pp=probable prime, c=composite)");
  printf("Elapsed time: %1.4lf seconds.\n", now - startTime);
  msgLog("", "sqrtTime: %1.1lf", now - startTime);

  mpz_clear(p); mpz_clear(q);

SS_DONE:
  return res;
}  
