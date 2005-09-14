/**************************************************************/
/* fbmisc.c                                                   */
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ggnfs.h"

#define ULL_NO_UL
#include "if.h"

/* The optimal value for this one is probably different than that in rels.c. */
/* The reason is that, when we are factoring things via factRel(), we know   */
/* that they are smooth. Here, on the other hand, they often will not be.    */
/* This means that more of the high-end time here is wasted time. So this    */
/* constant should probably always be smaller.                               */
#define FB_TRIAL_DIV_FRAC  0.5
#define MAX_TRIAL_DIV_SIZE 40000



/*********************************************************************/
void initFB(nfs_fb_t *FB)
/*********************************************************************/
{ int i;

  FB->name[0]=0;
  mpz_poly_init(FB->f);
  mpz_poly_init(FB->g);
  FB->skew = 0.0;
  mpz_init_set_ui(FB->n, 1);
  mpz_init_set_ui(FB->m, 1);
  mpz_init_set_ui(FB->y0, 1);
  mpz_init_set_ui(FB->y1, 1);
  mpz_init_set_ui(FB->disc, 0);
  mpz_init_set_ui(FB->cd, 1);
  mpz_init_set_ui(FB->knownDiv, 1);
  FB->rfb_size = FB->rLim = 0;
  FB->afb_size = FB->aLim = 0;
  FB->qcb_size = 0;
  FB->maxP_r = 1;
  FB->maxP_a = 1;
  FB->rfb = NULL;
  FB->rfb_log = NULL;
  FB->afb = NULL;
  FB->afb_log = NULL;
  FB->qcb = NULL;
  FB->maxLP = 0;
  FB->maxLPA = 0;
  for (i=0; i<MAXPOLYDEGREE; i++) {
    mpf_init2(FB->zeros[i].mpr, 128); 
    mpf_init2(FB->zeros[i].mpi, 128);
  }

}
        
/*********************************************************************/
int readPoly(FILE *fp, nfs_fb_t *FB)
/********************************************************************/
{ char token[128], value[512], thisLine[1024];
  u32_t i; 
  int cont=1 ,set_y=0, read[]={0,0,0,0,0,0,0};

  FB->f->degree = 0;
  for (i=0; i<=6; i++)
    mpz_set_ui(&FB->f->coef[i], 0);
  mpz_set_ui(FB->n, 0); mpz_set_ui(FB->m, 0); mpz_set_ui(FB->y1, 1);
  FB->skew = 100; /* A reasonable default for unknown cases. */
  while (cont) {
    readBinField(thisLine, 1023, fp);
    token[0]=value[0]=0;
    if (strlen(thisLine) && (thisLine[0] != '#')) {/* comment character. */
      for (i=0; i<strlen(thisLine)-1; i++) {
        if ((thisLine[i]==':') && (thisLine[i+1] != ' ')) {
          /* Insert a space. */
            memmove(thisLine+i+2, thisLine+i+1, strlen(thisLine)-i);
            thisLine[i+1]=' ';
            break;
        }
      }
      sscanf(thisLine, "%128s %512s", token, value);
      if (strncmp(token, "n:", 2)==0) {
        mpz_set_str(FB->n, value, 10);
      } else if (strncmp(token, "m:", 2)==0)  {
        mpz_set_str(FB->m, value, 10);
      } else if ((strncmp(token, "M", 1)==0)&&(strlen(token)==1)) {
        /* For (some) compatibility w/ Franke. */
        mpz_set_str(FB->m, value, 10);
      } else if (strncmp(token, "Y0:", 3)==0)  {
        mpz_set_str(FB->y0, value, 10);
        set_y |= 1;
      } else if (strncmp(token, "Y1:", 3)==0)  {
        mpz_set_str(FB->y1, value, 10);
        set_y |= 2;
      } else if (strncmp(token, "deg:", 4)==0)  {
        FB->f->degree = atoi(value);
      } else if (strncmp(token, "skew:", 5)==0) {
        FB->skew = atof(value);
      } else if (strncmp(token, "name:", 5)==0) {
        strncpy(FB->name, value, MAX_NAMESIZE);
      } else if ((token[0]=='c') && (token[1] >= '0') && (token[1] <= '6')) {
        mpz_set_str(&FB->f->coef[token[1]-'0'], value, 10);
        read[token[1]-'0']=1;
        FB->f->degree = MAX(FB->f->degree, (unsigned int)(token[1]-'0'));
      } else if ((token[0]=='X') && (token[1] >= '0') && (token[1] <= '6')) {
        /* For (some) compatibility w/ Franke. */
        mpz_set_str(&FB->f->coef[token[1]-'0'], value, 10);
        read[token[1]-'0']=1;
        FB->f->degree = MAX(FB->f->degree, (unsigned int)(token[1]-'0'));
      } else if (strncmp(token, "END_POLY",8)==0) {
        cont=0;
      } 
#ifdef _NO
      else {
        msgLog("", "Warning: readPoly() ignoring line:");
        msgLog("", thisLine);
      }
#endif
    }
    if (feof(fp)) cont=0;
  }
  if (mpz_sgn(FB->n)==0) {
    fprintf(stderr, "readPoly(): Did not find a valid 'n' value!\n");
    return -1;
  }
  if (set_y != 3){
    if (mpz_sgn(FB->m)==0) {
      fprintf(stderr, "readPoly(): Did not find a valid 'm' value!\n");
      return -1;
    }
    mpz_set(FB->y0, FB->m);
    mpz_neg(FB->y0, FB->y0);
    mpz_set_ui(FB->y1, 1);
  }
  if ((FB->f->degree < 2) || (FB->f->degree > 6)) {
    fprintf(stderr, "readPoly(): Poly degree %d invalid! Should be in [2,6]!\n",
          FB->f->degree);
    return -2;
  }
  return 0;
}
  
/*********************************************************************/
int writePoly(FILE *fp, nfs_fb_t *FB)
/********************************************************************/
{ int i;
  char str[1024], str2[1024];


  sprintf(str, "name: %s", FB->name); 
  writeBinField(fp, str); 
  sprintf(str, "n: %s", mpz_get_str(str2, 10, FB->n));
  writeBinField(fp, str);
  if(mpz_sgn(FB->m)){
    sprintf(str, "m: %s", mpz_get_str(str2, 10, FB->m));
    writeBinField(fp, str);
  }
  sprintf(str, "Y0: %s", mpz_get_str(str2, 10, FB->y0));
  writeBinField(fp, str);
  sprintf(str, "Y1: %s", mpz_get_str(str2, 10, FB->y1));
  writeBinField(fp, str);

  for (i=FB->f->degree; i>=0; i--) { 
    sprintf(str, "c%d: %s", i, mpz_get_str(str2, 10, &FB->f->coef[i]));
    writeBinField(fp, str);
  }
  sprintf(str, "skew: %1.1lf\n", FB->skew);
  writeBinField(fp, str);
  writeBinField(fp, "END_POLY");
  return 0;  
}

/* Bases will be chosen for logarithms with the goal of
   hitting these targets:
*/
#define RAT_LOG_TARGET 220
#define ALG_LOG_TARGET 200
/**************************************************/
int setLogs(nfs_fb_t *FB, s32 a0, s32 a1, s32 b0, s32 b1)
/**************************************************/
/* Initialize and compute values for              */
/* FB->*fb_log_base, FB->log_*lb and FB->*fb_log. */
/**************************************************/
{ s32   i, size;
  double t;
  mpz_t  tmp1, tmp2;

  if (!(FB->rfb_log = (int *)malloc(FB->rfb_size*sizeof(int)))) {
    fprintf(stderr, "setLogs(): Memory allocation error!\n");
    exit(-1);
  }
  if (!(FB->afb_log = (int *)malloc(FB->afb_size*sizeof(int)))) {
    fprintf(stderr, "setLogs(): Memory allocation error!\n");
    exit(-1);
  }
  mpz_init(tmp1); mpz_init(tmp2); 
  
  /* Decide on a base for the rational logs. */
  mpz_mul_si(tmp1, FB->y0, a1);
  mpz_abs(tmp1, tmp1);
  mpz_mul_si(tmp2, FB->y1, b1);
  mpz_abs(tmp2, tmp2);
  if (mpz_cmp(tmp1, tmp2)>0)
    t = mpz_sizeinbase(tmp1, 2)/M_LOG2E;
  else
    t = mpz_sizeinbase(tmp2, 2)/M_LOG2E;
  FB->rfb_log_base = exp(t/RAT_LOG_TARGET);
  FB->log_rlb = t/RAT_LOG_TARGET;

  /* Decide on a base for the algebraic logs. */
  /* Guess a rough max for the algebraic norms: */
  for (i=0; i<100; i++) {
    mpz_evalF(tmp2, a0+((rand()<<16)^rand())%(a1-a0), b1, FB->f);
    mpz_abs(tmp2, tmp2); 
    if (mpz_cmp(tmp2, tmp1) > 0) mpz_set(tmp1, tmp2);
  }
  t = mpz_sizeinbase(tmp1, 2)/M_LOG2E;
  FB->afb_log_base = exp(t/ALG_LOG_TARGET);
  FB->log_alb = t/ALG_LOG_TARGET;
  mpz_clear(tmp1); mpz_clear(tmp2);
  
  size = FB->rfb_size;
  for (i=0; i<size; i++) 
    FB->rfb_log[i] = fplog(FB->rfb[2*i], FB->log_rlb);
  
  size = FB->afb_size;
  for (i=0; i<size; i++) 
    FB->afb_log[i] = fplog(FB->afb[2*i], FB->log_alb);
  return 0;
}


/**************************************************************/
int isSmooth_rat(s32 a, s32 b, nfs_fb_t *FB)
/**************************************************************/
/* Input: a,b with b>0, and the factor base structure.        */
/* Return value:                                              */
/*  -1  :  a - b*m is not smooth over the RFB.                */
/*   k  :  a - b*m factors over the AFB with 'k' leftover     */
/*         prime factors below FB->maxP_r.                    */
/**************************************************************/
{ s32   i, fbSize;
  int    numFactors, numLarge;
  static mpz_t temp, temp2;
  static int initialized=0;
  u32   p[64], maxRFB;
  
  if (!(initialized)) {
    mpz_init(temp); mpz_init(temp2);
    initialized=1;
  }
  mpz_set_si(temp2, a);
  mpz_mul_si(temp, FB->m, b);
  mpz_sub(temp, temp2, temp);
  mpz_abs(temp, temp);  

  fbSize = (s32)(MIN(FB_TRIAL_DIV_FRAC*FB->rfb_size, MAX_TRIAL_DIV_SIZE));

/* Consider doing a mulmod32() here to check first. */
  for (i=0; i<2*fbSize; i+=2) {
    while (mpz_fdiv_ui(temp, FB->rfb[i])==0) 
      mpz_tdiv_q_ui(temp, temp, FB->rfb[i]);
  }
  if (mpz_cmp_ui(temp, 1)==0)
    return 0;
  numFactors = factor(p, temp, 1);
  if (numFactors <= 0)
    return -1; /* There was a factor larger than 2^32. */
  numLarge=0;
  maxRFB = FB->rfb[2*(FB->rfb_size-1)];
  for (i=0; i<numFactors; i++) {
    if (p[i] > FB->maxP_r) return -2;
    if (p[i] > maxRFB) numLarge++;
  }
  if (numLarge > FB->maxLP) return -3;
  return numLarge;
}

/**************************************************************/
int isSmooth_alg(s32 a, s32 b, nfs_fb_t *FB)
/**************************************************************/
/* Input: a,b with b>0, and the factor base structure.        */
/* Return value:                                              */
/*   -1 :  F(a,b) is not smooth over the AFB.                 */
/*    k :  F(a,b) factors over the AFB with 'k' leftover      */
/*         factors below FB->maxP_a.                          */
/**************************************************************/
{ s32   i, fbSize;
  u32  p[64], maxAFB;
  int    numFactors, numLarge;
  static int initialized=0;
  static mpz_t temp;
  
  if (!(initialized)) {
    mpz_init(temp);
    initialized=1;
  }
  mpz_evalF(temp, a, b, FB->f);
  mpz_abs(temp, temp);

  i=0;
  fbSize = (s32)(2*MIN(FB_TRIAL_DIV_FRAC*FB->afb_size, MAX_TRIAL_DIV_SIZE));
/* Again: consider adding a mulmod32() test here! */
  while (i<fbSize) {
    while (mpz_fdiv_ui(temp, FB->afb[i])==0)  {
      mpz_div_ui(temp, temp, FB->afb[i]);
    }
    i+=2;
  }
  
  if (mpz_cmp_ui(temp, 1)==0)
    return 0;
  numFactors = factor(p, temp, 1);
  if (numFactors <= 0)
    return -1; /* There was a factor larger than 2^32. */
  numLarge=0;
  maxAFB = FB->afb[2*(FB->afb_size-1)];
  for (i=0; i<numFactors; i++) {
    if (p[i] > FB->maxP_a) return -2;
    if (p[i] > maxAFB) numLarge++;
  }
  if (numLarge > FB->maxLPA) return -3;
  return numLarge;
}

/************************************************************/
int generateAFB(nfs_fb_t *FB, int verbose)
/************************************************************/
{ size_t  i; 
  u32  thisP;
  s32  zeros[MAXPOLYDEGREE], maxSize;
  u32  numZeros; 
  int d = FB->f->degree, cont;
  s32  total;
  char  str[128];    
  mpz_t cd;  
  u32  size=FB->afb_size, lim=FB->aLim;

  if (verbose) {
    printf("Generating AFB with norms upto %" PRIu32 "...\n", lim);
/*
    if (lim > 0) 
      printf("Generating AFB with norms upto %" PRId32 "...\n", lim);
    else
      printf("Generating AFB of size %" PRId32 "...\n", size);
*/
  }
  mpz_init_set(cd, &(FB->f->coef[d]));
  if (verbose)
    printf("Making algebraic factor base.\n");
  thisP = 2;
  total=0;
  if (lim > 0) 
    maxSize = getMaxP(1, lim) + d;
  else
    maxSize = size + d;
  
  if (!(FB->afb = (s32 *)malloc(maxSize*2*sizeof(s32)))) {
    fprintf(stderr, "generateAFB(): Memory allocation error!\n");
    exit(-1);
  }
  
  cont=1;
  while (cont) {
    if (verbose && ((total%10000)==0)) {
      sprintf(str, "Checking p=%" PRId32 "...(total=%" PRId32 ")", thisP, total);
      printf("%s",str); fflush(stdout);
      for (i=0; i<strlen(str); i++)
       printf("\b");
    }
#ifndef _OLD_ROOT_CODE
    numZeros = root_finder(zeros, (mpz_t*)FB->f->coef, d, thisP);
#else
    mpz_poly_modp(f_, FB->f, thisP);
    numZeros = poly_getZeros(zeros, f_, thisP);
#endif
    for (i=0; i<numZeros; i++) {
      FB->afb[2*total] = thisP;
      FB->afb[2*total+1] = zeros[i];
      total++;
    }
    if ((mpz_fdiv_ui(cd, thisP)==0) && (mpz_fdiv_ui(&(FB->f->coef[d-1]), thisP))) {
      /* This is a prime at infinity. */
      FB->afb[2*total] = thisP;
      FB->afb[2*total+1] = thisP;
      total++;
    }
    thisP = getNextPrime(thisP);
    if ((lim > 0) && (thisP > lim))
      cont=0;
/*
    else if ((lim <=0) && (total > size))
      cont=0;
*/
  }
  if (verbose) printf("\n");
    
  FB->afb_size = total;
  FB->aLim = FB->afb[2*(total-1)];
  mpz_clear(cd);
  return 0;
}


/************************************************************/
int generateQCB(nfs_fb_t *FB, int size)
/************************************************************/
/* There must already be an AFB in 'FB'.                    */
/************************************************************/
{ s32 i, thisP, zeros[MAXPOLYDEGREE], maxSize;
  u32 j;
  int  numZeros, usedZeros;
  s32 total;
  mpz_poly df;
  mpz_t    df_s, tmp, tmp2;
  poly_t   f_;
  
  
  total=0;
  maxSize = size + FB->f->degree;
  mpz_init(df_s);
  mpz_init(tmp);
  mpz_init(tmp2);
  for (i=0; i<=MAXPOLYDEGREE; i++)
    mpz_init(&(df->coef[i]));
  mpz_poly_diff(df, FB->f);
  
  if (!(FB->qcb = (s32 *)malloc(maxSize*2*sizeof(s32)))) {
    fprintf(stderr, "generateAFB(): Memory allocation error!\n");
    exit(-1);
  }

  /* The QCB will start right after the AFB. */
  /* For large prime variation, it should go after maxP_a. */
  thisP = 0x3FF00000;
  thisP = getNextPrime(thisP);
  
  while (total < size) {
    mpz_poly_modp(f_, FB->f, thisP);
    numZeros = poly_getZeros(zeros, f_, thisP);
    usedZeros = 0;
    for (i=0; i<numZeros; i++) {
      /* Check that f'(s) != 0. */
      mpz_set_ui(df_s, 0);
      mpz_set_ui(tmp, 1);
      for (j=0; j<=df->degree; j++) {
        mpz_mul(tmp2, tmp, &(df->coef[j]));
        mpz_add(df_s, df_s, tmp2);
        mpz_mod_ui(df_s, df_s, thisP);
        mpz_mul_ui(tmp, tmp, zeros[i]);
      }
      if (mpz_sgn(df_s) &&(total+usedZeros<size)) {
        FB->qcb[2*(total + usedZeros)] = thisP;
        FB->qcb[2*(total + usedZeros)+1] = zeros[i];
        usedZeros++;
      }
    }
    total += usedZeros;
    thisP = getNextPrime(thisP);
  }
  printf("done.\n");
    
  FB->qcb_size = total;
  mpz_clear(df_s);
  mpz_clear(tmp);
  mpz_clear(tmp2);
  for (i=0; i<=MAXPOLYDEGREE; i++)
    mpz_clear(&(df->coef[i]));
  return 0;
}


/************************************************************/
double getIPrimes(s32 *res, s32 numP, s32 upperBound, nfs_fb_t *FB)
/************************************************************/
/* Input: An empty array 'res' of size 'numP', and 'FB'.    */
/* Output: 'res' will contain the largest numP inert primes */
/*        below 'upperBound'. These are primes for which    */
/*        'f mod p' is irreducible and (p, c_d)=1.          */
/* Return value: The sum of the ln(p).                      */
/************************************************************/
{ s32   p, i, d;
  double sumLog=0.0;
  poly_t f;
        
  p = upperBound;
  d = FB->f->degree;
  i=0;
  while ((p>2) && (i<numP)) {
    p = getPrevPrime(p);
    mpz_poly_modp(f, FB->f, p);
    if (f->coef[d] && poly_irreducible_modp(f, p)) {
      res[i++] = p;
      sumLog += log((double)p);
    } 
  }
  return sumLog;              
    
}

/******************************************************/
int get_g(mpz_poly g, nfs_fb_t *FB) 
/******************************************************/
/* Compute the polynomial g = c_d^{d-1}f(x/c_d).      */
/* It is assumed that 'g' has been initialized.       */
/******************************************************/
{ int   i, d;
  mpz_t cdPow;

  mpz_init_set_ui(cdPow, 1);
  d = FB->f->degree;
  g->degree = d;
  mpz_set_ui(&(g->coef[d]), 1);
  for (i=d-1; i>=0; i--) { 
    mpz_mul(&(g->coef[i]), &(FB->f->coef[i]), cdPow);
    mpz_mul(cdPow, cdPow, &(FB->f->coef[d]));
  }  
  mpz_clear(cdPow);
  return 0;
}    

/*********************************************************************/
int loadFB(char *fName, nfs_fb_t *FB)
/********************************************************************/
{ size_t size;
  char token[128], value[512], thisLine[1024];
  int  cont=1, t;
  FILE *fp;

  if (!(fp = fopen(fName, "rb"))) {
    fprintf(stderr, "loadFB(): Could not open %s for read!\n", fName);
    return -1;
  }
  
  if (readPoly(fp, FB)) {
    fclose(fp);
    fprintf(stderr, "loadFB(): readPoly() reported an error!\n");
    return -1;
  }
  
  while (cont) {
    readBinField(thisLine, 1024, fp);      
    token[0]=value[0]=0;
    if (strlen(thisLine) && (thisLine[0] != '#')) {/* comment character. */
      sscanf(thisLine, "%128s %512s", token, value);
      if (strncmp(token, "npr:", 4)==0) 
        FB->maxLP = atoi(value);
      else if (strncmp(token, "mpr:", 4)==0) {
        t = atoi(value); t = MIN(t, 31);
        FB->maxP_r = (1<<t);
      } else if (strncmp(token, "npa:", 4)==0) 
        FB->maxLPA = atoi(value);
      else if (strncmp(token, "mpa:", 4)==0) {
        t = atoi(value); t = MIN(t, 31);
        FB->maxP_a = (1<<t);
      } else if (strncmp(token, "RFBsize:", 8)==0) 
        FB->rfb_size = atol(value);
      else if (strncmp(token, "AFBsize:", 8)==0) 
        FB->afb_size = atol(value);
      else if (strncmp(token, "END_HEADER",10)==0)
        cont=0;
    }
    if (feof(fp)) cont=0;
  }
  
  /********************************/
  /* get the rational factor base */
  /********************************/
  size = FB->rfb_size;
  if (!(FB->rfb = (s32 *)malloc(size*2*sizeof(s32)))) {
    fprintf(stderr, "loadfb(): Memory allocation error!\n");
    fclose(fp); return -1;
  }
  if (fread(FB->rfb, 2*sizeof(s32), size, fp) < size) {
    fprintf(stderr, "Error: factor base appears corrupt!\n");
    fclose(fp); return -1;
  }
  
  /*********************************/
  /* get the algebraic factor base */
  /*********************************/
  size = FB->afb_size;
  if (!(FB->afb = (s32 *)malloc(size*2*sizeof(s32)))) {
    fprintf(stderr, "loadfb(): Memory allocation error!\n");
    fclose(fp); return -1;
  }
  if (fread(FB->afb, 2*sizeof(s32), size, fp) < size) {
    fprintf(stderr, "Error: factor base appears corrupt!\n");
    fclose(fp); return -1;
  }
  FB->rLim = FB->rfb[2*(FB->rfb_size-1)];
  FB->aLim = FB->afb[2*(FB->afb_size-1)];
  fclose(fp);
  FB->qcb_size = 0;
  return 0;
}
  
  
    

/***********************************************/
int saveFB(char *fName, nfs_fb_t *FB)
/***********************************************/
{ s32  size;
  char str[1024];
  FILE *fp;
  
  if (!(fp = fopen(fName, "wb"))) {
    fprintf(stderr, "saveFB(): Error opening %s for write!\n", fName);
    return -1;
  }
  
  writePoly(fp, FB);
  
  sprintf(str, "npr: %d", FB->maxLP); writeBinField(fp, str);
  sprintf(str, "mpr: %d", (int)(0.5+log((double)FB->maxP_r)/M_LN2)); writeBinField(fp, str);
  sprintf(str, "npa: %d", FB->maxLPA); writeBinField(fp, str);
  sprintf(str, "mpa: %d", (int)(0.5+log((double)FB->maxP_a)/M_LN2)); writeBinField(fp, str);
  sprintf(str, "RFBsize: %" PRId32, FB->rfb_size); writeBinField(fp, str);
  sprintf(str, "AFBsize: %" PRId32, FB->afb_size); writeBinField(fp, str);
  sprintf(str, "END_HEADER"); writeBinField(fp, str);
  /*******************/
  /* output the RFB. */
  /*******************/

  size = FB->rfb_size;
  fwrite(FB->rfb, 2*sizeof(s32), size, fp);

  /*******************/
  /* output the AFB. */
  /*******************/
  size = FB->afb_size;
  fwrite(FB->afb, 2*sizeof(s32), size, fp);

  fclose(fp);
  return 0;
}
  
static __mpz_struct norms[MAX_PAIRS_PER_CALL];
static int normsInitialized=0;
/**************************************************************/
int isSmooth_rat_withInfo_par(relation_t *R, int numRels, nfs_fb_t *FB)
/**************************************************************/
/* Parallel version, for (a,b) pairs with the same b.         */
/**************************************************************/
{ s32   i, j, fbSize, locIndex, b, residue;
  u32   P;
  int    numFactors, numLarge, numGood=0;
  static mpz_t temp2;
  static int initialized=0;
  u32   p[64], maxRFB;
  
  fbSize = (s32)(MIN(FB_TRIAL_DIV_FRAC*FB->rfb_size, MAX_TRIAL_DIV_SIZE));
  if (!(initialized)) {
    mpz_init(temp2);
    if (!normsInitialized) {
      for (i=0; i<MAX_PAIRS_PER_CALL; i++)
        mpz_init(&norms[i]);
      normsInitialized=1;
    }
    initialized=1;
  }
  j=0;
  while ((j<numRels) && (R[j].b==0))
    j++;
  if (j >= numRels)
    return 0;
  b = R[j].b;
  mpz_mul_si(temp2, FB->m, b);
  for (j=0; j<numRels; j++) {
    if (R[j].b>0) {
      mpz_set_sll(&norms[j], R[j].a);
      mpz_sub(&norms[j], &norms[j], temp2);
      mpz_abs(&norms[j], &norms[j]);
      R[j].rFSize = 0;
    }
  }

  for (i=0; i<2*fbSize; i+=2) {
    P = FB->rfb[i];
    residue = mulmod32(b, FB->rfb[i+1], P);
    for (j=0; j<numRels; j++) {
      if ((R[j].a - residue)%P == 0) {
        R[j].rFactors[R[j].rFSize] = i/2;
        R[j].rFSize += 1;
        while (mpz_fdiv_ui(&norms[j], P)==0) 
          mpz_tdiv_q_ui(&norms[j], &norms[j], P);
      }
    }
  }
  maxRFB = FB->rfb[2*(FB->rfb_size-1)];
  for (j=0; j<numRels; j++) {
    R[j].p[0] = R[j].p[1] = 1;
    if (mpz_cmp_ui(&norms[j], 1) && (R[j].b>0)) {
      numFactors = factor(p, &norms[j], 1); /* I think the last arg could be a zero also. */
      if (numFactors <= 0) {
        R[j].b=0; /* There was a factor larger than 2^32. */
      } else  {
        numLarge=0;
        for (i=0; i<numFactors; i++) {
          P = p[i];
          if ((P > FB->maxP_r) || (P<0)) { /* Latter is for signed overflow checking. */
            R[j].b=0;
          } else if (P > maxRFB) {
            if (numLarge < FB->maxLP) {
              R[j].p[numLarge] = P;
              numLarge++;
            } else {
              R[j].b=0;
            }
          } else {
            /* It's in the RFB : find it. */
            locIndex = lookupRFB(P, FB);
            if (locIndex >= 0) { 
              R[j].rFactors[R[j].rFSize] = locIndex;
              R[j].rFSize += 1;
            } else { /* What the heck ?!?!#? */
              R[j].b=0;
            }
          }
        } /* for i=0..numFactors */
      }
    }
    if (R[j].b>0) numGood++;
  }
  while (j<numRels)
    R[j++].b=0;
  return numGood;
}

/***************************************************************/
int isSmooth_alg_withInfo_par(relation_t *R, int numRels, nfs_fb_t *FB)
/***************************************************************/
/* Parallel version of above, for (a,b) pairs with the same b. */
/***************************************************************/
{ s32   i, j, fbSize, r, locIndex, b, residue;
  u32  P;
  u32  p[64], maxAFB;
  int    numFactors, numLarge, numGood=0;
  static int initialized=0;
  static mpz_t temp;
  
  fbSize = (s32)(MIN(FB_TRIAL_DIV_FRAC*FB->afb_size, MAX_TRIAL_DIV_SIZE));
  if (!(initialized)) {
    mpz_init(temp);
    if (!normsInitialized) {
      for (i=0; i<MAX_PAIRS_PER_CALL; i++)
        mpz_init(&norms[i]);
      normsInitialized=1;
    }
    initialized=1;
  }
  j=0;
  while ((j<numRels) && (R[j].b==0))
    j++;
  if (j >= numRels)
    return 0;
  b = R[j].b;
  for (j=0; j<numRels; j++) {
    /* Consider making a new parallel evaluation function - it might be
       possible to do it much faster. But it might not have much runtime
       impact - I'll have to do another profile to find out.
    */
    if (R[j].b>0) {
      mpz_evalF(&norms[j], R[j].a, b, FB->f);
      mpz_abs(&norms[j], &norms[j]);
      R[j].aFSize = 0;
    }
  }
 
  for (i=0; i<2*fbSize; i+=2) {
    P = FB->afb[i];
    residue = mulmod32(b, FB->afb[i+1], P);
    for (j=0; j<numRels; j++) {
      if ((R[j].a - residue)%P == 0) {
        R[j].aFactors[R[j].aFSize] = i/2;
        R[j].aFSize += 1;
        while (mpz_fdiv_ui(&norms[j], P)==0) 
          mpz_tdiv_q_ui(&norms[j], &norms[j], P);
      }
    }
  }
  maxAFB = FB->afb[2*(FB->afb_size-1)];
  for (j=0; j<numRels; j++) {
    R[j].a_p[0] = R[j].a_p[1] = 1;
    if (mpz_cmp_ui(&norms[j], 1) && (R[j].b>0)) {
      numFactors = factor(p, &norms[j], 1); /* I think the last arg could be a zero also. */
      if (numFactors <= 0)
        R[j].b=0; /* There was a factor larger than 2^32. */
      else {
        numLarge=0;
        for (i=0; (i<numFactors)&&(R[j].b>0); i++) {
          P = p[i];
          if ((P > FB->maxP_a) || (P<0)) R[j].b=0;
          else if (P > maxAFB) {
            if (numLarge < FB->maxLPA) {
              R[j].a_p[numLarge] = P;
              R[j].a_r[numLarge] = mulmod32((P + (s32)(R[j].a%P))%P, inverseModP(b, P), P);
              numLarge++;
            } else R[j].b=0;
          } else {
            /* It's in the AFB : find it. */
            if (b%P==0) {
              r=P; /* prime @ infty. Don't do anything with it - the server will handle it. */
            } else {
              r = mulmod32( ((P+(s32)(R[j].a%P)) )%P, inverseModP(b, P), P);
              locIndex = lookupAFB(P, r, FB);
              if ((FB->afb[2*locIndex]==P) && (FB->afb[2*locIndex+1]==r)) {
                R[j].aFactors[R[j].aFSize] = locIndex;
                R[j].aFSize++;
              } else {
                printf("Strange error!?!?\n");
                R[j].b=0;
              }
            } /* for i=0..numFactors */
          }
        }
      }
    }
    if (R[j].b>0) numGood++;
  }
  while (j<numRels)
    R[j++].b=0;
  return numGood;
}
