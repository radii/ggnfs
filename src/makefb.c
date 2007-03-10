/**************************************************************/
/* makefb.c                                                   */
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
#include <gmp.h>
#include <time.h>
#include "ggnfs.h"

#define MINRFBSIZE 5
#define MAXRFBSIZE 5000000
#define MINAFBSIZE 5
#define MAXAFBSIZE 5000000

#define MIN_RLIM 11
#define MAX_RLIM 100000000
#define MIN_ALIM 11
#define MAX_ALIM 100000000

#define USAGE \
"[OPTIONS]\n"\
"--help            : show this help and exit\n"\
"-of <filename>    : save factor base to <filename>\n"\
"-if <filename>    : read n, m, f from <filename>\n"\
"-rl <limit>       : create a rational factor base upto <limit>\n"\
"-al <limit>       : create an algebraic factor base with norms upto <limit>\n"\
"-mpr <long>       : max large rational prime for large-prime variation\n"\
"-mpa <long>       : max large algebraic prime for large-prime variation\n"\
"-lpbr <int>       : max bits in large rational prime.\n"\
"-lpba <int>       : max bits in large algebraic prime.\n"\
"-2p               : Upto 2 large rat. and algebraic primes (default is 1).\n"\
"-3p               : Upto 3 large rat. and algebraic primes (default is 1).\n"\
"  The next two options are for backward-compatibility:\n"\
"-rs <size>        : create a rational factor base of size <size>\n"\
"-as <size>        : create an algebraic factor base of size <size>\n"


/*********************************************************************/
int makefb(nfs_fb_t *FB)
/*********************************************************************/
{ s32  size, *tmpPtr;
  u32  tmp_x;
  int   status, i, j;
  mpz_t eval, mpow, tmp, neg_y0;

  if (mpz_cmp_ui(FB->n, 1) <= 0) {
    fprintf(stderr, "makefb() error: invalid input FB->n = ");
    mpz_out_str(stderr, 10, FB->n);
    fprintf(stderr, "\n");
    return -1;
  }

  if (mpz_probab_prime_p(FB->n,10)) {
    fprintf(stderr, "makefb() error: input is probably prime FB->n = ");
    mpz_out_str(stderr, 10, FB->n);
    fprintf(stderr, "\n");
    return -1;
  }

  if (mpz_perfect_power_p(FB->n)) {
    fprintf(stderr, "makefb() error: input is a perfect power FB->n = ");
    mpz_out_str(stderr, 10, FB->n);
    fprintf(stderr, "\n");    
    return -1;
  }

  /* Verify the polynomial. */
  mpz_init(eval); mpz_init(mpow); mpz_init(tmp); mpz_init(neg_y0);
  mpz_set_ui(eval, 0);
  mpz_set_ui(mpow, 1);
  mpz_neg(neg_y0, FB->y0);
  
  for(i=0; i<= FB->f->degree; i++) {
    mpz_mul(tmp, mpow, &(FB->f->coef[i]));
      for(j=i; j < FB->f->degree; j++)
        mpz_mul(tmp, tmp, FB->y1);
    mpz_add(eval, eval, tmp);
    mpz_mul(mpow, mpow, neg_y0);
  }
  if (mpz_sgn(eval)==0) {
    printf("Error: Bad polynomial: f(m)=0.\n");
    printf("m = "); mpz_out_str(stdout, 10, FB->m); printf("\n");
    exit(-1);
  }
  mpz_mod(eval, eval, FB->n);
  if (mpz_sgn(eval)) {
    printf("Error: Bad polynomial: f(m) !=0 (mod n)\n");
    printf("m = "); mpz_out_str(stdout, 10, FB->m); printf("\n");
    exit(-1);
  }
  
  /********************************/
  /* get the rational factor base */
  /********************************/
  if (FB->rLim > 0) {
    size = getMaxP(1, FB->rLim);
    tmpPtr = malloc(size*sizeof(s32));
    size = pSieve(tmpPtr, size, 1, FB->rLim);
    FB->rfb_size = size;
  } else {
    size = FB->rfb_size;
    tmpPtr = getPList(&size);
  }
  if (size < FB->rfb_size) {
    fprintf(stderr, "Error getting rational primes!\n");
    return -1;
  }
  if (!(FB->rfb = (s32 *)malloc(size*2*sizeof(s32)))) {
    fprintf(stderr, "Memory allocation error!\n");
    exit(-1);
  }
  for (i=0; i<size; i++) {
    FB->rfb[2*i] = tmpPtr[i];
    if ((tmp_x = mpz_fdiv_ui(FB->y1, tmpPtr[i])))
      FB->rfb[2*i+1] = mulmod32(mpz_fdiv_ui(neg_y0, tmpPtr[i]), inverseModP(tmp_x, tmpPtr[i]), tmpPtr[i]);
    else FB->rfb[2*i+1] = tmpPtr[i];
  }
  free(tmpPtr);


  size = FB->afb_size;
  /*********************************/
  /* get the algebraic factor base */
  /*********************************/
  status = generateAFB(FB, 1);
  if (status) {
    fprintf(stderr, "Error getting algebraic factor base!\n");
    return -1;
  }
  
  return 0;
}


  
/*************************************************************************/
#ifdef _MAKEFB_STANDALONE
int main(int argC, char *args[])
{ char     ifname[MAXFNAMESIZE], ofname[MAXFNAMESIZE], str[1024];
  FILE    *ifp;
  int      i, status, twoP=0, threeP=0;
  time_t   start, stop;
  nfs_fb_t _FB, *FB;

  ifname[0]=ofname[0]=0;
  initFB(&_FB);
  FB = &_FB;
  
  /***************************/
  /* Parse command-line args */
  /***************************/
  for (i=1; i<argC; i++) {
    if (strcmp(args[i], "--help")==0) {
      printf("%s\n%s\n", args[0], USAGE);
      return 1;
    }
    else if (strcmp(args[i], "-if")==0) {
      if ((++i)<argC) {
        strcpy(ifname, args[i]);
      }
    }
    else if (strcmp(args[i], "-of")==0) {
      if ((++i)<argC) {
        strcpy(ofname, args[i]);
      }
    }
    else if (strcmp(args[i], "-as")==0) {
      if ((++i)<argC) {
        FB->afb_size = atol(args[i]);
      }
    }
    else if (strcmp(args[i], "-al")==0) {
      if ((++i)<argC) {
        FB->aLim = atol(args[i]);
      }
    }
    else if (strcmp(args[i], "-rs")==0) {
      if ((++i)<argC) {
        FB->rfb_size = atol(args[i]);
      }
    }
    else if (strcmp(args[i], "-rl")==0) {
      if ((++i)<argC) {
        FB->rLim = atol(args[i]);
      }
    }
    else if (strcmp(args[i], "-mpr")==0) {
      if ((++i)<argC) {
        FB->maxP_r = atol(args[i]);
      }
    }
    else if (strcmp(args[i], "-mpa")==0) {
      if ((++i)<argC) {
        FB->maxP_a = atol(args[i]);
      }
    } else if (strcmp(args[i], "-lpbr")==0) {
      if ((++i)<argC) {
        FB->maxP_r = (1<<atoi(args[i]));
      }
    } else if (strcmp(args[i], "-lpba")==0) {
      if ((++i)<argC) {
        FB->maxP_a = (1<<atoi(args[i]));
      }
    } else if (strcmp(args[i], "-2p")==0) {
      twoP=1;
    } else if (strcmp(args[i], "-3p")==0) {
      threeP=1;
    }
  }
  if (ifname[0]==0) { 
    printf("%s\n%s\n", args[0], USAGE);
    return 1;
  }
  if (FB->maxP_r) {
    if (threeP) FB->maxLP = 3;
    else if (twoP) FB->maxLP = 2;
    else FB->maxLP=1;
  }
  if (FB->maxP_a) {
    if (threeP) FB->maxLPA = 3;
    else if (twoP) FB->maxLPA = 2;
    else FB->maxLPA=1;
  }

  if (!(ifp = fopen(ifname, "r"))) {
    fprintf(stderr, "Error opening %s for read!\n", ifname);
    return -1;
  }
  if (readPoly(ifp, FB)) {
    fprintf(stderr, "readPoly() reported some error!\n");
    exit(-1);
  }
  fclose(ifp);

#else
int createFB(nfs_fb_t *FB, char *ofname)
{ char     str[1024];
  int      i, status;
  time_t   start, stop;

#endif

  
  /*********************/
  /* Verify parameters */
  /*********************/
  if ((FB->rfb_size < MINRFBSIZE)||(FB->rfb_size > MAXRFBSIZE)) {
    if ((FB->rLim < MIN_RLIM) || (FB->rLim > MAX_RLIM)) {
      fprintf(stderr, "Error! Invalid RFB size or limit.\n");
      return 0;
    }
  }
  if ((FB->afb_size < MINAFBSIZE)||(FB->afb_size > MAXAFBSIZE)) {
    if ((FB->aLim < MIN_ALIM) || (FB->aLim > MAX_ALIM)) {
      fprintf(stderr, "Error! Invalid AFB size or limit.\n");
      return 0;
    }
  }
  msgLog("", "GGNFS-%s : makefb", GGNFS_VERSION);
  { mpz_t g;

    mpz_init_set(g, &FB->f->coef[0]);
    for (i=1; i<=FB->f->degree; i++)
      mpz_gcd(g, g, &FB->f->coef[i]);
    mpz_abs(g,g);
    if (mpz_cmp_ui(g, 1)) {
      msgLog("", "Error: Polynomial coefficients should have gcd of 1!");
      printf("Error: Polynomial coefficients should have gcd of 1, but it is:\n");
      mpz_out_str(stdout, 10, g);
      printf("\nProbably you can just divide this factor out and try again?!?\n");
      exit(-1);
    }
  }
  printf("Making factor base..."); fflush(stdout);
  time(&start);
  status = makefb(FB);
  time(&stop);
  printf("Elapsed time: %ld seconds.\n", stop-start);
  if (status) 
    fprintf(stderr, "\a\aThere were some errors creating the factor base!\n");
  else {
    saveFB(ofname, FB);
    fprintf(stderr, "Factor base sucessfully created.\n");

    msgLog("", "name: %s", FB->name);
    mpz_get_str(str, 10, FB->n);
    msgLog("", "n=%s (%ld digits)", str, strlen(str));
    for (i=0; i<=FB->f->degree; i++) {
      mpz_get_str(str, 10, &FB->f->coef[i]);
      msgLog("", "c%d: %s", i, str);
    }
    msgLog("", "RFBsize: %" PRId32 " (upto %" PRId32 ")", FB->rfb_size, 
           FB->rfb[2*(FB->rfb_size-1)]);
    msgLog("", "AFBsize: %" PRId32 " (upto %" PRId32 ")", FB->afb_size,
           FB->afb[2*(FB->afb_size-1)]);
    msgLog("", "maxNumLargeRatPrimes: %d", FB->maxLP);
    msgLog("", "maxLargeRatPrime: %" PRId32, FB->maxP_r);
    msgLog("", "maxNumLargeAlgPrimes: %d", FB->maxLPA);
    msgLog("", "maxLargeAlgPrime: %" PRId32, FB->maxP_a);
    


  }
  return status;
}
  
  

