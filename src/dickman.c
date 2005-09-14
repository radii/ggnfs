/*********************************************************/
/* dickman.c                                             */
/* Copyright 2004, Chris Monico, Andrei Belenko.         */
/* Portions (everything but dickmanStrong()) are         */
/*   Copyright (C) 2001 Jens Franke, also subect to GPL. */
/*********************************************************/
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

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ggnfs.h"


/*
  Following are fast dickmann function calculation and table generation code
  taken from Franke's gnfs-mmqpg implementation and modified a little bit to suit GGNFS needs
*/



u32 dtab_max = 0;
double dickm_maxarg = 0;
u32 dickm_nstep = 0;

double *dickmann_table = NULL;

/*
  Initializes table required to compute dickmann function
  Should be called only once but before any call to dickman()
  Returns nonzero value on success or zero otherwise.
*/
int dickmann_init()
{ u32 last_arg= 31;
  u32 n_steps= 32; 
  double eps= 3e-3; /* Results in about 59K. */
  double *dtab, *old_dtab=NULL, x, y, nq;
  double f,maxf;
  u32 i,j;
  u32 n_rounds= 0;

  do {
    nq= 1.0/(2*n_steps);
    dtab= (double*)malloc((n_steps*last_arg + 1) * sizeof(double));
    if( dtab == NULL )
      return 0;
  
    for(i= 0;i<=n_steps;i++) 
      dtab[i]= 1.0;
    for(;i<=2*n_steps;i++) {
      x= i;
      x/= n_steps;
      dtab[i]= 1-log(x);
    }
    for(;i<=last_arg*n_steps;i++) {
      x= i;
      x/= n_steps;

      j= i-n_steps;
      y= dtab[j]/2;
    
      for(j++;j<i;j++)
        y+= dtab[j];
    
      y+= dtab[i-n_steps]/(x*12*n_steps);
      y-= dtab[i-2*n_steps]/((x-1)*12*n_steps);
      y/= n_steps;
      dtab[i]= y/(x-nq);
    }
    if(old_dtab) {
      for(i= 0,j= 0,maxf= 0;i<=n_steps*last_arg;i+= 2,j++) {
        f= fabs(dtab[i]-old_dtab[j])/dtab[i];
        if(maxf<f)
          maxf= f;
        if( (i > 0) && (i < n_steps*last_arg) ) {
          f= fabs((dtab[i+1]+dtab[i-1])/2-dtab[i])/dtab[i];
          if(maxf<f)
            maxf= f;
        }
      }
      if(maxf<eps) {
        //printf("Konvergenz nach %u Runden\n",n_rounds);
        
        // yeah, I've studied deutch when I was 13 :-) --A.B.
        //printf("Convergence at round %u\n", n_rounds);
    
        dickm_nstep = n_steps;
        dickm_maxarg = (double)last_arg;
        dtab_max = n_steps*(last_arg-2);
        
        if (!(dickmann_table = malloc( (dtab_max + 1)*sizeof(double) )))
          return 0;
        
        j = 0;
        dickmann_table[j++] = dtab[2*n_steps];
        for(i= 2*n_steps+1;i<=n_steps*last_arg;i++) 
          dickmann_table[j++] = dtab[i];
        //free(old_dtab);
        //free(dtab);
        return 1;
      }
      free(old_dtab);
    }
    old_dtab= dtab;
    n_steps*= 2;
    n_rounds++;
  } while(1);
  return 0;
}

int dickmann_cleanup()
{
  if( dickmann_table )
    free(dickmann_table);

  return 1;
}

/*
  This is function actually computes dickmann's function
*/
double dickman(double x)
{ double h, ii;
  u32    i;
  
  if( dickmann_table == NULL )
    dickmann_init();
  
  if(x<=1.0) return 1.0;
  if(x<=2.0) return 1-log(x);
  x=x-2.0;
  h=modf(x*dickm_nstep,&ii);
  assert(ii <= UINT_MAX);
  i=(u32)ii;
  
  if(i>=dtab_max) {
    if(x==dtab_max) {
      return dickmann_table[i];
    }
    printf("\n!!!dickman(%.2f) argument is out of range!!!\n", x);
  }

  return h*dickmann_table[i+1]+(1-h)*dickmann_table[i];
}



/* Code to compute values of Dickman's rho function at double precision.
   We use the technique described in Bach, Peralta, ``Asymptotic semismoothness
   probabilities,'', Math.Comp. 65,216. There, though, they attribute this
   method to Patterson and Rumsey.

   Although the above function is now the primary one, we should keep
   this one around in case we need more accurate numbers at a later stage
   of polynomial selection.
*/
#define DICKMAN_MAXTERMS 55
double dickmanStrong(double x, int numTerms)
{ int    K, k, i, j, nt;
  double xi, c[DICKMAN_MAXTERMS], newC[DICKMAN_MAXTERMS], sum, xiPow;

  nt = MAX(MIN(DICKMAN_MAXTERMS, numTerms), 5);
  K = (int)ceil(x);
  xi = K-x;
  /* Initialize */
  if (K >= 2) {
    c[0] = 1.0 - M_LN2;
    for (i=1; i<nt; i++)
      c[i] = 1.0/((double)i*pow(2.0, (double)i));
    k=2;
  } else {
    c[0]=1.0;
    for (i=1; i<nt; i++)
      c[i] = 0.0;
    k=1;
  }
  /* It's tempting to factor the k^{-i} out of the summation,
     but don't do it --- it will cause (probably severe) rounding
     errors! If speed becomes a real concern, we can build some tables
     with the Taylor coefficients upto some reasonable bound.
  */
  while (k < K) {
    for (i=1; i<nt; i++) {
      sum = 0.0;
      for (j=0; j<i; j++) 
        sum += c[j]/((double)i*pow((double)(k+1), (double)i-j));
      newC[i] = sum;
    }
    /* Finally, compute newC[0] and we're done. */
    sum=0.0;
    for (j=1; j<nt; j++)
      sum += newC[j]/(double)(j+1);
    newC[0] = sum/(double)(k);
    memcpy(c, newC, nt*sizeof(double));
    k++;
  }
  /* Now, we can compute rho(x) = rho(K - xi) = sum[ c[i] * \xi^{i} ]. */
  /* Note: there is an obvious typo in this eq. in Bach & Peralta:     */
  /* it wrongly says \xi^k.                                            */
  sum=0.0;
  xiPow=1.0;
  for (i=0; i<nt; i++) {
    sum += c[i]*xiPow;
    xiPow *= xi;
  }
  return sum;
}


#ifdef _TEST_MAIN
int main(int argC, char *args[])
{ double x, d1, d2;
  int    i;
  

}
#endif
