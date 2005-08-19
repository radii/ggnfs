/* emc.c -- Integer factorization using the Elliptic Curve Method
   See http://www.loria.fr/~zimmerma/records/ecmnet.html

  Copyright (C) 1998 Paul Zimmermann, INRIA Lorraine, zimmerma@loria.fr
  See http://www.loria.fr/~zimmerma/records/ecmnet.html

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.

Changes with respect to 2a:
- use base-2 division for factors of 2^k+/-1
- added division using pre-inversion (macro PREINVERT): no saving
- now prints a warning for probable prime input
- now checks if factors found are composite or not
- now checks for prime powers
Changes with respect to 2b:
- saved a factor of two in step 2 initialization, and a factor of two
  in memory needed by step 2
- changed B2 and m in step 2 to be double's --> no overflow any more,
  even on 32-bit machines.
- fixed bug for multiple-line input (thanks to Torbjorn and P. Leyland)
Changes with respect to 2c:
- added LARGE macro for large input, like Fermat numbers, to disable
  printing of input number, values of A and x
- now does primality and perfect-power tests only with CHECK=1 (default 0)
Changes with respect to 2d:
- no normalization by default in step 1
- no gcd after each prime in step 1 (only one at the end of step 1),
   unless GROUPORDER is defined (useful for computing group order)
- no gcd after each giant-step in phase 2 (only one at the end), unless
   GROUPORDER is defined
- in step 2, if NOGCD if defined, does no gcd (3 times more multiplications)
  ==> useful for Fermat numbers
- included changes from Alex Stuebinger for ANSIONLY compilers (08/04/98)
- LARGE=1 by default (does not print A=..., x=...)
- removed compare2 (always use base-2 division for divisors of 2^n+/-1)
Changes with respect to 2f:
- Change several variables from type double to just int.
- Change sieve table (pr) to have values 0 and 1, not '0' and '1'.
- Change pr to unsigned char (faster on alpha, sparc, etc).
- Avoid repeated %q in step2 (as well as step1).
- Sieve only on odd i's in pr[i] (step 1) and even i's in step 2
- Avoid expensive %q computations in step 2 (now only at initialization)
Changes with respect to 2g:
- combine k prime powers in step 1 and replace k muls by one inverse
To do:
- use Peter's trick to replace k gcdexts by one gcdext and (k-1) muls
  (both in step 1 and 2)
Changes with respect to 2h:
- removed PREINVERT stuff
- new add3 procedure: removed all mpz_set in multiply
- implemented Peter Montgomery's PRAC algorithm [4]
Changes with respect to 2i:
- implemented Peter Montgomery's MODMULN algorithm
Changes with respect to 2ia:
- completely rewritten step 2
- removed NOGCD stuff
Changes with respect to version 3a:
- added GNU General Public License
- removed REDC stuff
Changes in version 4 (April 1999):
- uses fast multipoint polynomial evaluation in step 2,
  together with O(n^1.59) Karatsuba multiplication
- uses k passes in step 2: cost [30+12*(k-1)]*M(n/2)
Changes in version 4a:
- cost of recip reduced to 9/2*M(n/2) instead of 6*M(n/2)
  i.e. cost of step 2 is now [57/2+12*(k-1)]*M(n/2)
- improved help
Changes in version 4b (November 1999):
- replaced use of reciprocals in polyeval by Burnickel/Ziegler algorithm
  i.e. cost of step 2 is now (4k+2)*K(n)
Changes in version 4c (December 1999):
- in step 2, call mpz_mod only *after* Karatsuba (i.e. 2*k calls to mpz_mod
  instead of k^1.59)

This version uses Montgomery's form (8) from [2] which avoids gcds:

        b*y^2*z = x^3 + a*x^2*z + x*z^2

References:
[1] "Speeding the Pollard and Elliptic Curve Methods of Factorization", by 
   Peter Montgomery, Math. of Comp. 48 (177), pages 243-264, January 1987.
[2] "Factorization of the tenth and eleventh Fermat numbers", by Richard Brent,
ftp://nimbus.anu.edu.au/pub/Brent/rpb161tr.dvi.gz
[3] Torbjorn Granlund, Peter L. Montgomery: Division by Invariant Integers 
    using Multiplication. PLDI 1994: 61-72, SIGPLAN Notices 29(6) (June 1994)
[4] "Evaluating recurrences of form X_{m+n} = f(X_m,X_n,X_{m-n}) via Lucas
     chains", Peter Montgomery, ftp.cwi.nl:/pub/pmontgom/Lucas.ps.gz

Examples (log and timing lines omitted):

% echo 137703491 | ecm 100 6
********** Factor found in step 1: 17389

(if compiled with -DGROUPORDER)
% echo 137703491 | ecm 100 13
********** Factor found in step 2: 7919

(bug found by T. Granlund in 1st version of ecm2f)
% echo 17061648125571273329563156588435816942778260706938821014533 | ecm 174000 585928442
********** Factor found in step 2: 4562371492227327125110177

From [2], page 15 (factorization of 55^126+1):
% echo 5394204444759808120647321820789847518754252780933425517607611172590240019087317088600360602042567541009369753816111824690753627535877960715703346991252857 | ecm 345551 805816989
********** Factor found in step 1: 25233450176615986500234063824208915571213

% ecm 314263 14152267 4677853 < F10.cofactor
Input number is 607820568181834328745927047401406785398975700821911559763928675076909152806525747797078707978021962487854849079350770968904705424125269800765765006449689562590686195386366153585734177565092347016126765195631310982002631912943551551593959032889971392442015624176361633631364310142874363629569
********** Factor found in step 2: 4659775785220018543264560743076778192897

# first Cunningham factor found by GMP-ECM (06 Dec 1997)
% echo 449590253344339769860648131841615148645295989319968106906219761704350259884936939123964073775456979170209297434164627098624602597663490109944575251386017 | ecm 1000000 63844855
********** Factor found in step 2: 241421225374647262615077397

# p48 found by Richard Brent on October 9, 1997
% echo 3923385745693995079670229419275984584311007321932374190635656246740175165573932140787529348954892963218868359081838772941945556717 | ecm 141667 876329474 150814537
********** Factor found in step 2: 662926550178509475639682769961460088456141816377

# p45 found by Richard Brent on October 24, 1997
% echo 89101594496537524661600025466303491594098940711325290746374420963129505171895306244425914080753573576861992127359576789001 | ecm 325001 877655087 1032299
********** Factor found in step 2: 122213491239590733375594767461662771175707001
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#include "if.h"

#if defined (__MINGW32__) || defined (MINGW32)
#ifndef _MSC_VER
#define _MSC_VER
#endif
#endif

#define PTR(x) ((x)->_mp_d)
/* #endif
*/
#include "prand.h"


#ifndef max
#define max(a,b) (((a)>(b)) ? (a) : (b))
#define min(a,b) (((a)<(b)) ? (a) : (b))
#endif


/* ANSI Prototypes */
extern int  ecm(__mpz_struct *,__mpz_struct *,__mpz_struct *,double,int,int);
extern int  step1(__mpz_struct *);
extern void initprimes(double,int );
extern void prac(unsigned int);
extern void add3(__mpz_struct *,__mpz_struct *,__mpz_struct *,__mpz_struct *,
		 __mpz_struct *,__mpz_struct *,__mpz_struct *,__mpz_struct *);
extern void duplicate(__mpz_struct *,__mpz_struct *,__mpz_struct *,__mpz_struct *);
extern int  step2(__mpz_struct *,__mpz_struct *,unsigned int,double,
		 __mpz_struct *,__mpz_struct *,int,unsigned int);
extern int  cputime(void);
extern int  isbase2(__mpz_struct *,double);
extern void mod2plus(__mpz_struct *,__mpz_struct *,__mpz_struct *);
extern void mod2minus(__mpz_struct *,__mpz_struct *,__mpz_struct *);
extern void mpz_mod_n(__mpz_struct *,__mpz_struct *,__mpz_struct *);
void (*mod)(mpz_t a,mpz_t b,mpz_t n)=NULL;
extern int  multiplyW(__mpz_struct *,__mpz_struct *,__mpz_struct *,
		     __mpz_struct *,__mpz_struct *,unsigned int,__mpz_struct *,
		     __mpz_struct *,__mpz_struct *,__mpz_struct *);
extern int  multiplyW2(__mpz_struct *,__mpz_struct *,__mpz_struct *,
		     __mpz_struct *,__mpz_struct *,__mpz_struct *,__mpz_struct *,
		     __mpz_struct *,__mpz_struct *,__mpz_struct *);
extern int  subW(__mpz_struct *,__mpz_struct *,__mpz_struct *,__mpz_struct *,
		__mpz_struct *,__mpz_struct *,__mpz_struct *,__mpz_struct *,
		__mpz_struct *,__mpz_struct *);
extern int  addWn(__mpz_struct *,mpz_t *,mpz_t *,
		 __mpz_struct *,mpz_t *,mpz_t *,int);
extern int  duplicateW(__mpz_struct *,__mpz_struct *,__mpz_struct *,
		      __mpz_struct *,__mpz_struct *,__mpz_struct *,
		      __mpz_struct *,__mpz_struct *,__mpz_struct *);
extern int  addW(__mpz_struct *,__mpz_struct *,__mpz_struct *,__mpz_struct *,
		__mpz_struct *,__mpz_struct *,__mpz_struct *,__mpz_struct *,
		__mpz_struct *,__mpz_struct *);
extern void polyeval(mpz_t *,mpz_t **,mpz_t *,unsigned int);
extern void buildF(mpz_t **,mpz_t *,unsigned int);
extern double default_B2(double);
extern void polymul(mpz_t*,mpz_t*,mpz_t*,unsigned int,mpz_t*);
extern void karatsuba(mpz_t*,mpz_t*,mpz_t*,unsigned int,mpz_t*);
extern void karatsuba1(mpz_t*,mpz_t*,mpz_t*,unsigned int,mpz_t*);
extern void buildG(mpz_t*,mpz_t*,unsigned int);
extern void polymulmod(mpz_t*,mpz_t*,mpz_t*,int);
extern int  initpoly(int,mpz_t*,mpz_t*,mpz_t,mpz_t,mpz_t,mpz_t,int,int);
extern void set_dickson(mpz_t*, int);
extern void RecursiveDivision(mpz_t *,mpz_t *,mpz_t *,unsigned int,mpz_t *);
extern void DivThreeLongHalvesByTwo(mpz_t *,mpz_t *,mpz_t *,
				    unsigned int, mpz_t *);

/********* global variables ***********/
static unsigned int bb, *prime, nbprimes, mul, gcdexts, base2=1, go=0;
static mp_size_t    sizen;
static int          B1, ispower2, use_dickson=0; 
static unsigned char *pr;
static mpz_t        a, b, n, initial_n, u, v, w, x, z, x1, z1, x2, z2;
static mpz_t        one,y,invn, x3, z3, x4, z4, *dick=NULL;
#ifdef MODMULN
static mp_limb_t    Nprim;
#endif
/**************************************/



/*******************************************************/
unsigned int nb_digits(mpz_t n)
{ unsigned int size; 
  char        *str;
 
  str = mpz_get_str(NULL,10,n);
  size = strlen(str);
  free(str);
  return size;
}

/*******************************************************/
void ecm_initGlobals()
{
  mpz_init(n);
  mpz_init(a); mpz_init(b); mpz_init(u); mpz_init(v);
  mpz_init(w); mpz_init(x); mpz_init(x1); mpz_init(z1); mpz_init(x2);
  mpz_init(z2); mpz_init(z); mpz_init(x3); mpz_init(z3);
  mpz_init(x4); mpz_init(z4); mpz_init(initial_n);
  mpz_init_set_ui(one,1); mpz_init(y); mpz_init(invn);
}
/*******************************************************/
void ecm_clearGlobals()
{
  mpz_clear(n);
  mpz_clear(a); mpz_clear(b); mpz_clear(u); mpz_clear(v);
  mpz_clear(w); mpz_clear(x); mpz_clear(x1); mpz_clear(z1); mpz_clear(x2);
  mpz_clear(z2); mpz_clear(z); mpz_clear(one); mpz_clear(y); mpz_clear(invn);
  mpz_clear(x3); mpz_clear(z3); mpz_clear(x4); mpz_clear(z4);
  mpz_clear(initial_n);
}

/************************************************************/
int bestD(double B2, int *fft_size, int *phi)
/************************************************************/
/* we must have 2D mod 6 = 0, B2<2*D*fft_size and           */
/* phi=(phi(2D)/2)<fft_size.                                */
/* The following optimized values come from Aiichi YAMASAKI,*/
/* Department of Mathematics, Kyoto University.             */
/************************************************************/
{
  if (B2<=336.0) { *fft_size=8; *phi=6; return(42/2); }
  else if (B2<=1440.0) { *fft_size=16; *phi=12; return(90/2); }
  else if (B2<=6720.0) { *fft_size=32; *phi=24; return(210/2); }
  else if (B2<=29568.0) { *fft_size=64; *phi=60; return(462/2); }
  else if (B2<=134400.0) { *fft_size=128; *phi=120; return(1050/2); }
  else if (B2<=591360.0) { *fft_size=256; *phi=240; return(2310/2); }
  else if (B2<=2365440.0) { *fft_size=512; *phi=480; return(4620/2); }
  else if (B2<=9461760.0) { *fft_size=1024; *phi=960; return(9240/2); }
  else if (B2<=39137280.0) { *fft_size=2048; *phi=2016; return(9555); }
  else if (B2<=160849920.0) { *fft_size=4096; *phi=3840; return(19635); }
  else if (B2<=648560640.0) { *fft_size=8192; *phi=8064; return(39585); }
  else if (B2<=2611445760.0) { *fft_size=16384; *phi=15840; return(79695); }
  else if (B2<=10824253440.0) { *fft_size=32768; *phi=31680; return(165165); }
  else if (B2<=45265059840.0) { *fft_size=65536; *phi=63360; return(345345); }
  else { printf("Error: too large B2\n"); exit(1); }
}

/*******************************************************/
int best_e(int K)
{ int e;

  if (K==256) e=2;
  else if (K==512) e=3;
  else if (K==1024) e=6;
  else if (K==2048) e=12;
  else if (K==4096) e=30;
  else if (K==8192) e=60;
  else if (K==16384) e=120;
  else e=1;
  return e;
}

/*******************************************************/
double default_B2 (double B1)
{ double c1, c2, B2, iter, oldB2, e; 
  int    K, phi, D;

  c1 = 9.0*B1/log(2.0); /* estimated number of modular mult. for step 1 */
  if (verbose) printf("estimated muls for step 1: %1.0f\n", c1);
  oldB2 = B2 = 100.0 * B1;
  /* the following gives the expected cost of step 2 with iter=6,
     assuming 2*D/phi(D) is about 2/log(B2/iter) */
  do {
    iter = 8.0; /* optimal value is iter=2 */
    D = bestD(B2/iter,&K,&phi);
    while (iter>0 && 2*(iter-1)*K*(double)D>=B2) iter--;
    e = (double) best_e(K);
    c2 = 6.0*e*((double)D/3.0) /* computation of nQx */
       + (4.0*iter+2.0)*pow((double)K, log(3.0)/log(2.0)) /* poly. stuff */
       + 6.0*e*iter*(double)K; /* computation of (2D)^m*Q */
    if (verbose) printf("estimated muls for B2=%u*%u*%u: %1.0f\n",
			(unsigned) iter, K-1, D, c2);
    if (c2<c1/2.0) { oldB2=B2; B2*=1.1; }
  } while (c2<c1/2.0);
  B2 = oldB2;
  return B2;
}

/********************************************************************/
int ecm(mpz_t p, mpz_t n, mpz_t s, double B2, int e, int iter)
/********************************************************************/
/* factors n and puts the result in p, s is the seed (0 -> random)  */
/* returns 0 iff no factor found                                    */
/* iter=0 means choice left to the program, otherwise imposed.      */
/********************************************************************/
{  unsigned int st,res,D,K,phi;

   mul=0; gcdexts=0; sizen=mpz_size(n); mpz_set(initial_n, n);
   mpz_set_ui(p,6); mpz_gcd(p,n,p); if (mpz_cmp(p,one)) return(1);
   /* now gcd(n,6)=1 */
   /* slower in step 1 than usual division for 2,568+ c120 (1.43)
      but faster for 2,671- c145 (1.40)
   */
   if ((base2 && (ispower2 = isbase2(n,1.4)))) {
     printf("recognized factor of 2^%d",(ispower2>0) ? ispower2 : -ispower2);
     if (ispower2>0) {
       mod=mod2plus;
       printf("+"); 
     }
     else {
       mod=mod2minus;
       printf("-");
     }
     printf("1, using special base-2 division\n");
     fflush(stdout);
   }
   else
#ifdef MODMULN
     { 
       mod=mpz_mod_n;
       mpz_set_ui(v,1); mpz_mul_2exp(v,v,mp_bits_per_limb); /* v=2^k */
       mpz_gcdext(z,u,NULL,n,v);
       /* z should be 1 since n is odd and v a power of 2 */
       if (mpz_cmp_ui(z,1)!=0) { fprintf(stderr,"gcd(n,R) is not 1\n"); exit(1);}
       mpz_neg(u,u); mpz_mod(u,u,v);
       Nprim=PTR(u)[0]; /* Nprim * n = -1 mod v=2^k, the word base */ 
     }
#else
   mod=(void*)mpz_mod;
#endif
   if (iter==0) iter=12; /* max number of iterations */
   D=bestD(B2/iter,&K,&phi);
   if (e==0) e=best_e(K);
   if (e==1) use_dickson=0; /* Dickson(1) = x^1 */
   if (use_dickson && dick==NULL) {
     dick=(mpz_t*) malloc((e/2+1)*sizeof(mpz_t));
     set_dickson(dick, e);
   }
   /* adjust iter if too large */
   if (B2>0.0) 
     while (2*(iter-1)*(K-1)*(double)D>=B2) 
       iter--;
   if (2*iter*(K-1)*(double)D<B2) 
     iter = ceil(B2/2.0/((double)K-1.0)/(double)D);
   printf("Using B1=%d, B2=%1.0f",B1,2*iter*(K-1)*(double)D);
   if (verbose) 
     printf(" (%u*%u*%u)",iter,K-1,2*D);
   if (use_dickson==0) 
     printf(", polynomial x^%u",e);
   else printf(", polynomial Dickson(%u)",e);
   if (mpz_cmp_ui(x,0)) /* start from given a and x instead of s */
     mpz_set_ui(z,1);
   else {
     /* generates a random starting point using (11) from [2], or take the 's' given */
     if (mpz_cmp_ui(s,0)) 
       mpz_set(u,s);
     else { /* generate a random sigma */
       mpz_set_ui(v,prand()); /* thanks to Conrad Curry, generates 31-bits */
       mpz_mod(u,v,n); 
     }
     printf(", sigma="); mpz_out_str(stdout,10,u);
     mpz_mul_ui(w,u,4); mpz_mod(v,w,n); /* v = (4*s) mod n */
     mpz_mul(x,u,u); mpz_sub_ui(w,x,5); mpz_mod(u,w,n); /* u = (s^2-5) mod n */
     mpz_mul(x,u,u); mpz_mul(w,x,u); mpz_mod(x,w,n); /* x = u^3 mod n */
     mpz_mul(z,v,v); mpz_mul(w,z,v); mpz_mod(z,w,n); /* z:=v^3 mod n */
     mpz_mul(b,x,v); mpz_mul_ui(w,b,4); mpz_mod(b,w,n); /* b = (4*x*v) mod n */
     mpz_sub(a,v,u); mpz_mul(w,a,a); mpz_mul(w,w,a); mpz_mod(w,w,n); /* w = (v-u)^3*/
     mpz_mul_ui(a,u,3); mpz_add(a,a,v); mpz_mul(w,w,a); mpz_mod(a,w,n);
     /* a = ((v-u)^3*(3*u+v)) mod n */
     mpz_gcdext(p,u,NULL,b,initial_n); gcdexts++; /* w = gcd(b,n) = u*b mod n */
     if (mpz_cmp(p,one)) 
       goto youpi;
     mpz_mul(a,a,u); mpz_sub_ui(a,a,2); 
     mpz_mod(a,a,initial_n); /* a = a/b-2 mod n */
   }
   printf("\n"); fflush(stdout);
   if (verbose || go) {
     printf("A="); mpz_out_str(stdout,10,a); printf("\n"); fflush(stdout);
   }
   mpz_add_ui(b,a,2);
   if (mpz_mod_ui(w,b,2)) mpz_add(b,b,initial_n); mpz_tdiv_q_2exp(b,b,1); /* b = b/2 */
   if (mpz_mod_ui(w,b,2)) mpz_add(b,b,initial_n); mpz_tdiv_q_2exp(b,b,1); /* b = b/2 */
   /* now b = (a+2)/4 mod n */
   mpz_gcdext(p,u,NULL,z,initial_n); gcdexts++; if (mpz_cmp(p,one)) goto youpi;
   mpz_mul(x,x,u); mpz_mod(x,x,initial_n);
   mpz_set_ui(z,1);
   if (verbose || go) {
     printf("starting point: x="); mpz_out_str(stdout,10,x); printf("\n");
     fflush(stdout);
   }
   /* Step 1 */
   st=cputime();
   res=step1(p);
   printf("Step 1 took %dms for %d muls, %d gcdexts\n",
		      cputime()-st,mul,gcdexts); fflush(stdout);
   if (res) {
     printf("           Factor found in step 1: "); mpz_out_str(stdout,10,p);
     printf("\n"); fflush(stdout); goto youpi;
   }
   mul=gcdexts=0;
   st=cputime();
   res=step2(p,n,B1,B2,x,a,(int)e,iter);
   printf("Step 2 took %dms for %d muls, %d gcdexts\n",
		      cputime()-st,mul,gcdexts); fflush(stdout);
   if (res) {
     printf("           Factor found in step 2: "); mpz_out_str(stdout,10,p);
     printf("\n"); fflush(stdout); goto youpi;
   }
   return(0);
 youpi:
   return(1);
}

#ifdef MODMULN
/* multiplies c by R^k modulo n where R=2^mp_bits_per_limb 
   n is supposed odd. Does not need to be efficient. */
void mod_mul2exp(c,n,k) mpz_t c,n; unsigned int k;
{
  mpz_mul_2exp(c,c,k*mp_bits_per_limb);
  mpz_mod(c,c,n);
}

/* divides c by R^k modulo n where R=2^mp_bits_per_limb
   n is supposed odd. Does not need to be efficient. */
void mod_div2exp(c,n,k) mpz_t c,n; unsigned int k;
{
  mpz_t invR,g,R;

  /* first computes the inverse of R mod n */
  mpz_init(invR); mpz_init(g); mpz_init(R);
  mpz_set_ui(R,1); mpz_mul_2exp(R,R,mp_bits_per_limb);
  mpz_gcdext(g,invR,NULL,R,n); /* g = 1 = invR*R mod n */
  while (k-->0) {
    mpz_mul(c,c,invR);
    mpz_mod(c,c,n);
  }
  mpz_clear(invR); mpz_clear(g); mpz_clear(R);
}
#endif

/****************************************************************/
int step1(mpz_t p)
/****************************************************************/
/* returns 0 iff no factor found, otherwise returns factor in p */
/****************************************************************/
{ unsigned int l, i, j, q, imax, lmax, pp;
	
#ifdef MODMULN
   if (ispower2==0) {
   /* multiply (x,z) by R^sizen */
   mod_mul2exp(x,n,sizen); mod_mul2exp(z,n,sizen);
   mod_mul2exp(b,n,sizen); /* for duplicate */
   _mpz_realloc(x,2*sizen+1); _mpz_realloc(z,2*sizen+1);
   _mpz_realloc(x1,2*sizen+1); _mpz_realloc(z1,2*sizen+1);
   _mpz_realloc(x2,2*sizen+1); _mpz_realloc(z2,2*sizen+1);
   _mpz_realloc(x3,2*sizen+1); _mpz_realloc(z3,2*sizen+1);
   _mpz_realloc(x4,2*sizen+1); _mpz_realloc(z4,2*sizen+1);
   _mpz_realloc(u,2*sizen+1); _mpz_realloc(v,2*sizen+1);
   _mpz_realloc(w,2*sizen+1); }
#endif
  /* treat the cases p=2 and p=3 separately */
 for (q=2;q<=B1;q*=2) duplicate(x,z,x,z);
  for (q=3;q<=B1;q*=3) { duplicate(x1,z1,x,z); add3(x,z,x,z,x1,z1,x,z); }
  lmax = B1/bb;
  for (l=0;l<=lmax;l++) {
    /* check range l*bb <= p < (l+1)*bb */
    if (l) { /* sieve primes, pr[i] corresponds to l*bb+i */
      for (i=0;i<bb;i++) pr[i]='1';
      for (j=1;j<=nbprimes;j++) {
	/* delete multiples of prime[j] */
	q=prime[j];
	i=(q-((l*bb)%q)) % q;
	for(;i<bb;i+=q) pr[i]='0';
      }
    }
    else {
      for (i=0;i<bb;i++) pr[i]='0';
      for (j=3;j<=nbprimes;j++) pr[prime[j]]='1';
    }
    imax = ((B1+1)<(l+1)*bb) ? B1+1-l*bb : bb;
    for (i=0;i<imax;i++)
      if (pr[i]=='1') {
	pp=l*bb+i; for (q=1;q<=B1/pp;q*=pp) prac(pp);
      }
  }
#ifdef MODMULN
   /* divide (x,z) by R^sizen before gcd */
   if (ispower2==0) { mod_div2exp(x,n,sizen); mod_div2exp(z,n,sizen); }
#endif
   mpz_gcdext(p,w,NULL,z,initial_n); gcdexts++; if (mpz_cmp(p,one)) return(1);
  /* normalizes z to 1 */
  mpz_mul(x,x,w); mpz_mod(x,x,initial_n); mpz_set_ui(z,1);
  return(0);
}

/*****************************************************/
void initprimes(B,b) double B; int b;
/*****************************************************/
/* initializes tables of primes up to max(sqrt(B),b) */
/*****************************************************/
{
  int i,j;

  i = (int)ceil(sqrt(B)+0.5);
  if (i>b) b=i;
  if (b%2) b++; /* ensures b is even for Step 1 */
  if (b<=(int)bb) return; /* already done */
  if (pr != NULL) free(pr);
  pr = (unsigned char*) malloc(b+1);
  /* compute primes up to b */
  for (i=2;i<=b;i++) pr[i]=1;
  j=2; do {
    for (i=j*j;i<=b;i+=j) pr[i]=0;
    while (pr[++j]==0);
  } while (j*j<=b);
  for (nbprimes=0,i=2;i<=b;i++) if (pr[i]!=0) nbprimes++;
  if (prime != NULL) free(prime);
  prime = (unsigned int*) malloc((nbprimes+1)*sizeof(int));
  for (j=0,i=2;i<=b;i++) if (pr[i]!=0) prime[++j]=i;
  bb=b;
}

#define ADD 6 /* number of multiplications in an addition */
#define DUP 5 /* number of multiplications in a duplicate */

#define START(d,v) ((d)/1.6180339887498948482-128.0+(v))


/*******************************************************/
unsigned int lucas_cost(n, v) unsigned int n; double v;
/*******************************************************/
/* returns the number of modular multiplications       */
/*******************************************************/
{ unsigned int c,d,e,r;

  d=n; r=(unsigned int)((double)d/v+0.5);
  if (r>=n) return(ADD*n);
  d=n-r; e=2*r-n; c=DUP+ADD; /* initial duplicate and final addition */
  while (d!=e) {
    if (d<e) { r=d; d=e; e=r; }
    if (4*d<=5*e && ((d+e)%3)==0) { /* condition 1 */
      r=(2*d-e)/3; e=(2*e-d)/3; d=r; c+=3*ADD; /* 3 additions */
    } else
    if (4*d<=5*e && (d-e)%6==0) { /* condition 2 */
      d=(d-e)/2; c+=ADD+DUP; /* one addition, one duplicate */
    } else
    if (d<=(4*e)) { /* condition 3 */
      d-=e; c+=ADD; /* one addition */
    } else
    if ((d+e)%2==0) { /* condition 4 */
      d=(d-e)/2; c+=ADD+DUP; /* one addition, one duplicate */
    } else
    if (d%2==0) { /* condition 5 */
      d/=2; c+=ADD+DUP; /* one addition, one duplicate */
    } else
    if (d%3==0) { /* condition 6 */
      d=d/3-e; c+=3*ADD+DUP; /* three additions, one duplicate */
    } else
    if ((d+e)%3==0) { /* condition 7 */
      d=(d-2*e)/3; c+=3*ADD+DUP; /* three additions, one duplicate */
    } else
    if ((d-e)%3==0) { /* condition 8 */
      d=(d-e)/3; c+=3*ADD+DUP; /* three additions, one duplicate */
    } else 
    if (e%2==0) { /* condition 9 */
      e/=2; c+=ADD+DUP; /* one addition, one duplicate */
    } else
      { printf("no condition qualifies for d=%u e=%u\n",d,e); exit(1); }
  }
  return(c);
}

#define NV 10

/***********************************************************************/
void prac(n) unsigned int n;
/***********************************************************************/
/* computes nP from P=(x:z) and puts the result in (x:z). Assumes n>2. */
/***********************************************************************/
{ unsigned int  d, e, r, i=0;
  __mpz_struct *xA, *zA, *xB, *zB, *xC, *zC, *xT, *zT, *xT2, *zT2, *t;
  static double v[10] = 
     {1.61803398875,1.72360679775,1.618347119656,1.617914406529,1.612429949509,
    1.632839806089,1.620181980807,1.580178728295,1.617214616534,1.38196601125};
   
   /* chooses the best value of v */
   for (d=0,r=ADD*n;d<NV;d++) {
     e=lucas_cost(n,v[d]);
     if (e<r) { r=e; i=d; }
   }
   d=n;
   r=(int)((double)d/v[i]+0.5);
   /* A=(x:z) B=(x1:z1) C=(x2:z2) T=T1=(x3:z3) T2=(x4:z4) */
   xA=x; zA=z; xB=x1; zB=z1; xC=x2; zC=z2; xT=x3; zT=z3; xT2=x4; zT2=z4;
   /* first iteration always begins by Condition 3, then a swap */
   d=n-r; e=2*r-n; 
   mpz_set(xB,xA); mpz_set(zB,zA); /* B=A */
   mpz_set(xC,xA); mpz_set(zC,zA); /* C=A */
   duplicate(xA,zA,xA,zA); /* A=2*A */
   while (d!=e) {
         if (d<e) { r=d; d=e; e=r; t=xA; xA=xB; xB=t; t=zA; zA=zB; zB=t; }
	 /* do the first line of Table 4 whose condition qualifies */
	 if (4*d<=5*e && ((d+e)%3)==0) { /* condition 1 */
	    r=(2*d-e)/3; e=(2*e-d)/3; d=r;
	    add3(xT,zT,xA,zA,xB,zB,xC,zC); /* T = f(A,B,C) */
	    add3(xT2,zT2,xT,zT,xA,zA,xB,zB); /* T2 = f(T,A,B) */
	    add3(xB,zB,xB,zB,xT,zT,xA,zA); /* B = f(B,T,A) */
	    t=xA; xA=xT2; xT2=t; t=zA; zA=zT2; zT2=t; /* swap A and T2 */
	  } else
	 if (4*d<=5*e && (d-e)%6==0) { /* condition 2 */
	   d=(d-e)/2; 
	   add3(xB,zB,xA,zA,xB,zB,xC,zC); /* B = f(A,B,C) */
	   duplicate(xA,zA,xA,zA); /* A = 2*A */
	 } else
	 if (d<=(4*e)) { /* condition 3 */
	   d-=e; 
	   add3(xT,zT,xB,zB,xA,zA,xC,zC); /* T = f(B,A,C) */
	   t=xB; xB=xT; xT=xC; xC=t;
	   t=zB; zB=zT; zT=zC; zC=t; /* circular permutation (B,T,C) */
	 } else
	 if ((d+e)%2==0) { /* condition 4 */
	   d=(d-e)/2; 
	   add3(xB,zB,xB,zB,xA,zA,xC,zC); /* B = f(B,A,C) */
	   duplicate(xA,zA,xA,zA); /* A = 2*A */
	 } else
	 if (d%2==0) { /* condition 5 */
	   d/=2; 
	   add3(xC,zC,xC,zC,xA,zA,xB,zB); /* C = f(C,A,B) */
	   duplicate(xA,zA,xA,zA); /* A = 2*A */
	 } else
	 if (d%3==0) { /* condition 6 */
	   d=d/3-e; 
	   duplicate(xT,zT,xA,zA); /* T1 = 2*A */
	   add3(xT2,zT2,xA,zA,xB,zB,xC,zC); /* T2 = f(A,B,C) */
	   add3(xA,zA,xT,zT,xA,zA,xA,zA); /* A = f(T1,A,A) */
	   add3(xT,zT,xT,zT,xT2,zT2,xC,zC); /* T1 = f(T1,T2,C) */
	   t=xC; xC=xB; xB=xT; xT=t;
	   t=zC; zC=zB; zB=zT; zT=t; /* circular permutation (C,B,T) */
	 } else
	 if ((d+e)%3==0) { /* condition 7 */
	   d=(d-2*e)/3; 
	   add3(xT,zT,xA,zA,xB,zB,xC,zC); /* T1 = f(A,B,C) */
	   add3(xB,zB,xT,zT,xA,zA,xB,zB); /* B = f(T1,A,B) */
	   duplicate(xT,zT,xA,zA); add3(xA,zA,xA,zA,xT,zT,xA,zA); /* A = 3*A */
	 } else
	 if ((d-e)%3==0) { /* condition 8 */
	   d=(d-e)/3; 
	   add3(xT,zT,xA,zA,xB,zB,xC,zC); /* T1 = f(A,B,C) */
	   add3(xC,zC,xC,zC,xA,zA,xB,zB); /* C = f(A,C,B) */
	   t=xB; xB=xT; xT=t; t=zB; zB=zT; zT=t; /* swap B and T */
	   duplicate(xT,zT,xA,zA);
	   add3(xA,zA,xA,zA,xT,zT,xA,zA); /* A = 3*A */
	 } else
	 if (e%2==0) { /* condition 9 */
	   e/=2; 
	   add3(xC,zC,xC,zC,xB,zB,xA,zA); /* C = f(C,B,A) */
	   duplicate(xB,zB,xB,zB); /* B = 2*B */
	 } else
	 { printf("no condition qualifies for d=%u e=%u\n",d,e); exit(1); }
       }
       add3(xA,zA,xA,zA,xB,zB,xC,zC);
#ifdef DEBUG
   if (d!=1) { printf("d!=1 at the end of PRAC\n"); exit(1); }
#endif
   if (x!=xA) { mpz_set(x,xA); mpz_set(z,zA); }
}

#define mpz_mulmod(a,b,c,n) { mpz_mul(a,b,c); mod(a,a,n); }

/* adds Q=(x2:z2) and R=(x1:z1) and puts the result in (x3:z3),
     using 5/6 mul, 6 add/sub and 6 mod. One assumes that Q-R=P or R-Q=P where P=(x:z).
     Uses the following global variables:
     - n : number to factor
     - x, z : coordinates of P
     - u, v, w : auxiliary variables
Modifies: x3, z3, u, v, w.
(x3,z3) may be identical to (x2,z2) and to (x,z)
*/
void add3(x3,z3,x2,z2,x1,z1,x,z) mpz_t x3,z3,x2,z2,x1,z1,x,z;
{
   mpz_sub(u,x2,z2); mpz_add(v,x1,z1);
   /*   u = x2-z2, v = x1+z1 */
   mpz_mulmod(u,u,v,n);
   /* u = (x2-z2)*(x1+z1) */
   mpz_add(w,x2,z2); mpz_sub(v,x1,z1);
   /* w = x2+z2, v = x1-z1 */
   mpz_mulmod(v,w,v,n);
   /* v = (x2+z2)*(x1-z1) */
   mpz_add(w,u,v); mpz_sub(v,u,v);
   /* w = 2*(x1*x2-z1*z2), v = 2*(x2*z1-x1*z2) */
   mpz_mulmod(w,w,w,n);
   /* w = 4*(x1*x2-z1*z2)^2 */
   mpz_mulmod(v,v,v,n);
   /* v = 4*(x2*z1-x1*z2)^2 */
   if (x==x3) {
     mpz_set(u,x);
     mpz_mulmod(w,w,z,n);
     mpz_mulmod(z3,u,v,n);
     mpz_set(x3,w);
   }
   else {
     mpz_mulmod(x3,w,z,n);
     /* x3 = 4*z*(x1*x2-z1*z2)^2 */
     mpz_mulmod(z3,x,v,n);
     /* z3 = 4*x*(x2*z1-x1*z2)^2 mod n */
   }
   mul += 6;
}

/* computes 2P=(x2:z2) from P=(x1:z1), with 5 mul, 4 add/sub, 5 mod.
     Uses the following global variables:
     - n : number to factor
     - b : (a+2)/4 mod n
     - u, v, w : auxiliary variables
Modifies: x2, z2, u, v, w
*/
void duplicate(x2,z2,x1,z1) mpz_t x2,z2,x1,z1;
{
   mpz_add(w,x1,z1); mpz_mulmod(u,w,w,n); /* u = (x1+z1)^2 mod n */
   mpz_sub(w,x1,z1); mpz_mulmod(v,w,w,n); /* v = (x1-z1)^2 mod n */
   mpz_mulmod(x2,u,v,n); /* x2 = (u*v) mod n */
   mpz_sub(w,u,v); /* w = u-v = 4*x1*z1 */
   mpz_mulmod(u,b,w,n);
   mpz_add(u,u,v); /* u = (v+b*w) mod n */
   mpz_mulmod(z2,w,u,n); /* z2 = (w*u) mod n */
   mul += 5;
}

unsigned int igcd(a,b) unsigned int a,b;
{
  unsigned int t;
  while (b!=0) {
    t=a; a=b; b=t%b;
  }
  return a;
}

/* returns k such that n=2^k */
int lg(n) int n;
{
  int k=0; 
  while (n>1) { 
    if (n%2!=0) { printf("Error: not a power of two\n"); exit(1); } 
    n/=2; k++; }
  return k;
}

void printpol(f,k,flag) mpz_t *f; unsigned int k,flag;
{
  int i;
  for (i=0;i<k;i++) {
    if (i>0 && mpz_cmp_ui(f[i],0)>=0) printf("+");
    mpz_out_str(stdout,10,f[i]);
    if (i>0) { printf("*x"); if (i>1) printf("^%d",i); }
  }
  if (flag) printf("+x^%d",k);
  printf("\n");
}

/* Step 2: improved standard continuation, cf [2] p. 7-8.
   - p: variable to put factor found
   - n: number to factor
   - B1,B2: bounds for step 1 and 2
   - x: x-coordinate of Q at the end of step 1 (z normalized to 1)
   - a: parameter from the curve b*y^2 = x^3 + a*x^2 + x used in step 1
   Returns 0 iff no factor found, otherwise puts factor in p.

   1. Computation of (6i+1)^e*Q: (D/3)*(6e) multiplications
   2. buildF(): 3^lg(K) multiplications
   4. Computation of (2mD)^e*Q: iter*K*(6e) multiplications
   5. buildG(): iter*3^lg(K) multiplications
   6. reducing G*H mod F: 3*(iter-1)*3^lg(K) multiplications
   7. polyeval(): 4*3^lg(K) multiplications
*/
int step2(p,n,B1,B2,x,a,e,iter) 
mpz_t p,n,x,a; unsigned int B1,iter; double B2; int e;
{
   mpz_t *T=NULL,**F=NULL,*nQx,*G,*lx,*ly,g,y,*tx,*ty; double m;
   int i,j,st,D,K,phi,residue,k; unsigned int res=0,nbit,extra,doit;
   unsigned long *Res=NULL, lgK;

   st=cputime();
   if (verbose || go) {
     printf("x="); mpz_out_str(stdout,10,x); printf("\n"); fflush(stdout);
   }
   /* faster for 2,951- c158 (1.82)
      but slower for 2,749- c123 (1.85) */
   if ((base2 && (ispower2 = isbase2(n,1.8)))) {
     if (ispower2>0) mod=mod2plus; else mod=mod2minus;
   }
   else mod=(void*)mpz_mod; /* can't use modmuln here */
   /* determines g,y such that g*y^2 = x^3 + a*x^2 + x */
   mpz_init(y); mpz_set_ui(y,1);
   mpz_init(g); mpz_add(g,x,a); mpz_mul(g,g,x); mpz_add_ui(g,g,1);
   mpz_mul(g,g,x); mod(g,g,n);
   /* change of coordinates x=g*X-a/3, y=g*Y to return to Weierstrass form */
   mpz_mul_ui(u,g,3); mpz_mul(u,u,g); mod(u,u,n);
   mpz_gcdext(p,v,NULL,u,initial_n);
   if (mpz_cmp_ui(p,1)!=0) return(1); /* v=1/(3g^2)*/
   mpz_mul_ui(x,x,3); mpz_add(x,x,a); mpz_mul(x,x,g);
   mpz_mul(x,x,v); mod(x,x,n); /* x = (x+a/3)/g = g*(3x+a)/(3g^2) */
   mpz_mul_ui(y,y,3); mpz_mul(y,y,g); mpz_mul(y,y,v); mod(y,y,n);
   mpz_mul(a,a,a); mpz_sub_ui(a,a,3); mpz_neg(a,a);
   mpz_mul(a,a,v); mod(a,a,n);
   D=bestD(B2/iter,&K,&phi);
   if (phi>=K) { printf("error: phi>=K\n"); exit(1); }
   initprimes(B2,2*D);
   /* with Q the point obtained by Step 1, we compute (6d+1)^e*Q for 0<=d<=D/3,
     (6i+1)^e*Q is stored in nQ[i], and 2DQ is stored in nQ[D/3+1]
     George Woltman noticed we don't need 2*d*Q for D/2 < d < D */
   nQx = (mpz_t*) malloc(K*sizeof(mpz_t));
   for (i=0;i<K;i++) mpz_init(nQx[i]);
   /* G[k] stores (2(k+1)D)^e*Q */
   G = (mpz_t*) malloc((K+1)*sizeof(mpz_t));
   for (k=0;k<=K;k++) mpz_init(G[k]);
   lx = (mpz_t*) malloc(2*(e+1)*sizeof(mpz_t));
   ly = (mpz_t*) malloc(2*(e+1)*sizeof(mpz_t));
   tx = lx+(e+1); ty = ly+(e+1);
   for (j=0;j<2*(e+1);j++) { mpz_init(lx[j]); mpz_init(ly[j]); }
   if ((res=initpoly(e,lx,ly,x,y,p,a,6,1))) goto youpi3;
   mpz_set(nQx[0], lx[0]);
   extra = K-1-phi; /* extra places to fill */
   j = D/3-phi; /* non prime residues mod 2D */
   if (go) { 
     Res = (unsigned long*) malloc(K*sizeof(unsigned int));
     Res[0]=1;
   }
   st=cputime();
   for (i=7,residue=1;i<=2*D;i+=6) {
     res=addWn(p,lx,ly,n,tx,ty,e); if (res) goto youpi3;
     if (igcd(i,D)==1 || j<=extra) {
       if (go) Res[residue]=i;
       mpz_set(nQx[residue++],lx[0]);
       if (residue==K-1) break;
     }
     else j--;
   }
   if (verbose) printf("computation of (6i+1)^%u Q for i=0 to %u took %dms\n", 
		       e, D/3-1, cputime()-st);
   if (residue>=K) { printf("error: residue>=K\n"); exit(1); }
   mpz_set_ui(nQx[K-1], 0);

   lgK = lg(K);
   F = (mpz_t**) malloc((lgK+1)*sizeof(mpz_t*));
   F[0] = nQx;
   T = (mpz_t*) malloc(5*K*sizeof(mpz_t)); /* for Karatsuba */
   for (i=0;i<5*K;i++) mpz_init(T[i]);
   if (!go) buildF(F,T,K);
   /* now F[i] contains polynomials of degree 2^i */

   /* Q <- (2D)^e*Q */
   if ((res=initpoly(e,lx,ly,x,y,p,a,2*D,2*D))) goto youpi2;
   /* now (lx[0],ly[0]) is (2D)^e*Q */
   for (nbit=0,doit=0,m=4.0*D;nbit<iter;nbit++) {
     st = cputime();
     for (k=0;k<K-1;m+=2.0*(double)D) {
       res = (e==1 && m==4.0*D) /* only case where lx[j]=lx[j+1] */
	 ? duplicateW(p,lx[0],ly[0],lx[0],ly[0],n,a,tx[0],ty[0])
	 : addWn(p,lx,ly,n,tx,ty,e);
       if (res) goto youpi2;
       /* now (lx[0],ly[0]) is m^e*Q */
       if (doit) mpz_set(G[k++],lx[0]); 
       else doit = (m+2.0*(double)D>(double)B1);
     }
     if (verbose)
       printf("computation of (mD)^%u Q for m=%1.0f to %1.0f took %dms\n",
                       e, m/D, m/D+2.0*K, cputime()-st);
     if (go) {
       printf("checking range %1.0f..%1.0f\n",m-(K-1)*2.0*(double)D,
	      m);
       /* tries to find a factor between nQx[0..K-2] and G[0..K-2] */
       for (i=0;i<K-1;i++) {
	 mpz_set_ui(g, 1);
	 for (k=0;k<K-1;k++) {
	   mpz_sub(w, nQx[i], G[k]); mpz_mul(g, g, w); mpz_mod(g, g, n);
	 }
	 mpz_gcd(p, g, initial_n); 
	 if (mpz_cmp_ui(p, 1)) { /* match between nQx[i] and some G[0..K-2] */
	   for (k=0;k<K-1;k++) {
	     mpz_sub(w, nQx[i], G[k]);
	     mpz_gcd(p, w, initial_n);
	     if (mpz_cmp_ui(p, 1)) {
	       printf("largest group order factor divides %1.0f^%u-%lu^%u\n",
		      m+(k-K+1)*2.0*(double)D,2*e,Res[i],2*e);
	       exit(1);
	     }
	   }
	 }
       }
     }
     else {
       mpz_set_ui(G[K-1],0);
       buildG(G,T+K,K);
       if (iter>1) { 
	 if (nbit==0) { /* copies  into T */
	   for (k=0;k<K;k++) mpz_set(T[k],G[k]);
	 }
	 else polymulmod(G,F[lgK],T,K);
       }
     }
   }
   if (iter>1) for (k=0;k<K;k++) mpz_set(G[k], T[k]);
   if (!go) {
     polyeval(G,F,T,K);
     /* now G[0]..G[K-1] contains the values of G(nQx[0])..G(nQx[K-1]) */
     mpz_set(g,G[0]);
     for (i=1;i<K;i++) mpz_mulmod(g,g,G[i],n);
     mpz_gcd(p,g,initial_n); if (mpz_cmp(p,one)) res=1;
   }
 youpi2:
   if (!go) for (i=1;i<lgK+1;i++) {
     for (j=0;j<=K;j++) mpz_clear(F[i][j]);
     free(F[i]); }
   free(F);
   for (i=0;i<5*K;i++) mpz_clear(T[i]); free(T);
 youpi3: 
   mpz_clear(g); /* thanks to Paul Leyland */
   mpz_clear(y);
   for (i=0;i<K;i++) mpz_clear(nQx[i]); free(nQx);
   for (k=0;k<=K;k++) mpz_clear(G[k]); free(G);
   for (i=0;i<2*(e+1);i++) { mpz_clear(lx[i]); mpz_clear(ly[i]); }
   free(lx); free(ly); if (go) free(Res);
   return(res);
 }

/* T <- G*T mod F */
void polymulmod(G,F,T,K) mpz_t *G,*F,*T; int K;
{
  int k; unsigned int st=0;

  if (verbose) st=cputime();
  /* assumes previous remainder is in T[0]..T[K-1] */
  karatsuba(T+K, T, G, K, T+3*K); /* T[K]..T[3K-2] <- T*G */
  mpz_set_ui(T[3*K-1], 0); /* as RecursiveDivision expects 2K coefficients */
  RecursiveDivision(T, T+K, F, K, T+3*K);
  for (k=0;k<K;k++) mpz_set(T[k], T[K+k]);
  if (verbose) printf("Reducing g*h mod f took %dms\n",cputime()-st);
}

/* Return user CPU time measured in milliseconds. Thanks to Torbjorn. */
#if defined (ANSIONLY) || defined (USG) || defined (__SVR4) || defined (_UNICOS) || defined(__hpux) || defined(_MSC_VER) || defined(__ppc__)
#include <time.h>

int cputime()
{
  return (int) ((double) clock () * 1000 / CLOCKS_PER_SEC);
}
#else
#include <sys/types.h>
#include <sys/resource.h>

int cputime ()
{ struct rusage rus;

  getrusage (0, &rus);
  return rus.ru_utime.tv_sec * 1000 + rus.ru_utime.tv_usec / 1000;
}
#endif

/**************************************************************/
int isbase2(mpz_t n, double threshold)
/**************************************************************/
/* returns +/-k if n is a factor of N=2^k+/-1 with            */
/* N<=n^threshold, 0 otherwise                                */
/**************************************************************/
{ unsigned int k,lo; 
  int          res=0; 
  mpz_t        u,w;

  mpz_init(u); mpz_init(w);
  lo=mpz_sizeinbase(n,2)-1;
  mpz_set_ui(u,1); mpz_mul_2exp(u,u,2*lo);
  mpz_mod(w,u,n); /* 2^(2lo) mod n = +/-2^(2lo-k) if m*n = 2^k+/-1 */
  k = mpz_sizeinbase(w,2)-1;
  /* try w = 2^k */
  mpz_set_ui(u,1); mpz_mul_2exp(u,u,k);
  if (mpz_cmp(w,u)==0) res=k-2*lo;
  else {
    /* try w = -2^k */
    mpz_neg(w,w); mpz_mod(w,w,n);
    k = mpz_sizeinbase(w,2)-1;
    mpz_set_ui(u,1); mpz_mul_2exp(u,u,k);
    if (mpz_cmp(w,u)==0) res=2*lo-k;
  }
  mpz_clear(u); mpz_clear(w);
  if (abs(res)>(int)(threshold*lo)) res=0;
  return(res);
}


void mod2plus(mpz_t a, mpz_t b, mpz_t n) 
/* N = 2^ispower2 + 1 */
{
  /* 2^k = -1 */
    mpz_tdiv_r_2exp(y,b,ispower2);
    mpz_tdiv_q_2exp(a,b,ispower2);
    mpz_sub(a,y,a);
    mpz_mod(a,a,n);
}

void mod2minus(mpz_t a, mpz_t b, mpz_t n)
/* N = 2^k - 1, ispower2<0 */
{
  /* 2^k = 1 */
    mpz_tdiv_r_2exp(y,b,-ispower2);
    mpz_tdiv_q_2exp(a,b,-ispower2);
    mpz_add(a,y,a);
    mpz_mod(a,a,n);
}

#ifdef MODMULN
static void mpn_incr(mp_ptr p, mp_limb_t incr)
{ mp_limb_t x;

  x = *p + incr;
  *p++ = x;
  if (x >= incr)
    return;
  while (++(*(p++)) == 0)
    ;
}

/*********************************************************/
void mpz_mod_n (c, a, n) mpz_t c,a,n;
/*********************************************************/
/* Computes c/R^nn mod n, where n are nn limbs           */
/* and c has space for size(c)+1 limbs.  n must be odd.  */
/*********************************************************/
{ mp_ptr    cp=PTR(c), np=PTR(n);
  mp_limb_t cy;
  mp_limb_t q;
  size_t    j,nn=sizen;

  for (j=ABS(SIZ(c));j<=2*nn;j++) cp[j] = 0;
  for (j = 0; j < nn; j++) {
    q = cp[0] * Nprim;
    cy = mpn_addmul_1 (cp, np, nn, q);
    mpn_incr (cp + nn, cy);
    cp++;
  }
  cp -= nn;
  if (cp[2*nn]) {
    cy = cp[2*nn] - mpn_sub_n (cp, cp + nn, np, nn);
    while (cy) cy -= mpn_sub_n (cp, cp, np, nn);
  }
  else 
    MPN_COPY (cp, cp + nn, nn);
  MPN_NORMALIZE (cp, nn);
  SIZ(c) = SIZ(c) < 0 ? -nn : nn;
}
#endif

/* a <- b^e */
#define poly1 mpz_ui_pow_ui

/* a = Dickson(e)(b) */
void poly2(a, b, e) mpz_t a, b; unsigned long int e;
{
  int i;

  if (use_dickson==0) mpz_pow_ui(a, b, e);
  else {
    mpz_set_ui(a,1);
    for (i=e/2-1;i>=0;i--) {
      mpz_mul(a,a,b); mpz_mul(a,a,b); mpz_add(a,a,dick[i]);
    }
    if (e==3) mpz_mul(a,a,b);
  }
}

/* puts in l[i], i=0..k the coefficients of Dickson(k,1,x) 
   if k is odd, l[0] contains coeff(1), l[1] contains coeff(3), ...
   if k is even l[0] contains coeff(0), l[1] contains coeff(2), ...
*/
void set_dickson(l, k) mpz_t *l; int k;
{
  int i,j; mpz_t *t;
  
  t = (mpz_t*) malloc((k/2+1)*sizeof(mpz_t));
  for (i=0;i<=k/2;i++) { mpz_init(t[i]); mpz_init(l[i]); }
  mpz_set_ui(l[0], 1); /* l = Dickson(1,1,x) = x */
  mpz_set_si(t[0], -2); mpz_set_ui(t[1], 1);
  for (i=3;i<=k;i++) {
    if (i%2) { /* puts result in l */
      for (j=0;j<i/2;j++) mpz_sub(l[j], t[j], l[j]);
      mpz_set_ui(l[i/2], 1);
    }
    else {
      mpz_neg(t[0], t[0]);
      for (j=1;j<i/2;j++) mpz_sub(t[j], l[j-1], t[j]);
      mpz_set_ui(t[i/2], 1);
    }
  }
  if (k%2==0) {
    for (i=0;i<=k/2;i++) mpz_set(l[i], t[i]);
  }
  for (i=0;i<=k/2;i++) mpz_clear(t[i]); free(t);
}

/* sets (lx[i]:ly[i]) to P_e(m*i+i0)*(x:y) for 0 <= i <= e */
int initpoly(e,lx,ly,x,y,p,a,m,i0) int e,m,i0; mpz_t *lx,*ly,x,y,p,a;
{
  int i,j,k,res=0,st,mul0,gcdext0; mpz_t tx,ty,*ll,q,r;

  st=cputime(); mul0=mul; gcdext0=gcdexts;
  mpz_init(tx); mpz_init(ty); mpz_init(q); mpz_init(r);
  ll = (mpz_t*) malloc((e+2)*sizeof(mpz_t));
  for (j=0;j<=e+1;j++) mpz_init(ll[j]);
  mpz_set_ui(tx, i0);
  for (i=0;i<=e;i++) {
    poly2(ll[i], tx, e); /* tx = m*i+i0 */
    mpz_add_ui(tx, tx, m);
  }
  for (i=1;i<=e;i++)
     for (j=e;j>=i;j--) mpz_sub(ll[j], ll[j], ll[j-1]);
  for (j=0;j<=e&&res==0;j++) {
    if (j>=1) {
      mpz_set(r, ll[j]); mpz_set_ui(lx[j], 0);
      for (k=j-1;k>=0&&res==0;k--) 
	if (mpz_cmp(r, ll[k])>=0) {
	/* compute difference wrt ll[j-1] */
	mpz_tdiv_qr(q, r, r, ll[k]);
        /* r = q*ll[k]+r' */
        res = multiplyW2(p,lx[j+1],ly[j+1],lx[k],ly[k],q,n,a,tx,ty);
	if (mpz_cmp_ui(lx[j], 0))
	  res = res || addW(p,lx[j],ly[j],lx[j+1],ly[j+1],lx[j],ly[j],n,tx,ty);
	else { mpz_set(lx[j], lx[j+1]); mpz_set(ly[j], ly[j+1]); }
      }
      if (mpz_cmp_ui(r, 0)) { 
	res = multiplyW2(p,lx[j+1],ly[j+1],x,y,r,n,a,tx,ty);
	res = res || addW(p,lx[j],ly[j],lx[j+1],ly[j+1],lx[j],ly[j],n,tx,ty);
      }
    }
    else res=multiplyW2(p,lx[j],ly[j],x,y,ll[j],n,a,tx,ty);
  }
  mpz_clear(tx); mpz_clear(ty); mpz_clear(q); mpz_clear(r);
  for (j=0;j<=e+1;j++) mpz_clear(ll[j]); free(ll);
  if (verbose) printf("initialization of table of differences took %dms, %d muls and %d gcdexts\n",cputime()-st,mul-mul0,gcdexts-gcdext0);
  return res;
}

/* (x1:y1) <- q*(x:y) where q is a small integer */
int multiplyW(p,x1,y1,x,y,q,n,a,u,v) mpz_t p,x1,y1,x,y,n,a,u,v; unsigned int q;
{
  unsigned int j,r,restore; mpz_t x2,y2;
  restore=(x1==x);
  if (restore) { mpz_init(x2); mpz_init(y2); x1=x2; y1=y2; }
  for (r=q,j=1;r!=1;r/=2,j<<=1);
  j >>= 1; 
  r=duplicateW(p,x1,y1,x,y,n,a,u,v);
  if (r) return(r);
  if (q&j) r=addW(p,x1,y1,x1,y1,x,y,n,u,v);
  if (r) return(r);
  j >>= 1;
  while (j!=0) {
    if (duplicateW(p,x1,y1,x1,y1,n,a,u,v)) return(1);
    if (q&j) if (addW(p,x1,y1,x1,y1,x,y,n,u,v)) return(1);
    j >>= 1;
  }
  if (restore) {   mpz_set(x,x1); mpz_set(y,y1);
		   mpz_clear(x2); mpz_clear(y2); }
  return(0);
}

#define getbit(x,i) (PTR(x)[i/mp_bits_per_limb] & ((mp_limb_t)1<<(i%mp_bits_per_limb)))

/* (x1:y1) <- q*(x:y) where q is a large integer */
int multiplyW2(p,x1,y1,x,y,q,n,a,u,v) mpz_t p,x1,y1,x,y,n,a,u,v,q;
{
  int j,neg=0;

  if (mpz_cmp_ui(q,0)<0) { neg=1; mpz_neg(q,q); }
  if (mpz_cmp_ui(q,1)==0) {
    mpz_set(x1,x); mpz_set(y1,y);
    goto exit_multiplyW2;
  }
  j=mpz_sizeinbase(q,2)-2;
  if (duplicateW(p,x1,y1,x,y,n,a,u,v)) return 1;
  if (getbit(q,j) && addW(p,x1,y1,x1,y1,x,y,n,u,v)) return 1;
  j--;
  while (j>=0) {
    if (duplicateW(p,x1,y1,x1,y1,n,a,u,v)) return 1;
    if (getbit(q,j) && addW(p,x1,y1,x1,y1,x,y,n,u,v)) return 1;
    j--;
  }
 exit_multiplyW2:
  if (neg) { mpz_neg(y1,y1); mpz_neg(q,q); }
  return 0;
}

/* (x,y) can be identical to (x1,y1) */
int duplicateW(p,x1,y1,x,y,n,a,u,v) mpz_t p,x1,y1,x,y,n,a,u,v;
{
  mpz_mul_ui(u,y,2);
  mpz_gcdext(p,v,NULL,u,initial_n);
  if (mpz_cmp_ui(p,1)!=0) return(1);
  mpz_mul(u,x,x); mpz_mul_ui(u,u,3); mpz_add(u,u,a); mod(u,u,n);
  mpz_mulmod(p,u,v,n);
  mpz_mul(u,p,p); mpz_mul_ui(v,x,2); mpz_sub(u,u,v); mod(u,u,n);
  mpz_sub(v,x,u); mpz_mul(v,v,p); mpz_sub(y1,v,y); mod(y1,y1,n);
  mpz_set(x1,u);
  mul+=4; gcdexts++;
  return(0);
}

/* performs the following loop with only one gcdext, using Montgomery's trick:
   for (j=0;j<e;j++) {
       res=addW(p,x[j],y[j],x[j],y[j],x[j+1],y[j+1],n,u[0],v[0]);
       if (res) return(1); }
   return(0);

   Uses one inversion and 6*e multiplications for e>1 (3 muls for e=1)
*/
int addWn(p,x,y,n,u,v,e) mpz_t p,*x,*y,n,*u,*v; int e;
{
  int j;
  mpz_sub(u[e-1],x[e],x[e-1]); mpz_set(v[e-1],u[e-1]);
  for (j=e-2;j>=0;j--) {
    mpz_sub(u[j],x[j+1],x[j]);
    mpz_mulmod(v[j],u[j],v[j+1],n); /* v[j] = u[j]*u[j+1]*...*u[e-1] */
    mul++;
  }
  mpz_gcdext(p,v[e],NULL,v[0],initial_n); if (mpz_cmp_ui(p,1)!=0) return(1);
  gcdexts++;
  for (j=0;j<e;j++) {
    /* loop invariant: v[e] = 1/(u[j]*u[j+1]*...*u[e-1]) */
    if (j!=e-1) {
      mpz_mulmod(v[j+1],v[j+1],v[e],n);
      /* restore v[e] for next loop and make u[j] free */
      mpz_mulmod(v[e],v[e],u[j],n); mul+=2; }
    /* now v[j+1] = 1/(x[j+1]-x[j]) mod n */
    mpz_sub(p,y[j+1],y[j]); mpz_mulmod(p,v[j+1],p,n);
    mpz_mul(u[j],p,p); mpz_sub(u[j],u[j],x[j]); 
    mpz_sub(x[j],u[j],x[j+1]); mod(x[j],x[j],n);
    mpz_sub(u[j],x[j+1],x[j]); mpz_mul(u[j],u[j],p); 
    mpz_sub(y[j],u[j],y[j+1]); mod(y[j],y[j],n);
    mul+=3;
  }
  return(0);
}

int addW(p,x,y,x1,y1,x2,y2,n,u,v) mpz_t p,x,y,x1,y1,x2,y2,n,u,v;
{
  mpz_sub(u,x2,x1);
  mpz_gcdext(p,v,NULL,u,initial_n); if (mpz_cmp_ui(p,1)!=0) return(1);
  mpz_sub(p,y2,y1); mpz_mulmod(p,v,p,n);
  mpz_mul(u,p,p); mpz_sub(u,u,x1); mpz_sub(v,u,x2); mod(v,v,n);
  mpz_sub(u,x1,v); mpz_mul(u,u,p); mpz_sub(y,u,y1); mod(y,y,n);
  mpz_set(x,v);
  mul+=3; gcdexts++;
  return(0);
}

/* (x,y) can be identical to (x1,y1) */
int subW(p,x,y,x1,y1,x2,y2,n,u,v) mpz_t p,x,y,x1,y1,x2,y2,n,u,v;
{
  mpz_sub(u,x1,x2);
  mpz_gcdext(p,v,NULL,u,initial_n); if (mpz_cmp_ui(p,1)!=0) return(1);
  mpz_add(p,y1,y2); mpz_mulmod(p,v,p,n);
  mpz_mul(u,p,p); mpz_sub(u,u,x1); mpz_sub(v,u,x2); mod(v,v,n);
  mpz_sub(u,x1,v); mpz_mul(u,u,p); mpz_sub(y,u,y1); mod(y,y,n);
  mpz_set(x,v);
  mul+=3; gcdexts++;
  return(0);
}

/* assuming a[0]..a[n-1] are in F[0]..F[n-1], for n a power of two,
   put the coefficients of products of 2^d consecutive x-a[i] in F[d]
   for d=1..lg(n) */
void buildF(F,T,n) mpz_t **F,*T; unsigned int n;
{
   unsigned int d,D,st=0,st1=0; int i;

   if (verbose) st=cputime(); d=0; D=1; while (D<n) {
      F[d+1]=(mpz_t*) malloc((n+1)*sizeof(mpz_t));
      for (i=0;i<=n;i++) mpz_init(F[d+1][i]);
      if (verbose) st1=cputime();
      for (i=0;i<n;i+=2*D) polymul(F[d+1]+i,F[d]+i,F[d]+i+D,D,T);
      if (verbose && D==n/2)
	printf("Product of two polynomials of degree %d took %dms\n",D,
	       cputime()-st1);
      D=2*D; d++; }
   if (verbose)
     printf("Building f from its roots took %dms\n",cputime()-st); 
}

/* puts in G the coefficients from (x-G[0])...(x-G[n-1])
   using 2*n cells in T */
void buildG(G,T,n) mpz_t *G,*T; unsigned int n;
{
   unsigned int st=0,d,D; int i;

   if (verbose) st=cputime();
   d=0; D=1; while (D<n) {
     if (d%2==0)
       for (i=0;i<n;i+=2*D) polymul(T+i,G+i,G+i+D,D,T+n);
     else
       for (i=0;i<n;i+=2*D) polymul(G+i,T+i,T+i+D,D,T+n);
     d++; D *= 2;
   }
   if (d%2) for (i=0;i<n;i++) mpz_set(G[i],T[i]);
   for (i=0;i<n-1;i++) mpz_set(G[i], G[i+1]); mpz_set_ui(G[n-1], 1);
   if (verbose) printf("Building g from its roots took %dms\n",cputime()-st);
}

/* algorithm POLYEVAL from section 3.7 of Peter Montgomery's dissertation.
Input: 
n - a power of two
G - a table (or array) of elements of R, G[i], 0<=i<n-1
a - a table (or array) of elements of R, a[i], 0<=i<n
Output: the sequence of values of G(a[i]) are stored in G[i] for 0<=i<n
*/
void polyeval(G,F,T,n) mpz_t *G,**F,*T; unsigned int n;
{
   unsigned int d,D,st=0,j; int i;

   if (verbose) st=cputime();
   D = n/2; d=lg(n)-1; while (D>=1) {
     for (i=n-2*D;i>=0;i-=2*D) {
       /* 3D copies is expensive: perhaps we could save some? */
       for (j=0;j<2*D;j++) mpz_set(T[j], G[i+j]);
       /* divide G[i..i+2D-1] by F[d][i..i+D-1] and put rem in place */
       RecursiveDivision(T+2*D, G+i, F[d]+i, D, T+3*D);
       /* divide by F[d][i+D..i+2D-1] and store remainder in G[i+D..i+2D-1] */
       RecursiveDivision(T+2*D, T, F[d]+D+i, D, T+3*D);
       for (j=0;j<D;j++) mpz_set(G[i+D+j], T[j]);
     }
     D /= 2; d--;
   } 
   if (verbose) printf("Evaluating g on roots of f took %dms\n",cputime()-st);
}

/* multiplies b[0]+b[1]*x+...+b[K-1]*x^(K-1) by c[0]+c[1]*x+...+c[K-1]*x^(K-1)
     and puts the result in a[0]+a[1]*x+...+a[2*K-2]*x^(2K-2)
     where K is a power of two. 
     t is an auxiliary storage of at least 2K coefficients.
*/
void karatsuba(a,b,c,K,t) mpz_t *a,*b,*c,*t; unsigned int K;
{
  int i;

  karatsuba1(a,b,c,K,t);
  for (i=0;i<2*K-1;i++) mpz_mod(a[i], a[i], n);
}

void karatsuba1(a,b,c,K,t) mpz_t *a,*b,*c,*t; unsigned int K;
{
   if (K==1) { mpz_mul(a[0],b[0],c[0]); mul++; }
   else { int i,k=K/2;
      for (i=0;i<k;i++) { 
         mpz_add(t[i],b[i],b[k+i]); mpz_add(t[k+i],c[i],c[k+i]); 
       }
      karatsuba1(t+K,t,t+k,k,a); /* puts (b0+b1)*(c0+c1) in t[K..2K-2] */
      karatsuba1(a,b,c,k,t); /* puts b0*c0 in a[0..K-2] */
      mpz_set_ui(a[K-1],0);
      karatsuba1(a+K,b+k,c+k,k,t); /* puts b1*c1 in a[K..2K-2] */
      /* a[K-1] = a[2K-1] = t[2K-1] = 0 */
      for (i=0;i<K-1;i++) {
         mpz_sub(t[K+i],t[K+i],a[i]); mpz_sub(t[K+i],t[K+i],a[K+i]);
      }
      for (i=0;i<K-1;i++) mpz_add(a[k+i],a[k+i],t[K+i]);
   }
}

/*
  divides a[0]+a[1]*x+...+a[2K-1]*x^(2K-1)
  by b[0]+b[1]*x+...+b[K-1]*x^(K-1)+x^K
  puts the quotient in q[0]+q[1]*x+...+q[K-1]*x^(K-1)
  and the remainder in a[0]+a[1]*x+...+a[K-1]*x^(K-1)
  K must be a power of two.
  Needs space for 2K coefficients in t.
*/
void RecursiveDivision(mpz_t *q, mpz_t *a, mpz_t *b, unsigned int K, mpz_t *t)
{
  if (K==1) { /* a0+a1*x = a1*(b0+x) + a0-a1*b0 */
    mpz_mulmod(q[0], a[1], b[0], n); mul++;
    mpz_sub(a[0], a[0], q[0]);
    mpz_set(q[0], a[1]);
  }
  else {
    DivThreeLongHalvesByTwo(q+K/2, a+K/2, b, K/2, t);
    DivThreeLongHalvesByTwo(q, a, b, K/2, t);
  }
}

/*
  divides a[0]+a[1]*x+...+a[2K-1]*x^(3K-1)
  by b[0]+b[1]*x+...+b[2K-1]*x^(2K-1)+x^(2K)
  puts the quotient in q[0]+q[1]*x+...+q[K-1]*x^(K-1)
  and the remainder in a[0]+a[1]*x+...+a[2K-1]*x^(2K-1)
  K must be a power of two.
  Needs space for 4K coefficients in t.
*/
void DivThreeLongHalvesByTwo(mpz_t *q, mpz_t *a, mpz_t *b, 
			     unsigned int K, mpz_t *t)
{
  if (K==1) { 
    /* a0+a1*x+a2*x^2 = a2*(b0+b1*x+x^2) + (a1-a2*b1)*x + (a0-a2*b0) */
    mpz_mulmod(q[0], a[2], b[1], n);
    mpz_sub(a[1], a[1], q[0]);
    mpz_mulmod(q[0], a[2], b[0], n); mul += 2;
    mpz_sub(a[0], a[0], q[0]);
    mpz_set(q[0], a[2]);
  }
  else { int i;
    RecursiveDivision(q, a+K, b+K, K, t);
    karatsuba1(t, q, b, K, t+2*K); /* needs 2K memory in t+2K */
    for (i=0;i<=2*K-2;i++) {
      mpz_sub(a[i], a[i], t[i]); mpz_mod(a[i], a[i], n);
    }
    /* we could also have a special version of karatsuba that directly
       subtract the result to a */
  }
}

/* multiplies b[0]+...+b[k-1]*x^(k-1)+x^k by c[0]+...+c[k-1]*x^(k-1)+x^k */
void polymul(a,b,c,k,t) mpz_t *a,*b,*c,*t; unsigned int k;
{
  unsigned int i;
  karatsuba1(a,b,c,k,t);
  for (i=0;i<k;i++) mpz_mod(a[i],a[i],n);
  for (i=k;i<2*k-1;i++) {
    mpz_add(a[i],a[i],b[i-k]); mpz_add(a[i],a[i],c[i-k]); mpz_mod(a[i],a[i],n);
  }
  mpz_add(a[2*k-1],b[k-1],c[k-1]); mpz_mod(a[2*k-1],a[2*k-1],n);
}



/*******************************************************/
int ecmFactor(mpz_t p, mpz_t N, s32 _B1, double B2, int iter, mpz_t s)
{ int r, e=32;
  static int initialized=0;

  if (!(initialized)) {
    ecm_initGlobals();
    initialized=1;
  }
  mpz_set(n, N);
  B1 = _B1;
  if (B1<0) { 
    printf("Error: negative B1\n");
    return -1;
  }

   /* initialize table of primes */
   bb=0; initprimes((double)B1,0);
   
   if (B2 <= (double)B1)
//       B2 = default_B2((double) B1);
       B2 = default_B2((double) B1)*4;
   mpz_set_ui(x,0);
   


     if ((r=ecm(p,n,s,B2,e,iter))) {
       if (mpz_cmp(p,n)) {
	 if (mpz_probab_prime_p(p,25))
	   printf("Found probable prime factor");
	 else printf("Found COMPOSITE factor");
	 printf(" of %u digits: ",nb_digits(p));
	 mpz_out_str(stdout,10,p); putchar('\n');
         if (mpz_probab_prime_p(p,25) && nb_digits(p)>=48) {
	   printf("Report your potential champion to Richard Brent <rpb@comlab.ox.ac.uk>\n");
	   printf("(see ftp://ftp.comlab.ox.ac.uk/pub/Documents/techpapers/Richard.Brent/champs.txt)\n");
	 }
	 mpz_divexact(n,n,p);
	 if (mpz_probab_prime_p(n,25)) printf("Probable prime");
	 else printf("Composite");
	 printf(" cofactor ");
	 mpz_out_str(stdout,10,n);
	 printf(" has %u digits",nb_digits(n));
	   }
       else printf("Found input number N");
       printf("\n"); fflush(stdout);
     }
   return r;
}

#ifdef _MY_MAIN
/*********************************************/
int main(int argc, char *argv[])
{ int    r,iter=0; 
  mpz_t  p,s; 
  long   B1;
  double B2;
  time_t now;
  clock_t cl;
  
  mpz_init(p);
  mpz_init(s);

  
  if (argc < 3) {
    printf("Usage: %s <B1> <N>.\n", argv[0]);
    exit(-1);
  }
  B1 = atoi(argv[1]);
  B2 = 0;
  mpz_set_str(n,argv[2], 10);
  mpz_set_ui(s, 0);

  cl = clock();
  time(&now);
  prandseed(now, 17*now + 41007*cl, cl);
  
  r = ecmFactor(p, n, B1, B2, iter, s);
  printf("r = %d\n", r);
  printf("n = "); mpz_out_str(stdout, 10, n); printf("\n");
  
}


#endif
