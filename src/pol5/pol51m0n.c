/* pol51m0n.c

   Copyright (C) 2005 T. Kleinjung.
   This file is part of pol5, distributed under the terms of the
   GNU General Public Licence and WITHOUT ANY WARRANTY.
                                                                                                               
  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  2111-1307, USA.
                                                                                                               
  CJM, 2/18/05: Hacked up for inclusion in GGNFS. It was originally written by
  T. Kleinjung and/or Jens Franke.
*/

/*
#define ZEIT
*/


#ifdef HAVE_ASM_INTEL
#define HAVE_ASM
#endif

#ifdef HAVE_ASM_ALPHA
#define HAVE_ASM
#endif

/* We need to write a floorl() for Cygwin. In the meantime: */
#if defined(__CYGWIN__) || defined(__MINGW32__) || defined(MINGW32)
#undef HAVE_FLOAT64
#endif

#include <unistd.h>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <stdio.h>
#include "gmp.h"
#include "if.h"
#include <limits.h>
#include "fnmatch.h"
#include <string.h>
#define START_MESSAGE \
"----------------------------------------------------\n"\
"|    pol51m0n GNFS polynomial selection program    |\n"\
"| This program is copyright (c) 2005, by Thorsten  |\n"\
"| Kleinjung and Jens Franke, and is subject to the |\n"\
"| terms of the GNU GPL version 2.                  |\n"\
"| This program is part of gnfs4linux.              |\n"\
"----------------------------------------------------\n"



#define  MULTIPLIER            60   /* 2*2*3*5 */
#define  P0_MAX             46300   /* <=2^15.5 */

extern int asm_hash1(uint);
extern int asm_hash2(uint);

mpz_t gmp_N, gmp_a5_begin, gmp_a5_end;
int compress;
char *input_line=NULL;
size_t input_line_alloc=0;
double norm_max, skewness_min, skewness_max, a3_max;
char *basename, *filename_data, *output_name;
FILE *outputfile;
mpz_t gmp_root;
mpz_t gmp_Na5, gmp_approx;
mpz_t gmp_help1, gmp_help2, gmp_help3, gmp_help4;
mpz_t gmp_a5, gmp_a4, gmp_a3, gmp_a2, gmp_a1, gmp_a0, gmp_m0;
double alpha, dbl_a5_min, dbl_a5_max, a3_b_help, dbl_N, dbl_a5;



#define NPR5      50

uint N_mod_pr[NPR5], p_N_inv[NPR5], p_N_mod_p2[NPR5], p_a5_mod_p2[NPR5];
uint N_inv_mod_pr[NPR5], a5_mod_pr[NPR5], p_minus5a5_mod_p[NPR5];
double p_log[NPR5], p_size_max;

#ifdef HAVE_FLOAT64
long double ld_p_inv[NPR5];
#endif

int p_ind[NPR5];
uint p_pr[NPR5], p_fr[NPR5][5];
uint p_inv_table[NPR5][NPR5];
uint p_mod_p2[NPR5+1];
int npr_in_p, npr_excess, npr_total;
uint pr_step[NPR5], pr_start[NPR5];
uchar *bitarray_5power[NPR5];
uchar ucmask[8]={ 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

mpz_t gmp_prod, gmp_d, gmp_disp;
mpz_t gmp_D0, gmp_D[NPR5][5];
mpz_t gmp_kappa0, gmp_kappa[NPR5][5];
mpz_t gmp_kappa_help[NPR5+1];
double dbl_kappa_help[NPR5+1];

uint64_t ull_kappa0, ull_kappa[NPR5][5], raw_ull_kappa[NPR5][5];
uint ul_d_help[NPR5][5];
uint64_t *s1, *s2, *s1sort, *s2sort;
int s1len, s2len;
uint64_t ull_bound, lambda0, raw_ull_bound;
volatile uint raw_bound asm("raw_bound");
uint64_t *stored_pairs, *raw_stored_pairs;
int store_len, nstore, raw_store_len, nraw_store;
int *raw_cand, nraw_cand;
volatile int*raw_cand_ptr asm ("raw_cand_ptr");
uint *raw_cand_hash;

uint64_t ull_kappa_p_inv[NPR5+1];
uint ull_kappa_help0[NPR5+1];

uint prep_p[NPR5], prep_5a5[NPR5], prep_N_mod_p2[NPR5];
uint prep_other_mod_p2[NPR5];
uint prep_inv_table[NPR5][NPR5];
uint64_t ull_p_inv[NPR5];


#define HASHSHIFT  52  /*dont change this (or disable/adapt the asm-functions)*/
#define NHASH   (1<<(64-HASHSHIFT))

#define HASHSHIFT32  (HASHSHIFT-32)

uint hashdata[32*NHASH+(NHASH>>2)];
volatile uint *hashdataptr asm ("hashdataptr");

uint *hashptr32[NHASH];
uint *s1sortl, *s2sortl;
uint *s11l asm("s11l");
uint *s12l asm("s12l");
uint *s21l asm("s21l");
uint *s22l asm("s22l");;

uint64_t *shelp;
uint s11len asm ("s11len");
uint s12len asm ("s12len");
uint s21len asm ("s21len");
uint s22len asm ("s22len");
int len1, len2, len11, len12, len21, len22;


int sort[NHASH];
uint64_t *hashptr[NHASH];

uint pr_mod5[]={
11, 31, 41, 61, 71, 101, 131, 151, 181, 191, 
211, 241, 251, 271, 281, 311, 331, 401, 421, 431, 
461, 491, 521, 541, 571, 601, 631, 641, 661, 691, 
701, 751, 761, 811, 821, 881, 911, 941, 971, 991,
1021, 1031, 1051, 1061, 1091, 1151, 1171, 1181, 1201, 1231
};  /* primes =1 mod 5, < 2^(31/3) */
int npr_mod5=50;
int success=0;

uint *p0_list, *p0_root, last_p0, p0_limit;
int p0_list_len, p0_list_len_max;
int data_not_copied;
mpz_t gmp_prod_1;
uint ul_d_help_1[NPR5][5], p0_inv_p[NPR5], p_inv_p0[NPR5];

uint smallprimes[]={
2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
233, 239, 241, 251
};
uint nsmallprimes=54;


/* new sorted knapsack: */

volatile uint hashpart_shift asm ("hashpart_shift");
volatile uint hash_shift asm ("hash_shift");
uint n_hash_parts;
uint *s12l_sort asm("s12l_sort");
uint *s12l_ind, s12l_sort_maxlen, s12l_sort_len;

uint *s22l_sort asm ("s22l_sort");
uint *s22l_ind, s22l_sort_maxlen, s22l_sort_len;
uint64_t *shelpsort;
int *indhelpsort, *indhelp;
uint *s11_begin asm ("s11_begin");
uint *s21_begin asm ("s21_begin");


/* statistics */
uint64_t stat_n_a5=0, stat_n_pr=0, stat_n_p=0, stat_n_p_p0=0, stat_n_raw=0;
uint64_t stat_n_polexpand=0, stat_n_survivors=0;

/* ----------------------------------------------- */

volatile uint modulo32 asm ("modulo32");

#ifdef HAVE_ASM
#ifdef HAVE_ASM_INTEL
static inline uint modmul32(uint x,uint y)
{
  uint res,clobber;
#ifdef _MSC_VER
  __asm
  {
		mov		eax,x
		mul		y
		div		modulo32
		mov		res,edx
  }
#else
  __asm__ __volatile__("mull %%ecx\n"
           "divl modulo32" : "=d" (res), "=a" (clobber) : "a" (x), "c" (y) :
           "cc");
#endif
  return res;
}
#else
#if 0
static inline uint modmul32(uint x,uint y)
{
  uint res;
  __asm__ __volatile__("mulq %1,%2,%0\n"
           "remq %0,%3,%0" : "=r" (res) :
           "r" (x), "r" (y), "r" (modulo32) :
           "0" );
  return res;
}
#else
static inline uint modmul32(uint x,uint y)
{
  uint64_t yy, rr;

  rr=(uint64_t)x; yy=(uint64_t)y;
  rr*=yy; rr%=((uint64_t)modulo32);
  return (uint)rr;
}
#endif
#endif
#else
static inline uint modmul32(uint x,uint y)
{
  uint res;
  uint64_t xx, yy, rr;

  xx=(uint64_t)x; yy=(uint64_t)y;
  rr=xx*yy; rr%=((uint64_t)modulo32);
  res=(uint)rr;
  return res;
}
#endif


uint powmod32(uint a, uint e)
{
  uint ex=e;
  uint aa, res;

  if (!ex) return 1;
  aa=a; res=1;
  while (1) {
    if (ex&1) res=modmul32(res,aa);
    ex>>=1;
    if (!ex) return res;
    aa=modmul32(aa,aa);
  }
  return 0; /* never reached */
}


void mpz_get_ull_64(uint64_t* targptr, mpz_t src)
{
  if (src->_mp_size == 0) *targptr=0;
  else if (src->_mp_size ==1) *targptr=(uint64_t)(src->_mp_d[0]);
  else {
    *targptr=(uint64_t)(src->_mp_d[1]);
    *targptr<<=32;
    *targptr+=(uint64_t)(src->_mp_d[0]);
  }
}


#ifdef HAVE_FLOAT64
long double mpz_get_ld(mpz_t src) /* long double has 64bit precision */
{
  long double res;
  int nb;

  if (mpz_sgn(src)<0) complain("mpz_get_ld: negative\n");
  nb=mpz_size(src)-1;
  if (nb<0) return 0.L;
  res=(long double)mpz_getlimbn(src,nb); nb--;
  if (nb<0) return res;
  res*=4294967296.L;
  res+=(long double)mpz_getlimbn(src,nb); nb--;
  if (nb<0) return res;
  res*=4294967296.L;
  res+=(long double)mpz_getlimbn(src,nb);
  while (nb>0) { res*=4294967296.L; nb--; }
  return res;
}
#endif

int five_pow(int n)
{
  int res, i;

  res=1;
  for (i=0; i<n; i++) res*=5;
  return res;
}

/* ----------------------------------------------- */

void get_options(int argc, char **argv)
{
  char c;

  basename=NULL;
  compress=0; npr_in_p=7;
  norm_max=1e+20;
  p0_limit=P0_MAX;
  mpz_init(gmp_a5_begin); mpz_init(gmp_a5_end);
  while ((c=getopt(argc,argv,"b:a:A:l:n:p:vz")) != (char)(-1)) {
    switch(c) {
    case 'b':
      basename=optarg;
      break;
    case 'n':
      if(sscanf(optarg,"%lf",&norm_max)!=1)
        complain("Bad argument to -n!\n");
      break;
    case 'p':
      numread(optarg,&npr_in_p);
      break;
    case 'a':
      mpz_set_str(gmp_a5_begin,optarg,10);
      mpz_mul_ui(gmp_a5_begin,gmp_a5_begin,1000000);
      mpz_add_ui(gmp_a5_begin,gmp_a5_begin,1);
      break;
    case 'A':
      mpz_set_str(gmp_a5_end,optarg,10);
      mpz_mul_ui(gmp_a5_end,gmp_a5_end,1000000);
      break;
    case 'l':
      numread(optarg,&p0_limit);
      if (p0_limit>P0_MAX) p0_limit=P0_MAX;
      if (p0_limit<1) p0_limit=1;
      break;
    case 'v':
      verbose++;
      break;
    case 'z':
      compress=1;
      break;
    default:
      fprintf(stderr,"Bad Option %c\n",(char)c);
      Schlendrian("");
    }
  }
  if (basename==NULL) complain("argument '-b basename' is necessary\n");
  asprintf(&filename_data,"%s.data",basename);
  if (npr_in_p<4) complain("npr_in_p (option -p) must be >=4\n");
}



void read_data()
{
  char *line;
  FILE *fi;

  mpz_init(gmp_N);
  if ((fi=fopen(filename_data,"r"))==NULL)
    complain("File not found: %s\n",filename_data);
  while (getline(&input_line,&input_line_alloc,fi)>0) {
    line=input_line;
    if (*line=='N') {
      mpz_set_str(gmp_N,line+2,10);
      continue;
    }
  }
  fclose(fi);
}


int find_output_name()
{
  struct stat statbuf;
  char *tmp_name;

  asprintf(&output_name,"%s.51.m",basename);
  if (stat(output_name,&statbuf)) {
    asprintf(&output_name,"%s.51.m.gz",basename);
    if (stat(output_name,&statbuf)) return 0;
  } else {
    asprintf(&tmp_name,"%s.51.m.gz",basename);
    if (!stat(tmp_name,&statbuf))
      complain("Both files %s and %s exist.\n",output_name,tmp_name);
    free(tmp_name);
  }
  return 1;
}


void open_outputfile()
{
  int exists=0;
  char *output_cmd;

  exists=find_output_name();
  if (!exists) {
    if (compress) {
      asprintf(&output_name,"%s.51.m.gz",basename);
    } else {
      asprintf(&output_name,"%s.51.m",basename);
    }
  }
  if (fnmatch("*.gz",output_name,0)) {
    if ((outputfile=fopen(output_name,"a"))==NULL) {
      fprintf(stderr,"%s ",output_name);
      complain("cannot open file %s",output_name);
    }
  } else {
    asprintf(&output_cmd,"gzip --best --stdout >> %s",output_name);
    if ((outputfile=popen(output_cmd,"w"))==NULL) {
      fprintf(stderr,"%s ",output_name);
      complain("cannot open file %s",output_name);
    }
    free(output_cmd);
  }
}


void close_outputfile()
{
  if (fnmatch("*.gz",output_name,0)) fclose(outputfile);
  else pclose(outputfile);
}

/* ----------------------------------------------------- */
/*
Computation of the fifth root of N/a5 given the root for the
previous a5. Probably this is no longer time-critical and might be
replaced by mpz_root.
*/

void gmp_root5(mpz_t targ, mpz_t approx) /* computes (N/a5)^0.2 */
{
  mpz_set(gmp_approx,approx);
  mpz_fdiv_q(gmp_Na5,gmp_N,gmp_a5);
  while (1) {
    mpz_mul(gmp_help1,gmp_approx,gmp_approx);
    mpz_mul(gmp_help1,gmp_help1,gmp_help1);
    mpz_mul(gmp_help2,gmp_help1,gmp_approx);
    mpz_mul_2exp(gmp_help2,gmp_help2,2);
    mpz_add(gmp_help2,gmp_help2,gmp_Na5);
                      /* now: gmp_help2=4*gmp_approx^5+gmp_Na5 */
    mpz_mul_ui(gmp_help1,gmp_help1,5);
    mpz_fdiv_q(gmp_help2,gmp_help2,gmp_help1);
    mpz_sub(gmp_help1,gmp_help2,gmp_approx);
    mpz_abs(gmp_help1,gmp_help1);
    mpz_set(gmp_approx,gmp_help2);
    if (mpz_cmp_ui(gmp_help1,2)<1) break;
  }
  mpz_set(targ,gmp_approx);
}


void compute_root(shift)
{
  mpz_set_ui(gmp_help1,shift);
  mpz_mul_ui(gmp_help1,gmp_help1,MULTIPLIER);
  mpz_add(gmp_a5,gmp_a5,gmp_help1);
  gmp_root5(gmp_root,gmp_root);
}

/* ----------------------------------------------- */
/*
The following functions check whether given (a5,p,d) is useful, i. e.
whether we can find a4,...,a0 and a skewness such that the de-skewed
sup-norm is less than normmax.
This is called very seldom and the following implementation is slow.
*/

double check_a2_eval(double b5, double b4, double b3, double b2, double sk) 
{
  double sk_sq=sqrt(sk), m, m0;

  m=fabs(b5*sk*sk*sk_sq);
  m0=fabs(b4*sk*sk_sq); if (m0>m) m=m0;
  m0=fabs(b3*sk_sq); if (m0>m) m=m0;
  m0=fabs(b2/sk_sq); if (m0>m) m=m0;
  return m/norm_max;
}

#define TRANSLATE k=(int)dk; dh=(double)k; d5=c5; d4=c4+5.*dh*c5; \
                  d3=c3+4.*dh*c4+10.*dh*dh*c5; \
                  d2=c2+3.*dh*c3+6.*dh*dh*c4+10.*dh*dh*dh*c5

int check_a2(double b5, double b4, double b3, double b2, double sk)
                  /* very bad!!!! */
{
  int i, k;
  double c2=b2, c3=b3, c4=b4, c5=b5, d2, d3, d4, d5;
  double s=sk, s0, ds, dk, dh, val, val0;

  ds=s/10.; dk=ds; val=check_a2_eval(c5,c4,c3,c2,s);
  for (i=0; i<1000; i++) {
/* skewness */
    if (ds>2.) {
      s0=s+ds; if (s0>skewness_max) s0=skewness_max;
      val0=check_a2_eval(c5,c4,c3,c2,s0);
      if (val0<val) { s=s0; ds*=1.4; val=val0; }
      else {
        s0=s-ds; if (s0<skewness_min) s0=skewness_min;
        val0=check_a2_eval(c5,c4,c3,c2,s0);
        if (val0<val) { s=s0; ds*=1.4; val=val0; }
        else ds/=1.4;
      }
      if (val<1.) return 1;
    }
/* translation */
    if (dk>2.) {
      TRANSLATE;
      val0=check_a2_eval(d5,d4,d3,d2,s);
      if (val0<val) { c2=d2; c3=d3; c4=d4; c5=d5; dk*=1.4; val=val0; }
      else {
        dk=-dk; TRANSLATE; dk=-dk;
        val0=check_a2_eval(d5,d4,d3,d2,s);
        if (val0<val) { c2=d2; c3=d3; c4=d4; c5=d5; dk*=1.4; val=val0; }
        else dk/=1.4;
      }
      if (val<1.) return 1;
    }
    if ((ds<10.) && (dk<10.)) break;
  }
  if (verbose>1) if (val<100.) printf(".");
  return 0;
}


int check_pol()
{
  double sk, s;
  double dbl_a5, dbl_a4, dbl_a3, dbl_a2;
  int ok;

  stat_n_polexpand++;
/* compute coefficients */
  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_help4);
  mpz_mul(gmp_help4,gmp_help4,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_a5);
  mpz_sub(gmp_help3,gmp_N,gmp_help4);
  mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_prod);
  if (mpz_sgn(gmp_help1)) complain("pol-dev5\n");

  if (!mpz_invert(gmp_help2,gmp_d,gmp_prod)) complain("pol-dev.i\n");
  mpz_mul(gmp_help4,gmp_help2,gmp_help2);
  mpz_mul(gmp_help4,gmp_help4,gmp_help4);
  mpz_mul(gmp_help4,gmp_help4,gmp_help3);
  mpz_fdiv_r(gmp_a4,gmp_help4,gmp_prod);
  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_help4);
  mpz_mul(gmp_help4,gmp_help4,gmp_a4);
  mpz_sub(gmp_help3,gmp_help3,gmp_help4);
  mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_prod);
  if (mpz_sgn(gmp_help1)) complain("pol-dev4\n");

  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_d);
  mpz_fdiv_q(gmp_a3,gmp_help3,gmp_help4);
  mpz_fdiv_q(gmp_a3,gmp_a3,gmp_prod);
  mpz_mul(gmp_a3,gmp_a3,gmp_prod);
  mpz_mul(gmp_help4,gmp_help2,gmp_help2);
  mpz_mul(gmp_help4,gmp_help4,gmp_help2);
  mpz_mul(gmp_help4,gmp_help4,gmp_help3);
  mpz_fdiv_r(gmp_help1,gmp_help4,gmp_prod);
  mpz_add(gmp_a3,gmp_a3,gmp_help1);
  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_a3);
  mpz_sub(gmp_help3,gmp_help3,gmp_help4);
  mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_prod);
  if (mpz_sgn(gmp_help1)) complain("pol-dev3\n");

  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_fdiv_q(gmp_a2,gmp_help3,gmp_help4);
  mpz_fdiv_q(gmp_a2,gmp_a2,gmp_prod);
  mpz_mul(gmp_a2,gmp_a2,gmp_prod);
  mpz_mul(gmp_help4,gmp_help2,gmp_help2);
  mpz_mul(gmp_help4,gmp_help4,gmp_help3);
  mpz_fdiv_r(gmp_help1,gmp_help4,gmp_prod);
  mpz_add(gmp_a2,gmp_a2,gmp_help1);
  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_a2);
  mpz_sub(gmp_help3,gmp_help3,gmp_help4);
  mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_prod);
  if (mpz_sgn(gmp_help1)) complain("pol-dev2\n");

  mpz_fdiv_q(gmp_a1,gmp_help3,gmp_d);
  mpz_fdiv_q(gmp_a1,gmp_a1,gmp_prod);
  mpz_mul(gmp_a1,gmp_a1,gmp_prod);
  mpz_mul(gmp_help4,gmp_help3,gmp_help2);
  mpz_fdiv_r(gmp_help1,gmp_help4,gmp_prod);
  mpz_add(gmp_a1,gmp_a1,gmp_help1);
  mpz_mul(gmp_help4,gmp_d,gmp_a1);
  mpz_sub(gmp_help3,gmp_help3,gmp_help4);
  mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_prod);
  if (mpz_sgn(gmp_help1)) complain("pol-dev1\n");

  mpz_set(gmp_a0,gmp_help3);

  mpz_fdiv_qr(gmp_help1,gmp_a3,gmp_a3,gmp_d);
  mpz_add(gmp_help2,gmp_a3,gmp_a3);
  if (mpz_cmp(gmp_d,gmp_help2)<0) {
    mpz_sub(gmp_a3,gmp_a3,gmp_d);
    mpz_add_ui(gmp_help1,gmp_help1,1);
  }
  mpz_mul(gmp_help1,gmp_help1,gmp_prod);
  mpz_add(gmp_a4,gmp_a4,gmp_help1);

  mpz_fdiv_qr(gmp_help1,gmp_a2,gmp_a2,gmp_d);
  mpz_add(gmp_help2,gmp_a2,gmp_a2);
  if (mpz_cmp(gmp_d,gmp_help2)<0) {
    mpz_sub(gmp_a2,gmp_a2,gmp_d);
    mpz_add_ui(gmp_help1,gmp_help1,1);
  }
  mpz_mul(gmp_help1,gmp_help1,gmp_prod);
  mpz_add(gmp_a3,gmp_a3,gmp_help1);

  if (verbose>2) {
    printf("polynomial:\n");
    mpz_out_str(stdout,10,gmp_a5); printf("\n");
    mpz_out_str(stdout,10,gmp_a4); printf("\n");
    mpz_out_str(stdout,10,gmp_a3); printf("\n");
    mpz_out_str(stdout,10,gmp_a2); printf("\n");
    mpz_out_str(stdout,10,gmp_a1); printf("\n");
    mpz_out_str(stdout,10,gmp_a0); printf("\n\n");
  }

  dbl_a5=mpz_get_d(gmp_a5);
  dbl_a4=mpz_get_d(gmp_a4);
  dbl_a3=mpz_get_d(gmp_a3);
  dbl_a2=mpz_get_d(gmp_a2);

  sk=skewness_max;
  s=pow(norm_max/fabs(dbl_a4),2./3.); if (s<sk) sk=s;
  s=pow(norm_max/fabs(dbl_a3),2.); if (s<sk) sk=s;

  if (verbose>2) {
    printf("quality: %f %f %f %f\n",dbl_a5/norm_max*sk*sk*sqrt(sk),dbl_a4/norm_max*sk*sqrt(sk),dbl_a3/norm_max*sqrt(sk),dbl_a2/norm_max/sqrt(sk));
  }
  if (fabs(dbl_a2/norm_max/sqrt(sk))>1.) {
#ifdef ZEIT
zeitA(19);
#endif
    ok=check_a2(dbl_a5,dbl_a4,dbl_a3,dbl_a2,sk);
#ifdef ZEIT
zeitB(19);
#endif
    if (!ok) return 0;
  } /* else pol is ok */

  stat_n_survivors++;
  return 1;
}

/* ----------------------------------------------- */
/*
Given p (product of small primes) the following functions search for 
suitable d such that we can find a polynomial of degree 5 with d/p as
zero and 'small' a5, a4 and a3.
With reasonable parameters for most p there does not exist a suitable
d. Therefore the search is first done in an approximate manner
(knapsack_raw) and them exact (knapsack_exact).
*/

#ifdef HAVE_FLOAT64
inline uint64_t dull(long double d)
{
  long double dh;

  dh=d-(long double)((int)d);
  dh*=18446744073709551616.L;
  return (uint64_t)dh;
}
#endif

/* ---------------------------------- */

#define DEG  5

void combine_raw_sort0(int len, int ind)
{
  int i, j, i0, i1, i2, ii, jj;
  uint64_t h, v0, v1, v2;
  int bt_len, bt_ind[DEG], bt_j[DEG];
  uint64_t bt_s[DEG];
  int crs_begin[DEG], crs_end[DEG];
  uint64_t *ruk;

  ruk=raw_ull_kappa[ind];
  memcpy(shelpsort,shelp,len*sizeof(uint64_t));
  memcpy(shelpsort+len,shelp,len*sizeof(uint64_t));
  memcpy(indhelpsort,indhelp,len*sizeof(int));
  memcpy(indhelpsort+len,indhelp,len*sizeof(int));
  for (j=0; j<DEG; j++) { /* binary search for begin */
    h=-ruk[j];
    i0=0; i2=len-1;
    v0=shelp[i0]; v2=shelp[i2];
    if ((v0>=h) || (v2<h)) i2=0;
    else
      while (1) {   /* v0<h<=v2 */
        i1=(i0+i2)/2; v1=shelp[i1];
        if (v1<h) { i0=i1; v0=v1; } else { i2=i1; v2=v1; }
        if (i0+1>=i2) break;
      }
    crs_begin[j]=i2;
    crs_end[j]=crs_begin[j]+len;
  }
/* initialize binary tree */
  bt_len=DEG;
  for (j=0; j<DEG; j++) {
    bt_s[j]=shelpsort[crs_begin[j]]+ruk[j];
    bt_ind[j]=j*len+indhelpsort[crs_begin[j]];
    bt_j[j]=j;
    crs_begin[j]++;
  }
  j=0;
  while (j<DEG-1) {  /* slow but DEG is small */
    if (bt_s[j+1]>=bt_s[j]) { j++; continue; }
    h=bt_s[j+1]; bt_s[j+1]=bt_s[j]; bt_s[j]=h;
    jj=bt_j[j+1]; bt_j[j+1]=bt_j[j]; bt_j[j]=jj;
    jj=bt_ind[j+1]; bt_ind[j+1]=bt_ind[j]; bt_ind[j]=jj;
    if (j) j--;
  }

  ii=0;
  while (ii<len*DEG) {
    if (!bt_len) complain("combine_raw_sort0.len=0\n");
    shelp[ii]=bt_s[0]; indhelp[ii]=bt_ind[0];
    j=bt_j[0];
    if (crs_begin[j]<crs_end[j]) {
      bt_s[0]=shelpsort[crs_begin[j]]+ruk[j];
      bt_ind[0]=j*len+indhelpsort[crs_begin[j]];
      crs_begin[j]++;  /* bt_j[0] already contains j */
    } else {
      bt_len--;
      bt_s[0]=bt_s[bt_len];
      bt_ind[0]=bt_ind[bt_len];
      bt_j[0]=bt_j[bt_len];
    }
/* sort binary tree */
    i=0;
    while (1) {
      if (2*i+1>=bt_len) break;
      j=2*i+1; if ((2*i+2<bt_len) && (bt_s[2*i+1]>bt_s[2*i+2])) j++;
      if (bt_s[i]<=bt_s[j]) break;
      h=bt_s[j]; bt_s[j]=bt_s[i]; bt_s[i]=h;
      jj=bt_j[j]; bt_j[j]=bt_j[i]; bt_j[i]=jj;
      jj=bt_ind[j]; bt_ind[j]=bt_ind[i]; bt_ind[i]=jj;
      i=j;
    }
    ii++;
  }
/* check */
/*
  for (i=1; i<len*DEG; i++) if (shelp[i-1]>shelp[i]) complain("combine_raw_sort0\n");
*/
}


void combine_raw_sort0last(uint *targ, uint *index, int len, int ind)
{
  int i, j, i0, i1, i2, ii, jj;
  uint64_t h, v0, v1, v2;
  int bt_len, bt_ind[DEG], bt_j[DEG];
  uint bt_s[DEG];
  int crs_begin[DEG], crs_end[DEG];
  uint64_t *ruk;

  ruk=raw_ull_kappa[ind];
  memcpy(shelpsort,shelp,len*sizeof(uint64_t));
  memcpy(shelpsort+len,shelp,len*sizeof(uint64_t));
  memcpy(indhelpsort,indhelp,len*sizeof(int));
  memcpy(indhelpsort+len,indhelp,len*sizeof(int));
  for (j=0; j<DEG; j++) { /* binary search for begin */
    h=-ruk[j];
    i0=0; i2=len-1;
    v0=shelp[i0]; v2=shelp[i2];
    if ((v0>=h) || (v2<h)) i2=0;
    else
      while (1) {   /* v0<h<=v2 */
        i1=(i0+i2)/2; v1=shelp[i1];
        if (v1<h) { i0=i1; v0=v1; } else { i2=i1; v2=v1; }
        if (i0+1>=i2) break;
      }
    crs_begin[j]=i2;
    crs_end[j]=crs_begin[j]+len;
  }
/* initialize binary tree */
  bt_len=DEG;
  for (j=0; j<DEG; j++) {
    bt_s[j]=(uint)((shelpsort[crs_begin[j]]+ruk[j])>>32);
    bt_ind[j]=j*len+indhelpsort[crs_begin[j]];
    bt_j[j]=j;
    crs_begin[j]++;
  }
  j=0;
  while (j<DEG-1) {  /* slow but DEG is small */
    if (bt_s[j+1]>=bt_s[j]) { j++; continue; }
    jj=bt_s[j+1]; bt_s[j+1]=bt_s[j]; bt_s[j]=jj;
    jj=bt_j[j+1]; bt_j[j+1]=bt_j[j]; bt_j[j]=jj;
    jj=bt_ind[j+1]; bt_ind[j+1]=bt_ind[j]; bt_ind[j]=jj;
    if (j) j--;
  }

#ifdef ZEIT
zeitA(22);
#endif
  ii=0;
  while (ii<len*DEG) {
    if (!bt_len) complain("combine_raw_sort0.len=0\n");
    targ[ii]=bt_s[0]; index[ii]=bt_ind[0];
    j=bt_j[0];
    if (crs_begin[j]<crs_end[j]) {
      bt_s[0]=(uint)((shelpsort[crs_begin[j]]+ruk[j])>>32);
      bt_ind[0]=j*len+indhelpsort[crs_begin[j]];
      crs_begin[j]++;  /* bt_j[0] already contains j */
    } else {
      bt_len--;
      bt_s[0]=bt_s[bt_len];
      bt_ind[0]=bt_ind[bt_len];
      bt_j[0]=bt_j[bt_len];
    }
/* sort binary tree */
    i=0;
    while (1) {
      if (2*i+1>=bt_len) break;
      j=2*i+1; if ((2*i+2<bt_len) && (bt_s[2*i+1]>bt_s[2*i+2])) j++;
      if (bt_s[i]<=bt_s[j]) break;
      jj=bt_s[j]; bt_s[j]=bt_s[i]; bt_s[i]=jj;
      jj=bt_j[j]; bt_j[j]=bt_j[i]; bt_j[i]=jj;
      jj=bt_ind[j]; bt_ind[j]=bt_ind[i]; bt_ind[i]=jj;
      i=j;
    }
    ii++;
  }
#ifdef ZEIT
zeitB(22);
#endif  
/* check */
/*
  for (i=1; i<len*DEG; i++) if (targ[i-1]>targ[i]) complain("combine_raw_sort0last\n");
*/
}


void combine_raw_sort(uint *targ, uint *index, int i0, int i1)
{
  int i, len;

  if (i1-i0<1) complain("combine_raw_sort\n");
  shelp[0]=0; indhelp[0]=0; len=1;  /* ??? ersten direkt sortieren ??? */
  for (i=i0; i<i1-1; i++) {
    combine_raw_sort0(len,i);
    len*=DEG;
  }
  combine_raw_sort0last(targ,index,len,i1-1); len*=DEG;
}


void raw_hash_sort_extend()
{
  int i;
  uint targ, i0, i1, i2, v0, v1, v2, *s;

  targ=1<<hashpart_shift;
  s=s12l_sort;
  i0=0; i2=s12len-1;
  if (s[0]<targ) {
    if (s[i2]<targ) i0=s12len;
    else {
      v0=s[i0]; v2=s[i2]; /* v0<v<=v2 */
      while (i0+1<i2) {
        i1=(i0+i2)/2; v1=s[i1];
        if (v1<targ) { i0=i1; v0=v1; } else { i2=i1; v2=v1; }
      }
      i0=i2;
    }
  }
/* append first i0 entries and targ */
  if (s12len+i0+1>s12l_sort_maxlen) {
    s12l_sort_maxlen=s12len+i0+1;
    s12l_sort=(uint *)xrealloc(s12l_sort,s12l_sort_maxlen*sizeof(uint));
    s12l_ind=(uint *)xrealloc(s12l_ind,s12l_sort_maxlen*sizeof(uint));
  }
  s12l_sort_len=s12len+i0+1;
  for (i=0; i<i0; i++) {
    s12l_sort[s12len+i]=s12l_sort[i];
    s12l_ind[s12len+i]=s12l_ind[i];
  }
  s12l_sort[s12len+i0]=targ;
  s12l_ind[s12len+i0]=len12;

  s=s22l_sort;
  i0=0; i2=s22len-1;
  if (s[0]<targ) {
    if (s[i2]<targ) i0=s22len;
    else {
      v0=s[i0]; v2=s[i2]; /* v0<v<=v2 */
      while (i0+1<i2) {
        i1=(i0+i2)/2; v1=s[i1];
        if (v1<targ) { i0=i1; v0=v1; } else { i2=i1; v2=v1; }
      }
      i0=i2;
    }
  }
/* append first i0 entries and targ */
  if (s22len+i0+1>s22l_sort_maxlen) {
    s22l_sort_maxlen=s22len+i0+1;
    s22l_sort=(uint *)xrealloc(s22l_sort,s22l_sort_maxlen*sizeof(uint));
    s22l_ind=(uint *)xrealloc(s22l_ind,s22l_sort_maxlen*sizeof(uint));
  }
  s22l_sort_len=s22len+i0+1;
  for (i=0; i<i0; i++) {
    s22l_sort[s22len+i]=s22l_sort[i];
    s22l_ind[s22len+i]=s22l_ind[i];
  }
  s22l_sort[s22len+i0]=targ;
  s22l_ind[s22len+i0]=s22len;
}


void raw_hash_sort_init(uint *sort, uint sortlen, uint *data, uint *begin, uint beginlen)
{
  int i;
  uint v, i0, i1, i2, v0, v1, v2;

  for (i=0; i<beginlen; i++) { /* search j min. s.th. data[i]+sort[j]>=2^32 */
    v=-data[i];
    i0=0; i2=sortlen-1;
    if (sort[0]<v) {
      if (sort[i2]<v) i0=sortlen;
      else {
        v0=sort[i0]; v2=sort[i2]; /* v0<v<=v2 */
        while (i0+1<i2) {
          i1=(i0+i2)/2; v1=sort[i1];
          if (v1<v) { i0=i1; v0=v1; } else { i2=i1; v2=v1; }
        }
        i0=i2;
      }
    }
    begin[i]=i0;
  }
}


int raw_hash_sort_1(uint hb)
{
  int i, j;
  uint ind, h, hsub, hp;
  uint add, *hash, mask;
  uchar *sort;

  mask=(1<<hashpart_shift)-1; mask=~mask;
  hsub=hb<<hashpart_shift;
  memset(hashdata,0,NHASH*sizeof(uchar));
  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
  for (i=0; i<s11len; i++) {
    add=s11l[i];
    j=s11_begin[i];
    while (1) {
      h=add+s12l_sort[j];
#if 0
      hp=h>>hashpart_shift;
      if (hp!=hb) break;
#else
      hp=h&mask;
      if (hp!=hsub) break;
#endif
      ind=(h-hsub)>>hash_shift;
      if (sort[ind]>=32) return 1;
      hash[sort[ind]*NHASH+ind]=h;
      sort[ind]++;
      j++;
    }
    if (j>=s12len) j-=s12len;
    s11_begin[i]=j;
  }
  return 0;
}


void raw_hash_sort_2(uint hb)
{
  int i, j, k;
  uint ind, h, hsub, hp;
  uint add, *hash;
  uchar *sort;

  hsub=hb<<hashpart_shift;
  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
  for (i=0; i<s21len; i++) {
    add=s21l[i];
    j=s21_begin[i];
    while (1) {
      h=add+s22l_sort[j];
      hp=h>>hashpart_shift;
      if (hp!=hb) break;
      ind=(h-hsub)>>hash_shift;
      for (k=0; k<sort[ind]; k++) {
        if (hash[k*NHASH+ind]-h<raw_bound) {
          if (j>=s22len) *raw_cand_ptr++=(i+((j-s22len)<<6));
          else *raw_cand_ptr++=(i+(j<<6)); /* assumes npr_in_p<25 */
          break;
        }
      }
      j++;
    }
    if (j>=s22len) j-=s22len;
    s21_begin[i]=j;
  }
}


int raw_hash_sort_3(uint hb)
{
  int k, i, j;
  uint h, ind, *hash, hsub;
  uchar *sort;

  nraw_cand=raw_cand_ptr-raw_cand;
  hsub=hb<<hashpart_shift;
  memset(hashdata,0,NHASH*sizeof(uchar));
  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
  for (k=0; k<nraw_cand; k++) {
    i=raw_cand[k]&0x0000003f; j=raw_cand[k]>>6;
/* if (j>=s22len) complain("too long\n");*/
    h=s21l[i]+s22l_sort[j];
    raw_cand_hash[k]=h;
    ind=(h-hsub)>>hash_shift;
    if (sort[ind]>=32) return 1;
    hash[sort[ind]*NHASH+ind]=h;
    sort[ind]++;
  }
  return 0;
}


void raw_store(int i1, int i2);

void raw_hash_sort_4(uint hb)
{
  int i, j, k, l, m, i1, i2;
  uint ind, h, h1, hh;
  uint add, *hash, hsub, hp;
  uchar *sort;
  uint64_t sum1, sum2;

  hsub=hb<<hashpart_shift;
  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
  for (i=0; i<s11len; i++) {
    add=s11l[i];
    j=s11_begin[i];
    if (j+s12len<s12l_sort_len) j+=s12len;
    while (1) {
      if (!j) break; j--;
      h=add+s12l_sort[j];
      hp=h>>hashpart_shift;
      if (hp!=hb) break;
      ind=(h-hsub)>>hash_shift;
      for (k=0; k<sort[ind]; k++) {
        if ((h-hash[k*NHASH+ind])<raw_bound) {
          h1=hash[k*NHASH+ind];
          for (l=0; l<nraw_cand; l++)
            if (h1==raw_cand_hash[l]) {
              i1=i+s12l_ind[j]*s11len;
              i2=(raw_cand[l]&0x0000003f)+s22l_ind[(raw_cand[l]>>6)]*s21len;
              sum1=0; sum2=0;
              for (m=0; m<len1; m++) {
                sum1+=raw_ull_kappa[m][i1%DEG];
                i1-=(i1%DEG); i1/=DEG;
              }
              for (m=0; m<len2; m++) {
                sum2+=raw_ull_kappa[len1+m][i2%DEG];
                i2-=(i2%DEG); i2/=DEG;
              }
#if 1
/* check */
              hh=(uint)((sum1+(8589934592ULL+raw_ull_bound))>>32);
              if (hh>h) {
                if (hh-h>3) complain("hash4.0\n");
              } else {
                if (h-hh>3) complain("%u %u %u %u %" PRIu64 " hash4.1\n",h,hh,h1,raw_bound,raw_ull_bound);
              }
              hh=(uint)(sum2>>32);
              if (hh>h1) {
                if (hh-h1>3) complain("%u %u %u %u hash4.2\n",h,hh,h1,raw_bound);
              } else {
                if (h1-hh>3) complain("hash4.3\n");
             }
#endif
              if (((sum1-sum2)<raw_ull_bound) || ((sum2-sum1)<raw_ull_bound)) {
                i1=i+s12l_ind[j]*s11len;
                i2=(raw_cand[l]&0x0000003f)+s22l_ind[(raw_cand[l]>>6)]*s21len;
                raw_store(i1,i2);
              }
            }
        }
      }
    }
  }
}

/* ------------------------------- */


void combine_raw0(int len, int ind)
{
  int i, j, disp;
  uint64_t add;

  if (raw_ull_kappa[ind][0]) complain("combine_raw0\n");
  disp=len;
  for (j=1; j<5; j++) {
    add=raw_ull_kappa[ind][j];
    for (i=0; i<len; i++) shelp[i+disp]=add+shelp[i];
    disp+=len;
  }
}


void combine_raw0last(uint *targ, int len, int ind)
{
  int i, j, disp;
  uint64_t add;

  if (raw_ull_kappa[ind][0]) complain("combine_raw0last\n");
  for (i=0; i<len; i++) targ[i]=(uint)(shelp[i]>>32);
  disp=len;
  for (j=1; j<5; j++) {
    add=raw_ull_kappa[ind][j];
    for (i=0; i<len; i++) targ[i+disp]=(uint)((add+shelp[i])>>32);
    disp+=len;
  }
}


void combine_raw(uint *targ, int i0, int i1)
{
  int i, len;

  if (i1-i0<1) complain("combine_raw\n");
  if (i1-i0==1) {
    for (i=0; i<5; i++) targ[i]=(uint)(raw_ull_kappa[i0][i]>>32);
    return;
  }
  for (i=0; i<5; i++) shelp[i]=raw_ull_kappa[i0][i];
  len=5;
  for (i=i0+1; i<i1-1; i++) {
    combine_raw0(len,i);
    len*=5;
  }
  combine_raw0last(targ,len,i1-1);
}


void raw_store(int i1, int i2)
{
  while (nraw_store+1>=raw_store_len) {
    raw_store_len+=16;
    raw_stored_pairs=(uint64_t *)xrealloc(raw_stored_pairs,2*raw_store_len*sizeof(uint64_t));
  }
  raw_stored_pairs[2*nraw_store]=i1; raw_stored_pairs[2*nraw_store+1]=i2;
  nraw_store++;
}

/* also have an asm-function */
int raw_hash_1()
{
  int i, j;
  uint ind, h;
  uint add, *hash;
  uchar *sort;

  memset(hashdata,0,NHASH*sizeof(uchar));
  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
  for (i=0; i<s11len; i++) {
    add=s11l[i];
    for (j=0; j<s12len; j++) {
      h=add+s12l[j];
      ind=h>>HASHSHIFT32;
      if (sort[ind]>=32) return 1;
      hash[sort[ind]*NHASH+ind]=h;
      sort[ind]++;
    }
  }
  return 0;
}


void raw_hash_2()
{
  int i, j, k;
  uint ind, h;
  uint add, *hash;
  uchar *sort;

  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
  for (i=0; i<s21len; i++) {
    add=s21l[i];
    for (j=0; j<s22len; j++) {
      h=add+s22l[j];
      ind=h>>HASHSHIFT32;
      for (k=0; k<sort[ind]; k++) {
        if (hash[k*NHASH+ind]-h<raw_bound) {
          *raw_cand_ptr++=(i+(j<<16)); /* assumes npr_in_p<25 */
          k=sort[ind];  /* avoid duplicates */
        }
      }
    }
  }
}


int raw_hash_3()
{
  int k, i, j;
  uint h, ind, *hash;
  uchar *sort;

  nraw_cand=raw_cand_ptr-raw_cand;
  memset(hashdata,0,NHASH*sizeof(uchar));
  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
  for (k=0; k<nraw_cand; k++) {
    i=raw_cand[k]&0x0000ffff; j=raw_cand[k]>>16;
    h=s21l[i]+s22l[j];
    raw_cand_hash[k]=h;
    ind=h>>HASHSHIFT32;
    if (sort[ind]>=32) return 1;
    hash[sort[ind]*NHASH+ind]=h;
    sort[ind]++;
  }
  return 0;
}


void raw_hash_4()
{
  int i, j, k, l, m, i1, i2;
  uint ind, h, h1, hh;
  uint add, *hash;
  uchar *sort;
  uint64_t sum1, sum2;

  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
  for (i=0; i<s11len; i++) {
    add=s11l[i];
    for (j=0; j<s12len; j++) {
      h=add+s12l[j];
      ind=h>>HASHSHIFT32;
      for (k=0; k<sort[ind]; k++) {
        if ((h-hash[k*NHASH+ind])<raw_bound) {
          h1=hash[k*NHASH+ind];
          for (l=0; l<nraw_cand; l++)
            if (h1==raw_cand_hash[l]) {
              i1=i+j*s11len;
              i2=(raw_cand[l]&0x0000ffff)+(raw_cand[l]>>16)*s21len;
              sum1=0; sum2=0;
              for (m=0; m<len1; m++) {
                sum1+=raw_ull_kappa[m][i1%5];
                i1-=(i1%5); i1/=5;
              }
              for (m=0; m<len2; m++) {
                sum2+=raw_ull_kappa[len1+m][i2%5];
                i2-=(i2%5); i2/=5;
              }
#if 1
/* check */
              hh=(uint)((sum1+(8589934592ULL+raw_ull_bound))>>32);
              if (hh>h) {
                if (hh-h>3) complain("hash4.0\n");
              } else {
                if (h-hh>3) complain("%u %u %u %u %" PRIu64 " hash4.1\n",h,hh,h1,raw_bound,raw_ull_bound);
              }
              hh=(uint)(sum2>>32);
              if (hh>h1) {
                if (hh-h1>3) complain("%u %u %u %u hash4.2\n",h,hh,h1,raw_bound);
              } else {
                if (h1-hh>3) complain("%u %u %u %u %u %u hash4.3\n",h,hh,h1,raw_bound,l,nraw_cand);
             }
#endif
              if (((sum1-sum2)<raw_ull_bound) || ((sum2-sum1)<raw_ull_bound)) {
                i1=i+j*s11len;
                i2=(raw_cand[l]&0x0000ffff)+(raw_cand[l]>>16)*s21len;
                raw_store(i1,i2);
              }
            }
        }
      }
    }
  }
}


/* ---------------------------------------------------- */

#ifdef HAVE_ASM_INTEL
static void ulladdmul(uint64_t *resptr, uint ulf, uint64_t *ullfptr)
{
#ifdef _MSC_VER
		__asm
		{	mov		esi,ullfptr
			mov		edi,resptr
			mov		ecx,ulf
			mov		eax,[esi]
			mul		ecx
			add		[edi],eax
			adc		[edi+4],edx
			mov		eax,[esi+4]
			mul		ecx
			add		[edi+4],eax
		}
#else
  __asm__ __volatile__ ("movl (%%esi),%%eax\n"
          "mull %%ecx\n"
          "addl %%eax,(%%edi)\n"
          "adcl %%edx,4(%%edi)\n"
          "movl 4(%%esi),%%eax\n"
          "mull %%ecx\n"
          "addl %%eax,4(%%edi)" : : "S" (ullfptr), "D" (resptr), "c" (ulf) : 
          "%edx", "%eax", "cc");
#endif
}
#elif defined HAVE_ASM_ALPHA
static void inline ulladdmul(uint64_t *resptr, uint ulf, uint64_t *ullfptr)
{
  uint64_t res, ullf, h;

  res=*resptr;
  h=(*ullfptr)*((uint64_t)ulf);
  res+=h;
  *resptr=res;
}
#else
static void ulladdmul(uint64_t *resptr, uint ulf, uint64_t *ullfptr)
{
  uint64_t res, h;

  res=*resptr;
  h=(*ullfptr)*((uint64_t)ulf);
  res+=h;
  *resptr=res;
}
#endif

#ifdef HAVE_ASM_INTEL
/* for intel this is not necessary since we have 64bit long double */
static void ull_mulh(uint64_t *resptr, uint64_t *ullf1, uint64_t *ullf2)
{
#ifdef _MSC_VER
	__asm
	{
		mov		esi,ullf1
		mov		edi,resptr
		mov		ecx,ullf2
		mov		eax,[esi+4]
        mov		ebx,[ecx+4]
        mul		ebx
        mov		[edi],eax
        mov		[edi+4],edx
        mov		eax,[esi]
        mul		ebx
        add		[edi],edx
        adc		[edi+4],0
        mov		ebx,[ecx]
        mov		eax,[esi+4]
        mul		ebx
        add		[edi],edx
        add		[edi+4],0
	}
#else
  __asm__ __volatile__ ("movl 4(%%esi),%%eax\n"
          "movl 4(%%ecx),%%ebx\n"
          "mull %%ebx\n"
          "movl %%eax,(%%edi)\n"
          "movl %%edx,4(%%edi)\n"
          "movl (%%esi),%%eax\n"
          "mull %%ebx\n"
          "addl %%edx,(%%edi)\n"
          "adcl $0,4(%%edi)\n"
          "movl (%%ecx),%%ebx\n"
          "movl 4(%%esi),%%eax\n"
          "mull %%ebx\n"
          "addl %%edx,(%%edi)\n"
          "addl $0,4(%%edi)" : : "S" (ullf1), "D" (resptr), "c" (ullf2) :
          "%edx", "%eax", "%ebx", "cc");
#endif
}
#elif defined HAVE_ASM_ALPHA
static void ull_mulh(uint64_t *resptr, uint64_t *ullf1, uint64_t *ullf2)
{
#if 0
  uint64_t res, f1, f2;
  uint64_t f10, f11, f20, f21, h;

  f1=*ullf1; f2=*ullf2;
  f10=f1&0xffffffffULL; f11=f1>>32;
  f20=f2&0xffffffffULL; f21=f2>>32;
  res=f11*f21;
  h=((f10*f21)>>32)+((f11*f20)>>32);
  res+=h;
  *resptr=res;
#else
  uint64_t res;
  __asm__ __volatile__("umulh %1,%2,%0" :
           "=r" (res) :
           "r" (*ullf1), "r" (*ullf2) :
           "0" );
  *resptr=res;
#endif
}
#else
static void ull_mulh(uint64_t *resptr, uint64_t *ullf1, uint64_t *ullf2)
{
  uint64_t res, f1, f2;
  uint64_t f10, f11, f20, f21, h;

  f1=*ullf1; f2=*ullf2;
  f10=f1&0xffffffffULL; f11=f1>>32;
  f20=f2&0xffffffffULL; f21=f2>>32;
  res=f11*f21;
  h=((f10*f21)>>32)+((f11*f20)>>32);
  res+=h;
  *resptr=res;
}
#endif


#define PREMUL

void init_knapsack_raw()
{
  int i, j, l;
  uint p, h, hh, inv, h6;
  uint dp, d0p, p2, hp, N_mod_p2, a5_mod_p2;
  double db;
#ifdef HAVE_FLOAT64
  long double dbl_prod;
  long double dbl_5a5p;
#else
  uint64_t ull_5a5p;
#endif

  mpz_set(gmp_m0,gmp_root);
/* preparation */
  mpz_set_ui(gmp_prod,1);
  for (j=0; j<npr_in_p; j++) mpz_mul_ui(gmp_prod,gmp_prod,p_pr[p_ind[j]]);
#ifdef HAVE_FLOAT64
  dbl_prod=mpz_get_ld(gmp_prod);
  dbl_5a5p=(5.L*mpz_get_ld(gmp_a5))/dbl_prod;
#else
  mpz_mul_ui(gmp_help1,gmp_a5,5);
  mpz_mul_2exp(gmp_help1,gmp_help1,64);
  mpz_fdiv_q(gmp_help4,gmp_help1,gmp_prod);
  if (mpz_sizeinbase(gmp_help4,2)>64) complain("ull5a5p-comp\n");
  mpz_get_ull_64(&ull_5a5p,gmp_help4);
#endif

  for (i=0; i<npr_in_p; i++)
    for (j=0; j<npr_in_p; j++)
      prep_inv_table[i][j]=p_inv_table[p_ind[i]][p_ind[j]];
  for (i=0; i<npr_in_p; i++) {
    p=p_pr[p_ind[i]];
    prep_p[i]=p;
    prep_5a5[i]=p_minus5a5_mod_p[p_ind[i]];
    prep_N_mod_p2[i]=p_N_mod_p2[p_ind[i]];
  }

#ifdef ZEIT
zeitA(7);
#endif
/* computation of D_{i,j} */
  p=1; db=0.;
  for (i=0; i<npr_in_p; i++) {
    p=prep_p[i];
    if (mpz_fdiv_q_ui(gmp_help1,gmp_prod,p)) complain("search.prod\n");
    p_mod_p2[i]=mpz_fdiv_ui(gmp_help1,p*p);
    inv=invert(p_mod_p2[i]%p,p);
    ull_kappa_p_inv[i]=ull_p_inv[p_ind[i]];

    modulo32=p;
    ull_kappa_help0[i]=modmul32(inv,inv);
    prep_5a5[i]=modmul32(prep_5a5[i],inv);
    h=modmul32(inv,p_fr[p_ind[i]][0]);
    ul_d_help[i][0]=h;
    mpz_mul_ui(gmp_D[i][0],gmp_help1,h);
    db+=(double)h/(double)p;

    for (j=1; j<5; j++) {
      ul_d_help[i][j]=modmul32(inv,p_fr[p_ind[i]][j]+(p-p_fr[p_ind[i]][0]));
    }
  }

#ifdef PREMUL
  for (l=0; l<npr_in_p; l++) {
    for (i=0; i<npr_in_p; i++) {
        prep_inv_table[i][l]=prep_5a5[l]*prep_inv_table[i][l];
    }
  }
#endif

  db+=0.5*(double)npr_in_p;
  mpz_set_ui(gmp_D0,0);
  for (i=0; i<npr_in_p; i++) mpz_add(gmp_D0,gmp_D0,gmp_D[i][0]);
  mpz_fdiv_r(gmp_help1,gmp_m0,gmp_prod);
  mpz_sub(gmp_disp,gmp_m0,gmp_help1);
  mpz_mul_si(gmp_help1,gmp_prod,(int)db);
  mpz_sub(gmp_disp,gmp_disp,gmp_help1);
  mpz_add(gmp_D0,gmp_D0,gmp_disp);


#ifdef ZEIT
zeitB(7); zeitA(8);
#endif
  for (i=0; i<npr_in_p; i++)
    for (j=0; j<5; j++)
      raw_ull_kappa[i][j]=0;

/* compute kappa_0 and kappa_{i,j} */
  ull_kappa0=0;
  for (i=0; i<npr_in_p; i++) {
    p=prep_p[i]; p2=p*p;
    d0p=mpz_fdiv_ui(gmp_D0,p2);
    a5_mod_p2=p_a5_mod_p2[p_ind[i]];
    N_mod_p2=p_N_mod_p2[p_ind[i]];

    modulo32=p2;
    hp=modmul32(d0p,d0p); hp=modmul32(hp,hp); hp=modmul32(hp,d0p);
    hp=modmul32(hp,a5_mod_p2);
    hp=p2-hp+N_mod_p2;
    if (hp%p) complain("%u %u %u neu-1\n",p,d0p,hp);
    h=hp/p;
    modulo32=p;
    inv=modmul32(a5_mod_p2,p_N_inv[p_ind[i]]*ull_kappa_help0[i]);
    h=modmul32(h,inv*d0p);
    ulladdmul(&ull_kappa0,h,ull_kappa_p_inv+i);

#ifdef ZEIT
zeitA(18);
#endif
    modulo32=p2;
    hh=modmul32(d0p,d0p); hh=modmul32(hh,d0p);
    h6=modmul32(hh,hh);  /* h6=d0p^6 */

    for (j=1; j<5; j++) {
      modulo32=p2;
      dp=modmul32(p_mod_p2[i],ul_d_help[i][j]);
      hh=d0p+dp; h=modmul32(hh,hh);
      h=modmul32(hh,h); h=modmul32(h,h); /* h=(d0p+dp)^6 */
      h=modmul32(h6+(p2-h),a5_mod_p2);
      h+=modmul32(N_mod_p2,dp);

      if (h%p) complain("%u %u %u neu1\n",p,d0p,dp);
      h/=p;
      modulo32=p; h=modmul32(h,inv);

      ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+i);
#if 1
    }
    for (l=0; l<i; l++) {
      modulo32=prep_p[l];
      for (j=1; j<5; j++) {
#ifndef PREMUL
        h=modmul32(prep_5a5[l]*ul_d_help[i][j],prep_inv_table[i][l]);
#else
        h=modmul32(ul_d_help[i][j],prep_inv_table[i][l]);
#endif
        ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+l);
      }
    }
    for (l=i+1; l<npr_in_p; l++) {
      modulo32=prep_p[l];
      for (j=1; j<5; j++) {
#ifndef PREMUL
        h=modmul32(prep_5a5[l]*ul_d_help[i][j],prep_inv_table[i][l]);
#else
        h=modmul32(ul_d_help[i][j],prep_inv_table[i][l]);
#endif
        ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+l);
      }
    }
#else
      for (l=0; l<i; l++) {
        modulo32=prep_p[l];
#ifndef PREMUL
        h=modmul32(prep_5a5[l]*ul_d_help[i][j],prep_inv_table[i][l]);
#else
        h=modmul32(ul_d_help[i][j],prep_inv_table[i][l]);
#endif
        ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+l);
      }
      for (l=i+1; l<npr_in_p; l++) {
        modulo32=prep_p[l];
#ifndef PREMUL
        h=modmul32(prep_5a5[l]*ul_d_help[i][j],prep_inv_table[i][l]);
#else
        h=modmul32(ul_d_help[i][j],prep_inv_table[i][l]);
#endif
        ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+l);
      }
    }
#endif
#ifdef ZEIT
zeitB(18);
#endif
  }

#ifdef ZEIT
zeitB(8); zeitA(9);
#endif




#ifdef ZEIT
zeitA(19);
#endif
#ifdef HAVE_FLOAT64
{ 
/* compute kappa0/prod + 5*a5*(D0-root)/prod^2 - 10*a5*(D0-root)^2/prod^2/D0 */ 
uint64_t ullh;
long double ddd, dde;
  mpz_sub(gmp_help2,gmp_root,gmp_D0);
  dde=-mpz_get_ld(gmp_help2);

  ddd=dde/dbl_prod;
  ddd*=dbl_5a5p;
  ullh=(uint64_t)(18446744073709551616.L*(ddd-floorl(ddd)));
  ddd*=2.*dde; ddd/=mpz_get_ld(gmp_root);
  ullh-=(uint64_t)(18446744073709551616.L*(ddd-floorl(ddd)));
  lambda0=ullh+ull_kappa0;
}
#else
{
uint64_t ullh1, ullh2;

  mpz_sub(gmp_help2,gmp_D0,gmp_root);  /* D0-root */
  mpz_mul_ui(gmp_help1,gmp_a5,5);      /* 5*a5 */
  mpz_mul(gmp_help1,gmp_help1,gmp_help2);  /* 5*a5*(D0-root) */
  mpz_mul(gmp_help3,gmp_prod,gmp_prod);    /* prod^2 */
  mpz_mul_2exp(gmp_help1,gmp_help1,64);
  mpz_fdiv_q(gmp_help4,gmp_help1,gmp_help3);
  if (mpz_sizeinbase(gmp_help4,2)>64) complain("lambda-comp\n");
  mpz_get_ull_64(&ullh1,gmp_help4);

#if 0
  mpz_add(gmp_help2,gmp_help2,gmp_help2);
  mpz_mul(gmp_help1,gmp_help1,gmp_help2);  /* 10*a5*(D0-root)^2 * 2^64 */
  mpz_mul(gmp_help3,gmp_help3,gmp_D0);
  mpz_fdiv_q(gmp_help4,gmp_help1,gmp_help3);
#else
  mpz_add(gmp_help2,gmp_help2,gmp_help2);  /* seems to be a bit faster */
  mpz_mul(gmp_help4,gmp_help4,gmp_help2);
  mpz_fdiv_q(gmp_help4,gmp_help4,gmp_D0);
#endif
  if (mpz_sizeinbase(gmp_help4,2)>64) complain("lambda-comp\n");
  mpz_get_ull_64(&ullh2,gmp_help4);
  lambda0=-ullh1+ullh2+ull_kappa0;
}
#endif
#ifdef ZEIT
zeitB(19);
#endif
#ifdef ZEIT
zeitA(19);
#endif

{
  uint64_t ullh;

  for (i=0; i<npr_in_p; i++) {
#ifdef HAVE_FLOAT64
    ullh=dull(dbl_5a5p*ld_p_inv[p_ind[i]]);
#else
    ull_mulh(&ullh,&(ull_p_inv[p_ind[i]]),&(ull_5a5p));
#endif
    for (j=1; j<5; j++) {
      ulladdmul(&(raw_ull_kappa[i][j]),ul_d_help[i][j],&ullh);
    }
  }
  for (j=0; j<5; j++) raw_ull_kappa[0][j]+=lambda0;
}
#ifdef ZEIT
zeitB(19);
#endif

  for (i=len1; i<npr_in_p; i++)
    for (j=0; j<5; j++) raw_ull_kappa[i][j]=-raw_ull_kappa[i][j];

  db=a3_max/mpz_get_d(gmp_m0);
  if ((db>=1.) || (db<0.))
    complain("a3-bound so high that all pols will pass: %f\n",db);
  db*=18446744073709551616.;
  raw_ull_bound=(uint64_t)db;
  raw_ull_bound+=(uint64_t)((npr_in_p+2)*npr_in_p*p_pr[p_ind[npr_in_p-1]]);
  if (verbose>3) printf("raw_ull_bound: %" PRIu64 "\n",raw_ull_bound);
  raw_bound=(uint)(raw_ull_bound>>31)+1+4;
#ifdef ZEIT
zeitB(9);
#endif
  data_not_copied=1;
}


void init_knapsack_raw_p0(uint p0, uint r0)
{
  int i, j, l;
  uint p, h, hh, h6, inv;
  uint dp, d0p, p2, hp, N_mod_p2, a5_mod_p2;
  double db;
#ifdef HAVE_FLOAT64
  long double dbl_prod;
  long double dbl_5a5p;
#else
  uint64_t ull_5a5p;
#endif
  uint p_minus5a5_mod_p0;

  if (data_not_copied) {
    mpz_set(gmp_prod_1,gmp_prod);
    for (i=0; i<npr_in_p; i++)
      for (j=1; j<5; j++)
        ul_d_help_1[i][j]=ul_d_help[i][j];

    data_not_copied=0;
  }

  if (last_p0!=p0) {
    for (i=0; i<npr_in_p; i++) {
      p0_inv_p[i]=invert(p0%p_pr[p_ind[i]],p_pr[p_ind[i]]);
      p_inv_p0[i]=invert(p_pr[p_ind[i]]%p0,p0);
    }
  }
  modulo32=p0;
  prep_p[npr_in_p]=p0;
  p_minus5a5_mod_p0=modmul32(mpz_fdiv_ui(gmp_a5,p0),5*p0-5);;
  ull_kappa_p_inv[npr_in_p]=18446744073709551615ULL/((uint64_t)p0);

  for (i=0; i<npr_in_p; i++)
    for (j=0; j<npr_in_p; j++)
      prep_inv_table[i][j]=p_inv_table[p_ind[i]][p_ind[j]];
  for (i=0; i<npr_in_p; i++) {
    prep_inv_table[i][npr_in_p]=p_inv_p0[i];
    prep_inv_table[npr_in_p][i]=p0_inv_p[i];
  }

  for (i=0; i<npr_in_p; i++) {
    prep_5a5[i]=p_minus5a5_mod_p[p_ind[i]];
    prep_N_mod_p2[i]=p_N_mod_p2[p_ind[i]];
  }
#ifdef ZEIT
zeitA(7);
#endif

/* compute D_{i,j} */
  mpz_mul_ui(gmp_prod,gmp_prod_1,p0);
  p=1;

#ifdef HAVE_FLOAT64
  dbl_prod=mpz_get_ld(gmp_prod);
  dbl_5a5p=(5.L*mpz_get_ld(gmp_a5))/dbl_prod;
#else
  mpz_mul_ui(gmp_help1,gmp_a5,5);
  mpz_mul_2exp(gmp_help1,gmp_help1,64);
  mpz_fdiv_q(gmp_help4,gmp_help1,gmp_prod);
  if (mpz_sizeinbase(gmp_help4,2)>64) complain("ull5a5p-comp\n");
  mpz_get_ull_64(&ull_5a5p,gmp_help4);
#endif

  h=mpz_fdiv_ui(gmp_prod_1,p0);
  inv=invert(h,p0);
  p_mod_p2[npr_in_p]=mpz_fdiv_ui(gmp_prod_1,p0*p0);
  modulo32=p0;
  ull_kappa_help0[npr_in_p]=modmul32(inv,inv);
  prep_5a5[npr_in_p]=modmul32(p_minus5a5_mod_p0,inv);
  h=modmul32(inv,r0);
  mpz_mul_ui(gmp_D0,gmp_prod_1,h);
  db=(double)h/(double)p0;

  for (i=0; i<npr_in_p; i++) {
    p=prep_p[i];
    if (mpz_fdiv_q_ui(gmp_help1,gmp_prod,p)) complain("search.prod\n");
    p_mod_p2[i]=mpz_fdiv_ui(gmp_help1,p*p);
    inv=invert(p_mod_p2[i]%p,p);

    modulo32=p;
    ull_kappa_help0[i]=modmul32(inv,inv);
    prep_5a5[i]=modmul32(prep_5a5[i],inv);
    h=modmul32(inv,p_fr[p_ind[i]][0]);
    ul_d_help[i][0]=h;
    mpz_mul_ui(gmp_D[i][0],gmp_help1,h);
    db+=(double)h/(double)p;

    for (j=1; j<5; j++) {
      ul_d_help[i][j]=modmul32(ul_d_help_1[i][j],p0_inv_p[i]);
    }
  }

#ifdef PREMUL
  for (l=0; l<npr_in_p; l++) {
    for (i=0; i<npr_in_p; i++) {
      prep_inv_table[i][l]=prep_5a5[l]*prep_inv_table[i][l];
    }
  }
  modulo32=p0;
  for (i=0; i<npr_in_p; i++) {
    prep_inv_table[i][npr_in_p]=modmul32(prep_5a5[npr_in_p],p_inv_p0[i]);
  }
#endif

  db+=0.5*(double)npr_in_p;
  for (i=0; i<npr_in_p; i++) mpz_add(gmp_D0,gmp_D0,gmp_D[i][0]);
  mpz_fdiv_r(gmp_help1,gmp_m0,gmp_prod);
  mpz_sub(gmp_disp,gmp_m0,gmp_help1);
  mpz_mul_si(gmp_help1,gmp_prod,(int)db);
  mpz_sub(gmp_disp,gmp_disp,gmp_help1);
  mpz_add(gmp_D0,gmp_D0,gmp_disp);

#ifdef ZEIT
zeitB(7); zeitA(8);
#endif
  for (i=0; i<npr_in_p; i++)
    for (j=0; j<5; j++)
      raw_ull_kappa[i][j]=0;

/* compute kappa_0 and kappa_{i,j} */
  ull_kappa0=0;

  p=p0; p2=p*p;
  d0p=mpz_fdiv_ui(gmp_D0,p2);
  a5_mod_p2=mpz_fdiv_ui(gmp_a5,p2);
  N_mod_p2=mpz_fdiv_ui(gmp_N,p2);

  modulo32=p2;
  hp=powmod32(d0p,5);
  hp=modmul32(hp,a5_mod_p2);
  hp=p2-hp+N_mod_p2;
  if (hp%p) complain("%u %u %u neu-1r\n",p,d0p,hp);
  h=hp/p;
  modulo32=p;
  inv=invert(N_mod_p2%p0,p0);
  inv=modmul32(inv,ull_kappa_help0[npr_in_p]);
  h=modmul32(h,inv);
  h=modmul32(h,a5_mod_p2);
  h=modmul32(h,d0p);
  ulladdmul(&ull_kappa0,h,ull_kappa_p_inv+npr_in_p);

  for (i=0; i<npr_in_p; i++) {
    p=prep_p[i]; p2=p*p;
    d0p=mpz_fdiv_ui(gmp_D0,p2);
    a5_mod_p2=p_a5_mod_p2[p_ind[i]];
    N_mod_p2=p_N_mod_p2[p_ind[i]];

    modulo32=p2;
    hp=modmul32(d0p,d0p); hp=modmul32(hp,hp); hp=modmul32(hp,d0p);
    hp=modmul32(hp,a5_mod_p2);
    hp=p2-hp+N_mod_p2;
    if (hp%p) complain("%u %u %u neu-1r\n",p,d0p,hp);
    h=hp/p;
    modulo32=p;
    inv=modmul32(a5_mod_p2,p_N_inv[p_ind[i]]*ull_kappa_help0[i]);
    h=modmul32(h,inv*d0p);
    ulladdmul(&ull_kappa0,h,ull_kappa_p_inv+i);

#ifdef ZEIT
zeitA(18);
#endif
    modulo32=p2;
    hh=modmul32(d0p,d0p); hh=modmul32(hh,d0p);
    h6=modmul32(hh,hh);  /* h6=d0p^6 */

    for (j=1; j<5; j++) {
      modulo32=p2;
      dp=modmul32(p_mod_p2[i],ul_d_help[i][j]);
      hh=d0p+dp; h=modmul32(hh,hh);
      h=modmul32(hh,h); h=modmul32(h,h); /* h=(d0p+dp)^6 */
      h=modmul32(h6+(p2-h),a5_mod_p2);
      h+=modmul32(N_mod_p2,dp);

      if (h%p) complain("%u %u %u neu1\n",p,d0p,dp);
      h/=p;
      modulo32=p; h=modmul32(h,inv);
      ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+i);

#ifndef PREMUL
      modulo32=p0;
      h=modmul32(prep_5a5[npr_in_p],ul_d_help[i][j]);
      h=modmul32(h,p_inv_p0[i]);
      ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+npr_in_p);

      for (l=0; l<i; l++) {
        modulo32=prep_p[l];
        h=modmul32(prep_5a5[l]*ul_d_help[i][j],prep_inv_table[i][l]);
        ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+l);
      }
      for (l=i+1; l<npr_in_p; l++) {
        modulo32=prep_p[l];
        h=modmul32(prep_5a5[l]*ul_d_help[i][j],prep_inv_table[i][l]);
        ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+l);
      }
    }
#else
    }
    for (l=0; l<i; l++) {
      modulo32=prep_p[l];
      for (j=1; j<5; j++) {
        h=modmul32(ul_d_help[i][j],prep_inv_table[i][l]);
        ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+l);
      }
    }
    for (l=i+1; l<npr_in_p+1; l++) {
      modulo32=prep_p[l];
      for (j=1; j<5; j++) {
        h=modmul32(ul_d_help[i][j],prep_inv_table[i][l]);
        ulladdmul(&(raw_ull_kappa[i][j]),h,ull_kappa_p_inv+l);
      }
    }
#endif

#ifdef ZEIT
zeitB(18);
#endif
  }

#ifdef ZEIT
zeitB(8); zeitA(9);
#endif

#ifdef ZEIT
zeitA(19);
#endif
#ifdef HAVE_FLOAT64
{ 
/* compute kappa0/prod - 5*a5*(D0-root)/prod^2 + 10*a5*(D0-root)^2/prod^2/D0 */ 
uint64_t ullh;
long double ddd, dde;
  mpz_sub(gmp_help2,gmp_root,gmp_D0);
  dde=-mpz_get_ld(gmp_help2);

  ddd=dde/dbl_prod;
  ddd*=dbl_5a5p;
  ullh=(uint64_t)(18446744073709551616.L*(ddd-floorl(ddd)));
  ddd*=2.*dde; ddd/=mpz_get_ld(gmp_root);
  ullh-=(uint64_t)(18446744073709551616.L*(ddd-floorl(ddd)));
  lambda0=ullh+ull_kappa0;
}
#else
{
uint64_t ullh1, ullh2;

  mpz_sub(gmp_help2,gmp_D0,gmp_root);  /* D0-root */
  mpz_mul_ui(gmp_help1,gmp_a5,5);      /* 5*a5 */
  mpz_mul(gmp_help1,gmp_help1,gmp_help2);  /* 5*a5*(D0-root) */
  mpz_mul(gmp_help3,gmp_prod,gmp_prod);    /* prod^2 */
  mpz_mul_2exp(gmp_help1,gmp_help1,64);
  mpz_fdiv_q(gmp_help4,gmp_help1,gmp_help3);
  if (mpz_sizeinbase(gmp_help4,2)>64) complain("lambda-comp\n");
  mpz_get_ull_64(&ullh1,gmp_help4);

#if 0
  mpz_add(gmp_help2,gmp_help2,gmp_help2);
  mpz_mul(gmp_help1,gmp_help1,gmp_help2);  /* 10*a5*(D0-root)^2 * 2^64 */
  mpz_mul(gmp_help3,gmp_help3,gmp_D0);
  mpz_fdiv_q(gmp_help4,gmp_help1,gmp_help3);
#else
  mpz_add(gmp_help2,gmp_help2,gmp_help2);  /* seems to be a bit faster */
  mpz_mul(gmp_help4,gmp_help4,gmp_help2);
  mpz_fdiv_q(gmp_help4,gmp_help4,gmp_D0);
#endif
  if (mpz_sizeinbase(gmp_help4,2)>64) complain("lambda-comp\n");
  mpz_get_ull_64(&ullh2,gmp_help4);
  lambda0=-ullh1+ullh2+ull_kappa0;
}
#endif
#ifdef ZEIT
zeitB(19);
#endif

{
  uint64_t ullh;

  for (i=0; i<npr_in_p; i++) {
#ifdef HAVE_FLOAT64
    ullh=dull(dbl_5a5p*ld_p_inv[p_ind[i]]);
#else
    ull_mulh(&ullh,&(ull_p_inv[p_ind[i]]),&(ull_5a5p));
#endif
    for (j=1; j<5; j++) {
      ulladdmul(&(raw_ull_kappa[i][j]),ul_d_help[i][j],&ullh);
    }
  }
  for (j=0; j<5; j++) raw_ull_kappa[0][j]+=lambda0;
}

  for (i=len1; i<npr_in_p; i++)
    for (j=0; j<5; j++) raw_ull_kappa[i][j]=-raw_ull_kappa[i][j];

  db=a3_max/mpz_get_d(gmp_m0);
  if ((db>=1.) || (db<0.))
    complain("a3-bound so high that all pols will pass: %f\n",db);
  db*=18446744073709551616.;
  raw_ull_bound=(uint64_t)db;
  raw_ull_bound+=(uint64_t)((npr_in_p+2)*npr_in_p*p_pr[p_ind[npr_in_p-1]]);
  if (verbose>3) printf("raw_ull_bound: %" PRIu64 "\n",raw_ull_bound);
  raw_bound=(uint)(raw_ull_bound>>31)+1+4;
#ifdef ZEIT
zeitB(9);
#endif
}

/*
Given n primes p_i such that a_5*x^5=N has 5 five solution modulo p_i
we want to compute all 5^n solutions of this congruence modulo the product
of the p_i and find an expansion N=a_5*d^5+a_4*p*d^4+a_3*p^2*d^3 where
a3 is 'small'. Let prod be the product of the p_i and fr_{ij} a fifth root
modulo p_i of N/a_5 (0<=j<5). Let d_{ij}=h_{ij}*prod/p_i such that
d_{ij}=fr_{ij} modulo p_i and 0<=d_{ij}<prod. Then the fifth roots
modulo prod of N/a5 are: sum_i d_{i,j_i}. Let D_{ij}=d_{ij}-d_{i0}.
We choose D_0 near (N/a_5)^0.2 such that D_0=sum_i d_{i0} modulo prod.
Then we can use d=D_0+sum_i D_{i,j_i} for the expansion mentioned above.
The value of a_4 corresponding to (j_i) can be written as
kappa_0+sum_i kappa_{i,j_i} modulo prod with kappa_{i0}=0. We consider 
these kappa's as integers between 0 and prod. They can be computed as follows:
kappa_{i,j}=-5*a_5*(h_{ij}-h_{i0})/p_i modulo p_k, k!=i and
kappa_{i,j}=a_5/N*(p_i/prod)*
           ((N*D_{i,j}-a_5*((d_0+D_{i,j})^6-D_0^6))/p_i) modulo p_i
and applying the chinese remainder theorem.
Now we can estimate the a_3's:
Let lambda_0 be such that:
  ((a_5*D_0^5-N)/prod+kappa_0*D_0^4)/prod=lambda_0*D_0^4
(we can replace the two D_0^4 by (N/a_5)^0.8 without loosing too much
precision).
For the other d with corresponding kappa we have:
((a_5*d^5-N)/prod+kappa*d^4)=
 (lambda_0+(kappa-kappa_0)/prod-5*a_5*(d-D_0)/prod^2+
  10*a_5*(d-D_0)^2/prod^2/d+...)*d^4
where ... are very small terms (of order n^2*a_5/d^2).
Since (kappa-kappa_0) and (d-D_0) depend 'linearly' on (j_i) this
approximation of a_3/d can be written as lambda_0+sum_i lambda_{i,j_i}
and we have to solve this knapsack problem.


The variables used are similar to those in the explanation above.
(In search_p_raw various optimizations are done so search_p is closer
to the description above.)
*/

int knapsack_raw()
{
  int j;
  uint h;
  int res;

#ifdef ZEIT
zeitA(10);
#endif
  for (j=0; j<5; j++) raw_ull_kappa[0][j]+=(8589934592ULL+raw_ull_bound);

  combine_raw(s11l,0,len11);
  combine_raw_sort(s12l_sort,s12l_ind,len11,len1);
  combine_raw(s21l,len1,len1+len21);
  combine_raw_sort(s22l_sort,s22l_ind,len1+len21,npr_in_p);

  for (j=0; j<5; j++) raw_ull_kappa[0][j]-=(8589934592ULL+raw_ull_bound);
#ifdef ZEIT
zeitB(10);
#endif

  nraw_store=0;
  raw_hash_sort_extend();
  raw_hash_sort_init(s12l_sort,s12len,s11l,s11_begin,s11len);
  raw_hash_sort_init(s22l_sort,s22len,s21l,s21_begin,s21len);

  for (h=0; h<n_hash_parts; h++) {
    raw_cand_ptr=raw_cand;
#ifdef ZEIT
zeitA(11);
#endif
#ifdef HAVE_ASM_INTEL
    res=asm_hash1(h);
#elif defined HAVE_ASM_ALPHA
#if 1
    res=asm_hash1(h);    
#else
    res=raw_hash_sort_1(h);
#endif
#else
    res=raw_hash_sort_1(h);
#endif
#ifdef ZEIT
zeitB(11);
#endif
    if (res) {
      printf("line too long in raw_hash_sort_1, part %u\n",h);
      continue;
    }
#ifdef ZEIT
zeitA(12);
#endif
#ifdef HAVE_ASM_INTEL
    asm_hash2(h);
#elif defined HAVE_ASM_ALPHA
#if 1
    asm_hash2(h);
#else
    raw_hash_sort_2(h);
#endif
#else
    raw_hash_sort_2(h);
#endif
#ifdef ZEIT
zeitB(12);
#endif
    if (raw_cand_ptr==raw_cand) continue;

#ifdef ZEIT
zeitA(13);
#endif
    res=raw_hash_sort_3(h);  /* aus Schleife entfernen ??? */
#ifdef ZEIT
zeitB(13);
#endif
    if (res) {
      printf("line too long in raw_hash_sort_3, part %u\n",h);
      continue;
    }
#ifdef ZEIT
zeitA(14);
#endif
    raw_hash_sort_4(h);
#ifdef ZEIT
zeitB(14);
#endif
  }
  return nraw_store;
}


int knapsack_raw_p0(uint p0)
{
  return knapsack_raw();
}

/* ----------------------------------------------- */
/*
The following is the exact version of the above 'raw' functions.
This is very time-critical and not optimized.
In fact this is quite old and should be revised.
*/

#define DEG 5

uint64_t renorm(mpz_t n)
{
  uint64_t res;

  mpz_mul_2exp(gmp_help1,n,64);
  mpz_fdiv_q(gmp_help1,gmp_help1,gmp_prod);
  if (mpz_sizeinbase(gmp_help1,2)>64) complain("renorm\n");
  mpz_get_ull_64(&res,gmp_help1);
  return res;
}


uint64_t renorm2(mpz_t n)
{
  uint64_t res;

  mpz_mul_2exp(gmp_help1,n,64);
  mpz_mul(gmp_help1,gmp_help1,gmp_a5);
  mpz_mul_ui(gmp_help1,gmp_help1,5);
  mpz_fdiv_q(gmp_help1,gmp_help1,gmp_prod);
  mpz_fdiv_q(gmp_help1,gmp_help1,gmp_prod);
  mpz_get_ull_64(&res,gmp_help1);
  return res;
}


void init_knapsack_exact()
{
  int i, j, l;
  uint p, h, hh, inv;
  uint64_t lambda0;
  uint dp, d0p, p2, hp, N_mod_p2, an_mod_p2;
  double db, dq;

/* compute D_{i,j} */
  p=1; db=0.;
  for (i=0; i<npr_in_p; i++) {
    p=p_pr[p_ind[i]];
    if (mpz_fdiv_q_ui(gmp_help1,gmp_prod,p)) complain("init_knapsack_exact.prod\n");
    h=mpz_fdiv_ui(gmp_help1,p);
    inv=invert(h,p);
    mpz_mul_ui(gmp_kappa_help[i],gmp_help1,inv);
    p_mod_p2[i]=mpz_fdiv_ui(gmp_help1,p*p);

    dbl_kappa_help[i]=(double)inv/(double)p;

    h=inv*p_fr[p_ind[i]][0]; h%=p;
    ul_d_help[i][0]=h;
    mpz_mul_ui(gmp_D[i][0],gmp_help1,h);
    db+=(double)h/(double)p;

    for (j=1; j<DEG; j++) {
      h=inv*p_fr[p_ind[i]][j]; h%=p;
      h+=(p-ul_d_help[i][0]);
      if (h>=p) h-=p;
      ul_d_help[i][j]=h;
      h+=ul_d_help[i][0];
      mpz_mul_ui(gmp_D[i][j],gmp_help1,h);
    }
  }
  db+=0.5*(double)npr_in_p;

  mpz_set_ui(gmp_D0,0);
  for (i=0; i<npr_in_p; i++) mpz_add(gmp_D0,gmp_D0,gmp_D[i][0]);
  mpz_fdiv_r(gmp_help1,gmp_m0,gmp_prod);
  mpz_sub(gmp_disp,gmp_m0,gmp_help1);
  mpz_mul_si(gmp_help1,gmp_prod,(int)db);
  mpz_sub(gmp_disp,gmp_disp,gmp_help1);
  mpz_add(gmp_D0,gmp_D0,gmp_disp);

/* compute kappa_0 and kappa_{i,j} */
  mpz_set_ui(gmp_help1,0);
  for (i=0; i<npr_in_p; i++) {
    mpz_set_ui(gmp_kappa[i][0],0);
    p=p_pr[p_ind[i]]; p2=p*p;
    d0p=mpz_fdiv_ui(gmp_D0,p2);
    an_mod_p2=p_a5_mod_p2[p_ind[i]];
    N_mod_p2=p_N_mod_p2[p_ind[i]];

    modulo32=p2;
    hp=powmod32(d0p,DEG);
    hp=modmul32(hp,an_mod_p2);
    hp=p2-hp+N_mod_p2;
    if (hp%p) complain("%u %u %u neu-1\n",p,d0p,hp);
    h=hp/p;
    modulo32=p;
    inv=modmul32(an_mod_p2,p_N_inv[p_ind[i]]);
    for (j=0; j<i; j++)
      inv=modmul32(inv,p_inv_table[p_ind[j]][p_ind[i]]);
    for (j=i+1; j<npr_in_p; j++)
      inv=modmul32(inv,p_inv_table[p_ind[j]][p_ind[i]]);
    h=modmul32(h,inv);
    h=modmul32(h,d0p);
    mpz_mul_ui(gmp_help2,gmp_kappa_help[i],h);
    mpz_add(gmp_help1,gmp_help1,gmp_help2);

    for (j=1; j<DEG; j++) {
      modulo32=p2;
      dp=modmul32(p_mod_p2[i],ul_d_help[i][j]);
      hh=d0p+dp;
      h=powmod32(hh,DEG+1);
      hh=powmod32(d0p,DEG+1);
      h=modmul32(hh+(p2-h),an_mod_p2);
      h+=modmul32(N_mod_p2,dp);

      if (h%p) complain("%u %u %u neu1\n",p,d0p,dp);
      h/=p;
      modulo32=p; h=modmul32(h,inv);
      mpz_mul_ui(gmp_help4,gmp_kappa_help[i],h);
      dq=(double)h*dbl_kappa_help[i];

      for (l=0; l<npr_in_p; l++)
        if (l!=i) {
          modulo32=p_pr[p_ind[l]];
          h=modmul32(p_minus5a5_mod_p[p_ind[l]]*ul_d_help[i][j],p_inv_table[p_ind[i]][p_ind[l]]);
          mpz_mul_ui(gmp_help3,gmp_kappa_help[l],h);
          mpz_add(gmp_help4,gmp_help4,gmp_help3);
          dq+=(((double)h)*dbl_kappa_help[l]);
        }

      h=(uint)(dq);
      mpz_mul_ui(gmp_help3,gmp_prod,h);
      mpz_sub(gmp_kappa[i][j],gmp_help4,gmp_help3);
      if ((mpz_sgn(gmp_kappa[i][j])<0)
           || (mpz_cmp(gmp_kappa[i][j],gmp_prod)>=0))
        mpz_fdiv_r(gmp_kappa[i][j],gmp_kappa[i][j],gmp_prod);
    }
  }
  mpz_fdiv_r(gmp_kappa0,gmp_help1,gmp_prod);

/* check */
  for (i=0; i<npr_in_p; i++) {
    for (j=1; j<DEG; j++) {
      mpz_sub(gmp_help2,gmp_D[i][j],gmp_D[i][0]);
      mpz_add(gmp_help1,gmp_D0,gmp_help2);
      mpz_pow_ui(gmp_help3,gmp_help1,DEG);
      mpz_mul(gmp_help3,gmp_help3,gmp_a5);
      mpz_sub(gmp_help3,gmp_help3,gmp_N);
      mpz_fdiv_qr(gmp_help3,gmp_help4,gmp_help3,gmp_prod);
      if (mpz_sgn(gmp_help4)) complain("neu2\n");
      mpz_mul(gmp_help3,gmp_help3,gmp_help1);
      mpz_mul(gmp_help3,gmp_help3,gmp_a5);
      mpz_add(gmp_help1,gmp_kappa0,gmp_kappa[i][j]);
      mpz_mul(gmp_help2,gmp_help1,gmp_N);
      mpz_add(gmp_help3,gmp_help3,gmp_help2);
      mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_prod);
      if (mpz_sgn(gmp_help1)) {
        mpz_out_str(stdout,10,gmp_prod); printf("\n");
        mpz_out_str(stdout,10,gmp_kappa0); printf("\n");
        mpz_out_str(stdout,10,gmp_kappa[i][j]); printf("\n");
        mpz_out_str(stdout,10,gmp_D0); printf("\n");
        mpz_out_str(stdout,10,gmp_D[i][0]); printf("\n");
        mpz_out_str(stdout,10,gmp_D[i][j]); printf("\n");
        mpz_out_str(stdout,10,gmp_N); printf("\n");
        mpz_out_str(stdout,10,gmp_help1); printf("\n");
        mpz_out_str(stdout,10,gmp_help4); printf("\n");
        complain("check at %d %d\n",i,j);
      }
    }
  }

  mpz_pow_ui(gmp_help1,gmp_D0,DEG-1);
  mpz_mul(gmp_help3,gmp_help1,gmp_D0);
  mpz_mul(gmp_help3,gmp_help3,gmp_a5);
  mpz_sub(gmp_help3,gmp_help3,gmp_N);
  mpz_fdiv_qr(gmp_help3,gmp_help4,gmp_help3,gmp_prod);
  if (mpz_sgn(gmp_help4)) complain("neu4\n");
  mpz_mul(gmp_help2,gmp_kappa0,gmp_help1);
  mpz_add(gmp_help3,gmp_help3,gmp_help2);
  mpz_fdiv_qr(gmp_help3,gmp_help4,gmp_help3,gmp_prod);
  if (mpz_sgn(gmp_help4)) complain("neu5\n");
  mpz_fdiv_r(gmp_help3,gmp_help3,gmp_help1);  /* modulo D0^(DEG-1) */
  mpz_mul_2exp(gmp_help3,gmp_help3,64);
  mpz_fdiv_q(gmp_help3,gmp_help3,gmp_help1);
  if (mpz_sizeinbase(gmp_help3,2)>64) complain("neu6\n");
  mpz_get_ull_64(&lambda0,gmp_help3);

/* renormalize kappa's to 2^64-range */
  ull_kappa0=renorm(gmp_kappa0);

  for (i=0; i<npr_in_p; i++) {
    ull_kappa[i][0]=0;
    for (j=1; j<DEG; j++) ull_kappa[i][j]=renorm(gmp_kappa[i][j])+renorm2(gmp_D[i][j])-renorm2(gmp_D[i][0]);
  }
  for (j=0; j<DEG; j++) ull_kappa[0][j]+=lambda0;

  db=a3_max/mpz_get_d(gmp_m0);
  if ((db>=1.) || (db<0.))
    complain("a_n-2-bound so high that all pols will pass: %f\n",db);
  db*=18446744073709551616.;
  ull_bound=(uint64_t)db;
  if (verbose>3) printf("ull_bound: %" PRIu64 "\n",ull_bound);
/*  bound=ull_bound+1;*/

  for (i=len1; i<npr_in_p; i++)
    for (j=0; j<DEG; j++) ull_kappa[i][j]=-ull_kappa[i][j];
}


void init_knapsack_exact_p0(uint p0, uint r0)
{
  int i, j, l;
  uint p, h, hh, inv;
  uint64_t lambda0;
  uint dp, d0p, p2, hp, N_mod_p2, an_mod_p2;
  double db, dq;

/* compute D_{i,j} */
  if (mpz_fdiv_q_ui(gmp_help1,gmp_prod,p0))
    complain("init_knapsack_exact.p0.prod\n");
  h=mpz_fdiv_ui(gmp_help1,p0);
  inv=invert(h,p0);
  mpz_mul_ui(gmp_kappa_help[npr_in_p],gmp_help1,inv);
  p_mod_p2[npr_in_p]=mpz_fdiv_ui(gmp_help1,p0*p0);
  dbl_kappa_help[npr_in_p]=(double)inv/(double)p0;

  p=1; db=(double)h/(double)p0;
  for (i=0; i<npr_in_p; i++) {
    p=p_pr[p_ind[i]];
    if (mpz_fdiv_q_ui(gmp_help1,gmp_prod,p)) complain("init_knapsack_exact.prod\n");
    h=mpz_fdiv_ui(gmp_help1,p);
    inv=invert(h,p);
    mpz_mul_ui(gmp_kappa_help[i],gmp_help1,inv);
    p_mod_p2[i]=mpz_fdiv_ui(gmp_help1,p*p);

    dbl_kappa_help[i]=(double)inv/(double)p;

    h=inv*p_fr[p_ind[i]][0]; h%=p;
    ul_d_help[i][0]=h;
    mpz_mul_ui(gmp_D[i][0],gmp_help1,h);
    db+=(double)h/(double)p;

    for (j=1; j<DEG; j++) {
      h=inv*p_fr[p_ind[i]][j]; h%=p;
      h+=(p-ul_d_help[i][0]);
      if (h>=p) h-=p;
      ul_d_help[i][j]=h;
      h+=ul_d_help[i][0];
      mpz_mul_ui(gmp_D[i][j],gmp_help1,h);
    }
  }
  db+=0.5*(double)npr_in_p;

  if (mpz_fdiv_q_ui(gmp_help1,gmp_prod,p0)) complain("search.prod.p0\n");
  h=mpz_fdiv_ui(gmp_help1,p0);
  h=invert(h,p0); h=(h*r0)%p0;
  mpz_mul_ui(gmp_D0,gmp_help1,h);
  mpz_set(gmp_help2,gmp_D0);
  for (i=0; i<npr_in_p; i++) mpz_add(gmp_D0,gmp_D0,gmp_D[i][0]);
  mpz_fdiv_r(gmp_help1,gmp_m0,gmp_prod);
  mpz_sub(gmp_disp,gmp_m0,gmp_help1);
  mpz_mul_si(gmp_help1,gmp_prod,(int)db);
  mpz_sub(gmp_disp,gmp_disp,gmp_help1);
  mpz_add(gmp_D0,gmp_D0,gmp_disp);
  mpz_add(gmp_disp,gmp_disp,gmp_help2);

/* compute kappa_0 and kappa_{i,j} */
  p=p0; p2=p*p;
  d0p=mpz_fdiv_ui(gmp_D0,p2);
  an_mod_p2=mpz_fdiv_ui(gmp_a5,p2);
  N_mod_p2=mpz_fdiv_ui(gmp_N,p2);

  modulo32=p2;
  hp=powmod32(d0p,DEG);
  hp=modmul32(hp,an_mod_p2);
  hp=p2-hp+N_mod_p2;
  if (hp%p) complain("%u %u %u neu-1re\n",p,d0p,hp);
  h=hp/p;
  modulo32=p;
  inv=N_mod_p2;
  for (j=0; j<npr_in_p; j++) inv=modmul32(inv,p_pr[p_ind[j]]);
  inv=invert(inv,p);
  h=modmul32(h,inv);
  h=modmul32(h,an_mod_p2);
  h=modmul32(h,d0p);
  mpz_mul_ui(gmp_help1,gmp_kappa_help[npr_in_p],h);

  for (i=0; i<npr_in_p; i++) {
    mpz_set_ui(gmp_kappa[i][0],0);
    p=p_pr[p_ind[i]]; p2=p*p;
    d0p=mpz_fdiv_ui(gmp_D0,p2);
    an_mod_p2=p_a5_mod_p2[p_ind[i]];
    N_mod_p2=p_N_mod_p2[p_ind[i]];

    modulo32=p2;
    hp=powmod32(d0p,DEG);
    hp=modmul32(hp,an_mod_p2);
    hp=p2-hp+N_mod_p2;
    if (hp%p) complain("%u %u %u neu-1re\n",p,d0p,hp);
    h=hp/p;
    modulo32=p;
    inv=modmul32(an_mod_p2,p_N_inv[p_ind[i]]);
    for (j=0; j<i; j++)
      inv=modmul32(inv,p_inv_table[p_ind[j]][p_ind[i]]);
    for (j=i+1; j<npr_in_p; j++)
      inv=modmul32(inv,p_inv_table[p_ind[j]][p_ind[i]]);
    inv=modmul32(inv,invert(p0%p,p));

    h=modmul32(h,inv);
    h=modmul32(h,d0p);
    mpz_mul_ui(gmp_help2,gmp_kappa_help[i],h);
    mpz_add(gmp_help1,gmp_help1,gmp_help2);

    for (j=1; j<DEG; j++) {
      modulo32=p2;
      dp=modmul32(p_mod_p2[i],ul_d_help[i][j]);
      hh=d0p+dp;
      h=powmod32(hh,DEG+1);
      hh=powmod32(d0p,DEG+1);
      h=modmul32(hh+(p2-h),an_mod_p2);
      h+=modmul32(N_mod_p2,dp);

      if (h%p) complain("%u %u %u neu1\n",p,d0p,dp);
      h/=p;
      modulo32=p; h=modmul32(h,inv);
      mpz_mul_ui(gmp_help4,gmp_kappa_help[i],h);
      dq=(double)h*dbl_kappa_help[i];

      modulo32=p0;
      h=mpz_fdiv_ui(gmp_a5,p0);
      h=modmul32(h,5*p0-5);
      h=modmul32(h,ul_d_help[i][j]);
      h=modmul32(h,invert(p_pr[p_ind[i]]%p0,p0));
      mpz_mul_ui(gmp_help3,gmp_kappa_help[npr_in_p],h);
      mpz_add(gmp_help4,gmp_help4,gmp_help3);
      dq+=(((double)h)*dbl_kappa_help[npr_in_p]);

      for (l=0; l<npr_in_p; l++)
        if (l!=i) {
          modulo32=p_pr[p_ind[l]];
          h=modmul32(p_minus5a5_mod_p[p_ind[l]]*ul_d_help[i][j],p_inv_table[p_ind[i]][p_ind[l]]);
          mpz_mul_ui(gmp_help3,gmp_kappa_help[l],h);
          mpz_add(gmp_help4,gmp_help4,gmp_help3);
          dq+=(((double)h)*dbl_kappa_help[l]);
        }

      h=(uint)(dq);
      mpz_mul_ui(gmp_help3,gmp_prod,h);
      mpz_sub(gmp_kappa[i][j],gmp_help4,gmp_help3);
      if ((mpz_sgn(gmp_kappa[i][j])<0)
           || (mpz_cmp(gmp_kappa[i][j],gmp_prod)>=0))
        mpz_fdiv_r(gmp_kappa[i][j],gmp_kappa[i][j],gmp_prod);
    }
  }
  mpz_fdiv_r(gmp_kappa0,gmp_help1,gmp_prod);

/* check */
  for (i=0; i<npr_in_p; i++) {
    for (j=1; j<DEG; j++) {
      mpz_sub(gmp_help2,gmp_D[i][j],gmp_D[i][0]);
      mpz_add(gmp_help1,gmp_D0,gmp_help2);
      mpz_pow_ui(gmp_help3,gmp_help1,DEG);
      mpz_mul(gmp_help3,gmp_help3,gmp_a5);
      mpz_sub(gmp_help3,gmp_help3,gmp_N);
      mpz_fdiv_qr(gmp_help3,gmp_help4,gmp_help3,gmp_prod);
      if (mpz_sgn(gmp_help4)) complain("neu2\n");
      mpz_mul(gmp_help3,gmp_help3,gmp_help1);
      mpz_mul(gmp_help3,gmp_help3,gmp_a5);
      mpz_add(gmp_help1,gmp_kappa0,gmp_kappa[i][j]);
      mpz_mul(gmp_help2,gmp_help1,gmp_N);
      mpz_add(gmp_help3,gmp_help3,gmp_help2);
      mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_prod);
      if (mpz_sgn(gmp_help1)) {
        mpz_out_str(stdout,10,gmp_a5); printf("\n");
        mpz_out_str(stdout,10,gmp_prod); printf("\n");
        mpz_out_str(stdout,10,gmp_kappa0); printf("\n");
        mpz_out_str(stdout,10,gmp_kappa[i][j]); printf("\n");
        mpz_out_str(stdout,10,gmp_D0); printf("\n");
        mpz_out_str(stdout,10,gmp_D[i][0]); printf("\n");
        mpz_out_str(stdout,10,gmp_D[i][j]); printf("\n");
        mpz_out_str(stdout,10,gmp_N); printf("\n");
        mpz_out_str(stdout,10,gmp_help1); printf("\n");
        mpz_out_str(stdout,10,gmp_help4); printf("\n");
        complain("check at %d %d   p0: %u\n",i,j,p0);
      }
    }
  }

  mpz_pow_ui(gmp_help1,gmp_D0,DEG-1);
  mpz_mul(gmp_help3,gmp_help1,gmp_D0);
  mpz_mul(gmp_help3,gmp_help3,gmp_a5);
  mpz_sub(gmp_help3,gmp_help3,gmp_N);
  mpz_fdiv_qr(gmp_help3,gmp_help4,gmp_help3,gmp_prod);
  if (mpz_sgn(gmp_help4)) complain("neu4\n");
  mpz_mul(gmp_help2,gmp_kappa0,gmp_help1);
  mpz_add(gmp_help3,gmp_help3,gmp_help2);
  mpz_fdiv_qr(gmp_help3,gmp_help4,gmp_help3,gmp_prod);
  if (mpz_sgn(gmp_help4)) complain("neu5\n");
  mpz_fdiv_r(gmp_help3,gmp_help3,gmp_help1);  /* modulo D0^(DEG-1) */
  mpz_mul_2exp(gmp_help3,gmp_help3,64);
  mpz_fdiv_q(gmp_help3,gmp_help3,gmp_help1);
  if (mpz_sizeinbase(gmp_help3,2)>64) complain("neu6\n");
  mpz_get_ull_64(&lambda0,gmp_help3);

/* renormalize kappa's to 2^64-range */
  ull_kappa0=renorm(gmp_kappa0);

  for (i=0; i<npr_in_p; i++) {
    ull_kappa[i][0]=0;
    for (j=1; j<DEG; j++) ull_kappa[i][j]=renorm(gmp_kappa[i][j])+renorm2(gmp_D[i][j])-renorm2(gmp_D[i][0]);
  }
  for (j=0; j<DEG; j++) ull_kappa[0][j]+=lambda0;

  db=a3_max/mpz_get_d(gmp_m0);
  if ((db>=1.) || (db<0.))
    complain("a_n-2-bound so high that all pols will pass: %f\n",db);
  db*=18446744073709551616.;
  ull_bound=(uint64_t)db;
  if (verbose>3) printf("ull_bound: %" PRIu64 "\n",ull_bound);
/*  bound=ull_bound+1;*/

  for (i=len1; i<npr_in_p; i++)
    for (j=0; j<DEG; j++) ull_kappa[i][j]=-ull_kappa[i][j];
}


void check_raw()
{
  int i, j, i1, i2;
  uint64_t sum1, sum2;

  i=0; /* remove duplicates: */
  while (i<nraw_store-1) {
    for (j=i+1; j<nraw_store; j++)
      if ((raw_stored_pairs[2*i]==raw_stored_pairs[2*j]) &&
          (raw_stored_pairs[2*i+1]==raw_stored_pairs[2*j+1])) break;
    if (j<nraw_store) {
      raw_stored_pairs[2*j]=raw_stored_pairs[2*nraw_store-2];
      raw_stored_pairs[2*j+1]=raw_stored_pairs[2*nraw_store-1];
      nraw_store--;
    }
    i++;
  }

  for (i=0; i<nraw_store; i++) {
    i1=raw_stored_pairs[2*i];
    if (i1>=s1len) complain("check_raw_stored_pairs.1 %d %d \n",i1,i);
    i2=raw_stored_pairs[2*i+1];
    if (i2>=s2len) complain("check_raw_stored_pairs.2 %d %d \n",i2,i);
/*printf("pair: %d %d\n",i1,i2);*/

    sum1=0; sum2=0;
    for (j=0; j<len1; j++) {
      sum1+=ull_kappa[j][i1%DEG];
      i1-=(i1%DEG); i1/=DEG;
    }
    for (j=len1; j<len1+len2; j++) {
      sum2+=ull_kappa[j][i2%DEG];
      i2-=(i2%DEG); i2/=DEG;
    }
    if ((sum1-sum2>ull_bound) && (sum2-sum1>ull_bound)) {
#if 1
      if (sum1>sum2)
        printf("%" PRIu64 " %" PRIu64 "! ",sum1-sum2,ull_bound);
      else
        printf("%" PRIu64 " %" PRIu64 "! ",sum2-sum1,ull_bound);
#else
      printf("!");
#endif
/*complain("debug\n");*/
      continue;
    }

    i1=raw_stored_pairs[2*i]; i2=raw_stored_pairs[2*i+1];
    mpz_set_ui(gmp_d,0); mpz_set_ui(gmp_help4,0);
    for (j=0; j<len1; j++) {
      mpz_add(gmp_d,gmp_d,gmp_D[j][i1%DEG]);
      mpz_add(gmp_help4,gmp_help4,gmp_kappa[j][i1%DEG]);
      i1-=(i1%DEG); i1/=DEG;
    }
    for (j=len1; j<len1+len2; j++) {
      mpz_add(gmp_d,gmp_d,gmp_D[j][i2%DEG]);
      mpz_add(gmp_help4,gmp_help4,gmp_kappa[j][i2%DEG]);
      i2-=(i2%DEG); i2/=DEG;
    }
    mpz_add(gmp_d,gmp_d,gmp_disp);
    mpz_add(gmp_help4,gmp_help4,gmp_kappa0);
/* CHECK */
    mpz_pow_ui(gmp_help3,gmp_d,DEG);
    mpz_mul(gmp_help3,gmp_help3,gmp_a5);
    mpz_sub(gmp_help3,gmp_help3,gmp_N);
    mpz_fdiv_qr(gmp_help3,gmp_help2,gmp_help3,gmp_prod);
    if (mpz_sgn(gmp_help2)) { printf("\n");
      mpz_out_str(stdout,10,gmp_d); printf("\n");
      mpz_out_str(stdout,10,gmp_prod); printf("\n");
      mpz_out_str(stdout,10,gmp_a5); printf("\n");
      mpz_out_str(stdout,10,gmp_help3); printf("\n");
      complain("check2.raw\n");
    }
    mpz_mul(gmp_help3,gmp_help3,gmp_d);
    mpz_mul(gmp_help3,gmp_help3,gmp_a5);
    mpz_mul(gmp_help1,gmp_N,gmp_help4);
    mpz_add(gmp_help3,gmp_help3,gmp_help1);
    mpz_fdiv_r(gmp_help3,gmp_help3,gmp_prod);

    if (mpz_sgn(gmp_help3)) {
      mpz_out_str(stdout,10,gmp_d); printf("\n");
      mpz_out_str(stdout,10,gmp_prod); printf("\n");
      mpz_out_str(stdout,10,gmp_help4); printf("\n");
      complain("check_raw_stored_pairs check\n");
    }

    if (check_pol()) {
      open_outputfile();
      mpz_out_str(outputfile,10,gmp_a5); fprintf(outputfile," ");
      mpz_out_str(outputfile,10,gmp_prod); fprintf(outputfile," ");
      mpz_out_str(outputfile,10,gmp_d); fprintf(outputfile,"\n");
      close_outputfile();
      success=1;
      if (verbose>1) {
        mpz_out_str(stdout,10,gmp_a5); printf(" ");
        mpz_out_str(stdout,10,gmp_d); printf(" ");
        mpz_out_str(stdout,10,gmp_prod); printf("\n");
        mpz_out_str(stdout,10,gmp_kappa0); printf(" ");
        mpz_out_str(stdout,10,gmp_D0); printf("\n");
      }
    }
  }
/*  if (verbose>2) printf("%d pairs checked\n",nraw_store);*/
}


void check_raw_p0(uint p0)
{
  return check_raw();
}


void check_stored_pairs()
{
  int i, j, i1, i2;

  i=0;
  while (i<nstore-1) {
    for (j=i+1; j<nstore; j++)
      if ((stored_pairs[2*i]==stored_pairs[2*j]) &&
          (stored_pairs[2*i+1]==stored_pairs[2*j+1])) break;
    if (j<nstore) {
      stored_pairs[2*j]=stored_pairs[2*nstore-2];
      stored_pairs[2*j+1]=stored_pairs[2*nstore-1];
      nstore--;
    }
    i++;
  }  /* slow algorithm since nstore usually is 1 */
  for (i=0; i<nstore; i++) {
/* TODO: there might be more than one i1 satisfying s1[i1]==stored_pairs[2*i] */
    for (i1=0; i1<s1len; i1++)
      if (stored_pairs[2*i]==s1[i1]) break;
    if (i1>=s1len) complain("check_stored_pairs.1\n");
    for (i2=0; i2<s2len; i2++)
      if (stored_pairs[2*i+1]==s2[i2]) break;
    if (i2>=s2len) complain("check_stored_pairs.2\n");

    mpz_set_ui(gmp_d,0); mpz_set_ui(gmp_help4,0);
    for (j=0; j<len1; j++) {
      mpz_add(gmp_d,gmp_d,gmp_D[j][i1%5]);
      mpz_add(gmp_help4,gmp_help4,gmp_kappa[j][i1%5]);
      i1-=(i1%5); i1/=5;
    }
    for (j=len1; j<len1+len2; j++) {
      mpz_add(gmp_d,gmp_d,gmp_D[j][i2%5]);
      mpz_add(gmp_help4,gmp_help4,gmp_kappa[j][i2%5]);
      i2-=(i2%5); i2/=5;
    }
    mpz_add(gmp_d,gmp_d,gmp_disp);
    mpz_add(gmp_help4,gmp_help4,gmp_kappa0);
/* CHECK */
    mpz_mul(gmp_help3,gmp_d,gmp_d);
    mpz_mul(gmp_help3,gmp_help3,gmp_help3);
    mpz_mul(gmp_help3,gmp_help3,gmp_d);
    mpz_mul(gmp_help3,gmp_help3,gmp_a5);
    mpz_sub(gmp_help3,gmp_help3,gmp_N);
    mpz_fdiv_qr(gmp_help3,gmp_help2,gmp_help3,gmp_prod);
    if (mpz_sgn(gmp_help2)) complain("check2\n");
    mpz_mul(gmp_help3,gmp_help3,gmp_d);
    mpz_mul(gmp_help3,gmp_help3,gmp_a5);
    mpz_mul(gmp_help1,gmp_N,gmp_help4);
    mpz_add(gmp_help3,gmp_help3,gmp_help1);
    mpz_fdiv_r(gmp_help3,gmp_help3,gmp_prod);

    if (mpz_sgn(gmp_help3)) {
      mpz_out_str(stdout,10,gmp_d); printf("\n");
      mpz_out_str(stdout,10,gmp_prod); printf("\n");
      mpz_out_str(stdout,10,gmp_help4); printf("\n");
      complain("check_stored_pairs check\n");
    }

    if (check_pol()) {
      open_outputfile();
      mpz_out_str(outputfile,10,gmp_a5); fprintf(outputfile," ");
      mpz_out_str(outputfile,10,gmp_prod); fprintf(outputfile," ");
      mpz_out_str(outputfile,10,gmp_d); fprintf(outputfile,"\n");
      close_outputfile();
      success=1;
      if (verbose>1) {
        mpz_out_str(stdout,10,gmp_a5); printf(" ");
        mpz_out_str(stdout,10,gmp_d); printf(" ");
        mpz_out_str(stdout,10,gmp_prod); printf("\n");
        mpz_out_str(stdout,10,gmp_kappa0); printf(" ");
        mpz_out_str(stdout,10,gmp_D0); printf("\n");
      }
    }
  }
  if (verbose>2) printf("%d pairs checked\n",nstore);
}


void store(int i1, int i2)
{
  while (nstore+1>=store_len) {
    store_len+=16;
    stored_pairs=(uint64_t *)xrealloc(stored_pairs,2*store_len*sizeof(uint64_t));
  }
  stored_pairs[2*nstore]=s1sort[i1]; stored_pairs[2*nstore+1]=s2sort[i2];
  nstore++;
}


void combine0(uint64_t *targ, int len, int ind)
{
  int i, j, disp;
  uint64_t add;

  if (ull_kappa[ind][0]) complain("combine0\n");
  disp=len;
  for (j=1; j<5; j++) {
    add=ull_kappa[ind][j];
    for (i=0; i<len; i++) targ[i+disp]=add+targ[i];
    disp+=len;
  }
}


void combinelast0(uint64_t *targ, int len, int ind)
{
  int i, j, disp;
  uint64_t add, h;

  if (ull_kappa[ind][0]) complain("combinelast0\n");
  memset(sort,0,NHASH*sizeof(int));
  disp=0;
  for (j=0; j<5; j++) {
    add=ull_kappa[ind][j];
    for (i=0; i<len; i++) {
      h=add+targ[i];
      targ[i+disp]=h;
      h>>=HASHSHIFT;
      sort[h]++;
    }
    disp+=len;
  }
}


void combine(uint64_t *targ, uint64_t *targsort, int i0, int i1)
{
  int i, len;
  uint64_t h, hh, *ptr;

  if (i1-i0<1) complain("combine\n");
  for (i=0; i<5; i++) targ[i]=ull_kappa[i0][i];
  len=5;
  for (i=i0+1; i<i1-1; i++) {
    combine0(targ,len,i);
    len*=5;
  }
  combinelast0(targ,len,i1-1); len*=5;
  hashptr[0]=targsort;
  for (i=0; i<NHASH-1; i++) hashptr[i+1]=hashptr[i]+sort[i];
  for (i=0; i<len; i++) {
    h=targ[i]; hh=h>>HASHSHIFT;
    *(hashptr[hh])=h; hashptr[hh]++;
  }
  for (i=0; i<NHASH; i++) {
    if (sort[i]>3) {
      ptr=hashptr[i]-sort[i];
      while (ptr<hashptr[i]-1) {
        if (ptr[0]>ptr[1]) {
          h=ptr[0]; ptr[0]=ptr[1]; ptr[1]=h;
          if (ptr+sort[i]>hashptr[i]) ptr--;
        } else ptr++;
      }
      continue;
    }
    if (sort[i]>1) {
      ptr=hashptr[i]-sort[i];
      if (sort[i]==2) {
        if (ptr[0]>ptr[1]) { h=ptr[0]; ptr[0]=ptr[1]; ptr[1]=h; }
      } else {
        if (ptr[0]>ptr[1]) { h=ptr[0]; ptr[0]=ptr[1]; ptr[1]=h; }
        if (ptr[1]>ptr[2]) {
          if (ptr[0]>ptr[2]) {
            h=ptr[2]; ptr[2]=ptr[1]; ptr[1]=ptr[0]; ptr[0]=h;
          } else {
            h=ptr[1]; ptr[1]=ptr[2]; ptr[2]=h;
          }
        }
      }
    }
  }
/*  for (i=0; i<len-1; i++) if (targsort[i]>targsort[i+1]) complain("combine\n");*/
}


void knapsack_exact()
{
  int i, i1, i2;
  uint64_t bound, v1, v2;

  combine(s1,s1sort,0,len1);
  combine(s2,s2sort,len1,len1+len2);

  bound=ull_bound+1;
  nstore=0;
  if ((s1sort[0]-s2sort[s2len-1])<bound) store(0,s2len-1);
  if ((-s1sort[s1len-1]+s2sort[0])<bound) store(s1len-1,0);
  i1=0; i2=0;
  while (1) {
    v1=s1sort[i1]; v2=s2sort[i2];
    if (v1<v2) {
      if ((v2-v1)<bound) {
        store(i1,i2);
        for (i=i2+1; i<s2len; i++)
          if ((s2sort[i]-v1)<bound) store(i1,i); else break;
      }
      if (i1<s1len-1) i1++;
      else break;
    } else {
      if ((v1-v2)<bound) {
        store(i1,i2);
        for (i=i1+1; i<s1len; i++)
          if ((s1sort[i]-v2)<bound) store(i,i2); else break;
      }
      if (i2<s2len-1) i2++;
      else break;
    }
  }
#ifdef ZEIT
zeitA(19);
#endif
  if (nstore) check_stored_pairs();
#ifdef ZEIT
zeitB(19);
#endif
}

/* ----------------------------------------------- */
/*
Preparations of auxilary factors.
*/

int coprime(uint a, uint b)
{
  uint r, bb=b, aa=a;

  while (bb) { r=aa%bb; aa=bb; bb=r; }
  return (aa==1);
}


void aux_factor_store(uint n, uint r)
{
  if (p0_list_len>=p0_list_len_max) {
    p0_list_len_max+=16;
    p0_list=(uint *)xrealloc(p0_list,p0_list_len_max*sizeof(uint));
    p0_root=(uint *)xrealloc(p0_root,p0_list_len_max*sizeof(uint));
    if (p0_list_len>=p0_list_len_max) complain("aux_factor_store\n");
  }
  p0_list[p0_list_len]=n;
  p0_root[p0_list_len]=r;
  p0_list_len++;
}


uint aux_factor_phi(uint n)  /* n<=2^16 */
{
  int i;
  uint phi, d, p;

  if (n>65536) complain("aux_factor_phi\n");
  d=n; phi=1;
  for (i=0; i<nsmallprimes; i++) {
    p=smallprimes[i];
    if (d%p==0) {
      phi*=(p-1); d/=p;
      while (d%p==0) { phi*=p; d/=p; }
      if (d==1) break;
    }
  }
  if (d>1) phi*=(d-1);
  return phi;
}


void aux_factor_roots(uint n)
{
  int i;
  uint h, r, a5, N, p;

  a5=mpz_fdiv_ui(gmp_a5,n);
  if (!coprime(a5,n)) return;
  N=mpz_fdiv_ui(gmp_N,n);
  if (n<20) {
    if (n==11) return;
    modulo32=n;
    for (r=1; r<n; r++) {
      h=powmod32(r,5); h=modmul32(h,a5);
      if (h==N) aux_factor_store(n,r);
    }
    return;
  }

  for (i=0; i<NPR5; i++) {
    p=pr_mod5[i];
    if (p>n) break;
    if (n%p==0) return;
  }

/* for deg=5 exactly one root exists */
  h=aux_factor_phi(n);
  if (h%5==0) return;
/*  if (h%5==0) complain("aux_factor_roots.phi %lu %lu\n",n,h);*/
  r=1+h; while (r%5) r+=h;
  r/=5;

  modulo32=n;
  h=invert(a5,n); h=modmul32(h,N);
  r=powmod32(h,r);
  aux_factor_store(n,r);
}


void init_aux_factors(double dbl_p0_max)
{
  uint n, p0m;

  p0_list[0]=1; p0_root[0]=0; p0_list_len=1;
  if (dbl_p0_max>(double)p0_limit) p0m=p0_limit;
  else p0m=(uint)rint(dbl_p0_max);
  if (p0m<2) return;
  for (n=2; n<=p0m; n++) {
    if (!coprime(n,MULTIPLIER)) continue;
    aux_factor_roots(n);
  }
}

/* ----------------------------------------------- */
/*
All the initialising functions for finding suitable a5's, i. e.
such that N=a5*x^5 mod p has many solutions modulo small p=1 (5)
*/

uint powmod(uint a, uint e, uint p)  /* assumes a<2^16 */
{
  uint ex=e;
  uint aa, res;

  if (!ex) return 1;
  aa=a; res=1;
  while (1) {
    if (ex&1) { res*=aa; res%=p; }
    ex>>=1;
    if (!ex) return res;
    aa*=aa; aa%=p;
  }
  return 0; /* never reached */
}


int is_5power(uint a, uint p) /* p!=5, 0 is not considered as fifth power */
{
  if (p%5!=1) return 0;
  if (powmod(a,(p-1)/5,p)==1) return 1;
  return 0;
}


int find_primes()
{
  int i, j, ind;
  uint h, p;
  double size;

  ind=0;
  for (i=0; i<npr_mod5; i++) {
    p=pr_mod5[i];
    pr_start[i]+=pr_step[i]; if (pr_start[i]>=p) pr_start[i]-=p;
    h=pr_start[i];
    if (bitarray_5power[i][h>>3]&ucmask[h&7]) {
      p_pr[ind]=p;
      p_log[ind]=log((double)p);
      ind++;
      if (ind==npr_in_p) break;
    }
  }
  i++;
  if (ind<npr_in_p) return 0;
  size=0.;
  for (j=0; j<npr_in_p; j++) size+=p_log[j];
  if (size>p_size_max) { /* just updating */
    for (; i<npr_mod5; i++) {
      p=pr_mod5[i];
      pr_start[i]+=pr_step[i]; if (pr_start[i]>=p) pr_start[i]-=p;
    }
    return 0;
  } /* else search for further candidates */
  for (; i<npr_mod5; i++) {
    p=pr_mod5[i];
    pr_start[i]+=pr_step[i]; if (pr_start[i]>=p) pr_start[i]-=p;
    h=pr_start[i];
    if (bitarray_5power[i][h>>3]&ucmask[h&7]) {
      p_pr[ind]=p;
      p_log[ind]=log((double)p);
      ind++;
    }
  }
  npr_total=ind; npr_excess=npr_total-npr_in_p;
  return 1;
}


void fifth_root(uint *ro, uint a, uint p) /* perhaps improve this */
{
  int ind;
  uint h, r;

  ind=0;
  for (r=1; r<p; r++) {
    h=r*r; h%=p; h=h*h; h%=p; h*=r; h%=p;
    if (h==a) {
      ro[ind++]=r;
      if (ind==5) break;
    }
  }
  if (ind<5) complain("cannot find enough fifth roots %d %d\n",a,p);
}


void init_all_pr()
{
  int i, j;
  uint h, p;

  for (i=0; i<npr_total; i++) {
    p=p_pr[i]; modulo32=p;
    p_N_mod_p2[i]=mpz_fdiv_ui(gmp_N,p*p);
    p_a5_mod_p2[i]=mpz_fdiv_ui(gmp_a5,p*p);
    p_minus5a5_mod_p[i]=modmul32(p-5,p_a5_mod_p2[i]);
    p_N_inv[i]=invert(p_N_mod_p2[i]%p,p);
    h=p_a5_mod_p2[i]%p;
    h=invert(h,p);
    h*=p_N_mod_p2[i]; h%=p;
#ifdef HAVE_FLOAT64
    ld_p_inv[i]=1.L/((long double)p);
#endif
    ull_p_inv[i]=18446744073709551615ULL/((uint64_t)p);
    fifth_root(p_fr[i],h,p);
  }
  for (i=0; i<npr_total; i++)
    for (j=0; j<npr_total; j++)
      if (i!=j)
        p_inv_table[i][j]=invert(p_pr[i]%p_pr[j],p_pr[j]);
}

/* ----------------------------------------------- */

void init_search()
{
  int i, j, len;

  hashdataptr=&(hashdata[0]); /* for asm-functions */
  mpz_init(gmp_help1);
  mpz_init(gmp_help2);
  mpz_init(gmp_help3);
  mpz_init(gmp_help4);
  mpz_init(gmp_approx);
  mpz_init(gmp_Na5);
  mpz_init(gmp_a5);
  mpz_init(gmp_a4);
  mpz_init(gmp_a3);
  mpz_init(gmp_a2);
  mpz_init(gmp_a1);
  mpz_init(gmp_a0);
  mpz_init(gmp_m0);

  mpz_init(gmp_root);
  mpz_set_ui(gmp_root,1); /* important for first root-computation */

  mpz_fdiv_q_ui(gmp_a5,gmp_a5_begin,MULTIPLIER);
  mpz_mul_ui(gmp_a5,gmp_a5,MULTIPLIER);
  if (mpz_cmp(gmp_a5,gmp_a5_begin)<0) mpz_add_ui(gmp_a5,gmp_a5,MULTIPLIER);

  dbl_N=mpz_get_d(gmp_N);
  dbl_a5=mpz_get_d(gmp_a5);
  skewness_min=pow(dbl_N/dbl_a5/pow(norm_max,5.),2./15.);
  p_size_max=norm_max/sqrt(skewness_min)/skewness_min; /* we need this in find_primes */

  dbl_a5_max=pow(norm_max,8.)/dbl_N;
  dbl_a5_max=sqrt(dbl_a5_max);
  dbl_a5_min=0.;
  for (i=0; i<npr_in_p; i++) dbl_a5_min+=log((double)pr_mod5[i]);
  dbl_a5_min*=5.;
  dbl_a5_min-=10.*log(norm_max);
  dbl_a5_min+=log(dbl_N);
  dbl_a5_min=exp(dbl_a5_min);
  dbl_a5_min=floor(dbl_a5_min)+1.;

  printf("Parameters: number: %.3e, norm_max: %.3e, a5: %.3e - %.3e\n",dbl_N,norm_max,dbl_a5,mpz_get_d(gmp_a5_end));
  printf("Accected a5-range: %f - %f \n",dbl_a5_min,dbl_a5_max);

  if (dbl_a5_min>mpz_get_d(gmp_a5))
    complain("a5_begin smaller than %f\n",dbl_a5_min);
  if (dbl_a5_max<mpz_get_d(gmp_a5_end))
    complain("a5_end bigger than %f\n",dbl_a5_max);

  printf("Searching in subinterval: ");
  mpz_out_str(stdout,10,gmp_a5); printf(" - ");
  mpz_out_str(stdout,10,gmp_a5_end); printf("\n");


  mpz_init(gmp_prod);
  mpz_init(gmp_prod_1);
  mpz_init(gmp_d);
  mpz_init(gmp_disp);
  mpz_init(gmp_D0);

  for (i=0; i<NPR5; i++) for (j=0; j<5; j++) mpz_init(gmp_D[i][j]);
  mpz_init(gmp_kappa0);
  for (i=0; i<NPR5; i++) for (j=0; j<5; j++) mpz_init(gmp_kappa[i][j]);
  for (i=0; i<NPR5+1; i++) mpz_init(gmp_kappa_help[i]);
  stored_pairs=NULL; store_len=0;
  raw_store_len=0; raw_stored_pairs=NULL;


  len1=npr_in_p/2; len2=npr_in_p-len1;  /* len1<=len2 */
  len11=1; if (len1>3) len11=2;
  len12=len1-len11;
  len21=1; if (len2>3) len21=2;
  len22=len2-len21;
  s11len=five_pow(len11); s12len=five_pow(len12);
  s21len=five_pow(len21); s22len=five_pow(len22);
  s1len=s11len*s12len; s2len=s21len*s22len;
  shelp=(uint64_t *)xmalloc(s22len*sizeof(uint64_t)); /* s22len is maximum */

  s11l=(uint *)xmalloc(s11len*sizeof(uint));
  s12l=(uint *)xmalloc(s12len*sizeof(uint));
  s21l=(uint *)xmalloc(s21len*sizeof(uint));
  s22l=(uint *)xmalloc(s22len*sizeof(uint));
  s1sortl=(uint *)xmalloc(s1len*sizeof(uint));
  s2sortl=(uint *)xmalloc(s2len*sizeof(uint));
  s1=(uint64_t *)xmalloc(s1len*sizeof(uint64_t));
  s2=(uint64_t *)xmalloc(s2len*sizeof(uint64_t));
  s1sort=(uint64_t *)xmalloc(s1len*sizeof(uint64_t));
  s2sort=(uint64_t *)xmalloc(s2len*sizeof(uint64_t));
  raw_cand=(int *)xmalloc(s2len*sizeof(int));
  raw_cand_hash=(uint *)xmalloc(s2len*sizeof(uint));

  hashpart_shift=1;
  while (s22len/(1<<hashpart_shift)>64) hashpart_shift++;  /* 32 or 64 or 128 seem to be best on PIII */
  hashpart_shift=32-hashpart_shift;
  hash_shift=hashpart_shift-(64-HASHSHIFT);
  n_hash_parts=1<<(32-hashpart_shift);
  if (verbose>1) printf("number of hashparts: %u\n",n_hash_parts);

  shelpsort=(uint64_t *)xmalloc(2*s22len/DEG*sizeof(uint64_t));
  indhelpsort=(int *)xmalloc(2*s22len/DEG*sizeof(int)); /* s22len>=s12len */
  indhelp=(int *)xmalloc(s22len*sizeof(int));
  s11_begin=(uint *)xmalloc(s11len*sizeof(uint));
  s21_begin=(uint *)xmalloc(s21len*sizeof(uint));

  s12l_sort_maxlen=s12len+(s12len/(1<<(31-hashpart_shift)));
  s22l_sort_maxlen=s22len+(s22len/(1<<(31-hashpart_shift)));
  s12l_sort_len=0; s22l_sort_len=0;
  s12l_sort=(uint *)xmalloc(s12l_sort_maxlen*sizeof(uint));
  s12l_ind=(uint *)xmalloc(s12l_sort_maxlen*sizeof(uint));
  s22l_sort=(uint *)xmalloc(s22l_sort_maxlen*sizeof(uint));
  s22l_ind=(uint *)xmalloc(s22l_sort_maxlen*sizeof(uint));

  mpz_sub_ui(gmp_a5,gmp_a5,MULTIPLIER);
  for (i=0; i<npr_mod5; i++) {
    modulo32=pr_mod5[i];
    N_mod_pr[i]=mpz_fdiv_ui(gmp_N,pr_mod5[i]);
    N_inv_mod_pr[i]=invert(N_mod_pr[i],pr_mod5[i]);
    pr_step[i]=modmul32(MULTIPLIER,N_inv_mod_pr[i]);
    a5_mod_pr[i]=mpz_fdiv_ui(gmp_a5,pr_mod5[i]);
    pr_start[i]=modmul32(a5_mod_pr[i],N_inv_mod_pr[i]);
  }

  len=0;
  for (i=0; i<npr_mod5; i++) len+=(pr_mod5[i]/8+1);
  bitarray_5power[0]=(uchar *)xmalloc(len*sizeof(uchar));
  memset(bitarray_5power[0],0,len*sizeof(uchar));
  for (i=1; i<npr_mod5; i++)
    bitarray_5power[i]=bitarray_5power[i-1]+(pr_mod5[i-1]/8+1);
  for (i=0; i<npr_mod5; i++)
    for (j=1; j<pr_mod5[i]; j++)
      if (is_5power(j,pr_mod5[i]))
        bitarray_5power[i][j>>3]|=ucmask[j&7];

  p0_list_len_max=256;
  p0_list=(uint *)xmalloc(p0_list_len_max*sizeof(uint));
  p0_root=(uint *)xmalloc(p0_list_len_max*sizeof(uint));
}


/* main loop */
void search_a5()
{
  int i, j, ind;
  double p_size, dp0_max;
  uint p0, p0_max, r0;
  int passed;
  uint shift=0;

  while (1) { stat_n_a5++;
    shift++;
    if (shift>=1000) {
 /* update some variables, not very precise at the beginning of a5-range */
      mpz_set_ui(gmp_help1,shift);
      mpz_mul_ui(gmp_help1,gmp_help1,MULTIPLIER);
      mpz_add(gmp_a5,gmp_a5,gmp_help1);
      if (mpz_cmp(gmp_a5,gmp_a5_end)>0) break;
      dbl_a5=mpz_get_d(gmp_a5);
      skewness_min=pow(dbl_N/dbl_a5/pow(norm_max,5.),2./15.);
      a3_max=norm_max/sqrt(skewness_min);
      p_size_max=a3_max/skewness_min;
      p_size_max=log(p_size_max);
      shift=0;
    }
    if (find_primes()) { stat_n_pr++;
#ifdef ZEIT
zeita(1);
#endif
      compute_root(shift); shift=0;
#ifdef ZEIT
zeitb(1);
#endif
      if (mpz_cmp(gmp_a5,gmp_a5_end)>0) break;
      if (verbose>3) { mpz_out_str(stdout,10,gmp_a5); printf(" \n"); }
#ifdef ZEIT
zeitA(16);
#endif
      init_all_pr();
#ifdef ZEIT
zeitB(16);
#endif
      p_size=0;
      dbl_a5=mpz_get_d(gmp_a5);
      skewness_max=pow(norm_max/dbl_a5,0.4);
      skewness_min=pow(dbl_N/dbl_a5/pow(norm_max,5.),2./15.);
      a3_max=norm_max/sqrt(skewness_min);
      p_size=0;
      p_size_max=a3_max/skewness_min;
      p_size_max=log(p_size_max);
      for (i=0; i<npr_in_p; i++) { p_ind[i]=i; p_size+=p_log[i]; }
#ifdef ZEIT
zeitA(20);
#endif
      init_aux_factors(exp(p_size_max-p_size));
#ifdef ZEIT
zeitB(20);
#endif
      while (1) {
#ifdef ZEIT
zeitA(2);
#endif
        dp0_max=exp(p_size_max-p_size);
        if (dp0_max>(double)(p0_limit)) p0_max=p0_limit;
        else p0_max=(uint)(rint(dp0_max)-1);
        if (p0_max<1) p0_max=1;
        p0=1;
        for (ind=0; ind<p0_list_len; ind++) {
          last_p0=p0;
          p0=p0_list[ind]; r0=p0_root[ind];
          if (p0>p0_max) break;

          if (verbose>3) {
            printf("search for: ");
            for (i=0; i<npr_in_p; i++) printf("%d ",p_pr[p_ind[i]]);
            if (p0!=1) printf("   %u (%u)",p0,r0);
            printf("\n");
          }
stat_n_p++;
#ifdef ZEIT
zeitA(3);
#endif
          if (p0==1) init_knapsack_raw();
          else { stat_n_p_p0++; init_knapsack_raw_p0(p0,r0); }
#ifdef ZEIT
zeitB(3); zeitA(4);
#endif
          passed=knapsack_raw();
#ifdef ZEIT
zeitB(4);
#endif
          if (passed>0) {
#ifdef ZEIT
zeitA(15);
#endif
            if (p0==1) init_knapsack_exact();
            else init_knapsack_exact_p0(p0,r0);
#ifdef ZEIT
zeitB(15); zeitA(6);
#endif
            check_raw();
#ifdef ZEIT
zeitB(6);
#endif
          } else if (passed!=0) {
#ifdef ZEIT
zeitA(15);
#endif
            if (p0==1) init_knapsack_exact();
            else init_knapsack_exact_p0(p0,r0);
#ifdef ZEIT
zeitB(15); zeitA(17);
#endif
            knapsack_exact();
#ifdef ZEIT
zeitB(17);
#endif
          }
        }
#ifdef ZEIT
zeitB(2); zeitA(5);
#endif
        i=npr_in_p-1;
        while (i>=0) {
          if (p_size>p_size_max) {
            i--;
            while ((i>=0) && (p_ind[i]+1==p_ind[i+1])) i--;
            if (i<0) break;
            p_ind[i]++; for (j=i+1; j<npr_in_p; j++) p_ind[j]=p_ind[j-1]+1;
            p_size=0; for (j=0; j<npr_in_p; j++) p_size+=p_log[p_ind[j]];
            if (p_size<=p_size_max) break;
            continue;
          }
          while ((i>=0) && (p_ind[i]>=i+npr_excess)) i--;
          if (i<0) break;
          p_ind[i]++; for (j=i+1; j<npr_in_p; j++) p_ind[j]=p_ind[j-1]+1;
          p_size=0; for (j=0; j<npr_in_p; j++) p_size+=p_log[p_ind[j]];
          if (p_size<=p_size_max) break;
        }
        if (i<0) break;
#ifdef ZEIT
zeitB(5);
#endif
      }
#ifdef ZEIT
zeitB(5);
#endif
    }
  }
}




int main(int argc, char **argv)
{
#ifdef ZEIT
initzeit(25); zeita(0);
#endif
  printf("%s\n", START_MESSAGE);

  setbuf(stdout,NULL);
  get_options(argc,argv);
  read_data();
  init_search();
  search_a5();
  if (verbose) {
    printf("Statistics:\n");
    printf("a5-values: %" PRIu64 ", suitable for %d primes: %" PRIu64 "\n",stat_n_a5,npr_in_p,stat_n_pr);
    printf("raw checks of (a5,p): %" PRIu64 " (%" PRIu64 "),  fine checks of (a5,p): %" PRIu64 "\n",stat_n_p,stat_n_p_p0,stat_n_raw);
    printf("polynomials computed: %" PRIu64 ",  survivors: %" PRIu64 "\n",stat_n_polexpand,stat_n_survivors);
    printf("Total number of checked polynomials: %" PRIu64 "\n",(uint64_t)(five_pow(npr_in_p))*stat_n_p);
  }
  if (success) printf("success\n");
#ifdef ZEIT
zeitb(0);
  printf("\nTiming:");
  printf("\nTotal            "); printzeit(0);
  printf("\n  Root             "); printzeit(1);
  printf("\n  Init primes      "); printzeit(16);
  printf("\n  Knapsack         "); printzeit(2);
  printf("\n    Init knaps. raw  "); printzeit(3);
  printf("\n      D's              "); printzeit(7);
  printf("\n      kappa's          "); printzeit(8); printzeit(18);
  printf("\n      rest             "); printzeit(9);
  printf("\n    Knapsack raw     "); printzeit(4);
  printf("\n      combine          "); printzeit(10);
  printf("\n      hash1/hash2      "); printzeit(11); printzeit(12);
  printf("\n      hash3/hash4      "); printzeit(13); printzeit(14);
  printf("\n    Init knaps. ex.  "); printzeit(15);
  printf("\n    Knapsack ex.       "); printzeit(17);
  printf("\n    Check pol        "); printzeit(6);
  printf("\n  Next p           "); printzeit(5);
  printf("\n");
#endif
  return 0;
}

