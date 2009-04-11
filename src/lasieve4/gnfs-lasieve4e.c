/* gnfs-lasieve4e.c
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

#include "lasieve.h"

#include <assert.h>
#include <stdio.h>
#include <sys/types.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <ctype.h>

#if defined (_MSC_VER) && defined (_DEBUG)
	#include <crtdbg.h>
#endif

#if defined (_MSC_VER) || defined (__MINGW32__) || defined (MINGW32)
#include "getopt.h"
#define popen  _popen
#define pclose _pclose
#define bzero(x,n)	memset(x,0,n)
#else
#include <fnmatch.h>
#endif
#ifdef LINUX
#include <endian.h>
#endif
#if !defined (_MSC_VER)
#include <sys/time.h>
#endif

#include <gmp.h>
#include <signal.h>
#include <setjmp.h>

#ifdef GGNFS_HOST_GENERIC
const u32_t schedule_primebounds[N_PRIMEBOUNDS]={0x100000,0x200000,0x400000,0x800000,0x1000000,0x2000000,UINT_MAX};
const u32_t schedule_sizebits[N_PRIMEBOUNDS]={20,21,22,23,24,25,32};
#endif

#define GCD_SIEVE_BOUND 10
#define MAX_TINY_2POW 4
#define L1_SIZE (1UL<<L1_BITS)
#define CANDIDATE_SEARCH_STEPS 128
#define MAX_LPFACTORS 3
#define RI_SIZE 2
#define SE_SIZE 2
#define SCHEDFBI_MAXSTEP 0x10000
#define USE_MEDSCHED
#define TINY_SIEVEBUFFER_SIZE 420
#define TINY_SIEVE_MIN 8
#define ALGEBRAIC_SIDE 0
#define RATIONAL_SIDE 1
#define NO_SIDE 2
#define USER_INTERRUPT 1
#define SCHED_PATHOLOGY 2

static float FB_bound[2], sieve_report_multiplier[2];
static size_t sieve_min[2], max_primebits[2], max_factorbits[2];
static u32_t *(FB[2]), *(proots[2]), FBsize[2];
static double *(tpoly_f[2]), *(tpoly_fr[2]);
static i16_t *rroots1, *rroots2;
static u16_t n_rroots1, n_rroots2;
static size_t rroots1_alloc = 0, rroots2_alloc = 0;
static unsigned char *plog_lb1, *plog_lb2;

u16_t get_plog_lb(unsigned char *,i16_t **,size_t *,double *,i32_t,double);

static u32_t first_spq, last_spq, sieve_count;
u32_t all_spq_done;
static mpz_t g_m, g_N, g_aux1, g_aux2, g_aux3, g_sr_a, g_sr_b;
static mpz_t *(g_poly[2]);
double *(poly_f[2]), poly_norm[2];
u32_t  g_poldeg[2], poldeg_max;
u32_t  keep_factorbase;
u32_t  g_resume;

static mpz_t rational_rest, algebraic_rest;
mpz_t factors[MAX_LPFACTORS];
static u32_t yield=0, n_mpqsfail[2]={0,0}, n_mpqsvain[2]={0,0};
static u32_t mpqs_clock=0;
static clock_t sieve_clock=0, sch_clock=0, td_clock=0, tdi_clock=0;
static u32_t cs_clock[2]={0,0}, Schedule_clock=0, medsched_clock=0;
static u32_t si_clock[2]={0,0}, s1_clock[2]={0,0};
static u32_t s2_clock[2]={0,0}, s3_clock[2]={0,0};
static u32_t tdsi_clock[2]={0,0}, tds1_clock[2]={0,0}, tds2_clock[2]={0,0};
static u32_t tds3_clock[2]={0,0}, tds4_clock[2]={0,0};
char  *base_name;
char  *input_line=NULL;
size_t input_line_alloc=0;
static u32_t ncand;
static u16_t *cand;
static unsigned char *fss_sv;
u32_t process_no;
char *sysload_cmd;
double sieveStartTime;
static u16_t short_output = 0;

int rho_factor(unsigned long *factors, mpz_t n);

/****************************************************/
static int tdcand_cmp(const void *x, const void *y)
/****************************************************/
{ return (int) (*((u16_t *) x)) - (int) (*((u16_t *) y)); 
}

typedef struct xFBstruct {
  u32_t p, pp, q, qq, r, l;
} *xFBptr;

static volatile xFBptr xFB[2];
static volatile u32_t xFBs[2];
static void xFBtranslate(u16_t *rop, xFBptr op);
static int xFBcmp(const void *, const void *);
static u32_t add_primepowers2xaFB(size_t * aFB_alloc_ptr,
                                  u32_t pp_bound, u32_t side, u32_t p,
                                  u32_t r);

i32_t a0, a1, b0, b1;

#if 0
u32_t I_bits;
#endif
u32_t J_bits, i_shift, n_I, n_J;
u32_t root_no;
float sigma;
static u32_t oddness_type;
static u32_t n_i, n_j, i_bits, j_bits;
u32_t spq_i, spq_j, spq_x;
u32_t fbi1[2];
u32_t fbis[2];
static u32_t j_per_strip, jps_bits, jps_mask, n_strips;
u32_t n_spq, n_spq_discard;
static pr32_struct special_q_ps;
u32_t special_q;
double special_q_log;

static struct schedule_struct {
  u16_t ***schedule;
  u32_t *fbi_bounds;
  u32_t n_pieces;
  unsigned char *schedlogs;
  u16_t n_strips, current_strip;
  size_t alloc, alloc1;
  u32_t *ri;
} *(schedules[2]);

u32_t n_schedules[2];
static u32_t *(LPri[2]), *(LPri1[2]);


static u32_t *(current_ij[2]);


#ifdef USE_MEDSCHED
static u16_t **(med_sched[2]);
static u32_t *(medsched_fbi_bounds[2]);
static unsigned char *(medsched_logs[2]);
static size_t medsched_alloc[2];
static u16_t n_medsched_pieces[2];
#endif

static unsigned char *sieve_interval = NULL, *(FB_logs[2]);
static unsigned char *tiny_sieve_buffer;

static double sieve_multiplier[2], FB_maxlog[2];

void   do_scheduling(struct schedule_struct *, u32_t, u32_t, u32_t);
static u16_t *(smallsieve_aux[2]), *(smallsieve_auxbound[2][5]);
static u16_t *(smallsieve_tinybound[2]);
static u16_t *(smallsieve_aux1[2]), *(smallsieve_aux1_ub_odd[2]);
static u16_t *(smallsieve_aux1_ub[2]), *(smallsieve_tinybound1[2]);
static u16_t *(smallsieve_aux2[2]), *(smallsieve_aux2_ub[2]);
static u16_t *(smallpsieve_aux[2]), *(smallpsieve_aux_ub_pow1[2]);
static u16_t *(smallpsieve_aux_ub_odd[2]), *(smallpsieve_aux_ub[2]);
static unsigned char *horizontal_sievesums;
static u16_t *(x2FB[2]), x2FBs[2];
static u16_t **(smalltdsieve_aux[2]);

#ifdef PREINVERT
static u32_t *(smalltd_pi[2]);
#endif

#ifdef GCD_SIEVE_BOUND
static u32_t np_gcd_sieve;
static unsigned char *gcd_sieve_buffer;
static void gcd_sieve(void);
#endif

u16_t **schedbuf;


#ifdef OFMT_CWI
static char ulong2cwi(u32_t);
#endif

void dumpsieve(u32_t j_offset, u32_t side);

u32_t *(td_buf[2]), **td_buf1;
size_t td_buf_alloc[2] = { 1024, 1024 };

static unsigned char tds_coll[UCHAR_MAX];
u32_t **tds_fbi = NULL;
u32_t **tds_fbi_curpos = NULL;

#ifndef TDFBI_ALLOC
#define TDFBI_ALLOC 256
static size_t tds_fbi_alloc = TDFBI_ALLOC;
#endif

static mpz_t td_rests[L1_SIZE];
static mpz_t large_factors[2], *(large_primes[2]);
static mpz_t FBb_sq[2];

/**************************************************/
void Usage()
/**************************************************/
{
  complain("Usage");
}

static u32_t n_prereports = 0, n_reports = 0, n_rep1 = 0, n_rep2 = 0;
static u32_t n_tdsurvivors[2] = { 0, 0 };
static FILE *g_ofile;
static char *g_ofile_name;

#ifdef STC_DEBUG
FILE *debugfile;
#endif

static u16_t special_q_side, first_td_side, first_sieve_side;
static u16_t first_psp_side, first_mpqs_side;
static u16_t cmdline_first_sieve_side = USHRT_MAX;
static u16_t cmdline_first_td_side = USHRT_MAX;


jmp_buf termination_jb;

static void terminate_sieving(int signo)
{
  longjmp(termination_jb, USER_INTERRUPT);
}

/******************************************************/
void getFB(int force_aFBcalc)
/******************************************************/
{ size_t FBS_alloc = 4096;
  u32_t  prime;
  pr32_struct ps;
  char  *afbname;
  FILE  *afbfile;
  u32_t  side, i, j, k, l=0, nr, x;
  u32_t  srfbs, safbs;
  double ld;

  initprime32(&ps);

  for (side = 0; side < 2; side++) {
    if (g_poldeg[side] == 1) {
      FB[side] = xmalloc(FBS_alloc * sizeof(u32_t));
      proots[side] = xmalloc(FBS_alloc * sizeof(u32_t));
      prime = firstprime32(&ps);
      for (prime = nextprime32(&ps), fbi1[side] = 0, FBsize[side] = 0;
           prime < FB_bound[side]; prime = nextprime32(&ps)) {
        x = mpz_fdiv_ui(g_poly[side][1], prime);
        if (x > 0) {
          modulo32 = prime;
          x = modmul32(modinv32(x), mpz_fdiv_ui(g_poly[side][0], prime));
          x = x > 0 ? prime - x : 0;
        } else
          x = prime;
        if (prime < L1_SIZE)
          fbi1[side] = FBsize[side];
        if (prime < n_i)
          fbis[side] = FBsize[side];
        if (FBsize[side] == FBS_alloc) {
          FBS_alloc *= 2;
          FB[side] = xrealloc(FB[side], FBS_alloc * sizeof(u32_t));
          proots[side] = xrealloc(proots[side], FBS_alloc * sizeof(u32_t));
        }
        proots[side][FBsize[side]] = x;
        FB[side][FBsize[side]++] = prime;
      }
      proots[side] = xrealloc(proots[side], FBsize[side] * sizeof(u32_t));
      FB[side] = xrealloc(FB[side], FBsize[side] * sizeof(u32_t));
      fbi1[side]++;
      fbis[side]++;
      if (fbi1[side] < fbis[side])
        fbi1[side] = fbis[side];
    } else {
      asprintf(&afbname, "%s.afb.%u", base_name, side);
      if (force_aFBcalc > 0 || (afbfile = fopen(afbname, "rb")) == NULL) {
        u32_t *root_buffer;
        size_t aFB_alloc;

        root_buffer = xmalloc(g_poldeg[side] * sizeof(*root_buffer));
        aFB_alloc = 4096;
        FB[side] = xmalloc(aFB_alloc * sizeof(**FB));
        proots[side] = xmalloc(aFB_alloc * sizeof(**proots));
        for (prime = firstprime32(&ps), FBsize[side] = 0;
             prime < FB_bound[side]; prime = nextprime32(&ps)) {

          nr = root_finder(root_buffer, g_poly[side], g_poldeg[side], prime);
          for (i = 0; i < nr; i++) {
            if (aFB_alloc <= FBsize[side]) {
              aFB_alloc *= 2;
              FB[side] = xrealloc(FB[side], aFB_alloc * sizeof(**FB));
              proots[side] =
                xrealloc(proots[side], aFB_alloc * sizeof(**proots));
            }
            FB[side][FBsize[side]] = prime;
            proots[side][FBsize[side]] = root_buffer[i];
            if (prime > 2)
              FBsize[side]++;
          }
        }
        FB[side] = xrealloc(FB[side], FBsize[side] * sizeof(**FB));
        proots[side] =
          xrealloc(proots[side], FBsize[side] * sizeof(**proots));
        free(root_buffer);

        if (keep_factorbase > 0) {
          if ((afbfile = fopen(afbname, "wb")) == NULL) {
            complain("Cannot open %s for output of aFB: %m\n", afbname);
          }
          write_u32(afbfile, &(FBsize[side]), 1);
          write_u32(afbfile, FB[side], FBsize[side]);
          write_u32(afbfile, proots[side], FBsize[side]);
          write_u32(afbfile, (void*)&(xFBs[side]), 1);
          fclose(afbfile);
        }

      } else {
        if (read_u32(afbfile, &(FBsize[side]), 1) != 1) {
          complain("Cannot read aFB size from %s: %m\n", afbname);
        }
        FB[side] = xmalloc(FBsize[side] * sizeof(u32_t));
        proots[side] = xmalloc(FBsize[side] * sizeof(u32_t));
        if (read_u32(afbfile, FB[side], FBsize[side]) != FBsize[side] ||
            read_u32(afbfile, proots[side], FBsize[side]) != FBsize[side]) {
          complain("Cannot read aFB from %s: %m\n", afbname);
        }
        if (read_u32(afbfile, (void*)&(xFBs[side]), 1) != 1) {
          complain("%s: Cannot read xFBsize\n", afbname);
        }
        fclose(afbfile);
      }

      for (j = 0, k = 0, l = 0; j < FBsize[side]; j++) {
        if (FB[side][j] < L1_SIZE)
          k = j;
        if (FB[side][j] < n_i)
          l = j;
        if (FB[side][j] > L1_SIZE && FB[side][j] > n_I)
          break;
      }
      if (FBsize[side] > 0) {
        if (k < l)
          k = l;
        fbis[side] = l + 1;
        fbi1[side] = k + 1;
      } else {
        fbis[side] = 0;
        fbi1[side] = 0;
      }
    }
  }

  for (i = 0, srfbs = 0; i < xFBs[1]; i++) {
    if (xFB[1][i].p == xFB[1][i].pp)
      srfbs++;
  }
  for (i = 0, safbs = 0; i < xFBs[0]; i++) {
    if (xFB[0][i].p == xFB[0][i].pp)
      safbs++;
  }
  logbook(0, "FBsize %u+%u (deg %u), %u+%u (deg %u)\n",
          FBsize[0], safbs, g_poldeg[0], FBsize[1], srfbs, g_poldeg[1]);
  free(afbname);

  for (i = 0; i < 2; i++) {
    tpoly_f[i] = xmalloc((1 + g_poldeg[i]) * sizeof(**tpoly_f));
    tpoly_fr[i] = xmalloc((1 + g_poldeg[i]) * sizeof(**tpoly_fr));
  }
  plog_lb1 = xmalloc(2 * n_J / CANDIDATE_SEARCH_STEPS);
  plog_lb2 = xmalloc(2 * n_J / CANDIDATE_SEARCH_STEPS);

  if (sieve_count == 0)
    exit(0);

  for (side = 0; side < 2; side++) {
    struct xFBstruct *s;
    u32_t *root_buffer;
    size_t xaFB_alloc = 0;

    FB_logs[side] = xmalloc(FBsize[side]);
    if (g_poldeg[side] == 1)
      sieve_multiplier[side] = (UCHAR_MAX-50)/log(n_J*fabs(mpz_get_d(g_poly[side][0])));
    else
      sieve_multiplier[side] = (UCHAR_MAX-50)/log(poly_norm[side]);

    root_buffer = xmalloc(g_poldeg[side]*sizeof(*root_buffer));
    prime = 2;
    nr = root_finder(root_buffer, g_poly[side], g_poldeg[side], prime);

    for (i = 0; i < nr; i++) {
      adjust_bufsize((void **) &(xFB[side]), &xaFB_alloc, 1+xFBs[side],16,sizeof(**xFB));
      s = xFB[side] + xFBs[side];
      s->p = prime;
      s->pp = prime;
      if (root_buffer[i] == prime) {
        s->qq = prime;
        s->q = 1;
        s->r = 1;
      } else {
        s->qq = 1;
        s->q = prime;
        s->r = root_buffer[i];
      }
      xFBs[side]++;
      add_primepowers2xaFB(&xaFB_alloc, n_I, side, 0, 0);
    }
    free(root_buffer);
    for (i = 0; i < FBsize[side]; i++) {
      u32_t l1;

      prime = FB[side][i];
      if (prime > n_I / prime)
        break;
      ld = log(prime);
      l1 = add_primepowers2xaFB(&xaFB_alloc, n_I, side, prime, proots[side][i]);
	  assert(rint(l1 * ld * sieve_multiplier[side]) <= UCHAR_MAX);
      FB_logs[side][i] = (unsigned char)rint(l1 * ld * sieve_multiplier[side]);
    }
    while (i < FBsize[side]) {
      ld = log(FB[side][i]);
      if (l > FB_maxlog[side])
        FB_maxlog[side] = ld;
	  assert(rint(sieve_multiplier[side] * ld) <= UCHAR_MAX);
      FB_logs[side][i++] = (unsigned char)rint(sieve_multiplier[side] * ld);
    }
    FB_maxlog[side] *= sieve_multiplier[side];
    qsort(xFB[side], xFBs[side], sizeof(*(xFB[side])), xFBcmp);
  }
}
/*******************************************************/
double sTime()
/*******************************************************/
#if !defined (_MSC_VER) && !defined (__MINGW32__) && !defined (MINGW32)
{ static struct  timeval  this_tv;
  static struct  timezone dumbTZ;
  double t;

  gettimeofday(&this_tv, &dumbTZ);
  t = this_tv.tv_sec + 0.000001*this_tv.tv_usec;
  return t;
}
#else
{
        return clock() / (double)CLOCKS_PER_SEC;
}
#endif

typedef unsigned long bc_t;
#define BC_ONES ((~0UL)/0xFFU)
#define BC_MASK (BC_ONES*0x80U)
inline void optsieve(uint32_t st1, uchar* i_o, uchar* i_max, size_t j) {
  // align i_o & i_max to 32-byte boundary
  for(;i_o<i_max && ((size_t)i_o & 0x1F);++i_o) {
    if (*i_o >= st1) {
	  assert((i_o - sieve_interval) <= USHRT_MAX);
      cand[ncand] = (u16_t)(i_o - sieve_interval);
      fss_sv[ncand++] = *i_o + horizontal_sievesums[j];
    }
  }
  while(i_o<i_max && ((size_t)i_max & 0x1F)) {
    if (*--i_max >= st1) {
	  assert((i_max - sieve_interval) <= USHRT_MAX);
      cand[ncand] = (u16_t)(i_max - sieve_interval);
      fss_sv[ncand++] = *i_max + horizontal_sievesums[j];
    }
  }
#ifndef HAVE_SSIMD
  if (st1 < 0x80) {
    bc_t bc, *i_oo;

    bc = BC_MASK - (BC_ONES * st1);
    for (i_oo = (bc_t *) i_o; i_oo < (bc_t *) i_max;
	 i_oo++) {
      bc_t v = *i_oo;

      if (((v & BC_MASK) | ((v + bc) & BC_MASK)) == 0) continue;
      for (i_o=(uchar *)i_oo; i_o<(uchar *)(i_oo+1); i_o++) {
	if (*i_o >= st1) {
	  assert((i_o - sieve_interval) <= USHRT_MAX);
	  cand[ncand] = (u16_t)(i_o - sieve_interval);
	  fss_sv[ncand++] = *i_o + horizontal_sievesums[j];
	}
      }
    }
  } else {
    bc_t *i_oo;

    for (i_oo=(bc_t *)i_o; i_oo<(bc_t*)i_max; i_oo++) {
      if ((*i_oo & BC_MASK) == 0) continue;
      for (i_o=(uchar *)i_oo; i_o<(uchar *)(i_oo+1); i_o++) {
	if (*i_o >= st1) {
	  assert((i_o - sieve_interval) <= USHRT_MAX);
	  cand[ncand] = (u16_t)(i_o - sieve_interval);
	  fss_sv[ncand++] = *i_o + horizontal_sievesums[j];
	}
      }
    }
  }

#else

  {
  uint64_t x;

  x = st1 - 1;
  x |= x << 8;
  x |= x << 16;
  x |= x << 32;
  while (i_o < i_max) {
#if defined(GGNFS_x86_32_MSCASM_MMX)
  __asm {
	  mov    	esi,[i_o]
	  mov    	edi,[i_max]
	  lea    	eax,[x]
	  movq	        mm7,[eax]
  l1:     movq	        mm1,[esi]
	  movq	        mm0,[esi+8]
	  pmaxub	mm1,[esi+16]
	  pmaxub	mm0,[esi+24]
	  pmaxub	mm1,mm7
	  pmaxub	mm0,mm1
	  pcmpeqb       mm0,mm7
	  pmovmskb      eax,mm0
	  cmp		eax,255
	  jnz		l2
	  lea		esi,[esi+32]
	  cmp		edi,esi
	  ja		l1
  l2:     mov		[i_o],esi
	  emms
  }
#elif defined(GGNFS_x86_64_ATTASM)
  asm volatile (
    "movq     (%%rax),%%mm7\n"
    ".align 32\n"
    "1:\n"
    "movq     (%%rsi),%%mm1\n"
    "movq     8(%%rsi),%%mm0\n"
    "pmaxub   16(%%rsi),%%mm1\n"
    "pmaxub   24(%%rsi),%%mm0\n"
    "pmaxub   %%mm7,%%mm1\n"
    "pmaxub   %%mm1,%%mm0\n"
    "pcmpeqb  %%mm7,%%mm0\n"
    "pmovmskb %%mm0,%%rax\n"
    "cmpq     $255,%%rax\n"
    "jnz      2f\n"
    "leaq     32(%%rsi),%%rsi\n"
    "cmpq     %%rsi,%%rdi\n"
    "ja       1b\n"
    "2:\n"
    "emms":"=S" (i_o):"a"(&x),
    "S"(i_o), "D"(i_max)
  );
#elif defined(GGNFS_x86_32_ATTASM_MMX)
  asm volatile (
    "movq     (%%eax),%%mm7\n"
    "1:\n"
    "movq     (%%esi),%%mm1\n"
    "movq     8(%%esi),%%mm0\n"
    "pmaxub   16(%%esi),%%mm1\n"
    "pmaxub   24(%%esi),%%mm0\n"
    "pmaxub   %%mm7,%%mm1\n"
    "pmaxub   %%mm1,%%mm0\n"
    "pcmpeqb  %%mm7,%%mm0\n"
    "pmovmskb %%mm0,%%eax\n"
    "cmpl     $255,%%eax\n"
    "jnz      2f\n"
    "leal     32(%%esi),%%esi\n"
    "cmpl     %%esi,%%edi\n"
    "ja       1b\n"
    "2:\n"
    "emms":"=S" (i_o):"a"(&x),
    "S"(i_o), "D"(i_max)
   );
#else
    #error Unsupported assembler model!
#endif
    if (i_o < i_max) {
      uchar *i_max2 = i_o + 32;

      while (i_o < i_max2) {
	if (*i_o >= st1) {
	  assert((i_o - sieve_interval) <= USHRT_MAX);
	  cand[ncand] = (u16_t)(i_o - sieve_interval);
	  fss_sv[ncand++] = *i_o + horizontal_sievesums[j];
	}
	i_o++;
      }
    }
  }
  }
#endif
}

/******************************************************************/
int lasieve()
/******************************************************************/
{ clock_t last_clock;
  u32_t  *r_ptr, s, i;
  u32_t   subsieve_nr, j_offset;
  u32_t   absa0, absa1, absb0, absb1;
  char    a0s, a1s;
  double  tStart, tNow, lastReport;

  initprime32(&special_q_ps);
  last_clock = clock();
  n_spq = 0;
  n_spq_discard = 0;
  r_ptr = xmalloc(poldeg_max * sizeof(*r_ptr));
  special_q = pr32_seek(&special_q_ps, first_spq);
  tStart = lastReport = sTime();
  for ( ; special_q < last_spq && special_q != 0; special_q = nextprime32(&special_q_ps)) {
    u32_t nr;

    special_q_log = log(special_q);
    if (cmdline_first_sieve_side == USHRT_MAX) {
      double nn[2];

      for (s = 0; s < 2; s++) {
        nn[s] = log(poly_norm[s] * (special_q_side == s ? 1 : special_q));
        nn[s] = nn[s] / (sieve_report_multiplier[s] * log(FB_bound[s]));
      }
      if (nn[0] < nn[1])
        first_sieve_side = 1;
      else
        first_sieve_side = 0;
    } else {
      first_sieve_side = cmdline_first_sieve_side;
      if (first_sieve_side >= 2)
        complain("First sieve side must not be %u\n",
                 (u32_t) first_sieve_side);
    }
    logbook(1, "First sieve side: %u\n", (u32_t) first_sieve_side);
    if (cmdline_first_td_side != USHRT_MAX)
      first_td_side = cmdline_first_td_side;
    else
      first_td_side = first_sieve_side;
    if (g_poldeg[special_q_side] > 1) {
      nr = root_finder(r_ptr, g_poly[special_q_side], g_poldeg[special_q_side], special_q);
      if (nr == 0)
        continue;
      if (r_ptr[nr - 1] == special_q) {
        nr--;
      }
    } else {
      u32_t x = mpz_fdiv_ui(g_poly[special_q_side][1], special_q);

      if (x == 0) {
        n_spq_discard++;
        continue;
      }
      modulo32 = special_q;
      x = modmul32(modinv32(x), mpz_fdiv_ui(g_poly[special_q_side][0], special_q));
      r_ptr[0] = x == 0 ? 0 : special_q - x;
      nr = 1;
    }

    for (root_no = 0; root_no < nr; root_no++) {
      clock_t new_clock;
      u32_t termination_condition;

      if ((termination_condition = setjmp(termination_jb)) != 0) {
        if (termination_condition == USER_INTERRUPT) {
          char *hn, *ofn;
          FILE *of;

          hn = xmalloc(100);
#if 0
          if (gethostname(hn, 99) == 0)
            asprintf(&ofn, "%s.%s.last_spq%d", base_name, hn, process_no);
          else
            asprintf(&ofn, "%s.unknown_host.last_spq%d", base_name,
                     process_no);
          free(hn);
#else
          asprintf(&ofn, ".last_spq%d", process_no);
#endif

          if ((of = fopen(ofn, "wb")) != 0) {
            fprintf(of, "%u\n", (unsigned int)special_q);
            fclose(of);
          }
          free(ofn);
          all_spq_done = 0;
          break;
        }
        else {
          char *cmd;

          asprintf(&cmd, "touch badsched.%s.%u.%u.%u", base_name,
                   special_q_side, special_q, r_ptr[root_no]);
          system(cmd);
          free(cmd);
          continue;
        }
      }
      if (sysload_cmd != NULL) {
        if (system(sysload_cmd) != 0) {
          longjmp(termination_jb, USER_INTERRUPT);
        }
      }
      n_spq++;
      reduce2(&a0, &b0, &a1, &b1, (i32_t) special_q, 0, (i32_t) r_ptr[root_no],
              1, sigma * sigma);

      if (b0 % ((i32_t) special_q) == 0 && b1 % ((i32_t) special_q) == 0) {
        i32_t x;

        x = a0 % ((i32_t) special_q);
        if (x < 0)
          x += (i32_t) special_q;
        spq_i = x;
        x = a1 % ((i32_t) special_q);
        if (x < 0)
          x += (i32_t) special_q;
        spq_j = x;
      } else {
        i32_t x;

        x = b0 % ((i32_t) special_q);
        if (x < 0)
          x += (i32_t) special_q;
        spq_i = x;
        x = b1 % ((i32_t) special_q);
        if (x < 0)
          x += (i32_t) special_q;
        spq_j = x;
      }
      modulo32 = special_q;
      spq_x = modmul32(spq_i, i_shift);

      if (i_bits == 0) {
        logbook(1, "(%d,%d) (%d,%d)\n", a0, b0, a1, b1);
        continue;
      }

#define GET_ABSSIG(abs,sig,arg) if(arg> 0) { abs= (u32_t)arg; sig= '+';} \
      else { abs= (u32_t)(-arg); sig= '-'; }


      GET_ABSSIG(absa0, a0s, a0);
      GET_ABSSIG(absa1, a1s, a1);
      absb0 = b0;
      absb1 = b1;

      for (s = 0; s < 2; s++) {
        u32_t fbi;
        u16_t *abuf;
        u16_t *ibuf;

        abuf = smallsieve_aux[s];
        ibuf = smallpsieve_aux[s];
        for (fbi = 0; fbi < fbis[s]; fbi++) {
          u32_t aa, bb;

          modulo32 = FB[s][fbi];

          aa = absa0 % FB[s][fbi];
          if (a0s == '-' && aa != 0)
            aa = FB[s][fbi] - aa;
          bb = absb0 % FB[s][fbi];
          if (proots[s][fbi] != FB[s][fbi]) {
            u32_t x;

            x = modsub32(aa, modmul32(proots[s][fbi], bb));
            if (x != 0) {
              aa = absa1 % FB[s][fbi];
              if (a1s == '-' && aa != 0)
                aa = FB[s][fbi] - aa;
              bb = absb1 % FB[s][fbi];
              x =
                modmul32(asm_modinv32(x),
                         modsub32(modmul32(proots[s][fbi], bb), aa));
              abuf[0] = (u16_t) (FB[s][fbi]);
              abuf[1] = (u16_t) x;
              abuf[2] = (u16_t) (FB_logs[s][fbi]);
              abuf += 4;
            } else {
              ibuf[0] = (u16_t) (FB[s][fbi]);
              ibuf[1] = (u16_t) (FB_logs[s][fbi]);
              ibuf += 3;
            }
          } else {

            if (bb != 0) {
              u32_t x;

              x = modulo32 - bb;
              bb = absb1 % FB[s][fbi];
              abuf[0] = (u16_t) (FB[s][fbi]);
              abuf[1] = (u16_t) (modmul32(asm_modinv32(x), bb));
              abuf[2] = (u16_t) (FB_logs[s][fbi]);
              abuf += 4;
            } else {
              ibuf[0] = (u16_t) (FB[s][fbi]);
              ibuf[1] = (u16_t) (FB_logs[s][fbi]);
              ibuf += 3;
            }
          }
        }
        smallsieve_auxbound[s][0] = abuf;
        smallpsieve_aux_ub_pow1[s] = ibuf;
      }


      for (s = 0; s < 2; s++) {
        u16_t *buf;
        u16_t *buf2;
        u16_t *ibuf;

        buf = smallsieve_aux1[s];
        buf2 = x2FB[s];
        ibuf = smallpsieve_aux_ub_pow1[s];
        for (i = 0; i < xFBs[s]; i++) {
          if (xFB[s][i].p == 2) {
            xFBtranslate(buf2, xFB[s] + i);
            buf2 += 4;
          } else {
            xFBtranslate(buf, xFB[s] + i);
            if (buf[0] == 1) {
              ibuf[1] = xFB[s][i].l;
              ibuf[0] = xFB[s][i].pp;
              ibuf += 3;
            } else
              buf += 6;
          }
        }
		assert(((buf2 - x2FB[s]) / 4) <= USHRT_MAX);
        x2FBs[s] = (u16_t)((buf2 - x2FB[s]) / 4);
        smallpsieve_aux_ub_odd[s] = ibuf;
        smallsieve_aux1_ub_odd[s] = buf;
      }

      for (s = 0; s < 2; s++) {
        u16_t *x;

        for (i = 0, x = smallsieve_aux[s];
             x < smallsieve_auxbound[s][0]; i++, x += 4) {
          u32_t k, r, pr;

          modulo32 = *x;
          r = x[1];
          pr = r;
          for (k = 0; k < j_per_strip; k++) {
            smalltdsieve_aux[s][k][i] = r;
            r = modadd32(r, pr);
          }
#ifdef PREINVERT

          {
            u32_t pinv;

            pinv = modulo32;
            pinv = 2 * pinv - pinv * pinv * modulo32;
            pinv = 2 * pinv - pinv * pinv * modulo32;
            pinv = 2 * pinv - pinv * pinv * modulo32;
#if 0
            pinv = 2 * pinv - pinv * pinv * modulo32;
#endif
            smalltd_pi[s][i] = 2 * pinv - pinv * pinv * modulo32;
          }
#endif
        }
      }

      for (s = 0; s < 2; s++) {
        u16_t *x, *xx, k, pbound, copy_buf[6];

        k = 0;
        pbound = TINY_SIEVE_MIN;
        for (x = smallsieve_aux[s]; x < smallsieve_auxbound[s][0];
             x += 4) {
          if (*x > pbound) {
            if (k == 0)
              smallsieve_tinybound[s] = x;
            else
              smallsieve_auxbound[s][5 - k] = x;
            k++;
            if (k < 5)
              pbound = n_i / (5 - k);
            else
              break;
          }
        }
        while (k < 5)
          smallsieve_auxbound[s][5 - (k++)] = x;
        for (x = (xx = smallsieve_aux1[s]);
             x < smallsieve_aux1_ub_odd[s]; x += 6) {
          if (x[0] < TINY_SIEVE_MIN) {
            if (x != xx) {
              memcpy(copy_buf, x, 6 * sizeof(*x));
              memcpy(x, xx, 6 * sizeof(*x));
              memcpy(xx, copy_buf, 6 * sizeof(*x));
            }
            xx += 6;
          }
        }
        smallsieve_tinybound1[s] = xx;
      }

      for (s = 0; s < 2; s++) {
#ifdef SCHEDULING_FUNCTION_CALCULATES_RI
        lasieve_setup(FB[s] + fbis[s], proots[s] + fbis[s],
                      fbi1[s] - fbis[s], a0, a1, b0, b1, LPri[s]);
#else
        lasieve_setup(FB[s] + fbis[s], proots[s] + fbis[s],
                      FBsize[s] - fbis[s], a0, a1, b0, b1, LPri[s]);
#endif
        LPri1[s] = LPri[s] + (fbi1[s] - fbis[s]) * RI_SIZE;
      }

      {
        u32_t k;

        for (i = 0; i < 2; i++) {
          tpol(tpoly_f[i], poly_f[i], g_poldeg[i], a0, a1, b0, b1);
        }

		if (first_sieve_side >= sizeof(tpoly_f)/sizeof(tpoly_f[0]))
		{
			assert(0);
			complain("Invalid first_sieve_side value!\n", first_sieve_side);
#ifdef _MSC_VER
			exit(1); // just to make PreFAST shutup
#endif
		}

        n_rroots1 = get_plog_lb(plog_lb1, &rroots1, &rroots1_alloc,
                                tpoly_f[first_sieve_side],
                                g_poldeg[first_sieve_side],
                                sieve_multiplier[first_sieve_side]);
        for (k = 0; k < 2; k++) {
          for (i = 0; i <= g_poldeg[k]; i++)
            tpoly_fr[k][g_poldeg[k] - i] = tpoly_f[k][i];
        }
        n_rroots2 = get_plog_lb(plog_lb2, &rroots2, &rroots2_alloc,
                                tpoly_fr[first_sieve_side],
                                g_poldeg[first_sieve_side],
                                sieve_multiplier[first_sieve_side]);
      }

      new_clock = clock();
      sch_clock += (clock_t)((1000.0 * (new_clock - last_clock)) / CLOCKS_PER_SEC);
      last_clock = new_clock;

      for (oddness_type = 1; oddness_type < 4; oddness_type++) {
        for (s = 0; s < 2; s++) {
          switch (oddness_type) {
              u16_t *x;

            case 1:
              for (x=smallsieve_aux[s]; x<smallsieve_auxbound[s][0]; x+=4) {
                u32_t p;

                p = x[0];
                x[3] = ((i_shift + p) / 2) % p;
              }
              for (x=smallsieve_aux1[s]; x<smallsieve_aux1_ub_odd[s]; x+=6) {
                u32_t p;

                p = x[0];
                x[4] = ((i_shift + p) / 2) % p;
                x[5] = 0;
              }
              for (x=smallpsieve_aux[s]; x<smallpsieve_aux_ub_odd[s]; x+=3)
                x[2] = 0;

              {
                u16_t *x, *y, *z;

                x = smallsieve_aux1_ub_odd[s];
                y = smallpsieve_aux_ub_odd[s];
                z = smallsieve_aux2[s];
                for (i = 0; i < (u32_t)(4 * x2FBs[s]); i += 4) {
                  u32_t p, pr, d, l;
                  u16_t **a;

                  d = x2FB[s][i + 1];
                  if (d == 1)
                    continue;
                  p = x2FB[s][i];
                  pr = x2FB[s][i + 2];
                  l = x2FB[s][i + 3];
                  if (p < 4) {
                    if (p == 1) {
                      *y = d / 2;
                      *(y + 2) = 0;
                    } else {
                      *y = d;
                      *(y + 2) = d / 2;
                    }
                    *(y + 1) = l;
                    y += 3;
                    continue;
                  }
                  p = p / 2;
                  if (p <= MAX_TINY_2POW)
                    a = &z;
                  else
                    a = &x;
                  **a = p;
                  *(1 + *a) = d;
                  *(2 + *a) = pr % p;
                  *(3 + *a) = l;
                  *(4 + *a) = ((i_shift + pr) / 2) % p;
                  *(5 + *a) = d / 2;
                  *a += 6;
                }
                smallsieve_aux1_ub[s] = x;
                smallpsieve_aux_ub[s] = y;
                smallsieve_aux2_ub[s] = z;
              }
              break;

            case 2:
              for (x=smallsieve_aux[s]; x<smallsieve_auxbound[s][0]; x+=4) {
                u32_t p, pr;

                p = x[0];
                pr = x[1];
                x[3] = (pr%2 == 0 ? ((i_shift+pr)/2)%p : ((i_shift+pr+p)/2)%p);
              }
              for (x=smallsieve_aux1[s]; x<smallsieve_aux1_ub_odd[s]; x+=6) {
                u32_t p, d, pr;

                p = x[0];
                d = x[1];
                pr = x[2];
                x[4] = (pr%2==0 ? ((i_shift+pr)/2)%p : ((i_shift+pr+p)/2)%p);
                x[5] = d / 2;
              }
              for (x=smallpsieve_aux[s]; x<smallpsieve_aux_ub_odd[s]; x+=3)
                x[2] = (x[0]) / 2;

              {
                u16_t *x, *y, *z;

                x = smallsieve_aux1_ub_odd[s];
                y = smallpsieve_aux_ub_odd[s];
                z = smallsieve_aux2[s];
                for (i = 0; i < (u32_t)(4 * x2FBs[s]); i += 4) {
                  u32_t p, pr, d, l;
                  u16_t **a;

                  d = x2FB[s][i + 1];
                  if (d != 1)
                    continue;
                  pr = x2FB[s][i + 2];
                  if (pr % 2 != 0)
                    continue;
                  p = x2FB[s][i];
                  l = x2FB[s][i + 3];
                  if (p < 4) {
                    if (p == 1) {
                      Schlendrian("Use 1=2^0 for sieving?\n");
                    }
                    *y = d;
                    *(y + 1) = l;
                    *(y + 2) = 0;
                    y += 3;
                    continue;
                  }
                  p = p / 2;
                  if (p <= MAX_TINY_2POW)
                    a = &z;
                  else
                    a = &x;
                  **a = p;
                  *(1 + *a) = d;
                  *(2 + *a) = pr % p;
                  *(3 + *a) = l;
                  *(4 + *a) = ((i_shift + pr) / 2) % p;
                  *(5 + *a) = 0;
                  *a += 6;
                }
                smallsieve_aux1_ub[s] = x;
                smallpsieve_aux_ub[s] = y;
                smallsieve_aux2_ub[s] = z;
              }
              break;

            case 3:
              for (x=smallsieve_aux[s]; x<smallsieve_auxbound[s][0]; x+=4) {
                u32_t p, pr;

                p = x[0];
                pr = x[1];
                x[3] = (pr%2 == 1 ? ((i_shift+pr)/2)%p : ((i_shift+pr+p)/2)%p);
              }
              for (x=smallsieve_aux1[s]; x<smallsieve_aux1_ub_odd[s]; x+=6) {
                u32_t p, d, pr;

                p = x[0];
                d = x[1];
                pr = x[2];
                x[4] = (pr%2 == 1 ? ((i_shift+pr)/2)%p : ((i_shift+pr+p)/2)%p);
                x[5] = d/2;
              }
              for (x=smallpsieve_aux[s]; x<smallpsieve_aux_ub_odd[s]; x+=3)
                x[2] = (x[0])/2;

              {
                u16_t *x, *y, *z;

                x = smallsieve_aux1_ub_odd[s];
                y = smallpsieve_aux_ub_odd[s];
                z = smallsieve_aux2[s];
                for (i = 0; i < (u32_t)(4 * x2FBs[s]); i += 4) {
                  u32_t p, pr, d, l;
                  u16_t **a;

                 d = x2FB[s][i + 1];
                  if (d != 1)
                    continue;
                  pr = x2FB[s][i + 2];
                  if (pr % 2 != 1)
                    continue;
                  p = x2FB[s][i];
                  l = x2FB[s][i + 3];
                  if (p < 4) {
                    if (p == 1) {
                      Schlendrian("Use 1=2^0 for sieving?\n");
                    }
                    *y = d;
                    *(y + 1) = l;
                    *(y + 2) = 0;
                    y += 3;
                    continue;
                  }
                  p = p / 2;
                  if (p <= MAX_TINY_2POW)
                    a = &z;
                  else
                    a = &x;
                  **a = p;
                  *(1 + *a) = d;
                  *(2 + *a) = pr % p;
                  *(3 + *a) = l;
                  *(4 + *a) = ((i_shift + pr) / 2) % p;
                  *(5 + *a) = 0;
                  *a += 6;
                }
                smallsieve_aux1_ub[s] = x;
                smallpsieve_aux_ub[s] = y;
                smallsieve_aux2_ub[s] = z;
              }
              break;
          }
        }
#ifdef GCD_SIEVE_BOUND
        for (i = 0; i < np_gcd_sieve; i++) {
          gcd_sieve_buffer[2 * i + 1] =
            (oddness_type / 2) * (gcd_sieve_buffer[2 * i] / 2);
        }
#endif

        j_offset = 0;
        for (s = 0; s < 2; s++) {
          for (i = 0; i < n_schedules[s]; i++) {
            u32_t ns;

            ns = schedules[s][i].n_strips;
            if (ns > n_strips)
              ns = n_strips;
            do_scheduling(schedules[s] + i, ns, oddness_type, s);
            schedules[s][i].current_strip = 0;
          }
        }
#ifdef GATHER_STAT
        Schedule_clock +=
          (1000.0 * (clock() - last_clock)) / CLOCKS_PER_SEC;
#endif

        last_clock = clock();
#ifdef ZSS_STAT
        nss += n_strips;
#endif
        for (subsieve_nr = 0; subsieve_nr < n_strips;
             subsieve_nr++, j_offset += j_per_strip) {
          u16_t s, stepno;

#ifdef USE_MEDSCHED

          for (s = 0; s < 2; s++) {
            u32_t ll, *sched, *ri;

            if (n_medsched_pieces[s] == 0)
              continue;
            for (ll = 0, sched = (u32_t *) med_sched[s][0], ri = LPri[s];
                 ll < n_medsched_pieces[s]; ll++) {
              ri = medsched(ri, current_ij[s] + medsched_fbi_bounds[s][ll],
                            current_ij[s] + medsched_fbi_bounds[s][ll + 1],
                            &sched, medsched_fbi_bounds[s][ll],
                            j_offset == 0 ? oddness_type : 0);
              med_sched[s][ll + 1] = (u16_t *) sched;
            }
          }

          {
            new_clock = clock();
            medsched_clock +=
              (clock_t)((1000.0 * (new_clock - last_clock)) / CLOCKS_PER_SEC);
            last_clock = new_clock;
          }
#endif
          for (s = first_sieve_side, stepno = 0; stepno < 2; stepno++, s = 1 - s) {
            clock_t clock_diff;

            {
              u32_t j;

              for (j = 0; j < j_per_strip; j++) {
                unsigned char *si_ub;

                bzero(tiny_sieve_buffer, TINY_SIEVEBUFFER_SIZE);
                si_ub = tiny_sieve_buffer + TINY_SIEVEBUFFER_SIZE;

                {
                  u16_t *x;

                  for (x = smallsieve_aux[s]; x < smallsieve_tinybound[s]; x += 4) {
                    u32_t p, r, pr;
                    unsigned char l, *si;

					assert(x[2] <= UCHAR_MAX);
                    p = x[0];
                    pr = x[1];
                    l = (unsigned char)x[2];
                    r = x[3];
                    si = tiny_sieve_buffer + r;
                    while (si < si_ub) {
                      *si += l;
                      si += p;
                    }
                    r = r + pr;
                    if (r >= p)
                      r = r - p;
                    x[3] = r;
                  }
                }

                {
                  u16_t *x;

                  for (x = smallsieve_aux2[s]; x < smallsieve_aux2_ub[s]; x += 6) {
                    u32_t p, r, pr, d, d0;
                    unsigned char l, *si;

					assert(x[3] <= UCHAR_MAX);
                    p = x[0];
                    d = x[1];
                    pr = x[2];
                    l = (unsigned char)x[3];
                    r = x[4];
                    d0 = x[5];
                    if (d0 > 0) {
                      x[5]--;
                      continue;
                    }
                    si = tiny_sieve_buffer + r;
                    while (si < si_ub) {
                      *si += l;
                      si += p;
                    }
                    r = r + pr;
                    if (r >= p)
                      r = r - p;
                    x[4] = r;
                    x[5] = d - 1;
                  }
                }

                {
                  u16_t *x;

                  for (x = smallsieve_aux1[s]; x < smallsieve_tinybound1[s]; x += 6) {
                    u32_t p, r, pr, d, d0;
                    unsigned char l, *si;

					assert(x[3] <= UCHAR_MAX);
                    p = x[0];
                    d = x[1];
                    pr = x[2];
                    l = (unsigned char)x[3];
                    r = x[4];

                    d0 = x[5];
                    if (d0 > 0) {
                      x[5]--;
                      continue;
                    }
                    si = tiny_sieve_buffer + r;
                    while (si < si_ub) {
                      *si += l;
                      si += p;
                    }
                    r = r + pr;
                    if (r >= p)
                      r = r - p;
                    x[4] = r;
                    x[5] = d - 1;
                  }
                }

                {
                  unsigned char *si;

                  si = sieve_interval + (j << i_bits);
                  si_ub = sieve_interval + ((j + 1) << i_bits);
                  while (si + TINY_SIEVEBUFFER_SIZE < si_ub) {
                    memcpy(si, tiny_sieve_buffer, TINY_SIEVEBUFFER_SIZE);
                    si += TINY_SIEVEBUFFER_SIZE;
                  }
                  memcpy(si, tiny_sieve_buffer, si_ub - si);
                }

              }
            }

#ifdef ZSS_STAT
            if (s == 1 && ncand == 0)
              nzss[0]++;
#endif
            new_clock = clock();
            clock_diff =
              (clock_t)((1000.0 * (new_clock - last_clock)) / CLOCKS_PER_SEC);
            si_clock[s] += clock_diff;
            sieve_clock += clock_diff;
            last_clock = new_clock;

#ifdef ASM_LINESIEVER
            slinie(smallsieve_tinybound[s], smallsieve_auxbound[s][4],
                   sieve_interval);
#else
            {
              u16_t *x;

              for (x = smallsieve_tinybound[s]; x < smallsieve_auxbound[s][4]; x += 4) {
                u32_t p, r, pr;
                unsigned char l, *y;

				assert(x[2] <= UCHAR_MAX);
                p = x[0];
                pr = x[1];
                l = (unsigned char)x[2];
                r = x[3];
                for (y = sieve_interval; y < sieve_interval + L1_SIZE;
                     y += n_i) {
                  unsigned char *yy, *yy_ub;

                  yy_ub = y + n_i - 3 * p;
                  for (yy = y + r; yy < yy_ub; yy = yy + 4 * p) {
                    *(yy) += l;
                    *(yy + p) += l;
                    *(yy + 2 * p) += l;
                    *(yy + 3 * p) += l;
                  }
                  while (yy < y + n_i) {
                    *(yy) += l;
                    yy += p;
                  }
                  r = r + pr;
                  if (r >= p)
                    r = r - p;
                }
#if 0
                x[3] = r;
#endif
              }
            }
#endif

#if 1
            {
              u16_t *x;

              for (x = smallsieve_auxbound[s][4]; x < smallsieve_auxbound[s][3]; x += 4) {
                u32_t p, r, pr;
                unsigned char l, *y;

				assert(x[2] <= UCHAR_MAX);
                p = x[0];
                pr = x[1];
                l = (unsigned char)x[2];
                r = x[3];
                for (y = sieve_interval; y < sieve_interval + L1_SIZE;
                     y += n_i) {
                  unsigned char *yy;

                  yy = y + r;
                  *(yy) += l;
                  *(yy + p) += l;
                  *(yy + 2 * p) += l;
                  yy += 3 * p;
                  if (yy < y + n_i)
                    *(yy) += l;
                  r = r + pr;
                  if (r >= p)
                    r = r - p;
                }
#if 0
                x[3] = r;
#endif
              }
            }
#endif

#if 1
            {
              u16_t *x;

              for (x = smallsieve_auxbound[s][3]; x < smallsieve_auxbound[s][2]; x += 4) {
                u32_t p, r, pr;
                unsigned char l, *y;

				assert(x[2] <= UCHAR_MAX);
                p = x[0];
                pr = x[1];
                l = (unsigned char)x[2];
                r = x[3];
                for (y = sieve_interval; y < sieve_interval + L1_SIZE;
                     y += n_i) {
                  unsigned char *yy;

                  yy = y + r;
                  *(yy) += l;
                  *(yy + p) += l;
                  yy += 2 * p;
                  if (yy < y + n_i)
                    *(yy) += l;
                  r = r + pr;
                  if (r >= p)
                    r = r - p;
                }
#if 0
                x[3] = r;
#endif
              }
            }
#endif

#if 1
            {
              u16_t *x;

              for (x = smallsieve_auxbound[s][2]; x < smallsieve_auxbound[s][1]; x += 4) {
                u32_t p, r, pr;
                unsigned char l, *y;

				assert(x[2] <= UCHAR_MAX);
                p = x[0];
                pr = x[1];
                l = (unsigned char)x[2];
                r = x[3];
                for (y = sieve_interval; y < sieve_interval + L1_SIZE;
                     y += n_i) {
                  unsigned char *yy;

                  yy = y + r;
                  *(yy) += l;
                  yy += p;
                  if (yy < y + n_i)
                    *(yy) += l;
                  r = r + pr;
                  if (r >= p)
                    r = r - p;
                }
#if 0
                x[3] = r;
#endif
              }
            }
#endif

#if 1
            {
              u16_t *x;

              for (x = smallsieve_auxbound[s][1]; x < smallsieve_auxbound[s][0]; x += 4) {
                u32_t p, r, pr;
                unsigned char l, *y;

				assert(x[2] <= UCHAR_MAX);
                p = x[0];
                pr = x[1];
                l = (unsigned char)x[2];
                r = x[3];
                for (y = sieve_interval; y < sieve_interval + L1_SIZE;
                     y += n_i) {
                  if (r < n_i)
                    *(y + r) += l;
                  r = r + pr;
                  if (r >= p)
                    r = r - p;
                }
#if 0
                x[3] = r;
#endif
              }
            }
#endif

#if 1
            {
              u16_t *x;

              for (x = smallsieve_tinybound1[s]; x < smallsieve_aux1_ub[s]; x += 6) {
                u32_t p, r, pr, d, d0;
                unsigned char l;

				assert(x[3] <= UCHAR_MAX);
                p = x[0];
                d = x[1];
                pr = x[2];
                l = (unsigned char)x[3];
                r = x[4];

                for (d0 = x[5]; d0 < j_per_strip; d0 += d) {
                  unsigned char *y, *yy, *yy_ub;

                  y = sieve_interval + (d0 << i_bits);
                  yy_ub = y + n_i - 3 * p;
                  for (yy = y + r; yy < yy_ub; yy = yy + 4 * p) {
                    *(yy) += l;
                    *(yy + p) += l;
                    *(yy + 2 * p) += l;
                    *(yy + 3 * p) += l;
                  }
                  while (yy < y + n_i) {
                    *(yy) += l;
                    yy += p;
                  }
                  r = r + pr;
                  if (r >= p)
                    r = r - p;
                }
                x[4] = r;
                x[5] = d0 - j_per_strip;
              }
            }
#endif

#if 1
            {
              u16_t *x;

              bzero(horizontal_sievesums, j_per_strip);
              for (x = smallpsieve_aux[s]; x < smallpsieve_aux_ub[s]; x += 3) {
                u32_t p, d;
                unsigned char l;

				assert(x[1] <= UCHAR_MAX);
                p = x[0];
                l = (unsigned char)x[1];
                d = x[2];
                while (d < j_per_strip) {
                  horizontal_sievesums[d] += l;
                  d += p;
                }
#if 0
                  x[2] = d - j_per_strip;
#endif
              }
            }
#else
            bzero(horizontal_sievesums, j_per_strip);
#endif

            new_clock = clock();
            clock_diff =
              (clock_t)((1000.0 * (new_clock - last_clock)) / CLOCKS_PER_SEC);
            s1_clock[s] += clock_diff;
            sieve_clock += clock_diff;
            last_clock = new_clock;
#ifdef GGNFS_BIGENDIAN
#define MEDSCHED_SI_OFFS 1
#else
#define MEDSCHED_SI_OFFS 0
#endif
#ifdef ASM_SCHEDSIEVE1
            schedsieve(medsched_logs[s], n_medsched_pieces[s],
                       med_sched[s], sieve_interval);
#else
            {
              u32_t l;

              for (l = 0; l < n_medsched_pieces[s]; l++) {
                unsigned char x;
                u16_t *schedule_ptr;

                x = medsched_logs[s][l];
                for (schedule_ptr = med_sched[s][l] + MEDSCHED_SI_OFFS;
                     schedule_ptr + 3 * SE_SIZE < med_sched[s][l + 1];
                     schedule_ptr += 4 * SE_SIZE) {
                  sieve_interval[*schedule_ptr] += x;
                  sieve_interval[*(schedule_ptr + SE_SIZE)] += x;
                  sieve_interval[*(schedule_ptr + 2 * SE_SIZE)] += x;
                  sieve_interval[*(schedule_ptr + 3 * SE_SIZE)] += x;
                }
                for (;
                     schedule_ptr < med_sched[s][l + 1];
                     schedule_ptr += SE_SIZE)
                  sieve_interval[*schedule_ptr] += x;
              }
            }
#endif

            new_clock = clock();
            clock_diff =
              (clock_t)((1000.0 * (new_clock - last_clock)) / CLOCKS_PER_SEC);
            s2_clock[s] += clock_diff;
#ifdef BADSCHED
            ncand = 0;
            continue;
#endif
            sieve_clock += clock_diff;
            last_clock = new_clock;
#ifdef GGNFS_BIGENDIAN
#define SCHED_SI_OFFS 1
#else
#define SCHED_SI_OFFS 0
#endif
            {
              u32_t j;

              for (j = 0; j < n_schedules[s]; j++) {
                if (schedules[s][j].current_strip == schedules[s][j].n_strips) {
                  u32_t ns;

                  ns = schedules[s][j].n_strips;
                  if (ns > n_strips - subsieve_nr)
                    ns = n_strips - subsieve_nr;
                  do_scheduling(schedules[s] + j, ns, 0, s);
                  schedules[s][j].current_strip = 0;
                }
              }
#ifdef GATHER_STAT
              new_clock = clock();
              Schedule_clock +=
                (1000.0 * (new_clock - last_clock)) / CLOCKS_PER_SEC;
              last_clock = new_clock;
#endif

              for (j = 0; j < n_schedules[s]; j++) {
#ifdef ASM_SCHEDSIEVE1
                u32_t k;

                k = schedules[s][j].current_strip;
                for (i = 0; i <= schedules[s][j].n_pieces; i++) {
                  schedbuf[i] = schedules[s][j].schedule[i][k];
                }
                schedsieve(schedules[s][j].schedlogs,
                           schedules[s][j].n_pieces, schedbuf,
                           sieve_interval);
#else
                u32_t l, k;

                k = schedules[s][j].current_strip;
                l = 0;
                while (l < schedules[s][j].n_pieces) {
                  unsigned char x;
                  u16_t *schedule_ptr, *sptr_ub;

                  x = schedules[s][j].schedlogs[l];
                  schedule_ptr =
                    schedules[s][j].schedule[l][k] + SCHED_SI_OFFS;
                  while (l < schedules[s][j].n_pieces)
                    if (schedules[s][j].schedlogs[++l] != x)
                      break;
                  sptr_ub = schedules[s][j].schedule[l][k];

#ifdef ASM_SCHEDSIEVE
                  schedsieve(x, sieve_interval, schedule_ptr, sptr_ub);
#else
                  while (schedule_ptr + 3 * SE_SIZE < sptr_ub) {
                    sieve_interval[*schedule_ptr] += x;
                    sieve_interval[*(schedule_ptr + SE_SIZE)] += x;
                    sieve_interval[*(schedule_ptr + 2 * SE_SIZE)] += x;
                    sieve_interval[*(schedule_ptr + 3 * SE_SIZE)] += x;
                    schedule_ptr += 4 * SE_SIZE;
                  }
                  while (schedule_ptr < sptr_ub) {
                    sieve_interval[*schedule_ptr] += x;
                    schedule_ptr += SE_SIZE;
                  }
#endif
                }
#endif
              }
            }

#if 0
            dumpsieve(j_offset, s);
#endif
            new_clock = clock();
            clock_diff =
              (clock_t)((1000.0 * (new_clock - last_clock)) / CLOCKS_PER_SEC);
            sieve_clock += clock_diff;
            s3_clock[s] += clock_diff;
            last_clock = new_clock;

            if (s == first_sieve_side) {
#ifdef GCD_SIEVE_BOUND
              gcd_sieve();
#endif

              {
                u32_t true_i_lb, true_i_ub;

                ncand = 0;
                true_i_ub = j_offset + j_per_strip;

                {
                  i32_t lt;

                  {
                    double a;

                    a = g_poldeg[s] * log(j_offset + 1);
                    if (special_q_side == s)
                      a -= special_q_log;
					
                    a = rint(a*sieve_multiplier[s] - sieve_report_multiplier[s]*FB_maxlog[s]);
					assert( (a >= INT_MIN) && (a <= INT_MAX) );
                    lt = (i32_t)a;
                  }
                  for (true_i_lb=0; true_i_lb<true_i_ub; true_i_lb+=CANDIDATE_SEARCH_STEPS) 
				  {
                    i32_t sieve_threshold;
                    u32_t lulb, luub;
                    size_t j;

                    lulb = (true_i_lb<<(J_bits-1))/(j_offset+j_per_strip);
                    lulb = (n_J+lulb)/(2*CANDIDATE_SEARCH_STEPS);
                    i = j_offset == 0 ? 1 : j_offset;
                    luub = ((true_i_lb+CANDIDATE_SEARCH_STEPS+i-1)<<(J_bits-1))/i;
                    luub = (n_J+luub+2*CANDIDATE_SEARCH_STEPS-1)/(2*CANDIDATE_SEARCH_STEPS);
                    if (luub > n_J / CANDIDATE_SEARCH_STEPS)
                      luub = n_J / CANDIDATE_SEARCH_STEPS;
                    for (i = lulb, sieve_threshold = INT_MAX; i < luub; i++)
                      if (plog_lb1[i] < sieve_threshold)
                        sieve_threshold = plog_lb1[i];
                    sieve_threshold += lt;
                    for (j = 0; j < j_per_strip; j++) {
                      unsigned char *i_o, *i_max, *i_maxx, st1;

                      i_o = sieve_interval + (j << i_bits) + i_shift / 2;
                      i_maxx = i_o + j + j_offset;
                      i_o += true_i_lb;
                      i_max = i_o + CANDIDATE_SEARCH_STEPS;
                      if (sieve_threshold <= (i32_t) horizontal_sievesums[j])
                        st1 = 1;
                      else
                        st1 = sieve_threshold - horizontal_sievesums[j];
                      optsieve(st1,i_o,MIN(i_max,i_maxx),j);
                    }

                    i = luub;
                    luub = n_J / CANDIDATE_SEARCH_STEPS - lulb;
                    lulb = n_J / CANDIDATE_SEARCH_STEPS - i;
                    for (i = lulb, sieve_threshold = INT_MAX; i < luub; i++)
                      if (plog_lb1[i] < sieve_threshold)
                        sieve_threshold = plog_lb1[i];
                    sieve_threshold += lt;
                    for (j = 0; j < j_per_strip; j++) {
                      unsigned char *i_o, *i_min, *i_max, st1;

                      i_o = sieve_interval + (j << i_bits) + i_shift / 2;
                      i_min = i_o - (j + j_offset);
                      i_max = i_o - true_i_lb;
                      i_o = i_max - CANDIDATE_SEARCH_STEPS;
                      if (sieve_threshold <= (i32_t) horizontal_sievesums[j])
                        st1 = 1;
                      else
                        st1 = sieve_threshold - horizontal_sievesums[j];
                      optsieve(st1,MAX(i_o,i_min),i_max,j);
                    }
                  }
                }

                {
                  u32_t next_lt_update;
                  i32_t lt=0;

                  for (true_i_lb=j_offset&(~(CANDIDATE_SEARCH_STEPS - 1)),
                       next_lt_update=0; true_i_lb<n_i/2; true_i_lb += CANDIDATE_SEARCH_STEPS) {
                    u32_t lulb, luub;
                    u32_t j;
                    i32_t sieve_threshold;

                    if (true_i_lb >= next_lt_update) {
                      double a;

                      a = g_poldeg[s] * log(true_i_lb + 1);
                      if (special_q_side == first_sieve_side)
                        a -= log(special_q);

					  a = rint(a*sieve_multiplier[s]-sieve_report_multiplier[s]*FB_maxlog[s]);
					  assert( (a >= INT_MIN) && (a <= INT_MAX) );
                      lt = (i32_t)a;
                      next_lt_update *= 2;
                      next_lt_update++;
                    }
                    lulb = (j_offset<<(J_bits-1))/(true_i_lb+CANDIDATE_SEARCH_STEPS);
                    lulb = (n_J + lulb) / (2 * CANDIDATE_SEARCH_STEPS);
                    i = true_i_lb == 0 ? 1 : true_i_lb;
                    luub = ((j_per_strip+j_offset+i-1)<<(J_bits-1))/i;
                    luub = (n_J+luub+2*CANDIDATE_SEARCH_STEPS-1)/(2*CANDIDATE_SEARCH_STEPS);
                    if (luub > n_J / CANDIDATE_SEARCH_STEPS)
                      luub = n_J / CANDIDATE_SEARCH_STEPS;
                    for (i = lulb, sieve_threshold = INT_MAX; i < luub; i++)
                      if (plog_lb2[i] < sieve_threshold)
                        sieve_threshold = plog_lb2[i];
                    sieve_threshold += lt;
                    for (j = 0; j < j_per_strip; j++) {
                      unsigned char *i_o, *i_min, *i_max, st1;

                      i_o = sieve_interval + (j << i_bits) + i_shift / 2;

                      i_min = i_o + j + j_offset;
                      i_o += true_i_lb;
                      i_max = i_o + CANDIDATE_SEARCH_STEPS;
                      if (sieve_threshold <= (i32_t) horizontal_sievesums[j])
                        st1 = 1;
                      else
                        st1 = sieve_threshold - horizontal_sievesums[j];
                      optsieve(st1,MAX(i_o,i_min),i_max,j);
                    }

                    i = luub;
                    luub = n_J / CANDIDATE_SEARCH_STEPS - lulb;
                    lulb = n_J / CANDIDATE_SEARCH_STEPS - i;
                    if (luub > n_J / CANDIDATE_SEARCH_STEPS)
                      luub = n_J / CANDIDATE_SEARCH_STEPS;
                    for (i = lulb, sieve_threshold = INT_MAX; i < luub; i++)
                      if (plog_lb2[i] < sieve_threshold)
                        sieve_threshold = plog_lb2[i];
                    sieve_threshold += lt;
                    for (j = 0; j < j_per_strip; j++) {
                      unsigned char *i_o, *i_max, *i_maxx, st1;

                      i_o = sieve_interval + (j << i_bits) + i_shift / 2;
                      i_maxx = i_o - (j + j_offset);
                      i_max = i_o - true_i_lb;
                      i_o = i_max - CANDIDATE_SEARCH_STEPS;
                      if (sieve_threshold <= (i32_t) horizontal_sievesums[j])
                        st1 = 1;
                      else
                        st1 = sieve_threshold - horizontal_sievesums[j];
                      optsieve(st1,i_o,MIN(i_max,i_maxx),j);
                    }
                  }
                }

              }

            } else {
              u32_t nc1;

              n_prereports += ncand;
              for (i = 0, nc1 = 0; i < ncand; i++) {
                u16_t st_i, t_j, ii, jj, j;
                double pvl;

                j = cand[i] >> i_bits;
                jj = j_offset + j;
                ii = cand[i] & (n_i - 1);
                st_i = 2 * ii + (oddness_type == 2 ? 0 : 1);
                t_j = 2 * jj + (oddness_type == 1 ? 0 : 1);
                pvl = log(fabs(rpol_eval(tpoly_f[s], g_poldeg[s], 
                              (double) st_i - (double) i_shift, (double) t_j)));
                if (special_q_side == s)
                  pvl -= special_q_log;
                pvl *= sieve_multiplier[s];
                pvl -= sieve_report_multiplier[s] * FB_maxlog[s];
                if ((double)(sieve_interval[cand[i]]+horizontal_sievesums[j]) >= pvl) {
                  fss_sv[nc1] = fss_sv[i];
                  cand[nc1++] = cand[i];
                }
              }
              ncand = nc1;
            }

            new_clock = clock();
            clock_diff = (clock_t)((1000.0 * (new_clock - last_clock)) / CLOCKS_PER_SEC);
            sieve_clock += clock_diff;
            cs_clock[s] += clock_diff;
            last_clock = new_clock;
          }

#ifndef NO_TDCODE
          {
            u32_t ci;
            u32_t nc1;
            u16_t side, tdstep;
            clock_t last_tdclock, newclock;

            for (ci = 0, nc1 = 0; ci < ncand; ci++) {
              u16_t strip_i, strip_j;
              u16_t st_i, true_j;
              u16_t s;
              double pvl;

              {
                u16_t jj;

                strip_j = cand[ci] >> i_bits;
                jj = j_offset + strip_j;
                strip_i = cand[ci] & (n_i - 1);
                st_i = 2 * strip_i + (oddness_type == 2 ? 0 : 1);
                true_j = 2 * jj + (oddness_type == 1 ? 0 : 1);
              }

              n_reports++;
              s = first_sieve_side;
#ifdef STC_DEBUG
              fprintf(debugfile, "%hu %hu\n", st_i, true_j);
#endif
              if (gcd32(st_i < i_shift ? i_shift - st_i : st_i - i_shift,true_j) != 1)
                continue;
              n_rep1++;
              pvl = log(fabs(rpol_eval(tpoly_f[s], g_poldeg[s],
                         (double) st_i-(double) i_shift, (double) true_j)));
              if (special_q_side == s)
                pvl -= special_q_log;
              pvl *= sieve_multiplier[s];
              if ((double) fss_sv[ci] + sieve_report_multiplier[s] * FB_maxlog[s] < pvl)
                continue;
              n_rep2++;

              modulo32 = special_q;
              if (modadd32(modmul32(st_i, spq_i), modmul32(true_j, spq_j)) == spq_x)
                continue;
              cand[nc1++] = cand[ci];
            }

#ifdef ZSS_STAT
            if (ncand == 0)
              nzss[1]++;
#endif
            last_tdclock = clock();
            tdi_clock += (clock_t)((1000.0 * (last_tdclock - last_clock)) / CLOCKS_PER_SEC);
            ncand = nc1;
            qsort(cand, ncand, sizeof(*cand), tdcand_cmp);
            td_buf1[0] = td_buf[first_td_side];
            for (side=first_td_side, tdstep=0; tdstep<2; side=1-side, tdstep++) {
#ifdef ZSS_STAT
              if (tdstep == 1 && ncand == 0)
                nzss[2]++;
#endif

              {
                size_t nfbp;
                u32_t p_bound;
                u16_t last_j, strip_i, strip_j;

                nfbp = 0;

                if (ncand > 0)
                  p_bound = (2 * n_i * j_per_strip) / (5 * ncand);
                else
                  p_bound = UINT_MAX;

                {
                  unsigned char ht, allcoll;

                  bzero(sieve_interval, L1_SIZE);
                  bzero(tds_coll, UCHAR_MAX - 1);
                  for (ci = 0, ht = 1, allcoll = 0; ci < ncand; ci++) {
                    unsigned char cht;

                    cht = sieve_interval[cand[ci]];
                    if (cht == 0) {
                      cht = ht;
                      if (ht < UCHAR_MAX)
                        ht++;
                      else {
                        ht = 1;
                        allcoll = 1;
                      }
                      tds_coll[cht - 1] = allcoll;
                      sieve_interval[cand[ci]] = cht;
                    } else {
                      tds_coll[cht - 1] = 1;
                    }
                    fss_sv[ci] = cht - 1;
                  }
                }

                {
                  u16_t *x, *y;

                  y = smalltdsieve_aux[side][j_per_strip - 1];
                  for (i=0, x=smallsieve_aux[side]; x<smallsieve_tinybound[side]; i++, x += 4) {
                    modulo32 = x[0];
                    x[3] = modsub32((u32_t) x[3], (u32_t) y[i]);
                  }
                }

                newclock = clock();
                tdsi_clock[side] += (u32_t)((1000.0 * (newclock - last_tdclock)) / CLOCKS_PER_SEC);
                last_tdclock = newclock;

                memcpy(tds_fbi_curpos, tds_fbi, UCHAR_MAX * sizeof(*tds_fbi));

#ifdef ASM_SCHEDTDSIEVE
                {
                  u32_t x, *(y[2]);

                  x = 0;
                  y[0] = med_sched[side][0];
                  y[1] = med_sched[side][n_medsched_pieces[side]];
                  schedtdsieve(&x, 1, y, sieve_interval, tds_fbi_curpos);
                }
#else
                {
                  u32_t l;

                  for (l = 0; l < n_medsched_pieces[side]; l++) {
                    u16_t *x, *x_ub;

                    x_ub = med_sched[side][l + 1];

                    for (x=med_sched[side][l]+MEDSCHED_SI_OFFS; x+6<x_ub; x += 8) {
                      unsigned char z;

                      if ((sieve_interval[*x] | sieve_interval[*(x + 2)] |
                           sieve_interval[*(x + 4)] | sieve_interval[*(x + 6)]) == 0) {
                        continue;
                      }
                      if ((z = sieve_interval[*x]) != 0)
                        *(tds_fbi_curpos[z - 1]++) = *(x + 1 - 2 * MEDSCHED_SI_OFFS);
                      if ((z = sieve_interval[*(x + 2)]) != 0)
                        *(tds_fbi_curpos[z - 1]++) = *(x + 3 - 2 * MEDSCHED_SI_OFFS);
                      if ((z = sieve_interval[*(x + 4)]) != 0)
                        *(tds_fbi_curpos[z - 1]++) = *(x + 5 - 2 * MEDSCHED_SI_OFFS);
                      if ((z = sieve_interval[*(x + 6)]) != 0)
                        *(tds_fbi_curpos[z - 1]++) = *(x + 7 - 2 * MEDSCHED_SI_OFFS);
                    }
                    while (x < x_ub) {
                      unsigned char z;

                      if ((z = sieve_interval[*x]) != 0)
                        *(tds_fbi_curpos[z - 1]++) = *(x + 1 - 2 * MEDSCHED_SI_OFFS);
                      x += 2;
                    }
                  }
                }
#endif
                newclock = clock();
                tds2_clock[side] += (u32_t)((1000.0 * (newclock - last_tdclock)) / CLOCKS_PER_SEC);
                last_tdclock = newclock;

                {
                  u32_t j;

                  for (j = 0; j < n_schedules[side]; j++) {
#ifdef ASM_SCHEDTDSIEVE
                    u32_t k;

                    k = schedules[side][j].current_strip++;
                    for (i = 0; i <= schedules[side][j].n_pieces; i++) {
                      schedbuf[i] = schedules[side][j].schedule[i][k];
                    }
                    schedtdsieve(schedules[side][j].fbi_bounds, schedules[side][j].n_pieces, 
                                 schedbuf, sieve_interval, tds_fbi_curpos);
#else
                    u32_t l, k;

                    k = schedules[side][j].current_strip++;
                    for (l = 0; l < schedules[side][j].n_pieces; l++) {
                      u16_t *x, *x_ub;
                      u32_t fbi_offset;

                      x_ub = schedules[side][j].schedule[l + 1][k];
                      fbi_offset = schedules[side][j].fbi_bounds[l];
                      for (x=schedules[side][j].schedule[l][k]+SCHED_SI_OFFS; x+6<x_ub; x+=8) {
                        unsigned char z;

                        if ((sieve_interval[*x] | sieve_interval[*(x + 2)]
                             | sieve_interval[*(x + 4)] | sieve_interval[*(x + 6)]) == 0) {
                          continue;
                        }
                        if ((z = sieve_interval[*x]) != 0)
                          *(tds_fbi_curpos[z-1]++) = fbi_offset + *(x+1-2*SCHED_SI_OFFS);
                        if ((z = sieve_interval[*(x + 2)]) != 0)
                          *(tds_fbi_curpos[z-1]++) = fbi_offset + *(x+3-2*SCHED_SI_OFFS);
                        if ((z = sieve_interval[*(x + 4)]) != 0)
                          *(tds_fbi_curpos[z-1]++) = fbi_offset + *(x+5-2*SCHED_SI_OFFS);
                        if ((z = sieve_interval[*(x + 6)]) != 0)
                          *(tds_fbi_curpos[z-1]++) = fbi_offset + *(x+7-2*SCHED_SI_OFFS);
                      }
                      while (x < x_ub) {
                        unsigned char z;

                        if ((z = sieve_interval[*x]) != 0)
                          *(tds_fbi_curpos[z-1]++) = fbi_offset + *(x+1-2*SCHED_SI_OFFS);
                        x += 2;
                      }
                    }
#endif
                  }
                }
                newclock = clock();
                tds3_clock[side] += (u32_t)((1000.0 * (newclock - last_tdclock)) / CLOCKS_PER_SEC);
                last_tdclock = newclock;

                for (i = 0; i < UCHAR_MAX && i < ncand; i++) {
                  u32_t *p;

                  for (p = tds_fbi[i]; p < tds_fbi_curpos[i]; p++)
                    *p = FB[side][*p];
                }

                {
                  u16_t *x;

                  for (x=smallsieve_auxbound[side][0]-4; x>=smallsieve_aux[side] && *x > p_bound;
                       x = x - 4) {
                    u32_t p, r, pr;
                    unsigned char *y;

                    p = x[0];
                    pr = x[1];
                    r = x[3];
                    modulo32 = p;
                    for (y=sieve_interval; y<sieve_interval+L1_SIZE; y+=n_i) {
                      unsigned char *yy, *yy_ub;

                      yy_ub = y + n_i;
                      yy = y + r;
                      while (yy < yy_ub) {
                        if (*yy != 0)
                          *(tds_fbi_curpos[*yy - 1]++) = p;
                        yy += p;
                      }
                      r = modadd32(r, pr);
                    }
                    x[3] = r;
                  }
                  newclock = clock();
                  tds1_clock[side] += (u32_t)((1000.0 * (newclock - last_tdclock)) / CLOCKS_PER_SEC);
                  last_tdclock = newclock;
                }

                last_j = 0;
                for (ci = 0, nc1 = 0; ci < ncand; ci++) {
                  u32_t *fbp_buf;
                  u32_t *fbp_ptr;
                  u16_t st_i, true_j;
                  i32_t true_i;

                  u32_t coll;

                  {
                    u16_t jj;

                    strip_j = cand[ci] >> i_bits;
                    jj = j_offset + strip_j;
                    strip_i = cand[ci] & (n_i - 1);
                    st_i = 2 * strip_i + (oddness_type == 2 ? 0 : 1);
                    true_j = 2 * jj + (oddness_type == 1 ? 0 : 1);
                  }

                  if (strip_j != last_j) {
                    u16_t j_step;

                    if (strip_j <= last_j)
                      Schlendrian("TD: Not sorted\n");
                    j_step = strip_j - last_j;
                    last_j = strip_j;

                    {
                      u16_t *x, *y;

                      y = smalltdsieve_aux[side][j_step - 1];
                      for (i=0, x=smallsieve_aux[side]; x<smallsieve_auxbound[side][0]; i++,x += 4) {
                        modulo32 = x[0];
                        if (modulo32 > p_bound)
                          break;
                        x[3] = modadd32((u32_t) x[3], (u32_t) y[i]);
                      }
                      for (x=smallpsieve_aux[side]; x<smallpsieve_aux_ub[side]; x += 3) {
                        modulo32 = x[0];
                        x[2] = modsub32(x[2], (j_step) % modulo32);
                      }
                    }

                  }
                  true_i = (i32_t) st_i - (i32_t) i_shift;

                  mpz_set_si(g_aux1, true_i);
                  mpz_mul_si(g_aux1, g_aux1, a0);
                  mpz_set_si(g_aux2, a1);
                  mpz_mul_ui(g_aux2, g_aux2, (u32_t) true_j);
                  mpz_add(g_sr_a, g_aux1, g_aux2);

                  mpz_set_si(g_aux1, true_i);
                  mpz_mul_si(g_aux1, g_aux1, b0);
                  mpz_set_si(g_aux2, b1);
                  mpz_mul_ui(g_aux2, g_aux2, (u32_t) true_j);
                  mpz_add(g_sr_b, g_aux1, g_aux2);

                  if (mpz_sgn(g_sr_b) < 0) {
                    mpz_neg(g_sr_b, g_sr_b);
                    mpz_neg(g_sr_a, g_sr_a);
                  }

                  {

                    i = 1;
                    mpz_set(g_aux2, g_sr_a);
                    mpz_set(g_aux1, g_poly[side][0]);
                    for (;;) {
                      mpz_mul(g_aux1, g_aux1, g_sr_b);
                      mpz_mul(g_aux3, g_aux2, g_poly[side][i]);
                      mpz_add(g_aux1, g_aux1, g_aux3);
                      if (++i > g_poldeg[side])
                        break;
                      mpz_mul(g_aux2, g_aux2, g_sr_a);
                    }
                  }

                  if (td_buf_alloc[side] < nfbp + mpz_sizeinbase(g_aux1, 2)) {
                    td_buf_alloc[side] += 1024;
                    while (td_buf_alloc[side] < nfbp + mpz_sizeinbase(g_aux1, 2)) {
                      td_buf_alloc[side] += 1024;
                    }
                    td_buf[side] = xrealloc(td_buf[side], td_buf_alloc[side] * sizeof(**td_buf));
                    if (side == first_td_side) {
                      u32_t *oldptr;

                      oldptr = td_buf1[0];
                      for (i = 0; i <= nc1; i++)
                        td_buf1[i] = td_buf[side] + (td_buf1[i] - oldptr);
                    }
                  }
                  if (side == first_td_side)
                    fbp_buf = td_buf1[nc1];
                  else
                    fbp_buf = td_buf[side];
                  fbp_ptr = fbp_buf;

                  {
                    u32_t *x, *y;

                    coll = tds_coll[fss_sv[ci]];
                    x = tds_fbi[fss_sv[ci]];
                    y = tds_fbi_curpos[fss_sv[ci]];

                    while (x < y) {
                      u32_t p;

                      p = *(x++);
                      if (coll != 0) {
                        while (mpz_fdiv_q_ui(g_aux2, g_aux1, p) == 0) {
                          mpz_set(g_aux1, g_aux2);
                          *(fbp_ptr++) = p;
                        }
                        continue;
                      }

                      if (mpz_fdiv_q_ui(g_aux1, g_aux1, p) != 0)
                        Schlendrian("TD sieve on side %u: %u does not divide\n",side, p);
                      *(fbp_ptr++) = p;
                      for (;;) {
                        if (mpz_fdiv_q_ui(g_aux2, g_aux1, p) != 0)
                          break;
                        mpz_set(g_aux1, g_aux2);
                        *(fbp_ptr++) = p;
                      }
                    }
                  }

                  {
                    u16_t *x;

#ifdef PREINVERT

                    {
                      u32_t *p_inv;

                      p_inv = smalltd_pi[side];
                      for (x=smallsieve_aux[side]; x<smallsieve_auxbound[side][0]
                               && *x <= p_bound; x += 4, p_inv++) {
                        modulo32 = *x;

                        if (((modsub32((u32_t)strip_i,(u32_t)(x[3]))*(*p_inv))&0xffff0000)==0) {
                          u32_t p;

                          p = *x;
                          if (mpz_fdiv_q_ui(g_aux1, g_aux1, p) != 0)
                            Schlendrian("TD small FB on side %u: %u does not divide\n",side, p);
                          *(fbp_ptr++) = p;
                          for (;;) {
                            if (mpz_fdiv_q_ui(g_aux2, g_aux1, p) != 0)
                              break;
                            mpz_set(g_aux1, g_aux2);
                            *(fbp_ptr++) = p;
                          }
                        }
                      }
                    }

#else
                    for (x=smallsieve_aux[side]; x<smallsieve_auxbound[side][0]
                            && *x <= p_bound; x += 4) {
                      u32_t p;

                      p = *x;
                      if (strip_i % p == x[3]) {
                        if (mpz_fdiv_q_ui(aux1, aux1, p) != 0)
                          Schlendrian("TD small FB on side %u: %u does not divide\n",side, p);
                        *(fbp_ptr++) = p;
                        for (;;) {
                          if (mpz_fdiv_q_ui(aux2, aux1, p) != 0)
                            break;
                          mpz_set(aux1, aux2);
                          *(fbp_ptr++) = p;
                        }
                      }
                    }
#endif
                    for (x=smallpsieve_aux[side]; x<smallpsieve_aux_ub_pow1[side]; x+=3) {
                      if (x[2] == 0) {
                        u32_t p;

                        p = *x;
                        if (coll && g_poldeg[side] > 0 && p > p_bound) {
                          while (mpz_fdiv_q_ui(g_aux2, g_aux1, p) == 0) {
                            mpz_set(g_aux1, g_aux2);
                            *(fbp_ptr++) = p;
                          }
                          continue;
                        }
                        if (mpz_fdiv_q_ui(g_aux1, g_aux1, p) != 0)
                          Schlendrian("TD projective FB on side %u: %u does not divide\n",
                                        side, p);
                        *(fbp_ptr++) = p;
                        for (;;) {
                          if (mpz_fdiv_q_ui(g_aux2, g_aux1, p) != 0)
                            break;
                          mpz_set(g_aux1, g_aux2);
                          *(fbp_ptr++) = p;
                        }
                      }
                    }
                  }

                  while (mpz_get_ui(g_aux1) % 2 == 0) {
                    mpz_fdiv_q_2exp(g_aux1, g_aux1, 1);
                    *(fbp_ptr++) = 2;
                  }

                  if (side == special_q_side) {
                    if (mpz_fdiv_q_ui(g_aux1, g_aux1, special_q) != 0) {
                      Schlendrian("Special q %u does not divide\n", special_q);
                    }
                    *(fbp_ptr++) = special_q;
                  }

                  if (mpz_sizeinbase(g_aux1, 2) <= max_factorbits[side]) {
                    n_tdsurvivors[side]++;
                    if (side == first_td_side) {
                      if (mpz_sgn(g_aux1) > 0)
                        mpz_set(td_rests[nc1], g_aux1);
                      else
                        mpz_neg(td_rests[nc1], g_aux1);
                      cand[nc1++] = cand[ci];
                      td_buf1[nc1] = fbp_ptr;
                      nfbp = fbp_ptr - td_buf[side];
                      continue;
                    }
#define KLEINJUNG_MPQS
                    {
                      u32_t s, *(fbp_buffers[2]), *(fbp_buffers_ub[2]);
                      i16_t need_mpqs[2];
                      size_t nlp[2];
                      clock_t cl;
                      u32_t ov;

                      s = first_td_side;
                      fbp_buffers[s] = td_buf1[ci];
                      fbp_buffers_ub[s] = td_buf1[ci + 1];
                      fbp_buffers[1 - s] = fbp_buf;
                      fbp_buffers_ub[1 - s] = fbp_ptr;
                      mpz_set(large_factors[s], td_rests[ci]);
                      if (mpz_sgn(g_aux1) > 0)
                        mpz_set(large_factors[1 - s], g_aux1);
                      else
                        mpz_neg(large_factors[1 - s], g_aux1);
#if 0
                      printf("%u %u\n",
                             mpz_sizeinbase(large_factors[0], 2),
                             mpz_sizeinbase(large_factors[1], 2));
#endif
                      for (s = 0; s < 2; s++) {
                        i16_t is_prime;
                        u16_t s1;

                        s1 = s ^ first_psp_side;
                        if (mpz_cmp_ui(large_factors[s1], 1) == 0) {
                          nlp[s1] = 0;
                          need_mpqs[s1] = 0;
                        } else {
                          if (mpz_cmp(large_factors[s1], FBb_sq[s1]) < 0)
                            is_prime = 1;
                          else
//                            is_prime = psp(large_factors[s1], 1);
                            is_prime = psp(large_factors[s1]);
                          if (is_prime == 1) {
                            if (mpz_sizeinbase(large_factors[s1], 2) >
                                max_primebits[s1])
                              break;
                            else {
                              mpz_set(large_primes[s1][0],
                                      large_factors[s1]);
                              need_mpqs[s1] = 0;
                              nlp[s1] = 1;
                            }
                          } else {
                            need_mpqs[s1] = 1;
                          }
                        }
                      }
                      if (s != 2)
                        continue;

                      cl = clock();
                      ov = verbose;
                      verbose = 0;
                      for (s = 0; s < 2; s++) {
                        u16_t s1;

                        s1 = s ^ first_mpqs_side;
                        if (need_mpqs[s1]) {
#if 0
#define KLEINJUNG_MPQS
#endif
						  long nf;
                          mpz_t *mf;
#ifdef KLEINJUNG_MPQS

                          if ((nf = mpqs_factor(large_factors[s1], max_primebits[s1], &mf)) < 0) {
#if 0
                            n_mpqsfail[s1]++;
                            break;
#else
                            goto attempt_mpqs4linux;
#endif
                          }
                          if (nf == 0) {
                            n_mpqsvain[s1]++;
                            break;
                          }
                          for (i = 0; i < nf; i++)
			    mpz_set(large_primes[s1][i], mf[i]);
                          nlp[s1] = nf;
                          continue;
#endif
                        attempt_mpqs4linux:
#if 0
                          nlp[s1] =
                            mpqs(large_primes[s1], large_factors[s1],
                                 NULL);
#endif
#define TRY_RHO_ON_FAILURES
#ifdef TRY_RHO_ON_FAILURES
                        { unsigned long small_factors[10];

                          nf = rho_factor(small_factors, large_factors[s1]);
                          for (i = 0; i < nf; i++) {
                            mpz_set_ui(large_primes[s1][i], small_factors[i]);
                            if (mpz_sizeinbase(large_primes[s1][i],2) > max_primebits[s1]) {
                              n_mpqsvain[s1]++;
                              break;
                            }
                          }
                          if ((i >= nf) && (nf >= 1))
                            nlp[s1] = nf;
                          else { nlp[s1]=0; }
                        }
#else
                          nlp[s1] = 0;
#endif
                          if (nlp[s1] == 0) {
                            if (ov > 1) {
                              fprintf(stderr, "mpqs failed for ");
                              mpz_out_str(stderr, 10, large_factors[s1]);
                              fprintf(stderr, "(a,b): ");
                              mpz_out_str(stderr, 10, g_sr_a);
                              fprintf(stderr, " ");
                              mpz_out_str(stderr, 10, g_sr_b);
                              fprintf(stderr, "\n");
                            }
                            n_mpqsfail[s1]++;
                            break;
                          }
#ifdef KLEINJUNG_MPQS
                          if (ov > 1) {
                            fprintf(stderr, "mpqs factored ");
                            mpz_out_str(stderr, 10, large_factors[s1]);
                            fprintf(stderr, "\n");
                          }
#endif
                          for (i = 0; i < nlp[s1]; i++) {
                            if (mpz_sizeinbase(large_primes[s1][i], 2) >
                                max_primebits[s1]) {
                              n_mpqsvain[s1]++;
printf("Too large!\n");
                              break;
                            }
                          }
                          if (i < nlp[s1])
                            break;
                        }
                      }
                      verbose = ov;
                      mpqs_clock += (clock_t)((1000.0 * (clock() - cl)) / CLOCKS_PER_SEC);
                      if (s != 2)
                        continue;
//                        fprintf(ofile, "W ");
#define OBASE 16
                      yield++;
                      mpz_out_str(g_ofile, 10, g_sr_a);
                      fprintf(g_ofile, ",");
                      mpz_out_str(g_ofile, 10, g_sr_b);

                      if (short_output == 0) /* Sten: added -s parameter to the command line. */
                      {
#ifdef _ORIG_OUTPUT_FORMAT
                         for (s = 0; s < 2; s++) {
                           u32_t *x;
                          
                           fprintf(ofile, "\n%c", 'X' + s);
                           for (i = 0; i < nlp[s]; i++) {
                             fprintf(ofile, " ");
                             mpz_out_str(ofile, OBASE, large_primes[s][i]);
                           }
                           for (x = fbp_buffers[s]; x < fbp_buffers_ub[s]; x++) {
                             fprintf(ofile, " %X", *x);
                           }
                         }
#else
                         { int numR=0;
                          u32_t *x;
    
                          fprintf(g_ofile, ":");
                          for (i = 0; i < nlp[1]; i++) { /* rational first. */
                            if (i>0) fprintf(g_ofile, ",");
                            mpz_out_str(g_ofile, OBASE, large_primes[1][i]);
                            numR++;
                          }
                          for (x = fbp_buffers[1]; x < fbp_buffers_ub[1];x++) {
                            if (numR>0) fprintf(g_ofile, ",%X", (unsigned int)*x);
                            else { fprintf(g_ofile, "%X", (unsigned int)*x); numR++;}
                          }
                        }
                        { int numA=0;
                          u32_t *x;
    
                          fprintf(g_ofile, ":");
                          for (i = 0; i < nlp[0]; i++) { /* algebraic next. */
                            if (i>0) fprintf(g_ofile, ",");
                            mpz_out_str(g_ofile, OBASE, large_primes[0][i]);
                            numA++;
                          }
                          for (x = fbp_buffers[0]; x < fbp_buffers_ub[0];x++) {
                            if (numA>0) fprintf(g_ofile, ",%X", (unsigned int)*x);
                            else { fprintf(g_ofile, "%X", (unsigned int)*x); numA++;}
                          }
                        }
#endif
                      } /* if (short_output == 0) */

                      fprintf(g_ofile, "\n");
                    }
                  } else
                    continue;
                }
                {
                  u16_t j_step;

                  j_step = j_per_strip - last_j;
                  {
                    u16_t *x, *y;

                    y = smalltdsieve_aux[side][j_step - 1];
                    for (i = 0, x = smallsieve_aux[side];
                         x < smallsieve_auxbound[side][0]; i++, x += 4) {
                      modulo32 = x[0];
                      if (modulo32 > p_bound)
                        break;
                      x[3] = modadd32((u32_t) x[3], (u32_t) y[i]);
                    }
                    for (x=smallpsieve_aux[side]; x<smallpsieve_aux_ub[side]; x+=3) {
                      modulo32 = x[0];
                      x[2] = modsub32(x[2], (j_step) % modulo32);
                    }
                  }
                }
                newclock = clock();
                tds4_clock[side] += (u32_t)((1000.0*(newclock-last_tdclock))/CLOCKS_PER_SEC);
                last_tdclock = newclock;
                ncand = nc1;
              }
            }
          }
          new_clock = clock();
          td_clock += (clock_t)((1000.0 * (new_clock - last_clock)) / CLOCKS_PER_SEC);
          last_clock = new_clock;
        }
      }
      last_clock = new_clock;
    }
    if (root_no < nr) {
      break;
    }
    tNow = sTime();
    if (tNow > lastReport + 5.0) {
      lastReport = sTime();
      fprintf(stderr, "\rtotal yield: %u, q=%u (%1.5lf sec/rel)", 
    	    (unsigned int)yield, (unsigned int)special_q, (tNow - tStart)/yield);
      fflush(stderr);
      { char *ofn;
        FILE *of;

        asprintf(&ofn, ".last_spq%d", process_no);
        if ((of = fopen(ofn, "wb")) != 0) {
          fprintf(of, "%u\n", (unsigned int)special_q);
          fclose(of);
        }
        free(ofn);
      }
    }
  }
  fprintf(stderr, "\rtotal yield: %u, q=%u (%1.5lf sec/rel)\n", 
	(unsigned int)yield, (unsigned int)special_q, (sTime() - tStart)/yield);
  free(r_ptr);
  return 0;
}


/**************************************************/
int parseJobFile(char *fName)
/**************************************************/
{ FILE *fp;
  char token[256], value[512], thisLine[1024];


  sieve_min[0] = sieve_min[1]=0;

  if (!(fp = fopen(fName, "rb"))) {
    printf("Error opening %s for read!\n", fName);
    return -1;
  }
  input_poly(g_N, g_poly, g_poldeg, g_poly + 1, g_poldeg + 1, g_m, fp);
  rewind(fp);
  while (!(feof(fp))) {
    thisLine[0] = 0;
    fgets(thisLine, 1023, fp);
    /* Special case: If there's a polynomial, handle it seperately: */
    if (strncmp(thisLine, "START_POLY", 10)==0) {
      while (!(feof(fp)) && strncmp(thisLine, "END_POLY", 8)) 
        fgets(thisLine, 1023, fp);
    } else  if ((sscanf(thisLine, "%255s %511s", token, value)==2) && 
                (thisLine[0] != '#')) {

	  token[sizeof(token)-1] = 0;
      if (strncmp(token, "skew:", 5)==0) {
        sigma = (float)atof(value);
      } else if (strncmp(token, "q0:", 3)==0) {
        first_spq = atol(value);
      } else if (strncmp(token, "qintsize:", 9)==0) {
        sieve_count = atol(value);
      } else if ((strncmp(token, "skip0:", 6)==0) ||
                 (strncmp(token, "askip:", 6)==0)) {
        sieve_min[0] = atol(value);
      } else if ((strncmp(token, "skip1:", 6)==0) ||
                 (strncmp(token, "rskip:", 6)==0)) {
        sieve_min[1] = atol(value);
      } else if ((strncmp(token, "lim0:", 5)==0) ||
                 (strncmp(token, "alim:", 5)==0)) {
        FB_bound[0] = (float)atol(value);
      } else if ((strncmp(token, "lim1:", 5)==0)||
                 (strncmp(token, "rlim:", 5)==0)) {
        FB_bound[1] = (float)atof(value);
      } else if ((strncmp(token, "lpb0:", 5)==0) ||
                 (strncmp(token, "lpba:", 5)==0)) {
        max_primebits[0] = atoi(value);
      } else if ((strncmp(token, "lpb1:", 5)==0) ||
                 (strncmp(token, "lpbr:", 5)==0)) {
        max_primebits[1] = atoi(value);
      } else if ((strncmp(token, "mfb0:", 5)==0) ||
                 (strncmp(token, "mfba:", 5)==0)) {
        max_factorbits[0] = atoi(value);
      } else if ((strncmp(token, "mfb1:", 5)==0) ||
                 (strncmp(token, "mfbr:", 5)==0)) {
	    value[sizeof(value)-1] = 0;
        max_factorbits[1] = atoi(value);
      } else if ((strncmp(token, "lambda0:", 8)==0) ||
                 (strncmp(token, "alambda:", 8)==0)) {
        sieve_report_multiplier[0] = (float)atof(value);
      } else if ((strncmp(token, "lambda1:", 8)==0) ||
                 (strncmp(token, "rlambda:", 8)==0)) {
        sieve_report_multiplier[1] = (float)atof(value);
      } 
#if 1
        else if (strncmp(token, "n:", 2)==0) { /* redundant, but prevents frustration with e.g. lbpr: */
      } else if (strncmp(token, "m:", 2)==0) {
      } else if ((token[0]=='c') && (token[1] >= '0') && (token[1] <= '8')) {
      } else if ((token[0]=='Y') && (token[1] >= '0') && (token[1] <= '8')) {
      } else if (strncmp(token, "name", 4)==0) {
      } else if (strncmp(token, "type", 4)==0) {
      } else if (strncmp(token, "lss:", 4)==0) {
      } else if (strncmp(token, "deg:", 4)==0) {
      } else {
        printf("Warning: Ignoring input line:\n%s\n", thisLine);
      }
#endif
    }
  }
  fclose(fp);
  return 0;

}

/**************************************************/
void logTotalTime()
/**************************************************/
{ double t=sTime()-sieveStartTime;
  FILE *fp=fopen("ggnfs.log", "a");

  fprintf(fp, "\tLatSieveTime: %ld\n", (long)t);
  fclose(fp);
}

/**************************************************/
int parse_q_from_line(char *buf) {
/**************************************************/
  char *p, *tmp, *next_field;
  u32_t q, q0, i, side;
  static int first=0;

  for(p=tmp=buf; *p && isspace(*p); p++);
  if(!*p) return 0; /* empty line, skip */

  side = (special_q_side == RATIONAL_SIDE) ? 0 : 1;
  for(i=0; *p; p++) {
    if(*p==':') {
      if(i++ == side) tmp = p; /* we will only scan this section for a q0 */
    } else if(!(*p=='-' || *p==',' || isspace(*p) || isxdigit(*p))) {
      if(first++ == 0) printf(" Warning! some corrupt lines in the original file\n");
      return -1;
    }
  }
  if(i!=2) {
    printf(" Warning: an incomplete line in the original file; if just a few, it's ok, they will be skipped\n");
    return -1;           /* must have two ':' some ',' and hexdigits */
  }
  
  q0 = first_spq;
  do {
    q = strtoul(tmp + 1, &next_field, 16);
    if(q >= first_spq && q < first_spq+sieve_count)
      q0 = q;
    tmp = next_field;
  } while(tmp[0] == ',' && isxdigit(tmp[1]));

  /* I've seen cases when q0 is not the last reported in the comma-separated list */
  /* However, the closer it is to the end of the line the more likely it was the true q0 */
  /* In 99% cases it is the last value, but we don't want to depend on that */

  if(q0 > first_spq && q0 < first_spq+sieve_count) {
    sieve_count -= (q0 - first_spq);
    first_spq = q0;
  }
  return 1;
}  

/**************************************************/
int main(int argc, char **argv)
/**************************************************/
{ u16_t zip_output, force_aFBcalc;
  u16_t catch_signals;
  u32_t s, i;

#if defined (_MSC_VER) && defined (_DEBUG)
  int tmpDbgFlag;
  tmpDbgFlag = _CrtSetDbgFlag(_CRTDBG_REPORT_FLAG);
/*
  tmpDbgFlag |= _CRTDBG_CHECK_ALWAYS_DF;
  tmpDbgFlag |= _CRTDBG_CHECK_CRT_DF;
  tmpDbgFlag |= _CRTDBG_DELAY_FREE_MEM_DF;
*/
  tmpDbgFlag |= _CRTDBG_LEAK_CHECK_DF;
  _CrtSetDbgFlag(tmpDbgFlag);
#endif

  n_spq = 0;
  n_spq_discard = 0;

  sieveStartTime = sTime();

#ifdef STC_DEBUG
  debugfile = fopen("rtdsdebug", "wb");
#endif
    mpz_init(g_N); mpz_init(g_m); 
    mpz_init(g_aux1); mpz_init(g_aux2); mpz_init(g_aux3);
    mpz_init(g_sr_a); mpz_init(g_sr_b);
    mpz_init(rational_rest); mpz_init(algebraic_rest);
    mpz_ull_init();

  {
    i32_t option;

    sieve_count = UINT_MAX;
    g_ofile_name = "spairs.out";
    zip_output = 0;
    special_q_side = NO_SIDE;
    sigma = 0;
    keep_factorbase = 0;
    g_resume = 0;
    base_name = NULL;
    first_spq = 0;
    sieve_count = 1;
    force_aFBcalc = 0;
    sysload_cmd = NULL;
    process_no = 0;
    catch_signals = 0;

    first_psp_side = 2;
    first_mpqs_side = 0;
    J_bits = UINT_MAX;

#define NumRead(x) if(sscanf(optarg, "%u" ,(unsigned int*)&x)!=1) Usage()
#define NumRead16(x) if(sscanf(optarg, "%hu" ,(unsigned short*)&x)!=1) Usage()

    while ((option =
            getopt(argc, argv, "FJ:L:M:N:P:RS:ab:c:f:i:kn:o:rst:vz")) != -1) {
      switch (option) {
        case 'R':
          g_resume = 1; break;
        case 'F':
          force_aFBcalc = 1; break;
        case 'J':
          NumRead(J_bits); break;
        case 'L':
          sysload_cmd = optarg; break;
        case 'M':
          NumRead16(first_mpqs_side); break;
        case 'P':
          NumRead16(first_psp_side); break;
        case 'S':
          if (sscanf(optarg, "%f", &sigma) != 1) {
            errprintf("Cannot read floating point number %s\n", optarg);
            Usage();
          }
          break;
        case 'a':
          if (special_q_side != NO_SIDE) {
            errprintf("Ignoring -a\n"); break;
          }
          special_q_side = ALGEBRAIC_SIDE; break;
        case 'b':
            base_name = optarg;
          break;
        case 'c':
          NumRead(sieve_count); break;
        case 'f':
          NumRead(first_spq); break;
        case 'i':
          if (sscanf(optarg, "%hu", &cmdline_first_sieve_side) != 1)
            complain("-i %s ???\n", optarg);
          break;
        case 'k':
          keep_factorbase = 1; break;
        case 'n':
          catch_signals = 1; //break; /* CJM: added `break'. */
        case 'N':
          NumRead(process_no); break;
        case 'o':
          g_ofile_name = optarg; break;
        case 'r':
          if (special_q_side != NO_SIDE) {
            errprintf("Ignoring -r\n"); break;
          }
          special_q_side = RATIONAL_SIDE; break;
        case 's':
          short_output = 1;
          break;
        case 't':
          if (sscanf(optarg, "%hu", &cmdline_first_td_side) != 1)
            complain("-t %s ???\n", optarg);
          break;
        case 'v':
          verbose++; break;
        case 'z':
          zip_output = 1; break;
      }
    }

#define LINE_BUF_SIZE 300

    if (g_resume != 0) {
      char buf[LINE_BUF_SIZE];
      int ret;
      
      if (zip_output != 0)
	complain("Cannot resume gzipped file. gunzip, and retry without -z\n");
      if (g_ofile_name == NULL)
	complain("Cannot resume without the file name\n");
      if (strcmp(g_ofile_name, "-") == 0)
	complain("Cannot resume using stdout\n");
      if ((g_ofile = fopen(g_ofile_name, "ab+")) == NULL)
	complain("Cannot open %s for append: %m\n", g_ofile_name);
      while(fgets(buf, LINE_BUF_SIZE, g_ofile)) {
	ret = parse_q_from_line(buf);
      }
      if(ret < 0) fprintf(g_ofile, "\n"); /* encapsulating the last incomplete line */
      printf(" Resuming with -f %d -c %d\n", first_spq, sieve_count);
    }

    if (J_bits == UINT_MAX)
      J_bits = I_bits - 1;
    if (first_psp_side == 2)
      first_psp_side = first_mpqs_side;
#ifndef I_bits
#error Must #define I_bits
#endif
    if (optind < argc && base_name == NULL) {
      base_name = argv[optind];
      optind++;
    }
    if (base_name == NULL)
      base_name = "gnfs";
    if (optind < argc)
      fprintf(stderr, "Ignoring %u trailing command line args\n",
              argc - optind);

    if (parseJobFile(base_name)) 
      complain("Bad job file: %s\ngiving up...\n", base_name);

    last_spq = first_spq + sieve_count;
    if (last_spq >= INT_MAX / 2) {
      complain("Cannot handle special q >= %d\n", INT_MAX / 2);
    }

    for (i = 0; i < 2; i++) {
      if (FB_bound[i] < 4 || sieve_report_multiplier[i] <= 0) {
        complain("Please set all bounds to reasonable values!\n");
      }
      if (max_primebits[i] > 33) {
        complain("Only large primes up to 33 bits are allowed.\n");
      }
    }
    if (sieve_count != 0) {
      if (sigma == 0)
        complain("Please set a skewness\n");
      if (special_q_side == NO_SIDE) {
        errprintf("Please use -a or -r\n");
        Usage();
      }

	  if (special_q_side >= sizeof(FB_bound)/sizeof(FB_bound[0]))
	  {
		  assert(0);
		  complain("Invalid special_q_side value!\n", special_q_side);
#ifdef _MSC_VER
		  exit(1); // just to make PreFAST shutup
#endif
	  }

      if (FB_bound[special_q_side] > first_spq) {
	FB_bound[special_q_side] = (float) first_spq-1;
	printf(" Warning:  lowering FB_bound to %u.\n",first_spq-1);
        /* complain("Special q lower bound %u below FB bound %g\n", 
	   first_spq, FB_bound[special_q_side]); */
      }
    }
    if (g_poldeg[0] < g_poldeg[1])
      poldeg_max = g_poldeg[1];
    else
      poldeg_max = g_poldeg[0];

    i_shift = 1 << (I_bits - 1);
    n_I = 1 << I_bits;
    n_i = n_I / 2;
    i_bits = I_bits - 1;
    n_J = 1 << J_bits;
    n_j = n_J / 2;
    j_bits = J_bits - 1;

    {
      u32_t  j;
      double x, y, z;

      x = sqrt(first_spq * sigma) * n_I;
      y = x / sigma;
      for (j = 0; j < 2; j++) {
        poly_f[j] = xmalloc((g_poldeg[j] + 1) * sizeof(*poly_f[j]));

        for (i = 0, z = 1, poly_norm[j] = 0; i <= g_poldeg[j]; i++) {
          poly_f[j][i] = mpz_get_d(g_poly[j][i]);
          poly_norm[j] = poly_norm[j] * y + fabs(poly_f[j][i]) * z;
          z *= x;
        }
      }
    }
  }


  siever_init();

  if (sieve_count != 0) {
    if (g_ofile_name == NULL) {
      if (zip_output == 0) {
        asprintf(&g_ofile_name, "%s.lasieve-%u.%u-%u", base_name,
                 special_q_side, first_spq, last_spq);
      } else {
        asprintf(&g_ofile_name,"gzip --best --stdout > %s.lasieve-%u.%u-%u.gz", 
                  base_name, special_q_side, first_spq, last_spq);
      }
    } else {
      if (strcmp(g_ofile_name, "-") == 0) {
        if (zip_output == 0) {
          g_ofile = stdout;
          g_ofile_name = "to stdout";
          goto done_opening_output;
        } else
          g_ofile_name = "gzip --best --stdout";
      } else {
          zip_output = 0;
      }
    }
    if (zip_output == 0) {
      if (g_resume != 0) {
        goto done_opening_output;
      }
      if ((g_ofile = fopen(g_ofile_name, "rb")) != NULL)
        complain(" Will not overwrite existing file %s for output; rename it, move it away, or use -R option (resume)\n", g_ofile_name);
      /* CJM: consider making this a "a". */
      if ((g_ofile = fopen(g_ofile_name, "wb")) == NULL)
        complain("Cannot open %s for output: %m\n", g_ofile_name);
    } else {
      if ((g_ofile = popen(g_ofile_name, "w")) == NULL)
        complain("Cannot exec %s for output: %m\n", g_ofile_name);
    }
  done_opening_output:
  /*    fprintf(ofile, "F 0 X %u 1\n", poldeg[0]); */
    ;
  }
  

  getFB(force_aFBcalc);


  sieve_interval = xvalloc(L1_SIZE);
  cand = xvalloc(L1_SIZE * sizeof(*cand));
  fss_sv = xvalloc(L1_SIZE);
  tiny_sieve_buffer = xmalloc(TINY_SIEVEBUFFER_SIZE);
  if (n_i > L1_SIZE)
    complain("Strip length %u exceeds L1 size %u\n", n_i, L1_SIZE);
  j_per_strip = L1_SIZE / n_i;
  jps_bits = L1_BITS - i_bits;
  jps_mask = j_per_strip - 1;
  if (j_per_strip != 1 << jps_bits)
    Schlendrian("Expected %u j per strip, calculated %u\n",j_per_strip,1<<jps_bits);
  n_strips = n_j >> (L1_BITS - i_bits);
  rec_info_init(n_i, n_j);


  horizontal_sievesums = xmalloc(j_per_strip * sizeof(*horizontal_sievesums));
  for (s = 0; s < 2; s++) {
    u32_t fbi;
    size_t maxent;

    smallsieve_aux[s] = xmalloc(4 * fbis[s] * sizeof(*(smallsieve_aux[s])));
#ifdef PREINVERT
    smalltd_pi[s] = xmalloc(fbis[s] * sizeof(*(smalltd_pi[s])));
#endif
    smalltdsieve_aux[s] = xmalloc(j_per_strip*sizeof(*(smalltdsieve_aux[s])));
    for (fbi = 0; fbi < j_per_strip; fbi++)
      smalltdsieve_aux[s][fbi] = xmalloc(fbis[s]*sizeof(**(smalltdsieve_aux[s])));
    smallsieve_aux1[s] = xmalloc(6*xFBs[s]*sizeof(*(smallsieve_aux1[s])));

    maxent = fbis[s];
    maxent += xFBs[s];
    smallpsieve_aux[s] = xmalloc(3*maxent*sizeof(*(smallpsieve_aux[s])));
    maxent = 0;
    for (fbi = 0; fbi < xFBs[s]; fbi++) {
      if (xFB[s][fbi].p == 2)
        maxent++;
    }
    smallsieve_aux2[s] = xmalloc(4*maxent*sizeof(*(smallsieve_aux2[s])));
    x2FB[s] = xmalloc(maxent * 6 * sizeof(*(x2FB[s])));
  }

#ifdef GCD_SIEVE_BOUND
  {
    u32_t p;

    firstprime32(&special_q_ps);
    np_gcd_sieve = 0;
    for (p=nextprime32(&special_q_ps); p<GCD_SIEVE_BOUND; p=nextprime32(&special_q_ps))
      np_gcd_sieve++;
    gcd_sieve_buffer = xmalloc(2 * np_gcd_sieve * sizeof(*gcd_sieve_buffer));

    firstprime32(&special_q_ps);
    i = 0;
    for (p=nextprime32(&special_q_ps); p<GCD_SIEVE_BOUND; p=nextprime32(&special_q_ps))
      gcd_sieve_buffer[2*i++] = p;
  }
#endif

  for (s = 0; s < 2; s++) {
    if (sieve_min[s] < TINY_SIEVE_MIN && sieve_min[s] != 0) {
      errprintf("Sieving with all primes on side %u since\n", s);
      errprintf("tiny sieve procedure is being used\n");
      sieve_min[s] = 0;
    }
    current_ij[s] = xmalloc(FBsize[s] * sizeof(*current_ij[s]));
    LPri[s] = xmalloc(FBsize[s] * sizeof(**LPri) * RI_SIZE);
  }

  {
    size_t total_alloc;
    u16_t *sched_buf;
    double pvl_max[2];

    total_alloc = 0;
    for (s = 0; s < 2; s++) {
      u32_t fbi_lb;

      if (sigma >= 1)
        pvl_max[s] = g_poldeg[s] * log(last_spq * sqrt(sigma));
      else
        pvl_max[s] = g_poldeg[s] * log(last_spq / sqrt(sigma));
      pvl_max[s] += log(poly_norm[s]);
      if (fbi1[s] >= FBsize[s] || i_bits + j_bits <= L1_BITS) {
        n_schedules[s] = 0;
        continue;
      }
      for (i = 0; i < N_PRIMEBOUNDS; i++)
        if (FB_bound[s] <= schedule_primebounds[i] || i_bits+j_bits <= schedule_sizebits[i])
          break;
      n_schedules[s] = i + 1;
      schedules[s] = xmalloc(n_schedules[s] * sizeof(**schedules));
      fbi_lb = fbi1[s];
      for (i = 0; i < n_schedules[s]; i++) {
  	    u32_t fbp_ub;
		u32_t fbp_lb;
        u32_t fbi, fbi_ub;
        u32_t sp_i;
        u32_t n, sl_i;
        u32_t ns;
        size_t allocate, all1;

        if (i == n_schedules[s] - 1)
		{
		  assert(FB_bound[s] < UINT_MAX);
          fbp_ub = (u32_t)FB_bound[s];
		}
        else
          fbp_ub = schedule_primebounds[i];
        if (i == 0)
          fbp_lb = FB[s][fbi1[s]];
        else
          fbp_lb = schedule_primebounds[i - 1];

        if (i_bits + j_bits < schedule_sizebits[i])
          ns = 1 << (i_bits + j_bits - L1_BITS);
        else
          ns = 1 << (schedule_sizebits[i] - L1_BITS);
        schedules[s][i].n_strips = ns;

#define SCHED_PAD 32
#define SCHED_TOL 1.2
		assert(rint(SCHED_PAD + SCHED_TOL * n_i * j_per_strip * log(log(fbp_ub) / log(fbp_lb))) <= ULONG_MAX);
        allocate = (size_t)rint(SCHED_PAD + SCHED_TOL * n_i * j_per_strip * log(log(fbp_ub) / log(fbp_lb)));
        allocate *= SE_SIZE;

		assert(((double)allocate + n_i * ceil(pvl_max[s] / log(fbp_lb)) * SE_SIZE) <= ULONG_MAX);
        all1 = allocate + (size_t)(n_i * ceil(pvl_max[s] / log(fbp_lb)) * SE_SIZE);
        schedules[s][i].alloc = allocate;
        schedules[s][i].alloc1 = all1;

        for (n = 0, fbi = fbi_lb; fbi < FBsize[s];) {
          u32_t fbi_ub1;

          fbi_ub1 = fbi + SCHEDFBI_MAXSTEP;
          if (fbi_ub1 >= FBsize[s])
            fbi_ub1 = FBsize[s];
          else {
            if (FB[s][fbi_ub1] > fbp_ub) {
              while (FB[s][fbi_ub1] > fbp_ub)
                fbi_ub1--;
              fbi_ub1++;
            }
          }
          if (FB_logs[s][fbi] == FB_logs[s][fbi_ub1 - 1]) {
            n++;
            fbi = fbi_ub1;
          } else {
            u32_t l;

            n += FB_logs[s][fbi_ub1 - 1] - FB_logs[s][fbi];
            fbi = fbi_ub1 - 1;
            l = FB_logs[s][fbi];
            while (FB_logs[s][fbi] == l)
              fbi--;
            fbi++;
          }
          if (fbi >= FBsize[s] || FB[s][fbi] > fbp_ub)
            break;
        }
        fbi_ub = fbi;
        schedules[s][i].n_pieces = n;
        n++;
        schedules[s][i].schedule =
          xmalloc(n * sizeof(*(schedules[s][i].schedule)));
        for (sl_i = 0; sl_i < n; sl_i++)
          schedules[s][i].schedule[sl_i] =
            xmalloc(ns * sizeof(**(schedules[s][i].schedule)));
        schedules[s][i].schedule[0][0] = (u16_t *) total_alloc;
        total_alloc += all1;
        for (sp_i = 1; sp_i < ns; sp_i++) {
          schedules[s][i].schedule[0][sp_i] = (u16_t *) total_alloc;
          total_alloc += allocate;
        }
        schedules[s][i].fbi_bounds =
          xmalloc(n * sizeof(*(schedules[s][i].fbi_bounds)));
        schedules[s][i].schedlogs = xmalloc(n);
        for (n = 0, fbi = fbi_lb; fbi < fbi_ub;) {
          u32_t fbi_ub1;

          fbi_ub1 = fbi + SCHEDFBI_MAXSTEP;
          if (fbi_ub1 > fbi_ub)
            fbi_ub1 = fbi_ub;
          if (FB_logs[s][fbi] == FB_logs[s][fbi_ub1 - 1]) {
            schedules[s][i].fbi_bounds[n++] = fbi;
            fbi = fbi_ub1;
          } else {
            u32_t l, lmax;

            lmax = FB_logs[s][fbi_ub1 - 1];
            for (l = FB_logs[s][fbi]; l < lmax; l++) {
              schedules[s][i].fbi_bounds[n++] = fbi;
              while (fbi < fbi_ub && FB_logs[s][fbi] == l)
                fbi++;
            }
          }
        }
        if (n != schedules[s][i].n_pieces)
          Schlendrian("Expected %u schedule pieces on side %u, have %u\n",
                      schedules[s][i].n_pieces, s, n);
        schedules[s][i].fbi_bounds[n++] = fbi;
        for (n = 0; n < schedules[s][i].n_pieces; n++)
          schedules[s][i].schedlogs[n] =
            FB_logs[s][schedules[s][i].fbi_bounds[n]];
        schedules[s][i].ri =
          LPri[s] + (schedules[s][i].fbi_bounds[0] - fbis[s]) * RI_SIZE;
        fbi_lb = fbi_ub;
      }
    }

    sched_buf = xmalloc((total_alloc + 65536 * SE_SIZE * j_per_strip) *
                        sizeof(**((**schedules).schedule)));
    for (s = 0; s < 2; s++) {

      for (i = 0; i < n_schedules[s]; i++) {
        u32_t sp_i;

        for (sp_i = 0; sp_i < schedules[s][i].n_strips; sp_i++)
          schedules[s][i].schedule[0][sp_i] =
            sched_buf + (size_t) (schedules[s][i].schedule[0][sp_i]);
      }
    }

#ifdef USE_MEDSCHED
    for (s = 0; s < 2; s++) {
      if (fbis[s] < fbi1[s]) {
        u32_t fbi;
        u32_t n;
        unsigned char oldlog;

        medsched_alloc[s] = j_per_strip * (fbi1[s] - fbis[s]) * SE_SIZE;

		assert(((double)medsched_alloc[s] + n_i * ceil(pvl_max[s] / log(n_i)) * SE_SIZE) <= ULONG_MAX);
        medsched_alloc[s] += (size_t)(n_i * ceil(pvl_max[s] / log(n_i)) * SE_SIZE);
        n_medsched_pieces[s] =
          1 + FB_logs[s][fbi1[s] - 1] - FB_logs[s][fbis[s]];
        med_sched[s] =
          xmalloc((1 + n_medsched_pieces[s]) * sizeof(**med_sched));
        med_sched[s][0] = xmalloc(medsched_alloc[s] * sizeof(***med_sched));

        medsched_fbi_bounds[s] =
          xmalloc((1 + n_medsched_pieces[s]) * sizeof(**medsched_fbi_bounds));
        medsched_logs[s] = xmalloc(n_medsched_pieces[s]);

        for (n = 0, fbi = fbis[s], oldlog = UCHAR_MAX; fbi < fbi1[s]; fbi++) {
          if (FB_logs[s][fbi] != oldlog) {
            medsched_fbi_bounds[s][n] = fbi;
            oldlog = FB_logs[s][fbi];
            medsched_logs[s][n++] = oldlog;
          }
        }
        if (n != n_medsched_pieces[s])
          Schlendrian("Expected %u medium schedule pieces on side %u, have %u\n",
                       n_medsched_pieces[s], s, n);
        medsched_fbi_bounds[s][n] = fbi;
      } else {

        n_medsched_pieces[s] = 0;
      }
    }
#endif
  }

  {
    size_t schedbuf_alloc;

    for (s = 0, schedbuf_alloc = 0; s < 2; s++) {
      for (i = 0; i < n_schedules[s]; i++)
        if (schedules[s][i].n_pieces > schedbuf_alloc)
          schedbuf_alloc = schedules[s][i].n_pieces;
    }
    schedbuf = xmalloc((1 + schedbuf_alloc) * sizeof(*schedbuf));
  }

  td_buf1 = xmalloc((1 + L1_SIZE) * sizeof(*td_buf1));
  td_buf[0] = xmalloc(td_buf_alloc[0] * sizeof(**td_buf));
  td_buf[1] = xmalloc(td_buf_alloc[1] * sizeof(**td_buf));

  if (tds_fbi == NULL) {
    tds_fbi = xmalloc(UCHAR_MAX * sizeof(*tds_fbi));
    tds_fbi_curpos = xmalloc(UCHAR_MAX * sizeof(*tds_fbi));
    for (i = 0; i < UCHAR_MAX; i++)
      tds_fbi[i] = xmalloc(tds_fbi_alloc * sizeof(**tds_fbi));
  }

  for (i = 0; i < L1_SIZE; i++) {
    mpz_init(td_rests[i]);
  }
  for (s = 0; s < 2; s++) {
    mpz_init(large_factors[s]);
    large_primes[s] =
      xmalloc(max_factorbits[s] * sizeof(*(large_primes[s])));

	assert((max_factorbits[s] >= 0) && (max_factorbits[s] <= INT_MAX));
    for (i = 0; i < (u32_t)max_factorbits[s]; i++) {
      mpz_init(large_primes[s][i]);
    }
    mpz_init_set_d(FBb_sq[s], FB_bound[s]);
    mpz_mul(FBb_sq[s], FBb_sq[s], FBb_sq[s]);
  }

  all_spq_done = 1;

  if (catch_signals != 0) {
    signal(SIGTERM, terminate_sieving);
    signal(SIGINT, terminate_sieving);
  }

  lasieve(); /* CJM, 6/17/04. */

  if (sieve_count != 0) {
    if (zip_output != 0)
      pclose(g_ofile);
    else
      fclose(g_ofile);
  }
  logbook(0, "%u Special q, %u reduction iterations\n", n_spq, n_iter);

  if (n_spq_discard > 0)
    logbook(0, "%u Special q discarded\n", n_spq_discard);

  {
    u32_t side;

    logbook(0, "reports: %u->%u->%u->%u->%u->%u\n",
            n_prereports, n_reports, n_rep1, n_rep2,
            n_tdsurvivors[first_td_side], n_tdsurvivors[1 - first_td_side]);
    logbook(0,
            "Number of relations with k rational and l algebraic primes for (k,l)=:\n");
    logbook(0, "\nTotal yield: %u\n", yield);
    if (n_mpqsfail[0] != 0 || n_mpqsfail[1] != 0 ||
        n_mpqsvain[0] != 0 || n_mpqsvain[1] != 0) {
      logbook(0, "%u/%u mpqs failures, %u/%u vain mpqs\n", n_mpqsfail[0],
              n_mpqsfail[1], n_mpqsvain[0], n_mpqsvain[1]);
    }
    logbook(0, "milliseconds total: Sieve %u Sched %u medsched %u\n",
            sieve_clock, Schedule_clock, medsched_clock);
    logbook(0, "TD %u (Init %u, MPQS %u) Sieve-Change %u\n",
            td_clock, tdi_clock, mpqs_clock, sch_clock);
    for (side = 0; side < 2; side++) {
      logbook(0,
              "TD side %u: init/small/medium/large/search: %u %u %u %u %u\n",
              side, tdsi_clock[side], tds1_clock[side], tds2_clock[side],
              tds3_clock[side], tds4_clock[side]);
      logbook(0, "sieve: init/small/medium/large/search: %u %u %u %u %u\n",
              si_clock[side], s1_clock[side], s2_clock[side], s3_clock[side],
              cs_clock[side]);
    }
#ifdef ZSS_STAT
    fprintf(stderr,
            "%u subsieves, zero: %u first sieve, %u second sieve %u first td\n",
            nss, nzss[0], nzss[1], nzss[2]);
#endif
  }

  logTotalTime();
  if (special_q >= last_spq && all_spq_done != 0)
    exit(0);
  exit(1);
}



/****************************************************************/
u16_t get_plog_lb(unsigned char *lb_array, i16_t ** root_array,
            size_t * root_array_alloc, double *poly, i32_t poldeg, double sm)
/****************************************************************/
{ u16_t nroots;
  u32_t k;

  for (k = 0, nroots = 0; k < n_J / CANDIDATE_SEARCH_STEPS; k++) {
    i32_t l, m, n;
    double pv_min;

    l = 2 * k * CANDIDATE_SEARCH_STEPS - n_J;
    n = l + 2 * CANDIDATE_SEARCH_STEPS;
    pv_min = -1;
    for (m = n; l < n; l = m, m = n) {
      double y;

      while ((y = rpol_lb(poly, poldeg, (double) l / (double) n_J,
             (double) m / (double) n_J)) < 0) {

        if (m > l + 1)
          m = (l + m) / 2;
        else {
          adjust_bufsize((void **) root_array, root_array_alloc,
                         1 + nroots, 1, sizeof(**root_array));
          (*root_array)[nroots++] = l;
          break;
        }
      }
      if (y > 0 && (y < pv_min || pv_min < 0))
        pv_min = y;
    }
    if (pv_min > 0)
	{
      assert(rint(sm * log(pv_min)) <= UCHAR_MAX);
      lb_array[k] = (unsigned char)rint(sm * log(pv_min));
	}
    else
      lb_array[k] = 0;
  }
  return nroots;
}

/****************************************************************/
void do_scheduling(struct schedule_struct *sched, u32_t ns, u32_t ot, u32_t s)
/****************************************************************/
{ u32_t ll, n1_j, *ri;

  n1_j = ns << (L1_BITS - i_bits);
  for (ll = 0, ri = sched->ri; ll < sched->n_pieces; ll++) {
    u32_t fbi_lb, fbi_ub, fbio;

    memcpy(sched->schedule[ll + 1], sched->schedule[ll], ns*sizeof(u16_t **));
    fbio = sched->fbi_bounds[ll];
    fbi_lb = fbio;
    fbi_ub = sched->fbi_bounds[ll + 1];
#ifdef SCHEDULING_FUNCTION_CALCULATES_RI
    if (ot == 1)
      lasieve_setup(FB[s] + fbi_lb, proots[s] + fbi_lb, fbi_ub - fbi_lb,
                    a0, a1, b0, b1, LPri[s] + (fbi_lb - fbis[s]) * RI_SIZE);
#endif
    ri = lasched(ri, current_ij[s] + fbi_lb, current_ij[s] + fbi_ub, n1_j, 
                 (u32_t **) (sched->schedule[ll + 1]), fbi_lb - fbio, ot);

    { u32_t k;

      for (k = 0; k < ns; k++)
        if (sched->schedule[ll + 1][k] >= sched->schedule[0][k] + sched->alloc) {
          if (k == 0 && sched->schedule[ll + 1][k] < sched->schedule[0][k] + sched->alloc1)
            continue;
          longjmp(termination_jb, SCHED_PATHOLOGY);
        }
    }

  }
}
#else
#define BADSCHED
#endif

#ifdef GCD_SIEVE_BOUND
/****************************************************************/
static void gcd_sieve() 
/****************************************************************/
{ u32_t i;

  for (i = 0; i < np_gcd_sieve; i++) {
    u32_t x, p;

    x = gcd_sieve_buffer[2 * i + 1];
    p = gcd_sieve_buffer[2 * i];
    while (x < j_per_strip) {
      unsigned char *z, *z_ub;

      z = sieve_interval + (x << i_bits);
      z_ub = z + n_i - 3 * p;
      z +=  oddness_type == 2 ? (n_i / 2) % p : ((n_i + p - 1) / 2) % p;
      while (z < z_ub) {
        *z = 0;
        *(z + p) = 0;
        z += 2 * p;
        *z = 0;
        *(z + p) = 0;
        z += 2 * p;
      }
      z_ub += 3 * p;
      while (z < z_ub) {
        *z = 0;
        z += p;
      }
      x = x + p;
    }
    gcd_sieve_buffer[2 * i + 1] = x - j_per_strip;
  }
}
#endif

/****************************************************************/
static void xFBtranslate(u16_t *rop, xFBptr op) 
/****************************************************************/
{ u32_t x, y, am, bm, rqq;

  modulo32 = op->pp;
  rop[3] = op->l;
  am = a1 > 0 ? ((u32_t) a1) % modulo32 : modulo32 - ((u32_t) (-a1)) % modulo32;
  if (am == modulo32)
    am = 0;
  bm = b1 > 0 ? ((u32_t) b1) % modulo32 : modulo32 - ((u32_t) (-b1)) % modulo32;
  if (bm == modulo32)
    bm = 0;
  x = modsub32(modmul32(op->qq, am), modmul32(op->r, bm));
  am = a0 > 0 ? ((u32_t) a0) % modulo32 : modulo32 - ((u32_t) (-a0)) % modulo32;
  if (am == modulo32)
    am = 0;
  bm = b0 > 0 ? ((u32_t) b0) % modulo32 : modulo32 - ((u32_t) (-b0)) % modulo32;
  if (bm == modulo32)
    bm = 0;
  y = modsub32(modmul32(op->r, bm), modmul32(op->qq, am));
  rqq = 1;
  if (y != 0) {
    while (y % (op->p) == 0) {
      y = y / (op->p);
      rqq *= op->p;
    }
  } else {
    rqq = op->pp;
  }
  modulo32 = modulo32 / rqq;
  rop[0] = modulo32;
  rop[1] = rqq;
  if (modulo32 > 1)
    rop[2] = modmul32(modinv32(y), x);
  else
    rop[2] = 0;
  rop[4] = op->l;
}

/****************************************************************/
static int xFBcmp(const void *opA, const void *opB) 
/****************************************************************/
{ xFBptr op1, op2;

  op1 = (xFBptr) opA;
  op2 = (xFBptr) opB;
  if (op1->pp < op2->pp)
    return -1;
  if (op1->pp == op2->pp)
    return 0;
  return 1;
}

/****************************************************************/
static u32_t add_primepowers2xaFB(size_t * xaFB_alloc_ptr, u32_t pp_bound,
                                   u32_t s, u32_t p, u32_t r) 
/****************************************************************/
{ u32_t a, b, q, qo, *rbuf, nr=0, *Ar, exponent, init_xFB;
  size_t rbuf_alloc;

  if (xFBs[s] == 0 && p == 0)
    Schlendrian("add_primepowers2xaFB on empty xaFB\n");

  rbuf_alloc = 0;
  Ar = xmalloc((1 + g_poldeg[s]) * sizeof(*Ar));

  if (p != 0) {
    init_xFB = 0;
    q = p;
    if (r == p) {
      a = 1;
      b = p;
    } else {
      a = r;
      b = 1;
    }
  } else {
    init_xFB = 1;
    q = xFB[s][xFBs[s] - 1].pp;
    p = xFB[s][xFBs[s] - 1].p;
    a = xFB[s][xFBs[s] - 1].r;
    b = xFB[s][xFBs[s] - 1].qq;
  }

  qo = q;
  exponent = 1;
  for (;;) {
    u32_t j;

    if (q > pp_bound / p)
      break;
    modulo32 = p * q;
    for (j = 0; j <= g_poldeg[s]; j++)
      Ar[j] = mpz_fdiv_ui(g_poly[s][j], modulo32);
    if (b == 1) {
      for (r = a, nr = 0; r < modulo32; r += qo) {
        u32_t pv;

        for (j = 1, pv = Ar[g_poldeg[s]]; j <= g_poldeg[s]; j++) {
          pv = modadd32(Ar[g_poldeg[s] - j], modmul32(pv, r));
        }
        if (pv == 0) {
          adjust_bufsize((void **) &rbuf, &rbuf_alloc, 1 + nr, 4, sizeof(*rbuf));
          rbuf[nr++] = r;
        } else if (pv % q != 0)
          Schlendrian("xFBgen: %u not a root mod %u\n", r, q);
      }
    } else {
      for (r = (modmul32(b, modinv32(a))) % qo, nr = 0; r < modulo32; r += qo) {
        u32_t pv;

        for (j = 1, pv = Ar[0]; j <= g_poldeg[s]; j++) {
          pv = modadd32(Ar[j], modmul32(pv, r));
        }
        if (pv == 0) {
          adjust_bufsize((void **) &rbuf, &rbuf_alloc, 1 + nr, 4, sizeof(*rbuf));
          rbuf[nr++] = r;
        } else if (pv % q != 0)
          Schlendrian("xFBgen: %u^{-1} not a root mod %u\n", r, q);
      }
    }

    if (qo * nr != modulo32)
      break;
    q = modulo32;
    exponent++;
  }
  if (init_xFB != 0)
  {
    assert(rint(sieve_multiplier[s]*log(q)) - rint(sieve_multiplier[s] * log(qo / p)) <= UINT_MAX);
    xFB[s][xFBs[s] - 1].l = (u32_t)(rint(sieve_multiplier[s]*log(q)) - rint(sieve_multiplier[s] * log(qo / p)));
  }

  if (q <= pp_bound / p) {
    u32_t j;

    for (j = 0; j < nr; j++) {
      xFBptr f;

      adjust_bufsize((void **) &(xFB[s]), xaFB_alloc_ptr, 1+xFBs[s], 16, sizeof(**xFB));
      f = xFB[s] + xFBs[s];
      f->p = p;
      f->pp = q * p;
      if (b == 1) {
        f->qq = 1;
        f->r = rbuf[j];
        f->q = f->pp;
      } else {
        modulo32 = (q * p) / b;
        rbuf[j] = rbuf[j] / b;
        if (rbuf[j] == 0) {
          f->qq = f->pp;
          f->q = 1;
          f->r = 1;
        } else {
          while (rbuf[j] % p == 0) {
            rbuf[j] = rbuf[j] / p;
            modulo32 = modulo32 / p;
          }
          f->qq = (f->pp) / modulo32;
          f->q = modulo32;
          f->r = modinv32(rbuf[j]);
        }
      }
      xFBs[s]++;
      add_primepowers2xaFB(xaFB_alloc_ptr, pp_bound, s, 0, 0);
    }
  }
  if (rbuf_alloc > 0)
    free(rbuf);
  free(Ar);
  return exponent;
}

#ifdef OFMT_CWI
/****************************************************************/
static char ulong2cwi(u32_t n) 
/****************************************************************/
{
  if (n < 10) return '0' + n;
  n = n - 10;
  if (n < 26) return 'A' + n;
  n = n - 26;
  if (n < 26) return 'a' + n;
  return '\0';
}
#endif

#ifdef DEBUG
/****************************************************************/
int mpout(mpz_t X) 
/****************************************************************/
{
  mpz_out_str(stdout, 10, X);
  puts("");
  return 1;
}
#endif

/****************************************************************/
void dumpsieve(u32_t j_offset, u32_t side) 
/****************************************************************/
{ FILE *ofile;
  char *ofn;

  asprintf(&ofn, "sdump4e.ot%u.j%u.s%u", oddness_type, j_offset, side);
  if ((ofile = fopen(ofn, "wb")) == NULL) {
    free(ofn);
    return;
  }
  fwrite(sieve_interval, 1, L1_SIZE, ofile);
  fclose(ofile);
  free(ofn);
  asprintf(&ofn, "hzsdump4e.ot%u.j%u.s%u", oddness_type, j_offset, side);
  if ((ofile = fopen(ofn, "wb")) == NULL) {
    free(ofn);
    return;
  }
  fwrite(horizontal_sievesums, 1, j_per_strip, ofile);
  fclose(ofile);
  free(ofn);
} 
