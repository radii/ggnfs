/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/
/* Written by T. Kleinjung, modified by J. Franke */
/* CJM, 11/30/04: This is a hack of generic/mpqs.c for easy compilation
   and inclusion in GGNFS.
*/
/*
#define MPQS_STAT
*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <memory.h>
#include <gmp.h>

#include "lasieve.h"

#ifdef MPQS_STAT
extern u64_t stat_asm_eval, stat_asm_td;
extern u32_t stat_asm_div, stat_final_mulmod;
#endif

#define MPQS_MAX_FBSIZE      512
#define MPQS_MIN_EXCESS      10
#define MPQS_MAX_ADIV_ALL    6
#define MPQS_SIEVELEN        (1<<L1_BITS)
#define MPQS_REL_ENTRIES     32
#define MPQS_TD_MAX_NDIV     (MPQS_REL_ENTRIES-5)
#define MPQS_MAX_NRELS       256  /* ??? */
#define MPQS_MAX_NPRELS      2048  /* ??? */
#define MPQS_MAX_NCRELS      256  /* ??? */
#define MPQS_MAX_NRELS_BUF   1024  /* ??? */
#define MPQS_MIN_RELBUFFER   16
#define MPQS_GAUSS_MAX       512
#define MPQS_MAX_FACTORS     5
#define MPQS_MAX_NPRIMES     1024
#define MPQS_HASH_OVERFLOW   65535
#define MPQS_MULT_NTESTPR    25       /* prime(MPQS_MULT_NTESTPR)<2^8 */
#define MPQS_MULT_NCAND      8        /* largest candidate<2^8 */
#define MPQS_FB_MAXPRIME     4096     /* used in sieve and final */

/* common with asm functions */
static u16_t mpqs_nFBk_1;                          /* = mpqs_nFBk + 1 */
static u16_t mpqs_td_begin, mpqs_sievebegin;
static u32_t mpqs_FB_inv_info[4*MPQS_MAX_NPRIMES];
u16_t mpqs_nFBk, mpqs_FBk[MPQS_MAX_NPRIMES];       /* list of primes p that divides kN */
unsigned char mpqs_FB_log[MPQS_MAX_FBSIZE];
u16_t mpqs_nFB, mpqs_nAdiv_total;
u16_t mpqs_Adiv_all[MPQS_MAX_ADIV_ALL];
u32_t mpqs_FBk_inv[3], mpqs_FB_A_inv[MPQS_MAX_ADIV_ALL];
/* end */

static int mpqs_isinit = 0;
static mpz_t mpqs_gmp_help, mpqs_dummy, mpqs_help;
static mpz_t mpqs_N, mpqs_kN, mpqs_factors[MPQS_MAX_FACTORS];
static unsigned char mpqs_f_status[MPQS_MAX_FACTORS];
static u32_t mpqs_multiplier, mpqs_nfactors;
static double mpqs_complexity;
static u16_t mpqs_prime_table[MPQS_MAX_NPRIMES];
static u16_t mpqs_prime_sqrthelp[MPQS_MAX_NPRIMES];
static u16_t mpqs_nAdiv, mpqs_accept;
static unsigned char mpqs_Adiv_log[MPQS_MAX_ADIV_ALL];;
static u16_t mpqs_pmax; /* maximum prime p in our factor base such that, (kN, p) == 1. */
static u16_t mpqs_Adiv[MPQS_MAX_ADIV_ALL];
static u16_t mpqs_Adiv_sqrt[MPQS_MAX_ADIV_ALL], mpqs_Adiv_all_sqrt[MPQS_MAX_ADIV_ALL];
static u16_t mpqs_Adiv_start1[MPQS_MAX_ADIV_ALL], mpqs_Adiv_start2[MPQS_MAX_ADIV_ALL];
static unsigned char mpqs_Adiv_active[MPQS_MAX_ADIV_ALL];
static unsigned char mpqs_A_mask[(1<<MPQS_MAX_ADIV_ALL)];
static u16_t mpqs_nA;
static i16_t mpqs_A_index, mpqs_B_index;
static u16_t mpqs_SI_add[MPQS_MAX_ADIV_ALL][MPQS_MAX_FBSIZE];
static u16_t mpqs_2A_all_inv[MPQS_MAX_FBSIZE];
static u16_t mpqs_Adiv_SI_add[MPQS_MAX_ADIV_ALL][MPQS_MAX_ADIV_ALL];
static i64_t mpqs_A, mpqs_B, mpqs_C, mpqs_2B, mpqs_Bi[MPQS_MAX_ADIV_ALL];
static u32_t mpqs_disp;
static u16_t mpqs_nsurvivors;
static u32_t mpqs_nrels, mpqs_nlp, mpqs_excess;
static u32_t mpqs_nsp; /* = 1 + mpqs_nFB + mpqs_nFBk */
static u16_t mpqs_nprels, mpqs_ncrels;
static u16_t **mpqs_relations, **mpqs_part_rels, **mpqs_comb_rels;
static u16_t **mpqs_rel_buffer;
static u16_t mpqs_rel_hash_table[256][16];
static u16_t mpqs_hash_table[128][16];
static u32_t mpqs_hash_index[128][16];
static u32_t mpqs_lp[128*15];
static u32_t **mpqs_gauss_row;
static unsigned char mpqs_sol[MPQS_GAUSS_MAX];
static i16_t mpqs_gauss_c[MPQS_GAUSS_MAX], mpqs_gauss_d[MPQS_GAUSS_MAX];
static u32_t mpqs_gauss_m, mpqs_gauss_n, mpqs_gauss_n32, mpqs_gauss_k;
static i16_t mpqs_exp[128*16+MPQS_MAX_FBSIZE+3+MPQS_MAX_ADIV_ALL];
static mpz_t mpqs_sq1, mpqs_sq2;
static double mpqs_kN_dbl;
static u64_t mpqs_kN_64, mpqs_A_inv_64;

static u16_t mpqs_FB[2*MPQS_MAX_FBSIZE];
static u16_t mpqs_FB_start[2*MPQS_MAX_FBSIZE];
static u32_t mpqs_sievelen;
static unsigned char *mpqs_sievearray;
u32_t mpqs_FB_inv[MPQS_MAX_FBSIZE];

unsigned char mpqs_256_inv_table[128]=
{
	1,   171, 205, 183, 57,  163, 197, 239, 
	241, 27,  61,  167, 41,  19,  53,  223, 
	225, 139, 173, 151, 25,  131, 165, 207, 
	209, 251, 29,  135, 9,   243, 21,  191, 
	193, 107, 141, 119, 249, 99,  133, 175, 
	177, 219, 253, 103, 233, 211, 245, 159, 
	161, 75,  109, 87,  217, 67,  101, 143, 
	145, 187, 221, 71,  201, 179, 213, 127, 
	129, 43,  77,  55,  185, 35,  69,  111, 
	113, 155, 189, 39,  169, 147, 181, 95, 
	97,  11,  45,  23,  153, 3,   37,  79, 
	81,  123, 157, 7,   137, 115, 149, 63, 
	65,  235, 13,  247, 121, 227, 5,   47, 
	49,  91,  125, 231, 105, 83,  117, 31, 
	33,  203, 237, 215, 89,  195, 229, 15, 
	17,  59,  93,  199, 73,  51,  85,  255 
};

#ifdef MPQS_STAT
u32_t stat_mpqs_nsieves, stat_mpqs_nsurvivors, stat_mpqs_ntrials, stat_mpqs_ndiv;
#endif

static u32_t mpqs_inv_32(u32_t a)
{
	u32_t inv, h;

	inv = (u32_t)mpqs_256_inv_table[(a & 0xff) >> 1];

	h = a * inv; 
	h &= 0xff00;
	h *= inv; 
	inv -= h;

	h = a * inv; 
	h &= 0xffff0000U;
	h *= inv;  
	inv -= h;

#ifdef _DEBUG
	if (inv * a != 1ULL) 
		complain("mpqs_inv_32() failed!");
#endif
	return inv;
}

static u64_t mpqs_inv_64(u64_t a)
{
  u64_t inv, h;

  inv = (u64_t)mpqs_256_inv_table[(a & 0xff) >> 1];
  
  h = a*inv; 
  h &= 0xff00;
  h *= inv; 
  inv -= h;
  
  h = a*inv; 
  h &= 0xffff0000U;
  h *= inv;  
  inv -= h;

  h = a*inv; 
  h &= 0xffffffff00000000ULL;
  h *= inv; 
  inv -= h;

#ifdef _DEBUG
  if (inv * a != 1ULL) 
	  complain("mpqs_inv_64 failed!");
#endif
  return inv;
}


static unsigned char mpqs_jacobi_tab2[8]={0,0,0,1,0,1,0,0};

static int mpqs_jacobi(u16_t a, u16_t b) /* always gcd(a,b)=1 and 0<a<b */
{
  int e, l, m, n, r;

  assert(0 < a);
  assert(a < b);

  m = (int)a; 
  n = (int)b;

  if (!(n & 1)) 
	  complain("mpqs_jacobi: b is even!\n");

  e = 0;
  
  if (!(m & 1)) 
  {
    l = mpqs_jacobi_tab2[n & 7];
    
	do 
	{
      m /= 2; 
	  e += l;
    } while (!(m & 1));
  }

  while (m > 1) 
  {
    r = n - m;

    if (r > 0) 
	{
      l = m & n;
      l >>= 1;
      e += l;
      n = m; 
	  m = r;
    } 
	else 
	{
      m =- r;
    }

    if (m & 1) 
		continue;
    
	l = mpqs_jacobi_tab2[n & 7];

    do 
	{
      m /= 2; 
	  e += l;
    } while (!(m & 1));
  }

  if (!m) 
	  complain("mpqs_jacobi(): something odd.\n");
 
  return 1 - 2*(e & 1);
}

/*
 *   return such 'e' that e*a = 1 mod p
 */
static u16_t mpqs_invert(u16_t a, u16_t p)
{
  u32_t v1 = 0;
  u32_t v2 = 1;
  u32_t b  = (u32_t)a;
  u32_t pp = (u32_t)p;
  u32_t q;

  while (b > 1) 
  {
    pp -= b; 
	v1 += v2;
    
	if (pp >= b) 
	{
      pp -= b; v1 += v2;
      
	  if (pp >= b) 
	  {
        q = pp / b; 
		pp %= b;
        v1 += q * v2;
      }
    }
    
	if (pp <= 1) 
	{ 
		v2 = p - v1; 
		break; 
	}

    b -= pp; 
	v2 += v1;
    
	if (b >= pp) 
	{
      b -= pp; 
	  v2 += v1;
      
	  if (b >= pp) 
	  {
        q = b / pp; 
		b %= pp;
        v2 += q * v1;
      }
    }
  }
/* BUGBUG: this check fails in the debug build! Something seems to be wrong. */
/*
#ifdef _DEBUG
  if ((v2*a - 1) % (u32_t)p) 
  {
	  assert(0);
	  complain("mpqs_invert(): something odd.");
  }
#endif
*/
  return v2;
}

/*
 *  Calculate a^e mod p
 */
static u16_t mpqs_powmod(u16_t a, u16_t e, u16_t p)
{
  u16_t ex = e;
  u32_t aa, res;

  if (!ex) 
	  return 1;

  aa = (u32_t)a;
  res = 1;
  
  while (1, 1) 
  {
    if (ex & 1) 
	{ 
		res *= aa; 
		res %= (u32_t)p; 
	}

    ex >>= 1;
    
	if (!ex) 
		return (u16_t)res;

    aa *= aa; 
	aa %= (u32_t)p;
  }

  assert(0); /* never reached */
  return 0;
}


static u16_t mpqs_sqrt_init(u16_t p)
{
  u16_t e, u, i, j;
  u32_t b;
  u32_t g = 0;

  if (p & 2) 
	  return 0;

  if (p & 4) 
	  return mpqs_powmod(2, (p - 1) >> 2, p);

  e = 0; 
  u = p - 1;

  while (!(u & 1)) 
  { 
	  e++; 
	  u >>= 1; 
  }

  for (i = 2; i < p; i++) 
  {
    g = mpqs_powmod(i, u, p); /* g = i^u mod p */

    if (g == 1) 
		continue;

    b = g;
    	
	for (j = 0; j < (e - 1); j++) 
	{ 
		b *= b; 
		b %= (u32_t)p; 
	}
    
	if (b == (u32_t)(p - 1)) 
		break;
  }

#ifdef _DEBUG
  if (i >= p) 
	  complain("mpqs_sqrt_init() something odd.\n");
#endif

  return g;
}

static u16_t mpqs_sqrt(u16_t a, u16_t p, u16_t help)
{
  u16_t e, u, i, l;
  u32_t b, g, r, k;

  assert(a != 0);
  assert(a < p);
  assert(p != 2);
  /* assert(isprime(p)); */
  assert(mpqs_jacobi(a, p) == 1);

  if (p & 2) 
	  return mpqs_powmod(a, (p + 1) >> 2, p);
  
  if (p & 4) 
  {
    b = mpqs_powmod(a, (p + 3) >> 3, p);
    g = b * b; g %= (u32_t)p;
    
	if (g == a) 
		return b;
    
	b *= (u32_t)help; 
	b %= (u32_t)p;
    
	return (u16_t)b;
  }

  e = 0; 
  u = p - 1;
  
  while (!(u & 1)) 
  { 
	  e++; 
	  u >>= 1; 
  }

  g = help;

#if 1
  r = mpqs_powmod(a, u >> 1, p);
  k = r*r; 
  k %= (u32_t)p;
  k *= a; 
  k %= (u32_t)p;
  r *= a; 
  r %= (u32_t)p;
#else
  r = mpqs_powmod(a, (u + 1)/2, p);
  k = mpqs_powmod(a, u, p);
#endif

  while (k != 1) 
  {
    b = k;
    
	for (i = 0; i < e && b != 1; i++) 
	{ 
		b *= b; 
		b %= (u32_t)p; 
	}
    
	for (l = i + 1; l < e; l++) 
	{ 
		g *= g; 
		g %= (u32_t)p; 
	}
    
	r *= g; 
	r %= (u32_t)p;
    g *= g; 
	g %= (u32_t)p;
    k *= g; 
	k %= (u32_t)p;
    e = i;
  }

  return (u16_t)r;
}

/* ---------------------------------------------------------- */

void mpqs_init()
{
  u32_t i, j;
  u32_t p, add, d;

  mpz_init(mpqs_N);
  mpz_init(mpqs_kN);
  mpz_init(mpqs_dummy);
  mpz_init(mpqs_gmp_help);
  mpz_init(mpqs_help);
  mpz_init(mpqs_sq1);
  mpz_init(mpqs_sq2);
  
  for (i = 0; i < MPQS_MAX_FACTORS; i++) 
	  mpz_init(mpqs_factors[i]);

  mpqs_sievearray = (unsigned char *)xmalloc(MPQS_SIEVELEN*sizeof(unsigned char));
  mpqs_relations = (u16_t **)xmalloc(MPQS_MAX_NRELS*sizeof(u16_t *));
  mpqs_relations[0] = (u16_t *)xmalloc(MPQS_MAX_NRELS*MPQS_REL_ENTRIES*sizeof(u16_t));
  
  for (i = 1; i < MPQS_MAX_NRELS; i++)
    mpqs_relations[i] = mpqs_relations[i - 1] + MPQS_REL_ENTRIES;

  mpqs_part_rels = (u16_t **)xmalloc(MPQS_MAX_NPRELS*sizeof(u16_t *));
  mpqs_part_rels[0] = (u16_t *)xmalloc(MPQS_MAX_NPRELS*MPQS_REL_ENTRIES*sizeof(u16_t));
  
  for (i = 1; i < MPQS_MAX_NPRELS; i++)
    mpqs_part_rels[i]=mpqs_part_rels[i-1] + MPQS_REL_ENTRIES;
  
  mpqs_comb_rels = (u16_t **)xmalloc(MPQS_MAX_NCRELS*sizeof(u16_t *));
  mpqs_comb_rels[0] = (u16_t *)xmalloc(MPQS_MAX_NCRELS*2*MPQS_REL_ENTRIES*sizeof(u16_t));
  
  for (i = 1; i < MPQS_MAX_NCRELS; i++)
    mpqs_comb_rels[i] = mpqs_comb_rels[i-1] + 2*MPQS_REL_ENTRIES;
  
  mpqs_rel_buffer = (u16_t **)xmalloc(MPQS_MAX_NRELS_BUF*sizeof(u16_t *));
  mpqs_rel_buffer[0] = (u16_t *)xmalloc(MPQS_MAX_NRELS_BUF*MPQS_REL_ENTRIES*sizeof(u16_t));
  
  for (i = 1; i < MPQS_MAX_NRELS_BUF; i++)
    mpqs_rel_buffer[i] = mpqs_rel_buffer[i-1] + MPQS_REL_ENTRIES;

  mpqs_gauss_row=(u32_t **)xmalloc(MPQS_GAUSS_MAX*sizeof(u32_t *));
  mpqs_gauss_row[0]=(u32_t *)xmalloc(MPQS_GAUSS_MAX*MPQS_GAUSS_MAX/32*sizeof(u32_t));

#if 0
  prime_table_init();
  for (i = 0; i < MPQS_MAX_NPRIMES; i++)
  {
    mpqs_prime_table[i] = (u16_t)get_next_prime();
    mpqs_prime_sqrthelp[i] = mpqs_sqrt_init(mpqs_prime_table[i]);
  }
#else
  mpqs_prime_table[0] = 2;
  mpqs_prime_table[1] = 3;
  mpqs_prime_table[2] = 5; 
  mpqs_prime_sqrthelp[2] = mpqs_sqrt_init(5);
  i = 3; 
  p = 7; 
  add = 4;
  
  /*
   *  Perform Eratosthenes sieve
   */
  while (i < MPQS_MAX_NPRIMES) 
  {
    for (j = 2; j < i; j++) 
	{
      d = (u32_t)mpqs_prime_table[j];
      
	  if (d*d >= p) 
		  break;
      
	  if (p%d == 0) 
		  break;
    }

    if (d*d > p) 
	{
      mpqs_prime_table[i] = (u16_t)p;
      mpqs_prime_sqrthelp[i] = mpqs_sqrt_init(p);
      i++;
    }

    p += add; 
	add = 6 - add;
  }
#endif

  mpqs_isinit = 1;
}


static u16_t mpqs_multiplier_cand[4][MPQS_MULT_NCAND]=
{
 { 1, 17, 33, 41, 57, 65, 73, 89 },
 { 3, 11, 19, 35, 43, 51, 59, 67 },
 { 5, 13, 21, 29, 37, 53, 61, 69 },
 { 7, 15, 23, 31, 39, 47, 55, 71 },
};

/*
 *  STEN: Choose multiplifier k such that kN = 1 mod 8. We also want to choose 
 *        such k, that gives us a factor base rich with small primes.
 */
static void mpqs_choose_multiplier()
{
  u16_t Nmod8, mult, p, mm;
  u16_t residue[MPQS_MULT_NTESTPR];
  double value[MPQS_MULT_NCAND], v, vp, vmin, dn;
  u32_t i, j;  

  /*
   *  Calc residue_i = N mod p_i
   */
  for (i = 0; i < MPQS_MULT_NTESTPR; i++)
    residue[i] = (u16_t)mpz_mod_ui(mpqs_dummy, mpqs_N, (u32_t)mpqs_prime_table[i]);
  
  Nmod8 = (u16_t)mpz_mod_ui(mpqs_dummy, mpqs_N, 8); /* Nmod8 = N mod 8      */
  dn = log(mpz_get_d(mpqs_N)) / log(2.);            /* dn = log2(N)         */
  dn = 100.9 - dn;                                  /* dn = 100.9 - log2(N) */
  
  if (dn > 7) 
	  mm = 128;
  else 
	  mm = (u16_t)(exp(dn*log(2)));                 /* mm = e^(100.9 - log2(N)) */
  
  /*
   * STEN:
   *
   * We now pick correct multiplifier candidate row and for all multiplifier
   * candidates in it do the following:
   *
   *   Calculate some characteristic 'v' of the multiplifier and store it to the
   *   value[j] array. We'll pick up later the multiplifier with the lowest 
   *   characteristic of all.
   *
   * The characteristic calculation process is done as following 
   * (some variation of the Knuth-Schroeppel function ??):
   * 
   *   For the first MPQS_MULT_NCAND primes p_i we calculate the
   *   jacobian function of the (k * N) mod p_i, where k - is our candidate
   *   multiplifier. If jacobian(k*N mod p_i, p_i) == 1 we substract 
   *   log(p_i)/p_i from the characteristic otherwise we add this value thus 
   *   making the characteristic somewhat less acceptable.
   *
   */
  for (j = 0; j < MPQS_MULT_NCAND; j++) 
  {
	/* we prefer kN = 1 mod 8 since 2 \in FB only under this condition */
    mult = mpqs_multiplier_cand[Nmod8/2][j];
   
	if (mult <= mm) 
	{
      v = log((double)mult) - 1.5*log(2);
      
	  for (i = 1; i < MPQS_MULT_NTESTPR; i++) 
	  {
        p = mpqs_prime_table[i]; 
		vp = log((double)p)/((double)p); /* vp = log(p)/p */
      
		if (mult % p) 
		{
          if (mpqs_jacobi((mult * residue[i])%p, p) == 1) 
			  v -= vp;
          else 
			  v += vp;
        }
      }
    } 
	else 
		v = 1000;
    
	value[j] = v;
  }

  vmin = value[0];
  
  for (j = 1; j < MPQS_MULT_NCAND; j++)
  {
    if (value[j] < vmin) 
		vmin = value[j];
  }

  for (j = 0; j < MPQS_MULT_NCAND; j++)
  {
    if (value[j]==vmin) 
		break;
  }

  if (j >= MPQS_MULT_NCAND) 
	  complain("mpqs_choose_multiplier.vmin\n");

  mpqs_complexity = (vmin + log(mpz_get_d(mpqs_N))) / log(2.);
  mpqs_multiplier = mpqs_multiplier_cand[Nmod8/2][j];

#ifdef MPQS_STAT
printf("%ld ",mpqs_multiplier);
#endif

  mpz_mul_ui(mpqs_kN, mpqs_N, mpqs_multiplier); /* kN = k*N */
  mpqs_kN_64 = mpz_get_ull(mpqs_kN);
}

/*
 * STEN: Array of MPQS parameters. The correct row is choosed depending 
 *       on mpqs_complexity value.
 *
 *       {nFB, sievebegin, nAdiv, nAdiv_total, accept, td_begin}
 *
 */
static u16_t mpqs_param[12][6] =
{
 { 40, 3, 2, 4, 11, 16},    
 { 50, 3, 2, 4, 12, 16},
 { 60, 3, 3, 4, 15, 16},
 { 70, 3, 3, 5, 14, 16},
 { 80, 3, 3, 5, 14, 16},
 { 90, 3, 3, 5, 15, 20},
 { 110, 3, 3, 5, 17, 20},
 { 120, 3, 3, 5, 19, 20},
 { 140, 3, 4, 6, 18, 30},
 { 140, 3, 4, 6, 20, 40},
 { 160, 3, 4, 6, 21, 50},
 { 180, 4, 4, 6, 23, 70}
};

/* 
 * STEN: Choose mpqs parameters.
 *
 * input:  mpqs_complexity
 * output: mpqs_nFB, mpqs_sievebegin, mpqs_nAdiv, mpqs_nAdiv_total, 
 *         mpqs_accept, mpqs_td_begin, 
 *         and 
 *         mpqs_sievelen, mpqs_disp, mpqs_nfactors, mpqs_nrels,
 *         mpqs_nlp, mpqs_nprels, mpqs_ncrels.
 *
 */
static void mpqs_choose_parameter(u16_t retry)
{
  u16_t n, r;

  n = (u16_t)(mpqs_complexity + 0.5); /* n = round(mpqs_complexity) */
  if (n > 96) n = 96;
  if (n < 51) n = 51;
  n = (n + 3) / 4;
  n -= 13;                           /* 0 <= n <= 11 */

  r = retry;
  
  if ((r) && (n < 9)) 
  { 
	  n++; 
	  r--; 
  }

  mpqs_nFB = mpqs_param[n][0];

  if (r) 
	  mpqs_nFB += r * mpqs_nFB / 4;

  if (mpqs_nFB >= MPQS_MAX_FBSIZE) 
	  mpqs_nFB = MPQS_MAX_FBSIZE - 1;

  mpqs_sievebegin = mpqs_param[n][1];
  
  if (r) 
	  mpqs_sievebegin += r;

  mpqs_nAdiv       = mpqs_param[n][2];
  mpqs_nAdiv_total = mpqs_param[n][3];
  mpqs_accept      = mpqs_param[n][4];
  mpqs_td_begin    = mpqs_param[n][5];

  mpqs_sievelen    = MPQS_SIEVELEN;  /* !!! */
  mpqs_disp        = mpqs_sievelen / 2;
  mpqs_nfactors    = 0; 
  mpqs_nrels       = 0; 
  mpqs_nlp         = 0; 
  mpqs_nprels      = 0; 
  mpqs_ncrels      = 0;

  memset(mpqs_hash_table, 0, sizeof(mpqs_hash_table));
  memset(mpqs_rel_hash_table, 0, sizeof(mpqs_rel_hash_table));
}

/*
 * Initialize and fill with primes mpqs_FB and mpqs_FBk factor bases.
 */
static void mpqs_generate_FB()
{
  u16_t *fb, p, rest, i, nfb;

  i   = 0; 
  nfb = 0; 
  fb  = mpqs_FB; 
  mpqs_nFBk = 0;
  
  /* 2 is always in factor base since we have choosen such k that kN = 1 mod 8 */
  *fb++=2; 
  *fb++=1; 
  nfb++;
  
  for (i = 1; i < MPQS_MAX_NPRIMES; i++) 
  {
    p = mpqs_prime_table[i];
    
	if (p > MPQS_FB_MAXPRIME) 
		break;
    
	rest = (u16_t)mpz_mod_ui(mpqs_dummy, mpqs_kN, p);  /* rest = kN mod p */
    
	if (rest) 
	{
      if (mpqs_jacobi(rest, p) == 1)
	  {
		/* we have found such p that (kN/p) == 1 */
        *fb++ = p;
        *fb++ = mpqs_sqrt(rest, p, mpqs_prime_sqrthelp[i]); /* ca 13M */
        nfb++;
       
		if (nfb >= mpqs_nFB) 
			break;
      }
    } 
	else 
	{
      /* we have found p that divides kN */
      if (mpqs_multiplier%(u32_t)p)
        complain("mpqs: N has small divisor: %u\n",p);
    
	  mpqs_FBk[mpqs_nFBk++] = p;
    }
  }

  assert(2*(nfb - 1) < sizeof(mpqs_FB)/sizeof(mpqs_FB[0]));
  mpqs_nFB       = nfb;
  mpqs_nFBk_1    = mpqs_nFBk + 1;
  mpqs_nsp       = 1 + mpqs_nFB + mpqs_nFBk;
  mpqs_pmax      = mpqs_FB[2*(nfb - 1)];
  mpqs_FB[2*nfb] = mpqs_sievelen;
}

static int mpqs_SI_init()
{
  double d;
  int  a, i, j, k, l, n;
  i64_t prod;
  u16_t inv, p, *fb = mpqs_FB;
  double A_div_log[MPQS_MAX_ADIV_ALL]; /* maximal number of A's */
  double v;

  d = mpz_get_d(mpqs_kN);              /* d = kN                                      */
  mpqs_kN_dbl = d;
  d *= 8.;                             /* d = 8 * kN                                  */ 
  d = sqrt(d);                         /* d = sqrt(8kN)                               */
  d /= (double)mpqs_sievelen;          /* d = sqrt(8kN) / sievelen                    */
  d = log(d);                          /* d = log(sqrt(8kN) / sievelen)               */
  d /= (double)mpqs_nAdiv;             /* d = log(sqrt(8kN) / sievelen) / nAdiv       */
  d = exp(d);                          /* d = e^(log(sqrt(8kN) / sievelen) / nAdiv)   */
  a = (long)d;                         /* a = [e^(log(sqrt(8kN) / sievelen) / nAdiv)] */

#ifdef _DEBUG
  if (a >= mpqs_pmax) 
  { 
	  printf("mpqs_SI_init err: %ld %ld ", a, mpqs_pmax); 
	  return -1; 
  }
#endif
  
  if (a >= mpqs_pmax) 
	  a = mpqs_pmax - 1;
  
  /*
   *  STEN: Search for the first prime in factor base that is greater then 
   *        'a' calculated above. After the loop i contains this prime's index.
   *  
   */
  for (i = 0; i < mpqs_nFB; i++)
  {
	  if (mpqs_FB[2*i] > a) 
		  break;
  }

  i -= mpqs_nAdiv_total/2;
  
  if (i < 1) 
	  i = 1; /* first prime >2 */
  
  if (mpqs_FB[2*i] < 3) 
	  return -2;
  
  if (i + mpqs_nAdiv_total > mpqs_nFB)
  {
	  i = mpqs_nFB - mpqs_nAdiv_total - 1;
  }
	  
  for (j = 0; j < mpqs_nAdiv_total; j++) 
  {
    mpqs_Adiv_all[j] = mpqs_FB[2*(i + j)];
    mpqs_Adiv_all_sqrt[j] = mpqs_FB[2*(i + j) + 1];
  }

  for (j = i + mpqs_nAdiv_total; j < mpqs_nFB; j++) 
  {
    mpqs_FB[2*(j - mpqs_nAdiv_total)] = mpqs_FB[2*j];
    mpqs_FB[2*(j - mpqs_nAdiv_total) + 1] = mpqs_FB[2*j + 1];
  }

  mpqs_nFB -= mpqs_nAdiv_total;
  
  if (i < mpqs_sievebegin) 
	  mpqs_sievebegin = (u16_t)i;

  p = mpqs_FB[2];
  mpqs_FB_inv[1] = mpqs_inv_32(p);
  mpqs_FB_inv_info[2] = p; 
  mpqs_FB_inv_info[3] = p;
  mpqs_FB_inv_info[6] = mpqs_FB_inv[1];
  mpqs_FB_inv_info[7] = mpqs_FB_inv[1];
  
  if (mpqs_td_begin & 1) 
	  complain("mpqs_td_begin is odd\n");
  
  for (j = 2; j < mpqs_td_begin; j += 2) 
  {
    p = mpqs_FB[2*j];
    mpqs_FB_inv[j] = mpqs_inv_32(p);
    mpqs_FB_inv_info[4*j] = p; 
	mpqs_FB_inv_info[4*j + 1] = p;
    mpqs_FB_inv_info[4*j + 4] = mpqs_FB_inv[j];
    mpqs_FB_inv_info[4*j + 5] = mpqs_FB_inv[j];
    
	p = mpqs_FB[2*(j + 1)];
    mpqs_FB_inv[j + 1] = mpqs_inv_32(p);
    mpqs_FB_inv_info[4*j + 2] = p; 
	mpqs_FB_inv_info[4*j+3] = p;
    mpqs_FB_inv_info[4*j + 6] = mpqs_FB_inv[j + 1];
    mpqs_FB_inv_info[4*j + 7] = mpqs_FB_inv[j + 1];
  }

  for (; j < mpqs_nFB; j++) 
	  mpqs_FB_inv[j] = mpqs_inv_32(mpqs_FB[2*j]);

  for (j = 0; j < mpqs_nFBk; j++) 
	  mpqs_FBk_inv[j] = mpqs_inv_32(mpqs_FBk[j]);
  
  for (j = 0; j < mpqs_nAdiv_total; j++)
    mpqs_FB_A_inv[j] = mpqs_inv_32(mpqs_Adiv_all[j]);

  /* compute log-approximations */
  d = mpz_get_d(mpqs_kN);                      /* d = kN                                                      */
  d /= 8.;                                     /* d = kN / 8                                                  */
  d = sqrt(d);                                 /* d = sqrt(kN / 8)                                            */
  d *= (double)mpqs_sievelen;                  /* d = sqrt(kN / 8) * sievelen                                 */
  d = log(d);                                  /* d = log(sqrt(kN / 8) * sievelen)                            */
  d -= log((double)((u64_t)1 << mpqs_accept)); /* d = log(sqrt(kN / 8) * sievelen / (1 << mpqs_accept))       */
  d = 128. / d;                                /* d = 128 / log(sqrt(kN / 8) * sievelen / (1 << mpqs_accept)) */

  /* precomute mpqs_FB_log[] */
  for (i = 0; i < mpqs_nFB; i++) 
  {
    p = mpqs_FB[2*i];
    mpqs_FB_log[i] = (unsigned char)(0.5 + d * log((double)p));
  }

  /* precomute mpqs_Adiv_log[] */
  for (i = 0; i < mpqs_nAdiv_total; i++) 
  {
    p = mpqs_Adiv_all[i];
    mpqs_Adiv_log[i] = (unsigned char)(0.5 + d*log((double)p));
  }

  /* compute mask for choice of A-divisors */
  n = 1 << mpqs_nAdiv_total;
  
  for (i = 0; i < n; i++)
  {
	  mpqs_A_mask[i] = 0;
  }

  for (i = 0; i < mpqs_nAdiv_total; i++) 
  {
    k = 1 << i;
    
	for (j = k; j < n; j += (2*k))
	{
      for (l = 0; l < k; l++)
	  {
		  mpqs_A_mask[j + l]++;
	  }
	}
  } 
  
  /* now mpqs_A_mask[i] contains the number of ones of i (binary) */
  j = 0;
  for (i = 0; i < n; i++)
  {
    if (mpqs_A_mask[i] == mpqs_nAdiv) 
	{
	  assert(i <= UCHAR_MAX);
      mpqs_A_mask[j] = (unsigned char)i; 
	  j++;
    }
  }

  for (i = 0; i < mpqs_nAdiv_total; i++)
  {
    A_div_log[i] = log((double)(mpqs_Adiv_all[i])) / log(2.);
  }

  d = log(mpqs_kN_dbl) / log(2.); /* d ~ number_of_bits(kN) */
  
  for (i = 0; i < j; i++) 
  {
    v = 0;
    
	for (k = 0; k < mpqs_nAdiv_total; k++)
	{
      if (mpqs_A_mask[i] & (1 << k)) 
		  v += A_div_log[k];
	}

    if (d - v > 63.95) 
	{
      j--; /* printf(":");*/
      mpqs_A_mask[i] = mpqs_A_mask[j];
      i--;
    }
  }

  assert(j <= USHRT_MAX);
  mpqs_nA = (u16_t)j; 
  mpqs_A_index = -1; 
  mpqs_B_index = (1 << mpqs_nAdiv) - 1;

  prod = 1; 
  d = 1.;
  
  for (i = 0; i < mpqs_nAdiv_total; i++) 
  {
    prod *= (i64_t)(mpqs_Adiv_all[i]);
    d *= (double)(mpqs_Adiv_all[i]);
  }
  
  /* ensure that |C| <= 2^64 */
  if (d > (double)LLONG_MAX) 
	  return -4;
  
  for (i = 1; i < mpqs_nFB; i++) 
  {
    p = fb[2*i];
    inv = (u16_t)(prod % (i64_t)p);
    inv = mpqs_invert(inv, p);
    inv <<= 1; 
	
	if (inv >= p) 
		inv -= p;
    
	mpqs_2A_all_inv[i] = inv;
  }

  mpqs_FB[2*mpqs_nFB] = mpqs_sievelen;
  return 1;
}


static int mpqs_next_pol()
{
  u32_t i, ind, j;
  i64_t add, bi, prod_other;
  u16_t p, s1, s2, inv, *fb=mpqs_FB, *fbs=mpqs_FB_start;
  i16_t sh;
  u32_t bb, cc, c1, c2;
  unsigned char mask;
  u16_t aa;

#define GCC_BUG
#ifndef GCC_BUG
  u32_t al, *ul;
#endif

  mpqs_B_index++;

  if (mpqs_B_index >= (1 << (mpqs_nAdiv-1))) 
  {
    mpqs_A_index++;
    
	if (mpqs_A_index >= mpqs_nA) 
		return 0;

    mask = mpqs_A_mask[mpqs_A_index];
    
	for (i = 0; i < mpqs_nAdiv_total; i++)
	{
      if (mask & (1 << i)) 
		  mpqs_Adiv_active[i] = 1;
      else 
		  mpqs_Adiv_active[i] = 0;
	}

    j = 0;
    
	for (i=0; i<mpqs_nAdiv_total; i++)
	{
      if (mpqs_Adiv_active[i]) 
	  {
        mpqs_Adiv[j] = mpqs_Adiv_all[i];
        mpqs_Adiv_sqrt[j] = mpqs_Adiv_all_sqrt[i];
        j++;
      }
	}

    mpqs_A = 1;
    
	for (i = 0; i < mpqs_nAdiv; i++) 
		mpqs_A *= (i64_t)(mpqs_Adiv[i]);

	/* compute B_i's */
    for (i = 0; i < mpqs_nAdiv; i++) 
	{
      p   = mpqs_Adiv[i];
      bi  = mpqs_A / (i64_t)(p);
      inv = (u16_t)(bi % (i64_t)p); /* bi>0 */
      inv = mpqs_invert(inv,p);
      bi *= (i64_t)inv; 
	  bi %= mpqs_A;
      bi *= (i64_t)mpqs_Adiv_sqrt[i]; 
	  bi %= mpqs_A;
      mpqs_Bi[i] = bi;
    }

    mpqs_B = 0;
    
	for (i = 0; i < mpqs_nAdiv; i++)
	{
		mpqs_B += mpqs_Bi[i];
	}

    prod_other = 1;
    
	for (i=0; i<mpqs_nAdiv_total; i++)
	{
      if (!mpqs_Adiv_active[i]) 
		  prod_other *= (i64_t)(mpqs_Adiv_all[i]);
	}

    for (i = 1; i < mpqs_nFB; i++) 
	{
      p = fb[2*i];
	  assert(p);

	  bb = (u32_t)(((i64_t)(mpqs_2A_all_inv[i]) * prod_other) % (u32_t)p);

      for (j = 0; j < mpqs_nAdiv; j++) 
	  {
        cc = bb * (u32_t)(mpqs_Bi[j] % (i64_t)p); 
		cc %= (u32_t)p;
        mpqs_SI_add[j][i] = (u16_t)cc;
      }

      if (bb & 1) 
		  bb += p;
      
	  bb >>= 1;  /* now bb=1/A mod p */
      cc = p - (u32_t)(mpqs_B % (i64_t)p);  /* mpqs_B > 0 */
      c1 = cc + fb[2*i+1]; 
	  c1 *= bb; 
	  c1 += mpqs_disp; 
	  c1 %= (u32_t)p;
      c2 = cc + (p - fb[2*i+1]); 
	  c2 *= bb; 
	  c2 += mpqs_disp;
	  c2 %= (u32_t)p;
      
	  fbs[2*i]  =(u16_t)c1; 
	  fbs[2*i+1]=(u16_t)c2;
    }

    mpqs_2B = 2 * mpqs_B;
    mpqs_A_inv_64 = mpqs_inv_64(mpqs_A);
    mpqs_C = (mpqs_B * mpqs_B - mpqs_kN_64) * mpqs_A_inv_64; /* C = (B^2 - kN) / A */

	// ***
    for (i=0; i<mpqs_nAdiv_total; i++)
	{
      if (!mpqs_Adiv_active[i]) 
	  {
        p = mpqs_Adiv_all[i];
        bb = (u32_t)(mpqs_A % (i64_t)p);
        bb = mpqs_invert((u16_t)bb, p);
        bb <<= 1;

		if (bb >= p) 
			bb -= p;

        for (j = 0; j < mpqs_nAdiv; j++) 
		{
          cc = bb * (u32_t)(mpqs_Bi[j] % (i64_t)p); 
		  cc %= (u32_t)p;
          mpqs_Adiv_SI_add[j][i] = (u16_t)cc;
        }

        if (bb & 1) 
			bb += p;

        bb >>= 1;  /* now bb = 1/A mod p */

        cc = p - (u32_t)(mpqs_B % (i64_t)p);  /* mpqs_B > 0 */
        c1 = cc + mpqs_Adiv_all_sqrt[i]; 
		c1 *= bb; 
		c1 += mpqs_disp; 
		c1 %= (u32_t)p;
        c2 = cc + (p - mpqs_Adiv_all_sqrt[i]); 
		c2 *= bb; 
		c2 += mpqs_disp; 
		c2 %= (u32_t)p;
        
		if (c1 < c2) 
		{
          mpqs_Adiv_start1[i] = (u16_t)c1; 
		  mpqs_Adiv_start2[i] = (u16_t)c2;
        } 
		else 
		{
          mpqs_Adiv_start1[i] = (u16_t)c2; 
		  mpqs_Adiv_start2[i] = (u16_t)c1;
        }
      } 
	  else 
	  {
        p = mpqs_Adiv_all[i];
        sh = (i16_t)(mpqs_2B % (i64_t)p);
        
		if (sh > 0) 
			bb = (u32_t)sh; 
		else 
			bb = (u32_t)((i16_t)p + sh);
        
		bb = (u32_t)mpqs_invert(bb, p);
        sh = (i16_t)(mpqs_C % (i64_t)p); 
		sh = -sh;
        
		if (sh > 0) 
			cc = (u32_t)sh; 
		else 
			cc = (u32_t)((i16_t)p + sh);

        cc *= bb; 
		cc %= (u32_t)p;  /* now p|(2B * cc + C) */
        cc += mpqs_disp; 
		cc %= (u32_t)p;
        mpqs_Adiv_start1[i] = (u16_t)cc;
      }
	}

    mpqs_B_index = 0;
    return 1;
  }

  ind = mpqs_B_index; 
  i = 0;
  
  while (1) 
  {
    if (ind & 1) 
		break;
    
	i++; 
	ind >>= 1;
  }

  ind >>= 1; 
  ind++;
  add = 2 * mpqs_Bi[i];
  
  if (ind & 1) 
  {
    mpqs_B -= add;
    
	for (j = 1; j < mpqs_nFB; j++) 
	{
      p = fb[2*j]; 
/* 
	  s1 = fbs[2*j];   
	  s2 = fbs[2*j + 1]; 
*/
      aa = mpqs_SI_add[i][j];
#ifndef GCC_BUG
      al = (((u32_t)aa) << 16) + (u32_t)aa;
      ul = (u32_t *)(fbs + 2*j);
      *ul += al;
#else
      fbs[2*j] += aa; 
	  fbs[2*j + 1] += aa;
#endif
      if (fbs[2*j] >= p) 
		  fbs[2*j] -= p;
      
	  if (fbs[2*j + 1] >= p) 
		  fbs[2*j + 1] -=p;
    }
  } 
  else 
  {
    mpqs_B += add;
    for (j = 1; j < mpqs_nFB; j++) 
	{
      p = fb[2*j]; /*s1=fbs[2*j]; s2=fbs[2*j+1];*/
      aa = p - mpqs_SI_add[i][j];
#ifndef GCC_BUG
      al = (((u32_t)aa) << 16) + (u32_t)aa;
      ul = (u32_t *)(fbs + 2*j);
      *ul += al;
#else
      fbs[2*j] += aa; 
	  fbs[2*j+1] += aa;
#endif
      if (fbs[2*j] >= p) 
		  fbs[2*j] -= p;
      
	  if (fbs[2*j+1] >= p) 
		  fbs[2*j+1] -= p;
    }
  }

  mpqs_2B = 2 * mpqs_B;
  mpqs_C = (mpqs_B * mpqs_B - mpqs_kN_64) * mpqs_A_inv_64; /* C = (B^2 - kN) / A */

  /* printf("mpqs_next_pol() pol: %Ld %Ld %Ld \n", mpqs_A, mpqs_B, mpqs_C); */
  
  for (j=0; j<mpqs_nAdiv_total; j++)
  {
    if (!mpqs_Adiv_active[j]) 
	{
      p  = mpqs_Adiv_all[j];
      s1 = mpqs_Adiv_start1[j]; 
	  s2 = mpqs_Adiv_start2[j];
      
	  if (ind & 1) 
	  {
        s1 += mpqs_Adiv_SI_add[i][j]; 	
		s2 += mpqs_Adiv_SI_add[i][j]; 
      } 
	  else 
	  {
        s1 += (p - mpqs_Adiv_SI_add[i][j]);        
		s2 += (p - mpqs_Adiv_SI_add[i][j]); 
      }

	  if (s1 >= p) s1 -= p;
	  if (s2 >= p) s2 -= p;

      mpqs_Adiv_start1[j] = s1; 
	  mpqs_Adiv_start2[j] = s2;
    } 
	else 
	{
      p = mpqs_Adiv_all[i];
      sh = (i16_t)(mpqs_2B % (i64_t)p);
      
	  if (sh > 0) 
		  bb = (u32_t)sh; 
	  else 
		  bb = (u32_t)((i16_t)p + sh);
      
	  bb = (u32_t)mpqs_invert((u16_t)bb, p);
      sh = (i16_t)(mpqs_C % (i64_t)p); 
	  sh = -sh;
      
	  if (sh > 0) 
		  cc = (u32_t)sh; 
	  else 
		  cc = (u32_t)((i16_t)p + sh);
      
	  cc *= bb; 
	  cc %= (u32_t)p;  /* now p|(2B * cc + C) */      
	  cc += mpqs_disp; 
	  cc %= (u32_t)p;
      mpqs_Adiv_start1[i] = (u16_t)cc;
    }
  }

  return 1;
}

static void mpqs_sieve()
{
  unsigned char *sv=mpqs_sievearray, *fblog=mpqs_FB_log, *svend;
  u16_t *fb=mpqs_FB, *fbs=mpqs_FB_start;
  u16_t p, s1, s2, pb;
  u32_t i;
  unsigned char lo;
  u32_t *ulsv, *ulsvend, mask;

  ulsv = (u32_t *)mpqs_sievearray;
  ulsvend = ulsv + mpqs_sievelen/4;
  
  if (mpqs_B & 1) 
  {
    mask = (u32_t)(mpqs_FB_log[0]); 
	mask *= 0x00040004;
  } 
  else 
  {
    mask = (u32_t)(mpqs_FB_log[0]); 
	mask *= 0x04000400;
  }

  while (ulsv < ulsvend) 
  {
    *ulsv++ = mask;
    *ulsv++ = mask;
    *ulsv++ = mask;
    *ulsv++ = mask;
  }

  pb = mpqs_sievelen >> 2;
  i = mpqs_sievebegin; 
  fb += 2*i; 
  fbs += 2*i;

#ifdef ASM_MPQS_SIEVER
  asm_sieve();
#else
  for (;;) 
  {
    p = *fb; 
	
	if (p > pb)
		break;
    
	lo = fblog[i]; 
	sv = mpqs_sievearray;
    fb += 2;
    s1 = *fbs++; 
	s2 = *fbs++; 
    svend = sv + mpqs_sievelen - 4*p;
    
	while (sv < svend) 
	{
      sv[s1] += lo; sv[s2] += lo; sv += p;
      sv[s1] += lo; sv[s2] += lo; sv += p;
      sv[s1] += lo; sv[s2] += lo; sv += p;
      sv[s1] += lo; sv[s2] += lo; sv += p;
    }

    svend += (p + p);
    
	if (sv < svend) 
	{
      sv[s1] += lo; sv[s2] += lo; sv += p;
      sv[s1] += lo; sv[s2] += lo; sv += p;
    }
    
	svend += p;
    if (sv < svend) 
	{
      sv[s1] += lo; sv[s2] += lo; sv += p;
    }

    svend += p;
    
	if ((sv + s1) < svend) sv[s1] += lo;
    if ((sv + s2) < svend) sv[s2] += lo;
    
	i++;
  }
#endif
  sv = mpqs_sievearray;
  
  for (i = 0; i < mpqs_nAdiv_total; i++)
  {
    if (!mpqs_Adiv_active[i]) 
	{
      p  = mpqs_Adiv_all[i]; 
	  lo = mpqs_Adiv_log[i];
      s1 = mpqs_Adiv_start1[i];
      s2 = mpqs_Adiv_start2[i];
      
	  while (s1 < mpqs_sievelen) 
	  { 
		  sv[s1] += lo; 
		  s1 += p; 
	  }
      
	  while (s2 < mpqs_sievelen) 
	  { 
		  sv[s2] += lo; 
		  s2 += p; 
	  }
    } 
	else 
	{
      p  = mpqs_Adiv_all[i]; 
	  lo = mpqs_Adiv_log[i];
      s1 = mpqs_Adiv_start1[i];
      while (s1 < mpqs_sievelen) 
	  { 
		  sv[s1] += lo; 
		  s1 += p; 
	  }
    }
  }
}

static u16_t mpqs_evaluate()
{
  u32_t *ulsv, *ulsvend, h;
  unsigned char *sv;
  u16_t **rels, buffer[256];
  u32_t i, nmaxsurv;

  mpqs_nsurvivors = 0;
  nmaxsurv = MPQS_MAX_NRELS - mpqs_nrels;
  
  if ((u32_t)(MPQS_MAX_NPRELS - mpqs_nprels) < nmaxsurv)
    nmaxsurv = MPQS_MAX_NPRELS - mpqs_nprels;
  
  if ((u32_t)(MPQS_MAX_NCRELS - mpqs_ncrels) < nmaxsurv)
    nmaxsurv = MPQS_MAX_NCRELS - mpqs_ncrels;
  
  if (nmaxsurv < MPQS_MIN_RELBUFFER) 
	  return 0;
  
  if (nmaxsurv > 255) 
	  nmaxsurv = 255;  /* 1 weniger als buffer-Länge !!! */
  
  ulsv = (u32_t *)mpqs_sievearray;
  ulsvend = ulsv + (mpqs_sievelen >> 2);  
  rels = mpqs_rel_buffer;

#ifndef ASM_MPQS_EVAL
  while (ulsv < ulsvend) 
  {
    h = *ulsv;
    
	if (h & 0x80808080) 
	{
      sv = (unsigned char *)ulsv;
      
	  for (i = 0; i < 4; i++)
	  {
        if (*sv++ & 0x80)
          buffer[mpqs_nsurvivors++] = (u16_t)(sv - mpqs_sievearray - 1);
	  }

      /* CJM: test added here. */
      if (mpqs_nsurvivors > nmaxsurv - 4) 
	  {
        while (ulsv<ulsvend) *ulsv++ = 0;
        break;
      }
    }

    *ulsv++ = 0;
    h = *ulsv;
    
	if (h & 0x80808080) 
	{
      sv = (unsigned char *)ulsv;
      
	  for (i = 0; i < 4; i++)
	  {
        if (*sv++ & 0x80)
          buffer[mpqs_nsurvivors++] = (u16_t)(sv - mpqs_sievearray - 1);
	  }

      if (mpqs_nsurvivors > nmaxsurv - 4) 
	  {
        while (ulsv < ulsvend) 
			*ulsv++ = 0;

        break;
      }
    }

    *ulsv++ = 0;
  }
#else
  mpqs_nsurvivors = asm_evaluate(ulsv, ulsvend, buffer, nmaxsurv); 
  /* printf("%ld ", mpqs_nsurvivors); */

#ifdef MPQS_STAT
  stat_asm_eval++;
#endif
#endif

  for (i = 0; i < mpqs_nsurvivors; i++) 
  {
    mpqs_sievearray[buffer[i]] = (unsigned char)(i + 1);
    rels[i][0] = buffer[i];
  }

  return mpqs_nsurvivors;
}

/* change hash_table size to 256*15 ??? */
static inline u16_t mpqs_hash_lp(u32_t lp)
{
  u32_t j;
  u32_t lp1;
  u16_t lphash, nh;

  lphash = (u16_t)(lp & 0x00ff) >> 1; 
  lp1 = lp >> 8;
  nh = mpqs_hash_table[lphash][0];
  
  if (nh >= 15) 
  { 
	  assert(0);
	  return MPQS_HASH_OVERFLOW; 
  }

  for (j = 1; j <= nh; j++)
  {
    if (mpqs_hash_table[lphash][j] == lp1) 
		break;
  }

  if (j > nh) 
  { /* new lp */
    mpqs_hash_table[lphash][j] = lp1;
    mpqs_hash_table[lphash][0]++;
    mpqs_nlp++;
    mpqs_hash_index[lphash][j] = mpqs_nprels;
    mpqs_lp[mpqs_nprels] = lp;
  }

  return mpqs_hash_index[lphash][j];
}

static inline int mpqs_hash_rel(i64_t axb)
{
  u32_t j;
  u32_t hash, nh;
  u16_t hash2;

  hash = (u32_t)(axb & 0xffffff);
  hash2 = (u16_t)(hash >> 8); 
  hash &= 0xff;
  nh = mpqs_rel_hash_table[hash][0];
  
  if (nh >= 15) 
	  return 1;
  
  for (j = 1; j <= nh; j++)
  {
    if (mpqs_rel_hash_table[hash][j] == hash2) 
		return 1;
  }

  mpqs_rel_hash_table[hash][nh + 1] = hash2;
  mpqs_rel_hash_table[hash][0]++;
  return 0;
}

static int mpqs_decompose()
{
  u16_t i, j;
  u16_t *fb=mpqs_FB, *fbs=mpqs_FB_start;
  unsigned char *sv, *svend, i1, i2;
  u16_t ind, p, s1, s2, nr, lpi, ii;
  u16_t **rels;
  i64_t axb, *llp;
  u64_t qx, pr;
  u16_t minus;
  double dbl_qx, dx;
  i16_t x;
  u32_t lp, ulqx;
  u32_t inv, ls1, ls2;
  u32_t ax, ay, az;

  rels = mpqs_rel_buffer - 1;
  
  for (i = 1; i <= mpqs_nsurvivors; i++) 
	  rels[i][4] = 0;

  i = mpqs_td_begin; 
  fb  += 2*i; 
  fbs += 2*i;
  
  for (; i < mpqs_nFB; i++) 
  {
    p = *fb; 
	fb += 2; 
	sv = mpqs_sievearray; 
	svend = sv + mpqs_sievelen - 2*p;
    s1 = *fbs++; 
	s2 = *fbs++;
    
	while (sv < svend) 
	{
      i1 = sv[s1]; 
	  i2 = sv[s2]; 
	  sv += p;
      
	  if (i1||i2) 
	  {
        if (i1) 
		{ 
			nr=rels[i1][4]; 
			rels[i1][nr + 5] = mpqs_nFBk_1 + i; 
			rels[i1][4]++; 
		}
        
		if (i2) 
		{ 
			nr=rels[i2][4]; 
			rels[i2][nr + 5] = mpqs_nFBk_1 + i; 
			rels[i2][4]++; 
		}      
	  }
      
	  i1 = sv[s1]; 
	  i2 = sv[s2]; 
	  sv += p;
      
	  if (i1 || i2) 
	  {
        if (i1) 
		{ 
			nr = rels[i1][4]; 
			rels[i1][nr + 5] = mpqs_nFBk_1 + i; 
			rels[i1][4]++; 
		}        
		
		if (i2) 
		{ 
			nr = rels[i2][4]; 
			rels[i2][nr + 5] = mpqs_nFBk_1 + i; 
			rels[i2][4]++; 
		}      
	  }
    }

    svend += p;
    
	while (sv < svend) 
	{
      i1 = sv[s1]; 
	  i2 = sv[s2]; 
	  sv += p;
      
	  if (i1 || i2) 
	  {
        if (i1) 
		{ 
			nr = rels[i1][4]; 
			rels[i1][nr + 5] = mpqs_nFBk_1 + i; 
			rels[i1][4]++; 
		}
        
		if (i2) 
		{ 
			nr = rels[i2][4]; 
			rels[i2][nr + 5] = mpqs_nFBk_1 + i; 
			rels[i2][4]++; 
		}
      }
    }

    svend += p;
    
	if ((sv + s1) < svend) 
	{
      i1 = sv[s1];
      
	  if (i1) 
	  { 
		  nr = rels[i1][4]; 
		  rels[i1][nr + 5] = mpqs_nFBk_1 + i; 
		  rels[i1][4]++; 
	  }
    }
    if ((sv + s2) < svend) 
	{
      i2 = sv[s2];
      
	  if (i2) 
	  { 
		  nr = rels[i2][4]; 
		  rels[i2][nr + 5] =mpqs_nFBk_1 + i; 
		  rels[i2][4]++; 
	  }
    }
  }
  
  for (i = 1; i <= mpqs_nsurvivors; i++) 
  {
    nr = rels[i][4];
    ind = rels[i][0];
    x = (i16_t)(ind) - mpqs_disp;

    axb = mpqs_A * (i64_t)x + mpqs_B;
    
	if (axb < 0) 
		axb =- axb;
    
	if (mpqs_hash_rel(axb)) 
		goto next;

    qx = mpqs_A * (u64_t)x + mpqs_2B;
    qx *= (u64_t)x;
    qx += (u64_t)mpqs_C;

/* das folgende besser durch ungefähre Berechnung der Nullstellen
   und Vergleiche ersetzen ?? */
    minus = 0; 
	dx = (double)x;
    dbl_qx = (double)mpqs_A * dx + (double)mpqs_2B;
    dbl_qx *= dx;
    dx = (double)mpqs_C; 
	
	if (dx > 0.) 
		dx -= (double)(ULLONG_MAX) + 1.;

	dbl_qx += dx;

	if (fabs(dbl_qx) >= ((double)(ULLONG_MAX) + 1.)) 
	{
		/*printf(";");*/ 
		goto next; 
	}

    if (dbl_qx < 0.) 
	{ 
		minus = 1; 
		qx = 0 - qx; 
	}

#ifdef ASM_MPQS_TD
#ifdef MPQS_STAT
stat_asm_td++;
#endif
if (asm_td(rels[i],minus,&qx,&ulqx)) {
  goto next;
}
/*
if (nr>=40) {
complain("%Lu %ld, %ld, ",qx,ulqx,nr);
}
*/
/*printf(",");*/
/*if (nr>27) {
  goto next;
}*/
/*printf("%lu ",ulqx);*/
/*      for (j=0; j<mpqs_nFBk; j++) {
        p=mpqs_FBk[j];
        if (ulqx%(u32_t)p==0) {
          complain("%Lu %lu %lu ",qx,ulqx,p);
        }
      }*/
/*      for (j=0; j<mpqs_nFBk; j++) {
        p=mpqs_FBk[j];
        if (ulqx%(u32_t)p==0) {
          complain("%Lu %lu %lu ",qx,ulqx,p);
        }
      }*/

#else
    for (j = 1; j < mpqs_td_begin; j++) 
	{
      p = mpqs_FB[2*j];
      inv = mpqs_FB_inv[j];
      ls1 = ind + (p - mpqs_FB_start[2*j]); 
	  ls2 = ind + (p - mpqs_FB_start[2*j+1]);
      ls1 *= inv; 
	  ls2 *= inv;
      ls1 &= 0xffff0000; 
	  ls2 &= 0xffff0000;
      
	  if (!ls1 || !ls2) 
	  {
          rels[i][5+nr] = mpqs_nFBk_1 + j; 
		  nr++; 
		  /*printf("%d.",p);*/
      }
    }

    pr = 1;
    
	for (j = 0; j < nr; j++) 
		pr *= (u64_t)(mpqs_FB[2*(rels[i][5+j] - mpqs_nFBk_1)]);
    
	rels[i][4] = nr;
    
	if (minus) 
	{
      rels[i][5 + rels[i][4]] = 0; rels[i][4]++;
    }
    
	while (!(qx & 1)) 
	{
      if (rels[i][4] < MPQS_TD_MAX_NDIV) 
	  {
        rels[i][5 + rels[i][4]] = mpqs_nFBk_1; 
		rels[i][4]++;
        qx >>= 1;
      } 
	  else 
	  {
        goto next;
      }
    }
    
	if (pr < ((u64_t)UINT_MAX + 1)) 
	{
      if (qx > (pr << 32)) 
	  {
        /* printf(","); */
        goto next;
      }
    }

	/* STEN: ensure that (u32_t)qx converstion is safe */
	if (qx > ((u64_t)UINT_MAX + 1))
		goto next;

    ax = (u32_t)qx; 
	ay = (u32_t)pr;
	az = ax * mpqs_inv_32(ay);

#ifdef _DEBUG
	if (qx % pr) 
		complain("mpqs_decompose: 'qx %% pr'!");
	
	if (qx/pr != az) 
		complain("mpqs_decompose: 'qx/pr != az'!");
#endif

	ulqx = az;

    for (j = 0; j < nr; j++) 
	{
      ii = rels[i][5 + j];
      p = mpqs_FB[2*(ii - mpqs_nFBk_1)];
      
	  while (ulqx % (u32_t)p == 0) 
	  {
        if (rels[i][4] < MPQS_TD_MAX_NDIV) 
		{
          rels[i][5 + rels[i][4]] = ii; 
		  rels[i][4]++;
          ulqx /= (u32_t)p;
        } 
		else 
		{
          goto next;
        }
      }
    }

    if (ulqx > 1)
	{
      for (j = 0; j < mpqs_nFBk; j++) 
	  {
        p = mpqs_FBk[j];
        while (ulqx % (u32_t)p == 0) 
		{
          if (rels[i][4] < MPQS_TD_MAX_NDIV) 
		  {
            rels[i][5 + rels[i][4]] = 1 + j; 
			rels[i][4]++;
            ulqx /= (u32_t)p;
          } 
		  else 
		  {
            goto next;
          }
        }

        if (ulqx == 1) 
			break;
      }
	}

    if (ulqx > 1)
	{
      for (j = 0; j < mpqs_nAdiv_total; j++) 
	  {
        p = mpqs_Adiv_all[j];
        
		while (ulqx % (u32_t)p == 0) 
		{
          if (rels[i][4] < MPQS_TD_MAX_NDIV) 
		  {
            rels[i][5 + rels[i][4]] = 1 + mpqs_nFB + mpqs_nFBk + j; 
			rels[i][4]++;
            ulqx /= (u32_t)p;
          } 
		  else 
		  {
            goto next;
          }
        }
        
		if (ulqx == 1) 
			break;
      }
	}
#endif

    if (rels[i][4] < MPQS_TD_MAX_NDIV - mpqs_nAdiv) 
	{
      nr = rels[i][4]; 
	  ind = 0;
      
	  for (j = 0; j < mpqs_nAdiv_total; j++)
	  {
        if (mpqs_Adiv_active[j]) 
		{
          rels[i][5 + nr + ind] = 1 + mpqs_nFB + mpqs_nFBk + j;
          ind++;
        }
	  }

      rels[i][4] += mpqs_nAdiv;
    } 
	else 
	{
      goto next;
    }

    if (ulqx == 1) 
	{
      llp = (i64_t *)(mpqs_relations[mpqs_nrels]);
      *llp = axb;
      nr = rels[i][4];
      
	  for (j = 0; j <= nr; j++)
	  {
        mpqs_relations[mpqs_nrels][j + 4] = rels[i][j + 4];
	  }

	  mpqs_nrels++;
      goto next;
    }

    if (ulqx < 0x100000) 
	{ 
      /*printf("%d ",x);*/  /* !!! */    
	  
	  if (ulqx <= mpqs_pmax) 
	  {
        printf("%d %" PRIu64 " %u ", x, qx, ulqx);
		return -1; /* CJM, 11/30/04. */
        complain("dec.2\n");
      }

      lp = ulqx;
      lpi = mpqs_hash_lp(lp);
      
	  if (lpi == MPQS_HASH_OVERFLOW) 
		  goto next;
      
	  if (lpi < mpqs_nprels) 
	  {
        if (mpqs_lp[lpi] != lp) 
			complain("lp");

        llp = (i64_t *)(mpqs_comb_rels[mpqs_ncrels]);
        *llp = axb;
        llp = (i64_t *)(mpqs_part_rels[lpi]);
        axb = *llp;
        llp = (i64_t *)(mpqs_comb_rels[mpqs_ncrels]); 
		llp++;
        *llp = axb;
        nr = rels[i][4];
        
		for (j = 0; j < nr; j++)
          mpqs_comb_rels[mpqs_ncrels][j + 9] = rels[i][j + 5];
        
		for (j=0; j<mpqs_part_rels[lpi][4]; j++)
          mpqs_comb_rels[mpqs_ncrels][j + 9 + nr] = mpqs_part_rels[lpi][j + 5];
        
		mpqs_comb_rels[mpqs_ncrels][8] = nr + 256*mpqs_part_rels[lpi][4];
        mpqs_ncrels++;
      } 
	  else 
	  {
        llp = (i64_t *)(mpqs_part_rels[mpqs_nprels]);
        *llp = axb;
        nr = rels[i][4];
        
		for (j = 0; j <= nr; j++)
          mpqs_part_rels[mpqs_nprels][j + 4] = rels[i][j + 4];
        
		mpqs_nprels++;
      }
    } 
	/*else printf("<");*/

next:
    ; /* Null statement to hush compiler warning. */
  }

  if (mpqs_nsp < (mpqs_nrels + mpqs_ncrels))
      mpqs_excess = mpqs_nrels + mpqs_ncrels-mpqs_nsp;
  else 
	  mpqs_excess = 0;

  return 0;
}


/*
order of the primes:
-1, primes in mpqs_FBk, primes in mpqs_FB, primes in mpqs_Adiv_total
*/

/*
static u32_t mpqs_gauss_mask0[]  =
{
  0x80000000UL, 0x40000000UL, 0x20000000UL, 0x10000000UL,
  0x08000000UL, 0x04000000UL, 0x02000000UL, 0x01000000UL,
  0x00800000UL, 0x00400000UL, 0x00200000UL, 0x00100000UL,
  0x00080000UL, 0x00040000UL, 0x00020000UL, 0x00010000UL,
  0x00008000UL, 0x00004000UL, 0x00002000UL, 0x00001000UL,
  0x00000800UL, 0x00000400UL, 0x00000200UL, 0x00000100UL,
  0x00000080UL, 0x00000040UL, 0x00000020UL, 0x00000010UL,
  0x00000008UL, 0x00000004UL, 0x00000002UL, 0x00000001UL
};
*/
static u32_t mpqs_gauss_mask[]  =
{
  0x00000001UL, 0x00000002UL, 0x00000004UL, 0x00000008UL,
  0x00000010UL, 0x00000020UL, 0x00000040UL, 0x00000080UL,
  0x00000100UL, 0x00000200UL, 0x00000400UL, 0x00000800UL,
  0x00001000UL, 0x00002000UL, 0x00004000UL, 0x00008000UL,
  0x00010000UL, 0x00020000UL, 0x00040000UL, 0x00080000UL,
  0x00100000UL, 0x00200000UL, 0x00400000UL, 0x00800000UL,
  0x01000000UL, 0x02000000UL, 0x04000000UL, 0x08000000UL,
  0x10000000UL, 0x20000000UL, 0x40000000UL, 0x80000000UL
};

#define PRUNING

static void mpqs_matrix_init()
{
  u32_t i, j;
  u16_t nr, ii, pi;

#ifdef PRUNING
  u16_t mpqs_gauss_wt[300];  /* to be improved !!!! */
  memset(mpqs_gauss_wt, 0, sizeof(mpqs_gauss_wt));
#endif

  mpqs_gauss_n = mpqs_nrels + mpqs_ncrels;
  
  if (mpqs_gauss_n >= MPQS_GAUSS_MAX) 
	  return;
  
  mpqs_gauss_m = mpqs_nsp;
  
  if (mpqs_gauss_m >= mpqs_gauss_n) 
	  complain("gauss: no excess\n");

  /* build matrix */
  mpqs_gauss_n32 = (mpqs_gauss_n + 31)/32;
  
  for (i = 1; i < mpqs_gauss_m; i++)
    mpqs_gauss_row[i] = mpqs_gauss_row[i-1] + mpqs_gauss_n32;
  
  for (i = 0; i < mpqs_gauss_m; i++)
  {
    for (j = 0; j<mpqs_gauss_n32; j++)
      mpqs_gauss_row[i][j] = 0;
  }

  for (i = 0; i < mpqs_nrels; i++) 
  {
    u32_t gauss_mask;

    nr = mpqs_relations[i][4]; 
	gauss_mask = mpqs_gauss_mask[i%32];
    
	for (j = 0; j < nr; j++) 
	{
      pi = mpqs_relations[i][5 + j];
      mpqs_gauss_row[pi][i/32] ^= gauss_mask;

#ifdef PRUNING
      mpqs_gauss_wt[pi]++;
#endif
    }
  }

  for (i = 0; i < mpqs_ncrels; i++) 
  {
    u32_t gauss_mask;

    nr = mpqs_comb_rels[i][8]; 
	nr = (nr & 255) + (nr >> 8);
    ii = i + mpqs_nrels; 
	gauss_mask = mpqs_gauss_mask[ii % 32];
    
	for (j = 0; j < nr; j++) 
	{
      pi = mpqs_comb_rels[i][9 + j];
      mpqs_gauss_row[pi][ii/32] ^= gauss_mask;

#ifdef PRUNING
      mpqs_gauss_wt[pi]++;
#endif
    }
  }

  for (i = 0; i < mpqs_gauss_m; i++) 
	  mpqs_gauss_c[i] = -1;
  
  mpqs_gauss_k = 0;
  
  for (i = 0; i < mpqs_gauss_n; i++) 
	  mpqs_gauss_d[i] = -1;

#ifdef PRUNING
/* simple pruning */
{
  u32_t k, k32;

  for (j = 0; j < mpqs_gauss_m; j++)
  {
    if (mpqs_gauss_wt[j] == 1) 
	{
      k32 = 0;
      
	  while (mpqs_gauss_row[j][k32] == 0) 
		  k32++;
      
	  for (k = 0; k < 32; k++)
        if (mpqs_gauss_row[j][k32] & mpqs_gauss_mask[k]) 
			break;
      
	  k += (k32 << 5);
      
	  if (k >= mpqs_gauss_n) 
		  complain("pruning\n");

      mpqs_gauss_d[k] = j; 
	  mpqs_gauss_c[j] = k;
	  /* No need to clear the k-th column since it will not appear in a solution. */
    }
  }
}
#endif
}

static int mpqs_matrix()
{
  u32_t i, j, l, t, k32;
  u32_t mask;
 
  /* solve matrix */
  while (mpqs_gauss_k < mpqs_gauss_n) 
  {
#ifdef PRUNING
    assert(mpqs_gauss_k < sizeof(mpqs_gauss_d)/sizeof(mpqs_gauss_d[0]));
    if (mpqs_gauss_d[mpqs_gauss_k] >= 0) 
	{ 
		mpqs_gauss_k++; 
		continue; 
	}
#endif

    mask = mpqs_gauss_mask[mpqs_gauss_k%32]; 
	k32 = mpqs_gauss_k / 32;
    j = 0;
    
	while ((j < mpqs_gauss_m) && (mpqs_gauss_c[j] >= 0 ||
		   ((mpqs_gauss_row[j][k32] & mask) == 0))) 
		j++;
    
	if (j == mpqs_gauss_m) 
	{ 
		assert(mpqs_gauss_k < sizeof(mpqs_gauss_d)/sizeof(mpqs_gauss_d[0]));
		mpqs_gauss_d[mpqs_gauss_k] = -1;
		break; 
	}
    
	mpqs_gauss_d[mpqs_gauss_k] = j; 
	mpqs_gauss_c[j] = mpqs_gauss_k;
    
	for (t = 0; t < j; t++)
      if (mpqs_gauss_row[t][k32] & mask)
        for (l = k32; l < mpqs_gauss_n32; l++)
          mpqs_gauss_row[t][l] ^= mpqs_gauss_row[j][l];
    
	for (t = j + 1; t < mpqs_gauss_m; t++)
	{
      if (mpqs_gauss_row[t][k32] & mask)
	  {
        for (l = k32; l < mpqs_gauss_n32; l++)
          mpqs_gauss_row[t][l] ^= mpqs_gauss_row[j][l];
	  }
	}

    mpqs_gauss_k++;
  }

  if (mpqs_gauss_k >= mpqs_gauss_n) 
  {
    return 0;
  }

  for (i = 0; i < mpqs_gauss_n; i++) 
	  mpqs_sol[i] = 0;

  if (mpqs_gauss_d[mpqs_gauss_k] != -1) 
	  complain("gauss1\n");

  for (i = 0; i < mpqs_gauss_k; i++)
  {
    if ((mpqs_gauss_d[i] != -1) && (mpqs_gauss_row[mpqs_gauss_d[i]][k32] & mask))
        mpqs_sol[i] = 1;
  }

  mpqs_sol[mpqs_gauss_k] = 1;
  mpqs_gauss_k++;
  return 1;
}

static int mpqs_final()
{
  u32_t j, k;
  u16_t nr, nr1, nr2;
  u32_t p, prod1, prod2;
  unsigned long int *up;
  int split;

  if (!mpqs_nfactors) 
  {
    mpqs_nfactors = 1;
    mpz_set(mpqs_factors[0],mpqs_N);
    mpqs_f_status[0]=0;
  }
  
  while (mpqs_matrix()) 
  {
#ifdef MPQS_STAT
stat_mpqs_ntrials++;
#endif
   
	for (j = 0; j < mpqs_nsp; j++) 
		mpqs_exp[j] = 0;
    
	for (j = 0; j < mpqs_nrels; j++)
	{
      if (mpqs_sol[j]) 
	  {
        nr = mpqs_relations[j][4];
        
		if (j & 1) /* 'random' distribution of full relations */
          for (k = 0; k < nr; k++) mpqs_exp[mpqs_relations[j][5 + k]]--;
        else
          for (k = 0; k < nr; k++) mpqs_exp[mpqs_relations[j][5 + k]]++;
      }
	}

    for (j = 0; j < mpqs_ncrels; j++)
	{
      if (mpqs_sol[mpqs_nrels + j]) 
	  {
        nr = mpqs_comb_rels[j][8];
        nr1 = nr & 255; 
		nr2 = nr >> 8;
        
		for (k = 0; k < nr1; k++) 
			mpqs_exp[mpqs_comb_rels[j][9 + k]]++;
        
		for (k = 0; k < nr2; k++) 
			mpqs_exp[mpqs_comb_rels[j][9 + nr1 + k]]--;
      }
	}

    for (j = 0; j < mpqs_nsp; j++)
	{
      if (mpqs_exp[j] & 1) 
		  complain("final.odd\n");
      else 
		  mpqs_exp[j] >>= 1;
	}

    mpz_set_ui(mpqs_sq1,1); 
    mpz_set_ui(mpqs_sq2,1); 
	prod1 = 1;
	prod2 = 1;

	for (j = 1; j < mpqs_nsp; j++)
	{
      if (mpqs_exp[j]) 
	  {
        if (j < (u32_t)(1 + mpqs_nFBk)) 
		{
          p = (u32_t)(mpqs_FBk[j - 1]);
        } 
		else 
		if (j < (u32_t)(1 + mpqs_nFB + mpqs_nFBk)) 
		{
          k = j - mpqs_nFBk - 1;
          p = (u32_t)(mpqs_FB[2*k]);
        } 
		else 
		{
          k = j - mpqs_nFB - mpqs_nFBk - 1;
          p = (u32_t)(mpqs_Adiv_all[k]);
        }

        if (mpqs_exp[j] > 0) 
		{
          for (k = 0; k < (u32_t)mpqs_exp[j]; k++) 
		  {
            prod1 *= p;
            if (prod1 & 0xfff00000) 
			{  
			  /* p<4096 */
              mpz_mul_ui(mpqs_sq1, mpqs_sq1, prod1);
              mpz_mod(mpqs_sq1, mpqs_sq1, mpqs_N);
#ifdef MPQS_STAT
stat_final_mulmod++;
#endif
              prod1 = 1;
            }
          }
        } 
		else 
		{
          for (k = 0; k < (u32_t)(0 - mpqs_exp[j]); k++) 
		  {
            prod2 *= p;
            if (prod2 & 0xfff00000) 
			{  
			  /* p < 4096 */
              mpz_mul_ui(mpqs_sq2, mpqs_sq2, prod2);
              mpz_mod(mpqs_sq2, mpqs_sq2, mpqs_N);
#ifdef MPQS_STAT
stat_final_mulmod++;
#endif
              prod2 = 1;
            }
          }
        }
      }
	}

    if (prod1 > 1) 
	{
      mpz_mul_ui(mpqs_sq1, mpqs_sq1, prod1);
      mpz_mod(mpqs_sq1, mpqs_sq1, mpqs_N);
    }
    
	if (prod2 > 1) 
	{
      mpz_mul_ui(mpqs_sq2, mpqs_sq2, prod2);
      mpz_mod(mpqs_sq2, mpqs_sq2, mpqs_N);
    }

    for (j = 0; j < mpqs_nrels; j++)
	{
      if (mpqs_sol[j]) 
	  {
        up = (unsigned long int *)(mpqs_relations[j]); 
		/* printf("Q: %Lu  ",*up); */
        mpz_set_ull(mpqs_gmp_help, *up);
        
		if (j & 1) 
		{
          mpz_mul(mpqs_sq1, mpqs_sq1, mpqs_gmp_help);
          mpz_mod(mpqs_sq1, mpqs_sq1, mpqs_N);
        } 
		else 
		{
          mpz_mul(mpqs_sq2, mpqs_sq2, mpqs_gmp_help);
          mpz_mod(mpqs_sq2, mpqs_sq2, mpqs_N);
        }
      }
	}

    for (j = 0; j < mpqs_ncrels; j++)
	{
      if (mpqs_sol[mpqs_nrels+j]) 
	  {
        up = (unsigned long int *)(mpqs_comb_rels[j]); 
		/* printf("Q: %Lu  ", *up); */
        mpz_set_ull(mpqs_gmp_help, *up);
        mpz_mul(mpqs_sq2, mpqs_sq2, mpqs_gmp_help);
        mpz_mod(mpqs_sq2, mpqs_sq2, mpqs_N);
        up++;
        mpz_set_ull(mpqs_gmp_help, *up);
        mpz_mul(mpqs_sq1, mpqs_sq1, mpqs_gmp_help);
        mpz_mod(mpqs_sq1, mpqs_sq1, mpqs_N);
      }
	}

    mpz_add(mpqs_gmp_help, mpqs_sq1, mpqs_sq2);
    mpz_gcd(mpqs_gmp_help, mpqs_gmp_help, mpqs_N);
    
	if (mpz_cmp_ui(mpqs_gmp_help, 1) == 0) 
		continue;
    
	if (mpz_cmp(mpqs_gmp_help, mpqs_N) == 0) 
		continue;

    for (j = 0; j < mpqs_nfactors; j++)
	{
      if (mpqs_f_status[j] == 0) 
	  {
        mpz_gcd(mpqs_sq1, mpqs_gmp_help, mpqs_factors[j]);
        
		if (mpz_cmp_ui(mpqs_sq1, 1) == 0)
			continue;
        
		if (mpz_cmp(mpqs_sq1, mpqs_factors[j]) == 0) 
			continue;
        
		mpz_set(mpqs_factors[mpqs_nfactors], mpqs_sq1);
        mpz_divexact(mpqs_factors[j], mpqs_factors[j], mpqs_sq1);
        mpqs_f_status[j] = psp(mpqs_factors[j]);
        mpqs_f_status[mpqs_nfactors] = psp(mpqs_factors[mpqs_nfactors]);
        mpqs_nfactors++;
        break;
      }
	}

    split = 1;
    
	for (j = 0; j < mpqs_nfactors; j++) 
	{
		if (mpqs_f_status[j] == 0) 
			split = 0;
	}

    if (split) 
		return 1;
    
	if (mpqs_nfactors >= MPQS_MAX_FACTORS)
      complain("final: too many factors\n");
  }

  return 0;
}

static long mpqs_factor0(mpz_t N, size_t max_bits, mpz_t **factors, u16_t retry)
{
  int ev;
  u32_t i;
  long err;

#ifdef MPQS_STAT
  stat_mpqs_nsieves = 0; 
  stat_mpqs_nsurvivors = 0;
  stat_mpqs_ntrials = 0; 
  stat_mpqs_ndiv = 0;
#endif

  if (!mpqs_isinit) 
	  mpqs_init();

  if (mpz_sizeinbase(N, 2) > 96) 
  {
    fprintf(stderr,"warning: mpqs called with >96 Bit\n");
    return -2;
  }

  mpz_set(mpqs_N, N);
  mpqs_choose_multiplier();
  mpqs_choose_parameter(retry);
  mpqs_generate_FB();
  
  if ((err = mpqs_SI_init()) < 1) 
  {
    printf("%ld ", err);
    fprintf(stderr,"warning: mpqs: error in self-initialization ");
    mpz_out_str(stderr, 10, N); 
	printf("\n");
    return -3;
  }

  while (1, 1) 
  { 
	/* printf("%ld, %ld; ", stat_mpqs_nsieves, mpqs_nrels); */

    if (!mpqs_next_pol()) 
	{
      if (retry) 
	  {
        fprintf(stderr,"warning: not enough polynomials in mpqs with ");
        mpz_out_str(stderr, 10, N); 
		fprintf(stderr,"\n");
      }
    
	  return -4;
    }

    mpqs_sieve();

#ifdef MPQS_STAT
stat_mpqs_nsieves++;
#endif
    
    ev = mpqs_evaluate();

#ifdef MPQS_STAT
stat_mpqs_nsurvivors += ev; /* printf("%d %Ld  ", ev, mpqs_C); */
#endif

	if (ev) 
	{
		{ 
			int mpqs_dec_res; /* CJM, 11/30/04. */

			if ((mpqs_dec_res = mpqs_decompose()))
				return -1;
		}

		if (mpqs_nsurvivors && (mpqs_excess > MPQS_MIN_EXCESS)) 
		{

        mpqs_matrix_init();
        if (mpqs_final()) 
		{
#ifdef MPQS_STAT
printf("Stat: %lu %lu %lu %lu %lu %lu %lu \n",stat_mpqs_nsieves,stat_mpqs_nsurvivors,stat_mpqs_ntrials,stat_mpqs_ndiv,mpqs_nrels,mpqs_nprels,mpqs_ncrels);
#endif
          for (i = 0; i < mpqs_nfactors; i++)
		  {
            if (mpz_sizeinbase(mpqs_factors[i], 2) > max_bits) 
				return 0;
		  }

          *factors = mpqs_factors;
          return mpqs_nfactors;
        }
      }
    }

    if (mpqs_nrels > MPQS_MAX_NRELS - MPQS_MIN_RELBUFFER) 
	{
      fprintf(stderr, "warning: too many relations in mpqs\n");
      return -5;
    }
  }

  return -1;
}

long mpqs_factor(mpz_t N, size_t max_bits, mpz_t **factors)
{
  long err;

  err = mpqs_factor0(N, max_bits, factors, 0);
  if (err == -4) err = mpqs_factor0(N, max_bits, factors, 1);
  if (err == -4) err = mpqs_factor0(N, max_bits, factors, 3);
  return err;
}

