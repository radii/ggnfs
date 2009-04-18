@* Configuration file for GNFS lattice siever.

Copyright (C) 2001 Jens Franke.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@(siever-config.h@>=
#define L1_BITS 15
#define ULONG_RI
#define HAVE_CMOV
#define HAVE_SSIMD
typedef unsigned long u32_t;
typedef long int i32_t;
typedef short int i16_t;
typedef unsigned short u16_t;
typedef unsigned long long u64_t;
typedef long long int i64_t;
#define U32_MAX 0xffffffff
#define I32_MAX INT_MAX
#define ULL_NO_UL

#define PREINVERT
#define HAVE_ASM_GETBC
#define ASM_SCHEDSIEVE
#define ASM_SCHEDTDSIEVE2
#define ASM_LINESIEVER
#define ASM_LINESIEVER3
#define ASM_LINESIEVER2
#define ASM_LINESIEVER1
#define ASM_SEARCH0
#define ASM_TDSLINIE1
#define ASM_TDSLINIE2
#define ASM_TDSLINIE3
#define ASM_TDSLINIE

#if 1
#define TDS_FB_PREFETCH(x) asm volatile ("prefetcht0 (%0)" : : "r" (x));
#endif

@ Use this to determine how many primes are treated by true trial division
before using the trial division sieve.
@(siever-config.h@>=
#if 0
#define SET_TDS_PBOUND(ni,jps,ncand) ( ncand==0 ? 0xffffffff : (ni*jps)/ncand)
#endif

@ Trial division using a modular inverse of the factor base prime appears
be give little or no advantage. The code needed to try this comes with
this distribution, however.
@c
#if 1
#define ASM_MPZ_TD
#endif

@
@(siever-config.h@>=
#if 1
void MMX_TdAllocate(int,size_t,size_t);
u16_t* MMX_TdInit(int,u16_t*,u16_t*,u32_t*,int);
void MMX_TdUpdate(int,int);
#define MMX_TD
#endif

@ Bounds for the primes treated by the various types of schedule.
@(siever-config.h@>=
#if 1
#define N_PRIMEBOUNDS 9
#else
#define N_PRIMEBOUNDS 1
#endif

@ Machine specific initializations.
@c
static void
siever_init(void)
{
  init_montgomery_multiplication();
}

@
@c
#if 0
const ulong schedule_primebounds[N_PRIMEBOUNDS]={0x20000,0x80000,
		0x200000,0x800000,ULONG_MAX};

const ulong schedule_sizebits[N_PRIMEBOUNDS]={19,20,22,24,32};
#else
#if 1
const ulong schedule_primebounds[N_PRIMEBOUNDS]=
	{0x100000,0x200000,0x400000,0x800000,0x1000000,0x2000000,0x4000000,
	0x8000000,ULONG_MAX};

const ulong schedule_sizebits[N_PRIMEBOUNDS]={20,21,22,23,24,25,26,27,32};
#else
const ulong schedule_primebounds[N_PRIMEBOUNDS]={ULONG_MAX,};

const ulong schedule_sizebits[N_PRIMEBOUNDS]={32};
#endif
#endif
