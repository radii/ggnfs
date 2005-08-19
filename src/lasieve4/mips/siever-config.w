@* Configuration file for GNFS lattice siever.


Copyright (C) 2002 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@(siever-config.h@>=
#define L1_BITS 15
#define ULONG_RI
#define ULL_NO_UL

#define PREINVERT
#define modinv32 asm_modinv32
#define ASM_LINESIEVER
#define HAVE_ASM_LASIEVE_SETUP
#define FLOAT_SETUP_BOUND1 0x4000000
#define FLOAT_SETUP_BOUND2 4503599627370496.0
#define ASM_SCHEDTDSIEVE
#define ASM_SCHEDSIEVE1
#define SCHEDULING_FUNCTION_CALCULATES_RI

@
@(siever-config.h@>=
#if 1
#define GGNFS_BIGENDIAN
#endif

@ Use this if your compile cannot handle inline funtions.
@(siever-config.h@>=
#if 0
#define inline
#endif

@
@(siever-config.h@>=
#if 1
#define NEED_ASPRINTF
#endif

@
@(siever-config.h@>=
#if 1
#define NEED_GETLINE
#endif

@ Machine specific initializations.
@c
static void
siever_init(void)
{
  return;
}

@ Bounds for the primes treated by the various types of schedule.
@(siever-config.h@>=
#if 0
#if 1
#define N_PRIMEBOUNDS 7
#else
#define N_PRIMEBOUNDS 1
#endif
#else
#define N_PRIMEBOUNDS 4
#endif

@
@c
#if 0
#if 0
const ulong schedule_primebounds[N_PRIMEBOUNDS]={0x20000,0x80000,
		0x200000,0x800000,ULONG_MAX};

const ulong schedule_sizebits[N_PRIMEBOUNDS]={19,20,22,24,32};
#else
#if 1
const ulong schedule_primebounds[N_PRIMEBOUNDS]=
	{0x100000,0x200000,0x400000,0x800000,0x1000000,0x2000000,ULONG_MAX};

const ulong schedule_sizebits[N_PRIMEBOUNDS]={20,21,22,23,24,25,32};
#else
const ulong schedule_primebounds[N_PRIMEBOUNDS]={ULONG_MAX,};

const ulong schedule_sizebits[N_PRIMEBOUNDS]={32};
#endif
#endif
#else
const ulong schedule_primebounds[N_PRIMEBOUNDS]={0x100000,0x200000,0x4000000,ULONG_MAX};

const ulong schedule_sizebits[N_PRIMEBOUNDS]={20,21,22,32};
#endif