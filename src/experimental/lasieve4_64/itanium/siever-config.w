@* Configuration file for GNFS lattice siever.

Copyright (C) 2001 Jens Franke.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.

@(siever-config.h@>=
#define L1_BITS 16
#define ULONG_RI
typedef unsigned u32_t;
typedef int i32_t;
typedef short int i16_t;
typedef unsigned short u16_t;
typedef unsigned long u64_t;
typedef long int i64_t;
#define U32_MAX 0xffffffff
#define I32_MAX INT_MAX
#define asm_modinv32(x) asm_modinv32a(x,modulo32)
#define HAVE_ASM_GETBC
#define ASM_SCHEDSIEVE
#define PREINVERT
#define MMX_TD
#define ASM_LINESIEVER
u16_t *MMX_TdInit(int,u16_t*,u16_t*,u32_t*,int);
void MMX_TdUpdate(int,int);
u32_t *MMX_Td(u32_t*,int,u16_t);

#if 1
#define HAVE_ASM_LASIEVE_SETUP
#define FLOAT_SETUP_BOUND1 0x80000000
#define FLOAT_SETUP_BOUND2 4503599627370496.0
#endif

@ Trial division using a modular inverse of the factor base prime appears
be give little or no advantage. The code needed to try this comes with
this distribution, however.
@c
#if 1
#define ASM_MPZ_TD
#endif

@ Bounds for the primes treated by the various types of schedule.
@(siever-config.h@>=
#define N_PRIMEBOUNDS 3

@ Machine specific initializations.
@c
static void
siever_init(void)
{
}

@
@c
const ulong schedule_primebounds[N_PRIMEBOUNDS]={0x200000,0x4000000,ULONG_MAX};

const ulong schedule_sizebits[N_PRIMEBOUNDS]={21,22,32};
