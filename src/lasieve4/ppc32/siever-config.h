/* Configuration file for GNFS lattice siever.

Copyright (C) 2002 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA. */

#include <stdint.h>
#include <gmp.h> 

#define L1_BITS 14
#define ULONG_RI
typedef long int i32_t;
typedef short int i16_t;
typedef unsigned short u16_t;
typedef long long int i64_t;
#define U32_MAX 0xffffffff
#define I32_MAX INT_MAX
#define ULL_NO_UL

typedef uint32_t u32_t;
typedef uint64_t u64_t;

#define PREINVERT
#define BIGENDIAN
#define NEED_GETLINE
#define N_PRIMEBOUNDS 7
#define GNFS_CS32

#include "32bit.h"

static inline void siever_init() { };

u32_t *medsched(u32_t*,u32_t*,u32_t*,u32_t**,u32_t,u32_t);
u32_t *lasched(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t,u32_t);
long mpqs_factor(mpz_t N, long max_bits, mpz_t **factors);
int psp(mpz_t n);
