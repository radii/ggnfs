/*
  Copyright (C) 2000 Jens Franke
  This file is part of mpqs4linux, distributed under the terms of the 
  GNU General Public Licence and WITHOUT ANY WARRANTY.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

volatile extern u32_t modulo32;
u32_t gcd32(u32_t x, u32_t y);
int jac32(u32_t x,u32_t y);
u32_t modpow32(u32_t x,u32_t a);
u32_t modsqrt32(u32_t x);
u32_t modinv32(u32_t x);

static inline u32_t modsq32(u32_t x)
{
  u32_t rv,clobber;

  asm ("mult %1,%1\n"
       "mfhi %1\n"
       "dsll %1,32\n"
       "mflo %0\n"
       "dsll %0,32\n"
       "dsrl %0,32\n"
       "or %0,%0,%1\n"
       "lw %1,modulo32\n"
       "dsll %1,32\n"
       "dsrl %1,32\n"
       "ddivu $0,%0,%1\n"
       "mfhi %0" : "=&r" (rv), "=r" (clobber) : "1" (x)  );
  return rv;
}

static inline u32_t modmul32(u32_t x,u32_t y)
{
  u32_t rv,clobber;

  asm ("mult %0,%1\n"
       "mfhi %1\n"
       "dsll %1,32\n"
       "mflo %0\n"
       "dsll %0,32\n"
       "dsrl %0,32\n"
       "or %0,%0,%1\n"
       "lw %1,modulo32\n"
       "dsll %1,32\n"
       "dsrl %1,32\n"
       "ddivu $0,%0,%1\n"
       "mfhi %0" : "=r" (rv), "=r" (clobber) : "0" (x), "1" (y)  );
  return rv;
}

static inline u32_t modsub32(u32_t minuend,u32_t subtrahend)
{
  /* This appears to save 1/3, compared to the generic version. */
  u32_t rv;

  asm ("sleu %0,%2,%1\n"
       "subu %0,%0,1\n"
       "and %0,%0,%3\n"
       "addu %0,%0,%1\n"
       "subu %0,%0,%2\n" : "=&r" (rv) :
       "r" (minuend), "r" (subtrahend) , "r" (modulo32) );
  return rv;
}

/* I am indebted to J. Leherbauer for pointing this out to me. */

static inline u32_t modadd32(u32_t x,u32_t y)
{
  return modsub32(x,modulo32-y);
}
