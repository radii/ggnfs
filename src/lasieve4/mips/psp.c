/*
Copyright (C) 2001 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>

extern ulong montgomery_inv_n;
extern ulong *montgomery_modulo_n;
extern ulong montgomery_modulo_R2[3];
extern ulong montgomery_64bit;

int psp(mpz_t n)
{
  if(mpz_size(n)<2) {
    if(mpz_sgn(n)!=0) return pt64(n[0]._mp_d[0]);
    return 0;
  }
  return mpz_probab_prime_p(n,1);
}

