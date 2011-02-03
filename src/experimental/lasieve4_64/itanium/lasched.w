@* Scheduling function for the lattice siever.
Copyright (C) 2002 Jens Franke, T. Kleinjung.
This file is part of gnfs4linux, distributed under the terms of the 
GNU General Public Licence and WITHOUT ANY WARRANTY.

You should have received a copy of the GNU General Public License along
with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.
@(lasched.h@>=
u32_t *lasched(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t,u32_t);

@
@c
#include <sys/types.h>
#include <limits.h>

#include "siever-config.h"
#include "lasched.h"
#include "../if.h"

#define L1_SIZE (1<<L1_BITS)
#define i_bits (I_bits-1)
#define n_i (1<<i_bits)

/* Amount by which to shift a larger number to provide space for a
   16-bit number. */

#define U16_SHIFT (CHAR_BIT*sizeof(u16_t))

u32_t *lasched0(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t *lasched1(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t *lasched2(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);
u32_t *lasched3(u32_t*,u32_t*,u32_t*,u32_t,u32_t**,u32_t);

u32_t *
lasched(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs,ot)
     u32_t *ri;
     u32_t *ij_ptr,*ij_ptr_ub,n1_j,**sched_ptr,fbi_offs,ot;
{
  u32_t ij,ij_ub;
  u32_t ot_mask,ot_tester;
  u32_t *ru;

  ij_ub=n1_j<<i_bits;
  ru=(u32_t*)ri;

  if(ot!=0) {
    if(ot==1) return lasched1(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
    if(ot==2) return lasched2(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
    if(ot==3) return lasched3(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);
    ot_tester=(ot&1)|((ot&2)<<(i_bits-1));
    ot_mask=n_i|1;
  } else return lasched0(ri,ij_ptr,ij_ptr_ub,n1_j,sched_ptr,fbi_offs);

  while ( ij_ptr < ij_ptr_ub) {
    u16_t a,b;

    a=n_i-(ru[0]&(n_i-1));
    b=n_i-(ru[1]&(n_i-1));
    if(ot == 0) ij=*ij_ptr;
    else @<Calculate first sieving event from |ri|@>@;
    while(ij<ij_ub) {
      u16_t i;

      *(sched_ptr[ij>>L1_BITS]++)=(fbi_offs<<U16_SHIFT) | (ij&(L1_SIZE-1));
      i=ij&(n_i-1);
      if(i<b) ij+=ru[0];
      if(i>=a) ij+=ru[1];
    }
    ru+=2;
    *(ij_ptr++)=ij-ij_ub;
    fbi_offs++;
  }
  return ru;
}

@
@<Calculate first sieving event from |ri|@>=
{
  ij=0;
  if( (ru[0]&ot_mask) == ot_tester ) ij=ru[0];
  else {
    if( (ru[1]&ot_mask) == (ot_tester^n_i) ) ij=ru[1];
    else {
      if((ru[0]&(n_i-1))<=(ru[1]&(n_i-1)) && ru[0]<=ru[1]) {
	/* This corresponds to the line
	   |if(b+c<=A && s<=t)| in recurrence6.w */
	if((ru[0]&(n_i-1))==(ru[1]&(n_i-1))) ij=ru[1]-ru[0];
	else ij=n_i;
	if(ot != 2)
	  Schlendrian("Exceptional situation for oddness type %u ?\n",
		      ot);
      }
      else ij=ru[0]+ru[1];
    }
  }
  ij=(ij+((~ot_tester)&n_i))/2;
}
