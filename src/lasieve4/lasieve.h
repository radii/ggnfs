/* lasieve.h
  Hacked up for inclusion in GGNFS by Chris Monico.

  Copyright (C) 2001 Jens Franke.
  This file is part of gnfs4linux, distributed under the terms of the 
  GNU General Public Licence and WITHOUT ANY WARRANTY.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#ifndef _LASIEVE_H
#define _LASIEVE_H

#include <gmp.h>
#include <stdlib.h>
#include <stdarg.h> 
#include <stdio.h> 

#include "asm/siever-config.h"
#include "if.h"

#ifdef _MSC_VER
double rint(double);
#endif

/* fbgen.h */
//u32_t root_finder(u32_t *root_buf,mpz_t *A,u32_t adeg,u32_t p);
extern u32_t * polr, polr_alloc;

/* input-poly.h */
void input_poly(mpz_t,mpz_t**,int32_t*,mpz_t**,int32_t*,mpz_t,FILE*);

/* lasieve-prepn.h */
void lasieve_setup(uint32_t*,uint32_t*,uint32_t,int32_t,int32_t,int32_t,int32_t,uint32_t*);

/* primgen32.h */
typedef struct{
  uint32_t Pind;
  uint32_t first_in_sieve;
  uint32_t nPrim;
  uint32_t Prime;
  unsigned char*PDiffs;
  uint32_t use_private;
  uint32_t PDiffs_allocated;
} pr32_struct;

void initprime32(pr32_struct*ps);
uint32_t firstprime32(pr32_struct*ps);
uint32_t nextprime32(pr32_struct*ps);
void clearprime32(pr32_struct*ps);
void p32_seek(pr32_struct*ps,uint32_t lb);

/* real-poly-aux.h */
void tpol(double *rop,double *op,uint32_t deg,int32_t x0,int32_t x1,int32_t y0,int32_t y1);
double rpol_eval(double *p,int32_t d,double x,double y);
double rpol_lb(double *pol,uint32_t poldeg,double a,double b);


/* recurrence6.h */
void rec_info_init(uint32_t A,uint32_t ub);
uint32_t get_recurrence_info(uint32_t*res_ptr,uint32_t p,uint32_t r);

/* redu2.h */
void reduce2(int32_t*,int32_t*,int32_t*,int32_t*,int32_t,int32_t,int32_t,int32_t,float);
extern int32_t n_iter;

/* primgen32.c */
uint32_t pr32_seek(pr32_struct * ps, uint32_t lb);

/* psp.c */
int psp(mpz_t);

/* mpqs.c */
long mpqs_factor(mpz_t N, long max_bits, mpz_t **factors);

#endif
