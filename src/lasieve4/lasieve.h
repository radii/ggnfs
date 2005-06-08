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

#ifdef _MSC_VER
double rint(double);
#endif

typedef unsigned long long ullong;

/* fbgen.h */
u32_t root_finder(u32_t *root_buf,mpz_t *A,u32_t adeg,u32_t p);
extern u32_t * polr,polr_alloc;

/* gmp-aux.h */
void adjust_mpz_bufsize(mpz_t**x,size_t*alloc_ptr,size_t size,size_t increment);
int string2mpz(mpz_t rop,char*x,int base);
#if __GNU_MP_VERSION < 3 || __GNU_MP_VERSION_MINOR == 0
#define NEED_MPZ_MUL_SI
void mpz_mul_si(mpz_t x,mpz_t y,long int z);
#endif

#define mpz_add_si(r,o,s) \
  ({ int _s; _s= (s); \
     _s>=0 ? mpz_add_ui(r,o,(ulong)(_s)) : mpz_sub_ui(r,o,(ulong)(-_s)); })

void mpz_set_ull(mpz_t targ,ullong src);
ullong mpz_get_ull(mpz_t src);
int mpz_cmp_ull(mpz_t op1,ullong op2);
#ifdef ULL_NO_UL
void mpz_ull_init();
void mpz_mul_ull(mpz_t rop,mpz_t op1,ullong op2);
int mpz_fits_sllong_p(mpz_t);
int mpz_fits_ullong_p(mpz_t);
long long int mpz_get_sll(mpz_t);
void mpz_set_sll(mpz_t,long long int);
unsigned long long mpz_get_ull(mpz_t);
void mpz_set_ull(mpz_t,unsigned long long);
#else
#define mpz_ull_init()
#define mpz_mul_ull mpz_mul_ui
#define mpz_fits_sllong_p mpz_fits_slong_p
#define mpz_fits_ullong_p mpz_fits_ulong_p
#define mpz_get_sll mpz_get_si
#define mpz_set_sll mpz_set_si
#define mpz_get_ull mpz_get_ui
#define mpz_set_ull mpz_set_ui
#endif

/* if.h */

void*xmalloc(size_t size);
void*xvalloc(size_t size);
void*xcalloc(size_t n,size_t s);
void*xrealloc(void*x,size_t size);
void complain(char*fmt,...);
#ifdef __ppc__
void Schlendrian(char*fmt,...);
#else
void Schlendrian(char*fmt,...) NAME("Schlendrian");
#endif
void logbook(int l,char*fmt,...);
int errprintf(char*fmt,...);
void adjust_bufsize(void**,size_t*,size_t,size_t,size_t);
extern int verbose;
extern FILE*logfile;
#if !defined( _MSC_VER ) && !defined(__MINGW32__) && defined( BIGENDIAN)
int write_i64(FILE*,i64_t*,size_t);
int write_u64(FILE*,u64_t*,size_t);
int write_i32(FILE*,i32_t*,size_t);
int write_u32(FILE*,u32_t*,size_t);
int read_i64(FILE*,i64_t*,size_t);
int read_u64(FILE*,u64_t*,size_t);
int read_i32(FILE*,i32_t*,size_t);
int read_u32(FILE*,u32_t*,size_t);
#else
#define write_i64(ofile,buffer,count) fwrite((void*)buffer,sizeof(i64_t),count,ofile)
#define write_u64(ofile,buffer,count) fwrite((void*)buffer,sizeof(u64_t),count,ofile)
#define write_u32(ofile,buffer,count) fwrite((void*)buffer,sizeof(u32_t),count,ofile)
#define write_i32(ofile,buffer,count) fwrite((void*)buffer,sizeof(i32_t),count,ofile)
#define read_i64(ofile,buffer,count) fread((void*)buffer,sizeof(i64_t),count,ofile)
#define read_u64(ofile,buffer,count) fread((void*)buffer,sizeof(u64_t),count,ofile)
#define read_u32(ofile,buffer,count) fread((void*)buffer,sizeof(u32_t),count,ofile)
#define read_i32(ofile,buffer,count) fread((void*)buffer,sizeof(i32_t),count,ofile)
#endif 
int yn_query(char*fmt,...);
ssize_t skip_blank_comments(char**,size_t*,FILE*);

#ifdef NEED_ASPRINTF
int asprintf(char**,const char*,...);
#endif

#ifdef NEED_GETLINE
ssize_t getline(char**,size_t*,FILE*);
#endif


/* input-poly.h */
void input_poly(mpz_t,mpz_t**,i32_t*,mpz_t**,i32_t*,mpz_t,FILE*);

/* lasieve-prepn.h */
void lasieve_setup(u32_t*,u32_t*,u32_t,i32_t,i32_t,i32_t,i32_t,u32_t*);

/* primgen32.h */
typedef struct{
  u32_t Pind;
  u32_t first_in_sieve;
  u32_t nPrim;
  u32_t Prime;
  unsigned char*PDiffs;
  u32_t use_private;
  u32_t PDiffs_allocated;
} pr32_struct;

void initprime32(pr32_struct*ps);
u32_t firstprime32(pr32_struct*ps);
u32_t nextprime32(pr32_struct*ps);
void clearprime32(pr32_struct*ps);
void p32_seek(pr32_struct*ps,u32_t lb);

/* real-poly-aux.h */
void tpol(double *rop,double *op,i32_t deg,i32_t x0,i32_t x1,i32_t y0,i32_t y1);
double rpol_eval(double *p,i32_t d,double x,double y);
double rpol_lb(double *pol,i32_t poldeg,double a,double b);


/* recurrence6.h */
void rec_info_init(u32_t A,u32_t ub);
u32_t get_recurrence_info(u32_t*res_ptr,u32_t p,u32_t r);

/* redu2.h */
void reduce2(i32_t*,i32_t*,i32_t*,i32_t*,i32_t,i32_t,i32_t,i32_t,float);
extern i32_t n_iter;

/* primgen32.c */
u32_t pr32_seek(pr32_struct * ps, u32_t lb);

/* These are supposed to be obtained by #define _GNU_SOURCE
   before including stdio.h, but that seems to break things,
   so here we go: 
*/
int asprintf(char **strp, const char *fmt, ...);
int vasprintf(char **strp, const char *fmt, va_list ap);


#endif
