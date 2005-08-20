/* if.h
  Hacked up for inclusion in GGNFS by Chris Monico.

  Copyright (C) 2001 Jens Franke.
  This file is part of gnfs4linux, distributed under the terms of the 
  GNU General Public Licence and WITHOUT ANY WARRANTY.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#ifndef _IF_H
#define _IF_H

#include <gmp.h>
#include <stdlib.h>
#include <stdarg.h> 
#include <stdio.h> 

#include "ggnfs.h"

void*xmalloc(size_t size);
void*xvalloc(size_t size);
void*xcalloc(size_t n,size_t s);
void*xrealloc(void*x,size_t size);
void complain(char*fmt,...);
void logbook(int l,char*fmt,...);
int errprintf(char*fmt,...);
void adjust_bufsize(void**,size_t*,size_t,size_t,size_t);
extern int verbose;
extern FILE*logfile;
void Schlendrian(char*fmt,...);
#if defined(GGNFS_BIGENDIAN)
int write_i64(FILE*,int64_t*,size_t);
int write_u64(FILE*,uint64_t*,size_t);
int write_i32(FILE*,int32_t*,size_t);
int write_u32(FILE*,uint32_t*,size_t);
int read_i64(FILE*,int64_t*,size_t);
int read_u64(FILE*,uint64_t*,size_t);
int read_i32(FILE*,int32_t*,size_t);
int read_u32(FILE*,uint32_t*,size_t);
#else
#define write_i64(ofile,buffer,count) fwrite((void*)buffer,sizeof(int64_t),count,ofile)
#define write_u64(ofile,buffer,count) fwrite((void*)buffer,sizeof(uint64_t),count,ofile)
#define write_u32(ofile,buffer,count) fwrite((void*)buffer,sizeof(uint32_t),count,ofile)
#define write_i32(ofile,buffer,count) fwrite((void*)buffer,sizeof(int32_t),count,ofile)
#define read_i64(ofile,buffer,count) fread((void*)buffer,sizeof(int64_t),count,ofile)
#define read_u64(ofile,buffer,count) fread((void*)buffer,sizeof(uint64_t),count,ofile)
#define read_u32(ofile,buffer,count) fread((void*)buffer,sizeof(uint32_t),count,ofile)
#define read_i32(ofile,buffer,count) fread((void*)buffer,sizeof(int32_t),count,ofile)
#endif 
int yn_query(char*fmt,...);
int skip_blank_comments(char**,size_t*,FILE*);

#ifdef NEED_ASPRINTF
int asprintf(char**,const char*,...);
#endif

//#ifdef NEED_GETLINE
size_t getline(char**,size_t*,FILE*);
//#endif


/* These are supposed to be obtained by #define _GNU_SOURCE
   before including stdio.h, but that seems to break things,
   so here we go: 
*/
int asprintf(char **strp, const char *fmt, ...);
int vasprintf(char **strp, const char *fmt, va_list ap);


/* gmp-aux.h */
void adjust_mpz_bufsize(mpz_t**x,size_t*alloc_ptr,size_t size,size_t increment);
int string2mpz(mpz_t rop,char*x,int base);

#if __GNU_MP_VERSION < 3 || __GNU_MP_VERSION_MINOR == 0
#define NEED_MPZ_MUL_SI
void mpz_mul_si(mpz_t x,mpz_t y,long int z);
#endif

#define mpz_add_si(r,o,s) \
  ({ int _s; _s= (s); \
     _s>=0 ? mpz_add_ui(r,o,(unsigned long)(_s)) : mpz_sub_ui(r,o,(unsigned long)(-_s)); })

#ifdef ULL_NO_UL
void mpz_ull_init();
void mpz_mul_ull(mpz_t rop,mpz_t op1,uint64_t op2);
int mpz_fits_sllong_p(mpz_t);
int mpz_fits_ullong_p(mpz_t);
int64_t mpz_get_sll(mpz_t);
void mpz_set_sll(mpz_t,int64_t);
uint64_t mpz_get_ull(mpz_t);
void mpz_set_ull(mpz_t,uint64_t);
int mpz_cmp_ull(mpz_t,uint64_t);
#else
#define mpz_ull_init()
#define mpz_mul_ull mpz_mul_ui
#define mpz_fits_sllong_p mpz_fits_slong_p
#define mpz_fits_ullong_p mpz_fits_ulong_p
#define mpz_get_sll mpz_get_si
#define mpz_set_sll mpz_set_si
#define mpz_get_ull mpz_get_ui
#define mpz_set_ull mpz_set_ui
#define mpz_cmp_ull mpz_cmp_ui
#endif

void numread(char *arg, unsigned *x);

#endif
