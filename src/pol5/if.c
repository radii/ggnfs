/* if.c

   Copyright (C) 2005 T. Kleinjung.
   This file is part of pol5, distributed under the terms of the
   GNU General Public Licence and WITHOUT ANY WARRANTY.
                                                                                                               
  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  2111-1307, USA.
                                                                                                               
  CJM, 2/18/05: Hacked up for inclusion in GGNFS. It was originally written by
  T. Kleinjung and/or Jens Franke.
*/

#include <unistd.h>
#include <time.h> 
#if !(defined(__MINGW32__)) && !(defined(_MSC_VER))
#include <sys/times.h> 
#endif
#include <stdio.h> 
#include <stdlib.h> 
#include <stdarg.h> 
#include <string.h> 
#include "gmp.h" 
#include <limits.h> 

#if defined( __CYGWIN__ ) || defined( _MSC_VER ) || defined(__MINGW32__) || defined(__FreeBSD__) || defined(__OpenBSD__) || defined(__NetBSD__)
#define NEED_GETLINE
#define NEED_ASPRINTF
#define valloc(n) malloc(n)
#define ulong unsigned long
#endif

#include "if.h"

static uint used_cols, ncol=80;
static uint log_clock=0;



#define ErrEx(s) { perror(s); exit(-1); }

#ifdef NEED_ASPRINTF
/****************************************************/
int vasprintf(char **ptr, const char *template, va_list ap)
/****************************************************/
{
  int n, size = 32;

  while (1) {
    *ptr = xmalloc(size);

    n = vsnprintf(*ptr, size, template, ap);
    if (1 + strlen(*ptr) < size)
      return n;
    size *= 2;
    free(*ptr);
  }
}

/****************************************************/
int asprintf(char **ptr, const char *template, ...)
/****************************************************/
{
  int rv;
  va_list ap;
                                                                                                               
  va_start(ap, template);
  rv = vasprintf(ptr, template, ap);
  va_end(ap);
  return rv;
}
#endif

#ifdef NEED_GETLINE
#define GETL_INCR 128
/****************************************************/
ssize_t getline(char **lineptr, size_t * n, FILE * stream)
/****************************************************/
{
  int rv;
                                                                                                 
  if (*n == 0) {
    *n = GETL_INCR;
    *lineptr = xmalloc(*n);
  }
                                                                                                 
  rv = 0;
  for (;;) {
    int m;
                                                                                                 
    m = *n - rv;
    if (fgets(*lineptr + rv, m - 1, stream) == NULL)
      break;
    rv = strlen(*lineptr);
    if (rv == 0 || (*lineptr)[rv - 1] == '\n')
      break;
    *n += GETL_INCR;
    *lineptr = xrealloc(*lineptr, *n);
  }
  return rv;
}
#endif


void *xmalloc(size_t size)
{
  char *x;
  if ((x=malloc(size))==NULL) ErrEx("xmalloc: ")
  else return x;
}


void *xvalloc(size_t size)
{
  char *x;
  if ((x=valloc(size))==NULL) ErrEx("xmalloc: ")
  else return x;
}


void *xcalloc(size_t n, size_t s)
{
  void *p;
  if ((p=calloc(n,s))!=NULL) return p;
  ErrEx("calloc");
}


void *xrealloc(void *x, size_t size)
{
  char *y;
  if ((y=realloc(x,size))==NULL) ErrEx("xrealloc: ")
  else return y;
}


FILE *xfopen (char *Fn, char *mode)
{
  FILE *x;
  if ((x=fopen(Fn,mode))==NULL) ErrEx("Open: ")
  return x;
}


void xfclose(FILE *file)
{
  if (fclose(file)) ErrEx("fclose")
}


int xfgetc(FILE *str)
{
  int i;
  if ((i=fgetc(str))!=EOF) return i;
  ErrEx("fgetc: ");
}


void complain(char *fmt,...)
{
  va_list arglist;
  va_start(arglist,fmt);
  vfprintf(stderr,fmt,arglist);
  exit(1);
}


void Schlendrian(char *fmt,...)
{
  va_list arglist;
  va_start(arglist,fmt);
  vfprintf(stderr,fmt,arglist);
  abort();
}


void numread(char *arg, unsigned *x)
{
  char*fmt;
  int offs;
  if (strlen(arg)==0) {
    *x=0;
    return;
  }
  if (*arg=='0') {
    if (strlen(arg)==1) {
      *x=0;
      return;
    }
    if (*(arg+1)=='x') {
      fmt="%X";
      offs=2;
    } else {
      fmt="%o";
      offs=1;
    }
  } else {
    offs=0;
    fmt="%u";
  }
  if (sscanf(arg+offs,fmt,x)!=1) {
    fprintf(stderr,"Bad integer value for some option!\n");
    /*Usage();*/ /* !!!!!!!!!!!!!!! */
  }
}


#define NEW_NUMBER -0x10000
#if 0
void logbook(int l, char *fmt,...)
{
  if (l==NEW_NUMBER) {
    used_cols=0;
    return;
  }
  if (l<verbose) {
    va_list arglist;
    char *output_str;
    uint sl;
    va_start(arglist,fmt);
    vasprintf(&output_str,fmt,arglist);
    sl= strlen(output_str);
    if(used_cols+sl>ncol) {
      fprintf(stderr,"\n");
      used_cols=0;
    }
    fputs(output_str,stderr);
    if (output_str[sl-1]=='\n') used_cols=0;
    else used_cols+=sl;
    free(output_str);
  }
}
#endif

void mpz_set_ull(mpz_t targ, unsigned long long src)
{
  mpz_set_ui(targ,
(ulong)((src&((unsigned long long)ULONG_MAX<<BITS_PER_ULONG))>>BITS_PER_ULONG));
  mpz_mul_2exp(targ,targ,BITS_PER_ULONG);
  mpz_add_ui(targ,targ,(ulong)(src&ULONG_MAX));
}


void mpz_get_ull(ull *targptr, mpz_t src)  /* preliminary */
{
  unsigned long long res;
  char *str;

  str=mpz_get_str(NULL,10,src);
  res=strtoull(str,NULL,10);
  free(str);
  *targptr=res;
}


void mpz_set_sll(mpz_t targ, long long src)
{
  if (src>=0) mpz_set_ull(targ,src);
  else { mpz_set_ull(targ,(ull)(-src)); mpz_neg(targ,targ); }
}

#if 0
void
mpz_mul_si(mpz_t x, mpz_t y, int i)
{
  uint m;
  if(i<0) {
    m=-i;
    mpz_mul_ui(x,y,m);
    mpz_neg(x,x);
  } else {
    m=i;
    mpz_mul_ui(x,y,m);
  }
}
#endif

