/* if.c
  By Jens Franke.
  6/13/04: Hacked up for use in GGNFS by Chris Monico.
                                                                                                                                                                                                             
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
                                                                                                                                                                                                             
  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#ifdef _MSC_VER
#pragma warning (disable: 4996) /* warning C4996: 'function' was declared deprecated */
#endif

#include "if.h"

#include <time.h>

#if defined (__MINGW32__) || defined (MINGW32)
#include <io.h>
#elif defined(_MSC_VER)
#include <io.h>
#include <unistd.h>
#else
#include <unistd.h>
#include <sys/times.h>
#endif

#include <string.h>
#include <gmp.h>
#include <limits.h>

int verbose = 0;
static size_t used_cols, ncol = 80;

FILE *logfile = NULL;
void complain(char * fmt, ...);

/****************************************************/
void *xmalloc(size_t size)
/****************************************************/
{
  char *x;

  if (size == 0)
    return NULL;
  if ((x = malloc(size)) == NULL)
    complain("xmalloc: %m");
  return x;
}

/****************************************************/
void *xvalloc(size_t size)
/****************************************************/
{
  char *x;

  if (size == 0)
    return NULL;
  if ((x = malloc(size)) == NULL)
    complain("xvalloc: %m");
  memset(x,0,size);
  return x;
}

/****************************************************/
void *xcalloc(size_t n, size_t s)
/****************************************************/
{
  void *p;

  if (n == 0 || s == 0)
    return NULL;
  if ((p = calloc(n, s)) == NULL)
    complain("calloc: %m");
  return p;
}

/****************************************************/
void *xrealloc(void *x, size_t size)
/****************************************************/
{
  char *y;

  if (size == 0) {
    if (x != NULL)
      free(x);
    return NULL;
  }
  if ((y = realloc(x, size)) == NULL && size != 0)
    complain("xrealloc: %m");
  return y;
}


/****************************************************/
void complain(char *fmt, ...)
/****************************************************/
{
  va_list arglist;

  va_start(arglist, fmt);
  vfprintf(stderr, fmt, arglist);
  if (logfile != NULL)
    vfprintf(logfile, fmt, arglist);
  exit(1);
}

/*****************************************************/
void Schlendrian(char *fmt, ...)
/****************************************************/
{
  va_list arglist;
  va_start(arglist, fmt);
  vfprintf(stderr, fmt, arglist);
  if (logfile != NULL)
    vfprintf(logfile, fmt, arglist);
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
/****************************************************/
void logbook(int l, char *fmt, ...)
/****************************************************/
{
  if (l == NEW_NUMBER) {
    used_cols = 0;
    return;
  }
  if (l < verbose) {
    va_list arglist;
    char *output_str;
    size_t sl;

    va_start(arglist, fmt);
    vasprintf(&output_str, fmt, arglist);
    sl = strlen(output_str);
    if (used_cols + sl > ncol) {
      fprintf(stderr, "\n");
      if (logfile != NULL)
        fprintf(logfile, "\n");
      used_cols = 0;
    }
    fputs(output_str, stderr);
    if (logfile != NULL)
      fputs(output_str, logfile);
    if (output_str[sl - 1] == '\n')
      used_cols = 0;
    else
      used_cols += sl;
    free(output_str);
  }
}

/****************************************************/
int errprintf(char *fmt, ...)
/****************************************************/
{
  va_list arglist;
  int res;

  va_start(arglist, fmt);
  if (logfile != NULL)
    vfprintf(logfile, fmt, arglist);
  res = vfprintf(stderr, fmt, arglist);
  return res;
}

/****************************************************/
void adjust_bufsize(void **buf, size_t * alloc, size_t req,
                    size_t incr, size_t item_size)
/****************************************************/
{
  if (req > *alloc) {
    size_t new_alloc;

    new_alloc = *alloc + incr * ((req + incr - 1 - *alloc) / incr);
    if (*alloc > 0)
      *buf = xrealloc(*buf, new_alloc * item_size);
    else
      *buf = xmalloc(new_alloc * item_size);
    *alloc = new_alloc;
  }
}

/****************************************************/
int yn_query(char *fmt, ...)
/****************************************************/
{
  va_list arglist;
  char answer[10];

  va_start(arglist, fmt);
  if (logfile != NULL)
    vfprintf(logfile, fmt, arglist);
  vfprintf(stderr, fmt, arglist);
#ifdef _MSC_VER
  if(_isatty(_fileno(stdin)) || _isatty(_fileno(stderr)))
#else
  if (!isatty(STDIN_FILENO) || !isatty(STDERR_FILENO))
#endif
	  return 0;

  fflush(stderr);
  while (scanf("%9s", answer) != 1 ||
         (strcasecmp(answer, "yes") != 0 && strcasecmp(answer, "no") != 0)) {
    fprintf(stderr, "Please answer yes or no!\n");
    vfprintf(stderr, fmt, arglist);
  }
  if (strcasecmp(answer, "yes") == 0)
    return 1;
  return 0;
}

/****************************************************/
int skip_blanks_comments(char **iline, size_t * iline_alloc, FILE * ifi)
/****************************************************/
{
  while (getline(iline, iline_alloc, ifi) > 0) {
    if (**iline != '#' && strspn(*iline, "\n\t ") < strlen(*iline))
      return 1;
  }
  return 0;
}

#ifdef GGNFS_BIGENDIAN

/****************************************************/
static u32_t bswap_32(u32_t x)
/****************************************************/
{
  return ((x & 0x000000ffUL) << 24) | ((x & 0x0000ff00UL) << 8) |
    ((x & 0x00ff0000UL) >> 8) | ((x & 0xff000000UL) >> 24);
}

/****************************************************/
static u64_t bswap_64(u64_t x)
/****************************************************/
{
  return ((x & 0xffULL) << 56) | ((x & 0xff00ULL) << 40) | ((x & 0xff0000ULL)
                                                            << 24) | ((x &
                                                                       0xff000000ULL)
                                                                      << 8) |
    ((x & 0xff00000000ULL) >> 8) | ((x & 0xff0000000000ULL) >> 24) |
    ((x & 0xff000000000000ULL) >> 40) | ((x & 0xff00000000000000ULL) >> 56);
}

/****************************************************/
int write_i64(FILE * ofile, i64_t * buffer, size_t count)
/****************************************************/
{
  size_t i;
  int res;

  for (i = 0; i < count; i++)
    buffer[i] = bswap_64(buffer[i]);
  res = fwrite(buffer, sizeof(*buffer), count, ofile);
  for (i = 0; i < count; i++)
    buffer[i] = bswap_64(buffer[i]);
  return res;
}

/****************************************************/
int write_u64(FILE * ofile, u64_t * buffer, size_t count)
/****************************************************/
{
  size_t i;
  int res;

  for (i = 0; i < count; i++)
    buffer[i] = bswap_64(buffer[i]);
  res = fwrite(buffer, sizeof(*buffer), count, ofile);
  for (i = 0; i < count; i++)
    buffer[i] = bswap_64(buffer[i]);
  return res;
}

/****************************************************/
int write_i32(FILE * ofile, i32_t * buffer, size_t count)
/****************************************************/
{
  size_t i;
  int res;

  for (i = 0; i < count; i++)
    buffer[i] = bswap_32(buffer[i]);
  res = fwrite(buffer, sizeof(*buffer), count, ofile);
  for (i = 0; i < count; i++)
    buffer[i] = bswap_32(buffer[i]);
  return res;
}

/****************************************************/
int write_u32(FILE * ofile, u32_t * buffer, size_t count)
/****************************************************/
{
  size_t i;
  int res;

  for (i = 0; i < count; i++)
    buffer[i] = bswap_32(buffer[i]);
  res = fwrite(buffer, sizeof(*buffer), count, ofile);
  for (i = 0; i < count; i++)
    buffer[i] = bswap_32(buffer[i]);
  return res;
}

/****************************************************/
int read_i64(FILE * ifile, i64_t * buffer, size_t count)
/****************************************************/
{
  size_t i;
  int res;

  res = fread(buffer, sizeof(*buffer), count, ifile);
  for (i = 0; i < count; i++)
    buffer[i] = bswap_64(buffer[i]);
  return res;
}

/****************************************************/
int read_u64(FILE * ifile, u64_t * buffer, size_t count)
/****************************************************/
{
  size_t i;
  int res;

  res = fread(buffer, sizeof(*buffer), count, ifile);
  for (i = 0; i < count; i++)
    buffer[i] = bswap_64(buffer[i]);
  return res;
}

/****************************************************/
int read_i32(FILE * ifile, i32_t * buffer, size_t count)
/****************************************************/
{
  size_t i;
  int res;

  res = fread(buffer, sizeof(*buffer), count, ifile);
  for (i = 0; i < count; i++)
    buffer[i] = bswap_32(buffer[i]);
  return res;
}

/****************************************************/
int read_u32(FILE * ifile, u32_t * buffer, size_t count)
/****************************************************/
{
  size_t i;
  int res;

  res = fread(buffer, sizeof(*buffer), count, ifile);
  for (i = 0; i < count; i++)
    buffer[i] = bswap_32(buffer[i]);
  return res;
}

#endif

#ifdef GGNFS_GNU_MISSING
#define GETL_INCR 128
/****************************************************/
ssize_t getline(char **lineptr, size_t * n, FILE * stream)
/****************************************************/
{
  ssize_t rv = 0;

  if (*n == 0) {
    *n = GETL_INCR;
    *lineptr = xmalloc(*n);
  }

  for (;;) {
    int m;

    m = (int)(*n - rv);
    if (fgets(*lineptr + rv, m - 1, stream) == NULL)
      break;
    rv = (ssize_t)strlen(*lineptr);
    if (rv == 0 || (*lineptr)[rv - 1] == '\n')
      break;
    *n += GETL_INCR;
    *lineptr = xrealloc(*lineptr, *n);
  }
  return rv;
}

/****************************************************/
int vasprintf(char **ptr, const char *template, va_list ap)
/****************************************************/
{
  int n;
  size_t size = 32;

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
