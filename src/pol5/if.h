/* if.h

   Copyright (C) 2005 T. Kleinjung.
   This file is part of pol5, distributed under the terms of the
   GNU General Public Licence and WITHOUT ANY WARRANTY.
                                                                                                               
  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  2111-1307, USA.
                                                                                                               
  CJM, 2/18/05: Hacked up for inclusion in GGNFS. It was originally written by
  T. Kleinjung and/or Jens Franke. This is ugly - I hacked this into a general
  purpose header file with prototypes for things commonly used across several
  files. It should really be renamed to something like pol.h probably.
*/

#include <unistd.h>
#include <stdarg.h> 
#include <stdio.h> 
#include <gmp.h>
#if !defined(__FreeBSD__) && !defined(__OpenBSD__) && !defined(__NetBSD__)
#include <malloc.h>
#endif
#include "defs.h"

#define BITS_PER_ULONG 32
/*
#define valloc(n)    _aligned_malloc((n), 4096)
#define vfree(p)     _aligned_free(p)
*/
#define valloc(n) malloc(n)
#define vfree(n)  free(n)

void *xmalloc(size_t size);
void *xvalloc(size_t size);
void *xcalloc(size_t n, size_t s);
void *xrealloc(void *x, size_t size);
void xfscanf(FILE *file, char *fmt_string,...);
FILE *xfopen(char *Fn, char *mode);
void xfclose(FILE *file);
int xfgetc(FILE *str);
void complain(char *fmt,...);
void Schlendrian(char *fmt,...);
void logbook(int l, char *fmt,...);
void numread(char *arg, unsigned *x);
void mpz_set_ull(mpz_t targ, unsigned long long src);
void mpz_set_sll(mpz_t targ, long long src);


ssize_t getline(char **lineptr, size_t * n, FILE * stream);
int asprintf(char **ptr, const char *template, ...);

/* assess.c */
uint invert(uint a, uint p);
void murphy_en(double *me, int deg0, double *dbl_coeff0, int deg1, double *dbl_coeff1, 
               double alpha0, double alpha1, double skewness, int nsm);
void murphy_e(double *me, int deg0, double *dbl_coeff0, int deg1, double *dbl_coeff1, 
              double alpha0, double alpha1, double skewness);
int compute_alpha(double *alpha, int deg, uint **coeffmod, mpz_t *gmp_coeff, double alpha_targ);
void compute_alpha_exact(double *alpha, int deg, uint **coeffmod, mpz_t *gmp_coeff, uint pb);
void init_assess(double b0, double b1, double area, int pb);


/* roots.c */
int find_optima(int *deg, double **coeff, double skewness, double **optima);

/* primes.c */
void prime_table_init();
uint get_next_prime();

