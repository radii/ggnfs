/* primes.c

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

#include <fcntl.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "if.h"

#define SMALLPRIMELIMIT   65536
#define NSMALLPRIMES      6542
#define TABLELEN          32768         /* SMALLPRIMELIMIT=2*TABLELEN */

#define uchar   unsigned char

ushort table_small_primes[NSMALLPRIMES];
ushort *table_small_ptr_end, *prime_small_ptr;
uchar table_primes[TABLELEN];
uchar *table_ptr_end, *prime_ptr;
uint table_displacement;


void prime_table_init()
{
  ushort pr[]={ 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };
  uchar tmp[SMALLPRIMELIMIT];
  int ind, i, p, x;

  table_small_ptr_end=table_small_primes+NSMALLPRIMES;
  table_ptr_end=table_primes+TABLELEN;
  
  for (ind=0; ind<11; ind++) table_small_primes[ind]=pr[ind];
  memset(tmp,1,SMALLPRIMELIMIT);
  tmp[1]=0;
  for (i=0; i<11; i++) {
    p=(int)(pr[i]); x=0;
    while (x<SMALLPRIMELIMIT) { tmp[x]=0; x+=p; }
  }
  for (i=32; i<256; i++)
    if (tmp[i]) table_small_primes[ind++]=i;
  for (i=11; i<ind; i++) {
    p=(int)(table_small_primes[i]); x=p;
    while (x<SMALLPRIMELIMIT) { tmp[x]=0; x+=p; }
  }
  for (i=256; i<SMALLPRIMELIMIT; i++)
    if (tmp[i]) {
      table_small_primes[ind++]=i;
      if (ind>NSMALLPRIMES) Schlendrian("prime_table_init");
    }
  if (ind!=NSMALLPRIMES) Schlendrian("prime_table_init !=");
  prime_small_ptr=table_small_primes-1;
  prime_ptr=NULL;
}


void prime_table_zero()
{
  prime_small_ptr=table_small_primes-1;
  prime_ptr=NULL;
}


void compute_prime_table()
{
  uint i, p, x;

  table_displacement+=2*TABLELEN;
  memset(table_primes,1,TABLELEN);
  for (i=1; i<NSMALLPRIMES; i++) {
    p=(uint)(table_small_primes[i]);
    x=(table_displacement+1)%p;
    if (x) x=p-x;
    if (x&1) x+=p; x/=2;
    while (x<TABLELEN) { table_primes[x]=0; x+=p; }
  }
  prime_ptr=table_primes;
}


void set_prime_table(uint p)
{
  if (p<=SMALLPRIMELIMIT) {
    prime_ptr=NULL;
    prime_small_ptr=table_small_primes;
    while (prime_small_ptr<table_small_ptr_end) {
      if ((uint)(*prime_small_ptr)>=p) {
        prime_small_ptr--;
        return;
      }
      prime_small_ptr++;
    }
    return;
  } 
  table_displacement=p-2*TABLELEN;
  if (table_displacement&1) table_displacement--;
  compute_prime_table();
  prime_ptr--;
}


uint get_next_prime()
{
  if (prime_ptr==NULL) {
    prime_small_ptr++;
    if (prime_small_ptr<table_small_ptr_end)
      return ((uint)(*prime_small_ptr));
    table_displacement=0;
    compute_prime_table();
    prime_ptr--;
  }
  prime_ptr++;
  while (1) {
    while (prime_ptr<table_ptr_end) {
      if (*prime_ptr)
        return (1+table_displacement+2*(prime_ptr-table_primes));
      prime_ptr++;
    }
    compute_prime_table();
  }
/* never reached */
  return 0;
}

