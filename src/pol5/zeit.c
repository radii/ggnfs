/* zeit.c

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

#ifdef HAVE_ASM_INTEL
#define HAVE_ASM
#endif

#ifdef HAVE_ASM_ALPHA
#define HAVE_ASM
#endif


#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"

clock_t *zeitcounter;
ull *asmzeitcounter;
double *zeitsum;
uint zeitcounteranz;
#ifdef HAVE_ASM_ALPHA
uint *asmzeitcounter_begin;
#endif

extern void *xmalloc(size_t size);

#ifdef HAVE_ASM_INTEL
static inline void asmgetclock(ull *clptr)
{
#ifdef _MSC_VER
	__asm
	{
		rdtsc
		mov		ecx,clptr
		mov		[ecx],eax
		mov		[ecx+4],edx
	}
#else
  asm("rdtsc\n"
      "movl %%eax,(%%ecx)\n"
      "movl %%edx,4(%%ecx)" : : "c" (clptr) :
      "%eax", "%edx");
#endif
}
#elif defined HAVE_ASM_ALPHA
static void asmgetclock(ull *clptr)
{
  asm_getclock(clptr);
}

#endif

void zeita(int i)
{ 

#ifdef HAVE_ASM_INTEL
  ull asmcl;

  zeitcounter[i]=clock();
  asmgetclock(&asmcl);
  asmzeitcounter[i]-=asmcl;
#elif defined HAVE_ASM_ALPHA
  ull asmcl;

  zeitcounter[i]=clock();
  asm_getclock(&asmcl);
  asmzeitcounter_begin[i]=(uint)asmcl;
#endif
}

void zeitA(int i)
{
#ifdef HAVE_ASM_INTEL
  ull asmcl;

  asmgetclock(&asmcl);
  asmzeitcounter[i]-=asmcl;
#elif defined HAVE_ASM_ALPHA
  ull asmcl;

  asm_getclock(&asmcl);
  asmzeitcounter_begin[i]=(uint)asmcl;
#endif
}

void zeitb(int i)
{
#if defined(HAVE_ASM_INTEL) || defined(HAVE_ASM_ALPHA)
  ull asmcl;
#endif

  zeitsum[i]+=(double)(clock()-zeitcounter[i]);
#ifdef HAVE_ASM_INTEL
  asmgetclock(&asmcl);
  asmzeitcounter[i]+=asmcl;
#elif defined HAVE_ASM_ALPHA
  asm_getclock(&asmcl);
  asmzeitcounter[i]+=(ull)((uint)(asmcl)-asmzeitcounter_begin[i]);
#endif
}

void zeitB(int i)
{
#ifdef HAVE_ASM_INTEL
  ull asmcl;

  asmgetclock(&asmcl);
  asmzeitcounter[i]+=asmcl;
#elif defined HAVE_ASM_ALPHA
  ull asmcl;

  asm_getclock(&asmcl);
  asmzeitcounter[i]+=(ull)((uint)(asmcl)-asmzeitcounter_begin[i]);
#endif
}

void initzeit(int i)
{
  zeitcounter=(clock_t *)xmalloc(i*sizeof(clock_t));
  memset(zeitcounter,0,i*sizeof(clock_t));
  zeitsum=(double *)xmalloc(i*sizeof(double));
  memset(zeitsum,0,i*sizeof(double));
  asmzeitcounter=(ull *)xmalloc(i*sizeof(ull));
  memset(asmzeitcounter,0,i*sizeof(ull));
#ifdef HAVE_ASM_ALPHA
  asmzeitcounter_begin=(uint *)xmalloc(i*sizeof(uint));
#endif
  zeitcounteranz=i;
}

void printzeit(int i)
{
  printf("%.2fs (%llu)  ",zeitsum[i]/CLOCKS_PER_SEC,asmzeitcounter[i]);
}

void printzeit2(int i1, int i2)
{
  printf("%.2fs (%llu)  ",(zeitsum[i1]+zeitsum[i2])/CLOCKS_PER_SEC,asmzeitcounter[i1]+asmzeitcounter[i2]);
}

void printzeitall()
{
  int i;

  printf("\nZeit: ");
  for (i=0; i<zeitcounteranz; i++)
    printf("%d: %.3fs  ",i,zeitsum[i]/CLOCKS_PER_SEC);
  printf("\n");
  printf("\nCycles: ");
  for (i=0; i<zeitcounteranz; i++)
    printf("%d: %llu  ",i,asmzeitcounter[i]);
  printf("\n");
}

