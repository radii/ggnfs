#include <sys/types.h>
#include <limits.h>

#include "siever-config.h"
#include "32bit.h"
#include "../if.h"

/*
 * Auxilliary information for trial division using MMX instructions,
 * stored in the form (root0,root1,root2,root3,prime0,
 *                     prime0,prime1,prime2,prime3,
 *                     mi0,mi1,mi2,mi3),
 * where the mi are modular Inverses to the primes.
 */

static u16_t *(MMX_TdAux[2]),*(MMX_TdBound[2]);

/*
 * The projective roots used to update MMX_TdAux when the line is changed.
 */

static u16_t **(MMX_TdPr[2]);

static size_t MMX_TdAlloc[2]={0,0};

static int jps=0;
static int MMX_TdNl[2];

/* Read-Ahead safety. */
#define RAS 4

void
MMX_TdAllocate(int jps_arg,size_t s0,size_t s1)
{
  int side;

  MMX_TdAlloc[0]=4*((s0+3)/4);
  MMX_TdAlloc[1]=4*((s1+3)/4);
  jps=jps_arg;

  for(side=0;side<2;side++) {
    int i;

    MMX_TdAux[side]=xmalloc(3*(MMX_TdAlloc[side]+RAS)*sizeof(*MMX_TdAux));
    MMX_TdPr[side]=xmalloc(jps*sizeof(*(MMX_TdPr[side])));
    for(i=0;i<jps;i++)
      MMX_TdPr[side][i]=xmalloc((MMX_TdAlloc[side]+RAS)*sizeof(**(MMX_TdPr[side])));
  }
}

u16_t*
MMX_TdInit(int side,u16_t *x,u16_t *x_ub,u32_t *pbound_ptr,
	   int initialize)
{
  u16_t *y,*z,*u;
  u32_t p_bound;

  if(initialize==1) {
    int i;
    if(x_ub>x+4*MMX_TdAlloc[side])
      Schlendrian("Buffer overflow in MMX_TdInit\n");
    z=MMX_TdAux[side]+4;
    y=x;
    i=0;
    while(y+16<x_ub) {
      int j;
      for(j=0;j<4;j++,y+=4,z++,i++) {
	u16_t mi;	u32_t r,rr;
	int k;

	modulo32=*y;
	mi=*y;
	*z=mi;
	mi=2*mi-mi*mi*modulo32;
	mi=2*mi-mi*mi*modulo32;
	mi=2*mi-mi*mi*modulo32;
	*(z+4)=mi;
	r=y[1];
	rr=r;
	for(k=0;k<jps;k++) {
	  MMX_TdPr[side][k][i]=rr;
	  rr=modadd32(rr,r);
	}
      }
      z+=8;
    }
  }
  z=MMX_TdAux[side];
  u=MMX_TdPr[side][jps-1];
  p_bound=*pbound_ptr;
  for(y=x;y<x_ub-16;y=y+16,z+=12,u+=4) {
    if(z[4]>p_bound) break;
    modulo32=y[0];
    z[0]=modsub32(0,y[3]);
    y[3]=modadd32(y[3],u[0]);
    modulo32=y[4];
    z[1]=modsub32(0,y[7]);
    y[7]=modadd32(y[7],u[1]);
    modulo32=y[8];
    z[2]=modsub32(0,y[11]);
    y[11]=modadd32(y[11],u[2]);
    modulo32=y[12];
    z[3]=modsub32(0,y[15]);
    y[15]=modadd32(y[15],u[3]);
  }
  *pbound_ptr=z[-5];
  MMX_TdBound[side]=z;
  MMX_TdNl[side]=(z-MMX_TdAux[side])/12;
  return y;
}

void asm_TdUpdate(u16_t*,size_t,u16_t*);

void
MMX_TdUpdate(int side,int j_step)
{
#if 0
  u16_t *x,*y;

  y=MMX_TdPr[side][j_step-1];
  for(x=MMX_TdAux[side];x<MMX_TdBound[side];x+=12,y+=4) {
    modulo32=x[4];
    x[0]=modsub32(x[0],y[0]);
    modulo32=x[5];
    x[1]=modsub32(x[1],y[1]);
    modulo32=x[6];
    x[2]=modsub32(x[2],y[2]);
    modulo32=x[7];
    x[3]=modsub32(x[3],y[3]);
  }
#else
  asm_TdUpdate(MMX_TdAux[side],
	       (MMX_TdBound[side]-MMX_TdAux[side])/4,MMX_TdPr[side][j_step-1]);
#endif
}

u32_t *asm_MMX_Td(u32_t*,u64_t,u16_t*,u16_t*);

u64_t MMX_TdNloop=0;

u32_t *
MMX_Td(u32_t *pbuf,int side,u16_t strip_i)
{
#if 1
  u64_t Strip_i=strip_i;
  MMX_TdNloop+=(MMX_TdBound[side]-MMX_TdAux[side])/4;
  return asm_MMX_Td(pbuf,Strip_i|(Strip_i<<16)|(Strip_i<<32)|(Strip_i<<48),
		    MMX_TdAux[side],MMX_TdNl[side]);
#else
  u16_t* x;

  for(x=MMX_TdAux[side];x<MMX_TdBound[side];x+=8) {
    int i;

    for(i=0;i<4;i++,x++) {
      u16_t t;

      modulo32=x[4];
    
      t=strip_i+x[0];
      t*=x[8];
      if(((modulo32*(u32_t)t)&0xffff0000)==0)
	*(pbuf++)=modulo32;
    }
  }
  return pbuf;
#endif
}
