/**************************************************************/
/* blanczos64.c                                               */
/* Copyright 2004, Chris Monico.                              */
/**************************************************************/
/*  This file is part of GGNFS.
*
*   GGNFS is free software; you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation; either version 2 of the License, or
*   (at your option) any later version.
*
*   GGNFS is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with GGNFS; if not, write to the Free Software
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include "prand.h"
#include "ggnfs.h"

/*************************************************************/ 
/* Implementation of the Block Lanczos algorithm for finding */
/* vectors in the kernel of a large sparse matrix.           */
/* The implementation here is as descibed in:                */
/* ``A Block Lanczos Algorithm for Finding Dependencies over */
/*   GF(2)'', Peter Montgomery.                              */
/* The paper I have does not have a publication source on it,*/
/* so it's probably a preprint.                              */
/*************************************************************/ 


#ifdef __GNUC__
#define ALIGNED16 __attribute__ ((aligned (16)))
#else
#define ALIGNED16
#endif

#if defined(__MINGW32__) || defined(__CYGWIN__)
#define ASM_UL "_"
#else
#define ASM_UL ""
#endif

#if defined( __GNUC__ )
#define	USE_MMX_GCC
#elif defined( _MSC_VER )
#define	USE_MMX_MSC
#endif

#if defined( USE_MMX_GCC )
#if defined(__MINGW32__) || defined(MINGW32)
#define malloc_aligned64(p,a,n) (!((p) = (u64*)__mingw_aligned_malloc((n)*sizeof(u64),(a))) || ((uintptr_t)(p)%(a)))
#define free_aligned64(x) __mingw_aligned_free(x)
#else
#define malloc_aligned64(p,a,n) (!((p) = (u64*)memalign((a),(n)*sizeof(u64))) || ((uintptr_t)(p)%(a)))
#define free_aligned64(x) free(x)
#endif
#define XEMMS asm("emms")


#define XOR64(_a, _b) {       asm("movq (%0), %%mm0\n" \
                                  "pxor (%1), %%mm0\n" \
                                  "movq %%mm0, (%0)\n" \
                                  : : "r"(_a), "r"(_b) : "memory" ); } 
#elif defined( USE_MMX_MSC )
#define	malloc_aligned64(p,a,n)  (!((p) = (u64*)_aligned_malloc((n) * sizeof(u64), (a))) || ((uintptr_t)(p)%(a)))
#define free_aligned64(x)      _aligned_free(x)
#define	XEMMS                  __asm	emms
#define XOR64(_a, _b)          __asm {				\
                               __asm	mov	 eax,_a		\
                               __asm	mov	 edx,_b		\
                               __asm	movq mm0,[eax]		\
                               __asm	pxor mm0,[edx]		\
                               __asm	movq [eax],mm0	}
#else
#define	malloc_aligned64(p,a,n)	(!((p) = (u64*)malloc((n) * sizeof(u64))) || ((uintptr_t)(p)%(a)))
#define free_aligned64(x)	        free(x)
#define	XEMMS
#define XOR64(_a, _b)           { *(_a) ^= *(_b); }
#endif

static const u64 bit64[] ALIGNED16 ={
0x0000000000000001ULL,0x0000000000000002ULL,0x0000000000000004ULL,0x0000000000000008ULL,
0x0000000000000010ULL,0x0000000000000020ULL,0x0000000000000040ULL,0x0000000000000080ULL,
0x0000000000000100ULL,0x0000000000000200ULL,0x0000000000000400ULL,0x0000000000000800ULL,
0x0000000000001000ULL,0x0000000000002000ULL,0x0000000000004000ULL,0x0000000000008000ULL,
0x0000000000010000ULL,0x0000000000020000ULL,0x0000000000040000ULL,0x0000000000080000ULL,
0x0000000000100000ULL,0x0000000000200000ULL,0x0000000000400000ULL,0x0000000000800000ULL,
0x0000000001000000ULL,0x0000000002000000ULL,0x0000000004000000ULL,0x0000000008000000ULL,
0x0000000010000000ULL,0x0000000020000000ULL,0x0000000040000000ULL,0x0000000080000000ULL,
0x0000000100000000ULL,0x0000000200000000ULL,0x0000000400000000ULL,0x0000000800000000ULL,
0x0000001000000000ULL,0x0000002000000000ULL,0x0000004000000000ULL,0x0000008000000000ULL,
0x0000010000000000ULL,0x0000020000000000ULL,0x0000040000000000ULL,0x0000080000000000ULL,
0x0000100000000000ULL,0x0000200000000000ULL,0x0000400000000000ULL,0x0000800000000000ULL,
0x0001000000000000ULL,0x0002000000000000ULL,0x0004000000000000ULL,0x0008000000000000ULL,
0x0010000000000000ULL,0x0020000000000000ULL,0x0040000000000000ULL,0x0080000000000000ULL,
0x0100000000000000ULL,0x0200000000000000ULL,0x0400000000000000ULL,0x0800000000000000ULL,
0x1000000000000000ULL,0x2000000000000000ULL,0x4000000000000000ULL,0x8000000000000000ULL};

#ifdef BIT
#undef BIT
#endif
#ifdef BIT64
#undef BIT64
#endif
#define BIT64(_i) bit64[(_i)]



/***** Prototypes for locally used functions. *****/
void MultB64(u64 *Product, u64 *x, void *P);
void MultB_T64(u64 *Product, u64 *x, void *P);
void multT(u64 *res, u64 *A, u64 *B, s32 n);
void multS(u64 *D, int *S);
void mult64x64(u64 *res, u64 *A, u64 *B);
void preMult(u64 *A, u64 *B);
void multnx64(u64 *C_n, u64 *A_n, u64 *B, s32 n);
void addmultnx64(u64 *C_n, u64 *A_n, u64 *B, s32 n);
void getW_S(u64 *Wi, int *Si, u64 *T, int *Si_1);
int  isZeroV(u64 *A, s32 size);
int  doColumnOps(u64 *A, u64 *B, s32 n);
int  doColumnOps64(u64 *A, u64 *B, s32 n);
/**************************************************/


void seedBlockLanczos(s32 seed)
{ prandseed(seed, 712*seed + 21283, seed^0xF3C91D1A);
}

/**************************************************/
int testMult(nfs_sparse_mat_t *P)
/**************************************************/
/* Test the matrix multiplication routines.       */
/**************************************************/
{ u64 *u, *x, *y, *z, t;
  int  i, fail=0, totalFailures=0;
  long j, n=P->numCols;

  u = (u64 *)malloc(n*sizeof(u64));
  x = (u64 *)malloc(n*sizeof(u64));
  y = (u64 *)malloc(n*sizeof(u64));
  z = (u64 *)malloc(n*sizeof(u64));
  if (!(x&&y&&z)) {
    printf("testMult(): Memory allocation error!\n");
    exit(-1);
  }
  memset(u, 0x00, n*sizeof(u64));
  memset(x, 0x00, n*sizeof(u64));
  memset(y, 0x00, n*sizeof(u64));
  memset(z, 0x00, n*sizeof(u64));
  for (i=0; i<30; i++) {
    for (j=0; j<n; j++) {
      t = prand();
      t = (t<<32)^prand();
      u[j]=t;
      x[j] ^= t;
    }
    MultB64(y, u, (void *)P);
    for (j=0; j<n; j++)
      z[j] ^= y[j];
  }
  MultB64(y, x, (void *)P);
  for (j=0; j<n; j++) {
    if (y[j] != z[j]) {
      printf("product MultB64() failure at column %ld!\n", j);
      fail=1;
    }
  }
  if (!fail) 
    printf("First matrix product test passed.\n");

  totalFailures += fail;
  fail=0;
  memset(u, 0x00, n*sizeof(u64));
  memset(x, 0x00, n*sizeof(u64));
  memset(y, 0x00, n*sizeof(u64));
  memset(z, 0x00, n*sizeof(u64));
  for (i=0; i<30; i++) {
    for (j=0; j<n; j++) {
      t = prand();
      t = (t<<32)^prand();
      u[j]=t;
      x[j] ^= t;
    }
    MultB_T64(y, u, (void *)P);
    for (j=0; j<n; j++)
      z[j] ^= y[j];
  }
  MultB_T64(y, x, (void *)P);
  for (j=0; j<n; j++) {
    if (y[j] != z[j]) {
      printf("product MultB_T64() failure at column %ld!\n", j);
      fail=1;
    }
  }

  if (!fail) 
    printf("Second matrix product test passed.\n");

  totalFailures += fail;
  fail=0;


  return totalFailures;
}



void MultB64(u64 *Product, u64 *x, void *P) {
  nfs_sparse_mat_t *M = (nfs_sparse_mat_t *)P;
  memset(Product, 0, M->numCols * sizeof(u64)); 
  {
    int i;
    for (i = 0; i < M->numDenseBlocks; i++) {
      multT(Product + M->denseBlockIndex[i], M->denseBlocks[i], x, M->numCols);
    }
  }
#if defined(L2_CACHE_SIZE) && (L2_CACHE_SIZE > 0)
// L2_CACHE_SIZE has to be a power of 2.
// MULTB64_PAGESIZE is a half of L2 cache size.
#define MULTB64_PAGESIZE (L2_CACHE_SIZE * 1024 / 2 / sizeof(u64))
#if L2_CACHE_SIZE == 256
#define MULTB64_PAGEMASK "-16384"
#elif L2_CACHE_SIZE == 512
#define MULTB64_PAGEMASK "-32768"
#elif L2_CACHE_SIZE == 1024
#define MULTB64_PAGEMASK "-65536"
#else
#error 1
#endif
  {
    s32 n = M->numCols;
    u32 *cEntry = M->cEntry;
    s32 *cIndex = M->cIndex;
    u32 pagestart;
    for (pagestart = 0; pagestart < n; pagestart += MULTB64_PAGESIZE) {
      asm volatile("\
	movl	%0, %%esi			#cEntry		\n\
	movl	%1, %%edi			#Product	\n\
	xorl	%%ecx, %%ecx			#i		\n\
	xorl	%%edx, %%edx			#p		\n\
1:								\n\
	movl	%2, %%eax			#cIndex		\n\
	movl	4(%%eax,%%ecx,4), %%ebx		#s=cIndex[i+1]	\n\
	subl	(%%eax,%%ecx,4), %%ebx		#s-=cIndex[i]	\n\
	jle	6f						\n\
	movl	%3, %%eax			#x		\n\
	movq	(%%eax,%%ecx,8), %%mm1		#t=x[i]		\n\
	andl	$-16, %%ebx			#s&=-16		\n\
	jz	3f						\n\
	addl	%%edx, %%ebx			#s+=p		\n\
2:								\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	andl	$"MULTB64_PAGEMASK", %%eax	#&pagemask	\n\
	cmpl	%5, %%eax			#==pagestart	\n\
	jne	7f						\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	movq	(%%edi,%%eax,8), %%mm0		#Product[p[0]]	\n\
	pxor	%%mm1, %%mm0			#^=t		\n\
	movq	%%mm0, (%%edi,%%eax,8)				\n\
7:								\n\
	.set	n, 1				#p[1]..p[15]	\n\
	.rept	15						\n\
		movl	4*n(%%esi,%%edx,4), %%eax		\n\
		andl	$"MULTB64_PAGEMASK", %%eax		\n\
		cmpl	%5, %%eax				\n\
		jne	7f					\n\
		movl	4*n(%%esi,%%edx,4), %%eax		\n\
		movq	(%%edi,%%eax,8), %%mm0			\n\
		pxor	%%mm1, %%mm0				\n\
		movq	%%mm0, (%%edi,%%eax,8)			\n\
7:								\n\
		.set	n, n+1					\n\
	.endr							\n\
	addl	$16, %%edx			#p+=16		\n\
	cmpl	%%ebx, %%edx			#p<s		\n\
	jl	2b						\n\
3:								\n\
	movl	%2, %%eax			#cIndex		\n\
	movl	4(%%eax,%%ecx,4), %%ebx		#s=cIndex[i+1]	\n\
	cmpl	%%ebx, %%edx					\n\
	jge	5f						\n\
4:								\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	andl	$"MULTB64_PAGEMASK", %%eax	#&pagemask	\n\
	cmpl	%5, %%eax			#==pagestart	\n\
	jne	7f						\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	movq	(%%edi,%%eax,8), %%mm0		#Product[p[0]]	\n\
	pxor	%%mm1, %%mm0			#^=t		\n\
	movq	%%mm0, (%%edi,%%eax,8)				\n\
7:								\n\
	addl	$1, %%edx			#p++		\n\
	cmpl	%%ebx, %%edx			#p<s		\n\
	jl	4b						\n\
5:								\n\
6:								\n\
	addl	$1, %%ecx			#i++		\n\
	cmpl	%4, %%ecx			#i<n		\n\
	jl	1b						\n\
	emms" : : "m"(cEntry), "m"(Product), "m"(cIndex), "m"(x), "m"(n),
                   "m"(pagestart) :
                   "%eax", "%ebx", "%ecx", "%edx", "%esi", "%edi");
    }
  }
#else
  {
    s32 n = M->numCols;
    u32 *cEntry = M->cEntry;
    s32 *cIndex = M->cIndex;
    asm volatile("\
	movl	%0, %%esi			#cEntry		\n\
	movl	%1, %%edi			#Product	\n\
	xorl	%%ecx, %%ecx			#i		\n\
	xorl	%%edx, %%edx			#p		\n\
1:								\n\
	movl	%2, %%eax			#cIndex		\n\
	movl	4(%%eax,%%ecx,4), %%ebx		#s=cIndex[i+1]	\n\
	subl	(%%eax,%%ecx,4), %%ebx		#s-=cIndex[i]	\n\
	jle	6f						\n\
	movl	%3, %%eax			#x		\n\
	movq	(%%eax,%%ecx,8), %%mm1		#t=x[i]		\n\
	andl	$-16, %%ebx			#s&=-16		\n\
	jz	3f						\n\
	addl	%%edx, %%ebx			#s+=p		\n\
2:								\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	movq	%%mm1, %%mm0					\n\
	pxor	(%%edi,%%eax,8), %%mm0		#Product[p[0]]	\n\
	movq	%%mm0, (%%edi,%%eax,8)		#^=t		\n\
	.set	n, 1				#p[1]..p[15]	\n\
	.rept	15						\n\
		movl	4*n(%%esi,%%edx,4), %%eax		\n\
		movq	%%mm1, %%mm0				\n\
		pxor	(%%edi,%%eax,8), %%mm0			\n\
		movq	%%mm0, (%%edi,%%eax,8)			\n\
		.set	n, n+1					\n\
	.endr							\n\
	addl	$16, %%edx			#p+=16		\n\
	cmpl	%%ebx, %%edx			#p<s		\n\
	jl	2b						\n\
3:								\n\
	movl	%2, %%eax			#cIndex		\n\
	movl	4(%%eax,%%ecx,4), %%ebx		#s=cIndex[i+1]	\n\
	cmpl	%%ebx, %%edx					\n\
	jge	5f						\n\
4:								\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	movq	%%mm1, %%mm0					\n\
	pxor	(%%edi,%%eax,8), %%mm0		#Product[p[0]]	\n\
	movq	%%mm0, (%%edi,%%eax,8)		#^=t		\n\
	addl	$1, %%edx			#p++		\n\
	cmpl	%%ebx, %%edx			#p<s		\n\
	jl	4b						\n\
5:								\n\
6:								\n\
	addl	$1, %%ecx			#i++		\n\
	cmpl	%4, %%ecx			#i<n		\n\
	jl	1b						\n\
	emms" : : "m"(cEntry), "m"(Product), "m"(cIndex), "m"(x), "m"(n) :
                 "%eax", "%ebx", "%ecx", "%edx", "%esi", "%edi");
  }
#endif
}

void MultB_T64(u64 *Product, u64 *x, void *P) {
  nfs_sparse_mat_t *M = (nfs_sparse_mat_t *)P;
  memset(Product, 0, M->numCols * sizeof(u64));
  {
    int i;
    for (i = 0; i < M->numDenseBlocks; i++) {
      addmultnx64(Product, M->denseBlocks[i], x + M->denseBlockIndex[i], M->numCols);
    }
  }
#if defined(L2_CACHE_SIZE) && (L2_CACHE_SIZE > 0)
// L2_CACHE_SIZE has to be a power of 2.
// MULTB_T64_PAGESIZE is a half of L2 cache size.
#define MULTB_T64_PAGESIZE (L2_CACHE_SIZE * 1024 / 2 / sizeof(u64))
#if L2_CACHE_SIZE == 256
#define MULTB_T64_PAGEMASK "-16384"
#elif L2_CACHE_SIZE == 512
#define MULTB_T64_PAGEMASK "-32768"
#elif L2_CACHE_SIZE == 1024
#define MULTB_T64_PAGEMASK "-65536"
#else
#error 1
#endif
  {
    s32 n = M->numCols;
    u32 *cEntry = M->cEntry;
    s32 *cIndex = M->cIndex;
    u32 pagestart;
    for (pagestart = 0; pagestart < n; pagestart += MULTB_T64_PAGESIZE) {
      asm volatile("\
	movl	%0, %%esi			#cEntry		\n\
	movl	%1, %%edi			#x		\n\
	xorl	%%ecx, %%ecx			#i		\n\
	xorl	%%edx, %%edx			#p		\n\
1:								\n\
	movl	%2, %%eax			#cIndex		\n\
	movl	4(%%eax,%%ecx,4), %%ebx		#s=cIndex[i+1]	\n\
	subl	(%%eax,%%ecx,4), %%ebx		#s-=cIndex[i]	\n\
	jle	6f						\n\
	movl	%3, %%eax			#Product	\n\
	movq	(%%eax,%%ecx,8), %%mm0		#t=Product[i]	\n\
	andl	$-16, %%ebx			#s&=-16		\n\
	jz	3f						\n\
	addl	%%edx, %%ebx			#s+=p		\n\
2:								\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	andl	$"MULTB_T64_PAGEMASK", %%eax	#&pagemask	\n\
	cmpl	%5, %%eax			#==pagestart	\n\
	jne	7f						\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	pxor	(%%edi,%%eax,8), %%mm0		#t^=x[p[0]]	\n\
7:								\n\
	.set	n, 1				#p[1]..p[15]	\n\
	.rept	15						\n\
		movl	4*n(%%esi,%%edx,4), %%eax		\n\
		andl	$"MULTB_T64_PAGEMASK", %%eax		\n\
		cmpl	%5, %%eax				\n\
		jne	7f					\n\
		movl	4*n(%%esi,%%edx,4), %%eax		\n\
		pxor	(%%edi,%%eax,8), %%mm0			\n\
7:								\n\
		.set	n, n+1					\n\
	.endr							\n\
	addl	$16, %%edx			#p+=16		\n\
	cmpl	%%ebx, %%edx			#p<s		\n\
	jl	2b						\n\
3:								\n\
	movl	%2, %%eax			#cIndex		\n\
	movl	4(%%eax,%%ecx,4), %%ebx		#s=cIndex[i+1]	\n\
	cmpl	%%ebx, %%edx					\n\
	jge	5f						\n\
4:								\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	andl	$"MULTB_T64_PAGEMASK", %%eax	#&pagemask	\n\
	cmpl	%5, %%eax			#==pagestart	\n\
	jne	7f						\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	pxor	(%%edi,%%eax,8), %%mm0		#t^=x[p[0]]	\n\
7:								\n\
	addl	$1, %%edx			#p++		\n\
	cmpl	%%ebx, %%edx			#p<s		\n\
	jl	4b						\n\
5:								\n\
	movl	%3, %%eax			#Product	\n\
	movq	%%mm0, (%%eax,%%ecx,8)		#Product[i]=t	\n\
6:								\n\
	addl	$1, %%ecx			#i++		\n\
	cmpl	%4, %%ecx			#i<n		\n\
	jl	1b						\n\
	emms" : : "m"(cEntry), "m"(x), "m"(cIndex), "m"(Product), "m"(n),
                   "m"(pagestart) :
                   "%eax", "%ebx", "%ecx", "%edx", "%esi", "%edi");
    }
  }
#else
  {
    s32 n = M->numCols;
    u32 *cEntry = M->cEntry;
    s32 *cIndex = M->cIndex;
    asm volatile("\
	movl	%0, %%esi			#cEntry		\n\
	movl	%1, %%edi			#x		\n\
	xorl	%%ecx, %%ecx			#i		\n\
	xorl	%%edx, %%edx			#p		\n\
1:								\n\
	movl	%2, %%eax			#cIndex		\n\
	movl	4(%%eax,%%ecx,4), %%ebx		#s=cIndex[i+1]	\n\
	subl	(%%eax,%%ecx,4), %%ebx		#s-=cIndex[i]	\n\
	jle	6f						\n\
	movl	%3, %%eax			#Product	\n\
	movq	(%%eax,%%ecx,8), %%mm0		#t=Product[i]	\n\
	andl	$-16, %%ebx			#s&=-16		\n\
	jz	3f						\n\
	addl	%%edx, %%ebx			#s+=p		\n\
2:								\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	pxor	(%%edi,%%eax,8), %%mm0		#t^=x[p[0]]	\n\
	.set	n, 1				#p[1]..p[15]	\n\
	.rept	15						\n\
		movl	4*n(%%esi,%%edx,4), %%eax		\n\
		pxor	(%%edi,%%eax,8), %%mm0			\n\
		.set	n, n+1					\n\
	.endr							\n\
	addl	$16, %%edx			#p+=16		\n\
	cmpl	%%ebx, %%edx			#p<s		\n\
	jl	2b						\n\
3:								\n\
	movl	%2, %%eax			#cIndex		\n\
	movl	4(%%eax,%%ecx,4), %%ebx		#s=cIndex[i+1]	\n\
	cmpl	%%ebx, %%edx			#p<s		\n\
	jge	5f						\n\
4:								\n\
	movl	(%%esi,%%edx,4), %%eax		#p[0]		\n\
	pxor	(%%edi,%%eax,8), %%mm0		#t^=x[p[0]]	\n\
	addl	$1, %%edx			#p++		\n\
	cmpl	%%ebx, %%edx			#p<s		\n\
	jl	4b						\n\
5:								\n\
	movl	%3, %%eax			#Product	\n\
	movq	%%mm0, (%%eax,%%ecx,8)		#Product[i]=t	\n\
6:								\n\
	addl	$1, %%ecx			#i++		\n\
	cmpl	%4, %%ecx			#i<n		\n\
	jl	1b						\n\
	emms" : : "m"(cEntry), "m"(x), "m"(cIndex), "m"(Product), "m"(n) :
                 "%eax", "%ebx", "%ecx", "%edx", "%esi", "%edi");
  }
#endif
}

/**********************************************************************/
int blockLanczos64(u64 *deps, MAT_MULT_FUNC_PTR64 MultB, 
                  MAT_MULT_FUNC_PTR64 MultB_T, void *P, s32 n)
/**********************************************************************/
{ u64 *Y=NULL, *X=NULL, *Vi=NULL, *Vi_1=NULL, *Vi_2=NULL, *tmp_n=NULL, *tmp2_n=NULL;
  u64 *V0=NULL, *tmp3_n=NULL, *Z=NULL, *AZ=NULL;
  u64 D[64] ALIGNED16, E[64] ALIGNED16, F[64] ALIGNED16, Wi[64] ALIGNED16;
  u64 Wi_1[64] ALIGNED16, Wi_2[64] ALIGNED16, T[64] ALIGNED16, T_1[64] ALIGNED16;
  u64 tmp[64] ALIGNED16;
  u64 U[64] ALIGNED16, U_1[64] ALIGNED16, tmp2[64] ALIGNED16;
  int  Si[64], Si_1[64];
  u64 i, j, m, mask, isZero, r1,r2;
  u32  iterations;
  int  errs=0, numDeps=-1, cont, s;
  double startTime, now, estTotal;

  if (testMult((nfs_sparse_mat_t *)P)) {
    printf("Self test reported some errors! Stopping...\n");
    exit(-1);
  }

  /* Memory allocation: */
  if (malloc_aligned64(Y, 16, n)) errs++;
  if (malloc_aligned64(X, 16, n)) errs++;
  if (malloc_aligned64(Vi, 16, n)) errs++;
  if (malloc_aligned64(V0, 16, n)) errs++;
  if (malloc_aligned64(Vi_1, 16, n)) errs++;
  if (malloc_aligned64(Vi_2, 16, n)) errs++;
  if (malloc_aligned64(tmp_n, 16, n)) errs++;
  if (malloc_aligned64(tmp2_n, 16, n)) errs++;
  if (malloc_aligned64(tmp3_n, 16, n)) errs++;
  if (errs) {
    fprintf(stderr, "blanczos(): Memory allocation error!\n");
    goto SHORT_CIRC_STOP;
  }
  
  /******************************************************************/
  /* Throughout, 'A' means the matrix A := (B^T)B. In fact, all the */
  /* notation is directly from Montgomery's paper, except that I    */
  /* don't bother indexing the D_i, E_i, F_i. They are simply D,E,F.*/
  /******************************************************************/
  /* Initialization: */
  for (j=0; j<n; j++) {
    r1 = prand(); r2 = prand(); 
    Y[j] = r1^(r2<<32);
    X[j] = Vi_1[j] = Vi_2[j] = tmp_n[j] = tmp2_n[j] = tmp3_n[j] = 0;
  }
  for (i=0; i<64; i++) {
    Wi[i] = Wi_1[i] = Wi_2[i] = 0;
    T[i] = T_1[i] = U[i] = U_1[i] = tmp[i] = 0;
    Si[i] = Si_1[i] = i; /* Si and Si_1 are both I_64. */
  }
  MultB(tmp_n, Y, P); 
  MultB_T(V0, tmp_n, P);
  memcpy(Vi, V0, n*sizeof(u64));

  /* Prime 'T' for the loop, so that we always have */
  /* T = (Vi^T)A(Vi) and T_1 = (Vi_1^T)A(Vi_1).     */
  
  MultB(tmp2_n, Vi, P);
  MultB_T(tmp_n, tmp2_n, P);  /* tmp_n <-- A*Vi */
  multT(T, Vi, tmp_n, n);      /* T <-- (Vi^T)(tmp_n) = (Vi^T)A(Vi) */
  
  cont = 1;
  i = 0;
  getW_S(Wi, Si, T, Si_1); /* Compute W0 and S0. */
  /* Initialize X <-- (V0)(W0)(V0^T)(V0). */
  multT(tmp, V0, V0, n);    /* tmp <-- (V0^T)(V0). */
  mult64x64(tmp2, Wi, tmp); /* tmp2 <-- (W0)(tmp) = (W0)(V0^T)(V0). */
  multnx64(X, V0, tmp2, n); /* X <-- V0(tmp2). */
  iterations = 0;
  startTime = sTime();
  do {
    /* Iteration step. */
    iterations++;

    /********** Compute D_{i+1}. **********/
    /* tmp_n = A*Vi from initialization, or the previous iteration. */

    multT(U, tmp_n, tmp_n, n);  /* U <-- (tmp_n)^T(tmp_n) = (Vi^T)(A^2)(Vi) */
    multS(U, Si);            /* U <-- (Vi^T)(A^2)(Vi)(Si)(Si^T). */
    memcpy(D, U, 64*sizeof(u64)); /* D <-- U. */
    for (j=0; j<64; j++)
      D[j] ^= T[j]; /* D <-- D + (Vi^T)A(Vi). */
    preMult(D, Wi); /* D <-- (Wi)D. */
    for (j=0; j<64; j++)
      D[j] ^= BIT64(j); /* D <-- D + I_{64}. */

    /********** Compute E_{i+1}. **********/
    mult64x64(E, Wi_1, T); /* E <-- (Wi_1)(Vi^T)A(Vi). */
    multS(E, Si);          /* E <-- E(Si)(Si^T).       */

    /**************** Compute F_{i+1}. *************/
    /* Note: We should, at this point, always have */
    /* T_1 = (Vi_1^T)A(Vi_1) and                   */
    /* U_1 = (Vi_1^T)A^2(Vi_1)(Si_1)(Si_1^T)       */
    /***********************************************/
    for (j=0; j<64; j++)
      F[j] = U_1[j] ^ T_1[j];
    multS(F, Si);

    mult64x64(tmp, T_1, Wi_1); /* tmp <-- (Vi_1^T)A(Vi_1)(Wi_1). */
    for (j=0; j<64; j++)
      tmp[j] ^= BIT64(j); /* tmp <-- tmp + I_64. */
    preMult(tmp, Wi_2);
    preMult(F, tmp);
    /*** Done computing 'F'. */

    /* Finally, do the actual iteration step, putting V_{i+1} in 'tmp_n'. */
    /* Right now, we still have tmp_n = A*Vi, so we should first use it.  */
    
    /* Do tmp_n <-- tmp_n*(Si)(Si^T). */
    mask = 0x00000000;
    for (j=0; j<64; j++) {
      s = Si[j];
      if ((s>=0) && (s<64))
        mask |= BIT64(s);
    }
    for (j=0; j<n; j++) 
      tmp_n[j] &= mask;
    
    addmultnx64(tmp_n, Vi, D, n);   /* tmp_n <-- tmp_n + (Vi)D.        */

    addmultnx64(tmp_n, Vi_1, E, n); /* tmp_n <-- tmp_n + (Vi_1)E */

    addmultnx64(tmp_n, Vi_2, F, n); /* tmp_n <-- tmp_n + (Vi_2)F */
    /*** Done. We now have tmp_n = V_{i+1}. ***/
   

    i++;
    memcpy(Vi_2, Vi_1, n*sizeof(u64));
    memcpy(Vi_1, Vi, n*sizeof(u64));
    memcpy(Vi, tmp_n, n*sizeof(u64));
    memcpy(Wi_2, Wi_1, 64*sizeof(u64));
    memcpy(Wi_1, Wi, 64*sizeof(u64));
    memcpy(T_1, T, 64*sizeof(u64));
    memcpy(U_1, U, 64*sizeof(u64));
    memcpy(Si_1, Si, 64*sizeof(int));
    
    /******** My `step 3'. ********/
    /* Is (Vi^T)(A)(Vi) == 0 ? */
    MultB(tmp2_n, Vi, P);
    MultB_T(tmp_n, tmp2_n, P);  /* tmp_n <-- A*Vi       */

    multT(T, Vi, tmp_n, n); /* T <-- (Vi^T)A(Vi). */
    getW_S(Wi, Si, T, Si_1);
    
    if (!(isZeroV(T, 64))) {
      /* X <-- X + (Vi)(Wi)(Vi^T)(V_0) */
      multT(tmp, Vi, V0, n); /* tmp <-- (Vi^T)(V0). */
      preMult(tmp, Wi);   /* tmp <-- (Wi)(tmp).  */
      
      addmultnx64(X, Vi, tmp, n);   /* X <-- X + (Vi)(tmp)    */
    } else {
      cont=0;
      m = i;
    }
    now = sTime();
    estTotal = ((double)1.02*n/(iterations*64.0))*(now-startTime);
    printTmp("Lanczos: Estimate %1.1lf%% complete (%1.1lf seconds / %1.1lf seconds)...",
              (double)100.0*64.0*iterations/(1.02*n), now-startTime, estTotal);  
    if ((double)100.0*64.0*iterations/n > 250) {
      fprintf(stderr, "Some error has occurred: Lanczos is not converging!\n");
      fprintf(stderr, "Number of iterations is %" PRIu32 ".\n", iterations);
      /* Add some debugging stuff here! */
      fprintf(stderr, "Terminating...\n");
      exit(-1);
    }
  } while (cont);
  printf("\nBlock Lanczos used %" PRIu32 " iterations.\n", iterations);

          
  if (malloc_aligned64(Z, 16, 2*n)) Z=NULL;
  if (malloc_aligned64(AZ, 16, 2*n)) AZ=NULL;
/*
  Z = (u64 *)malloc(2*n*sizeof(u64));
  AZ = (u64 *)malloc(2*n*sizeof(u64));
*/
  if (!(Z&&AZ)) {
    fprintf(stderr, "blanczos(): Memory allocation error!\n");
    goto SHORT_CIRC_STOP;
  }
  
  if (isZeroV(Vi, n)) {
    /* Then <X+Y> < ker(A). Later, we will expect AX=0, so do X <-- X+Y. */
    printf("After Block Lanczos iteration, Vm=0.\n");

    for (i=0; i<n; i++)
      X[i] ^= Y[i];
  } else {
    printf("After Block Lanczos iteration, Vm is nonzero. Finishing...\n");
    /* We need more memory (for convenience), so free off what we don't need. */
    free_aligned64(V0);   free_aligned64(Vi_1); 
    free_aligned64(Vi_2); free_aligned64(tmp3_n);
    V0 = Vi_1 = Vi_2 = tmp3_n = NULL;

    /* Construct Z=[ X+Y | Vi] and compute AZ=[ A(X+Y) | AVi ] */
    /* X <-- X+Y, for convenience.                     */
    for (j=0; j<n; j++) {
      X[j] ^= Y[j];
      Z[2*j] = X[j];
      Z[2*j+1] = Vi[j];
    }
    /* Now, compute AZ = [ AX | AVi ]. */
    MultB(tmp2_n, X, P);  
    MultB_T(tmp_n, tmp2_n, P); /* tmp_n <-- AX       */
    for (j=0; j<n; j++) 
      AZ[2*j] = tmp_n[j];
    MultB(tmp2_n, Vi, P);  
    MultB_T(tmp_n, tmp2_n, P); /* tmp_n <-- AVi       */
    for (j=0; j<n; j++) 
      AZ[2*j+1] = tmp_n[j];
    /* Now, AZ should have small rank, and we should be  */
    /* able to do column ops to make zero vectors. Doing */
    /* the same ops to 'Z', we get vectors in ker(A).    */
    doColumnOps(Z, AZ, n);
    /* Now, look for zero columns in AZ. If a column is zero, */
    /* copy the corresponding column of Z into 'X'.           */
    for (i=0; i<n; i++)
      X[i] = 0;
    numDeps=0;
    for (i=0; (i<64) && (numDeps < 64); i++) {
      j=0;
      while ((j<n) && ((AZ[2*j + i/64]&BIT64(i%64))==0))
        j++;
      if (j==n) {
        /* Copy column i of Z into column 'numDeps' of X. */
        for (j=0; j<n; j++) {
          if (Z[2*j + i/64]&BIT64(i%64))
            X[j] ^= BIT64(numDeps);
        }
        numDeps++;
      }
    }
    printf("Found %d dependencies for A=(B^T)B.\n", numDeps);
  }
  j=0;
  while ((j<n)&&(X[j]==0))
    j++;
  if (j==n)
    printf("Probable error: The matrix X is identically zero!!!\n");


  printf("Getting dependencies for original matrix, B...\n");
  /***************************************************/
  /* At this point, we should have AX = (B^T)BX = 0. */
  /* Do simultaneous column ops to BX and X, putting */
  /* BX in RREF, to get vectors in ker(B).           */
  /***************************************************/

  for (i=0; i<n; i++)
    deps[i] = 0;
  numDeps=0;

  MultB(tmp_n, X, P);  /* tmp_n <-- BX. */
  doColumnOps64(X, tmp_n, n);
  /* We only want 32 of the dependencies. */
  for (i=0; i<32; i++) {
    for (j=0, isZero=1; j<n; j++)
      if (tmp_n[j]&BIT64(i)) {
        isZero=0;
      }
    if (isZero) {
      for (j=0, isZero=1; j<n; j++) {
        if (X[j]&BIT64(i)) {
          deps[j] ^= BIT64(numDeps);
          isZero=0;
        }
      }
      if (!(isZero))
        numDeps++;
    }
  }

  if (numDeps) {
    printf("Found %d dependencies for 'B'. Verifying...\n", numDeps);
  } else {
    printf("Some error occurred: all dependencies found seem to be trivial!\n");
    printf("Is Rank(A) = Rank((B^T)B) too small?\n");
    goto SHORT_CIRC_STOP;
  }

  MultB(tmp_n, deps, P);
  for (i=0, isZero=1; (i<n)&&(isZero); i++)
    if (tmp_n[i])
      isZero = 0;
  if (!(isZero))
    printf("Some error occurred: Final product (B)(deps) is nonzero (i=%" PRIu64 ")!\n", i);
  else
    printf("Verified.\n"); 

SHORT_CIRC_STOP:
  
  if (Y != NULL)      free_aligned64(Y);
  if (X != NULL)      free_aligned64(X);
  if (Vi != NULL)     free_aligned64(Vi);
  if (V0 != NULL)     free_aligned64(V0);
  if (Vi_1 != NULL)   free_aligned64(Vi_1);
  if (Vi_2 != NULL)   free_aligned64(Vi_2);
  if (tmp_n != NULL)  free_aligned64(tmp_n);
  if (tmp2_n != NULL) free_aligned64(tmp2_n);
  if (tmp3_n != NULL) free_aligned64(tmp3_n);
  if (Z != NULL)      free_aligned64(Z);
  if (AZ != NULL)     free_aligned64(AZ);
  return numDeps;
}


u64 mult_w[2048] ALIGNED16;

void multT(u64 *c, u64 *a, u64 *b, s32 n) {
  memset(mult_w, 0, sizeof(u64) * 256 * 8);
  asm volatile("\
	movl	%0, %%esi					\n\
	movl	%1, %%edi					\n\
	movl	%2, %%ecx					\n\
	leal	(%%esi,%%ecx,8), %%esi				\n\
	leal	(%%edi,%%ecx,8), %%edi				\n\
	negl	%%ecx						\n\
1:								\n\
	movq	(%%edi,%%ecx,8), %%mm1				\n\
	movl	(%%esi,%%ecx,8), %%ebx				\n\
	movzbl	%%bl, %%eax					\n\
	movzbl	%%bh, %%edx					\n\
	movq	%%mm1, %%mm0					\n\
	pxor	"ASM_UL"mult_w(,%%eax,8), %%mm0			\n\
	movq	%%mm0, "ASM_UL"mult_w(,%%eax,8)			\n\
	shrl	$16, %%ebx					\n\
	movq	%%mm1, %%mm0					\n\
	pxor	"ASM_UL"mult_w+8*256*1(,%%edx,8), %%mm0		\n\
	movq	%%mm0, "ASM_UL"mult_w+8*256*1(,%%edx,8)		\n\
	movzbl	%%bl, %%eax					\n\
	movl	4(%%esi,%%ecx,8), %%edx				\n\
	shrl	$8, %%ebx					\n\
	movq	%%mm1, %%mm0					\n\
	pxor	"ASM_UL"mult_w+8*256*2(,%%eax,8), %%mm0		\n\
	movq	%%mm0, "ASM_UL"mult_w+8*256*2(,%%eax,8)		\n\
	movzbl	%%dl, %%eax					\n\
	movq	%%mm1, %%mm0					\n\
	pxor	"ASM_UL"mult_w+8*256*3(,%%ebx,8), %%mm0		\n\
	movq	%%mm0, "ASM_UL"mult_w+8*256*3(,%%ebx,8)		\n\
	movzbl	%%dh, %%ebx					\n\
	shrl	$16, %%edx					\n\
	movq	%%mm1, %%mm0					\n\
	pxor	"ASM_UL"mult_w+8*256*4(,%%eax,8), %%mm0		\n\
	movq	%%mm0, "ASM_UL"mult_w+8*256*4(,%%eax,8)		\n\
	movzbl	%%dl, %%eax					\n\
	movq	%%mm1, %%mm0					\n\
	pxor	"ASM_UL"mult_w+8*256*5(,%%ebx,8), %%mm0		\n\
	movq	%%mm0, "ASM_UL"mult_w+8*256*5(,%%ebx,8)		\n\
	shrl	$8, %%edx					\n\
	movq	%%mm1, %%mm0					\n\
	pxor	"ASM_UL"mult_w+8*256*6(,%%eax,8), %%mm0		\n\
	movq	%%mm0, "ASM_UL"mult_w+8*256*6(,%%eax,8)		\n\
	movq	%%mm1, %%mm0					\n\
	pxor	"ASM_UL"mult_w+8*256*7(,%%edx,8), %%mm0		\n\
	movq	%%mm0, "ASM_UL"mult_w+8*256*7(,%%edx,8)		\n\
	addl	$1, %%ecx					\n\
	jnz	1b						\n\
	emms" : : "m"(a), "m"(b), "m"(n) :
               "%eax", "%ebx", "%ecx", "%edx", "%esi", "%edi");
  {
    int i;
    for (i = 0; i < 64; i += 8) {
      asm volatile("\
	movq	2040(%0), %%mm0		#r0 = a[255]	255->254	\n\
	pxor	2032(%0), %%mm0		#r0 ^= a[254]			\n\
	movq	%%mm0, 2032(%0)		#a[254] = r0			\n\
	pxor	2016(%0), %%mm0		#r0 ^= a[252]	254->252	\n\
	pxor	2024(%0), %%mm0		#r0 ^= a[253]	253->252	\n\
	movq	%%mm0, 2016(%0)		#a[252] = r0			\n\
	movq	2008(%0), %%mm1		#r1 = a[251]	251->250	\n\
	pxor	2000(%0), %%mm1		#r1 ^= a[250]			\n\
	movq	%%mm1, 2000(%0)		#a[250] = r1			\n\
	pxor	1984(%0), %%mm0		#r0 ^= a[248]	252->248	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	250->248	\n\
	pxor	1992(%0), %%mm0		#r0 ^= a[249]	249->248	\n\
	movq	%%mm0, 1984(%0)		#a[248] = r0			\n\
	movq	1976(%0), %%mm1		#r1 = a[247]	247->246	\n\
	pxor	1968(%0), %%mm1		#r1 ^= a[246]			\n\
	movq	%%mm1, 1968(%0)		#a[246] = r1			\n\
	pxor	1952(%0), %%mm1		#r1 ^= a[244]	246->244	\n\
	pxor	1960(%0), %%mm1		#r1 ^= a[245]	245->244	\n\
	movq	%%mm1, 1952(%0)		#a[244] = r1			\n\
	movq	1944(%0), %%mm2		#r2 = a[243]	243->242	\n\
	pxor	1936(%0), %%mm2		#r2 ^= a[242]			\n\
	movq	%%mm2, 1936(%0)		#a[242] = r2			\n\
	pxor	1920(%0), %%mm0		#r0 ^= a[240]	248->240	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	244->240	\n\
	pxor	%%mm2, %%mm0		#r0 ^= r2	242->240	\n\
	pxor	1928(%0), %%mm0		#r0 ^= a[241]	241->240	\n\
	movq	%%mm0, 1920(%0)		#a[240] = r0			\n\
	movq	1912(%0), %%mm1		#r1 = a[239]	239->238	\n\
	pxor	1904(%0), %%mm1		#r1 ^= a[238]			\n\
	movq	%%mm1, 1904(%0)		#a[238] = r1			\n\
	pxor	1888(%0), %%mm1		#r1 ^= a[236]	238->236	\n\
	pxor	1896(%0), %%mm1		#r1 ^= a[237]	237->236	\n\
	movq	%%mm1, 1888(%0)		#a[236] = r1			\n\
	movq	1880(%0), %%mm2		#r2 = a[235]	235->234	\n\
	pxor	1872(%0), %%mm2		#r2 ^= a[234]			\n\
	movq	%%mm2, 1872(%0)		#a[234] = r2			\n\
	pxor	1856(%0), %%mm1		#r1 ^= a[232]	236->232	\n\
	pxor	%%mm2, %%mm1		#r1 ^= r2	234->232	\n\
	pxor	1864(%0), %%mm1		#r1 ^= a[233]	233->232	\n\
	movq	%%mm1, 1856(%0)		#a[232] = r1			\n\
	movq	1848(%0), %%mm2		#r2 = a[231]	231->230	\n\
	pxor	1840(%0), %%mm2		#r2 ^= a[230]			\n\
	movq	%%mm2, 1840(%0)		#a[230] = r2			\n\
	pxor	1824(%0), %%mm2		#r2 ^= a[228]	230->228	\n\
	pxor	1832(%0), %%mm2		#r2 ^= a[229]	229->228	\n\
	movq	%%mm2, 1824(%0)		#a[228] = r2			\n\
	movq	1816(%0), %%mm3		#r3 = a[227]	227->226	\n\
	pxor	1808(%0), %%mm3		#r3 ^= a[226]			\n\
	movq	%%mm3, 1808(%0)		#a[226] = r3			\n\
	pxor	1792(%0), %%mm0		#r0 ^= a[224]	240->224	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	232->224	\n\
	pxor	%%mm2, %%mm0		#r0 ^= r2	228->224	\n\
	pxor	%%mm3, %%mm0		#r0 ^= r3	226->224	\n\
	pxor	1800(%0), %%mm0		#r0 ^= a[225]	225->224	\n\
	movq	%%mm0, 1792(%0)		#a[224] = r0			\n\
	movq	1784(%0), %%mm1		#r1 = a[223]	223->222	\n\
	pxor	1776(%0), %%mm1		#r1 ^= a[222]			\n\
	movq	%%mm1, 1776(%0)		#a[222] = r1			\n\
	pxor	1760(%0), %%mm1		#r1 ^= a[220]	222->220	\n\
	pxor	1768(%0), %%mm1		#r1 ^= a[221]	221->220	\n\
	movq	%%mm1, 1760(%0)		#a[220] = r1			\n\
	movq	1752(%0), %%mm2		#r2 = a[219]	219->218	\n\
	pxor	1744(%0), %%mm2		#r2 ^= a[218]			\n\
	movq	%%mm2, 1744(%0)		#a[218] = r2			\n\
	pxor	1728(%0), %%mm1		#r1 ^= a[216]	220->216	\n\
	pxor	%%mm2, %%mm1		#r1 ^= r2	218->216	\n\
	pxor	1736(%0), %%mm1		#r1 ^= a[217]	217->216	\n\
	movq	%%mm1, 1728(%0)		#a[216] = r1			\n\
	movq	1720(%0), %%mm2		#r2 = a[215]	215->214	\n\
	pxor	1712(%0), %%mm2		#r2 ^= a[214]			\n\
	movq	%%mm2, 1712(%0)		#a[214] = r2			\n\
	pxor	1696(%0), %%mm2		#r2 ^= a[212]	214->212	\n\
	pxor	1704(%0), %%mm2		#r2 ^= a[213]	213->212	\n\
	movq	%%mm2, 1696(%0)		#a[212] = r2			\n\
	movq	1688(%0), %%mm3		#r3 = a[211]	211->210	\n\
	pxor	1680(%0), %%mm3		#r3 ^= a[210]			\n\
	movq	%%mm3, 1680(%0)		#a[210] = r3			\n\
	pxor	1664(%0), %%mm1		#r1 ^= a[208]	216->208	\n\
	pxor	%%mm2, %%mm1		#r1 ^= r2	212->208	\n\
	pxor	%%mm3, %%mm1		#r1 ^= r3	210->208	\n\
	pxor	1672(%0), %%mm1		#r1 ^= a[209]	209->208	\n\
	movq	%%mm1, 1664(%0)		#a[208] = r1			\n\
	movq	1656(%0), %%mm2		#r2 = a[207]	207->206	\n\
	pxor	1648(%0), %%mm2		#r2 ^= a[206]			\n\
	movq	%%mm2, 1648(%0)		#a[206] = r2			\n\
	pxor	1632(%0), %%mm2		#r2 ^= a[204]	206->204	\n\
	pxor	1640(%0), %%mm2		#r2 ^= a[205]	205->204	\n\
	movq	%%mm2, 1632(%0)		#a[204] = r2			\n\
	movq	1624(%0), %%mm3		#r3 = a[203]	203->202	\n\
	pxor	1616(%0), %%mm3		#r3 ^= a[202]			\n\
	movq	%%mm3, 1616(%0)		#a[202] = r3			\n\
	pxor	1600(%0), %%mm2		#r2 ^= a[200]	204->200	\n\
	pxor	%%mm3, %%mm2		#r2 ^= r3	202->200	\n\
	pxor	1608(%0), %%mm2		#r2 ^= a[201]	201->200	\n\
	movq	%%mm2, 1600(%0)		#a[200] = r2			\n\
	movq	1592(%0), %%mm3		#r3 = a[199]	199->198	\n\
	pxor	1584(%0), %%mm3		#r3 ^= a[198]			\n\
	movq	%%mm3, 1584(%0)		#a[198] = r3			\n\
	pxor	1568(%0), %%mm3		#r3 ^= a[196]	198->196	\n\
	pxor	1576(%0), %%mm3		#r3 ^= a[197]	197->196	\n\
	movq	%%mm3, 1568(%0)		#a[196] = r3			\n\
	movq	1560(%0), %%mm4		#r4 = a[195]	195->194	\n\
	pxor	1552(%0), %%mm4		#r4 ^= a[194]			\n\
	movq	%%mm4, 1552(%0)		#a[194] = r4			\n\
	pxor	1536(%0), %%mm0		#r0 ^= a[192]	224->192	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	208->192	\n\
	pxor	%%mm2, %%mm0		#r0 ^= r2	200->192	\n\
	pxor	%%mm3, %%mm0		#r0 ^= r3	196->192	\n\
	pxor	%%mm4, %%mm0		#r0 ^= r4	194->192	\n\
	movq	%%mm0, 1536(%0)		#a[192] = r0			\n\
	pxor	1544(%0), %%mm0		#r0 ^= a[193]	193->192	\n\
	movq	%%mm0, 1536(%0)		#a[192] = r0			\n\
	movq	1528(%0), %%mm0		#r0 = a[191]	191->190	\n\
	pxor	1520(%0), %%mm0		#r0 ^= a[190]			\n\
	movq	%%mm0, 1520(%0)		#a[190] = r0			\n\
	pxor	1504(%0), %%mm0		#r0 ^= a[188]	190->188	\n\
	pxor	1512(%0), %%mm0		#r0 ^= a[189]	189->188	\n\
	movq	%%mm0, 1504(%0)		#a[188] = r0			\n\
	movq	1496(%0), %%mm1		#r1 = a[187]	187->186	\n\
	pxor	1488(%0), %%mm1		#r1 ^= a[186]			\n\
	movq	%%mm1, 1488(%0)		#a[186] = r1			\n\
	pxor	1472(%0), %%mm0		#r0 ^= a[184]	188->184	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	186->184	\n\
	pxor	1480(%0), %%mm0		#r0 ^= a[185]	185->184	\n\
	movq	%%mm0, 1472(%0)		#a[184] = r0			\n\
	movq	1464(%0), %%mm1		#r1 = a[183]	183->182	\n\
	pxor	1456(%0), %%mm1		#r1 ^= a[182]			\n\
	movq	%%mm1, 1456(%0)		#a[182] = r1			\n\
	pxor	1440(%0), %%mm1		#r1 ^= a[180]	182->180	\n\
	pxor	1448(%0), %%mm1		#r1 ^= a[181]	181->180	\n\
	movq	%%mm1, 1440(%0)		#a[180] = r1			\n\
	movq	1432(%0), %%mm2		#r2 = a[179]	179->178	\n\
	pxor	1424(%0), %%mm2		#r2 ^= a[178]			\n\
	movq	%%mm2, 1424(%0)		#a[178] = r2			\n\
	pxor	1408(%0), %%mm0		#r0 ^= a[176]	184->176	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	180->176	\n\
	pxor	%%mm2, %%mm0		#r0 ^= r2	178->176	\n\
	pxor	1416(%0), %%mm0		#r0 ^= a[177]	177->176	\n\
	movq	%%mm0, 1408(%0)		#a[176] = r0			\n\
	movq	1400(%0), %%mm1		#r1 = a[175]	175->174	\n\
	pxor	1392(%0), %%mm1		#r1 ^= a[174]			\n\
	movq	%%mm1, 1392(%0)		#a[174] = r1			\n\
	pxor	1376(%0), %%mm1		#r1 ^= a[172]	174->172	\n\
	pxor	1384(%0), %%mm1		#r1 ^= a[173]	173->172	\n\
	movq	%%mm1, 1376(%0)		#a[172] = r1			\n\
	movq	1368(%0), %%mm2		#r2 = a[171]	171->170	\n\
	pxor	1360(%0), %%mm2		#r2 ^= a[170]			\n\
	movq	%%mm2, 1360(%0)		#a[170] = r2			\n\
	pxor	1344(%0), %%mm1		#r1 ^= a[168]	172->168	\n\
	pxor	%%mm2, %%mm1		#r1 ^= r2	170->168	\n\
	pxor	1352(%0), %%mm1		#r1 ^= a[169]	169->168	\n\
	movq	%%mm1, 1344(%0)		#a[168] = r1			\n\
	movq	1336(%0), %%mm2		#r2 = a[167]	167->166	\n\
	pxor	1328(%0), %%mm2		#r2 ^= a[166]			\n\
	movq	%%mm2, 1328(%0)		#a[166] = r2			\n\
	pxor	1312(%0), %%mm2		#r2 ^= a[164]	166->164	\n\
	pxor	1320(%0), %%mm2		#r2 ^= a[165]	165->164	\n\
	movq	%%mm2, 1312(%0)		#a[164] = r2			\n\
	movq	1304(%0), %%mm3		#r3 = a[163]	163->162	\n\
	pxor	1296(%0), %%mm3		#r3 ^= a[162]			\n\
	movq	%%mm3, 1296(%0)		#a[162] = r3			\n\
	pxor	1280(%0), %%mm0		#r0 ^= a[160]	176->160	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	168->160	\n\
	pxor	%%mm2, %%mm0		#r0 ^= r2	164->160	\n\
	pxor	%%mm3, %%mm0		#r0 ^= r3	162->160	\n\
	movq	%%mm0, 1280(%0)		#a[160] = r0			\n\
	pxor	1288(%0), %%mm0		#r0 ^= a[161]	161->160	\n\
	movq	%%mm0, 1280(%0)		#a[160] = r0			\n\
	movq	1272(%0), %%mm0		#r0 = a[159]	159->158	\n\
	pxor	1264(%0), %%mm0		#r0 ^= a[158]			\n\
	movq	%%mm0, 1264(%0)		#a[158] = r0			\n\
	pxor	1248(%0), %%mm0		#r0 ^= a[156]	158->156	\n\
	pxor	1256(%0), %%mm0		#r0 ^= a[157]	157->156	\n\
	movq	%%mm0, 1248(%0)		#a[156] = r0			\n\
	movq	1240(%0), %%mm1		#r1 = a[155]	155->154	\n\
	pxor	1232(%0), %%mm1		#r1 ^= a[154]			\n\
	movq	%%mm1, 1232(%0)		#a[154] = r1			\n\
	pxor	1216(%0), %%mm0		#r0 ^= a[152]	156->152	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	154->152	\n\
	pxor	1224(%0), %%mm0		#r0 ^= a[153]	153->152	\n\
	movq	%%mm0, 1216(%0)		#a[152] = r0			\n\
	movq	1208(%0), %%mm1		#r1 = a[151]	151->150	\n\
	pxor	1200(%0), %%mm1		#r1 ^= a[150]			\n\
	movq	%%mm1, 1200(%0)		#a[150] = r1			\n\
	pxor	1184(%0), %%mm1		#r1 ^= a[148]	150->148	\n\
	pxor	1192(%0), %%mm1		#r1 ^= a[149]	149->148	\n\
	movq	%%mm1, 1184(%0)		#a[148] = r1			\n\
	movq	1176(%0), %%mm2		#r2 = a[147]	147->146	\n\
	pxor	1168(%0), %%mm2		#r2 ^= a[146]			\n\
	movq	%%mm2, 1168(%0)		#a[146] = r2			\n\
	pxor	1152(%0), %%mm0		#r0 ^= a[144]	152->144	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	148->144	\n\
	pxor	%%mm2, %%mm0		#r0 ^= r2	146->144	\n\
	movq	%%mm0, 1152(%0)		#a[144] = r0			\n\
	pxor	1160(%0), %%mm0		#r0 ^= a[145]	145->144	\n\
	movq	%%mm0, 1152(%0)		#a[144] = r0			\n\
	movq	1144(%0), %%mm0		#r0 = a[143]	143->142	\n\
	pxor	1136(%0), %%mm0		#r0 ^= a[142]			\n\
	movq	%%mm0, 1136(%0)		#a[142] = r0			\n\
	pxor	1120(%0), %%mm0		#r0 ^= a[140]	142->140	\n\
	pxor	1128(%0), %%mm0		#r0 ^= a[141]	141->140	\n\
	movq	%%mm0, 1120(%0)		#a[140] = r0			\n\
	movq	1112(%0), %%mm1		#r1 = a[139]	139->138	\n\
	pxor	1104(%0), %%mm1		#r1 ^= a[138]			\n\
	movq	%%mm1, 1104(%0)		#a[138] = r1			\n\
	pxor	1088(%0), %%mm0		#r0 ^= a[136]	140->136	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	138->136	\n\
	movq	%%mm0, 1088(%0)		#a[136] = r0			\n\
	pxor	1096(%0), %%mm0		#r0 ^= a[137]	137->136	\n\
	movq	%%mm0, 1088(%0)		#a[136] = r0			\n\
	movq	1080(%0), %%mm0		#r0 = a[135]	135->134	\n\
	pxor	1072(%0), %%mm0		#r0 ^= a[134]			\n\
	movq	%%mm0, 1072(%0)		#a[134] = r0			\n\
	pxor	1056(%0), %%mm0		#r0 ^= a[132]	134->132	\n\
	movq	%%mm0, 1056(%0)		#a[132] = r0			\n\
	pxor	1064(%0), %%mm0		#r0 ^= a[133]	133->132	\n\
	movq	%%mm0, 1056(%0)		#a[132] = r0			\n\
	movq	1048(%0), %%mm0		#r0 = a[131]	131->130	\n\
	pxor	1040(%0), %%mm0		#r0 ^= a[130]			\n\
	movq	%%mm0, 1040(%0)		#a[130] = r0			\n\
	movq	1016(%0), %%mm0		#r0 = a[127]	127->126	\n\
	pxor	1008(%0), %%mm0		#r0 ^= a[126]			\n\
	movq	%%mm0, 1008(%0)		#a[126] = r0			\n\
	pxor	992(%0), %%mm0		#r0 ^= a[124]	126->124	\n\
	pxor	1000(%0), %%mm0		#r0 ^= a[125]	125->124	\n\
	movq	%%mm0, 992(%0)		#a[124] = r0			\n\
	movq	984(%0), %%mm1		#r1 = a[123]	123->122	\n\
	pxor	976(%0), %%mm1		#r1 ^= a[122]			\n\
	movq	%%mm1, 976(%0)		#a[122] = r1			\n\
	pxor	960(%0), %%mm0		#r0 ^= a[120]	124->120	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	122->120	\n\
	pxor	968(%0), %%mm0		#r0 ^= a[121]	121->120	\n\
	movq	%%mm0, 960(%0)		#a[120] = r0			\n\
	movq	952(%0), %%mm1		#r1 = a[119]	119->118	\n\
	pxor	944(%0), %%mm1		#r1 ^= a[118]			\n\
	movq	%%mm1, 944(%0)		#a[118] = r1			\n\
	pxor	928(%0), %%mm1		#r1 ^= a[116]	118->116	\n\
	pxor	936(%0), %%mm1		#r1 ^= a[117]	117->116	\n\
	movq	%%mm1, 928(%0)		#a[116] = r1			\n\
	movq	920(%0), %%mm2		#r2 = a[115]	115->114	\n\
	pxor	912(%0), %%mm2		#r2 ^= a[114]			\n\
	movq	%%mm2, 912(%0)		#a[114] = r2			\n\
	pxor	896(%0), %%mm0		#r0 ^= a[112]	120->112	\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	116->112	\n\
	pxor	%%mm2, %%mm0		#r0 ^= r2	114->112	\n\
	pxor	904(%0), %%mm0		#r0 ^= a[113]	113->112	\n\
	movq	%%mm0, 896(%0)		#a[112] = r0			\n\
	movq	888(%0), %%mm1		#r1 = a[111]	111->110	\n\
	pxor	880(%0), %%mm1		#r1 ^= a[110]			\n\
	movq	%%mm1, 880(%0)		#a[110] = r1			\n\
	pxor	864(%0), %%mm1		#r1 ^= a[108]	110->108	\n\
	pxor	872(%0), %%mm1		#r1 ^= a[109]	109->108	\n\
	movq	%%mm1, 864(%0)		#a[108] = r1			\n\
	movq	856(%0), %%mm2		#r2 = a[107]	107->106	\n\
	pxor	848(%0), %%mm2		#r2 ^= a[106]			\n\
	movq	%%mm2, 848(%0)		#a[106] = r2			\n\
	pxor	832(%0), %%mm1		#r1 ^= a[104]	108->104	\n\
	pxor	%%mm2, %%mm1		#r1 ^= r2	106->104	\n\
	pxor	840(%0), %%mm1		#r1 ^= a[105]	105->104	\n\
	movq	%%mm1, 832(%0)		#a[104] = r1			\n\
	movq	824(%0), %%mm2		#r2 = a[103]	103->102	\n\
	pxor	816(%0), %%mm2		#r2 ^= a[102]			\n\
	movq	%%mm2, 816(%0)		#a[102] = r2			\n\
	pxor	800(%0), %%mm2		#r2 ^= a[100]	102->100	\n\
	pxor	808(%0), %%mm2		#r2 ^= a[101]	101->100	\n\
	movq	%%mm2, 800(%0)		#a[100] = r2			\n\
	movq	792(%0), %%mm3		#r3 = a[99]	99->98		\n\
	pxor	784(%0), %%mm3		#r3 ^= a[98]			\n\
	movq	%%mm3, 784(%0)		#a[98] = r3			\n\
	pxor	768(%0), %%mm0		#r0 ^= a[96]	112->96		\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	104->96		\n\
	pxor	%%mm2, %%mm0		#r0 ^= r2	100->96		\n\
	pxor	%%mm3, %%mm0		#r0 ^= r3	98->96		\n\
	movq	%%mm0, 768(%0)		#a[96] = r0			\n\
	pxor	776(%0), %%mm0		#r0 ^= a[97]	97->96		\n\
	movq	%%mm0, 768(%0)		#a[96] = r0			\n\
	movq	760(%0), %%mm0		#r0 = a[95]	95->94		\n\
	pxor	752(%0), %%mm0		#r0 ^= a[94]			\n\
	movq	%%mm0, 752(%0)		#a[94] = r0			\n\
	pxor	736(%0), %%mm0		#r0 ^= a[92]	94->92		\n\
	pxor	744(%0), %%mm0		#r0 ^= a[93]	93->92		\n\
	movq	%%mm0, 736(%0)		#a[92] = r0			\n\
	movq	728(%0), %%mm1		#r1 = a[91]	91->90		\n\
	pxor	720(%0), %%mm1		#r1 ^= a[90]			\n\
	movq	%%mm1, 720(%0)		#a[90] = r1			\n\
	pxor	704(%0), %%mm0		#r0 ^= a[88]	92->88		\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	90->88		\n\
	pxor	712(%0), %%mm0		#r0 ^= a[89]	89->88		\n\
	movq	%%mm0, 704(%0)		#a[88] = r0			\n\
	movq	696(%0), %%mm1		#r1 = a[87]	87->86		\n\
	pxor	688(%0), %%mm1		#r1 ^= a[86]			\n\
	movq	%%mm1, 688(%0)		#a[86] = r1			\n\
	pxor	672(%0), %%mm1		#r1 ^= a[84]	86->84		\n\
	pxor	680(%0), %%mm1		#r1 ^= a[85]	85->84		\n\
	movq	%%mm1, 672(%0)		#a[84] = r1			\n\
	movq	664(%0), %%mm2		#r2 = a[83]	83->82		\n\
	pxor	656(%0), %%mm2		#r2 ^= a[82]			\n\
	movq	%%mm2, 656(%0)		#a[82] = r2			\n\
	pxor	640(%0), %%mm0		#r0 ^= a[80]	88->80		\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	84->80		\n\
	pxor	%%mm2, %%mm0		#r0 ^= r2	82->80		\n\
	movq	%%mm0, 640(%0)		#a[80] = r0			\n\
	pxor	648(%0), %%mm0		#r0 ^= a[81]	81->80		\n\
	movq	%%mm0, 640(%0)		#a[80] = r0			\n\
	movq	632(%0), %%mm0		#r0 = a[79]	79->78		\n\
	pxor	624(%0), %%mm0		#r0 ^= a[78]			\n\
	movq	%%mm0, 624(%0)		#a[78] = r0			\n\
	pxor	608(%0), %%mm0		#r0 ^= a[76]	78->76		\n\
	pxor	616(%0), %%mm0		#r0 ^= a[77]	77->76		\n\
	movq	%%mm0, 608(%0)		#a[76] = r0			\n\
	movq	600(%0), %%mm1		#r1 = a[75]	75->74		\n\
	pxor	592(%0), %%mm1		#r1 ^= a[74]			\n\
	movq	%%mm1, 592(%0)		#a[74] = r1			\n\
	pxor	576(%0), %%mm0		#r0 ^= a[72]	76->72		\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	74->72		\n\
	movq	%%mm0, 576(%0)		#a[72] = r0			\n\
	pxor	584(%0), %%mm0		#r0 ^= a[73]	73->72		\n\
	movq	%%mm0, 576(%0)		#a[72] = r0			\n\
	movq	568(%0), %%mm0		#r0 = a[71]	71->70		\n\
	pxor	560(%0), %%mm0		#r0 ^= a[70]			\n\
	movq	%%mm0, 560(%0)		#a[70] = r0			\n\
	pxor	544(%0), %%mm0		#r0 ^= a[68]	70->68		\n\
	movq	%%mm0, 544(%0)		#a[68] = r0			\n\
	pxor	552(%0), %%mm0		#r0 ^= a[69]	69->68		\n\
	movq	%%mm0, 544(%0)		#a[68] = r0			\n\
	movq	536(%0), %%mm0		#r0 = a[67]	67->66		\n\
	pxor	528(%0), %%mm0		#r0 ^= a[66]			\n\
	movq	%%mm0, 528(%0)		#a[66] = r0			\n\
	movq	504(%0), %%mm0		#r0 = a[63]	63->62		\n\
	pxor	496(%0), %%mm0		#r0 ^= a[62]			\n\
	movq	%%mm0, 496(%0)		#a[62] = r0			\n\
	pxor	480(%0), %%mm0		#r0 ^= a[60]	62->60		\n\
	pxor	488(%0), %%mm0		#r0 ^= a[61]	61->60		\n\
	movq	%%mm0, 480(%0)		#a[60] = r0			\n\
	movq	472(%0), %%mm1		#r1 = a[59]	59->58		\n\
	pxor	464(%0), %%mm1		#r1 ^= a[58]			\n\
	movq	%%mm1, 464(%0)		#a[58] = r1			\n\
	pxor	448(%0), %%mm0		#r0 ^= a[56]	60->56		\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	58->56		\n\
	pxor	456(%0), %%mm0		#r0 ^= a[57]	57->56		\n\
	movq	%%mm0, 448(%0)		#a[56] = r0			\n\
	movq	440(%0), %%mm1		#r1 = a[55]	55->54		\n\
	pxor	432(%0), %%mm1		#r1 ^= a[54]			\n\
	movq	%%mm1, 432(%0)		#a[54] = r1			\n\
	pxor	416(%0), %%mm1		#r1 ^= a[52]	54->52		\n\
	pxor	424(%0), %%mm1		#r1 ^= a[53]	53->52		\n\
	movq	%%mm1, 416(%0)		#a[52] = r1			\n\
	movq	408(%0), %%mm2		#r2 = a[51]	51->50		\n\
	pxor	400(%0), %%mm2		#r2 ^= a[50]			\n\
	movq	%%mm2, 400(%0)		#a[50] = r2			\n\
	pxor	384(%0), %%mm0		#r0 ^= a[48]	56->48		\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	52->48		\n\
	pxor	%%mm2, %%mm0		#r0 ^= r2	50->48		\n\
	movq	%%mm0, 384(%0)		#a[48] = r0			\n\
	pxor	392(%0), %%mm0		#r0 ^= a[49]	49->48		\n\
	movq	%%mm0, 384(%0)		#a[48] = r0			\n\
	movq	376(%0), %%mm0		#r0 = a[47]	47->46		\n\
	pxor	368(%0), %%mm0		#r0 ^= a[46]			\n\
	movq	%%mm0, 368(%0)		#a[46] = r0			\n\
	pxor	352(%0), %%mm0		#r0 ^= a[44]	46->44		\n\
	pxor	360(%0), %%mm0		#r0 ^= a[45]	45->44		\n\
	movq	%%mm0, 352(%0)		#a[44] = r0			\n\
	movq	344(%0), %%mm1		#r1 = a[43]	43->42		\n\
	pxor	336(%0), %%mm1		#r1 ^= a[42]			\n\
	movq	%%mm1, 336(%0)		#a[42] = r1			\n\
	pxor	320(%0), %%mm0		#r0 ^= a[40]	44->40		\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	42->40		\n\
	movq	%%mm0, 320(%0)		#a[40] = r0			\n\
	pxor	328(%0), %%mm0		#r0 ^= a[41]	41->40		\n\
	movq	%%mm0, 320(%0)		#a[40] = r0			\n\
	movq	312(%0), %%mm0		#r0 = a[39]	39->38		\n\
	pxor	304(%0), %%mm0		#r0 ^= a[38]			\n\
	movq	%%mm0, 304(%0)		#a[38] = r0			\n\
	pxor	288(%0), %%mm0		#r0 ^= a[36]	38->36		\n\
	movq	%%mm0, 288(%0)		#a[36] = r0			\n\
	pxor	296(%0), %%mm0		#r0 ^= a[37]	37->36		\n\
	movq	%%mm0, 288(%0)		#a[36] = r0			\n\
	movq	280(%0), %%mm0		#r0 = a[35]	35->34		\n\
	pxor	272(%0), %%mm0		#r0 ^= a[34]			\n\
	movq	%%mm0, 272(%0)		#a[34] = r0			\n\
	movq	248(%0), %%mm0		#r0 = a[31]	31->30		\n\
	pxor	240(%0), %%mm0		#r0 ^= a[30]			\n\
	movq	%%mm0, 240(%0)		#a[30] = r0			\n\
	pxor	224(%0), %%mm0		#r0 ^= a[28]	30->28		\n\
	pxor	232(%0), %%mm0		#r0 ^= a[29]	29->28		\n\
	movq	%%mm0, 224(%0)		#a[28] = r0			\n\
	movq	216(%0), %%mm1		#r1 = a[27]	27->26		\n\
	pxor	208(%0), %%mm1		#r1 ^= a[26]			\n\
	movq	%%mm1, 208(%0)		#a[26] = r1			\n\
	pxor	192(%0), %%mm0		#r0 ^= a[24]	28->24		\n\
	pxor	%%mm1, %%mm0		#r0 ^= r1	26->24		\n\
	movq	%%mm0, 192(%0)		#a[24] = r0			\n\
	pxor	200(%0), %%mm0		#r0 ^= a[25]	25->24		\n\
	movq	%%mm0, 192(%0)		#a[24] = r0			\n\
	movq	184(%0), %%mm0		#r0 = a[23]	23->22		\n\
	pxor	176(%0), %%mm0		#r0 ^= a[22]			\n\
	movq	%%mm0, 176(%0)		#a[22] = r0			\n\
	pxor	160(%0), %%mm0		#r0 ^= a[20]	22->20		\n\
	movq	%%mm0, 160(%0)		#a[20] = r0			\n\
	pxor	168(%0), %%mm0		#r0 ^= a[21]	21->20		\n\
	movq	%%mm0, 160(%0)		#a[20] = r0			\n\
	movq	152(%0), %%mm0		#r0 = a[19]	19->18		\n\
	pxor	144(%0), %%mm0		#r0 ^= a[18]			\n\
	movq	%%mm0, 144(%0)		#a[18] = r0			\n\
	movq	120(%0), %%mm0		#r0 = a[15]	15->14		\n\
	pxor	112(%0), %%mm0		#r0 ^= a[14]			\n\
	movq	%%mm0, 112(%0)		#a[14] = r0			\n\
	pxor	96(%0), %%mm0		#r0 ^= a[12]	14->12		\n\
	movq	%%mm0, 96(%0)		#a[12] = r0			\n\
	pxor	104(%0), %%mm0		#r0 ^= a[13]	13->12		\n\
	movq	%%mm0, 96(%0)		#a[12] = r0			\n\
	movq	88(%0), %%mm0		#r0 = a[11]	11->10		\n\
	pxor	80(%0), %%mm0		#r0 ^= a[10]			\n\
	movq	%%mm0, 80(%0)		#a[10] = r0			\n\
	movq	56(%0), %%mm0		#r0 = a[7]	7->6		\n\
	pxor	48(%0), %%mm0		#r0 ^= a[6]			\n\
	movq	%%mm0, 48(%0)		#a[6] = r0			\n\
	movq	8(%0), %%mm0		#r0 = a[1]			\n\
	pxor	24(%0), %%mm0		#r0 ^= a[3]			\n\
	pxor	40(%0), %%mm0		#r0 ^= a[5]			\n\
	pxor	56(%0), %%mm0		#r0 ^= a[7]			\n\
	pxor	72(%0), %%mm0		#r0 ^= a[9]			\n\
	pxor	88(%0), %%mm0		#r0 ^= a[11]			\n\
	pxor	104(%0), %%mm0		#r0 ^= a[13]			\n\
	pxor	120(%0), %%mm0		#r0 ^= a[15]			\n\
	pxor	136(%0), %%mm0		#r0 ^= a[17]			\n\
	pxor	152(%0), %%mm0		#r0 ^= a[19]			\n\
	pxor	168(%0), %%mm0		#r0 ^= a[21]			\n\
	pxor	184(%0), %%mm0		#r0 ^= a[23]			\n\
	pxor	200(%0), %%mm0		#r0 ^= a[25]			\n\
	pxor	216(%0), %%mm0		#r0 ^= a[27]			\n\
	pxor	232(%0), %%mm0		#r0 ^= a[29]			\n\
	pxor	248(%0), %%mm0		#r0 ^= a[31]			\n\
	pxor	264(%0), %%mm0		#r0 ^= a[33]			\n\
	pxor	280(%0), %%mm0		#r0 ^= a[35]			\n\
	pxor	296(%0), %%mm0		#r0 ^= a[37]			\n\
	pxor	312(%0), %%mm0		#r0 ^= a[39]			\n\
	pxor	328(%0), %%mm0		#r0 ^= a[41]			\n\
	pxor	344(%0), %%mm0		#r0 ^= a[43]			\n\
	pxor	360(%0), %%mm0		#r0 ^= a[45]			\n\
	pxor	376(%0), %%mm0		#r0 ^= a[47]			\n\
	pxor	392(%0), %%mm0		#r0 ^= a[49]			\n\
	pxor	408(%0), %%mm0		#r0 ^= a[51]			\n\
	pxor	424(%0), %%mm0		#r0 ^= a[53]			\n\
	pxor	440(%0), %%mm0		#r0 ^= a[55]			\n\
	pxor	456(%0), %%mm0		#r0 ^= a[57]			\n\
	pxor	472(%0), %%mm0		#r0 ^= a[59]			\n\
	pxor	488(%0), %%mm0		#r0 ^= a[61]			\n\
	pxor	504(%0), %%mm0		#r0 ^= a[63]			\n\
	pxor	520(%0), %%mm0		#r0 ^= a[65]			\n\
	pxor	536(%0), %%mm0		#r0 ^= a[67]			\n\
	pxor	552(%0), %%mm0		#r0 ^= a[69]			\n\
	pxor	568(%0), %%mm0		#r0 ^= a[71]			\n\
	pxor	584(%0), %%mm0		#r0 ^= a[73]			\n\
	pxor	600(%0), %%mm0		#r0 ^= a[75]			\n\
	pxor	616(%0), %%mm0		#r0 ^= a[77]			\n\
	pxor	632(%0), %%mm0		#r0 ^= a[79]			\n\
	pxor	648(%0), %%mm0		#r0 ^= a[81]			\n\
	pxor	664(%0), %%mm0		#r0 ^= a[83]			\n\
	pxor	680(%0), %%mm0		#r0 ^= a[85]			\n\
	pxor	696(%0), %%mm0		#r0 ^= a[87]			\n\
	pxor	712(%0), %%mm0		#r0 ^= a[89]			\n\
	pxor	728(%0), %%mm0		#r0 ^= a[91]			\n\
	pxor	744(%0), %%mm0		#r0 ^= a[93]			\n\
	pxor	760(%0), %%mm0		#r0 ^= a[95]			\n\
	pxor	776(%0), %%mm0		#r0 ^= a[97]			\n\
	pxor	792(%0), %%mm0		#r0 ^= a[99]			\n\
	pxor	808(%0), %%mm0		#r0 ^= a[101]			\n\
	pxor	824(%0), %%mm0		#r0 ^= a[103]			\n\
	pxor	840(%0), %%mm0		#r0 ^= a[105]			\n\
	pxor	856(%0), %%mm0		#r0 ^= a[107]			\n\
	pxor	872(%0), %%mm0		#r0 ^= a[109]			\n\
	pxor	888(%0), %%mm0		#r0 ^= a[111]			\n\
	pxor	904(%0), %%mm0		#r0 ^= a[113]			\n\
	pxor	920(%0), %%mm0		#r0 ^= a[115]			\n\
	pxor	936(%0), %%mm0		#r0 ^= a[117]			\n\
	pxor	952(%0), %%mm0		#r0 ^= a[119]			\n\
	pxor	968(%0), %%mm0		#r0 ^= a[121]			\n\
	pxor	984(%0), %%mm0		#r0 ^= a[123]			\n\
	pxor	1000(%0), %%mm0		#r0 ^= a[125]			\n\
	pxor	1016(%0), %%mm0		#r0 ^= a[127]			\n\
	pxor	1032(%0), %%mm0		#r0 ^= a[129]			\n\
	pxor	1048(%0), %%mm0		#r0 ^= a[131]			\n\
	pxor	1064(%0), %%mm0		#r0 ^= a[133]			\n\
	pxor	1080(%0), %%mm0		#r0 ^= a[135]			\n\
	pxor	1096(%0), %%mm0		#r0 ^= a[137]			\n\
	pxor	1112(%0), %%mm0		#r0 ^= a[139]			\n\
	pxor	1128(%0), %%mm0		#r0 ^= a[141]			\n\
	pxor	1144(%0), %%mm0		#r0 ^= a[143]			\n\
	pxor	1160(%0), %%mm0		#r0 ^= a[145]			\n\
	pxor	1176(%0), %%mm0		#r0 ^= a[147]			\n\
	pxor	1192(%0), %%mm0		#r0 ^= a[149]			\n\
	pxor	1208(%0), %%mm0		#r0 ^= a[151]			\n\
	pxor	1224(%0), %%mm0		#r0 ^= a[153]			\n\
	pxor	1240(%0), %%mm0		#r0 ^= a[155]			\n\
	pxor	1256(%0), %%mm0		#r0 ^= a[157]			\n\
	pxor	1272(%0), %%mm0		#r0 ^= a[159]			\n\
	pxor	1288(%0), %%mm0		#r0 ^= a[161]			\n\
	pxor	1304(%0), %%mm0		#r0 ^= a[163]			\n\
	pxor	1320(%0), %%mm0		#r0 ^= a[165]			\n\
	pxor	1336(%0), %%mm0		#r0 ^= a[167]			\n\
	pxor	1352(%0), %%mm0		#r0 ^= a[169]			\n\
	pxor	1368(%0), %%mm0		#r0 ^= a[171]			\n\
	pxor	1384(%0), %%mm0		#r0 ^= a[173]			\n\
	pxor	1400(%0), %%mm0		#r0 ^= a[175]			\n\
	pxor	1416(%0), %%mm0		#r0 ^= a[177]			\n\
	pxor	1432(%0), %%mm0		#r0 ^= a[179]			\n\
	pxor	1448(%0), %%mm0		#r0 ^= a[181]			\n\
	pxor	1464(%0), %%mm0		#r0 ^= a[183]			\n\
	pxor	1480(%0), %%mm0		#r0 ^= a[185]			\n\
	pxor	1496(%0), %%mm0		#r0 ^= a[187]			\n\
	pxor	1512(%0), %%mm0		#r0 ^= a[189]			\n\
	pxor	1528(%0), %%mm0		#r0 ^= a[191]			\n\
	pxor	1544(%0), %%mm0		#r0 ^= a[193]			\n\
	pxor	1560(%0), %%mm0		#r0 ^= a[195]			\n\
	pxor	1576(%0), %%mm0		#r0 ^= a[197]			\n\
	pxor	1592(%0), %%mm0		#r0 ^= a[199]			\n\
	pxor	1608(%0), %%mm0		#r0 ^= a[201]			\n\
	pxor	1624(%0), %%mm0		#r0 ^= a[203]			\n\
	pxor	1640(%0), %%mm0		#r0 ^= a[205]			\n\
	pxor	1656(%0), %%mm0		#r0 ^= a[207]			\n\
	pxor	1672(%0), %%mm0		#r0 ^= a[209]			\n\
	pxor	1688(%0), %%mm0		#r0 ^= a[211]			\n\
	pxor	1704(%0), %%mm0		#r0 ^= a[213]			\n\
	pxor	1720(%0), %%mm0		#r0 ^= a[215]			\n\
	pxor	1736(%0), %%mm0		#r0 ^= a[217]			\n\
	pxor	1752(%0), %%mm0		#r0 ^= a[219]			\n\
	pxor	1768(%0), %%mm0		#r0 ^= a[221]			\n\
	pxor	1784(%0), %%mm0		#r0 ^= a[223]			\n\
	pxor	1800(%0), %%mm0		#r0 ^= a[225]			\n\
	pxor	1816(%0), %%mm0		#r0 ^= a[227]			\n\
	pxor	1832(%0), %%mm0		#r0 ^= a[229]			\n\
	pxor	1848(%0), %%mm0		#r0 ^= a[231]			\n\
	pxor	1864(%0), %%mm0		#r0 ^= a[233]			\n\
	pxor	1880(%0), %%mm0		#r0 ^= a[235]			\n\
	pxor	1896(%0), %%mm0		#r0 ^= a[237]			\n\
	pxor	1912(%0), %%mm0		#r0 ^= a[239]			\n\
	pxor	1928(%0), %%mm0		#r0 ^= a[241]			\n\
	pxor	1944(%0), %%mm0		#r0 ^= a[243]			\n\
	pxor	1960(%0), %%mm0		#r0 ^= a[245]			\n\
	pxor	1976(%0), %%mm0		#r0 ^= a[247]			\n\
	pxor	1992(%0), %%mm0		#r0 ^= a[249]			\n\
	pxor	2008(%0), %%mm0		#r0 ^= a[251]			\n\
	pxor	2024(%0), %%mm0		#r0 ^= a[253]			\n\
	pxor	2040(%0), %%mm0		#r0 ^= a[255]			\n\
	movq	%%mm0, 0(%1)		#b[0] = r0			\n\
	movq	16(%0), %%mm0		#r0 = a[2]			\n\
	pxor	24(%0), %%mm0		#r0 ^= a[3]			\n\
	pxor	48(%0), %%mm0		#r0 ^= a[6]			\n\
	pxor	80(%0), %%mm0		#r0 ^= a[10]			\n\
	pxor	112(%0), %%mm0		#r0 ^= a[14]			\n\
	pxor	144(%0), %%mm0		#r0 ^= a[18]			\n\
	pxor	176(%0), %%mm0		#r0 ^= a[22]			\n\
	pxor	208(%0), %%mm0		#r0 ^= a[26]			\n\
	pxor	240(%0), %%mm0		#r0 ^= a[30]			\n\
	pxor	272(%0), %%mm0		#r0 ^= a[34]			\n\
	pxor	304(%0), %%mm0		#r0 ^= a[38]			\n\
	pxor	336(%0), %%mm0		#r0 ^= a[42]			\n\
	pxor	368(%0), %%mm0		#r0 ^= a[46]			\n\
	pxor	400(%0), %%mm0		#r0 ^= a[50]			\n\
	pxor	432(%0), %%mm0		#r0 ^= a[54]			\n\
	pxor	464(%0), %%mm0		#r0 ^= a[58]			\n\
	pxor	496(%0), %%mm0		#r0 ^= a[62]			\n\
	pxor	528(%0), %%mm0		#r0 ^= a[66]			\n\
	pxor	560(%0), %%mm0		#r0 ^= a[70]			\n\
	pxor	592(%0), %%mm0		#r0 ^= a[74]			\n\
	pxor	624(%0), %%mm0		#r0 ^= a[78]			\n\
	pxor	656(%0), %%mm0		#r0 ^= a[82]			\n\
	pxor	688(%0), %%mm0		#r0 ^= a[86]			\n\
	pxor	720(%0), %%mm0		#r0 ^= a[90]			\n\
	pxor	752(%0), %%mm0		#r0 ^= a[94]			\n\
	pxor	784(%0), %%mm0		#r0 ^= a[98]			\n\
	pxor	816(%0), %%mm0		#r0 ^= a[102]			\n\
	pxor	848(%0), %%mm0		#r0 ^= a[106]			\n\
	pxor	880(%0), %%mm0		#r0 ^= a[110]			\n\
	pxor	912(%0), %%mm0		#r0 ^= a[114]			\n\
	pxor	944(%0), %%mm0		#r0 ^= a[118]			\n\
	pxor	976(%0), %%mm0		#r0 ^= a[122]			\n\
	pxor	1008(%0), %%mm0		#r0 ^= a[126]			\n\
	pxor	1040(%0), %%mm0		#r0 ^= a[130]			\n\
	pxor	1072(%0), %%mm0		#r0 ^= a[134]			\n\
	pxor	1104(%0), %%mm0		#r0 ^= a[138]			\n\
	pxor	1136(%0), %%mm0		#r0 ^= a[142]			\n\
	pxor	1168(%0), %%mm0		#r0 ^= a[146]			\n\
	pxor	1200(%0), %%mm0		#r0 ^= a[150]			\n\
	pxor	1232(%0), %%mm0		#r0 ^= a[154]			\n\
	pxor	1264(%0), %%mm0		#r0 ^= a[158]			\n\
	pxor	1296(%0), %%mm0		#r0 ^= a[162]			\n\
	pxor	1328(%0), %%mm0		#r0 ^= a[166]			\n\
	pxor	1360(%0), %%mm0		#r0 ^= a[170]			\n\
	pxor	1392(%0), %%mm0		#r0 ^= a[174]			\n\
	pxor	1424(%0), %%mm0		#r0 ^= a[178]			\n\
	pxor	1456(%0), %%mm0		#r0 ^= a[182]			\n\
	pxor	1488(%0), %%mm0		#r0 ^= a[186]			\n\
	pxor	1520(%0), %%mm0		#r0 ^= a[190]			\n\
	pxor	1552(%0), %%mm0		#r0 ^= a[194]			\n\
	pxor	1584(%0), %%mm0		#r0 ^= a[198]			\n\
	pxor	1616(%0), %%mm0		#r0 ^= a[202]			\n\
	pxor	1648(%0), %%mm0		#r0 ^= a[206]			\n\
	pxor	1680(%0), %%mm0		#r0 ^= a[210]			\n\
	pxor	1712(%0), %%mm0		#r0 ^= a[214]			\n\
	pxor	1744(%0), %%mm0		#r0 ^= a[218]			\n\
	pxor	1776(%0), %%mm0		#r0 ^= a[222]			\n\
	pxor	1808(%0), %%mm0		#r0 ^= a[226]			\n\
	pxor	1840(%0), %%mm0		#r0 ^= a[230]			\n\
	pxor	1872(%0), %%mm0		#r0 ^= a[234]			\n\
	pxor	1904(%0), %%mm0		#r0 ^= a[238]			\n\
	pxor	1936(%0), %%mm0		#r0 ^= a[242]			\n\
	pxor	1968(%0), %%mm0		#r0 ^= a[246]			\n\
	pxor	2000(%0), %%mm0		#r0 ^= a[250]			\n\
	pxor	2032(%0), %%mm0		#r0 ^= a[254]			\n\
	movq	%%mm0, 8(%1)		#b[1] = r0			\n\
	movq	32(%0), %%mm0		#r0 = a[4]			\n\
	pxor	40(%0), %%mm0		#r0 ^= a[5]			\n\
	pxor	48(%0), %%mm0		#r0 ^= a[6]			\n\
	pxor	96(%0), %%mm0		#r0 ^= a[12]			\n\
	pxor	160(%0), %%mm0		#r0 ^= a[20]			\n\
	pxor	224(%0), %%mm0		#r0 ^= a[28]			\n\
	pxor	288(%0), %%mm0		#r0 ^= a[36]			\n\
	pxor	352(%0), %%mm0		#r0 ^= a[44]			\n\
	pxor	416(%0), %%mm0		#r0 ^= a[52]			\n\
	pxor	480(%0), %%mm0		#r0 ^= a[60]			\n\
	pxor	544(%0), %%mm0		#r0 ^= a[68]			\n\
	pxor	608(%0), %%mm0		#r0 ^= a[76]			\n\
	pxor	672(%0), %%mm0		#r0 ^= a[84]			\n\
	pxor	736(%0), %%mm0		#r0 ^= a[92]			\n\
	pxor	800(%0), %%mm0		#r0 ^= a[100]			\n\
	pxor	864(%0), %%mm0		#r0 ^= a[108]			\n\
	pxor	928(%0), %%mm0		#r0 ^= a[116]			\n\
	pxor	992(%0), %%mm0		#r0 ^= a[124]			\n\
	pxor	1056(%0), %%mm0		#r0 ^= a[132]			\n\
	pxor	1120(%0), %%mm0		#r0 ^= a[140]			\n\
	pxor	1184(%0), %%mm0		#r0 ^= a[148]			\n\
	pxor	1248(%0), %%mm0		#r0 ^= a[156]			\n\
	pxor	1312(%0), %%mm0		#r0 ^= a[164]			\n\
	pxor	1376(%0), %%mm0		#r0 ^= a[172]			\n\
	pxor	1440(%0), %%mm0		#r0 ^= a[180]			\n\
	pxor	1504(%0), %%mm0		#r0 ^= a[188]			\n\
	pxor	1568(%0), %%mm0		#r0 ^= a[196]			\n\
	pxor	1632(%0), %%mm0		#r0 ^= a[204]			\n\
	pxor	1696(%0), %%mm0		#r0 ^= a[212]			\n\
	pxor	1760(%0), %%mm0		#r0 ^= a[220]			\n\
	pxor	1824(%0), %%mm0		#r0 ^= a[228]			\n\
	pxor	1888(%0), %%mm0		#r0 ^= a[236]			\n\
	pxor	1952(%0), %%mm0		#r0 ^= a[244]			\n\
	pxor	2016(%0), %%mm0		#r0 ^= a[252]			\n\
	movq	%%mm0, 16(%1)		#b[2] = r0			\n\
	movq	64(%0), %%mm0		#r0 = a[8]			\n\
	pxor	72(%0), %%mm0		#r0 ^= a[9]			\n\
	pxor	80(%0), %%mm0		#r0 ^= a[10]			\n\
	pxor	96(%0), %%mm0		#r0 ^= a[12]			\n\
	pxor	192(%0), %%mm0		#r0 ^= a[24]			\n\
	pxor	320(%0), %%mm0		#r0 ^= a[40]			\n\
	pxor	448(%0), %%mm0		#r0 ^= a[56]			\n\
	pxor	576(%0), %%mm0		#r0 ^= a[72]			\n\
	pxor	704(%0), %%mm0		#r0 ^= a[88]			\n\
	pxor	832(%0), %%mm0		#r0 ^= a[104]			\n\
	pxor	960(%0), %%mm0		#r0 ^= a[120]			\n\
	pxor	1088(%0), %%mm0		#r0 ^= a[136]			\n\
	pxor	1216(%0), %%mm0		#r0 ^= a[152]			\n\
	pxor	1344(%0), %%mm0		#r0 ^= a[168]			\n\
	pxor	1472(%0), %%mm0		#r0 ^= a[184]			\n\
	pxor	1600(%0), %%mm0		#r0 ^= a[200]			\n\
	pxor	1728(%0), %%mm0		#r0 ^= a[216]			\n\
	pxor	1856(%0), %%mm0		#r0 ^= a[232]			\n\
	pxor	1984(%0), %%mm0		#r0 ^= a[248]			\n\
	movq	%%mm0, 24(%1)		#b[3] = r0			\n\
	movq	128(%0), %%mm0		#r0 = a[16]			\n\
	pxor	136(%0), %%mm0		#r0 ^= a[17]			\n\
	pxor	144(%0), %%mm0		#r0 ^= a[18]			\n\
	pxor	160(%0), %%mm0		#r0 ^= a[20]			\n\
	pxor	192(%0), %%mm0		#r0 ^= a[24]			\n\
	pxor	384(%0), %%mm0		#r0 ^= a[48]			\n\
	pxor	640(%0), %%mm0		#r0 ^= a[80]			\n\
	pxor	896(%0), %%mm0		#r0 ^= a[112]			\n\
	pxor	1152(%0), %%mm0		#r0 ^= a[144]			\n\
	pxor	1408(%0), %%mm0		#r0 ^= a[176]			\n\
	pxor	1664(%0), %%mm0		#r0 ^= a[208]			\n\
	pxor	1920(%0), %%mm0		#r0 ^= a[240]			\n\
	movq	%%mm0, 32(%1)		#b[4] = r0			\n\
	movq	256(%0), %%mm0		#r0 = a[32]			\n\
	pxor	264(%0), %%mm0		#r0 ^= a[33]			\n\
	pxor	272(%0), %%mm0		#r0 ^= a[34]			\n\
	pxor	288(%0), %%mm0		#r0 ^= a[36]			\n\
	pxor	320(%0), %%mm0		#r0 ^= a[40]			\n\
	pxor	384(%0), %%mm0		#r0 ^= a[48]			\n\
	pxor	768(%0), %%mm0		#r0 ^= a[96]			\n\
	pxor	1280(%0), %%mm0		#r0 ^= a[160]			\n\
	pxor	1792(%0), %%mm0		#r0 ^= a[224]			\n\
	movq	%%mm0, 40(%1)		#b[5] = r0			\n\
	movq	512(%0), %%mm0		#r0 = a[64]			\n\
	pxor	520(%0), %%mm0		#r0 ^= a[65]			\n\
	pxor	528(%0), %%mm0		#r0 ^= a[66]			\n\
	pxor	544(%0), %%mm0		#r0 ^= a[68]			\n\
	pxor	576(%0), %%mm0		#r0 ^= a[72]			\n\
	pxor	640(%0), %%mm0		#r0 ^= a[80]			\n\
	pxor	768(%0), %%mm0		#r0 ^= a[96]			\n\
	pxor	1536(%0), %%mm0		#r0 ^= a[192]			\n\
	movq	%%mm0, 48(%1)		#b[6] = r0			\n\
	movq	1024(%0), %%mm0		#r0 = a[128]			\n\
	pxor	1032(%0), %%mm0		#r0 ^= a[129]			\n\
	pxor	1040(%0), %%mm0		#r0 ^= a[130]			\n\
	pxor	1056(%0), %%mm0		#r0 ^= a[132]			\n\
	pxor	1088(%0), %%mm0		#r0 ^= a[136]			\n\
	pxor	1152(%0), %%mm0		#r0 ^= a[144]			\n\
	pxor	1280(%0), %%mm0		#r0 ^= a[160]			\n\
	pxor	1536(%0), %%mm0		#r0 ^= a[192]			\n\
	movq	%%mm0, 56(%1)		#b[7] = r0			\n\
	emms" : : "r"(mult_w + (i << 5)), "r"(c + i));
    }
  }
}

  
/*********************************************/
void multS(u64 *D, int *S)
/*********************************************/
/* Input: D is a 64x64 matrix,               */
/*        'S' is a subset of the columns,    */
/*        {0,1,2,...,31}.                    */
/* Output: Columns of 'D' not in 'S' are     */
/*         zeroed, and the others unchanged. */
/*********************************************/
{ int i, s;
  u64 mask;

  mask = 0x00000000;
  for (i=0; i<64; i++) {
    s = S[i];
    if ((s>=0) && (s<64))
      mask |= BIT64(s);
  }
  for (i=0; i<64; i++) 
    D[i] &= mask;
}

void mult64x64(u64 *c, u64 *a, u64 *b) {
  asm volatile("\
	movl	%0, %%esi			#a		\n\
	movl	%1, %%edx			#b		\n\
	movl	%2, %%edi			#c		\n\
	leal	8*64(%%esi), %%esi				\n\
	leal	8*64(%%edi), %%edi				\n\
	movl	$-64, %%ecx			#i=0		\n\
1:								\n\
	movl	(%%esi,%%ecx,8), %%eax		#t=(u32)a[i]	\n\
	pxor	%%mm0, %%mm0			#u=0		\n\
	testb	$1, %%al					\n\
	jz	2f						\n\
	movq	(%%edx), %%mm0			#b[0]		\n\
2:								\n\
	.set	n, 1						\n\
	.rept	6				#b[1]..b[6]	\n\
		testb	$1<<n, %%al				\n\
		jz	2f					\n\
		pxor	8*n(%%edx), %%mm0			\n\
2:								\n\
		.set	n, n+1					\n\
	.endr							\n\
	testb	%%al, %%al					\n\
	jns	2f						\n\
	pxor	8*7(%%edx), %%mm0		#b[7]		\n\
2:								\n\
	.set	n, 8						\n\
	.rept	7				#b[8]..b[14]	\n\
		testb	$1<<(n-8), %%ah				\n\
		jz	2f					\n\
		pxor	8*n(%%edx), %%mm0			\n\
2:								\n\
		.set	n, n+1					\n\
	.endr							\n\
	testb	%%ah, %%ah					\n\
	jns	2f						\n\
	pxor	8*15(%%edx), %%mm0		#b[15]		\n\
2:								\n\
	.set	n, 16						\n\
	.rept	15				#b[16]..b[30]	\n\
		testl	$1<<n, %%eax				\n\
		jz	2f					\n\
		pxor	8*n(%%edx), %%mm0			\n\
2:								\n\
		.set	n, n+1					\n\
	.endr							\n\
	testl	%%eax, %%eax					\n\
	jns	2f						\n\
	pxor	8*31(%%edx), %%mm0		#b[31]		\n\
2:								\n\
	movl	4(%%esi,%%ecx,8), %%eax		#t=a[i]>>32	\n\
	.set	n, 32						\n\
	.rept	7				#b[32]..b[38]	\n\
		testb	$1<<(n-32), %%al			\n\
		jz	2f					\n\
		pxor	8*n(%%edx), %%mm0			\n\
2:								\n\
		.set	n, n+1					\n\
	.endr							\n\
	testb	%%al, %%al					\n\
	jns	2f						\n\
	pxor	8*39(%%edx), %%mm0		#b[39]		\n\
2:								\n\
	.set	n, 40						\n\
	.rept	7				#b[40]..b[46]	\n\
		testb	$1<<(n-40), %%ah			\n\
		jz	2f					\n\
		pxor	8*n(%%edx), %%mm0			\n\
2:								\n\
		.set	n, n+1					\n\
	.endr							\n\
	testb	%%ah, %%ah					\n\
	jns	2f						\n\
	pxor	8*47(%%edx), %%mm0		#b[47]		\n\
2:								\n\
	.set	n, 48						\n\
	.rept	15				#b[48]..b[62]	\n\
		testl	$1<<(n-32), %%eax			\n\
		jz	2f					\n\
		pxor	8*n(%%edx), %%mm0			\n\
2:								\n\
		.set	n, n+1					\n\
	.endr							\n\
	testl	%%eax, %%eax					\n\
	jns	2f						\n\
	pxor	8*63(%%edx), %%mm0		#b[63]		\n\
2:								\n\
	movq	%%mm0, (%%edi,%%ecx,8)		#c[i]=u		\n\
	addl	$1, %%ecx			#i++		\n\
	jnz	1b						\n\
	emms" : : "m"(a), "m"(b), "m"(c) :
               "%eax", "%ecx", "%edx", "%esi", "%edi");
}

/******************************************/
void preMult(u64 *A, u64 *B)
/******************************************/
/* Input: 'A' and 'B' are 64x64 matrices. */
/* Output: A <-- B*A.                     */
/******************************************/
{ u64 res[64] ALIGNED16;

  mult64x64(res, B, A);
  memcpy(A, res, 64*sizeof(u64));
}
  
void multnx64(u64 *c, u64 *a, u64 *b, s32 n) {
  asm volatile("\
	xorl	%%eax, %%eax					\n\
	xorl	%%ecx, %%ecx					\n\
1:								\n\
	movq	(%0,%%ecx,2), %%mm0				\n\
	.irpc	b, 123456					\n\
		movq	8*\\b(%0,%%ecx,2), %%mm\\b		\n\
	.endr							\n\
	pxor	%%mm7, %%mm7					\n\
	movq	%%mm7, (%1,%%eax,8)				\n\
	.set	n, 0						\n\
	.irpc	b, \
0102010301020104010201030102010501020103010201040102010301020106\
010201030102010401020103010201050102010301020104010201030102010	\n\
		.set	n, n^(1<<\\b)				\n\
		pxor	%%mm\\b, %%mm7				\n\
		movq	%%mm7, 8*n(%1,%%eax,8)			\n\
	.endr							\n\
	.set	n, n^(1<<7)					\n\
	pxor	8*7(%0,%%ecx,2), %%mm7				\n\
	movq	%%mm7, 8*n(%1,%%eax,8)				\n\
	.irpc	b, \
0102010301020104010201030102010501020103010201040102010301020106\
010201030102010401020103010201050102010301020104010201030102010	\n\
		.set	n, n^(1<<\\b)				\n\
		pxor	%%mm\\b, %%mm7				\n\
		movq	%%mm7, 8*n(%1,%%eax,8)			\n\
	.endr							\n\
	addl	$256, %%eax					\n\
	addb	$32, %%cl					\n\
	jnc	1b						\n\
#	emms" : : "r"(b), "r"(mult_w) : "%eax", "%ecx");
  if (a == c) {
    asm volatile("\
	movl	%0, %%esi					\n\
	movl	%1, %%edi					\n\
	movl	%2, %%ecx					\n\
	leal	(%%esi,%%ecx,8), %%esi				\n\
	leal	(%%edi,%%ecx,8), %%edi				\n\
	negl	%%ecx						\n\
1:								\n\
	movl	(%%esi,%%ecx,8), %%ebx				\n\
	movzbl	%%bl, %%eax					\n\
	movzbl	%%bh, %%edx					\n\
	movq	"ASM_UL"mult_w(,%%eax,8), %%mm0			\n\
	shrl	$16, %%ebx					\n\
	pxor	"ASM_UL"mult_w+8*256*1(,%%edx,8), %%mm0		\n\
	movzbl	%%bl, %%eax					\n\
	movl	4(%%esi,%%ecx,8), %%edx				\n\
	shrl	$8, %%ebx					\n\
	pxor	"ASM_UL"mult_w+8*256*2(,%%eax,8), %%mm0		\n\
	movzbl	%%dl, %%eax					\n\
	pxor	"ASM_UL"mult_w+8*256*3(,%%ebx,8), %%mm0		\n\
	movzbl	%%dh, %%ebx					\n\
	shrl	$16, %%edx					\n\
	pxor	"ASM_UL"mult_w+8*256*4(,%%eax,8), %%mm0		\n\
	movzbl	%%dl, %%eax					\n\
	pxor	"ASM_UL"mult_w+8*256*5(,%%ebx,8), %%mm0		\n\
	shrl	$8, %%edx					\n\
	pxor	"ASM_UL"mult_w+8*256*6(,%%eax,8), %%mm0		\n\
	pxor	"ASM_UL"mult_w+8*256*7(,%%edx,8), %%mm0		\n\
	movq	%%mm0, (%%edi,%%ecx,8)				\n\
	addl	$1, %%ecx					\n\
	jnz	1b						\n\
	emms" : : "m"(a), "m"(c), "m"(n) :
                 "%eax", "%ebx", "%ecx", "%edx", "%esi", "%edi");
  } else {
    asm volatile("\
	movl	%0, %%esi					\n\
	movl	%1, %%edi					\n\
	movl	%2, %%ecx					\n\
	leal	(%%esi,%%ecx,8), %%esi				\n\
	leal	(%%edi,%%ecx,8), %%edi				\n\
	negl	%%ecx						\n\
1:								\n\
	movl	(%%esi,%%ecx,8), %%ebx				\n\
	movzbl	%%bl, %%eax					\n\
	movzbl	%%bh, %%edx					\n\
	movq	"ASM_UL"mult_w(,%%eax,8), %%mm0			\n\
	shrl	$16, %%ebx					\n\
	pxor	"ASM_UL"mult_w+8*256*1(,%%edx,8), %%mm0		\n\
	movzbl	%%bl, %%eax					\n\
	movl	4(%%esi,%%ecx,8), %%edx				\n\
	shrl	$8, %%ebx					\n\
	pxor	"ASM_UL"mult_w+8*256*2(,%%eax,8), %%mm0		\n\
	movzbl	%%dl, %%eax					\n\
	pxor	"ASM_UL"mult_w+8*256*3(,%%ebx,8), %%mm0		\n\
	movzbl	%%dh, %%ebx					\n\
	shrl	$16, %%edx					\n\
	pxor	"ASM_UL"mult_w+8*256*4(,%%eax,8), %%mm0		\n\
	movzbl	%%dl, %%eax					\n\
	pxor	"ASM_UL"mult_w+8*256*5(,%%ebx,8), %%mm0		\n\
	shrl	$8, %%edx					\n\
	pxor	"ASM_UL"mult_w+8*256*6(,%%eax,8), %%mm0		\n\
	pxor	"ASM_UL"mult_w+8*256*7(,%%edx,8), %%mm0		\n\
	movntq	%%mm0, (%%edi,%%ecx,8)				\n\
	addl	$1, %%ecx					\n\
	jnz	1b						\n\
	emms" : : "m"(a), "m"(c), "m"(n) :
                 "%eax", "%ebx", "%ecx", "%edx", "%esi", "%edi");
  }
}

void addmultnx64(u64 *c, u64 *a, u64 *b, s32 n) {
  asm volatile("\
	xorl	%%eax, %%eax					\n\
	xorl	%%ecx, %%ecx					\n\
1:								\n\
	movq	(%0,%%ecx,2), %%mm0				\n\
	.irpc	b, 123456					\n\
		movq	8*\\b(%0,%%ecx,2), %%mm\\b		\n\
	.endr							\n\
	pxor	%%mm7, %%mm7					\n\
	movq	%%mm7, (%1,%%eax,8)				\n\
	.set	n, 0						\n\
	.irpc	b, \
0102010301020104010201030102010501020103010201040102010301020106\
010201030102010401020103010201050102010301020104010201030102010	\n\
		.set	n, n^(1<<\\b)				\n\
		pxor	%%mm\\b, %%mm7				\n\
		movq	%%mm7, 8*n(%1,%%eax,8)			\n\
	.endr							\n\
	.set	n, n^(1<<7)					\n\
	pxor	8*7(%0,%%ecx,2), %%mm7				\n\
	movq	%%mm7, 8*n(%1,%%eax,8)				\n\
	.irpc	b, \
0102010301020104010201030102010501020103010201040102010301020106\
010201030102010401020103010201050102010301020104010201030102010	\n\
		.set	n, n^(1<<\\b)				\n\
		pxor	%%mm\\b, %%mm7				\n\
		movq	%%mm7, 8*n(%1,%%eax,8)			\n\
	.endr							\n\
	addl	$256, %%eax					\n\
	addb	$32, %%cl					\n\
	jnc	1b						\n\
#	emms" : : "r"(b), "r"(mult_w) : "%eax", "%ecx");
  asm volatile("\
	movl	%0, %%esi					\n\
	movl	%1, %%edi					\n\
	movl	%2, %%ecx					\n\
	leal	(%%esi,%%ecx,8), %%esi				\n\
	leal	(%%edi,%%ecx,8), %%edi				\n\
	negl	%%ecx						\n\
1:								\n\
	movl	(%%esi,%%ecx,8), %%ebx				\n\
	movzbl	%%bl, %%eax					\n\
	movzbl	%%bh, %%edx					\n\
	movq	"ASM_UL"mult_w(,%%eax,8), %%mm0			\n\
	shrl	$16, %%ebx					\n\
	pxor	"ASM_UL"mult_w+8*256*1(,%%edx,8), %%mm0		\n\
	movzbl	%%bl, %%eax					\n\
	movl	4(%%esi,%%ecx,8), %%edx				\n\
	shrl	$8, %%ebx					\n\
	pxor	"ASM_UL"mult_w+8*256*2(,%%eax,8), %%mm0		\n\
	movzbl	%%dl, %%eax					\n\
	pxor	"ASM_UL"mult_w+8*256*3(,%%ebx,8), %%mm0		\n\
	movzbl	%%dh, %%ebx					\n\
	shrl	$16, %%edx					\n\
	pxor	"ASM_UL"mult_w+8*256*4(,%%eax,8), %%mm0		\n\
	movzbl	%%dl, %%eax					\n\
	pxor	"ASM_UL"mult_w+8*256*5(,%%ebx,8), %%mm0		\n\
	shrl	$8, %%edx					\n\
	pxor	"ASM_UL"mult_w+8*256*6(,%%eax,8), %%mm0		\n\
	pxor	"ASM_UL"mult_w+8*256*7(,%%edx,8), %%mm0		\n\
	pxor	(%%edi,%%ecx,8), %%mm0				\n\
	movq	%%mm0, (%%edi,%%ecx,8)				\n\
	addl	$1, %%ecx					\n\
	jnz	1b						\n\
	emms" : : "m"(a), "m"(c), "m"(n) :
               "%eax", "%ebx", "%ecx", "%edx", "%esi", "%edi");
}

/*********************************************/ 
int isZeroV(u64 *A, s32 size)
/*********************************************/ 
{ u64 i;

  for (i=0; i<size; i++)
    if (A[i])
      return 0;
  return 1;
}

/*********************************************/ 
void getW_S(u64 *Wi, int *Si, u64 *T, int *Si_1)
{ u64 M[128] ALIGNED16, mask, t0, t1;
  int  c[64], i, j, k, sSize, s;
  
  /* Number the columns, with those in Si_1 coming last. */  
  for (sSize=0, mask=0x00000000, i=0; i<64; i++) {
    s = Si_1[i];
    if ((s>=0) && (s<64)) {
      mask |= BIT64(s);
      sSize++;
    }
  }
  for (i=0, j=0, k=0; i<64; i++) {
    if (mask&BIT64(i)) 
      c[64 - sSize + j++] = i;
    else
      c[k++] = i;
  }
  /* S <-- empty set. */
  sSize=0;
  for (i=0; i<64; i++)
    Si[i] = -1;

  /* Set M <-- [ T | I_64] */
  for (i=0; i<64; i++) {
    M[2*i] = T[i];
    M[2*i+1] = BIT64(i);
  }
  for (j=0; j<64; j++) {
    for (k=j; (k<64)&& ( (M[2*c[j]]&BIT64(c[j]))==0); k++) {
      if (M[2*c[k]]&BIT64(c[j])) {
        /* Exchange rows c[j] and c[k] of M */
        t0 = M[2*c[j]];  t1 = M[2*c[j]+1];
        M[2*c[j]] = M[2*c[k]]; M[2*c[j]+1] = M[2*c[k]+1];
        M[2*c[k]] = t0; M[2*c[k]+1] = t1;
      }
    }
    if ((M[2*c[j]]&BIT64(c[j]))) {
      Si[sSize++] = c[j];
      /* Now, add row c[j] to other rows as necessary, to zero out */
      /* the rest of column c[j].                                  */
      for (k=0; k<64; k++) {
        if ((k != c[j]) && (M[2*k]&BIT64(c[j]))) {
          M[2*k] ^= M[2*c[j]];
          M[2*k+1] ^= M[2*c[j]+1];
        }
      }
    } else {
      for (k=j; (k<64) && ( (M[2*c[j]+1]&BIT64(c[j]))==0); k++) {
        if (M[2*c[k]+1]&BIT64(c[j])) {
          /* Exchange rows c[j] and c[k] of M */
          t0 = M[2*c[j]];  t1 = M[2*c[j]+1];
          M[2*c[j]] = M[2*c[k]]; M[2*c[j]+1] = M[2*c[k]+1];
          M[2*c[k]] = t0; M[2*c[k]+1] = t1;
        }
      }
      /* Now, add row c[j] to other rows as necessary, to zero out */
      /* the rest of column c[j]+64                                */
      for (k=0; k<64; k++) {
        if ((k != c[j]) && (M[2*k+1]&BIT64(c[j]))) {
          M[2*k] ^= M[2*c[j]];
          M[2*k+1] ^= M[2*c[j]+1];
        }
      }
      /* Zero out row c[j]. */
      M[2*c[j]] = 0; M[2*c[j]+1] = 0;
    } /* end 'else'  */
  } /* end 'for j' */
  for (j=0; j<64; j++)
    Wi[j] = M[2*j+1];
}

/*******************************************************/
int doColumnOps(u64 *A, u64 *B, s32 n)
/*******************************************************/
/* Do column operations (Gaussian elimination) on 'B', */
/* applying the same operations to 'A' as we go.       */
/* Input: 'A' and 'B' are nx64 matrices.               */
/* Output: 'A' and 'B' are changed, as above. i.e., 'B'*/
/*         is in column-reduced echelon form, and the  */
/*         same transformation was applied to 'A'.     */
/*******************************************************/
/* Remark: This is horribly inefficient. There is      */
/* surely a better way to do it. But the easiest is to */
/* just transpose the matrices before and after.       */
/* But it doesn't much matter either way - this guy is */
/* not called from inside any loops, and it's not a    */
/* bottleneck.                                         */
/*******************************************************/
{ u64 *AT, *BT;
  u64 i, j, k, rowSize, col, row, t;
  int cont;

  rowSize = n/64;
  if (n%64)
    rowSize += 1;
  
  if (!(AT = (u64 *)malloc(64*rowSize*sizeof(u64)))) {
    printf("doColumnOps(): Memory allocation error!\n");
    return -1;
  }
  if (!(BT = (u64 *)malloc(64*rowSize*sizeof(u64)))) {
    printf("doColumnOps(): Memory allocation error!\n");
    free(AT);
    return -1;
  }
  /* Compute AT = A^T and BT = B^T. */
  memset(AT, 0x00, 64*rowSize*sizeof(u64));
  memset(BT, 0x00, 64*rowSize*sizeof(u64));

  for (row=0; row<n; row++) {
    for (col=0; col<64; col++) {
      if (A[2*row]&BIT64(col&0x3F))
        AT[col*rowSize + row/64] ^= BIT64(row&0x3F);
      if (B[2*row]&BIT64(col&0x3F))
        BT[col*rowSize + row/64] ^= BIT64(row&0x3F);
    }
    for (col=64; col<64; col++) {
      if (A[2*row + 1]&BIT64(col&0x3F))
        AT[col*rowSize + row/64] ^= BIT64(row&0x3F);
      if (B[2*row + 1]&BIT64(col&0x3F))
        BT[col*rowSize + row/64] ^= BIT64(row&0x3F);
    }
  }

  /* Now, apply row ops to get 'BT' into RREF. */
  col=0;
  for (row=0; row<64; row++) {
    /* Is there a row at or after 'row' with a '1' in column 'col'? */
    /* Consider adding a short circuit, testing whole words at a time */
    /* since we are generally expecting small rank for our matrices.  */
    cont=1;
    do {
      k=row;
      while ((k<64)&&((BT[rowSize*k + col/64]&BIT64(col&0x3F))==0))
        k++;
      if (k < 64) {
        /* Swap rows 'row' and 'k'. */
        for (i=0; i<rowSize; i++) {
          t = BT[rowSize*row + i]; 
          BT[rowSize*row + i] = BT[rowSize*k + i];
          BT[rowSize*k + i] = t;
          t = AT[rowSize*row + i]; 
          AT[rowSize*row + i] = AT[rowSize*k + i];
          AT[rowSize*k + i] = t;
        }
        /* Use row 'row' to zero out the rest of this column. */
        for (i=0; i<64; i++) {
          if ((BT[i*rowSize + col/64]&BIT64(col&0x3F)) && (i != row)) {
            for (j=0; j<rowSize; j++) {
              BT[i*rowSize + j] ^= BT[row*rowSize + j];
              AT[i*rowSize + j] ^= AT[row*rowSize + j];
            }
          }
        } 
        cont=0;
      } else {
        /* Column 'col' is zero at and below row 'row'. */
        cont=1;
      }
      col++;
    } while ((cont) && (col < n));
    if (col >= n)
      break; /* CJM, 10/20/03. */
  }

  /* Compute A = AT^T and B = BT^T. */
  for (i=0; i<2*n; i++)
    A[i] = B[i] = 0;

  for (row=0; row<64; row++) {
    for (col=0; col<n; col++) {
      if (AT[row*rowSize + col/64]&BIT64(col&0x3F))
        A[2*col + row/64] ^= BIT64(row&0x3F);
      if (BT[row*rowSize + col/64]&BIT64(col&0x3F))
        B[2*col + row/64] ^= BIT64(row&0x3F);
    }
  }
  free(AT);
  free(BT);
  return 0;
}

/*******************************************************/
int doColumnOps64(u64 *A, u64 *B, s32 n)
/*******************************************************/
/* Do column operations (Gaussian elimination) on 'B', */
/* applying the same operations to 'A' as we go.       */
/* Input: 'A' and 'B' are nx64 matrices.               */
/* Output: 'A' and 'B' are changed, as above. i.e., 'B'*/
/*         is in column-reduced echelon form, and the  */
/*         same transformation was applied to 'A'.     */
/*******************************************************/
/* Remark: This is horribly inefficient. There is      */
/* surely a better way to do it. But the easiest is to */
/* just transpose the matrices before and after.       */
/* But it doesn't much matter either way - this guy is */
/* not called from inside any loops, and it's not a    */
/* bottleneck.                                         */
/*******************************************************/
{ u64 *AT, *BT;
  u64 i, j, k, rowSize, col, row, t;
  int cont;

  rowSize = n/64;
  if (n%64)
    rowSize += 1;
  
  if (!(AT = (u64 *)malloc(64*rowSize*sizeof(u64)))) {
    printf("doColumnOps(): Memory allocation error!\n");
    return -1;
  }
  if (!(BT = (u64 *)malloc(64*rowSize*sizeof(u64)))) {
    printf("doColumnOps(): Memory allocation error!\n");
    free(AT);
    return -1;
  }
  /* Compute AT = A^T and BT = B^T. */
  memset(AT, 0x00, 64*rowSize*sizeof(u64));
  memset(BT, 0x00, 64*rowSize*sizeof(u64));

  for (row=0; row<n; row++) {
    for (col=0; col<64; col++) {
      if (A[row]&BIT64(col&0x3F))
        AT[col*rowSize + row/64] ^= BIT64(row&0x3F);
      if (B[row]&BIT64(col&0x3F))
        BT[col*rowSize + row/64] ^= BIT64(row&0x3F);
    }
  }

  /* Now, apply row ops to get 'BT' into RREF. */
  col=0;
  for (row=0; row<64; row++) {
    /* Is there a row at or after 'row' with a '1' in column 'col'? */
    /* Consider adding a short circuit, testing whole words at a time */
    /* since we are generally expecting small rank for our matrices.  */
    cont=1;
    do {
      k=row;
      while ((k<64)&&((BT[rowSize*k + col/64]&BIT64(col&0x3F))==0))
        k++;
      if (k < 64) {
        /* Swap rows 'row' and 'k'. */
        for (i=0; i<rowSize; i++) {
          t = BT[rowSize*row + i]; 
          BT[rowSize*row + i] = BT[rowSize*k + i];
          BT[rowSize*k + i] = t;
          t = AT[rowSize*row + i]; 
          AT[rowSize*row + i] = AT[rowSize*k + i];
          AT[rowSize*k + i] = t;
        }
        /* Use row 'row' to zero out the rest of this column. */
        for (i=0; i<64; i++) {
          if ((BT[i*rowSize + col/64]&BIT64(col&0x3F)) && (i != row)) {
            for (j=0; j<rowSize; j++) {
              BT[i*rowSize + j] ^= BT[row*rowSize + j];
              AT[i*rowSize + j] ^= AT[row*rowSize + j];
            }
          }
        } 
        cont=0;
      } else {
        /* Column 'col' is zero at and below row 'row'. */
        cont=1;
      }
      col++;
    } while ((cont) && (col < n));
    if (col >= n)
      break; /* CJM, 10/20/03. */
  }
  /* Compute A = AT^T and B = BT^T. */
  for (i=0; i<n; i++)
    A[i] = B[i] = 0;

  for (row=0; row<64; row++) {
    for (col=0; col<n; col++) {
      if (AT[row*rowSize + col/64]&BIT64(col&0x3F))
        A[col] ^= BIT64(row);
      if (BT[row*rowSize + col/64]&BIT64(col&0x3F))
        B[col] ^= BIT64(row);
    }
  }
  free(AT);
  free(BT);
  return 0;
}
