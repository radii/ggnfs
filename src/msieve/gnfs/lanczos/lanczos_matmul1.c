/*----------------------------------------------------------------------
Copyright 2007, Jason Papadopoulos

This file is part of GGNFS.

GGNFS is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GGNFS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GGNFS; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

$Id: lanczos_matmul1.c,v 1.1 2007-12-16 03:54:47 jasonp_sf Exp $
----------------------------------------------------------------------*/

#include "lanczos.h"

	/* code for handling matrix multiplies when the
	   matrix is in packed format */

/*-------------------------------------------------------------------*/
static void mul_one_med_block(packed_block_t *curr_block,
			uint64 *curr_col, uint64 *curr_b) {

	uint16 *entries = curr_block->med_entries;

	while (1) {
		uint64 accum;
#if (defined(__GNUC__) || defined(__ICL) ) && defined(__x86_64__)
		uint64 i = 0;
		uint64 row = entries[0];
		uint64 count = entries[1];
#else
		uint32 i = 0;
		uint32 row = entries[0];
		uint32 count = entries[1];
#endif

		if (count == 0)
			break;

		/* Unlike the sparse blocks, medium-dense blocks
		   have enough entries that they can be stored in
		   row-major order, with many entries in each row.
		   One iteration of the while loop handles an entire
		   row at a time */

		/* curr_col and curr_b are both cached, so we have to
		   minimize the number of memory accesses and calculate
		   pointers as early as possible */

#if (defined(__GNUC__) || defined(__ICL)) && defined(__i386__) && \
	defined(NDEBUG) && defined(HAS_MMX)

	#define _txor(k)				\
		"movzwl %%ax, %%edx		\n\t"	\
		"pxor (%2,%%edx,8), %0		\n\t"	\
		"shrl $16, %%eax		\n\t"	\
		"pxor (%2,%%eax,8), %%mm0	\n\t"	\
		"movl 2*(2+4+(" #k "))(%3,%4,2), %%eax \n\t"	\
		"movzwl %%cx, %%edx		\n\t"	\
		"pxor (%2,%%edx,8), %0		\n\t"	\
		"shrl $16, %%ecx		\n\t"	\
		"pxor (%2,%%ecx,8), %%mm0	\n\t"	\
		"movl 2*(2+6+(" #k "))(%3,%4,2), %%ecx \n\t"

	asm volatile(
		"movl 2*(2+0)(%3,%4,2), %%eax	\n\t"
		"movl 2*(2+2)(%3,%4,2), %%ecx	\n\t"
		"pxor %0, %0			\n\t"
		"pxor %%mm0, %%mm0		\n\t"
		"cmpl $0, %5			\n\t"
		"je 1f				\n\t"
		".p2align 4,,7			\n\t"
		"0:				\n\t"

		_txor(0) _txor(4) _txor(8) _txor(12)

		"addl $16, %4			\n\t"
		"cmpl %5, %4			\n\t"
		"jne 0b				\n\t"
		"pxor %%mm0, %0			\n\t"
		"1:				\n\t"

		:"=y"(accum), "=r"(i)
		:"r"(curr_col), "r"(entries), "1"(i), 
			"g"(count & (uint32)(~15))
		:"%eax", "%ecx", "%edx", "%mm0");

	#undef _txor

#elif (defined(__GNUC__) || defined(__ICL)) && defined(__x86_64__)

	#define _txor(k)				\
		"movzwq %%ax, %%rdx		\n\t"	\
		"xorq (%2,%%rdx,8), %0		\n\t"	\
		"shrq $16, %%rax		\n\t"	\
		"xorq (%2,%%rax,8), %%rsi	\n\t"	\
		"movl 2*(2+4+(" #k "))(%3,%4,2), %%eax \n\t"	\
		"movzwq %%cx, %%rdx		\n\t"	\
		"xorq (%2,%%rdx,8), %0		\n\t"	\
		"shrq $16, %%rcx		\n\t"	\
		"xorq (%2,%%rcx,8), %%rsi	\n\t"	\
		"movl 2*(2+6+(" #k "))(%3,%4,2), %%ecx \n\t"

	asm volatile(
		"movl 2*(2+0)(%3,%4,2), %%eax	\n\t"
		"movl 2*(2+2)(%3,%4,2), %%ecx	\n\t"
		"xorq %0, %0			\n\t"
		"xorq %%rsi, %%rsi		\n\t"
		"cmpq $0, %5			\n\t"
		"je 1f				\n\t"
		".p2align 4,,7			\n\t"
		"0:				\n\t"

		_txor(0) _txor(4) _txor(8) _txor(12)

		"addq $16, %4			\n\t"
		"cmpq %5, %4			\n\t"
		"jne 0b				\n\t"
		"xorq %%rsi, %0			\n\t"
		"1:				\n\t"

		:"=r"(accum), "=r"(i)
		:"r"(curr_col), "r"(entries), "1"(i), 
			"r"(count & (uint64)(~15))
		:"%rax", "%rcx", "%rdx", "%rsi");

	#undef _txor

#elif defined(_MSC_VER) && !defined(_WIN64) && \
	defined(NDEBUG) && defined(HAS_MMX)

	#define _txor(k)				\
	    __asm movzx edx, ax				\
	    __asm pxor	mm1, [esi+edx*8]		\
	    __asm shr 	eax, 16				\
	    __asm pxor 	mm0, [esi+eax*8]		\
	    __asm mov 	eax, [2*(k+2+4)+ebx+edi*2]	\
	    __asm movzx edx, cx				\
	    __asm pxor 	mm1, [esi+edx*8]		\
	    __asm shr	ecx, 16				\
	    __asm pxor 	mm0, [esi+ecx*8]		\
	    __asm mov 	ecx, [2*(k+2+6)+ebx+edi*2]
			
	__asm
	{	
		push ebp
		push ebx
		push esi
		push edi
		mov esi, curr_col
		mov ebx, entries
		mov edi, i
		mov ebp, count
		and ebp, ~15
		mov eax, [2*(2+0)+ebx+edi*2]
		mov ecx, [2*(2+2)+ebx+edi*2]
		pxor mm1, mm1
		pxor mm0, mm0
		cmp ebp, 0
		je L1
		align 16
	L0:	_txor(0) _txor(4) _txor(8) _txor(12)
		add edi, 16
		cmp edi, ebp
		jne L0
		pxor mm1, mm0
	L1:	movq accum, mm1
		mov i, edi
		pop edi
		pop esi
		pop ebx
		pop ebp
	}

	#undef _txor

#else
	accum = 0;
	for (i = 0; i < (count & (uint32)(~15)); i += 16) {
		accum ^= curr_col[entries[i+2+0]] ^
		         curr_col[entries[i+2+1]] ^
		         curr_col[entries[i+2+2]] ^
		         curr_col[entries[i+2+3]] ^
		         curr_col[entries[i+2+4]] ^
		         curr_col[entries[i+2+5]] ^
		         curr_col[entries[i+2+6]] ^
		         curr_col[entries[i+2+7]] ^
		         curr_col[entries[i+2+8]] ^
		         curr_col[entries[i+2+9]] ^
		         curr_col[entries[i+2+10]] ^
		         curr_col[entries[i+2+11]] ^
		         curr_col[entries[i+2+12]] ^
		         curr_col[entries[i+2+13]] ^
		         curr_col[entries[i+2+14]] ^
		         curr_col[entries[i+2+15]];
	}

#endif
		for (; i < count; i++)
			accum ^= curr_col[entries[i+2]];
		curr_b[row] ^= accum;
		entries += count + 2;
	}
}

/*-------------------------------------------------------------------*/
static void mul_one_block(packed_block_t *curr_block,
			uint64 *curr_col, uint64 *curr_b) {

	uint32 i = 0; 
	uint32 j = 0;
	uint32 k;
	uint32 num_entries = curr_block->num_entries;
	entry_idx_t *entries = curr_block->entries;

	/* unroll by 16, i.e. the number of matrix elements
	   in one cache line (usually). For 32-bit x86, we get
	   a huge performance boost by using either SSE or MMX
	   registers; not because they intrinsically are faster,
	   but because using them cuts the number of memory
	   operations in half, allowing the processor to buffer
	   more xor operations. Also replace two 16-bit loads
	   with a single 32-bit load and extra arithmetic to
	   unpack the array indices */

#if (defined(__GNUC__) || defined(__ICL)) && \
	defined(__i386__) && defined(HAS_MMX)

	#define _txor(x)				\
		"movl 4*" #x "(%1,%4,4), %%eax  \n\t"	\
		"movzwl %%ax, %%ecx             \n\t"	\
		"movq (%2,%%ecx,8), %%mm0       \n\t"	\
		"shrl $16, %%eax                \n\t"	\
		"pxor (%3,%%eax,8), %%mm0       \n\t"	\
		"movq %%mm0, (%2,%%ecx,8)       \n\t"

	asm volatile(
		"cmpl $0, %5			\n\t"
		"je 1f				\n\t"
		".p2align 4,,7			\n\t"
		"0:				\n\t"

		_txor( 0) _txor( 1) _txor( 2) _txor( 3)
		_txor( 4) _txor( 5) _txor( 6) _txor( 7)
		_txor( 8) _txor( 9) _txor(10) _txor(11)
		_txor(12) _txor(13) _txor(14) _txor(15)

		"addl $16, %4			\n\t"
		"cmpl %5, %4			\n\t"
		"jne 0b				\n\t"
		"1:				\n\t"

		:"=r"(i)
		:"r"(entries), "r"(curr_b), "r"(curr_col), 
		 "0"(i), "g"(num_entries & (uint32)(~15))
		:"%eax", "%ecx", "%mm0", "memory");

#elif defined(_MSC_VER) && !defined(_WIN64) && \
	defined(NDEBUG) && defined(HAS_MMX)

	#define _txor(x)				\
		__asm mov	eax, [4*x+edi+esi*4]	\
		__asm movzx ecx, ax			\
		__asm movq 	mm0, [ebx+ecx*8]	\
		__asm shr 	eax, 16			\
		__asm pxor 	mm0, [ebp+eax*8]	\
		__asm movq 	[ebx+ecx*8], mm0
	
	__asm
	{
		push ebp
		push ebx
		push esi
		push edi
		mov esi, i
		mov edi, entries
		mov ebx, curr_b
		mov ebp, curr_col
		mov edx, num_entries
		and edx, ~15
		cmp edx, 0
		je L1
		align 16
	L0:	_txor( 0) _txor( 1) _txor( 2) _txor( 3)
		_txor( 4) _txor( 5) _txor( 6) _txor( 7)
		_txor( 8) _txor( 9) _txor(10) _txor(11)
		_txor(12) _txor(13) _txor(14) _txor(15)
		add esi,16
		cmp esi,edx
		jne L0
	L1:	mov i, esi
		pop edi
		pop esi
		pop ebx
		pop ebp
	}

#else
	#define _txor(x) curr_b[entries[i+x].row_off] ^= \
				 curr_col[entries[i+x].col_off]

	for (i = 0; i < (num_entries & (uint32)(~15)); i += 16) {
		#ifdef MANUAL_PREFETCH
		PREFETCH(entries + i + 48);
		#endif

		_txor( 0); _txor( 1); _txor( 2); _txor( 3);
		_txor( 4); _txor( 5); _txor( 6); _txor( 7);
		_txor( 8); _txor( 9); _txor(10); _txor(11);
		_txor(12); _txor(13); _txor(14); _txor(15);
	}
#endif
	#undef _txor

	for (; i < num_entries; i++) {
		j = entries[i].row_off;
		k = entries[i].col_off;
		curr_b[j] ^= curr_col[k];
	}
}

/*-------------------------------------------------------------------*/
void mul_packed_core(thread_data_t *t) {

	uint64 *x = t->x;
	uint64 *b = t->b;
	uint32 i;

	/* proceed block by block. We assume that blocks access
	   the matrix in row-major order; when computing b = A*x
	   this will write to the same block of b repeatedly, and
	   will read from all of x. This reduces the number of
	   dirty cache writebacks, improving performance slightly */

	for (i = 0; i < t->num_blocks; i++) {
		packed_block_t *curr_block = t->blocks + i;
		if (curr_block->med_entries)
			mul_one_med_block(curr_block, 
					x + curr_block->start_col,
					b + curr_block->start_row);
		else
			mul_one_block(curr_block, 
					x + curr_block->start_col,
					b + curr_block->start_row);
	}

	/* multiply the densest few rows by x (in batches of 64 rows) */

	for (i = 0; i < (t->num_dense_rows + 63) / 64; i++) {
		mul_64xN_Nx64(t->dense_blocks[i], 
				x + t->blocks[0].start_col, 
				b + 64 * i, t->ncols);
	}
}
