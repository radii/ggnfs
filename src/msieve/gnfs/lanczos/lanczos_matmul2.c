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
----------------------------------------------------------------------*/

#include "lanczos.h"

	/* code for handling transposed matrix multiplies 
	   when the matrix is in packed format */

/*-------------------------------------------------------------------*/
static void mul_trans_one_med_block(packed_block_t *curr_block,
			uint64 *curr_row, uint64 *curr_b) {

	uint16 *entries = curr_block->med_entries;

	while (1) {
		uint64 t;
#if (defined(__GNUC__) || defined(__ICL)) && defined(__x86_64__)
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

		t = curr_row[row];

		/* Unlike the sparse blocks, medium-dense blocks
		   have enough entries that they can be stored in
		   row-major order, with many entries in each row.
		   One iteration of the while loop handles an entire
		   row at a time */

		/* curr_row and curr_b are both cached, so we have to
		   minimize the number of memory accesses and calculate
		   pointers as early as possible */

#if (defined(__GNUC__) || defined(__ICL)) && defined(__i386__) && \
	defined(NDEBUG) && defined(HAS_MMX)

	#define _txor(k)				\
		"movl 2*(2+2+(" #k "))(%2,%3,2), %%ecx	\n\t"	\
		"movzwl %%ax, %%edx      	\n\t"	\
		"shrl $16, %%eax         	\n\t"	\
		"movq (%1,%%edx,8), %%mm0	\n\t"	\
		"movq (%1,%%eax,8), %%mm1	\n\t"	\
		"pxor %5, %%mm0			\n\t"	\
		"pxor %5, %%mm1			\n\t"	\
		"movq %%mm0, (%1,%%edx,8)	\n\t"	\
		"movq %%mm1, (%1,%%eax,8)	\n\t"	\
		"movl 2*(2+4+(" #k "))(%2,%3,2), %%eax	\n\t"	\
		"movzwl %%cx, %%edx      	\n\t"	\
		"shrl $16, %%ecx         	\n\t"	\
		"movq (%1,%%edx,8), %%mm0	\n\t"	\
		"movq (%1,%%ecx,8), %%mm1	\n\t"	\
		"pxor %5, %%mm0			\n\t"	\
		"pxor %5, %%mm1			\n\t"	\
		"movq %%mm0, (%1,%%edx,8)	\n\t"	\
		"movq %%mm1, (%1,%%ecx,8)	\n\t"

	asm volatile(
		"movl 2*(2+0)(%2,%3,2), %%eax	\n\t"
		"cmpl $0, %4			\n\t"
		"je 1f				\n\t"
		".p2align 4,,7			\n\t"
		"0:				\n\t"

		_txor(0) _txor(4) _txor(8) _txor(12)

		"addl $16, %3			\n\t"
		"cmpl %4, %3			\n\t"
		"jne 0b				\n\t"
		"1:				\n\t"
		:"=r"(i)
		:"r"(curr_b), "r"(entries), "0"(i), 
			"g"(count & (uint32)(~15)), "y"(t)
		:"%eax", "%ecx", "%edx", "%mm0", "%mm1", "memory");

	#undef _txor

#elif (defined(__GNUC__) || defined(__ICL)) && defined(__x86_64__)

	#define _txor(k)				\
		"movzwq %%r8w, %%r9          	\n\t"	\
		"xorq %5, (%1,%%r9,8)         	\n\t"	\
		"shrq $16, %%r8              	\n\t"	\
		"xorq %5, (%1,%%r8,8)         	\n\t"	\
		"movl 2*(2+8+" #k ")(%2,%3,2), %%r8d	\n\t"	\
		"movzwq %%r10w, %%r11          	\n\t"	\
		"xorq %5, (%1,%%r11,8)         	\n\t"	\
		"shrq $16, %%r10              	\n\t"	\
		"xorq %5, (%1,%%r10,8)         	\n\t"	\
		"movl 2*(2+10+" #k ")(%2,%3,2), %%r10d	\n\t"	\
		"movzwq %%r12w, %%r13          	\n\t"	\
		"xorq %5, (%1,%%r13,8)         	\n\t"	\
		"shrq $16, %%r12              	\n\t"	\
		"xorq %5, (%1,%%r12,8)         	\n\t"	\
		"movl 2*(2+12+" #k ")(%2,%3,2), %%r12d	\n\t"	\
		"movzwq %%r14w, %%r15          	\n\t"	\
		"xorq %5, (%1,%%r15,8)         	\n\t"	\
		"shrq $16, %%r14              	\n\t"	\
		"xorq %5, (%1,%%r14,8)         	\n\t"	\
		"movl 2*(2+14+" #k ")(%2,%3,2), %%r14d	\n\t"

	asm volatile(
		"movl 2*(2+0)(%2,%3,2), %%r8d	\n\t"
		"movl 2*(2+2)(%2,%3,2), %%r10d	\n\t"
		"movl 2*(2+4)(%2,%3,2), %%r12d	\n\t"
		"movl 2*(2+6)(%2,%3,2), %%r14d	\n\t"
		"cmpq $0, %4			\n\t"
		"je 1f				\n\t"
		".p2align 4,,7			\n\t"
		"0:				\n\t"

		_txor(0) _txor(8)

		"addq $16, %3			\n\t"
		"cmpq %4, %3			\n\t"
		"jne 0b				\n\t"
		"1:				\n\t"
		:"=r"(i)
		:"r"(curr_b), "r"(entries), "0"(i), 
			"r"(count & (uint64)(~15)), "r"(t)
		:"%r8", "%r9", "%r10", "%r11", 
		 "%r12", "%r13", "%r14", "%r15", "memory");

	#undef _txor

#elif defined(_MSC_VER) && !defined(_WIN64) && \
	defined(NDEBUG) && defined(HAS_MMX)

	#define _txor(k)				\
		__asm mov ecx, [2*(2+2+k)+ebx+esi*2]	\
		__asm movzx edx, ax			\
		__asm shr  eax, 16			\
		__asm movq mm0, [edi+edx*8]		\
		__asm movq mm1, [edi+eax*8]		\
		__asm pxor mm0, mm2			\
		__asm pxor mm1, mm2			\
		__asm movq [edi+edx*8], mm0		\
		__asm movq [edi+eax*8], mm1		\
		__asm mov eax, [2*(2+4+k)+ebx+esi*2]	\
		__asm movzx edx, cx			\
		__asm shr ecx, 16			\
		__asm movq mm0, [edi+edx*8]		\
		__asm movq mm1, [edi+ecx*8]		\
		__asm pxor mm0, mm2			\
		__asm pxor mm1, mm2			\
		__asm movq [edi+edx*8], mm0		\
		__asm movq [edi+ecx*8], mm1

	__asm
	{	
		push ebp
		push ebx
		push esi
		push edi
		mov esi, i
		mov edi, curr_b
		mov ebx, entries
		mov ebp, count
		and ebp, ~15
		movq mm2, t
		mov eax,[2*(2+0)+ebx+esi*2]
		cmp ebp, 0
		je L1
		align 16
	L0:	_txor(0) _txor(4) _txor(8) _txor(12)
		add esi, 16
		cmp esi, ebp
		jne L0
	L1: 	mov i, esi
		pop edi
		pop esi
		pop ebx
		pop ebp
	}

	#undef _txor

#else
		for (i = 0; i < (count & (uint32)(~15)); i += 16) {
			curr_b[entries[i+2+ 0]] ^= t;
			curr_b[entries[i+2+ 1]] ^= t;
			curr_b[entries[i+2+ 2]] ^= t;
			curr_b[entries[i+2+ 3]] ^= t;
			curr_b[entries[i+2+ 4]] ^= t;
			curr_b[entries[i+2+ 5]] ^= t;
			curr_b[entries[i+2+ 6]] ^= t;
			curr_b[entries[i+2+ 7]] ^= t;
			curr_b[entries[i+2+ 8]] ^= t;
			curr_b[entries[i+2+ 9]] ^= t;
			curr_b[entries[i+2+10]] ^= t;
			curr_b[entries[i+2+11]] ^= t;
			curr_b[entries[i+2+12]] ^= t;
			curr_b[entries[i+2+13]] ^= t;
			curr_b[entries[i+2+14]] ^= t;
			curr_b[entries[i+2+15]] ^= t;
		}
#endif
		for (; i < count; i++)
			curr_b[entries[i+2]] ^= t;
		entries += count + 2;
	}
}

/*-------------------------------------------------------------------*/
static void mul_trans_one_block(packed_block_t *curr_block,
				uint64 *curr_row, uint64 *curr_b) {

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
	   more xor operations. Also convert two 16-bit reads into
	   a single 32-bit read with unpacking arithmetic */

#if (defined(__GNUC__) || defined(__ICL)) && \
	defined(__i386__) && defined(HAS_MMX)

	#define _txor(x)				\
		"movl 4*" #x "(%1,%4,4), %%eax    \n\t"	\
		"movzwl %%ax, %%ecx               \n\t"	\
		"movq (%3,%%ecx,8), %%mm0         \n\t"	\
		"shrl $16, %%eax                  \n\t"	\
		"pxor (%2,%%eax,8), %%mm0         \n\t"	\
		"movq %%mm0, (%2,%%eax,8)         \n\t"

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
		:"r"(entries), "r"(curr_b), "r"(curr_row), 
		 "0"(i), "g"(num_entries & (uint32)(~15))
		:"%eax", "%ecx", "%mm0", "memory");

#elif defined(_MSC_VER) && !defined(_WIN64) && \
	defined(NDEBUG) && defined(HAS_MMX)

	#define _txor(x)			\
		__asm mov eax,[4*x+edi+esi*4]	\
		__asm movzx ecx, ax		\
		__asm movq mm0, [ebp+ecx*8]	\
		__asm shr eax, 16		\
		__asm pxor mm0, [ebx+eax*8]	\
		__asm movq [ebx+eax*8], mm0

	__asm
	{
		push ebp
		push ebx
		push esi
		push edi
		mov esi, i
		mov edi, entries
		mov ebx, curr_b
		mov ebp, curr_row
		mov edx, num_entries
		and edx, ~15
		cmp edx, 0
		je L1
		align 16
	L0:	_txor( 0) _txor( 1) _txor( 2) _txor( 3)
		_txor( 4) _txor( 5) _txor( 6) _txor( 7)
		_txor( 8) _txor( 9) _txor(10) _txor(11)
		_txor(12) _txor(13) _txor(14) _txor(15)
		add esi, 16
		cmp esi, edx
		jne L0
	L1:	mov i, esi
		pop edi
		pop esi
		pop ebx
		pop ebp
	}

#else
	#define _txor(x) curr_b[entries[i+x].col_off] ^= \
				 curr_row[entries[i+x].row_off]	

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
		curr_b[k] ^= curr_row[j];
	}
}

/*-------------------------------------------------------------------*/
void mul_trans_packed_core(thread_data_t *t) {

	uint64 *x = t->x;
	uint64 *b = t->b;
	uint32 i;

	/* you would think that doing the matrix multiply
	   in column-major order would be faster, since this would
	   also minimize dirty cache writes. Except that it's slower;
	   it appears that row-major matrix access minimizes some
	   kind of thrashing somewhere */

	for (i = 0; i < t->num_blocks; i++) {
		packed_block_t *curr_block = t->blocks + i;
		if (curr_block->med_entries)
			mul_trans_one_med_block(curr_block, 
					x + curr_block->start_row,
					b + curr_block->start_col);
		else
			mul_trans_one_block(curr_block, 
					x + curr_block->start_row,
					b + curr_block->start_col);
	}

	/* multiply the densest few rows by x (in batches of 64 rows) */

	for (i = 0; i < (t->num_dense_rows + 63) / 64; i++) {
		mul_Nx64_64x64_acc(t->dense_blocks[i], x + 64 * i, 
				   b + t->blocks[0].start_col, t->ncols);
	}
}
