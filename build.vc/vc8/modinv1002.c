/* modinv1002.s
   Written by T. Kleinjung and/or Jens Franke.
   6/13/04: Hacked up for use in GGNFS by Chris Monico.
                                                                                                       
   Copyright (C) 2002 Jens Franke, T.Kleinjung
   This file is part of gnfs4linux, distributed under the terms of the
   GNU General Public Licence and WITHOUT ANY WARRANTY.
                                                                                                       
   You should have received a copy of the GNU General Public License along
   with this program; see the file COPYING.  If not, write to the Free
   Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
   2111-1307, USA.
                                                                                                       
                                                                                                       
	The function we want to write.
	ulong asm_modinv32(x)
	Modular inverse of x modulo modulo32
	and y satisfies 0<y<x.

	Number of trial subtractions before doing a division
	
	%eax and %edx are used in the division.
*/

volatile extern unsigned long modulo32;

unsigned long asm_modinv32(unsigned long x)
{	unsigned long res;

	__asm
	{
		mov		edi,[x]	
		mov		esi,[modulo32]
		test	edi,edi
		jz		badargs
		cmp		esi,edi
		jbe		badargs
	
		xor		ebx,ebx
		xor		ecx,ecx
		inc		ecx
		cmp		edi,1
		jbe		have_inverse2
divide:
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		sub		esi,edi
		add		ebx,ecx
		cmp		esi,edi
		jb		xlarger
		mov		eax,esi
		xor		edx,edx
		div		edi
		mov		esi,edx
		mul		ecx
		add		ebx,eax
xlarger:
		cmp		esi,1
		jbe		have_inverse1
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		sub		edi,esi
		add		ecx,ebx
		cmp		edi,esi
		jb		ylarger
		mov		eax,edi
		xor		edx,edx
		div		esi
		mov		edi,edx
		mul		ebx
		add		ecx,eax
ylarger:
		cmp		edi,1
		ja		divide
have_inverse2:
		jne		badargs
		mov		[res],ecx
		}
		return res;

have_inverse1:
	__asm
	{
		jne		badargs
		mov		eax,[modulo32]
		sub		eax,ebx
		mov		[res],eax
	}
	return res;

badargs:
	printf("\nError in asm_modinv32");
	return 0;
}
