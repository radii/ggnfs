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


	.section	.rodata
	.error_string:
	.string "Bad args to asm_modinv32\n"
	.section	.data
	.iy:
	.short	
	.text
	.align	4
	.globl asm_modinv32
asm_modinv32:
	pushl %ebx
	pushl %esi
	pushl %edi
	# Get args from stack. Test their sanity.
	movl 16(%esp),%edi
	movl modulo32,%esi
	testl %edi,%edi
	jz badargs
	cmpl %edi,%esi
	jbe badargs
	# Set xc and yc to their initial values.
	xorl %ebx,%ebx
	xorl %ecx,%ecx
	incl %ecx
	cmpl $1,%edi
	jbe have_inverse2
divide:

	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	subl %edi,%esi
	addl %ecx,%ebx
	cmpl %edi,%esi
	jb xlarger
	movl %esi,%eax
	xorl %edx,%edx
	divl %edi
	movl %edx,%esi
	mull %ecx
	addl %eax,%ebx

xlarger:
	cmpl $1,%esi
	jbe have_inverse1
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	subl %esi,%edi
	addl %ebx,%ecx
	cmpl %esi,%edi
	jb ylarger
	movl %edi,%eax
	xorl %edx,%edx
	divl %esi
	movl %edx,%edi
	mull %ebx
	addl %eax,%ecx
ylarger:
	cmpl $1,%edi
	ja divide
have_inverse2:
	jne badargs
	movl %ecx,%eax
	popl %edi
	popl %esi
	popl %ebx
	ret
have_inverse1:
	jne badargs
	movl modulo32,%eax
	subl %ebx,%eax
	popl %edi
	popl %esi
	popl %ebx
	ret
badargs:
	pushl $.error_string
	call Schlendrian
