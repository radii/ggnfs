# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# The function we want to write.
# ulong asm_modinv32(x)
# Modular inverse of x modulo modulo32
# and y satisfies 0<y<x.

x		.ASSIGNC	"%edi"
y		.ASSIGNC	"%esi"
xc		.ASSIGNC	"%ecx"
yc		.ASSIGNC	"%ebx"
# Number of trial subtractions before doing a division
nts		.ASSIGNA	26

# %eax and %edx are used in the division.

	.MACRO s1 k kc l lc
	subl \k,\l
	addl \kc,\lc
	.ENDM

	.MACRO ts1 k kc l lc
	cmpl \k,\l
	setbe %al
	decl %eax
	movl %eax,%edx
	andl \k,%edx
	subl %edx,\l
	andl \kc,%eax
	addl %eax,\lc
	.ENDM

.section	.rodata
.error_string:
	.string "Bad args to asm_modinv32\n"
.section	.data
.iy:
	.dw
.text
	.align 4
.globl asm_modinv32
	.type	 asm_modinv32,@function	
asm_modinv32:
	pushl %ebx
	pushl %esi
	pushl %edi
# Get args from stack. Test their sanity.
	movl 16(%esp),\&x
	movl modulo32,\&y
	testl \&x,\&x
	jz badargs
	cmpl \&x,\&y
	jbe badargs
# Set xc and yc to their initial values.
	xorl \&yc,\&yc
	xorl \&xc,\&xc
	incl \&xc
	cmpl $1,\&x
	jbe have_inverse2
divide:

ntsc	.ASSIGNA 0
.AWHILE	\&ntsc LE \&nts
	s1 \&x \&xc \&y \&yc
	cmpl \&x,\&y
	jb xlarger
ntsc	.ASSIGNA \&ntsc + 1
.AENDW
	movl \&y,%eax
	xorl %edx,%edx
	divl \&x
	movl %edx,\&y
	mull \&xc
	addl %eax,\&yc

xlarger:
	cmpl $1,\&y
	jbe have_inverse1
ntsc	.ASSIGNA 0
.AWHILE	\&ntsc LE \&nts
	s1 \&y \&yc \&x \&xc
	cmpl \&y,\&x
	jb ylarger
ntsc	.ASSIGNA \&ntsc + 1
.AENDW
        movl \&x,%eax
	xorl %edx,%edx
	divl \&y
	movl %edx,\&x
	mull \&yc
	addl %eax,\&xc
ylarger:
	cmpl $1,\&x
	ja divide
have_inverse2:
	jne badargs
	movl \&xc,%eax
	popl %edi
	popl %esi
	popl %ebx
	ret
have_inverse1:
	jne badargs
	movl modulo32,%eax
	subl \&yc,%eax
	popl %edi
	popl %esi
	popl %ebx
	ret
badargs:
	pushl $.error_string
	call Schlendrian
	.END
