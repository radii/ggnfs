# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# The function we want to write.
# ulong asm_getbc(x)
# Modular inverse of x modulo modulo32
# and y satisfies 0<y<x.

x		.ASSIGNC	"%edi"
y		.ASSIGNC	"%esi"
xc		.ASSIGNC	"%ecx"
yc		.ASSIGNC	"%ebx"
A		.ASSIGNC	"%ebp"

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

.text
	.align 4
.globl asm_getbc
	.type	 asm_getbc,@function	
asm_getbc:
	pushl %ebx
	pushl %esi
	pushl %edi
	pushl %ebp
# Get args from stack. Test their sanity.
# Set xc and yc to their initial values.
	movl 20(%esp),\&x
	xorl \&xc,\&xc
	movl 28(%esp),\&A
	xorl \&yc,\&yc
	movl 24(%esp),\&y
	incl \&xc
	cmpl \&x,\&A
	ja   have_bs
divide:
.AREPEAT \&nts
	s1 \&x \&xc \&y \&yc
	cmpl \&x,\&y
	jb test_y
.AENDR
	movl \&y,%eax
	xorl %edx,%edx
	divl \&x
	movl %edx,\&y
	mull \&xc
	addl %eax,\&yc
test_y:
	cmpl \&y,\&A
	ja have_ct
ntsc	.ASSIGNA 0
.AREPEAT \&nts
	s1 \&y \&yc \&x \&xc
	cmpl \&y,\&x
	jb test_x
.AENDR
        movl \&x,%eax
	xorl %edx,%edx
	divl \&y
	movl %edx,\&x
	mull \&yc
	addl %eax,\&xc
test_x:
	cmpl \&x,\&A
	jbe divide
have_bs:
.AREPEAT \&nts
	s1 \&x \&xc \&y \&yc
	cmpl \&y,\&A
	ja have_bsct
.AENDR
	movl \&y,%eax
	xorl %edx,%edx
	subl \&A,%eax
	divl \&x
	incl %eax
	movl %eax,\&A
	mull \&x
	subl %eax,\&y
	movl \&A,%eax
	mull \&xc
	addl %eax,\&yc
	jmp have_bsct
have_ct:
.AREPEAT \&nts
	s1 \&y \&yc \&x \&xc
	cmpl \&x,\&A
	ja have_bsct
.AENDR
	movl \&x,%eax
	xorl %edx,%edx
	subl \&A,%eax
	divl \&y
	incl %eax
	movl %eax,\&A
	mull \&y
	subl %eax,\&x
	movl \&A,%eax
	mull \&yc
	addl %eax,\&xc
have_bsct:
	movl 32(%esp),%eax
	movl 36(%esp),%edx
	movl 40(%esp),%ebp
	movl \&x,(%eax)
	movl 44(%esp),%eax
	movl \&xc,(%edx)
	movl \&y,(%ebp)
	movl \&yc,(%eax)
	popl %ebp
	popl %edi
	popl %esi
	popl %ebx
	ret
