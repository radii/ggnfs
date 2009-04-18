# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# Written by T. Kleinjung

sv           .ASSIGNC  "16(%esp)"
svend        .ASSIGNC  "20(%esp)"
buffer       .ASSIGNC  "24(%esp)"
nmax         .ASSIGNC  "28(%esp)"


	.align 16
.globl asm_evaluate
# asm_evaluate(sievebegin,sieveend,buffer,nmax):
# scans sievearray between sievebegin and sieveend for entries
# >127 and stores them in buffer (2 Bytes), stores at most nmax
# entries, returns number of stored entries

# edx: counts entries found so far
# esi: points to location of array we are investigating
# edi: points at end of array (=sieveend)
# mm7: 0
# mm0-3:
# ebx, ecx:

# Modified by J. Franke

	.type	 asm_evaluate,@function
asm_evaluate:
	pushl %edi
	pushl %esi
	pushl %ebx


	movl \&sv,%esi
	movl $0x80808080,%eax
	movl \&svend,%edi
	movd %eax,%mm6
	movl \&buffer,%edx
	movq %mm6,%mm1
	psllq $32,%mm6
	movl \&nmax,%ecx
	pxor %mm7,%mm7
	por %mm1,%mm6
	leal (%edx,%ecx,1),%ecx
	movl %ecx,\&nmax	
	jmp entry32

loop32:
	leal 32(%esi),%esi
	movq %mm7,-32(%esi)
	movq %mm7,-24(%esi)
	movq %mm7,-16(%esi)
	movq %mm7,-8(%esi)

entry32:
	cmpl %edi,%esi
	jz end

	movq (%esi),%mm0
	movq 8(%esi),%mm1
	movq 16(%esi),%mm2
	movq 24(%esi),%mm3

	por %mm0,%mm1
	por %mm2,%mm3
	por %mm1,%mm3
	pand %mm6,%mm3
	pmovmskb %mm3,%ecx
	testl %ecx,%ecx
	jz loop32

	movq 8(%esi),%mm1
	pand %mm6,%mm0
	movq 24(%esi),%mm3
	pand %mm6,%mm2
	pmovmskb %mm0,%eax
	pand %mm6,%mm1
	pmovmskb %mm2,%ecx
	pand %mm6,%mm3
	sall $16,%ecx
	pmovmskb %mm1,%ebx
	orl %ecx,%eax
	pmovmskb %mm3,%ecx
	sall $8,%ebx
	sall $24,%ecx
	orl %ebx,%eax
	subl \&sv,%esi
	orl %ecx,%eax
	xorl %ebx,%ebx
loop1:	
	bsfl %eax,%ecx
	addl %ecx,%ebx
	addl %ebx,%esi
	shrl %cl,%eax
	movw %si,(%edx)
	leal 2(%edx),%edx 
	subl %ebx,%esi
	incl %ebx
	shrl $1,%eax
	cmpl %edx,\&nmax
	jbe buffer_full
	testl %eax,%eax
	jnz loop1
	addl \&sv,%esi
	jmp loop32

buffer_full:
	addl \&sv,%esi
loop032:
	cmpl %edi,%esi
	jz end

	leal 32(%esi),%esi
	movq %mm7,-32(%esi)
	movq %mm7,-24(%esi)
	movq %mm7,-16(%esi)
	movq %mm7,-8(%esi)
	jmp loop032

end:
	movl %edx,%eax
	emms
	subl \&buffer,%eax
	popl %ebx
	popl %esi
	popl %edi
	shrl $1,%eax
	ret

	.align 16
.globl asm_evaluate16

	.type	 asm_evaluate16,@function
asm_evaluate16:
	pushl %edi
	pushl %esi
	pushl %ebx


	movl \&sv,%esi
	movl \&svend,%edi
	xorl %edx,%edx
	pxor %mm7,%mm7

	jmp entry32a

loop32a:
	leal 16(%esi),%esi
	movq %mm7,-16(%esi)
	movq %mm7,-8(%esi)

entry32a:
	cmpl %edi,%esi
	jz enda

	movq (%esi),%mm0
	movq 8(%esi),%mm1

	por %mm0,%mm1

	movd %mm1,%ebx
	psrlq $32,%mm1
	movd %mm1,%ecx
	orl %ebx,%ecx
	and $0x80808080,%ecx
	jz loop32a

	movl $16,%ecx
loop1a:
	movb (%esi),%al
	movb $0,(%esi)
	leal 1(%esi),%esi
	andb $128,%al
	jnz founda
entry1a:
	decl %ecx
	jnz loop1a
	jmp entry32a

founda:
	movl %esi,%eax
	decl %eax
	subl \&sv,%eax
	movl \&buffer,%ebx
	movl %eax,(%ebx,%edx,2)
	incl %edx
	cmpl %edx,\&nmax
	ja entry1a

loop0a:
	movb $0,(%esi)
	leal 1(%esi),%esi
	decl %ecx
	jnz loop0a

loop032a:
	cmpl %esi,%edi
	jz enda

	leal 16(%esi),%esi
	movq %mm7,-16(%esi)
	movq %mm7,-8(%esi)
	jmp loop032a

enda:
	emms
	movl %edx,%eax
	popl %ebx
	popl %esi
	popl %edi
	ret
	.END


