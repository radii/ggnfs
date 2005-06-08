	#  Copyright (C) 2002 Jens Franke, T. Kleinjung.
	#  This file is part of gnfs4linux, distributed under the terms of the
	#  GNU General Public Licence and WITHOUT ANY WARRANTY.
	#
	#  You should have received a copy of the GNU General Public License along
	#  with this program; see the file COPYING.  If not, write to the Free
	#  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
	#  2111-1307, USA.

	# root-sieve



	.text


	.align	4
	.globl asm_root_sieve
	.globl _asm_root_sieve
asm_root_sieve:
_asm_root_sieve:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp

	movl 24(%esp),%ebx
	movl 28(%esp),%ecx
	leal (%ebx,%ecx,4),%edi        # ptrend
	movl 32(%esp),%ebp
	movl 20(%esp),%edx
	movl (%edx),%esi
	movl 36(%esp),%ecx
	subl 28(%esp),%ecx
	leal (%ebp,%ecx,4),%ecx        # sv[end]
loop1:
	cmpl %edi,%esi
	jnc loop2
	movl (%esi),%eax
	leal 4(%esi),%esi
	addl %eax,(%ebp)
	leal 4(%ebp),%ebp
	jmp loop1

loop2:
	movl 28(%esp),%eax
	andl $1,%eax
	jz loop2_even
loop2_odd:
	cmpl %ecx,%ebp
	jnc loopend
	movl %ebx,%esi
	movl (%esi),%eax
	leal 4(%esi),%esi
	addl %eax,(%ebp)
	leal 4(%ebp),%ebp
innerloop_odd:
	cmpl %edi,%esi
	jnc loop2_odd
	movl (%esi),%eax
	movl 4(%esi),%edx
	addl %eax,(%ebp)
	leal 4(%ebp),%ebp
	leal 8(%esi),%esi
	addl %edx,(%ebp)
	leal 4(%ebp),%ebp
	jmp innerloop_odd

loop2_even:
	cmpl %ecx,%ebp
	jnc loopend
	movl %ebx,%esi
innerloop_even:
	cmpl %edi,%esi
	jnc loop2_even
	movl (%esi),%eax
	leal 4(%esi),%esi
	addl %eax,(%ebp)
	leal 4(%ebp),%ebp
	movl (%esi),%eax
	leal 4(%esi),%esi
	addl %eax,(%ebp)
	leal 4(%ebp),%ebp
	jmp innerloop_even

loopend:
	movl %ebx,%esi
	movl 28(%esp),%ebx
	leal (%ecx,%ebx,4),%ecx
loop3:
	cmpl %ecx,%ebp
	jnc end
	movl (%esi),%eax
	leal 4(%esi),%esi
	addl %eax,(%ebp)
	leal 4(%ebp),%ebp
	jmp loop3

end:
	movl 20(%esp),%edx
	movl %esi,(%edx)
	#	emms

	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret



	.align	4
	.globl asm_root_sieve8
	.globl _asm_root_sieve8
asm_root_sieve8:
_asm_root_sieve8:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp

	movl 24(%esp),%ebx
	movl 28(%esp),%ecx
	leal (%ebx,%ecx,8),%edi        # ptrend
	movl 32(%esp),%ebp
	movl 20(%esp),%edx
	movl (%edx),%esi
	movl 36(%esp),%ecx
	subl 28(%esp),%ecx
	leal (%ebp,%ecx,8),%ecx        # sv[end]
loop1_8:
	cmpl %edi,%esi
	jnc loop2_8
	movq (%esi),%mm0
	movq (%ebp),%mm1
	leal 8(%esi),%esi
	paddusw %mm0,%mm1
	movq %mm1,(%ebp)
	leal 8(%ebp),%ebp
	jmp loop1_8

loop2_8:
	movl 28(%esp),%eax
	andl $1,%eax
	jz loop2_even_8
loop2_odd_8:
	cmpl %ecx,%ebp
	jnc loopend_8
	movl %ebx,%esi
	movq (%esi),%mm0
	movq (%ebp),%mm1
	leal 8(%esi),%esi
	paddw %mm0,%mm1
	movq %mm1,(%ebp)
	leal 8(%ebp),%ebp
innerloop_odd_8:
	cmpl %edi,%esi
	jnc loop2_odd_8
	movq (%esi),%mm0
	movq (%ebp),%mm1
	leal 8(%esi),%esi
	paddw %mm0,%mm1
	movq (%esi),%mm2
	movq %mm1,(%ebp)
	leal 8(%ebp),%ebp
	movq (%ebp),%mm3
	leal 8(%esi),%esi
	paddw %mm2,%mm3
	movq %mm3,(%ebp)
	leal 8(%ebp),%ebp
	jmp innerloop_odd_8

loop2_even_8:
	cmpl %ecx,%ebp
	jnc loopend_8
	movl %ebx,%esi
innerloop_even_8:
	cmpl %edi,%esi
	jnc loop2_even_8
	movq (%esi),%mm0
	movq (%ebp),%mm1
	leal 8(%esi),%esi
	paddw %mm0,%mm1
	movq %mm1,(%ebp)
	leal 8(%ebp),%ebp
	movq (%esi),%mm0
	movq (%ebp),%mm1
	leal 8(%esi),%esi
	paddw %mm0,%mm1
	movq %mm1,(%ebp)
	leal 8(%ebp),%ebp
	jmp innerloop_even_8

loopend_8:
	movl %ebx,%esi
	movl 28(%esp),%ebx
	leal (%ecx,%ebx,8),%ecx
loop3_8:
	cmpl %ecx,%ebp
	jnc end_8
	movq (%esi),%mm0
	movq (%ebp),%mm1
	leal 8(%esi),%esi
	paddw %mm0,%mm1
	movq %mm1,(%ebp)
	leal 8(%ebp),%ebp
	jmp loop3_8

end_8:
	movl 20(%esp),%edx
	movl %esi,(%edx)
	emms

	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret

	.IF 0
lohnt:	sich momentan nicht !



	.align	4
	.globl asm_root_eval8
	.globl _asm_root_eval8
asm_root_eval8:
_asm_root_eval8:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp


	...
	emms

	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret
	.ENDIF
