/* asm_mul128.s
   Written by T. Kleinjung
   6/13/04: Hacked up for use in GGNFS by Chris Monico.

   Copyright (C) 2002 Jens Franke, T.Kleinjung
   This file is part of gnfs4linux, distributed under the terms of the
   GNU General Public Licence and WITHOUT ANY WARRANTY.

   You should have received a copy of the GNU General Public License along
   with this program; see the file COPYING.  If not, write to the Free
   Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
   2111-1307, USA.
*/



	.comm _montgomery_modulo_n,4
	.comm _montgomery_inv_n,4

	.text


	# asm_zero128(a): a=0
	.align	4
	.globl asm_zero128
asm_zero128:
	movl 4(%esp),%edx
	xorl %eax,%eax
	movl %eax,(%edx)
	movl %eax,4(%edx)
	movl %eax,8(%edx)
	movl %eax,12(%edx)
	ret

	# asm_copy128(b,a): b=a
	.align	4
	.globl asm_copy128
asm_copy128:
	movl 4(%esp),%edx
	movl 8(%esp),%ecx
	movl (%ecx),%eax
	movl %eax,(%edx)
	movl 4(%ecx),%eax
	movl %eax,4(%edx)
	movl 8(%ecx),%eax
	movl %eax,8(%edx)
	movl 12(%ecx),%eax
	movl %eax,12(%edx)
	ret

	# asm_sub_n128(b,a): b-=a  mod 2^128
	.align	4
	.globl asm_sub_n128
asm_sub_n128:
	movl 4(%esp),%edx
	movl 8(%esp),%ecx
	movl (%ecx),%eax
	subl %eax,(%edx)
	movl 4(%ecx),%eax
	sbbl %eax,4(%edx)
	movl 8(%ecx),%eax
	sbbl %eax,8(%edx)
	movl 12(%ecx),%eax
	sbbl %eax,12(%edx)
	ret

	# asm_sub128(c,a,b): c=a-b mod N
	.align	4
	.globl asm_sub128
asm_sub128:
	pushl %esi
	movl 12(%esp),%edx
	movl 16(%esp),%esi
	movl 8(%esp),%ecx
	movl (%edx),%eax
	subl (%esi),%eax
	movl %eax,(%ecx)
	movl 4(%edx),%eax
	sbbl 4(%esi),%eax
	movl %eax,4(%ecx)
	movl 8(%edx),%eax
	sbbl 8(%esi),%eax
	movl %eax,8(%ecx)
	movl 12(%edx),%eax
	sbbl 12(%esi),%eax
	movl %eax,12(%ecx)
	jnc sub_end
	movl _montgomery_modulo_n,%edx
	movl (%edx),%eax
	addl %eax,(%ecx)
	movl 4(%edx),%eax
	adcl %eax,4(%ecx)
	movl 8(%edx),%eax
	adcl %eax,8(%ecx)
	movl 12(%edx),%eax
	adcl %eax,12(%ecx)
sub_end:
	popl %esi
	ret


	# asm_half128(a): a/=2 mod N
	.align	4
	.globl asm_half128
asm_half128:
	movl 4(%esp),%ecx
	movl (%ecx),%eax
	testl $1,%eax
	jnz half_odd
	# a is even
	movl 12(%ecx),%eax
	shrl $1,%eax
	movl %eax,12(%ecx)
	movl 8(%ecx),%eax
	shrl $1,%eax
	movl %eax,8(%ecx)
	movl 4(%ecx),%eax
	rcrl $1,%eax
	movl %eax,4(%ecx)
	movl (%ecx),%eax
	rcrl $1,%eax
	movl %eax,(%ecx)
	ret
	# a is odd, compute (a+N)/2
half_odd:
	pushl %esi
	movl _montgomery_modulo_n,%esi
	movl (%esi),%eax
	addl %eax,(%ecx)
	movl 4(%esi),%eax
	adcl %eax,4(%ecx)
	movl 8(%esi),%eax
	adcl %eax,8(%ecx)
	movl 12(%esi),%eax
	adcl 12(%ecx),%eax
	rcrl $1,%eax
	movl %eax,12(%ecx)
	rcrl $1,8(%ecx)
	rcrl $1,4(%ecx)
	rcrl $1,(%ecx)

	popl %esi
	ret

	# asm_diff128(c,a,b): c=|a-b|
	.align	4
	.globl asm_diff128
asm_diff128:
	pushl %esi
	pushl %edi
	pushl %ebx
	movl 20(%esp),%edi
	movl 24(%esp),%esi
	movl 16(%esp),%ebx

	movl 12(%esi),%eax
	cmpl 12(%edi),%eax
	jc b_smaller_a
	jnz a_smaller_b
	movl 8(%esi),%eax
	cmpl 8(%edi),%eax
	jc b_smaller_a
	jnz a_smaller_b
	movl 4(%esi),%eax
	cmpl 4(%edi),%eax
	jc b_smaller_a
	jnz a_smaller_b
	movl (%esi),%eax
	cmpl (%edi),%eax
	jc b_smaller_a
	subl (%edi),%eax
	movl %eax,(%ebx)
	xorl %eax,%eax
	movl %eax,4(%ebx)
	movl %eax,8(%ebx)
	jmp diff_end
a_smaller_b:
	movl (%esi),%eax
	subl (%edi),%eax
	movl %eax,(%ebx)
	movl 4(%esi),%eax
	sbbl 4(%edi),%eax
	movl %eax,4(%ebx)
	movl 8(%esi),%eax
	sbbl 8(%edi),%eax
	movl %eax,8(%ebx)
	movl 12(%esi),%eax
	sbbl 12(%edi),%eax
	movl %eax,12(%ebx)
	jmp diff_end
b_smaller_a:
	movl (%edi),%eax
	subl (%esi),%eax
	movl %eax,(%ebx)
	movl 4(%edi),%eax
	sbbl 4(%esi),%eax
	movl %eax,4(%ebx)
	movl 8(%edi),%eax
	sbbl 8(%esi),%eax
	movl %eax,8(%ebx)
	movl 12(%edi),%eax
	sbbl 12(%esi),%eax
	movl %eax,12(%ebx)
diff_end:
	popl %ebx
	popl %edi
	popl %esi
	ret

	# asm_add128(b,a): b+=a mod N
	.align	4
	.globl asm_add128
asm_add128:
	pushl %esi
	pushl %edi
	movl 16(%esp),%esi
	movl 12(%esp),%edi
	movl (%esi),%eax
	addl %eax,(%edi)
	movl 4(%esi),%eax
	adcl %eax,4(%edi)
	movl 8(%esi),%eax
	adcl %eax,8(%edi)
	movl 12(%esi),%eax
	adcl %eax,12(%edi)

	movl _montgomery_modulo_n,%esi
	jc sub
	movl 12(%esi),%eax
	cmpl 12(%edi),%eax
	jc sub
	jnz add_end
	movl 8(%esi),%eax
	cmpl 8(%edi),%eax
	jc sub
	jnz add_end
	movl 4(%esi),%eax
	cmpl 4(%edi),%eax
	jc sub
	jnz add_end
	movl (%esi),%eax
	cmpl (%edi),%eax
	jnc add_end

sub:
	movl (%esi),%eax
	subl %eax,(%edi)
	movl 4(%esi),%eax
	sbbl %eax,4(%edi)
	movl 8(%esi),%eax
	sbbl %eax,8(%edi)
	movl 12(%esi),%eax
	sbbl %eax,12(%edi)

	#	jnc add_end
	#	movl _montgomery_modulo_n,%esi
	#	movl (%esi),%eax
	#	subl %eax,(%edi)
	#	movl 4(%esi),%eax
	#	sbbl %eax,4(%edi)
	#	movl 8(%esi),%eax
	#	sbbl %eax,8(%edi)
add_end:
	popl %edi
	popl %esi
	ret

	# asm_add128_ui(b,a): b+=a mod N, a is ulong
	.align	4
	.globl asm_add128_ui
asm_add128_ui:
	pushl %esi
	pushl %edi
	movl 12(%esp),%edi
	movl 16(%esp),%eax
	addl %eax,(%edi)
	adcl $0,4(%edi)
	adcl $0,8(%edi)
	adcl $0,12(%edi)
	jnc add_ui_end
	movl _montgomery_modulo_n,%esi
	movl (%esi),%eax
	subl %eax,(%edi)
	movl 4(%esi),%eax
	sbbl %eax,4(%edi)
	movl 8(%esi),%eax
	sbbl %eax,8(%edi)
	movl 12(%esi),%eax
	sbbl %eax,12(%edi)
add_ui_end:
	popl %edi
	popl %esi
	ret

	.IF 0
	.align	4
	.globl asm_cmp128a

asm_cmp128a:
	pushl %esi
	movl 8(%esp),%esi
	movl 12(%esp),%edx
	xorl %eax,%eax
	movl (%esi),%ecx
	subl (%edx),%ecx
	orl %ecx,%eax
	movl 4(%esi),%ecx
	subl 4(%edx),%ecx
	orl %ecx,%eax
	movl 8(%esi),%ecx
	subl 8(%edx),%ecx
	orl %ecx,%eax
	popl %esi
	ret
	.ENDIF

	# asm_mulm128(c,a,b): c=a*b mont-mod N
	.align	4
	.globl asm_mulm128

asm_mulm128:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp
	subl $60,%esp

	# begin of multiplication

	xorl %ebx,%ebx
	movl %ebx,4(%esp)
	movl %ebx,8(%esp)
	movl %ebx,12(%esp)
	movl %ebx,16(%esp)
	movl %ebx,20(%esp)
	movl %ebx,24(%esp)
	movl %ebx,28(%esp)
	movl %ebx,32(%esp)

	movl 84(%esp),%edi
	movl 88(%esp),%esi

	movl (%edi),%ecx
	movl (%esi),%eax
	mull %ecx
	movl %eax,4(%esp)
	movl %edx,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	movl %eax,8(%esp)
	movl $0,%ebx
	adcl %edx,%ebx

	movl 8(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	movl %eax,12(%esp)
	movl $0,%ebx
	adcl %edx,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	movl %eax,16(%esp)
	movl $0,%ebx
	adcl %edx,%ebx
	movl %ebx,20(%esp)
	movl $0,%ebx

	movl 4(%edi),%ecx
	movl (%esi),%eax
	mull %ecx
	addl %eax,8(%esp)
	adcl %edx,%ebx
	addl %ebx,12(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %eax,12(%esp)
	adcl %edx,%ebx
	addl %ebx,16(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	addl %eax,16(%esp)
	adcl %edx,%ebx
	addl %ebx,20(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %eax,20(%esp)
	adcl %edx,%ebx
	movl %ebx,24(%esp)
	movl $0,%ebx

	movl 8(%edi),%ecx
	movl (%esi),%eax
	mull %ecx
	addl %eax,12(%esp)
	adcl %edx,%ebx
	addl %ebx,16(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %eax,16(%esp)
	adcl %edx,%ebx
	addl %ebx,20(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	addl %eax,20(%esp)
	adcl %edx,%ebx
	addl %ebx,24(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %eax,24(%esp)
	adcl %edx,%ebx
	movl %ebx,28(%esp)
	movl $0,%ebx

	movl 12(%edi),%ecx
	movl (%esi),%eax
	mull %ecx
	addl %eax,16(%esp)
	adcl %edx,%ebx
	addl %ebx,20(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %eax,20(%esp)
	adcl %edx,%ebx
	addl %ebx,24(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	addl %eax,24(%esp)
	adcl %edx,%ebx
	addl %ebx,28(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %eax,28(%esp)
	adcl %edx,%ebx
	movl %ebx,32(%esp)
	movl $0,%ebx
	#		jmp debug1
	# end of multiplication

	# begin of reduction

	movl _montgomery_inv_n,%eax
	movl 4(%esp),%ecx
	mull %ecx
	movl %eax,%ecx
	movl _montgomery_modulo_n,%esi

	movl (%esi),%eax
	mull %ecx
	xorl %ebx,%ebx
	addl %eax,4(%esp)
	adcl %edx,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	adcl $0,%edx
	addl %eax,8(%esp)
	adcl %edx,12(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	#	addl %ebx,%eax
	addl %ebx,%edx       # no carry
	addl %eax,12(%esp)
	adcl %edx,16(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %ebx,%edx
	addl %eax,16(%esp)
	adcl %edx,20(%esp)
	adcl $0,24(%esp)
	adcl $0,28(%esp)
	adcl $0,32(%esp)
	#		jmp debug1
	movl _montgomery_inv_n,%eax
	movl 8(%esp),%ecx
	mull %ecx
	movl %eax,%ecx

	movl (%esi),%eax
	mull %ecx
	xorl %ebx,%ebx
	addl %eax,8(%esp)
	adcl %edx,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	adcl $0,%edx
	addl %eax,12(%esp)
	adcl %edx,16(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	#	addl %ebx,%eax
	addl %ebx,%edx        # no carry possible
	addl %eax,16(%esp)
	adcl %edx,20(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %ebx,%edx
	addl %eax,20(%esp)
	adcl %edx,24(%esp)
	adcl $0,28(%esp)
	adcl $0,32(%esp)

	movl _montgomery_inv_n,%eax
	movl 12(%esp),%ecx
	mull %ecx
	movl %eax,%ecx

	movl (%esi),%eax
	mull %ecx
	xorl %ebx,%ebx
	addl %eax,12(%esp)
	adcl %edx,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	adcl $0,%edx
	addl %eax,16(%esp)
	adcl %edx,20(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	#	addl %ebx,%eax
	addl %ebx,%edx            #no carry possible
	addl %eax,20(%esp)
	adcl %edx,24(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %ebx,%edx
	addl %eax,24(%esp)
	adcl %edx,28(%esp)
	adcl $0,32(%esp)

	movl _montgomery_inv_n,%eax
	movl 16(%esp),%ecx
	mull %ecx
	movl %eax,%ecx

	movl (%esi),%eax
	mull %ecx
	xorl %ebx,%ebx
	addl %eax,16(%esp)
	adcl %edx,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	adcl $0,%edx
	addl %eax,20(%esp)
	adcl %edx,24(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	#	addl %ebx,%eax
	addl %ebx,%edx    # no carry
	addl %eax,24(%esp)
	adcl %edx,28(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %ebx,%edx
	addl %eax,28(%esp)
	adcl %edx,32(%esp)

	movl 80(%esp),%edi
	movl 20(%esp),%eax
	movl %eax,(%edi)
	movl 24(%esp),%eax
	movl %eax,4(%edi)
	movl 28(%esp),%eax
	movl %eax,8(%edi)
	movl 32(%esp),%eax
	movl %eax,12(%edi)

	jc subtract
	cmpl 12(%esi),%eax
	jc ende
	jnz subtract
	movl 8(%edi),%eax
	cmpl 8(%esi),%eax
	jc ende
	jnz subtract
	movl 4(%edi),%eax
	cmpl 4(%esi),%eax
	jc ende
	jnz subtract
	movl (%edi),%eax
	cmpl (%esi),%eax
	jc ende

subtract:
	movl (%esi),%eax
	subl %eax,(%edi)
	movl 4(%esi),%eax
	sbbl %eax,4(%edi)
	movl 8(%esi),%eax
	sbbl %eax,8(%edi)
	movl 12(%esi),%eax
	sbbl %eax,12(%edi)
	jmp ende
debug1:
	movl 80(%esp),%edi
	movl 4(%esp),%eax
	movl %eax,(%edi)
	movl 8(%esp),%eax
	movl %eax,4(%edi)
	movl 12(%esp),%eax
	movl %eax,8(%edi)
	movl 16(%esp),%eax
	movl %eax,12(%edi)
	jmp ende
debug2:
	movl 80(%esp),%edi
	movl 20(%esp),%eax
	movl %eax,(%edi)
	movl 24(%esp),%eax
	movl %eax,4(%edi)
	movl 28(%esp),%eax
	movl %eax,8(%edi)
	movl 32(%esp),%eax
	movl %eax,12(%edi)
	jmp ende

ende:
	addl $60,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret



	# asm_sqm128(b,a): b=a*a mont-mod N
	.align	4
	.globl asm_sqm128

asm_sqm128:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp
	subl $80,%esp

	# begin of squaring

	# gemischte produkte
	movl $0,%ebx
	movl %ebx,4(%esp)
	movl %ebx,8(%esp)
	movl %ebx,12(%esp)
	movl %ebx,16(%esp)
	movl %ebx,20(%esp)
	movl %ebx,24(%esp)
	movl %ebx,28(%esp)
	movl %ebx,32(%esp)

	movl 104(%esp),%ebp
	movl (%ebp),%ecx    # a0
	movl 4(%ebp),%eax
	mull %ecx
	movl %eax,8(%esp)     # a0*a1
	movl %edx,%ebx

	movl 8(%ebp),%eax
	mull %ecx
	xorl %edi,%edi
	addl %ebx,%eax
	movl %eax,12(%esp)    # a0*a2
	adcl %edx,%edi

	movl 12(%ebp),%eax
	mull %ecx
	xorl %ebx,%ebx
	addl %edi,%eax
	movl %eax,16(%esp)    # a0*a3
	adcl %edx,%ebx
	movl %ebx,20(%esp)    # no carry !

	movl 4(%ebp),%ecx    # a1
	movl 8(%ebp),%eax
	mull %ecx
	xorl %ebx,%ebx
	movl %ebx,32(%esp)
	addl %eax,16(%esp)    # a1*a2
	adcl %edx,20(%esp)
	adcl $0,%ebx

	movl 12(%ebp),%eax
	mull %ecx
	addl %ebx,%edx
	#	adcl $0,%edx
	addl %eax,20(%esp)    # a1*a3
	adcl %edx,24(%esp)    # no carry !

	movl 8(%ebp),%ecx     # a2
	movl 12(%ebp),%eax
	mull %ecx
	addl %eax,24(%esp)    # a2*a3
	adcl $0,%edx
	movl %edx,28(%esp)    # no carry !
	# gemischte produkte ende
	# jmp debug2sq
	# jetzt *2:
	shll $1,8(%esp)
	rcll $1,12(%esp)
	rcll $1,16(%esp)
	rcll $1,20(%esp)
	rcll $1,24(%esp)
	rcll $1,28(%esp)
	adcl $0,32(%esp)
	# *2 ende

	# jetzt a_i^2 dazuaddieren
	movl (%ebp),%eax
	mull %eax
	movl %eax,4(%esp)
	movl %edx,40(%esp)

	movl 4(%ebp),%eax
	mull %eax
	movl %eax,44(%esp)
	movl %edx,48(%esp)

	movl 8(%ebp),%eax
	mull %eax
	movl %eax,52(%esp)
	movl %edx,56(%esp)

	movl 12(%ebp),%eax
	mull %eax

	movl 40(%esp),%ebx
	addl %ebx,8(%esp)
	movl 44(%esp),%ecx
	adcl %ecx,12(%esp)
	movl 48(%esp),%ebx
	adcl %ebx,16(%esp)
	movl 52(%esp),%ecx
	adcl %ecx,20(%esp)
	movl 56(%esp),%ebx
	adcl %ebx,24(%esp)
	adcl %eax,28(%esp)
	adcl %edx,32(%esp)
	# end of squaring
	#	jmp debug2sq
	# begin of reduction

	movl _montgomery_inv_n,%eax
	movl 4(%esp),%ecx
	mull %ecx
	movl _montgomery_modulo_n,%esi
	movl %eax,%ecx

	movl (%esi),%eax
	mull %ecx
	xorl %ebx,%ebx
	addl %eax,4(%esp)
	adcl %edx,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	adcl $0,%edx
	addl %eax,8(%esp)
	adcl %edx,12(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	#	addl %ebx,%eax
	addl %ebx,%edx      # no carry
	addl %eax,12(%esp)
	adcl %edx,16(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %ebx,%edx
	addl %eax,16(%esp)
	adcl %edx,20(%esp)
	adcl $0,24(%esp)
	adcl $0,28(%esp)
	adcl $0,32(%esp)

	movl _montgomery_inv_n,%eax
	movl 8(%esp),%ecx
	mull %ecx
	movl %eax,%ecx

	movl (%esi),%eax
	mull %ecx
	xorl %ebx,%ebx
	addl %eax,8(%esp)
	adcl %edx,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	adcl $0,%edx
	addl %eax,12(%esp)
	adcl %edx,16(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	#	addl %ebx,%eax
	addl %ebx,%edx         # no carry
	addl %eax,16(%esp)
	adcl %edx,20(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %ebx,%edx
	addl %eax,20(%esp)
	adcl %edx,24(%esp)
	adcl $0,28(%esp)
	adcl $0,32(%esp)

	movl _montgomery_inv_n,%eax
	movl 12(%esp),%ecx
	mull %ecx
	movl %eax,%ecx

	movl (%esi),%eax
	mull %ecx
	xorl %ebx,%ebx
	addl %eax,12(%esp)
	adcl %edx,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	adcl $0,%edx
	addl %eax,16(%esp)
	adcl %edx,20(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	#	addl %ebx,%eax
	addl %ebx,%edx         # no carry
	addl %eax,20(%esp)
	adcl %edx,24(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %ebx,%edx
	addl %eax,24(%esp)
	adcl %edx,28(%esp)
	adcl $0,32(%esp)

	movl _montgomery_inv_n,%eax
	movl 16(%esp),%ecx
	mull %ecx
	movl %eax,%ecx

	movl (%esi),%eax
	mull %ecx
	xorl %ebx,%ebx
	addl %eax,16(%esp)
	adcl %edx,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %ebx,%eax
	adcl $0,%edx
	addl %eax,20(%esp)
	adcl %edx,24(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 8(%esi),%eax
	mull %ecx
	#	addl %ebx,%eax
	addl %ebx,%edx     # no carry
	addl %eax,24(%esp)
	adcl %edx,28(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 12(%esi),%eax
	mull %ecx
	addl %ebx,%edx
	addl %eax,28(%esp)
	adcl %edx,32(%esp)

	movl 100(%esp),%ebp
	movl 20(%esp),%eax
	movl %eax,(%ebp)
	movl 24(%esp),%eax
	movl %eax,4(%ebp)
	movl 28(%esp),%eax
	movl %eax,8(%ebp)
	movl 32(%esp),%eax
	movl %eax,12(%ebp)

	jc subtractsq
	cmpl 12(%esi),%eax
	jc endesq
	jnz subtractsq
	movl 8(%ebp),%eax
	cmpl 8(%esi),%eax
	jc endesq
	jnz subtractsq
	movl 4(%ebp),%eax
	cmpl 4(%esi),%eax
	jc endesq
	jnz subtractsq
	movl (%ebp),%eax
	cmpl (%esi),%eax
	jc endesq

subtractsq:
	movl (%esi),%eax
	subl %eax,(%ebp)
	movl 4(%esi),%eax
	sbbl %eax,4(%ebp)
	movl 8(%esi),%eax
	sbbl %eax,8(%ebp)
	movl 12(%esi),%eax
	sbbl %eax,12(%ebp)
	jmp endesq
debug1sq:
	movl 100(%esp),%ebp
	movl 4(%esp),%eax
	movl %eax,(%ebp)
	movl 8(%esp),%eax
	movl %eax,4(%ebp)
	movl 12(%esp),%eax
	movl %eax,8(%ebp)
	movl 16(%esp),%eax
	movl %eax,12(%ebp)
	jmp endesq
debug2sq:
	movl 100(%esp),%ebp
	movl 20(%esp),%eax
	movl %eax,(%ebp)
	movl 24(%esp),%eax
	movl %eax,4(%ebp)
	movl 28(%esp),%eax
	movl %eax,8(%ebp)
	movl 32(%esp),%eax
	movl %eax,12(%ebp)

endesq:
	addl $80,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret
