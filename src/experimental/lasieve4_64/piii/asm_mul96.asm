# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# Written by T. Kleinjung

# 96 Bit


f1 			.ASSIGNC "64(%esp)"
f2			.ASSIGNC "68(%esp)"
prod			.ASSIGNC "60(%esp)"


.comm montgomery_modulo_n,4
.comm montgomery_inv_n,4

.text


# asm_zero96(a): a=0
	.align 4
.globl asm_zero96
	.type    asm_zero96,@function
asm_zero96:
	movl 4(%esp),%edx
	xorl %eax,%eax
	movl %eax,(%edx)
	movl %eax,4(%edx)
	movl %eax,8(%edx)
	ret

# asm_copy96(b,a): b=a
	.align 4
.globl asm_copy96
	.type    asm_copy96,@function
asm_copy96:
	movl 4(%esp),%edx
	movl 8(%esp),%ecx
	movl (%ecx),%eax
	movl %eax,(%edx)
	movl 4(%ecx),%eax
	movl %eax,4(%edx)
	movl 8(%ecx),%eax
	movl %eax,8(%edx)
	ret

# asm_sub_n96(b,a): b-=a  mod 2^96
	.align 4
.globl asm_sub_n96
	.type    asm_sub_n96,@function
asm_sub_n96:
	movl 4(%esp),%edx
	movl 8(%esp),%ecx
	movl (%ecx),%eax
	subl %eax,(%edx)
	movl 4(%ecx),%eax
	sbbl %eax,4(%edx)
	movl 8(%ecx),%eax
	sbbl %eax,8(%edx)
	ret

# asm_sub96(c,a,b): c=a-b mod N
	.align 4
.globl asm_sub96
	.type    asm_sub96,@function
asm_sub96:
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
	jnc sub_end
	movl montgomery_modulo_n,%edx
	movl (%edx),%eax
	addl %eax,(%ecx)
	movl 4(%edx),%eax
	adcl %eax,4(%ecx)
	movl 8(%edx),%eax
	adcl %eax,8(%ecx)
sub_end:
	popl %esi
	ret


# asm_half96(a): a/=2 mod N
	.align 4
.globl asm_half96
	.type    asm_half96,@function
asm_half96:
	movl 4(%esp),%ecx
	movl (%ecx),%eax
	testl $1,%eax
	jnz half_odd
# a is even
#	movl 8(%ecx),%eax
#	shrl $1,%eax
#	movl %eax,8(%ecx)
#	movl 4(%ecx),%eax
#	rcrl $1,%eax
#	movl %eax,4(%ecx)
#	movl (%ecx),%eax
#	rcrl $1,%eax
#	movl %eax,(%ecx)
	shrl $1,8(%ecx)
	rcrl $1,4(%ecx)
	rcrl $1,(%ecx)
	ret
# a is odd, compute (a+N)/2
half_odd:
	pushl %esi
	movl montgomery_modulo_n,%esi
	movl (%esi),%eax
	addl %eax,(%ecx)
	movl 4(%esi),%eax
	adcl %eax,4(%ecx)
	movl 8(%esi),%eax
	adcl 8(%ecx),%eax
	rcrl $1,%eax
	movl %eax,8(%ecx)
	rcrl $1,4(%ecx)
	rcrl $1,(%ecx)

	popl %esi
	ret

# asm_diff96(c,a,b): c=|a-b|
	.align 4
.globl asm_diff96
	.type    asm_diff96,@function
asm_diff96:
	pushl %esi
	pushl %edi
	pushl %ebx
	movl 20(%esp),%edi
	movl 24(%esp),%esi
	movl 16(%esp),%ebx
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
diff_end:
	popl %ebx
	popl %edi
	popl %esi
	ret

# asm_add96(b,a): b+=a mod N
	.align 4
.globl asm_add96
	.type    asm_add96,@function
asm_add96:
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

	movl montgomery_modulo_n,%esi
	jc sub
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

#	jnc add_end
#	movl montgomery_modulo_n,%esi
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

# asm_add96_ui(b,a): b+=a mod N, a is ulong
	.align 4
.globl asm_add96_ui
	.type    asm_add96_ui,@function
asm_add96_ui:
	pushl %esi
	pushl %edi
	movl 12(%esp),%edi
	movl 16(%esp),%eax
	addl %eax,(%edi)
	adcl $0,4(%edi)
	adcl $0,8(%edi)
	jnc add_ui_end
	movl montgomery_modulo_n,%esi
	movl (%esi),%eax
	subl %eax,(%edi)
	movl 4(%esi),%eax
	sbbl %eax,4(%edi)
	movl 8(%esi),%eax
	sbbl %eax,8(%edi)
add_ui_end:
	popl %edi
	popl %esi
	ret

.IF 0
	.align 4
.globl asm_cmp96a
	.type    asm_cmp96a,@function

asm_cmp96a:
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

# asm_mulm96(c,a,b): c=a*b mont-mod N
	.align 4
.globl asm_mulm96
	.type	 asm_mulm96,@function

asm_mulm96:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp
	subl $40,%esp

# begin of multiplication

	xorl %ebx,%ebx
	movl %ebx,4(%esp)
	movl %ebx,8(%esp)
	movl %ebx,12(%esp)
	movl %ebx,16(%esp)
	movl %ebx,20(%esp)
	movl %ebx,24(%esp)

	movl \&f1,%edi
	movl \&f2,%esi

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
	movl %ebx,16(%esp)
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
	movl %ebx,20(%esp)
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
	movl %ebx,24(%esp)
	movl $0,%ebx

# end of multiplication

# begin of reduction

	movl montgomery_inv_n,%eax
	movl 4(%esp),%ecx
	mull %ecx
	movl %eax,%ecx
	movl montgomery_modulo_n,%esi

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
	addl %ebx,%edx
	addl %eax,12(%esp)
	adcl %edx,16(%esp)
	adcl $0,20(%esp)
	adcl $0,24(%esp)

	movl montgomery_inv_n,%eax
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
	addl %ebx,%edx
	addl %eax,16(%esp)
	adcl %edx,20(%esp)
	adcl $0,24(%esp)

	movl montgomery_inv_n,%eax
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
	addl %ebx,%edx
	addl %eax,20(%esp)
	adcl %edx,24(%esp)

	movl \&prod,%edi
	movl 16(%esp),%eax
	movl %eax,(%edi)
	movl 20(%esp),%eax
	movl %eax,4(%edi)
	movl 24(%esp),%eax
	movl %eax,8(%edi)

	jc subtract
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
#	jmp ende
#debug:
#	movl \&prod,%edi
#	movl 4(%esp),%eax
#	movl %eax,(%edi)
#	movl 8(%esp),%eax
#	movl %eax,4(%edi)
#	movl 12(%esp),%eax
#	movl %eax,8(%edi)
#	movl 16(%esp),%eax
#	movl %eax,12(%edi)
ende:
	addl $40,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret

f                       .ASSIGNC "84(%esp)"
sq                      .ASSIGNC "80(%esp)"
tmp0                    .ASSIGNC "52(%esp)"
tmp1                    .ASSIGNC "56(%esp)"


# asm_sqm96(b,a): b=a*a mont-mod N
	.align 4
.globl asm_sqm96
	.type	 asm_sqm96,@function

asm_sqm96:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp
	subl $60,%esp

# begin of squaring

.IF 0
	movl \&f,%ebp
	movl (%ebp),%ecx       # a0
	movl 4(%ebp),%eax
	mull %ecx              # a0*a1
	movl 8(%ebp),%esi      # a2
	xorl %ebx,%ebx
	movl %ebx,24(%esp)
#
	movl %eax,8(%esp)
	movl %edx,12(%esp)

	movl %ecx,%eax
	mull %esi              # a0*a2
	movl 4(%ebp),%edi      # a1
#
	movl %eax,%ebx
	movl %edx,16(%esp)

	movl %edi,%eax
	mull %esi              # a1*a2
	addl %ebx,12(%esp)
	adcl $0,16(%esp)
#
	movl %eax,%ebx
	movl %edx,20(%esp)

	movl %ecx,%eax
	mull %ecx               # a0*a0
	addl %ebx,16(%esp)
	adcl $0,20(%esp)
	movl 8(%esp),%ebx
	movl 12(%esp),%ecx
#
	movl %eax,4(%esp)
	movl %edx,8(%esp)

	movl %edi,%eax
	mull %edi               # a1*a1
	addl %ebx,%ebx
	adcl %ecx,%ecx
	sbbl %edi,%edi          # store carry of multiplication by 2
	addl %ebx,8(%esp)
	adcl $0,%ecx
	movl %ecx,12(%esp)
	sbbl %ebx,%ebx          # store carry of addition
#
	movl %eax,36(%esp)
	movl %edx,40(%esp)

	movl %esi,%eax
	mull %esi               # a2*a2
	addl $1,%edi            # set carry
	movl 16(%esp),%ecx
	adcl %ecx,16(%esp)
	movl 20(%esp),%ecx
	adcl %ecx,20(%esp)
	movl $0,%ecx
	adcl %ecx,24(%esp)

	movl %eax,44(%esp)
	movl %edx,48(%esp)

	addl $1,%ebx
	adcl $0,40(%esp)
	adcl $0,44(%esp)
	adcl $0,48(%esp)


	movl 36(%esp),%eax
	addl %eax,12(%esp)
	movl 40(%esp),%eax
	adcl %eax,16(%esp)
	movl 44(%esp),%eax
	adcl %eax,20(%esp)
	movl 48(%esp),%eax
	adcl %eax,24(%esp)
.ENDIF

.IF 0
# fehlerhaft !!
	movl \&f,%ebp
	movl (%ebp),%ecx        # a0
	movl 4(%ebp),%eax       # a1
	mull %ecx               # a0*a1
	movl 8(%ebp),%esi       # a2
#
	movl %eax,8(%esp)
	movl %edx,12(%esp)

	movl %ecx,%eax
	mull %esi               # a0*a2
	movl 4(%ebp),%edi       # a1
#
	movl %eax,%ebx
	movl %edx,16(%esp)

	movl %edi,%eax
	mull %esi               # a1*a2
	addl %ebx,12(%esp)
	adcl $0,16(%esp)
#
	movl %eax,%ebx
	movl %edx,20(%esp)
	movl %ecx,%eax
	mull %ecx               # a0*a0
	addl %ebx,16(%esp)
	adcl $0,20(%esp)
	movl 8(%esp),%ebx
	movl 12(%esp),%ecx
#
	movl %eax,4(%esp)
	movl %edx,8(%esp)

	movl %edi,%eax
	mull %edi               # a1*a1
	addl %ebx,%ebx
	adcl %ecx,%ecx
	sbbl %edi,%edi          # store carry of multiplication by 2
	addl %ebx,8(%esp)
	adcl $0,%ecx
	movl %ecx,12(%esp)
	sbbl %ebx,%ebx          # store carry of addition
#
	movl %eax,\&tmp0
	movl %edx,\&tmp1

	movl %esi,%eax
	mull %esi               # a2*a2
	addl $1,%edi            # set carry
	rcll $1,16(%esp)
	movl $0,%ecx
	rcll $1,20(%esp)
	adcl $0,%ecx            # store carry
	addl $1,%ebx            # set carry of addition
	movl \&tmp0,%ebx
	adcl $0,%ebx
	movl \&tmp1,%edi
	adcl $0,%edi            # no carry
	addl %ebx,12(%esp)
	adcl %edi,16(%esp)
#
	adcl $0,%eax
	adcl $0,%edx
	addl $1,%ecx
	adcl $0,%edx
	addl 20(%esp),%eax
	adcl $0,%edx
	movl %eax,20(%esp)
	movl %edx,24(%esp)
.ENDIF


.IF 1
	movl \&f,%ebp
	movl (%ebp),%ecx
	movl 4(%ebp),%eax
	mull %ecx
	movl 8(%ebp),%esi      # a2
	xorl %ebx,%ebx
	movl %ebx,24(%esp)
	movl %eax,8(%esp)
	movl %edx,%ebx

	movl %ecx,%eax
	mull %esi
	xorl %edi,%edi
	addl %ebx,%eax
	movl %eax,12(%esp)
	adcl %edx,%edi

	movl 4(%ebp),%eax
	mull %esi
	xorl %ecx,%ecx
	addl %edi,%eax
	adcl $0,%edx

	shll $1,8(%esp)
	rcll $1,12(%esp)
	rcll $1,%eax
	movl %eax,16(%esp)
	rcll $1,%edx
	movl %edx,20(%esp)
	adcl %ecx,24(%esp)

	movl (%ebp),%eax
	mull %eax
	movl %eax,4(%esp)
	movl %edx,32(%esp)

	movl 4(%ebp),%eax
	mull %eax
	movl %eax,36(%esp)
	movl %edx,40(%esp)

	movl 8(%ebp),%eax
	mull %eax

	movl 32(%esp),%ebx
	addl %ebx,8(%esp)
	movl 36(%esp),%ecx
	adcl %ecx,12(%esp)
	movl 40(%esp),%ebx
	adcl %ebx,16(%esp)
	adcl %eax,20(%esp)
	adcl %edx,24(%esp)
.ENDIF
# end of squaring

# begin of reduction

	movl montgomery_inv_n,%eax
	movl 4(%esp),%ecx
	mull %ecx
	movl montgomery_modulo_n,%esi
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
	addl %ebx,%edx
	addl %eax,12(%esp)
	adcl %edx,16(%esp)
	adcl $0,20(%esp)
	adcl $0,24(%esp)

	movl montgomery_inv_n,%eax
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
	addl %ebx,%edx
	addl %eax,16(%esp)
	adcl %edx,20(%esp)
	adcl $0,24(%esp)

	movl montgomery_inv_n,%eax
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
	addl %ebx,%edx
	addl %eax,20(%esp)
	adcl %edx,24(%esp)

	movl \&sq,%ebp
	movl 16(%esp),%eax
	movl %eax,(%ebp)
	movl 20(%esp),%eax
	movl %eax,4(%ebp)
	movl 24(%esp),%eax
	movl %eax,8(%ebp)

	jc subtractsq
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
#	jmp ende
#debug:
#	movl \&prod,%ebp
#	movl 4(%esp),%eax
#	movl %eax,(%ebp)
#	movl 8(%esp),%eax
#	movl %eax,4(%ebp)
#	movl 12(%esp),%eax
#	movl %eax,8(%ebp)
#	movl 16(%esp),%eax
#	movl %eax,12(%ebp)
endesq:
	addl $60,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret
	.END
