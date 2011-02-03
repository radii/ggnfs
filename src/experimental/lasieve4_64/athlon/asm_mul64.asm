# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# Written by T. Kleinjung

# 64 Bit


f1 			.ASSIGNC "56(%esp)"
f2			.ASSIGNC "60(%esp)"
prod			.ASSIGNC "52(%esp)"

.comm montgomery_modulo_n,4
.comm montgomery_inv_n,4

.text

# asm_zero64(a): a=0
	.align 4
.globl asm_zero64
	.type    asm_zero64,@function
asm_zero64:
	movl 4(%esp),%edx
	xorl %eax,%eax
	movl %eax,(%edx)
	movl %eax,4(%edx)
	ret

# asm_copy64(b,a): b=a
	.align 4
.globl asm_copy64
	.type    asm_copy64,@function
asm_copy64:
	movl 4(%esp),%edx
	movl 8(%esp),%ecx
	movl (%ecx),%eax
	movl %eax,(%edx)
	movl 4(%ecx),%eax
	movl %eax,4(%edx)
	ret


# asm_sub64(c,a,b): c=a-b mod N
	.align 4
.globl asm_sub64
	.type    asm_sub64,@function
asm_sub64:
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
	jnc sub_end
	movl montgomery_modulo_n,%edx
	movl (%edx),%eax
	addl %eax,(%ecx)
	movl 4(%edx),%eax
	adcl %eax,4(%ecx)
sub_end:
	popl %esi
	ret


# asm_half64(a): a/=2 mod N
	.align 4
.globl asm_half64
	.type    asm_half64,@function
asm_half64:
	movl 4(%esp),%ecx
	movl (%ecx),%eax
	testl $1,%eax
	jnz half_odd
# a is even
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
	movl montgomery_modulo_n,%esi
	movl (%esi),%eax
	addl (%ecx),%eax
	movl %eax,(%ecx)
	movl 4(%esi),%eax
	adcl 4(%ecx),%eax
	movl %eax,4(%ecx)

	rcrl $1,4(%ecx)
	rcrl $1,(%ecx)

	popl %esi
	ret


# asm_sub_n64(b,a): b-=a  mod 2^64
        .align 4
.globl asm_sub_n64
        .type    asm_sub_n64,@function
asm_sub_n64:
        movl 4(%esp),%edx
        movl 8(%esp),%ecx
        movl (%ecx),%eax
        subl %eax,(%edx)
        movl 4(%ecx),%eax
        sbbl %eax,4(%edx)
        ret


# asm_diff64(c,a,b): c=|a-b|
        .align 4
.globl asm_diff64
        .type    asm_diff64,@function
asm_diff64:
	pushl %esi
	pushl %edi
	pushl %ebx
	movl 20(%esp),%edi
	movl 24(%esp),%esi
	movl 16(%esp),%ebx
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
	jmp diff_end
a_smaller_b:
	movl (%esi),%eax
	subl (%edi),%eax
	movl %eax,(%ebx)
	movl 4(%esi),%eax
	sbbl 4(%edi),%eax
	movl %eax,4(%ebx)
	jmp diff_end
b_smaller_a:
	movl (%edi),%eax
	subl (%esi),%eax
	movl %eax,(%ebx)
	movl 4(%edi),%eax
	sbbl 4(%esi),%eax
	movl %eax,4(%ebx)
diff_end:
	popl %ebx
	popl %edi
	popl %esi
	ret

# asm_add64(b,a): b+=a mod N
        .align 4
.globl asm_add64
asm_add64:
	pushl %esi
	pushl %edi
	movl 16(%esp),%esi
	movl 12(%esp),%edi
	movl (%esi),%eax
	addl %eax,(%edi)
	movl 4(%esi),%eax
	adcl %eax,4(%edi)
	movl montgomery_modulo_n,%esi
	jc sub
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
add_end:
	popl %edi
	popl %esi
	ret

# asm_add64_ui(b,a): b+=a mod N, a is ulong
        .align 4
.globl asm_add64_ui
asm_add64_ui:
        pushl %esi
	pushl %edi
        movl 12(%esp),%edi
        movl 16(%esp),%eax
        addl %eax,(%edi)
        adcl $0,4(%edi)
        jnc add_ui_end
        movl montgomery_modulo_n,%esi
        movl (%esi),%eax
        subl %eax,(%edi)
        movl 4(%esi),%eax
        sbbl %eax,4(%edi)
add_ui_end:
	popl %edi
        popl %esi
        ret

# asm_gcd64:



# asm_mulm64(): ??? mod N
	.align 4
.globl asm_mulm64
	.type	 asm_mulm64,@function

asm_mulm64:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp
	subl $32,%esp
# begin
	xorl %ebx,%ebx
	movl %ebx,4(%esp)
	movl %ebx,8(%esp)
	movl %ebx,12(%esp)
	movl %ebx,16(%esp)

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

	movl 4(%edi),%ecx
	movl (%esi),%eax
	mull %ecx
	addl %eax,8(%esp)
	adcl %edx,%ebx
	movl %ebx,12(%esp)
	movl $0,%ebx
	adcl $0,%ebx

	movl 4(%esi),%eax
	mull %ecx
	addl %eax,12(%esp)
	adcl %edx,%ebx
	movl %ebx,16(%esp)

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
	adcl $0,16(%esp)

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

	movl \&prod,%edi
	movl 12(%esp),%eax
	movl %eax,(%edi)
	movl 16(%esp),%eax
	movl %eax,4(%edi)

	jc subtract
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
	addl $32,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret

f                       .ASSIGNC "56(%esp)"
sq                      .ASSIGNC "52(%esp)"

	.align 4
.globl asm_sqm64
	.type	 asm_sqm64,@function

asm_sqm64:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp
	subl $32,%esp

# begin of squaring

	xorl %ebx,%ebx
	movl %ebx,4(%esp)
	movl %ebx,8(%esp)
	movl %ebx,12(%esp)
	movl %ebx,16(%esp)

	movl \&f,%edi
	movl (%edi),%eax
	mull %eax
	movl %eax,4(%esp)
	movl %edx,8(%esp)

	movl 4(%edi),%eax
	mull %eax
	movl %eax,12(%esp)
	movl %edx,16(%esp)

	movl 4(%edi),%ecx
	movl (%edi),%eax
	mull %ecx
	addl %eax,%eax
	adcl %edx,%edx
	movl $0,%ebx
	adcl $0,%ebx

	addl %eax,8(%esp)
	adcl %edx,12(%esp)
	adcl %ebx,16(%esp)

# end of squaring

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
	adcl $0,16(%esp)

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

	movl \&sq,%edi
	movl 12(%esp),%eax
	movl %eax,(%edi)
	movl 16(%esp),%eax
	movl %eax,4(%edi)

	jc subtractsq
	cmpl 4(%esi),%eax
	jc endesq
	jnz subtractsq
	movl (%edi),%eax
	cmpl (%esi),%eax
	jc endesq

subtractsq:
        movl (%esi),%eax
        subl %eax,(%edi)
        movl 4(%esi),%eax
        sbbl %eax,4(%edi)
#	jmp endesq
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
endesq:
	addl $32,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret
	.END
