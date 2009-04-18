# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# Written by T. Kleinjung

relptr       .ASSIGNC  "40(%esp)"
minus        .ASSIGNC  "44(%esp)"
qx           .ASSIGNC  "48(%esp)"
ulqx         .ASSIGNC  "52(%esp)"


.comm mpqs_nFBk_1,2
.comm mpqs_td_begin,2
.comm mpqs_sievebegin,2
.comm mpqs_FB_inv_info,4
.comm mpqs_FB_start,4
.comm mpqs_256_inv_table,4
.comm mpqs_FB_inv,4
.comm stat_asm_div,4


	.align 16
.globl asm_td
	.type    asm_td,@function
# mm0: p,p
# mm1: inv,inv
# mm2: s1,s2
# mm3: ind,ind
# mm4: computing
# mm5: 0

asm_td:
	pushl %edi
	pushl %esi
	pushl %ebx
	pushl %ebp
	subl $20,%esp

	movl \&relptr,%ebp    # rel[i]
	movl (%ebp),%eax      # ind
	andl $0x0000ffff,%eax
	movl %eax,%ebx
	shll $16,%eax
	orl %ebx,%eax
	movd %eax,%mm3        # ind,ind
	psllq $32,%mm3
	movd %eax,%mm5
	paddd %mm5,%mm3        # ind,ind,ind,ind
	pxor %mm5,%mm5

#	movl $0,%edx
	movzwl 8(%ebp),%edx      # nr

#	movl $0,%ecx
	movzwl mpqs_td_begin,%ecx

	movl $mpqs_FB_inv_info,%esi
	movl $mpqs_FB_start,%edi

# prime mpqs_FB[1]
	movd 4(%esi),%mm0
	movd 12(%esi),%mm1
	movd 4(%edi),%mm2
	movq %mm0,%mm4
	psubw %mm2,%mm4
	paddw %mm3,%mm4        # ind+p-s1,ind+p-s2
	pmullw %mm1,%mm4
	pmulhw %mm0,%mm4
	pcmpeqw %mm5,%mm4
	movd %mm4,%eax
	orl $0,%eax
	jz loop2a
# found divisor mpqs_FB[1]
	movw mpqs_nFBk_1,%ax
	incw %ax
	movw %ax,10(%ebp,%edx,2)
	incl %edx

loop2a:
	leal 16(%esi),%esi
	leal 8(%edi),%edi
	movq (%esi),%mm0
	movq 8(%esi),%mm1
	movq (%edi),%mm2

loop2:
	subl $2,%ecx
	jz prod

	movq %mm0,%mm4
	psubw %mm2,%mm4
	paddw %mm3,%mm4        # ind+p-s1,ind+p-s2,ind+P-S1,ind+P-S2
	pmullw %mm1,%mm4
	movq 8(%edi),%mm2
	leal 16(%esi),%esi
	leal 8(%edi),%edi
	pmulhw %mm0,%mm4
	movq (%esi),%mm0
	movq 8(%esi),%mm1
	pcmpeqw %mm5,%mm4
	pmovmskb %mm4,%eax
	testl %eax,%eax
	jz loop2
	movw mpqs_nFBk_1,%bx
	addw mpqs_td_begin,%bx
	subw %cx,%bx
# found divisor
	testl $15,%eax
	jz testsecond

	movw %bx,10(%ebp,%edx,2)
	incl %edx
testsecond:
	testl $240,%eax
	jz loop2
	incw %ebx
	movw %bx,10(%ebp,%edx,2)
	incl %edx
	jmp loop2

prod:
	movl $mpqs_FB,%esi
#	movl $0,%ebx
	movzwl mpqs_nFBk_1,%ebx
	addl %ebx,%ebx
	addl %ebx,%ebx
	subl %ebx,%esi
	movl $0,%ecx
	movl %edx,4(%esp)              # nr
	xorl %edx,%edx
	movl $1,%eax
	movl $0,%ebx
	movl %ebx,%edi
prodloop:
	addl %edi,%edx
	cmpl 4(%esp),%ecx
	jnc prodend
	movzwl 10(%ebp,%ecx,2),%ebx
	movl %eax,%edi
	movzwl (%esi,%ebx,4),%ebx
	movl %edx,%eax
	mull %ebx
	xchg %eax,%edi
	incl %ecx
	mull %ebx
	jmp prodloop
	
prodend:
	movl %eax,12(%esp)
	movl %edx,16(%esp)

	movl 4(%esp),%ebx          # nr
	movl \&minus,%eax
	testl $1,%eax
	jz positive
	movw $0,10(%ebp,%ebx,2)
	incl %ebx

positive:
	movl \&qx,%edi
	movl (%edi),%eax
	movl 4(%edi),%edx
posloop:
	testl $0x00000001,%eax
	jnz odd
	incl %ebx
	cmpl $27,%ebx
	jnc gotonext
	movw mpqs_nFBk_1,%cx
	movw %cx,8(%ebp,%ebx,2)
	shrl $1,%edx
	rcrl $1,%eax
	jmp posloop

odd:
	movw %bx,8(%ebp)

	cmpl $0,16(%esp)
	jnz division
	cmpl 12(%esp),%edx
	jnc gotonext

division:
	movl %eax,8(%esp)         # ax
	movl 12(%esp),%ebx         # ay
	movl %ebx,%edx
	movl $0,%ecx
	andl $0x000000ff,%edx
	shrl $1,%edx
	movl $mpqs_256_inv_table,%edi
	movb (%edi,%edx),%cl         # inv
	movl %ebx,%eax
	mull %ecx
	andl $0x0000ff00,%eax
	mull %ecx
	subl %eax,%ecx
	movl %ebx,%eax
	mull %ecx
	andl $0xffff0000,%eax
	mull %ecx
	subl %eax,%ecx

	movl 8(%esp),%eax
	mull %ecx

# trial divison of sieved primes
	movl $mpqs_FB_inv,%edi
#	movl $0,%ecx
	movzwl mpqs_nFBk_1,%ecx
	addl %ecx,%ecx
	addl %ecx,%ecx
	subl %ecx,%edi
	movl $0,%ecx
	movl 4(%esp),%edx
	movl %edx,(%esp)
	movl $0,%ebx
	movl %eax,8(%esp)
tdloop:
	movl 8(%esp),%eax
	cmpl $0,(%esp)
	jz tdend
	movl (%esp),%edx
	movzwl 8(%ebp,%edx,2),%ebx  # bx: ii
	decl (%esp)
	movzwl (%esi,%ebx,4),%ecx   # cx: p
	movl (%edi,%ebx,4),%edx  # edx: inv
	movl %edx,12(%esp)
divloop:
	mull %edx
	movl %eax,16(%esp)        # rr
	mull %ecx
	testl %edx,%edx
	jnz tdloop
	movw 8(%ebp),%dx
	cmpw $27,%dx
	jnc gotonext
	movw %bx,10(%ebp,%edx,2)
	incw 8(%ebp)
	movl 16(%esp),%eax
	movl 12(%esp),%edx
	movl %eax,8(%esp)
	jmp divloop


tdend:
# trial division of mpqs_FBk-primes
	cmpl $1,%eax
	jz end

	movl $mpqs_FBk_inv,%edi
#	movl $0,%ebx
	movl $mpqs_FBk,%esi
	xorl %ecx,%ecx
	movzwl mpqs_nFBk,%ebx
	incl %ebx
tdloopk:
	decl %ebx
	movl 8(%esp),%eax
	cmpl $0,%ebx
	jz tdendk
	movzwl -2(%esi,%ebx,2),%ecx   # cx: p
	movl -4(%edi,%ebx,4),%edx  # edx: inv
	movl %edx,12(%esp)
divloopk:
	mull %edx
	movl %eax,16(%esp)        # rr
	mull %ecx
	testl %edx,%edx
	jnz tdloopk
	movw 8(%ebp),%dx
	cmpw $27,%dx
	jnc gotonext
	movw %bx,10(%ebp,%edx,2)
	incw 8(%ebp)
	movl 16(%esp),%eax
	movl 12(%esp),%edx
	movl %eax,8(%esp)
	jmp divloopk

tdendk:
# trial division of mpqs_FB_Adiv-primes
	cmpl $1,%eax
	jz end

	movl $mpqs_FB_A_inv,%edi
	movl $0,%ecx
	movl $mpqs_Adiv_all,%esi
#	movl $0,%ebx
        movw mpqs_nFB,%cx
        addw mpqs_nFBk,%cx
	movzwl mpqs_nAdiv_total,%ebx
	incl %ebx
	movl %ecx,(%esp)           # mpqs_nFB+mpqs_nFBk
#	xorl %ecx,%ecx
tdloopa:
	decl %ebx
	movl 8(%esp),%eax
	cmpl $0,%ebx
	jz tdenda
	movzwl -2(%esi,%ebx,2),%ecx   # cx: p
	movl -4(%edi,%ebx,4),%edx  # edx: inv
	movl %edx,12(%esp)
divloopa:
	mull %edx
	movl %eax,16(%esp)        # rr
	mull %ecx
	testl %edx,%edx
	jnz tdloopa
	movw 8(%ebp),%dx
	cmpw $27,%dx
	jnc gotonext
	addl (%esp),%ebx
	movw %bx,10(%ebp,%edx,2)
	incw 8(%ebp)
	subl (%esp),%ebx
	movl 16(%esp),%eax
	movl 12(%esp),%edx
	movl %eax,8(%esp)
	jmp divloopa

tdenda:


end:
	movl 8(%esp),%eax
	movl \&ulqx,%edi
	movl %eax,(%edi)

	xorl %eax,%eax
	emms
	addl $20,%esp
	popl %ebp
	popl %ebx
	popl %esi
	popl %edi
	ret

gotonext:
	movl $1,%eax
	emms
	addl $20,%esp
	popl %ebp
	popl %ebx
	popl %esi
	popl %edi
	ret
# only used for debugging
dbg:
	addl $40,4(%esp)
	movl \&ulqx,%edi
	movl 8(%esp),%edx
	movl %edx,(%edi)
	jmp end
	.END
