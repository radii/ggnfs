# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# Written by T. Kleinjung

.comm mpqs_nFBk_1,2
.comm mpqs_td_begin,2
.comm mpqs_sievebegin,2
.comm mpqs_FB_inv_info,4
.comm mpqs_FB,4
.comm mpqs_FB_start,4
.comm mpqs_sievearray,4
.comm mpqs_sievelen,4



	.align 16
.globl asm_sieve
	.type    asm_sieve,@function
# 4(%esp)  mpqs_FB-ptr
# 8(%esp)  mpqs_FB_start-ptr
# 12(%esp) mpqs_FB_log-ptr
# 16(%esp) mpqs_sievelen/4


asm_sieve:
	pushl %edi
	pushl %esi
	pushl %ebx
	pushl %ebp

	subl $20,%esp

	movl $mpqs_FB,%eax
#	movl $0,%ebx
	movzwl mpqs_sievebegin,%ebx
	addl %ebx,%ebx
	addl %ebx,%ebx
	addl %ebx,%eax
	movl %eax,4(%esp)

	movl $mpqs_FB_start,%eax
	addl %ebx,%eax
	movl %eax,8(%esp)

	movl $mpqs_FB_log,%eax
#	movl $0,%ebx
	movzwl mpqs_sievebegin,%ebx
	addl %ebx,%eax
	movl %eax,12(%esp)

	movl mpqs_sievelen,%eax
	shrl $2,%eax
	movl %eax,16(%esp)

#	xorl %edx,%edx
#	xorl %ebx,%ebx
#	xorl %ecx,%ecx
#	xorl %eax,%eax

mainloop:
	movl 4(%esp),%esi
	movzwl (%esi),%ecx
	leal 4(%esi),%esi
	movl %esi,4(%esp)    # p 
	cmpl %ecx,16(%esp)
	jc end

#	movl $0,%ebx
	movl 8(%esp),%esi
	movzwl (%esi),%ebx
	movzwl 2(%esi),%edx
	leal 4(%esi),%esi
	movl %esi,8(%esp)    # s1, s2

	movl 12(%esp),%esi
	movb (%esi),%al
	leal 1(%esi),%esi
	movl %esi,12(%esp)   # lo

	movl mpqs_sievelen,%esi
	movl %ecx,%edi
	shll $2,%edi
	subl %edi,%esi
	movl mpqs_sievearray,%edi
	addl %edi,%esi

loop4:
	addb %al,(%edi,%ebx)
	addb %al,(%edi,%edx)
	leal (%edi,%ecx),%edi
	addb %al,(%edi,%ebx)
	addb %al,(%edi,%edx)
	leal (%edi,%ecx),%edi
	addb %al,(%edi,%ebx)
	addb %al,(%edi,%edx)
	leal (%edi,%ecx),%edi
	addb %al,(%edi,%ebx)
	addb %al,(%edi,%edx)
	leal (%edi,%ecx),%edi
	cmpl %esi,%edi
	jc loop4

	addl %ecx,%esi
	addl %ecx,%esi
	cmpl %esi,%edi
	jnc check

	addb %al,(%edi,%ebx)
	addb %al,(%edi,%edx)
	leal (%edi,%ecx),%edi
	addb %al,(%edi,%ebx)
	addb %al,(%edi,%edx)
	leal (%edi,%ecx),%edi
check:
	addl %ecx,%esi
	cmpl %esi,%edi
	jnc check1
	addb %al,(%edi,%ebx)
	addb %al,(%edi,%edx)
	addl %ecx,%edi
check1:
	addl %ecx,%esi
	addl %edi,%ebx
	cmpl %esi,%ebx
	jnc check2
	addb %al,(%ebx)
check2:
	addl %edx,%edi
	cmpl %esi,%edi
	jnc loopend
	addb %al,(%edi)
loopend:
	jmp mainloop


end:
        addl $20,%esp
	popl %ebp
	popl %ebx
	popl %esi
	popl %edi
	ret

        .align 16
.globl asm_sieve1
        .type    asm_sieve1,@function

asm_sieve1:
	pushl %edi
	pushl %esi
	pushl %ebx
	pushl %ebp

	subl $20,%esp

	movl $mpqs_FB,%eax
#	movl $0,%ebx
	movzwl mpqs_sievebegin,%ebx
	leal (%eax,%ebx,4),%eax
	movl %eax,4(%esp)

	movl $mpqs_FB_start,%eax
        leal (%eax,%ebx,4),%eax
	movl %eax,8(%esp)

	movl $mpqs_FB_log,%eax
#	movl $0,%ebx
	movzwl mpqs_sievebegin,%ebx
	addl %ebx,%eax
	movl %eax,12(%esp)

	movl mpqs_sievelen,%eax
	shrl $2,%eax
	movl %eax,16(%esp)

#	xorl %edx,%edx
#	xorl %ebx,%ebx
#	xorl %ecx,%ecx
	xorl %eax,%eax

mainloopa:
	movl 4(%esp),%esi
	movzwl (%esi),%ecx
	leal 4(%esi),%esi
	movl %esi,4(%esp)    # p 
	cmpl %ecx,16(%esp)
	jc enda

	movl $0,%ebx
	movl $0,%edx
	movl 8(%esp),%esi
	movzwl (%esi),%ebx
	movzwl 2(%esi),%edx
	leal 4(%esi),%esi
	movl %esi,8(%esp)    # s1, s2

	cmpl %ebx,%edx
	jnc noxch
	xchg %ebx,%edx 

noxch:                # now ebx<=edx

	movl 12(%esp),%esi
	movb (%esi),%al
	leal 1(%esi),%esi
	movl %esi,12(%esp)   # lo
	movb %al,%ah

	movl mpqs_sievelen,%esi
	movl %ecx,%edi
	leal (%edi,%ecx,2),%edi
	subl %edi,%esi
	movl mpqs_sievearray,%edi
	addl %edi,%esi

	addl %edi,%ebx
	addl %edi,%edx

loop4a:
	addb %al,(%ebx)
	addb %ah,(%edx)
	leal (%ebx,%ecx,2),%ebp
	leal (%edx,%ecx,2),%edi
	addb %al,(%ebx,%ecx)
	addb %ah,(%edx,%ecx)

	addb %al,(%ebp)
	addb %ah,(%edi)
	leal (%ebp,%ecx,2),%ebx
	leal (%edi,%ecx,2),%edx
	addb %al,(%ebp,%ecx)
	addb %ah,(%edi,%ecx)
	cmpl %esi,%edx
	jc loop4a

	leal (%esi,%ecx,2),%esi
	cmpl %esi,%edx
	jnc checka

	addb %al,(%ebx)
	addb %ah,(%edx)
	addb %al,(%ebx,%ecx)
	addb %ah,(%edx,%ecx)
	leal (%ebx,%ecx,2),%ebx
	leal (%edx,%ecx,2),%edx

checka:
	addl %ecx,%esi
	cmpl %esi,%edx
	jnc check1a
	addb %al,(%ebx)
	addb %ah,(%edx)
	leal (%ebx,%ecx),%ebx

check1a:
	cmpl %esi,%ebx
	jnc check2a
	addb %al,(%ebx)
check2a:
loopenda:
	jmp mainloopa


enda:
        addl $20,%esp
	popl %ebp
	popl %ebx
	popl %esi
	popl %edi
	ret
	.END
