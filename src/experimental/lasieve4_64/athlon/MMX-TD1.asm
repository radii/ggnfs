define(aux_ptr,%esi)dnl
define(prime_buffer,%edi)dnl
define(aux_ptr_ub,%ebp)dnl
define(aux0,%eax)dnl
define(aux1,%ebx)dnl
dnl %ecx is used directly as it has special relation to shift operations
define(aux2,%edx)dnl
dnl 4 saved regs + return address
define(prime_buffer_arg,20(%esp))dnl
define(strip_i_arg,24(%esp))dnl
define(aux_ptr_arg,28(%esp))dnl
define(aux_ptr_ub_arg,32(%esp))dnl
function_head(asm_MMX_Td)
	pushl %ebx
	pushl %esi
	pushl %edi
	pushl %ebp
	movl aux_ptr_arg,aux_ptr
	movl strip_i_arg,aux0
	movl strip_i_arg,aux1
	movq (aux_ptr),%mm2
	shll $16,aux1
	movq 8(aux_ptr),%mm0
	orl aux1,aux0
	movq 16(aux_ptr),%mm1
	leal 24(aux_ptr),aux_ptr
	movd aux0,%mm3
	movl aux_ptr_ub_arg,aux_ptr_ub
	psllq $32,%mm3
	movd aux0,%mm4
	leal 24(aux_ptr_ub),aux_ptr_ub
	por %mm4,%mm3
	pxor %mm4,%mm4
	paddw %mm3,%mm2
	cmpl aux_ptr,aux_ptr_ub
	movl prime_buffer_arg,prime_buffer
MMX_TdLoop:
	cmpl aux_ptr,aux_ptr_ub
	pmullw %mm2,%mm1
	movq (aux_ptr),%mm2
	jbe MMX_TdEnde
	pmulhw %mm1,%mm0
	movq 16(aux_ptr),%mm1
	pcmpeqw %mm4,%mm0
	paddw %mm3,%mm2
	pmovmskb %mm0,aux0
	movq 8(aux_ptr),%mm0
	leal 24(aux_ptr),aux_ptr
	testl aux0,aux0
	jz MMX_TdLoop
	andl $0x55,aux0
	xorl aux1,aux1
MMX_TdLoop1:
	bsfl aux0,%ecx
	shrl %cl,aux0
	shrl $1,%ecx
	addl %ecx,aux1
	movzwl -40(aux_ptr,aux1,2),aux2
	shrl $2,aux0
	movl aux2,(prime_buffer)
	incl aux1
	testl aux0,aux0
	leal 4(prime_buffer),prime_buffer
	jz MMX_TdLoop
	jmp MMX_TdLoop1
.align 4
MMX_TdEnde:
	emms
	movl prime_buffer,%eax
	popl %ebp
	popl %edi
	popl %esi
	popl %ebx
	ret

define(`auxptr_arg',4(%esp))dnl
define(`auxptr_ub_arg',8(%esp))dnl
define(`uptr_arg',12(%esp))dnl
define(`auxptr',%eax)dnl
define(`auxptr_ub',%ecx)dnl
define(`uptr',%edx)dnl
function_head(asm_TdUpdate)
	movl auxptr_arg,auxptr
	movl auxptr_ub_arg,auxptr_ub
	cmpl auxptr,auxptr_ub_arg
	movl uptr_arg,uptr
	leal -24(auxptr_ub),auxptr_ub
	movq (auxptr),%mm2
	movq (uptr),%mm0
	leal 8(uptr),uptr
	movq 8(auxptr),%mm1
	movq %mm0,%mm3
	jbe TdUpdateEnde
TdUpdateLoop:
	pcmpgtw %mm2,%mm3
	psubw %mm0,%mm2
	cmpl auxptr,auxptr_ub
	movq (uptr),%mm0
	pand %mm1,%mm3
	movq 32(auxptr),%mm1
	leal 8(uptr),uptr
	paddw %mm2,%mm3
	movq 24(auxptr),%mm2
	movq %mm3,(auxptr)
	movq %mm0,%mm3
	leal 24(auxptr),auxptr
	ja TdUpdateLoop
TdUpdateEnde:
	emms
	ret
