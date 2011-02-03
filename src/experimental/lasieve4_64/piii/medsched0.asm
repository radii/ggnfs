dnl Offset of first arg rel to esp, after all regs are safed and space for
dnl the auto vars is created
dnl 4 saved regs,return address
dnl
define(ri_stack_offs,20)dnl
define(ij_ptr_stack_offs,eval(ri_stack_offs+4))dnl
define(ij_ptr_ub_stack_offs,eval(ri_stack_offs+8))dnl
define(sched_ptr_stack_offs,eval(ri_stack_offs+12))dnl
define(fbi_offs_stack_offs,eval(ri_stack_offs+16))dnl
dnl
dnl Now, create aliases for the arguments
dnl
define(ri_arg,ri_stack_offs()(%esp))dnl
define(ij_ptr,ij_ptr_stack_offs()(%esp))dnl
define(ij_ptr_ub,ij_ptr_ub_stack_offs()(%esp))dnl
define(sched_ptr_arg,sched_ptr_stack_offs()(%esp))dnl
define(fbi_offs_arg,fbi_offs_stack_offs()(%esp))dnl
dnl
dnl and for the auxilliary regs
dnl
define(ri,%eax)dnl
define(ij,%ebx)dnl
define(aux1,%ecx)dnl
define(aux2,%edx)dnl
define(_a,%esi)dnl
define(sched,%ebp)dnl
define(_b,%edi)dnl
function_head(medsched0)
	pushl %ebp
	pushl %edi
	pushl %esi
	pushl %ebx
	movl fbi_offs_arg,aux1
	movl sched_ptr_arg,aux2
	shll $16,aux1
	movl aux1,fbi_offs_arg
	movl (aux2),sched
	movl ri_arg,ri
	movl ij_ptr,aux1
	movl ij_ptr_ub,aux2
	cmpl aux1,aux2
	movl (aux1),ij
	movl (ri),aux1
	movl 4(ri),aux2
	jbe medsched0_end
medsched0_fbi_loop:
	prefetcht0 128(ri)
	movl $n_i,_a
	andl $n_i_mask,aux1
	movl _a,_b
	andl $n_i_mask,aux2
	subl aux1,_a
	subl aux2,_b
	xorl aux2,aux2
	cmpl $l1_size,ij
	movl ij,aux1
	jae medsched0_ijloop_end
medsched0_ij_loop:
	andl $n_i_mask,aux1
	orl  fbi_offs_arg,ij
	cmpl aux1,_b
	movl ij,(sched)
	leal 4(sched),sched
	cmoval (ri),aux2
	andl $l1_mask,ij
	addl aux2,ij
	xorl aux2,aux2
	cmpl aux1,_a
	cmovbel 4(ri),aux2
	addl aux2,ij
	xorl aux2,aux2
	cmpl $l1_size,ij
	movl ij,aux1
	jb medsched0_ij_loop
medsched0_ijloop_end:
	movl ij_ptr,aux1
	movl ij_ptr_ub,aux2
	prefetcht0 128(aux1)
	subl $l1_size,ij
	leal 8(ri),ri
	movl ij,(aux1)
	leal 4(aux1),aux1
	addl $65536,fbi_offs_arg
	cmpl aux1,aux2
	movl (aux1),ij
	movl aux1,ij_ptr
	movl (ri),aux1
	movl 4(ri),aux2
	ja medsched0_fbi_loop
medsched0_end:
	movl sched_ptr_arg,aux2
### Save end of sched array
	movl sched,(aux2)
	popl %ebx
	popl %esi
	popl %edi
	popl %ebp
	ret
