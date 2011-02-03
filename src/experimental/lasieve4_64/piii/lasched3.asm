dnl Offset of first arg rel to esp, after all regs are safed and space for
dnl the auto vars is created
dnl 2 auto ulongs,4 saved regs,return address
dnl
define(ri_stack_offs,24)dnl
define(ij_ptr_stack_offs,eval(ri_stack_offs+4))dnl
define(ij_ptr_ub_stack_offs,eval(ri_stack_offs+8))dnl
define(n1j_stack_offs,eval(ri_stack_offs+12))dnl
define(sched_ptr_stack_offs,eval(ri_stack_offs+16))dnl
define(fbi_offs_stack_offs,eval(ri_stack_offs+20))dnl
dnl
dnl Now, create aliases for the arguments
dnl
define(ri_arg,ri_stack_offs()(%esp))dnl
define(ij_ptr,ij_ptr_stack_offs()(%esp))dnl
define(ij_ptr_ub,ij_ptr_ub_stack_offs()(%esp))dnl
define(ij_ub,n1j_stack_offs()(%esp))dnl
define(sched_ptr_arg,sched_ptr_stack_offs()(%esp))dnl
define(fbi_offs_arg,fbi_offs_stack_offs()(%esp))dnl
dnl
dnl and the auto variables
dnl
define(_a,(%esp))dnl
define(_b,2(%esp))dnl
dnl
dnl and for the auxilliary regs
dnl
define(ri,%eax)dnl
define(ij,%ebx)dnl
define(aux1,%ecx)dnl
define(aux1w,%cx)dnl
define(aux2,%edx)dnl
define(aux2w,%dx)dnl
define(aux3,%esi)dnl
define(sched,%ebp)dnl
define(aux4,%edi)dnl
define(aux4w,%di)dnl
define(ot_mask,eval(n_i|1))dnl
define(ot_tester1,ot_mask)dnl
define(ot_tester2,eval(n_i^ot_tester1))dnl
function_head(lasched3)
	pushl %ebp
	pushl %edi
	pushl %esi
	pushl %ebx
	subl $4,%esp
	movl fbi_offs_arg,aux1
	movl ij_ub,aux2
	shll $n_i_bits,aux2
	shll $16,aux1
	movl aux1,fbi_offs_arg
	movl sched_ptr_arg,sched
	movl ri_arg,ri
	movl aux2,ij_ub
	movl ij_ptr,aux1
	movl ij_ptr_ub,aux2
	cmpl aux1,aux2
	movl (aux1),ij
	movl (ri),aux3
	movl 4(ri),aux2
	jbe lasched3_end
lasched3_fbi_loop:
	prefetcht1 32(ri)
### First, do the same thing as the module @<Calculate first sieving event...@>
### in generic/lasched.w
	movl aux3,aux1
	movl aux2,aux4
	andl $ot_mask,aux1
	andl $ot_mask,aux4
	xorl ij,ij
	cmpl $ot_tester1,aux1
	movl $n_i,aux1
	cmovnel aux2,ij
	cmpl $ot_tester2,aux4
	movl $0,aux4
	cmovnel aux3,aux4
	addl aux4,ij
	movl aux1,aux4
	andl $n_i_mask,aux3
	andl $n_i_mask,aux2
	subl aux3,aux1
	subl aux2,aux4
	shrl $1,ij
	movw aux1w,_a
	movw aux4w,_b
	xorl aux3,aux3
	cmpl ij,ij_ub
	movl aux4,aux2
	movl ij,aux4
	movl ij,aux1
	jbe lasched3_ijloop_end
lasched3_ij_loop:
	andl $n_i_mask,aux4
	shrl $l1_bits,aux1
#	cmpw aux4w,_b
	cmpw aux4w,aux2w
	movl ij,aux2
	leal (sched,aux1,4),aux1
	cmoval (ri),aux3
	andl $l1_mask,aux2
	orl  fbi_offs_arg,aux2
	addl aux3,ij
	movl (aux1),aux3
	addl $4,(aux1)
	xorl aux1,aux1
	cmpw aux4w,_a
	movl aux2,(aux3)
	movl ij_ub,aux2
	cmovbel 4(ri),aux1
	xorl aux3,aux3
	addl aux1,ij
	cmpl ij,aux2
	movw _b,aux2w
	movl ij,aux4
	movl ij,aux1
	ja lasched3_ij_loop
lasched3_ijloop_end:
	movl ij_ptr,aux1
	movl ij_ptr_ub,aux2
	prefetcht1 32(aux1)
	subl ij_ub,ij
	leal 8(ri),ri
	movl ij,(aux1)
	leal 4(aux1),aux1
	addl $65536,fbi_offs_arg
	cmpl aux1,aux2
	movl (aux1),ij
	movl (ri),aux3
	movl 4(ri),aux2
	movl aux1,ij_ptr
	ja lasched3_fbi_loop
lasched3_end:
### Undo 2 static ulongs
	addl $4,%esp
	popl %ebx
	popl %esi
	popl %edi
	popl %ebp
	ret
