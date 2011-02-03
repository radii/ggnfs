dnl tdslinie3(aux_ptr,aux_ptr_ub,sieve_interval,tds_buffer)
dnl We save %ebx,%esi,%edi,%ebp and also have one auto variable and the
dnl return address on the stack. Therefore, the stack offset of the
dnl first arg is 24.
define(aux_ptr_arg,24(%esp))dnl
define(aux_ptr_ub,28(%esp))dnl
define(sieve_interval,32(%esp))dnl
define(tds_buffer_arg,36(%esp))dnl
dnl Now, the registers which we are going to use
define(sieve_ptr,%edi)dnl
define(sieve_ptr_ub,%ebp)dnl
define(root,%edx)dnl
define(prime,%ecx)dnl
define(tds_buffer,%eax)dnl
define(sv0,%bh)dnl
dnl The bx-register may also be used for auxilliary 32-bit values if sv1
dnl is not used
define(auxreg,%ebx)dnl
define(aux_ptr,%esi)dnl
dnl Offset of the various things from this pointer
define(prime_src,(%esi))dnl
define(proot_src,2(%esi))dnl
define(root_src,6(%esi))dnl
dnl We store the int difference projective_root-prime here:
define(proot,(%esp))dnl
dnl This macro is taken from the GNU info documentation of m4.
define(`forloop',
       `pushdef(`$1', `$2')_forloop(`$1', `$2', `$3', `$4')popdef(`$1')')dnl
define(`_forloop',
       `$4`'ifelse($1, `$3', ,
     		 `define(`$1', incr($1))_forloop(`$1', `$2', `$3', `$4')')')dnl
function_head(tdslinie3)
	pushl %ebx
	pushl %esi
	pushl %edi
	pushl %ebp
	subl $4,%esp
	movl aux_ptr_arg,aux_ptr
	movl $-8,auxreg
	cmpl aux_ptr,aux_ptr_ub
	jbe tdslinie3_ende
	addl auxreg,aux_ptr_ub
	movl tds_buffer_arg,tds_buffer
tdslinie3_fbi_loop:
	movzwl proot_src,auxreg
	movzwl prime_src,prime
	movzwl root_src,root
	subl prime,auxreg
	movl auxreg,proot	
	movl sieve_interval,sieve_ptr_ub
forloop(`i',1,j_per_strip,`
	movl root,sieve_ptr
	xorl auxreg,auxreg
	addl proot,root
	leal (sieve_ptr_ub,sieve_ptr),sieve_ptr
	cmovncl prime,auxreg
	add $n_i,sieve_ptr_ub
	add auxreg,root
	movb (sieve_ptr),sv0
	subl prime,sieve_ptr_ub
	orb (sieve_ptr,prime),sv0
	leal (sieve_ptr,prime,2),sieve_ptr
	orb (sieve_ptr),sv0
	cmpl sieve_ptr,sieve_ptr_ub
	jbe tdslinie3_t2_`'i
	orb (sieve_ptr,prime),sv0
tdslinie3_t2_`'i:
	testb sv0,sv0
	leal (sieve_ptr_ub,prime),sieve_ptr_ub
	jz tdslinie3_next_j`'i
	subl prime,sieve_ptr
	subl prime,sieve_ptr
tdslinie3_tloop`'i:
	movzbl (sieve_ptr),auxreg
	testl auxreg,auxreg
	leal (sieve_ptr,prime),sieve_ptr
	jz tdslinie3_tdloop`'i`'a
	decl auxreg
	pushl tds_buffer
	leal (tds_buffer,auxreg,4),auxreg
	movl (auxreg),tds_buffer
	movl prime,(tds_buffer)
	leal 4(tds_buffer),tds_buffer
	movl tds_buffer,(auxreg)
	popl tds_buffer
tdslinie3_tdloop`'i`'a:
	cmpl sieve_ptr,sieve_ptr_ub
	ja tdslinie3_tloop`'i
tdslinie3_next_j`'i:
')
	cmpl aux_ptr,aux_ptr_ub
	movw root,root_src
	leal 8(aux_ptr),aux_ptr
	ja tdslinie3_fbi_loop
tdslinie3_ende:
	addl $4,%esp
	popl %ebp
	popl %edi
	popl %esi
	popl %ebx
	ret
