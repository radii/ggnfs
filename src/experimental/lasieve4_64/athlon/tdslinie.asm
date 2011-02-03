dnl tdslinie(aux_ptr,aux_ptr_ub,sieve_interval,tds_buffer)
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
function_head(tdslinie)
	pushl %ebx
	pushl %esi
	pushl %edi
	pushl %ebp
	subl $4,%esp
	movl aux_ptr_arg,aux_ptr
	movl $-8,auxreg
	cmpl aux_ptr,aux_ptr_ub
	jbe tdslinie_ende
	addl auxreg,aux_ptr_ub
	movl tds_buffer_arg,tds_buffer
tdslinie_fbi_loop:
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
	leal (prime,prime,4),auxreg
	subl auxreg,sieve_ptr_ub
tdslinie_loop`'i:
	movb (sieve_ptr),sv0
	orb (sieve_ptr,prime),sv0
	leal (sieve_ptr,prime,2),sieve_ptr
	orb (sieve_ptr),sv0
	orb (sieve_ptr,prime),sv0
	testb sv0,sv0
	jz tdslinie_looptest`'i
	pushl sieve_ptr_ub
	pushl sieve_ptr
	leal (sieve_ptr,prime,2),sieve_ptr_ub
	subl prime,sieve_ptr
	subl prime,sieve_ptr
tdslinie_tloop`'i`'a:
	movzbl (sieve_ptr),auxreg
	testl auxreg,auxreg
	leal (sieve_ptr,prime),sieve_ptr
	jz tdslinie_tloop`'i`'b
	decl auxreg
	pushl tds_buffer
	leal (tds_buffer,auxreg,4),auxreg
	movl (auxreg),tds_buffer
	movl prime,(tds_buffer)
	leal 4(tds_buffer),tds_buffer
	movl tds_buffer,(auxreg)
	popl tds_buffer
tdslinie_tloop`'i`'b:
	cmpl sieve_ptr,sieve_ptr_ub
	ja tdslinie_tloop`'i`'a
	popl sieve_ptr
	popl sieve_ptr_ub
tdslinie_looptest`'i:
	cmpl sieve_ptr,sieve_ptr_ub
	leal (sieve_ptr,prime,2),sieve_ptr
	ja tdslinie_loop`'i
	leal (sieve_ptr_ub,prime,4),sieve_ptr_ub
	xorb sv0,sv0
	pushl sieve_ptr
	cmpl sieve_ptr,sieve_ptr_ub
	jbe tdslinie_lasttesta`'i
	movb (sieve_ptr),sv0
	orb (sieve_ptr,prime),sv0
	leal (sieve_ptr,prime,2),sieve_ptr
tdslinie_lasttesta`'i:
	leal (sieve_ptr_ub,prime),sieve_ptr_ub
	cmpl sieve_ptr,sieve_ptr_ub
	jbe tdslinie_lasttestb`'i
	orb (sieve_ptr),sv0
tdslinie_lasttestb`'i:
	popl sieve_ptr
	testb sv0,sv0
	jz tdslinie_next_j`'i
tdslinie_tloop`'i`'c:
	movzbl (sieve_ptr),auxreg
	testl auxreg,auxreg
	leal (sieve_ptr,prime),sieve_ptr
	jz tdslinie_tloop`'i`'d
	decl auxreg
	pushl tds_buffer
	leal (tds_buffer,auxreg,4),auxreg
	movl (auxreg),tds_buffer
	movl prime,(tds_buffer)
	leal 4(tds_buffer),tds_buffer
	movl tds_buffer,(auxreg)
	popl tds_buffer
tdslinie_tloop`'i`'d:
	cmpl sieve_ptr,sieve_ptr_ub
	ja tdslinie_tloop`'i`'c
tdslinie_next_j`'i:
')
	cmpl aux_ptr,aux_ptr_ub
	movw root,root_src
	leal 8(aux_ptr),aux_ptr
	ja tdslinie_fbi_loop
tdslinie_ende:
	addl $4,%esp
	popl %ebp
	popl %edi
	popl %esi
	popl %ebx
	ret
