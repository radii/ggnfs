define(prime,%ebx)dnl
define(modular_inverse,%ebp)dnl
define(limb_ptr,%esi)dnl
define(limb_ptr_ub,%edi)dnl
define(carryl,%cl)dnl
define(carry,%ecx)dnl
dnl
dnl Now, define our arguments. We have 4 saved regs and the return address.
dnl Therefore, our first argument is:
dnl
define(prime_arg,20(%esp))dnl
define(mi_arg,24(%esp))dnl
define(limb_ptr_arg,28(%esp))dnl
define(nlimbs_arg,32(%esp))dnl
dnl
function_head(mpz_asm_td)
	pushl %ebx
	pushl %esi
	pushl %edi
	pushl %ebp
	movl nlimbs_arg,%edx
	movl limb_ptr_arg,limb_ptr
	testl %edx,%edx
	movl prime_arg,prime
	movl (limb_ptr),%eax
	jz mpz_asm_td_ret
	xorl carry,carry
	cmpl $2,%edx
	leal (limb_ptr,%edx,4),limb_ptr_ub
	movl mi_arg,modular_inverse
	leal -4(limb_ptr_ub),limb_ptr_ub
	jb mpz_asm_td_last_limb
mpz_asm_td_loop:
	mull modular_inverse
	movl %eax,(limb_ptr)
	mull prime
	leal 4(limb_ptr),limb_ptr
	movl (limb_ptr),%eax
	addl carry,%edx
	subl %edx,%eax
	setb carryl
	cmpl limb_ptr,limb_ptr_ub
	ja mpz_asm_td_loop
mpz_asm_td_last_limb:
	mull modular_inverse
	movl %eax,(limb_ptr)
	mull prime
	addl carry,%edx
mpz_asm_td_ret:
	movl %edx,%eax
	popl %ebp
	popl %edi
	popl %esi
	popl %ebx
	ret
