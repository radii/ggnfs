dnl return value = carry
define(carry,r8)dnl
define(limb,r9)dnl
define(new_limb,r10)dnl
define(limb_dest,r11)dnl
define(bound,r14)dnl
define(new_limb,r15)dnl
dnl
define(flimb,f8)dnl
define(fprime,f9)dnl
define(fmodular_inverse,f10)dnl
dnl
define(prime,r32)dnl
define(modular_inverse,r33)dnl
define(limb_src,r34)dnl
define(limb_ptr_ub,r35)dnl
function_head(mpz_asm_td)
### At the Moment, limb_ptr_ub is in fact number of limbs
	alloc r31=ar.pfs,4,0,0,0
	setf.sig fmodular_inverse=modular_inverse
	cmp.ltu p8,p9=r0,limb_ptr_ub
	mov limb_dest=limb_src
	setf.sig fprime=prime
	cmp.eq p6,p7=r0,r0
	shladd limb_ptr_ub=limb_ptr_ub,3,limb_src ;;
	(p8)ld8 carry=[limb_src],8 
	(p9)mov carry=r0
	(p9)br.spnt b0 ;;
mpz_asm_td_loop:
	setf.sig flimb=carry
	cmp.ltu p8,p9=limb_src,limb_ptr_ub ;;
	(p8)ld8 limb=[limb_src],8 
	xma.l flimb=flimb,fmodular_inverse,f0 ;;
	getf.sig new_limb=flimb
	xma.hu flimb=flimb,fprime,f0 ;;
	getf.sig carry=flimb
	st8 [limb_dest]=new_limb,8 ;;
	(p7)add carry=1,carry ;;
	cmp.leu p6,p7=carry,limb
	(p8)sub carry=limb,carry
	(p8)br.dptk mpz_asm_td_loop;;
	br.ret.sptk.many b0
