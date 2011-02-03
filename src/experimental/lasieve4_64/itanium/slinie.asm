define(aux_ptr,r32)dnl
define(aux_ptr_ub,r33)dnl
define(sieve_interval,r34)dnl
define(prime,r2)dnl
define(root,r3)dnl
define(proot,r8)dnl
define(sieve_base,r9)dnl
define(sieve_base2,r10)dnl
define(sieve_ptr_ub,r11)dnl
define(sieve_ptr0,r14)dnl
define(sieve_ptr1,r15)dnl
define(sieve_ptr2,r16)dnl
define(sieve_ptr3,r17)dnl
define(sv0,r18)dnl
define(sv1,r19)dnl
define(sv2,r20)dnl
define(sv3,r21)dnl
define(slog,r22)dnl
define(n_i_reg,r23)dnl
define(prime4,r24)dnl
function_head(slinie)
	alloc r31=ar.pfs,3,0,0,0
	mov n_i_reg=1 ;;
	shl n_i_reg=n_i_reg,n_i_bits ;;
slinie_fbi_loop:
	ld2 prime=[aux_ptr],2
	mov sieve_base=sieve_interval
	mov sieve_ptr_ub=sieve_interval;;
	ld2 proot=[aux_ptr],2
	sub sieve_ptr_ub=sieve_ptr_ub,prime
	shl prime4=prime,2
	shladd sieve_base2=prime,1,sieve_interval ;;
	cmp.eq p6,p7=proot,r0
	ld2 slog=[aux_ptr],2 ;;
	(p7)sub proot=prime,proot
	ld2 root=[aux_ptr],2 ;;
	sub sieve_ptr_ub=sieve_ptr_ub,prime4
	cmp.ltu p12,p13=aux_ptr,aux_ptr_ub ;;
forloop(`i',1,j_per_strip,`
	add sieve_ptr0=sieve_base,root
	add sieve_ptr2=sieve_base2,root
	cmp.ltu p6,p7=root,proot
	add sieve_ptr_ub=sieve_ptr_ub,n_i_reg
	add sieve_base=sieve_base,n_i_reg
	add sieve_base2=sieve_base2,n_i_reg
	sub root=root,proot ;;
	(p6)add root=root,prime ;;
slinie_loop`'i:
	ld1 sv0=[sieve_ptr0]
	ld1 sv2=[sieve_ptr2]
	cmp.ltu p6,p7=sieve_ptr2,sieve_ptr_ub
	add sieve_ptr1=sieve_ptr0,prime
	add sieve_ptr3=sieve_ptr2,prime ;;
	ld1 sv1=[sieve_ptr1]
	ld1 sv3=[sieve_ptr3] ;;
	add sv0=sv0,slog
	add sv2=sv2,slog ;;
	add sv1=sv1,slog
	add sv3=sv3,slog ;;
	st1 [sieve_ptr0]=sv0
	st1 [sieve_ptr2]=sv2
	add sieve_ptr0=sieve_ptr0,prime4
	st1 [sieve_ptr1]=sv1
	st1 [sieve_ptr3]=sv3
	add sieve_ptr2=sieve_ptr2,prime4
	(p6)br.dptk slinie_loop`'i ;;
	add sieve_ptr1=sieve_ptr0,prime
	cmp.ltu p6,p7=sieve_ptr0,sieve_base
	cmp.ltu p8,p9=sieve_ptr2,sieve_base ;;
	(p6)ld1 sv0=[sieve_ptr0]
	(p8)ld1 sv2=[sieve_ptr2]
	cmp.ltu p10,p11=sieve_ptr1,sieve_base ;;
	(p10)ld1 sv1=[sieve_ptr1] ;;
	(p6)add sv0=sv0,slog
	(p8)add sv2=sv2,slog ;;
	(p6)st1 [sieve_ptr0]=sv0
	(p8)st1 [sieve_ptr2]=sv2 ;;
	(p10)add sv1=sv1,slog ;;
	(p10)st1 [sieve_ptr1]=sv1')
	(p12)br.dptk slinie_fbi_loop
	br.ret.sptk b0
