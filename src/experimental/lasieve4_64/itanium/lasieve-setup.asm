define(FB_ptr,r32)dnl
define(proot_ptr,r33)dnl
define(FB_ptr_ub,r34)dnl
define(A0,r35)dnl
define(A1,r36)dnl
define(B0,r37)dnl
define(B1,r38)dnl
define(ri,r39)dnl
define(prime,r40)dnl
define(pr64,r41)dnl
define(saved_pfs,r42)dnl
define(saved_b0,r43)dnl
define(saved_cc,r44)dnl
define(out1,r45)dnl
define(out2,r46)dnl
define(out3,r47)dnl
dnl
define(auxx1,r8)dnl
define(auxx2,r9)dnl
dnl
define(proot,f16)dnl
define(inv,f17)dnl
define(fprime,f18)dnl
define(aux2,f19)dnl
define(minus_prime,f20)dnl
define(A0x,f21)dnl
define(A1x,f22)dnl
define(B0x,f23)dnl
define(B1x,f24)dnl
dnl
define(err,f6)dnl
define(q1,f7)dnl
define(q2,f8)dnl
define(aux1,f9)dnl
define(aux3,f10)dnl
define(aux4,f11)dnl
function_head(asm_lasieve_setup)
	alloc saved_pfs=ar.pfs,8,5,3,0
	mov saved_cc=pr
### At the beginning, FB_pr_ub holds the FB size
	cmp.gtu p6,p7=FB_ptr_ub,r0
	mov saved_b0=b0
	mov auxx2=ar.fpsr
	mov auxx1=1
	sxt4 B0=B0
	sxt4 B1=B1
	sxt4 A0=A0
	sxt4 A1=A1
	add r12=-80,r12;;
	shl auxx1=auxx1,25;;
	xor auxx2=auxx1,auxx2;;
	mov ar.fpsr=auxx2
	sub B0=r0,B0
	sub A1=r0,A1
	setf.sig A0x=A0
	setf.sig B1x=B1
	(p7)br.spnt lasieve_setup_ende
	mov auxx1=r12
	shladd FB_ptr_ub=FB_ptr_ub,2,FB_ptr
	ld4 prime=[FB_ptr],4 ;;
	setf.sig B0x=B0
	setf.sig A1x=A1
	stf.spill [auxx1]=fprime,16
	ld4 pr64=[proot_ptr],4;;
	stf.spill [auxx1]=inv,16
	sub auxx2=r0,prime
	setf.sig fprime=prime;;
	stf.spill [auxx1]=proot,16
	fcvt.xf fprime=fprime
	setf.sig proot=pr64;;
	stf.spill [auxx1]=aux2,16
### Inv 1
	frcpa.s1 inv,p6=f1,fprime;;
	cmp.eq p4,p6=pr64,prime
	stf.spill [auxx1]=minus_prime,16
	setf.sig minus_prime=auxx2
### Inv 2
	fnma.s1 err=fprime,inv,f1;;
lasieve_setup_loop:
	cmp.ltu p5,p6=FB_ptr,FB_ptr_ub
	xma.l aux1=B0x,proot,A0x
	(p4)br.spnt proot_infinity
### Load next projective root.
	ld4 pr64=[proot_ptr],4
### Inv 3
	fma.s1 inv=inv,err,inv
### Inv 4
	fma.s1 err=err,err,f0
	xma.l aux2=B1x,proot,A1x ;;
lasieve_setup_have_aux:
	fcvt.xf aux3=aux1
### Inv 5
	fma.s1 inv=inv,err,inv
	fcvt.xf aux4=aux2 ;;
### Ende Inversion
	fma.s1 q1=aux3,inv,f0
	fma.s1 q2=aux4,inv,f0 ;;
	fnma.s1 aux3=q1,fprime,aux3
	fnma.s1 aux4=q2,fprime,aux4 ;;
	fma.s1 q1=aux3,inv,q1
	fma.s1 q2=aux4,inv,q2 ;;
	fcvt.fx.s1 q1=q1
	fcvt.fx.s1 q2=q2 ;;
	xma.l aux1=q1,minus_prime,aux1
	xma.l aux2=q2,minus_prime,aux2;;
	getf.sig out1=aux1 ;;
	cmp.lt p8,p9=out1,r0
	cmp.eq p3,p4=0,out1 ;;
	(p5)setf.sig proot=pr64
	(p8)add out1=out1,prime
	mov out2=prime
	(p4)br.call.sptk.many b0=asm_modinv32a#;;
	mov out2=prime
	mov out1=ri
	(p4)setf.sig aux1=r8 ;;
	(p5)ld4 prime=[FB_ptr],4
	(p4)xma.l aux1=aux1,aux2,f0;;
	(p3)mov out3=out2
	(p4)fcvt.xf aux2=aux1 ;;
	(p4)fma.s1 q2=aux2,inv,f0 ;;
	(p4)fnma.s1 aux2=q2,fprime,aux2
	(p5)setf.sig fprime=prime
	(p5)sub auxx1=r0,prime ;;
	(p4)fma.s1 q2=aux2,inv,q2 ;;
	(p5)fcvt.xf fprime=fprime
	(p4)fcvt.fx.s1 q1=q2 ;;
	(p4)xma.l aux1=q1,minus_prime,aux1 ;;
	(p4)getf.sig out3=aux1
	cmp.eq p4,p6=prime,pr64;;
### Inv 1
	(p5)frcpa.s1 inv,p7=f1,fprime
	setf.sig minus_prime=auxx1
	cmp.lt p8,p9=out3,r0;;
	(p8)add out3=out3,out2
	br.call.sptk.many b0=get_recurrence_info#
	add ri=8,ri
### Inv 2
	fnma.s1 err=fprime,inv,f1
	(p5)br.sptk lasieve_setup_loop
lasieve_setup_ende:
	mov pr=saved_cc,-1
	mov b0=saved_b0
	mov ar.pfs=saved_pfs
	mov auxx1=r12;;
	add r12=80,r12
	ldf.fill fprime=[auxx1],16 ;;
	ldf.fill inv=[auxx1],16 ;;
	ldf.fill proot=[auxx1],16 ;;
	ldf.fill aux2=[auxx1],16 ;;
	ldf.fill minus_prime=[auxx1],16 ;;
	br.ret.sptk b0

proot_infinity:
	mov aux1=B0x
	mov aux2=B1x
	### Load next projective root.
	ld4 pr64=[proot_ptr],4
### Inv 3
	fma.s1 inv=inv,err,inv
### Inv 4
	fma.s1 err=err,err,f0
	br.sptk lasieve_setup_have_aux ;;
