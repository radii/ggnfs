define(aux_src,r32)dnl
define(lcnt,r33)dnl
define(shift_src,r34)dnl
dnl
define(primes,r2)dnl
define(aux,r3)dnl
define(rootshift,r8)dnl
define(res,r9)dnl
define(saved_cc,r10)dnl
define(saved_lc,r11)dnl
define(saved_ec,r14)dnl
define(aux_dest,r15)dnl
define(carry,r16)dnl
define(diff,r17)dnl
function_head(asm_TdUpdate)
	alloc r31=ar.pfs,3,0,0,0
	mov saved_ec=ar.ec
	mov saved_lc=ar.lc
	mov saved_cc=pr
	mov aux_dest=aux_src
	mov ar.lc=lcnt
	mov ar.ec = 2
	mov pr.rot=1<<16 ;;
MMX_TdUpdateLoop:
	(p17)pcmp2.gt carry=aux,rootshift
	(p17)psub2 diff=aux,rootshift
	(p16)ld8 aux=[aux_src],8 ;;
	(p16)ld8 rootshift=[shift_src],8
	(p18)st8 [aux_dest]=res,24
	(p17)andcm carry=primes,carry ;;
	(p16)ld8 primes=[aux_src],16
	(p17)padd2 res=diff,carry
	br.ctop.dptk MMX_TdUpdateLoop ;;
	mov ar.ec=saved_ec
	mov ar.lc=saved_lc
	mov pr=saved_cc,-1
	mov ar.pfs=r31 ;;
	br.ret.sptk b0

define(aux_ptr,r2)dnl
define(auxreg,r3)dnl
define(pbuf,r8)dnl
define(aux_ptr1,r9)dnl
define(`saved_cc',r10)dnl
define(`saved_lc',r11)dnl
define(`saved_ec',r14)dnl
define(sieve_i,r15)dnl
define(mask16,r16)dnl
function_head(asm_MMX_Td)
	alloc r31=ar.pfs,4,28,0,32
	mov aux_ptr=r34
	mov saved_cc=pr ;;
	mov mask16=1 ;;
	mov saved_lc=ar.lc
	mov saved_ec=ar.ec
	mov ar.lc=r35
	pcmp2.gt mask16=mask16,r0
	add aux_ptr1=16,aux_ptr
	mov ar.ec=6
	mov sieve_i=r33
	mov pbuf=r32 
	mov pr.rot=1<<16 ;;
asm_MMX_TdLoop:
	(p16)ld8 r32=[aux_ptr],8
	(p16)ld8 r40=[aux_ptr1],24
	(p18)pmpyshr2.u r34=r34,r42,0 ;;
	(p16)ld8 r48=[aux_ptr],16
	(p17)padd2 r33=r33,sieve_i
	(p61)br.dpnt MMX_TdStoreDivisor
asm_MMX_TdLoopEntry1:
	(p19)pmpyshr2.u r35=r35,r51,16
	(p20)pcmp2.eq r36=r36,r0
	(p21)cmp.eq.unc p6,p60=r37,r0
	br.ctop.sptk asm_MMX_TdLoop ;;
	mov ar.ec=saved_ec
	mov ar.lc=saved_lc
	mov pr=saved_cc,-1 ;;
	br.ret.sptk.many b0

MMX_TdStoreDivisor:
	add auxreg=-1,r38 ;;
	xor auxreg=auxreg,r38 ;;
	popcnt auxreg=auxreg ;;
	add auxreg=-1,auxreg ;;
	shr.u r54=r54,auxreg
	add auxreg=16,auxreg ;;
	shr.u r38=r38,auxreg
	and auxreg=mask16,r54
	shr.u r54=r54,16 ;;
	st4 [pbuf]=auxreg,4
	cmp.eq p6,p7=r38,r0 ;;
	(p7)br.dpnt MMX_TdStoreDivisor
	br.sptk asm_MMX_TdLoopEntry1
