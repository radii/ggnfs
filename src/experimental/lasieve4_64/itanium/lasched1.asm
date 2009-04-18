define(ri_ptr,r8)dnl
define(ij_ptr,r9)dnl
define(ij_ptr_ub,r10)dnl
define(ij_ub,r11)dnl
define(sched,r14)dnl
define(fbi_offs,r15)dnl
define(aux0,r16)dnl
define(sched_index,r17)dnl
define(sched_ptr,r18)dnl
define(sched_elem,r19)dnl
define(sched_entry,r20)dnl
define(la,r21)dnl
define(lb,r22)dnl
define(li,r23)dnl
define(n_i_reg,r24)dnl
define(n_i_mask_reg,r25)dnl
define(reg65536,r26)dnl
define(l1_mask_reg,r27)dnl
define(advanced_ij_ptr,r29)
define(advanced_ri_ptr,r30)
define(ij,r33)dnl
define(aux1,r32)dnl
define(ri0,r37)dnl
define(next_ri0,r36)dnl
define(ri1,r39)dnl
define(next_ri1,r38)dnl
dnl
define(ot_mask_reg,r2)dnl
define(ot_tester1,1)dnl
define(ot_tester2,ot_mask_reg)dnl
function_head(lasched1)
	alloc r31=ar.pfs,6,2,0,8
	mov ij_ptr_ub=r34
	mov ij_ptr=r33;;
	mov ri_ptr=r32
	shl ij_ub=r35,n_i_bits
	cmp.ltu p12,p13=ij_ptr,ij_ptr_ub;;
	mov sched=r36
	shl fbi_offs=r37,16
	add advanced_ij_ptr=128,ij_ptr
	add advanced_ri_ptr=128,ri_ptr
	(p12)ld4 ri0=[ri_ptr],4
	(p12)ld4 ij=[ij_ptr]
	(p13)br.ret.spnt b0;;
	ld4 ri1=[ri_ptr],4
	sub aux1=ij_ptr_ub,ij_ptr
	mov reg65536=65536;;
	mov l1_mask_reg=l1_mask
	shr.u aux1=aux1,2;;
	sub aux1=aux1,r0,1
	mov ar.ec=r0;;
	mov ar.lc=aux1
	mov n_i_reg=n_i;;
	or ot_mask_reg=1,n_i_reg
	mov n_i_mask_reg=n_i_mask;;
lasched1_fbi_loop:
#	lfetch.nt1 [advanced_ij_ptr],4
	mov ij=n_i
	and aux0=ri0,ot_mask_reg
	and aux1=ri1,ot_mask_reg;;
	cmp.eq p8,p9=ot_tester1,aux0
	cmp.eq p10,p11=ot_tester2,aux1
	and la=n_i_mask_reg,ri0;;
	(p9)add ij=ij,ri1
	sub la=n_i_reg,la
	and lb=n_i_mask_reg,ri1;;
	(p11)add ij=ij,ri0
	ld4 next_ri0=[ri_ptr],4;;
	shr.u ij=ij,1
	ld4 next_ri1=[ri_ptr],4
	sub lb=n_i_reg,lb;;
	cmp.ltu p6,p7=ij,ij_ub
	shr.u sched_index=ij,l1_bits
	and li=ij,n_i_mask_reg;;
	shladd sched_elem=sched_index,3,sched;;
	lfetch.nt1 [advanced_ri_ptr],8
	(p7)br.dpnt lasched1_ij_loop_ende
lasched1_ij_loop:
	ld8 sched_ptr=[sched_elem]
	cmp.ltu p8,p9=li,lb
	cmp.geu p10,p11=li,la;;
	and sched_entry=ij,l1_mask_reg
	(p8)add ij=ij,ri0;;
	or sched_entry=sched_entry,fbi_offs
	(p10)add ij=ij,ri1;;
	st4 [sched_ptr]=sched_entry,4
	shr.u sched_index=ij,l1_bits
	cmp.ltu p6,p7=ij,ij_ub;;
	st8 [sched_elem]=sched_ptr
	and li=ij,n_i_mask_reg
	shladd sched_elem=sched_index,3,sched;;
	(p6)br.dptk lasched1_ij_loop
lasched1_ij_loop_ende:
	sub ij=ij,ij_ub
	add fbi_offs=fbi_offs,reg65536;;
	st4 [ij_ptr]=ij,4
	br.ctop.sptk lasched1_fbi_loop;;
	add ri_ptr=-8,ri_ptr
	mov ar.pfs=r31
	br.ret.sptk b0
