define(ri_ptr,r8)dnl
define(ij_ptr,r9)dnl
define(ij_ptr_ub,r10)dnl
define(ij_ub,r11)dnl
define(sched,r14)dnl
define(fbi_offs,r15)dnl
define(next_ij_ptr,r16)dnl
define(sched_ptr,r17)dnl
define(la,r18)dnl
define(lb,r19)dnl
define(li,r20)dnl
define(n_i_reg,r21)dnl
define(n_i_mask_reg,r22)dnl
define(reg65536,r23)dnl
define(l1_mask_reg,r24)dnl
define(advanced_ij_ptr,r25)
define(advanced_ri_ptr,r26)
define(sched_entry,r27)
define(ij,r33)dnl
define(next_ij,r32)dnl
define(ri0,r37)dnl
define(next_ri0,r36)dnl
define(ri1,r39)dnl
define(next_ri1,r38)dnl
function_head(medsched0)
	alloc r31=ar.pfs,5,3,0,8
	mov ij_ptr_ub=r34
	mov ij_ptr=r33;;
	mov ri_ptr=r32
	mov ij_ub=l1_size
	cmp.ltu p12,p13=ij_ptr,ij_ptr_ub;;
	mov sched_ptr=r35
	ld8 sched=[r35]
	shl fbi_offs=r36,16
	add next_ij_ptr=4,ij_ptr;;
	add advanced_ij_ptr=128,ij_ptr
	add advanced_ri_ptr=128,ri_ptr
	(p12)ld4 ri0=[ri_ptr],4
	(p12)ld4 ij=[ij_ptr]
	(p13)br.ret.spnt b0;;
	ld4 ri1=[ri_ptr],4
	sub next_ij=ij_ptr_ub,ij_ptr
	mov reg65536=65536;;
	mov l1_mask_reg=l1_mask
	shr.u next_ij=next_ij,2;;
	sub next_ij=next_ij,r0,1
	mov ar.ec=r0;;
	mov ar.lc=next_ij
	mov n_i_reg=n_i
	mov n_i_mask_reg=n_i_mask;;
medsched0_fbi_loop:
	lfetch.nt1 [advanced_ij_ptr],4
	cmp.ltu p6,p7=ij,ij_ub
	and la=n_i_mask_reg,ri0;;
	ld4 next_ij=[next_ij_ptr],4
	sub la=n_i_reg,la
	and lb=n_i_mask_reg,ri1;;
	ld4 next_ri0=[ri_ptr],4
	sub lb=n_i_reg,lb;;
	ld4 next_ri1=[ri_ptr],4
	and li=ij,n_i_mask_reg;;
	lfetch.nt1 [advanced_ri_ptr],8
	(p7)br.dpnt medsched0_ij_loop_ende
medsched0_ij_loop:
	cmp.ltu p8,p9=li,lb
	cmp.geu p10,p11=li,la;;
	or sched_entry=ij,fbi_offs
	(p8)add ij=ij,ri0;;
	(p10)add ij=ij,ri1
	st4 [sched]=sched_entry,4;;
	cmp.ltu p6,p7=ij,ij_ub
	and li=ij,n_i_mask_reg;;
	(p6)br.dptk medsched0_ij_loop
medsched0_ij_loop_ende:
	sub ij=ij,ij_ub
	add fbi_offs=fbi_offs,reg65536;;
	st4 [ij_ptr]=ij,4
	br.ctop.sptk medsched0_fbi_loop;;
	add ri_ptr=-8,ri_ptr
	st8 [sched_ptr]=sched
	mov ar.pfs=r31
	br.ret.sptk b0
