define(sieve_interval,32(%esp))dnl
define(hzs_ptr_arg,36(%esp))dnl
define(hzs_ptr_ub,40(%esp))dnl
define(srb_ptr,44(%esp))dnl
define(srb_ptr_ub,48(%esp))dnl
define(cand_array_arg,52(%esp))dnl
define(csv_arg,56(%esp))dnl
dnl
dnl Create three stack variable
dnl
define(si_ptr1,(%esp))dnl
define(si_ptr_ub,4(%esp))dnl
define(hzs_ptr_var,8(%esp))dnl
dnl
define(si_ptr,%esi)dnl
define(si_ptrw,%si)dnl
define(hzs_ptr,%edi)dnl
define(srb,%ah)dnl
define(srb1,%al)dnl
define(auxreg1,%eax)dnl
define(auxreg2_byte,%cl)dnl
define(auxreg2,%ecx)dnl
define(auxreg,%edx)dnl
define(cand_array,%ebx)dnl
define(si_index,%ebp)dnl
define(cs_steps,128)dnl
define(cs_2steps,eval(2*cs_steps))dnl
define(n_i_minus_cssteps,eval(n_i-cs_steps))dnl
function_head(lasieve_search0)
	pushl %ebx
	pushl %esi
	pushl %edi
	pushl %ebp
	subl $12,%esp
	movl cand_array_arg,cand_array
	movl sieve_interval,si_ptr
	movl hzs_ptr_arg,hzs_ptr
	movl si_ptr,auxreg
	movl srb_ptr,auxreg2
	leal cs_steps`'(auxreg),auxreg
	movl cand_array_arg,cand_array
	movl auxreg,si_ptr_ub
	movl si_ptr,si_ptr1
search_j_loop:
	pxor %mm1,%mm1
	movb (auxreg2),srb
	pxor %mm7,%mm7
	movb (hzs_ptr),srb1
	incb srb1
	leal 1(hzs_ptr),hzs_ptr
	subb srb1,srb
	jb store_many_sievereports
	movzbl srb,auxreg
	movd auxreg,%mm1
	movd auxreg,%mm7
	psllq $8,%mm7
	por %mm7,%mm1
	xorl auxreg,auxreg
	movq %mm1,%mm7
	psllq $16,%mm1
	por %mm1,%mm7
	movq %mm7,%mm1
	psllq $32,%mm7
	por %mm1,%mm7
search_i_loop:
	movq   (si_ptr),%mm0
	movq  8(si_ptr),%mm1
	movq %mm7,%mm4
	movq 16(si_ptr),%mm2
	pmaxub %mm0,%mm4
	movq 24(si_ptr),%mm3
	movq %mm2,%mm5
	pmaxub %mm1,%mm4
	pmaxub %mm3,%mm5
	leal 32(si_ptr),si_ptr
	pmaxub %mm4,%mm5
	pcmpeqb %mm7,%mm5
	pmovmskb %mm5,auxreg
	cmpl $255,auxreg
	jne search_sievereport
i_loop_entry2:
	cmpl si_ptr,si_ptr_ub
	ja search_i_loop
i_loop_ende:
	movl si_ptr_ub,auxreg
	cmpl hzs_ptr,hzs_ptr_ub
	movl srb_ptr,auxreg2
	leal n_i`'(auxreg),auxreg
	leal n_i_minus_cssteps`'(si_ptr),si_ptr
	movl auxreg,si_ptr_ub
	ja search_j_loop
cs_ret2l0:
	movl hzs_ptr_arg,hzs_ptr
	leal 1(auxreg2),auxreg2
	movl si_ptr1,si_ptr
	cmpl auxreg2,srb_ptr_ub
	movl auxreg2,srb_ptr
	leal cs_2steps`'(si_ptr),auxreg
	leal cs_steps`'(si_ptr),si_ptr
	movl auxreg,si_ptr_ub
	movl si_ptr,si_ptr1
	ja search_j_loop
lss_ende:
	movl cand_array,%eax
	subl cand_array_arg,%eax
	addl $12,%esp
	popl %ebp
	emms
	popl %edi
	shrl $1,%eax
	popl %esi
	popl %ebx
	ret

.align 4
search_sievereport:
	pmaxub %mm7,%mm0
	pmaxub %mm7,%mm1
	movl hzs_ptr,hzs_ptr_var
	pmaxub %mm7,%mm2
	pmaxub %mm7,%mm3
	pcmpeqb %mm7,%mm0
	pcmpeqb %mm7,%mm1
	pcmpeqb %mm7,%mm2
	pcmpeqb %mm7,%mm3
	pmovmskb %mm0,auxreg1
	movl csv_arg,hzs_ptr
#	pand %mm6,%mm1
	pmovmskb %mm2,auxreg2
#	pand %mm6,%mm3
	sall $16,auxreg2
	pmovmskb %mm1,auxreg
	orl auxreg2,auxreg1
	pmovmskb %mm3,auxreg2
	sall $8,auxreg
	sall $24,auxreg2
	orl auxreg,auxreg1
	xorl auxreg,auxreg
	orl auxreg2,auxreg1
	leal -32(si_ptr),si_ptr
	xorl $0xffffffff,auxreg1
loop1:	
	bsfl auxreg1,auxreg2
	addl auxreg2,auxreg
	subl sieve_interval,si_ptr
	addl auxreg,si_ptr
	shrl auxreg2_byte,auxreg1
	movw si_ptrw,(cand_array)
	addl sieve_interval,si_ptr
	leal 2(cand_array),cand_array 
	movb (si_ptr),auxreg2_byte
	subl auxreg,si_ptr
	incl auxreg
	shrl $1,auxreg1
	testl auxreg1,auxreg1
	movb auxreg2_byte,(hzs_ptr)
	leal 1(hzs_ptr),hzs_ptr
	jnz loop1
	leal 32(si_ptr),si_ptr
	movl hzs_ptr,csv_arg
#	addl sieve_interval,si_ptr
	movl hzs_ptr_var,hzs_ptr
	jmp i_loop_entry2

.align 4
store_many_sievereports:
	subl sieve_interval,si_ptr
	movl hzs_ptr,hzs_ptr_var
	movl si_ptr,auxreg
	shll $16,auxreg
	addl si_ptr,auxreg
	addl $0x00010000,auxreg
	movd auxreg,%mm0
	addl $0x00020002,auxreg
	movd auxreg,%mm1
	addl $0x00020002,auxreg
	movl $0x00080008,auxreg1
	psllq $32,%mm1
	movd auxreg,%mm2
	addl $0x00020002,auxreg
	movd auxreg1,%mm3
	por %mm1,%mm0
	addl sieve_interval,si_ptr
	movd auxreg,%mm1
	movq %mm3,%mm4
	movl csv_arg,hzs_ptr
	psllq $32,%mm1
	psllq $32,%mm4
	por %mm2,%mm1
	por %mm4,%mm3
many_sr_loop:
	movq (si_ptr),%mm2
	movq 8(si_ptr),%mm4
	movq %mm0,(cand_array)
	leal 16(si_ptr),si_ptr
	paddw %mm3,%mm0
	movq %mm1,8(cand_array)
	cmpl si_ptr,si_ptr_ub
	paddw %mm3,%mm1
	movq %mm2,(hzs_ptr)
	movq %mm0,16(cand_array)
	paddw %mm3,%mm0
	movq %mm4,8(hzs_ptr)
	leal 16(hzs_ptr),hzs_ptr
	movq %mm1,24(cand_array)
	leal 32(cand_array),cand_array
	paddw %mm3,%mm1
	ja many_sr_loop
	movl hzs_ptr,csv_arg
	movl hzs_ptr_var,hzs_ptr
	jmp i_loop_ende
