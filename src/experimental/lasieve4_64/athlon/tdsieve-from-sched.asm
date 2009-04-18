define(sched,%esi)dnl
define(sched_ub,%ebp)dnl
define(sieve_interval,%edi)dnl
define(sched_tds_buffer,%eax)dnl
define(si0,%ebx)dnl
define(si1,%ecx)dnl
define(svo,%dh)dnl
dnl 4 saved regs and return address, therefore first arg has offs 20
define(sched_arg,20(%esp))dnl
define(sched_ub_arg,24(%esp))dnl
define(sieve_interval_arg,28(%esp))dnl
define(stds_buf_arg,32(%esp))dnl
define(stds_buf_ub,36(%esp))dnl
function_head(tdsieve_sched2buf)
	pushl %ebx
	pushl %esi
	pushl %edi
	pushl %ebp
	movl sched_ub_arg,sched_ub
	movl sched_arg,si0
	leal -12(sched_ub),sched_ub
	movl (si0),sched
	movl sieve_interval_arg,sieve_interval
	cmpl sched,sched_ub
	movl stds_buf_arg,sched_tds_buffer
	jbe sched_tds_loop_ende
sched_tds_loop:
	prefetcht0 128(sched)
	movzwl (sched),si0
	movzwl 4(sched),si1
	movb (sieve_interval,si0),svo
	movzwl 8(sched),si0
	orb (sieve_interval,si1),svo
	movzwl 12(sched),si1
	orb (sieve_interval,si0),svo
	orb (sieve_interval,si1),svo
	testb svo,svo
	jnz stds_search_divisor
sched_tds_loop_entry1:
	leal 16(sched),sched
	cmpl sched,sched_ub
	ja sched_tds_loop
sched_tds_loop_ende:
	leal 8(sched_ub),sched_ub
	cmpl sched,sched_ub
	leal 4(sched_ub),sched_ub
	movzwl (sched),si0
	jbe sched_tds_last_test
	movb (sieve_interval,si0),svo
	movzwl 4(sched),si1
	testb svo,svo
	movb (sieve_interval,si1),svo
	jz sched_tds_2last_test
	movl sched,(sched_tds_buffer)
	leal 4(sched_tds_buffer),sched_tds_buffer
sched_tds_2last_test:
	testb svo,svo
	movzwl 8(sched),si0
	leal 4(sched),si1
	leal 8(sched),sched
	jz sched_tds_last_test
	movl si1,(sched_tds_buffer)
	leal 4(sched_tds_buffer),sched_tds_buffer
sched_tds_last_test:
	cmpl sched,sched_ub
	jbe sched_tds_return
	movb (sieve_interval,si0),svo
	testb svo,svo
	jz sched_tds_ende
	movl sched,(sched_tds_buffer)
	leal 4(sched_tds_buffer),sched_tds_buffer
sched_tds_ende:
	leal 4(sched),sched
sched_tds_return:
	movl sched_arg,si0
dnl Note sched_tds_buffer is in the register containing the return value
	movl sched,(si0)
	popl %ebp
	popl %edi
	popl %esi
	popl %ebx
	ret

.align 4
stds_search_divisor:
	movzwl (sched),si0
	movb (sieve_interval,si0),svo
	movzwl 4(sched),si1
	testb svo,svo
	jz stds_search2
	movl sched,(sched_tds_buffer)
	leal 4(sched_tds_buffer),sched_tds_buffer
stds_search2:
	movb (sieve_interval,si1),svo
	movzwl 8(sched),si0
	testb svo,svo
	leal 4(sched),si1
	jz tds_search3
	movl si1,(sched_tds_buffer)
	leal 4(sched_tds_buffer),sched_tds_buffer
tds_search3:
	movb (sieve_interval,si0),svo
	leal 8(sched),si0
	movzwl 12(sched),si1
	testb svo,svo
        jz tds_search4
	movl si0,(sched_tds_buffer)
	leal 4(sched_tds_buffer),sched_tds_buffer
tds_search4:
	movb (sieve_interval,si1),svo
	leal 12(sched),si1
	testb svo,svo
	jz sched_tds_test_buffer
	movl si1,(sched_tds_buffer)
        leal 4(sched_tds_buffer),sched_tds_buffer
sched_tds_test_buffer:
	cmpl sched_tds_buffer,stds_buf_ub
	ja sched_tds_loop_entry1
	leal 16(sched),sched
	jmp sched_tds_return
