	# knapsack-functions for HASHSHIFT=52


	.comm s11len,4
	.comm s12len,4
	.comm s21len,4
	.comm s22len,4
	.comm s11l,4
	.comm s12l,4
	.comm s21l,4
	.comm s22l,4
	.comm hashdataptr,4
	.comm raw_bound,4
	.comm raw_cand_ptr,4

	.text


	# asm_hash1(): hashes the values s11l[i]+s12l[j], 0<=i<s11len, 0<=j<s12len
	.align	4
	.globl asm_hash1
	.globl _asm_hash1
asm_hash1:
_asm_hash1:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp
	subl $8,%esp
	movl hashdataptr,%ebp
	# memset(hashdata,0,NHASH*sizeof(uchar)):
	pxor %mm0,%mm0
	movl $256,%ecx
	movl %ebp,%edi
zeroloop:
	movq %mm0,(%edi)
	movq %mm0,8(%edi)
	leal 16(%edi),%edi
	decl %ecx
	jnz zeroloop
	emms

	movl s11l,%esi
	movl %esi,(%esp)
	movl s11len,%eax
	leal (%esi,%eax,4),%esi
	movl %esi,4(%esp)

outerloop1:
	movl (%esp),%esi
	cmpl 4(%esp),%esi
	jnc hash1_end
	movl (%esi),%edx
	leal 4(%esi),%esi
	movl %esi,(%esp)

	movl s12l,%edi
	movl %edi,%esi
	movl s12len,%eax
	leal (%esi,%eax,4),%esi    # *s12l+s12len
innerloop1:
	cmpl %esi,%edi
	jnc outerloop1
	movl (%edi),%eax      # s12l[j]
	leal 4(%edi),%edi
	addl %edx,%eax        # s11l[i]+s12l[j]
	movl %eax,%ebx
	shrl $20,%ebx
	movzbl (%ebp,%ebx),%ecx
	cmpl $32,%ecx
	jnc return1_1          # too many entries for this hash value, return 1
	shll $12,%ecx
	incb (%ebp,%ebx)
	addl %ebx,%ecx
	movl %eax,4096(%ebp,%ecx,4)

	# we use the sloppy variant which will miss approximately
	# raw_bound out of 2^20 solutions; since raw_bound is small (often 5)
	# this is much faster
	# for using the exact variant comment out the following instruction
	jmp innerloop1

	movl %eax,%ecx
	subl $raw_bound,%ecx
	shrl $20,%ecx
	cmpl %ebx,%ecx
	jz innerloop1

	movzbl (%ebp,%ecx),%ebx
	cmpl $32,%ebx
	jnc return1_1          # too many entries for this hash value, return 1
	shll $12,%ebx
	incb (%ebp,%ecx)
	addl %ecx,%ebx
	movl %eax,4096(%ebp,%ebx,4)
	jmp innerloop1

hash1_end:
	xorl %eax,%eax
	addl $8,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret
return1_1:
	movl $1,%eax
	addl $8,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret


	# asm_hash2():
	.align	4
	.globl asm_hash2
	.globl _asm_hash2
asm_hash2:
_asm_hash2:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp
	subl $8,%esp
	movl hashdataptr,%ebp

	movl s21l,%esi
	movl %esi,(%esp)
	movl s21len,%eax
	leal (%esi,%eax,4),%esi
	movl %esi,4(%esp)

outerloop2:
	movl (%esp),%esi
	cmpl 4(%esp),%esi
	jnc hash2_end
	movl (%esi),%edx
	leal 4(%esi),%esi
	movl %esi,(%esp)

	movl s22l,%edi
	movl %edi,%esi
	movl s22len,%eax
	leal (%esi,%eax,4),%esi    # *s22l+s22len
innerloop2:
	cmpl %esi,%edi
	jnc outerloop2
	movl (%edi),%ebx      # s22l[j]
	leal 4(%edi),%edi
	addl %edx,%ebx        # s21l[i]+s22l[j]
	movl %ebx,%eax
	shrl $20,%eax
	movzbl (%ebp,%eax),%ecx
	testl %ecx,%ecx
	jz innerloop2

	decl %ecx
	shll $12,%ecx
	addl %eax,%ecx
testloop:
	movl 4096(%ebp,%ecx,4),%eax
	subl %ebx,%eax
	cmpl raw_bound,%eax
	jnc nostore          # jc for a-b<raw_bound, jna for a-b<=raw_bound
	movl %edi,%ebx
	subl s22l,%ebx
	subl $4,%ebx
	shll $14,%ebx         # j<<16
	movl (%esp),%eax
	subl s21l,%eax
	subl $4,%eax
	shrl $2,%eax          # i
	addl %ebx,%eax
	movl raw_cand_ptr,%ebx
	movl %eax,(%ebx)
	leal 4(%ebx),%ebx
	movl %ebx,raw_cand_ptr
nostore:
	subl $4096,%ecx
	jnc testloop
	jmp innerloop2

hash2_end:
	addl $8,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret
