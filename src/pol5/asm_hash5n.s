	# knapsack-functions for HASHSHIFT=52


	.comm s11len,4
	.comm s12len,4
	.comm s21len,4
	.comm s22len,4
	.comm s11l,4
	.comm s12l,4
	.comm s21l,4
	.comm s22l,4
	.comm s12l_sort,4
	.comm s22l_sort,4
	.comm s21_begin,4
	.comm hashdataptr,4
	.comm raw_bound,4
	.comm raw_cand_ptr,4
	.comm s11_begin,4
	.comm hashpart_shift,4
	.comm hash_shift,4

	.text


	# asm_hash1(hb): hashes the values s11l[i]+s12l[j], 0<=i<s11len, 0<=j<s12len
	#                for which (s11l[i]+s12l[j])>>hashpart_shift=hb
	#
	#  hsub=hb<<hashpart_shift;
	#  memset(hashdata,0,NHASH*sizeof(uchar));
	#  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
	#  for (i=0; i<s11len; i++) {
	#    add=s11l[i];
	#    j=s11_begin[i];
	#    while (1) {
	#      h=add+s12l_sort[j];
	#      hp=h>>hashpart_shift;
	#      if (hp!=hb) break;
	#      ind=(h-hsub)>>hash_shift;
	#      if (sort[ind]>=32) return 1;
	#      hash[sort[ind]*NHASH+ind]=h;
	#      sort[ind]++;
	#      j++;
	#    }
	#    if (j>=s12len) j-=s12len;
	#    s11_begin[i]=j;
	#  }
	#  return 0;


	.align	4
	.globl asm_hash1
	.globl _asm_hash1
asm_hash1:
_asm_hash1:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp
	subl $20,%esp
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
	movl s11_begin,%esi
	movl %esi,8(%esp)
	movl hashpart_shift,%ecx
	movl 40(%esp),%eax
	shll %cl,%eax
	movl %eax,16(%esp)         # hb*2^hashpart_shift
	movl $1,%eax
	shll %cl,%eax
	negl %eax
	movl %eax,12(%esp)         # 2^32-2^hashpart_shift
	movl hash_shift,%ecx

outerloop1:
	movl (%esp),%esi
	cmpl 4(%esp),%esi
	jnc hash1_end
	movl (%esi),%edx         # s11l[i]
	leal 4(%esi),%esi
	movl %esi,(%esp)

	movl 8(%esp),%esi
	movl (%esi),%ebx         # j=s11_begin[i]
	leal 4(%esi),%esi
	movl %esi,8(%esp)

	movl s12l_sort,%edi
	leal (%edi,%ebx,4),%edi    # *s12l_sort+j
	#.align 16
innerloop1:
	movl (%edi),%ebx      # s12l_sort[j]
	addl %edx,%ebx        # h=s11l[i]+s12l_sort[j]
	movl %ebx,%eax
	subl 16(%esp),%ebx
	movl %ebx,%esi
	andl 12(%esp),%esi
	jnz innerloop1_end
	shrl %cl,%ebx      # ind=(h-hsub)>>hash_shift
	leal 4(%edi),%edi
	movzbl (%ebp,%ebx),%esi
	cmpl $32,%esi
	jnc return1_1          # too many entries for this hash value, return 1
	shll $12,%esi
	incb (%ebp,%ebx)
	addl %ebx,%esi
	movl %eax,4096(%ebp,%esi,4)
	jmp innerloop1

innerloop1_end:
	movl s12l_sort,%edx
	subl %edx,%edi
	shrl $2,%edi             # j
	movl s12len,%ebx
	xorl %eax,%eax
	subl %ebx,%edi
	cmovcl %ebx,%eax
	addl %eax,%edi           # if (j>=s12len) j-=s12len

	movl 8(%esp),%esi
	movl %edi,-4(%esi)       # s11_begin[i]=j
	jmp outerloop1

hash1_end:
	xorl %eax,%eax
	addl $20,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret
return1_1:
	movl $1,%eax
	addl $20,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret


	# asm_hash2(hb):
	#
	#  hsub=hb<<hashpart_shift;
	#  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
	#  for (i=0; i<s21len; i++) {
	#    add=s21l[i];
	#    j=s21_begin[i];
	#    while (1) {
	#      h=add+s22l_sort[j];
	#      hp=h>>hashpart_shift;
	#      if (hp!=hb) break;
	#      ind=(h-hsub)>>hash_shift;
	#      for (k=0; k<sort[ind]; k++) {
	#        if (hash[k*NHASH+ind]-h<raw_bound) {
	#          if (j>=s22len) *raw_cand_ptr++=(i+((j-s22len)<<6));
	#          else *raw_cand_ptr++=(i+(j<<6)); /* assumes npr_in_p<25 */
	#          break;
	#        }
	#      }
	#      j++;
	#    }
	#    if (j>=s22len) j-=s22len;
	#    s21_begin[i]=j;
	#  }


	.align	4
	.globl asm_hash2
	.globl _asm_hash2
_asm_hash2:
asm_hash2:
	pushl %esi
	pushl %edi
	pushl %ebx
	pushl %ebp
	subl $24,%esp
	movl hashdataptr,%ebp

	movl s21l,%esi
	movl %esi,(%esp)
	movl s21len,%eax
	leal (%esi,%eax,4),%esi
	movl %esi,4(%esp)
	movl s21_begin,%esi
	movl %esi,8(%esp)
	movl hashpart_shift,%ecx
	movl 44(%esp),%eax
	shll %cl,%eax
	movl %eax,20(%esp)         # hb*2^hashpart_shift
	movl $1,%eax
	shll %cl,%eax
	negl %eax
	movl %eax,16(%esp)         # 2^32-2^hashpart_shift
	movl hash_shift,%ecx

outerloop2:
	movl (%esp),%esi
	cmpl 4(%esp),%esi
	jnc hash2_end
	movl (%esi),%edx        # s21l[i]
	leal 4(%esi),%esi
	movl %esi,(%esp)

	movl 8(%esp),%esi
	movl (%esi),%ebx         # j=s21_begin[i]
	leal 4(%esi),%esi
	movl %esi,8(%esp)

	movl s22l_sort,%edi
	leal (%edi,%ebx,4),%edi    # *s22l_sort+j
innerloop2:
	movl (%edi),%ebx      # s22l_sort[j]
	addl %edx,%ebx        # h=s21l[i]+s22l_sort[j]
	movl %ebx,%eax
	subl 20(%esp),%ebx
	movl %ebx,%esi
	andl 16(%esp),%esi
	jnz innerloop2_end
	shrl %cl,%ebx         # ind=(h-hsub)>>hash_shift
	leal 4(%edi),%edi
	movzbl (%ebp,%ebx),%esi   # sort[ind]
	testl %esi,%esi
	jz innerloop2

	decl %esi
	shll $12,%esi
	addl %ebx,%esi
testloop:
	movl 4096(%ebp,%esi,4),%ebx
	subl %eax,%ebx
	cmpl raw_bound,%ebx
	jc store          # jc for b-a<raw_bound, jna for b-a<=raw_bound
	subl $4096,%esi
	jnc testloop
	jmp innerloop2
store:
	# the following part is used rarely -> not optimized
	movl %eax,12(%esp)     # save this, free registers: eax, ebx, esi
	movl s22l_sort,%ebx
	movl %edi,%esi
	subl %ebx,%esi
	shrl $2,%esi             # j+1
	decl %esi                # j
	movl s22len,%ebx
	xorl %eax,%eax
	subl %ebx,%esi
	cmovcl %ebx,%eax
	addl %eax,%esi           # if (j>=s22len) j-=s22len
	shll $6,%esi
	movl (%esp),%eax
	subl s21l,%eax
	shrl $2,%eax
	decl %eax
	addl %eax,%esi           # i+(j<<6) resp. i+((j-s22len)<<6)

	movl raw_cand_ptr,%ebx
	movl %esi,(%ebx)
	leal 4(%ebx),%ebx
	movl %ebx,raw_cand_ptr   # *raw_cand_ptr++=...

	movl 12(%esp),%eax
	jmp innerloop2        # break

innerloop2_end:
	movl s22l_sort,%edx
	subl %edx,%edi
	shrl $2,%edi             # j
	movl s22len,%ebx
	xorl %eax,%eax
	subl %ebx,%edi
	cmovcl %ebx,%eax
	addl %eax,%edi           # if (j>=s22len) j-=s22len

	movl 8(%esp),%esi
	movl %edi,-4(%esi)       # s21_begin[i]=j
	jmp outerloop2

hash2_end:
	addl $24,%esp
	popl %ebp
	popl %ebx
	popl %edi
	popl %esi
	ret
