
	extern _s12l_sort
	extern _s21_begin
	extern _s22l_sort
	
	section .bss

	extern _s11len
	extern _s12len
	extern _s21len
	extern _s22len
	extern _s11l
	extern _s12l
	extern _s21l
	extern _s22l
	extern _hashdataptr
	extern _raw_bound
	extern _raw_cand_ptr
	extern _s11_begin
	extern _hashpart_shift
	extern _hash_shift

	extern _hashdata
	
	section .text
	align 16
	global _asm_hash1

_asm_hash1:
	push       esi
	push       edi
	push       ebx
	push       ebp
	sub        esp, 20
	lea		   ebp, [_hashdata]
	
	pxor       mm0, mm0
	mov        ecx, 256
	mov        edi, ebp
zeroloop:
	movq       [edi], mm0
	movq       [edi+8], mm0
	lea        edi, [edi+16]
	dec        ecx
	jnz        zeroloop
	emms
	
	mov        esi, dword [_s11l]
	mov        dword [esp], esi
	mov        eax, dword [_s11len]
	lea        esi, [esi+eax*4]
	mov        dword [esp+4], esi
	mov        esi, dword [_s11_begin]
	mov        dword [esp+8], esi
	mov        ecx, dword [_hashpart_shift]
	mov        eax, dword [esp+40]
	shl        eax, cl
	mov        dword [esp+16], eax
	mov        eax, 1
	shl        eax, cl
	neg        eax
	mov        dword [esp+12], eax
	mov        ecx, dword [_hash_shift]
	
outerloop1:
	mov        esi, dword [esp]
	cmp        esi, dword [esp+4]
	jnc        hash1_end
	mov        edx, dword [esi]
	lea        esi, [esi+4]
	mov        dword [esp], esi
	
	mov        esi, dword [esp+8]
	mov        ebx, dword [esi]
	lea        esi, [esi+4]
	mov        dword [esp+8], esi
	
	mov        edi, dword [_s12l_sort]
	lea        edi, [edi+ebx*4]
	
innerloop1:
	mov        ebx, dword [edi]
	add        ebx, edx
	mov        eax, ebx
	sub        ebx, dword [esp+16]
	mov        esi, ebx
	and        esi, dword [esp+12]
	jnz        innerloop1_end
	shr        ebx, cl
	lea        edi, [edi+4]
	movzx      esi, byte [ebp+ebx]
	cmp        esi, 32
	jnc        return1_1          
	shl        esi, 12
	inc        byte [ebp+ebx]
	add        esi, ebx
	mov        dword [ebp+esi*4+4096], eax
	jmp        innerloop1
	
innerloop1_end:
	mov        edx, dword [_s12l_sort]
	sub        edi, edx
	shr        edi, 2
	mov        ebx, dword [_s12len]
	xor        eax, eax
	sub        edi, ebx
	cmovc      eax, ebx
	add        edi, eax
	
	mov        esi, dword [esp+8]
	mov        dword [esi-4], edi
	jmp        outerloop1
	
hash1_end:
	xor        eax, eax
	add        esp, 20
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret
return1_1:
	mov        eax, 1
	add        esp, 20
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret

	align 16
	global _asm_hash2

_asm_hash2:
	push       esi
	push       edi
	push       ebx
	push       ebp
	sub        esp, 24
	lea        ebp, [_hashdata]
	
	mov        esi, dword [_s21l]
	mov        dword [esp], esi
	mov        eax, dword [_s21len]
	lea        esi, [esi+eax*4]
	mov        dword [esp+4], esi
	mov        esi, dword [_s21_begin]
	mov        dword [esp+8], esi
	mov        ecx, dword [_hashpart_shift]
	mov        eax, dword [esp+44]
	shl        eax, cl
	mov        dword [esp+20], eax
	mov        eax, 1
	shl        eax, cl
	neg        eax
	mov        dword [esp+16], eax
	mov        ecx, dword [_hash_shift]
	
outerloop2:
	mov        esi, dword [esp]
	cmp        esi, dword [esp+4]
	jnc        hash2_end
	mov        edx, dword [esi]
	lea        esi, [esi+4]
	mov        dword [esp], esi
	
	mov        esi, dword [esp+8]
	mov        ebx, dword [esi]
	lea        esi, [esi+4]
	mov        dword [esp+8], esi
	
	mov        edi, dword [_s22l_sort]
	lea        edi, [edi+ebx*4]
innerloop2:
	mov        ebx, dword [edi]
	add        ebx, edx
	mov        eax, ebx
	sub        ebx, dword [esp+20]
	mov        esi, ebx
	and        esi, dword [esp+16]
	jnz        innerloop2_end
	shr        ebx, cl
	lea        edi, [edi+4]
	movzx      esi, byte [ebp+ebx]
	test       esi, esi
	jz         innerloop2
	
	dec        esi
	shl        esi, 12
	add        esi, ebx
testloop:
	mov        ebx, dword [ebp+esi*4+4096]
	sub        ebx, eax
	cmp        ebx, dword [_raw_bound]
	jc         store          
	sub        esi, 4096
	jnc        testloop
	jmp        innerloop2
store:
	
	mov        dword [esp+12], eax
	mov        ebx, dword [_s22l_sort]
	mov        esi, edi
	sub        esi, ebx
	shr        esi, 2
	dec        esi
	mov        ebx, dword [_s22len]
	xor        eax, eax
	sub        esi, ebx
	cmovc      eax, ebx
	add        esi, eax
	shl        esi, 6
	mov        eax, dword [esp]
	sub        eax, dword [_s21l]
	shr        eax, 2
	dec        eax
	add        esi, eax
	
	mov        ebx, dword [_raw_cand_ptr]
	mov        dword [ebx], esi
	lea        ebx, [ebx+4]
	mov        dword [_raw_cand_ptr   ], ebx
	
	mov        eax, dword [esp+12]
	jmp        innerloop2        
	
innerloop2_end:
	mov        edx, dword [_s22l_sort]
	sub        edi, edx
	shr        edi, 2
	mov        ebx, dword [_s22len]
	xor        eax, eax
	sub        edi, ebx
	cmovc      eax, ebx
	add        edi, eax
	
	mov        esi, dword [esp+8]
	mov        dword [esi-4], edi
	jmp        outerloop2
	
hash2_end:
	add        esp, 24
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret
	          
	end
