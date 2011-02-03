
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
	
;	extern _hashdata

	section .text
	
	align 16
	global _asm_hash1

_asm_hash1:
	push       esi
	push       edi
	push       ebx
	push       ebp
	sub        esp, 8
	lea		   ebp, [_hashdataptr]     ; BRG: was _hashdata
	
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
	
outerloop1:
	mov        esi, dword [esp]
	cmp        esi, dword [esp+4]
	jnc        hash1_end
	mov        edx, dword [esi]
	lea        esi, [esi+4]
	mov        dword [esp], esi
	
	mov        edi, dword [_s12l]
	mov        esi, edi
	mov        eax, dword [_s12len]
	lea        esi, [esi+eax*4]
innerloop1:
	cmp        edi, esi
	jnc        outerloop1
	mov        eax, dword [edi]
	lea        edi, [edi+4]
	add        eax, edx
	mov        ebx, eax
	shr        ebx, 20
	movzx      ecx, byte [ebp+ebx]
	cmp        ecx, 32
	jnc        return1_1          
	shl        ecx, 12
	inc        byte [ebp+ebx]
	add        ecx, ebx
	mov        dword [ebp+ecx*4+4096], eax
	
	jmp        innerloop1
	
	mov        ecx, eax
	sub        ecx, dword _raw_bound
	shr        ecx, 20
	cmp        ecx, ebx
	jz         innerloop1
	
	movzx      ebx, byte [ebp+ecx]
	cmp        ebx, 32
	jnc        return1_1          
	shl        ebx, 12
	inc        byte [ebp+ecx]
	add        ebx, ecx
	mov        dword [ebp+ebx*4+4096], eax
	jmp        innerloop1
	
hash1_end:
	xor        eax, eax
	add        esp, 8
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret
return1_1:
	mov        eax, 1
	add        esp, 8
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
	sub        esp, 8
	lea		   ebp, [_hashdataptr]     ; BRG: was _hashdata
	
	mov        esi, dword [_s21l]
	mov        dword [esp], esi
	mov        eax, dword [_s21len]
	lea        esi, [esi+eax*4]
	mov        dword [esp+4], esi
	
outerloop2:
	mov        esi, dword [esp]
	cmp        esi, dword [esp+4]
	jnc        hash2_end
	mov        edx, dword [esi]
	lea        esi, [esi+4]
	mov        dword [esp], esi
	
	mov        edi, dword [_s22l]
	mov        esi, edi
	mov        eax, dword [_s22len]
	lea        esi, [esi+eax*4]
innerloop2:
	cmp        edi, esi
	jnc        outerloop2
	mov        ebx, dword [edi]
	lea        edi, [edi+4]
	add        ebx, edx
	mov        eax, ebx
	shr        eax, 20
	movzx      ecx, byte [ebp+eax]
	test       ecx, ecx
	jz         innerloop2
	
	dec        ecx
	shl        ecx, 12
	add        ecx, eax
testloop:
	mov        eax, dword [ebp+ecx*4+4096]
	sub        eax, ebx
	cmp        eax, dword [_raw_bound]
	jnc        nostore          
	mov        ebx, edi
	sub        ebx, dword [_s22l]
	sub        ebx, 4
	shl        ebx, 14
	mov        eax, dword [esp]
	sub        eax, dword [_s21l]
	sub        eax, 4
	shr        eax, 2
	add        eax, ebx
	mov        ebx, dword [_raw_cand_ptr]
	mov        dword [ebx], eax
	lea        ebx, [ebx+4]
	mov        dword [_raw_cand_ptr], ebx
nostore:
	sub        ecx, 4096
	jnc        testloop
	jmp        innerloop2
	
hash2_end:
	add        esp, 8
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret
	          
	end
