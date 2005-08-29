
%if 0
		
	section .bss

	global _montgomery_modulo_n
	global _montgomery_inv_n

_montgomery_modulo_n	resb 4
_montgomery_inv_n		resb 4

%else
	extern	_montgomery_modulo_n
	extern	_montgomery_inv_n

%endif
	
	section .text
	
	align 16
	global _asm_zero64
	
_asm_zero64:
	mov        edx, dword [esp+4]
	xor        eax, eax
	mov        dword [edx], eax
	mov        dword [edx+4], eax
	ret
	
	align 16
	global _asm_copy64

_asm_copy64:
	mov        edx, dword [esp+4]
	mov        ecx, dword [esp+8]
	mov        eax, dword [ecx]
	mov        dword [edx], eax
	mov        eax, dword [ecx+4]
	mov        dword [edx+4], eax
	ret
	
	align 16
	global _asm_sub64

_asm_sub64:
	push       esi
	mov        edx, dword [esp+12]
	mov        esi, dword [esp+16]
	mov        ecx, dword [esp+8]
	mov        eax, dword [edx]
	sub        eax, dword [esi]
	mov        dword [ecx], eax
	mov        eax, dword [edx+4]
	sbb        eax, dword [esi+4]
	mov        dword [ecx+4], eax
	jnc        sub_end
	mov        edx, dword [_montgomery_modulo_n]
	mov        eax, dword [edx]
	add        dword [ecx], eax
	mov        eax, dword [edx+4]
	adc        dword [ecx+4], eax
sub_end:
	pop        esi
	ret
	
	align 16
	global _asm_half64

_asm_half64:
	mov        ecx, dword [esp+4]
	mov        eax, dword [ecx]
	test       eax, 1
	jnz        half_odd
	mov        eax, dword [ecx+4]
	rcr        eax, 1
	mov        dword [ecx+4], eax
	mov        eax, dword [ecx]
	rcr        eax, 1
	mov        dword [ecx], eax
	ret
	
half_odd:
	push       esi
	mov        esi, dword [_montgomery_modulo_n]
	mov        eax, dword [esi]
	add        eax, dword [ecx]
	mov        dword [ecx], eax
	mov        eax, dword [esi+4]
	adc        eax, dword [ecx+4]
	mov        dword [ecx+4], eax
	rcr        dword [ecx+4], 1
	rcr        dword [ecx], 1
	pop        esi
	ret
	
	align 16
	global _asm_sub_n64

_asm_sub_n64:
	mov        edx, dword [esp+4]
	mov        ecx, dword [esp+8]
	mov        eax, dword [ecx]
	sub        dword [edx], eax
	mov        eax, dword [ecx+4]
	sbb        dword [edx+4], eax
	ret
	
	align 16
	global _asm_diff64

_asm_diff64:
	push       esi
	push       edi
	push       ebx
	mov        edi, dword [esp+20]
	mov        esi, dword [esp+24]
	mov        ebx, dword [esp+16]
	mov        eax, dword [esi+4]
	cmp        eax, dword [edi+4]
	jc         b_smaller_a
	jnz        a_smaller_b
	mov        eax, dword [esi]
	cmp        eax, dword [edi]
	jc         b_smaller_a
	sub        eax, dword [edi]
	mov        dword [ebx], eax
	xor        eax, eax
	mov        dword [ebx+4], eax
	jmp        diff_end
a_smaller_b:
	mov        eax, dword [esi]
	sub        eax, dword [edi]
	mov        dword [ebx], eax
	mov        eax, dword [esi+4]
	sbb        eax, dword [edi+4]
	mov        dword [ebx+4], eax
	jmp        diff_end
b_smaller_a:
	mov        eax, dword [edi]
	sub        eax, dword [esi]
	mov        dword [ebx], eax
	mov        eax, dword [edi+4]
	sbb        eax, dword [esi+4]
	mov        dword [ebx+4], eax
diff_end:
	pop        ebx
	pop        edi
	pop        esi
	ret
	
	align 16
	global _asm_add64

_asm_add64:
	push       esi
	push       edi
	mov        esi, dword [esp+16]
	mov        edi, dword [esp+12]
	mov        eax, dword [esi]
	add        dword [edi], eax
	mov        eax, dword [esi+4]
	adc        dword [edi+4], eax
	mov        esi, dword [_montgomery_modulo_n]
	jc         lsub
	mov        eax, dword [esi+4]
	cmp        eax, dword [edi+4]
	jc         lsub
	jnz        add_end
	mov        eax, dword [esi]
	cmp        eax, dword [edi]
	jnc        add_end
lsub:
	mov        eax, dword [esi]
	sub        dword [edi], eax
	mov        eax, dword [esi+4]
	sbb        dword [edi+4], eax
add_end:
	pop        edi
	pop        esi
	ret
	
	align 16
	global _asm_add64_ui

_asm_add64_ui:
	push       esi
	push       edi
	mov        edi, dword [esp+12]
	mov        eax, dword [esp+16]
	add        dword [edi], eax
	adc        dword [edi+4], 0
	jnc        add_ui_end
	mov        esi, dword [_montgomery_modulo_n]
	mov        eax, dword [esi]
	sub        dword [edi], eax
	mov        eax, dword [esi+4]
	sbb        dword [edi+4], eax
add_ui_end:
	pop        edi
	pop        esi
	ret
	
	align 16
	global _asm_mulm64

_asm_mulm64:
	push       esi
	push       edi
	push       ebx
	push       ebp
	sub        esp, 32
	xor        ebx, ebx
	mov        dword [esp+4], ebx
	mov        dword [esp+8], ebx
	mov        dword [esp+12], ebx
	mov        dword [esp+16], ebx
	mov        edi, dword [esp+56]
	mov        esi, dword [esp+60]
	mov        ecx, dword [edi]
	mov        eax, dword [esi]
	mul        ecx
	mov        dword [esp+4], eax
	mov        ebx, edx
	mov        eax, dword [esi+4]
	mul        ecx
	add        eax, ebx
	mov        dword [esp+8], eax
	mov        ebx, 0
	adc        ebx, edx
	mov        ecx, dword [edi+4]
	mov        eax, dword [esi]
	mul        ecx
	add        dword [esp+8], eax
	adc        ebx, edx
	mov        dword [esp+12], ebx
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+4]
	mul        ecx
	add        dword [esp+12], eax
	adc        ebx, edx
	mov        dword [esp+16], ebx
	mov        eax, dword [_montgomery_inv_n]
	mov        ecx, dword [esp+4]
	mul        ecx
	mov        ecx, eax
	mov        esi, dword [_montgomery_modulo_n]
	mov        eax, dword [esi]
	mul        ecx
	xor        ebx, ebx
	add        dword [esp+4], eax
	adc        ebx, edx
	mov        eax, dword [esi+4]
	mul        ecx
	add        eax, ebx
	adc        edx, 0
	add        dword [esp+8], eax
	adc        dword [esp+12], edx
	adc        dword [esp+16], 0
	mov        eax, dword [_montgomery_inv_n]
	mov        ecx, dword [esp+8]
	mul        ecx
	mov        ecx, eax
	mov        eax, dword [esi]
	mul        ecx
	xor        ebx, ebx
	add        dword [esp+8], eax
	adc        ebx, edx
	mov        eax, dword [esi+4]
	mul        ecx
	add        eax, ebx
	adc        edx, 0
	add        dword [esp+12], eax
	adc        dword [esp+16], edx
	mov        edi, dword [esp+52]
	mov        eax, dword [esp+12]
	mov        dword [edi], eax
	mov        eax, dword [esp+16]
	mov        dword [edi+4], eax
	jc         subtract
	cmp        eax, dword [esi+4]
	jc         ende
	jnz        subtract
	mov        eax, dword [edi]
	cmp        eax, dword [esi]
	jc         ende
subtract:
	mov        eax, dword [esi]
	sub        dword [edi], eax
	mov        eax, dword [esi+4]
	sbb        dword [edi+4], eax
ende:
	add        esp, 32
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret
	
	align 16
	global _asm_sqm64

_asm_sqm64:
	push       esi
	push       edi
	push       ebx
	push       ebp
	sub        esp, 32
	xor        ebx, ebx
	mov        dword [esp+4], ebx
	mov        dword [esp+8], ebx
	mov        dword [esp+12], ebx
	mov        dword [esp+16], ebx
	mov        edi, dword [esp+56]
	mov        eax, dword [edi]
	mul        eax
	mov        dword [esp+4], eax
	mov        dword [esp+8], edx
	mov        eax, dword [edi+4]
	mul        eax
	mov        dword [esp+12], eax
	mov        dword [esp+16], edx
	mov        ecx, dword [edi+4]
	mov        eax, dword [edi]
	mul        ecx
	add        eax, eax
	adc        edx, edx
	mov        ebx, 0
	adc        ebx, 0
	add        dword [esp+8], eax
	adc        dword [esp+12], edx
	adc        dword [esp+16], ebx
	mov        eax, dword [_montgomery_inv_n]
	mov        ecx, dword [esp+4]
	mul        ecx
	mov        ecx, eax
	mov        esi, dword [_montgomery_modulo_n]
	mov        eax, dword [esi]
	mul        ecx
	xor        ebx, ebx
	add        dword [esp+4], eax
	adc        ebx, edx
	mov        eax, dword [esi+4]
	mul        ecx
	add        eax, ebx
	adc        edx, 0
	add        dword [esp+8], eax
	adc        dword [esp+12], edx
	adc        dword [esp+16], 0
	mov        eax, dword [_montgomery_inv_n]
	mov        ecx, dword [esp+8]
	mul        ecx
	mov        ecx, eax
	mov        eax, dword [esi]
	mul        ecx
	xor        ebx, ebx
	add        dword [esp+8], eax
	adc        ebx, edx
	mov        eax, dword [esi+4]
	mul        ecx
	add        eax, ebx
	adc        edx, 0
	add        dword [esp+12], eax
	adc        dword [esp+16], edx
	mov        edi, dword [esp+52]
	mov        eax, dword [esp+12]
	mov        dword [edi], eax
	mov        eax, dword [esp+16]
	mov        dword [edi+4], eax
	jc         subtractsq
	cmp        eax, dword [esi+4]
	jc         endesq
	jnz        subtractsq
	mov        eax, dword [edi]
	cmp        eax, dword [esi]
	jc         endesq
subtractsq:
	mov        eax, dword [esi]
	sub        dword [edi], eax
	mov        eax, dword [esi+4]
	sbb        dword [edi+4], eax
endesq:
	add        esp, 32
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret
	          
	end
	