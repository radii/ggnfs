
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
	global _asm_zero96

_asm_zero96:
	mov        edx, dword [esp+4]
	xor        eax, eax
	mov        dword [edx], eax
	mov        dword [edx+4], eax
	mov        dword [edx+8], eax
	ret
	
	align 16
	global _asm_copy96

_asm_copy96:
	mov        edx, dword [esp+4]
	mov        ecx, dword [esp+8]
	mov        eax, dword [ecx]
	mov        dword [edx], eax
	mov        eax, dword [ecx+4]
	mov        dword [edx+4], eax
	mov        eax, dword [ecx+8]
	mov        dword [edx+8], eax
	ret
	
	align 16
	global _asm_sub_n96

_asm_sub_n96:
	mov        edx, dword [esp+4]
	mov        ecx, dword [esp+8]
	mov        eax, dword [ecx]
	sub        dword [edx], eax
	mov        eax, dword [ecx+4]
	sbb        dword [edx+4], eax
	mov        eax, dword [ecx+8]
	sbb        dword [edx+8], eax
	ret
	
	align 16
	global _asm_sub96

_asm_sub96:
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
	mov        eax, dword [edx+8]
	sbb        eax, dword [esi+8]
	mov        dword [ecx+8], eax
	jnc        sub_end
	mov        edx, dword [_montgomery_modulo_n]
	mov        eax, dword [edx]
	add        dword [ecx], eax
	mov        eax, dword [edx+4]
	adc        dword [ecx+4], eax
	mov        eax, dword [edx+8]
	adc        dword [ecx+8], eax
sub_end:
	pop        esi
	ret
	
	align 16
	global _asm_half96

_asm_half96:
	mov        ecx, dword [esp+4]
	mov        eax, dword [ecx]
	test       eax, 1
	jnz        half_odd	
	shr        dword [ecx+8], 1
	rcr        dword [ecx+4], 1
	rcr        dword [ecx], 1
	ret
half_odd:
	push       esi
	mov        esi, dword [_montgomery_modulo_n]
	mov        eax, dword [esi]
	add        dword [ecx], eax
	mov        eax, dword [esi+4]
	adc        dword [ecx+4], eax
	mov        eax, dword [esi+8]
	adc        eax, dword [ecx+8]
	rcr        eax, 1
	mov        dword [ecx+8], eax
	rcr        dword [ecx+4], 1
	rcr        dword [ecx], 1
	pop        esi
	ret

	align 16
	global _asm_diff96

_asm_diff96:
	push       esi
	push       edi
	push       ebx
	mov        edi, dword [esp+20]
	mov        esi, dword [esp+24]
	mov        ebx, dword [esp+16]
	mov        eax, dword [esi+8]
	cmp        eax, dword [edi+8]
	jc         b_smaller_a
	jnz        a_smaller_b
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
	mov        dword [ebx+8], eax
	jmp        diff_end
a_smaller_b:
	mov        eax, dword [esi]
	sub        eax, dword [edi]
	mov        dword [ebx], eax
	mov        eax, dword [esi+4]
	sbb        eax, dword [edi+4]
	mov        dword [ebx+4], eax
	mov        eax, dword [esi+8]
	sbb        eax, dword [edi+8]
	mov        dword [ebx+8], eax
	jmp        diff_end
b_smaller_a:
	mov        eax, dword [edi]
	sub        eax, dword [esi]
	mov        dword [ebx], eax
	mov        eax, dword [edi+4]
	sbb        eax, dword [esi+4]
	mov        dword [ebx+4], eax
	mov        eax, dword [edi+8]
	sbb        eax, dword [esi+8]
	mov        dword [ebx+8], eax
diff_end:
	pop        ebx
	pop        edi
	pop        esi
	ret
	
	align 16
	global _asm_add96

_asm_add96:
	push       esi
	push       edi
	mov        esi, dword [esp+16]
	mov        edi, dword [esp+12]
	mov        eax, dword [esi]
	add        dword [edi], eax
	mov        eax, dword [esi+4]
	adc        dword [edi+4], eax
	mov        eax, dword [esi+8]
	adc        dword [edi+8], eax
	mov        esi, dword [_montgomery_modulo_n]
	jc         lsub
	mov        eax, dword [esi+8]
	cmp        eax, dword [edi+8]
	jc         lsub
	jnz        add_end
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
	mov        eax, dword [esi+8]
	sbb        dword [edi+8], eax
add_end:
	pop        edi
	pop        esi
	ret
	
	align 16
	global _asm_add96_ui

_asm_add96_ui:
	push       esi
	push       edi
	mov        edi, dword [esp+12]
	mov        eax, dword [esp+16]
	add        dword [edi], eax
	adc        dword [edi+4], 0
	adc        dword [edi+8], 0
	jnc        add_ui_end
	mov        esi, dword [_montgomery_modulo_n]
	mov        eax, dword [esi]
	sub        dword [edi], eax
	mov        eax, dword [esi+4]
	sbb        dword [edi+4], eax
	mov        eax, dword [esi+8]
	sbb        dword [edi+8], eax
add_ui_end:
	pop        edi
	pop        esi
	ret
	
%if 0
	align 16
	global _asm_cmp96a

_asm_cmp96a:
	push       esi
	mov        esi, dword [esp+8]
	mov        edx, dword [esp+12]
	xor        eax, eax
	mov        ecx, dword [esi]
	sub        ecx, dword [edx]
	or         eax, ecx
	mov        ecx, dword [esi+4]
	sub        ecx, dword [edx+4]
	or         eax, ecx
	mov        ecx, dword [esi+8]
	sub        ecx, dword [edx+8]
	or         eax, ecx
	pop        esi
	ret
%endif

	align 16
	global _asm_mulm96

_asm_mulm96:
	push       esi
	push       edi
	push       ebx
	push       ebp
	sub        esp, 40	
	xor        ebx, ebx
	mov        dword [esp+4], ebx
	mov        dword [esp+8], ebx
	mov        dword [esp+12], ebx
	mov        dword [esp+16], ebx
	mov        dword [esp+20], ebx
	mov        dword [esp+24], ebx
	mov        edi, dword [esp+64]
	mov        esi, dword [esp+68]
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
	mov        eax, dword [esi+8]
	mul        ecx
	add        eax, ebx
	mov        dword [esp+12], eax
	mov        ebx, 0
	adc        ebx, edx
	mov        dword [esp+16], ebx
	mov        ebx, 0
	mov        ecx, dword [edi+4]
	mov        eax, dword [esi]
	mul        ecx
	add        dword [esp+8], eax
	adc        ebx, edx
	add        dword [esp+12], ebx
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+4]
	mul        ecx
	add        dword [esp+12], eax
	adc        ebx, edx
	add        dword [esp+16], ebx
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+8]
	mul        ecx
	add        dword [esp+16], eax
	adc        ebx, edx
	mov        dword [esp+20], ebx
	mov        ebx, 0
	mov        ecx, dword [edi+8]
	mov        eax, dword [esi]
	mul        ecx
	add        dword [esp+12], eax
	adc        ebx, edx
	add        dword [esp+16], ebx
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+4]
	mul        ecx
	add        dword [esp+16], eax
	adc        ebx, edx
	add        dword [esp+20], ebx
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+8]
	mul        ecx
	add        dword [esp+20], eax
	adc        ebx, edx
	mov        dword [esp+24], ebx
	mov        ebx, 0
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
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+8]
	mul        ecx
	add        edx, ebx
	add        dword [esp+12], eax
	adc        dword [esp+16], edx
	adc        dword [esp+20], 0
	adc        dword [esp+24], 0
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
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+8]
	mul        ecx
	add        edx, ebx
	add        dword [esp+16], eax
	adc        dword [esp+20], edx
	adc        dword [esp+24], 0
	mov        eax, dword [_montgomery_inv_n]
	mov        ecx, dword [esp+12]
	mul        ecx
	mov        ecx, eax
	mov        eax, dword [esi]
	mul        ecx
	xor        ebx, ebx
	add        dword [esp+12], eax
	adc        ebx, edx
	mov        eax, dword [esi+4]
	mul        ecx
	add        eax, ebx
	adc        edx, 0
	add        dword [esp+16], eax
	adc        dword [esp+20], edx
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+8]
	mul        ecx
	add        edx, ebx
	add        dword [esp+20], eax
	adc        dword [esp+24], edx
	mov        edi, dword [esp+60]
	mov        eax, dword [esp+16]
	mov        dword [edi], eax
	mov        eax, dword [esp+20]
	mov        dword [edi+4], eax
	mov        eax, dword [esp+24]
	mov        dword [edi+8], eax
	jc         subtract
	cmp        eax, dword [esi+8]
	jc         ende
	jnz        subtract
	mov        eax, dword [edi+4]
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
	mov        eax, dword [esi+8]
	sbb        dword [edi+8], eax
ende:
	add        esp, 40
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret
	
	align 16
	global _asm_sqm96

_asm_sqm96:
	push       esi
	push       edi
	push       ebx
	push       ebp
	sub        esp, 60
%if 0
	mov        ebp, dword [esp+84]
	mov        ecx, dword [ebp]
	mov        eax, dword [ebp+4]
	mul        ecx
	mov        esi, dword [ebp+8]
	xor        ebx, ebx
	mov        dword [esp+24], ebx
	mov        dword [esp+8], eax
	mov        dword [esp+12], edx
	mov        eax, ecx
	mul        esi
	mov        edi, dword [ebp+4]
	mov        ebx, eax
	mov        dword [esp+16], edx
	mov        eax, edi
	mul        esi
	add        dword [esp+12], ebx
	adc        dword [esp+16], 0
	mov        ebx, eax
	mov        dword [esp+20], edx
	mov        eax, ecx
	mul        ecx
	add        dword [esp+16], ebx
	adc        dword [esp+20], 0
	mov        ebx, dword [esp+8]
	mov        ecx, dword [esp+12]
	mov        dword [esp+4], eax
	mov        dword [esp+8], edx
	mov        eax, edi
	mul        edi
	add        ebx, ebx
	adc        ecx, ecx
	sbb        edi, edi
	add        dword [esp+8], ebx
	adc        ecx, 0
	mov        dword [esp+12], ecx
	sbb        ebx, ebx
	mov        dword [esp+36], eax
	mov        dword [esp+40], edx
	mov        eax, esi
	mul        esi
	add        edi, 1
	mov        ecx, dword [esp+16]
	adc        dword [esp+16], ecx
	mov        ecx, dword [esp+20]
	adc        dword [esp+20], ecx
	mov        ecx, 0
	adc        dword [esp+24], ecx
	mov        dword [esp+44], eax
	mov        dword [esp+48], edx
	add        ebx, 1
	adc        dword [esp+40], 0
	adc        dword [esp+44], 0
	adc        dword [esp+48], 0
	mov        eax, dword [esp+36]
	add        dword [esp+12], eax
	mov        eax, dword [esp+40]
	adc        dword [esp+16], eax
	mov        eax, dword [esp+44]
	adc        dword [esp+20], eax
	mov        eax, dword [esp+48]
	adc        dword [esp+24], eax
%endif

%if 0
	mov        ebp, dword [esp+84]
	mov        ecx, dword [ebp]
	mov        eax, dword [ebp+4]
	mul        ecx
	mov        esi, dword [ebp+8]
	mov        dword [esp+8], eax
	mov        dword [esp+12], edx
	mov        eax, ecx
	mul        esi
	mov        edi, dword [ebp+4]
	mov        ebx, eax
	mov        dword [esp+16], edx
	mov        eax, edi
	mul        esi
	add        dword [esp+12], ebx
	adc        dword [esp+16], 0
	mov        ebx, eax
	mov        dword [esp+20], edx
	mov        eax, ecx
	mul        ecx
	add        dword [esp+16], ebx
	adc        dword [esp+20], 0
	mov        ebx, dword [esp+8]
	mov        ecx, dword [esp+12]
	mov        dword [esp+4], eax
	mov        dword [esp+8], edx
	mov        eax, edi
	mul        edi
	add        ebx, ebx
	adc        ecx, ecx
	sbb        edi, edi
	add        dword [esp+8], ebx
	adc        ecx, 0
	mov        dword [esp+12], ecx
	sbb        ebx, ebx
	mov        dword [esp+52], eax
	mov        dword [esp+56], edx
	mov        eax, esi
	mul        esi
	add        edi, 1
	rcl        dword [esp+16], 1
	mov        ecx, 0
	rcl        dword [esp+20], 1
	adc        ecx, 0
	add        ebx, 1
	mov        ebx, dword [esp+52]
	adc        ebx, 0
	mov        edi, dword [esp+56]
	adc        edi, 0
	add        dword [esp+12], ebx
	adc        dword [esp+16], edi
	adc        eax, 0
	adc        edx, 0
	add        ecx, 1
	adc        edx, 0
	add        eax, dword [esp+20]
	adc        edx, 0
	mov        dword [esp+20], eax
	mov        dword [esp+24], edx
%endif

%if 1
	mov        ebp, dword [esp+84]
	mov        ecx, dword [ebp]
	mov        eax, dword [ebp+4]
	mul        ecx
	mov        esi, dword [ebp+8]
	xor        ebx, ebx
	mov        dword [esp+24], ebx
	mov        dword [esp+8], eax
	mov        ebx, edx
	mov        eax, ecx
	mul        esi
	xor        edi, edi
	add        eax, ebx
	mov        dword [esp+12], eax
	adc        edi, edx
	mov        eax, dword [ebp+4]
	mul        esi
	xor        ecx, ecx
	add        eax, edi
	adc        edx, 0
	shl        dword [esp+8], 1
	rcl        dword [esp+12], 1
	rcl        eax, 1
	mov        dword [esp+16], eax
	rcl        edx, 1
	mov        dword [esp+20], edx
	adc        dword [esp+24], ecx
	mov        eax, dword [ebp]
	mul        eax
	mov        dword [esp+4], eax
	mov        dword [esp+32], edx
	mov        eax, dword [ebp+4]
	mul        eax
	mov        dword [esp+36], eax
	mov        dword [esp+40], edx
	mov        eax, dword [ebp+8]
	mul        eax
	mov        ebx, dword [esp+32]
	add        dword [esp+8], ebx
	mov        ecx, dword [esp+36]
	adc        dword [esp+12], ecx
	mov        ebx, dword [esp+40]
	adc        dword [esp+16], ebx
	adc        dword [esp+20], eax
	adc        dword [esp+24], edx
%endif
	mov        eax, dword [_montgomery_inv_n]
	mov        ecx, dword [esp+4]
	mul        ecx
	mov        esi, dword [_montgomery_modulo_n]
	mov        ecx, eax
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
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+8]
	mul        ecx
	add        edx, ebx
	add        dword [esp+12], eax
	adc        dword [esp+16], edx
	adc        dword [esp+20], 0
	adc        dword [esp+24], 0
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
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+8]
	mul        ecx
	add        edx, ebx
	add        dword [esp+16], eax
	adc        dword [esp+20], edx
	adc        dword [esp+24], 0
	mov        eax, dword [_montgomery_inv_n]
	mov        ecx, dword [esp+12]
	mul        ecx
	mov        ecx, eax
	mov        eax, dword [esi]
	mul        ecx
	xor        ebx, ebx
	add        dword [esp+12], eax
	adc        ebx, edx
	mov        eax, dword [esi+4]
	mul        ecx
	add        eax, ebx
	adc        edx, 0
	add        dword [esp+16], eax
	adc        dword [esp+20], edx
	mov        ebx, 0
	adc        ebx, 0
	mov        eax, dword [esi+8]
	mul        ecx
	add        edx, ebx
	add        dword [esp+20], eax
	adc        dword [esp+24], edx
	mov        ebp, dword [esp+80]
	mov        eax, dword [esp+16]
	mov        dword [ebp], eax
	mov        eax, dword [esp+20]
	mov        dword [ebp+4], eax
	mov        eax, dword [esp+24]
	mov        dword [ebp+8], eax
	jc         subtractsq
	cmp        eax, dword [esi+8]
	jc         endesq
	jnz        subtractsq
	mov        eax, dword [ebp+4]
	cmp        eax, dword [esi+4]
	jc         endesq
	jnz        subtractsq
	mov        eax, dword [ebp]
	cmp        eax, dword [esi]
	jc         endesq
subtractsq:
	mov        eax, dword [esi]
	sub        dword [ebp], eax
	mov        eax, dword [esi+4]
	sbb        dword [ebp+4], eax
	mov        eax, dword [esi+8]
	sbb        dword [ebp+8], eax
endesq:
	add        esp, 60
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret
	          
	end
