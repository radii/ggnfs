
	section .text
	align 16
	global _asm_root_sieve

_asm_root_sieve:
	push       esi
	push       edi
	push       ebx
	push       ebp	
	mov        ebx, dword [esp+24]
	mov        ecx, dword [esp+28]
	lea        edi, [ebx+ecx*4]
	mov        ebp, dword [esp+32]
	mov        edx, dword [esp+20]
	mov        esi, dword [edx]
	mov        ecx, dword [esp+36]
	sub        ecx, dword [esp+28]
	lea        ecx, [ebp+ecx*4]
loop1:
	cmp        esi, edi
	jnc        loop2
	mov        eax, dword [esi]
	lea        esi, [esi+4]
	add        dword [ebp], eax
	lea        ebp, [ebp+4]
	jmp        loop1
	
loop2:
	mov        eax, dword [esp+28]
	and        eax, 1
	jz         loop2_even
loop2_odd:
	cmp        ebp, ecx
	jnc        loopend
	mov        esi, ebx
	mov        eax, dword [esi]
	lea        esi, [esi+4]
	add        dword [ebp], eax
	lea        ebp, [ebp+4]
innerloop_odd:
	cmp        esi, edi
	jnc        loop2_odd
	mov        eax, dword [esi]
	mov        edx, dword [esi+4]
	add        dword [ebp], eax
	lea        ebp, [ebp+4]
	lea        esi, [esi+8]
	add        dword [ebp], edx
	lea        ebp, [ebp+4]
	jmp        innerloop_odd
	
loop2_even:
	cmp        ebp, ecx
	jnc        loopend
	mov        esi, ebx
innerloop_even:
	cmp        esi, edi
	jnc        loop2_even
	mov        eax, dword [esi]
	lea        esi, [esi+4]
	add        dword [ebp], eax
	lea        ebp, [ebp+4]
	mov        eax, dword [esi]
	lea        esi, [esi+4]
	add        dword [ebp], eax
	lea        ebp, [ebp+4]
	jmp        innerloop_even
	
loopend:
	mov        esi, ebx
	mov        ebx, dword [esp+28]
	lea        ecx, [ecx+ebx*4]
loop3:
	cmp        ebp, ecx
	jnc        end1
	mov        eax, dword [esi]
	lea        esi, [esi+4]
	add        dword [ebp], eax
	lea        ebp, [ebp+4]
	jmp        loop3
	
end1:
	mov        edx, dword [esp+20]
	mov        dword [edx], esi
	
	
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret
	
	align 16
	global _asm_root_sieve8

_asm_root_sieve8:
	push       esi
	push       edi
	push       ebx
	push       ebp
	
	mov        ebx, dword [esp+24]
	mov        ecx, dword [esp+28]
	lea        edi, [ebx+ecx*8]
	mov        ebp, dword [esp+32]
	mov        edx, dword [esp+20]
	mov        esi, dword [edx]
	mov        ecx, dword [esp+36]
	sub        ecx, dword [esp+28]
	lea        ecx, [ebp+ecx*8]
loop1_8:
	cmp        esi, edi
	jnc        loop2_8
	movq       mm0, [esi]
	movq       mm1, [ebp]
	lea        esi, [esi+8]
	paddusw    mm1, mm0
	movq       [ebp], mm1
	lea        ebp, [ebp+8]
	jmp        loop1_8
	
loop2_8:
	mov        eax, dword [esp+28]
	and        eax, 1
	jz         loop2_even_8
loop2_odd_8:
	cmp        ebp, ecx
	jnc        loopend_8
	mov        esi, ebx
	movq       mm0, [esi]
	movq       mm1, [ebp]
	lea        esi, [esi+8]
	paddw      mm1, mm0
	movq       [ebp], mm1
	lea        ebp, [ebp+8]
innerloop_odd_8:
	cmp        esi, edi
	jnc        loop2_odd_8
	movq       mm0, [esi]
	movq       mm1, [ebp]
	lea        esi, [esi+8]
	paddw      mm1, mm0
	movq       mm2, [esi]
	movq       [ebp], mm1
	lea        ebp, [ebp+8]
	movq       mm3, [ebp]
	lea        esi, [esi+8]
	paddw      mm3, mm2
	movq       [ebp], mm3
	lea        ebp, [ebp+8]
	jmp        innerloop_odd_8
	
loop2_even_8:
	cmp        ebp, ecx
	jnc        loopend_8
	mov        esi, ebx
innerloop_even_8:
	cmp        esi, edi
	jnc        loop2_even_8
	movq       mm0, [esi]
	movq       mm1, [ebp]
	lea        esi, [esi+8]
	paddw      mm1, mm0
	movq       [ebp], mm1
	lea        ebp, [ebp+8]
	movq       mm0, [esi]
	movq       mm1, [ebp]
	lea        esi, [esi+8]
	paddw      mm1, mm0
	movq       [ebp], mm1
	lea        ebp, [ebp+8]
	jmp        innerloop_even_8
	
loopend_8:
	mov        esi, ebx
	mov        ebx, dword [esp+28]
	lea        ecx, [ecx+ebx*8]
loop3_8:
	cmp        ebp, ecx
	jnc        end_8
	movq       mm0, [esi]
	movq       mm1, [ebp]
	lea        esi, [esi+8]
	paddw      mm1, mm0
	movq       [ebp], mm1
	lea        ebp, [ebp+8]
	jmp        loop3_8
	
end_8:
	mov        edx, dword [esp+20]
	mov        dword [edx], esi
	emms
	
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret

%if 0
	
lohnt:
	
	align 16
	global _asm_root_eval8

_asm_root_eval8:
	push       esi
	push       edi
	push       ebx
	push       ebp
	emms
	pop        ebp
	pop        ebx
	pop        edi
	pop        esi
	ret

%endif
          
	end
