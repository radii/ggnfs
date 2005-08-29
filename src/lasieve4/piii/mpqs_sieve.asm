
	extern _mpqs_FB_log

%if 1

	section .bss

	global _mpqs_nFBk_1
	global _mpqs_td_begin
	global _mpqs_sievebegin
	global _mpqs_FB_inv_info
	global _mpqs_FB
	global _mpqs_FB_start
	global _mpqs_sievearray
	global _mpqs_sievelen

_mpqs_nFBk_1		resb 2
_mpqs_td_begin		resb 2
_mpqs_sievebegin	resb 2
_mpqs_FB_inv_info	resb 4
_mpqs_FB			resb 4
_mpqs_FB_start		resb 4
_mpqs_sievearray	resb 4
_mpqs_sievelen		resb 4

%else

	extern _mpqs_nFBk_1
	extern _mpqs_td_begin
	extern _mpqs_sievebegin
	extern _mpqs_FB_inv_info
	extern _mpqs_FB
	extern _mpqs_FB_start
	extern _mpqs_sievearray
	extern _mpqs_sievelen

%endif
	
	section .text
	
	align 256
	global _asm_sieve
	
_asm_sieve:
	push       edi
	push       esi
	push       ebx
	push       ebp
	sub        esp, 20
	mov        eax, dword _mpqs_FB
	movzx      ebx, word [_mpqs_sievebegin]
	add        ebx, ebx
	add        ebx, ebx
	add        eax, ebx
	mov        dword [esp+4], eax
	mov        eax, dword _mpqs_FB_start
	add        eax, ebx
	mov        dword [esp+8], eax
	mov        eax, dword _mpqs_FB_log
	movzx      ebx, word [_mpqs_sievebegin]
	add        eax, ebx
	mov        dword [esp+12], eax
	mov        eax, dword [_mpqs_sievelen]
	shr        eax, 2
	mov        dword [esp+16], eax
mainloop:
	mov        esi, dword [esp+4]
	movzx      ecx, word [esi]
	lea        esi, [esi+4]
	mov        dword [esp+4], esi
	cmp        dword [esp+16], ecx
	jc         lend
	mov        esi, dword [esp+8]
	movzx      ebx, word [esi]
	movzx      edx, word [esi+2]
	lea        esi, [esi+4]
	mov        dword [esp+8], esi
	mov        esi, dword [esp+12]
	mov        al, byte [esi]
	lea        esi, [esi+1]
	mov        dword [esp+12], esi
	mov        esi, dword [_mpqs_sievelen]
	mov        edi, ecx
	shl        edi, 2
	sub        esi, edi
	mov        edi, dword [_mpqs_sievearray]
	add        esi, edi
loop4:
	add        byte [edi+ebx], al
	add        byte [edi+edx], al
	lea        edi, [edi+ecx]
	add        byte [edi+ebx], al
	add        byte [edi+edx], al
	lea        edi, [edi+ecx]
	add        byte [edi+ebx], al
	add        byte [edi+edx], al
	lea        edi, [edi+ecx]
	add        byte [edi+ebx], al
	add        byte [edi+edx], al
	lea        edi, [edi+ecx]
	cmp        edi, esi
	jc         loop4
	add        esi, ecx
	add        esi, ecx
	cmp        edi, esi
	jnc        check
	add        byte [edi+ebx], al
	add        byte [edi+edx], al
	lea        edi, [edi+ecx]
	add        byte [edi+ebx], al
	add        byte [edi+edx], al
	lea        edi, [edi+ecx]
check:
	add        esi, ecx
	cmp        edi, esi
	jnc        check1
	add        byte [edi+ebx], al
	add        byte [edi+edx], al
	add        edi, ecx
check1:
	add        esi, ecx
	add        ebx, edi
	cmp        ebx, esi
	jnc        check2
	add        byte [ebx], al
check2:
	add        edi, edx
	cmp        edi, esi
	jnc        loopend
	add        byte [edi], al
loopend:
	jmp        mainloop
lend:
	add        esp, 20
	pop        ebp
	pop        ebx
	pop        esi
	pop        edi
	ret
	
	align 256
	global _asm_sieve1

_asm_sieve1:
	push       edi
	push       esi
	push       ebx
	push       ebp
	sub        esp, 20
	mov        eax, dword _mpqs_FB
	movzx      ebx, word [_mpqs_sievebegin]
	lea        eax, [eax+ebx*4]
	mov        dword [esp+4], eax	
	mov        eax, dword _mpqs_FB_start
	lea        eax, [eax+ebx*4]
	mov        dword [esp+8], eax
	mov        eax, dword _mpqs_FB_log
	movzx      ebx, word [_mpqs_sievebegin]
	add        eax, ebx
	mov        dword [esp+12], eax
	mov        eax, dword [_mpqs_sievelen]
	shr        eax, 2
	mov        dword [esp+16], eax
	xor        eax, eax
mainloopa:
	mov        esi, dword [esp+4]
	movzx      ecx, word [esi]
	lea        esi, [esi+4]
	mov        dword [esp+4], esi
	cmp        dword [esp+16], ecx
	jc         enda
	mov        ebx, 0
	mov        edx, 0
	mov        esi, dword [esp+8]
	movzx      ebx, word [esi]
	movzx      edx, word [esi+2]
	lea        esi, [esi+4]
	mov        dword [esp+8], esi
	cmp        edx, ebx
	jnc        noxch
	xchg       edx, ebx
noxch:
	mov        esi, dword [esp+12]
	mov        al, byte [esi]
	lea        esi, [esi+1]
	mov        dword [esp+12], esi
	mov        ah, al
	mov        esi, dword [_mpqs_sievelen]
	mov        edi, ecx
	lea        edi, [edi+ecx*2]
	sub        esi, edi
	mov        edi, dword [_mpqs_sievearray]
	add        esi, edi
	add        ebx, edi
	add        edx, edi
loop4a:
	add        byte [ebx], al
	add        byte [edx], ah
	lea        ebp, [ebx+ecx*2]
	lea        edi, [edx+ecx*2]
	add        byte [ebx+ecx], al
	add        byte [edx+ecx], ah
	add        byte [ebp], al
	add        byte [edi], ah
	lea        ebx, [ebp+ecx*2]
	lea        edx, [edi+ecx*2]
	add        byte [ebp+ecx], al
	add        byte [edi+ecx], ah
	cmp        edx, esi
	jc         loop4a
	lea        esi, [esi+ecx*2]
	cmp        edx, esi
	jnc        checka
	add        byte [ebx], al
	add        byte [edx], ah
	add        byte [ebx+ecx], al
	add        byte [edx+ecx], ah
	lea        ebx, [ebx+ecx*2]
	lea        edx, [edx+ecx*2]
checka:
	add        esi, ecx
	cmp        edx, esi
	jnc        check1a
	add        byte [ebx], al
	add        byte [edx], ah
	lea        ebx, [ebx+ecx]
check1a:
	cmp        ebx, esi
	jnc        check2a
	add        byte [ebx], al
check2a:
loopenda:
	jmp        mainloopa
enda:
	add        esp, 20
	pop        ebp
	pop        ebx
	pop        esi
	pop        edi
	ret
	          
	end
