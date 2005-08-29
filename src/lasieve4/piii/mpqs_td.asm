
	extern _mpqs_FB
	extern _mpqs_FBk_inv
	extern _mpqs_FBk
	extern _mpqs_nFBk
	extern _mpqs_FB_A_inv
	extern _mpqs_Adiv_all
	extern _mpqs_nFB
	extern _mpqs_nAdiv_total
	extern _mpqs_256_inv_table	
	extern	_mpqs_FB_inv

%if 0

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
	global _asm_td

_asm_td:
	push       edi
	push       esi
	push       ebx
	push       ebp
	sub        esp, 20	
	mov        ebp, dword [esp+40]
	mov        eax, dword [ebp]
	and        eax, 00000ffffh
	mov        ebx, eax
	shl        eax, 16
	or         eax, ebx
	movd       mm3, eax
	psllq      mm3, 32
	movd       mm5, eax
	paddd      mm3, mm5
	pxor       mm5, mm5
	movzx      edx, word [ebp+8]
	movzx      ecx, word [_mpqs_td_begin]
	mov        esi, dword _mpqs_FB_inv_info
	mov        edi, dword _mpqs_FB_start
	movd       mm0, [esi+4]
	movd       mm1, [esi+12]
	movd       mm2, [edi+4]
	movq       mm4, mm0
	psubw      mm4, mm2
	paddw      mm4, mm3
	pmullw     mm4, mm1
	pmulhw     mm4, mm0
	pcmpeqw    mm4, mm5
	movd       eax, mm4
	or         eax, 0
	jz         loop2a
	mov        ax, word [_mpqs_nFBk_1]
	inc        ax
	mov        word [ebp+edx*2+10], ax
	inc        edx
loop2a:
	lea        esi, [esi+16]
	lea        edi, [edi+8]
	movq       mm0, [esi]
	movq       mm1, [esi+8]
	movq       mm2, [edi]
loop2:
	sub        ecx, 2
	jz         prod
	movq       mm4, mm0
	psubw      mm4, mm2
	paddw      mm4, mm3
	pmullw     mm4, mm1
	movq       mm2, [edi+8]
	lea        esi, [esi+16]
	lea        edi, [edi+8]
	pmulhw     mm4, mm0
	movq       mm0, [esi]
	movq       mm1, [esi+8]
	pcmpeqw    mm4, mm5
	movd       ebx, mm4
	psrlq      mm4, 32
	movd       eax, mm4
	or         eax, ebx
	jz         loop2
	movd       eax, mm4
	or         eax, 0
	jz         testsecond
	mov        ax, word [_mpqs_nFBk_1]
	add        ax, word [_mpqs_td_begin]
	sub        ax, cx
	inc        ax
	mov        word [ebp+edx*2+10], ax
	inc        edx
testsecond:
	or         ebx, 0
	jz         loop2
	mov        ax, word [_mpqs_nFBk_1]
	add        ax, word [_mpqs_td_begin]
	sub        ax, cx
	mov        word [ebp+edx*2+10], ax
	inc        edx
	jmp        loop2
prod:
	mov        esi, dword _mpqs_FB
	movzx      ebx, word [_mpqs_nFBk_1]
	add        ebx, ebx
	add        ebx, ebx
	sub        esi, ebx
	mov        ecx, 0
	mov        dword [esp+4], edx
	xor        edx, edx
	mov        eax, 1
	mov        ebx, 0
	mov        edi, ebx
prodloop:
	add        edx, edi
	cmp        ecx, dword [esp+4]
	jnc        prodend
	movzx      ebx, word [ebp+ecx*2+10]
	mov        edi, eax
	movzx      ebx, word [esi+ebx*4]
	mov        eax, edx
	mul        ebx
	xchg       edi, eax
	inc        ecx
	mul        ebx
	jmp        prodloop
prodend:
	mov        dword [esp+12], eax
	mov        dword [esp+16], edx
	mov        ebx, dword [esp+4]
	mov        eax, dword [esp+44]
	test       eax, 1
	jz         positive
	mov        word [ebp+ebx*2+10], 0
	inc        ebx
positive:
	mov        edi, dword [esp+48]
	mov        eax, dword [edi]
	mov        edx, dword [edi+4]
posloop:
	test       eax, 000000001h
	jnz        odd
	inc        ebx
	cmp        ebx, 27
	jnc        gotonext
	mov        cx, word [_mpqs_nFBk_1]
	mov        word [ebp+ebx*2+8], cx
	shr        edx, 1
	rcr        eax, 1
	jmp        posloop
odd:
	mov        word [ebp+8], bx
	cmp        dword [esp+16], 0
	jnz        division
	cmp        edx, dword [esp+12]
	jnc        gotonext
division:
	mov        dword [esp+8], eax
	mov        ebx, dword [esp+12]
	mov        edx, ebx
	mov        ecx, 0
	and        edx, 0000000ffh
	shr        edx, 1
	mov        edi, dword _mpqs_256_inv_table
	mov        cl, byte [edi+edx]
	mov        eax, ebx
	mul        ecx
	and        eax, 00000ff00h
	mul        ecx
	sub        ecx, eax
	mov        eax, ebx
	mul        ecx
	and        eax, 0ffff0000h
	mul        ecx
	sub        ecx, eax
	mov        eax, dword [esp+8]
	mul        ecx	
	mov        edi, dword [_mpqs_FB_inv]
	movzx      ecx, word [_mpqs_nFBk_1]
	add        ecx, ecx
	add        ecx, ecx
	sub        edi, ecx
	mov        ecx, 0
	mov        edx, dword [esp+4]
	mov        dword [esp], edx
	mov        ebx, 0
	mov        dword [esp+8], eax
tdloop:
	mov        eax, dword [esp+8]
	cmp        dword [esp], 0
	jz         tdend
	mov        edx, dword [esp]
	movzx      ebx, word [ebp+edx*2+8]
	dec        dword [esp]
	movzx      ecx, word [esi+ebx*4]
	mov        edx, dword [edi+ebx*4]
	mov        dword [esp+12], edx
divloop:
	mul        edx
	mov        dword [esp+16], eax
	mul        ecx
	test       edx, edx
	jnz        tdloop
	mov        dx, word [ebp+8]
	cmp        dx, 27
	jnc        gotonext
	mov        word [ebp+edx*2+10], bx
	inc        word [ebp+8]
	mov        eax, dword [esp+16]
	mov        edx, dword [esp+12]
	mov        dword [esp+8], eax
	jmp        divloop
tdend:
	cmp        eax, 1
	jz         lend
	mov        edi, dword _mpqs_FBk_inv
	mov        esi, dword _mpqs_FBk
	xor        ecx, ecx
	movzx      ebx, word [_mpqs_nFBk]
	inc        ebx
tdloopk:
	dec        ebx
	mov        eax, dword [esp+8]
	cmp        ebx, 0
	jz         tdendk
	movzx      ecx, word [esi+ebx*2-2]
	mov        edx, dword [edi+ebx*4-4]
	mov        dword [esp+12], edx
divloopk:
	mul        edx
	mov        dword [esp+16], eax
	mul        ecx
	test       edx, edx
	jnz        tdloopk
	mov        dx, word [ebp+8]
	cmp        dx, 27
	jnc        gotonext
	mov        word [ebp+edx*2+10], bx
	inc        word [ebp+8]
	mov        eax, dword [esp+16]
	mov        edx, dword [esp+12]
	mov        dword [esp+8], eax
	jmp        divloopk
tdendk:
	cmp        eax, 1
	jz         lend
	mov        edi, dword _mpqs_FB_A_inv
	mov        ecx, 0
	mov        esi, dword _mpqs_Adiv_all
	mov        cx, word [_mpqs_nFB]
	add        cx, word [_mpqs_nFBk]
	movzx      ebx, word [_mpqs_nAdiv_total]
	inc        ebx
	mov        dword [esp], ecx
tdloopa:
	dec        ebx
	mov        eax, dword [esp+8]
	cmp        ebx, 0
	jz         tdenda
	movzx      ecx, word [esi+ebx*2-2]
	mov        edx, dword [edi+ebx*4-4]
	mov        dword [esp+12], edx
divloopa:
	mul        edx
	mov        dword [esp+16], eax
	mul        ecx
	test       edx, edx
	jnz        tdloopa
	mov        dx, word [ebp+8]
	cmp        dx, 27
	jnc        gotonext
	add        ebx, dword [esp]
	mov        word [ebp+edx*2+10], bx
	inc        word [ebp+8]
	sub        ebx, dword [esp]
	mov        eax, dword [esp+16]
	mov        edx, dword [esp+12]
	mov        dword [esp+8], eax
	jmp        divloopa
tdenda:
lend:
	mov        eax, dword [esp+8]
	mov        edi, dword [esp+52]
	mov        dword [edi], eax
	xor        eax, eax
	emms
	add        esp, 20
	pop        ebp
	pop        ebx
	pop        esi
	pop        edi
	ret
gotonext:
	mov        eax, 1
	emms
	add        esp, 20
	pop        ebp
	pop        ebx
	pop        esi
	pop        edi
	ret
dbg:
	add        dword [esp+4], 40
	mov        edi, dword [esp+52]
	mov        edx, dword [esp+8]
	mov        dword [edi], edx
	jmp        end
	          
	end
