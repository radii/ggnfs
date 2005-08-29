
	section .text
		
	align 256
	global _asm_evaluate

_asm_evaluate:
	push       edi
	push       esi
	push       ebx
	mov        esi, dword [esp+16]
	mov        edi, dword [esp+20]
	xor        edx, edx
	pxor       mm7, mm7
	jmp        entry32
	
loop32:
	lea        esi, [esi+32]
	movq       [esi-32], mm7
	movq       [esi-24], mm7
	movq       [esi-16], mm7
	movq       [esi-8], mm7
entry32:
	cmp        esi, edi
	jz         lend
	movq       mm0, [esi]
	movq       mm1, [esi+8]
	movq       mm2, [esi+16]
	movq       mm3, [esi+24]
	por        mm1, mm0
	por        mm3, mm2
	por        mm3, mm1
	movd       ebx, mm3
	psrlq      mm3, 32
	movd       ecx, mm3
	or         ecx, ebx
	and        ecx, 080808080h
	jz         loop32
	mov        ecx, 32
loop1:
	mov        al, byte [esi]
	mov        byte [esi], 0
	lea        esi, [esi+1]
	and        al, 128
	jnz        found
entry1:
	dec        ecx
	jnz        loop1
	jmp        entry32
found:
	mov        eax, esi
	dec        eax
	sub        eax, dword [esp+16]
	mov        ebx, dword [esp+24]
	mov        dword [ebx+edx*2], eax
	inc        edx
	cmp        dword [esp+28], edx
	ja         entry1
	jmp        loop0entry
loop0:
	mov        byte [esi], 0
	lea        esi, [esi+1]
loop0entry:
	dec        ecx
	jnz        loop0
loop032:
	cmp        esi, edi
	jz         lend
	lea        esi, [esi+32]
	movq       [esi-32], mm7
	movq       [esi-24], mm7
	movq       [esi-16], mm7
	movq       [esi-8], mm7
	jmp        loop032	
lend:
	emms
	mov        eax, edx
	pop        ebx
	pop        esi
	pop        edi
	ret
	
	align 256
	global _asm_evaluate16
	
_asm_evaluate16:
	push       edi
	push       esi
	push       ebx
	mov        esi, dword [esp+16]
	mov        edi, dword [esp+20]
	xor        edx, edx
	pxor       mm7, mm7
	jmp        entry32a
loop32a:
	lea        esi, [esi+16]
	movq       [esi-16], mm7
	movq       [esi-8], mm7
entry32a:
	cmp        esi, edi
	jz         enda
	movq       mm0, [esi]
	movq       mm1, [esi+8]
	por        mm1, mm0
	movd       ebx, mm1
	psrlq      mm1, 32
	movd       ecx, mm1
	or         ecx, ebx
	and        ecx, 080808080h
	jz         loop32a
	mov        ecx, 16
loop1a:
	mov        al, byte [esi]
	mov        byte [esi], 0
	lea        esi, [esi+1]
	and        al, 128
	jnz        founda
entry1a:
	dec        ecx
	jnz        loop1a
	jmp        entry32a
founda:
	mov        eax, esi
	dec        eax
	sub        eax, dword [esp+16]
	mov        ebx, dword [esp+24]
	mov        dword [ebx+edx*2], eax
	inc        edx
	cmp        dword [esp+28], edx
	ja         entry1a
loop0a:
	mov        byte [esi], 0
	lea        esi, [esi+1]
	dec        ecx
	jnz        loop0a
loop032a:
	cmp        edi, esi
	jz         enda
	lea        esi, [esi+16]
	movq       [esi-16], mm7
	movq       [esi-8], mm7
	jmp        loop032a
enda:
	emms
	mov        eax, edx
	pop        ebx
	pop        esi
	pop        edi
	ret
	
	end
	          
