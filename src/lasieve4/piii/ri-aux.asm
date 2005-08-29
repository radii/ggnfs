	
	section .text
	
	align 16
	global _asm_getbc

_asm_getbc:
	push       ebx
	push       esi
	push       edi
	push       ebp
	mov        edi, dword [esp+20]
	xor        ecx, ecx
	mov        ebp, dword [esp+28]
	xor        ebx, ebx
	mov        esi, dword [esp+24]
	inc        ecx
	cmp        ebp, edi
	ja         have_bs
divide:
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         test_y
	mov        eax, esi
	xor        edx, edx
	div        edi
	mov        esi, edx
	mul        ecx
	add        ebx, eax
test_y:
	cmp        ebp, esi
	ja         have_ct
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         test_x
	mov        eax, edi
	xor        edx, edx
	div        esi
	mov        edi, edx
	mul        ebx
	add        ecx, eax
test_x:
	cmp        ebp, edi
	jbe        divide
have_bs:
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	sub        esi, edi
	add        ebx, ecx
	cmp        ebp, esi
	ja         have_bsct
	mov        eax, esi
	xor        edx, edx
	sub        eax, ebp
	div        edi
	inc        eax
	mov        ebp, eax
	mul        edi
	sub        esi, eax
	mov        eax, ebp
	mul        ecx
	add        ebx, eax
	jmp        have_bsct
have_ct:
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	sub        edi, esi
	add        ecx, ebx
	cmp        ebp, edi
	ja         have_bsct
	mov        eax, edi
	xor        edx, edx
	sub        eax, ebp
	div        esi
	inc        eax
	mov        ebp, eax
	mul        esi
	sub        edi, eax
	mov        eax, ebp
	mul        ebx
	add        ecx, eax
have_bsct:
	mov        eax, dword [esp+32]
	mov        edx, dword [esp+36]
	mov        ebp, dword [esp+40]
	mov        dword [eax], edi
	mov        eax, dword [esp+44]
	mov        dword [edx], ecx
	mov        dword [ebp], esi
	mov        dword [eax], ebx
	pop        ebp
	pop        edi
	pop        esi
	pop        ebx
	ret
	          
	end
