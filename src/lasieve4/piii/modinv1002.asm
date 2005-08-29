
	extern _modulo32
	extern _Schlendrian
	
	section .data

_error_string:
	db	"Bad args to asm_modinv32\n",0
	
	section .text
	
	align 16
	global _asm_modinv32
	
_asm_modinv32:
	push       ebx
	push       esi
	push       edi
	mov        edi, dword [esp+16]
	mov        esi, dword [_modulo32]
	test       edi, edi
	jz         badargs
	cmp        esi, edi
	jbe        badargs
	xor        ebx, ebx
	xor        ecx, ecx
	inc        ecx
	cmp        edi, 1
	jbe        have_inverse2
divide:
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	sub        esi, edi
	add        ebx, ecx
	cmp        esi, edi
	jb         xlarger
	mov        eax, esi
	xor        edx, edx
	div        edi
	mov        esi, edx
	mul        ecx
	add        ebx, eax
xlarger:
	cmp        esi, 1
	jbe        have_inverse1
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	sub        edi, esi
	add        ecx, ebx
	cmp        edi, esi
	jb         ylarger
	mov        eax, edi
	xor        edx, edx
	div        esi
	mov        edi, edx
	mul        ebx
	add        ecx, eax
ylarger:
	cmp        edi, 1
	ja         divide
have_inverse2:
	jne        badargs
	mov        eax, ecx
	pop        edi
	pop        esi
	pop        ebx
	ret
have_inverse1:
	jne        badargs
	mov        eax, dword [_modulo32]
	sub        eax, ebx
	pop        edi
	pop        esi
	pop        ebx
	ret
badargs:
	push       dword _error_string
	call       _Schlendrian
	
	end

	          
