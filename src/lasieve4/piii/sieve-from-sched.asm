
	section .text
	
	global _schedsieve

_schedsieve:
	push       esi
	push       edi
	push       ebx
	mov        edx, dword [esp+28]
	mov        al, byte [esp+16]
	mov        edi, dword [esp+20]
	sub        edx, 12
	mov        esi, dword [esp+24]
	cmp        edx, esi
	movzx      ebx, word [esi]
	jbe        fat_loop_end
fat_loop:
	prefetcht1 [esi+128]	
	movzx      ecx, word [esi+4]
	add        byte [edi+ebx], al
	movzx      ebx, word [esi+8]
	add        byte [edi+ecx], al
	movzx      ecx, word [esi+12]
	lea        esi, [esi+16]
	add        byte [edi+ebx], al
	movzx      ebx, word [esi]
	add        byte [edi+ecx], al
	cmp        edx, esi
	ja         fat_loop
fat_loop_end:
	add        edx, 12
	cmp        edx, esi
	movzx      ecx, word [esi+4]
	lea        esi, [esi+4]
	jbe        schedsieve_end
	add        byte [edi+ebx], al
	cmp        edx, esi
	movzx      ebx, word [esi+4]
	lea        esi, [esi+4]
	jbe        schedsieve_end
	add        byte [edi+ecx], al
	cmp        edx, esi
	jbe        schedsieve_end
	add        byte [edi+ebx], al
schedsieve_end:
	pop        ebx
	pop        edi
	pop        esi
	ret
	          
	end
