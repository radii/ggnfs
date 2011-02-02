
; u32_t *function(u32_t*, u32_t*, u32_t*, u32_t,  u32_t**,    u32_t)
; GAS    rax         rdi,    rsi,    rdx,   rcx,       r8,       r9
; YASM   rax         rcx,    rdx,     r8,    r9, [rsp+40], [rsp+48] 

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi, rbx, r12, r13, r14, r15

%macro expand 1
        FRAME_PROC lasched%1, 0, reg_save_list
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx,  r8
        mov     ecx, r9d
        mov      r8, [rsp+stack_use+40]
        mov     r9d, dword [rsp+stack_use+48]

        shl     rcx, n_i_bits
        shr     r9, 16
        cmp     rdx, rsi
        lea     rdx, [rdx-4]
        mov     r14d, [rdi]
        mov     r15d, [rdi+4]
        jbe     %%5
%%1:    mov     ebx, r14d
        mov     r10d, r15d
        neg     r14d
        and     ebx, (n_i|1)
        neg     r15d
        cmp     ebx, ((%1&1)|((%1&2)*(1<<(n_i_bits-1))))
        mov     r11d, 0
        cmovne  r11d, [rdi+4]
        and     r10d, (n_i|1)
        and     r14d, n_i_mask
        cmp     r10d, (n_i^((%1&1)|((%1&2)*(1<<(n_i_bits-1)))))
        mov     eax, 0
        cmovne  eax, [rdi]
        and     r15d, n_i_mask
        add     eax, r11d
%if %1 == 1  
        add     eax, n_i
%endif  
        shr     eax, 1
%if %1 == 2  
        cmp     r15d, r14d
        jbe     %%6
%%2:  
%endif  
        cmp     ecx, eax
        mov     ebx, eax
        mov     r10d, eax
        prefetcht0 [rdi+32]
        jbe     %%4
%%3:    mov     r11d, eax
        and     rbx, n_i_mask
        shr     r10, l1_bits
        and     r11d, l1_mask
        cmp     r14d, ebx
        mov     r12d, 0
        mov     r13d, 0
        cmovbe  r12d, [rdi+4]
        or      r11d, r9d
        cmp     r15, rbx
        mov     rbx, [r8+r10*8]
        cmova   r13d, [rdi]
        add     eax, r12d
        mov     [rbx], r11d
        lea     rbx, [rbx+4]
        add     eax, r13d
        mov     [r8+r10*8], rbx
        cmp     ecx, eax
        mov     ebx, eax
        mov     r10d, eax
        ja      %%3
%%4:    prefetcht0 [rsi+32]
        sub     eax, ecx
        lea     rdi, [rdi+8]
        mov     [rsi], eax
        add     r9d, 65536
        cmp     rdx, rsi
        lea     rsi, [rsi+4]
        mov     r14d, [rdi]
        mov     r15d, [rdi+4]
        ja      %%1
%%5:    mov     rax, rdi
        jmp     %%7
                
%if %1 == 2  
%%6:    mov     ebx, [rdi]
        cmp     [rdi+4], ebx
        jb      %%2
        cmp     r15d, r14d
        mov     eax, n_i
        mov     r10d, 0
        cmove   eax, [rdi+4]
        cmove   r10d, ebx
        sub     eax, r10d
        shr     eax, 1
        jmp     %%2
%endif
%%7:    END_PROC reg_save_list
%endmacro
  
    expand 0
    expand 1
    expand 2
    expand 3

    end
