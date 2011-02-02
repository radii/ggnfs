
; u32_t *function(u32_t*, u32_t*, u32_t*, u32_t,  u32_t**,    u32_t)
; GAS    rax         rdi,    rsi,    rdx,   rcx,       r8,       r9
; YASM   rax         rcx,    rdx,     r8,    r9, [rsp+40], [rsp+48] 

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi, rbx, rbp, r12, r13, r14, r15

; %1 = ot %2 = nt_sched
%macro expand 2
%if %2 == 1
%define name lasched%1nt
%else
%define name  lasched%1
%endif
    	FRAME_PROC name, 1, reg_save_list
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx,  r8
        mov     ecx, r9d
        mov      r8, [rsp+stack_use+40]
        mov     r9d, dword [rsp+stack_use+48]

        shl     rcx, n_i_bits
        mov     eax, [rsi]
        cmp     rdx, rsi
        lea     rdx, [rdx-4]
        jbe     %%4
        shr     r9, 16
        mov     [rsp], rdx
%%1:    mov     r14d, [rdi]
        mov     r15d, [rdi+4]
        cmp     ecx, eax
        mov     ebx, eax
        mov     r10d, eax
        prefetch [rdi+128]
        jbe     %%3
        mov     edx, r14d
        neg     r14d
        mov     ebp, r15d
        neg     r15d
        and     r14d, n_i_mask
        and     r15d, n_i_mask
%%2:    mov     r11d, eax
        and     rbx, n_i_mask
        shr     r10, l1_bits
        and     r11d, l1_mask
        cmp     r14d, ebx
        mov     r12d, 0
        mov     r13d, 0
        cmovbe  r12d, ebp
        or      r11d, r9d
        cmp     r15, rbx
        mov     rbx, [r8+r10*8]
        cmova   r13d, edx
        add     eax, r12d
%if %2 == 1  
        movnti  [rbx], r11d
%else   
        mov     [rbx], r11d
%endif  
        lea     rbx, [rbx+4]
        add     eax, r13d
        mov     [r8+r10*8], rbx
        cmp     ecx, eax
        mov     ebx, eax
        mov     r10d, eax
        ja      %%2
%%3:    prefetch [rsi+64]
        sub     eax, ecx
        lea     rdi, [rdi+8]
        mov     [rsi], eax
        add     r9d, 65536
        cmp     [rsp], rsi
        mov     eax, [rsi+4]
        lea     rsi, [rsi+4]
        ja      %%1
%%4:    mov     rax, rdi
        END_PROC reg_save_list
%endmacro

        ; %1 = ot %2 = nt_sched
        expand 0, 0
        expand 0, 1

    end