
; u32_t  *function(u32_t*, u32_t*, u32_t*, u32_t**,    u32_t)
; GAS    rax          rdi,    rsi,    rdx,     rcx,       r8
; YASM   rax          rcx,    rdx,     r8,      r9, [rsp+40] 

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi, rbx, r12

        FRAME_PROC medsched0, 0, reg_save_list 
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx,  r8
        mov     rcx,  r9
        mov     r8d, dword [rsp+stack_use+40]
        
        shl     r8d, 16
        cmp     rdx, rsi
        lea     rdx, [rdx-4]
        mov     r12, [rcx]
        jbe     .4
.1:     mov     r9d, [rsi]
        mov     eax, [rdi]
        prefetcht0 [rdi+128]
        cmp     r9d, l1_size
        mov     ebx, [rdi+4]
        jae     .3
        neg     eax
        neg     ebx
        mov     r10d, r9d
        xor     r11d, r11d
        and     eax, n_i_mask
        and     ebx, n_i_mask
.2:     and     r10d, n_i_mask
        or      r9d, r8d
        cmp     eax, r10d
        mov     [r12], r9d
        lea     r12, [r12+4]
        cmovbe  r11d, [rdi+4]
        and     r9d, l1_mask
        cmp     ebx, r10d
        mov     r10d, 0
        cmova   r10d, [rdi]
        add     r9d, r11d
        add     r9d, r10d
        xor     r11d, r11d
        cmp     r9d, l1_size
        mov     r10d, r9d
        jb      .2
.3:     prefetcht0 [rsi+128]
        sub     r9d, l1_size
        add     r8d, 65536
        cmp     rdx, rsi
        mov     [rsi], r9d
        lea     rdi, [rdi+8]
        lea     rsi, [rsi+4]
        ja      .1
.4:     mov     [rcx], r12
        pop     r12
        mov     rax, rdi
        END_PROC reg_save_list

        end
