
;  slinie3(aux_ptr,aux_ptr_ub,sieve_interval)  
;  Now, the registers which we are going to use  
;  The bx-register may also be used for auxilliary 32-bit values if sv1  
;  is not used  
;  Offset of the various things from this pointer  
;  We store the int difference projective_root-prime here:  
        
; void slinie (u16_t*, u16_t*, uchar*)
; GAS    rax      rdi,    rsi,    rdx
; YASM   rax      rcx,    rdx,     r8

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi, rbx, r12

%if j_per_strip != 1
%define j_per_strip_minus1 j_per_strip - 1

        FRAME_PROC slinie3, 0, reg_save_list
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx,  r8

        cmp     rsi, rdi
        jbe     .3
        sub     rsi, 8
.1:     movzx   rbx, word[rdi+2]
        movzx   r9, word[rdi]
        movzx   r8, word[rdi+6]
        sub     rbx, r9
        mov     r10b, [rdi+4]
        mov     r12, rbx
        mov     rax, rdx

%assign i 0
%rep    j_per_strip_minus1
%assign i i + 1
%define label(x) .L %+ x
        mov     rcx, r8
        xor     rbx, rbx
        add     r8, r12
        lea     rcx, [rax+rcx]
        cmovnc  rbx, r9
        add     rax, n_i
        add     r8, rbx
        mov     r11b, [rcx+r9]
        sub     rax, r9
        add     [rcx], r10b
        add     r11b, r10b
        mov     [rcx+r9], r11b
        lea     rcx, [rcx+r9*2]
        mov     r11b, [rcx]
        add     r11b, r10b
        cmp     rax, rcx
        mov     [rcx], r11b
        lea     rax, [rax+r9]
        jbe     label(i)
        add     [rcx+r9], r10b
label(i):  
%endrep   
        mov     rcx, r8
        lea     rcx, [rax+rcx]
        add     rax, n_i
        mov     r11b, [rcx+r9]
        sub     rax, r9
        add     [rcx], r10b
        add     r11b, r10b
        mov     [rcx+r9], r11b
        lea     rcx, [rcx+r9*2]
        mov     r11b, [rcx]
        add     r11b, r10b
        cmp     rax, rcx
        mov     [rcx], r11b
        lea     rax, [rax+r9]
        jbe     .2
        add     [rcx+r9], r10b
.2:     cmp     rsi, rdi
        lea     rdi, [rdi+8]
        ja      .1
.3:     END_PROC reg_save_list

%else
%define j_per_strip_minus1  1
        
        FRAME_PROC slinie3, 0, reg_save_list

        cmp     rsi, rdi
        jbe     .3
        sub     rsi, 8
.1:     movzx   rbx, word[rdi+2]
        movzx   r9, word[rdi]
        movzx   r8, word[rdi+6]
        sub     rbx, r9
        mov     r10b, [rdi+4]
        mov     r12, rbx
        mov     rax, rdx
        mov     rcx, r8
        lea     rcx, [rax+rcx]
        add     rax, n_i
        mov     r11b, [rcx+r9]
        sub     rax, r9
        add     [rcx], r10b
        add     r11b, r10b
        mov     [rcx+r9], r11b
        lea     rcx, [rcx+r9*2]
        mov     r11b, [rcx]
        add     r11b, r10b
        cmp     rax, rcx
        mov     [rcx], r11b
        lea     rax, [rax+r9]
        jbe     .2
        add     [rcx+r9], r10b
.2:     cmp     rsi, rdi
        lea     rdi, [rdi+8]
        ja      .1
.3:     END_PROC reg_save_list

%endif
        end
