
;  slinie(aux_ptr,aux_ptr_ub,sieve_interval)  
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

        FRAME_PROC slinie, 0, reg_save_list
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx,  r8

        cmp     rsi, rdi
        push    rbx
        push    r12
        jbe     .2
        sub     rsi, 8
.1:  
        movzx   rbx, word[rdi+2]
        movzx   r9, word[rdi]
        movzx   r8, word[rdi+6]
        sub     rbx, r9
        mov     r10b, [rdi+4]
        mov     r12, rbx
        mov     rax, rdx

%assign i 0
%rep    j_per_strip
%assign i i + 1
%define lab1(x) .L %+ x
%define lab2(x) .M %+ x
%define lab3(x) .N %+ x
        mov     rcx, r8
        xor     rbx, rbx
        add     r8, r12
        lea     rcx, [rax+rcx]
        cmovnc  rbx, r9
        add     rax, n_i
        add     r8, rbx
        lea     rbx, [r9+r9*4]
        sub     rax, rbx
lab1(i):add     [rcx+r9], r10b
        add     [rcx], r10b
        lea     rcx, [rcx+r9*2]
        add     [rcx+r9], r10b
        add     [rcx], r10b
        cmp     rax, rcx
        lea     rcx, [rcx+r9*2]
        ja      lab1(i)
        lea     rax, [rax+r9*4]
        cmp     rax, rcx
        jbe     lab2(i)
        add     [rcx], r10b
        add     [rcx+r9], r10b
        lea     rcx, [rcx+r9*2]
lab2(i):  
        lea     rax, [rax+r9]
        cmp     rax, rcx
        jbe     lab3(i)
        add     [rcx], r10b
lab3(i):  
%endrep  
        cmp     rsi, rdi
        lea     rdi, [rdi+8]
        ja      .1
.2:     END_PROC reg_save_list
        
        end
