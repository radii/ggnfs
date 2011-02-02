
;  tdslinie3(aux_ptr,aux_ptr_ub,sieve_interval,tds_buffer)  
;  Now, the registers which we are going to use  
;  The bx-register may also be used for auxilliary 32-bit values if sv1  
;  is not used  
;  Offset of the various things from this pointer  
;  We store the int difference projective_root-prime here:  

; void function(u16_t, u16_t*, uchar*, u32_t**);
; GAS    rax      rdi,    rsi,    rdx,     rcx
; YASM   rax      rcx,    rdx,     r8,      r9
      
%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi, rbx, r12

        FRAME_PROC tdslinie3, 0, reg_save_list
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx,  r8
        mov     rcx,  r9

        cmp     rsi, rdi
        jbe     .2
        sub     rsi, 8
.1:     movzx   r12, word[rdi+2]
        movzx   r11, word[rdi]
        movzx   r10, word[rdi+6]
        sub     r12, r11
        mov     rbx, r12
        mov     r9, rdx
%assign i 0
%rep    j_per_strip
%assign i i + 1
%define lb1(x) .L %+ x 
%define lb2(x) .M %+ x 
%define lb3(x) .N %+ x 
%define lb4(x) .O %+ x 
        mov     r8, r10
        xor     r12, r12
        add     r10, rbx
        lea     r8, [r9+r8]
        cmovnc  r12, r11
        add     r9, n_i
        add     r10, r12
        mov     al, [r8]
        sub     r9, r11
        or      al, [r8+r11]
        lea     r8, [r8+r11*2]
        or      al, [r8]
        cmp     r9, r8
        jbe     lb1(i)
        or      al, [r8+r11]
lb1(i): test    al, al
        lea     r9, [r9+r11]
        jz      lb4(i)
        sub     r8, r11
        sub     r8, r11
lb2(i): movzx   r12, byte[r8]
        test    r12, r12
        lea     r8, [r8+r11]
        jz      lb3(i)
        dec     r12
        lea     r12, [rcx+r12*8]
        mov     rax, [r12]
        mov     [rax], r11d
        lea     rax, [rax+4]
        mov     [r12], rax
lb3(i): cmp     r9, r8
        ja      lb2(i)
lb4(i):  
%endrep    
        cmp     rsi, rdi
        mov     [rdi+6], r10
        lea     rdi, [rdi+8]
        ja      .1
.2:     END_PROC reg_save_list

        end
