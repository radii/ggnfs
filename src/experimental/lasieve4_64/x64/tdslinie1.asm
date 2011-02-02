
;  tdslinie1(aux_ptr,aux_ptr_ub,sieve_interval,tds_buffer)  
;  We save %ebx,%esi,%edi,%ebp and also have one auto variable and the  
;  return address on the stack. Therefore, the stack offset of the  
;  first arg is 24.  
;  Now, the registers which we are going to use  
;  Offset of the various things from this pointer  
;  We store the int difference projective_root-prime here:  
;  This macro is taken from the GNU info documentation of m4.  

; void function(u16_t, u16_t*, uchar*, u32_t**);
; GAS    rax      rdi,    rsi,    rdx,     rcx
; YASM   rax      rcx,    rdx,     r8,      r9

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi, rbx, r12

        FRAME_PROC tdslinie1, 0, reg_save_list
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx,  r8
        mov     rcx,  r9

        cmp     rsi, rdi
        push    rbx
        push    r12
        jbe     .3
        sub     rsi, 8
        lea     rcx, [rcx-8]
.1:     movzx   r12, word[rdi+2]
        movzx   r11, word[rdi]
        movzx   r10, word[rdi+6]
        sub     r12, r11
        mov     rbx, r12
        mov     r9, rdx
        xor     al, al

%assign i 0
%rep    j_per_strip
%assign i i + 1
        mov     r8, r10
        xor     r12, r12
        add     r10, rbx
        lea     r8, [r9+r8]
        cmovnc  r12, r11
        add     r9, n_i
        or      al, [r8]
        add     r10, r12
        sub     r9, r11
        xor     r12, r12
        cmp     r9, r8
        lea     r9, [r9+r11]
        cmova   r12, [r8+r11]
        or      al, r12b
%endrep
        test    al, al
        jnz     .4
.2:     cmp     rsi, rdi
        mov     [rdi+6], r10
        lea     rdi, [rdi+8]
        ja      .1
.3:     EXIT_PROC reg_save_list
 
        align   8
.4:     movzx   r12, word[rdi+2]
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
        mov     r8, r10
        xor     r12, r12
        add     r10, rbx
        lea     r8, [r9+r8]
        cmovnc  r12, r11
        add     r9, n_i
        add     r10, r12
        sub     r9, r11
        mov     al, [r8]
        cmp     r9, r8
        lea     r9, [r9+r11]
        jbe     lb1(i)
        or      al, [r8+r11]
lb1(i): test    al, al
        jz      lb3(i)
        movzx   r12, byte[r8]
        test    r12, r12
        lea     r8, [r8+r11]
        jz      lb2(i)
        mov     rax, [rcx+r12*8]
        mov     [rax], r11d
        lea     rax, [rax+4]
        mov     [rcx+r12*8], rax
lb2(i): cmp     r9, r8
        jbe     lb3(i)
        movzx   r12, byte[r8]
        test    r12, r12
        jz      lb3(i)
        mov     r8, [rcx+r12*8]
        mov     [r8], r11d
        lea     r8, [r8+4]
        mov     [rcx+r12*8], r8
lb3(i):  
%endrep
        jmp     .2

        END_PROC reg_save_list

        end
