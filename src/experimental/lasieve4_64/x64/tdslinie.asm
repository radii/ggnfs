

;  tdslinie(aux_ptr,aux_ptr_ub,sieve_interval,tds_buffer)  
;  We save %ebx,%esi,%edi,%ebp and also have one auto variable and the  
;  return address on the stack. Therefore, the stack offset of the  
;  first arg is 24.  
;  Now, the registers which we are going to use  
;  The ax-register may also be used for auxilliary 32-bit values if sv1  
;  is not used  
;  Offset of the various things from this pointer  
;  We store the int difference projective_root-prime here:  
;  This macro is taken from the GNU info documentation of m4.  

; void function(u16_t, u16_t*, uchar*, u32_t**);
; GAS    rax      rdi,    rsi,    rdx,     rcx
; YASM   rax      rcx,    rdx,     r8,      r9

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi, rbx, r12, r13, r14, r15

        FRAME_PROC tdslinie, 0, reg_save_list
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx,  r8
        mov     rcx,  r9

        cmp     rsi, rdi
        jbe     .2
        sub     rsi, 8
.1:     movzx   r15, word[rdi+2]
        movzx   r11, word[rdi]
        movzx   r10, word[rdi+6]
        sub     r15, r11
        mov     rbx, r15
        mov     r9, rdx

%assign i 0
%rep    j_per_strip
%assign i i + 1
%define lb1(x) .L %+ x 
%define lb2(x) .M %+ x 
%define lb3(x) .N %+ x 
%define lb4(x) .O %+ x 
%define lb5(x) .P %+ x 
%define lb6(x) .Q %+ x 
%define lb7(x) .R %+ x 
%define lb8(x) .S %+ x 
%define lb9(x) .T %+ x 
        mov     r8, r10
        xor     r15, r15
        add     r10, rbx
        lea     r8, [r9+r8]
        cmovnc  r15, r11
        add     r9, n_i
        add     r10, r15
        lea     r15, [r11+r11*4]
        sub     r9, r15
lb1(i): mov     al, [r8]
        or      al, [r8+r11]
        lea     r8, [r8+r11*2]
        or      al, [r8]
        or      al, [r8+r11]
        test    al, al
        jz      lb4(i)
        lea     r13, [r8+r11*2]
        neg     r11
        lea     r12, [r8+r11*2]
        neg     r11
lb2(i): movzx   r15, byte[r12]
        test    r15, r15
        lea     r12, [r12+r11]
        jz      lb3(i)
        dec     r15
        lea     r15, [rcx+r15*8]
        mov     r14, [r15]
        mov     [r14], r11d
        lea     r14, [r14+4]
        mov     [r15], r14
lb3(i): cmp     r13, r12
        ja      lb2(i)
lb4(i): cmp     r9, r8
        lea     r8, [r8+r11*2]
        ja      lb1(i)
        lea     r9, [r9+r11*4]
        xor     al, al
        cmp     r9, r8
        mov     r12, r8
        jbe     lb5(i)
        mov     al, [r8]
        or      al, [r8+r11]
        lea     r8, [r8+r11*2]
lb5(i): lea     r9, [r9+r11]
        cmp     r9, r8
        jbe     lb6(i)
        or      al, [r8]
lb6(i): test    al, al
        jz      lb9(i)
lb7(i): movzx   r15, byte[r12]
        test    r15, r15
        lea     r12, [r12+r11]
        jz      lb8(i)
        dec     r15
        lea     r15, [rcx+r15*8]
        mov     r14, [r15]
        mov     [r14], r11d
        lea     r14, [r14+4]
        mov     [r15], r14
lb8(i): cmp     r9, r12
        ja      lb7(i)
lb9(i):  
%endrep    
        cmp     rsi, rdi
        mov     [rdi+6], r10
        lea     rdi, [rdi+8]
        ja      .1
.2:     END_PROC reg_save_list

        end
