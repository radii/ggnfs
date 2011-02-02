
;  auxreg will be destroyed by squaring.  
;  other auxregs which may be used outside the squaring  
;  are %rax and %rdx  
;  auxreg2 is preserved by squaring.  

; u64_t   pt64(u64_t);
; GAS    rax     rdi
; YASM   rax     rcx
       
%include "ls-defs.inc"

        bits 64
        text
    
        extern  mpqs_256_inv_table

%define reg_save_list rsi, rdi

modsq64:mov     rax, r11
        mul     r11
        mov     r11, rdx
        mul     r8
        mul     rdi
        xor     rax, rax
        sub     r11, rdx
        cmovb   rax, rdi
        add     r11, rax
        ret     
        
        FRAME_PROC pt64, 0, reg_save_list
        mov     rdi, rcx
         
        mov     rdx, 1
        xor     rax, rax
        div     rdi
        mov     rcx, rdi
        mov     r9, rdi
        and     rcx, 255
        sub     r9, 1
        shr     rcx, 1
        lea     r8, [rel mpqs_256_inv_table]
        movzx   r8, byte[rcx+r8]
        mov     rsi, rdx
.1:  
%rep    3
        mov     rax, r8
        mul     r8
        shl     r8, 1
        mul     rdi
        sub     r8, rax
%endrep

.2:     bsf     rcx, r9
        shr     r9, cl
        mov     r10, rcx
.3:     bsr     rcx, r9
        mov     rax, 1
        shl     rax, cl
        mov     rcx, rax        
        mov     r11, rsi
        jmp     .5
.4:     call    modsq64
        test    r9, rcx
        jz      .6
.5:     xor     rax, rax
        shl     r11, 1
        cmovc   rax, rdi
        cmp     r11, rdi
        cmovae  rax, rdi
        sub     r11, rax
.6:     shr     rcx, 1
        jnz     .4
        cmp     rsi, r11
        mov     rcx, rdi
        je      .9
        sub     rcx, rsi
.7:     cmp     rcx, r11
        je      .9
        dec     r10
        jz      .8
        call    modsq64
        jmp     .7
        
        align   8
.8:     xor     rax, rax
        EXIT_PROC reg_save_list
        
        align   8
.9:     mov     rax, 1
        END_PROC reg_save_list

        end
