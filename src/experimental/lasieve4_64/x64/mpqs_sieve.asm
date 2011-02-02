
; void asm_sieve(void);

%include "ls-defs.inc"

        bits 64
        text
    
        extern  mpqs_FB
        extern  mpqs_FB_log
        extern  mpqs_FB_start
        extern  mpqs_sievearray
        extern  mpqs_sievebegin
        extern  mpqs_sievelen

%define reg_save_list rsi, rdi, rbx, r12, r13, r14, r15

        FRAME_PROC asm_sieve, 0, reg_save_list 

        movzx   rax, word[rip+mpqs_sievebegin]
        movzx   r10, word[rip+mpqs_sievelen]
        mov     rcx, mpqs_FB
        mov     r8, mpqs_FB_start
        mov     r9, mpqs_FB_log
        lea     rcx, [rcx+rax*4]
        lea     r8, [r8+rax*4]
        lea     r9, [r9+rax]
        mov     r11, r10
        shr     r10, 2
        mov     r13, [rip+mpqs_sievearray]
        mov     r15, r13
        add     r15, r11
        
        align   16
.1:     movzx   r14, word[rcx]    ; p 
        cmp     r10, r14
        lea     rcx, [rcx+4]
        jc      .4
        movzx   rsi, word[r8]
        movzx   rdi, word[r8+2]
        lea     r8, [r8+4]
        cmp     rdi, rsi    ; sort %rsi %rdi: 
        mov     r12, rsi
        mov     rdx, rdi
        mov     al, [r9]
        lea     r9, [r9+1]
        cmovc   rdi, r12
        cmovc   rsi, rdx
        add     rsi, r13
        add     rdi, r13
        mov     rbx, r15
        lea     r12, [r14+r14*2]
        lea     rdx, [r14+r14]
        sub     rbx, r12    
.2:     add     [rdi], al
        add     [rdi+r14], al
        add     [rsi], al
        add     [rsi+r14], al
        add     [rdi+rdx], al
        add     [rdi+r12], al
        add     [rsi+rdx], al
        add     [rsi+r12], al
        lea     rdi, [rdi+rdx*2]
        lea     rsi, [rsi+rdx*2]
        cmp     rdi, rbx
        jc      .2
        add     rbx, rdx
        cmp     rdi, rbx
        jnc     .3
        add     [rsi], al
        add     [rsi+r14], al
        add     [rdi], al
        add     [rdi+r14], al
        add     rsi, rdx
        add     rdi, rdx
.3:     cmp     rsi, r15
        cmovnc  rsi, r15
        cmp     rdi, r15
        cmovnc  rdi, r15
        add     [rsi], al
        add     rsi, r14
        add     [rdi], al
        cmp     rsi, r15
        cmovnc  rsi, r15
        add     [rsi], al
        jmp     .1
.4:     mov     rax, r11
        xor     rdx, rdx    ; %rdx is not used at the moment 
        mov     r10, 3
        div     r10
        mov     r10, rax    ; %r10=%r11/3 in this part 
        lea     rcx, [rcx-4]   
.5:     movzx   r14, word[rcx]    ; p 
        cmp     r10, r14
        lea     rcx, [rcx+4]
        jc      .6
        movzx   rsi, word[r8]
        movzx   rdi, word[r8+2]
        lea     r8, [r8+4]
        add     rsi, r13
        add     rdi, r13
        mov     al, [r9]
        lea     r9, [r9+1]
        lea     r12, [r14+r14*2]
        lea     rdx, [r14+r14]
        add     [rsi], al
        add     [rsi+r14], al
        add     [rdi], al
        add     [rdi+r14], al
        add     [rsi+rdx], al
        add     [rdi+rdx], al
        add     rsi, r12
        add     rdi, r12
        cmp     rsi, r15
        cmovnc  rsi, r15
        cmp     rdi, r15
        cmovnc  rdi, r15
        add     [rsi], al
        add     [rdi], al
        jmp     .5
.6:     mov     r10, r11
        shr     r10, 1
        lea     rcx, [rcx-4]    
.7:     movzx   r14, word[rcx]    ; p 
        cmp     r10, r14
        lea     rcx, [rcx+4]
        jc      .8    
        movzx   rsi, word[r8]
        movzx   rdi, word[r8+2]
        lea     r8, [r8+4]
        add     rsi, r13
        add     rdi, r13
        mov     al, [r9]
        lea     r9, [r9+1]
        lea     rdx, [r14+r14]
        add     [rsi], al
        add     [rsi+r14], al
        add     [rdi], al
        add     [rdi+r14], al
        add     rsi, rdx
        add     rdi, rdx
        cmp     rsi, r15
        cmovnc  rsi, r15
        cmp     rdi, r15
        cmovnc  rdi, r15
        add     [rsi], al
        add     [rdi], al
        jmp     .7
.8:     mov     r10, r11
        sub     r10, 1
        lea     rcx, [rcx-4]     
.9:     movzx   r14, word[rcx]    ; p 
        cmp     r10, r14
        lea     rcx, [rcx+4]
        jc      .10
        movzx   rsi, word[r8]
        movzx   rdi, word[r8+2]
        lea     r8, [r8+4]
        add     rsi, r13
        add     rdi, r13
        mov     al, [r9]
        lea     r9, [r9+1]
        add     [rsi], al
        add     [rdi], al
        add     rsi, r14
        add     rdi, r14
        cmp     rsi, r15
        cmovnc  rsi, r15
        cmp     rdi, r15
        cmovnc  rdi, r15
        add     [rsi], al
        add     [rdi], al
        jmp     .9
.10:    mov     r10, 0xffff
        lea     rcx, [rcx-4]   
.11     movzx   r14, word[rcx]    ; p 
        cmp     r10, r14
        lea     rcx, [rcx+4]
        jz      .12
        movzx   rsi, word[r8]
        movzx   rdi, word[r8+2]
        lea     r8, [r8+4]
        add     rsi, r13
        add     rdi, r13
        mov     al, [r9]
        lea     r9, [r9+1]
        cmp     rsi, r15
        cmovnc  rsi, r15
        cmp     rdi, r15
        cmovnc  rdi, r15
        add     [rsi], al
        add     [rdi], al
        jmp     .11
.12:    END_PROC reg_save_list
        
        end
