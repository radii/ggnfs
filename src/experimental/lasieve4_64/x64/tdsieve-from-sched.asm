
; u16_t **tdsieve_sched2buf(u16_t**, u16_t*, unsigned char*, u16_t**,  u16_t**);
; GAS    rax                    rdi,    rsi,            rdx,     rcx,       r8
; YASM   rax                    rcx,    rdx,             r8,      r9,  [rsp+40]

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi

        FRAME_PROC tdsieve_sched2buf, 0, reg_save_list 
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx, r8
        mov     rcx, r9
        mov      r8, qword [rsp+stack_use+40]

        mov     r9, [rdi]
        lea     rsi, [rsi-12]
        cmp     rsi, r9
        jbe     .3
.1:     prefetcht0 [r9+128]
        movzx   r10, word[r9]
        movzx   r11, word[r9+4]
        mov     al, [rdx+r10]
        movzx   r10, word[r9+8]
        or      al, [rdx+r11]
        movzx   r11, word[r9+12]
        or      al, [rdx+r10]
        or      al, [rdx+r11]
        test    al, al
        jnz     .8
.2:     lea     r9, [r9+16]
        cmp     rsi, r9
        ja      .1
.3:     lea     rsi, [rsi+8]
        cmp     rsi, r9
        lea     rsi, [rsi+4]
        movzx   r10, word[r9]
        jbe     .5
        mov     al, [rdx+r10]
        movzx   r11, word[r9+4]
        test    al, al
        mov     al, [rdx+r11]
        jz      .4
        mov     [rcx], r9
        lea     rcx, [rcx+8]
.4:     test    al, al
        movzx   r10, word[r9+8]
        lea     r11, [r9+4]
        lea     r9, [r9+8]
        jz      .5
        mov     [rcx], r11
        lea     rcx, [rcx+8]
.5:     cmp     rsi, r9
        jbe     .7
        mov     al, [rdx+r10]
        test    al, al
        jz      .6
        mov     [rcx], r9
        lea     rcx, [rcx+8]
.6:     lea     r9, [r9+4]
.7:     mov     rax, rcx
        mov     [rdi], r9
        EXIT_PROC reg_save_list

        align   4
.8:     movzx   r10, word[r9]
        mov     al, [rdx+r10]
        movzx   r11, word[r9+4]
        test    al, al
        jz      .9
        mov     [rcx], r9
        lea     rcx, [rcx+8]
.9:     mov     al, [rdx+r11]
        movzx   r10, word[r9+8]
        test    al, al
        lea     r11, [r9+4]
        jz      .10
        mov     [rcx], r11
        lea     rcx, [rcx+8]
.10:     mov     al, [rdx+r10]
        lea     r10, [r9+8]
        movzx   r11, word[r9+12]
        test    al, al
        jz      .11
        mov     [rcx], r10
        lea     rcx, [rcx+8]
.11:    mov     al, [rdx+r11]
        lea     r11, [r9+12]
        test    al, al
        jz      .12
        mov     [rcx], r11
        lea     rcx, [rcx+8]
.12:  
        cmp     r8, rcx
        ja      .2
        lea     r9, [r9+16]
        jmp     .7
        END_PROC reg_save_list

        end
