
; lasieve_search0(uchar*, uchar*, uchar*, uchar*,   uchar*,   u16_t*, uchar*);
; GAS    rax         rdi,     rsi,   rdx,    rcx,       r8,       r9, [rsp+8]
; YASM   rax         rcx,     rdx,    r8,     r9, [rsp+40], [rsp+48], [rsp+56],

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi, rbx, r12, r13, r14, r15

        FRAME_PROC lasieve_search0, 0, reg_save_list
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx, r8
        mov     rcx, r9
        mov     r8, qword [rsp+stack_use+40]
        mov     r9, qword [rsp+stack_use+48]
        mov     r14, qword [rsp+stack_use+56]

        mov     r10, rdi
        lea     r15, [rdi+128]
        mov     r11, rcx
        mov     r13, r10
.1:     mov     al, [r11]
        mov     bl, [rsi]
        inc     bl
        lea     rsi, [rsi+1]
        sub     al, bl
        jb      .9
        movzx   rcx, al
        movd    xmm7, rcx
        punpcklbw xmm7, xmm7
        xor     rbx, rbx
        punpcklwd xmm7, xmm7
        pshufd  xmm7, xmm7, 0
.2:     movdqa  xmm0, [r10]
        movdqa  xmm1, [r10+16]
        movdqa  xmm4, xmm7
        movdqa  xmm2, [r10+32]
        pmaxub  xmm4, xmm0
        movdqa  xmm3, [r10+48]
        movdqa  xmm5, xmm2
        pmaxub  xmm4, xmm1
        pmaxub  xmm5, xmm3
        lea     r10, [r10+64]
        pmaxub  xmm5, xmm4
        pcmpeqb xmm5, xmm7
        pmovmskb rbx, xmm5
        cmp     rbx, 65535
        jne     .7
.3:     cmp     r15, r10
        ja      .2
.4:     cmp     rdx, rsi
        lea     r15, [r15+n_i]
        lea     r10, [r10+(n_i-128)]
        ja      .1
.5:     lea     rsi, [rdx+-j_per_strip]
        lea     r11, [r11+1]
        mov     r10, r13
        cmp     r8, r11
        lea     r15, [r10+256]
        lea     r10, [r10+128]
        mov     r13, r10
        ja      .1
.6:     mov     rax, r14
        sub     rax, [rsp+48]
        emms    
        EXIT_PROC reg_save_list
        
        align   4
.7:     pmaxub  xmm0, xmm7
        pmaxub  xmm1, xmm7
        pmaxub  xmm2, xmm7
        pmaxub  xmm3, xmm7
        pcmpeqb xmm0, xmm7
        pcmpeqb xmm1, xmm7
        pcmpeqb xmm2, xmm7
        pcmpeqb xmm3, xmm7
        pmovmskb rax, xmm0
        pmovmskb rcx, xmm2
        sal     rcx, 32
        pmovmskb rbx, xmm1
        or      rax, rcx
        pmovmskb rcx, xmm3
        sal     rbx, 16
        sal     rcx, 48
        or      rax, rbx
        xor     rbx, rbx
        or      rax, rcx
        lea     r10, [r10-64]
        xor     rax, 0xffffffffffffffff
        sub     r10, rdi
.8:     bsf     rcx, rax
        add     rbx, rcx
        add     r10, rbx
        shr     rax, cl
        mov     [r9], r10w
        lea     r9, [r9+2]
        mov     cl, [rdi+r10]
        sub     r10, rbx
        inc     rbx
        shr     rax, 1
        test    rax, rax
        mov     [r14], cl
        lea     r14, [r14+1]
        jnz     .8
        lea     r10, [r10+64]
        lea     r10, [rdi+r10]
        jmp     .3
        
        align   4
.9:     sub     r10, rdi
        mov     rbx, r10
        shl     rbx, 16
        add     rbx, r10
        add     rbx, 0x00010000
        movd    xmm0, rbx
        add     rbx, 0x00020002
        movd    xmm1, rbx
        add     rbx, 0x00020002
        mov     rax, 0x00080008
        psllq   xmm1, 32
        movd    xmm2, rbx
        add     rbx, 0x00020002
        movd    xmm3, rax
        por     xmm0, xmm1
        add     r10, rdi
        movd    xmm1, rbx
        movq    xmm4, xmm3
        psllq   xmm1, 32
        psllq   xmm4, 32
        por     xmm1, xmm2
        por     xmm3, xmm4
.10:    movq    xmm2, [r10]
        movq    xmm4, [r10+8]
        movq    [r9], xmm0
        lea     r10, [r10+16]
        paddw   xmm0, xmm3
        movq    [r9+8], xmm1
        cmp     r15, r10
        paddw   xmm1, xmm3
        movq    [r14], xmm2
        movq    [r9+16], xmm0
        paddw   xmm0, xmm3
        movq    [r14+8], xmm4
        lea     r14, [r14+16]
        movq    [r9+24], xmm1
        lea     r9, [r9+32]
        paddw   xmm1, xmm3
        ja      .10
        jmp     .4

        END_PROC reg_save_list

        end
