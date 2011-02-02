
; u32_t *function(u32_t*,  u32_t, u16_t*, u16_t*)
; GAS    rax         rdi,    esi,    rdx,    rcx
; YASM   rax         rcx,    edx,     r8,     r9 

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi

        FRAME_PROC asm_MMX_Td8, 0, reg_save_list
        mov     rdi, rcx
        mov     esi, edx
        mov     rdx, r8
        mov     ecx, r9d

        mov     r8d, esi
        movdqa  xmm2, [rdx]
        shl     r8d, 16
        movdqa  xmm0, [rdx+16]
        or      esi, r8d
        movdqa  xmm1, [rdx+32]
        lea     rdx, [rdx+48]
        movd    xmm3, rsi
        mov     rax, rcx
        lea     rax, [rax+24]
        pshufd  xmm3, xmm3, 0
        pxor    xmm4, xmm4
        paddw   xmm2, xmm3
        movd    r10d, xmm3
.1:     cmp     rax, rdx
        pmullw  xmm1, xmm2
        movdqa  xmm2, [rdx]
        jbe     .3
        pmulhw  xmm0, xmm1
        movdqa  xmm1, [rdx+32]
        pcmpeqw xmm0, xmm4
        paddw   xmm2, xmm3
        pmovmskb esi, xmm0
        movdqa  xmm0, [rdx+16]
        lea     rdx, [rdx+48]
        test    esi, esi
        jz      .1
        and     esi, 0x5555
        xor     r8, r8
.2:     bsf     rcx, rsi
        shr     esi, cl
        shr     rcx, 1
        add     r8, rcx
        movzx   r9d, word[rdx+r8*2-80]
        shr     esi, 2
        mov     [rdi], r9d
        inc     r8d
        test    esi, esi
        lea     rdi, [rdi+4]
        jz      .1
        jmp     .2

        align   4
.3:     mov     rax, rdi
        END_PROC reg_save_list
        
        FRAME_PROC asm_TdUpdate8, 0, reg_save_list
        mov     rdi, rcx
        mov     esi, edx
        mov     rdx, r8
        mov     ecx, r9d

        mov     rax, rdi
        mov     rcx, rsi
        cmp     rsi, rax
        mov     r8, rdx
        lea     rcx, [rcx-48]
        movdqa  xmm2, [rax]
        movdqa  xmm0, [r8]
        lea     r8, [r8+16]
        movdqa  xmm1, [rax+16]
        movdqa  xmm3, xmm0
        jbe     .5
.4:     pcmpgtw xmm3, xmm2
        psubw   xmm2, xmm0
        cmp     rcx, rax
        movdqa  xmm0, [r8]
        pand    xmm3, xmm1
        movdqa  xmm1, [rax+64]
        lea     r8, [r8+16]
        paddw   xmm3, xmm2
        movdqa  xmm2, [rax+48]
        movdqa  [rax], xmm3
        movdqa  xmm3, xmm0
        lea     rax, [rax+48]
        ja      .4
.5:     END_PROC reg_save_list

        end
