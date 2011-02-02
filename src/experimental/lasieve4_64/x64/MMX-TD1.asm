
; u32_t *function(u32_t*,  u32_t, u16_t*, u16_t*)
; GAS    rax         rdi,    esi,    rdx,    rcx
; YASM   rax         rcx,    edx,     r8,     r9 

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi

        FRAME_PROC asm_MMX_Td4, 0, reg_save_list
        mov     rdi, rcx
        mov     esi, edx
        mov     rdx, r8
        mov     ecx, r9d

        mov     r8d, esi
        movq    mm2, [rdx]
        shl     r8d, 16
        movq    mm0, [rdx+8]
        or      esi, r8d
        movq    mm1, [rdx+16]
        lea     rdx, [rdx+24]
        movd    mm3, esi
        mov     rax, rcx
        psllq   mm3, 32
        movd    mm4, esi
        lea     rax, [rax+24]
        por     mm3, mm4
        pxor    mm4, mm4
        paddw   mm2, mm3
.1:     cmp     rax, rdx
        pmullw  mm1, mm2
        movq    mm2, [rdx]
        jbe     .3
        pmulhw  mm0, mm1
        movq    mm1, [rdx+16]
        pcmpeqw mm0, mm4
        paddw   mm2, mm3
        pmovmskb esi, mm0
        movq    mm0, [rdx+8]
        lea     rdx, [rdx+24]
        test    esi, esi
        jz      .1
        and     esi, 0x55
        xor     r8, r8
.2:     bsf     rcx, rsi
        shr     esi, cl
        shr     rcx, 1
        add     r8, rcx
        movzx   r9d, word[rdx+r8*2-40]
        shr     esi, 2
        mov     [rdi], r9d
        inc     r8d
        test    esi, esi
        lea     rdi, [rdi+4]
        jz      .1
        jmp     .2

        align   4
.3:     emms    
        mov     rax, rdi
        END_PROC reg_save_list

        FRAME_PROC asm_TdUpdate4, 0, reg_save_list
        mov     rdi, rcx
        mov     esi, edx
        mov     rdx, r8
        mov     ecx, r9d

        mov     rax, rdi
        mov     rcx, rsi
        cmp     rsi, rax
        mov     r8, rdx
        lea     rcx, [rcx-24]
        movq    mm2, [rax]
        movq    mm0, [r8]
        lea     r8, [r8+8]
        movq    mm1, [rax+8]
        movq    mm3, mm0
        jbe     .5
.4:     pcmpgtw mm3, mm2
        psubw   mm2, mm0
        cmp     rcx, rax
        movq    mm0, [r8]
        pand    mm3, mm1
        movq    mm1, [rax+32]
        lea     r8, [r8+8]
        paddw   mm3, mm2
        movq    mm2, [rax+24]
        movq    [rax], mm3
        movq    mm3, mm0
        lea     rax, [rax+24]
        ja      .4
.5:     emms
        END_PROC reg_save_list
        
        end
