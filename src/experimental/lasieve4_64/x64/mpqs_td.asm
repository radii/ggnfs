; Copyright (C) 2004 Jens Franke, T.Kleinjung  
; This file is part of gnfs4linux, distributed under the terms of the  
; GNU General Public Licence and WITHOUT ANY WARRANTY.  
        
; You should have received a copy of the GNU General Public License along  
; with this program; see the file COPYING.  If not, write to the Free  
; Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
; 02111-1307, USA.  
        
; Written by T. Kleinjung  
; Modifications by J. Franke  

; u32_t function(u16_t**, u16_t, u64_t*, u32_t*)
; GAS    rax         rdi,   esi,    rdx,    rcx
; YASM   rax         rcx,   edx,     r8      r9

%include "ls-defs.inc"

        bits 64
        text

        extern  abort
        extern  mpqs_256_inv_table
        extern  mpqs_Adiv_all
        extern  mpqs_FB
        extern  mpqs_FB_A_inv
        extern  mpqs_FB_inv
        extern  mpqs_FB_inv_info
        extern  mpqs_FB_start
        extern  mpqs_FBk
        extern  mpqs_FBk_inv
        extern  mpqs_nAdiv_total
        extern  mpqs_nFB
        extern  mpqs_nFBk
        extern  mpqs_nFBk_1
        extern  mpqs_td_begin

%define reg_save_list rsi, rdi, rbx, r12, r13

        FRAME_PROC asm_td, 0, reg_save_list 
        mov     rdi, rcx
        mov     esi, edx
        mov     rdx, r8
        mov     rcx, r9

        mov     r8, rdx
        movzx   rax, word[rdi]    ; ind 
        mov     rcx, rax
        shl     rax, 16
        or      rax, rcx    ; ind ind 
        mov     rcx, rax
        shl     rax, 32
        or      rax, rcx    ; ind ind ind ind 
        movd    xmm5, rax
        movd    xmm3, rax
        pslldq  xmm5, 8
        paddw   xmm3, xmm5    ; ind ind ind ind ind ind ind ind 
        pxor    xmm5, xmm5
        
        movzx   r9, word[rdi+8]    ; %r9         
        movzx   rcx, word[rip+mpqs_td_begin]
        mov     rbx, mpqs_FB_inv_info
        mov     r11, mpqs_FB_start
        mov     r13w, [rip+mpqs_nFBk_1]
        add     r13w, 4
        movaps  xmm0, [rbx]
        movaps  xmm1, [rbx+16]
        movaps  xmm2, [r11]
        movaps  xmm4, xmm0
        psubw   xmm4, xmm2
        paddw   xmm4, xmm3    ; ind+p-s1 ind+p-s2 ... 
        pmullw  xmm4, xmm1
        movaps  xmm2, [r11+16]
        lea     rbx, [rbx+32]
        lea     r11, [r11+16]
        pmulhuw xmm4, xmm0
        movaps  xmm0, [rbx]
        movaps  xmm1, [rbx+16]
        psubw   xmm2, xmm3
        pcmpeqw xmm4, xmm5
        pmovmskb r12d, xmm4
        and     r12d, 0xfff0    ; only considering primes %r9 1 2 3 
        movaps  xmm4, xmm0
        psubw   xmm4, xmm2
        jz      .1
        mov     eax, r12d
        sub     r13w, 3    ; %r13w=mpqs_nFBk_1+1 
        shr     eax, 2
        or      r12d, eax        
        shr     r12d, 5
        mov     [rdi+r9*2+10], r13w
        adc     r9, 0
        inc     r13w    ; %r13w=mpqs_nFBk_1+2 
        shr     r12d, 4
        mov     [rdi+r9*2+10], r13w
        adc     r9, 0
        inc     r13w    ; %r13w=mpqs_nFBk_1+3 
        shr     r12d, 4
        mov     [rdi+r9*2+10], r13w
        adc     r9, 0
        inc     r13w    ; %r13w=mpqs_nFBk_1+4 
.1:     add     r13w, 4
        sub     rcx, 4
        lea     rbx, [rbx+32]
        lea     r11, [r11+16]
        jz      .2    
        pmullw  xmm4, xmm1
        movaps  xmm2, [r11]
        psubw   xmm2, xmm3
        pmulhuw xmm4, xmm0
        movaps  xmm0, [rbx]
        movaps  xmm1, [rbx+16]
        pcmpeqw xmm4, xmm5
        pmovmskb r12d, xmm4
        test    r12d, r12d
        movaps  xmm4, xmm0
        psubw   xmm4, xmm2
        jz      .1
        mov     eax, r12d
        sub     r13w, 4
        shr     eax, 2
        or      r12d, eax    ; significant bits at position 0 4 8 12         
        shr     r12d, 1
        mov     [rdi+r9*2+10], r13w
        adc     r9, 0
        inc     r13w
        shr     r12d, 4
        mov     [rdi+r9*2+10], r13w
        adc     r9, 0
        inc     r13w
        shr     r12d, 4
        mov     [rdi+r9*2+10], r13w
        adc     r9, 0
        inc     r13w
        shr     r12d, 4
        mov     [rdi+r9*2+10], r13w
        adc     r9, 0
        inc     r13w
        jmp     .1
.2:   
        mov     r10, mpqs_FB
        movzx   rdx, word[rip+mpqs_nFBk_1]
        add     rdx, rdx
        add     rdx, rdx
        sub     r10, rdx
        mov     rcx, 0
        mov     rax, 1
.3:     cmp     rcx, r9
        jnc     .4
        movzx   rdx, word[rdi+rcx*2+10]
        movzx   rdx, word[r10+rdx*4]
        mul     rdx
        inc     rcx
        jmp     .3
.4:     mov     rbx, rax
        test    rsi, 1
        mov     rsi, r9
        jz      .5
        mov     word[rdi+r9*2+10], 0
        inc     r9  
.5:     test    r8, 0x00000001
        jnz     .6
        inc     r9
        cmp     r9, 27
        jnc     .16
        mov     cx, [rip+mpqs_nFBk_1]
        mov     [rdi+r9*2+8], cx
        shr     r8, 1
        jmp     .5  
.6:     mov     r11, rbx
        mov     rdx, rbx
        mov     r13, rbx
        mov     rcx, rbx
        and     rdx, 0xff
        shr     r13, 32
        shl     rcx, 32
        shr     rdx, 1
        inc     rcx
        mov     eax, ebx
        cmp     r8, rcx
        movzx   r11d, byte[rdx+mpqs_256_inv_table]
        adc     r13, 0
        mul     r11d
        test    r13, r13
        jz      .16
        and     eax, 0xff00
        mul     r11d
        sub     r11d, eax
        mov     eax, ebx
        mul     r11d
        and     eax, 0xffff0000
        mul     r11d
        sub     r11d, eax
        mov     eax, r8d
        mul     r11d
        mov     r8d, eax
        mov     r11, mpqs_FB_inv
        movzx   rcx,  word[rip+mpqs_nFBk_1]
        add     rcx, rcx
        add     rcx, rcx
        sub     r11, rcx
        mov     rcx, 0
.7:     test    rsi, rsi
        mov     eax, r8d
        jz      .9
        movzx   r13, word[rdi+rsi*2+8]    ; bx: ii 
        dec     rsi
        movzx   ecx, word[r10+r13*4]    ; cx: p 
        mov     ebx, [r11+r13*4]    ; edx: inv 
.8:     mul     ebx
        mov     r12d, eax
        mul     ecx
        test    edx, edx
        jnz     .7
        cmp     r9, 27
        jnc     .16
        mov     eax, r12d
        mov     [rdi+r9*2+10], r13w
        inc     r9
        mov     r8d, r12d
        jmp     .8     
.9:     cmp     r8d, 1
        jz      .15  
        mov     r11, mpqs_FBk_inv
        mov     r10, mpqs_FBk
        movzx   rsi, word[rip+mpqs_nFBk]
        inc     rsi
.10:    dec     rsi
        mov     eax, r8d
        cmp     rsi, 0
        jz      .12
        movzx   ecx, word[r10+rsi*2-2]    ; cx: p 
        mov     ebx, [r11+rsi*4-4]    ; %rbx: inv 
.11:    mul     ebx
        mov     r12d, eax    ; rr 
        mul     ecx
        test    edx, edx
        jnz     .10
        cmp     r9, 27
        jnc     .16
        mov     eax, r12d
        mov     [rdi+r9*2+10], rsi
        inc     r9
        mov     r8d, r12d
        jmp     .11    
.12:    cmp     eax, 1
        jz      .15  
        mov     r11, mpqs_FB_A_inv
        mov     r10, mpqs_Adiv_all
        mov     r13w, [rip+mpqs_nFB]
        add     r13w, [rip+mpqs_nFBk]
        movzx   rsi, word[mpqs_nAdiv_total]
        inc     rsi
.13:    dec     rsi
        mov     eax, r8d
        cmp     rsi, 0
        jz      .15
        movzx   ecx, word[r10+rsi*2-2]    ; cx: p 
        mov     ebx, [r11+rsi*4-4]    ; %rbx: inv 
.14:    mul     ebx
        mov     r12d, eax    ; rr 
        mul     ecx
        test    edx, edx
        jnz     .13
        mov     eax, r12d
        cmp     r9, 27
        jnc     .16
        add     r13w, si
        mov     [rdi+r9*2+10], r13w
        inc     r9
        sub     r13w, si
        mov     r8d, r12d
        jmp     .14    
.15:    mov     [rdi+8], r9w
        mov     eax, r8d
        emms
        EXIT_PROC reg_save_list 
.16:    xor     rax, rax
        emms    
        EXIT_PROC reg_save_list
.17:    call    abort
        END_PROC reg_save_list

        end
