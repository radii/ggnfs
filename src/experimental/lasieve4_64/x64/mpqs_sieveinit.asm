

; Copyright (C) 2002 Jens Franke, T.Kleinjung  
; This file is part of gnfs4linux, distributed under the terms of the  
; GNU General Public Licence and WITHOUT ANY WARRANTY.  
        
; You should have received a copy of the GNU General Public License along  
; with this program; see the file COPYING.  If not, write to the Free  
; Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
; 02111-1307, USA.  

; void asm_function(unsigned char*, u32_t, ushort*, u64_t*, unsigned char*,    u32_t)
; GAS    rax                   rdi,   esi,     rdx,    rcx,             r8,      r9d
; YASM   rax                   rcx,   edx,      r8,     r9,       [rsp+40], [rsp+48]

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi, rbx, r12

        FRAME_PROC asm_sieve_init, 0, reg_save_list 
        mov     rdi, rcx
        mov     esi, edx
        mov     rdx, r8
        mov     rcx, r9
        mov      r8, [rsp+stack_use+40]
        mov      r9d, [rsp+stack_use+48]

        mov     rbx, rdx
        mov     rcx, [rcx]
        mov     r11, r8
        mov     r12, r8
        mov     r8, r9
        shl     r8, 4
        add     r12, r8    ; tinyend 
        movzx   rdx, word[rbx]    ; this is 0 
.1:     movzx   r10, word[rbx+4]
        sub     r10, rdx    ; length in byte 
        shr     r10, 4    ; length in 32 byte 
        movzx   rax, word[rbx+2]
        lea     rbx, [rbx+4]
        mul     rcx
        movd    mm2, rax

        align   16
.2:     movq    mm0, [r11]
        movq    mm1, [r11+8]
        lea     r11, [r11+16]
        xor     rax, rax
        cmp     r11, r12
        paddb   mm0, mm2
        paddb   mm1, mm2
        cmovnc  rax, r8
        sub     r11, rax
        dec     r10    ; length>0 
        movq    [rdi], mm0
        movq    [rdi+8], mm1
        lea     rdi, [rdi+16]
        jnz     .2
        movzx   rdx, word[rbx]
        cmp     rdx, rsi
        jnz     .1
        emms    
        END_PROC reg_save_list
               
        FRAME_PROC asm_sieve_init16, 0, reg_save_list 
        mov     rdi, rcx
        mov     esi, edx
        mov     rdx, r8
        mov     rcx, r9
        mov      r8, [rsp+stack_use+40]
        mov      r9d, [rsp+stack_use+48]

        mov     rbx, rdx
        mov     rcx, [rcx]
        mov     r11, r8
        mov     r12, r8
        mov     r8, r9
        shl     r8, 4
        add     r12, r8    ; tinyend 
        movzx   rdx, word[rbx]    ; this is 0 
.3:  
        movzx   r10, word[rbx+4]
        sub     r10, rdx    ; length in byte 
        shr     r10, 5    ; length in 32 byte 
        movzx   rax, word[rbx+2]
        lea     rbx, [rbx+4]
        mul     rcx
        movd    mm2, rax

        align   16
.4:     movq    mm0, [r11]
        movq    mm1, [r11+8]
        movq    mm4, [r11+16]
        movq    mm5, [r11+24]
        lea     r11, [r11+32]
        paddb   mm0, mm2
        paddb   mm1, mm2
        cmp     r11, r12
        mov     rax, 0
        cmovnc  rax, r8
        paddb   mm4, mm2
        paddb   mm5, mm2
        sub     r11, rax
        dec     r10    ; length>0 
        movq    [rdi], mm0
        movq    [rdi+8], mm1
        movq    [rdi+16], mm4
        movq    [rdi+24], mm5
        lea     rdi, [rdi+32]     
        jnz     .4
        movzx   rdx, word[rbx]
        cmp     rdx, rsi
        jnz     .3
        emms    
        END_PROC reg_save_list

        end
         