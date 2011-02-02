
; Copyright (C) 2001 Jens Franke, T. Kleinjung.  
; This file is part of gnfs4linux, distributed under the terms of the  
; GNU General Public Licence and WITHOUT ANY WARRANTY.  
;         
; You should have received a copy of the GNU General Public License along  
; with this program; see the file COPYING.  If not, write to the Free  
; Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
; 02111-1307, USA.        

; void schedsieve(unsigned char*, u16_t, u16_t*, unsigned char*);
; GAS    rax                 rdi,   esi,    rdx,            rcx
; YASM   rax                 rcx,   edx,     r8,             r9

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi

        FRAME_PROC schedsieve, 0, reg_save_list 
        mov     rdi, rcx
        mov     esi, edx
        mov     rdx, r8
        mov     rcx, r9

        mov     eax, edi
        sub     rcx, 12
        cmp     rcx, rdx
        movzx   r8, word[rdx]
        movzx   r9, word[rdx+4]
        movzx   r10, word[rdx+8]
        movzx   r11, word[rdx+12]
        lea     rdx, [rdx+16]
        jbe     .2
.1:     prefetch [rdx+128]
        add     [rsi+r8], al
        movzx   r8, word[rdx]
        add     [rsi+r9], al
        movzx   r9, word[rdx+4]
        add     [rsi+r10], al
        movzx   r10, word[rdx+8]
        add     [rsi+r11], al
        cmp     rcx, rdx
        movzx   r11, word[rdx+12]
        lea     rdx, [rdx+16]
        ja      .1
.2:     add     rcx, 28
        cmp     rcx, rdx
        lea     rdx, [rdx+4]
        jbe     .3
        add     [rsi+r8], al
        cmp     rcx, rdx
        lea     rdx, [rdx+4]
        jbe     .3
        add     [rsi+r9], al
        cmp     rcx, rdx
        jbe     .3
        add     [rsi+r10], al
.3:     END_PROC reg_save_list

        end
