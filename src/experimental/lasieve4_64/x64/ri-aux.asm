

;  Copyright (C) 2002 Jens Franke, T.Kleinjung  
;  This file is part of gnfs4linux, distributed under the terms of the  
;  GNU General Public Licence and WITHOUT ANY WARRANTY.  
        
;  You should have received a copy of the GNU General Public License along  
;  with this program; see the file COPYING.  If not, write to the Free  
;  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
;  02111-1307, USA.  
        
;  The function we want to write.  
;  ulong asm_getbc(x)  
;  Modular inverse of x modulo modulo32  
;  and y satisfies 0<y<x.  
            
;  Number of trial subtractions before doing a division  
        
; void asm_getbc(u32_t, u32_t, u32_t, u32_t*,   u32_t*,    u32_t*,  u32_t*);
; GAS    rax       edi,   esi,   edx,    rcx,       r8,        r9,  [rsp+8]
; YASM   rax       ecx,   edx,   r8d,     r9, [rsp+40],  [rsp+48], [rsp+56]

%include "ls-defs.inc"

        bits 64
        text
    
        extern  mpqs_256_inv_table

%define reg_save_list rsi, rdi, rbx

        FRAME_PROC asm_getbc, 0, reg_save_list
        mov     edi, ecx
        mov     esi, edx
        mov     edx, r8d
        mov     rcx,  r9
        mov      r8, [rsp+stack_use+40]
        mov      r9, [rsp+stack_use+48]

        xor     r10d, r10d
        mov     ebx, edx
        xor     r11d, r11d
        inc     r10d
        cmp     ebx, edi
        ja      .4
.1:  
%rep    15
        sub     esi, edi
        add     r11d, r10d
        cmp     esi, edi
        jb      .2
%endrep
        mov     eax, esi
        xor     edx, edx
        div     edi
        mov     esi, edx
        mul     r10d
        add     r11d, eax     
.2:     cmp     ebx, esi
        ja      .5
        
%rep    15
        sub     edi, esi
        add     r10d, r11d
        cmp     edi, esi
        jb      .3
%endrep
        mov     eax, edi
        xor     edx, edx
        div     esi
        mov     edi, edx
        mul     r11d
        add     r10d, eax
.3:     cmp     ebx, edi
        jbe     .1
.4: 
%rep    15
        sub     esi, edi
        add     r11d, r10d
        cmp     ebx, esi
        ja      .6
%endrep
        mov     eax, esi
        xor     edx, edx
        sub     eax, ebx
        div     edi
        inc     eax
        mov     ebx, eax
        mul     edi
        sub     esi, eax
        mov     eax, ebx
        mul     r10d
        add     r11d, eax
        jmp     .6
.5:
%rep    15
        sub     edi, esi
        add     r10d, r11d
        cmp     ebx, edi
        ja      .6
%endrep
        mov     eax, edi
        xor     edx, edx
        sub     eax, ebx
        div     esi
        inc     eax
        mov     ebx, eax
        mul     esi
        sub     edi, eax
        mov     eax, ebx
        mul     r11d
        add     r10d, eax
.6:     mov     rax,[rsp+stack_use+56]
        mov     [rcx], edi
        mov     [r8], r10d
        mov     [r9], esi
        mov     [rax], r11d
        END_PROC reg_save_list
        
        end
