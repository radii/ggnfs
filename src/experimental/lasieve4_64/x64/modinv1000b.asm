
;  Copyright (C) 2004 Jens Franke, T.Kleinjung  
;  This file is part of gnfs4linux, distributed under the terms of the  
;  GNU General Public Licence and WITHOUT ANY WARRANTY.  
        
;  You should have received a copy of the GNU General Public License along  
;  with this program; see the file COPYING.  If not, write to the Free  
;  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
;  02111-1307, USA.  
        
;  The function we want to write.  
;  ulong asm_modinv32b(x,y)  
;  Modular inverse of x modulo y  
;  when x satisfies 0<x<y.  
;  If x==0 or x==y, 0 is returned  
;  All args outside these bounds,  
;  and all other non-coprime args, cause abort  
        
;  Number of trial subtractions before doing a division  

; unsigned long asm_modinv32b(u32_t, u32_t)
; GAS    rax                    edi,   esi
; YASM   rax                    ecx,   edx

%include "ls-defs.inc"

        bits 64
        text

        extern  abort
        
        LEAF_PROC asm_modinv32b
        mov     r10d, ecx
        mov     r11d, edx
          
        test    r10d, r10d
        mov     r9d, r11d
        jz      .8
        cmp     r11d, r10d
        jbe     .6
        xor     r8d, r8d
        xor     ecx, ecx
        inc     ecx
        cmp     r10d, 1
        jbe     .4
.1:  
%rep    15
        sub     r11d, r10d
        add     r8d, ecx
        cmp     r11d, r10d
        jb      .2
%endrep      
        mov     eax, r11d
        xor     edx, edx
        div     r10d
        mov     r11d, edx
        mul     ecx
        add     r8d, eax      
.2:     cmp     r11d, 1
        jbe     .5
%rep    15
        sub     r10d, r11d
        add     ecx, r8d
        cmp     r10d, r11d
        jb      .3
%endrep
        mov     eax, r10d
        xor     edx, edx
        div     r11d
        mov     r10d, edx
        mul     r8d
        add     ecx, eax
.3:     cmp     r10d, 1
        ja      .1
.4:     jne     .7
        mov     eax, ecx
        ret         
.5:     jne     .7
        mov     eax, r9d
        sub     eax, r8d
        ret     
.6:     je      .8
.7:     call    abort
.8:     xor     eax, eax
        ret     

        end
   