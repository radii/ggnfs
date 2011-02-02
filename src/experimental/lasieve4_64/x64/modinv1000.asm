
;  Copyright (C) 2004 Jens Franke, T.Kleinjung  
;  This file is part of gnfs4linux, distributed under the terms of the  
;  GNU General Public Licence and WITHOUT ANY WARRANTY.  
        
;  You should have received a copy of the GNU General Public License along  
;  with this program; see the file COPYING.  If not, write to the Free  
;  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
;  02111-1307, USA.  
        
;  The function we want to write.  
;  ulong asm_modinv32(x)  
;  Modular inverse of x modulo modulo32  
;  and y satisfies 0<y<x.  

;  Number of trial subtractions before doing a division

; unsigned long asm_modinv32b(u32_t)
; GAS    rax                    edi
; YASM   rax                    ecx

%include "ls-defs.inc"

        bits 64
        text
    
        extern  abort
        extern  modulo32

        LEAF_PROC asm_modinv32
        mov     r10d, ecx
        mov     r9d, [rel modulo32]
        test    r10d, r10d
        jz      .6
        cmp     r9d, r10d
        jbe     .6
        xor     r8d, r8d
        xor     ecx, ecx
        inc     ecx
        cmp     r10d, 1
        jbe     .4
.1:  
%rep    15
        sub     r9d, r10d
        add     r8d, ecx
        cmp     r9d, r10d
        jb      .2
%endrep      
        mov     eax, r9d
        xor     edx, edx
        div     r10d
        mov     r9d, edx
        mul     ecx
        add     r8d, eax  
.2:     cmp     r9d, 1
        jbe     .5
%rep    15
        sub     r10d, r9d
        add     ecx, r8d
        cmp     r10d, r9d
        jb      .3
%endrep
        mov     eax, r10d
        xor     edx, edx
        div     r9d
        mov     r10d, edx
        mul     r8d
        add     ecx, eax
.3:     cmp     r10d, 1
        ja      .1
.4:     jne     .6
        mov     eax, ecx
        ret
.5:     jne     .6
        mov     eax, [rel modulo32]
        sub     eax, r8d
        ret
.6:     call    abort
 
        end
           