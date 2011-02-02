

;  Copyright (C) 2002,2004 Jens Franke, T.Kleinjung  
;  This file is part of gnfs4linux, distributed under the terms of the  
;  GNU General Public Licence and WITHOUT ANY WARRANTY.  
        
;  You should have received a copy of the GNU General Public License along  
;  with this program; see the file COPYING.  If not, write to the Free  
;  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
;  02111-1307, USA.  
        
;  Written by T. Kleinjung  
;  Modified by J. Franke  
;  Athlon64 version by J. Franke  
        
;  asm_evaluate(sievebegin,sieve.13,buffer,nmax):  
;  scans sievearray between sievebegin and sieve.13 for entries  
;  >127 and stores them in buffer (2 Bytes), stores at most nmax  
;  entries, returns number of stored entries  
        
;  edx: counts entries found so far  
;  esi: points to location of array we are investigating  
;  edi: points at .13 of array (=sieve.13)  
;  mm7: 0  
;  mm0-3:  
;  ebx, ecx:  

; asm_evaluate(unsigned long*, unsigned long*, unsignhed short*,  long);
; GAS    rax              rdi,            rsi,              rdx,   ecx
; YASM   rax              rcx,            rdx,               r8,   r9d
        
%include "ls-defs.inc"

        bits 64
        text
 
 %define reg_save_list  rsi, rdi

        FRAME_PROC asm_evaluate0, 0, reg_save_list  
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx, r8
        mov     ecx, r9d

        mov     r8, rdi
        mov     rax, rdx
        pxor    mm7, mm7
        lea     r9, [rax+rcx*2]
        jmp     .2   
.1:     lea     r8, [r8+32]
        movq    [r8-32], mm7
        movq    [r8-24], mm7
        movq    [r8-16], mm7
        movq    [r8-8], mm7
.2:     cmp     r8, rsi
        jz      .6   
        movq    mm0, [r8]
        movq    mm1, [r8+8]
        movq    mm2, [r8+16]
        movq    mm3, [r8+24]
        por     mm1, mm0
        por     mm3, mm2
        por     mm3, mm1
        pmovmskb r11, mm3
        test    r11, r11
        jz      .1
        movq    mm1, [r8+8]
        movq    mm3, [r8+24]
        pmovmskb r10, mm0
        pmovmskb r11, mm2
        sal     r11, 16
        pmovmskb rcx, mm1
        or      r10, r11
        pmovmskb r11, mm3
        sal     rcx, 8
        sal     r11, 24
        or      r10, rcx
        sub     r8, rdi
        or      r10, r11
        xor     r11, r11
.3:     bsf     rcx, r10
        add     r11, rcx
        inc     rcx
        add     r8, r11
        shr     r10, cl
        mov     [rax], r8w
        lea     rax, [rax+2]
        sub     r8, r11
        inc     r11
        cmp     r9, rax
        jbe     .4
        test    r10, r10
        jnz     .3
        add     r8, rdi
        jmp     .1    
.4:     add     r8, rdi
.5:     cmp     r8, rsi
        jz      .6   
        lea     r8, [r8+32]
        movq    [r8-32], mm7
        movq    [r8-24], mm7
        movq    [r8-16], mm7
        movq    [r8-8], mm7
        jmp     .5
.6:     sub     rax, rdx
        emms    
        shr     rax, 1
        END_PROC reg_save_list
        
        FRAME_PROC asm_evaluate, 0, reg_save_list
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx, r8
        mov     ecx, r9d
          
        mov     r8, rdi
        mov     rax, rdx
        lea     r9, [rax+rcx*2]
        sub     r8, 32        
.11:    lea     r8, [r8+32]
        cmp     r8, rsi
        jz      .13   
        movq    mm0, [r8]
        movq    mm1, [r8+8]
        movq    mm2, [r8+16]
        movq    mm3, [r8+24]
        por     mm1, mm0
        por     mm3, mm2
        por     mm3, mm1
        pmovmskb r11, mm3
        test    r11, r11
        jz      .11
        movq    mm1, [r8+8]
        movq    mm3, [r8+24]
        pmovmskb r10, mm0
        pmovmskb r11, mm2
        sal     r11, 16
        pmovmskb rcx, mm1
        or      r10, r11
        pmovmskb r11, mm3
        sal     rcx, 8
        sal     r11, 24
        or      r10, rcx
        sub     r8, rdi
        or      r10, r11
        xor     r11, r11
.12:    bsf     rcx, r10
        add     r11, rcx
        inc     rcx
        add     r8, r11
        shr     r10, cl
        mov     [rax], r8w
        lea     rax, [rax+2]
        sub     r8, r11
        inc     r11
        cmp     r9, rax
        jbe     .13
        test    r10, r10
        jnz     .12
        add     r8, rdi
        jmp     .11   
.13:    sub     rax, rdx
        emms    
        shr     rax, 1
        END_PROC reg_save_list
        
        FRAME_PROC asm_evaluate0_xmm, 0, reg_save_list  
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx, r8
        mov     ecx, r9d

        mov     r8, rdi
        mov     rax, rdx
        xorps   xmm7, xmm7
        lea     r9, [rax+rcx*2]
        jmp     .22    
.21:    lea     r8, [r8+32]
.22:    cmp     r8, rsi
        jz      .26   
        movaps  xmm0, [r8]
        movaps  [r8], xmm7
        movaps  xmm1, [r8+16]
        movaps  [r8+16], xmm7
        movaps  xmm2, xmm1
        orps    xmm1, xmm0
        pmovmskb r11, xmm1
        test    r11, r11
        jz      .21
        pmovmskb r10, xmm0
        pmovmskb r11, xmm2
        shl     r11, 16
        or      r10, r11
        sub     r8, rdi
        xor     r11, r11
.23:    bsf     rcx, r10
        add     r11, rcx
        inc     rcx
        add     r8, r11
        shr     r10, cl
        mov     [rax], r8w
        lea     rax, [rax+2]
        sub     r8, r11
        inc     r11
        cmp     r9, rax
        jbe     .24
        test    r10, r10
        jnz     .23
        add     r8, rdi
        jmp     .21       
.24:    add     r8, rdi
.25:    cmp     r8, rsi
        jz      .26  
        lea     r8, [r8+32]
        movaps  [r8-32], xmm7
        movaps  [r8-16], xmm7
        jmp     .25
.26:    sub     rax, rdx
        emms    
        shr     rax, 1
        END_PROC reg_save_list
        
        FRAME_PROC asm_evaluate_xmm, 0, reg_save_list  
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx, r8
        mov     ecx, r9d

        mov     r8, rdi
        mov     rax, rdx
        lea     r9, [rax+rcx*2]
        sub     r8, 32
.31:    lea     r8, [r8+32]
        cmp     r8, rsi
        jz      .33  
        movaps  xmm0, [r8]
        movaps  xmm1, [r8+16]
        movaps  xmm2, xmm0
        orps    xmm0, xmm1
        pmovmskb r11, xmm0
        test    r11, r11
        jz      .31
        pmovmskb r11, xmm1
        pmovmskb r10, xmm2
        shl     r11, 16
        or      r10, r11
        sub     r8, rdi
        xor     r11, r11
.32:    bsf     rcx, r10
        add     r11, rcx
        inc     rcx
        add     r8, r11
        shr     r10, cl
        mov     [rax], r8w
        lea     rax, [rax+2]
        sub     r8, r11
        inc     r11
        cmp     r9, rax
        jbe     .33
        test    r10, r10
        jnz     .32
        add     r8, rdi
        jmp     .31  
.33:    sub     rax, rdx
        emms    
        shr     rax, 1
        END_PROC reg_save_list

        end

