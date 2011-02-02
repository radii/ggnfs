; Copyright (C) 2002 Jens Franke, T.Kleinjung  
; This file is part of gnfs4linux, distributed under the terms of the  
; GNU General Public Licence and WITHOUT ANY WARRANTY.  
        
; You should have received a copy of the GNU General Public License along  
; with this program; see the file COPYING.  If not, write to the Free  
; Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
; 02111-1307, USA.  

; void asm_gauss(void);

%include "ls-defs.inc"

        bits 64
        text
    
        extern  mpqs_gauss_c
        extern  mpqs_gauss_col
        extern  mpqs_gauss_d
        extern  mpqs_gauss_j
        extern  mpqs_gauss_k
        extern  mpqs_gauss_m
        extern  mpqs_gauss_mat
        extern  mpqs_gauss_n32
        extern  mpqs_gauss_row

%define reg_save_list rsi, rdi

        FRAME_PROC asm_gauss, 0, reg_save_list

.1:     mov     eax, [rip+mpqs_gauss_k]
        dec     eax
        mov     r8d, 1
        mov     [rip+mpqs_gauss_k], eax
        jl      .14
        mov     ecx, eax
        and     ecx, 31
        shl     r8d, cl
        mov     ecx, eax
        shr     ecx, 5    ; %r11d 
        mov     r11d, ecx      
        mov     ecx, [rip+mpqs_gauss_j]
        mov     rdi, mpqs_gauss_c
        
        mov     rsi, [rip+mpqs_gauss_row]
        mov     eax, r11d
        mov     edx, [rip+mpqs_gauss_n32]
        mov     rsi, [rsi+rcx*8]
        sub     rax, rdx
        dec     ecx
        lea     rsi, [rsi+rax*4]
        mov     eax, [rip+mpqs_gauss_n32]
        
        align   16
.2:     inc     ecx
.3:     cmp     ecx, [rip+mpqs_gauss_m]
        jnc     .15
        lea     rsi, [rsi+rax*4]
        mov     edx, [rsi]
        and     edx, r8d
        jz      .2        
.4:     cmp     ecx, [rip+mpqs_gauss_j]
        jz      .7
        mov     rsi, [rip+mpqs_gauss_row]
        mov     rax, [rsi+rcx*8]    ; mpqs_gauss_row[j] 
        mov     edx, [rip+mpqs_gauss_j]
        mov     r9, [rsi+rdx*8]    ; mpqs_gauss_row[mpqs_gauss_j] 
        xor     rdx, rdx
.5:     movq    mm0, [rax+rdx*4]
        movq    mm1, [r9+rdx*4]
        movq    [r9+rdx*4], mm0
        movq    [rax+rdx*4], mm1
        add     rdx, 2
        cmp     edx, [rip+mpqs_gauss_n32]
        jc      .5    
.6:     mov     ecx, [rip+mpqs_gauss_j]
.7:     mov     rsi, mpqs_gauss_d
        mov     eax, [rip+mpqs_gauss_k]
        mov     [rsi+rax*2], cx
        mov     [rdi+rcx*2], ax
        inc     dword[rip+mpqs_gauss_j]
        mov     r10d, ecx
        mov     rsi, [rip+mpqs_gauss_mat]
        mov     rdi, mpqs_gauss_col
        mov     edx, r11d
        test    ecx, ecx
        mov     eax, [rip+mpqs_gauss_n32]
        lea     rsi, [rsi+rdx*4]
        mov     rdx, 2
        jz      .9
        xor     ecx, ecx
        
        align   16
.8:     mov     r9d, [rsi]
        and     r9d, r8d
        mov     r9, 0
        mov     [rdi], cx
        cmovnz  r9, rdx
        inc     ecx
        add     rdi, r9
        cmp     ecx, r10d
        lea     rsi, [rsi+rax*4]
        jc      .8
.9:     inc     ecx
        cmp     ecx, [rip+mpqs_gauss_m]
        lea     rsi, [rsi+rax*4]
        jnc     .11
        
        align   16
.10:    mov     r9d, [rsi]
        and     r9d, r8d
        mov     r9, 0
        mov     [rdi], cx
        cmovnz  r9, rdx
        inc     rcx
        add     rdi, r9
        cmp     ecx, [rip+mpqs_gauss_m]
        lea     rsi, [rsi+rax*4]
        jc      .10
.11:    mov     r8, mpqs_gauss_col
        cmp     rdi, r8
        jz      .1
        mov     ecx, r10d
        mov     rsi, [rip+mpqs_gauss_row]
        mov     rsi, [rsi+rcx*8]   
        mov     edx, r11d
        cmp     edx, 2
        jc      .16
        cmp     edx, 4
        jc      .18
        cmp     edx, 6
        jc      .20

        align   16
.12:    movzx   r9, word[r8]
        mov     rax, [rip+mpqs_gauss_row]
        mov     rax, [rax+r9*8]
        xor     ecx, ecx
.13:    mov     edx, [rsi+rcx*4]
        xor     [rax+rcx*4], edx
        inc     ecx
        cmp     r11d, ecx
        jnc     .13  
        lea     r8, [r8+2]
        cmp     r8, rdi
        jc      .12
        jmp     .1
.14:    mov     eax, r10d
        emms    
        EXIT_PROC reg_save_list
        
.15:    mov     rdi, mpqs_gauss_d
        mov     eax, [rip+mpqs_gauss_k]
        mov     r9w, -1
        mov     [rdi+rax*2], r9w
        mov     r10d, ecx
        jmp     .14        
.16     movzx   r9, word[r8]
        movq    mm0, [rsi]
        mov     rsi, [rip+mpqs_gauss_row]
        mov     rax, [rsi+r9*8]
.17:    movq    mm1, [rax]
        lea     r8, [r8+2]
        movzx   r9, word[r8]
        pxor    mm1, mm0
        cmp     r8, rdi
        movq    [rax], mm1
        mov     rax, [rsi+r9*8]
        jc      .17
        jmp     .1
.18:    movq    mm0, [rsi]
        movzx   r9, word[r8]
        movq    mm2, [rsi+8]
        mov     rsi, [rip+mpqs_gauss_row]
        mov     rax, [rsi+r9*8]
.19:    movq    mm1, [rax]
        lea     r8, [r8+2]
        movq    mm3, [rax+8]
        pxor    mm1, mm0
        movzx   r9, word[r8]
        pxor    mm3, mm2
        cmp     r8, rdi
        movq    [rax], mm1
        movq    [rax+8], mm3
        mov     rax, [rsi+r9*8]
        jc      .19
        jmp     .1
.20     movq    mm0, [rsi]
        movzx   r9, word[r8]
        movq    mm2, [rsi+8]
        movq    mm4, [rsi+16]
        mov     rsi, [rip+mpqs_gauss_row]
        mov     rax, [rsi+r9*8]
.21:    movq    mm1, [rax]
        lea     r8, [r8+2]
        movq    mm3, [rax+8]
        pxor    mm1, mm0
        movq    mm5, [rax+16]
        movzx   r9, word[r8]
        pxor    mm3, mm2
        movq    [rax], mm1
        pxor    mm5, mm4
        movq    [rax+8], mm3
        cmp     r8, rdi
        movq    [rax+16], mm5
        mov     rax, [rsi+r9*8]
        jc      .21
        jmp     .1
        END_PROC reg_save_list

        end

        