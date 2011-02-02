
;  Copyright (C) 2004 Jens Franke, T.Kleinjung  
;  This file is part of gnfs4linux, distributed under the terms of the  
;  GNU General Public Licence and WITHOUT ANY WARRANTY.  
        
;  You should have received a copy of the GNU General Public License along  
;  with this program; see the file COPYING.  If not, write to the Free  
;  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
;  02111-1307, USA.  
                
;  Hilfsvariablen  
        
;  The pair of factor base primes  
;  The pair of corresponding projective roots  
;  Inverse of the pair of factor base primes  
;  hat a0 im ersten und dritten dword zu stehen  
;  Same for a1, etc  
        
;  Use this comment for lines contributing to the calculation of x  
;  %eax and %edx are used in the division.  
        
;  In general we have to calculate  
;   ( +- a_1 + rb_1 ) / (+- a_0 - rb_0 ) modulo p.  
;  The signs of a_0 and a_1 depend on the lattice; for each choice we have  
;  a function asm_lasieve_mm_setup_mp_ , i=0,1,2,3.  
;       
;  The calculation is done as follows:  
;  Let p be fixed and R(a)=a/2^32 mod p.  
;  We first compute num=R(+- a_1 + rb_1) and den=R(R(+- a_0 - rb_0)).  
;  For den!=0 the result is R(num*den^-1).  
;  If two successive prime ideals of the factor base lie over the same prime p  
;  we try to save one inversion using a trick of Montgomery (rarely one of the  
;  denominators is zero; in this case we do the two calculations seperately).  
;  R(a) is calculated as in Montgomery multiplication.  
        
        
;  case a0>=0, a1>=0:	i=0  
;  case a0>=0, a1<0:	i=1  
;  case a0<0, a1>=0:	i=2  
;  case a0<0, a1<0:	i=3  
;  asm_lasieve_mm_setup_mp_(FB,proots,fbsz,absa0,b0_ul,absa1,b1_ul,ri_ptr)  

; u32_t *function(u32_t*, u32_t*, size_t, u32_t,    u32_t,    u32_t,    u32_t,   u32_t*)
; GAS    rax         rdi,    rsi,    rdx,   rcx,       r8,       r9,  [rsp+8], [rsp+16]
; YASM   rax         rcx,    rdx,     r8,    r9, [rsp+40], [rsp+48], [rsp+56], [rsp+64] 

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi, rbx, rbp, r12, r13, r14, r15

        extern  asm_modinv32b
        extern  get_recurrence_info
        extern  mpqs_256_inv_table

%macro expand 1
        FRAME_PROC asm_lasieve_mm_setup%1, 0, reg_save_list
        mov     rdi, rcx
        mov     rsi, rdx
        mov     edx, r8d
        mov     ecx, r9d
        mov     r8d, dword [rsp+stack_use+40]
        mov     r9d, dword [rsp+stack_use+48]

        pxor    xmm3, xmm3
        shr     rdx, 1
        mov     r12, rdi
        pxor    xmm5, xmm5
        pxor    xmm4, xmm4
        adc     rdx, 0
        pxor    xmm6, xmm6
        movd    xmm3, rcx
        mov     rbp, rsi
        movd    xmm5, r9
        movd    xmm4, r8
        lea     r14, [r12+rdx*8]
        movd    xmm6, dword [rsp+stack_use+56]
        pshufd  xmm3, xmm3, 0x88
        pshufd  xmm5, xmm5, 0x88
        mov     r13, [rsp+stack_use+64]
        pshufd  xmm4, xmm4, 0x88
        pshufd  xmm6, xmm6, 0x88
        mov     rax, [r12]
        mov     rdx, rax
        movd    xmm0, rax
        movq    xmm1, [rbp]
        shr     rdx, 32
        and     rax, 0x7fffffff
        shr     rax, 1
        shr     rdx, 1
        and     rax, 0x7f
        and     rdx, 0x7f
        movzx   rax, byte[rax+mpqs_256_inv_table]
        movzx   rdx, byte[rdx+mpqs_256_inv_table]
        pxor    xmm2, xmm2
        pinsrw  xmm2, rax, 0
        pinsrw  xmm2, rdx, 4
        pshufd  xmm0, xmm0, 0x98
        movdqa  xmm11, xmm2
        pmuludq xmm2, xmm0
        pmuludq xmm2, xmm11
        pslld   xmm11, 1
        psubd   xmm11, xmm2
        movdqa  xmm2, xmm11
        pmuludq xmm11, xmm0
        pmuludq xmm11, xmm2
        pshufd  xmm12, xmm1, 0x98
        pslld   xmm2, 1
%%1:    movdqa  xmm8, xmm5
        mov     ebx, [r12]
        mov     r15d, [r12+4]    
        movdqa  xmm1, xmm12
        lea     r12, [r12+8]
        pcmpeqd xmm12, xmm0
%if (%1 - 1) * (%1 - 2) == 0  
        pcmpeqd xmm13, xmm13
%endif         
        lea     rbp, [rbp+8]
        pshufd  xmm12, xmm12, 0xa0        
        psubd   xmm2, xmm11
%if (%1 - 1) * (%1 - 2) == 0  
        paddd   xmm13, xmm0
%endif  
%if %1 * (%1 - 1) == 0  
        movdqa  xmm7, xmm0
%else   
        movdqa  xmm7, xmm1
%endif  
        
%if (%1 - 1) * (%1 - 2) == 0  
        pmuludq xmm8, xmm13
%endif  
        
%if %1 * (%1 - 2) == 0  
        movdqa  xmm9, xmm0
%else   
        movdqa  xmm9, xmm1
%endif  
        pand    xmm8, xmm12
        movdqa  xmm10, xmm6
%if %1 * (%1 - 1) == 0  
        psubd   xmm7, xmm1
%endif  
        
%if %1 * (%1 - 2) == 0  
        psubd   xmm9, xmm1
%endif  
        pmuludq xmm7, xmm5
        movdqa  xmm11, xmm12
        paddq   xmm7, xmm3
        pand    xmm10, xmm12
        pandn   xmm11, xmm7
        movq    xmm1, [rbp]
        por     xmm11, xmm8
        movdqa  xmm8, xmm11
        pmuludq xmm11, xmm2
        pxor    xmm7, xmm7
        pmuludq xmm11, xmm0
        psubq   xmm8, xmm11
        pmuludq xmm9, xmm6
        psrldq  xmm8, 4
        pcmpgtd xmm7, xmm8
        pand    xmm7, xmm0
        paddd   xmm8, xmm7        
        movdqa  xmm7, xmm8
        pmuludq xmm8, xmm2
        paddq   xmm9, xmm4
        pmuludq xmm8, xmm0
        pandn   xmm12, xmm9
        psubq   xmm7, xmm8
        psrldq  xmm7, 4
        por     xmm12, xmm10
        paddd   xmm7, xmm0
%%2:    movd    edi, xmm7
        movdqa  xmm9, xmm12
        psrldq  xmm7, 8
        mov     esi, ebx
        pmuludq xmm12, xmm2
        call    asm_modinv32b
        movd    edi, xmm7
        pmuludq xmm12, xmm0
        movd    xmm11, eax
        mov     esi, r15d
        pxor    xmm10, xmm10
        call    asm_modinv32b
        movd    xmm8, eax
        psubq   xmm9, xmm12
        mov     rax, [r12]
        pxor    xmm12, xmm12
        pshufd  xmm8, xmm8, 0x8a
        psrldq  xmm9, 4
        mov     rdx, rax
%if %1 * (%1 - 3) == 0  
        movdqa  xmm7, xmm0
        pcmpgtd xmm10, xmm9
        por     xmm8, xmm11
        pcmpeqd xmm11, xmm11
        pand    xmm10, xmm0
        pcmpeqd xmm12, xmm8
        psubd   xmm7, xmm8
%else   
        movdqa  xmm7, xmm11
        pcmpgtd xmm10, xmm9
        por     xmm7, xmm8
        pand    xmm10, xmm0
        pcmpeqd xmm11, xmm11
        pcmpeqd xmm12, xmm7
%endif  
        paddd   xmm9, xmm10
        pxor    xmm10, xmm10        
        shr     rdx, 33
        shr     rax, 1
        pmuludq xmm7, xmm9
        and     rax, 0x7f
        and     rdx, 0x7f
        movdqa  xmm11, xmm7
        pshufd  xmm8, xmm8, 0x98
        mov     rdi, r13
        pmuludq xmm11, xmm2
        pxor    xmm2, xmm2
        movzx   rax, byte[rax+mpqs_256_inv_table]
        movzx   rdx, byte[rdx+mpqs_256_inv_table]
        pmuludq xmm11, xmm0
        pinsrw  xmm2, rax, 0
        pinsrw  xmm2, rdx, 4
        psubq   xmm7, xmm11
        mov     esi, ebx
        movdqa  xmm11, xmm2
        psrldq  xmm7, 4
        pmuludq xmm2, xmm8
        pcmpgtd xmm10, xmm7
        pand    xmm10, xmm0
        paddd   xmm7, xmm10
        pmuludq xmm2, xmm11
        pslld   xmm11, 1
        pand    xmm7, xmm12
        pandn   xmm12, xmm0
        psubd   xmm11, xmm2
        por     xmm7, xmm12        
%%3:    movd    edx, xmm7
        psrldq  xmm7, 8
        movdqa  xmm2, xmm11
        pmuludq xmm11, xmm8
        call    get_recurrence_info
        mov     esi, r15d
        lea     rdi, [r13+8]
        movd    edx, xmm7
        pmuludq xmm11, xmm2
        movdqa  xmm0, xmm8
        call    get_recurrence_info
        cmp     r14, r12
        pslld   xmm2, 1
        pshufd  xmm12, xmm1, 0x98
        lea     r13, [r13+16]
        ja      %%1
        END_PROC reg_save_list
%endmacro

        expand 0
        expand 1
        expand 2
        expand 3

        end
