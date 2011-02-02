
;  Copyright (C) 2002 Jens Franke, T.Kleinjung  
;  This file is part of gnfs4linux, distributed under the terms of the  
;  GNU General Public Licence and WITHOUT ANY WARRANTY.  
        
;  You should have received a copy of the GNU General Public License along  
;  with this program; see the file COPYING.  If not, write to the Free  
;  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  
;  02111-1307, USA.  
             
;  for (i=0; i<mpqs_nFB; i++) {  
;    p=fb[2*i];  
;    mmi=mpqs_FB_mm_inv[i];  
;    pi=fbs[2*i]; bbb=fbs[2*i+1];  
;    cc=bbb;  
;    if (cc&1) cc+=p; cc>>=1;  
;    cc=mpqs_FB_disp[i]+(p-cc); if (cc>=p) cc-=p;  
;    cc1=fb[2*i+1];  
;    h32=cc1*pi;  
;    MMREDUCE; cc1=h32;  
;    if (cc1>=p) cc1-=p;  
;    cc2=p-cc1;  
;    cc1+=cc; if (cc1>=p) cc1-=p;  
;    cc2+=cc; if (cc2>=p) cc2-=p;  
;    fbs[2*i]=(ushort)cc1; fbs[2*i+1]=(ushort)cc2;  
;  }    

%include "ls-defs.inc"

        bits 64
        text
    
        extern  mpqs_FB_disp
        extern  mpqs_FB_mm_inv
        extern  mpqs_FB_np_p
        extern  mpqs_FB_np_px
        extern  mpqs_FB_start

        LEAF_PROC asm_next_pol11_xmm
;       mov     rcx, rdi 
        mov     rax, 0x00010001
        mov     r9, rax
        shl     rax, 32
        or      rax, r9
        movd    xmm7, rax
        pslldq  xmm7, 8
        movd    xmm6, rax
        paddw   xmm7, xmm6    ; 0x00010001000100010001000100010001 
        mov     r8, mpqs_FB_disp
        mov     r9, mpqs_FB_np_px
        mov     r10, mpqs_FB_mm_inv
        mov     r11, mpqs_FB_start
.1:     movaps  xmm0, [r9]    ; p 
        movaps  xmm1, [r9+16]    ; sqrt 
        movaps  xmm2, [r10]    ; mm_inv 
        movaps  xmm3, [r11]    ; pi 
        movaps  xmm4, [r11+16]    ; cc         
        movaps  xmm6, xmm3
        pmullw  xmm6, xmm1    ; low sqrt*cc1 
        pmulhuw xmm3, xmm1    ; high sqrt*cc1 
        pmullw  xmm6, xmm2    ; h=low(sqrt*cc1)*mm_inv 
        pxor    xmm5, xmm5
        pcmpeqw xmm5, xmm6    ; 0xffff iff h=0 
        pand    xmm5, xmm7    ; 0x0001 iff h=0 
        pmulhuw xmm6, xmm0    ; high h*p 
        paddw   xmm3, xmm7
        psubw   xmm3, xmm5    ; carry 
        paddw   xmm3, xmm6    ; res 
        movaps  xmm6, xmm3    ; if >=p subtract p 
        paddw   xmm6, xmm7    ; res+1 
        pcmpgtw xmm6, xmm0    ; 0xffff iff res>=p 
        pand    xmm6, xmm0    ; p iff res>=p 
        psubw   xmm3, xmm6    ; res mod p = cc1     
        movaps  xmm5, xmm4
        pand    xmm5, xmm7    ; 0x0001 iff cc odd 
        pcmpeqw xmm5, xmm7    ; 0xffff iff cc odd 
        pand    xmm5, xmm0    ; p iff cc odd 
        paddw   xmm5, xmm4
        psrlw   xmm5, 1    ; cc/2 mod p = cc 
        movaps  xmm1, [r8]    ; disp 
        movaps  xmm4, xmm0
        psubw   xmm4, xmm5    ; p-cc 
        paddw   xmm4, xmm1
        movaps  xmm6, xmm4    ; if >=p subtract p 
        paddw   xmm6, xmm7    ; res+1 
        pcmpgtw xmm6, xmm0    ; 0xffff iff res>=p 
        pand    xmm6, xmm0    ; p iff res>=p 
        psubw   xmm4, xmm6    ; res mod p = cc 
        movaps  xmm2, xmm0
        psubw   xmm2, xmm3    ; p-cc1 = cc2 
        paddw   xmm3, xmm4
        movaps  xmm6, xmm3    ; if >=p subtract p 
        paddw   xmm6, xmm7    ; res+1 
        pcmpgtw xmm6, xmm0    ; 0xffff iff res>=p 
        pand    xmm6, xmm0    ; p iff res>=p 
        psubw   xmm3, xmm6    ; res mod p = cc1 
        paddw   xmm2, xmm4
        movaps  xmm5, xmm2    ; if >=p subtract p 
        paddw   xmm5, xmm7    ; res+1 
        pcmpgtw xmm5, xmm0    ; 0xffff iff res>=p 
        pand    xmm5, xmm0    ; p iff res>=p 
        psubw   xmm2, xmm5    ; res mod p = cc2 
        movaps  xmm4, xmm3
        punpcklwd xmm3, xmm2
        punpckhwd xmm4, xmm2
        movaps  [r11], xmm3
        movaps  [r11+16], xmm4
        dec     rcx
        lea     r11, [r11+32]
        lea     r8, [r8+16]
        lea     r9, [r9+32]
        lea     r10, [r10+16]
        jnz     .1
        emms    
        ret     
        
;  for (i=1; i<mpqs_nFB; i++) {  
;    p=fb[2*i];  
;    mmi=mpqs_FB_mm_inv[i];  
;    cc=invptr[i];  
;    cc*=bimul; h32=cc;  
;    MMREDUCE; cc=h32; if (cc>=p) cc-=p;  
;    ropptr[i]=(ushort)cc;  
;    pi=fbs[2*i];  
;    pi*=invptr[i]; h32=pi;  
;    MMREDUCE; fbs[2*i]=(ushort)h32;  
;    bbb=fbs[2*i+1]+cc; if (bbb>=p) bbb-=p; fbs[2*i+1]=bbb;  
;  }    

;       asm_next_pol10_xmm(%rcx,*invptr,*%r8,bimul) 
 
        LEAF_PROC asm_next_pol10_xmm
;       mov     r8, rdx             ; p3 ok
        mov     rax, 0x00010001
        mov     rdx, rax
        shl     rax, 32
        or      rax, rdx
        movd    xmm7, rax
        pslldq  xmm7, 8
        movd    xmm6, rax
        paddw   xmm7, xmm6          ; 0x00010001000100010001000100010001 
        mul     r9                  ; p4 ***
        movd    xmm1, rax
        pslldq  xmm1, 8
        movd    xmm5, rax
        paddw   xmm1, xmm5     
;       mov     rcx, rdi            ; p1 ok
        mov     r9, mpqs_FB_np_px
        mov     r10, mpqs_FB_mm_inv
        mov     r11, mpqs_FB_start
;       mov     rdx, rsi            ; p2 ok
.2:     movaps  xmm0, [r9]    ; p 
        movaps  xmm2, [r10]    ; mm_inv 
        movaps  xmm3, [r11]    ; pi 
        movaps  xmm4, [rdx]    ; cc=invptr[i] 
        movaps  xmm6, xmm3
        pmullw  xmm6, xmm4    ; low pi*cc 
        pmulhuw xmm3, xmm4    ; high pi*cc 
        pmullw  xmm6, xmm2    ; h=low(pi*cc)*mm_inv 
        pxor    xmm5, xmm5
        pcmpeqw xmm5, xmm6    ; 0xffff iff h=0 
        pand    xmm5, xmm7    ; 0x0001 iff h=0 
        pmulhuw xmm6, xmm0    ; high h*p 
        paddw   xmm3, xmm7
        psubw   xmm3, xmm5    ; carry 
        paddw   xmm3, xmm6    ; res 
        movaps  [r11], xmm3
        movaps  xmm6, xmm4
        pmullw  xmm6, xmm1    ; low pi*cc 
        pmulhuw xmm4, xmm1    ; high pi*cc 
        pmullw  xmm6, xmm2    ; h=low(pi*cc)*mm_inv 
        pxor    xmm5, xmm5
        pcmpeqw xmm5, xmm6    ; 0xffff iff h=0 
        pand    xmm5, xmm7    ; 0x0001 iff h=0 
        pmulhuw xmm6, xmm0    ; high h*p 
        paddw   xmm4, xmm7
        psubw   xmm4, xmm5    ; carry 
        paddw   xmm4, xmm6    ; res 
        movaps  xmm6, xmm4    ; if >=p subtract p 
        paddw   xmm6, xmm7    ; res+1 
        pcmpgtw xmm6, xmm0    ; 0xffff iff res>=p 
        pand    xmm6, xmm0    ; p iff res>=p 
        psubw   xmm4, xmm6    ; res mod p 
        movaps  [r8], xmm4
        movaps  xmm3, [r11+16]    ; fbs[2*i+1] 
        paddw   xmm3, xmm4    ; fbs[2*i+1]+cc 
        movaps  xmm6, xmm3    ; if >=p subtract p 
        paddw   xmm6, xmm7    ; res+1 
        pcmpgtw xmm6, xmm0    ; 0xffff iff res>=p 
        pand    xmm6, xmm0    ; p iff res>=p 
        psubw   xmm3, xmm6    ; res mod p 
        movaps  [r11+16], xmm3
        dec     rcx
        lea     r11, [r11+32]
        lea     r8, [r8+16]
        lea     rdx, [rdx+16]
        lea     r9, [r9+32]
        lea     r10, [r10+16]
        jnz     .2
        emms    
        ret     
        
;  asm_next_pol3plus(len,*SI_add)  

        LEAF_PROC asm_next_pol3plus
        mov    r8d, ecx
        mov     r9, rdx

        mov     rcx, mpqs_FB_start
        mov     rdx, mpqs_FB_np_p
        mov     rax, 0x00010001
        movd    mm7, rax
        psllq   mm7, 32
        movd    mm6, rax
        paddw   mm7, mm6    ; 0x0001000100010001       
.3:     movq    mm2, [rdx]
        movq    mm4, [r9]
        movq    mm0, mm2
        movq    mm1, mm2
        punpcklwd mm0, mm2    ; p0 p0 p1 p1 
        punpckhwd mm1, mm2    ; p2 p2 p3 p3 
        movq    mm2, mm4
        punpcklwd mm2, mm4    ; a0 a0 a1 a1 
        movq    mm3, mm4
        punpckhwd mm3, mm4    ; a2 a2 a3 a3         
        movq    mm4, [rcx]
        paddw   mm4, mm2
        movq    mm2, mm7
        paddw   mm2, mm4
        pcmpgtw mm2, mm0
        pand    mm2, mm0
        psubw   mm4, mm2
        movq    [rcx], mm4
        movq    mm5, [rcx+8]
        paddw   mm5, mm3
        movq    mm3, mm7
        paddw   mm3, mm5
        pcmpgtw mm3, mm1
        pand    mm3, mm1
        psubw   mm5, mm3
        movq    [rcx+8], mm5       
        lea     rdx, [rdx+16]
        dec     r8
        lea     r9, [r9+8]
        lea     rcx, [rcx+16]
        jnz     .3
        emms    
        ret     
        
;  asm_next_pol3plus_xmm(len,*SI_add)  

        LEAF_PROC asm_next_pol3plus_xmm
        mov    r8d, ecx
        mov     r9, rdx

        mov     rcx, mpqs_FB_start
        mov     rdx, mpqs_FB_np_px
        mov     rax, 0x00010001
        mov     r9, rax
        shl     rax, 32
        or      rax, r9
        movd    xmm7, rax
        pslldq  xmm7, 8
        movd    xmm6, rax
        paddw   xmm7, xmm6    ; 0x00010001000100010001000100010001  
.4:     movaps  xmm2, [rdx]
        movaps  xmm4, [r9]
        movaps  xmm0, xmm2
        movaps  xmm1, xmm2
        punpcklwd xmm0, xmm2    ; p0 p0 p1 p1 
        punpckhwd xmm1, xmm2    ; p2 p2 p3 p3 
        movaps  xmm2, xmm4
        punpcklwd xmm2, xmm4    ; a0 a0 a1 a1 
        movaps  xmm3, xmm4
        punpckhwd xmm3, xmm4    ; a2 a2 a3 a3 
        movaps  xmm4, [rcx]
        paddw   xmm4, xmm2
        movaps  xmm2, xmm7
        paddw   xmm2, xmm4
        pcmpgtw xmm2, xmm0
        pand    xmm2, xmm0
        psubw   xmm4, xmm2
        movaps  [rcx], xmm4
        movaps  xmm5, [rcx+16]
        paddw   xmm5, xmm3
        movaps  xmm3, xmm7
        paddw   xmm3, xmm5
        pcmpgtw xmm3, xmm1
        pand    xmm3, xmm1
        psubw   xmm5, xmm3
        movaps  [rcx+16], xmm5
        lea     rdx, [rdx+32]
        dec     r8
        lea     r9, [r9+16]
        lea     rcx, [rcx+32]
        jnz     .4
        emms    
        ret     
        
        
;  asm_next_pol3minus(len,*SI_add)  

        LEAF_PROC asm_next_pol3minus
        mov    r8d, ecx
        mov     r9, rdx

        mov     rcx, mpqs_FB_start
        mov     rdx, mpqs_FB_np_p
        mov     rax, 0x00010001
        movd    mm7, rax
        psllq   mm7, 32
        movd    mm6, rax
        paddw   mm7, mm6    ; 0x0001000100010001    
.5:     movq    mm2, [rdx]
        movq    mm4, [r9]
        movq    mm0, mm2
        movq    mm1, mm2
        punpcklwd mm0, mm2    ; p0 p0 p1 p1 
        punpckhwd mm1, mm2    ; p2 p2 p3 p3 
        movq    mm2, mm4
        punpcklwd mm2, mm4    ; a0 a0 a1 a1 
        movq    mm3, mm4
        punpckhwd mm3, mm4    ; a2 a2 a3 a3      
        movq    mm4, [rcx]
        psubw   mm2, mm0
        psubw   mm4, mm2
        movq    mm2, mm7
        paddw   mm2, mm4
        pcmpgtw mm2, mm0
        pand    mm2, mm0
        psubw   mm4, mm2
        movq    [rcx], mm4
        movq    mm5, [rcx+8]
        psubw   mm3, mm1
        psubw   mm5, mm3
        movq    mm3, mm7
        paddw   mm3, mm5
        pcmpgtw mm3, mm1
        pand    mm3, mm1
        psubw   mm5, mm3
        movq    [rcx+8], mm5
        lea     rdx, [rdx+16]
        dec     r8
        lea     r9, [r9+8]
        lea     rcx, [rcx+16]
        jnz     .5
        emms    
        ret     
        
;  asm_next_pol3minus_xmm(len,*SI_add)  

        LEAF_PROC asm_next_pol3minus_xmm
        mov    r8d, ecx
        mov     r9, rdx

        mov     rcx, mpqs_FB_start
        mov     rdx, mpqs_FB_np_px
        mov     rax, 0x00010001
        mov     r9, rax
        shl     rax, 32
        or      rax, r9
        movd    xmm7, rax
        pslldq  xmm7, 8
        movd    xmm6, rax
        paddw   xmm7, xmm6    ; 0x00010001000100010001000100010001 
.6:     movaps  xmm2, [rdx]
        movaps  xmm4, [r9]
        movaps  xmm0, xmm2
        movaps  xmm1, xmm2
        punpcklwd xmm0, xmm2    ; p0 p0 p1 p1 
        punpckhwd xmm1, xmm2    ; p2 p2 p3 p3 
        movaps  xmm2, xmm4
        punpcklwd xmm2, xmm4    ; a0 a0 a1 a1 
        movaps  xmm3, xmm4
        punpckhwd xmm3, xmm4    ; a2 a2 a3 a3     
        movaps  xmm4, [rcx]
        psubw   xmm2, xmm0
        psubw   xmm4, xmm2
        movaps  xmm2, xmm7
        paddw   xmm2, xmm4
        pcmpgtw xmm2, xmm0
        pand    xmm2, xmm0
        psubw   xmm4, xmm2
        movaps  [rcx], xmm4
        movaps  xmm5, [rcx+16]
        psubw   xmm3, xmm1
        psubw   xmm5, xmm3
        movaps  xmm3, xmm7
        paddw   xmm3, xmm5
        pcmpgtw xmm3, xmm1
        pand    xmm3, xmm1
        psubw   xmm5, xmm3
        movaps  [rcx+16], xmm5
        lea     rdx, [rdx+32]
        dec     r8
        lea     r9, [r9+16]
        lea     rcx, [rcx+32]
        jnz     .6
        emms    
        ret     
        
        end
