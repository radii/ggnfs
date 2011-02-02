
;  Now, define our arguments. We have 4 saved regs and the return address.  
;  Therefore, our first argument is:  

; mp_limb_t mpz_asm_td(mp_limb_t, mp_limb_t, mp_limb_t*, mp_size_t);
; GAS    rax                 rdi,       rsi,        rdx,        rcx
; YASM   rax                 rcx,       rdx,         r8,         r9
 

%include "ls-defs.inc"

        bits 64
        text
    
%define reg_save_list rsi, rdi

        FRAME_PROC mpz_asm_td, 0, reg_save_list 
        mov     rdi, rcx
        mov     rsi, rdx
        mov     rdx, r8
        mov     rcx, r9
 
        mov     r10, rdx
        xor     r9, r9
        lea     r8, [rdx+rcx*8]
        mov     rax, [rdx]
        cmp     rcx, 2
        lea     r8, [r8-8]
        jb      .2
.1:     mul     rsi
        mov     [r10], rax
        mul     rdi
        lea     r10, [r10+8]
        mov     rax, [r10]
        add     rdx, r9
        xor     r9, r9
        sub     rax, rdx
        adc     r9, r9
        cmp     r8, r10
        ja      .1
.2:     mul     rsi
        mov     [r10], rax
        mul     rdi
        add     rdx, r9
        mov     rax, rdx
        END_PROC reg_save_list

        end
