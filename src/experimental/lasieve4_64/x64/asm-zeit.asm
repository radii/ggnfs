
%include "ls-defs.inc"

        bits 64
        text
    
        extern  asmzeitcounter

        LEAF_PROC asmgetclock
        rdtsc   
        shl     rdx, 32
        or      rax, rdx
        ret     
        
        LEAF_PROC zeitA
        rdtsc   
        shl     rdx, 32
        or      rax, rdx
        mov     rdx, [rip+asmzeitcounter]
        sub     [rdx+rdi*8], rax
        ret     
        
        LEAF_PROC zeitB
        rdtsc   
        shl     rdx, 32
        or      rax, rdx
        mov     rdx, [rip+asmzeitcounter]
        add     [rdx+rdi*8], rax
        ret     
    
        end
