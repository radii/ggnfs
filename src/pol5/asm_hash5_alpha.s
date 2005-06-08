 #  Copyright (C) 2003 Jens Franke, T. Kleinjung.
 #  This file is part of gnfs4linux, distributed under the terms of the
 #  GNU General Public Licence and WITHOUT ANY WARRANTY.
 #
 #  You should have received a copy of the GNU General Public License along
 #  with this program; see the file COPYING.  If not, write to the Free
 #  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 #  02111-1307, USA.

 # knapsack-functions for HASHSHIFT=52

        .text

        .globl  s11len
        .comm   s11len,4
        .globl  s12len
        .comm   s12len,4
        .globl  s21len
        .comm   s21len,4
        .globl  s22len
        .comm   s22len,4
        .globl  s11l
        .comm   s11l,8
        .globl  s12l
        .comm   s12l,8
        .globl  s21l
        .comm   s21l,8
        .globl  s22l
        .comm   s22l,8
        .globl  hashdataptr
        .comm   hashdataptr,8
        .globl  raw_bound
        .comm   raw_bound,4
        .globl  raw_cand_ptr
        .comm   raw_cand_ptr,8

        .globl  modulo32
        .comm   modulo32,4

        .set noreorder
        .align 3
        .globl asm_modulo31
        .ent asm_modulo31
asm_modulo31:
        .frame $30,0,$26,0
        .prologue 0
	mulq $16,$17,$0
	ldl $1,modulo32
	remq $0,$1,$0
	ret     $31,($26),1
	.end    asm_modulo31

        .set noreorder
        .align 3
        .globl asm_getclock
        .ent asm_getclock
asm_getclock:
        .frame $30,0,$26,0
        .prologue 0
	rpcc $0
	stq $31,0($16)
	srl $0,32,$1
	addq $0,$1,$0
	stl $0,0($16)
	ret     $31,($26),1
	.end    asm_getclock


 #  memset(hashdata,0,NHASH*sizeof(uchar));
 #  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
 #  for (i=0; i<s11len; i++) {
 #    add=s11l[i];
 #    for (j=0; j<s12len; j++) {
 #      h=add+s12l[j];
 #      ind=h>>HASHSHIFT32;
 #      if (sort[ind]>=32) return 1;
 #      hash[sort[ind]*NHASH+ind]=h;
 #      sort[ind]++;
 #    }
 #  }
 #  return 0;

        .set noreorder
        .align 3
        .globl asm_hash1
        .ent asm_hash1
asm_hash1:    # hashes the values s11l[i]+s12l[j], 0<=i<s11len, 0<=j<s12len
        .frame $30,0,$26,0
        .prologue 0
 # memset(hashdata,0,NHASH*sizeof(uchar)):
	ldq $7,hashdataptr
	ldiq $1,256
	addq $31,$7,$8
zeroloop:
	subq $1,1,$1
	stq $31,0($8)
	stq $31,8($8)
	addq $8,16,$8
	bgt $1,zeroloop

	ldiq $16,1
	sll $16,32,$16
	subq $16,1,$16      # 0x00000000ffffffff in $16
	ldq $6,s11l
	ldl $1,s11len
	s4addq $1,$6,$5     # s11-end
	ldq $4,s12l
	ldl $1,s12len
	s4addq $1,$4,$3     # s12-end

outerloop1:
	cmpult $6,$5,$1
	beq $1,hash1_end    # if s11>=s11-end break
	ldq $4,s12l
	ldl $2,0($6)        # add=*s11
	addq $6,4,$6        # s11++

innerloop1:
	cmpult $4,$3,$1
	beq $1,outerloop1   # if s12>=s12-end break
	ldl $17,0($4)       # h=*s12
	addq $4,4,$4        # s12++
	addq $17,$2,$17     # h+=add
	and $17,$16,$17
	sra $17,20,$18      # ind=h>>20
	addq $18,$7,$19     # sort+ind
	ldbu $20,0($19)     # sort[ind]
	cmplt $20,32,$1     # sort[ind]<32 ?
	beq $1,return1_1
	sll $20,12,$21      # sort[ind]*NHASH
	addq $20,1,$20
	addq $21,$18,$21
	stb $20,0($19)      # sort[ind]++
	s4addq $21,$8,$21
	stl $17,0($21)      # hash[sort[ind]*NHASH+ind]=h
 # we use the sloppy variant which will miss approximately
 # raw_bound out of 2^20 solutions; since raw_bound is small (often 5)
 # this is much faster
        br innerloop1

hash1_end: # return 0
	mov 0,$0
	ret     $31,($26),1
return1_1: # return 1 on error
	mov 1,$0
	ret     $31,($26),1
	.end    asm_hash1



 #  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
 #  for (i=0; i<s21len; i++) {
 #    add=s21l[i];
 #    for (j=0; j<s22len; j++) {
 #      h=add+s22l[j];
 #      ind=h>>HASHSHIFT32;
 #      for (k=0; k<sort[ind]; k++) {
 #        if (hash[k*NHASH+ind]-h<raw_bound) {
 #          *raw_cand_ptr++=(i+(j<<16)); /* assumes npr_in_p<25 */
 #          k=sort[ind];  /* avoid duplicates */
 #        }
 #      }
 #    }
 #  }


        .set noreorder
        .align 3
        .globl asm_hash2
        .ent asm_hash2
asm_hash2:
        .frame $30,0,$26,0
        .prologue 0
	ldl $0,raw_bound
	ldq $7,hashdataptr
	ldiq $1,1
	sll $1,12,$1
	addq $1,$7,$8       # hashdata

	ldiq $16,1
	sll $16,32,$16
	subq $16,1,$16      # 0x00000000ffffffff in $16
	and $0,$16,$0
	ldq $6,s21l
	ldl $1,s21len
	s4addq $1,$6,$5     # s21-end
	ldq $4,s22l
	ldl $1,s22len
	s4addq $1,$4,$3     # s22-end

outerloop2:
	cmpult $6,$5,$1
	beq $1,hash2_end    # if s21>=s21-end break
	ldq $4,s22l
	ldl $2,0($6)        # add=*s21
	addq $6,4,$6        # s21++

innerloop2:
	cmpult $4,$3,$1
	beq $1,outerloop2   # if s22>=s22-end break
	ldl $17,0($4)       # h=*s22
	addq $4,4,$4        # s22++
	addq $17,$2,$17     # h+=add
	and $17,$16,$17
	sra $17,20,$18      # ind=h>>20
	addq $18,$7,$19     # sort+ind
	ldbu $20,0($19)     # sort[ind]
	beq $20,innerloop2
	s4addq $18,$8,$21   # hashdata+ind
	sll $20,12,$20      # NHASH*sort[ind]
	s4addq $20,$21,$20  # end
	mov 1,$19
	sll $19,12,$19      # NHASH

testloop:
	ldl $18,0($21)      # hash[k*NHASH+ind]
	subq $18,$17,$18    # hash[k*NHASH+ind]-h
	and $18,$16,$18
	cmpult $18,$0,$1
	beq $1,nostore      # if hash[k*NHASH+ind]-h>=raw_bound goto nostore
 # store i+(j<<16), probably not time critical
 # $17, ..., $21 free 
	ldq $17,s21l
	subq $6,$17,$17     # 4*(i+1)
	ldq $18,s22l
	subq $4,$18,$18     # 4*(j+1)
	subq $17,4,$17
	subq $18,4,$18
	sra $17,2,$17       # i
	sll $18,14,$18      # j<<16
	ldq $19,raw_cand_ptr
	addq $18,$17,$18    # i+(j<<16)
	stq $18,0($19)
	addq $19,4,$19
	stq $19,raw_cand_ptr
	br innerloop2       # avoid duplicates

nostore:
	s4addq $19,$21,$21
	cmpult $21,$20,$1
	beq $1,innerloop2
	br testloop

hash2_end:
        ret     $31,($26),1
        .end    asm_hash2

