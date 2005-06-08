
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
        .globl  s11_begin
        .comm   s11_begin,8
        .globl  hashpart_shift
        .comm   hashpart_shift,4
        .globl  hash_shift
        .comm   hash_shift,4



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



 # asm_hash1(hb): hashes the values s11l[i]+s12l[j], 0<=i<s11len, 0<=j<s12len
 #                for which (s11l[i]+s12l[j])>>hashpart_shift=hb
 #
 #  mask=(1<<hashpart_shift)-1; mask=~mask;
 #  hsub=hb<<hashpart_shift;
 #  memset(hashdata,0,NHASH*sizeof(uchar));
 #  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
 #  for (i=0; i<s11len; i++) {
 #    add=s11l[i];
 #    j=s11_begin[i];
 #    while (1) {
 #      h=add+s12l_sort[j];
 #      hp=h&mask;
 #      if (hp!=hsub) break;
 #      ind=(h-hsub)>>hash_shift;
 #      if (sort[ind]>=32) return 1;
 #      hash[sort[ind]*NHASH+ind]=h;
 #      sort[ind]++;
 #      j++;
 #    }
 #    if (j>=s12len) j-=s12len;
 #    s11_begin[i]=j;
 #  }
 #  return 0;

 # register usage
 # 0   sort[ind] et al
 # 1
 # 2   s12l_sort+j
 # 3   s12len
 # 4   s11_begin+i
 # 5   s11l+s11len
 # 6   s11l+i
 # 7   sort
 # 8   hash
 # 9   add
 # 10  s12l_sort
 # 11  hash_shift
 # 16  32 bit mask
 # 17  h
 # 18
 # 19  ind
 # 20  hsub
 # 21  mask


        .set noreorder
        .align 3
        .globl asm_hash1
        .ent asm_hash1
asm_hash1:    # hashes the values s11l[i]+s12l[j], 0<=i<s11len, 0<=j<s12len
        .frame $30,0,$26,0
        .prologue 0
	lda $30,-80($30)
	stq $9,8($30)
	stq $10,16($30)
	stq $11,24($30)
 # memset(hashdata,0,NHASH*sizeof(uchar)):
	ldq $7,hashdataptr
	ldiq $1,256         # 256=2^(64-HASHSHIFT)/16
	addq $31,$7,$8
zeroloop:
	subq $1,1,$1
	stq $31,0($8)
	stq $31,8($8)
	addq $8,16,$8
	bgt $1,zeroloop

	ldl $1,hashpart_shift
	sll $16,$1,$20      # hsub=hb<<hashpart_shift
	ldiq $21,1
	sll $21,$1,$21      # 2^hashpart_shift
	negq $21,$21        # mask
	ldiq $16,1
	sll $16,32,$16
	subq $16,1,$16      # 0x00000000ffffffff in $16
	ldq $6,s11l
	ldl $1,s11len
	ldl $3,s12len
	s4addq $1,$6,$5     # s11-end
	ldq $4,s11_begin
	ldl $3,s12len
	ldq $10,s12l_sort
	ldl $11,hash_shift

outerloop1:
	cmpult $6,$5,$1
	beq $1,hash1_end    # if s11>=s11-end break
	ldl $0,0($4)        # s11_begin[i]
	s4addq $0,$10,$2    # s12l_sort+j
	ldl $9,0($6)        # add=s11l[i]
	addq $6,4,$6        # s11++
	addq $4,4,$4        # s11_begin++

innerloop1:
	ldl $17,0($2)       # h=s12l_sort[j]
	addq $2,4,$2        # s12_sort++
	addq $17,$9,$17     # h+=add
	and $17,$16,$17
	and $17,$21,$1      # hp
	subq $1,$20,$1      # hp-=hsub
	bne $1,innerloop1end
	subq $17,$20,$19    # h-hsub
	sra $19,$11,$19     # ind
	addq $19,$7,$18     # sort+ind
	ldbu $0,0($18)      # sort[ind]
	cmplt $0,32,$1      # sort[ind]<32 ?
	beq $1,return1_1
	sll $0,12,$1        # sort[ind]*NHASH: 12=64-HASHSHIFT
	addq $0,1,$0
	addq $1,$19,$1
	stb $0,0($18)      # sort[ind]++
	s4addq $1,$8,$1
	stl $17,0($1)       # hash[sort[ind]*NHASH+ind]=h
	br innerloop1

innerloop1end:
	subq $2,$10,$2      # 4*(j+1)
	subq $2,4,$2        # 4*j
	sra $2,2,$2         # j
	mov 0,$0
	cmpult $2,$3,$1
	cmoveq $1,$3,$0
	subq $2,$0,$2       # if j>=s12len: j-=s12len
	stl $2,-4($4)       # s11_begin[i]=j
	br outerloop1

hash1_end: # return 0
	mov 0,$0
	ldq $9,8($30)
	ldq $10,16($30)
	ldq $11,24($30)
	lda $30,80($30)
	ret     $31,($26),1
return1_1: # return 1 on error
	mov 1,$0
	ldq $9,8($30)
	ldq $10,16($30)
	ldq $11,24($30)
	lda $30,80($30)
	ret     $31,($26),1
	.end    asm_hash1




 # asm_hash2(hb):
 #
 #  mask=(1<<hashpart_shift)-1; mask=~mask;
 #  hsub=hb<<hashpart_shift;
 #  sort=(uchar *)hashdata; hash=hashdata+(NHASH>>2);
 #  for (i=0; i<s21len; i++) {
 #    add=s21l[i];
 #    j=s21_begin[i];
 #    while (1) {
 #      h=add+s22l_sort[j];
 #      hp=h&mask;
 #      if (hp!=hsub) break;
 #      ind=(h-hsub)>>hash_shift;
 #      for (k=0; k<sort[ind]; k++) {
 #        if (hash[k*NHASH+ind]-h<raw_bound) {
 #          if (j>=s22len) *raw_cand_ptr++=(i+((j-s22len)<<6));
 #          else *raw_cand_ptr++=(i+(j<<6)); /* assumes npr_in_p<25 */
 #          break;
 #        }
 #      }
 #      j++;
 #    }
 #    if (j>=s22len) j-=s22len;
 #    s21_begin[i]=j;
 #  }


 # register usage
 # 0   sort[ind] et al
 # 1
 # 2   s22l_sort+j
 # 3   s22len
 # 4   s21_begin+i
 # 5   s21l+s11len
 # 6   s21l+i
 # 7   sort
 # 8   hash
 # 9   add
 # 10  s22l_sort
 # 11  hash_shift
 # 12  raw_bound
 # 16  32 bit mask
 # 17  h
 # 18
 # 19  ind
 # 20  hsub
 # 21  mask


        .set noreorder
        .align 3
        .globl asm_hash2
        .ent asm_hash2
asm_hash2:
        .frame $30,0,$26,0
        .prologue 0
	lda $30,-80($30)
	stq $9,8($30)
	stq $10,16($30)
	stq $11,24($30)
	stq $12,32($30)
	ldl $12,raw_bound
	ldq $7,hashdataptr
	ldiq $1,1
	sll $1,12,$1        # 12=64-HASHSHIFT
	addq $1,$7,$8       # hashdata

	ldl $1,hashpart_shift
	sll $16,$1,$20      # hsub=hb<<hashpart_shift
	ldiq $21,1
	sll $21,$1,$21      # 2^hashpart_shift
	negq $21,$21        # mask
	ldiq $16,1
	sll $16,32,$16
	subq $16,1,$16      # 0x00000000ffffffff in $16
	ldq $6,s21l
	ldl $1,s21len
	ldl $3,s22len
	s4addq $1,$6,$5     # s21-end
	ldq $4,s21_begin
	ldl $3,s22len
	ldq $10,s22l_sort
	ldl $11,hash_shift

outerloop2:
	cmpult $6,$5,$1
	beq $1,hash2_end    # if s21>=s21-end break
	ldl $0,0($4)        # s21_begin[i]
	s4addq $0,$10,$2    # s22l_sort+j
	ldl $9,0($6)        # add=s21l[i]
	addq $6,4,$6        # s21++
	addq $4,4,$4        # s21_begin++

innerloop2:
	ldl $17,0($2)       # h=s22l_sort[j]
	addq $2,4,$2        # s22_sort++
	addq $17,$9,$17     # h+=add
	and $17,$16,$17
	and $17,$21,$1      # hp
	subq $1,$20,$1      # hp-=hsub
	bne $1,innerloop2end
	subq $17,$20,$19    # h-hsub
	sra $19,$11,$19     # ind
	addq $19,$7,$18     # sort+ind
	ldbu $0,0($18)      # sort[ind]

	beq $0,innerloop2
	s4addq $19,$8,$19   # hashdata+ind
	sll $0,12,$0        # NHASH*sort[ind]: 12=64-HASHSHIFT
	s4addq $0,$19,$0    # end
	mov 1,$18
	sll $18,12,$18      # NHASH: 12=64-HASHSHIFT

testloop:
	ldl $1,0($19)      # hash[k*NHASH+ind]
	subq $1,$17,$1     # hash[k*NHASH+ind]-h
	and $1,$16,$1
	cmpult $1,$12,$1
	beq $1,nostore      # if hash[k*NHASH+ind]-h>=raw_bound goto nostore
 # store i+((j mod s22len)<<6), probably not time critical
 # $1, $18, $19, $0 free
	subq $2,$10,$0      # 4*(j+1)
	subq $0,4,$0        # 4*j
	sra $0,2,$0         # j
	mov 0,$18
	cmpult $0,$3,$1
	cmoveq $1,$3,$18
	subq $0,$18,$0      # if j>=s22len: j-=s22len
	sll $0,6,$0
	ldq $1,s21l
	subq $6,$1,$1       # 4*(i+1)
	subq $1,4,$1
	sra $1,2,$1         # i
	addq $0,$1,$0       # i+((j mod s22len)<<6)
	ldq $19,raw_cand_ptr
	stl $0,0($19)
	addq $19,4,$19
	stq $19,raw_cand_ptr
	br innerloop2

nostore:
	s4addq $18,$19,$19
	cmpult $19,$0,$1
	bne $1,testloop
	br innerloop2

innerloop2end:
	subq $2,$10,$2      # 4*(j+1)
	subq $2,4,$2        # 4*j
	sra $2,2,$2         # j
	mov 0,$0
	cmpult $2,$3,$1
	cmoveq $1,$3,$0
	subq $2,$0,$2       # if j>=s22len: j-=s22len
	stl $2,-4($4)       # s21_begin[i]=j
	br outerloop2

hash2_end:
	ldq $9,8($30)
	ldq $10,16($30)
	ldq $11,24($30)
	ldq $12,32($30)
	lda $30,80($30)
	ret     $31,($26),1
	.end    asm_hash2

