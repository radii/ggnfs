# Copyright (C) 2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the 
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# The function we want to write.
# ulong asm_mi32(x,y)
# Modular inverse of x modulo y, where y is odd and the gcd is 1
# and y satisfies 0<y<x.

x		.ASSIGNC	"$4"
y		.ASSIGNC	"$5"
xc		.ASSIGNC	"$6"
yc		.ASSIGNC	"$7"
aux		.ASSIGNC	"$2"
nts		.ASSIGNA	16

	.MACRO euclid_step k kc l lc
	divu $0,\l,\k
	mflo \&aux
	mfhi \l
	multu \&aux,\kc
	mflo \&aux
	addu \lc,\&aux
	.ENDM

	.MACRO subtra k kc l lc subtlabel need_first_diff=1
.AIF	\need_first_diff GT 0
	subu \l,\k
.AENDI
	bltu \l,\k,\subtlabel
	addu \lc,\kc	
	.ENDM

# %eax and %edx are used in the division.

.section	.rodata
.error_string:
#	.SDATA "Bad args to asm_modinv32\n"
		.byte	66,97,100,32,97,114,103,115,32,116,111,32,97,115,109,95,109,111,100,105,110,118,51,50,10
	.align 4
.errmsg:
	.word .error_string
.text
	.align 4
	.globl asm_modinv32
	.ent asm_modinv32

asm_modinv32:
	.set noreorder
	.frame $sp, 0, $31
	lw \&y,modulo32
	bleu \&x,1,badargs
	li \&xc,1
	bgeu \&x,\&y,badargs
	li \&yc,0
	subu \&y,\&x
divide:
ntsc	.ASSIGNA 0
.AWHILE	\&ntsc LT \&nts
	subtra \&x \&xc \&y \&yc test_y \&ntsc
ntsc	.ASSIGNA \&ntsc + 1
.AENDW
	euclid_step \&x \&xc \&y \&yc
test_y:
	bleu \&y,1,have_inverse1
ntsc	.ASSIGNA 0
.AWHILE	\&ntsc LT \&nts
	subtra \&y \&yc \&x \&xc test_x
ntsc	.ASSIGNA \&ntsc + 1
.AENDW
	euclid_step \&y \&yc \&x \&xc
test_x:	
	bgtu \&x,1,divide
	subu \&y,\&x
have_inverse2:
	bne \&x,1,badargs
	move $2,\&xc
	.set	nomacro
	j $31
	.align 4
	.set	macro
have_inverse1:
	lw $2,modulo32
	bne \&y,1,badargs
	subu $2,\&yc
ret_inv1:
	.set	nomacro
	j $31
	nop
badargs:
	beq \&x,1,ret_inv1
	li	$2,1 
	lw $4,.errmsg
	la	$25,Schlendrian
	jal	$31,$25
	.END
