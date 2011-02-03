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
A		.ASSIGNC	"$6"
b		.ASSIGNC	"$7"
s		.ASSIGNC	"$8"
c		.ASSIGNC	"$9"
t		.ASSIGNC	"$10"
xc		.ASSIGNC	"$11"
yc		.ASSIGNC	"$12"
aux		.ASSIGNC	"$13"
aux2		.ASSIGNC	"$14"
nts		.ASSIGNA	8

	.MACRO euclid_step k kc l lc
	divu $0,\l,\k
	mflo \&aux
	mfhi \l
	multu \&aux,\kc
	mflo \&aux
	addu \lc,\&aux
	.ENDM

	.MACRO subtra k kc l lc subtlabel
	subu \l,\k
	bltu \l,\k,\subtlabel
	addu \lc,\kc	
	.ENDM

# %eax and %edx are used in the division.

.section	.rodata
.error_string:
	.SDATA "Bad args to asm_ri_aux\n"
.text
	.align 4
	.globl asm_getbc
	.ent asm_getbc

asm_getbc:
	.set noreorder
#	.set	nomacro
#	.set	noat
	li \&xc,1
	bltu \&x,\&A,have_bs
	li \&yc,0
divide:
.AREPEAT	\&nts
	subtra \&x \&xc \&y \&yc test_y
.AENDR
	euclid_step \&x \&xc \&y \&yc
test_y:
	bltu \&y,\&A,have_ct
	nop
.AREPEAT	\&nts
	subtra \&y \&yc \&x \&xc test_x
.AENDR
	euclid_step \&y \&yc \&x \&xc
test_x:	
	bgeu \&x,\&A,divide
	nop
have_bs:
.AREPEAT \&nts
	subu \&y,\&x
	bltu \&y,\&A,have_bsct
	addu \&yc,\&xc
.AENDR
	subu \&aux,\&y,\&A
	divu $0,\&aux,\&x
	mflo \&aux
	addu \&aux,1
	multu \&aux,\&x
	mflo \&aux2
	multu \&aux,\&xc
	subu \&y,\&y,\&aux2
	mflo \&aux
	addu \&yc,\&yc,\&aux
	j have_bsct
	nop	# Delay slot
	.align 4
have_ct:
.AREPEAT \&nts
	subu \&x,\&y
	bltu \&x,\&A,have_bsct
	addu \&xc,\&yc
.AENDR
	subu \&aux,\&x,\&A
	divu $0,\&aux,\&y
	mflo \&aux
	addu \&aux,1
	multu \&aux,\&y
	mflo \&aux2
	multu \&aux,\&yc
	subu \&x,\&x,\&aux2
	mflo \&aux
	addu \&xc,\&xc,\&aux
have_bsct:
	sw \&x,(\&b)
	sw \&y,(\&c)
	sw \&xc,(\&s)
	sw \&yc,(\&t)
	j $31
	nop
	\(.end asm_getbc)
	.END

