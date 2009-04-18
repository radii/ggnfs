# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

## Let R=2^{64}, N an odd number, and let M be the modular inverse of
## -N modulo R. Our first task is to calculate M from N.
N	.ASSIGNC "$4"
M	.ASSIGNC "$5"
nt	.ASSIGNC "$2"
aux1	.ASSIGNC "$12"
aux2	.ASSIGNC "$13"
aux3	.ASSIGNC "$14"

## We use the following fact: If b is inverse to a mod N, then
## 2*b-b*b*a is inverse to a mod N^2

	.MACRO hensel_step
	dmultu \&M,\&M
	dsll \&M,\&M,1
	mflo \&aux1
	dmultu \&N,\&aux1
	mflo \&aux1
	dsubu \&M,\&aux1
	.ENDM

## Do the intitialization as explained above.
## Note that every odd number is inverse to itself modulo 8.
## Five Hensel steps improve this to 2^{96}, which is more than
## needed.

	.MACRO mm64_init
	move \&M,\&N
.AREPEAT 5
	hensel_step
.AENDR
	dnegu \&M
	.ENDM


## The following macro squares its argument, which must not be a member
## of { \&aux1; \&aux2; \&N; \&M }, the first two of which are clobbered.

	.MACRO mm64_sq arg
	dmultu \arg,\arg
# Throughout this explanation, let Hi and  Lo denot the hi and Lo
# after this dmultu
	mfhi \&aux1
	mflo \arg
# Result will be the modular addition of Hi with the modular product p
# of Lo and the modular inverse of 2^64. This equals modsub(p,N-Hi),
# N being the modulus of the modsub.
	dsubu \&aux2,\&N,\&aux1
	dmultu \arg,\&M
# Actually, the following chain of instructions will calculate p-1 instead
# of p, unless Lo equals zero. Therefore, decrement \&aux2
# (which currently holds N-Hi) unless \arg is zero.
	sne \&aux1,\arg,0
	mflo \arg
	dmultu \arg,\&N
	dsubu \&aux2,\&aux1
	mfhi \arg
# Now \arg holds p-1 or p, and \&aux2 holds N-Hi-1 or N-Hi, where the first
# case occurs iff Lo is not zero.
# The following five instruction subtract \&aux2 from \arg modulo \&N
	modsub \arg \&aux2 \&aux1
	.ENDM

## We have used the following macro, which calculates
## \arg1 - \arg2 (modulo \&N), clobbering its third argument.

	.MACRO modsub arg1 arg2 clobber
	sleu \clobber,\arg2,\arg1
	dsubu \clobber,1
	dsubu \arg1,\arg2
	and \clobber,\clobber,\&N
	daddu \arg1,\clobber
	.ENDM

	.MACRO modsub1 arg1 arg2 clobber rop
	sleu \clobber,\arg2,\arg1
	dsubu \clobber,1
	dsubu \rop,\arg1,\arg2
	and \clobber,\clobber,\&N
	daddu \rop,\clobber
	.ENDM

# It will also be necessary to calculate the Montgomery representation
# of 1.

one	.ASSIGNC "$6"

	.MACRO mm64_one
	dnegu \&one,\&N
	ddivu \&one,\&N
	mfhi \&one
	.ENDM

##\( Let e be an exponent in \&ex, which is >= 2^k but < 2^{k+1}.\)
##\(Let us also assume that 2^k is in \&exm. The\)
##\(following macro calculates 2^e [mod N] in \&p2, destroying \&exm\)

ex	.ASSIGNC "$7"
exm	.ASSIGNC "$8"
p2	.ASSIGNC "$9"

	.MACRO power
	.globl pt64_power
pt64_power:
	dsubu \&aux1,\&N,\&one
	modsub1 \&one \&aux1 \&aux2 \&p2
	dsrl \&exm,\&exm,1
	beqz \&exm,pt64_have_power
	move \&aux3,\&ex # Branch delay slot
	.globl pt64_power_loop
pt64_power_loop:
	mm64_sq \&p2
	and \&aux3,\&exm
	dsrl \&exm,\&exm,1
	seq \&aux1,\&aux3,0
	dsubu \&aux1,1
	move \&aux2,\&N
	and \&aux1,\&p2
	dsubu \&aux2,\&aux1
	modsub \&p2 \&aux2 \&aux1
	bnez \&exm,pt64_power_loop
	move \&aux3,\&ex # Delay slot
	.globl pt64_have_power
pt64_have_power:
	.ENDM

	.MACRO get_msb rop op
	and \&aux1,\op,0xffffffff
	dsrl \&aux2,\op,32
	sne \rop,\&aux2,0
	movn \&aux1,\&aux2,\&aux2
	sll \rop,5
	dsrl \&aux2,\&aux1,16
	and \&aux1,\&aux1,0xffff
	sne \&aux3,\&aux2,0
	movn \&aux1,\&aux2,\&aux2
	sll \&aux3,4
	dsrl \&aux2,\&aux1,8
	and \&aux1,\&aux1,0xff
	addu \rop,\&aux3
	sne \&aux3,\&aux2,0
	movn \&aux1,\&aux2,\&aux2
	sll \&aux3,3
	dsrl \&aux2,\&aux1,4
	and \&aux1,\&aux1,0xf
	addu \rop,\&aux3
	sne \&aux3,\&aux2,0
	movn \&aux1,\&aux2,\&aux2
	sll \&aux3,2
	dsrl \&aux2,\&aux1,2
	and \&aux1,\&aux1,0x3
	addu \rop,\&aux3
	sne \&aux3,\&aux2,0
	movn \&aux1,\&aux2,\&aux2
	sll \&aux3,1
	dsrl \&aux2,\&aux1,1
	addu \rop,\&aux3
	sne \&aux3,\&aux2,0
	addu \rop,\&aux3
	.ENDM

	.MACRO get_lsb rop op
	dsrl \&aux2,\op,32
	and \&aux1,\op,0xffffffff
	seq \rop,\&aux1,0
	movz \&aux1,\&aux2,\&aux1
	sll \rop,5
	dsrl \&aux2,\&aux1,16
	and \&aux1,\&aux1,0xffff
	seq \&aux3,\&aux1,0
	movz \&aux1,\&aux2,\&aux1
	sll \&aux3,4
	dsrl \&aux2,\&aux1,8
	and \&aux1,\&aux1,0xff
	addu \rop,\&aux3
	seq \&aux3,\&aux1,0
	movz \&aux1,\&aux2,\&aux1
	sll \&aux3,3
	dsrl \&aux2,\&aux1,4
	and \&aux1,\&aux1,0xf
	addu \rop,\&aux3
	seq \&aux3,\&aux1,0
	movz \&aux1,\&aux2,\&aux1
	sll \&aux3,2
	dsrl \&aux2,\&aux1,2
	and \&aux1,\&aux1,0x3
	addu \rop,\&aux3
	seq \&aux3,\&aux1,0
	movz \&aux1,\&aux2,\&aux1
	sll \&aux3,1
	and \&aux1,\&aux1,0x1
	addu \rop,\&aux3
	seq \&aux3,\&aux1,0
	addu \rop,\&aux3
	.ENDM

	.file 1 "pt64.asm"
	.text
	.align 4
	.globl pt64
	.ent pt64

pt64:
	.set noreorder
	dsubu \&ex,\&N,1
##\(At the moment, we are free to use \&M as auxilliary  variables.\)
	get_lsb \&nt,\&ex
	get_msb \&M,\&ex
	.globl pt64_have_exponents
pt64_have_exponents:
	subu \&M,\&nt
	li \&exm,1
	beqz \&nt,no_prime # even input
	dsrl \&ex,\&nt
	dsll \&exm,\&M
	mm64_init
	mm64_one
	power
	.globl pt64_first_test
pt64_first_test:
	beq \&p2,\&one,prime
	move \&exm,\&N # Delay slot
	dsubu \&exm,\&one
	subu \&nt,1
### Now, \&exm holds the Montgomery representation of -1.
	.globl pt64_test_loop
pt64_test_loop:
	beq \&p2,\&exm,prime
	mm64_sq \&p2
	bnez \&nt,pt64_test_loop
	subu \&nt,1 # Delay slot
no_prime:
	.align 4
	j $31
	li $2,0
prime:
	.align 4
	j $31
	li $2,1
	.END
