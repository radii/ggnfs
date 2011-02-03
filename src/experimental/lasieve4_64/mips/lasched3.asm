# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

ri_ptr		.ASSIGNC	"$4"
ij_ptr		.ASSIGNC	"$5"
ij_ptr_ub	.ASSIGNC	"$6"
ij_ub		.ASSIGNC	"$7"
sched_ptr	.ASSIGNC	"$8"
fbi_offs	.ASSIGNC	"$9"
lavec0		.ASSIGNC	"$2"
lavec1		.ASSIGNC	"$3"
labd0		.ASSIGNC	"$10"
labd1		.ASSIGNC	"$11"
aux1		.ASSIGNC	"$12"
aux2		.ASSIGNC	"$13"
aux3		.ASSIGNC	"$14"
aux4		.ASSIGNC	"$15"
aux5		.ASSIGNC	"$16"
ij		.ASSIGNC	"$17"
ij1		.ASSIGNC	"$1"

# Uses ula2.asm
	.MACRO lasched_head
.text
	.section .text
	.align 4
	.globl lasched
	.ent lasched
	.set noreorder
	.set nomacro
	.set noat
lasched:
	.frame	$sp, 16, $31
	sltu	\&aux1,		\&ij_ptr,	\&ij_ptr_ub
	sltu	\&aux2,		$10,		2
	subu	$sp, 16		
	sd	$16,		 0($sp)
	sd	$17,		 8($sp)
	beqz	\&aux1,		lasched_end0
## The loop over factor base elements will put the incrementation of
## \&ij_ptr into the delay slot of its branch. Therefore:
	subu	\&ij_ptr_ub,	4	# Delay slot
	beqz	\&aux2,		lasched23
	sltu	\&aux2,		$10,	3
	bnez	$10,		lasched1
	nop
	.globl lasched0
lasched0:
oddness_type	.ASSIGNA	0
	lasched_body
	.globl lasched1
lasched1:
oddness_type	.ASSIGNA	1
	lasched_body
	.globl lasched23
lasched23:
	beqz	\&aux2,		lasched3
	nop
	.globl lasched2
lasched2:
oddness_type	.ASSIGNA	2
	lasched_body
	.globl lasched3
lasched3:
oddness_type	.ASSIGNA	3
	lasched_body
	.ENDM

	.MACRO lasched_body
	sll	\&fbi_offs,	16

## The program gives us an upper bound for j, not the combination of i and
## j. Calculate what we need.
	sll	\&ij_ub,	\&i_bits

	.globl lasched_fbi_loop\&oddness_type
## Load them for first iteration of the outer loop.
	lwu	\&lavec0,	0(\&ri_ptr)
# Afterwards, this is done while updating \(\&ij_ptr)
	lwu	\&lavec1,	4(\&ri_ptr)
lasched_fbi_loop\&oddness_type:
.AIF	\&oddness_type	EQ		0
	lwu	\&ij,		0(\&ij_ptr)
	sltu	\&aux3,		\&ij,		\&ij_ub	#Bedingter Sprung unten
.AENDI
	pref	5,		128(\&ij_ptr)
.AIF	\&longer_ri		GT		0
	lhu	\&labd0,	6(\&ri_ptr)
	lhu	\&labd1,	8(\&ri_ptr)
ri_size		.ASSIGNC	"12"
.AELSE
	li	\&labd0,	\&n_i
	li	\&labd1,	\&n_i
	and	\&aux1,		\&lavec0,	\&n_i_mask
	and	\&aux2,		\&lavec1,	\&n_i_mask
	subu	\&labd0,	\&aux1
	subu	\&labd1,	\&aux2
.AIF	\&oddness_type	GT		0
	get_ij_from_ri lasched_ot2a	lasched_ot2b
.AENDI
ri_size		.ASSIGNC	"8"
.AENDI
	pref	4,		256(\&ri_ptr)
	beqz	\&aux3,		lasched_riloop_end\&oddness_type
	addu	\&ri_ptr,	\&ri_size
### From the second loop onward, the following instruction will
### run in the delay slot.
	srl	\&aux3,		\&ij,		\&l1_bits
	.globl lasched_ri_loop\&oddness_type
lasched_ri_loop\&oddness_type:
	ula_next \&ij1 \&ij \&aux1 \&aux2
	sll	\&aux3,		\&aux3,		2 # 4 Bytes/Sched entry
	and	\&aux4,		\&ij,		\&l1_mask
	ula_n2 \&ij1 \&ij \&aux1 \&aux2
	addu	\&aux3,		\&sched_ptr
	ula_n3 \&ij1 \&ij \&aux1 \&aux2
	lwu	\&aux5,		0(\&aux3)
	ula_n4 \&ij1 \&ij \&aux1 \&aux2
	ula_n5 \&ij1 \&ij \&aux1 \&aux2
	or	\&aux4,		\&fbi_offs
	ula_n6 \&ij1 \&ij \&aux1 \&aux2
	sw	\&aux4,		0(\&aux5)
	pref	5,		128(\&aux5)
	ula_n7 \&ij1 \&ij \&aux1 \&aux2
	addu	\&aux5,		4
	ula_n8 \&ij1 \&ij \&aux1 \&aux2
	sw	\&aux5,		0(\&aux3)
	sltu	\&aux3,		\&ij1,		\&ij_ub
	ula_next \&ij \&ij1 \&aux1 \&aux2
	beqz	\&aux3,		lasched_riloop_end\&oddness_type
	srl	\&aux3,		\&ij1,		\&l1_bits	# Delay slot
	sll	\&aux3,		\&aux3,		2 # 4 Bytes/Sched entry
	ula_n2 \&ij \&ij1 \&aux1 \&aux2
	and	\&aux4,		\&ij1,		\&l1_mask
	ula_n3 \&ij \&ij1 \&aux1 \&aux2
	addu	\&aux3,		\&sched_ptr
	ula_n4 \&ij \&ij1 \&aux1 \&aux2
	lwu	\&aux5,		0(\&aux3)
	ula_n5 \&ij \&ij1 \&aux1 \&aux2
	or	\&aux4,		\&fbi_offs
	ula_n6 \&ij \&ij1 \&aux1 \&aux2
	sw	\&aux4,		0(\&aux5)
	pref	5,		128(\&aux5)
	ula_n7 \&ij \&ij1 \&aux1 \&aux2
	addu	\&aux5,		4
	ula_n8 \&ij \&ij1 \&aux1 \&aux2
	sw	\&aux5,		0(\&aux3)
	sgtu	\&aux3,		\&ij_ub,	\&ij
	bnez	\&aux3,		lasched_ri_loop\&oddness_type
	srl	\&aux3,		\&ij,		\&l1_bits	# Delay slot
	.globl lasched_riloop_end\&oddness_type
lasched_riloop_end\&oddness_type:
### Increment what has to be incremented. For the pointer to the
### recurrence info, this has already been done in the dealay slot of
### the first test for the upper bound on \&ij.
### Also, load first recurrence info vector for next iteration of fbi loop
	lwu	\&lavec0,	0(\&ri_ptr)
.AIF	\&lasched_updates	GT	0
### If we reached this point by the non-taken branch of the last beqz, or 
### without ever passing through lasched_ri_loop, then everything is OK.
### But otherwise (we jumped to this point from the middle of lasched_ri_loop),
### it is necessary to replace \&ij by \&ij1.
	sltu	\&aux1,		\&ij,	\&ij_ub
	movn	\&ij,		\&ij1,	\&aux1
	li	\&aux2,		0x10000
### Now \&ij holds the correct value. Subtract \&ij_ub and store the
### difference while loading second rec info vector for next iteration.
	lwu	\&lavec1,	4(\&ri_ptr)
	sltu	\&aux3,		\&ij_ptr,	\&ij_ptr_ub
	subu	\&ij,		\&ij_ub
	addu	\&fbi_offs,	\&aux2
	sw	\&ij,		0(\&ij_ptr)
.AELSE
	li	\&aux1,		0x10000
	sltu	\&aux3,		\&ij_ptr,	\&ij_ptr_ub
### Second vector, next iteration of fbi loop
	lwu	\&lavec1,	4(\&ri_ptr)
	addu	\&fbi_offs,	\&aux1
.AENDI
	bnez	\&aux3		lasched_fbi_loop\&oddness_type
	addu	\&ij_ptr,	4 # Delay slot
	.globl lasched_end\&oddness_type
lasched_end\&oddness_type:
	ld	$16,		0($sp)
	ld	$17,		8($sp)
	addu	$sp,16
	j $31
	move	$2,		\&ri_ptr
	.ENDM

	.MACRO get_ij_from_ri label1 label2
ot_mask		.ASSIGNA	\&n_i + 1
ot_testa	.ASSIGNA	(\&oddness_type/2)*\&n_i + (\&oddness_type & 1)
ot_testb	.ASSIGNA	(1-\&oddness_type/2)*\&n_i+(\&oddness_type & 1)
.AIF	\&oddness_type	EQ		2
#the following macro jumps to label1 if this is a generic case.
#otherwise, it treats the degenerate case and jumps to label2.
	get_ij_from_ri_degenerate \label1 \label2
	.globl	\label1
\label1:
.AENDI
	get_ij_from_ri_generic
.AIF	\&oddness_type	EQ		2
	.globl	\label2
\label2:
.AENDI
.AIF	\&oddness_type	EQ		1
	addu	\&ij,	\&ij,	\&n_i
.AENDI
	srl	\&ij,	\&ij,	1
	sltu	\&aux3,		\&ij,		\&ij_ub	#Bedingter Sprung unten
	.ENDM

	
	.MACRO get_ij_from_ri_degenerate label1 label2
	sgtu	\&ij1,		\&aux1,		\&aux2
	sgtu	\&aux3,		\&lavec0,	\&lavec1
	bgtzl	\&ij1,		\label1
	nop
	sltu	\&ij1,		\&aux1,		\&aux2
	bgtz	\&aux3,		\label1
	subu	\&aux3,		\&lavec1,	\&lavec0
	li	\&aux4,		\&n_i
	movz	\&ij,		\&aux3,		\&ij1
	b	\label2
	movn	\&ij,		\&aux4,		\&ij1
	.ENDM

	.MACRO	get_ij_from_ri_generic
	and	\&aux1,		\&lavec0,	\&ot_mask
	and	\&aux2,		\&lavec1,	\&ot_mask
	move	\&ij,		$0
	subu	\&aux1,		\&aux1,		\&ot_testa
	subu	\&aux2,		\&aux2,		\&ot_testb
	movz	\&ij,		\&lavec0,	\&aux1
	movz	\&ij,		\&lavec1,	\&aux2
	movz	\&aux2,		\&lavec1,	\&ij
	movn	\&aux2,		$0,		\&ij
	movz	\&ij,		\&lavec0,	\&ij
	addu	\&ij,		\&ij,		\&aux2
	.ENDM

### Now, call the macro which generates the lasched function.
	lasched_head
	\(.end lasched)
