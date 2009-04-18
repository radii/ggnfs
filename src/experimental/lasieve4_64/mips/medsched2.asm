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
sched_ptr	.ASSIGNC	"$7"
fbi_offs	.ASSIGNC	"$8"
lavec0		.ASSIGNC	"$2"
lavec1		.ASSIGNC	"$3"
labd0		.ASSIGNC	"$9"
labd1		.ASSIGNC	"$10"
aux5		.ASSIGNC	"$16"
ij		.ASSIGNC	"$11"
sched		.ASSIGNC	"$12"
aux1		.ASSIGNC	"$13"
aux2		.ASSIGNC	"$14"
aux3		.ASSIGNC	"$15"
ij1		.ASSIGNC	"$1"
ij_ub		.ASSIGNC	"$16"

.text

	.MACRO medsched_head
.text
	.section .text
	.align 4
	.globl medsched
	.ent medsched
	.set noreorder
	.set nomacro
	.set noat
medsched:
	.frame	$sp, 16, $31
	subu	$sp,		$sp,		16	
	sltu	\&aux2,		$9,		2
	sd	$16,		0($sp)
	sd	$3,		8($sp)
	sltu	\&aux1,		\&ij_ptr,	\&ij_ptr_ub
	lwu	\&sched,	0(\&sched_ptr)
	beqz	\&aux1,		medsched_end0
	sll	\&fbi_offs,	\&fbi_offs,	16	# Delay slot
## The loop over factor base elements will put the incrementation of
## \&ij_ptr into the delay slot of its branch. Therefore:
	subu	\&ij_ptr_ub,	4
	beqz	\&aux2,		medsched23
	sltu	\&aux2,		$9,		3
	bnez	$9,		medsched1
	nop
	.globl medsched0
medsched0:
oddness_type	.ASSIGNA	0
	medsched_body
	.globl medsched1
medsched1:
oddness_type	.ASSIGNA	1
	medsched_body
	.globl medsched23
medsched23:
	beqz	\&aux2,		medsched3
	nop
	.globl medsched2
medsched2:
oddness_type	.ASSIGNA	2
	medsched_body
	.globl medsched3
medsched3:
oddness_type	.ASSIGNA	3
	medsched_body
	.ENDM
	
	.MACRO medsched_body
	li	\&ij_ub,	\&l1_size
	.globl medsched_fbi_loop\&oddness_type
medsched_fbi_loop\&oddness_type:
.AIF	\&oddness_type	EQ		0
	lwu	\&ij,		0(\&ij_ptr)
.AENDI
	lwu	\&lavec0,	0(\&ri_ptr)
	pref	5,		128(\&ij_ptr)
	lwu	\&lavec1,	4(\&ri_ptr)
.AIF	\&longer_ri	GT	 0
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
	get_ij_from_ri medsched_ot2a	medsched_ot2b
.AENDI
ri_size		.ASSIGNC	"8"
.AENDI
	sltu	\&aux3,		\&ij,		\&ij_ub
	pref	4,		128(\&ri_ptr)
	beqz	\&aux3,		medsched_riloop_end\&oddness_type
	addu	\&ri_ptr,	\&ri_size
	ula_next \&ij1 \&ij \&aux1 \&aux2
### From the second loop onward, the following instruction will
### run in the delay slot.
	or	\&aux3,		\&ij,		\&fbi_offs
	.globl medsched_ri_loop\&oddness_type
medsched_ri_loop\&oddness_type:
	ula_n2 \&ij1 \&ij \&aux1 \&aux2	
	ula_n3 \&ij1 \&ij \&aux1 \&aux2
	ula_n5 \&ij1 \&ij \&aux1 \&aux2
	sw	\&aux3,		0(\&sched)
	ula_n4 \&ij1 \&ij \&aux1 \&aux2
	ula_n6 \&ij1 \&ij \&aux1 \&aux2
	ula_n7 \&ij1 \&ij \&aux1 \&aux2
	addu   \&sched,	4
	ula_n8 \&ij1 \&ij \&aux1 \&aux2
	pref	5,	256(\&sched)
	sltu	\&aux3,		\&ij1,		\&ij_ub
	ula_next \&ij \&ij1 \&aux1 \&aux2
	beqz	\&aux3,		medsched_riloop_end\&oddness_type
	or	\&aux3,		\&ij1,		\&fbi_offs	# Delay slot
	ula_n2 \&ij \&ij1 \&aux1 \&aux2
	ula_n3 \&ij \&ij1 \&aux1 \&aux2
	ula_n5 \&ij \&ij1 \&aux1 \&aux2
	sw	\&aux3,		0(\&sched)
	ula_n4 \&ij \&ij1 \&aux1 \&aux2
	ula_n6 \&ij \&ij1 \&aux1 \&aux2
	ula_n7 \&ij \&ij1 \&aux1 \&aux2
	addu	\&sched,	4
	ula_n8 \&ij \&ij1 \&aux1 \&aux2
	sgtu	\&aux3,		\&ij_ub,	\&ij
	ula_next \&ij1 \&ij \&aux1 \&aux2
	bnez	\&aux3,		medsched_ri_loop\&oddness_type
	or	\&aux3,		\&ij,		\&fbi_offs	# Delay slot
	.globl medsched_riloop_end\&oddness_type
medsched_riloop_end\&oddness_type:
### Increment what has to be incremented. For the pointer to the
### recurrence info, this has already been done in the dealay slot of
### the first test for the upper bound on \&ij.

### If we reached this point by the non-taken branch of the last beqz, or 
### without ever passing through medsched_ri_loop, then everything is OK.
### But otherwise (we reached this point from the middle of medsched_ri_loop),
### it is necessary to replace \&ij by \&ij1.
	sltu	\&aux1,		\&ij,	\&ij_ub
	movn	\&ij,		\&ij1,	\&aux1
	li	\&aux2,		0x10000
### Now \&ij holds the correct value. Subtract \&ij_ub and store the
### difference.
	sltu	\&aux3,		\&ij_ptr,	\&ij_ptr_ub
	subu	\&ij,		\&ij_ub
	addu	\&fbi_offs,	\&aux2
	sw	\&ij,		0(\&ij_ptr)
	bnez	\&aux3		medsched_fbi_loop\&oddness_type
	addu	\&ij_ptr,	4 # Delay slot
	.globl medsched_end\&oddness_type
medsched_end\&oddness_type:
	sw	\&sched,	(\&sched_ptr)
	ld	$16,		0($sp)
	ld	$3,		8($sp)
	addu	$sp,		16
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
	.ENDM

	
	.MACRO get_ij_from_ri_degenerate label1 label2
	sgtu	\&ij1,		\&aux1,		\&aux2
	sgtu	\&aux3,		\&lavec0,	\&lavec1
	bgtzl	\&ij1,		\label1
	nop
	sltu	\&ij1,		\&aux1,		\&aux2
	bgtz	\&aux3,		\label1
	subu	\&aux3,		\&lavec1,	\&lavec0
	li	\&aux2,		\&n_i
	movz	\&ij,		\&aux3,		\&ij1
	b	\label2
	movn	\&ij,		\&aux2,		\&ij1
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

	medsched_head
	\(.end medsched)
	.END
