# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

fbi_offs_ptr	.ASSIGNC	 "$4"
fbi_offs_ptr_ub	.ASSIGNC	 "$5"
sched_ptr	.ASSIGNC	 "$6"
sieve_interval	.ASSIGNC	 "$7"
tds_fbi_curpos	.ASSIGNC	 "$8"
fbi_offset	.ASSIGNC	 "$9"
sched		.ASSIGNC	"$10"
sched_ub	.ASSIGNC	"$11"
aux1		.ASSIGNC	"$12"
aux2		.ASSIGNC	"$13"
aux3		.ASSIGNC	"$14"
aux4		.ASSIGNC	"$15"
aux5		.ASSIGNC	"$16"
aux6		.ASSIGNC	"$17"
aux7		.ASSIGNC	"$18"
aux8		.ASSIGNC	"$19"


	.section .text
	.align 4
	.globl schedtdsieve
	.ent schedtdsieve
	.set noreorder
	.set nomacro
	.set noat

schedtdsieve:
### The caller gives us the number of schedule pieces as the
### second argument. Transform this into a bound on fbi_offs_ptr
	sll	\&fbi_offs_ptr_ub,	\&fbi_offs_ptr_ub,		2
	subu	$sp,		$sp,	48
	addu	\&fbi_offs_ptr_ub,	\&fbi_offs_ptr_ub,	\&fbi_offs_ptr
	sd	$31,		($sp)
	sd	$16,		8($sp)
	sd	$17,		16($sp)
	sd	$18,		24($sp)
	sd	$19,		32($sp)
	sltu	\&aux1,		\&fbi_offs_ptr,	\&fbi_offs_ptr_ub
	lw	\&sched_ub,	(\&sched_ptr)
	addu	\&sched_ptr,	\&sched_ptr,	4
	beqz	\&aux1,		schedtdsieve_end
	move	\&sched,	\&sched_ub			# Delay slot
	.globl schedtdsieve_fbi_offsloop
schedtdsieve_fbi_offs_loop:
	lw	\&sched_ub,	(\&sched_ptr)
	addu	\&sched_ptr,	\&sched_ptr,	4
	lw	\&fbi_offset,	(\&fbi_offs_ptr)
	addu	\&sched,	\&sched,	2
	addu	\&sched_ub,	\&sched_ub,	4
	lhu	\&aux1,		(\&sched)
	addu	\&sched,	\&sched,	4
	addu	\&fbi_offs_ptr,	\&fbi_offs_ptr,	4
	lhu	\&aux4,		(\&sched)
	addu	\&sched,	\&sched,	4
	addu	\&aux1,		\&aux1,		\&sieve_interval
	lbu	\&aux6,		(\&aux1)
	lhu	\&aux7,		(\&sched)
	addu	\&sched,	\&sched,	4
	lhu	\&aux8,		(\&sched)
	addu	\&sched,	\&sched,	4
	sltu	\&aux1,		\&sched,	\&sched_ub
	move	\&aux3,		$0
	beqz	\&aux1,		schedtdsieve_sieveloop_end
	nop							# Delay slot
	.globl schedtdsieve_sieveloop
schedtdsieve_sieveloop:
	pref	4,		256(\&sched)
	lhu	\&aux1,		(\&sched)
	addu	\&aux4,		\&aux4,		\&sieve_interval
	lbu	\&aux5,		(\&aux4)
	addu	\&sched,	\&sched,	4
	or	\&aux3,		\&aux3,		\&aux6
	lhu	\&aux4,		(\&sched)
	addu	\&aux7,		\&aux7,		\&sieve_interval
	lbu	\&aux6,		(\&aux7)
	addu	\&sched,	\&sched,	4
	or	\&aux3,		\&aux3,		\&aux5
	lhu	\&aux7,		(\&sched)
	addu	\&aux8,		\&aux8,		\&sieve_interval
	lbu	\&aux5,		(\&aux8)
	addu	\&sched,	\&sched,	4
	or	\&aux3,		\&aux3,		\&aux6
	lhu	\&aux8,		(\&sched)
	addu	\&aux1,		\&aux1,		\&sieve_interval
	lbu	\&aux6,		(\&aux1)
	addu	\&sched,	\&sched,	4
	or	\&aux3,		\&aux3,		\&aux5
	bnez	\&aux3,		schedtdsieve_test
	sltu	\&aux5,		\&sched,	\&sched_ub	# Delay slot
	bnez	\&aux5,		schedtdsieve_sieveloop
	nop							# Delay slot
	.globl schedtdsieve_sieveloop_end
schedtdsieve_sieveloop_end:
	subu	\&sched_ub,	\&sched_ub,	4
	subu	\&sched,	\&sched,	16
.AREPEAT 3
	sltu	\&aux1,		\&sched,	\&sched_ub
	beqz	\&aux1,		schedtdsieve_next_fbi_offs
	lhu	\&aux1,		(\&sched)
	addu	\&aux1,		\&aux1,		\&sieve_interval
	lbu	\&aux2,		(\&aux1)
	subu	\&aux2,		\&aux2,		1
	bgezal	\&aux2,		schedtdsieve_store
	addu	\&sched,	\&sched,	4		# Delay slot
.AENDR
	.globl schedtdsieve_next_fbi_offs
schedtdsieve_next_fbi_offs:
	sltu	\&aux1,		\&fbi_offs_ptr,	\&fbi_offs_ptr_ub
	bnez	\&aux1,		schedtdsieve_fbi_offs_loop
	move	\&sched,	\&sched_ub			# Delay slot
	.globl schedtdsieve_end
schedtdsieve_end:
	ld	$31,	($sp)
	ld	$16,		8($sp)
	ld	$17,		16($sp)
	ld	$18,		24($sp)
	ld	$19,		32($sp)
	addu $sp,	$sp,	48
	j $31
	nop

	.globl schedtdsieve_test
schedtdsieve_test:
	subu	\&sched,	\&sched,	32
.AREPEAT 4
	lhu	\&aux1,		(\&sched)
	addu	\&aux1,		\&aux1,		\&sieve_interval
	lbu	\&aux2,		(\&aux1)
	subu	\&aux2,		\&aux2,		1
	bgezal	\&aux2,		schedtdsieve_store
	addu	\&sched,	\&sched,	4		# Delay slot
.AENDR
	addu	\&sched,	\&sched,	 16
	sltu	\&aux1,		\&sched,	\&sched_ub
	bnez	\&aux1,		schedtdsieve_sieveloop
	move	\&aux3,		$0				# Delay slot
	b	schedtdsieve_sieveloop_end
	nop							# Delay slot
	.globl schedtdsieve_store
schedtdsieve_store:
	subu	\&sched,	\&sched,	6
	lhu	\&aux3,		(\&sched)
	sll	\&aux2,		\&aux2,		2
	addu	\&aux2,		\&aux2,		\&tds_fbi_curpos
	lw	\&aux1,		(\&aux2)
	addu	\&aux3,		\&aux3,		\&fbi_offset
	sw	\&aux3,		(\&aux1)
	addu	\&aux1,		\&aux1,		4
	sw	\&aux1,		(\&aux2)
	jr $31
	addu	\&sched,	\&sched,	6		# Delay slot
	\(.end schedtdsieve)
	.END
