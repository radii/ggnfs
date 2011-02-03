# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

log_ptr		.ASSIGNC	 "$4"
log_ptr_ub	.ASSIGNC	 "$5"
sched_ptr	.ASSIGNC	 "$6"
sieve_interval	.ASSIGNC	 "$7"
slog		.ASSIGNC	 "$8"
sched		.ASSIGNC	 "$9"
sched_ub	.ASSIGNC	"$10"
aux1		.ASSIGNC	"$11"
aux2		.ASSIGNC	"$12"
aux3		.ASSIGNC	"$13"

sched_ub_safety	.ASSIGNC	"8"
sched_ub_s1	.ASSIGNC	"12"

	.section .text
	.align 4
	.globl schedsieve
	.ent schedsieve
	.set noreorder
	.set nomacro
	.set noat

schedsieve:
### The caller gives us the number of schedule pieces as the
### second argument. Transform this into a bound on log_ptr
	addu	\&log_ptr_ub,	\&log_ptr_ub,	\&log_ptr
	sltu	\&aux1,		\&log_ptr,	\&log_ptr_ub
	lw	\&sched_ub,	(\&sched_ptr)
	addu	\&sched_ptr,	\&sched_ptr,	4
	beqz	\&aux1,		schedsieve_end
	move	\&sched,	\&sched_ub	# Delay slot
	.globl schedsieve_logloop
schedsieve_logloop:
	addu	\&sched,	\&sched,	2
	lhu	\&aux1,		(\&sched)
	addu	\&sched,	\&sched,	4
	lw	\&sched_ub,	(\&sched_ptr)
	addu	\&sched_ptr,	\&sched_ptr,	4
	subu	\&sched_ub,	\&sched_ub,	\&sched_ub_safety
	lb	\&slog,		(\&log_ptr)
	sltu	\&aux3,		\&sched,	\&sched_ub
	addu	\&log_ptr,	\&log_ptr,	1
	beqz	\&aux3,		schedsieve_sieveloop_end
	addu	\&aux1,		\&aux1,		\&sieve_interval
	.globl schedsieve_sieveloop
schedsieve_sieveloop:
	pref	4,		256(\&sched)
.AREPEAT 2
	lb	\&aux2,		(\&aux1)
	lhu	\&aux3,		(\&sched)
	addu	\&sched,	\&sched,	4
	addu	\&aux2,		\&aux2,		\&slog
	addu	\&aux3,		\&aux3,		\&sieve_interval
	sb	\&aux2,		(\&aux1)
	lb	\&aux2,		(\&aux3)
	lhu	\&aux1,		(\&sched)
	addu	\&sched,	\&sched,	4
	addu	\&aux2,		\&aux2,		\&slog
	addu	\&aux1,		\&aux1,		\&sieve_interval
	sb	\&aux2,		(\&aux3)
.AENDR
	sltu	\&aux3,		\&sched,	\&sched_ub
	bnez	\&aux3,		schedsieve_sieveloop
	nop
	.globl schedsieve_sieveloop_end
schedsieve_sieveloop_end:
	addu	\&sched_ub,	\&sched_ub,	\&sched_ub_s1
.AREPEAT 3
	sltu	\&aux3,		\&sched,	\&sched_ub
	beqz	\&aux3,		schedsieve_nextlog
	lb	\&aux2,		(\&aux1)
	addu	\&aux2,		\&aux2,		\&slog
	sb	\&aux2,		(\&aux1)
	lhu	\&aux1,		(\&sched)
	addu	\&sched,	\&sched,	4
	addu	\&aux1,		\&aux1,		\&sieve_interval
.AENDR
	.globl schedsieve_nextlog
schedsieve_nextlog:
	sltu	\&aux1,		\&log_ptr,	\&log_ptr_ub
	subu	\&sched_ub,	\&sched_ub,	4
	bnez	\&aux1,		schedsieve_logloop
	move	\&sched,	\&sched_ub	# Delay slot
	.globl schedsieve_end
schedsieve_end:
	j $31
	nop
	\(.end schedsieve)
	.END
