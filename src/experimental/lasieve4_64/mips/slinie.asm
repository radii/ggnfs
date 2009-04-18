# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

prime_ptr		.ASSIGNC	"$4"
prime_ptr_ub		.ASSIGNC	"$5"
sieve_interval		.ASSIGNC	"$6"
addr0			.ASSIGNC	"$7"
addr1			.ASSIGNC	"$8"
addr2			.ASSIGNC	"$9"
sv0			.ASSIGNC	"$10"
sv1			.ASSIGNC	"$11"
sv2			.ASSIGNC	"$12"
prime			.ASSIGNC	"$2"
r			.ASSIGNC	"$3"
slog			.ASSIGNC	"$14"
linectr			.ASSIGNC	"$15"
aux0			.ASSIGNC	"$16"
aux1			.ASSIGNC	"$17"
prime7			.ASSIGNC	"$18"
n_i_minus_prime7	.ASSIGNC	"$19"
line_ctr		.ASSIGNC	"$20"
addr_ub			.ASSIGNC	"$21"
pr			.ASSIGNC	"$1"

	.section .text
	.align 4
	.globl slinie
	.ent slinie
	.set noreorder
	.set nomacro
	.set noat
slinie:
	.frame	$sp, 16, $31
	subu	$sp,		$sp,		48
	sd	$16,		  ($sp)
	sd	$17,		 8($sp)
	sd	$18,		16($sp)
	sd	$19,		24($sp)
	sd	$20,		32($sp)
	sd	$21,		40($sp)
	.globl slinie_fbi_loop
slinie_fbi_loop:
	lhu	\&prime,		(\&prime_ptr)
	sltu	\&slog,		\&prime_ptr,	\&prime_ptr_ub
	addu	\&prime_ptr,	\&prime_ptr,	2
	lhu	\&pr,		(\&prime_ptr)
	addu	\&prime_ptr,	\&prime_ptr,	2
	addu	\&aux0,		\&prime,	\&prime
	beqz	\&slog,		slinie_end
### pr is the projective root. Negate it, since we will carry out modular
### subtractions.
	move	\&sv0,		\&pr
	lhu	\&slog,		(\&prime_ptr)
	addu	\&prime_ptr,	\&prime_ptr,	2
	addu	\&aux1,		\&aux0,		\&aux0
	addu	\&aux0,		\&aux0,		\&prime
	lhu	\&r,		(\&prime_ptr)
	move	\&addr_ub,	\&sieve_interval
	addu	\&prime7,	\&aux1,		\&aux0
	li	\&n_i_minus_prime7,		\&n_i
	li	\&line_ctr,	\&j_per_strip
	subu	\&pr,		\&prime,	\&pr
	subu	\&n_i_minus_prime7,	\&n_i_minus_prime7,	\&prime7
	addu	\&prime_ptr,	\&prime_ptr,	2
	movz	\&pr,		\&sv0,		\&sv0
	.globl slinie_strip_loop
slinie_strip_loop:
	addu	\&addr0,	\&addr_ub,	\&r
	addu	\&addr_ub,	\&addr_ub,	\&n_i_minus_prime7
	beqz	\&line_ctr,	slinie_fbi_loop
	subu	\&line_ctr,	\&line_ctr,		1	# Delay  slot
	.globl slinie_line_loop
slinie_line_loop:
	lb	\&sv0,		(\&addr0)
	addu	\&addr1,	\&addr0,	\&prime
	sltu	\&aux0,		\&addr0,	\&addr_ub
	lb	\&sv1,		(\&addr1)
	addu	\&addr2,	\&addr1,	\&prime
	lb	\&sv2,		(\&addr2)
	addu	\&sv0,		\&sv0,		\&slog
	sb	\&sv0,		(\&addr0)
	addu	\&sv1,		\&sv1,		\&slog
	addu	\&sv2,		\&sv2,		\&slog
	sb	\&sv1,		(\&addr1)
	addu	\&addr1,	\&addr2,	\&prime
	lb	\&sv1,		(\&addr1)
	sb	\&sv2,		(\&addr2)
	addu	\&addr0,	\&addr1,	\&prime
	addu	\&sv1,		\&sv1,		\&slog
	sb	\&sv1,		(\&addr1)
	bnez	\&aux0,		slinie_line_loop
	nop

	addu	\&addr_ub,	\&addr_ub,	\&prime7
### Note that we have negated pr, hence we need a modular subtraction
	sltu	\&aux1,		\&r,		\&pr
	subu	\&r,		\&r,		\&pr
	sltu	\&aux0,		\&addr0,	\&addr_ub
	movn	\&aux1,		\&prime,	\&aux1
	beqz	\&aux0,		slinie_strip_loop
#	b slinie_strip_loop
	addu	\&r,		\&r,		\&aux1
	.globl	slinie_nochmal
slinie_nochmal:
	lb	\&sv0,		(\&addr0)
	addu	\&addr1,	\&addr0,	\&prime
	addu	\&sv0,		\&sv0,		\&slog
	sltu	\&aux0,		\&addr1,	\&addr_ub
	sb	\&sv0,		(\&addr0)
	beqz	\&aux0,		slinie_strip_loop
	nop
	.globl slinie_noch_ein_2tes_mal
slinie_noch_ein_2tes_mal:
	lb	\&sv1,		(\&addr1)
	addu	\&addr0,	\&addr1,	\&prime
	addu	\&sv1,		\&sv1,		\&slog
	sltu	\&aux0,		\&addr0,	\&addr_ub
	sb	\&sv1,		(\&addr1)
	beqz	\&aux0,		slinie_strip_loop
	nop

	.globl	slinie_zum_dritten_und_letzten
slinie_zum_dritten_und_letzten:
	lb	\&sv0,		(\&addr0)
	addu	\&sv0,		\&sv0,		\&slog
	sb	\&sv0,		(\&addr0)
	b			slinie_strip_loop
	nop

	.globl slinie_end
slinie_end:
	ld	$16,		  ($sp)
	ld	$17,		 8($sp)
	ld	$18,		16($sp)
	ld	$19,		24($sp)
	ld	$20,		32($sp)
	ld	$21,		40($sp)
	jr	$31
	addu	$sp,		$sp,		48
	\(.end slinie)
	.END
