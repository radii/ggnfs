# Copyright (C) 2002 Jens Franke, T.Kleinjung
# This file is part of gnfs4linux, distributed under the terms of the
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

FBptr		.ASSIGNC	"$16"
prptr		.ASSIGNC	"$17"
FBptr_ub	.ASSIGNC	"$18"
a0		.ASSIGNC	"$7"
a1		.ASSIGNC	"$8"
b0		.ASSIGNC	"$9"
b1		.ASSIGNC	"$10"
ri_ptr		.ASSIGNC	"$19"
prime		.ASSIGNC	"$20"
aux0		.ASSIGNC	 "$1"
aux1		.ASSIGNC	 "$2"
prime_f		.ASSIGNC	 "$f0"
pinv_f		.ASSIGNC	 "$f1"
x_f		.ASSIGNC	 "$f2"
y_f		.ASSIGNC	 "$f3"
aux0_f		.ASSIGNC	 "$f4"
aux1_f		.ASSIGNC	 "$f5"
a0f		.ASSIGNC	 "$f6"
a1f		.ASSIGNC	 "$f7"
b0f		.ASSIGNC	 "$f8"
b1f		.ASSIGNC	 "$f9"
proot_f		.ASSIGNC	"$f10"

.text
	.globl asm_lasieve_setup
	.ent asm_lasieve_setup
	.set noat
	.set noreorder
	.set nomacro
asm_lasieve_setup:
	.frame	$sp, 64, $31
	subu		$sp,		$sp,		64
	sd		$31,		  ($sp)
	sd		$16,		 8($sp)
	sd		$17,		16($sp)
	sd		$18,		24($sp)
	sd		$19,		32($sp)
	sd		$20,		40($sp)
	sd		$gp,		48($sp)
	lui	\&aux0,%hi(%neg(%gp_rel(asm_lasieve_setup)))
	addiu	\&aux0,\&aux0,%lo(%neg(%gp_rel(asm_lasieve_setup)))
	daddu	$gp,\&aux0,$25
	move		\&FBptr,	 $4
	move		\&prptr,	 $5
	move		\&FBptr_ub,	 $6
	move		\&ri_ptr	$11
	sll		\&FBptr_ub,	\&FBptr_ub,	2
	mtc1		\&a0,		\&a0f
	mtc1		\&a1,		\&a1f
	mtc1		\&b0,		\&b0f
	mtc1		\&b1,		\&b1f
	cvt.d.w		\&a0f,		\&a0f
	cvt.d.w		\&a1f,		\&a1f
	cvt.d.w		\&b0f,		\&b0f
	cvt.d.w		\&b1f,		\&b1f
	addu		\&FBptr_ub,	\&FBptr_ub,	\&FBptr
	sltu		\&aux0,		\&FBptr,	\&FBptr_ub
	beqz		\&aux0,		asm_lasieve_setup_end
	nop
	.globl asm_lasieve_setup_loop
asm_lasieve_setup_loop:
	lw		\&prime,	(\&FBptr)
	mtc1		\&prime,	\&prime_f
	pref		4,		256(\&FBptr)
	lw		\&aux0,		(\&prptr)
	mtc1		\&aux0,		\&proot_f
	sltu		\&aux1,		\&aux0,		\&prime
	lw		\&aux0,		%got_disp(modulo32)($gp)
	cvt.d.l		\&prime_f,	\&prime_f
	recip.d		\&pinv_f,	\&prime_f
	pref		4,		256(\&prptr)
	bnezl		\&aux1,		asm_lasieve_setup_affine_root
	cvt.d.l		\&proot_f,	\&proot_f	# Delay slot
	.globl asm_lasieve_setup_ab_infty
asm_lasieve_setup_ab_infty:
	mov.d		\&x_f,		\&b0f
	sw		\&prime,	(\&aux0)	# prime->modulo32
	neg.d		\&y_f,		\&b1f
	lw		$25,		%call16(asm_modinv32)($gp)
	mtc1		$0,		\&aux0_f
	b		asm_lasieve_setup_divide
	c.eq.d		\&aux0_f,	\&x_f		# Delay slot
asm_lasieve_setup_affine_root:
	mul.d		\&x_f,		\&b0f,		\&proot_f
	pref		5,		128(\&ri_ptr)
	mul.d		\&y_f,		\&b1f,		\&proot_f
	sub.d		\&x_f,		\&x_f,		\&a0f
	sub.d		\&y_f,		\&a1f,		\&y_f
	mul.d		\&aux0_f,	\&x_f,		\&pinv_f
	sw		\&prime,	(\&aux0)	# prime->modulo32
	mul.d		\&aux1_f,	\&y_f,		\&pinv_f
	floor.w.d	\&aux0_f,	\&aux0_f
	floor.w.d	\&aux1_f,	\&aux1_f
	cvt.d.w		\&aux0_f,	\&aux0_f
	cvt.d.w		\&aux1_f,	\&aux1_f
	lw		$25,		%call16(asm_modinv32)($gp)
	mul.d		\&aux0_f,	\&aux0_f,	\&prime_f
	mul.d		\&aux1_f,	\&aux1_f,	\&prime_f
	c.eq.d		$fcc0,		\&x_f,		\&aux0_f
	sub.d		\&x_f,		\&x_f,		\&aux0_f
	sub.d		\&y_f,		\&y_f,		\&aux1_f
asm_lasieve_setup_divide:
	cvt.w.d		\&x_f,		\&x_f
	bc1fl		$fcc0,		asm_lasieve_setup_call_modinv
	mfc1		$4,		\&x_f		# Delay slot
asm_lasieve_setup_infinity:
	lw		$25,		%call16(get_recurrence_info)($gp)
	move		$4,		\&ri_ptr
	move		$5,		\&prime
	move		$6,		\&prime
	b		asm_lasieve_setup_call_get_ri
	nop
	.globl asm_lasieve_setup_call_modinv
asm_lasieve_setup_call_modinv:
	jalr		$25
	nop
	mtc1		$2,		\&x_f
	cvt.d.w		\&x_f,		\&x_f
	mul.d		\&x_f,		\&x_f,		\&y_f
	move		$5,		\&prime
	mul.d		\&aux0_f,	\&x_f,		\&pinv_f
	floor.w.d	\&aux0_f,	\&aux0_f
	move		$4,		\&ri_ptr
	cvt.d.w		\&aux0_f,	\&aux0_f
	mul.d		\&aux0_f,	\&aux0_f,	\&prime_f
	lw		$25,		%call16(get_recurrence_info)($gp)
	sub.d		\&x_f,		\&x_f,		\&aux0_f
	cvt.w.d		\&x_f,		\&x_f
	mfc1		$6,		\&x_f
	.globl asm_lasieve_setup_call_get_ri
asm_lasieve_setup_call_get_ri:
	jalr		$25
	addu		\&FBptr,	\&FBptr,	4	#Delay slot
	addu		\&prptr,	\&prptr,	4
	sltu		\&aux0,		\&FBptr,	\&FBptr_ub
	bnez		\&aux0,		asm_lasieve_setup_loop
	addu		\&ri_ptr,	\&ri_ptr,	8
	.globl	asm_lasieve_setup_end
asm_lasieve_setup_end:
	ld		$31,		  ($sp)
	ld		$16,		 8($sp)
	ld		$17,		16($sp)
	ld		$18,		24($sp)
	ld		$19		32($sp)
	ld		$20,		40($sp)
	ld		$gp,		48($sp)
	jr		$31
	addu		$sp,		$sp,		64
	\(.end asm_lasieve_setup)
	.END
