# Copyright (C) 2002 Jens Franke
# This file is part of gnfs4linux, distributed under the terms of the 
# GNU General Public Licence and WITHOUT ANY WARRANTY.

# You should have received a copy of the GNU General Public License along
# with this program; see the file COPYING.  If not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.

# Universal lattice sieving macros
# Switch to the next sieving event. All  four arguments must be
# different.

	.MACRO ula_next rop op aux1 aux2
	and \aux1,\op,\&n_i_mask
	.ENDM

	.MACRO ula_n2 rop op aux1 aux2
	sltu \rop,\aux1,\&labd1
	.ENDM

	.MACRO ula_n3 rop op aux1 aux2
	sltu \aux2,\aux1,\&labd0
	.ENDM

	.MACRO ula_n4 rop op aux1 aux2
	subu \aux2,\aux2,1
	.ENDM

	.MACRO ula_n5 rop op aux1 aux2
	movn \rop,\&lavec0,\rop
	.ENDM

	.MACRO ula_n6 rop op aux1 aux2
	and \aux2,\aux2,\&lavec1
	.ENDM

	.MACRO ula_n7 rop op aux1 aux2
	addu \rop,\rop,\op
	.ENDM

	.MACRO ula_n8 rop op aux1 aux2
	addu \rop,\rop,\aux2
	.ENDM

