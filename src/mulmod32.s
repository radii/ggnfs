/**************************************************************/
/* mulmod32.s                                                 */
/* Copyright 2004, Chris Monico.                              */
/**************************************************************/
/*  This file is part of GGNFS.
*
*   GGNFS is free software; you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation; either version 2 of the License, or
*   (at your option) any later version.
*
*   GGNFS is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with GGNFS; if not, write to the Free Software
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

.text
	.align 4
.globl mulmod32
.globl _mulmod32
mulmod32:
_mulmod32:
	
	movl 4(%esp),%eax   # EAX <-- op1
	imull 8(%esp)       # EDX:EAX <-- EAX*op2

	idivl 12(%esp)       # EAX <-- EDX:EAX / modulus
	                    # EDX <-- EDX:EAX % modulus
	
	movl %edx,%eax      # return value goes into EAX
	ret

	
