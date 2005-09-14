;*************************************************************
;* mulmod32.s - adapted from mulmod32.s by Brian Gladman                                                 
;* Copyright 2003, Chris Monico.                              
;*************************************************************
;*  This file is part of GGNFS.
;*
;*   GGNFS is free software; you can redistribute it and/or modify
;*   it under the terms of the GNU General Public License as published by
;*   the Free Software Foundation; either version 2 of the License, or
;*   (at your option) any later version.
;*
;*   GGNFS is distributed in the hope that it will be useful,
;*   but WITHOUT ANY WARRANTY; without even the implied warranty of
;*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;*   GNU General Public License for more details.
;*
;*   You should have received a copy of the GNU General Public License
;*   along with GGNFS; if not, write to the Free Software
;*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    section .text ; use32
	align	4
	global	_mulmod32
	
_mulmod32:
	mov		eax,[esp+4]
	imul	dword [esp+8]
	idiv	dword [esp+12]
	mov		eax,edx
	ret

	
