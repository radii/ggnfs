/**************************************************************/
/* intutils.h                                                 */
/* Copyright 2004, Chris Monico.                              */
/* Initial revision by Sten.                                  */
/* Several common integer functions are moved here.           */
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
#ifndef __INTUTILS_H__
#define __INTUTILS_H__
#include "ggnfs.h"

#if defined (__cplusplus)
extern "C" {
#endif

/********************************************/
/* Compare arrays of two s32 values.        */
/********************************************/
int cmp2S32s(const void *a, const void *b);

/*********************************************************************/
/* Given a list of integers: list[0], ..., list[size-1], reduce it   */
/* so it contains only those integers that appeared an odd number    */
/* of times. That is, remove the elements that occur an even number  */
/* of times.                                                         */
/* Return value: The number of elements in the reduced list.         */
/*********************************************************************/
int removeEvens(s32 *list, s32 size);

/*************************************************/
/* Sort an array of s32s and remove duplicates.  */
/*************************************************/
s32 sortRMDups(s32 *L, s32 size);

/**********************************************************/
/* Sort an array of pairs of s32s and remove duplicates. */
/**********************************************************/
s32 sortRMDups2(s32 *L, s32 size);

#if defined (__cplusplus)
};
#endif

#endif /* __INTUTILS_H__ */
