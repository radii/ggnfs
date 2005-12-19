/**************************************************************/
/* intutils.c                                                 */
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
#include "intutils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifdef _MSC_VER
#pragma warning (disable: 4996) /* warning C4996: 'function' was declared deprecated */
#endif

/********************************************/
/* Compare arrays of two s32 values.        */
/********************************************/
int cmp2S32s(const void *a, const void *b)
{ 
	if (((s32 *)a)[0] < ((s32 *)b)[0]) return -1;
	if (((s32 *)a)[0] > ((s32 *)b)[0]) return 1;
	if (((s32 *)a)[1] < ((s32 *)b)[1]) return -1;
	if (((s32 *)a)[1] > ((s32 *)b)[1]) return 1;
	
	return 0;
}

/*********************************************************************/
/* Given a list of integers: list[0], ..., list[size-1], reduce it   */
/* so it contains only those integers that appeared an odd number    */
/* of times. That is, remove the elements that occur an even number  */
/* of times.                                                         */
/* Return value: The number of elements in the reduced list.         */
/*********************************************************************/
/* This is a crappy way to do it, but it's good for now.             */
int removeEvens(s32 *list, s32 size)
{ 
	s32 j, k, unique;

	if (size <= 1) 
		return size;

	qsort(list, size, sizeof(s32), cmpS32s);

	j = k = unique = 0;

	while (j < size) 
	{
		k = 1;
	
		while (((j+k) < size) && (list[j] == list[j+k]))
			k++;
	
		if (k & 0x01) 
		{
		    /* It occurs an odd number of times, so keep it. */
			list[unique++] = list[j];
		}
	
		j += k;
	}

	return unique;
}

/*************************************************/
/* Sort an array of s32s and remove duplicates.  */
/*************************************************/
s32 sortRMDups(s32 *L, s32 size)
{ 
	s32 i, unique;

	if (size <= 1) 
		return size;

	qsort(L, size, sizeof(s32), cmpS32s);

	i = 1; 
	unique = 0;

	while (i < size) 
	{
		if (L[i] == L[unique])
		{
			i++;
		}
	    else 
		{
			unique++;
			L[unique] = L[i];
			i++; 
		}
	}

	return unique + 1;
}

/**********************************************************/
/* Sort an array of pairs of s32s and remove duplicates. */
/**********************************************************/
s32 sortRMDups2(s32 *L, s32 size)
{ 
	s32 i, unique;

	if (size <= 1) 
		return size;

	qsort(L, size, 2*sizeof(s32), cmp2S32s);

	i = 1; 
	unique = 0;

	while (i < size) 
	{
	
		if ((L[2*i] == L[2*unique]) && (L[2*i+1] == L[2*unique+1]))
		{
			i++;
		}
		else 
		{ 
			unique++;
			L[2*unique] = L[2*i];
			L[2*unique+1] = L[2*i+1];
			i++;
		}
	}

	return unique + 1;
}
