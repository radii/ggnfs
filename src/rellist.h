/**************************************************************/
/* rellist.h                                                  */
/* Copyright 2004, Chris Monico.                              */
/* Initial revision by Sten.                                  */
/* Several common relation functions are moved here.          */
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
#ifndef __RELLIST_H__
#define __RELLIST_H__
#include "ggnfs.h"

#if defined (__cplusplus)
extern "C" {
#endif

/* These are in s32's. */
#define IO_BUFFER_SIZE  2*1024*1024

/******************************************************/
/* Allocate 'RL' so it can hold the largest of the    */
/* processed relation files.                          */
/******************************************************/
int allocateRelList(multi_file_t *prelF, rel_list *RL);

/*********************************************************************/
/* Allocate for and read in the specified relation file. Caller is   */
/* obviously responsible for freeing the memory when done!           */
/*********************************************************************/
rel_list *getRelList(multi_file_t *prelF, int index);

/*********************************************************************/
/* Clear relations list.                                             */
/*********************************************************************/
void clearRelList(rel_list *RL);

/************************************************************************/
/* removeFrac should be a fraction in [0,1). This function will remove  */
/* the heaviest removeFrac relations from the processed relation files  */
/* file, appending them in siever-output format to the file appendName. */
/************************************************************************/
void pruneRelLists(multi_file_t *prelF, char *appendName, double removeFrac, nfs_fb_t *FB);

#if defined (__cplusplus)
};
#endif

#endif /* __RELLIST_H__ */
