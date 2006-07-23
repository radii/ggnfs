/**************************************************************/
/* rellist.c                                                  */
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
#include "rellist.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifdef _MSC_VER
#pragma warning (disable: 4996) /* warning C4996: 'function' was declared deprecated */
#endif

/******************************************************/
/* Allocate 'RL' so it can hold the largest of the    */
/* processed relation files.                          */
/******************************************************/
int allocateRelList(multi_file_t *prelF, rel_list *RL)
{ 
	off_t maxSize;
	char prelName[512];
	int  i;
	struct stat fileInfo;

	maxSize = 0;

	for (i = 0; i < prelF->numFiles; i++) 
	{
	
		sprintf(prelName, "%s.%d", prelF->prefix, i);
	
		if (stat(prelName, &fileInfo) == 0) 
			maxSize = MAX(maxSize, fileInfo.st_size);
	}

	RL->numRels = 0;
	RL->maxDataSize = 1000 + maxSize/sizeof(s32);

	if (!(RL->relData = (s32 *)lxmalloc(RL->maxDataSize * sizeof(s32), 0))) 
	{
		fprintf(stderr, "Error allocating %" PRIu32 "MB for processed relation files!\n",
			(u32)(RL->maxDataSize * sizeof(s32)/1048576) );
		
		fprintf(stderr, "Try decreasing DEFAULT_MAX_FILESIZE and re-running.\n");
		exit(-1);
	}

	/* Again: it's a safe bet that any relation needs at least 20 s32s, so: */
	RL->maxRels = (u32)RL->maxDataSize/20;

	if (!(RL->relIndex = (s32 *)lxmalloc(RL->maxRels * sizeof(s32), 0))) 
	{
		fprintf(stderr, "Error allocating %" PRIu32 "MB for relation pointers!\n",
			(u32)(RL->maxRels * sizeof(s32)/1048756) );
	
		free(RL->relData);
		exit(-1);
	}
	
	return 0;
}

/*********************************************************************/
/* Allocate for and read in the specified relation file. Caller is   */
/* obviously responsible for freeing the memory when done!           */
/*********************************************************************/
rel_list *getRelList(multi_file_t *prelF, int index)
{ 
	rel_list *RL;
	FILE     *fp;
	char      fName[256];
	struct stat fileInfo;

	RL = (rel_list *)lxmalloc(sizeof(rel_list), 1);
	RL->maxDataSize = 0;
	sprintf(fName, "%s.%d", prelF->prefix, index);

	if (stat(fName, &fileInfo)) 
	{
		printf("Could not stat file %s!\n", fName);
		free(RL); return NULL;
	}

	RL->maxDataSize = 4096 + fileInfo.st_size/sizeof(s32);

	if ((fp = fopen(fName, "rb"))) 
	{
		readRaw32(&RL->maxRels, fp);
		fclose(fp);
	}

	RL->maxRels += 5;

	/* Now allocate for the relations. */
	if (!(RL->relData = (s32 *)lxmalloc(RL->maxDataSize * sizeof(s32), 0))) 
	{
		fprintf(stderr, "Error allocating %" PRIu32 "MB for reading relation list!\n",
			(u32)(RL->maxDataSize * sizeof(s32)/1048576) );
		free(RL); 
		return NULL;
	}

	if (!(RL->relIndex = (u32 *)lxmalloc(RL->maxRels * sizeof(s32), 0))) 
	{
		fprintf(stderr, "Error allocating %" PRIu32 "MB for relation pointers!\n", 
			(u32)(RL->maxRels * sizeof(s32)/1048576) );
		free(RL->relData); 
		free(RL);
		return NULL;
	}

	RL->numRels = 0;
	readRelList(RL, fName);
	return RL;
}

/*********************************************************************/
/* Clear relations list.                                             */
/*********************************************************************/
void clearRelList(rel_list *RL)
{
	if (RL->relData != NULL) 
		free(RL->relData);

	if (RL->relIndex != NULL)
		free(RL->relIndex);

	RL->relData = RL->relIndex = NULL;
	RL->maxDataSize = RL->maxRels = 0;
}

/************************************************************************/
/* removeFrac should be a fraction in [0,1). This function will remove  */
/* the heaviest removeFrac relations from the processed relation files  */
/* file, appending them in siever-output format to the file appendName. */
/************************************************************************/
void pruneRelLists(multi_file_t *prelF, char *appendName, double removeFrac, nfs_fb_t *FB, int short_form)
#define MAX_PR_SIZE 2048
{ 
	rel_list *RL;
	int       i;
	u32       j;
	long      numBySize[MAX_PR_SIZE], numRemove, t;
	s32       size;  
	char      outputStr[1024], str[128];
	s32       bufMax, bufIndex, *buf;
	relation_t R;
	FILE     *afp, *rfp;

	/* Prep the output buffers: */
	bufMax = IO_BUFFER_SIZE;
	if (!(buf = (s32 *)malloc(bufMax*sizeof(s32)))) 
	{
		printf("pruneRelLists() : Memory allocation error for buf!\n");
		exit(-1);
	}

	afp = fopen(appendName, "a");
	
	for (i = 0; i < prelF->numFiles; i++) 
	{
		printf("Pruning rel file %d/%d...", i+1, prelF->numFiles);
		fflush(stdout);
		RL = getRelList(prelF, i);
		memset(numBySize, 0x00, MAX_PR_SIZE*sizeof(long));

		for (j = 0; j < RL->numRels; j++) 
		{
			size = RL->relIndex[j+1] - RL->relIndex[j];
			size = MIN(size, MAX_PR_SIZE);
			numBySize[size] += 1;
		}
		
		numRemove = (long)((removeFrac*RL->numRels)) - 1;
		
		printf("from %" PRId32 " to %ld relations...",RL->numRels, RL->numRels-numRemove);
		
		/*
		 * We proceed by modifying numBySize so that it contains the
		 * number of relations we will keep with the given size.
		 */
		t = numRemove;
		
		for (j = MAX_PR_SIZE-1; (j >=0 ) && (t > 0); j--) 
		{
			if (t >= numBySize[j]) 
			{
				t -= numBySize[j];
				numBySize[j] = 0;
			} 
			else 
			{
				numBySize[j] -= t;
				t = 0;
			}
		}

		/* Now, go through the relations and decide which to keep and which to dump. */
		rfp = fopen(".tmp", "w");
		size = RL->numRels - numRemove;
		writeRaw32(rfp, &size);
		bufIndex = 0;
		
		for (j = 0; j < RL->numRels; j++) 
		{
			size = RL->relIndex[j+1] - RL->relIndex[j];
			if (numBySize[size] <= 0) 
			{
				/* We need to dump and remove this relation. */
				dataConvertToRel(&R, &RL->relData[RL->relIndex[j]]);
				makeOutputLine(outputStr, &R, FB, short_form);
				fprintf(afp, "%s\n", outputStr);
			} 
			else
			{
				numBySize[size] -= 1;
				/* And re-record this relation. */
				if (bufIndex + 1024 < bufMax) 
				{
					memcpy(buf+bufIndex, &RL->relData[RL->relIndex[j]], size*sizeof(s32));
					bufIndex += size;
				} 
				else 
				{
					fwrite(buf, sizeof(s32), bufIndex, rfp);
					bufIndex = 0;
					memcpy(buf+bufIndex, &RL->relData[RL->relIndex[j]], size*sizeof(s32));
					bufIndex += size;
				}
			}
		}

		fwrite(buf, sizeof(s32), bufIndex, rfp);
		fclose(rfp);
		sprintf(str, "%s.%d", prelF->prefix, i);
		remove(str);
		rename(".tmp", str);
		clearRelList(RL);
	}

	free(buf);
}
