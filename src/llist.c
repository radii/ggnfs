/**************************************************************/
/* llist.c                                                    */
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

/* Routines for efficient operation on (large) arrays
   with a variable number of s32s in each entry.
*/

#ifdef _MSC_VER
#pragma warning (disable: 4996) /* warning C4996: 'function' was declared deprecated */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "ggnfs.h"


#define MAXFIELDENTRIES 2048


/*************************************************************/
int ll_init(llist_t *L, s32 maxDataSize, s32 maxFields)
/*************************************************************/
{
  L->data = L->index = NULL;
  L->maxDataSize = L->maxFields = 0;
  L->numFields = 0;
  if (!(L->data = (s32 *)lxmalloc(maxDataSize*sizeof(s32),0))) {
    fprintf(stderr, "ll_init() Memory allocation error (%" PRIu32 "Mb).\n", 
            (u32)(maxDataSize*sizeof(s32)/(1024*1024)) );
    return -1;
  }
  if (!(L->index = (s32 *)lxmalloc((maxFields+1)*sizeof(s32),0))) {
    fprintf(stderr, "ll_init() Memory allocation error (%" PRIu32 "Mb).\n", 
            (u32)(maxFields*sizeof(s32)/(1024*1024)) );
    free(L->data);
    L->data=NULL;
    return -1;
  }
  L->maxDataSize = maxDataSize;
  L->maxFields = maxFields;
  L->index[0]=0;
  return 0;
}

/*************************************************************/
void ll_clear(llist_t *L)
/*************************************************************/
{
  if (L->data) free(L->data);
  if (L->index) free(L->index);
  L->data = L->index = NULL;
  L->maxDataSize = 0;
  L->maxFields = 0;
  L->numFields = 0;
}

/*************************************************************/
void ll_resize(llist_t *L, s32 maxDataSize)
/*************************************************************/
{ s32 *tmp;

  tmp = lxrealloc(L->data, maxDataSize*sizeof(s32),0);
  if (tmp != NULL) {
    L->data = tmp;
  } else {
    printf("ll_resize() Memory (re-)allocation error!\n");
  }
  L->maxDataSize = maxDataSize;
}

#define RESIZE 1024*1024*4
/*********************************************************/
int ll_appendField(llist_t *L, s32 *entries, s32 numEntries)
/*********************************************************/
{ s32 cSize = L->index[L->numFields];
  s32 *tmpPtr, i;

  if (cSize + numEntries + 2 >= L->maxDataSize) {
    tmpPtr = (s32 *)lxrealloc(L->data, (L->maxDataSize + numEntries + RESIZE)*sizeof(s32),0);
    if (tmpPtr != NULL) L->data = tmpPtr;
    else {
      fprintf(stderr, "ll_appendField() Memory (re-)allocation error!\n");
      return -1;
    }
    L->maxDataSize += numEntries + RESIZE;
  }
  if (L->numFields + 3 >= L->maxFields) {
    tmpPtr = (s32 *)lxrealloc(L->index, (L->maxFields + RESIZE/4)*sizeof(s32),0);
    if (tmpPtr != NULL) L->index = tmpPtr;
    else {
      fprintf(stderr, "ll_appendField() Memory (re-)allocation error!\n");
      return -1;
    }
    L->maxFields += RESIZE/4;
  }
  for (i=0; i<numEntries; i++)
    L->data[cSize+i] = entries[i];
  L->numFields += 1;
  L->index[L->numFields] = cSize+numEntries;
  return 0;
}
  

/*********************************************************/
int ll_verify(llist_t *L)
/*********************************************************/
/* Check that the llist_t appears to be valid.           */
/*********************************************************/
{ s32 i, s0, s1;

  if ((L->numFields < 0) || (L->numFields > L->maxFields)) {
    printf("ll_verify() L->numFields=%" PRId32 " vs. L->maxFields=%" PRId32 ".\n",
           L->numFields, L->maxFields);
    return -1;
  }

  for (i=0; i<L->numFields; i++) {
    s0 = L->index[i];
    s1 = L->index[i+1];
    if (s1 > L->maxDataSize) {
      printf("ll_verify() L->index[%" PRId32 "]=%" PRId32 " vs. L->maxDataSize=%" PRId32 ".\n",
             i+1, s1, L->maxDataSize);
      return -1;
    }
    if (s1 < s0) {
      printf("ll_verify() L->index[%" PRId32 "]=%" PRId32 " vs. L->index[%" PRId32 "+1]=%" PRId32 ".\n",
             i, s0, i, s1);
      return -1;
    }
  }
  return 0;
}
      


/*********************************************************/
int ll_deleteFields(llist_t *L, s32 *fields, s32 numFields)
/*********************************************************/
/* Delete the specified fields from the list. The        */
/* remaining fields will be shifted down to fill in the  */
/* gaps.                                                 */
/* `fields' should be sorted and unique!!!!!!!           */
/*********************************************************/
{ s32 destField, c0, numMove, dataSize, shift, j, m=numFields;
  int  i;
                                                                                   
  destField=fields[0];

  /* Sanity check: */
  if ((L->numFields > L->maxFields) || (L->index[L->numFields] > L->maxDataSize)) {
    printf("ll_deleteFields(): 'L' appears corrupt:\n");
    exit(-1);
  }                                                                                   
  for (i=0; i<m; i++) {
    c0 = fields[i]+1;
    numMove = (i<(m-1)) ? (fields[i+1]-c0):(L->numFields - c0);
    /* Now, move the numMove fields beginning at c0  */
    /* to destField.                                 */
    if (numMove > 0) {
      /* First the data. */
      dataSize = L->index[c0+numMove] - L->index[c0];
      memmove(&L->data[L->index[destField]], &L->data[L->index[c0]], dataSize*sizeof(s32));
      /* Now fix the index fields. */
      shift = L->index[c0] - L->index[destField];
      for (j=0; j<=numMove; j++)
        L->index[c0+j] -= shift;
      for (j=0; j<=numMove; j++)
        L->index[destField+j] = L->index[c0+j];
      destField += numMove;
    }
  }
  L->numFields -= m;
  return 0;
}

/************************************************************/
int ll_numCommonEntries(llist_t *L, s32 field0, s32 field1)
/************************************************************/
/* return the number of common entries in these two fields. */
/* The fields are assumed to have a fairly small number of  */
/* entries (i.e., say, less than 20 or so).                 */
/************************************************************/
{ s32 i, j, num=0;

  for (i=L->index[field0]; i<L->index[field0+1]; i++)
    for (j=L->index[field1]; j<L->index[field1+1]; j++)
      if (L->data[i]==L->data[j])
        num++;
  return num;
}
  


/*********************************************************/
int ll_catFields(llist_t *L, lpair_t *pairs, s32 numPairs, int mod2)
/*********************************************************/
/* For each f1=pairs[i].x, f2=pairs[i].y, concatenate the*/
/* data from field f2 onto field f1. If `mod2' is        */
/* nonzero, pairs of entries will be assumed to cancel,  */
/* so there is at most one copy of each resulting data   */
/* element in the field f1.                              */
/*   This is tricky, because some of the data may shift  */
/* in one direction while some shifts in the opposite.   */
/* But we handle this with the elegance it deserves.     */
/*    It is assumed that the list of pairs is sorted     */
/* ascending on the 'x' field, and the 'x' fields are    */
/* unique! (i.e., they should be strictly increasing).   */
/*********************************************************/
{ s32     i, j, c0, c1, shift, maxShift, indexShift;
  s32     newEntries[MAXFIELDENTRIES], numNewEntries=0, numOldEntries;
  s32     sIndex, dIndex, moveSize;
  s32    *tmpPtr;

  if (numPairs <= 0) return 0;

  /* Compute the distance that the fields must be shifted. */
  shift=maxShift=0;
  for (i=0; i<numPairs; i++) {
    c0 = pairs[i].x;
    c1 = pairs[i].y;
    numNewEntries=0;
    for (j=L->index[c0]; (j<L->index[c0+1]) && (numNewEntries < MAXFIELDENTRIES); j++)
      newEntries[numNewEntries++] = L->data[j];
    for (j=L->index[c1]; (j<L->index[c1+1]) && (numNewEntries < MAXFIELDENTRIES); j++)
      newEntries[numNewEntries++] = L->data[j];
    if (numNewEntries >= MAXFIELDENTRIES) {
      printf("MAXFIELDENTRIES exceeded (i=%" PRId32 ", c0=%" PRId32 ", c1=%" PRId32 ".). Ignoring...\n",i,c0,c1);
      return -1;
    }
    if (mod2)
      numNewEntries = removeS32Pairs(newEntries, numNewEntries);
    shift +=  numNewEntries - (L->index[c0+1] - L->index[c0]);
    if (labs(shift) > maxShift) 
      maxShift = labs(shift);
  }
  maxShift+=1;

  /* The idea is this: Shift everything to the right by the maxShift.
     Then, we can work by simply shifting everything left
     one block at a time and not have to worry about handling
     situations where some data is overwritten. The cost is nothing:
     each byte of data will be moved twice, instead of once. But moving
     it only once would come at a huge coding complexity cost, so this
     is probably even better.
  */
  if (L->index[L->numFields] + maxShift >= L->maxDataSize) {
    /* We need to do a reallocation: */
    tmpPtr = realloc(L->data, (L->index[L->numFields] + maxShift + 10)*sizeof(s32));
    if (tmpPtr == NULL) {
      printf("ll_catFields(): memory reallocation error!\n");
      printf("%" PRId32 " s32's requested.\n", (L->index[L->numFields] + maxShift + 10));
      printf("Old size was L->maxDataSize=%" PRId32 "\n", L->maxDataSize);
      printf("L->numFields = %" PRId32 ", maxShift = %" PRId32 "\n", L->numFields, maxShift);
      printf("L->index[%" PRId32 "] = %" PRId32 "\n", L->numFields, L->index[L->numFields]);
      printf("numPairs=%" PRId32 ", numNewEntries=%" PRId32 "\n", numPairs, numNewEntries);
#ifdef MALLOC_REPORTING
      printf("malloc useage is about %d MB.", mallocReport());
#endif
      free(pairs); return -1;
    }
    L->data = tmpPtr;
    L->maxDataSize = L->index[L->numFields] + maxShift + 10;
  }
  /* Move everything right by `maxShift'. Of course, there's no need 
     to move the first block (i.e., fields left of the first field
     to be changed.)
  */
  sIndex = L->index[pairs[0].x+1];
  dIndex = sIndex + maxShift;
  moveSize = L->index[L->numFields] - sIndex;
  memmove(&L->data[dIndex], &L->data[sIndex], moveSize*sizeof(s32));

  /* Now, do the field additions, moving blocks back to the left as we go: */
  indexShift = 0;
  for (i=0; i<numPairs; i++) {
    /* Right now, the fields in [0, pairs[i].x] are in proper position,
       starting at L->data[0], and they all have correct L->index values.
       The fields in [pairs[i].x + 1, ..., L->numFields) are all shifted
       in memory by `maxShift', and their corresponding L->index values
       are off by that amount as well.
    */

    c0 = pairs[i].x; c1 = pairs[i].y;
    /* Compute the new entries for column c0: */
    numNewEntries=0;
    for (j=L->index[c0]; j<L->index[c0+1]+indexShift; j++)
      newEntries[numNewEntries++] = L->data[j];
    numOldEntries = numNewEntries;
    if (c1 <=c0) {
      for (j=L->index[c1]; j<L->index[c1+1]; j++)
        newEntries[numNewEntries++] = L->data[j];
    } else {
      for (j=L->index[c1]; j<L->index[c1+1]; j++)
        newEntries[numNewEntries++] = L->data[j+maxShift];
    }
    if (mod2)
      numNewEntries = removeS32Pairs(newEntries, numNewEntries);
    /* Place the new entries: */
    for (j=0; j<numNewEntries; j++)
      L->data[L->index[c0]+j] = newEntries[j];

    /* Move the next block of data and adjust it's indicies. */
    indexShift += numNewEntries - numOldEntries;
    sIndex = L->index[c0+1] + maxShift;
    dIndex = L->index[c0]+numNewEntries;
    if (i<(numPairs-1)) {
      moveSize = L->index[pairs[i+1].x+1] - sIndex + maxShift;
      memmove(&L->data[dIndex], &L->data[sIndex], moveSize*sizeof(s32));
      for (j=c0+1; j<=pairs[i+1].x; j++) 
        L->index[j] += indexShift;
    } else {
      moveSize = L->index[L->numFields] - sIndex + maxShift;
      memmove(&L->data[dIndex], &L->data[sIndex], moveSize*sizeof(s32));
      for (j=c0+1; j<=L->numFields; j++) 
        L->index[j] += indexShift;
    }
  }  
  return 0;
}


/************************************************************/
/* These are for functions which need access to a llist_t,  */
/* but for some reason or another, it's not convenient to   */
/* pass it as a pointer to the function.                    */
/************************************************************/
static llist_t *llistPtr=NULL;
void ll_setGlobal(llist_t *L)
{ llistPtr = L; }
/****************************************************/
int ll_cmpFieldSize(const void *A, const void *B)
/****************************************************/
/* Uses the global list set by ll_setGlobal().      */
/****************************************************/
{ s32 *a=(s32 *)A, *b=(s32 *)B;
  s32 s1, s2;

  s1 = llistPtr->index[*a+1] - llistPtr->index[*a];
  s2 = llistPtr->index[*b+1] - llistPtr->index[*b];
  if (s1 < s2) return -1;
  if (s1 > s2) return 1;
  return 0;
}

/***************************************************************/
int ll_getsortOnFieldSize(s32 *fields, llist_t *L)
/***************************************************************/
/* `fields' is allocated for at least L->numFields s32s.       */
/* On output, fields[0], fields[1],... will be indicies to the */
/* fields of 'L', sorted in ascending order according to the   */
/* number of entries in each field.                            */
/***************************************************************/
/* This function is only used for sorting relations on the     */
/* number of large primes right now. If it is used for         */
/* something else later, this might need to be increased:      */
/***************************************************************/
#define MAX_LIST_ENTRIES 10
{ s32 i;
  s32 counts[MAX_LIST_ENTRIES];
  s32 c_idx;

  memset(counts,0,sizeof(s32)*MAX_LIST_ENTRIES);
  for (i=0; i<L->numFields; ++i) {
    c_idx=L->index[i+1]-L->index[i];
    if (c_idx<MAX_LIST_ENTRIES) {
      ++counts[c_idx+1];
    } else {
      printf("ll_getsortOnFieldSize: Severe: MAX_LIST_ENTRIES too small.\n");
      return -1;
    }
  }
  for (i=1; i<MAX_LIST_ENTRIES; ++i)
    counts[i]+=counts[i-1];

  for (i=0; i<L->numFields; ++i) {
    c_idx=L->index[i+1]-L->index[i];
    if (c_idx<MAX_LIST_ENTRIES) {
      fields[counts[c_idx]++]=i;
    } else {
      printf("ll_getsortOnFieldSize: Severe: MAX_LIST_ENTRIES too small.\n");
      return -1;
    }
  }
  return 0;
}

/***************************************************************/
int ll_write(char *fname, llist_t *C)
/***************************************************************/
{ FILE *fp;

  if (!(fp = fopen(fname, "wb"))) {
    fprintf(stderr, "Error opening %s for write!\n", fname);
    return -1;
  }
  fwrite(&C->numFields, sizeof(s32), 1, fp);
  fwrite(C->index, sizeof(s32), C->numFields+1, fp);
  fwrite(C->data, sizeof(s32), C->index[C->numFields], fp);
  fclose(fp);
  return 0;
}

/***************************************************************/
int ll_read(llist_t *C, char *fname)
/***************************************************************/
{ FILE *fp;

  if (!(fp = fopen(fname, "rb"))) {
    fprintf(stderr, "Error opening %s for read!\n", fname);
    return -1;
  }
  fread(&C->numFields, sizeof(s32), 1, fp);
  C->maxFields = C->numFields;

  C->index = (s32 *)lxmalloc(sizeof(s32)*(C->numFields+1),1);
  fread(C->index, sizeof(s32), C->numFields+1, fp);

  C->maxDataSize = C->index[C->numFields];
  C->data = (s32 *)lxmalloc(sizeof(s32)*(C->maxDataSize),1);
  fread(C->data, sizeof(s32), C->index[C->numFields], fp);
  fclose(fp);
  return 0;
}



