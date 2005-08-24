/**************************************************************/
/* matstuff.c                                                 */
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifndef _MSC_VER 
#include <sys/time.h>
#endif
#include "ggnfs.h"

/* The matrix is assumed to be setup so that each column corresponds
   to a prime from a factor base (or the QCB). Thus, the columns
   correspond to the relations.
*/

/* Max weight of a column: */
#define MAXCOLWT MAX_ROWS_IN_COL

/* We need a structure for recalling the changes we make
   during pruning, so that when we find some dependencies
   for the pruned matrix, we can recover corresponding
   dependencies on the original matrix. This structure is 
   such that columns
       origCol[colIndex[j]], ..., origCol[colIndex[j+1]-1]
   is the list of original columns contributing to column j
   of the new matrix.
*/


/* A quick way to test if some particular row is in one of
   the dense blocks.
*/
s32 *isDenseRow=NULL;
#define ISDENSEROW(_j) (isDenseRow[_j/32]&BIT(_j%32))

/*********************************************************/
int mat_verify(nfs_sparse_mat_t *M)
/*********************************************************/
/* Check that M appears to be valid.                     */
/*********************************************************/
{ s32 i, s0, s1;

  if (M->numCols <= 0) {
    printf("mat_verify() M->numCols=%" PRId32 ".\n", M->numCols);
    return -1;
  }

  for (i=0; i<M->numCols; i++) {
    s0 = M->cIndex[i];
    s1 = M->cIndex[i+1];
    if (s1 > M->maxDataSize) {
      printf("mat_verify() M->cIndex[%" PRId32 "]=%" PRId32 " vs. M->maxDataSize=%" PRId32 ".\n",
             i+1, s1, M->maxDataSize);
      return -1;
    }
    if (s1 < s0) {
      printf("mat_verify() M->cIndex[%" PRId32 "]=%" PRId32 " vs. M->cIndex[%" PRId32 "+1]=%" PRId32 ".\n",
             i, s0, i, s1);
      return -1;
    }
  }
  return 0;
}

/*********************************************************/
void setDenseRows(nfs_sparse_mat_t *M)
/*********************************************************/
{ s32 i, j, n, offset;

  /* Setup the strucure for recognizing dense rows. */
  n = 1+M->numRows/32;
  if (isDenseRow != NULL) free(isDenseRow);
  isDenseRow = (s32 *)calloc(n, sizeof(s32));
  for (i=0; i<M->numDenseBlocks; i++) {
    offset = M->denseBlockIndex[i];
    for (j=0; j<64; j++) 
      isDenseRow[(offset+j)/32] |= BIT((offset+j)%32);
  }
}

/*********************************************************/
int cm_init(llist_t *C, nfs_sparse_mat_t *M)
/*********************************************************/
{ s32 i, maxDataSize, maxCols;

  maxDataSize = 1.5*M->numCols;
  maxCols = M->numCols;
  if (ll_init(C, maxDataSize, maxCols))
    exit(-1);
  for (i=0; i<M->numCols; i++) {
    C->data[i] = i; C->index[i] = i;
  }
  C->index[i] = i; /* Da Terminator. */
  C->numFields = M->numCols;
  return 0;
}
  
/*********************************************************/
int cm_removeCols(llist_t *C, s32 *cols, s32 numCols)
/*********************************************************/
/* Adjust C to compensate for deleting the given columns */
/* from the matrix.                                      */
/* `cols' should be sorted, but we will sort it just in  */
/* case.                                                 */
/*********************************************************/
{ 
  return ll_deleteFields(C, cols, numCols);
}

/*********************************************************/
int cm_addCols_par(llist_t *C, lpair_t *pairs, s32 numCols)
/*********************************************************/
/* In the new matrix, add column `srcCol' to column      */
/* destCol.                                              */
/* The values of destCol must be distinct!!!!            */
/* This is tricky, because some of the data may shift    */
/* in one direction while some shifts in the opposite.   */
/* But we handle this with the elegance it deserves.     */
/*********************************************************/
{ 
  return ll_catFields(C, pairs, numCols, 1);
}


/*********************************************************/
int removeRows(nfs_sparse_mat_t *M, s32 *r, int m)
/*********************************************************/
/* Remove rows r[0], r[1],..., r[m-1] from the matrix M. */
/* These rows must already be empty or the matrix will   */
/* be corrupted! (i.e., the rows numbers r[0],... must   */
/* not appear anywhere in M->cEntry[].                   */
/*********************************************************/
{ s32 j, i, rT, rj;
  int  good;
  s32 *rowMap;

  rowMap=(s32 *)malloc(M->numRows*sizeof(s32));
  for (i=0; i<M->numRows; i++)
    rowMap[i]=i;

  rT = M->numRows-64;
  for (j=0; j<m; j++) {

    rj = r[j];
    if (ISDENSEROW(rj)) {
      fprintf(stderr, "removeRows() Error: attempt to delete a dense row!\n");
      free(rowMap); return -1;
    }
    rT--;
    do {
      good=1;
      for (i=0; i<m; i++) {
        if (rT==r[i]) {
          good=0;
          rT--;
        }
      }
    } while (!good);
    /* This shouldn't happen unless an awful lot of rows were deleted,
       or the matrix had some sparse rows, other then QCB, near the end
       of the matrix. We should think about fixing this code to handle
       that case as well, but it can wait.
    */
    if (ISDENSEROW(rT)) return -1;
    if (rT <= rj)
      break;
    /* Replace all occurrences of rT with rj. */
    rowMap[rT]=rj;
  }

  /* Finally, apply the map. */
  for (i=M->cIndex[M->numCols]-1; i>=0; i--)
    M->cEntry[i] = rowMap[M->cEntry[i]]; 

  /* We will definitely need to adjust the indicies for the QCB blocks: */
  M->denseBlockIndex[0] -= m;
  /* There shouldn't be any others dense blocks whose indicies need to be
     adjusted. If there were, the code above would have bailed out.
  */
  M->numRows -= m;
  setDenseRows(M);
  free(rowMap);
  return 0;
}

/******************************************************/
int removeCols(nfs_sparse_mat_t *M, llist_t *C, s32 *cols, s32 m)
/******************************************************/
{ s32 destCol, c0, numMove, dataSize, shift, j;
  int  i;

  if (m <= 0) return 0;
  qsort(cols, m, sizeof(s32), cmpS32s);
  for (i=1; i<m; i++)
    if (cols[i-1]==cols[i]) {
      printf("Error: removeCols(): duplicate columns!\n");
      exit(-1);
    }

  destCol=cols[0];
  for (i=0; i<m; i++) {
    c0 = cols[i]+1;
    numMove = (i<(m-1)) ? (cols[i+1]-c0):(M->numCols - c0);
    /* Now, move the numMove columns beginning at c0 */
    /* to destCol.                                   */
    if (numMove > 0) {
      /* First the data. */
      dataSize = M->cIndex[c0+numMove] - M->cIndex[c0];
      memmove(&M->cEntry[M->cIndex[destCol]], &M->cEntry[M->cIndex[c0]], dataSize*sizeof(s32));
      /* Now fix the cIndex fields. */
      shift = M->cIndex[c0] - M->cIndex[destCol];
      for (j=0; j<=numMove; j++)
        M->cIndex[c0+j] -= shift;
      for (j=0; j<=numMove; j++)
        M->cIndex[destCol+j] = M->cIndex[c0+j];

      /* Move the dense words. */
      for (j=0; j<M->numDenseBlocks; j++)
        memmove(&M->denseBlocks[j][destCol], &M->denseBlocks[j][c0], numMove*sizeof(u64));
      destCol += numMove;
    }
  }
  /* Update the colMap. */
  cm_removeCols(C, cols, m);
  M->numCols -= m;
  return 0;
}


/*********************************************************/
int addCols_par(nfs_sparse_mat_t *M, llist_t *C, s32 *destCol, s32 *srcCol, s32 _numCols)
/*********************************************************/
/* In the new matrix, add column `srcCol' to column      */
/* destCol.                                              */
/* The values of destCol must be distinct!!!!            */
/* This is tricky, because some of the data may shift    */
/* in one direction while some shifts in the opposite.   */
/* But we handle this with the elegance it deserves.     */
/* We will take some liberties here, and the caller may  */
/* not get all the additions asked for. In particular,   */
/* we make sure that the destCol's are distinct and not  */
/* appearing as any srcCol.                              */
/*********************************************************/
{ lpair_t *pairs;
  s32     i, j, k, c0, c1, shift, maxShift, indexShift, numCols=_numCols;
  s32     newEntries[MAXCOLWT], numNewEntries, numOldEntries;
  s32     sIndex, dIndex, moveSize;
  s32    *tmpPtr;

  if (numCols <= 0) return 0;

  /* Setup and sort the information: */
  pairs = (lpair_t *)malloc((numCols+1)*sizeof(lpair_t));
  for (i=0; i<numCols; i++) {
    pairs[i].x = destCol[i];
    pairs[i].y = srcCol[i];
  }
  pairs[i].x = M->numCols-1; /* For convenience - makes the bottom loop cleaner. */
  qsort(pairs, numCols, sizeof(lpair_t), cmp_lpair_t);

  /* We will take some liberties here to make the list
     of destCol's unique. In fact, we need somewhat more: to avoid
     any ambiguity, we need to make sure that any column appearing
     as a destination does not appear as a source (since we are doing this
     in parallel, it would very much complicate things!)
  */
  tmpPtr = (s32 *)malloc(numCols*sizeof(s32));
  memcpy(tmpPtr, srcCol, numCols*sizeof(s32));
  qsort(tmpPtr, numCols, sizeof(s32), cmpS32s);
  k=0;
  for (j=0,i=1; i<numCols; i++) {
    if (pairs[i].x != pairs[j].x) {
      while ((k<numCols) && (tmpPtr[k] < pairs[i].x))
        k++;
      if ((k>=numCols) || (pairs[i].x != tmpPtr[k])) {
        pairs[++j].x = pairs[i].x;
        pairs[j].y = pairs[i].y;
      }
    }
  }
  numCols = j+1;
  free(tmpPtr);

  /* Compute the distance that the columns must be shifted. */
  shift=maxShift=0;
  for (i=0; i<numCols; i++) {
    c0 = pairs[i].x;
    c1 = pairs[i].y;
    numNewEntries=0;
    for (j=M->cIndex[c0]; (j<M->cIndex[c0+1]) && (numNewEntries < MAXCOLWT); j++)
      newEntries[numNewEntries++] = M->cEntry[j];
    for (j=M->cIndex[c1]; (j<M->cIndex[c1+1])  && (numNewEntries < MAXCOLWT); j++)
      newEntries[numNewEntries++] = M->cEntry[j];
    if (numNewEntries >= MAXCOLWT) {
      printf("MAXCOLWT exceeded (i=%" PRId32 ", c0=%" PRId32 ", c1=%" PRId32 ".). Ignoring...\n",i,c0,c1);
      printf("         cIndex[c0]=%" PRId32 ", cIndex[c0+1] = %" PRId32 ".\n", M->cIndex[c0], M->cIndex[c0+1]);
      printf("         cIndex[c1]=%" PRId32 ", cIndex[c1+1] = %" PRId32 ".\n", M->cIndex[c1], M->cIndex[c1+1]);
      return -1;
    }
    numNewEntries = removeS32Pairs(newEntries, numNewEntries);
    shift +=  numNewEntries - (M->cIndex[c0+1] - M->cIndex[c0]);
    maxShift = MAX(labs(shift), maxShift);
    /* pairs[i].w is how much the data in columns
         [pairs[i].x+1, ..., pairs[i+1].x]
       must be shifted by from its current location.
    */
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
  if (M->cIndex[M->numCols] + maxShift >= M->maxDataSize) {
    /* We need to do a reallocation: */
    tmpPtr = realloc(M->cEntry, (M->cIndex[M->numCols] + maxShift + 10)*sizeof(s32));
    if (tmpPtr == NULL) {
      printf("addCols_par(): memory reallocation error!\n");
      free(pairs); return -1;
    }
    M->cEntry = tmpPtr;
    M->maxDataSize = M->cIndex[M->numCols] + maxShift + 10;
  }
  /* Move everything right by `maxShift'. Of course, there's no need 
     to move the first block (i.e., columns left of the first column
     to be changed.)
  */
  sIndex = M->cIndex[pairs[0].x+1];
  dIndex = sIndex + maxShift;
  moveSize = M->cIndex[M->numCols] - sIndex;
  memmove(&M->cEntry[dIndex], &M->cEntry[sIndex], moveSize*sizeof(s32));

  /* Now, do the column additions, moving blocks back to the left as we go: */
  indexShift = 0;
  for (i=0; i<numCols; i++) {
    c0 = pairs[i].x; c1 = pairs[i].y;
    /* At the top of this loop, columns [0, pairs[i].x] are in proper
       position and their cIndex values are correct. 
       The remaining columns, [pairs[i].x+1, M->numCols),
       are all shifted right by maxShift;
    */
    /* Do the easy part first: dense entries: */
    for (j=0; j<M->numDenseBlocks; j++)
      M->denseBlocks[j][c0] ^= M->denseBlocks[j][c1];

    /* Right now, the columns in [0, pairs[i].x] are in proper position,
       starting at M->cEntry[0], and they all have correct cIndex values.
       The columns in [pairs[i].x + 1, ..., M->numCols) are all shifted
       in memory by `maxShift', and they're corresponding cIndex values
       are off by that amount as well.
    */

    /* Compute the new entries for column c0: */
    numNewEntries=0;
    for (j=M->cIndex[c0]; j<M->cIndex[c0+1]+indexShift/*aha!*/; j++)
      newEntries[numNewEntries++] = M->cEntry[j];
    numOldEntries = numNewEntries;
    if (c1 <=c0) {
      for (j=M->cIndex[c1]; j<M->cIndex[c1+1]; j++)
        newEntries[numNewEntries++] = M->cEntry[j];
    } else {
      for (j=M->cIndex[c1]; j<M->cIndex[c1+1]; j++)
        newEntries[numNewEntries++] = M->cEntry[j+maxShift];
    }
    numNewEntries = removeS32Pairs(newEntries, numNewEntries);
    /* Place the new entries: */
    for (j=0; j<numNewEntries; j++)
      M->cEntry[M->cIndex[c0]+j] = newEntries[j];

    /* Move the next block of data and adjust it's indicies. */
    indexShift += (numNewEntries - numOldEntries);
    sIndex = M->cIndex[c0+1] + maxShift;
    dIndex = M->cIndex[c0]+numNewEntries;
    if (i<(numCols-1)) {
      moveSize = M->cIndex[pairs[i+1].x+1] - sIndex + maxShift;
      memmove(&M->cEntry[dIndex], &M->cEntry[sIndex], moveSize*sizeof(s32));
      for (j=c0+1; j<=pairs[i+1].x; j++) 
        M->cIndex[j] += indexShift;
    } else {
      moveSize = M->cIndex[M->numCols] - sIndex + maxShift;
      memmove(&M->cEntry[dIndex], &M->cEntry[sIndex], moveSize*sizeof(s32));
      for (j=c0+1; j<=M->numCols; j++) 
        M->cIndex[j] += indexShift;
    }
  }  
  /* Adjust the colmap: */
  cm_addCols_par(C, pairs, numCols);

  free(pairs);
  return 0;
}

/******************************************************/
int removeSingletons(nfs_sparse_mat_t *M, llist_t *C)
/******************************************************/
/* A `singleton' is the column corresponding to a row */
/* which has exactly 1 nonzero entry. Such columns    */
/* are of no use anyway.                              */
/******************************************************/
{ s32    i, j, singletons[1024], *rowMember;
  s32    totalSingletons, r;
  s32    numSingletons;

  rowMember = (s32 *)malloc(M->numRows*sizeof(s32));

  /* In this way, we'll never be able to remove column zero, but oh well. */
  memset(rowMember, 0x00, M->numRows*sizeof(s32));
  for (j=0; j<M->numCols; j++) {
    for (i=M->cIndex[j]; i<M->cIndex[j+1]; i++) {
      r = M->cEntry[i];
      if ((rowMember[r]==0)&&(r>0)) 
        rowMember[r] = j;
      else
        rowMember[r]=-1;
    }
  }

  totalSingletons=0;
  numSingletons=0;
  for (i=0; (i<(M->numRows-64)) && (numSingletons<1024); i++) {
    if (rowMember[i]>0) {
      singletons[numSingletons++] = rowMember[i];
      for (j=0; j<(numSingletons-1); j++) {
        if (singletons[j]==singletons[numSingletons-1]) {
          numSingletons--; break;
        }
      }        
    }
  }
  if (numSingletons) 
    removeCols(M, C, singletons, numSingletons);
  free(rowMember);
  return 0;
}

/******************************************************/
int removeSingletons2(nfs_sparse_mat_t *M, llist_t *C)
/******************************************************/
/* A `singleton' is the column corresponding to a row */
/* which has exactly 1 nonzero entry. Such columns    */
/* are of no use anyway.                              */
/******************************************************/
/* For the life of me, I can't figure out why, but this */
/* function seems to break sometimes. It's supposed to  */
/* be a cleaner version of removeSingletons(), and it   */
/* does seem to work sometimes. But other times, it is  */
/* just horribly wrong, and deletes non-singleton cols. */
{ s32    i, j, *singletons, *rowWt;
  s32    numSingletons, maxSingletons, r;

  rowWt = (s32 *)calloc(M->numRows, sizeof(s32));

  for (j=0; j<M->numCols; j++) {
    for (i=M->cIndex[j]; i<M->cIndex[j+1]; i++) {
      r = M->cEntry[i];
      rowWt[r] += 1;
    }
  }
  maxSingletons=0;
  for (i=0; i<(M->numRows-64); i++)
    if (rowWt[i] <= 1)
      maxSingletons++;
  if (maxSingletons <= 0) {
    free(rowWt);
    return 0;
  }

  singletons = (s32 *)malloc(maxSingletons*sizeof(s32));
  numSingletons=0; /* There could wind up being strictly less than maxSingletons, if
                      some columns have more than one of these singly occurring primes. */
  for (j=0; j<M->numCols; j++) {
    for (i=M->cIndex[j]; i<M->cIndex[j+1]; i++) {
      r = M->cEntry[i];
      if ((rowWt[r] == 1)&& (numSingletons < maxSingletons)) {
        singletons[numSingletons++] = j;
//        i = M->cIndex[j+1];
        break;
      }
    }
  }

  if (numSingletons) 
    removeCols(M, C, singletons, numSingletons);
  free(rowWt);
  free(singletons);
  return 0;
}

/******************************************************/
int removeEmptyRows(nfs_sparse_mat_t *M)
/******************************************************/
{ s32    i, j, *rowWeight, emptyRows[1024];
  s32    totalEmpty, c;
  int     numEmpty;

  rowWeight = (s32 *)malloc(M->numRows*sizeof(s32));

  memset(rowWeight, 0x00, M->numRows*sizeof(s32));
  for (j=M->cIndex[M->numCols]-1; j>=0; j--) {
    c = M->cEntry[j];
    if ((c>=0) && (c<M->numRows))
      rowWeight[c] = 1; /* Not weight anymore, but whatever. */
  }
  numEmpty=0; totalEmpty=0;

  for (i=0; (i<(M->numRows-64)) && (numEmpty<1024); i++) {
    if ((rowWeight[i]==0)&&(!(ISDENSEROW(i)))) {
      totalEmpty++;
      emptyRows[numEmpty++] = i;
    }  
  }
  if (numEmpty)
    removeRows(M, emptyRows, numEmpty);

  free(rowWeight);
  return 0;
}

/****************************************************/
int cmpCW(const void *a, const void *b)
{ s32 *A=(s32 *)a, *B=(s32 *)b;
                                                                                
  if (A[1]<B[1]) return 1;
  if (A[1]>B[1]) return -1;
  return 0;
}

/****************************************************/
int cmpCW2(const void *a, const void *b)
{ s32 *A=(s32 *)a, *B=(s32 *)b;
                                                                                
  if (A[1]<B[1]) return -1;
  if (A[1]>B[1]) return 1;
  return 0;
}

/***************************************************/
int removeHeavyColumns(nfs_sparse_mat_t *M, llist_t *C, s32 numC)
/***************************************************/
{ s32 *cwt, cols[2048];
  s32  c, numRemove;
                                                                                
  cwt = (s32 *)malloc(2*M->numCols*sizeof(s32));
  for (c=0; c<M->numCols; c++) {
    cwt[2*c]=c;
    cwt[2*c+1] = M->cIndex[c+1]-M->cIndex[c];
  }
  qsort(cwt, M->numCols, 2*sizeof(s32), cmpCW);
  for (numRemove=0; (numRemove<numC)&&(numRemove<2048); numRemove++)
    cols[numRemove] = cwt[2*numRemove];
  removeCols(M, C, cols, numRemove);
  free(cwt);
  return 0;
}

/***************************************************/
int removeHeavyColumnsByRows(nfs_sparse_mat_t *M, llist_t *C, s32 minC, s32 maxC, s32 maxMWt)
/***************************************************/
/* Same as above, but a little smarter: Remove     */
/* some heavy rows which have some entries from    */
/* low weight rows.                                */
/* We remove at least minC, at most maxC, stopping */
/* as soon as the total weight of the matrix is    */
/* maxMWt or we hit maxC.                          */
/* If maxMWt==0, it will be effectively ignored.   */
/***************************************************/
{ s32    *cwt, cols[2048], rmWt=0;
  s32    c, numRemove, x;
  s32    i, j, *rowWeight;
  s32    rowWtCt[512], ct, mWt = M->cIndex[M->numCols];
  int     maxWt;

  if (!(cwt = (s32 *)malloc(2*M->numCols*sizeof(s32)))) {
    printf("removeHeavyColumnsByRows() memory allocation error!\n");
    printf("( %" PRIu32 " bytes requested).\n", (u32)(2*sizeof(s32)*M->numCols) );
    return -1;
  }
  for (c=0; c<M->numCols; c++) {
    cwt[2*c]=c;
    cwt[2*c+1] = M->cIndex[c+1]-M->cIndex[c];
  }
  qsort(cwt, M->numCols, 2*sizeof(s32), cmpCW);

  if (!(rowWeight = (s32 *)malloc(M->numRows*sizeof(s32)))) {
    free(cwt);
    printf("removeHeavyColumnsByRows() memory allocation error!\n");
    printf("( %" PRIu32 " bytes requested).\n", (u32)(sizeof(s32)*M->numRows) );
    return -1;
  }

  memset(rowWeight, 0x00, M->numRows*sizeof(s32));
  for (j=M->cIndex[M->numCols]-1; j>=0; j--) {
    x = M->cEntry[j];
    if ((x>=0) && (x<M->numRows)) {
      rowWeight[x] += 1; 
    }
  }
  for (i=0; i<512; i++)
    rowWtCt[i]=0;

  for (i=0; i<(M->numRows-64); i++) {
    if (rowWeight[i] < 512)
      rowWtCt[rowWeight[i]] += 1; 
  }
  ct=0; i=0;
  while ((i<512) && (ct < maxC))
    ct += rowWtCt[i++];
  maxWt = i;

  /* Now scan the cwt list, which is sorted on column weight,
     looking for some with entries from rows with weight <= maxWt.
  */
  numRemove=0;
  for (i=0; (i<M->numCols) && (numRemove < maxC) && (numRemove < 2048); i++) {
    c = cwt[2*i];
    for (j=M->cIndex[c]; j<M->cIndex[c+1]; j++) {
      if (rowWeight[M->cEntry[j]] <= maxWt) {
        cols[numRemove++] = c;
        rmWt += M->cIndex[c+1] - M->cIndex[c];
        j = M->cIndex[c+1];
      }
    }
    if (((mWt - rmWt) < maxMWt) && (numRemove > minC)) break;
  }
  free(rowWeight);
  free(cwt);
  removeCols(M, C, cols, MIN(numRemove, maxC));

  return 0;
}

/***************************************************/
int cmpList1(const void *a, const void *b)
/***************************************************/
{ s32 *A=(s32 *)a, *B=(s32 *)b;

  if (A[1] < B[1]) return -1;
  if (A[1] > B[1]) return 1;
  return 0;
}

/***************************************************/
int addLightColumnsByRows(nfs_sparse_mat_t *M, llist_t *C, s32 numC)
/***************************************************/
/* Look for low weight columns which have entries  */
/* from a low weight row. Add that column to the   */
/* others.                                         */
/***************************************************/
{ s32    *cwt, *tmpPtr;
  s32    c, ct;
  s32    i, j, *rowWeight;
  s32    rowWtCt[512], numAdds, size;
  int     maxWt, minWt;
  s32    *list, m, *dest, *src;

  if (numC <= 0) return 0;
  cwt = (s32 *)malloc(2*M->numCols*sizeof(s32));
  for (c=0; c<M->numCols; c++) {
    cwt[2*c]=c;
    cwt[2*c+1] = M->cIndex[c+1]-M->cIndex[c];
  }
  qsort(cwt, M->numCols, 2*sizeof(s32), cmpCW2);

  rowWeight = (s32 *)calloc(M->numRows,sizeof(s32));
  for (j=M->cIndex[M->numCols]-1; j>=0; j--)
    rowWeight[M->cEntry[j]] += 1; 

  for (i=0; i<512; i++)
    rowWtCt[i]=0;

  for (i=0; i<(M->numRows-64); i++) {
    if (rowWeight[i] < 512)
      rowWtCt[rowWeight[i]] += 1; 
  }


  ct=0; i=0;
  while ((i<512) && (ct < numC))
    ct += rowWtCt[i++];
  maxWt = i;


  /* Now scan the cwt list, which is sorted on column weight,
     looking for some with entries from rows with weight <= maxWt.
  */
  numAdds=0;
  size = 4*numC;
  list = (s32 *)malloc(3*size*sizeof(s32));
  for (i=0; (i<M->numCols); i++) {
    c = cwt[2*i];
    for (j=M->cIndex[c]; j<M->cIndex[c+1]; j++) {
      if (rowWeight[M->cEntry[j]] <= maxWt) {
        /* This is a column we will use. What we need to do
           now is keep track of all columns with an entry
           in this row. When we're done, we'll take the lightest
           one (i.e., the first one we found, since they're sorted),
           and add it to the others.
        */
        list[3*numAdds] = c; 
        list[3*numAdds+1] = M->cEntry[j];
        list[3*numAdds+2] = M->cIndex[c+1] - M->cIndex[c];
        numAdds++;
        if (numAdds >= size) {
          size += numC;
          tmpPtr = realloc(list, (size+1)*3*sizeof(s32));
          if (tmpPtr != NULL) list=tmpPtr;
          else {
            printf("addLightColumnsbyRow() memory reallocation error!\n");
          }
        }
        break;
      }
    }
  }
  if (numAdds <= 0) {
    free(cwt);
    free(rowWeight);
    free(list);
    return 0;
  }

  /* Now we need to pair them off: */
  /* First, sort them by entry so we can see which should
     be added to which.
  */
  qsort(list, numAdds, 3*sizeof(s32), cmpList1);
  numAdds = MIN(numAdds, numC);

  dest = (s32 *)malloc(numAdds*sizeof(s32));
  src = (s32 *)malloc(numAdds*sizeof(s32));
  i=0;
  m=0;
  while (i<numAdds) {
    j=i;
    minWt=list[3*i+2]; c=i;
    while ((j<numAdds) && (list[3*j+1]==list[3*i+1])) {
      if (list[3*j+2] < minWt) {
        minWt = list[3*j+2];
        c = j;
      }
      j++;
    }
    
    /* Now, column 'c' is the one which we should add to the others,
       and columns list[3*i,3*(i+1),..., 3*(j-1)] are the ones with the same entry.
    */
    while (i<j) {
      if (i != c) {
        dest[m] = list[3*i];
        src[m] = list[3*c];
        m++;
      }
      i++;
    }
  }
  free(list);
  free(rowWeight);
  free(cwt);
  addCols_par(M, C, dest, src, m);

  free(dest); free(src);

  return 0;
}





typedef struct {
  s32 c, w;
} cw_t;
typedef struct {
  s32 c0, c1, c2, w;
} cw_t3;
/*******************************************/
int cmp_cw_t(const void *a, const void *b)
/*******************************************/
{ cw_t *A=(cw_t *)a, *B=(cw_t *)b;

  if (A->w < B->w) return 1;
  if (A->w > B->w) return -1;
  return 0;
}
/*******************************************/
int cmp_cw_t3(const void *a, const void *b)
/*******************************************/
{ cw_t3 *A=(cw_t3 *)a, *B=(cw_t3 *)b;

  if (A->w < B->w) return -1;
  if (A->w > B->w) return 1;
  return 0;
}
/******************************************************/
int removeDoubles(nfs_sparse_mat_t *M, llist_t *C, s32 numRemove)
/******************************************************/
/* A `singleton' is the column corresponding to a row */
/* which has exactly 1 nonzero entry. A `double' is   */
/* the obvious generalization.                        */
/******************************************************/
{ s32    i, j, D[1024], *rowMember, *rm2, k;
  s32    totalD, r, num2=0, c;
  int     numD;
  cw_t    *CW;

  rowMember = (s32 *)malloc(M->numRows*sizeof(s32));
  rm2 = (s32 *)malloc(M->numRows*sizeof(s32));
  memset(rowMember, 0x00, M->numRows*sizeof(s32));
  memset(rm2, 0x00, M->numRows*sizeof(s32));
  for (j=0; j<M->numCols; j++) {
    for (i=M->cIndex[j]; i<M->cIndex[j+1]; i++) {
      r = M->cEntry[i];
      if (r>0) {
        if (rowMember[r]==0) {
          if (rm2[r]==0) {
            rm2[r]=j;
          } else {
            rowMember[r] = j;
            num2++;
          }
        } else {
          if (rowMember[r] > 0)
            num2--;
          rowMember[r]=-1;
        }
      }
    }
  }
  /* Run through the double's and grab the column weights, so
     we can sort on that and remove the heaviest ones.
  */
  if (num2 <= 0) {
    free(rowMember); free(rm2);
    return 0;
  }
  k=0;
  CW = (cw_t *)malloc(num2*sizeof(cw_t));
  for (j=0; (j<M->numRows-64) && (k<num2); j++) {
    c = rowMember[j];
    if (c>0) {
      CW[k].c = c;
      CW[k++].w = M->cEntry[c+1]-M->cEntry[c];
    }
  }
  qsort(CW, k, sizeof(cw_t), cmp_cw_t);


  totalD=0;
  numD=0;
  for (i=0; (i<k) && (numD<1024) && (numD < numRemove); i++) {
    c = CW[i].c;
    if (c>0) {
      D[numD++] = c;
      for (j=0; j<(numD-1); j++) {
        if (D[j]==D[numD-1]) {
          numD--; break;
        }
      }        
    }
  }
  if (numD) 
    removeCols(M, C, D, numD);

  free(CW);
  free(rowMember);
  free(rm2);
  return 0;
}

/******************************************************/
int combineDoubles(nfs_sparse_mat_t *M, llist_t *C)
/******************************************************/
/* A `singleton' is the column corresponding to a row */
/* which has exactly 1 nonzero entry. We will combine */
/* (low weight) columns corresponding to `doubles'.   */
/* The result is that we lower the dimensions of the  */
/* matrix 1x1 without changing the weight (at least,  */
/* that's what will happen the next time somebody     */
/* calls removeSingletons() and removeEmptyRows(). )  */
/******************************************************/
{ s32    i, j, *rowMember, *rm2, k;
  s32    r, num2=0, c0, c1;
  s32    *src, *dest, size, *tmpPtr;

  rowMember = (s32 *)malloc(M->numRows*sizeof(s32));
  rm2 = (s32 *)malloc(M->numRows*sizeof(s32));
  memset(rowMember, 0x00, M->numRows*sizeof(s32));
  memset(rm2, 0x00, M->numRows*sizeof(s32));
  for (j=0; j<M->numCols; j++) {
    for (i=M->cIndex[j]; i<M->cIndex[j+1]; i++) {
      r = M->cEntry[i];
      if (r>0) {
        if (rowMember[r]==0) {
          if (rm2[r]==0) {
            rm2[r]=j;
          } else {
            rowMember[r] = j;
            num2++;
          }
        } else {
          if (rowMember[r] > 0)
            num2--;
          rowMember[r]=-1;
        }
      }
    }
  }

  size=4096;
  src = (s32 *)malloc(size*sizeof(s32));
  dest = (s32 *)malloc(size*sizeof(s32));
  if ((src==NULL)||(dest==NULL)) {
    printf("combineDoubles() memory allocation error!\n");
    exit(-1);
  }
  k=0;
  for (j=0; (j<M->numRows-64); j++) {
    c0 = rowMember[j];
    c1 = rm2[j];
    if (c0>0) {
      dest[k]=c0;
      src[k]=c1;
      if (++k >= (size-1)) {
        size += 512;
        tmpPtr = realloc(src, (size+1)*sizeof(s32));
        if (tmpPtr==NULL) {
          printf("Severe memory re-allocation error!\n");
          exit(-1);
        } else src=tmpPtr;
        tmpPtr = realloc(dest, (size+1)*sizeof(s32));
        if (tmpPtr==NULL) {
          printf("Severe memory re-allocation error!\n");
          exit(-1);
        } else dest=tmpPtr;
      }
    }
  }
  free(rowMember);
  free(rm2);

  if (k>0) {
    addCols_par(M, C, dest, src, k);
  }
  free(src); free(dest);
  return 0;
}

/****************************************************/
s32 matrixWeight(nfs_sparse_mat_t *M)
/****************************************************/
{ int  cwt[256];
  s32 i, j, k, w;
  u64  e;
                                                                                                            
  for (i=0; i<256; i++) {
    w=0;
    for (j=0;j<8;j++)
      if (i&BIT(j)) w++;
    cwt[i]=w;
  }
  w=0;
  for (k=0; k<M->numDenseBlocks; k++) {
    for (i=0; i<M->numCols; i++) {
      e = M->denseBlocks[k][i];
      w += cwt[e&0x00000000000000FFULL] + 
           cwt[(e&0x000000000000FF00ULL)>>8] +
           cwt[(e&0x0000000000FF0000ULL)>>16] + 
           cwt[(e&0x00000000FF000000ULL)>>24] +
           cwt[(e&0x000000FF00000000ULL)>>32] +
           cwt[(e&0x0000FF0000000000ULL)>>40] +
           cwt[(e&0x00FF000000000000ULL)>>48] +
           cwt[(e&0xFF00000000000000ULL)>>56];
    }
  }
  return w + M->cIndex[M->numCols];
}

/***************************************************/
int pruneMatrix(nfs_sparse_mat_t *M, s32 minExtraCols, double wtFactor,
                llist_t *C)
/***************************************************/
/* minExtraCols is the minimum number of extra     */
/* columns to leave after pruning. That is, if M'  */
/* is the pruned matrix we return, then            */
/* (cols M') - (rows M') will be >= minExtraCols.  */
/***************************************************/
{ char str[256];
  s32 extraCols, lastR, lastC, origWt;
  double wt;

  wt = MIN(0.95, MAX(0.05, wtFactor));
  sprintf(str, "Initial matrix is %" PRId32 " x %" PRId32 " with sparse part having weight %" PRId32 ".", 
          M->numRows, M->numCols, M->cIndex[M->numCols]);
  msgLog(NULL, "Pruning matrix with wt=%1.3lf", wt);
  origWt = M->cIndex[M->numCols];
  printf("%s\n",str); msgLog("", str);
  sprintf(str, "(total weight is %" PRId32 ")", matrixWeight(M));
  printf("%s\n",str); msgLog("", str);

  setDenseRows(M);
  mat_verify(M);

  removeSingletons(M, C);
  removeEmptyRows(M);
  do {
    mat_verify(M);
    lastR = M->numRows; lastC = M->numCols;
    extraCols = M->numCols - M->numRows - minExtraCols;
    combineDoubles(M, C);
    printTmp("Current matrix is %ld x %ld with weight %ld.          ", 
              M->numRows, M->numCols, M->cIndex[M->numCols]);
    removeSingletons(M, C);
    removeEmptyRows(M);
    extraCols = M->numCols - M->numRows - minExtraCols;
    printTmp("Current matrix is %ld x %ld with weight %ld..          ", 
              M->numRows, M->numCols, M->cIndex[M->numCols]);
    if (extraCols > 0) {
      if (extraCols < 200) {
        removeHeavyColumns(M, C, extraCols);
      } else {
        removeHeavyColumnsByRows(M, C, (s32)(0.3*wt*extraCols), (s32)(0.6*wt*extraCols), 0);
        removeHeavyColumns(M, C, (s32)(0.3*(1.0-wt)*extraCols));
        addLightColumnsByRows(M, C, (s32)(0.3*(1.0-wt)*extraCols));
      }
      printTmp("Current matrix is %ld x %ld with weight %ld...          ", 
                M->numRows, M->numCols, M->cIndex[M->numCols]);
    }
    removeSingletons(M, C);
    removeEmptyRows(M);
    extraCols = M->numCols - M->numRows - minExtraCols;

    printTmp("Current matrix is %ld x %ld with weight %ld....          ", 
              M->numRows, M->numCols, M->cIndex[M->numCols]);
  } while (extraCols > 0);

  mat_verify(M);
  sprintf(str, "Matrix pruned to %" PRId32 " x %" PRId32 " with weight %" PRId32 ".",
        M->numRows, M->numCols, M->cIndex[M->numCols]);
  printf("%s\n", str); msgLog("", str);
  mat_verify(M);
  return 0;
}


/***************************************************/
int getDependencies(nfs_sparse_mat_t *M, llist_t *C, s32 *deps)
/***************************************************/
{ double blstart, blstop, difficulty;
  int    res;
  s32   i, j, origC = M->numCols;
  u64   *tmpDeps;

  difficulty = (M->numCols/64.0)*(M->cIndex[M->numCols] + M->numCols*M->numDenseBlocks);
  difficulty /= 1000000.0;
  printf("Matrix difficulty is about %1.2lf\n", difficulty);
  printf("Doing block Lanczos...\n");
  blstart = sTime();
  tmpDeps = (u64 *)malloc(M->numCols*sizeof(u64));
  res = blockLanczos64(tmpDeps, MultB64, MultB_T64, (void *)M, M->numCols);
  blstop = sTime();
  printf("Returned %d. Block Lanczos took %1.2lf seconds.\n", res, blstop-blstart);
  msgLog("", "BLanczosTime: %1.1lf", blstop-blstart);
  if (res < 0) return res;

  /******************************************************/
  /* Now, we need to translate the dependencies back to */
  /* deps for the original matrix.                      */
  /******************************************************/
  memset(deps, 0x00, origC*sizeof(s32));

  for (i=0; i<M->numCols; i++) {
    for (j=C->index[i]; j<C->index[i+1]; j++)
      deps[C->data[j]] ^= (s32)(tmpDeps[i]&0xFFFFFFFF);
  }

  free(tmpDeps);  
  return 0;
}


/***************************************************/
int writeSparseMat(char *fname, nfs_sparse_mat_t *M)
/***************************************************/
{ long  i;
  FILE *fp;

  if (!(fp = fopen(fname, "wb"))) {
    fprintf(stderr, "Could not open %s for write!\n", fname);
    return -1;
  }
  fwrite(&M->numRows, sizeof(s32), 1, fp);
  fwrite(&M->numCols, sizeof(s32), 1, fp);
  fwrite(&M->maxDataSize, sizeof(s32), 1, fp);
  fwrite(&M->numDenseBlocks, sizeof(s32), 1, fp);
  fwrite(M->denseBlockIndex, sizeof(s32), M->numDenseBlocks, fp);
  fwrite(M->cIndex, sizeof(s32), M->numCols+1, fp);
  fwrite(M->cEntry, sizeof(s32), M->cIndex[M->numCols], fp);
  for (i=0; i<M->numDenseBlocks; i++)
    fwrite(M->denseBlocks[i], sizeof(u64), M->numCols, fp);
  fclose(fp);
  return 0;
}

/***************************************************/
int readSparseMat(nfs_sparse_mat_t *M, char *fname)
/***************************************************/
{ long  i;
  FILE *fp;

  if (!(fp = fopen(fname, "rb"))) {
    fprintf(stderr, "Could not open %s for read!\n", fname);
    return -1;
  }
  fread(&M->numRows, sizeof(s32), 1, fp);
  fread(&M->numCols, sizeof(s32), 1, fp);
  fread(&M->maxDataSize, sizeof(s32), 1, fp);
  fread(&M->numDenseBlocks, sizeof(s32), 1, fp);
  fread(M->denseBlockIndex, sizeof(s32), M->numDenseBlocks, fp);
  M->cIndex = (s32 *)lxmalloc((M->numCols+1)*sizeof(s32), 1);
  fread(M->cIndex, sizeof(s32), M->numCols+1, fp);
  M->cEntry = (s32 *)lxmalloc((M->cIndex[M->numCols]+1)*sizeof(s32),1);
  fread(M->cEntry, sizeof(s32), M->cIndex[M->numCols], fp);
  for (i=0; i<M->numDenseBlocks; i++) {
    M->denseBlocks[i] = (u64 *)lxmalloc(M->numCols*sizeof(u64),1);
    fread(M->denseBlocks[i], sizeof(u64), M->numCols, fp);
  }
  fclose(fp);
  return 0;
}

