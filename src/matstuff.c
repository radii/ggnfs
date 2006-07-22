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
  M->maxDataSize = M->cIndex[M->numCols]+1;
  M->cEntry = (s32 *)lxmalloc((M->cIndex[M->numCols]+1)*sizeof(s32),1);
  fread(M->cEntry, sizeof(s32), M->cIndex[M->numCols], fp);
  for (i=0; i<M->numDenseBlocks; i++) {
    M->denseBlocks[i] = (u64 *)lxmalloc(M->numCols*sizeof(u64),1);
    fread(M->denseBlocks[i], sizeof(u64), M->numCols, fp);
  }
  fclose(fp);
  return 0;
}


/*********************************************************/
static int cmp2L1(const void *a, const void *b)
/*********************************************************/
{ s32 *A=(s32 *)a, *B=(s32 *)b;

  if (A[0] < B[0]) return -1;
  if (A[1] > B[1]) return 1;
  return 0;
}

/*********************************************************/
static int colsAreEqual(nfs_sparse_mat_t *M, s32 c0, s32 c1)
/*********************************************************/
{ int i, j, found;
  s32 w0, w1, e;

  w0 = M->cIndex[c0+1]-M->cIndex[c0];
  w1 = M->cIndex[c1+1]-M->cIndex[c1];
  if (w0 != w1) return 0;
  for (i=M->cIndex[c0]; i<M->cIndex[c0+1]; i++) {
    e = M->cEntry[i];
    for (j=M->cIndex[c1], found=0; j<M->cIndex[c1+1]; j++)
      if (M->cEntry[j]==e) found=1;
    if (!(found)) return 0;
  }
  for (i=2; i<M->numDenseBlocks; i++) {
    if (M->denseBlocks[i][c0] != M->denseBlocks[i][c1])
      return 0;
  }
  return 1; 
}

/*********************************************************/
int checkMat(nfs_sparse_mat_t *M, s32 *delCols, s32 *numDel)
/*********************************************************/
/* Do a sanity check on the matrix: Make sure there are  */
/* no all zero columns and such.                         */
/*********************************************************/
{ s32 c, i, w;
  int nz, warn=0;
  s32 *colHash, h, hi;

  *numDel=0;
  colHash = (s32 *)lxmalloc(2*M->numCols*sizeof(s32),1);
  for (c=0; c<M->numCols; c++) {
    w = M->cIndex[c+1] - M->cIndex[c];
    if (w == 0) {
      /* Skip QCB & sign entries. */
      for (i=2, nz=0; i<M->numDenseBlocks; i++) {
        if (M->denseBlocks[i][c])
          nz=1;
      }
      if (nz==0) {
        printf("Warning: column %" PRId32 " is all zero!\n", c);
        if (*numDel < 2048) 
          delCols[(*numDel)++]=c;
        warn=1;
      }
    }
    h = 0x00000000;
    for (i=M->cIndex[c]; i<M->cIndex[c+1]; i++) {
      hi = NFS_HASH(0, M->cEntry[i], 0x8FFFFFFF);
      h ^= hi;
    }
    for (i=NUM_QCB_BLOCKS; i<M->numDenseBlocks; i++) 
      h ^= (M->denseBlocks[i][c] & 0xFFFFFFFF) ^ 
           ((M->denseBlocks[i][c] & 0xFFFFFFFF00000000ULL) >> 32);
    colHash[2*c] = h;
    colHash[2*c+1]=c;
  }
  qsort(colHash, M->numCols, 2*sizeof(s32), cmp2L1);
  for (i=0; i<(M->numCols-1); i++) {
    if (colHash[2*i]==colHash[2*i+2]) {
      if (colsAreEqual(M, colHash[2*i+1], colHash[2*i+3])) {
        printf("Bad matrix: column %" PRId32 " = column %" PRId32 "!\n", 
               colHash[2*i+1], colHash[2*i+3]);
        if (*numDel < 2048) 
          delCols[(*numDel)++]=colHash[2*i+1];
        
        warn=2;
      }
    }
  }

  if (warn) {
    printf("checkMat() did not like something about the matrix:\n");
    printf("This is probably a sign that something has gone horribly wrong\n");
    printf("in the matrix construction (matbuild).\n");
    if (*numDel < 2048) {
      printf("However, the number of bad columns is only %" PRId32 ",\n", *numDel);
      printf("so we will (probably) try delete them and attempt to continue.\n");
    }
    else warn = 2048;
  }
  free(colHash);
  return warn;
}

