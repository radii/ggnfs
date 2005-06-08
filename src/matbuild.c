/**************************************************************/
/* matbuild.c                                                 */
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
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>


#ifndef _MSC_VER 
#include <sys/time.h>
#endif
#include "ggnfs.h"




#define MAX_LPMEM_ALLOC 256000000
#define MAX_SPAIRS_ALLOC 12000000

/* These are in s32's. */
#define IO_BUFFER_SIZE  2*1024*1024

#define MAX_PBUF_RAM 32000000
#define DEFAULT_QCB_SIZE 62
#define DEFAULT_SEED 1
#define DEFAULT_NUM_FILES 4
#define TMP_FILE "tmpdata.000"
#define DEFAULT_DEPNAME "deps"
#define DEFAULT_COLNAME "cols"
#define DEFAULT_LPI_NAME "lpindex"
#define DEFAULT_PRELPREFIX "rels.bin"
#define DEFAULT_COLNAME "cols"
#define DEFAULT_WT_FACTOR 0.7

#define USAGE " -fb <fname> -prel <fname> -newrel <fname> [-qs <qcb size>] [-v]\n"\
"-fb <fname>           : Factor base.\n"\
"-prel <file prefix>   : File name prefix for input of processed relations.\n"\
"-minff <int>          : Minimum number of FF's (prevent R-S wt. reduction and\n"\
"                        writing of the column files if there are fewer than this).\n"\
"-maxrelsinff <int>    : Max relation-set weight.\n"

#define START_MSG \
"\n"\
" __________________________________________________________ \n"\
"|        This is the matbuild program for GGNFS.           |\n"\
"| Version: %-15s                                 |\n"\
"| This program is copyright 2004, Chris Monico, and subject|\n"\
"| to the terms of the GNU General Public License version 2.|\n"\
"|__________________________________________________________|\n"

/***********************************************************/
/* This structure is for columns that are being processed. */
/* Before processing, we know only the relations that      */
/* comprise each column. After processing, we know all the */
/* nonzero rows of the matrix (i.e., the primes appearing  */
/* with odd exponent and the QCB (sign bit is in QCB)).    */
/* During processing, we will know some of each. That is,  */
/* we might have already replaced some relations with their*/
/* factorizations, but not yet all of them. These data are */
/* what will be stored in the file cols.np.                */
/***********************************************************/
typedef struct {
  s32 numRels;
  s32 numPrimes;
  s32 Rels[MAX_RELS_IN_FF];
  s32 QCB[2];
  s32 rows[MAX_ROWS_IN_COL];
} column_t;



#define CC_OFF  0
#define CC_ON   1
#define CC_AUTO 2

/***** Globals *****/
int  verbose=0, discFact=1, cycleCount=CC_AUTO;
s32  initialFF=0, initialRelations=0, finalFF=0;
s32  totalLargePrimes=0;
s32  minFF;
long relsNumLP[8]={0,0,0,0,0,0,0,0};
s32  delCols[2048], numDel=0;

int cmp2S32s(const void *a, const void *b);


/******************************************************/
void readColSF(column_t *C, FILE *fp)
/******************************************************/
{ 
  fread(&C->numPrimes, sizeof(s32), 1, fp);
  fread(&C->QCB, sizeof(s32), 2, fp);
  if (C->numPrimes < MAX_ROWS_IN_COL) 
    fread(&C->rows, sizeof(s32), C->numPrimes, fp);
  else {
    printf("readColSF() : MAX_ROWS_IN_COL exceeded! Cannot continue!\n");
    exit(-1);
  }
}

/*********************************************************/
int cmp2L1(const void *a, const void *b)
/*********************************************************/
{ s32 *A=(s32 *)a, *B=(s32 *)b;

  if (A[0] < B[0]) return -1;
  if (A[1] > B[1]) return 1;
  return 0;
}

/*********************************************************/
int colsAreEqual(nfs_sparse_mat_t *M, s32 c0, s32 c1)
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
int checkMat(nfs_sparse_mat_t *M)
/*********************************************************/
/* Do a sanity check on the matrix: Make sure there are  */
/* no all zero columns and such.                         */
/*********************************************************/
{ s32 c, i, w;
  int nz, warn=0;
  s32 *colHash, h, hi;

  if (!(colHash = (s32 *)lxmalloc(2*M->numCols*sizeof(s32),0))) {
    printf("checkMat() Memory allocation error!\n");
    printf("** Cannot do sanity check on the matrix!\n");
    return -1;
  }
  for (c=0; c<M->numCols; c++) {
    w = M->cIndex[c+1] - M->cIndex[c];
    if (w == 0) {
      /* Skip QCB & sign entries. */
      for (i=2, nz=0; i<M->numDenseBlocks; i++) {
        if (M->denseBlocks[i][c])
          nz=1;
      }
      if (nz==0) {
        printf("Warning: column %ld is all zero!\n", c);
        if (numDel < 2048) 
          delCols[numDel++]=c;
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
        printf("Bad matrix: column %ld = column %ld!\n", 
               colHash[2*i+1], colHash[2*i+3]);
        if (numDel < 2048) 
          delCols[colHash[2*i+3]]=2*i+1;
        
        warn=2;
      }
    }
  }

  if (warn) {
    printf("checkMat() did not like something about the matrix:\n");
    printf("This is probably a sign that something has gone horribly wrong\n");
    printf("in the matrix construction (procrels).\n");
    if (numDel < 2048) {
      printf("However, the number of bad columns is only %ld,\n", numDel);
      printf("so we will delete them and attempt to continue.\n");
    }
  }
  free(colHash);
  return warn;
}

/*********************************************************/
s32 loadMat(nfs_sparse_mat_t *M, char *colName)
/*********************************************************/
{ FILE  *fp;
  s32   i, j, k, fileSize, index, nR, r, *rwt, blockWt;
  struct stat fileInfo;
  column_t     C;
  int    bSize;
  double dW;


  if (stat(colName, &fileInfo)) {
    printf("loadMat() Could not stat %s!\n", colName);
    return -1;
  }
  fileSize = fileInfo.st_size;

  if (!(fp = fopen(colName, "rb"))) {
    fprintf(stderr, "loadMat() Unable to open %s for read!\n", colName);
    return -1;
  }
  fread(&M->numCols, sizeof(s32), 1, fp);
  /* Scan the file once to find out how many dense rows there are. */
  rwt = (s32 *)lxmalloc(M->numCols*sizeof(s32),1);
  memset(rwt, 0x00, M->numCols*sizeof(s32));

  nR = 0;
  for (j=0; j<M->numCols ; j++) {
    C.numPrimes = 0;
    readColSF(&C, fp);
    for (i=0; i<C.numPrimes; i++) {
      r = C.rows[i];
      rwt[r] += 1;
      nR = MAX(nR, r);
    }
  }
  M->numRows = nR + 1 + 64; /* '1' b/c we count from zero. '64' for the additional QCB/sign bits. */
  M->numDenseBlocks = 1;
  M->denseBlockIndex[0] = M->numRows - 64;
  bSize=64;

  /* Condition for a block of rows to be considered dense: */
  dW=2.0;

  printf("Matrix scanned: it should be %ld x %ld.\n", M->numRows, M->numCols);
  i=0;
  /* This could be made a bit slicker, but whatever. */
  while ((i<nR-bSize) && (M->numDenseBlocks < MAX_DENSE_BLOCKS)) {
    /* Is the block of rows [i, i+bSize] dense? */
    for (j=0, blockWt=0; j<bSize; j++) 
      blockWt += rwt[i+j]; 
    if (blockWt > dW*M->numCols) {
      /* Check: Is this particular row fairly dense? */
      if (rwt[i] > 0.005*M->numCols) {
        /* Okay - this is a dense block. */
        M->denseBlockIndex[M->numDenseBlocks] = i;
        M->numDenseBlocks += 1;
        /* rwt is reused below, to check which entries are from dense blocks. */
        for (j=0; j<bSize; j++)
          rwt[i+j]=1; 
        i += bSize;
      } else {
        rwt[i++]=0;
      }
    } else {
      rwt[i++]=0;
    }
  }
  printf("Found %ld dense blocks. Re-reading matrix...\n", M->numDenseBlocks);
printf("The dense blocks consist of the following sets of rows:\n");
for (k=0; k<M->numDenseBlocks; k++) 
printf("[%ld, %ld]\n", M->denseBlockIndex[k], M->denseBlockIndex[k]+bSize-1);

  rewind(fp);
  fread(&M->numCols, sizeof(s32), 1, fp);

  /* We could do a little better, by discounting for the QCB and sign entries. */
  M->maxDataSize = 256 + fileSize/sizeof(s32);
  M->cEntry = (s32 *)lxmalloc(M->maxDataSize*sizeof(s32),1);
  M->cIndex = (s32 *)lxmalloc((M->numCols+1)*sizeof(s32),1);
  for (i=0; i<M->numDenseBlocks; i++) {
    if (!(M->denseBlocks[i] = (u64 *)lxcalloc((M->numCols+1)*sizeof(u64),0))) {
      fprintf(stderr, "loadMat() Error allocating %ld bytes for the QCB entries!\n",
              (M->numCols+1)*sizeof(u64));
      free(M->cIndex); free(M->cEntry); fclose(fp); return -1;
    }
  }

  M->cIndex[0] = 0; index=0;
  for (j=0; (j<M->numCols) && ((index+50) < M->maxDataSize); j++) {
    C.numPrimes = 0;
    readColSF(&C, fp);
    for (i=0; i<C.numPrimes; i++) {
      r = C.rows[i];
      /* Is this entry in a dense block? */
      if (rwt[r]==1) {
        /* Find the block (it won't be QCB, though). */
        for (k=64/bSize; k<M->numDenseBlocks; k++) {
          if ((r>=M->denseBlockIndex[k]) && (r < bSize+M->denseBlockIndex[k])) {
            M->denseBlocks[k][j] ^= BIT64(r-M->denseBlockIndex[k]);
            break;
          }
        }
      } else {
        M->cEntry[index] = r;
        index++;
      }
    }
    M->denseBlocks[0][j] = (C.QCB[0]^(((u64)C.QCB[1])<<32))|0x8000000000000000ULL;

    /* We have enough room for one extra index, and we use it */
    /* to determine size, so this is ok even on the last pass */
    /* through the loop:                                      */
    M->cIndex[j+1] = index;
  }
  fclose(fp);
  free(rwt);
  return M->numRows;
}    


/*********************************************************************/
rel_list *getRelList(multi_file_t *prelF, int index)
/*********************************************************************/
/* Allocate for and read in the specified relation file. Caller is   */
/* obviously responsible for freeing the memory when done!           */
/*********************************************************************/
{ rel_list *RL;
  FILE     *fp;
  char      fName[256];
  struct stat fileInfo;

  RL = (rel_list *)lxmalloc(sizeof(rel_list),1);
  RL->maxDataSize = 0;
  sprintf(fName, "%s.%d", prelF->prefix, index);
  
  if (stat(fName, &fileInfo)) {
    printf("Could not stat file %s!\n", fName);
    free(RL); return NULL;
  }
  RL->maxDataSize = 4096 + fileInfo.st_size/sizeof(s32);
  if ((fp = fopen(fName, "rb"))) {
    readRawS32(&RL->maxRels, fp);
    fclose(fp);
  }
  RL->maxRels += 5;
  /* Now allocate for the relations. */
  if (!(RL->relData = (s32 *)lxmalloc(RL->maxDataSize * sizeof(s32),0))) {
    fprintf(stderr, "Error allocating %ldMB for reading relation list!\n",
            RL->maxDataSize * sizeof(s32)/1048576);
    free(RL); return NULL;
  }
  if (!(RL->relIndex = (s32 *)lxmalloc(RL->maxRels * sizeof(s32),0))) {
    fprintf(stderr, "Error allocating %ldMB for relation pointers!\n", 
            RL->maxRels * sizeof(s32)/1048576);
    free(RL->relData); free(RL);
    return NULL;
  }
  RL->numRels = 0;
  readRelList(RL, fName);
  return RL;
}

/*********************************************************************/
void clearRelList(rel_list *RL)
/*********************************************************************/
{
  if (RL->relData != NULL) 
    free(RL->relData);
  if (RL->relIndex != NULL)
    free(RL->relIndex);
  RL->relData = RL->relIndex = NULL;
  RL->maxDataSize = RL->maxRels = 0;
}

/*********************************************************************/
int removeEvens(s32 *list, s32 size)
/*********************************************************************/
/* Given a list of integers: list[0], ..., list[size-1], reduce it   */
/* so it contains only those integers that appeared an odd number    */
/* of times. That is, remove the elements that occur an even number  */
/* of times.                                                         */
/* Return value: The number of elements in the reduced list.         */
/*********************************************************************/
/* This is a crappy way to do it, but it's good for now.             */
{ s32 j, k, unique;

  if (size <= 1) return size;
  qsort(list, size, sizeof(s32), cmpS32s);
  j=k=unique=0;
  while (j<size) {
    k=1;
    while (((j+k)<size) && (list[j] == list[j+k]))
      k++;
    if (k&0x01) {
      /* It occurrs an odd number of times, so keep it. */
      list[unique++] = list[j];
    }
    j += k;
  }
  return unique;
}

/*************************************************/
s32 sortRMDups(s32 *L, s32 size)
/*************************************************/
/* Sort an array of s32s and remove duplicates. */
/*************************************************/
{ s32 i, unique;

  if (size <= 1) return size;
  qsort(L, size, sizeof(s32), cmpS32s);
  i=1; unique=0;
  while (i<size) {
    if (L[i] == L[unique]) i++;
    else {
      ++unique;
      L[unique] = L[i];
      i++; 
    }
  }
  return unique+1;
}

/**********************************************************/
s32 sortRMDups2(s32 *L, s32 size)
/**********************************************************/
/* Sort an array of pairs of s32s and remove duplicates. */
/**********************************************************/
{ s32 i, unique;

  if (size <= 1) return size;
  qsort(L, size, 2*sizeof(s32), cmp2S32s);
  i=1; unique=0;
  while (i<size) {
    if ((L[2*i] == L[2*unique]) && (L[2*i+1]==L[2*unique+1])) 
      i++;
    else { 
      unique++;
      L[2*unique] = L[2*i];
      L[2*unique+1] = L[2*i+1];
      i++;
    }
  }
  return unique+1;
}



#define LP_LIST_INC_SIZE 1048576
/******************************************************************************/
llist_t *getLPList(multi_file_t *prelF)
/******************************************************************************/
/* Get the list of (indexed) primes corresponding to each relation. The return*/
/* value is the structure to be passed to combParts which will then combine   */
/* the partial relations to form full relation-sets.                          */
/******************************************************************************/
/* This is done as follows:                                                   */
/* (1) Make one pass through each prel file to index the large primes.        */
/* (2) Make a second pass through each prel file, grabbing the large primes   */
/*     occuring in it and looking up the index.                               */
/******************************************************************************/
{ s32 i, j, k, p, r, numRels=0, totalLP=0;
  s32 *lR=NULL, *lA=NULL, lRSize, lRMax, lASize, lAMax;
  s32 *buf, bufSize, bufIndex, bufMax;
  s32 *loc, key[2], index;
  u32 sF, lrpi, lapi;
  FILE *fp;
  int  numLR, numLA;
  rel_list *RL;
  llist_t *P;

  lRSize = lRMax=0;
  lASize = lAMax=0;
  for (i=0; i<prelF->numFiles; i++) {
    printf("Loading processed file %ld/%d...", i+1, prelF->numFiles);
    fflush(stdout);
    RL = getRelList(prelF, i);
    printf("Done. Processing...\n");
    numRels += RL->numRels;
    for (j=0; j<RL->numRels; j++) {
      sF = RL->relData[RL->relIndex[j]];
      numLR = GETNUMLRP(sF);
      numLA = GETNUMLAP(sF);
      lrpi = 3 + 2*(GETNUMRFB(sF) + GETNUMAFB(sF) + GETNUMSPB(sF)) + 2;
      lapi = 3 + 2*(GETNUMRFB(sF) + GETNUMAFB(sF) + GETNUMSPB(sF)) + 2 + numLR;
      totalLP += numLR + numLA;
      for (k=0; k<numLR; k++) {
        p = RL->relData[RL->relIndex[j] + lrpi + k];
        /* So p is a large rational prime in this relation. */
        if (lRSize + 8 > lRMax) {
          lRMax += LP_LIST_INC_SIZE;
          lR = lxrealloc(lR, lRMax*sizeof(s32),0);
          if (lR == NULL) {
            printf("getLPList() : Memory allocation error for lR!\n");
            exit(-1); 
          }
        }
        lR[lRSize++] = p;
      }
      for (k=0; k<numLA; k++) {
        p = RL->relData[RL->relIndex[j] + lapi + 2*k];
        r = RL->relData[RL->relIndex[j] + lapi + 2*k + 1];
        /* So (p,r) is a large algebraic prime in this relation. */
        if (lASize + 16 > lAMax) {
          lAMax += LP_LIST_INC_SIZE;
          lA = lxrealloc(lA, 2*lAMax*sizeof(s32),0);
          if (lA == NULL) {
            printf("getLPList() : Memory allocation error for lA!\n");
            exit(-1); 
          }
        }
        lA[2*lASize] = p;
        lA[2*lASize+1] = r;
        lASize++;
      }
    }
    /* All done with this relation list. */
    clearRelList(RL);
    free(RL);
    /* Now sort what we've got and remove duplicates: */
    printf("Sorting and filtering LR..."); fflush(stdout);
    lRSize = sortRMDups(lR, lRSize);
    printf("Done.\nSorting and filtering LA..."); fflush(stdout);
    lASize = sortRMDups2(lA, lASize);
    printf("Done.\n");
    printf("Found %ld distinct large rprimes and %ld large aprimes so far.\n",lRSize, lASize); 
  }

  bufMax = IO_BUFFER_SIZE;
  if (!(buf = (s32 *)lxmalloc(bufMax*sizeof(s32),0))) {
    printf("getLPList(): Memory allocation error for buf!\n");
    exit(-1);
  }
  bufSize=0;


  printf("There are %ld large primes versus %ld relations.\n", 
          lRSize+lASize, numRels);
  initialRelations = numRels;


#define NEED_LPA_DUMP
#ifdef NEED_LPA_DUMP
  printf("Dumping large primes to lpindex.L...\n");
  bufSize=0;
  if (!(fp = fopen("lpindex.L", "wb"))) {
    printf("Error opening lpindex.L for write!\n");
    exit(-1);
  }
  for (i=0; i<lRSize; i++) {
    if (bufSize + 128 > bufMax) {
      /* Flush the buffer. */
      fwrite(buf, sizeof(s32), bufSize, fp);
      bufSize=0;
    }
    buf[bufSize++]=lR[i];
    buf[bufSize++]=-1;
    buf[bufSize++]=i;
    buf[bufSize++]=2; /* Anything is ok, I think. */
  }
  for (i=0; i<lASize; i++) {
    if (bufSize + 128 > bufMax) {
      /* Flush the buffer. */
      fwrite(buf, sizeof(s32), bufSize, fp);
      bufSize=0;
    }
    buf[bufSize++]=lA[2*i];
    buf[bufSize++]=lA[2*i+1];
    buf[bufSize++]=lRSize+i;
    buf[bufSize++]=3; /* Anything is ok, I think. */
  }
  /* Flush the buffer: */
  fwrite(buf, sizeof(s32), bufSize, fp);
  bufSize=0;
  fclose(fp);

//#define DEBUG_LPCODE
#ifdef  DEBUG_LPCODE
  bufSize=0;
  if (!(fp = fopen("lpindex.txt", "w"))) {
    printf("Error opening lpindex.txt for write!\n");
    exit(-1);
  }
  for (i=0; i<lRSize; i++) {
    fprintf(fp, "%ld\n", lR[i]);
  }
  for (i=0; i<lASize; i++) {
    fprintf(fp, "%ld, %ld\n", lA[2*i], lA[2*i+1]);
  }
  fclose(fp);
#endif

#endif 

  
  /***********************************************************************/
  /* Okay - we now have a list of the distinct large primes occurring.   */
  /* Furthermore, they now have an implicit indexing given, since they   */
  /* are sorted. If we had the RAM to spare, we could compute the index  */
  /* of a given prime much more quickly than by searching, but these     */
  /* lists are already consuming quite a bit of RAM. It's not too bad,   */
  /* though. If there are 20M distinct primes on each side, then we need */
  /* about 24c operations per lookup, for some (modest) constant c.      */
  /***********************************************************************/
  if (!(fp = fopen(TMP_FILE, "wb"))) {
    printf("Error opening temp file %s for write!\n", TMP_FILE);
    exit(-1);
  }


  for (i=0; i<prelF->numFiles; i++) {
    printf("Loading processed file %ld/%d...", i+1, prelF->numFiles);
    fflush(stdout);
    RL = getRelList(prelF, i);
    printf("Done. Processing...\n");
    for (j=0; j<RL->numRels; j++) {
      if (bufSize + 128 > bufMax) {
        /* Flush the output buffer. */
        fwrite(buf, sizeof(s32), bufSize, fp);
        bufSize=0;
      }
      sF = RL->relData[RL->relIndex[j]];
      numLR = GETNUMLRP(sF);
      numLA = GETNUMLAP(sF);
      lrpi = 3 + 2*(GETNUMRFB(sF) + GETNUMAFB(sF) + GETNUMSPB(sF)) + 2;
      lapi = 3 + 2*(GETNUMRFB(sF) + GETNUMAFB(sF) + GETNUMSPB(sF)) + 2 + numLR;
      buf[bufSize++]=numLR+numLA;
      for (k=0; k<numLR; k++) {
        p = RL->relData[RL->relIndex[j] + lrpi + k];
        /* So p is a large rational prime in this relation: get it's index. */
        loc = bsearch(&p, lR, lRSize, sizeof(s32), cmpS32s);
        if (loc==NULL) {
          printf("Warning: Could not find large rational prime %ld in lR!\n", p);
          index = BAD_LP_INDEX; /* See the note at top of file. */
        } else {
          index = loc-lR;
        }
        buf[bufSize++]=index;
      }
      for (k=0; k<numLA; k++) {
        p = RL->relData[RL->relIndex[j] + lapi + 2*k];
        r = RL->relData[RL->relIndex[j] + lapi + 2*k + 1];
        /* So (p,r) is a large algebraic prime in this relation. */
        key[0]=p; key[1]=r;
        loc = bsearch(key, lA, lASize, 2*sizeof(s32), cmp2S32s);
        if (loc==NULL) {
          printf("Warning: Could not find large alg prime (%ld,%ld) in lR!\n",p,r);
          index = BAD_LP_INDEX;
        } else {
          index = lRSize + (loc - lA)/2;
        }
        buf[bufSize++]=index;
      }
    }
    clearRelList(RL);
    free(RL);
  }
  /* Flush the buffer: */
  fwrite(buf, sizeof(s32), bufSize, fp);
  bufSize=0;
  fclose(fp);
  totalLargePrimes = lRSize + lASize;
  free(lR);
  free(lA);

  P = (llist_t *)lxmalloc(sizeof(llist_t),1);
  if (ll_init(P, totalLP, numRels)) {
    printf("getLPList(): Memory allocation error doing ll_init()!\n");
    exit(-1);
  }
  if (!(fp = fopen(TMP_FILE, "rb"))) {
    printf("Error re-opening temp file!\n");
    exit(-1);
  }
  bufSize=bufIndex=0;
  while (!(feof(fp)) || (bufIndex < bufSize)) {
    if ((bufIndex + 128 > bufSize)&& !(feof(fp))) {
      /* Read the next block into the buffer: */
      memmove(buf, buf+bufIndex, (bufSize-bufIndex)*sizeof(s32));
      bufSize -= bufIndex;
      bufIndex = 0;
      bufSize += fread(buf+bufSize, sizeof(s32), (bufMax-bufSize), fp);
    }
    ll_appendField(P, buf+bufIndex+1, buf[bufIndex]);
    bufIndex += 1+buf[bufIndex];
  }
  remove(TMP_FILE);
  free(buf);
  fclose(fp);
  return P;
}

/******************************************************************************/
s32 doRowOps3(llist_t **P, llist_t *R, multi_file_t *prelF, int maxRelsInFF)
/******************************************************************************/
{ s32 numFull, numLP[10], j;
  llist_t *PL;
  s32 numLargeP;

  PL = getLPList(prelF);
  numLargeP = totalLargePrimes;
  printf("----------------------------\n");
  printf("There are %ld large primes versus %ld relations.\n", 
          numLargeP, initialRelations);
  msgLog(NULL, "largePrimes: %ld , relations: %ld", numLargeP, initialRelations);
  printf("----------------------------\n");

  *P = PL;
  for (j=0; j<10; j++)
    numLP[j]=0;
  for (j=0; j<PL->numFields; j++) 
    numLP[PL->index[j+1]-PL->index[j]] += 1;
  printf("Num large primes  |  Relations\n");
  printf("------------------------------\n");
  j=0;
  do {
    printf("                %ld | %ld\n", j, numLP[j]);
    j++;
  } while ((j<10) && (numLP[j]>0));
  printf("------------------------------\n");


  numFull = combParts(R, PL, maxRelsInFF, minFF);
  return numFull;
}



/*********************************************************************/
s32 getCols(char *colName, multi_file_t *prelF, multi_file_t *lpF, nfs_fb_t *FB,
             s32 minFull, int maxRelsInFF)
/*********************************************************************/
{ s32         i, j, k, l, R0, R1;
  s32         numFulls, rIndex, numFF, tPP;
  s32         aOffset, spOffset;
  char         fName[64], prelName[64], indexName[64];
  FILE         *fp, *ofp, *ofp2;
  column_t     C;
  rel_list    *RL;
  relation_t   R;
  mpz_t       tmp;
  llist_t     *P, Rl;
  s32       *buf=NULL, bufSize, bufMax, bufIndex;
  s32       *buf2=NULL, bufSize2, bufMax2, bufIndex2;

  mpz_init(tmp);
  tPP = approxPi_x(FB->maxP_r) + approxPi_x(FB->maxP_a);
  tPP = tPP - FB->rfb_size - FB->afb_size;
  printf("Max # of large primes is approximately %ld.\n", tPP);
  numFulls = doRowOps3(&P, &Rl, prelF, maxRelsInFF);

  if (numFulls < minFull) {
    ll_clear(P); ll_clear(&Rl); free(P);
    return numFulls;
  }
  /**************************************************************/
  /* Convert the full relations into column_t data.             */
  /* This is done in several passes. First, create the column   */
  /* file, with only the relations. Later, we will go through   */
  /* and replace the relations with their corresponding primes. */
  /**************************************************************/
  C.numPrimes = 0; 
  C.QCB[0] = C.QCB[1] = 0x00000000;
  sprintf(indexName, "%s.index", colName);
  if (!(ofp = fopen(colName, "wb"))) {
    fprintf(stderr, "getCols() Error opening %s for write!\n", colName);
    exit(-1);
  }
  if (!(ofp2 = fopen(indexName, "wb"))) {
    fprintf(stderr, "getCols() Error opening %s for write!\n", indexName);
    fclose(ofp);
    exit(-1);
  }
  /* Prep the output buffers: */
  bufMax = IO_BUFFER_SIZE;
  if (!(buf = (s32 *)lxmalloc(bufMax*sizeof(s32),0))) {
    printf("getCols() : Memory allocation error for buf!\n");
    exit(-1);
  }
  bufMax2 = IO_BUFFER_SIZE;
  if (!(buf2 = (s32 *)lxmalloc(bufMax2*sizeof(s32),0))) {
    printf("getCols() : Memory allocation error for buf2!\n");
    exit(-1);
  }
  bufSize = bufSize2=0;
  /* Rl should contain just full relation-sets now. */
  for (i=0,numFF=0; i<Rl.numFields; i++) {
    if ((Rl.index[i+1]-Rl.index[i]>0) && (P->index[i+1]==P->index[i])) {
      /* This is a full relation, coming from several (a,b) pairs. */
      C.numRels = Rl.index[i+1] - Rl.index[i];
      for (j=0; j<C.numRels; j++) {
        C.Rels[j] = Rl.data[Rl.index[i]+j];
#ifdef _DEBUG
        if ((C.Rels[j] < 0) || (C.Rels[j] > initialRelations)) {
          printf("Error: C.Rels[%ld] = %ld versus %ld relations!\n",j,C.Rels[j],initialRelations);
          exit(-1);
        }
        if (C.numRels > MAX_RELS_IN_FF) {
          printf("Error: relation-set %ld (i=%ld) has too many relations (%ld)!\n",
                  numFF, i, C.numRels);
          exit(-1);
        }
#endif
      }
      /* Write the column out in LF format. */
      if (bufSize +2048 > bufMax) {
        fwrite(buf, sizeof(s32), bufSize, ofp);
        bufSize=0;
      }
      buf[bufSize++] = C.numRels;
      buf[bufSize++] = C.numPrimes;
      memcpy(buf+bufSize, &C.Rels, C.numRels*sizeof(s32));
      bufSize += C.numRels;
      buf[bufSize++] = C.QCB[0];
      buf[bufSize++] = C.QCB[1];
      memcpy(buf+bufSize, &C.rows, C.numPrimes*sizeof(s32));
      bufSize += C.numPrimes;
      /* And the column index info: */
      if (bufSize2 + 2048 > bufMax2) {
        fwrite(buf2, sizeof(s32), bufSize2, ofp2);
        bufSize2=0;
      }
      buf2[bufSize2++] = C.numRels;
      memcpy(buf2+bufSize2, &C.Rels, C.numRels*sizeof(s32));
      bufSize2 += C.numRels;
      numFF++;
    }
  }
  /* Flush buffers. */
  fwrite(buf, sizeof(s32), bufSize, ofp);
  bufSize=0;
  fwrite(buf2, sizeof(s32), bufSize2, ofp2);
  bufSize2=0;

  fclose(ofp); fclose(ofp2);
  /* We won't be needing these anymore. */
  ll_clear(&Rl);
  ll_clear(P); free(P);


  printf("After re-scanning files and building column indicies, numFF=%ld.\n", numFF);
  bufSize = bufIndex = 0;
  bufSize2 = bufIndex2 = 0;
  /*******************************************************************************/
  /* We now have a completely unprocessed column file. It has only the relations */
  /* occurring in each column. We need to go through and replace each relation   */
  /* with the corresponding primes, QCB bits.                                    */
  /*******************************************************************************/

  aOffset = FB->rfb_size;
  spOffset = aOffset + FB->afb_size;

  printf("Creating %ld matrix columns...\n", numFF);
  strcpy(fName, TMP_FILE);
  R0=R1=0;
  for (i=0; i<prelF->numFiles; i++) {
    sprintf(prelName, "%s.%ld", prelF->prefix, i);
    RL = getRelList(prelF, i);
    R1 = R0 + RL->numRels;
    printf("Re-read %ld relations from %s :  [%ld, %ld).\n", RL->numRels, prelName, R0, R1);
    /* Now, we have in RAM the relations numbered [R0, R1). */
    if (!(fp = fopen(colName, "rb"))) {
      fprintf(stderr, "getCols() Error opening %s for read!\n", colName);
      exit(-1);
    }
    if (!(ofp = fopen(fName, "wb"))) {
      fprintf(stderr, "getCols() Error opening %s for write!\n", fName);
      fclose(fp);
      exit(-1);
    }
    /* Read the first block of data into the input buffer: */
    bufSize2 = fread(buf2, sizeof(s32), bufMax2, fp);
    bufIndex2 = 0;


    for (j=0; j<numFF; j++) {
      /* Grab the next col in LF format from the input buffer. */
      if ((bufIndex2 + 8192 > bufSize2) && !(feof(fp))) {
        memmove(buf2, buf2+bufIndex2, (bufSize2-bufIndex2)*sizeof(s32));
        bufSize2 -= bufIndex2;
        bufIndex2=0;
        bufSize2 += fread(buf2+bufSize2, sizeof(s32), bufMax2-bufSize2, fp);
      }
      C.numRels = buf2[bufIndex2++];
      C.numPrimes = buf2[bufIndex2++];
      memcpy(&C.Rels, buf2+bufIndex2, C.numRels*sizeof(s32));
      bufIndex2 += C.numRels;
      C.QCB[0] = buf2[bufIndex2++];
      C.QCB[1] = buf2[bufIndex2++];
      memcpy(&C.rows, buf2+bufIndex2, C.numPrimes*sizeof(s32));
      bufIndex2 += C.numPrimes;
      if (bufIndex2 > bufSize2) {
        printf("getCols(): Unknown error: data is probably corrupt!\n");
        exit(-1);
      }
      for (k=0; k<C.numRels; k++) {
        if ((C.Rels[k]>= R0) && (C.Rels[k] < R1)) {
          /* We have the factorization of this relation in RAM. */
          rIndex = C.Rels[k] - R0;
          for (l=k; l<C.numRels; l++)
            C.Rels[l] = C.Rels[l+1];
          C.numRels -= 1;
          /****************************************************************/        
          /* This is horribly inefficient. If we get ambitious, we should */
          /* replace this with code that pulls the factors directly out   */
          /* of the data. But for now, it'll do.                          */
          /****************************************************************/        
          dataConvertToRel(&R, &RL->relData[RL->relIndex[rIndex]]);
          /* Note: These '%2's cannot be `&0x01' b/c negative exponents can occur!   */
          if (C.numPrimes + R.rFSize + R.aFSize + R.spSize >= MAX_ROWS_IN_COL) {
            fprintf(stderr, "getCols(): MAX_ROWS_IN_COL exceeded!\n");
            exit(-1);
          }
          for (l=0; l<R.rFSize; l++) 
            if (R.rExps[l]%2) {
              C.rows[C.numPrimes] = R.rFactors[l]; C.numPrimes += 1;
            }
          for (l=0; l<R.aFSize; l++) 
            if (R.aExps[l]%2) {
              C.rows[C.numPrimes] = aOffset+ R.aFactors[l]; C.numPrimes += 1;
            }
          for (l=0; l<R.spSize; l++) 
            if (R.spExps[l]%2) {
              C.rows[C.numPrimes] = spOffset+ R.spFactors[l]; C.numPrimes += 1;
            }
          C.QCB[0] ^= R.qcbBits[0]; C.QCB[1] ^= R.qcbBits[1];
          C.numPrimes = removeEvens(C.rows, C.numPrimes);
          k--; /* Since we removed this relation and shifted everything down,
                  we need this so that after the increment, we look at the
                  same position (which will be the next relation index, if there
                  is another.)
               */
        } 
#ifdef _DEBUG
        else if ((C.Rels[k] < 0) || (C.Rels[k] > initialRelations)) {
          printf("Error: intermediate col file is corrupt!\n");
          printf("For relation-set %ld:\n", j);
          printf("C.numRels = %ld\n", C.numRels);
          printf("C.Rels[%ld] = %ld\n", k, C.Rels[k]);
          exit(-1);
        }
#endif
      }
      /* Write out the column in LF format. */
      if (bufSize + 2048 > bufMax) {
        fwrite(buf, sizeof(s32), bufSize, ofp);
        bufSize = 0;
      }
      buf[bufSize++] = C.numRels;
      buf[bufSize++] = C.numPrimes;
      memcpy(buf+bufSize, &C.Rels, C.numRels*sizeof(s32));
      bufSize += C.numRels;
      buf[bufSize++] = C.QCB[0];
      buf[bufSize++] = C.QCB[1];
      memcpy(buf+bufSize, &C.rows, C.numPrimes*sizeof(s32));
      bufSize += C.numPrimes;
    }
    R0 += RL->numRels;
    /* Flush the output buffer. */
    fwrite(buf, sizeof(s32), bufSize, ofp);
    bufSize=0;
    fclose(fp); fclose(ofp);
    remove(colName); rename(fName, colName);
    clearRelList(RL);
    free(RL);
  }

  /* Finally, go back through and clean up the column file - now that none */
  /* of the entries have any un-resolved relations left, all we need is    */
  /* the prime-indexing fields.                                            */
  if (!(fp = fopen(colName, "rb"))) {
    fprintf(stderr, "getCols() Error opening %s for read!\n", colName);
    exit(-1);
  }
  strcpy(fName, TMP_FILE);
  if (!(ofp = fopen(fName, "wb"))) {
    fprintf(stderr, "getCols() Error opening %s for write!\n", TMP_FILE);
    fclose(fp);
    exit(-1);
  }

  bufSize=0;
  bufSize2=bufIndex2=0;
  fwrite(&numFF, sizeof(s32), 1, ofp);
  /* Read the first block. */
  bufSize2 = fread(buf2, sizeof(s32), bufMax2, fp);

  for (j=0; j<numFF; j++) {
    /* Grab the next col in LF format from the input buffer. */
    if ((bufIndex2 + 1024 > bufSize2) && !(feof(fp))) {
      memmove(buf2, buf2+bufIndex2, (bufSize2-bufIndex2)*sizeof(s32));
      bufSize2 -= bufIndex2;
      bufIndex2=0;
      bufSize2 += fread(buf2+bufSize2, sizeof(s32), bufMax2-bufSize2, fp);
    }
    C.numRels = buf2[bufIndex2++];
    C.numPrimes = buf2[bufIndex2++];
    memcpy(&C.Rels, buf2+bufIndex2, C.numRels*sizeof(s32));
    bufIndex2 += C.numRels;
    C.QCB[0] = buf2[bufIndex2++];
    C.QCB[1] = buf2[bufIndex2++];
    memcpy(&C.rows, buf2+bufIndex2, C.numPrimes*sizeof(s32));
    bufIndex2 += C.numPrimes;
    if (bufIndex2 > bufSize2) {
      printf("getCols(): Unknown error: data is probably corrupt!\n");
      exit(-1);
    }

    if (C.numRels > 0) {
      printf("Error: relation-set %ld still has %ld unconverted relations!\n",
              j,C.numRels);
      printf("They are: ");
      for (i=0; i<C.numRels; i++) 
        printf("%ld ", C.Rels[i]);
      printf("\n");
      exit(-1);
    }
    /* Write the column in SF format to the output buffer. */
    buf[bufSize++] = C.numPrimes;
    buf[bufSize++] = C.QCB[0];
    buf[bufSize++] = C.QCB[1];
    memcpy(buf+bufSize, &C.rows, C.numPrimes*sizeof(s32));
    bufSize += C.numPrimes;
    if (bufSize + 1024 > bufMax) {
      fwrite(buf, sizeof(s32), bufSize, ofp);
      bufSize=0;
    }
  }
  fwrite(buf, sizeof(s32), bufSize, ofp);
  fclose(fp); fclose(ofp);
  remove(colName); rename(fName, colName);
  free(buf); free(buf2);
  mpz_clear(tmp);
  return numFF;
}

/******************************************************/
int count_prelF(multi_file_t *prelF)
/******************************************************/
/* Check to see how many files there are, and whether */
/* or not this needs to be increased. If so, increase.*/
/* Well, the increase will only be done if takeAction */
/* is nonzero - otherwise, do nothing but count them. */
/******************************************************/
{ int    i, cont;
  struct stat fileInfo;
  char   fName[512];
  s32   maxSize=0;

  i=0;
  do {
    cont=0;
    sprintf(fName, "%s.%d", prelF->prefix, i);
    printf("Checking file %s ...\n", fName);
    if (stat(fName, &fileInfo)==0) {
      maxSize = MAX(maxSize, fileInfo.st_size);
      cont=1;
      i++;
    }
  } while (cont);
/* CJM, 129/04 : Consider making this MAX(i, DEFAULT_NUM_FILES); */
  prelF->numFiles = MAX(i, 1);
  return 0;
}


/******************************************************/
int allocateRL(multi_file_t *prelF, rel_list *RL)
/******************************************************/
/* Allocate 'RL' so it can hold the largest of the    */
/* processed relation files.                          */
/******************************************************/
{ s32 maxSize;
  char prelName[512];
  int  i;
  struct stat fileInfo;

  maxSize = 0;
  for (i=0; i<prelF->numFiles; i++) {
    sprintf(prelName, "%s.%d", prelF->prefix, i);
    if (stat(prelName, &fileInfo)==0) 
      maxSize = MAX(maxSize, fileInfo.st_size);
  }

  RL->numRels = 0;
  RL->maxDataSize = 1000 + maxSize/sizeof(s32);
  if (!(RL->relData = (s32 *)lxmalloc(RL->maxDataSize * sizeof(s32),0))) {
    fprintf(stderr, "Error allocating %ldMB for processed relation files!\n",
            RL->maxDataSize * sizeof(s32)/1048576);
    fprintf(stderr, "Try decreasing DEFAULT_MAX_FILESIZE and re-running.\n");
    exit(-1);
  }
  /* Again: it's a safe bet that any relation needs at least 20 s32s, so: */
  RL->maxRels = (u32)RL->maxDataSize/20;
  if (!(RL->relIndex = (s32 *)lxmalloc(RL->maxRels * sizeof(s32),0))) {
    fprintf(stderr, "Error allocating %ldMB for relation pointers!\n",
            RL->maxRels * sizeof(s32)/1048756);
    free(RL->relData);
    exit(-1);
  }
  return 0;
}

/******************************************************/
void clearRL(rel_list *RL)
{
  if (RL->relData) free(RL->relData);
  if (RL->relIndex) free(RL->relIndex);
  RL->maxDataSize = RL->numRels = 0;
  RL->relData = RL->relIndex = NULL;
}

/********************************************/
int cmp2S32s(const void *a, const void *b)
/********************************************/
{ s32 *A = (s32 *)a, *B = (s32 *)b;

  if (A[0] < B[0]) return -1;
  if (A[0] > B[0]) return 1;
  if (A[1] < B[1]) return -1;
  if (A[1] > B[1]) return 1;
  return 0;
}

/****************************************************/
int buildSparseMat(char *colName, double wtFactor)
/****************************************************/
{ llist_t C;
  nfs_sparse_mat_t M;
  long minExtraCols;

  printf("Loading matrix into RAM...\n");
  loadMat(&M, colName);
  printf("Matrix loaded: it is %ld x %ld.\n", M.numRows, M.numCols);
  if (M.numCols < (M.numRows + 64)) {
    printf("More columns needed (current = %ld, min = %ld)\n",
           M.numCols, M.numRows+64);
    free(M.cEntry); free(M.cIndex);
    return 0;
  }
  if (checkMat(&M)) {
    printf("checkMat() returned some error! Terminating...\n");
    exit(-1);
  }
  /**************************************************************/
  /* This is to map new columns onto old ones, so we can delete */
  /* some, and still be able to recover dependencies from the   */
  /* original matrix.                                           */
  /**************************************************************/
  cm_init(&C, &M);

  /* First, delete any columns which we are required to remove (probably
     because they were somehow bad.
  */
  removeCols(&M, &C, delCols, numDel);

  minExtraCols = 64 + 0.005*M.numRows;
  pruneMatrix(&M, minExtraCols, wtFactor, &C);
  /* Sanity check: */
  if (M.numCols != C.numFields) {
    fprintf(stderr, "Error: M.numCols = %ld != %ld = C.numFields.\n",
            M.numCols, C.numFields);
    exit(-1);
  } else {
    printf("Sanity check: M.numCols = %ld = C.numFields. passed.\n",
            M.numCols);
  }
  writeSparseMat("spmat", &M);
  ll_write("sp-index", &C);
  return 0;
}


/****************************************************/
int main(int argC, char *args[])
/****************************************************/
{ char       fbName[64], prelName[40], newRelName[64], depName[64], colName[64];
  char       str[128], line[128];
  int        i, qcbSize = DEFAULT_QCB_SIZE, seed=DEFAULT_SEED, retVal=0;
  int        maxRelsInFF=MAX_RELS_IN_FF;
  double     startTime, stopTime, wtFactor=DEFAULT_WT_FACTOR;
  s32        totalRels, relsInFile;
  nfs_fb_t   FB;
  multi_file_t prelF, lpF;
  FILE      *fp;

lxmalloc(4000000,0);
  prelF.numFiles = DEFAULT_NUM_FILES;
  lpF.numFiles = 0;
  prelF.prefix[0] = lpF.prefix[0]=0;
  fbName[0] = newRelName[0] = 0;
  strcpy(depName, DEFAULT_DEPNAME);
  strcpy(colName, DEFAULT_COLNAME);
  strcpy(lpF.prefix, DEFAULT_LPI_NAME);
  strcpy(prelF.prefix, DEFAULT_PRELPREFIX);
  line[0]=0;
  printf(START_MSG, GGNFS_VERSION);
  minFF=0;
  for (i=1; i<argC; i++) {
    if (strcmp(args[i], "-fb")==0) {
      if ((++i) < argC) 
        strncpy(fbName, args[i], 64);
    } else if (strcmp(args[i], "-prel")==0) {
      if ((++i) < argC) 
        strncpy(prelF.prefix, args[i], 32);
    } else if (strcmp(args[i], "-qs")==0) {
      if ((++i) < argC)
        qcbSize = atoi(args[i]);
    } else if (strcmp(args[i], "-minff")==0) {
      if ((++i) < argC)
        minFF = atoi(args[i]);
    } else if (strcmp(args[i], "-seed")==0) {
      if ((++i) < argC) {
        seed = atoi(args[i]);
      }
    } else if (strcmp(args[i], "-v")==0) {
      verbose++;
    } else if (strcmp(args[i], "-maxrelsinff")==0) {
      if ((++i) < argC) 
        maxRelsInFF = atoi(args[i]);
    }
  }
  maxRelsInFF=MIN(MAX_RELS_IN_FF,maxRelsInFF);
  srand(seed);
  startTime = sTime();

  if ((fbName[0]==0) || (prelF.prefix[0]==0)) {
    printf("USAGE: %s %s\n", args[0], USAGE);
    exit(0);
  }
  msgLog("", "GGNFS-%s : matbuild", GGNFS_VERSION);
  sprintf(prelName, "%s.%d", prelF.prefix, 0);
  initFB(&FB);
  if (loadFB(fbName, &FB)) {
    printf("Could not load FB from %s!\n", fbName);
    exit(-1);
  }

  if (minFF < FB.rfb_size + FB.afb_size + 64 + 32)
    minFF = FB.rfb_size + FB.afb_size + 64 + 32;
  if (verbose)
    printf("Getting QCB of size %d...\n", qcbSize);

  /* Very strange! Why is this being done here?? */
  generateQCB(&FB, qcbSize); 

  count_prelF(&prelF);

  totalRels=0;
  for (i=0; i<prelF.numFiles; i++) {
    sprintf(prelName, "%s.%d", prelF.prefix, i);
    if ((fp = fopen(prelName, "rb"))) {
      rewind(fp);
      fread(&relsInFile, sizeof(s32), 1, fp);
      fclose(fp);
    }
    totalRels += relsInFile;
  }

  finalFF = getCols("cols", &prelF, &lpF, &FB, minFF, maxRelsInFF);
  msgLog("", "Heap stats for matbuild run, after cycle-building");
  logHeapStats();
  if (finalFF > 0)
    msgLog("", "rels:%ld, initialFF:%ld, finalFF:%ld", 
           initialRelations, initialFF, finalFF);
  if (finalFF < minFF) {
    printf("More columns needed (current = %ld, min = %ld)\n",
           finalFF, minFF);
    exit(0);
  }

  /* Ok - build and prune the sparse matrix. */
  buildSparseMat("cols", wtFactor);

  /* Write some header information that will eventually be needed
     in the `deps' file.
  */
  if (!(fp = fopen("depinf", "wb"))) {
    fprintf(stderr, "Error opening %s for write!\n", "depinf");
  } else {
    sprintf(str, "NUMCOLS: %8.8lx", finalFF); writeBinField(fp, str);
    sprintf(str, "COLNAME: %s.index", colName); writeBinField(fp, str);
    sprintf(str, "MAXRELS: %8.8lx", totalRels); writeBinField(fp, str);
    sprintf(str, "RELPREFIX: %s", prelF.prefix); writeBinField(fp, str);
    sprintf(str, "RELFILES: %x", prelF.numFiles); writeBinField(fp, str);
    sprintf(str, "LPFPREFIX: %s", lpF.prefix); writeBinField(fp, str);
    sprintf(str, "LPFFILES: %x", lpF.numFiles); writeBinField(fp, str);
    sprintf(str, "END_HEADER"); writeBinField(fp, str);
    fclose(fp);
  }
  printf("`depinf' written. You can now run matsolve.\n");
  msgLog("", "depinf file written. Run matsolve.");

  stopTime = sTime();
  printf("Total elapsed time: %1.2lf seconds.\n", stopTime-startTime);
  msgLog("", "Heap stats for matbuild run:");
  logHeapStats();

  return retVal;
}  






