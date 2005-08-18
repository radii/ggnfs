/**************************************************************/
/* matsolve.c                                                 */
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
#include <sys/stat.h>


#ifndef _MSC_VER 
#include <sys/time.h>
#endif
#include "ggnfs.h"

#define DEFAULT_SEED 1
#define DEFAULT_DEPNAME "deps"
#define DEFAULT_COLNAME "cols"
#define DEFAULT_WT_FACTOR 0.7


#define USAGE "[OPTIONS]\n"\
"-v           : verbose.\n"\
"-seed <int>  : Set the seed for the PRNG.\n"\
"--help       : Show this help and quit.\n"

#define START_MSG \
"\n"\
" __________________________________________________________ \n"\
"|        This is the matsolve program for GGNFS.           |\n"\
"| Version: %-25s                            |\n"\
"| This program is copyright 2004, Chris Monico, and subject|\n"\
"| to the terms of the GNU General Public License version 2.|\n"\
"|__________________________________________________________|\n"


/***********************************************************/
/* This structure is for columns that are being processed. */
/* Before processing, we know only the relations that      */
/* comprise each column. After processing, we know all the */
/* nonzero rows of the matrix (i.e., the primes appearing  */
/* with odd exponent and the QCB values and sign value).   */
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
  char sign;
  s32 rows[MAX_ROWS_IN_COL];
} column_t;




/***** Globals *****/
int verbose=0;
s32 delCols[2048], numDel=0;


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

  if (!(colHash = (s32 *)malloc(2*M->numCols*sizeof(s32)))) {
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
        printf("Warning: column %" PRId32 " is all zero!\n", c);
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
        printf("Bad matrix: column %" PRId32 " = column %" PRId32 "!\n", 
               colHash[2*i+1], colHash[2*i+3]);
        if (numDel < 2048) 
          delCols[colHash[2*i+3]]=c;
        
        warn=2;
      }
    }
  }

  if (warn) {
    printf("checkMat() did not like something about the matrix:\n");
    printf("This is probably a sign that something has gone horribly wrong\n");
    printf("in the matrix construction (matbuild).\n");
    if (numDel < 2048) {
      printf("However, the number of bad columns is only %" PRId32 ",\n", numDel);
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
  rwt = (s32 *)malloc(M->numCols*sizeof(s32));
  if (rwt == NULL) {
    fprintf(stderr, "loadMat(): Memory allocation error for rwt! (%" PRIu32 " bytes)\n",
            (u32)(M->numCols*sizeof(s32)) );
    exit(-1);
  }
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

  printf("Matrix scanned: it should be %" PRId32 " x %" PRId32 ".\n", M->numRows, M->numCols);
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
  printf("Found %" PRId32 " dense blocks. Re-reading matrix...\n", M->numDenseBlocks);
  printf("The dense blocks consist of the following sets of rows:\n");
  for (k=0; k<M->numDenseBlocks; k++) 
    printf("[%" PRId32 ", %" PRId32 "]\n", M->denseBlockIndex[k], M->denseBlockIndex[k]+bSize-1);

  rewind(fp);
  fread(&M->numCols, sizeof(s32), 1, fp);

  /* We could do a little better, by discounting for the QCB and sign entries. */
  M->maxDataSize = 256 + fileSize/sizeof(s32);
  if (!(M->cEntry = (s32 *)malloc(M->maxDataSize*sizeof(s32)))) {
    fprintf(stderr, "loadMat() Error allocating %" PRIu32 " bytes for the sparse matrix!\n",
            (u32)(M->maxDataSize*sizeof(s32)) );
    fclose(fp); return -1;
  }
  if (!(M->cIndex = (s32 *)malloc((M->numCols+1)*sizeof(s32)))) {
    fprintf(stderr, "loadMat() Error allocating %" PRIu32 " bytes for the sparse matrix indicies!\n",
            (u32)((M->numCols+1)*sizeof(s32)) );
    free(M->cEntry); fclose(fp); return -1;
  }
  for (i=0; i<M->numDenseBlocks; i++) {
    if (!(M->denseBlocks[i] = (u64 *)calloc((M->numCols+1),sizeof(u64)))) {
      fprintf(stderr, "loadMat() Error allocating %" PRIu32 " bytes for the QCB entries!\n",
              (u32)((M->numCols+1)*sizeof(u64)) );
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


/****************************************************/
int main(int argC, char *args[])
/****************************************************/
{ char       colName[64], depName[64], str[1024];
  double     startTime, stopTime;
  s32       *deps, origC, seed=DEFAULT_SEED;
  struct stat fileInfo;
  nfs_sparse_mat_t M;
  llist_t    C;
  int        i;
  FILE      *fp, *ifp;

  strcpy(colName, DEFAULT_COLNAME);
  strcpy(depName, DEFAULT_DEPNAME);
  printf(START_MSG, GGNFS_VERSION);
  seed=time(0);
  /* This probably shouldn't be needed, but whatever. */
  seed = ((seed % 1001)*seed) ^ (171*seed);

  for (i=1; i<argC; i++) {
    if (strcmp(args[i], "-v")==0) {
      verbose++;
    } else if (strcmp(args[i], "-seed")==0) {
      if ((++i) < argC) {
        seed = atol(args[i]);
      }
    } else if (strcmp(args[i], "--help")==0) {
      printf("USAGE: %s %s\n", args[0], USAGE);
      exit(0);
    }
  }
  srand(seed);
  if (stat("depinf", &fileInfo)) {
    printf("Could not stat depinf! Are you trying to run %s to soon?\n", args[0]);
    return -1;
  }
  seedBlockLanczos(seed);
  startTime = sTime();
  msgLog("", "GGNFS-%s : matsolve (seed=%" PRId32 ")", GGNFS_VERSION, seed);
  printf("Using PRNG seed=%" PRId32 ".\n", seed);


  readSparseMat(&M, "spmat");
  ll_read(&C, "sp-index");
  printf("Verifying column map..."); fflush(stdout);
  ll_verify(&C);
  printf("done.\n");


  printf("Matrix loaded: it is %" PRId32 " x %" PRId32 ".\n", M.numRows, M.numCols);
  if (M.numCols < (M.numRows + 64)) {
    printf("More columns needed (current = %" PRId32 ", min = %" PRId32 ")\n",
           M.numCols, M.numRows+64);
    free(M.cEntry); free(M.cIndex);
    exit(-1);
  }
  if (checkMat(&M)) {
    printf("checkMat() returned some error! Terminating...\n");
    exit(-1);
  }

  /* We need to know how many columns there were in the original, unpruned
     matrix, so we know how much memory to allocate for the dependencies.
  */
  if (!(ifp = fopen("depinf", "rb"))) {
    fprintf(stderr, "Error opening depinf for read!\n");
    exit(-1);
  }
  readBinField(str, 1024, ifp);
  while (!(feof(ifp)) && strncmp(str, "END_HEADER",10)) {
    if (strncmp(str, "NUMCOLS: ", 9)==0) {
      sscanf(&str[9], "%" SCNx32, &origC);
    }
    readBinField(str, 1024, ifp);
  }
  fclose(ifp); 
  printf("Original matrix had %" PRId32 " columns.\n", origC);

  if (!(deps = (s32 *)malloc(origC*sizeof(s32)))) {
    printf("Could not allocate %" PRIu32 " bytes for the dependencies.\n", (u32)(origC*sizeof(s32)) );
    free(M.cEntry); free(M.cIndex); return -1;
  }

  if (getDependencies(&M, &C, deps) == 0) {
    if (!(ifp = fopen("depinf", "rb"))) {
      fprintf(stderr, "Error opening depinf for read!\n");
      exit(-1);
    }
    printf("Writing dependencies to file %s.\n", depName);
    if (!(fp = fopen(depName, "wb"))) {
      fprintf(stderr, "Error opening %s for write!\n", depName);
      fclose(ifp);
    } else {
      /* Get the header information from depinf. */
      readBinField(str, 1024, ifp);
      while (!(feof(ifp)) && strncmp(str, "END_HEADER",10)) {
        writeBinField(fp, str);
        readBinField(str,1024,ifp);
      }
      if (strncmp(str, "END_HEADER",10)) {
        fprintf(stderr, "Error: depinf is corrupt!\n");
        fclose(ifp); fclose(fp); exit(-1);
      }
      writeBinField(fp, str);
      fclose(ifp);
      fwrite(deps, sizeof(s32), origC, fp);
      fclose(fp);
    }
  }

  stopTime = sTime();
  printf("Total elapsed time: %1.2lf seconds.\n", stopTime-startTime);
  msgLog("", "Heap stats for matsolve run:");
  logHeapStats();



  free(M.cEntry); free(M.cIndex); free(deps);
  return 0;
}  
