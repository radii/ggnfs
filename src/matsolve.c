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

#include "if.h"

#if !defined(_MSC_VER)
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
"-save <int>  : Interval (in minutes) between save files.\n"\
"-test        : Do not solve matrix; use it to test multiply operations.\n"\
"               This can help expose hardware problems or miscompilations.\n"\
"--help       : Show this help and quit.\n"

#define START_MSG \
"\n"\
" __________________________________________________________ \n"\
"|        This is the matsolve program for GGNFS.           |\n"\
"| Version: %-25s                       |\n"\
"| This program is copyright 2004, Chris Monico, and subject|\n"\
"| to the terms of the GNU General Public License version 2.|\n"\
"|__________________________________________________________|\n"





/***** Globals *****/
s32 delCols[2048], numDel=0;


/****************************************************/
int main(int argC, char *args[])
/****************************************************/
{ char       colName[64], depName[64], str[1024];
  double     startTime, stopTime;
  s32       *deps, origC, seed=DEFAULT_SEED;
  long       testMode=0;
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
    } else if (strcmp(args[i], "-save")==0) {
      if ((++i) < argC) {
        matsave_interval = 60 * atoi(args[i]);
      }
    } else if (strcmp(args[i], "-test")==0) {
      testMode = 1;
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
  if (checkMat(&M, delCols, &numDel)) {
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

  if (getDependencies(&M, &C, deps, origC, testMode) == 0) {
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
