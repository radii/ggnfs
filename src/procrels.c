/**************************************************************/
/* procrels.c                                                 */
/* Copyright 2004, Chris Monico.                              */
/* Many fixes and improvements due to KAMADA Makoto.          */
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

/* There is a hack in this file which needs to be properly fixed.
   For some reason, people are getting the occasional
   ``Could not find large xxxx prime ...'' message with exit.
   The temporary workaround is to use the value BAD_LP_INDEX
   to index that particular prime. Then, when combparts goes through
   it's ``throw out singletons'' phase, it will also keep an eye
   out for primes that are tagged with BAD_LP_INDEX and throw them
   out as well. This is a serious hack and needs to be fixed.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "if.h"

/* After reading this many new relations, they will be added to the
   larger hash table, to reduce sort/searches which would otherwise
   slow down relation processing when processing many new relaitons.
   Keep this extraordinarily large for now, because I think it's
   broken. But it should be fixed and backed down to 100K or so.
*/
#define MAX_AB_EXTRA_ENTRIES 80000


#if !defined(_MSC_VER)
#include <sys/time.h>
#endif
#include "ggnfs.h"
#include "prand.h"
#include "rellist.h"
#include "intutils.h"


/* I need to figure out what is roughly optimal
   for a machine with, say 512MB of RAM. This is
   the max file size for individual processed relation
   files. But beware - don't make it too small - some
   of the internal algorithms may be
   quadratic in the number of files! On the other hand,
   if this is too large, you'll run into VM and everything
   will slow to a crawl or outright crash.
    Anyway, if you run into RAM problems running this program,
   this is something you can try to decrease. Just don't get
   carried away.
*/
#define DEFAULT_MAX_FILESIZE 128000000

#define MAX_LPMEM_ALLOC 256000000
#define MAX_SPAIRS_ALLOC 12000000

/* One bit per entry in the (a,b) hash table. Higher means fewer collisions
   (and so, faster processing).
*/
#define AB_HASH_RAM 2*33554432

#define MAX_PBUF_RAM 32000000
#define DEFAULT_QCB_SIZE 62
#define DEFAULT_SEED 1
#define DEFAULT_NUM_FILES 4
#define TMP_FILE "tmpdata.000"
#define DEFAULT_DEPNAME "deps"
#define DEFAULT_COLNAME "cols"
#define DEFAULT_LPI_NAME "lpindex"
#define DEFAULT_PRELPREFIX "rels.bin"

#define USAGE " -fb <fname> -prel <fname> -newrel <fname> [-qs <qcb size>] [-v]\n"\
"-fb <fname>           : Factor base.\n"\
"-prel <file prefix>   : File name prefix for input/output of processed relations.\n"\
"-newrel <fname>       : File name for new, unprocessed relations.\n"\
"-nodfactor            : Don't try to factor the discriminant again.\n"\
"-cc <on | off | auto> : Cycle count: on, off, auto.\n"\
"-minff <int>          : Minimum number of FF's (prevent R-S wt. reduction and\n"\
"                        writing of the column files if there are fewer than this).\n"\
"-dump                 : Dump processed relations into siever-output formatted\n"\
"-s                    : Use short relations format (only a,b) with -dump and -prune\n"\
"-maxrelsinff <int>    : Max relation-set weight.\n"\
"-speedtest            : Do nothing but report a number representing the relative speed\n"\
"                        of this machine.\n"\
"-nolpcount            : Don't count large primes.\n"\
"-prune <float>        : EXPERIMENTAL! Remove the heaviest <float> fraction of processed\n"\
"                        relations (and dump them in siever-output format, just in case).\n"\
"                        ASCII files, then quit.\n"

#define START_MSG \
"\n"\
" __________________________________________________________ \n"\
"|        This is the procrels program for GGNFS.           |\n"\
"| Version: %-25s                       |\n"\
"| This program is copyright 2004, Chris Monico, and subject|\n"\
"| to the terms of the GNU General Public License version 2.|\n"\
"|__________________________________________________________|\n"



#define CC_OFF  0
#define CC_ON   1
#define CC_AUTO 2

/***** Globals *****/
int  discFact=1, cycleCount=CC_AUTO;
int  dump_short_mode=0; /* Sten: dump only a,b when -s parameter specified.  */
s32  initialFF=0, initialRelations=0, finalFF=0;
s32  totalLargePrimes=0;
s32  minFF;
long relsNumLP[8]={0,0,0,0,0,0,0,0};

/********************************************************************************/
long countLP(multi_file_t *prelF)
/********************************************************************************/
/* Find out how many distinct large primes are contained in the processed rels. */
/********************************************************************************/
{ s32 i, j, k, p, r, numRels=0, totalLP=0;
  s32 *lR=NULL, *lA=NULL, lRSize, lRMax, lASize, lAMax;
  u32 sF, lrpi, lapi;
  int  numLR, numLA;
  rel_list *RL;

  lRSize = lRMax=0;
  lASize = lAMax=0;
  for (i=0; i<prelF->numFiles; i++) {
    printf("Loading processed file %" PRId32 "/%d...", i+1, prelF->numFiles);
    fflush(stdout);
    RL = getRelList(prelF, i);
    printf("Done. Processing...\n");
    numRels += RL->numRels;
    for (j=0; j<RL->numRels; j++) {
      sF = RL->relData[RL->relIndex[j]];
      numLR = GETNUMLRP(sF);
      numLA = GETNUMLAP(sF);
      lrpi = 4 + 2*(GETNUMRFB(sF) + GETNUMAFB(sF) + GETNUMSPB(sF)) + 2;
      lapi = 4 + 2*(GETNUMRFB(sF) + GETNUMAFB(sF) + GETNUMSPB(sF)) + 2 + numLR;
      totalLP += numLR + numLA;
      for (k=0; k<numLR; k++) {
        p = RL->relData[RL->relIndex[j] + lrpi + k];
        /* So p is a large rational prime in this relation. */
        if (lRSize + 8 > lRMax) {
          lRMax += 65536;
          lR = realloc(lR, lRMax*sizeof(s32));
          if (lR == NULL) {
            printf("countLP() : Memory allocation error for lR!\n");
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
          lAMax += 65536;
          lA = realloc(lA, 2*lAMax*sizeof(s32));
          if (lA == NULL) {
            printf("countLP() : Memory allocation error for lA!\n");
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
    printf("Found %" PRId32 " distinct large rprimes and %" PRId32 " large aprimes so far.\n",lRSize, lASize); 
  }

  printf("There are %" PRId32 " large primes versus %" PRId32 " relations.\n", 
          lRSize+lASize, initialRelations);
  msgLog(NULL, "largePrimes: %" PRId32 " , relations: %" PRId32,
         lRSize+lASize, initialRelations);
  /* Now, check: should we do a cycle count? */
  free(lA); free(lR);
  return lRSize + lASize;
}

/******************************************************/
int set_prelF(multi_file_t *prelF, off_t maxFileSize, int takeAction)
/******************************************************/
/* Check to see how many files there are, and whether */
/* or not this needs to be increased. If so, increase.*/
/* Well, the increase will only be done if takeAction */
/* is nonzero - otherwise, do nothing but count them. */
/******************************************************/
{ int    i, cont, newFiles;
  struct stat fileInfo;
  char   fName[512], newName[512];
  off_t  maxSize=0;
  s32    where, k;
  rel_list   *RL;
  s32   bufSize, b, size, relsInFile;
  s64   a;
  s32  *newData[256], newDataIndex[256], newRels[256];
  int    fileno;
  char   prelname[256];
  FILE  *ofp;

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
  printf("Largest prel file size is %" PRIu64 " versus max allowed of %" PRIu64 ".\n", 
          (u64)maxSize, (u64)maxFileSize);
  if ((maxSize < maxFileSize)||(takeAction==0))
    return 0;

  /* We need to increase the number of files, which means
     moving the existing data around.
  */
  if (prelF->numFiles < 4)
    newFiles = 4;
  else newFiles = prelF->numFiles+4;

  printf("Increasing prelF->numFiles from %d to %d...\n",
          prelF->numFiles, newFiles);

  bufSize = MAX_PBUF_RAM/(MAX(1,newFiles)*sizeof(s32));
  for (i=0; i<newFiles; i++) {
    if (!(newData[i] = (s32 *)malloc(bufSize*sizeof(s32)))) {
      printf("Mem. allocation error for newData!\n");
      exit(-1);
    }
    newDataIndex[i]=0;
    newRels[i]=0;
  }

  for (i=0; i<prelF->numFiles; i++) {
    RL = getRelList(prelF, i);
    printf("Read %" PRId32 " relations from %s.%d\n", RL->numRels, prelF->prefix, i);

    for (k=0; k<RL->numRels; k++) {
      where = RL->relIndex[k];
      a = *((s64*) &(RL->relData[where+1]) );
      b = RL->relData[where+3];
      fileno = NFS_HASH(b, b, newFiles);
      size = RL->relIndex[k+1]-where;
      memcpy(&newData[fileno][newDataIndex[fileno]], &RL->relData[where], size*sizeof(s32));
      newDataIndex[fileno] += size;
      newRels[fileno] += 1;
      if (bufSize - newDataIndex[fileno]  < 500) {
        /* Append to file and clear. */
        /* It should be possible to do this with one fopen(), but I don't
           know about portability of doing it that way. */
        sprintf(prelname, "%s.%d_", prelF->prefix, fileno);
        if ((ofp = fopen(prelname, "rb"))) {
          rewind(ofp);
          fread(&relsInFile, sizeof(s32), 1, ofp);
          fclose(ofp);
        } else {
          relsInFile=0;
          /* And create the file. */
          ofp = fopen(prelname, "wb"); 
          fclose(ofp);
        }
        relsInFile += newRels[fileno];
        if ((ofp = fopen(prelname, "r+b"))) {
          rewind(ofp);
          fwrite(&relsInFile, sizeof(s32), 1, ofp);
          fseek(ofp, 0, SEEK_END);
          fclose(ofp);
        } 
        if ((ofp = fopen(prelname, "ab"))) {
          fwrite(newData[fileno], sizeof(s32), newDataIndex[fileno], ofp);
          fclose(ofp);
        }
        newDataIndex[fileno]=0;
        newRels[fileno]=0;
      }
    }
    /* Delete the old i-th file. */
    sprintf(prelname, "%s.%d", prelF->prefix, i);
    remove(prelname);
    clearRelList(RL);
    free(RL);
  }
  for (fileno=0; fileno<newFiles; fileno++) {
    sprintf(prelname, "%s.%d_", prelF->prefix, fileno);
    if ((ofp = fopen(prelname, "rb"))) {
      rewind(ofp);
      fread(&relsInFile, sizeof(s32), 1, ofp);
      fclose(ofp);
    } else {
      relsInFile=0;
      /* And create the file. */
      ofp = fopen(prelname, "wb"); 
      fclose(ofp);
    }
    relsInFile += newRels[fileno];
    if ((ofp = fopen(prelname, "r+b"))) {
      rewind(ofp);
      fwrite(&relsInFile, sizeof(s32), 1, ofp);
      fseek(ofp, 0, SEEK_END);
      fclose(ofp);
    } 
    if ((ofp = fopen(prelname, "ab"))) {
      fwrite(newData[fileno], sizeof(s32), newDataIndex[fileno], ofp);
      fclose(ofp);
    }
    newDataIndex[fileno]=0;
    newRels[fileno]=0;
  }
  prelF->numFiles = newFiles;
  for (i=0; i<prelF->numFiles; i++) {
    /* rename the files. */
    sprintf(fName, "%s.%d_", prelF->prefix, i);
    sprintf(newName, "%s.%d", prelF->prefix, i);
    rename(fName, newName);
  }
  for (i=0; i<newFiles; i++)
    free(newData[i]);

  return 0;
}

/* Notes: Avoiding duplicate (a,b) pairs is actually a much bigger
   problem than it may seem at first. After much toiling over possible
   ways to handle this little predicament, I have arrived at a reasonable
   (though ugly) method to achieve reasonable speed, with reasonable memory,
   while completely avoiding the possibility of accidentally throwing out
   non-duplicate relations. Here's how it goes:
   (1) Create two hash tables, abHash0 and abHash1, and two empty lists
       abList and abExtra.
   (2) Read in all processed relations. For each such relation, apply HASH0
       and see if it's value is already hashed out in abHash0. If not, hash
       it out. If the position was already hashed out, check the other
       table, abHash1, using the second hash function HASH1.
       In either case, add the relation to the list abList (it takes exactly
       two s32s to do this).
   (3) After all relations have been read in, sort the abList list (this
       takes about 15 seconds for 20M relations).
   (4) For each new relation to be processed, do the following:
       (i) compute HASH0, and see if the value is already hashed out of
           the table abHash0. If not, we know it's not a duplicate. Hash it out
           of both tables, add it to the list abExtra, and stop.
       (ii) If the HASH0 value was already hashed out of the table, try HASH1.
           If this value is not in the table, we know the relation is not
           a duplicate, so hash it out of abHash1, add it to the list abExtra,
           and stop.
       (iii) If the hash value collided in both tables abHash0 and abHash1,
           do a binary search on the original list, abList, to see if it's there.
           If so, it's a duplicate so throw it out.
       (iv) If the relation was not in the original list, check the list where
           the new relations are living, abExtra. If it's not there either,
           it's not a duplicate so add it to abExtra and stop. (In fact, we also
           maintain a hash table of the newly processed relations to prevent
           excess sorting and searching; actual sorting and searching of the 
           newly processed relations is a very last resort (hehe - no pun intended)).
  Note: A semi-obvious alternative at first is to use some binary trees to
  help. However, we really do need to minimize the amount of memory used in
  this process, and a btree just uses too much (it would double the amount
  of RAM needed). Furthermore, it would still need to be used in conjunction
  with some sort of hashing technique to achieve the speed we've got here.
  Notice that if everything is set up right, we will only need to do binary
  searches on rare occassions. For most relations, we will see from the hash
  tables that they are unique.
    This brings up another point: Why two hash tables? We cannot afford to
  use the astronomical amount of RAM that would be needed to make collisions
  in a single table unlikely. The method we've used here insures that the
  only relations which will make it to the second hash table are those that
  would collide with another earlier relation in the first table. This makes
  collisions in the second table far less likely than collisions would have
  been in a single table of twice the size (i.e., the birthday problem).
  It takes a bit of thought to see that this actually works. For example, we
  rely heavily on the assumption that any relations which have already been
  processed are unique; the whole method would die a horrible death if that
  assumption were false.
*/
/* These data are shared by the next three functions. */
u32 *abHash0=NULL, *abHash1=NULL, *abExtraHash=NULL;
s32 abHashSize=0, abHashWords=0;
s32 *abList=NULL, abListSize=0, abListMax=0;
s32 *abExtra=NULL, abExtraSize=0, abExtraMax=0;
int   abExtraSorted=1;
#define REL_HASH0(_a, _b, _s) (((u32)((_a)*314159265)+((_b)*577215664))%(_s))
#define REL_HASH1(_a, _b, _s) (((u32)((_a)*592653141)+((_b)*271828182))%(_s))
#define REL_HASH2(_a, _b, _s) (((u32)((_a)*531415926)+((_b)*828182271))%(_s))

/*****************************************************************/
s32 makeABLookup(multi_file_t *prelF)
/*****************************************************************/
/* Lookup all processed relations and build some structures for  */
/* checkAB() to be able to lookup an (a,b) pair to see if it's a */
/* duplicate.                                                    */
/*****************************************************************/
{ int  i;
  char prelname[256];
  s32  r, loc, b;
  s64  a;
  u32  s;
  u32  h0, h1;
  rel_list RL;

  abHashWords = AB_HASH_RAM/(3*4);
  abHashSize = 32*abHashWords;
  if (!(abHash0 = (s32 *)malloc(abHashWords*sizeof(u32)))) {
    printf("Memory allocation error for abHash0!\n");
    exit(-1);
  }
  if (!(abHash1 = (s32 *)malloc(abHashWords*sizeof(u32)))) {
    printf("Memory allocation error for abHash1!\n");
    exit(-1);
  }
  if (!(abExtraHash = (s32 *)malloc(abHashWords*sizeof(u32)))) {
    printf("Memory allocation error for abExtraHash!\n");
    exit(-1);
  }
  abListMax = 32768;
  abList = malloc(2*abListMax*sizeof(s32));
  abListSize=0;
  

  allocateRelList(prelF, &RL);

  memset(abHash0, 0x00, abHashWords*sizeof(s32));
  memset(abHash1, 0x00, abHashWords*sizeof(s32));
  memset(abExtraHash, 0x00, abHashWords*sizeof(s32));
  printf("Building (a,b) hash table..."); fflush(stdout);
  for (i=0; i<prelF->numFiles; i++) {
    printf("%d..", i); fflush(stdout);
    sprintf(prelname, "%s.%d", prelF->prefix, i);
    if (readRelList(&RL, prelname)) {
      printf("makeABList() Failed to open %s for read!\n", prelname);
      break;
    }
    /* Now go through the list and hash out each (a,b) pair. */
    for (r=0; r<RL.numRels; r++) {
      loc = RL.relIndex[r];
      s = RL.relData[loc];
      relsNumLP[GETNUMLRP(s)+GETNUMLAP(s)] += 1;
      a = *( (s64*)&(RL.relData[loc+1]) );
      b = RL.relData[loc+3];
      h0 = REL_HASH0(a, b, abHashSize);
      h1 = REL_HASH1(a, b, abHashSize);
      abHash0[h0/32] |= BIT(h0&0x0000001F);
      abHash1[h1/32] |= BIT(h1&0x0000001F);
      if (abListSize +1 > abListMax) {
        abListMax += 1048576;
        abList = realloc(abList, 2*abListMax*sizeof(s32));
        if (abList==NULL) {
          printf("Fatal memory allocation error for abList!\n");
          exit(-1);
        }
      }
      abList[2*abListSize]=a;
      abList[2*abListSize+1]=b;
      abListSize++;
    }
  }
  printf("\n");
  clearRelList(&RL);
  /* Sort the list. */
  printf("makeABLookup() : Sorting abList..."); fflush(stdout);
  qsort(abList, abListSize, 2*sizeof(s32), cmp2S32s);
  printf("Done.\n");
  return abListSize;
}

/*****************************************************************/
void clearABLookup()
/*****************************************************************/
/* Free up memory used for the (a,b) lookup tables.              */
/*****************************************************************/
{
  if (abHash0) free(abHash0);
  if (abHash1) free(abHash1);
  if (abExtraHash) free(abExtraHash);
  if (abList) free(abList);
  if (abExtra) free(abExtra);

  abHash0 = abHash1 = abList = abExtra = NULL;
  abHashSize = abHashWords = abListSize = abListMax = 0;
  abExtraSize = abExtraMax = 0;
}

s32 sortOps=0;
/*****************************************************************/
int checkAB(s64 a, s32 b)
/*****************************************************************/
/* Is this an already-processed (a,b) pair? If not, we add it to */
/* the list of processed pairs and return 0. Otherwise, return   */
/* some nonzero value.                                           */
/* This can be done better: I will make some improvements so that
   it doesn't slow down as much when nearing the end of a large
   input file.
*/
{ u32 h0, h1, h2;
  s32 key[2], *loc;

  h0 = REL_HASH0(a, b, abHashSize);
  h1 = REL_HASH1(a, b, abHashSize);
  if ((abHash0[h0/32]&BIT(h0&0x0000001F))==0) {
    abHash0[h0/32] |= BIT(h0&0x0000001F);
    abHash1[h1/32] |= BIT(h1&0x0000001F);
    if (abExtraSize + 1 >= abExtraMax) {
      abExtraMax += 32768;
      abExtra = realloc(abExtra, abExtraMax*2*sizeof(s32));
      if (abExtra == NULL) {
        printf("Memory (re-)allocation error for abExtra!\n");
        exit(-1);
      }
    }
    abExtra[2*abExtraSize] = a;
    abExtra[2*abExtraSize+1] = b;
    abExtraSorted=0;
    abExtraSize++;
    h2 = REL_HASH2(a, b, abHashSize);
    abExtraHash[h2/32] |= BIT(h2&0x0000001F);
    return 0;
  }
  if ((abHash1[h1/32]&BIT(h1&0x0000001F))==0) {
    abHash1[h1/32] |= BIT(h1&0x0000001F);
    if (abExtraSize + 1 >= abExtraMax) {
      abExtraMax += 32768;
      abExtra = realloc(abExtra, abExtraMax*2*sizeof(s32));
      if (abExtra == NULL) {
        printf("Memory (re-)allocation error for abExtra!\n");
        exit(-1);
      }
    }
    abExtra[2*abExtraSize] = a;
    abExtra[2*abExtraSize+1] = b;
    abExtraSorted=0;
    abExtraSize++;
    h2 = REL_HASH2(a, b, abHashSize);
    abExtraHash[h2/32] |= BIT(h2&0x0000001F);
    return 0;
  }
  /* Before we sort and search, make sure the size isn't getting too large: */
  if (abExtraSize > MAX_AB_EXTRA_ENTRIES) {
    /* Clear it, and add them to the big list, to reduce collisions. */
    abListMax = abListSize + abExtraSize;
    abList = realloc(abList, 2*abListMax*sizeof(s32));
    if (abList==NULL) {
      printf("Fatal memory allocation error for abList!\n");
      exit(-1);
    }
    memcpy(abList+2*abListSize, abExtra, abExtraSize*2*sizeof(s32));
    abListSize += abExtraSize;
    abExtraSize=0;
    qsort(abList, abListSize, 2*sizeof(s32), cmp2S32s);
    /* Reset the abExtraHash to empty. */
    memset(abExtraHash, 0x00, abHashWords*sizeof(s32));
  }

  /* Ok - we need to check both lists: */
  key[0]=a; key[1]=b;
  if (abListSize) {
    loc = bsearch(key, abList, abListSize, 2*sizeof(s32), cmp2S32s);
    if (loc) return -1;
  }

  h2 = REL_HASH2(a, b, abHashSize);
  if ((abExtraHash[h2/32]&BIT(h2&0x0000001F))==0) {
    abExtraHash[h2/32] |= BIT(h2&0x0000001F);
    if (abExtraSize + 1 >= abExtraMax) {
      abExtraMax += 32768;
      abExtra = realloc(abExtra, abExtraMax*2*sizeof(s32));
      if (abList == NULL) {
        printf("Memory (re-)allocation error for abList!\n");
        exit(-1);
      }
    }
    abExtra[2*abExtraSize] = a;
    abExtra[2*abExtraSize+1] = b;
    abExtraSorted=0;
    abExtraSize++;
    return 0;
  }


  /* Ok - we are having some seriously bad luck, or there actually were two duplicates
     in the new relation file. Let's check: (we might think about doing this with a
     binary tree instead).
  */
  qsort(abExtra, abExtraSize, 2*sizeof(s32), cmp2S32s);
  sortOps++;
  abExtraSorted=1;
  loc = bsearch(key, abExtra, abExtraSize, 2*sizeof(s32), cmp2S32s);
  if (loc) return -1;
  if (abExtraSize + 1 >= abExtraMax) {
    abExtraMax += 32768;
    abExtra = realloc(abExtra, abExtraMax*2*sizeof(s32));
    if (abList == NULL) {
      printf("Memory (re-)allocation error for abList!\n");
      exit(-1);
    }
  }
  abExtra[2*abExtraSize] = a;
  abExtra[2*abExtraSize+1] = b;
  abExtraSorted=0;
  return 0;
}


/*****************************************************************/
s32 getABHash(u32 *hash0, u32 *hash1, s32 hashSize, multi_file_t *prelF)
/*****************************************************************/
/* 'hash' is a hash table, one bit per entry. We will initialize */
/* it to zero and then turn on each bit corresponding to an      */
/* existing relation.                                            */
/*****************************************************************/
{ int i;
  char prelname[256];
  s32 r, loc, a, b, total=0;
  u32  h;
  rel_list RL;
  

  allocateRelList(prelF, &RL);

  memset(hash0, 0x00, hashSize/32);
  memset(hash1, 0x00, hashSize/32);
  printf("Building (a,b) hash table..."); fflush(stdout);
  for (i=0; i<prelF->numFiles; i++) {
    printf("%d..", i); fflush(stdout);
    sprintf(prelname, "%s.%d", prelF->prefix, i);
    if (readRelList(&RL, prelname))
      break;
    /* Now go through the list and hash out each (a,b) pair from both tables. */
    for (r=0; r<RL.numRels; r++) {
      loc = RL.relIndex[r];
      a = RL.relData[loc+1];
      b = RL.relData[loc+2];
      h = REL_HASH0(a, b, hashSize);
      hash0[h/32] |= BIT(h&0x0000001F);
      h = REL_HASH1(a, b, hashSize);
      hash1[h/32] |= BIT(h&0x0000001F);
      total++;
    }
  }
  printf("\n");
  clearRelList(&RL);
  return total;
}


/***************************************************************/
s32 addNewRelations5(multi_file_t *prelF, char *fName,  nf_t *N)
/***************************************************************/
/* Read new relations from fName and add them to the processed */
/* relation files. This function is very different from its    */
/* predecessor, addNewRelations4() in that it does everything  */
/* (and does it more efficiently).                             */
/* NOT DONE YET! */
/***************************************************************/
{ s32        numNew=0, numRead=0, total=0;
  relation_t R;
  FILE      *fp, *ofp;
  int        factRes, i, j;
  char       thisLine[512];
  double     startTime, now;
  s32        nextReportNumRead = 10000, collisions=0;
  s32        fSize, fBlockSize, fRemainSize, fTotalRead=0;
  u32        s;
  unsigned char *fData, *fPos, *fLimit = NULL, *fEol = NULL;
  unsigned char *fWarningTrack=NULL;
  static char xdigit[256] = {
    -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0x00+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0x10+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0x20+ */
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, -1, -1, -1, -1, -1, /* 0x30+ */
    -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0x40+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0x50+ */
    -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0x60+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0x70+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0x80+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0x90+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0xa0+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0xb0+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0xc0+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0xd0+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0xe0+ */
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, /* 0xf0+ */
  };
  nfs_fb_t *FB = N->FB;
  s32 maxRFB = FB->rfb[2 * (FB->rfb_size - 1)];
  s32 maxAFB = FB->afb[2 * (FB->afb_size - 1)];
  s32 bufSize;
  s32 *newData[256], newDataIndex[256], newRels[256], relsInFile;
  int  fileno;
  char prelname[256];

  /* If there is nothing to add, we still must return the
     total number of relations.
  */
  if (!(fp = fopen(fName, "r"))) {
    for (i=0; i<prelF->numFiles; i++) {
      sprintf(prelname, "%s.%d", prelF->prefix, i);
      if ((ofp = fopen(prelname, "rb"))) {
        rewind(ofp);
        fread(&relsInFile, sizeof(s32), 1, ofp);
        fclose(ofp);
      }
      total += relsInFile;
    }
    return total;
  }

  total = makeABLookup(prelF);
  printf("Before processing new relations, there are %" PRId32 " total.\n", total);

  /* Set up to prepare for the new data: */
  bufSize = MAX_PBUF_RAM/(MAX(1,prelF->numFiles)*sizeof(s32));
  for (i=0; i<prelF->numFiles; i++) {
    if (!(newData[i] = (s32 *)malloc(bufSize*sizeof(s32)))) {
      printf("Mem. allocation error for newData!\n");
      exit(-1);
    }
    newDataIndex[i]=0;
    newRels[i]=0;
  }
  if (!(fp = fopen(fName, "r"))) {
    fprintf(stderr, "Error opening %s for read.\n", fName);
    return 0;
  }
  fseek(fp, 0, SEEK_END);
  fSize = ftell(fp);
  if (fSize == 0) {
    fclose(fp);
    return 0;
  }
  fseek(fp, 0, SEEK_SET);
  fBlockSize = MIN(fSize, MAX_SPAIRS_ALLOC);
  if ((fData = (unsigned char *)malloc(fBlockSize + 1)) != NULL) {
    fTotalRead = fread(fData, 1, fBlockSize, fp);
    fLimit = fWarningTrack = fData + fBlockSize;
    if (fBlockSize >= MAX_SPAIRS_ALLOC)
      fWarningTrack -= 512; /* Where we stop and read the next block from file. */
    *fLimit = '\0';
    fEol = fData - 1;
  }
  startTime = sTime();
  for (;;) {
    s32 r, k;
    s64 p;
    int c, m;
    int short_form = 0; /* 0 - we have only a,b.
                           1 - we have a,b and all large factors. */

    if (fData != NULL) {
      fPos = fEol + 1;
      if (fPos >= fLimit) {
        break;
      } else if ((fPos > fWarningTrack)&&(fTotalRead < fSize)) {
        /* Read the next block of data from file. */
        fRemainSize = fLimit-fPos;
        memmove(fData, fPos, fRemainSize*sizeof(char));
        fBlockSize = fRemainSize + fread(fData+fRemainSize, 1, fBlockSize-fRemainSize,fp);
        fTotalRead += fBlockSize;
        fPos = fData;
//        fEol = fData - 1;
        fLimit = fData + fBlockSize;
        *fLimit = '\0';
      }
    } else {
      if (feof(fp)) {
        break;
      }
      thisLine[0] = '\0';
      fgets(thisLine, 512, fp);
      fPos = thisLine;
    }
    for (fEol = fPos; *fEol && *fEol != '\n'; fEol++) { /* search end-of-line */
      ;
    }
    if (*fPos == '#') { /* comment line */
      continue;
    }
    *fEol = '\0';

    /* expand and optimize parseOutputLine() */
    /* %ld */
    fPos += (m = *fPos == '-');
    if ((unsigned int)(p = *fPos - '0') >= 10) 
    {
      continue;
    }
    
    fPos++;
    
    while ((unsigned int)(c = *fPos - '0') < 10) 
    {
      p = (((p << 2) + p) << 1) + c; /* p = p * 10 + c */
      fPos++;
    }
    
    R.a = m ? -p : p;
    
    /* , */
    fPos++;
    
    /* %ld */
    fPos += (m = *fPos == '-');
    
    if ((unsigned int)(p = *fPos - '0') >= 10) 
    {
      continue;
    }
    
    fPos++;
    
    while ((unsigned int)(c = *fPos - '0') < 10)
    {
      p = (((p << 2) + p) << 1) + c; /* p = p * 10 + c */
      fPos++;
    }

    R.b = m ? -p : p;

    /* : */
    while (fPos < fEol && (*fPos != ':')) 
    {
      fPos++;
    }
    
    short_form = (fPos == fEol);
    fPos++;

    if (short_form == 0)
    {
        /* This is long form of the relation.  */
        /* %lx,%lx,... */
        for (j=0; j<FB->maxLP; j++)
        {
          R.p[j] = 1;
        }

        R.rFSize = 0;
        m = 0;
        
        while (fPos < fEol && *fPos != ':') 
        {
          if ((p = xdigit[*fPos++]) < 0) 
          {
            continue;
          }
          
          while ((c = xdigit[*fPos]) >= 0) 
          {
            p = (p << 4) + c;
            fPos++;
          }
          
          k = lookupRFB(p, FB);
          
          if (k >= 0)
          {
            R.rFactors[R.rFSize++] = k;
          }
          else if ((p > maxRFB) && (m < FB->maxLP)) 
          {
            R.p[m++] = p;
          }
        }

        /* : */
        fPos += *fPos == ':';
        
        /* %lx,%lx,... */
        
        for (j=0; j<FB->maxLPA; j++)
        {
          R.a_p[j] = R.a_r[j] = 1;
        }

        R.aFSize = 0;
        m = 0;
        
        while (fPos < fEol && *fPos != ':') 
        {
          if ((p = xdigit[*fPos++]) < 0) 
          {
            continue;
          }
          
          while ((c = xdigit[*fPos]) >= 0) 
          {
            p = (p << 4) + c;
            fPos++;
          }
          
          if (R.b % p) 
          { 
            /* p is always non-zero? */
            r = mulmod32(p + (R.a % p), inverseModP(R.b, p), p);
            k = lookupAFB(p, r, FB);
          
            if (k >= 0) 
            {
              R.aFactors[R.aFSize++] = k;
            } 
            else 
            if ((p > maxAFB) && (m < FB->maxLPA)) 
            {
              R.a_p[m] = p; 
              R.a_r[m++] = r;
            }
          }
        }
    } /* if (*fPos == ':')  */

    numRead++;

//    printf("Read (%" PRId64 ", %ld) from file\n", R.a, R.b );

    if (checkAB(R.a, R.b)==0) {
      fileno = NFS_HASH(R.b, R.b, prelF->numFiles);
      /* Sten: we smartly choose here if this is short format or long and thus if
               we should try to factor relation completely or only partly. */
      factRes = (short_form ? factRel(&R, N) : completePartialRelFact(&R, N, CLIENT_SKIP_R_PRIMES, CLIENT_SKIP_A_PRIMES));
      if (factRes == 0) {
        k = newDataIndex[fileno];
        newDataIndex[fileno] += relConvertToData(&newData[fileno][newDataIndex[fileno]], &R);
        s = newData[fileno][k];
        relsNumLP[GETNUMLRP(s)+GETNUMLAP(s)] += 1;

        newRels[fileno] += 1;
        if (bufSize - newDataIndex[fileno]  < 500) {
          /* Append to file and clear. */
          /* It should be possible to do this with one fopen(), but I don't
             know about portability of doing it that way. */
          sprintf(prelname, "%s.%d", prelF->prefix, fileno);
          if ((ofp = fopen(prelname, "rb"))) {
            rewind(ofp);
            fread(&relsInFile, sizeof(s32), 1, ofp);
            fclose(ofp);
          } else {
            relsInFile=0;
            /* And create the file. */
            ofp = fopen(prelname, "wb"); 
            fclose(ofp);
          }
          relsInFile += newRels[fileno];
          if ((ofp = fopen(prelname, "r+b"))) {
            rewind(ofp);
            fwrite(&relsInFile, sizeof(s32), 1, ofp);
            fseek(ofp, 0, SEEK_END);
            fclose(ofp);
          } 
          if ((ofp = fopen(prelname, "ab"))) {
            fwrite(newData[fileno], sizeof(s32), newDataIndex[fileno], ofp);
            fclose(ofp);
          }
          newDataIndex[fileno]=0;
          newRels[fileno]=0;
        }
         
        numNew++;
      } else {
#ifdef _DEBUG
        printf("Relation (%" PRId64 ", %ld) bad : return value %d.\n", R.a, R.b, factRes);
#endif
        ;
      }
    } else {
      collisions++;
    }
    if (numRead >= nextReportNumRead) {
      nextReportNumRead += 10000;
      now = sTime();
      printTmp("Status: processed %ld relations from %s... (at %1.2lf rels/sec)",
                 numRead, fName,
                 /* considering now == startTime */
                 now != startTime ? (double)numRead / (now - startTime) : 0.0);
    }
  }
  fclose(fp);

  /* Dump any remaining relations to their files. */
  for (i=0; i<prelF->numFiles; i++) {
    if (newRels[i] > 0) {
      sprintf(prelname, "%s.%d", prelF->prefix, i);
      if ((ofp = fopen(prelname, "rb"))) {
        rewind(ofp);
        fread(&relsInFile, sizeof(s32), 1, ofp);
        fclose(ofp);
      } else { 
        relsInFile=0;  
        /* Create the file. */
        ofp = fopen(prelname, "wb"); 
        fclose(ofp);
      }
      relsInFile += newRels[i];
      if ((ofp = fopen(prelname, "r+b"))) {
        rewind(ofp);
        fwrite(&relsInFile, sizeof(s32), 1, ofp);
        fseek(ofp, 0, SEEK_END);
        fclose(ofp);
      } 
      if ((ofp = fopen(prelname, "ab"))) {
        fwrite(newData[i], sizeof(s32), newDataIndex[i], ofp);
        fclose(ofp);
      }
      newDataIndex[i]=0;
      newRels[i]=0;
    }
  }
  printf("\n");

  if (fData != NULL) 
    free(fData);
  for (i=0; i<prelF->numFiles; i++) 
    free(newData[i]);
  clearABLookup();
  printf("   abExtra was sorted %" PRId32 " times.\n", sortOps);
  msgLog("", "There were %" PRId32 "/%" PRId32 " duplicates.",
         collisions, numRead);
  total += numNew;
  return total;
}


#define MAX_DUMP_PER_FILE 250000
/****************************************************/
s32 dumpPairs(char *fName, multi_file_t *prelF, nfs_fb_t *FB)
/****************************************************/
/* Dump all the relations into a text file with the */
/* siever output format (useful for debugging, or   */
/* changing factor base sizes mid-factorization).   */
/****************************************************/
{ int  i, fileNum=0;
  char prelName[512], outStr[1024], outName[256];
  FILE *ofp;
  relation_t R;
  s32 j, numRels, total=0;
  rel_list *RL;


  sprintf(outName, "%s.%d", fName, fileNum);
  if (!(ofp = fopen(outName, "wb"))) {
    printf("Error opening %s for write!\n", fName);
    exit(-1);
  }
  for (i=0; i<prelF->numFiles; i++) {
    sprintf(prelName, "%s.%d", prelF->prefix, i);
    RL = getRelList(prelF, i);
    numRels = RL->numRels;
    printf("Dumping %" PRId32 " relations from %s...\n", numRels, prelName);
    for (j=0; j<numRels; j++) {
      dataConvertToRel(&R, &RL->relData[RL->relIndex[j]]);
      makeOutputLine(outStr, &R, FB, dump_short_mode);
      fprintf(ofp, "%s\n", outStr);
      if ((++total % MAX_DUMP_PER_FILE)==0) {
        fclose(ofp);
        fileNum++;
        sprintf(outName, "%s.%d", fName, fileNum);
        if (!(ofp = fopen(outName, "wb"))) {
          printf("Error opening %s for write!\n", fName);
          exit(-1);
        }
      }
    }
    clearRelList(RL);
    free(RL);
  }
  fclose(ofp);
  printf("Dumped %" PRId32 " relations to %s.\n", total, fName);
  return total;
}

/****************************************************/
int fsingleVerbose(nf_t *N, char *line)
/****************************************************/
{ relation_t R;
  int res, i, n;

  sscanf(line, "%" SCNd64 ",%" SCNd32, &R.a, &R.b);
  printf("Attempting to factor relation (%" PRId64 ", %" PRId32 ")\n", R.a, R.b);
  if (R.b <= 0) {
    printf("Error: 'b' should be positive!\n");
    return -1;
  }
  res = factRel(&R, N);
  printf("res=%d\n", res);
  printf("RFB:\n");
  for (i=0; i<R.rFSize; i++) {
    printf("(%" PRId32 ")^%d ", N->FB->rfb[2*R.rFactors[i]], R.rExps[i]);
    if (i%6==5) printf("\n");
  }
  printf("Large:\n");
  for (i=0; i<MAX_LARGE_RAT_PRIMES; i++) {
    if (R.p[i]>1)
      printf("%" PRId32 " ", R.p[i]);
  }
  
  printf("\nAFB:\n");
  for (i=0; i<R.aFSize; i++) {
    printf("(%" PRId32 ",%" PRId32 ")^%d ", N->FB->afb[2*R.aFactors[i]], 
           N->FB->afb[2*R.aFactors[i]+1], R.aExps[i]);
    if (i%6==5) printf("\n");
  }
  printf("Large:\n");
  for (i=0; i<MAX_LARGE_ALG_PRIMES; i++) {
    if (R.a_p[i]>1)
      printf("(%" PRId32 ", %" PRId32 ") ", R.a_p[i], R.a_r[i]);
  }
  printf("Special primes:\n");
  for (i=0; i<R.spSize; i++) {
    n = R.spFactors[i];
    printf("  exponent: %d, ideal %d:", R.spExps[i], n);
    mpz_out_str(stdout, 10, N->sPrimes[n].p);
    printf(", ");
    mpz_poly_print(stdout, "", N->sPrimes[n].alpha);
    printf("\n");
  }
  printf("QCB entries: %8.8" PRIx32 " %8.8" PRIx32 "\n", R.qcbBits[0], R.qcbBits[1]);
  return 0;
}
  


/****************************************************/
int main(int argC, char *args[])
/****************************************************/
{ char       fbName[64], prelName[40], newRelName[64], depName[64], colName[64];
  char       tmpStr[1024], line[128];
  int        i, qcbSize = DEFAULT_QCB_SIZE, seed=DEFAULT_SEED, retVal=0, dump=0;
  int        fr=0, maxRelsInFF=MAX_RELS_IN_FF, doCountLP=1;
  double     startTime, rStart, rStop, pruneFrac=0.0;
  off_t      oldSize, newSize, maxSize;
  s32        totalRels, numNewRels;
  struct stat fileInfo;
  nf_t       N;
  mpz_fact_t D;
  multi_file_t prelF, lpF;
  FILE            *fp;

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
    } else if (strcmp(args[i], "-newrel")==0) {
      if ((++i) < argC) 
        strncpy(newRelName, args[i], 64);
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
    } else if (strcmp(args[i], "-nodfactor")==0) {
      discFact=0;
    } else if (strcmp(args[i], "-cc")==0) {
      if ((++i) < argC) {
        if (strcmp(args[i], "off")==0)
          cycleCount=CC_OFF;
        else if (strcmp(args[i], "on")==0) 
          cycleCount=CC_ON;
        else if (strcmp(args[i], "auto")==0)
          cycleCount=CC_AUTO;
      }
    } else if (strcmp(args[i], "-v")==0) {
      verbose++;
    } else if (strcmp(args[i], "-dump")==0) {
        dump=1;
    } else if (strcmp(args[i], "-s")==0) {
        dump_short_mode=1;
    } else if (strcmp(args[i], "-fr")==0) {
      fr=1;
      if ((++i)<argC) 
        strncpy(line, args[i],128);
    } else if (strcmp(args[i], "-maxrelsinff")==0) {
      if ((++i) < argC) 
        maxRelsInFF = atoi(args[i]);
    } else if (strcmp(args[i], "-prune")==0) {
      if ((++i) < argC) 
        pruneFrac = atof(args[i]);
    } else if (strcmp(args[i], "-nolpcount")==0) {
      doCountLP=0;
    } else if (strcmp(args[i], "-speedtest")==0) {
      u32 a,b[1024],c=rand();
      double start=sTime(), now;
      for (a=2; a<364500000; a++) {
        b[a&400] = (b[c&400]*a)+c; c = (c+a)%101;
        b[(a+c)&400] = b[c&400]*c; c *= (b[a&400]+a)%100001;
        c = prand();
      }
      now=sTime();
      printf("b[a%%1024] = %8.8" PRIx32 "\n", b[a%1024]); /* So the compiler doesn't remove the loop above! */
      printf("timeunit: %1.3lf\n",10.0/(now-start));
      exit(0);
    }
  }
  maxRelsInFF=MIN(MAX_RELS_IN_FF,maxRelsInFF);
  srand(seed);
  startTime = sTime();

  if ((fbName[0]==0) || (prelF.prefix[0]==0)) {
    printf("USAGE: %s %s\n", args[0], USAGE);
    exit(0);
  }
  msgLog("", "GGNFS-%s : procrels", GGNFS_VERSION);
  sprintf(prelName, "%s.%d", prelF.prefix, 0);
  if (stat(prelName, &fileInfo)==0) {
    printf("It appears this is not the first run. Setting discFact=-1.\n");
    discFact=0;
  }
  initNF(&N);
  N.FB = (nfs_fb_t *)malloc(sizeof(nfs_fb_t));
  initFB(N.FB);
  mpz_fact_init(&D);
  if (loadFB(fbName, N.FB)) {
    printf("Could not load FB from %s!\n", fbName);
    exit(-1);
  }

  if ((pruneFrac < 0.0) || (pruneFrac > 0.9)) {
    printf("unreasonable value of pruneFrac supplied! Should be in (0,1)!\n");
    return 0;
  }
  if (pruneFrac > 0.0000000001) {
    set_prelF(&prelF, DEFAULT_MAX_FILESIZE, 0);
    pruneRelLists(&prelF, "spairs.dump", pruneFrac, N.FB, dump_short_mode);
    return 0;
  }
  if (minFF < N.FB->rfb_size + N.FB->afb_size + 64 + 32)
    minFF = N.FB->rfb_size + N.FB->afb_size + 64 + 32;
  if (verbose)
    printf("Getting QCB of size %d...\n", qcbSize);
  generateQCB(N.FB, qcbSize); 
  if (verbose)
    printf("Determining needed info about the number field...\n");
  mpz_poly_cp(N.f, N.FB->f);
  get_g(N.T, N.FB);
  mpz_poly_discrim(D.N, N.T);
  mpz_fact_factorEasy(&D, D.N, discFact);
  printf("Monic polynomial: T="); mpz_poly_print(stdout, "", N.T);
  /* If it the discrim. didn't completely factor above, there's
     no sense trying again! */
  getIntegralBasis(&N, &D, 0); 

  printf("Obtained integral basis:\nW = \n");
  mpz_mat_print(stdout, N.W);
  printf("denominator = "); mpz_out_str(stdout, 10, N.W_d); printf("\n");

  /* This is for debugging: factor a single relation from the command line
     and exit.
  */
  if (fr && line[0]) {
    fsingleVerbose(&N, line);
    exit(0);
  }

  set_prelF(&prelF, DEFAULT_MAX_FILESIZE, 1);

  maxSize = 0;
  for (i=0; i<prelF.numFiles; i++) {
    sprintf(prelName, "%s.%d", prelF.prefix, i);
    if (stat(prelName, &fileInfo)) {
      printf("Warning: Could not stat processed file %s. Is this the first run?.\n", prelName);
      oldSize = 0;
    } else {
      oldSize = fileInfo.st_size;
      printf("Existing file %s %1.2lfMB.\n", prelName, (double)oldSize/(1024.0*1024.0));
    }
    maxSize = MAX(maxSize, oldSize);
  }
  if (dump) {  
    dumpPairs("spairs.dump", &prelF, N.FB);
    exit(0);
  }

  numNewRels=0;
  if (stat(newRelName, &fileInfo)) {
    printf("Warning: Could not stat new file %s. No new relations will be processed.\n", newRelName);
    newSize = 0;
  } else {
    newSize = fileInfo.st_size;
    printf("     New file is %1.5lfMB.\n", (double)newSize/(1024.0*1024.0));
    if ((fp = fopen(newRelName, "r"))) {
      while (!(feof(fp))) {
        fgets(tmpStr, 1023, fp);
        numNewRels++;
      }
      fclose(fp);
    }
    /* Each rel ends with <nl>, so that when the last rel is read from
     * the file it STILL does not indicate EOF.  The "while" loop therefore
     * continues.
     */
    printf("     New file has %" PRId32 " relations.\n", numNewRels-1);
  }

  totalRels = 0;
  rStart = sTime();
  totalRels = addNewRelations5(&prelF, newRelName, &N);
  rStop = sTime();
  msgLog("", "RelProcTime: %1.1lf", rStop-rStart);


  printf("--------------------------------------\n");
  printf("There are now a total of %" PRId32 " unique relations in %d files.\n",
          totalRels, prelF.numFiles);
  initialRelations = totalRels;
  if (totalRels <= 0) {
    printf("Warning: no valid relations were processed. Either you have\n");
    printf("chosen very bad parameters for this factorization, or your usage\n");
    printf("was wrong.\n");
    exit(-1);
  }
  printf("# l.p. | Relations\n");
  printf("--------------------------\n");
  for (i=0; i<8; i++) {
    if (relsNumLP[i]>0) {
      printf("     %d | %ld\n", i, relsNumLP[i]);
    }
  }
  if (doCountLP)
    countLP(&prelF);
  return retVal;
}  


