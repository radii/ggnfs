/**************************************************************/
/* combparts.c                                                */
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
/* Routines for combining partials to form full relations.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#ifndef _MSC_VER        
#include <sys/time.h>
#endif

#include "ggnfs.h"

#define RS_WT_REDUCTION

#define LARGE_BUFFER 100000
#define FIELDSIZE(_L, _i) (_L->index[(_i)+1] - _L->index[_i])

#define BITCOUNT(_c) (bitcount[_c&0x0000FFFF] + bitcount[(_c&0xFFFF0000)>>16])

/* Globals: */
unsigned char *counter=NULL;
s32 counterSize;
int  *bitcount=NULL;

s32 mkUniqueS32s(s32 *x, s32 size);

/*******************************************/
void init_bitcount()
/*******************************************/
{ int i, j, c;

  if (bitcount == NULL)
    bitcount = (int *)malloc(65536*sizeof(int));
  for (i=0; i<65536; i++) {
    for (c=0, j=0; j<16; j++) 
      if (i&BIT(j)) c++;
    bitcount[i] = c;
  }
}


/*************************************************************/
void sortByNumLP(s32 *fieldsBySize, llist_t *P, char *str)
/*************************************************************/
{ s32 i, size;
  s32 numBySize[8];

  ll_getsortOnFieldSize(fieldsBySize, P);
  if (str != NULL) {
    printf(str);
    for (i=0; i<8; i++) /* Overkill - 4 or 6 would do. */
      numBySize[i]=0;
    i=0;
    while (i<P->numFields) {
      if (fieldsBySize[i] > P->numFields) 
        printf("?!?WTF?!? fieldsBySize[%" PRId32 "]=%" PRId32 " vs. P->numFields=%" PRId32 ".\n",
               i,fieldsBySize[i], P->numFields);
      size = P->index[fieldsBySize[i]+1] - P->index[fieldsBySize[i]];
      if (size > 6) 
        printf("?!Error!? Relation has more than 6 large primes!\n");
      else if (size < 0) {
        printf("P has been corrupted!\n");
        printf("P->numFields = %" PRId32 ", P->index[%" PRId32 "]=%" PRId32 ", P->index[%" PRId32 "]=%" PRId32 ".\n",
               P->numFields, fieldsBySize[i], P->index[fieldsBySize[i]],
               fieldsBySize[i]+1, P->index[fieldsBySize[i]+1]);
      } else {
        numBySize[size] += 1;
      }
      i++;
    }
    for (i=0; i<7; i++)
      printf("There are %" PRId32 " relations with %" PRId32 " large primes.\n", numBySize[i], i);
  }
}

/*************************************************************/
s32 getMaxEntry(llist_t *P)
/*************************************************************/
{ s32 size, i, N=0;

  size = P->index[P->numFields];
if (size > P->maxDataSize) {
  printf("Error: P is broken!\n");
}
  for (i=0; i<size; i++)
    if ((P->data[i] > N) && (P->data[i] != (s32)BAD_LP_INDEX)) 
      N=P->data[i];
  return N;
}
/*************************************************************/
int mkLT(llist_t *L, llist_t *P, s32 P0, s32 P1)
/*************************************************************/
/* Make a reverse-lookup table for entries of P in the range */
/* [P0, P1].                                                 */
/*************************************************************/
{ s32     i, j, size, cum, t, loc, pnum;

  /* Make one quick pass through to see how large the table will be,
     and how large each entry will be.
  */
  if (L->maxFields < P1-P0+5) {
    L->maxFields = P1-P0+5;
    free(L->index);
    L->index = (s32 *)lxmalloc((L->maxFields)*sizeof(s32),1);
  } 
  memset(L->index, 0x00, (L->maxFields)*sizeof(s32));
  size = 0;
  for (i=0; i<P->numFields; i++) {
    for (j=P->index[i]; j<P->index[i+1]; j++) {
      pnum = P->data[j];
      if ((pnum>=P0) && (pnum<=P1)) {
        L->index[P->data[j]-P0] += 1;
        size++;
      }
    }
  }
  if (size == 0)
    printf("mkLT: There appear to be no large primes in the specified range.\n");

  if (L->maxDataSize < size+10) {
    L->maxDataSize = size+10;
    free(L->data);
    L->data = (s32 *)lxmalloc(L->maxDataSize*sizeof(s32),1);
  }
  /* Now each L->index field has the number of relations that
     will be asociated with that prime. Go through and make
     it cumulative, so it will point to the start of relations
     for the corresponding prime.
  */
  size = P1-P0+1;
  cum=0;
  for (i=0; i<size; i++) {
    t = L->index[i];
    L->index[i] = cum;
    cum += t;
  }
  L->index[size] = cum; /* terminator. */
  L->numFields = size;

  /* Initialize the fields to all hold 0xFFFFFF, so we can tell where
     we have put stuff and where we haven't.
  */
  memset(L->data, 0xFF, L->maxDataSize*sizeof(s32));

  /* Finally, put the data where it belongs. */
  for (i=0; i<P->numFields; i++) {
    for (j=P->index[i]; j<P->index[i+1]; j++) {
      pnum = P->data[j];
      if ((pnum>=P0) && (pnum<=P1)) {
        loc = L->index[pnum-P0];
        while (L->data[loc] != (s32)0xFFFFFFFF)
          loc++;
/* CJM, 4/13/05: Just looking for fishy places where the code could
   be going bad, this looks like a possibility - see if anyone
   gets this message.
*/
        if (loc >= L->maxDataSize) {
          printf("** mkLT() is broken! Fix me!\n");
        } else {
          L->data[loc] = i;
        }
      }
    }
  }
  return 0;
}



/*************************************************************/
int removeLPSingletons(llist_t *R, llist_t *P)
/*************************************************************/
/* Remove relations having a large prime which does not      */
/* appear in any of the other relations. This is done as     */
/* follows: (1) Determine an upper bound for the number of   */
/* large primes appearing.                                   */
/* (2) Make an array having 4 bits for each possible large   */
/*     prime.                                                */
/* (3) Make a pass through all relations. For each large     */
/*     prime encountered, set the corresponding counter to   */
/*     MIN(15, counter+1).                                   */
/* Finally, make one more pass through the data again looking*/
/* at the large primes. Any relations having a large prime   */
/* whose counter is == 1 will be deleted.                    */
/* So the time cost is about 2*(number of relations), and    */
/* the space cost is (max. # of large primes)/2 bytes.       */
/*  Upon completion, 'P' will be reallocated to its new size */
/* plus LARGE_BUFFER.                                        */
/*  The counter field is a global variable, as it will be    */
/* useful to some other functions. If anybody else ever      */
/* `free's it, they had better set it back to NULL as well!  */
/*************************************************************/
{ s32     Np, i, j, h, size;
  s32     *removeList, rLSize, numRemove;
  s32     *tmpPtr;
  unsigned char ct=0;


  Np=getMaxEntry(P)+4; /* a little padding, just in case. */
  /* We need four bits per prime in the hash table: */
  counterSize = 1 + Np/2;
  if (counter != NULL) free(counter);
  counter = (unsigned char *)lxmalloc(counterSize*sizeof(unsigned char),1);
  memset(counter, 0x00, counterSize*sizeof(unsigned char));

  size = P->index[P->numFields];
  for (i=0; i<size; i++) {
    h = P->data[i];
    if (h != (s32)BAD_LP_INDEX) {
      if (h%2==0) {
        ct = (counter[h/2]&0x0F);
        if (ct < 2)
          counter[h/2] += 0x01;
      } else {
          ct = (counter[h/2]&0xF0)>>4;
        if (ct < 2) 
          counter[h/2] += 0x10;
      }
    }
  }

  numRemove=0;
  rLSize = 32768;
  removeList = (s32 *)lxmalloc(rLSize*sizeof(s32),1);
  /* Now make a pass through looking for primes which appeared only once. */
  for (j=0; j<P->numFields; j++) {
    for (i=P->index[j]; i<P->index[j+1]; i++) {
      h = P->data[i];
      if (h != (s32)BAD_LP_INDEX) {
        if (h%2==0) 
          ct = (counter[h/2]&0x0F);
        else
          ct = (counter[h/2]&0xF0)>>4;
      }
      if ((ct == 1) || (h==(s32)BAD_LP_INDEX)) {
        removeList[numRemove++] = j;
        if (numRemove >= (rLSize-1)) {
          rLSize += 4*32768;
          tmpPtr = lxrealloc(removeList, rLSize*sizeof(s32), 0);
          if (tmpPtr==NULL) {
            fprintf(stderr, "removeLPSingletons() Memory (re-)allocation error!\n");
            exit(-1);
          }
          removeList = tmpPtr;
        }
        break;
      }
    }
  }
  /* 'removeList' should already be sorted and unique, but whatever. */
  numRemove = mkUniqueS32s(removeList, numRemove);
  printf("Deleting %" PRId32 " singleton large primes.\n", numRemove);

  ll_deleteFields(P, removeList, numRemove);
  ll_resize(P, P->index[P->numFields] + LARGE_BUFFER);
  ll_deleteFields(R, removeList, numRemove);
//  printf("There are %" PRId32 " relations remaining.\n", R->numFields);
  free(removeList);
  return (int)numRemove;
}

/*************************************************************/
s32 mkUniqueS32s(s32 *x, s32 size)
/*************************************************************/
{ s32 i, j;

  if (size<=1) return size;
  qsort(x, size, sizeof(s32), cmpS32s);
  for (i=0, j=1; j<size; j++) 
    if (x[j] != x[i])
      x[++i] = x[j];
  return i+1;
}

/*************************************************************/
s32 mkUniquePairs(lpair_t *pairs, s32 numPairs)
/*************************************************************/
/* Make sure there are no two pairs with the same x value.   */
/*************************************************************/
{ s32 i, j, k;

  if (numPairs<=1) return numPairs;
  qsort(pairs, numPairs, sizeof(lpair_t), cmp_lpair_t);
  k=0;
  for (j=0,i=1; i<numPairs; i++) {
    if (pairs[i].x != pairs[j].x) {
      if (pairs[i].x != pairs[i].y) {
        pairs[++j].x = pairs[i].x;
        pairs[j].y = pairs[i].y;
      } else 
        printf("mkUniquePairs(): Warning: somebody tried to add a column to itself!\n");
    }
  }
  return j+1;
}



#define MAX_MERGE_OPS 300000
/*************************************************************/
int merge(llist_t *R, llist_t *P, llist_t *revP, s32 P0, s32 P1, int level)
/*************************************************************/
/* This does the actual combining of relations to eliminate  */
/* unbalanced large primes.                                  */
/*************************************************************/
{ s32      i, j, k, pnum, relnum, numAdds, c1, ck;
  s32      rwt[1024], minrwt, pivot;
  lpair_t *pairs;
  s32     *removeList, numRemove=0;

  removeList = (s32 *)lxmalloc(revP->index[revP->numFields]*sizeof(s32),1);
  pairs = (lpair_t *)lxcalloc(MAX_MERGE_OPS*sizeof(lpair_t),1);
  numAdds=0;
  for (i=0; (i<revP->numFields)&&(numAdds < MAX_MERGE_OPS); i++) {
    pnum = P0+i;
    k = FIELDSIZE(revP, i);
    if (k>1) {
      /* This large prime appears in k relation-sets. */
      k = MIN(k, 1024); /* k should never be anywhere near 1024! */
      /* How many relations are contributing to each relation set? */
//      minrwt=1024;
      for (j=0; j<k; j++) {
        relnum = revP->data[revP->index[i]+j];
        rwt[j] = FIELDSIZE(R, relnum);
//        minrwt = MIN(minrwt, rwt[j]);
      }
#if 0
      /* How many large primes are there in each relation set? */
      minpwt=64;
      for (j=0; j<k; j++) {
        relnum = revP->data[revP->index[i]+j];
        pwt[j] = FIELDSIZE(P, relnum);
        minpwt = MIN(minpwt, pwt[j]);
      }
      pivot=-1;
      for (j=0; j<k; j++) {
        if ((pwt[j]==minpwt) && ((pivot==-1) || (rwt[j]<rwt[pivot])))
          pivot = j;
      }
      if (pivot >= 0) {
        /* Check to see if we can use it. */
        if (minpwt > 1) pivot=-1;
        else if (rwt[pivot] > level) pivot=-1;
      }
#else
/* This could be much smarter. It should be possible
   to make relation-sets with lower total weight than now,
   and perhaps even to generate significantly more of them
   (since the poor weighting choices limit how many can be converted
   now).
    At the very least, we can temporarily ignore the weights,
   but go back after all merges are done and look for combinations
   of relation-sets which can be used to reduce the weights after-the-fact.

*/
      pivot=-1;
      minrwt=1024;
      for (j=0; j<k; j++) {
        relnum = revP->data[revP->index[i]+j];
        if (FIELDSIZE(P, relnum)==1) {
          if (FIELDSIZE(R, relnum) < minrwt) {
            pivot = j;
            minrwt = FIELDSIZE(R, relnum);
          }
        }
      }
      if (minrwt > level) pivot = -1; /* Too heavy - can't use it. */
#endif
      if (pivot >= 0) {
        ck = revP->data[revP->index[i]+pivot];
        /* Ok - see what additions we can do. */
        for (j=0; j<k; j++) {
          if (j != pivot) {
            c1 = revP->data[revP->index[i]+j];
//            if ((rwt[pivot]+rwt[j] <= level)&&(numAdds<MAX_MERGE_OPS)) {
//              pairs[numAdds].x = c1; pairs[numAdds++].y = ck;
//            }
            if (rwt[pivot]+rwt[j] <= level) {
              if (numAdds<MAX_MERGE_OPS) {
                pairs[numAdds].x = c1; pairs[numAdds++].y = ck;
              }
            } else {
              removeList[numRemove++] = revP->data[revP->index[i]+j];
            }
          }
        }
      }
    }
  }
  numAdds = mkUniquePairs(pairs, numAdds);
  if (numAdds > 0) {
    printf("Doing %" PRId32 " additions...\n", numAdds);
    ll_catFields(P, pairs, numAdds, 1);
    ll_catFields(R, pairs, numAdds, 1);
  }
  numRemove = mkUniqueS32s(removeList, numRemove);
  ll_deleteFields(P, removeList, numRemove);
  ll_resize(P, P->index[P->numFields] + LARGE_BUFFER);
  ll_deleteFields(R, removeList, numRemove);

  free(removeList);
  free(pairs);
  return 0;
}


/* The maximum # of `s32s' we can allocate for a reverse-lookup
   table. The larger this is, the more memory makePass will use,
   but it will also be faster. When I'm less tired, I'll figure
   the runtime and put a more precise comment here.
*/
#define MAX_PART_DATASIZE 2000000
/*************************************************************/
s32 makePass(llist_t *R, llist_t *P)
/*************************************************************/
{ s32 *fieldsBySize, startNumLP[7];
  s32  i, k, bd, pMax, numParts, part, partSize, P0, P1;
  s32  lastSize, numFulls, s0, s1;
  llist_t revP;

  /* First do some sanity checks. */
  fieldsBySize = (s32 *)lxmalloc(P->numFields*sizeof(s32),1);
  printf("Before sortByNumLP()... Doing ll_verify(P)...\n");
  if (ll_verify(P)) {
    printf("ll_verify() reported an error!\n");
    exit(-1);
  }
  printf("ll_verify() reports that 'P' appears to be intact.\n");
  sortByNumLP(fieldsBySize, P, "makePass:\n");
  printf("After sortByNumLP()... Doing ll_verify(P)...\n");
  if (ll_verify(P)) {
    printf("ll_verify() reported an error!\n");
    exit(-1);
  }
  printf("ll_verify() reports that 'P' appears to be intact.\n");
  free(fieldsBySize);

  /* Remove any singletons. But the function removeLPSingletons()
     won't take the liberty of re-running itself after it removes
     the singletons. As a result, removing some singletons could 
     create some more. So repeat until it's stable.
       This will also set up `counter', which counts the number of
     times a given large prime appears in some relation-set. 
  */
  s0 = R->numFields;
  do {
    lastSize = P->numFields;
    removeLPSingletons(R, P); 
  } while (P->numFields < lastSize);
  s1 = R->numFields;
  printf("Total: %" PRId32 " singletons deleted.\n", s0-s1);

  /* Sort the relation-sets on the number of large primes
     appearing. We don't actually sort 'P', but rather obtain
     a list of fields of P which is sorted in ascending order.
  */
  fieldsBySize = (s32 *)lxmalloc(P->numFields*sizeof(s32),1);
  sortByNumLP(fieldsBySize, P, "makePass:\n");
  /* Find the start of fields with 1 LP, 2 LP,... */
  startNumLP[0]=0;
  bd = P->numFields;
  for (i=1; i<=6; i++) {
    k=startNumLP[i-1];
    while ((k < bd) && (FIELDSIZE(P,fieldsBySize[k])<i))
      k++;
    startNumLP[i]=k;
  }

  /*****************************************************************/
  /* Now do as follows:                                            */
  /* (1) Find the maximum large prime appearing, maxP.             */
  /* (2) Partition the possible primes [0, maxP] into subsets      */
  /*     of manageable size.                                       */
  /* (3) For each partition of the primes, create a reverse-lookup */
  /*     table, so we can look the prime up to see in which        */
  /*     relation-sets it appears.                                 */
  /* (4) Use the reverse-lookup and counter tables to decide what  */
  /*     column additions to make, if any.                         */
  /*****************************************************************/
  pMax=getMaxEntry(P);
  numParts = 1 + P->index[P->numFields]/MAX_PART_DATASIZE;
  P0 = 0;
  P1 = pMax/numParts;
  partSize = P1+1;
  part=0;
  ll_init(&revP, P1+100, P->maxDataSize);
  do {
    printf("Doing merge on chunk %" PRId32 "/%" PRId32 "  (P0=%" PRId32 ", P1=%" PRId32 ")...\n", part+1, numParts, P0, P1);
    /* Make the reverse-lookup table. */
    mkLT(&revP, P, P0, P1);
    merge(R, P, &revP, P0, P1, (s32)(1.5*MAX_RELS_IN_FF));

    /* Free up the reverse-lookup table. */
    if (++part < numParts) {
      P0 = P1+1;
      P1 = MIN(P0+partSize, pMax);
    }
  } while (part < numParts);
  ll_clear(&revP); 
  free(fieldsBySize);
  /* Count the number of full relations: */
  numFulls=0;
  for (i=0; i<P->numFields; i++)
    if (FIELDSIZE(P, i)==0)
      numFulls++;
  return numFulls;
}
/*************************************************************/
int checkR(llist_t *R)
/*************************************************************/
{ s32 i, j, k;

  for (i=0; i<R->numFields; i++) {
    if (R->index[i+1]-R->index[i] == 0) {
      printf("Error: R field %" PRId32 " is empty!\n", i);
      exit(-1);
    }
    for (j=R->index[i]; j<R->index[i+1]; j++) {
      for (k=j+1; k<R->index[i+1]; k++) {
        if (R->data[j]==R->data[k]) {
          printf("Error: R field %" PRId32 " has duplicate entries!\n", i);
          exit(-1);
        }
      }
    }
  }
  return 0;
}

/*************************************************************/
s32 keepFulls(llist_t *R, llist_t *P)
/*************************************************************/
/* Remove any relation-sets which still have large primes.   */
/*************************************************************/
{ s32 i, *rem=NULL, numR, maxR;

  numR=0;
  maxR=8192;
  rem = (s32 *)lxmalloc(maxR*sizeof(s32),1);
  for (i=0; i<R->numFields; i++) {
    if (FIELDSIZE(P, i)>0) {
      if (numR >= maxR) {
        maxR += 4*8192;
        rem = lxrealloc(rem, maxR*sizeof(s32),0);
        if (rem==NULL) {
          printf("keepFulls(): Memory re-allocation error for remove!\n");
          exit(-1);
        }
      }
      rem[numR++] = i;
    }
  }
  if (numR) {
    ll_deleteFields(P, rem, numR);
    ll_deleteFields(R, rem, numR);
    free(rem);
  }
  return R->numFields;
}


/*************************************************************/
int combinedSize(llist_t *R, s32 a, s32 b)
/*************************************************************/
/* Figure how many relation-sets would be left if we added   */
/* sets: a <-- a + b.                                        */
/* Return value: (combined size) - (original size),          */
/* or 128 if we could not compute it (i.e., 128 is large     */
/* enough that the caller will see that and decide that it   */
/* would not be a good idea to perform the addition).        */
/*************************************************************/
/* It is assumed that the entries of each field of R are     */
/* sorted ascending!                                         */
/*************************************************************/
{ s32 i, j, r;
  s32 sa, sb, *A, *B;
  int  sizeDiff=0;

  sa = FIELDSIZE(R, a);
  sb = FIELDSIZE(R, b);
  if (sb > 2*sa) return 128;
  A = R->data+R->index[a];
  B = R->data+R->index[b];
  for (i=j=0; i<sa; i++) {
    r = A[i];
    while ((j<sb) && (B[j] < r)) {
      j++;
      sizeDiff++;
    }
    if ((j<sb) && (B[j] == r)) {
      sizeDiff--;
      j++;
    }
  }
  sizeDiff += (sb - j);
  return sizeDiff;
}

/*************************************************************/
s32 reduceRelSets(llist_t *R, llist_t *P)
/*************************************************************/
/* Attempt to use some full relation-sets to reduce the      */
/* the weights of other relation sets.                       */
/* For example, if we have the two LP-free relation sets:    */
/*   (12, 1002, 17332, 19282, 20004)                         */
/*   (1032, 17332, 19282, 20004, 21920)                      */
/* (and this certainly could happen, by construction), we can*/
/* add the second to the first to obtain relation-sets with  */
/* less total weight:                                        */
/*   (12, 1002, 1032, 21920)                                 */
/*   (1032, 17332, 19282, 20004, 21920)                      */
/* The relation-sets should have already been put through    */
/* keepFulls(), so that only LP-free relation-sets remain.   */
/*************************************************************/
/* The idea: Use mkLT to make a reverse-lookup table on the  */
/* relation-sets. Then, go through the relation-sets one at  */
/* a time, and reverse lookup each relation occurring. For   */
/* each relation occurring, check to see if adding it would  */
/* reduce the total weight by calling combinedSize(). If one */
/* of these additions would reduce the size, then by all     */
/* means do it!                                              */
/*************************************************************/
{ llist_t revR;
  s32    r0, r1, rMax, part, numParts, partSize;
  s32    i, j, a, b, r, numAdds, totDiff=0, initW, s;
  s32   *hash;
  int     bestDiff, bdb, d, c1, c2;
  lpair_t *pairs;

  initW = R->index[R->numFields];
  pairs = (lpair_t *)lxcalloc(MAX_MERGE_OPS*sizeof(lpair_t),1);
  numAdds=0;

  printf("Attempting to reduce weight of relation sets.\n");
  printf("Initial weight is: %" PRId32 "\n", R->index[R->numFields]);
  printf("Sorting relation-sets..."); fflush(stdout);
  for (i=0; i<R->numFields; i++) {
    s = FIELDSIZE(R, i);
    qsort(R->data+R->index[i], s, sizeof(s32), cmpS32s);
  }
  printf("Done.\n");

  rMax=getMaxEntry(R);
  numParts = 1 + R->index[R->numFields]/MAX_PART_DATASIZE;
  r0 = 0;
  r1 = rMax/numParts;
  partSize = r1+1;
  part=0;

  if (!(hash = (s32 *)lxmalloc((R->numFields)*2*sizeof(s32),0))) {
    printf("reduceRelSets() : Memory allocation error for hash!\n");
    exit(-1);
  }
  memset(hash, 0x00, 2*R->numFields*sizeof(s32));
  for (a=0; a<R->numFields; a++) {
    for (i=R->index[a]; i<R->index[a+1]; i++) {
      r = R->data[i]%64;
      hash[2*a+(r/32)] |= BIT(r&0x1F);
    }
  }
  init_bitcount();

  ll_init(&revR, r1-r0+100, R->maxDataSize);
  do {
    printf("Making lookup table for chunk %" PRId32 " / %" PRId32 ": [%" PRId32 ", %" PRId32 ")...", 
            part+1, numParts, r0, r1);
    fflush(stdout);
    mkLT(&revR, R, r0, r1);
    printf("Done.\n");
    for (a=0; (a<R->numFields) && (numAdds < MAX_MERGE_OPS); a++) {
      bestDiff = 128; bdb=0;
      if (a%10000 == 0) 
        printTmp("Done examining %ld of %ld relation sets...", a, R->numFields);
      if (FIELDSIZE(R,a) > 6) {
        for (i=R->index[a]; i<R->index[a+1]; i++) {
          r = R->data[i];
          if ((r >= r0) && (r < r1)) {
            for (j=revR.index[r-r0]; j<revR.index[r-r0+1]; j++) {
              b = revR.data[j];
              if (a != b) {
                c1 = hash[2*a]&hash[2*b];
                c2 = hash[2*a+1]&hash[2*b+1];
                if (BITCOUNT(c1) + BITCOUNT(c2) > 2) {
                  d = MIN(bestDiff, combinedSize(R, a, b));
                  if (d < bestDiff) {
                    bestDiff = d;
                    bdb = b;
                  }
                }
              }
            }
          }
        }
      }
      if ((bestDiff < 0) && (numAdds < MAX_MERGE_OPS)) {
        pairs[numAdds].x = a;
        pairs[numAdds++].y = bdb;
        totDiff += bestDiff;
      }
    }
    part++;
    numAdds = mkUniquePairs(pairs, numAdds);
    /* I think mkUniquePairs() really doesn't do quite enough. We should
       also make sure that we don't have a circular addition like:
         x <-- x+y
         y <-- y+x
       But it seems to be okay, so I'll let it go for now. If this is a problem,
       it will probably manifest as a matrix that looks okay and solves okay,
       but something odd happens at the square root step.
    */
    if (numAdds > 0) {
      printf("Doing %" PRId32 " additions to reduce relation-set weight...\n", numAdds);
      ll_catFields(R, pairs, numAdds, 1);
      printf("Current weight is: %" PRId32 "\n", R->index[R->numFields]);
      numAdds=0;
    }
    r0 += partSize;
    r1 =  MIN(r1+partSize, rMax);
    
  } while ((part < numParts) && (numAdds < MAX_MERGE_OPS));
  ll_clear(&revR); 
  free(hash);
  free(pairs);
//  free(bitcount);
  printf("\nfinal weight is: %" PRId32 ".\n", R->index[R->numFields]);
  return R->index[R->numFields] - initW;
}

/*************************************************************/
s32 removeHeavyRelSets(llist_t *R, llist_t *P, int maxRelsInRS)
/*************************************************************/
/* Remove any relation-set having more than maxRelsInRS      */
/* relations.                                                */
/* And dump a list of the number of relation-sets with each  */
/* weight to stdout. This is to help the user if he/she needs*/
/* to re-run with a different maxRelsInFF.                   */
/*************************************************************/
{ s32 i, *remove=NULL, numR, maxR, fSize;
  u32 byWt[10*MAX_RELS_IN_FF], cum, cwt;

  numR=0;
  maxR=8192;
  remove = (s32 *)lxmalloc(maxR*sizeof(s32), 1);
  for (i=0; i<10*MAX_RELS_IN_FF; i++)
    byWt[i]=0;
  for (i=0; i<R->numFields; i++) {
    fSize = FIELDSIZE(R,i);
    byWt[fSize]+=1;
    if (fSize>maxRelsInRS) {
      if (numR >= maxR) {
        maxR += 4*8192;
        remove = lxrealloc(remove, maxR*sizeof(s32), 1);
      }
      remove[numR++] = i;
    }
  }
  printf("Before deleting relation sets heavier than wt %d, there were:\n", maxRelsInRS);
  printf("Wt  |  # R-S   | Cum. R-S | Cum. wt.\n");
  printf("---------------------------------\n");
  for (i=0,cum=0,cwt=0; i<10*MAX_RELS_IN_FF; i++) {
    cum += byWt[i];
    cwt += i*byWt[i];
    if (byWt[i]>0) {
      printf("%3" PRId32 " |%10" PRIu32 "|%10" PRIu32 "|%" PRIu32 "\n", i, byWt[i],cum,cwt);
    }
  }
  printf("--------------------------------------------------------------\n");
  printf("Wt = Weight\n# R-S = Number of relation-sets with this weight.\n");
  printf("Cum. R-S = # of relation-sets with at most this weight.\n");
  printf("Cum. wt. = cumulative weight of relation-sets upto here.\n");
  printf("--------------------------------------------------------------\n");
    
  if (numR) {
    ll_deleteFields(P, remove, numR);
    ll_deleteFields(R, remove, numR);
    free(remove);
  }
  return R->numFields;
}

    



/*************************************************************/
s32 combParts(llist_t *R, llist_t *P, int maxRelsInFF, s32 minFF)
/*************************************************************/
/*Input:
     'P' is a list containing the large primes. It should not
     actually contain the primes, but rather a unique integer associated
     with each large prime (i.e., algebraic and rational primes have
     disjoint sets of indicies).
     On input, it should be indexed by relation number.
     'R' should point to an uninitialized llist_t.
  Output:
     'R' is a list of full relations. Each field contains
     a list of relation numbers contributing to the 
     corresponding full relation. In Cavallar's nomenclature, 'R'
     will be a collection of relation sets.
     'P' will be modified in the process.
  Return value: The number of full relations built.

  This will be a rather memory intensive function, so be sure to
  free up any unneeded memory before calling us! Everything will
  be done in RAM instead of with disk I/O. So, for example, with
  50M relations having an average of 3 large primes each, you will
  be looking at a 'P' structure with 800MB (assuming 4 byte s32s).
  Then, tack onto that the RAM for the 'R' structure and miscellaneous,
  and you could easily be around 1.5BGB of RAM! 

  This function is unwritten at the moment - I want to put a little
  thought into it so I get it right this time and get something near
  an optimal algorithm which is at the same time efficient.
*/
/**************************************************************/ 
{ s32 i, wt0, wt1;
  s32 lastFull, full;
  int  pass=0;
  double shrink;

  /* We can reallocate if necessary. So if this seems to be way
     too much memory, decrease it at will.
  */
  printf("combParts() Doing ll_verify(P)...\n");
  if (ll_verify(P)) {
    printf("ll_verify() reported an error!\n");
    exit(-1);
  }
  printf("ll_verify() reports that 'P' appears to be intact.\n");

/* CJM, 3/21/05 : Changed from P->maxDataSize/4 in an effort
   to combat some of the large realloc's which seem to be causing
   problems. I didn't think too much about it before changing it,
   so it could be overkill, or could still be too small.
*/
  if (ll_init(R, P->maxDataSize, P->maxFields)) {
    fprintf(stderr, "combParts: ll_init() reports severe error!\n");
    exit(-1);
  }
  for (i=0; i<P->numFields; i++) {
    ll_appendField(R, &i, 1);
  }

  full = -1; /* Force at least two passes. */
  do {
    printf("  pass %d...\n", ++pass);
    lastFull = full;
    full = makePass(R, P);
    checkR(R);
    printf("* There are now %" PRId32 " full relations.\n", full);
  } while (lastFull < full);

  /* Drop any relation-sets still containing a large prime: */
  keepFulls(R, P); 
  printf("After keepFulls(), R->numFields = %" PRId32 "\n", R->numFields);

  /* Don't bother with the weight reduction unless we're close
     to having enough relations.
  */
  if (full < minFF) {
    return full;
  }
#ifdef RS_WT_REDUCTION
  printf("Reducing the weight of relation sets. This is painfully\n");
  printf("slow at the moment, but it's worth it.\n");
  do {
    wt0 = R->index[R->numFields];
    reduceRelSets(R, P);
    wt1 = R->index[R->numFields];
    msgLog("", "reduceRelSets dropped relation-set weight from %" PRId32 " to %" PRId32 ".",
           wt0, wt1);
    shrink = (double)(wt0-wt1)/wt0;
  }  while (shrink > 0.15);
#endif
  full = removeHeavyRelSets(R, P, maxRelsInFF);
  msgLog("", "After removing heavy rel-sets, weight is %" PRId32 ".", R->index[R->numFields]);
  printf("After removing heavy rel-sets, weight is %" PRId32 ".\n", R->index[R->numFields]);
  if (ll_verify(R)) {
    printf("ll_verify() reported an error for R!\n");
    exit(-1);
  }
  if (ll_verify(P)) {
    printf("ll_verify() reported an error for P!\n");
    exit(-1);
  }
  return full;
}
