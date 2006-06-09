/**************************************************************
 * combparts_tpie.cpp
 * Routines for combining partials to form full relations.
 * Copyright 2004-2005 Chris Monico, Anton Korobeynikov
 **************************************************************/

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
#include "ggnfs.h"
#include "llist_tpie.h"
#include "templates.h"

#define LARGE_BUFFER 200000
#define BIT(_i) ((0x00000001)<<(_i))
#define BITCOUNT(_c) (bitcount[_c&0x0000FFFF] + bitcount[(_c&0xFFFF0000)>>16])
#define MAX_RELS_IN_FF 48
#define MAX_RELS_IN_COL 128
#define MAX_PART_DATASIZE 2000000
#define MAX_MERGE_OPS 300000
#define DUMP_RELS

/* Globals: */
unsigned char *counter=NULL;
s32 counterSize;
int  *bitcount=NULL;

bool checkRelSets(LSList&R, LSList& dataList);

void init_bitnumFields() {
	int i, j, c;

	bitcount = (int *)malloc(65536*sizeof(int));
	for (i=0; i<65536; i++) {
		for (c=0, j=0; j<16; j++)  {
			if (i&BIT(j)) {
				c++;
			}
		}
		bitcount[i] = c;
	}
}

s32 mkUniqueS32s(s32 *x, s32 size) {
	s32 i, j;

	if (size<=1) {
		return size;
	}
	qsort(x, size, sizeof(s32), cmpS32s);
	for (i=0, j=1; j<size; j++)  {
		if (x[j] != x[i]) {
			x[++i] = x[j];
		}
	}
	return i+1;
}

/*************************************************************
 * Make sure there are no two pairs with the same x value.
 *************************************************************/
s32 mkUniquePairs(lpair_t *pairs, s32 numPairs) {
	s32 i, j;

	if (numPairs<=1) return numPairs;
	qsort(pairs, numPairs, sizeof(lpair_t), cmp_lpair_t);
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

/*************************************************************/
void sortByNumLP(s32 *fieldsBySize, LSList& L, char *str) {
	s32 i, size;
	s32 numBySize[8];

	L.getSortOnFieldSize(fieldsBySize);
	if (str != NULL) {
		printf(str);
		for (i=0; i<8; i++) { /* Overkill - 4 or 6 would do. */
			numBySize[i]=0;
		}
		i=0;
		while (i<(s32)L.numFields()) {
			if (fieldsBySize[i] > (s32)L.numFields()) {
				printf("?!?WTF?!? fieldsBySize[%" PRId32 "]=%" PRId32 " vs. numFields=%" PRId32 ".\n",
						i,fieldsBySize[i], (s32)L.numFields());
			}
			size = L.fieldSize(fieldsBySize[i]);
			if (size > 6)  {
				printf("?!Error!? Relation has more than 6 large primes!\n");
			}
			else if (size < 0) {
				printf("LList has been corrupted! Size < 0\n");
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
LookupTable* mkLT(LSList& P, s32 P0, s32 P1) {
	s32 i, j, size, pnum;

	/* Make one quick pass through to see how large the table will be,
	 * and how large each entry will be.
	 */

	s32* size_buffer = new s32[5+P1-P0];
	memset(size_buffer,0,(5+P1-P0)*sizeof(s32));

	size = 0;
	for (i=0; i<(s32)P.numFields(); ++i) {
		LSListField* lEntry=P.fetchField(i);
		for (j=0; j<(s32)lEntry->numEntries(); ++j) {
			pnum = lEntry->at(j);
			if ((pnum>=P0) && (pnum<=P1)) {
				size_buffer[pnum-P0] += 1;
				size++;
			}
		}
		P.releaseField(lEntry);
	}
	if (size == 0) {
		printf("mkLT: There appear to be no large primes in the specified range.\n");
	}

	for (i=1; i<P1-P0+1; ++i) {
		size_buffer[i]+=size_buffer[i-1];
	}


	LookupTable *L = new LookupTable(size_buffer,P1-P0+1);
	delete[] size_buffer;

	for (i=0; i<(s32)P.numFields(); ++i) {
		LSListField* lEntry=P.fetchField(i);
		for (j=0; j<(s32)lEntry->numEntries(); ++j) {
			pnum = lEntry->at(j);
			if ((pnum>=P0) && (pnum<=P1)) {
				L->appendEntry(pnum-P0,i);
			}
		}
		P.releaseField(lEntry);
	}
	return L;
}

/*************************************************************
 * Remove relations having a large prime which does not
 * appear in any of the other relations. This is done as
 * follows: (1) Determine an upper bound for the number of
 * large primes appearing.
 * (2) Make an array having 4 bits for each possible large
 * prime.
 * (3) Make a pass through all relations. For each large
 * prime encountered, set the corresponding counter to
 * MIN(15, counter+1).
 * Finally, make one more pass through the data again looking
 * at the large primes. Any relations having a large prime
 * whose counter is == 1 will be deleted.
 * So the time cost is about 2*(number of relations), and
 * the space cost is (max. # of large primes)/2 bytes.
 * Upon completion, 'P' will be reallocated to its new size
 * plus LARGE_BUFFER.
 * The counter field is a global variable, as it will be
 * useful to some other functions. If anybody else ever
 * `free's it, they had better set it back to NULL as well!
 **************************************************************/
int removeLPSingletons(LSList& R, LSList& P) {
	s32 Np, i, j, h;
	s32 numRemove;
	s32 rLSize, *removeList, *tmpPtr;
	unsigned char ct=0;

	Np=P.getMaxEntry()+4; /* a little padding, just in case. */
	/* We need four bits per prime in the hash table: */
	counterSize = 1 + Np/2;
	if (counter != NULL) {
		free(counter);
	}
	counter = (unsigned char *)calloc(counterSize,sizeof(unsigned char));

	for (i=0; i<(s32)P.numFields(); ++i) {
		LSListField* lEntry=P.fetchField(i);
		for (j=0; j<(s32)lEntry->numEntries(); ++j) {
			h = lEntry->at(j);
			if (h != (s32)BAD_LP_INDEX) {
				if (h%2==0) {
					ct = (counter[h/2]&0x0F);
					if (ct < 2) {
						counter[h/2] += 0x01; 
					}
				} else {
					ct = (counter[h/2]&0xF0)>>4;
					if (ct < 2) {
						counter[h/2] += 0x10;
					}
				}
			}
		}
		P.releaseField(lEntry);
	}

	numRemove=0;
	/* Now make a pass through looking for primes which appeared only once. */

	rLSize = 32768;
	// FIXME!
	removeList = (s32 *)calloc(rLSize,sizeof(s32));

	for (i=0; i<(s32)P.numFields(); ++i) {
		LSListField* pEntry=P.fetchField(i);
		for (j=0; j<(s32)pEntry->numEntries(); ++j) {
			h = pEntry->at(j);
			if (h != (s32)BAD_LP_INDEX) {
				if (h%2==0) {
					ct = (counter[h/2]&0x0F);
				}
				else {
					ct = (counter[h/2]&0xF0)>>4;
				}
			}
			if ((ct == 1) || (h==(s32)BAD_LP_INDEX)) {
				removeList[numRemove++] = i;
				if (numRemove >= (rLSize-1)) {
					rLSize += 4*32768;
					tmpPtr = (s32*)realloc(removeList, rLSize*sizeof(s32));
					if (tmpPtr==NULL) {
						fprintf(stderr, "removeLPSingletons() Memory (re-)allocation error!\n");
						exit(-1);
					}
					removeList = tmpPtr;
				}
				break;
			}
		}
		P.releaseField(pEntry);
	}
	numRemove = mkUniqueS32s(removeList, numRemove);
	printf("Deleting %" PRId32 " singleton large primes.\n", numRemove);
	P.deleteFields(removeList,numRemove);
	R.deleteFields(removeList,numRemove);
	free(removeList);

	// printf("There are %" PRId32 " relations remaining.", (s32)R.numFields());
	return (int)numRemove;
}


/*************************************************************/
int merge(LSList& R, LSList& P, LookupTable& revP, s32 P0, s32 P1, int level) {
	s32 j, k, pnum, relnum, numAdds, c1, ck;
	s32 rwt[1024], minrwt, pivot;
	lpair_t *pairs;

	pairs = (lpair_t *)calloc(MAX_MERGE_OPS, sizeof(lpair_t)); 
	numAdds=0;
	for (u32 i=0; (i<revP.numFields())&&(numAdds < MAX_MERGE_OPS); ++i) {
		pnum = P0+i;
		LSListField* revpEntry = revP.fetchField(i);
		k = revpEntry->numEntries();
		if (k>1) {
			/* This large prime appears in k relation-sets. */
			k = (k<1024 ? k: 1024); /* k should never be anywhere near 1024! */
			/* How many relations are contributing to each relation set? */
			minrwt=1024;
			for (j=0; j<k; ++j) {
				relnum = revpEntry->at(j);
				rwt[j] = R.fieldSize(relnum);
				minrwt = (minrwt < rwt[j] ? minrwt : rwt[j]);
			}
			pivot=-1;
			minrwt=1024;
			for (j=0; j<k; j++) {
				relnum = revpEntry->at(j);
				if (P.fieldSize(relnum)==1) {
					s32 rfSize = R.fieldSize(relnum);
					if (rfSize < minrwt) {
						pivot = j;
						minrwt = rfSize;
					}
				}
			}
			if (minrwt > level) {
				pivot = -1; /* Too heavy - can't use it. */
			}
			if (pivot >= 0) {
				ck = revpEntry->at(pivot);
				/* Ok - see what additions we can do. */
				for (j=0; j<k; ++j) {
					if (j != pivot) {
						c1 = revpEntry->at(j);
						if ((rwt[pivot]+rwt[j] <= level)&&(numAdds<MAX_MERGE_OPS)) {
							pairs[numAdds].x = c1; pairs[numAdds++].y = ck;
						} // if ((rwt[pivot]...
					} // if (j != pivot) ..
				} // for (j ...
			} //if (pivot ...
		} // if (k>1)
		revP.releaseField(revpEntry);
	} // for (i ..
	numAdds = mkUniquePairs(pairs, numAdds);
	if (numAdds > 0) {
		printf("Doing %" PRId32 " additions...\n", numAdds);
		P.catFields(pairs, numAdds, 1);
		R.catFields(pairs, numAdds, 1);
	}
	free(pairs);
	return 0;
}

s32 keepFulls(LSList& R, LSList& P) {
	s32 *rem=NULL, numR, maxR;

	numR=0;
	maxR=8192;
	if (!(rem = (s32 *)malloc(maxR*sizeof(s32)))) {
		printf("keepFulls(): Memory allocation error for remove!\n");
		exit(-1);
	}
	for (u32 i=0; i<R.numFields(); i++) {
		if (P.fieldSize(i)>0) {
			if (numR >= maxR) {
				maxR += 8192;
				rem = (s32*)realloc(rem, maxR*sizeof(s32));
				if (rem==NULL) {
					printf("keepFulls(): Memory re-allocation error for remove!\n");
					exit(-1);
				}
			}
			rem[numR++] = i;
		}
	}
	if (numR) {
		numR = mkUniqueS32s(rem, numR);
		P.deleteFields(rem, numR);
		R.deleteFields(rem, numR);
		free(rem);
	}
	return R.numFields();
}

/*************************************************************
 * Figure how many relation-sets would be left if we added
 * sets: a <-- a + b.
 * Return value: (combined size) - (original size),
 * or 128 if we could not compute it (i.e., 128 is large
 * enough that the caller will see that and decide that it
 * would not be a good idea to perform the addition).
 * It is assumed that the entries of each field of R are
 * sorted ascending!
 **************************************************************/
int combinedSize(LSList& R, s32 a, s32 b) {
	s32 i, j, r;
	s32 sa, sb;
	int	sizeDiff=0;

	sa = R.fieldSize(a);
	sb = R.fieldSize(b);
	if (sb > 2*sa) return 128;
	LSListField* A = R.fetchField(a);
	LSListField* B = R.fetchField(b);
	for (i=j=0; i<sa; i++) {
		r = A->at(i);
		while ((j<sb) && (B->at(j) < r)) {
			j++;
			sizeDiff++;
		}
		if ((j<sb) && (B->at(j) == r)) {
			sizeDiff--;
			j++;
		}
	}
	sizeDiff += (sb - j);
	R.releaseField(A);
	R.releaseField(B);
	return sizeDiff;
}

/*************************************************************
 * Remove any relation-set having more than maxRelsInRS
 * relations.
 * And dump a list of the number of relation-sets with each
 * weight to stdout. This is to help the user if he/she needs
 * to re-run with a different maxRelsInFF.
 **************************************************************/
s32 removeHeavyRelSets(LSList& R, LSList& P, int maxRelsInRS) {
	s32 *remove=NULL, numR, maxR, fSize;
	u32 byWt[10*MAX_RELS_IN_FF];

	numR=0;
	maxR=8192;
	if (!(remove = (s32 *)malloc(maxR*sizeof(s32)))) {
		printf("removeHeavyRelSets(): Memory allocation error for remove!\n");
		exit(-1);
	}
	for (u32 i=0; i<10*MAX_RELS_IN_FF; i++)
		byWt[i]=0;
	for (u32 i=0; i<R.numFields(); i++) {
		fSize = R.fieldSize(i);
		byWt[fSize]+=1;
		if (fSize>maxRelsInRS) {
			if (numR >= maxR) {
				maxR += 8192;
				remove = (s32*)realloc(remove, maxR*sizeof(s32));
				if (remove==NULL) {
					printf("removeHeavyRelSets(): Memory re-allocation error for remove!\n");
					exit(-1);
				}
			}
			remove[numR++] = i;
		}
	}
	printf("Before deleting relation sets heavier than wt %d, there were:\n", maxRelsInRS);
	printf("Wt	|	# R-S	 | Cum. R-S | Cum. wt.\n");
	printf("---------------------------------\n");
	for (u32 i=0,cum=0,cwt=0; i<10*MAX_RELS_IN_FF; i++) {
		cum += byWt[i];
		cwt += i*byWt[i];
		if (byWt[i]>0) {
			printf("%3" PRId32 " |%10" PRId32 "|%10" PRId32 "|%" PRId32 "\n", i, byWt[i],cum,cwt);
		}
	}
	printf("--------------------------------------------------------------\n");
	printf("Wt = Weight\n# R-S = Number of relation-sets with this weight.\n");
	printf("Cum. R-S = # of relation-sets with at most this weight.\n");
	printf("Cum. wt. = cumulative weight of relation-sets upto here.\n");
	printf("--------------------------------------------------------------\n");

	if (numR) {
		numR = mkUniqueS32s(remove, numR);
		P.deleteFields(remove, numR);
		R.deleteFields(remove, numR);
		free(remove);
	}
	return R.numFields();
}

int checkR(LSList& R) {
	for (u32 i=0; i<R.numFields(); ++i) {
		LSListField* cField = R.fetchField(i);
		if (cField->numEntries() == 0) {
			printf("Error: R field %" PRId32 " is empty!\n", i);
			exit(-1);
		}
		for (u32 j=0; j<cField->numEntries(); ++j) {
			for (u32 k=j+1; k<cField->numEntries(); ++k) {
				if (cField->at(j)==cField->at(k)) {
					printf("Error: R field %" PRId32 " has duplicate entries!\n", i);
					exit(-1);
				}
			}
		}
		R.releaseField(cField);
	}
	return 0;
}

/*************************************************************/
s32 makePass(LSList& R, LSList& P) {
	s32	*fieldsBySize, startNumLP[7];
	s32	i, k, bd, pMax, numParts, part, partSize, P0, P1;
	s32	numFulls, s0, s1;
	u32	lastSize;
	LookupTable* revP;

	/* First do some sanity checks. */
	if (!(fieldsBySize = (s32 *)malloc(P.numFields()*sizeof(s32)))) {
		fprintf(stderr, "makePass() Severe memory allocation error!\n");
		exit(-1);
	}
	sortByNumLP(fieldsBySize, P, "makePass:\n");
	free(fieldsBySize);

	/* Remove any singletons. But this function won't take the liberty
	 * of re-running itself after it removes the singletons. As a result,
	 * removing some singletons could create some more. So repeat until
	 * it's stable.
	 * This will also set up `counter', which counts the number of
	 * times a given large prime appears in some relation-set. 
	 */
	s0 = R.numFields();
	do {
		lastSize = P.numFields();
		removeLPSingletons(R, P); 
	} while (P.numFields() < lastSize);
	s1 = R.numFields();
	printf("Total: %" PRId32 " singletons deleted.\n", s0-s1);

	/* Sort the relation-sets on the number of large primes
	 * appearing. We don't actually sort 'P', but rather obtain
	 * a list of fields of P which is sorted in ascending order.
	 */
	if (!(fieldsBySize = (s32 *)malloc(P.numFields()*sizeof(s32)))) {
		fprintf(stderr, "makePass() Severe memory allocation error!\n");
		exit(-1);
	}
	sortByNumLP(fieldsBySize, P, "makePass:\n");
	/* Find the start of fields with 1 LP, 2 LP,... */
	startNumLP[0]=0;
	bd = P.numFields();
	for (i=1; i<=6; i++) {
		k=startNumLP[i-1];
		while ((k < bd) && (P.fieldSize(fieldsBySize[k])<i))
			k++;
		startNumLP[i]=k;
	}

	/*****************************************************************
	 * Now do as follows:
	 * (1) Find the maximum large prime appearing, maxP.
	 * (2) Partition the possible primes [0, maxP] into subsets
	 * of manageable size.
	 * (3) For each partition of the primes, create a reverse-lookup
	 * table, so we can look the prime up to see in which
	 * relation-sets it appears.
	 * (4) Use the reverse-lookup and counter tables to decide what
	 * column additions to make, if any.
	 *****************************************************************/
	pMax=P.getMaxEntry();
	printf("Current weight is: %" PRId32 "				 \n",P.Weight());
	numParts = 1 + P.Weight()/MAX_PART_DATASIZE;
	P0 = 0;
	P1 = pMax/numParts;
	partSize = P1+1;
	part=0;
	do {
		printf("Doing merge on chunk %" PRId32 "/%" PRId32 "  (P0=%" PRId32 ", P1=%" PRId32 ")...\n", part+1, numParts, P0, P1);
		/* Make the reverse-lookup table. */
		revP = mkLT(P, P0, P1);
		merge(R, P, *revP, P0, P1, (s32)(1.5*MAX_RELS_IN_FF));

		/* Free up the reverse-lookup table. */
		delete revP;
		if (++part < numParts) {
			P0 = P1+1;
			P1 = (P0+partSize < pMax ? P0+partSize : pMax);
		}
	} while (part < numParts);
	free(fieldsBySize);
	/* Count the number of full relations: */
	numFulls=0;
	for (u32 i=0; i<P.numFields(); i++) {
		if (P.fieldSize(i)==0) {
			numFulls++;
		}
	}
	return numFulls;
}

/*************************************************************
 * Attempt to use some full relation-sets to reduce the
 * the weights of other relation sets.
 * For example, if we have the two LP-free relation sets:
 * (12, 1002, 17332, 19282, 20004)
 * (1032, 17332, 19282, 20004, 21920)
 * (and this certainly could happen, by construction), we can
 * add the second to the first to obtain relation-sets with
 * less total weight:
 * (12, 1002, 1032, 21920)
 * (1032, 17332, 19282, 20004, 21920)
 * The relation-sets should have already been put through
 * keepFulls(), so that only LP-free relation-sets remain.	 
 * The idea: Use mkLT to make a reverse-lookup table on the
 * relation-sets. Then, go through the relation-sets one at
 * a time, and reverse lookup each relation occurring. For
 * each relation occurring, check to see if adding it would
 * reduce the total weight by calling combinedSize(). If one
 * of these additions would reduce the size, then by all
 * means do it!
 **************************************************************/

s32 reduceRelSets(LSList& R, LSList& P) { 
	s32	r0, r1, rMax, part, numParts, partSize;
	s32	b, r, numAdds, totDiff=0, initW;
	s32	 *hash;
	int	 bestDiff, bdb, d, c1, c2;
	lpair_t *pairs;
	LookupTable* revR;

	initW = R.Weight();
	pairs = (lpair_t *)calloc(MAX_MERGE_OPS, sizeof(lpair_t));
	numAdds=0;

	printf("Attempting to reduce weight of relation sets.\n");
	printf("Initial weight is: %" PRId32 "\n", initW);
	printf("Sorting relation-sets..."); fflush(stdout);
	for (u32 i=0; i<R.numFields(); ++i) {
		R.sortEntries(i);
	}
	printf("Done.\n");

	rMax=R.getMaxEntry();
	numParts = 1 + R.Weight()/MAX_PART_DATASIZE;
	r0 = 0;
	r1 = rMax/numParts;
	partSize = r1+1;
	part=0;

	if (!(hash = (s32 *)malloc((R.numFields())*2*sizeof(s32)))) {
		printf("reduceRelSets() : Memory allocation error for hash!\n");
		exit(-1);
	}
	memset(hash, 0x00, 2*R.numFields()*sizeof(s32));
	for (u32 a=0; a<R.numFields(); ++a) {
		LSListField* cField = R.fetchField(a);
		for (u32 i=0; i<cField->numEntries(); ++i) {
			r = (cField->at(i))%64;
			hash[2*a+(r/32)] |= BIT(r&0x1F);
		}
		R.releaseField(cField);
	}
	init_bitnumFields();

	do {
		printf("Making lookup table for chunk %" PRId32 " / %" PRId32 ": [%" PRId32 ", %" PRId32 ")...",
				part+1, numParts, r0, r1); 
		fflush(stdout);
		revR = mkLT(R, r0, r1);
		printf("Done.\n");
		for (u32 a=0; (a<R.numFields()) && (numAdds < MAX_MERGE_OPS); a++) {
			bestDiff = 128; bdb=0;
			if (a%10000 == 0) {
				printTmp("Done examining %" PRId32 " of %" PRId32 " relation sets...", a, (s32)R.numFields());
			}
			if (R.fieldSize(a) > 6) {
				LSListField* cField = R.fetchField(a);
				for (u32 i=0; i<cField->numEntries(); ++i) {
					r = cField->at(i);
					if ((r >= r0) && (r < r1)) {
						LSListField* revrField = revR->fetchField(r-r0);
						for (u32 j=0; j<revrField->numEntries(); ++j) {
							b = revrField->at(j);
							if ((s32)a != b) {
								c1 = hash[2*a]&hash[2*b];
								c2 = hash[2*a+1]&hash[2*b+1];
								if (BITCOUNT(c1) + BITCOUNT(c2) > 2) {
									s32 cs = combinedSize(R, a, b);
									d = (bestDiff < cs ? bestDiff : cs);
									if (d < bestDiff) {
										bestDiff = d;
										bdb = b;
									}
								}
							}			
						}
						revR->releaseField(revrField);
					}
				}
				R.releaseField(cField);
			}
			if ((bestDiff < 0) && (numAdds < MAX_MERGE_OPS)) {
				pairs[numAdds].x = a;
				pairs[numAdds++].y = bdb;
				totDiff += bestDiff;
			}
		}
		delete revR;
		part++;
		numAdds = mkUniquePairs(pairs, numAdds);
		/* I think (CJM) mkUniquePairs() really doesn't do quite enough. We should
		 * also make sure that we don't have a circular addition like:
		 * x <-- x+y
		 * y <-- y+x
		 * But it seems to be okay, so I'll let it go for now. If this is a problem,
		 * it will probably manifest as a matrix that looks okay and solves okay,
		 * but something odd happens at the square root step.
		 * ******
		 * The real life is much more complicated ;) The original implementation using
		 * in-memory llist_t is ummuned from this sort of error due to original implementation.
		 * of ll_catFields(). Unfortunately, while utilizing LSList, we cannot algorithm
		 * 'in the spirit' of ll_catFields() due to huge amount of seek operation. That's why,
		 * we need to add one more check to LSList::catFields(). I hope, that mkUniquePairs() will
		 * be rewritten in some efficient way to handle circular addition(s).
		 */
		if (numAdds > 0) {
			printf("Doing %" PRId32 " additions to reduce relation-set weight...\n", numAdds);
			R.catFields(pairs, numAdds, 1);
			checkR(R);
			printf("Current weight is: %" PRId32 "\n", R.Weight());
			numAdds=0;
		}
		r0 += partSize;
		r1 = (r1+partSize < rMax ? r1+partSize : rMax);

	} while ((part < numParts) && (numAdds < MAX_MERGE_OPS));
	free(hash);
	free(pairs);
	free(bitcount);
	s32 finalW = R.Weight();
	printf("\nfinal weight is: %" PRId32 ".\n", finalW);
	return finalW - initW;
}

/*************************************************************
 * 	Input:
 * 	'P' is a list containing the large primes. It should not
 * 	actually contain the primes, but rather a unique integer 
 * 	associated with each large prime (i.e., algebraic and
 * 	rational primes have disjoint sets of indicies).
 * 	On input, it should be indexed by relation number.
 * 	'R' should point to an uninitialized LSList.
 * 	Output:
 * 	'R' is a list of full relations. Each field contains
 * 	a list of relation numbers contributing to the 
 * 	corresponding full relation. In Cavallar's nomenclature, 'R'
 * 	will be a collection of relation sets.
 * 	'P' will be modified in the process.
 * 	Return value: The number of full relations built.
 **************************************************************/ 
s32 combParts_new(LSList& P, LSList& R, int maxRelsInFF, s32 minFF) {
	s32 wt0, wt1;
	s32 lastFull, full;
	int	pass=0;
	double shrink;

	/* 
	 * We can reallocate if necessary. So if this seems to be way
	 * too much memory, decrease it at will.
	*/
	printf("combParts() Doing ll_verify(P)...\n");
	printf("ll_verify() reports that 'P' appears to be intact.\n");

	for (u32 i=0; i<P.numFields(); ++i)
		R.appendField((s32*)&i, 1);

	full = -1; /* Force at least two passes. */
	do {
		printf("** pass %d...\n", ++pass);
		lastFull = full;
		full = makePass(R, P);
		checkR(R);
		printf("* There are now %" PRId32 " full relations.\n", full);
	} while (lastFull < full);

	/* Drop any relation-sets still containing a large prime: */
	keepFulls(R, P); 
	printf("After keepFulls(), R->numFields = %" PRId32 "\n", (s32)R.numFields());

	/* 
	 * Don't bother with the weight reduction unless we're close
	 * to having enough relations.
	 */
	if (full < minFF) {
		return full;
	}

	printf("Reducing the weight of relation sets. This is painfully\n");
	printf("slow at the moment, but it's worth it.\n");
	do {
		wt0 = R.Weight();
		reduceRelSets(R, P);
		wt1 = R.Weight();
		msgLog("", "reduceRelSets dropped relation-set weight from %" PRId32 " to %" PRId32 ".",
				wt0, wt1);
		shrink = (double)(wt0-wt1)/wt0;
	} while (shrink > 0.10);
	full = removeHeavyRelSets(R, P, maxRelsInFF);
	msgLog("", "After removing heavy rel-sets, weight is %" PRId32 ".", R.Weight());
	//printf("After removing heavy rel-sets, weight is %" PRId32 ".\n", R.Weight());

	return full;
}

bool checkRelSets(LSList&R, LSList& dataList) {
	checkR(R);
	for (u32 i=0; i<R.numFields(); ++i) {
		if (i%5000==0) {
			printTmp("Checking relation sets : [ %3.3d%% done]", 100*i/R.numFields());
		}
		LSListField* F = R.fetchField(i);
		s32 c=0;
		for (u32 j=0; j<F->numEntries();++j) {
			LSListField* cF = dataList.fetchField(F->at(j));
			for (u32 k=0; k<cF->numEntries();++k) {
				c^=cF->at(k);
			}
			dataList.releaseField(cF);
		}
		R.releaseField(F);
		if (c) {
			printf("Error: relation set no. %" PRId32 " appears to be corrupted!\n",i);
			return false;
		}
	}
	return true;
}

s32 combParts_tpie(llist_t *lR, llist_t *lP, int maxRelsInFF, s32 minFF) {
	MM_manager.ignore_memory_limit();

	LSList P(lP);
	LSList R(lR);
	LSList cmpP(lP);

	double t1=sTime();
	s32 numFulls=combParts_new(P,R,maxRelsInFF,minFF);
	if (checkRelSets(R,cmpP)) { 
		printf("checkRelSets() reports, that R has good chances to be fine ;)\n");
	} else {
		printf("checkRelSets() failed. Exiting...\n");
	}
	ll_clear(lP);
	P.toLList_t(lP);
	R.toLList_t(lR);

	double t2=sTime();
	msgLog("","Elapsed time for cycle-counting: %lf seconds.",t2-t1);
	cout	<< "Stream statistics (global):\n"
			<< "\tREAD ITEM:	"
			<< AMI_STREAM<s32>::gstats().get(ITEM_READ) << "\n"
			<< "\tWRITE ITEM:	 "
			<< AMI_STREAM<s32>::gstats().get(ITEM_WRITE) << "\n"
			<< "\tSEEK ITEM:	"
			<< AMI_STREAM<s32>::gstats().get(ITEM_SEEK) << "\n"
			<< "\tREAD BLOCK:	 "
			<< AMI_STREAM<s32>::gstats().get(BLOCK_READ) << "\n"
			<< "\tWRITE BLOCK:	"
			<< AMI_STREAM<s32>::gstats().get(BLOCK_WRITE) << "\n";
	return numFulls;
}
