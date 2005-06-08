/**************************************************************/
/* Classical sieve for the NFS.                               */
/**************************************************************/
/* Copyright 2005, Jason Papadopoulos                         */
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
#include "ggnfs.h"

extern int clForceStop;

/* sieving takes place in L1-cache-size blocks */
#define BLOCK_SIZE 65536
#define BLOCK_HASH(x) ((x) / BLOCK_SIZE)

/* the factor base updates to a given block are also
   explicitly listed out, to speedup up trial factoring
   and updating the sieve interval. Each sieve update
   is encoded into the following structure */
typedef struct {
  s32 p;       /* prime */
  u16 h;       /* offset within block */
  u8 logp;     /* log of prime */
} update_t;

/* because sieving takes place a block at a time, there
   is no need to tie up memory for blocks that have already
   been sieved. Hence the updates for a particular block
   are packed into a linked list of script_t structures,
   which are reallocated to other hash blocks later */
#define MAX_UPDATES 1000

typedef struct script_t {
  struct script_t *next;
  s32 num_used;
  double align_me;
  update_t list[MAX_UPDATES];
} script_t;

/* the blocks are indexed via a hashtable. Each block
   also maintains a cache of update structures that are
   packed together to conserve TLB resources */

typedef struct {
  script_t *script;  /* linked list for this block */
  u16 cache_off;     /* offset into the cache array where block cache starts */
  u16 cache_used;    /* number of entries currently used in block's cache */
} hashtable_t;

/* If the sieve values are expected to be small, an extra
   level of hashing is used that slows down sieving but massively
   speeds up trial factoring. */

#define SUBBLOCK_SIZE 2048
#define NUM_SUBBLOCKS (BLOCK_SIZE/SUBBLOCK_SIZE)
#define SUBBLOCK_HASH(x) ((x)/SUBBLOCK_SIZE)
#define SUBBLOCK_USE_THRESHOLD 15
#define SUBBLOCK_NUM_AVG_B 10

typedef struct {
  update_t *list;
  s32 num_alloc;
  s32 num_used;
} subblock_hash_t;

/* Primes below this value are not sieved */
#define TINY_R_CUTOFF 256
#define TINY_A_CUTOFF 256

/* how often the log targets for each sieve are
   recalculated */
#define A_NORM_INTERVAL 1024
#define R_NORM_INTERVAL 8192

/* fudge factor (base 2) for not sieving with small primes */
#define SMALL_PRIME_FF 30

#define POSITIVE 0
#define NEGATIVE 1

/* structure used for the rational and algebraic factor
   bases during sieving */
typedef struct {
  s32 p;        /* factor base prime */
  s32 r;        /* value where F(r)%p=0 (algebraic FB) or
                   (k1*r-k2)%p=0 (rational FB) */
  s32 root;     /* value where F(root,b)%p=0 (algebraic FB) or
                   (k1*root-k2*b)%p=0 (rational FB) */
  s32 next;     /* next sieve offset to update for this p */
  u8 logp;      /* value used to update the sieve */
} sieve_fb_t;
  
/* main structure controlling the sieving process */
typedef struct {
  FILE *outfp;

  script_t *mempool;

  s32 sieve_size;       /* size of (positive half of) sieve interval */
  s32 blocks;           /* number of blocks this is divided into */
  s32 block_window;     /* number of blocks covered by the hashtable */
  s32 cache_max;        /* maximum size of any block's cache */

  s32 num_tf_attempts[SUBBLOCK_NUM_AVG_B];
  s32 use_subblock;

  sieve_fb_t *rfb;            /* rational factor base */
  s32 rfb_size;               /* number of entries in RFB */
  s32 tiny_rfb_size;          /* number of entries in RFB */
  s32 med_rfb_size;           /* number of entries in RFB */
  u8 *rblock;                 /* piece of RFB (one block worth) */
  update_t *rcache;           /* concatenated array of all RFB caches */
  hashtable_t *rtable;        /* hashtable for RFB */
  s32 cutoff1_r;              /* fudge factor for sieve values */
  s32 cutoff2_r;              /* fudge factor for trial factoring */
  mpz_t LP1_max_r;            /* single large prime max bound */
  mpz_t LP2_min_r;            /* double large prime min bound */
  mpz_t LP2_max_r;            /* double large prime max bound */
  sieve_fb_t rfb_inf[64];     /* rational primes at infinity */
  s32 rfb_inf_size;           /* number of rational primes at infinity */
  u8 logTarget_r[2*BLOCK_SIZE/R_NORM_INTERVAL];
  subblock_hash_t subblock_r[NUM_SUBBLOCKS];

  /* the fields for the algebraic FB are completely
     equivalent to the above */
  sieve_fb_t *afb;
  s32 afb_size;
  s32 tiny_afb_size;
  s32 med_afb_size;
  u8 *ablock;
  update_t *acache;
  hashtable_t *atable;
  s32 cutoff1_a;
  s32 cutoff2_a;
  mpz_t LP1_max_a;
  mpz_t LP2_min_a;
  mpz_t LP2_max_a;
  sieve_fb_t afb_inf[64];
  s32 afb_inf_size;
  u8 logTarget_a[2*BLOCK_SIZE/A_NORM_INTERVAL];
  subblock_hash_t subblock_a[NUM_SUBBLOCKS];

} sieve_conf_t;
  
#if defined(__GNUC__)
#define PREFETCH(x) __builtin_prefetch(x)
#else
#define PREFETCH(x) /* nothing */
#endif

int doSieve(nfs_fb_t *FB, sieve_conf_t *conf, s32 b0, s32 b1);
int do1SieveEven(nfs_fb_t *FB, sieve_conf_t *conf, s32 b);
int do1SieveOdd(nfs_fb_t *FB, sieve_conf_t *conf, s32 b);
void fillHashtableEven(sieve_conf_t *conf, s32 b, s32 sign);
void fillHashtableOdd(sieve_conf_t *conf, s32 b, s32 sign);
int fillSieve(nfs_fb_t *FB, sieve_conf_t *conf, 
              s32 b, s32 block, s32 first_block, s32 sign);
int doFactoringOdd(nfs_fb_t *FB, sieve_conf_t *conf, 
                 s32 b, s32 block, s32 first_block, s32 sign);
int doFactoringEven(nfs_fb_t *FB, sieve_conf_t *conf, 
                 s32 b, s32 block, s32 first_block, s32 sign);
int do1Factoring(nfs_fb_t *FB, sieve_conf_t *conf, 
                 s32 off, s32 b, s32 rbits, s32 abits,
                 s32 block, s32 first_block, s32 sign);
void add2Hashtable(update_t *cache, s32 cache_max, 
                   hashtable_t *hash, s32 p, s32 h, u8 logp, 
                   script_t **mempool);
void add2Subblock(subblock_hash_t *entry, s32 p, s32 h);
void flushCache(update_t *cache, hashtable_t *hash, script_t **mempool);

void printRelation(sieve_conf_t *conf, char *buf);
void flushSavefile(sieve_conf_t *conf);
void clMakeOutputLine(char *str, relation_t *R);

/***********************************************************/     
int clSieve(nfs_sieve_job_t *J)
/***********************************************************
* Driver for classical sieve                               *
************************************************************/
{ s32   i, j, k;
  nfs_fb_t *FB=&J->FB;
  FILE *fp;
  sieve_conf_t conf;
  int rels;

  if (J->a0 != -(J->a1)) {
    fprintf(stderr, "fatal error: 'a' range of [%d,%d] "
            "not symmetric\n", (int)J->a0, (int)J->a1);
    return 0;
  }
  if (J->b1<J->b0 || J->b0<=0) {
    fprintf(stderr, "fatal error: invalid 'b' "
              "range of [%d,%d]\n", (int)J->b0, (int)J->b1);
    return 0;
  }
              
  if ((fp = fopen(J->outName, "a"))) {
    fprintf(fp, "# GGNFS version %s classical siever output.\n", GGNFS_VERSION);
    fprintf(fp, "# number: %s\n", FB->name);
    conf.outfp = fp;
  }

  setLogs(FB, J->a0, J->a1, J->b0, J->b1);
  FB->rfb_lambda = J->rfb_lambda;
  FB->afb_lambda = J->afb_lambda;

  printf("Doing classical sieve on a=[%d, %d] from b=%d to b=%d.\n",
          (int)J->a0, (int)J->a1, (int)J->b0, (int)J->b1);
  
  /**********************************
   * set up the rational FB         *
   **********************************/
  conf.rfb = (sieve_fb_t *)malloc(FB->rfb_size*sizeof(sieve_fb_t));

  for (i=j=k=0; i<FB->rfb_size; i++) {
    s32 p = FB->rfb[2*i];
    s32 r = FB->rfb[2*i+1];
    u8 logp = FB->rfb_log[i];
    s32 root;

    if (p==r) {
      conf.rfb_inf[j].p = p;
      conf.rfb_inf[j].logp = logp;
      j++;
      continue;
    }

    MULMOD32(root, J->b0, r, p);
    conf.rfb[k].p = p;
    conf.rfb[k].r = r;
    conf.rfb[k].root = root;
    conf.rfb[k].next = root;
    conf.rfb[k].logp = logp;
    k++;
  }
  conf.rfb_size = k;
  conf.rfb_inf_size = j;

  /*************************************************
   * find the offset of the max unsieved RFB prime *
   *************************************************/
  for (i=0; i<conf.rfb_size; i++) {
    if (conf.rfb[i].p>TINY_R_CUTOFF)
      break;
  }
  conf.tiny_rfb_size = i+1;

  for (i=0; i<conf.rfb_size; i++) {
    if (conf.rfb[i].p>BLOCK_SIZE)
      break;
  }
  conf.med_rfb_size = i+1;

  /************************************
   * fill in RFB large prime cutoffs  *
   ************************************/
  mpz_init_set_si(conf.LP1_max_r, FB->maxP_r);
  mpz_init_set_si(conf.LP2_max_r, FB->maxP_r);
  mpz_init_set_si(conf.LP2_min_r, FB->rLim);
  mpz_mul(conf.LP2_min_r, conf.LP2_min_r, conf.LP2_min_r);
  if (FB->maxLP==2) {
    if (FB->MFB_r>0)
      mpz_set_d(conf.LP2_max_r, pow(2.0, (double)FB->MFB_r));
    else
      mpz_mul_si(conf.LP2_max_r, conf.LP2_max_r, 
               (s32)(pow((double)FB->maxP_r, 0.95)+0.5));
  }

  /***************************************
   * fill in the trial factoring cutoffs *
   ***************************************/
  conf.cutoff1_r = (s32)(FB->rfb_lambda*fplog(FB->rLim, FB->log_rlb))+4;
  conf.cutoff2_r = fplog_mpz(conf.LP2_max_r, FB->log_rlb);
  conf.cutoff2_r = MAX(conf.cutoff1_r, conf.cutoff2_r);
  conf.cutoff1_r = conf.cutoff2_r + (s32)(SMALL_PRIME_FF*2/
                                          FB->rfb_log_base+0.5);
  
  /**********************************
   * set up the algebraic FB        *
   **********************************/
  conf.afb = (sieve_fb_t *)malloc(FB->afb_size * sizeof(sieve_fb_t));

  for (i=j=k=0; i<FB->afb_size; i++) {
    s32 p = FB->afb[2*i];
    s32 r = FB->afb[2*i+1];
    u8 logp = FB->afb_log[i];
    s32 root;

    if (p==r) {
      /* careful: AFB primes at infinity can occur
         more than once */
      if (j==0 || p!=conf.afb_inf[j-1].p) {
        conf.afb_inf[j].p = p;
        conf.afb_inf[j].logp = logp;
        j++;
      }
      continue;
    }

    MULMOD32(root, J->b0, r, p);
    conf.afb[k].p = p;
    conf.afb[k].r = r;
    conf.afb[k].root = root;
    conf.afb[k].next = root;
    conf.afb[k].logp = logp;
    k++;
  }
  conf.afb_size = k;
  conf.afb_inf_size = j;

  for (i=1,j=0; i<conf.afb_size; i++) {
    if (conf.afb[i].p!=conf.afb[i-1].p) {
      if (conf.afb[i].p>TINY_A_CUTOFF)
        break;
    }
  }
  conf.tiny_afb_size = i+1;

  for (i=1,j=0; i<conf.afb_size; i++) {
    if (conf.afb[i].p!=conf.afb[i-1].p) {
      if (conf.afb[i].p>BLOCK_SIZE)
        break;
    }
  }
  conf.med_afb_size = i+1;

  mpz_init_set_si(conf.LP1_max_a, FB->maxP_a);
  mpz_init_set_si(conf.LP2_max_a, FB->maxP_a);
  mpz_init_set_si(conf.LP2_min_a, FB->aLim);
  mpz_mul(conf.LP2_min_a, conf.LP2_min_a, conf.LP2_min_a);

  if (FB->maxLPA == 2) {
    if (FB->MFB_a>0)
      mpz_set_d(conf.LP2_max_a, pow(2.0, (double)FB->MFB_a));
    else
      mpz_mul_si(conf.LP2_max_a, conf.LP2_max_a, 
               (s32)(pow((double)FB->maxP_a, 0.95) + 0.5));
  }
  conf.cutoff1_a = (s32)(FB->afb_lambda*fplog(FB->aLim, FB->log_alb))+4;
  conf.cutoff2_a = fplog_mpz(conf.LP2_max_a, FB->log_alb);
  conf.cutoff2_a = MAX(conf.cutoff1_a, conf.cutoff2_a);
  conf.cutoff1_a = conf.cutoff2_a + (s32)(SMALL_PRIME_FF * 2 /
                                          FB->afb_log_base + 0.5);

  /************************
   * Allocate sieve stuff *
   ************************/
  conf.rblock = (u8 *)malloc(BLOCK_SIZE*sizeof(u8));
  conf.ablock = (u8 *)malloc(BLOCK_SIZE*sizeof(u8));
  conf.sieve_size = J->a1;
  conf.blocks = (J->a1+BLOCK_SIZE-1)/BLOCK_SIZE;

  /********************************************
   * The hashtable must contain enough blocks *
   * so that incrementing a value in block 0  *
   * by the largest possible prime will not   *
   * put that value in the last block         *
   ********************************************/

  conf.block_window = (MAX(FB->rLim, FB->aLim)+(BLOCK_SIZE-1))/
                                             BLOCK_SIZE + 1;
  conf.block_window = MIN(conf.blocks, conf.block_window);

  /********************************************
   * The cache for each block is large enough *
   * so that the hashtable and all the caches *
   * combined fit in BLOCK_SIZE bytes. These  *
   * two structures are packed contiguously so*
   * they won't conflict with each other. The *
   * cache arrays are each aligned to a 16-   *
   * byte boundary                            *
   ********************************************/
  conf.cache_max = (BLOCK_SIZE-conf.block_window*sizeof(hashtable_t))/
                                   (conf.block_window*sizeof(update_t));
  conf.cache_max = MIN(conf.cache_max, MAX_UPDATES);
  conf.rtable = (hashtable_t *)malloc(BLOCK_SIZE+16);
  conf.atable = (hashtable_t *)malloc(BLOCK_SIZE+16);
  conf.rcache = (update_t *)(conf.rtable+conf.block_window);
  conf.rcache = (update_t *)((u8 *)(conf.rcache)+16-((s32)(conf.rcache)%16));
  conf.acache = (update_t *)(conf.atable+conf.block_window);
  conf.acache = (update_t *)((u8 *)(conf.acache)+16-((s32)(conf.acache)%16));
  conf.mempool = NULL;

  for (i=0; i<conf.block_window; i++) {
    conf.rtable[i].script = NULL;
    conf.rtable[i].cache_off = i*conf.cache_max;
    conf.rtable[i].cache_used = 0;
    conf.atable[i].script = NULL;
    conf.atable[i].cache_off = i*conf.cache_max;
    conf.atable[i].cache_used = 0;
  }

  /*******************************************
   * prepare to allcate a second level of    *
   * hashtables, for use when a lot of trial *
   * factoring is expected. By default the   *
   * second hashtable starts off unused,     *
   * unless the first b value is small.      *
   * The actual sieving code will turn       *
   * this feature on and off as conditions   *
   * warrant it.                             *
   *******************************************/
  memset(conf.subblock_r, 0, sizeof(conf.subblock_r));
  memset(conf.subblock_a, 0, sizeof(conf.subblock_a));
  memset(conf.num_tf_attempts, 0, sizeof(conf.num_tf_attempts));
  conf.use_subblock = 0;
  if (J->b0<=10000) {
    conf.use_subblock = 1;
    for (i=0; i<NUM_SUBBLOCKS; i++) {
      conf.subblock_r[i].list = (update_t *)
                          malloc(1000*sizeof(update_t));
      conf.subblock_r[i].num_alloc = 1000;
      conf.subblock_a[i].list = (update_t *)
                          malloc(1000*sizeof(update_t));
      conf.subblock_a[i].num_alloc = 1000;
    }
  }

  msgLog("", "");
  msgLog("", "hashtable: %d bins of size %d", conf.block_window, BLOCK_SIZE);
  msgLog("", "hashtable cache: %d entries per bin", conf.cache_max);
  msgLog("", "Rational factor base:");
  msgLog("", "base of logs: %4.3f", FB->rfb_log_base);
  msgLog("", "factor base entries: %d (%3.1f MB)", conf.rfb_size,
              (double)conf.rfb_size*sizeof(sieve_fb_t)/1e6);
  msgLog("", "maximum factor base prime: %d", FB->rLim);
  msgLog("", "primes at infinity: %d", conf.rfb_inf_size);
  msgLog("", "hashed RFB entries: %d (%3.1f%%, max=%d)", 
               conf.rfb_size-conf.med_rfb_size, 
               100.0*(conf.rfb_size-conf.med_rfb_size)/conf.rfb_size,
               conf.rfb[conf.rfb_size-1].p);
  msgLog("", "sieved RFB entries: %d (%3.3f%%, max=%d)", 
               conf.med_rfb_size-conf.tiny_rfb_size, 
               100.0*(conf.med_rfb_size-conf.tiny_rfb_size)/conf.rfb_size,
               conf.rfb[conf.med_rfb_size-1].p);
  msgLog("", "unsieved RFB entries: %d", conf.tiny_rfb_size);
  msgLog("", "large prime cutoff: %d bits", mpz_sizeinbase(conf.LP1_max_r,2)-1);
  msgLog("", "trial factoring cutoff: %d bits", 
               (s32)(conf.cutoff2_r*M_LOG2E*FB->log_rlb)-1);
  msgLog("", "2-large prime cutoff: %d-%d bits", 
              mpz_sizeinbase(conf.LP2_min_r,2)-1,
              mpz_sizeinbase(conf.LP2_max_r,2)-1);

  msgLog("", "Algebraic factor base:");
  msgLog("", "base of logs: %4.3f", FB->afb_log_base);
  msgLog("", "factor base entries: %d (%3.1f MB)", conf.afb_size,
              (double)conf.afb_size*sizeof(sieve_fb_t)/1e6);
  msgLog("", "maximum factor base prime: %d", FB->aLim);
  msgLog("", "primes at infinity: %d", conf.afb_inf_size);
  msgLog("", "hashed AFB entries: %d (%3.1f%%, max=%d)", 
               conf.afb_size-conf.med_afb_size, 
               100.0*(conf.afb_size-conf.med_afb_size)/conf.afb_size,
               conf.afb[conf.afb_size-1].p);
  msgLog("", "sieved AFB entries: %d (%3.3f%%, max=%d)", 
               conf.med_afb_size-conf.tiny_afb_size, 
               100.0*(conf.med_afb_size-conf.tiny_afb_size)/conf.afb_size,
               conf.afb[conf.med_afb_size-1].p);
  msgLog("", "unsieved AFB entries: %d", conf.tiny_afb_size);
  msgLog("", "large prime cutoff: %d bits", mpz_sizeinbase(conf.LP1_max_a,2)-1);
  msgLog("", "trial factoring cutoff: %d bits", 
               (s32)(conf.cutoff2_a*M_LOG2E*FB->log_alb)-1);
  msgLog("", "2-large prime cutoff: %d-%d bits", 
              mpz_sizeinbase(conf.LP2_min_a,2)-1,
              mpz_sizeinbase(conf.LP2_max_a,2)-1);
  msgLog("", "");

  /**********************
   * do all the work :) *
   **********************/
  rels = doSieve(FB, &conf, J->b0, J->b1);

  /************
   * clean up *
   ************/
  for (i=0; i<NUM_SUBBLOCKS; i++) {
    free(conf.subblock_r[i].list);
    free(conf.subblock_a[i].list);
  }
  while (conf.mempool!=NULL) {
    script_t *next = conf.mempool->next;
    free(conf.mempool);
    conf.mempool = next;
  }
  free(conf.rtable);
  free(conf.atable);
  free(conf.rblock);
  free(conf.ablock);
  free(conf.rfb);
  free(conf.afb);
  fclose(conf.outfp);
  return rels;
}

/***********************************************************/     
int doSieve(nfs_fb_t *FB, sieve_conf_t *conf, s32 b0, s32 b1)
/***********************************************************     
 * Sieve lines for b = b0 to b1                            *
 ***********************************************************/
{ s32   i, j, k;
  s32 rfb_size = conf->rfb_size;
  s32 afb_size = conf->afb_size;
  sieve_fb_t *rfb = conf->rfb;
  sieve_fb_t *afb = conf->afb;
  int old_rels = 0, rels = 0;
  double startTime = sTime();
  double currTime;

  fflush(stdout);

  for (i=b0; i<=b1; i++) {

    for (j=0; j<SUBBLOCK_NUM_AVG_B-1; j++)
      conf->num_tf_attempts[j] = conf->num_tf_attempts[j+1];
    conf->num_tf_attempts[j] = 0;

    /************************************
     * do the sieving for the current b *
     ************************************/
    if (i%2)
      rels += do1SieveOdd(FB, conf, i);
    else
      rels += do1SieveEven(FB, conf, i);

    /**************************************
     * update the roots for each FB prime *
     **************************************/
    for (j=0; j<rfb_size; j++) {
      s32 p = rfb[j].p;
      s32 root = rfb[j].root + rfb[j].r;
      if (root>=p)
        root-=p;
      rfb[j].root = root;
      rfb[j].next = root;
    }
    for (j=0; j<afb_size; j++) {
      s32 p = afb[j].p;
      s32 root = afb[j].root + afb[j].r;
      if (root>=p)
        root-=p;
      afb[j].root = root;
      afb[j].next = root;
    }

    /*******************************************
     * Determine if the next line should slow  *
     * down the sieving stage to make trial    *
     * factoring much faster. The indication   *
     * that this is necessary comes from the   *
     * number of trial factoring attempts per  *
     * hash block, averaged over a few lines   *
     *******************************************/

    for (j=k=0; j<SUBBLOCK_NUM_AVG_B; j++)
      k += conf->num_tf_attempts[j];

    conf->use_subblock = 0;
    if (k>SUBBLOCK_USE_THRESHOLD*SUBBLOCK_NUM_AVG_B*conf->blocks) {
      if (conf->subblock_r[0].list==NULL) {
        for (j=0; j<NUM_SUBBLOCKS; j++) {
          conf->subblock_r[j].list = (update_t *)
                              malloc(1000*sizeof(update_t));
          conf->subblock_r[j].num_alloc = 1000;
          conf->subblock_a[j].list = (update_t *)
                              malloc(1000*sizeof(update_t));
          conf->subblock_a[j].num_alloc = 1000;
        }
      }
      conf->use_subblock = 1;
    }

    /*******************************************
     * print a progress message if appropriate *
     *******************************************/
    if (rels>old_rels+20 || clForceStop || i==b1) {
      currTime = sTime();
      printf("TR=%d b=%d newR=%d (%3.3f sec/rel)\n", (int)rels,
               (int)i, (int)(rels-old_rels), (currTime-startTime)/rels);
      fflush(stdout);
      old_rels = rels;
    }
    if (clForceStop) {
      i++; /* compensate for the i-1 below in the log message. */
      break;
    }
  }
  msgLog("", "Classical sieved [%ld, %ld]x[%ld, %ld]",
              -conf->sieve_size, conf->sieve_size, b0, i-1);
  flushSavefile(conf);
  return rels;
}

/***********************************************************/     
int do1SieveOdd(nfs_fb_t *FB, sieve_conf_t *conf, s32 b)
/**************************************
 * Do the sieving for odd values of b *
 **************************************/
{ s32   i, j;
  int rels = 0;

  fillHashtableOdd(conf, b, POSITIVE);
  for (i=0; i<conf->blocks; i+=conf->block_window) {
    s32 num_blocks = MIN(conf->block_window, conf->blocks-i);
    for (j=0; j<num_blocks; j++) {
      rels += fillSieve(FB, conf, b, j, i, POSITIVE);
    }
  }

  for (i=conf->tiny_rfb_size; i<conf->med_rfb_size; i++)
    conf->rfb[i].next = conf->rfb[i].p - conf->rfb[i].root;
  for (i=conf->tiny_afb_size; i<conf->med_afb_size; i++)
    conf->afb[i].next = conf->afb[i].p - conf->afb[i].root;

  fillHashtableOdd(conf, b, NEGATIVE);
  for (i=0; i<conf->blocks; i+=conf->block_window) {
    s32 num_blocks = MIN(conf->block_window, conf->blocks-i);
    for (j=0; j<num_blocks; j++) {
      rels += fillSieve(FB, conf, b, j, i, NEGATIVE);
    }
  }

  return rels;
}

/***********************************************************/     
int do1SieveEven(nfs_fb_t *FB, sieve_conf_t *conf, s32 b)
/***************************************
 * Do the sieving for even values of b *
 ***************************************/
{ s32   i, j;
  int rels = 0;

  for (i=conf->tiny_rfb_size; i<conf->med_rfb_size; i++) {
    sieve_fb_t *fb = conf->rfb + i;
    s32 p = fb->p;
    s32 next = fb->next;
    if (next%2==0)
      next += p;
    fb->next = next/2;
  }
  for (i=conf->tiny_afb_size; i<conf->med_afb_size; i++) {
    sieve_fb_t *fb = conf->afb + i;
    s32 p = fb->p;
    s32 next = fb->next;
    if (next%2==0)
      next += p;
    fb->next = next/2;
  }

  fillHashtableEven(conf, b, POSITIVE);
  for (i=0; i<(conf->blocks+1)/2; i+=conf->block_window) {
    s32 num_blocks = MIN(conf->block_window, (conf->blocks+1)/2-i);
    for (j=0; j<num_blocks; j++) {
      rels += fillSieve(FB, conf, b, j, i, POSITIVE);
    }
  }

  for (i=conf->tiny_rfb_size; i<conf->med_rfb_size; i++) {
    sieve_fb_t *fb = conf->rfb + i;
    s32 p = fb->p;
    s32 next = p - fb->root;
    if (next%2==0)
      next += p;
    fb->next = next/2;
  }
  for (i=conf->tiny_afb_size; i<conf->med_afb_size; i++) {
    sieve_fb_t *fb = conf->afb + i;
    s32 p = fb->p;
    s32 next = p - fb->root;
    if (next%2==0)
      next += p;
    fb->next = next/2;
  }

  fillHashtableEven(conf, b, NEGATIVE);
  for (i=0; i<(conf->blocks+1)/2; i+=conf->block_window) {
    s32 num_blocks = MIN(conf->block_window, (conf->blocks+1)/2-i);
    for (j=0; j<num_blocks; j++) {
      rels += fillSieve(FB, conf, b, j, i, NEGATIVE);
    }
  }

  return rels;
}


/***********************************************************/     
void fillHashtableOdd(sieve_conf_t *conf, s32 b, s32 sign)
{ s32 i;
  s32 sieve_size = conf->sieve_size;
  s32 rfb_size = conf->rfb_size;
  s32 afb_size = conf->afb_size;
  s32 med_rfb_size = conf->med_rfb_size;
  s32 med_afb_size = conf->med_afb_size;
  sieve_fb_t *rfb = conf->rfb;
  sieve_fb_t *afb = conf->afb;
  hashtable_t *rhash = conf->rtable;
  hashtable_t *ahash = conf->atable;

  /*************************************
   * Note that only the first sieve    *
   * update for each factor base prime *
   * is put into the hashtable         *
   *************************************/
  if (sign==POSITIVE) {
    for (i=med_rfb_size; i<rfb_size; i++) {
      sieve_fb_t *fb = rfb + i;
      s32 p = fb->p;
      s32 off = fb->next;
      u8 logp = fb->logp;

      if (off<sieve_size)
        add2Hashtable(conf->rcache, conf->cache_max,
                      rhash+BLOCK_HASH(off), 
                      p, off, logp, &conf->mempool);
    }
    for (i=med_afb_size; i<afb_size; i++) {
      sieve_fb_t *fb = afb + i;
      s32 p = fb->p;
      s32 off = fb->next;
      u8 logp = fb->logp;
      if (off<sieve_size)
        add2Hashtable(conf->acache, conf->cache_max,
                      ahash+BLOCK_HASH(off), 
                      p, off, logp, &conf->mempool);
    }
  }
  else {
    for (i=med_rfb_size; i<rfb_size; i++) {
      sieve_fb_t *fb = rfb + i;
      s32 p = fb->p;
      s32 off = p - fb->next;
      u8 logp = fb->logp;
      if (off<sieve_size)
        add2Hashtable(conf->rcache, conf->cache_max,
                      rhash+BLOCK_HASH(off), 
                      p, off, logp, &conf->mempool);
    }
    for (i=med_afb_size; i<afb_size; i++) {
      sieve_fb_t *fb = afb + i;
      s32 p = fb->p;
      s32 off = p - fb->next;
      u8 logp = fb->logp;
      if (off<sieve_size)
        add2Hashtable(conf->acache, conf->cache_max,
                      ahash+BLOCK_HASH(off), 
                      p, off, logp, &conf->mempool);
    }
  }
}

/***********************************************************/     
void fillHashtableEven(sieve_conf_t *conf, s32 b, s32 sign)
{ s32 i;
  s32 sieve_size = conf->sieve_size/2;
  s32 rfb_size = conf->rfb_size;
  s32 afb_size = conf->afb_size;
  s32 med_rfb_size = conf->med_rfb_size;
  s32 med_afb_size = conf->med_afb_size;
  sieve_fb_t *rfb = conf->rfb;
  sieve_fb_t *afb = conf->afb;
  hashtable_t *rhash = conf->rtable;
  hashtable_t *ahash = conf->atable;

  /********************************************
   * for even values of b, the sieve interval *
   * is compressed, so that updates are by    *
   * p and not 2p                             *
   ********************************************/
  if (sign==POSITIVE) {
    for (i=med_rfb_size; i<rfb_size; i++) {
      sieve_fb_t *fb = rfb + i;
      s32 p = fb->p;
      s32 off = fb->next;
      u8 logp = fb->logp;
      if (off%2==0)
        off += p;
      off = off/2;
      if (off<sieve_size)
        add2Hashtable(conf->rcache, conf->cache_max,
                      rhash+BLOCK_HASH(off), 
                      p, off, logp, &conf->mempool);
    }
    for (i=med_afb_size; i<afb_size; i++) {
      sieve_fb_t *fb = afb + i;
      s32 p = fb->p;
      s32 off = fb->next;
      u8 logp = fb->logp;
      if (off%2==0)
        off += p;
      off = off/2;
      if (off<sieve_size)
        add2Hashtable(conf->acache, conf->cache_max,
                      ahash+BLOCK_HASH(off), 
                      p, off, logp, &conf->mempool);
    }
  }
  else {
    for (i=med_rfb_size; i<rfb_size; i++) {
      sieve_fb_t *fb = rfb + i;
      s32 p = fb->p;
      s32 off = p - fb->next;
      u8 logp = fb->logp;
      if (off%2==0)
        off += p;
      off = off/2;
      if (off<sieve_size)
        add2Hashtable(conf->rcache, conf->cache_max,
                      rhash+BLOCK_HASH(off), 
                      p, off, logp, &conf->mempool);
    }
    for (i=med_afb_size; i<afb_size; i++) {
      sieve_fb_t *fb = afb + i;
      s32 p = fb->p;
      s32 off = p - fb->next;
      u8 logp = fb->logp;
      if (off%2==0)
        off += p;
      off = off/2;
      if (off<sieve_size)
        add2Hashtable(conf->acache, conf->cache_max,
                      ahash+BLOCK_HASH(off), 
                      p, off, logp, &conf->mempool);
    }
  }
}

/***********************************************************/     
int fillSieve(nfs_fb_t *FB, sieve_conf_t *conf, 
              s32 b, s32 block, s32 first_block, s32 sign)
/***********************************************************
 * Combine all the hashtable data with the sieving for the *
 * small factor base primes, then look for smooth relations*
 ***********************************************************/
{ s32 i, j, k;
  s32 tiny_rfb_size = conf->tiny_rfb_size;
  s32 med_rfb_size = conf->med_rfb_size;
  s32 tiny_afb_size = conf->tiny_afb_size;
  s32 med_afb_size = conf->med_afb_size;
  s32 block_start = (first_block+block)*BLOCK_SIZE;
  s32 last_block = conf->block_window;
  s32 sieve_size = (b&1)? conf->sieve_size : conf->sieve_size/2;
  hashtable_t *rhash = conf->rtable;
  hashtable_t *ahash = conf->atable;
  subblock_hash_t *subblock_r = conf->subblock_r;
  subblock_hash_t *subblock_a = conf->subblock_a;
  u8 *rblock = conf->rblock;
  u8 *ablock = conf->ablock;
  script_t *script;
  int rels = 0;

  /**********************************
   * clear the sieve block, pulling *
   * it into cache                  *
   **********************************/
  memset(rblock, 0, BLOCK_SIZE);

  if (!conf->use_subblock) {
    /*******************************************
     * sieve with the small factor base primes *
     *******************************************/
    for (i=tiny_rfb_size; i<med_rfb_size; i++) {
      sieve_fb_t *fb = conf->rfb + i;
      s32 p = fb->p;
      u8 logp = fb->logp;
      s32 off = fb->next;
  
      while (off<BLOCK_SIZE) {
        rblock[off] += logp;
        off += p;
      }
      fb->next = off-BLOCK_SIZE;
    }
  
    /**************************************
     * Also empty the block hash bin into *
     * the sieve interval. Since the block*
     * size is small, most updates will be*
     * cache hits and so will be fast     *
     **************************************/
    flushCache(conf->rcache, rhash+block, &conf->mempool);
    script = rhash[block].script;
    while (script!=NULL) {
      update_t *updates = script->list;
      k = script->num_used;
      for (j=0; j<k; j++) {
        s32 h = updates[j].h;
        u8 logp = updates[j].logp;
        if (j%8==0)
          PREFETCH(updates+j+16);
        rblock[h] += logp;
      }
      script = script->next;
    }
  }
  else {

    /*************************************
     * When b is small, add all of the   *
     * required updates into a hashtable *
     * whose block size is much smaller  *
     *************************************/
    for (i=0; i<NUM_SUBBLOCKS; i++)
      subblock_r[i].num_used = 0;

    for (i=tiny_rfb_size; i<med_rfb_size; i++) {
      sieve_fb_t *fb = conf->rfb + i;
      s32 p = fb->p;
      u8 logp = fb->logp;
      s32 off = fb->next;
  
      while (off<BLOCK_SIZE) {
        rblock[off] += logp;
        add2Subblock(subblock_r+SUBBLOCK_HASH(off), p, off);
        off += p;
      }
      fb->next = off-BLOCK_SIZE;
    }
  
    flushCache(conf->rcache, rhash+block, &conf->mempool);
    script = rhash[block].script;
    while (script!=NULL) {
      update_t *updates = script->list;
      k = script->num_used;
      for (j=0; j<k; j++) {
        s32 p = updates[j].p;
        s32 h = updates[j].h;
        u8 logp = updates[j].logp;
        if (j%8==0)
          PREFETCH(updates+j+16);
        rblock[h] += logp;
        add2Subblock(subblock_r+SUBBLOCK_HASH(h), p, h);
      }
      script = script->next;
    }
  }

  /****************************************
   * repeat for the algebraic factor base *
   ****************************************/
  memset(ablock, 0, BLOCK_SIZE);

  if (!conf->use_subblock) {
    for (i=tiny_afb_size; i<med_afb_size; i++) {
      sieve_fb_t *fb = conf->afb + i;
      s32 p = fb->p;
      u8 logp = fb->logp;
      s32 off = fb->next;
  
      while (off<BLOCK_SIZE) {
        ablock[off] += logp;
        off += p;
      }
      fb->next = off-BLOCK_SIZE;
    }
  
    flushCache(conf->acache, ahash+block, &conf->mempool);
    script = ahash[block].script;
    while (script!=NULL) {
      update_t *updates = script->list;
      k = script->num_used;
      for (j=0; j<k; j++) {
        s32 h = updates[j].h;
        u8 logp = updates[j].logp;
        if (j%8==0)
          PREFETCH(updates+j+16);
        ablock[h] += logp;
      }
      script = script->next;
    }
  }
  else {
    for (i=0; i<NUM_SUBBLOCKS; i++)
      subblock_a[i].num_used = 0;

    for (i=tiny_afb_size; i<med_afb_size; i++) {
      sieve_fb_t *fb = conf->afb + i;
      s32 p = fb->p;
      u8 logp = fb->logp;
      s32 off = fb->next;
  
      while (off<BLOCK_SIZE) {
        ablock[off] += logp;
        add2Subblock(subblock_a+SUBBLOCK_HASH(off), p, off);
        off += p;
      }
      fb->next = off-BLOCK_SIZE;
    }
  
    flushCache(conf->acache, ahash+block, &conf->mempool);
    script = ahash[block].script;
    while (script!=NULL) {
      update_t *updates = script->list;
      k = script->num_used;
      for (j=0; j<k; j++) {
        s32 p = updates[j].p;
        s32 h = updates[j].h;
        u8 logp = updates[j].logp;
        if (j%8==0)
          PREFETCH(updates+j+16);
        ablock[h] += logp;
        add2Subblock(subblock_a+SUBBLOCK_HASH(h), p, h);
      }
      script = script->next;
    }
  }

  /********************************
   * look for smooth sieve values *
   ********************************/
  if (b%2)
    rels += doFactoringOdd(FB, conf, b, block, first_block, sign);
  else
    rels += doFactoringEven(FB, conf, b, block, first_block, sign);

  /****************************************
   * The hash bins for the current block  *
   * are not needed any longer; increment *
   * each of the entries using their      *
   * associated prime, add the entry      *
   * back into some other hash bin, and   *
   * when the entire script is finished   *
   * return its memory to the pool for    *
   * other hash bins to use               *
   ****************************************/
  script = rhash[block].script;
  while (script!=NULL) {
    update_t *updates = script->list;
    script_t *next_script = script->next;
    k = script->num_used;
    for (j=0; j<k; j++) {
      s32 p = updates[j].p;
      u32 h = block_start + updates[j].h + p;
      u8 logp = updates[j].logp;
      /**************************************
       * pick the next hash bin that entry  *
       * j will go into. Hash bins are al-  *
       * located in a circular manner; it   *
       * doesn't matter which bin is chosen,*
       * as long as the all the updated     *
       * factor base primes that need to be *
       * together end up there              *
       **************************************/
      s32 new_block = BLOCK_HASH(h)-first_block;
      if (new_block>=last_block)
        new_block -= last_block;
      if (j%8==0)
        PREFETCH(updates+j+16);
      if (h < sieve_size)
        add2Hashtable(conf->rcache, conf->cache_max,
                      rhash+new_block,
                      p, h, logp, &conf->mempool);
    }
    script->next = conf->mempool;
    conf->mempool = script;
    script = next_script;
  }
  rhash[block].script = NULL;

  script = ahash[block].script;
  while (script!=NULL) {
    update_t *updates = script->list;
    script_t *next_script = script->next;
    k = script->num_used;
    for (j=0; j<k; j++) {
      s32 p = updates[j].p;
      u32 h = block_start + updates[j].h + p;
      u8 logp = updates[j].logp;
      s32 new_block = BLOCK_HASH(h)-first_block;
      if (new_block>=last_block)
        new_block -= last_block;
      if (j%8==0)
        PREFETCH(updates+j+16);
      if (h < sieve_size)
        add2Hashtable(conf->acache, conf->cache_max,
                      ahash+new_block,
                      p, h, logp, &conf->mempool);
    }
    script->next = conf->mempool;
    conf->mempool = script;
    script = next_script;
  }
  ahash[block].script = NULL;

  return rels;
}

/***********************************************************/     
int doFactoringOdd(nfs_fb_t *FB, sieve_conf_t *conf, 
                 s32 b, s32 block, s32 first_block, s32 sign)
/***********************************************************
 * Scan one sieve block looking for values to trial factor *
 ***********************************************************/
{ s32 i, j;
  u8 logTarget_r, logTarget_a;
  u8 *rblock = conf->rblock;
  u8 *ablock = conf->ablock;
  s32 blocksize = MIN(R_NORM_INTERVAL, A_NORM_INTERVAL);
  s32 cutoff1_r, cutoff1L_r, cutoff1R_r;
  s32 cutoff1_a, cutoff1L_a, cutoff1R_a;
  s32 block_start = (first_block+block)*BLOCK_SIZE;
  int rels = 0;

  static s32 initialized = 0;
  static mpz_t tmp;
  
  if (!initialized) {
    initialized = 1;
    mpz_init(tmp);
  }

  /**************************************
   * begin the progression of norm      *
   * sizes for use in the sieve interval*
   **************************************/
  if (sign==POSITIVE) {
    cutoff1L_a = fplog_evalF(block_start, b, FB);
    mpz_mul_si(tmp, FB->y1, block_start);
    mpz_addmul_ui(tmp, FB->y0, b);
    cutoff1L_r = fplog_mpz(tmp, FB->log_rlb);
  }
  else {
    cutoff1L_a = fplog_evalF(-block_start, b, FB);
    mpz_mul_si(tmp, FB->y1, -block_start);
    mpz_addmul_ui(tmp, FB->y0, b);
    cutoff1L_r = fplog_mpz(tmp, FB->log_rlb);
  }

  /******************************************
   * Work out the size of the rational      *
   * norms in block i. This is a little     *
   * wasteful if the rational poly is monic *
   ******************************************/
  for (j=0; j<BLOCK_SIZE; j+=R_NORM_INTERVAL) {
    if (sign==POSITIVE) {
      mpz_mul_si(tmp, FB->y1, block_start+j+R_NORM_INTERVAL);
      mpz_addmul_ui(tmp, FB->y0, b);
      cutoff1R_r = fplog_mpz(tmp, FB->log_rlb);
    }
    else {
      mpz_mul_si(tmp, FB->y1, -(block_start+j+R_NORM_INTERVAL));
      mpz_addmul_ui(tmp, FB->y0, b);
      cutoff1R_r = fplog_mpz(tmp, FB->log_rlb);
    }

    cutoff1_r = (cutoff1L_r+cutoff1R_r)/2 - conf->cutoff1_r;
    if (cutoff1_r<0)
      cutoff1_r = 0;
    conf->logTarget_r[j/R_NORM_INTERVAL] = cutoff1_r;
    cutoff1L_r = cutoff1R_r;
  }

  /******************************************************
   * Repeat the procedure for the algebraic factor base * 
   ******************************************************/

  for (j=0; j<BLOCK_SIZE; j+=A_NORM_INTERVAL) {
    if (sign==POSITIVE)
      cutoff1R_a = fplog_evalF(block_start+j+A_NORM_INTERVAL, b, FB);
    else
      cutoff1R_a = fplog_evalF(-(block_start+j+A_NORM_INTERVAL), b, FB);

    cutoff1_a = (cutoff1L_a+cutoff1R_a)/2 - conf->cutoff1_a;
    if (cutoff1_a<0)
      cutoff1_a = 0;
    conf->logTarget_a[j/A_NORM_INTERVAL] = cutoff1_a;
    cutoff1L_a = cutoff1R_a;
  }

  /*************************************************
   * Note that it would be much better to compare  *
   * blocks of bytes using some kind of multimedia *
   * packed unsigned compare, but the only way to  *
   * do this portably would be to decrement sieve  *
   * values instead of incrementing them, then to  *
   * test in parallel for sign bits that are set.  *
   * Unfortunately, this reduces the dynamic range *
   * that's allowed for log values and makes the   *
   * process slower overall.                       *
   *************************************************/
  for (i=0; i<BLOCK_SIZE; i+=blocksize) {
    logTarget_r = conf->logTarget_r[i/R_NORM_INTERVAL];
    logTarget_a = conf->logTarget_a[i/A_NORM_INTERVAL];

    for (j=0; j<blocksize; j++) {
      s32 off = i+j;
      s32 rbits = rblock[off];
      s32 abits = ablock[off];

      if (abits>logTarget_a &&
          rbits>logTarget_r &&
          gcd(block_start+off, b)==1) {
        rels += do1Factoring(FB, conf, off, b, 
                            rbits, abits, block, 
                            first_block, sign);
      }
    }
  }
  return rels;
}

/***********************************************************/     
int doFactoringEven(nfs_fb_t *FB, sieve_conf_t *conf, 
                 s32 b, s32 block, s32 first_block, s32 sign)
/***********************************************************
 * Scan one sieve block looking for values to trial factor *
 * (version specialized for even b)                        *
 ***********************************************************/
{ s32 i, j;
  u8 logTarget_r, logTarget_a;
  u8 *rblock = conf->rblock;
  u8 *ablock = conf->ablock;
  s32 blocksize = MIN(R_NORM_INTERVAL/2, A_NORM_INTERVAL/2);
  s32 cutoff1_r, cutoff1L_r, cutoff1R_r;
  s32 cutoff1_a, cutoff1L_a, cutoff1R_a;
  s32 block_start = (first_block+block)*BLOCK_SIZE;
  int rels = 0;

  /*************************************
   * the sieve interval for even b is  *
   * compressed, so that BLOCK_SIZE    *
   * entries count for 2*BLOCK_SIZE    *
   * sieve values (half of which are   *
   * ignored). That means array offset *
   * x corresponds to sieve offset     *
   * (2*x+1)                           *
   *************************************/
  static s32 initialized = 0;
  static mpz_t tmp;
  
  if (!initialized) {
    initialized = 1;
    mpz_init(tmp);
  }

  if (sign==POSITIVE) {
    cutoff1L_a = fplog_evalF(2*block_start, b, FB);
    mpz_mul_si(tmp, FB->y1, 2*block_start);
    mpz_addmul_ui(tmp, FB->y0, b);
    cutoff1L_r = fplog_mpz(tmp, FB->log_rlb);
  }
  else {
    cutoff1L_a = fplog_evalF(-2*block_start, b, FB);
    mpz_mul_si(tmp, FB->y1, -2*block_start);
    mpz_addmul_ui(tmp, FB->y0, b);
    cutoff1L_r = fplog_mpz(tmp, FB->log_rlb);
  }

  for (j=0; j<BLOCK_SIZE; j+=R_NORM_INTERVAL/2) {
    if (sign==POSITIVE) {
      mpz_mul_si(tmp, FB->y1, 2*(block_start+j+R_NORM_INTERVAL/2));
      mpz_addmul_ui(tmp, FB->y0, b);
      cutoff1R_r = fplog_mpz(tmp, FB->log_rlb);
    }
    else {
      mpz_mul_si(tmp, FB->y1, -2*(block_start+j+R_NORM_INTERVAL/2));
      mpz_addmul_ui(tmp, FB->y0, b);
      cutoff1R_r = fplog_mpz(tmp, FB->log_rlb);
    }

    cutoff1_r = (cutoff1L_r+cutoff1R_r)/2 - conf->cutoff1_r;
    if (cutoff1_r<0)
      cutoff1_r = 0;
    conf->logTarget_r[j/(R_NORM_INTERVAL/2)] = cutoff1_r;
    cutoff1L_r = cutoff1R_r;
  }

  for (j=0; j<BLOCK_SIZE; j+=A_NORM_INTERVAL/2) {
    if (sign==POSITIVE)
      cutoff1R_a = fplog_evalF(2*(block_start+j+A_NORM_INTERVAL/2), b, FB);
    else
      cutoff1R_a = fplog_evalF(-2*(block_start+j+A_NORM_INTERVAL/2), b, FB);

    cutoff1_a = (cutoff1L_a+cutoff1R_a)/2 - conf->cutoff1_a;
    if (cutoff1_a<0)
      cutoff1_a = 0;
    conf->logTarget_a[j/(A_NORM_INTERVAL/2)] = cutoff1_a;
    cutoff1L_a = cutoff1R_a;
  }

  for (i=0; i<BLOCK_SIZE; i+=blocksize) {
    logTarget_r = conf->logTarget_r[i/(R_NORM_INTERVAL/2)];
    logTarget_a = conf->logTarget_a[i/(A_NORM_INTERVAL/2)];

    for (j=0; j<blocksize; j++) {
      s32 off = i+j;
      s32 rbits = rblock[off];
      s32 abits = ablock[off];

      if (abits>logTarget_a &&
          rbits>logTarget_r &&
          gcd(2*(block_start+off)+1, b)==1) {
        rels += do1Factoring(FB, conf, off, b, 
                            rbits, abits, block, 
                            first_block, sign);
      }
    }
  }
  return rels;
}

/***********************************************************/     
int do1Factoring(nfs_fb_t *FB, sieve_conf_t *conf, 
                 s32 off, s32 b, s32 rbits, s32 abits,
                 s32 block, s32 first_block, s32 sign)
/***********************************************************
 * Do trial factoring to test a single sieve value for     *
 * smoothness                                              *
 ***********************************************************/
{ s32 i, j, k, rem = 0;
  hashtable_t *rhash = conf->rtable+block;
  hashtable_t *ahash = conf->atable+block;
  subblock_hash_t *subblock_r = conf->subblock_r+SUBBLOCK_HASH(off);
  subblock_hash_t *subblock_a = conf->subblock_a+SUBBLOCK_HASH(off);
  script_t *script;
  s32 cutoff2_r, cutoff2_a;
  static mpz_t rRes, aRes;
  static mpz_t quot, base, exponent;
  static s32 initialized = 0;
  relation_t Rel;
  s32 a = (block+first_block)*BLOCK_SIZE+off;
  s32 tiny_rfb_size = conf->tiny_rfb_size;
  s32 med_rfb_size = conf->med_rfb_size;
  s32 rfb_inf_size = conf->rfb_inf_size;
  s32 tiny_afb_size = conf->tiny_afb_size;
  s32 med_afb_size = conf->med_afb_size;
  s32 afb_inf_size = conf->afb_inf_size;
  s32 last_p;
  s32 rfactors = 0, afactors = 0;
  char printbuf[256];

  if (!initialized) {
    mpz_init(quot);
    mpz_init(rRes);
    mpz_init(aRes);
    mpz_init_set_ui(base, 2);
    mpz_init(exponent);
    initialized = 1;
  }

  /*********************************************************
   * Note that even when b is large, trial factoring       *
   * can take a *lot* of time, especially when sieving     *
   * is fast. Any number of tricks is justified to         *
   * accelerate this phase.                                *
   * The high-level view is that we can use sieve          *
   * information to reduce the work involved, and should   *
   * start with the fast parts of the process to weed      *
   * out unsuitable values early. We should also alternate *
   * between the rational and algebraic halves of the      *
   * trial factoring process, to determine as soon as      *
   * possible when one side or the other is not going      *
   * to be smooth                                          *
   *********************************************************/

  if (b%2==0)
    a = 2*a+1;
  if (sign==POSITIVE)
    Rel.a = a;
  else
    Rel.a = -a;
  Rel.b = b;
  mpz_mul_si(rRes, FB->y1, Rel.a);
  mpz_addmul_ui(rRes, FB->y0, b);
  mpz_evalF(aRes, Rel.a, b, FB->f);
  mpz_abs(rRes, rRes);
  mpz_abs(aRes, aRes);

  /**************************************************
   * Use the actual sieve values to determine the   *
   * cutoff for continuing the trial factoring; this*
   * treats unusually small sieve values more fairly*
   **************************************************/
  cutoff2_r = fplog_mpz(rRes, FB->log_rlb) - conf->cutoff2_r;
  if (cutoff2_r<0)
    cutoff2_r = 0;
  cutoff2_a = fplog_mpz(aRes, FB->log_alb) - conf->cutoff2_a;
  if (cutoff2_a<0)
    cutoff2_a = 0;

  /******************************************
   * Begin by folding in the unsieved AFB   *
   * entries, since the algebraic norm is   *
   * usually harder to reach. Update the    *
   * log value too, since this will account *
   * for repeated small factors             *
   ******************************************/
  last_p = -1;
  for (i=0; i<tiny_afb_size; i++) {
    sieve_fb_t *fb = conf->afb + i;
    s32 p = fb->p;
    s32 root = fb->root;
    u32 logp = fb->logp;

    if (p!=last_p)
      rem = a%p;     /* a%p doesn't change if p doesn't change */
    last_p = p;
    if (sign==NEGATIVE) {
      if (root)
        root = p-root;
    }
    if (rem==root) {
      Rel.aFactors[afactors++] = p;
      j = mpz_fdiv_q_ui(quot, aRes, p);
      do {
        abits += logp;
        mpz_swap(quot, aRes);
        j = mpz_fdiv_q_ui(quot, aRes, p);
      } while(j==0);
    }
  }

  /*****************************************
   * continue with the algebraic primes at *
   * infinity, then see if this value meets*
   * the cutoff                            *
   *****************************************/
  for (i=0; i<afb_inf_size; i++) {
    sieve_fb_t *fb = conf->afb_inf + i;
    s32 p = fb->p;
    u32 logp = fb->logp;

    if (b%p==0) {
      j = mpz_fdiv_q_ui(quot, aRes, p);
      if (j==0) {
        Rel.aFactors[afactors++] = p;
        do {
          abits += logp;
          mpz_swap(quot, aRes);
          j = mpz_fdiv_q_ui(quot, aRes, p);
        } while(j==0);
      }
    } 
  }

  if (abits<=cutoff2_a)
    return 0;
      
  /*************************************
   * repeat with the unsieved rational *
   * factor base primes                *
   *************************************/
  for (i=0; i<tiny_rfb_size; i++) {
    sieve_fb_t *fb = conf->rfb + i;
    s32 p = fb->p;
    s32 root = fb->root;
    u32 logp = fb->logp;

    rem = a%p;
    if (sign==NEGATIVE) {
      if (root)
        root = p-root;
    }
    if (rem==root) {
      Rel.rFactors[rfactors++] = p;
      j = mpz_fdiv_q_ui(quot, rRes, p);
      do {
        rbits += logp;
        mpz_swap(quot, rRes);
        j = mpz_fdiv_q_ui(quot, rRes, p);
      } while(j==0);
    }
  }

  for (i=0; i<rfb_inf_size; i++) {
    sieve_fb_t *fb = conf->rfb_inf + i;
    s32 p = fb->p;
    u32 logp = fb->logp;

    if (b%p==0) {
      j = mpz_fdiv_q_ui(quot, rRes, p);
      if (j==0) {
        Rel.rFactors[rfactors++] = p;
        do {
          rbits += logp;
          mpz_swap(quot, rRes);
          j = mpz_fdiv_q_ui(quot, rRes, p);
        } while(j==0);
      }
    } 
  }

  if (rbits<=cutoff2_r)
    return 0;
      
  conf->num_tf_attempts[SUBBLOCK_NUM_AVG_B-1]++;

  if (!conf->use_subblock) {
    /*********************************
     * Now trial factor the small    *
     * algebraic factor base primes. * 
     *********************************/
    last_p = -1;
    for (i=tiny_afb_size; i<med_afb_size; i++) {
      sieve_fb_t *fb = conf->afb + i;
      s32 p = fb->p;
      s32 root = fb->root;
  
      if (p!=last_p)
        rem = a%p;
      last_p = p;
      if (sign==NEGATIVE) {
        if (root)
          root = p-root;
      }
      if (rem==root) {
        j = mpz_fdiv_q_ui(quot, aRes, p);
        if (j==0) {
          Rel.aFactors[afactors++] = p;
          do {
            mpz_swap(quot, aRes);
            j = mpz_fdiv_q_ui(quot, aRes, p);
          } while(j==0);
        }
      }
    }
  
    /*******************************************
     * For the medium and large AFB primes,    *
     * walk through the appropriate hash bin   *
     * looking for hash offsets that match 'a'.*
     * Since each hash entry also contains the *
     * prime involved, this is equivalent to   *
     * the trial division step above. Not only *
     * is this division-free, the number of    *
     * entries to check is much smaller than   *
     * the factor base size.                   *
     *
     * Since the hashtable entries were halved *
     * if b is even, the sieve offset 'off' is *
     * always directly comparable to the values*
     * out of the hashtable                    *
     *******************************************/
    script = ahash->script;
    while (script!=NULL) {
      update_t *updates = script->list;
      k = script->num_used;
      for (i=0; i<k; i++) {
        if (i%8==0)
          PREFETCH(updates+i+16);
        if (updates[i].h==off) {
          s32 p = updates[i].p;
          j = mpz_fdiv_q_ui(quot, aRes, p);
          if (j==0) {
            Rel.aFactors[afactors++] = p;
            do {
              mpz_swap(quot, aRes);
              j = mpz_fdiv_q_ui(quot, aRes, p);
            } while(j==0);
          }
        }
      }
      script = script->next;
    }
  }
  else {

    /***********************************
     * when a large number of trial    *
     * factorings are expected, the    *
     * explicit remainder operations   *
     * above will take forever. Use the* 
     * subblock table instead, which   *
     * involves much less work         *
     ***********************************/
    update_t *updates = subblock_a->list;
    k = subblock_a->num_used;
    for (i=0; i<k; i++) {
      if (i%8==0)
        PREFETCH(updates+i+16);
      if (updates[i].h==off) {
        s32 p = updates[i].p;
        j = mpz_fdiv_q_ui(quot, aRes, p);
        if (j==0) {
          Rel.aFactors[afactors++] = p;
          do {
            mpz_swap(quot, aRes);
            j = mpz_fdiv_q_ui(quot, aRes, p);
          } while(j==0);
        }
      }
    }
  }

  /********************************************
   * Do not attempt to factor the leftover    *
   * residue, but do test for primality (since*
   * this is much cheaper than trial factoring*
   * the rational residue).                   *
   * Use a base-2 pseudoprime test, since     *
   * mpz_probable_prime does trial division   *
   * which would just waste time in this ap-  *
   * plication.                               *
   ********************************************/
  if (mpz_cmp(aRes, conf->LP2_max_a)>0)
    return 0;
  if (mpz_cmp(aRes, conf->LP2_min_a)>0) {
    mpz_sub_ui(exponent, aRes, 1);
    mpz_powm(exponent, base, exponent, aRes);
    if (mpz_cmp_ui(exponent, 1)==0)
      return 0;
  }
  else if (mpz_cmp(aRes, conf->LP1_max_a)>0) {
    return 0;
  }

  /**********************************
   * Do the trial factoring for the *
   * rational sieve value           *
   **********************************/
  if (!conf->use_subblock) {
    for (i=tiny_rfb_size; i<med_rfb_size; i++) {
      sieve_fb_t *fb = conf->rfb + i;
      s32 p = fb->p;
      s32 root = fb->root;
  
      rem = a%p;
      if (sign==NEGATIVE) {
        if (root)
          root = p-root;
      }
      if (rem==root) {
        Rel.rFactors[rfactors++] = p;
        j = mpz_fdiv_q_ui(quot, rRes, p);
        do {
          mpz_swap(quot, rRes);
          j = mpz_fdiv_q_ui(quot, rRes, p);
        } while(j==0);
      }
    }
  
    script = rhash->script;
    while (script!=NULL) {
      update_t *updates = script->list;
      k = script->num_used;
      for (i=0; i<k; i++) {
        if (i%8==0)
          PREFETCH(updates+i+16);
        if (updates[i].h==off) {
          s32 p = updates[i].p;
          Rel.rFactors[rfactors++] = p;
          j = mpz_fdiv_q_ui(quot, rRes, p);
          do {
            mpz_swap(quot, rRes);
            j = mpz_fdiv_q_ui(quot, rRes, p);
          } while(j==0);
        }
      }
      script = script->next;
    }
  }
  else {
    update_t *updates = subblock_r->list;
    k = subblock_r->num_used;
    for (i=0; i<k; i++) {
      if (i%8==0)
        PREFETCH(updates+i+16);
      if (updates[i].h==off) {
        s32 p = updates[i].p;
        j = mpz_fdiv_q_ui(quot, rRes, p);
        if (j==0) {
          Rel.rFactors[rfactors++] = p;
          do {
            mpz_swap(quot, rRes);
            j = mpz_fdiv_q_ui(quot, rRes, p);
          } while(j==0);
        }
      }
    }
  }

  /*******************************************
   * Test the rational residue for primality *
   *******************************************/
  if (mpz_cmp(rRes, conf->LP2_max_r)>0)
    return 0;
  if (mpz_cmp(rRes, conf->LP2_min_r)>0) {
    mpz_sub_ui(exponent, rRes, 1);
    mpz_powm(exponent, base, exponent, rRes);
    if (mpz_cmp_ui(exponent, 1)==0)
      return 0;
  }
  else if (mpz_cmp(rRes, conf->LP1_max_r)>0) {
    return 0;
  }

  Rel.rFSize = rfactors;
  Rel.aFSize = afactors;

  /**********************************
   * Find the large primes, if any, *
   * using SQUFOF                   *
   **********************************/
  if (mpz_cmp(aRes, conf->LP1_max_a)<=0) {
    Rel.a_p[0] = Rel.a_p[1] = 1; 
    Rel.a_p[2] = mpz_get_si(aRes);
  }
  else {
    i = squfof(aRes);
    if (i<=1)
      return 0;
    mpz_fdiv_q_ui(aRes, aRes, i);
    if (mpz_cmp_si(conf->LP1_max_a, i)<=0 ||
        mpz_cmp(conf->LP1_max_a, aRes)<=0) {
      return 0;
    }
    Rel.a_p[0] = 1; Rel.a_p[1] = i;
    Rel.a_p[2] = mpz_get_si(aRes);
  }

  if (mpz_cmp(rRes, conf->LP1_max_r)<=0) {
    Rel.p[0] = Rel.p[1] = 1; 
    Rel.p[2] = mpz_get_si(rRes);
  }
  else {
    i = squfof(rRes);
    if (i<=1)
      return 0;
    mpz_fdiv_q_ui(rRes, rRes, i);
    if (mpz_cmp_si(conf->LP1_max_r, i)<=0 ||
        mpz_cmp(conf->LP1_max_r, rRes)<=0) {
      return 0;
    }
    Rel.p[0] = 1; Rel.p[1] = i;
    Rel.p[2] = mpz_get_si(rRes);
  }

  /*******************************
   * Yay! Another relation found *
   *******************************/
  clMakeOutputLine(printbuf, &Rel);
  printRelation(conf, printbuf);
  return 1;
}

/***********************************************************/     
void add2Hashtable(update_t *cache, s32 cache_max, 
                   hashtable_t *hash, s32 p, s32 h, u8 logp, 
                   script_t **mempool)
{
  update_t *cache_base = cache+hash->cache_off;
  update_t new_update;

  new_update.p = p;
  new_update.h = h & (BLOCK_SIZE-1);
  new_update.logp = logp;
  
  /**********************************************
   * new hash entries never update the hashtable*
   * directly; instead they always go into the  *
   * cache. This conserves (processor) cache,   *
   * is much more TLB-friendly and involves     *
   * fewer safety checks in the common case     *
   **********************************************/
  if (hash->cache_used>=cache_max)
    flushCache(cache, hash, mempool);
  cache_base[hash->cache_used++] = new_update;
}

/***********************************************************/     
void flushCache(update_t *cache, hashtable_t *hash, script_t **mempool)
/***********************************************************
 * Perform a block copy of the cache for one hashtable     *
 * bin, copying to the linked list of entries for that bin.* 
 ***********************************************************/
{
  update_t *cache_base = cache+hash->cache_off;
  update_t *update_base;
  script_t *script = hash->script;
  s32 used = hash->cache_used;
  s32 i;

  /*****************************************
   * Allocate a new group of sieve updates *
   * if necessary                          *
   *****************************************/
  if (script==NULL || script->num_used+used >= MAX_UPDATES) {
    script_t *new_script = *mempool;
    if (new_script==NULL)
      new_script = (script_t *)malloc(sizeof(script_t));
    else
      *mempool = new_script->next;

    new_script->next = hash->script;
    new_script->num_used = 0;
    script = hash->script = new_script;
  }
  update_base = script->list+script->num_used;

  /****************************************************
   * we don't want to bring the destination into      *
   * processor cache, since that would kick out       *
   * our cache of hashtable updates. Use non-temporal *
   * stores if they're available; they help a *lot*   *
   * on the athlon, but they're slower on the K8      *
   ****************************************************/
#if defined(__GNUC__) && defined(__i386__)
  for (i=0; i<(used & ~7); i+=8) {
    asm("movq 0(%1,%2,8), %%mm0 \n\t"
        "movq 8(%1,%2,8), %%mm1 \n\t"
        "movq 16(%1,%2,8), %%mm2 \n\t"
        "movq 24(%1,%2,8), %%mm3 \n\t"
        "movq 32(%1,%2,8), %%mm4 \n\t"
        "movq 40(%1,%2,8), %%mm5 \n\t"
        "movq 48(%1,%2,8), %%mm6 \n\t"
        "movq 56(%1,%2,8), %%mm7 \n\t"
        "movntq %%mm0, 0(%0,%2,8) \n\t"
        "movntq %%mm1, 8(%0,%2,8) \n\t"
        "movntq %%mm2, 16(%0,%2,8) \n\t"
        "movntq %%mm3, 24(%0,%2,8) \n\t"
        "movntq %%mm4, 32(%0,%2,8) \n\t"
        "movntq %%mm5, 40(%0,%2,8) \n\t"
        "movntq %%mm6, 48(%0,%2,8) \n\t"
        "movntq %%mm7, 56(%0,%2,8)"
        ::"r"(update_base), "r"(cache_base), "r"(i): "memory");
  }

  for (; i<used; i++) {
    asm("movq 0(%1,%2,8), %%mm0 \n\t"
        "movntq %%mm0, 0(%0,%2,8) "
        ::"r"(update_base), "r"(cache_base), "r"(i): "memory");
  }
  asm("femms");

#else
  for (i=0; i<(used & ~7); i+=8) {
    PREFETCH(update_base+16);
    {
      update_t t0 = cache_base[i];
      update_t t1 = cache_base[i+1];
      update_t t2 = cache_base[i+2];
      update_t t3 = cache_base[i+3];
      update_t t4 = cache_base[i+4];
      update_t t5 = cache_base[i+5];
      update_t t6 = cache_base[i+6];
      update_t t7 = cache_base[i+7];
      update_base[i] = t0;
      update_base[i+1] = t1;
      update_base[i+2] = t2;
      update_base[i+3] = t3;
      update_base[i+4] = t4;
      update_base[i+5] = t5;
      update_base[i+6] = t6;
      update_base[i+7] = t7;
    }
  }
  for (; i<used; i++)
    update_base[i] = cache_base[i];
#endif

  hash->cache_used = 0;
  script->num_used += used;
}

/***********************************************************/     
void add2Subblock(subblock_hash_t *entry, s32 p, s32 h)
{
  update_t new_update;

  new_update.p = p;
  new_update.h = h & (BLOCK_SIZE-1);
  
  if (entry->num_used==entry->num_alloc) {
    entry->num_alloc *= 2;
    entry->list = (update_t *)realloc(entry->list, 
                     entry->num_alloc*sizeof(update_t));
  }
  entry->list[entry->num_used++] = new_update;
}

/***********************************************************/     
#define SAVEFILE_BUF_SIZE 65536
char savefile_buf[SAVEFILE_BUF_SIZE];
int savefile_buf_off = 0;

/***********************************************
 * Relations are buffered because I've noticed *
 * in other programs that high-speed formatted *
 * output sometimes corrupts itself            *
 ***********************************************/
void printRelation(sieve_conf_t *conf, char *buf) 
{
  if (savefile_buf_off + strlen(buf) + 1 >= SAVEFILE_BUF_SIZE) {
    flushSavefile(conf);
  }

  savefile_buf_off += sprintf(savefile_buf + 
        savefile_buf_off, "%s\n", buf);
}

void flushSavefile(sieve_conf_t *conf) 
{
  fprintf(conf->outfp, "%s", savefile_buf);
  fflush(conf->outfp);
  savefile_buf_off = 0;
}

/**************************************************************/
void clMakeOutputLine(char *str, relation_t *R)
/*************************************************************
 * Custom version of makeOutputLine the assumes a relation_t *
 * contains prime factors and not offsets into the factor    *
 * bases. This allows conf->rfb_size and conf->afb_size to   *
 * differ from their counterparts in an nfs_fb_t structure   *
 *************************************************************/
{ int i, numR=0, numA=0;
  char s[128];

  sprintf(str, "%ld,%ld:", R->a, R->b);
  for (i=0; i<R->rFSize; i++) {
    if (numR==0)
      sprintf(s, "%lx", R->rFactors[i]);
    else sprintf(s, ",%lx", R->rFactors[i]);
    strcat(str, s);
    numR++;
  }
  for (i=0; i<MAX_LARGE_RAT_PRIMES; i++) {
    if (R->p[i] > 1) {
      if (numR>0)
        sprintf(s, ",%lx", R->p[i]);
      else
        sprintf(s, "%lx", R->p[i]);
      strcat(str, s);
      numR++;
    }
  }
  strcat(str, ":");
  for (i=0; i<R->aFSize; i++) {
    if (numA==0)
      sprintf(s, "%lx", R->aFactors[i]);
    else sprintf(s, ",%lx", R->aFactors[i]);
    strcat(str, s);
    numA++;
  }
  for (i=0; i<MAX_LARGE_ALG_PRIMES; i++) {
    if (R->a_p[i] > 1) {
      if (numA>0)
        sprintf(s, ",%lx", R->a_p[i]);
      else
        sprintf(s, "%lx", R->a_p[i]);
      strcat(str, s);
      numA++;
    }
  }
} 
