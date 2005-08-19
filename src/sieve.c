/**************************************************************/
/* sieve.c                                                    */
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
#include <signal.h>
#include <time.h>
#include <gmp.h>
#include "ggnfs.h"

int clForceStop = 0;   /* referenced in clsieve.c */

#define _CATCH_SIGNALS
#define DEFAULT_OUTNAME "spairs.out"
#define DEFAULT_REPORTINTERVAL 5.0
#define SPEED_SMOOTH_B 10

/* These are (probably unreasonable) defaults! They should absolutely
   not be relied on as anything but perhaps a starting point!
  I am in a hurry at the moment.
*/
#define MAX_DEFAULTS 6
s32 defaultRLIM[MAX_DEFAULTS]={0,1000,2000,  4000,  8000,   15000};
s32 defaultALIM[MAX_DEFAULTS]={0,1000,2000,  4000,  8000,   15000};
s32 defaultMPR[MAX_DEFAULTS]={0,  0,  0,20000, 20000, 500000};
s32 defaultMPA[MAX_DEFAULTS]={0,  0,  0,20000, 20000, 500000};
int  defaultLP[MAX_DEFAULTS] ={0,  0,  1,    1,     1,      2};
int  defaultLPA[MAX_DEFAULTS]={0,  0,  1,    1,     1,      2};
double defaultRLambda[MAX_DEFAULTS]={1.0,1.1, 1.3, 1.4,    1.5,     1.6};
double defaultALambda[MAX_DEFAULTS]={1.0,1.1, 1.3, 1.4,    1.5,     1.6};
#define DEFAULT_FB_NAME "default.fb"

#define USAGE \
"[OPTIONS]\n"\
"--help              : show this help and exit\n"\
"-fb <filename>      : use factor base in file <filename>\n"\
"-j  <filename>      : use job file <filename>\n"

#define START_MSG \
"\n"\
" __________________________________________________________ \n"\
"|         This is the sieve program for GGNFS.             |\n"\
"| Version: %-25s                       |\n"\
"| This program is copyright 2004, Chris Monico, and subject|\n"\
"| to the terms of the GNU General Public License version 2.|\n"\
"|__________________________________________________________|\n"

/**************************************************/
void catch_term(int sig)
/**************************************************/
{ 
  printf("Term signal caught. Doing clean shutdown...\n");
  fflush(stdout);
  clForceStop=1;
}

/*******************************************************/
void setDependentDefaults(nfs_sieve_job_t *job)
/*******************************************************/
/* `job' has at least a number to be factored and the  */
/* `nDigits' field filled in. This function will set   */
/* any other fields which haven't yet been initialized */
/* to some default values depending on `nDigits'.      */
/*******************************************************/
{ int n, w=0;
  char warn[]="** Warning: using default value for %s!";

  n = job->nDigits/10;
  n = MIN(n, MAX_DEFAULTS-1);

  if ((job->FB.rfb_size == -1)&&(job->FB.rLim == -1)) {
    job->FB.rLim = defaultRLIM[n]; w++;
    msgLog("", warn, "rlim");
  }
  if ((job->FB.afb_size == -1)&&(job->FB.aLim == -1)) {
    job->FB.aLim = defaultALIM[n]; w++;
    msgLog("", warn, "alim");
  }
  if (job->FB.maxP_r == -1) {
    job->FB.maxP_r = defaultMPR[n]; w++;
    msgLog("", warn, "mpr");
  }
  if (job->FB.maxP_a == -1) {
    job->FB.maxP_a = defaultMPA[n]; w++;
    msgLog("", warn, "mpa");
  }
  if (job->FB.maxLP == -1)
    job->FB.maxLP = defaultLP[n];
  if (job->FB.maxLPA == -1)
    job->FB.maxLPA = defaultLPA[n];

  if (job->FB.rfb_lambda < 0 ) {
    job->FB.rfb_lambda = job->rfb_lambda = defaultRLambda[n]; w++;
    msgLog("", warn, "rlambda");
  }
  if (job->FB.afb_lambda < 0) {
    job->FB.afb_lambda = job->afb_lambda = defaultALambda[n]; w++;
    msgLog("", warn, "alambda");
  }
  if (w) 
    msgLog("", "** Warning: Unreliable default values used for %d important fields!", w);

}

/*******************************************************/
int parseJobFile(nfs_sieve_job_t *job, char *fName)
/*******************************************************/
{ FILE *fp;
  char token[256], value[256], thisLine[1024];
  int  t;

  if (!(fp = fopen(fName, "r"))) {
    printf("Error opening %s for read!\n", fName);
    return -1;
  }

  readPoly(fp, &job->FB);	/* read n, m, poly, etc. */
  job->nDigits = mpz_sizeinbase(job->FB.n, 10);
  rewind(fp); /* Back up to the beginning of the file. */
  
  /* Initialize fields to empty, or defaults which are independent of
     the number being factored. */
  job->type = SIEVE_CLASSICAL;
  job->a0 = job->a1 = job->b0 = job->b1 = 0;
  job->qIndex0 = job->qIndex1 = 0;
  job->sieveArea = 0.0;
  job->FB.rfb_size = job->FB.afb_size = -1;
  job->FB.rLim = job->FB.aLim = -1;
  job->FB.maxP_r = job->FB.maxP_a = -1;
  job->FB.maxLP = job->FB.maxLPA = -1;
  job->FB.MFB_r = job->FB.MFB_a = -1;
  job->rfb_lambda = job->afb_lambda = -1;
  strcpy(job->outName, DEFAULT_OUTNAME);

  while (!(feof(fp))) {
    thisLine[0] = 0;
    fgets(thisLine, 1023, fp);
    if ((sscanf(thisLine, "%255s %255s", token, value)==2) && (thisLine[0] != '#')) {
      if (strncmp(token, "a0:",3)==0) {
        job->a0 = atol(value);
      } else if (strncmp(token, "a1:",3)==0) {
        job->a1 = atol(value);
      } else if (strncmp(token, "b0:",3)==0) {
        job->b0 = atol(value);
      } else if (strncmp(token, "b1:",3)==0) {
        job->b1 = atol(value);
      } else if (strncmp(token, "q0:",3)==0) {
        job->qIndex0 = atol(value);
      } else if (strncmp(token, "q1:",3)==0) {
        job->qIndex1 = atol(value);
      } else if (strncmp(token, "sieveArea:",10)==0) {
        job->sieveArea = atof(value);
      } else if (strncmp(token, "rlambda:",8)==0) {
        job->rfb_lambda = atof(value);
      } else if (strncmp(token, "alambda:",8)==0) {
        job->afb_lambda = atof(value);
      } else if (strncmp(token, "fbname:", 7)==0) {
        strncpy(job->fbName, value, MAXFNAMESIZE);
      } else if (strncmp(token, "rfb_size:", 9)==0) {
        job->FB.rfb_size = atol(value);
      } else if (strncmp(token, "afb_size:", 9)==0) {
        job->FB.afb_size = atol(value);
      } else if (strncmp(token, "rlim:", 5)==0) {
        job->FB.rLim = atol(value);
      } else if (strncmp(token, "alim:", 5)==0) {
        job->FB.aLim = atol(value);
      } else if (strncmp(token, "lpbr:", 5)==0) {
        t = atoi(value); t = MIN(t, 31);
        job->FB.maxP_r = (1<<t);
      } else if (strncmp(token, "lpba:", 5)==0) {
        t = atoi(value); t = MIN(t, 31);
        job->FB.maxP_a = (1<<t);
      } else if (strncmp(token, "mfbr:", 5)==0) {
        job->FB.MFB_r = atoi(value);
      } else if (strncmp(token, "mfba:", 5)==0) {
        job->FB.MFB_a = atoi(value);
      } else if (strncmp(token, "npr:", 4)==0) {
        job->FB.maxLP = atoi(value);
      } else if (strncmp(token, "npa:", 4)==0) {
        job->FB.maxLPA = atoi(value);
      }
    }
  }
  fclose(fp);
  return 0;
}

/***************************************************/
int tryMakeFB(nfs_sieve_job_t *J)
/***************************************************/
/* Try to create a factor base from the info in J. */
/* Return 0 on success, nonzero on failure.        */
/***************************************************/
{ int res;

  /* setDependentDefaults() will try to fill in missing fields with
     (probably bad!) default values. It will log warnings
     when defaults values are used in important fields.
  */
  setDependentDefaults(J);

  res = createFB(&J->FB, J->fbName);
  return res;
}

/**************************************************/
int main(int argC, char *args[])
{ char            fbname[MAXFNAMESIZE], jobfile[MAXFNAMESIZE];
  int             i, fbLoaded=0;
  int             numRelations;
  nfs_sieve_job_t Job;
  double          now;
  double          startTime;
  
  fbname[0]=jobfile[0]=0;
  forceStop=0;
  strncpy(Job.fbName, DEFAULT_FB_NAME, MAXFNAMESIZE);

  printf(START_MSG, GGNFS_VERSION);

  /* Register signal handler for clean shutdowns. */

#ifdef _CATCH_SIGNALS
  signal(SIGTERM, catch_term);
  signal(SIGINT, catch_term);
#ifndef _MSC_VER
  signal(SIGQUIT, catch_term);
#endif
#endif

  /***************************/
  /* Parse command-line args */
  /***************************/
  for (i=1; i<argC; i++) {
    if (strcmp(args[i], "-fb")==0) {
      if ((++i)<argC) {
        strcpy(fbname, args[i]);
      }
    }
    else if (strcmp(args[i], "-j")==0) {
      if ((++i)<argC) {
        strcpy(jobfile, args[i]);
      }
    }
    else if (strcmp(args[i], "--help")==0) {
      printf("USAGE: %s %s\n", args[0], USAGE);
      return 0;
    }
  }
  if (!(jobfile[0])) {
    printf("USAGE: %s %s\n", args[0], USAGE);
    return 0;
  }
  msgLog("", "GGNFS-%s : sieve", GGNFS_VERSION);

  /***********************************/
  /* Initialize and get factor base. */
  /***********************************/
  initFB(&Job.FB);
  
  /******************/
  /* Get parameters */
  /******************/
  if (parseJobFile(&Job, jobfile)) {
    printf("Some error occured parsing job file. Cannot continue.\n");
    exit(-1);
  }
  if (loadFB(fbname, &Job.FB)) {
    fbLoaded=0;
    fprintf(stderr, "Could not load factor base from %s. "
		    "Attempting to create...\n", fbname);
    if (tryMakeFB(&Job)) {
      fprintf(stderr, "Factor base creation failed! Aborting.\n");
      exit(-1);
    }
  }

  startTime = sTime();
  numRelations = clSieve(&Job);
  now = sTime();
  msgLog("", "Found %d relations in %1.1lf sec", 
                      numRelations, now-startTime);
  printf("Found %d relations in %1.1lf sec\n", numRelations, now-startTime);
  return clForceStop;
}
