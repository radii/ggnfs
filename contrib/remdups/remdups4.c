/* remdups4.c
  By Greg Childers
  This is a standalone one-pass duplicate relations remover;
      only syntax (a,b:<csv-of-factors>:<csv-of-factors>) is checked and incomplete lines are discarded;
      validity of relations is not tested (nor is polynomial needed).
  Hashing of (a,b) values was changed to accomodate for gnfs projects with huge skews (S.Batalov).

  This version is a filter (stdin, stdout):
    you may redirect stdout to /dev/null and/or
    you may use many input relation files (or pipe from zcat, bzcat).
    However, this needs porting for Windows (use remdups.c instead or get CygWin).

  No makefile needed:
      cc -O3 -o remdups4 remdups4.c

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned long uint32;
typedef long int32;
typedef unsigned long long uint64;

#define MDVAL 327673

int strccnt(const char *s, int c)
{
	const unsigned char *us = (const unsigned char *) s;
	const unsigned char uc = c;
	int n = 0;
	if (!uc) return 1;
	while (*us) if (*us++ == uc) n++;
	return n;
}
int main(int argc, char **argv) {
	char buf[512];
        int mm,nn,av;
	FILE *badfile;
	uint64 ** arra;
        int n[MDVAL]={0};
	int numbad=0, numdups=0, numuniq=0,numskip=0;
	int DIM=1000;

        if (argc == 2) {
	  DIM=atoi(argv[1]); 
        } else {
	  fprintf(stderr,"\nusage: cat relations.file(s) | %s table_size > out_file \n"
		  "\t table_size is a number (1000-3000 recommended)\n\n", argv[0]);
	  exit(-1);
	}

	badfile = fopen("badrels.txt", "a");
	if (badfile == NULL) {
		fprintf(stderr,"cannot open badfile\n");
		exit(-1);
	}
	
	if ((DIM<50)||(DIM>100000)) { 
		fprintf(stderr,"DIM should be between 50 and 100000!\n");
		exit(1);
	}

	/* initialize arrays */
	arra = (uint64**)malloc(MDVAL * sizeof(uint64 *));
	for(mm = 0; mm < MDVAL; mm++){
#define CHUNK 2048
	  if((mm % CHUNK) == 0) {
	    arra[mm] = (uint64*)malloc(CHUNK * DIM * sizeof(uint64));
	    if(!arra[mm]) {
	      fprintf(stderr, "out of memory\n");
	      exit(1);
	    }
	    /* memset(arra[mm],0,CHUNK * DIM * sizeof(uint64)); */  /* unnecessary */
	  } else {
	    arra[mm] = arra[mm-1] + DIM;
	  }
	}

	while (fgets(buf, sizeof(buf), stdin)) {
		char *tmp, *field_end;
		uint64 a;
		int32 i, j, p;

		if (buf[0] == '#') {
			printf("%s", buf);
			continue;
		}
		 
		if (buf[0] == 'N') {
			printf("%s", buf);
			continue;
		}
		
		for(tmp = buf; *tmp && isspace(*tmp); tmp++);

		/* Hash used to be in a and b bins; it worked well for SNFS */
		/* However, for gnfs, the bins were very shallow */
		/* New hash value a is a nonsensical base-11 hybrid of both a and b -SB 2009 */
		if(*tmp=='-') {a=10; tmp++;} else a=0;
		for(         ; *tmp ; tmp++) {
		  if (isdigit(*tmp)) a=11*a+(*tmp-'0');
		  else if(*tmp==',') a=11*a+10;
		  else {
		    if(*tmp==':') {
		      if ((tmp-2>buf && tmp[-2]==',' && tmp[-1]=='0') || strccnt(tmp+1,':')==1)
		      	break;
		    }
		    fprintf(badfile, "%s", buf);
		    numbad++;
		    goto skip_;
		  }
		}
		p=a%MDVAL;
		for (i=0;i<n[p];i++)
		  if (a==arra[p][i]) { numdups++; goto skip_; }
		if (n[p]<DIM) arra[p][n[p]++]=a; /* Not quitting after a whole lot of work! Just stop extending the bin --SB. */
		else numskip++;
		numuniq++;
		if(numuniq % 500000 == 0)
			fprintf(stderr,"\r %.1fM relns \r", numuniq/1000000.0);
		for(tmp=buf;*tmp;tmp++) *tmp=tolower(*tmp); /* will compress better */
		printf("%s", buf);
	skip_:  ;
	}
		
        fprintf(stderr,"Found %d unique, %d duplicate, and %d bad relations.\n",numuniq, numdups, numbad);
	for (av=0,mm=1,nn=n[0];mm<MDVAL;mm++) {if (n[mm]>nn) nn=n[mm]; av+=n[mm];}
	fprintf(stderr,"Largest dimension used: %d of %d\n",nn,DIM);
	fprintf(stderr,"Average dimension used: %.1f of %d\n",((double)av)/MDVAL,DIM);
	if(nn>=DIM) fprintf(stderr,"*** Some redundant relations may have been retained (increase DIM)\n");
	if(numskip) fprintf(stderr,"*** %d (quasi-unique) relations were not hashed\n",numskip);
	if(numbad) fprintf(badfile, "\n");	/* usually last reln line is truncated and ends up in the bad bin. We don't want them to stick together */
	fclose(badfile);
	return 0;
}
