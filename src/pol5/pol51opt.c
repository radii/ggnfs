/* pol51opt.c

   Copyright (C) 2005 T. Kleinjung.
   This file is part of pol5, distributed under the terms of the
   GNU General Public Licence and WITHOUT ANY WARRANTY.
                                                                                                               
  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  2111-1307, USA.
                                                                                                               
  CJM, 2/18/05: Hacked up for inclusion in GGNFS. It was originally written by
  T. Kleinjung and/or Jens Franke.
*/


/*
#define ZEIT
*/

#include "ggnfs.h"

#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include "gmp.h"
#include "if.h"
#include <limits.h>
#include "fnmatch.h"
#include <string.h>
#if defined (_MSC_VER) || defined (__MINGW32__) || defined (MINGW32)
#include "getopt.h"
#include "rint.h"
#endif

#define START_MESSAGE \
"----------------------------------------------------\n"\
"|    pol51opt GNFS polynomial selection program    |\n"\
"| This program is copyright (c) 2005, by Thorsten  |\n"\
"| Kleinjung and Jens Franke, and is subject to the |\n"\
"| terms of the GNU GPL version 2.                  |\n"\
"| This program is part of gnfs4linux.              |\n"\
"----------------------------------------------------\n"


#define  MAX_PRIME_PROJ       100
#define  MAX_PRIME_AFF        200
#define  NPROJ_PRIMES          25
#define  NAFF_PRIMES           46
#define  MAX_X       1000000  /* !!! */
#define  MAX_Y           100  /* !!! */
#define  SIEVELEN           8192

unsigned int primes[46]={
 2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
 179, 181, 191, 193, 197, 199
};

mpz_t gmp_N;
int compress;
char *input_line=NULL;
size_t input_line_alloc=0;
char *base_name, *filename_data, *output_name, *m_name;
FILE *outputfile, *m_file;
mpz_t gmp_a[6], gmp_b[6], gmp_help1, gmp_help2, gmp_help3, gmp_help4;
mpz_t gmp_lina[2], gmp_linb[2], gmp_p, gmp_d, gmp_m, gmp_mb;
double dbl_a[6], dbl_m, dbl_p, dbl_d, dbl_b[6], dbl_mb;
double skewness, sk_b, e_value;
double Emax, alpha_proj;

double max_norm_1, max_norm_2, min_e, log_max_norm_2;
double pol_norm;
int xmin, xmax, ymin, ymax;
double bound0, bound1, area;
unsigned int p_bound;

/* statistics */
int64_t lanz0=0, lanz1=0, lanz2=0, lanz3=0, lanz4=0, lanz5=0;

/* ----------------------------------------------- */

int divmod(int a, int b, int pk)  /* a/b mod pk , a, b < 2^15 */
{
  int inv;

  inv=invert(b,pk); if (!inv) Schlendrian("divmod\n");
  inv=inv%pk; if (inv<0) inv+=pk;
  inv*=(a%pk);
  inv=inv%pk; if (inv<0) inv+=pk;
  return inv;
}

/* ----------------------------------------------- */

void get_options(int argc, char **argv)
{
  char c;

  base_name=NULL;
  compress=0;
  max_norm_1=1e20; max_norm_2=1e18; min_e=0.;
  p_bound=2000;
  bound0=1e7; bound1=5e6; area=1e16;
  while ((c=getopt(argc,argv,"A:b:e:F:f:n:N:P:vz")) != (char)(-1)) {
    switch(c) {
    case 'b':
      base_name=optarg;
      break;
    case 'A':
      if(sscanf(optarg,"%lf",&area)!=1)
        complain("Bad argument to -A!\n");
      break;
    case 'F':
      if(sscanf(optarg,"%lf",&bound0)!=1)
        complain("Bad argument to -F!\n");
      break;
    case 'f':
      if(sscanf(optarg,"%lf",&bound1)!=1)
        complain("Bad argument to -f!\n");
      break;
    case 'e':
      if(sscanf(optarg,"%lf",&min_e)!=1)
        complain("Bad argument to -e!\n");
      break;
    case 'n':
      if(sscanf(optarg,"%lf",&max_norm_1)!=1)
        complain("Bad argument to -n!\n");
      break;
    case 'N':
      if(sscanf(optarg,"%lf",&max_norm_2)!=1)
        complain("Bad argument to -N!\n");
      break;
    case 'P':
      if(sscanf(optarg,"%u",&p_bound)!=1)
        complain("Bad argument to -P!\n");
      break;
    case 'v':
      verbose++;
      break;
    case 'z':
      compress=1;
      break;
    default:
      fprintf(stderr,"Bad Option %c\n",(char)c);
      Schlendrian("");
    }
  }
  if (base_name==NULL) complain("argument '-b base_name' is necessary\n");
  asprintf(&filename_data,"%s.data",base_name);
  log_max_norm_2=log(max_norm_2);
}



void read_data()
{
  char *line;
  int read;
  FILE *fi;

  read=0;
  mpz_init(gmp_N);
  if ((fi=fopen(filename_data,"r"))==NULL)
    complain("File not found: %s\n",filename_data);
  while (getline(&input_line,&input_line_alloc,fi)>0) {
    line=input_line;
    if (*line=='N') {
      mpz_set_str(gmp_N,line+2,10); read|=1;
      continue;
    }
  }
  fclose(fi);
  if (read!=1) complain("Not enough data in file: %s\n",filename_data);
}


int find_m_name()
{
  struct stat statbuf;
  char *tmp_name;

  asprintf(&m_name,"%s.51.m",base_name);
  if (stat(m_name,&statbuf)) {
    asprintf(&m_name,"%s.51.m.gz",base_name);
    if (stat(m_name,&statbuf)) return 0;
  } else {
    asprintf(&tmp_name,"%s.51.m.gz",base_name);
    if (!stat(tmp_name,&statbuf))
      complain("Both files %s and %s exist.\n",m_name,tmp_name);
    free(tmp_name);
  }
  return 1;
}


void open_m()
{
  char *input_cmd;

  if (fnmatch("*.gz",m_name,0)) {
    if ((m_file=fopen(m_name,"r"))==NULL) {
      fprintf(stderr,"%s ",m_name);
      complain("cannot open file %s",m_name);
    }
  } else {
    asprintf(&input_cmd,"zcat %s",m_name);
    if ((m_file=popen(input_cmd,"r"))==NULL) {
      fprintf(stderr,"%s ",m_name);
      complain("cannot open file %s",m_name);
    }
    free(input_cmd);
  }
}



void close_m()
{
  if (fnmatch("*.gz",m_name,0)) fclose(m_file);
  else pclose(m_file);
}


void open_outputfile()
{
  char *output_cmd;

  output_cmd=NULL;
  if (compress!=0) {
    asprintf(&output_name,"%s.cand.gz",base_name);
    asprintf(&output_cmd,"gzip --best --stdout >> %s",output_name);
    if ((outputfile=popen(output_cmd,"w"))==NULL) {
      fprintf(stderr,"%s ",output_name);
      complain("cannot open file %s",output_name);
    }
/*    if (stat(output_name,&statbuf)==0) complain("Output file exists!\n");*/
    free(output_cmd);
  } else {
    asprintf(&output_name,"%s.cand",base_name);
    if ((outputfile=fopen(output_name,"a"))==NULL) {
      fprintf(stderr,"%s ",output_name);
      complain("cannot open file %s",output_name);
    }
/*    if (stat(output_name,&statbuf)==0) complain("Output file exists!\n");*/
  }
  setbuf(outputfile,NULL);
}


void close_outputfile()
{
  if (fnmatch("*.gz",output_name,0)) fclose(outputfile);
  else pclose(outputfile);
}


int read_a5pd()
{
  char *line, *end;

  if (getline(&input_line,&input_line_alloc,m_file)<=0) return 0;
  line=input_line;
  end=strchr(line,' '); if (end==NULL) return -1;
  *end=0;
  if (mpz_set_str(gmp_a[5],line,10)) return -1;
  line=end+1;
  end=strchr(line,' '); if (end==NULL) return -1;
  *end=0;
  if (mpz_set_str(gmp_p,line,10)) return -1;
  line=end+1;
  if (mpz_set_str(gmp_d,line,10)) return -1;
  return 1;
}

/* ---------------------------------------------------------- */

double ifs(double *coeff, double skewness)  /* degree=5 */
{
  double sq[6], d, s, res;
  int i, j, k;

  for (i=0; i<6; i++) sq[i]=coeff[i]*coeff[i];
  for (i=0; i<4; i++) {
    d=2*coeff[i]; k=i+1;
    for (j=i+2; j<6; j+=2, k++) sq[k]+=d*coeff[j];
  }
  s=skewness; res=0.;
  res+=s*sq[3]/35.;
  res+=sq[2]/s/35.;
  s*=(skewness*skewness);
  res+=s*sq[4]/27.;
  res+=sq[1]/s/27.;
  s*=(skewness*skewness);
  res+=s*sq[5]/11.;
  res+=sq[0]/s/11.;
  return res;
}


#define COMPUTE_IFS   s2=sc; res=0.; \
  res+=s2*sq[3]/35.; res+=sq[2]/s2/35.;   \
  s2*=(sc*sc);   \
  res+=s2*sq[4]/27.; res+=sq[1]/s2/27.;   \
  s2*=(sc*sc);   \
  res+=s2*sq[5]/11.; res+=sq[0]/s2/11.


double find_best_skewness(double *coeff, double s0)  /* degree=5 */
{
  double sq[6], d, s, sc, s2, res, ds, v;
  int i, j, k;

  for (i=0; i<6; i++) sq[i]=coeff[i]*coeff[i];
  for (i=0; i<4; i++) {
    d=2*coeff[i]; k=i+1;
    for (j=i+2; j<6; j+=2, k++) sq[k]+=d*coeff[j];
  }
  ds=10.;
  s=s0; sc=s; COMPUTE_IFS; v=res;
  while (ds>2.) {
    sc=s+ds; COMPUTE_IFS;
    if (res<v) { s=sc; v=res; ds*=1.1; continue; }
    if (ds<s) {
      sc=s-ds; COMPUTE_IFS;
      if (res<v) { s=sc; v=res; ds*=1.1; continue; }
    }
    ds*=0.9;
  }
  return s;
}

/* -------------------- translations --------------------------- */

void translate_dbl(double *dtarg, double *dsrc, int k)
{
  int i;
  double d, dk;

  for (i=0; i<6; i++) dtarg[i]=dsrc[i];
  dk=(double)(-k);
  for (i=0; i<5; i++) dtarg[i]+=(double)(i+1)*dsrc[i+1]*dk;
  d=dk*dk;
  dtarg[0]+=d*dsrc[2]; dtarg[1]+=3*d*dsrc[3];
  dtarg[2]+=6*d*dsrc[4]; dtarg[3]+=10*d*dsrc[5];
  d*=dk;
  dtarg[0]+=d*dsrc[3]; dtarg[1]+=4*d*dsrc[4]; dtarg[2]+=10*d*dsrc[5];
  d*=dk;
  dtarg[0]+=d*dsrc[4]; dtarg[1]+=5*d*dsrc[5];
  d*=dk;
  dtarg[0]+=d*dsrc[5];
}


void translate_gmp(mpz_t *gmp_z, mpz_t *lin, mpz_t gmp_mz, int k)
{

/*printf("trans %d\n",k);
mpz_out_str(stdout,10,gmp_p); printf(" ");
mpz_out_str(stdout,10,gmp_d); printf("\n");
for (ii=0; ii<6; ii++) { mpz_out_str(stdout,10,gmp_z[5-ii]); printf(" "); }
printf("\n");
*/
  if (k>=0) { mpz_set_ui(gmp_help1,k); mpz_neg(gmp_help1,gmp_help1); }
  else mpz_set_ui(gmp_help1,-k);  /* help1=-k */
  mpz_set(gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[1],gmp_help3);
  mpz_add(gmp_z[0],gmp_z[0],gmp_help2);
  mpz_mul(gmp_help3,gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[2],gmp_help3);
  mpz_add(gmp_z[0],gmp_z[0],gmp_help2);
  mpz_mul(gmp_help3,gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[3],gmp_help3);
  mpz_add(gmp_z[0],gmp_z[0],gmp_help2);
  mpz_mul(gmp_help3,gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[4],gmp_help3);
  mpz_add(gmp_z[0],gmp_z[0],gmp_help2);
  mpz_mul(gmp_help3,gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[5],gmp_help3);
  mpz_add(gmp_z[0],gmp_z[0],gmp_help2);
/* a0<-a0-a1*k+a2*k^2-a3*k^3+a4*k^4-a5*k^5 */

  mpz_set(gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[2],gmp_help3);
  mpz_mul_ui(gmp_help2,gmp_help2,2);
  mpz_add(gmp_z[1],gmp_z[1],gmp_help2);
  mpz_mul(gmp_help3,gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[3],gmp_help3);
  mpz_mul_ui(gmp_help2,gmp_help2,3);
  mpz_add(gmp_z[1],gmp_z[1],gmp_help2);
  mpz_mul(gmp_help3,gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[4],gmp_help3);
  mpz_mul_ui(gmp_help2,gmp_help2,4);
  mpz_add(gmp_z[1],gmp_z[1],gmp_help2);
  mpz_mul(gmp_help3,gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[5],gmp_help3);
  mpz_mul_ui(gmp_help2,gmp_help2,5);
  mpz_add(gmp_z[1],gmp_z[1],gmp_help2);
/* a1<-a1-2*a2*k+3*a3*k^2-4*a4*k^3+5*a5*k^4 */

  mpz_set(gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[3],gmp_help3);
  mpz_mul_ui(gmp_help2,gmp_help2,3);
  mpz_add(gmp_z[2],gmp_z[2],gmp_help2);
  mpz_mul(gmp_help3,gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[4],gmp_help3);
  mpz_mul_ui(gmp_help2,gmp_help2,6);
  mpz_add(gmp_z[2],gmp_z[2],gmp_help2);
  mpz_mul(gmp_help3,gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[5],gmp_help3);
  mpz_mul_ui(gmp_help2,gmp_help2,10);
  mpz_add(gmp_z[2],gmp_z[2],gmp_help2);
/* a2<-a2-3*a3*k+6*a4*k^2-10*a5*k^3 */

  mpz_set(gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[4],gmp_help3);
  mpz_mul_ui(gmp_help2,gmp_help2,4);
  mpz_add(gmp_z[3],gmp_z[3],gmp_help2);
  mpz_mul(gmp_help3,gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[5],gmp_help3);
  mpz_mul_ui(gmp_help2,gmp_help2,10);
  mpz_add(gmp_z[3],gmp_z[3],gmp_help2);
/* a3<-a3-4*a4*k+10*a5*k^2 */

  mpz_set(gmp_help3,gmp_help1);
  mpz_mul(gmp_help2,gmp_z[5],gmp_help3);
  mpz_mul_ui(gmp_help2,gmp_help2,5);
  mpz_add(gmp_z[4],gmp_z[4],gmp_help2);
/* a4<-a4-5*a5*k */

  mpz_sub(gmp_mz,gmp_mz,gmp_help1);
/* m<-m+k */

  mpz_mul(gmp_help2,lin[1],gmp_help1);
  mpz_add(lin[0],lin[0],gmp_help2);
/* lin0<-lin0-lin1*k */

/*  mpz_set(gmp_p,lin[1]);
  mpz_neg(gmp_d,lin[0]);*/

/*printf("result: \n");
mpz_out_str(stdout,10,gmp_p); printf(" ");
mpz_out_str(stdout,10,gmp_d); printf("\n");
for (ii=0; ii<6; ii++) { mpz_out_str(stdout,10,gmp_z[5-ii]); printf(" "); }
printf("\n");
*/
}

/* ----------------------------------------------- */

int pol_expand()
{
/* compute coefficients */
  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_help4);
  mpz_mul(gmp_help4,gmp_help4,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_a[5]);
  mpz_sub(gmp_help3,gmp_N,gmp_help4);
  mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_p);
  if (mpz_sgn(gmp_help1)) return 0;

  if (mpz_cmp_ui(gmp_p,1)==0) mpz_set_ui(gmp_help2,1);
  else {
    if (!mpz_invert(gmp_help2,gmp_d,gmp_p)) return 0;
  }
  mpz_mul(gmp_help4,gmp_help2,gmp_help2);
  mpz_mul(gmp_help4,gmp_help4,gmp_help4);
  mpz_mul(gmp_help4,gmp_help4,gmp_help3);
  mpz_fdiv_r(gmp_a[4],gmp_help4,gmp_p);
  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_help4);
  mpz_mul(gmp_help4,gmp_help4,gmp_a[4]);
  mpz_sub(gmp_help3,gmp_help3,gmp_help4);
  mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_p);
  if (mpz_sgn(gmp_help1)) return 0;

  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_d);
  mpz_fdiv_q(gmp_a[3],gmp_help3,gmp_help4);
  mpz_fdiv_q(gmp_a[3],gmp_a[3],gmp_p);
  mpz_mul(gmp_a[3],gmp_a[3],gmp_p);
  mpz_mul(gmp_help4,gmp_help2,gmp_help2);
  mpz_mul(gmp_help4,gmp_help4,gmp_help2);
  mpz_mul(gmp_help4,gmp_help4,gmp_help3);
  mpz_fdiv_r(gmp_help1,gmp_help4,gmp_p);
  mpz_add(gmp_a[3],gmp_a[3],gmp_help1);
  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_a[3]);
  mpz_sub(gmp_help3,gmp_help3,gmp_help4);
  mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_p);
  if (mpz_sgn(gmp_help1)) return 0;

  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_fdiv_q(gmp_a[2],gmp_help3,gmp_help4);
  mpz_fdiv_q(gmp_a[2],gmp_a[2],gmp_p);
  mpz_mul(gmp_a[2],gmp_a[2],gmp_p);
  mpz_mul(gmp_help4,gmp_help2,gmp_help2);
  mpz_mul(gmp_help4,gmp_help4,gmp_help3);
  mpz_fdiv_r(gmp_help1,gmp_help4,gmp_p);
  mpz_add(gmp_a[2],gmp_a[2],gmp_help1);
  mpz_mul(gmp_help4,gmp_d,gmp_d);
  mpz_mul(gmp_help4,gmp_help4,gmp_a[2]);
  mpz_sub(gmp_help3,gmp_help3,gmp_help4);
  mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_p);
  if (mpz_sgn(gmp_help1)) return 0;

  mpz_fdiv_q(gmp_a[1],gmp_help3,gmp_d);
  mpz_fdiv_q(gmp_a[1],gmp_a[1],gmp_p);
  mpz_mul(gmp_a[1],gmp_a[1],gmp_p);
  mpz_mul(gmp_help4,gmp_help3,gmp_help2);
  mpz_fdiv_r(gmp_help1,gmp_help4,gmp_p);
  mpz_add(gmp_a[1],gmp_a[1],gmp_help1);
  mpz_mul(gmp_help4,gmp_d,gmp_a[1]);
  mpz_sub(gmp_help3,gmp_help3,gmp_help4);
  mpz_fdiv_qr(gmp_help3,gmp_help1,gmp_help3,gmp_p);
  if (mpz_sgn(gmp_help1)) return 0;

  mpz_set(gmp_a[0],gmp_help3);

  mpz_fdiv_qr(gmp_help1,gmp_a[3],gmp_a[3],gmp_d);
  mpz_add(gmp_help2,gmp_a[3],gmp_a[3]);
  if (mpz_cmp(gmp_d,gmp_help2)<0) {
    mpz_sub(gmp_a[3],gmp_a[3],gmp_d);
    mpz_add_ui(gmp_help1,gmp_help1,1);
  }
  mpz_mul(gmp_help1,gmp_help1,gmp_p);
  mpz_add(gmp_a[4],gmp_a[4],gmp_help1);

  mpz_fdiv_qr(gmp_help1,gmp_a[2],gmp_a[2],gmp_d);
  mpz_add(gmp_help2,gmp_a[2],gmp_a[2]);
  if (mpz_cmp(gmp_d,gmp_help2)<0) {
    mpz_sub(gmp_a[2],gmp_a[2],gmp_d);
    mpz_add_ui(gmp_help1,gmp_help1,1);
  }
  mpz_mul(gmp_help1,gmp_help1,gmp_p);
  mpz_add(gmp_a[3],gmp_a[3],gmp_help1);

  mpz_set(gmp_lina[1],gmp_p);
  mpz_neg(gmp_lina[0],gmp_d);
  mpz_invert(gmp_m,gmp_p,gmp_N);
  mpz_mul(gmp_m,gmp_m,gmp_d);
  mpz_mod(gmp_m,gmp_m,gmp_N);

  if (verbose>2) {
    int i;

    printf("pol-expand\npol0: ");
    for (i=5; i>=0; i--) { mpz_out_str(stdout,10,gmp_a[i]); printf(" "); }
    printf("\npol1: ");
    mpz_out_str(stdout,10,gmp_lina[1]); printf(" ");
    mpz_out_str(stdout,10,gmp_lina[0]); printf("\n\n");
  }

  return 1;
}


void optimize_1()
{
  int dk, i, niter;
  int64_t di1, di0;
  double dbl_a0[6];
  double value, v0;
  double s, s0, ds, d;

  dk=16;
  niter=0;
  for (i=0; i<6; i++) dbl_a[i]=mpz_get_d(gmp_a[i]);
  dbl_d=mpz_get_d(gmp_d);
  s=1.;
  if (dbl_a[4]!=0)
    s=sqrt(fabs(dbl_a[2]/dbl_a[4]));
  s0=1.;
  if (dbl_a[3]!=0)
    s0=fabs(dbl_a[2]/dbl_a[3]);
  if (s0>s)
    s=s0;
  ds=2.;
  di1=(int64_t)(dbl_a[1]/dbl_d);
  di0=(int64_t)(dbl_a[0]/dbl_d);
  value=ifs(dbl_a,s);
  while (1) {
    if (niter>10000) {
      fprintf(stderr,"too many iterations in optimize_1\n");
      break;
    }
    if ((dk<2) && (di1<2) && (di0<2) && (ds<0.001)) break;
/* skewness */
    s0=s*(1.+ds); v0=ifs(dbl_a,s0);
    if (v0<value) {
      s=s0; value=v0; ds*=1.1;
    } else {
      s0=s/(1.+ds); v0=ifs(dbl_a,s0);
      if (v0<value) {
        s=s0; value=v0; ds*=1.1;
      } else ds/=1.1;
    }
/*printf("v=%e, s=%f, S  ",value,s); for (ii=0; ii<6; ii++) printf("%e ",dbl_a[5-ii]); printf("\n");
for (ii=0; ii<6; ii++) { mpz_out_str(stdout,10,gmp_a[5-ii]); printf(" "); }
printf("\n");*/
/* translation */
    translate_dbl(dbl_a0,dbl_a,dk); v0=ifs(dbl_a0,s);
    if (v0<value) {
      translate_gmp(gmp_a,gmp_lina,gmp_m,dk);
      value=v0; dk=(int)(1+1.1*(double)dk);
    } else {
      translate_dbl(dbl_a0,dbl_a,-dk); v0=ifs(dbl_a0,s);
      if (v0<value) {
        translate_gmp(gmp_a,gmp_lina,gmp_m,-dk);
        value=v0; dk=(int)(1+1.1*(double)dk);
      } else {
        dk=(int)(0.9*(double)dk-1);
        if (dk<1) dk=1;
      }
    }
  mpz_set(gmp_p,gmp_lina[1]);
  mpz_neg(gmp_d,gmp_lina[0]);
    for (i=0; i<6; i++) dbl_a[i]=mpz_get_d(gmp_a[i]);
    dbl_d=mpz_get_d(gmp_d); dbl_p=mpz_get_d(gmp_p);
/*printf("v=%e, s=%f, T  ",value,s); for (ii=0; ii<6; ii++) printf("%e ",dbl_a[5-ii]); printf("\n");
for (ii=0; ii<6; ii++) { mpz_out_str(stdout,10,gmp_a[5-ii]); printf(" "); }
printf("\n");*/
/* "rotation" with i1*x*(px-m) */
    d=(double)di1;
    dbl_a[2]+=d*dbl_p; dbl_a[1]-=d*dbl_d; v0=ifs(dbl_a,s);
    if (v0<value) {
      mpz_set_sll(gmp_help1,di1);
      mpz_mul(gmp_help2,gmp_help1,gmp_p);
      mpz_add(gmp_a[2],gmp_a[2],gmp_help2);
      mpz_mul(gmp_help1,gmp_help1,gmp_d);
      mpz_sub(gmp_a[1],gmp_a[1],gmp_help1);
      value=v0; di1=(int64_t)(1+1.1*(double)di1);
    } else {
      dbl_a[2]-=2*d*dbl_p; dbl_a[1]+=2*d*dbl_d; v0=ifs(dbl_a,s);
      if (v0<value) {
        mpz_set_sll(gmp_help1,-di1);
        mpz_mul(gmp_help2,gmp_help1,gmp_p);
        mpz_add(gmp_a[2],gmp_a[2],gmp_help2);
        mpz_mul(gmp_help1,gmp_help1,gmp_d);
        mpz_sub(gmp_a[1],gmp_a[1],gmp_help1);
        value=v0; di1=(int64_t)(1+1.1*(double)di1);
      } else {
        dbl_a[2]+=d*dbl_p; dbl_a[1]-=d*dbl_d; /* set it back */
        di1=(int64_t)(0.9*(double)di1-1);
        if (di1<1) di1=1;
      }
    }
/*printf("v=%e, s=%f, R  ",value,s); for (ii=0; ii<6; ii++) printf("%e ",dbl_a[5-ii]); printf("\n");
for (ii=0; ii<6; ii++) { mpz_out_str(stdout,10,gmp_a[5-ii]); printf(" "); }
printf("\n");*/
/* "rotation" with i0*(px-m) */
    d=(double)di0;
    dbl_a[1]+=d*dbl_p; dbl_a[0]-=d*dbl_d; v0=ifs(dbl_a,s);
    if (v0<value) {
      mpz_set_sll(gmp_help1,di0);
      mpz_mul(gmp_help2,gmp_help1,gmp_p);
      mpz_add(gmp_a[1],gmp_a[1],gmp_help2);
      mpz_mul(gmp_help1,gmp_help1,gmp_d);
      mpz_sub(gmp_a[0],gmp_a[0],gmp_help1);
      value=v0; di0=(int64_t)(1+1.1*(double)di0);
    } else {
      dbl_a[1]-=2*d*dbl_p; dbl_a[0]+=2*d*dbl_d; v0=ifs(dbl_a,s);
      if (v0<value) {
        mpz_set_sll(gmp_help1,-di0);
        mpz_mul(gmp_help2,gmp_help1,gmp_p);
        mpz_add(gmp_a[1],gmp_a[1],gmp_help2);
        mpz_mul(gmp_help1,gmp_help1,gmp_d);
        mpz_sub(gmp_a[0],gmp_a[0],gmp_help1);
        value=v0; di0=(int64_t)(1+1.1*(double)di0);
      } else {
        dbl_a[1]+=d*dbl_p; dbl_a[0]-=d*dbl_d; /* set it back */
        di0=(int64_t)(0.9*(double)di0-1);
        if (di0<1) di0=1;
      }
    }
    niter++;
/*printf("v=%e, s=%f, R  ",value,s); for (ii=0; ii<6; ii++) printf("%e ",dbl_a[5-ii]); printf("\n");
for (ii=0; ii<6; ii++) { mpz_out_str(stdout,10,gmp_a[5-ii]); printf(" "); }
printf("\n");*/
  }
#if 0
{
printf("pol=");
for (i=0; i<6; i++) { mpz_out_str(stdout,10,gmp_a[5-i]); printf("*x^%ld+ ",5-i); }
printf("0 \n");
}
{
for (i=0; i<6; i++) printf("%f ",dbl_a[5-i]); printf("\n");
}
#endif
  skewness=s;
  pol_norm=sqrt(ifs(dbl_a,s));

  if (verbose>2) {
    int i;

    printf("optimize 1\nskewness: %.2f norm: %.4e\npol0: ",skewness,pol_norm);
    for (i=5; i>=0; i--) { mpz_out_str(stdout,10,gmp_a[i]); printf(" "); }
    printf("\npol1: ");
    mpz_out_str(stdout,10,gmp_lina[1]); printf(" ");
    mpz_out_str(stdout,10,gmp_lina[0]); printf("\n\n");
  }
}


void optimize_2(double *norm_ptr)
{
  int dk, i, niter;
  double dbl_b0[6];
  double value, v0;
  double s, s0, ds;

  dk=16; ds=2.;
  niter=0;
  for (i=0; i<6; i++) dbl_b[i]=mpz_get_d(gmp_b[i]);
  s=skewness;
  value=ifs(dbl_b,s);
  while (1) {
    if (niter>10000) {
      fprintf(stderr,"too many iterations in optimize_2\n");
      break;
    }
    if ((dk<2) && (ds<1.001)) break;
/* skewness */
    s0=s*ds; v0=ifs(dbl_b,s0);
    if (v0<value) {
      s=s0; value=v0; ds*=1.1;
    } else {
      s0=s/ds; v0=ifs(dbl_b,s0);
      if (v0<value) {
        s=s0; value=v0; ds*=1.1;
      } else ds=1.+(ds-1.)/1.1;
    }
/* translation */
    translate_dbl(dbl_b0,dbl_b,dk); v0=ifs(dbl_b0,s);
    if (v0<value) {
      translate_gmp(gmp_b,gmp_linb,gmp_mb,dk);
      value=v0; dk=(int)(1+1.1*(double)dk);
    } else {
      translate_dbl(dbl_b0,dbl_b,-dk); v0=ifs(dbl_b0,s);
      if (v0<value) {
        translate_gmp(gmp_b,gmp_linb,gmp_mb,-dk);
        value=v0; dk=(int)(1+1.1*(double)dk);
      } else {
        dk=(int)(0.9*(double)dk-1);
        if (dk<1) dk=1;
      }
    }
    for (i=0; i<6; i++) dbl_b[i]=mpz_get_d(gmp_b[i]);
    niter++;
  }
  sk_b=s;
  *norm_ptr=sqrt(value);

  if (verbose>3) {
    int i;

    printf("optimize 2\nskewness: %.2f norm: %.4e\npol0: ",sk_b,*norm_ptr);
    for (i=5; i>=0; i--) { mpz_out_str(stdout,10,gmp_b[i]); printf(" "); }
    printf("\npol1: ");
    mpz_out_str(stdout,10,gmp_linb[1]); printf(" ");
    mpz_out_str(stdout,10,gmp_linb[0]); printf("\n\n");
  }
}


void optimize_3(double *norm_ptr, double *eptr, double *alphaptr)
{
  int dk, i, niter;
  double dbl_b0[6];
  double e, new_e, alpha;
  double s, s0, ds;
  double dbl_linb[2];

  alpha=*alphaptr;
  dk=16; ds=2.;
  niter=0;
  for (i=0; i<6; i++) dbl_b[i]=mpz_get_d(gmp_b[i]);
/*  s=skewness;*/
s=sk_b;
  dbl_linb[1]=mpz_get_d(gmp_linb[1]);
  dbl_linb[0]=mpz_get_d(gmp_linb[0]);
  murphy_e(&e,5,dbl_b,1,dbl_linb,alpha,0.,s);
  if (e<min_e) {
    sk_b=s;
    *eptr=e;
    *norm_ptr=sqrt(ifs(dbl_b,s));
    return;
  }
  while (1) {
    if (niter>10000) {
      fprintf(stderr,"too many iterations in optimize_3\n");
      break;
    }
    if ((dk<8) && (ds<1.001)) break;
/* skewness */
    s0=s*ds; murphy_e(&new_e,5,dbl_b,1,dbl_linb,alpha,0.,s0);
    if (new_e>e) {
      s=s0; e=new_e; ds*=1.1;
    } else {
      s0=s/ds; murphy_e(&new_e,5,dbl_b,1,dbl_linb,alpha,0.,s0); 
      if (new_e>e) {
        s=s0; e=new_e; ds*=1.1;
      } else ds=1.+(ds-1.)/1.1;
    }
/* translation */
    translate_dbl(dbl_b0,dbl_b,dk); dbl_linb[0]-=dk*dbl_linb[1];
    murphy_e(&new_e,5,dbl_b0,1,dbl_linb,alpha,0.,s);
    dbl_linb[0]+=dk*dbl_linb[1];
    if (new_e>e) {
      translate_gmp(gmp_b,gmp_linb,gmp_mb,dk);
      dbl_linb[0]-=dk*dbl_linb[1];
      e=new_e; dk=(int)(1+1.1*(double)dk);
    } else {
      dbl_linb[0]+=dk*dbl_linb[1];
      translate_dbl(dbl_b0,dbl_b,-dk);
      murphy_e(&new_e,5,dbl_b0,1,dbl_linb,alpha,0.,s);
      dbl_linb[0]-=dk*dbl_linb[1];
      if (new_e>e) {
        translate_gmp(gmp_b,gmp_linb,gmp_mb,-dk);
        dbl_linb[0]+=dk*dbl_linb[1];
        e=new_e; dk=(int)(1+1.1*(double)dk);
      } else {
        dk=(int)(0.9*(double)dk-1);
        if (dk<1) dk=1;
      }
    }
    for (i=0; i<6; i++) dbl_b[i]=mpz_get_d(gmp_b[i]);
    niter++;
  }
  sk_b=s;
  murphy_en(&e,5,dbl_b0,1,dbl_linb,alpha,0.,s,10000);
  if (e>=min_e) {
    if (p_bound<2000) {
      if (verbose>2) printf("(%f,%g) -> ",alpha,e);
      compute_alpha_exact(&alpha,5,NULL,gmp_b,2000);
      murphy_en(&e,5,dbl_b0,1,dbl_linb,alpha,0.,s,10000);
      if (verbose>2) printf("(%f,%g)\n",alpha,e);
    }
  }

  *alphaptr=alpha;
  *eptr=e;
  *norm_ptr=sqrt(ifs(dbl_b,s));

  if (verbose>2) {
    int i;

    printf("optimize 3\nskewness: %.2f norm: %.4e Murphy_E: %.3e\npol0: ",sk_b,*norm_ptr,*eptr);
    for (i=5; i>=0; i--) { mpz_out_str(stdout,10,gmp_b[i]); printf(" "); }
    printf("\npol1: ");
    mpz_out_str(stdout,10,gmp_linb[1]); printf(" ");
    mpz_out_str(stdout,10,gmp_linb[0]); printf("\n\n");
  }
}


/* ----------------------------------------------- */


typedef struct {
  int primepower;     /* p^k */
  int rooti;          /* zero =i mod p^k */
  int J;              /* f(i)+J*(i-m)=0 mod p^k */
  int step;           /* stepwidth p^(k-l) */
  int start;          /* sievearray[start] has zero =i mod p^k */
  int lineinc;        /* increment for change y->y+1 */
  unsigned short v;            /* scaled value */
  double value;        /* log(p)/(p^(k-1)*(p+1)) */
} primelist;

primelist pl[MAX_PRIME_AFF*MAX_PRIME_AFF];
int primelistlen;
double limit, tval, st_alpha;
unsigned short sievearray[SIEVELEN], cut;
int sievelen, nsubsieves;

int prep_len, prep_p_len[NAFF_PRIMES];
unsigned int prep_p_begin[NAFF_PRIMES][MAX_PRIME_AFF*2];
unsigned int *prep_p[NAFF_PRIMES];


void compute_proj_alpha()
{
  unsigned int i, j, k;
  unsigned int p, p2, p3;
  unsigned int w, b3, b4, b5;
  double value, dp, dl;
  double table[MAX_PRIME_PROJ];

  value=0.;
  for (k=0; k<NPROJ_PRIMES; k++) {
    p=primes[k];
    if (mpz_mod_ui(gmp_help1,gmp_a[5],p)==0) {
      dp=(double)p;
      value-=log(dp)/(dp+1.);
      if (MAX_PRIME_PROJ/p/p<p) {  /* consider only p^2 */
        p2=p*p;
        if (mpz_mod_ui(gmp_help1,gmp_a[5],p2)==0) {
          if (mpz_mod_ui(gmp_help1,gmp_a[4],p)==0)
            value-=log(dp)*dp/(dp*dp-1);
          else
            value-=log(dp)/(dp*dp-1);
        } else {
          if (mpz_mod_ui(gmp_help1,gmp_a[4],p)==0) value-=log(dp)/(dp*dp-1);
        }
      } else {
        p2=p*p; p3=p*p2;
        for (i=0; i<p2; i++) table[i]=0.;
        dl=log(dp)/((dp*dp-1.)*dp*dp*dp);
        for (i=0; i<p; i++) table[i*p]=dl;
        table[0]=log(dp)*(2.+1./(dp-1.))/((dp*dp-1.)*dp*dp*dp);
        b5=mpz_mod_ui(gmp_help1,gmp_a[5],p3); b5/=p;
        b4=mpz_mod_ui(gmp_help1,gmp_a[4],p2);
        b3=mpz_mod_ui(gmp_help1,gmp_a[3],p2);
/* (a5/p)*x^2 + a4*x*(y/p) + a3*p*(y/p)^2, b5=(a5/p),b4=a4,b3=a3,j=y/p,i=x */
        for (j=0; j<p2; j++)
          for (i=0; i<p2; i++)
            if (i%p) {
              w=b5*i*i+b4*i*j+b3*p*j*j; w=w%p2;
              value-=table[w];
            }
      }
    }
  }
  alpha_proj=value;
}

/* ----------------------------------------------- */


void write_polynomial_51(FILE *fi, int deg, mpz_t *coeff1, mpz_t *coeff2, mpz_t m, double skewness, double norm, double alpha, double murphy_e)
{
  int i;

  fprintf(fi,"BEGIN POLY ");
  fprintf(fi,"#skewness %.2f norm %.2e alpha %.2f Murphy_E %.2e\n",skewness,norm,alpha,murphy_e);
  for (i=0; i<deg+1; i++) {
    fprintf(fi,"X%d ",deg-i);
    mpz_out_str(fi,10,coeff1[deg-i]);
    fprintf(fi,"\n");
  }
  if (mpz_cmp_ui(coeff2[1],1)) {
    fprintf(fi,"Y1 ");
    mpz_out_str(fi,10,coeff2[1]);
    fprintf(fi,"\nY0 ");
    mpz_out_str(fi,10,coeff2[0]);
    fprintf(fi,"\n");
  }
  fprintf(fi,"M ");
  mpz_out_str(fi,10,m);
  fprintf(fi,"\nEND POLY\n");
}


void check(int x, int y)
{
  int i;
  double alpha, norm, murphye, alpha_max;

#ifdef ZEIT
zeita(6);
#endif
  if (verbose>3) printf("check at i=%d, j=%d\n",x,y);

  for (i=0; i<6; i++) mpz_set(gmp_b[i],gmp_a[i]);
  mpz_set_si(gmp_help1,y);
  mpz_mul(gmp_help2,gmp_help1,gmp_p);
  mpz_add(gmp_b[2],gmp_b[2],gmp_help2);
  mpz_mul(gmp_help1,gmp_help1,gmp_d);
  mpz_sub(gmp_b[1],gmp_b[1],gmp_help1);
  mpz_set_si(gmp_help1,x);
  mpz_mul(gmp_help2,gmp_help1,gmp_p);
  mpz_add(gmp_b[1],gmp_b[1],gmp_help2);
  mpz_mul(gmp_help1,gmp_help1,gmp_d);
  mpz_sub(gmp_b[0],gmp_b[0],gmp_help1);
  mpz_set(gmp_mb,gmp_m);
  for (i=0; i<2; i++) mpz_set(gmp_linb[i],gmp_lina[i]);
#ifdef ZEIT
zeitb(6);
zeita(8);
#endif
  optimize_2(&norm);
#ifdef ZEIT
zeitb(8);
zeita(7);
#endif
  alpha_max=log(max_norm_2/norm);
  if (compute_alpha(&alpha,5,NULL,gmp_b,alpha_max)) {
#ifdef ZEIT
    zeitb(7);
#endif
    if (verbose>2) printf("failed\n");
    return;
  }
  if (verbose>2) printf("alpha: %.3f\n",alpha);
#ifdef ZEIT
zeitb(7);
#endif

  if (norm*exp(alpha)>max_norm_2) {
    if (verbose>2) printf("failed\n");
    return;
  }

  if (verbose>2) printf("i: %d, j: %d\n",x,y);
#ifdef ZEIT
zeita(9);
#endif
  optimize_3(&norm,&murphye,&alpha);
#ifdef ZEIT
zeitb(9);
#endif
  if (murphye<min_e) {
    if (verbose>1) printf("f");
    return;
  }
  if (verbose>1) printf("P");
  write_polynomial_51(outputfile,5,gmp_b,gmp_linb,gmp_mb,sk_b,norm,alpha,murphye);
}

/* ----------------------------------------------- */

double compute_standard_alpha()
{
  int p, pk, k;
  double dp, dpk, dlp, dlog, standardalpha;

  standardalpha=0.;
  for (k=0; k<NAFF_PRIMES; k++) {
    p=primes[k];
    dp=(double)p; dlp=(double)(log(dp));
    pk=p;
    while (pk<MAX_PRIME_AFF) {
      dpk=(double)pk; dlog=dlp/dpk*dp/(dp+1.);
      standardalpha-=dlog;
      pk*=p;
    }
  }
  return standardalpha;
}


void compute_sieving_area()    /* sehr grob !!! */
{
  int x0, x1, xx, y0, y1, yy;
  double da1, da0, dd;
  double v, v0;

  v=max_norm_1*exp(-alpha_proj);

  da1=dbl_a[1]; da0=dbl_a[0]; dd=-dbl_d;
  y0=0; y1=MAX_Y;
  dbl_a[1]=da1-((double)y1)*dd;
  v0=sqrt(ifs(dbl_a,find_best_skewness(dbl_a,skewness)));
  if (v0>v) {
    while (y1-y0>1) {
      yy=(y0+y1)/2;
      dbl_a[1]=da1-((double)yy)*dd;
      v0=sqrt(ifs(dbl_a,find_best_skewness(dbl_a,skewness)));
      if (v0>v) y1=yy; else y0=yy;
    }
  }
  ymin=-y1;
  y0=0; y1=MAX_Y;
  dbl_a[1]=da1+((double)y1)*dd;
  v0=sqrt(ifs(dbl_a,find_best_skewness(dbl_a,skewness)));
  if (v0>v) {
    while (y1-y0>1) {
      yy=(y0+y1)/2;
      dbl_a[1]=da1+((double)yy)*dd;
      v0=sqrt(ifs(dbl_a,find_best_skewness(dbl_a,skewness)));
      if (v0>v) y1=yy; else y0=yy;
    }
  }
  ymax=y1;
  dbl_a[1]=da1;

  x0=0; x1=MAX_X;
  dbl_a[0]=da0-((double)x1)*dd;
  v0=sqrt(ifs(dbl_a,find_best_skewness(dbl_a,skewness)));
  if (v0>v) {
    while (x1-x0>1) {
      xx=(x0+x1)/2;
      dbl_a[0]=da0-((double)xx)*dd;
      v0=sqrt(ifs(dbl_a,find_best_skewness(dbl_a,skewness)));
      if (v0>v) x1=xx; else x0=xx;
    }
  }
  xmin=-x1;
  x0=0; x1=MAX_X;
  dbl_a[0]=da0+((double)x1)*dd;
  v0=sqrt(ifs(dbl_a,find_best_skewness(dbl_a,skewness)));
  if (v0>v) {
    while (x1-x0>1) {
      xx=(x0+x1)/2;
      dbl_a[0]=da0+((double)xx)*dd;
      v0=sqrt(ifs(dbl_a,find_best_skewness(dbl_a,skewness)));
      if (v0>v) x1=xx; else x0=xx;
    }
  }
  xmax=x1;
  dbl_a[0]=da0;
}


void initsieve()
{
  int exc;
  int p;
  int i, j, co[6], md, mp, len, pk, fi, mi, l, k;
  int sieveevaldiff[6], fis[6];
  double dp, dpk, dlp, dlog;
  int J, step, help;

for (i=0; i<6; i++) { mpz_out_str(stdout,10,gmp_a[5-i]); printf(" "); }
printf("\n");
mpz_out_str(stdout,10,gmp_p); printf(" ");
mpz_out_str(stdout,10,gmp_d); printf("\n\n");

  sievelen=SIEVELEN;
  if (2*(xmax-xmin)<sievelen) {
    sievelen/=2;
    exc=sievelen-(xmax-xmin); exc/=2;
    xmin-=exc;
    xmax=sievelen-1+xmin;
    nsubsieves=1;
  } else {
    nsubsieves=(xmax-xmin)/sievelen+1;
    if ((nsubsieves&1)==0) nsubsieves++;  /* odd number */
    exc=nsubsieves*sievelen-(xmax-xmin); exc/=2;
    xmin-=exc;
    xmax=nsubsieves*sievelen-1+xmin;
  }
  len=0;
  tval=0;
  for (k=0; k<NAFF_PRIMES; k++) {
    p=primes[k];
    dp=(double)p; dlp=(double)(log(dp));
    pk=p;
    while (pk<MAX_PRIME_AFF) {
      dpk=(double)pk; dlog=dlp/dpk*dp/(dp+1.);
      limit+=dlog;
      for (i=0; i<=5; i++) co[i]=mpz_mod_ui(gmp_help1,gmp_a[i],pk);
      md=mpz_mod_ui(gmp_help1,gmp_d,pk);
      mp=mpz_mod_ui(gmp_help1,gmp_p,pk);
      if (pk<8) {
        for (i=0; i<pk; i++) {
          fi=0; for (j=0; j<=5; j++) { fi*=i; fi+=co[5-j]; fi%=pk; }
          mi=md-((i*mp)%pk); if (mi<0) mi+=pk;
          if (mi==0) {
            if (fi==0) {
              limit-=dlog;   /* p^k | fi+j*(i*P-d) for all j */
              tval+=dlog;
            }
          } else {
            l=0; step=pk;
            while ((mi%p==0) && (fi%p==0)) {
              l++; mi/=p; fi/=p; step/=p;
              if (l>32) Schlendrian("initsieve\n");
            }
            if (mi%p) {
              J=divmod(fi,mi,step);
              pl[len].primepower=pk;
              pl[len].rooti=i;
              pl[len].J=J;
              pl[len].step=step;
              help=(J-ymin*i-xmin)%step;
              if (help<0) help+=step;
              pl[len].start=help;
              help=(-i+nsubsieves*sievelen)%step;
              if (help<0) help+=step;
              pl[len].lineinc=help;
              pl[len].value=dlog;
              len++;
            } /* otherwise no solutions */
          }
        }
      } else {
        for (i=0; i<6; i++) {
          fi=0; for (j=0; j<=5; j++) { fi*=i; fi+=co[5-j]; fi%=pk; }
          fis[i]=fi;
        }
        for (i=0; i<6; i++) {
          sieveevaldiff[i]=fis[0];
          for (j=0; j<5-i; j++) fis[j]=fis[j+1]-fis[j];
        }
        for (i=0; i<pk; i++) {
          mi=md-((i*mp)%pk); if (mi<0) mi+=pk;
          fi=sieveevaldiff[0];
          for (j=0; j<5; j++) sieveevaldiff[j]+=sieveevaldiff[j+1];
          for (j=0; j<5; j++) sieveevaldiff[j]%=pk;
          if (mi==0) {
            if (fi==0) limit-=dlog;   /* p^k | fi+j*(i-m) for all j */
          } else {
            l=0; step=pk;
            while ((mi%p==0) && (fi%p==0)) {
              l++; mi/=p; fi/=p; step/=p;
              if (l>32) Schlendrian("initsieve\n");
            }
            if (mi%p) {
              J=divmod(fi,mi,step);
              pl[len].primepower=pk;
              pl[len].rooti=i;
              pl[len].J=J;
              pl[len].step=step;
              help=(J-ymin*i-xmin)%step;
              if (help<0) help+=step;
              pl[len].start=help;
              help=(-i+nsubsieves*sievelen)%step;
              if (help<0) help+=step;
              pl[len].lineinc=help;
              pl[len].value=dlog;
              len++;
            } /* otherwise no solutions */
          }
        }
      }
      pk*=p;
    }
  }
/*
  for (i=0; i<len; i++) {
    printf("%d %d %d %d %d %d %f\n",pl[i].primepower,pl[i].rooti,pl[i].J,pl[i].step,pl[i].start,pl[i].lineinc,pl[i].value);
  }
*/
/* provisorisch */
  for (i=0; i<len; i++) {
    pl[i].v=(unsigned short)(pl[i].value*1000.+0.5);
    if (pl[i].v==0) pl[i].v=1;
  }
  cut=(unsigned short)(limit*1000.+0.5);
  if (cut==0) cut=1;

  primelistlen=len;
}


void prepare_sieve()
{
  int l, i0, i1, i, k, p, pk, ii, st;
  unsigned short *us_ptr;

  l=0; i0=0;
  for (k=0; k<NAFF_PRIMES; k++) {
    p=primes[k];
    for (i=i0; i<primelistlen; i++) if (pl[i].primepower==p) break;
    if (i==primelistlen) continue;  /* p not in pl-list */
    i0=i;
    for (i=i0+1; i<primelistlen; i++) if (pl[i].primepower%p) break;
    i1=i;   /* [i0,i1[: p-part of pl-list */
    pk=p;
    for (i=i0; i<i1; i++)
      if (pl[i].primepower>pk)
        pk=pl[i].primepower;         /* maximal primepower */

    for (i=0; i<pk; i++) prep_p_begin[l][i]=0;
    us_ptr=(unsigned short *)(prep_p_begin[l]);
    for (i=i0; i<i1; i++) {
      st=pl[i].step; ii=pl[i].start;
      while (ii<2*pk) {
        us_ptr[ii]+=pl[i].v; ii+=st;
      }
    }

#if 0
/* compress data, a bit slower */
    while (pk>p) {
      st=pk/p;
      for (i=0; i<pk-st; i++)
        if (prep_p_begin[l][i]!=prep_p_begin[l][i+st]) break;
      if (i==pk-st) pk=st; else break;
    }
#endif
    for (i=0; i<pk; i++) prep_p_begin[l][i+pk]=prep_p_begin[l][i];

    prep_p[l]=prep_p_begin[l];
    prep_p_len[l]=pk;
    l++;
  }
  prep_len=l;

#if 0
  for (i=0; i<20; i++) {
printf("%d: (%d,%d) %d %d %d %d %u %f\n",i,pl[i].primepower,pl[i].rooti,pl[i].J,pl[i].step,pl[i].start,pl[i].lineinc,pl[i].v,pl[i].value);
  }

printf("len:%d\n",l);
for (i=0; i<l; i++) {
  printf("  p: %d;    ",prep_p_len[i]);
  for (k=0; k<prep_p_len[i]; k++) {
    printf("%u ",prep_p_begin[i][k]);
  }
  printf("\n");
}
complain("debug\n");
#endif
}


void finish_sieve(int nsubsieves, int sievelen) /* provisorisch */
{
  int i, ind, add;

  for (i=0; i<primelistlen; i++) {
    ind=pl[i].start; add=pl[i].step;
    ind+=add;
    ind-=((sievelen*nsubsieves)%add);
    ind%=add;
    pl[i].start=ind;
  }
}

/* CJM, 2/28/05, prototype added: */
void asm_root_sieve8(unsigned int **p1, unsigned int *p2, int l1, unsigned int *p4, int l2);


void sieve_new(int len)
{
  int   i, len2;
  unsigned int *ul_sv;
#ifndef HAVE_ASM_INTEL
  int ind, end;
  unsigned int *ptr, *ptrbegin, *ptrend;
#endif

  len2=len/2; ul_sv=(unsigned int *)sievearray;
  memset(sievearray,0,len*sizeof(*sievearray));
  for (i=0; i<prep_len; i++) {
#ifdef HAVE_ASM_INTEL
#if 0
    asm_root_sieve(&(prep_p[i]),prep_p_begin[i],prep_p_len[i],ul_sv,len2);
#else
    asm_root_sieve8(&(prep_p[i]),prep_p_begin[i],prep_p_len[i],ul_sv,len/4);
#endif
#else
    ptr=prep_p[i]; ptrbegin=prep_p_begin[i];
    ptrend=ptrbegin+prep_p_len[i];
    ind=0; end=len2-prep_p_len[i];
    while (ptr<ptrend) ul_sv[ind++]+=*ptr++;
    while (ind<end) {
      for (ptr=ptrbegin; ptr<ptrend; ptr++, ind++) ul_sv[ind]+=*ptr;
    }
    for (ptr=ptrbegin; ind<len2; ptr++, ind++) ul_sv[ind]+=*ptr;
    prep_p[i]=ptr;
#endif
  }
}


void advancesieve()
{
  int i, s, add;

  for (i=0; i<primelistlen; i++) {
    s=pl[i].start+pl[i].lineinc;
    add=pl[i].step;
    if (s>=add) s-=add;
    pl[i].start=s;
  }
}


void do_sieving()
{
  int x, y, i, j;
  double lim0;
  double sk, v, dbl_sv[6];

#ifdef ZEIT
zeita(1);
zeita(5);
#endif
  compute_sieving_area(); /*ymin=-22; ymax=-22; xmin=-326415; xmax=xmin+4095;*/
  lim0=log(pol_norm);
  limit=lim0-log_max_norm_2+alpha_proj;
  initsieve();
  if (verbose>1) printf("area: [%d,%d]x[%d,%d]\n",xmin,xmax,ymin,ymax);
  if (verbose>2) printf("limit: %f ",limit);
  if (limit<0.) {   /* accept any (x,y) in sieving area */
#ifdef ZEIT
zeitb(5);
zeita(3);
#endif
    for (y=ymin; y<=ymax; y++)
      for (x=xmin; x<=xmax; x++)
        check(x,y);
#ifdef ZEIT
zeitb(3);
zeitb(1);
#endif
    return;
  }
  dbl_p=mpz_get_d(gmp_p); dbl_d=mpz_get_d(gmp_d);
  for (i=0; i<=5; i++) dbl_sv[i]=mpz_get_d(gmp_a[i]);
  dbl_sv[2]+=((double)ymin)*dbl_p;
  dbl_sv[1]-=((double)ymin)*dbl_d;
  sk=skewness;

printf("norm: %.3e  alpha_proj: %.3f  skewness: %.3f\n",pol_norm,alpha_proj,sk);
printf("limit: %f, lim0: %f   cut: %d\n",limit,lim0,(int)cut);

#ifdef ZEIT
zeitb(5);
#endif
  for (y=ymin; y<=ymax; y++) {
    sk=find_best_skewness(dbl_sv,sk);
    v=ifs(dbl_sv,sk);
    dbl_sv[2]+=dbl_p; dbl_sv[1]-=dbl_d;
    v=log(v)/2.;
    cut=(unsigned short)((limit+v-lim0)*1000+0.5);
    if (cut==0) cut=1;
    if (verbose>2) printf("y: %d, lim: %f cut: %d,  sk: %.2f\n",y,v,(int)cut,sk);

    x=xmin;
#ifdef ZEIT
zeita(5);
#endif
    prepare_sieve();
#ifdef ZEIT
zeitb(5);
#endif
    for (i=0; i<nsubsieves; i++) {
#ifdef ZEIT
zeita(2);
#endif
      sieve_new(sievelen);
#ifdef ZEIT
zeitb(2);
zeita(3);
#endif
      if (verbose>2) {
        for (j=0; j<sievelen; j++)
          if (sievearray[j]>=cut) {
            printf("sv: %d, cut: %d  at x: %d  ",sievearray[j],(int)cut,x+j);
            check(x+j,y);
          }
      } else {
        for (j=0; j<sievelen; j++) if (sievearray[j]>=cut) check(x+j,y);
      }
      x+=sievelen;
#ifdef ZEIT
zeitb(3);
#endif
    }
    finish_sieve(nsubsieves,sievelen);
    if (y<ymax) advancesieve();
    if (verbose>1) printf(".");
  }
#ifdef ZEIT
zeitb(1);
#endif
}

/* ----------------------------------------------- */

void init_optimize()
{
  int i;

  mpz_init(gmp_m);
  mpz_init(gmp_mb);
  mpz_init(gmp_p);
  mpz_init(gmp_d);
  mpz_init(gmp_lina[0]); mpz_init(gmp_lina[1]);
  mpz_init(gmp_linb[0]); mpz_init(gmp_linb[1]);
  for (i=0; i<6; i++) mpz_init(gmp_a[i]);
  for (i=0; i<6; i++) mpz_init(gmp_b[i]);
  mpz_init(gmp_help1);
  mpz_init(gmp_help2);
  mpz_init(gmp_help3);
  mpz_init(gmp_help4);
  st_alpha=compute_standard_alpha();
  init_assess(bound0,bound1,area,p_bound);
}


void optimize()
{
  int err;

  find_m_name();
  open_m();
  while (1) {
    err=read_a5pd();
    if (!err) break;
    if (err<0) {
      printf("ignoring line: %s\n",input_line);
      continue;
    }
#if 0
    printf("line: %s\n",input_line);
#endif
    if (!pol_expand()) {
      mpz_out_str(stdout,10,gmp_a[5]);
      fprintf(stderr,"expand failed\n");
      continue;
    }
#ifdef ZEIT
zeita(4);
#endif
    optimize_1();
    compute_proj_alpha();
#ifdef ZEIT
zeitb(4);
#endif
    if (verbose>2) printf("proj_alpha: %.3f\n",alpha_proj);
    if (pol_norm*exp(alpha_proj)>max_norm_1) continue;
if (0)
{
int i;
printf("pol=");
for (i=0; i<6; i++) { mpz_out_str(stdout,10,gmp_a[5-i]); printf("*x^%d+ ",5-i); }
printf("0 \n");
printf("m= "); mpz_out_str(stdout,10,gmp_m); printf("\n");
printf("s=%f\n",skewness);
printf("value=%f (%f)\n",e_value-alpha_proj,alpha_proj);
}
    do_sieving();
  }
  close_m();
}



int main(int argc, char **argv)
{
#ifdef ZEIT
initzeit(25); zeita(0);
#endif
  printf("%s\n", START_MESSAGE);

  setbuf(stdout,NULL);
  get_options(argc,argv);
  read_data();
  open_outputfile();
  init_optimize();
  optimize();
  close_outputfile();

#ifdef ZEIT
zeitb(0);
  printf("\nTiming:");
  printf("\nTotal            "); printzeit(0);
  printf("\n  1. Optimization  "); printzeit(4);
  printf("\n  Sieve all        "); printzeit(1);
  printf("\n    Init sieve       "); printzeit(5);
  printf("\n    Sieve            "); printzeit(2);
  printf("\n    eval             "); printzeit(3);
  printf("\n      gmp/alpha        "); printzeit(6); printzeit(7);
  printf("\n      2. Optimization  "); printzeit(8);
  printf("\n      3. Optimization  "); printzeit(9);
  printf("\n        polroots         "); printzeit(10);
  printf("\n        murphy-e sum     "); printzeit(11);
  printf("\n");
#endif
  return 0;
}

