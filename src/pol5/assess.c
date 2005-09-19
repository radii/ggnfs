/* assess.c 

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

#include <unistd.h>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <time.h>
#include <sys/timeb.h>
#include <stdlib.h>
#include <stdio.h>
#include "gmp.h"
#include "if.h"
#include <limits.h>
#include "fnmatch.h"
#include <string.h>


#define NSUMMANDS        1000
/*
#define LATTICE
*/
#define NSMALLPRIMES      6542

extern int verbose;

static double assess_bound0, assess_bound1, assess_area;
static double *assess_optima;

static unsigned int assess_primes[NSMALLPRIMES];
static unsigned int assess_prime_bound;
static double assess_alpha_max, assess_rat;
static unsigned int *assess_mod, *assess_root, *assess_coeffmod;
static int assess_mod_len, assess_root_len, assess_coeffmod_len;

/* ---------------------------------------------------------- */

unsigned int invert(unsigned int a, unsigned int p)  /* 0<b<p */
{
  unsigned int v1=0, v2=1, q, b=a, oldp=p;

  if (a==0) complain("cannot invert 0\n");
  while (b>1) {
    p-=b; v1+=v2;
    if (p>=b) {
      p-=b; v1+=v2;
      if (p>=b) {
        q=p/b;
        v1+=q*v2;
        q*=b;
        p-=q;
      }
    }
    if (p<=1) { v2=oldp-v1; break; }
    b-=p; v2+=v1;
    if (b>=p) {
      b-=p; v2+=v1;
      if (b>=p) {
        q=b/p;
        v2+=q*v1;
        q*=p;
        b-=q;
      }
    }
  }
  if ((v2*a-1)%oldp) complain("invert %u %u %u \n",v2,a,oldp);
  return v2;
}


void pol_mul_lin(unsigned int p, int deg, unsigned int *pol, unsigned int *modpol, unsigned int a)
{
  int i;
  unsigned int coeff, pol2[7];

  pol2[0]=0;
  for (i=0; i<=deg; i++) pol2[i+1]=pol[i];
  for (i=0; i<=deg; i++) pol2[i]+=((pol[i]*a)%p);
  if (pol2[deg+1]) {
    coeff=pol2[deg+1];
    for (i=0; i<=deg; i++) pol2[i]+=((coeff*modpol[i])%p);
  }
  for (i=0; i<=deg; i++) pol[i]=pol2[i]%p;
}


void pol_sqr(unsigned int p, int deg, unsigned int *pol, unsigned int *modpol)
{
  int i, j;
  unsigned int coeff, prod[12];

  for (i=0; i<=2*deg; i++) prod[i]=0;
  for (i=0; i<=deg; i++) prod[2*i]=(pol[i]*pol[i])%p;
  for (i=0; i<deg; i++)
    for (j=i+1; j<=deg; j++)
      prod[i+j]+=(((pol[i]*pol[j])%p)<<1);
  for (i=deg; i>=0; i--)
    if (prod[i+deg]) {
      coeff=prod[i+deg];
      for (j=0; j<=deg; j++) prod[i+j]+=((coeff*modpol[j])%p);
    }
  for (i=0; i<=deg; i++) pol[i]=prod[i]%p;
}


void pol_exp(unsigned int p, int deg, unsigned int *pol, unsigned int *modpol, unsigned int a, unsigned int ex)
{
  int i, n;
  unsigned int ma, pol2[6], inv;

  if (modpol[deg]!=p-1) {
    if (modpol[deg]==0) complain("pol_exp\n");
    inv=invert(p-modpol[deg],p);
    for (i=0; i<=deg; i++) { modpol[i]*=inv; modpol[i]%=p; }
  }
  pol2[0]=1; for (i=1; i<=deg; i++) pol2[i]=0;
  n=31; ma=1UL<<31;
  while (n>0) {
    if (ex&ma) break;
    ma>>=1; n--;
  }
  for (; n>=0; n--) {
    pol_sqr(p,deg,pol2,modpol);
    if (ex&ma) pol_mul_lin(p,deg,pol2,modpol,a);
    ma>>=1;
  }
  for (i=0; i<=deg; i++) pol[i]=pol2[i];
}


int pol_gcd(unsigned int p, int deg, unsigned int *res , unsigned int *pol1, unsigned int *pol)
{
  int i, d2, d3, diff;
  unsigned int h, leadinv, pol2[6], pol3[6];

  for (i=0; i<=deg; i++) { pol2[i]=pol[i]; pol3[i]=pol1[i]; }
  d2=deg; d3=deg;
  while (d2>=0) if (pol2[d2]) break; else d2--;
  while (d3>=0) if (pol3[d3]) break; else d3--;
  if (d2==-1) { for (i=0; i<=d3; i++) res[i]=pol3[i]; return d3; }
  if (d3==-1) { for (i=0; i<=d2; i++) res[i]=pol2[i]; return d2; }
  h=p-pol2[d2]; leadinv=invert(h,p);
  for (i=0; i<=d2; i++) pol2[i]=((leadinv*pol2[i])%p);
  h=p-pol3[d3]; leadinv=invert(h,p);
  for (i=0; i<=d3; i++) pol3[i]=((leadinv*pol3[i])%p);
  while (1) {
    if (d3==-1) { for (i=0; i<=d2; i++) res[i]=pol2[i]; return d2; }
    if (d2==-1) { for (i=0; i<=d3; i++) res[i]=pol3[i]; return d3; }
    if (d3>=d2) {
      diff=d3-d2;
      for (i=0; i<=d2; i++) {
        pol3[i+diff]+=p; pol3[i+diff]-=pol2[i];
        pol3[i+diff]%=p;
      }
      while (d3>=0) if (!pol3[d3]) d3--; else break;
      if (d3==-1) continue;
      if (pol3[d3]!=p-1) {
        h=p-pol3[d3]; leadinv=invert(h,p);
        for (i=0; i<=d3; i++) pol3[i]=((leadinv*pol3[i])%p);
      }
    } else {  /* d3<d2 */
      diff=d2-d3;
      for (i=0; i<=d3; i++) {
        pol2[i+diff]+=p; pol2[i+diff]-=pol3[i];
        pol2[i+diff]%=p;
      }
      while (d2>=0) if (!pol2[d2]) d2--; else break;
      if (d2==-1) continue;
      if (pol2[d2]!=p-1) {
        h=p-pol2[d2]; leadinv=invert(h,p);
        for (i=0; i<=d2; i++) pol2[i]=((leadinv*pol2[i])%p);
      }
    }
  }
  Schlendrian("pol_gcd\n");
  return 0;
}


void pol_div_lin(unsigned int p, int *deg, unsigned int *coeff, unsigned int r)
{
  int i, d;
  unsigned int h, res;

  d=*deg; h=coeff[d];
  for (i=d-1; i>=0; i--) {
    res=h; h*=r; h+=coeff[i]; h%=p;
    coeff[i]=res;
  }
  if (h) complain("pol_div: not divisible\n");
  *deg=d-1;
}


unsigned int find_root(unsigned int p, int deg, unsigned int *coeff)
{
  unsigned int rpol[6], rpol1[6], rpol2[6], res, inv, a, p2;
  int i, rd, d;

  for (i=0; i<=deg; i++) rpol[i]=coeff[i];
  if (p==2) {
    if (!rpol[0]) return 0;
    a=0;
    for (i=0; i<=deg; i++) a+=rpol[i];
    if ((a&1)==0) return 1;
    return 2; /* no root */
  }
  rd=deg; p2=(p-1)/2;
  if (rd>1) {
    pol_exp(p,rd,rpol1,rpol,0,p);
    if (rpol1[1]) rpol1[1]--; else rpol1[1]=p-1;
    d=0; for (i=0; i<=rd; i++) if (rpol1[i]) d=1;
    if (d) {
      d=pol_gcd(p,rd,rpol,rpol,rpol1);
      rd=d;
    } /* else rpol has rd simple roots */
    if (rd==0) return p; /* no root */
  }
  while (rd>1) {
    a=((unsigned int)rand())%p;
    pol_exp(p,rd,rpol1,rpol,a,p2);
    if (rpol1[0]) rpol1[0]--; else rpol1[0]=p-1;
    d=pol_gcd(p,rd,rpol2,rpol1,rpol);
    if (d) {
      rd=d; for (i=0; i<=rd; i++) rpol[i]=rpol2[i];
    }
  }
  if (rd==0) complain("find_root: degree=0?\n");
  inv=invert(rpol[1],p);
  res=inv*rpol[0]; res%=p; if (res) res=p-res;
  return res;
}


unsigned int pol_eval(unsigned int p, int deg, unsigned int *coeff, unsigned int r)
{
  int i;
  unsigned int res;

  res=0;
  for (i=deg; i>=0; i--) { res*=r; res+=coeff[i]; res%=p; }
  return res;
}


int find_pol_roots_homog(unsigned int p, int deg, unsigned int *polmod, unsigned int *roots)
{
  unsigned int r;
  int n, d, i;

  if (assess_mod_len<deg+1) {
    assess_mod_len=deg+1;
    assess_mod=(unsigned int *)xrealloc(assess_mod,assess_mod_len*sizeof(unsigned int));
  }
  for (i=0; i<=deg; i++) assess_mod[i]=polmod[i]%p;
  d=deg; n=0; 
  while ((d) && (assess_mod[d]==0)) { d--; roots[n++]=p; }
  if (d==0) {
    if (assess_mod[0]==0) complain("find_pol_roots_homog: zero pol\n");
    return n;
  }
  while (d) {
    r=find_root(p,d,assess_mod);
    if (r==p) break;
    if (pol_eval(p,d,assess_mod,r)) complain("find_pol_roots_homog\n");
    do {
      pol_div_lin(p,&d,assess_mod,r); roots[n++]=r;
    } while (pol_eval(p,d,assess_mod,r)==0);
  }
  return n;
}

/* ---------------------------------------------------------- */

#define MAX_BF      65535

double brute_force(unsigned int p, int deg, mpz_t *gmp_coeff, unsigned int r)
{
  unsigned int qmax, h, q, q0;
  unsigned int lifts[32], powers[32];
  int ind, indmax, j;
  double dp, lo, res;

  qmax=1; indmax=0;
  while (qmax<MAX_BF/p) { qmax*=p; powers[indmax++]=qmax; }
  if (assess_mod_len<deg+1) {
    assess_mod_len=deg+1;
    assess_mod=(unsigned int *)xrealloc(assess_mod,assess_mod_len*sizeof(unsigned int));
  }
  for (j=0; j<=deg; j++) assess_mod[j]=mpz_fdiv_ui(gmp_coeff[j],qmax);
  dp=(double)p; lo=log(dp)*dp/(dp+1);
  ind=0; lifts[0]=r; res=0.;
  while (1) {
    q=powers[ind]; res+=lo/((double)q);
    if (ind<indmax-1) { /* lift root */
      h=lifts[ind];
      while (h<q*p) {
        if (pol_eval(q*p,deg,assess_mod,h)==0) break;
        h+=q;
      }
      if (h<q*p) {
        ind++; lifts[ind]=h;
        continue;
      }
    } /* find next lift or go back */
    while (ind>=0) {
      q=powers[ind]; q0=q/p;
      h=lifts[ind]+q0;
      while (h<q) {
        if (pol_eval(q,deg,assess_mod,h)==0) break;
        h+=q0;
      }
      lifts[ind]=h;
      if (h<q) break; /* found another zero */
      ind--;
    }
    if (ind>0) continue;
    break; /* finished */
  }
  return res;
}


double brute_force_proj(unsigned int p, int deg, mpz_t *gmp_coeff)
{
  unsigned int qmax, h, q, q0;
  unsigned int lifts[32], powers[32];
  int ind, indmax, j;
  double dp, lo, res;

  qmax=1; indmax=0;
  while (qmax<MAX_BF/p) { qmax*=p; powers[indmax++]=qmax; }
  if (assess_mod_len<deg+1) {
    assess_mod_len=deg+1;
    assess_mod=(unsigned int *)xrealloc(assess_mod,assess_mod_len*sizeof(unsigned int));
  }
  for (j=0; j<=deg; j++) assess_mod[j]=mpz_fdiv_ui(gmp_coeff[deg-j],qmax);
  dp=(double)p; lo=log(dp)*dp/(dp+1);
  ind=0; lifts[0]=0; res=0.;
  while (1) {
    q=powers[ind]; res+=lo/((double)q);
    if (ind<indmax-1) { /* lift root */
      h=lifts[ind];
      while (h<q*p) {
        if (pol_eval(q*p,deg,assess_mod,h)==0) break;
        h+=q;
      }
      if (h<q*p) {
        ind++; lifts[ind]=h;
        continue;
      } /* did not find lift */
    } /* find next lift or go back */
    while (ind>=0) {
      q=powers[ind]; q0=q/p;
      h=lifts[ind]+q0;
      while (h<q) {
        if (pol_eval(q,deg,assess_mod,h)==0) break;
        h+=q0;
      }
      lifts[ind]=h;
      if (h<q) break; /* found another zero */
      ind--;
    }
    if (ind>0) continue;
    break; /* finished */
  }
  return res;
}

/* ---------------------------------------------------------- */

void compute_coeffmod(unsigned int p, int deg, mpz_t *coeff)
{
  int i;

  if (assess_coeffmod_len<deg+1) {
    assess_coeffmod_len=deg+1;
    assess_coeffmod=(unsigned int *)xrealloc(assess_coeffmod,assess_coeffmod_len*sizeof(unsigned int));
  }
  for (i=0; i<=deg; i++) assess_coeffmod[i]=mpz_fdiv_ui(coeff[i],p);
}


/* returns 1 if alpha>alpha_targ */
int compute_alpha(double *alpha, int deg, unsigned int **coeffmod, mpz_t *gmp_coeff, double alpha_targ)
{
  double al, dp;
  unsigned int p, r;
  int i, j, n;
  double al_prev, max_rest, rat_rest;

  al=0.; max_rest=assess_alpha_max; rat_rest=assess_rat;
  if (assess_root_len<deg+1) {
    assess_root_len=deg+1;
    assess_root=(unsigned int *)xrealloc(assess_root,assess_root_len*sizeof(unsigned int));
  }
  for (i=0; i<NSMALLPRIMES; i++) {
    al_prev=al;
    p=assess_primes[i];
    if (p>assess_prime_bound) break;
    if (coeffmod==NULL) {
      compute_coeffmod(p,deg,gmp_coeff);
      n=find_pol_roots_homog(p,deg,assess_coeffmod,assess_root);
    } else n=find_pol_roots_homog(p,deg,coeffmod[i],assess_root);
    j=0; dp=(double)p;
    while (j<n) {
      r=assess_root[j];
      if ((j+1<n) && (r==assess_root[j+1])) { /* multiple zero */
        if (r==p) al-=brute_force_proj(p,deg,gmp_coeff);
        else al-=brute_force(p,deg,gmp_coeff,r);
      } else {
        al-=dp*log(dp)/(dp*dp-1);
      }
      j++;
      while ((j<n) && (assess_root[j]==r)) j++;
    }
    al+=log(dp)/(dp-1);
    if (verbose>3) {
      printf("%u   %.3f %.3f:  ",p,al,al-al_prev);
      for (j=0; j<n; j++) printf("%u ",assess_root[j]);
      printf("\n");
    }
    if (p>100) {
      max_rest-=dp*log(dp)/(dp*dp-1);
      rat_rest-=log(dp)/(dp-1);
      if (al+rat_rest-max_rest*((double)deg)>alpha_targ) {
        return 1;
      }
    }
  }
  if (verbose>3) printf("\n");
  *alpha=al;
  return 0;
}


void compute_alpha_exact(double *alpha, int deg, unsigned int **coeffmod, mpz_t *gmp_coeff, unsigned int pb)
{
  double al, dp;
  unsigned int p, r;
  int i, j, n;
  double al_prev;

  al=0.;
  if (assess_root_len<deg+1) {
    assess_root_len=deg+1;
    assess_root=(unsigned int *)xrealloc(assess_root,assess_root_len*sizeof(unsigned int));
  }
  for (i=0; i<NSMALLPRIMES; i++) {
    al_prev=al;
    p=assess_primes[i];
    if (p>pb) break;
    if (coeffmod==NULL) {
      compute_coeffmod(p,deg,gmp_coeff);
      n=find_pol_roots_homog(p,deg,assess_coeffmod,assess_root);
    } else n=find_pol_roots_homog(p,deg,coeffmod[i],assess_root);
    j=0; dp=(double)p;
    while (j<n) {
      r=assess_root[j];
      if ((j+1<n) && (r==assess_root[j+1])) { /* multiple zero */
        if (r==p) al-=brute_force_proj(p,deg,gmp_coeff);
        else al-=brute_force(p,deg,gmp_coeff,r);
      } else {
        al-=dp*log(dp)/(dp*dp-1);
      }
      j++;
      while ((j<n) && (assess_root[j]==r)) j++;
    }
    al+=log(dp)/(dp-1);
    if (verbose>3) {
      printf("%u   %.3f %.3f:  ",p,al,al-al_prev);
      for (j=0; j<n; j++) printf("%u ",assess_root[j]);
      printf("\n");
    }
  }
  if (verbose>3) printf("\n");
  *alpha=al;
}

/* ---------------------------------------------------------- */

int compare_doubles(const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

#include "dickman.tab"
#define EULER_C   0.577215664901532

double dickman(double x)
{
  int i;
  double h;

  if (x<=1.) return 1.;
  if (x>=20.) return 0.;
  i=(int)(20*x);
  h=20*x-(double)i;
  return (h*rho_table[i]+(1-h)*rho_table[i-1]);
}


double prob(double r, double b)
{
  double u, logr;

  logr=log(fabs(r));
  u=logr/log(b);
  return dickman(u);
}


void murphy_en(double *me, int deg0, double *dbl_coeff0, int deg1, double *dbl_coeff1, double alpha0, double alpha1, double skewness, int nsm)
{
  int i, j;
  double x, y, sx, sy, theta;
  double v0, v1, xp, e, al0, al1;
  double e0, left, right, theta_left, theta_right, theta_len, theta_inc;
  int interval, nop, nsum;
  int deg[2];
  double *dbl_coeff[2];

  sx=sqrt(assess_area*skewness); sy=sx/skewness;
#ifdef LATTICE
  al0=exp(alpha0)/assess_bound0; /* assuming that special-q in lattice sieving
                                  is of magnitude assess_bound0 */
#else
  al0=exp(alpha0);
#endif
  al1=exp(alpha1);
  deg[0]=deg0; deg[1]=deg1;
  dbl_coeff[0]=dbl_coeff0; dbl_coeff[1]=dbl_coeff1;
  nop=find_optima(deg,dbl_coeff,skewness,&assess_optima);
  qsort(assess_optima,nop,sizeof(double),compare_doubles);

  e=0;
  for (interval=1; interval<nop; interval++) {
    left=assess_optima[interval-1]; right=assess_optima[interval];
    theta_left=atan(left); theta_right=atan(right);
    theta_len=theta_right-theta_left;
    nsum=(int)(((double)nsm)/M_PI*theta_len);
    if (nsum<10) nsum=10;
    theta_inc=theta_len/(double)(nsum);
    e0=0; theta=theta_left+theta_inc/2.;
    for (i=0; i<nsum; i++) {
      theta+=theta_inc;
      x=sx*sin(theta); y=sy*cos(theta);
      v0=dbl_coeff0[0]; xp=1.;
      for (j=1; j<=deg0; j++) { xp*=x; v0*=y; v0+=xp*dbl_coeff0[j]; } 
      v1=dbl_coeff1[0]; xp=1.;
      for (j=1; j<=deg1; j++) { xp*=x; v1*=y; v1+=xp*dbl_coeff1[j]; }
      e0+=prob(al0*v0,assess_bound0)*prob(al1*v1,assess_bound1);
    }
    e0/=nsum;
    e+=(e0*theta_len);
  }
  e/=M_PI;
  *me=e;
}


void murphy_e(double *me, int deg0, double *dbl_coeff0, int deg1, double *dbl_coeff1, double alpha0, double alpha1, double skewness)
{
  murphy_en(me,deg0,dbl_coeff0,deg1,dbl_coeff1,alpha0,alpha1,skewness,NSUMMANDS);

}

/* ----------------------------------------------- */

void init_assess(double b0, double b1, double area, unsigned int pb)
{
  int i;
  double dp;

  assess_bound0=b0;
  assess_bound1=b1;
  assess_area=area;
  assess_prime_bound=pb;
  prime_table_init();
  for (i=0; i<NSMALLPRIMES; i++)
    assess_primes[i]=get_next_prime();
  if (assess_primes[NSMALLPRIMES-1]!=65521) Schlendrian("init_assess\n");
  assess_alpha_max=0.; assess_rat=0.;
  for (i=0; i<NSMALLPRIMES; i++) {
    if (assess_primes[i]<100) continue;
    if (assess_primes[i]>assess_prime_bound) break;
    dp=(double)(assess_primes[i]);
    assess_alpha_max+=dp*log(dp)/(dp*dp-1);
    assess_rat+=log(dp)/(dp-1);
  }

  assess_mod_len=0; assess_mod=NULL;
  assess_root_len=0; assess_root=NULL;
  assess_coeffmod_len=0; assess_coeffmod=NULL;
  assess_optima=(double *)xmalloc(12*sizeof(double)); /* degree 5+1 !! */
}

