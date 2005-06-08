/* roots.c

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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "if.h"

extern void *xmalloc(size_t size);


static double **roots_pol_deriv, *roots_derivroots;
double *roots_optima;
static int roots_noptima, roots_deg;

static void init_find_roots(int deg, double *polcoeff)
{
  int i, j;

  roots_deg=deg;
  roots_optima=(double *)xmalloc(2*deg*sizeof(double));
  roots_noptima=0;
  roots_pol_deriv=(double **)xmalloc(deg*sizeof(double *));
  for (i=0; i<deg; i++)
    roots_pol_deriv[i]=(double *)xmalloc((deg+1)*sizeof(double));
  for (j=0; j<=deg; j++) roots_pol_deriv[0][j]=polcoeff[j];
  for (i=1; i<deg; i++) {
    for (j=0; j<deg; j++)
      roots_pol_deriv[i][j]=roots_pol_deriv[i-1][j+1]*(double)(j+1);
    roots_pol_deriv[i][deg]=0.;
  }
  roots_derivroots=(double *)xmalloc((deg+2)*sizeof(double));
}

static void exit_find_roots()
{
  int i;

  for (i=0; i<roots_deg; i++) free(roots_pol_deriv[i]);
  free(roots_pol_deriv);
  free(roots_derivroots);
}


static double roots_eval(int deriv, double x)
{
  double res, xn;
  int i;

  res=roots_pol_deriv[deriv][0]; xn=1.;
  for (i=1; i<=roots_deg-deriv; i++) {
    xn*=x;
    res+=roots_pol_deriv[deriv][i]*xn;
  }
  return res;
}


static int roots_find(int deriv, double d0, double d1, double *dr)
{
  double p0, p1, db, de, pm;
  volatile double dm; /* important for some compilers */

  p0=roots_eval(deriv,d0); p1=roots_eval(deriv,d1);
  if (p1==0.) { *dr=d1; return 1; }
  if (p0==0.) return 0;
  if (p0*p1>0.) return 1;
  db=d0; de=d1;
  while (de>db) {
    dm=(db+de)/2.;
    if ((dm==db) || (dm==de)) break;
    pm=roots_eval(deriv,dm);
    if (pm==0.) break;
    if (p0*pm>0.) { db=dm; p0=pm; } else { de=dm; p1=pm; }
  }
  *dr=dm;
  return 1;
}

#define DBL_MIN   -1e20
#define DBL_MAX    1e20

int find_roots(int deg, double *polcoeff, double **op)
{
  int n, deriv, i, j, nroots;
  double d, dd;

  if (deg<1) return 0;
  if (deg==1) {
    d=-polcoeff[0]/polcoeff[1];
    roots_optima=(double *)xmalloc(sizeof(double));
    *op=roots_optima;
    roots_optima[0]=d;
    return 1;
  }
  init_find_roots(deg,polcoeff);
  roots_derivroots[0]=DBL_MIN; roots_derivroots[1]=DBL_MAX;
  n=0;
  for (deriv=deg-1; deriv>=0; deriv--) {
    j=1; d=roots_derivroots[0];
    for (i=1; i<n+2; i++) {
      dd=roots_derivroots[i];
      if (roots_find(deriv,d,dd,roots_derivroots+j)) j++;
      d=dd;
    }
    if (roots_derivroots[j-1]<DBL_MAX) { roots_derivroots[j]=DBL_MAX; j++; }
    n=j-2;
  }
  for (i=0; i<n; i++) roots_optima[i]=roots_derivroots[i+1];
  nroots=n;

  i=0;
  while (i<nroots-1) {
    if (roots_optima[i]==roots_optima[i+1]) {
      for (j=i+1; j<nroots-1; j++) roots_optima[j]=roots_optima[j+1];
      nroots--;
      continue;
    }
    if (roots_optima[i]>roots_optima[i+1]) {
      d=roots_optima[i]; roots_optima[i]=roots_optima[i+1];
      roots_optima[i+1]=d;
      if (i>0) i--;
      continue;
    }
    i++;
  }
  exit_find_roots();
  *op=roots_optima;
  return nroots;
}


int find_optima(int *deg, double **coeff, double skewness, double **optima)
{
  int n, nop, s;
  double *op, *res;

/* zeros */
  nop=find_roots(deg[0],coeff[0],&res);
  if (!nop) res=xmalloc(sizeof(double));
  n=find_roots(deg[1],coeff[1],&op);
  if (n) {
    res=(double *)xrealloc(res,(n+nop)*sizeof(double));
    memcpy(res+nop,op,n*sizeof(double));
    free(op);
    nop+=n;
  }
  for (n=0; n<nop; n++) res[n]/=skewness;

/* extrema */
  for (s=0; s<2; s++) {
    double *deriv, *c, sp;
    int d, i;

    d=deg[s];
    deriv=(double *)xmalloc((d+1)*sizeof(double));
    c=(double *)xmalloc((d+1)*sizeof(double));
    sp=1.;
    for (i=0; i<=d; i++) { c[i]=coeff[s][i]*sp; sp*=skewness; }
    for (i=0; i<d; i++) deriv[i]=-c[i+1]*(double)(i+1);
    deriv[d]=0.;
    for (i=0; i<d; i++) deriv[i+1]+=c[i]*(double)(d-i);
    while (d) if (deriv[d]) break; else d--;
    n=find_roots(d,deriv,&op);
    if (n) {
      res=(double *)xrealloc(res,(n+nop)*sizeof(double));
      memcpy(res+nop,op,n*sizeof(double));
      free(op);
      nop+=n;
    }
    free(deriv); free(c);
  }
  nop+=2;
  res=(double *)xrealloc(res,nop*sizeof(double));
  res[nop-2]=-1e20;
  res[nop-1]=1e20;
  *optima=res;
  return nop;
}


