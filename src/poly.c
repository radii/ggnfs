/**************************************************************/
/* poly.c                                                     */
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

/********************************************************/
/* Functions for performing various tasks with our      */
/* number field defining polynomial, f.                 */
/* In particular, there are functions here to do        */
/* various types of evaluation, find the complex zeros, */
/* and find the base-m representation.                  */
/********************************************************/
/* Throughout,                                          */
/* 'f' is the polynomial f(x)                           */
/* 'F' is the homogenized version: F(x,y) = (y^d)f(x/y).*/
/* 'g' is the monic version, g(x) = c_d^{d-1}f(x/c_d).  */
/********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include "ggnfs.h"
#include "prand.h"

/******************************************************/
void poly_cp(poly_t dest, poly_t src)
/******************************************************/
{ int i;

  dest->degree = src->degree;
  for (i=dest->degree; i>=0; i--)
    dest->coef[i] = src->coef[i];
}

/******************************************************/
int poly_fixDeg(poly_t op)
/******************************************************/
{ int i = op->degree;

  while ((i>0) && (op->coef[i]==0))
    i--;
  op->degree = i;
  return i;
}
  
/******************************************************/    
int poly_Jacobi(poly_t f, poly_t g, s32 p)
/******************************************************/    
/* Jacobi symbol (f/g). If g is irreducible, this is  */
/* the quadratic character of 'f' in the finite field */
/* F_p[t]/<g>.                                        */
/* Bach & Shallit, Ch.6, p.142                        */
/******************************************************/    
{ poly_t u, v;
  s32   c;
  int    i, s;
  mpz_t  A, B;
  
  poly_cp(v, g);
  poly_mod(u, f, v, p);
  
  c = u->coef[u->degree];
  mpz_init_set_ui(A, c); mpz_init_set_ui(B, p);
  s = mpz_jacobi(A, B);
  mpz_clear(A); mpz_clear(B);
  if (v->degree%2==0)
    s = 1;
  if (u->degree==0)
    return s;

  c = inverseModP(c, p);
  for (i=0; i<=u->degree; i++)
    u->coef[i] = mulmod32(u->coef[i], c, p);
  if ( (((p-1)/2)&0x01) && (u->degree&0x01) && (v->degree&0x01))
    s = -s;
  return s*poly_Jacobi(v, u, p);
}

        
/***********************************************************************/
INLINE int poly_mul_modp(poly_t res, poly_t op1, poly_t op2, s32 p) 
/***********************************************************************/
/* Input: Polynomials 'op1' and 'op2' and a prime integer 'p'.         */
/* Output: Their product modulo 'p'.                                   */
/* res <-- op1*op2 (mod p)                                             */
/* Caveats:                                                            */
/* p must be < 2^{31} and res[] must have enough space allocated       */
/* for the result.                                                     */
/***********************************************************************/
{ int i, j, dd = op1->degree + op2->degree;
  s32 templ, tRes[2*MAXPOLYSIZE+1];

  for (i=dd; i>=0; i--)
    tRes[i] = 0;

  for (i=0; i<= op1->degree; i++) {
    for (j=0; j<= op2->degree; j++) {
      templ = mulmod32(op1->coef[i], op2->coef[j], p);
      tRes[i+j] = (tRes[i+j] + templ)%p;
    }
  }
  res->degree = dd;
  for (i=0; i<=dd; i++)
    res->coef[i] = tRes[i];
  poly_fixDeg(res);
  return 0;
}

/*********************************************************************/
INLINE int poly_mulmodpp(poly_t res, poly_t op1, poly_t op2, poly_t mod, s32 p)
/*********************************************************************/
/* Compute res <-- op1*op2 mod <p, mod(x)>                           */
/*********************************************************************/
{ poly_t tmp;
  

  if (mod->degree == 0) {
    res->coef[0]=0;
    res->degree = 0;
    return 1;
  }
  poly_mul_modp(tmp, op1, op2, p);
  poly_mod(res, tmp, mod, p);
  
  return 0;
}
  


/*********************************************************************/
int poly_mod(poly_t res, poly_t op, poly_t _mod, s32 p)
/*********************************************************************/
/* Input: Polynomials 'op' and 'mod' and a prime integer 'p'.        */
/* Output: The reduction of 'op' modulo 'mod'.                       */
/*   res <-- op (mod < mod, p> )                                     */
/*********************************************************************/
{ int  i;
  s32 lci, temp2l, lcl;
  poly_t tmp, mod;


  if (_mod->degree == 0) {
    res->coef[0]=0;
    res->degree = 0;
    return 1;
  }
  
  poly_cp(mod, _mod);
  poly_cp(tmp, op);

  /*********************/
  /* make "mod" monic  */
  /*********************/
  lci = mod->coef[mod->degree];
  if (lci != 1) {
    lci = inverseModP(lci, p);
    for (i=mod->degree; i>=0; i--) {
      mod->coef[i] = mulmod32(lci, mod->coef[i], p);
    }
  }
   
  while (tmp->degree >= mod->degree) {
    /************************************************/
    /* tmp <-- tmp - lc*mod*x^{*deg_tmp - deg_mod}  */
    /* where lc = tmp[*deg_tmp]                     */
    /************************************************/
    lcl = tmp->coef[tmp->degree];
    
    for (i=mod->degree; i>=0; i--) {
      temp2l = mulmod32(lcl, mod->coef[i], p);
      tmp->coef[tmp->degree -(mod->degree-i)] = 
          (p+tmp->coef[tmp->degree - (mod->degree-i)] - (temp2l%p))%p;
    }
    poly_fixDeg(tmp);
  }
  poly_cp(res, tmp);
  return 1;
}

/******************************************************************/
int poly_xpow_modpp(poly_t res, s32 n, poly_t mod, s32 p)
/******************************************************************/
/* Input: An exponent 'n', a polynomial 'mod', and a prime p.     */
/* Output: x^n mod <mod, p>.                                      */
/******************************************************************/
{ poly_t g;

  g->coef[0] = 0;
  g->coef[1] = 1;
  g->degree=1;

  return poly_pow_modpp(res, g, n, mod, p);
}


/******************************************************************/
int poly_pow_modpp(poly_t res, poly_t f, s32 n, poly_t mod, s32 p)
/******************************************************************/
/* res <-- f^n mod <mod, p>.                                      */
/******************************************************************/
{ poly_t g, h, temp;
  s32   remain;


  /* g <-- f, h <-- 1 */
  poly_cp(g, f);
  h->coef[0] = 1; h->degree = 0;

  remain=n;
  while (remain) {
    if (remain&0x01) {
      poly_mul_modp(temp, h, g, p);
      poly_mod(h, temp, mod, p);
    }
    remain = (remain>>1);
    poly_mul_modp(temp, g, g, p);
    poly_mod(g, temp, mod, p);
  }
  poly_cp(res, h);
  return res->degree;
}



/******************************************************************/
int poly_powmpz_mod(poly_t res, poly_t f, mpz_t n, poly_t mod, s32 p)
/******************************************************************/
/* res <-- f^n mod <mod, p>.                                      */
/******************************************************************/
{ poly_t g, h, temp;
  mpz_t  remain;

  mpz_init_set(remain, n);

  /* g <-- f, h <-- 1 */
  poly_cp(g, f);
  h->coef[0] = 1; h->degree = 0;
  
  while (mpz_sgn(remain)) {
    if (mpz_odd_p(remain)) {
      poly_mul_modp(temp, h, g, p);
      poly_mod(h, temp, mod, p);
    }
    mpz_div_2exp(remain, remain, 1);
    poly_mul_modp(temp, g, g, p);
    poly_mod(g, temp, mod, p);
  }
  
  poly_cp(res, h);
  mpz_clear(remain);
  return res->degree;
}

	

/******************************************************/    
int poly_gcd(poly_t g, poly_t h, s32 p)
/******************************************************/    
/* g <-- gcd(g,h)                                     */
/* g and h will both be modified.                     */
/******************************************************/    
{ poly_t r;

  /*************************************/
  /* Make sure the degrees are correct */
  /* or horrible things will happen.   */
  /*************************************/
  poly_fixDeg(h);
  poly_fixDeg(g);
    
  while ((h->degree>0) || (h->coef[h->degree])) {

    poly_mod(r, g, h, p);
    poly_cp(g, h);
    poly_cp(h, r);
  }
  poly_fixDeg(g);
  if (g->degree==0)
    g->coef[0] = 1;
  return 1;
}


/******************************************************/    
int poly_div(poly_t _q, poly_t _r, poly_t x, poly_t y, s32 p)
/******************************************************/    
/* Get q,r so that x = qy + r, deg(r) < deg(y).       */
/* i.e., compute x/y.                                 */
/******************************************************/    
{ int    i, j, n;
  s32   c, lr, ly_inv;
  poly_t q, r;	
  
  q->degree = MAX(x->degree - y->degree, 0);
  for (i=0; i<=q->degree; i++)
    q->coef[i]=0;
  poly_cp(r, x);

  poly_fixDeg(r);

  i=q->degree;
  ly_inv = inverseModP(y->coef[y->degree], p);
  while ((r->degree >= y->degree) && !((r->degree==0)&&(r->coef[0]==0))) {
    lr = r->coef[r->degree];
    c = mulmod32(ly_inv, lr, p);
    n = r->degree - y->degree;
    q->coef[n] = c;
    for (j=r->degree; j>=(r->degree - y->degree); j--) 
      r->coef[j] = (p+(r->coef[j] - mulmod32(y->coef[j-n], c, p)))%p;
    poly_fixDeg(r);
  }

  /* Why wasn't this here before? :*/
  poly_cp(_q, q); poly_cp(_r, r);

  return 0;
}
  


/******************************************************/    
int poly_inv(poly_t inv, poly_t f, poly_t m, s32 p)
/******************************************************/    
/* s <-- f^{-1} mod <p, m(x)>.                        */
/******************************************************/    
{ poly_t r, g, h, t, s1, s2, t1, t2, q, s;
  int    i;
  s32   c;
  
  poly_cp(g, f);
  poly_cp(h, m);

  /*************************************/
  /* Make sure the degrees are correct */
  /* or horrible things will happen.   */
  /*************************************/
  poly_fixDeg(h);
  poly_fixDeg(g);

  if ((h->degree==0) && (h->coef[0]==0)) {
    inv->coef[0]=1;
    inv->degree = 0;
    return 0;
  }
  
  s2->coef[0]=1; s2->degree=0;
  s1->coef[0]=0; s1->degree=0;
  t2->coef[0]=0; t2->degree=0;
  t1->coef[0]=1; t1->degree=0;

  while ((h->degree>0) || (h->coef[0])) {
    poly_div(q, r, g, h, p);

    /* s <-- s2 - q*s1. */ 
    poly_mul_modp(s, q, s1, p);   
    for (i=0; i<=s->degree; i++)
      s->coef[i] = (p-s->coef[i])%p;
    for (i=0; i<=s2->degree; i++)
      s->coef[i] = (s->coef[i] + s2->coef[i])%p;

    /* t <-- t2 - q*t1 */
    poly_mul_modp(t, q, t1, p);   
    for (i=0; i<=t->degree; i++)
      t->coef[i] = (p-t->coef[i])%p;
    for (i=0; i<=t2->degree; i++)
      t->coef[i] = (t->coef[i] + t2->coef[i])%p;
    
    /* g <-- h, h <-- r */
    poly_cp(g, h);
    poly_cp(h, r);
    /* s2 <-- s1, s1 <-- s */
    poly_cp(s2, s1);
    poly_cp(s1, s);
    /* t2 <-- t1, t1 <-- t */
    poly_cp(t2, t1);
    poly_cp(t1, t);
  }
  if (g->degree >0) {
    fprintf(stderr, "poly_inv(): Error - polynomial was not invertible!\n");
    exit(-1);
  }
  c = inverseModP(g->coef[0], p);

  inv->degree = s2->degree;
  for (i=0; i<=s2->degree; i++)
    inv->coef[i] = mulmod32(c, s2->coef[i], p);
  return 0;
}

	  
  

/***********************************************************/
int get_zeros_rec(s32 *zeros, int *numZeros, poly_t g, s32 p)
/***********************************************************/
/* Do not call this function directly! Use poly_getZeros() */
/* instead!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      */
/* get the zeros of a poly, g, that is known to split      */
/* completely over Z/pZ                                    */
/* degree of g must be >= 1                                */
/***********************************************************/
{ s32   b, w, w2, wl;
  poly_t temp, temp2, temp3, temp4, x_b;
  int    i, j;

  /*******************************************/
  /* termination condition for the recursion */
  /*******************************************/
  if (g->degree==1) {
    w = inverseModP(g->coef[1], p);
    w2 = mulmod32(w, g->coef[0], p);
    w2 = (p-w2)%p;  
    zeros[*numZeros] = w2;
    *numZeros += 1;
    return 1;
  }

CHOOSE_B: 
  b = prand()%p;
  x_b->coef[1]=1;
  x_b->coef[0]= (p-b)%p;
  x_b->degree = 1;
  for (i=MAXPOLYDEGREE; i>=0; i--)
    temp2->coef[i]=0;

  /**************************************/
  /*      Compute temp2 <-- g(x-b)      */
  /**************************************/
  temp->coef[0]=1;
  temp->degree=0;
  
  /* At the i-th iteration, temp = (x-b)^i */
  for (i=0; i<= g->degree; i++) {
    /**************************************************************/
    /*     temp2 <-- temp2 + g[i]*(x-b)^i = temp2 + g[i]*temp     */
    /**************************************************************/
    for (j=temp->degree; j>=0; j--) {
      wl = mulmod32(temp->coef[j], g->coef[i], p);
      temp2->coef[j] = (temp2->coef[j] + wl)%p;
    }
    poly_mul_modp(temp3, temp, x_b, p);
    poly_cp(temp, temp3);
  }
  temp2->degree = g->degree;
  
  /***********************************/
  /* If x | temp2, then -b is a root */
  /***********************************/
  if (temp2->coef[0]==0) {

    zeros[*numZeros] = (p-b)%p;
    *numZeros += 1;
    for (i=0; i<temp2->degree; i++)
      temp2->coef[i] = temp2->coef[i+1];
    if (--temp2->degree > 0) {
      i=*numZeros;
      get_zeros_rec(zeros, numZeros, temp2, p);
      for (j=i; j<*numZeros; j++)
        zeros[j] = (p+zeros[j]-b)%p;
    }
    return 1;
  }
  /********************************************/
  /* x does not divide temp2, so try to split */
  /* temp2 by computing:                      */
  /* gcd(temp2, x^{(p-1)/2} + 1) and          */
  /* gcd(temp2, x^{(p-1)/2} - 1) and          */
  /********************************************/
  poly_xpow_modpp(temp, ((p-1)>>1) , temp2, p);
  for (i=temp3->degree=temp->degree; i>=0; i--)
    temp3->coef[i] = temp->coef[i];
  temp->coef[0] = (temp->coef[0] + 1)%p;
  temp3->coef[0] =(temp3->coef[0]-1)%p;

  poly_cp(temp4, temp2);

  poly_gcd(temp, temp2, p);
  poly_gcd(temp3, temp4, p);

  if ((temp->degree==0)||(temp3->degree==0))  {
    goto CHOOSE_B;
  }

  /**********************************************/
  /* If we're here, we sucessfully split g(x-b) */
  /**********************************************/
  i = *numZeros;
  get_zeros_rec(zeros, numZeros, temp, p);
  for (j=i; j<*numZeros; j++)
    zeros[j] = (p+zeros[j]-b)%p;
  i=*numZeros;
  get_zeros_rec(zeros, numZeros, temp3, p);
  for (j=i; j<*numZeros; j++)
    zeros[j] = (p+zeros[j]-b)%p;
  return *numZeros;
}


/******************************************************/    
int poly_getZeros(s32 *zeros, poly_t _f, s32 p)
/******************************************************/    
/* It is assumed that f_big has degree > 1            */
/******************************************************/    
{ poly_t g, x_pow, f;
  int    numZeros;
  static int initialized=0;
  static mpz_t temp;
  
  if (!(initialized)) {
    mpz_init(temp);
    initialized=1;
  }
  poly_cp(f, _f);
    
  /*************************/
  /* Compute gcd(f, x^p-x) */
  /*************************/
  poly_xpow_modpp(x_pow, p, f, p);

  poly_cp(g, x_pow);
  if (g->degree<=0) { /* If it's 0 or a constant, that'll change in a sec. */
    g->degree = 1;
    g->coef[1] = p-1;
  } else {
    g->coef[1] = (p+g->coef[1] -1)%p;
  }
  /******************************************************************/
  /* It could happen that we just subtracted away the leading term, */
  /* or created a new leading term, so we must find the degree.     */
  /******************************************************************/
  poly_fixDeg(g);

  
  poly_gcd(g, f, p);
  /*************************************************************/
  /* Now g is the product of the linear factors of f_big mod p */
  /*************************************************************/

  numZeros=0;
  /**********************************/
  /* take care of the trivial cases */
  /**********************************/
  if (g->degree < 1) {  
    return 0;
  }
  if (g->degree == 1) {
    zeros[0] = mulmod32( (p-g->coef[0])%p, inverseModP(g->coef[1], p), p);
    numZeros=1;
    return 1;
  }
  
  get_zeros_rec(zeros, &numZeros, g, p);
  return numZeros;
}

/******************************************************/    
int poly_isSimpleZero(poly_t f, s32 p, s32 r)
/******************************************************/    
/* Assuming f(r)==0 (mod p), return f'(r) mod p.      */
/******************************************************/    
{ int    i;
  s32   y, xPow;
  poly_t dF;

  dF->degree=f->degree-1;
  for (i=0; i<f->degree; i++)
    dF->coef[i] = mulmod32((i+1), f->coef[i+1], p);
  poly_fixDeg(dF);
  xPow=1;
  y=0;
  for (i=0; i<dF->degree; i++) {
    y = (y + mulmod32(xPow, dF->coef[i], p))%p;
    xPow = mulmod32(xPow, r, p);
  }
  return y;
}

	

/**************************************************************************/
INLINE int poly_coprime(poly_t f, poly_t g, s32 p)
/**************************************************************************/
/* Return 1 if gcd(f,g)=1, where f and g are polys mod p.                 */
/**************************************************************************/
{ poly_t tmp1, tmp2;

  poly_cp(tmp1, f);
  poly_cp(tmp2, g);
 
  poly_gcd(tmp1, tmp2, p);
  if (tmp1->degree > 0)
    return 0;
  return 1;
}

/**************************************************************************/
int poly_irreducible_modp(poly_t _f, s32 p)
/**************************************************************************/
/* Return value: 1, if '_f' is irreducible mod p. 0 otherwise.            */
/**************************************************************************/
{ poly_t xPow, tmp, f;
  int    i;

  poly_cp(f, _f);
  
  /* First, compute x^p mod f. */
  poly_xpow_modpp(xPow, p, f, p);
  
  /* Verify that (x^p - x, f) = 1. */
  tmp->degree = (xPow->degree > 1) ? xPow->degree : 1;
  tmp->coef[1]=0;
  for (i=0; i<= xPow->degree; i++)
    tmp->coef[i] = xPow->coef[i];
  tmp->coef[1] -= 1;
  if (tmp->coef[1] < 0)
    tmp->coef[1] += p;
  if (!(poly_coprime(tmp, f, p)))
    return 0;

  /* Now, compute x^{p^{deg_f}} mod f. */
  for (i=2; i<=f->degree; i++)
    poly_pow_modpp(xPow, xPow, p, f, p);
  
  /* Check that x^{p^{deg_f}} == x (mod f). */
  if (xPow->degree != 1)
    return 0;
  if (xPow->coef[1] != 1)
    return 0;
  return 1;
}

