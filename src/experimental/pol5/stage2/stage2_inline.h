#ifndef _STAGE1_INLINE_H_
#define _STAGE1_INLINE_H_

#include "ggnfs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------*/
static INLINE unsigned int
invert(unsigned int a, unsigned int p)
{				/* 0<b<p */
	unsigned int v1 = 0, v2 = 1, q, b = a, oldp = p;

	if (a == 0)
		complain("cannot invert 0\n");
	while (b > 1) {
		p -= b;
		v1 += v2;
		if (p >= b) {
			p -= b;
			v1 += v2;
			if (p >= b) {
				q = p / b;
				v1 += q * v2;
				q *= b;
				p -= q;
			}
		}
		if (p <= 1) {
			v2 = oldp - v1;
			break;
		}
		b -= p;
		v2 += v1;
		if (b >= p) {
			b -= p;
			v2 += v1;
			if (b >= p) {
				q = b / p;
				v2 += q * v1;
				q *= p;
				b -= q;
			}
		}
	}
	if ((v2 * a - 1) % oldp)
		complain("invert %u %u %u \n", v2, a, oldp);
	return v2;
}

/*-------------------------------------------------------------------------*/
static INLINE double
ifs(double *coeff, double skewness)
{				/* degree=5 */
	double sq[6], d, s, res;
	int i, j, k;

	for (i = 0; i < 6; i++)
		sq[i] = coeff[i] * coeff[i];
	for (i = 0; i < 4; i++) {
		d = 2 * coeff[i];
		k = i + 1;
		for (j = i + 2; j < 6; j += 2, k++)
			sq[k] += d * coeff[j];
	}
	s = skewness;
	res = 0.;
	res += s * sq[3] / 35.;
	res += sq[2] / s / 35.;
	s *= (skewness * skewness);
	res += s * sq[4] / 27.;
	res += sq[1] / s / 27.;
	s *= (skewness * skewness);
	res += s * sq[5] / 11.;
	res += sq[0] / s / 11.;
	return res;
}


/*-------------------------------------------------------------------------*/
#define COMPUTE_IFS   s2=sc; res=0.; \
  res+=s2*sq[3]/35.; res+=sq[2]/s2/35.;   \
  s2*=(sc*sc);   \
  res+=s2*sq[4]/27.; res+=sq[1]/s2/27.;   \
  s2*=(sc*sc);   \
  res+=s2*sq[5]/11.; res+=sq[0]/s2/11.


static INLINE double
find_best_skewness(double *coeff, double s0)
{				/* degree=5 */
	double sq[6], d, s, sc, s2, res, ds, v;
	int i, j, k;

	for (i = 0; i < 6; i++)
		sq[i] = coeff[i] * coeff[i];
	for (i = 0; i < 4; i++) {
		d = 2 * coeff[i];
		k = i + 1;
		for (j = i + 2; j < 6; j += 2, k++)
			sq[k] += d * coeff[j];
	}
	ds = 10.;
	s = s0;
	sc = s;
	COMPUTE_IFS;
	v = res;
	while (ds > 2.) {
		sc = s + ds;
		COMPUTE_IFS;
		if (res < v) {
			s = sc;
			v = res;
			ds *= 1.1;
			continue;
		}
		if (ds < s) {
			sc = s - ds;
			COMPUTE_IFS;
			if (res < v) {
				s = sc;
				v = res;
				ds *= 1.1;
				continue;
			}
		}
		ds *= 0.9;
	}
	return s;
}

#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_INLINE_H_ */
