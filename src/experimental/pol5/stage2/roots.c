#include "stage2_impl.h"

#define DBL_MIN   -1e20
#define DBL_MAX    1e20

/*-----------------------------------------------------------------------*/
static double
poly_eval(int deg, double *poly, double x)
{
	double res, xn;
	int i;

	res = poly[0];
	xn = 1.;
	for (i = 1; i <= deg; i++) {
		xn *= x;
		res += poly[i] * xn;
	}
	return res;
}

/*-----------------------------------------------------------------------*/
static int
find_roots_core(int deg, double *poly, double d0, double d1, double *dr)
{
	double p0, p1, db, de, pm;
	volatile double dm;	/* important for some compilers */

	p0 = poly_eval(deg, poly, d0);
	p1 = poly_eval(deg, poly, d1);
	if (p1 == 0.) {
		*dr = d1;
		return 1;
	}
	if (p0 == 0.)
		return 0;
	if (p0 * p1 > 0.)
		return 1;
	db = d0;
	de = d1;
	while (de > db) {
		dm = (db + de) / 2.;
		if ((dm == db) || (dm == de))
			break;
		pm = poly_eval(deg, poly, dm);
		if (pm == 0.)
			break;
		if (p0 * pm > 0.) {
			db = dm;
			p0 = pm;
		}
		else {
			de = dm;
			p1 = pm;
		}
	}
	*dr = dm;
	return 1;
}

/*-----------------------------------------------------------------------*/
static int
find_roots(int deg, double *polcoeff, double *roots)
{
	int n, deriv, i, j, nroots;
	double d, dd;
	double pol_deriv[MAX_POLY_DEGREE][MAX_POLY_DEGREE + 1];
	double derivroots[MAX_POLY_DEGREE + 2];

	if (deg == 0) {
		return 0;
	}
	else if (deg == 1) {
		roots[0] = -polcoeff[0] / polcoeff[1];
		return 1;
	}

	for (i = 0; i <= deg; i++) {
		pol_deriv[0][i] = polcoeff[i];
	}
	for (i = 1; i < deg; i++) {
		for (j = 0; j < deg; j++)
			pol_deriv[i][j] = pol_deriv[i - 1][j + 1] * (j + 1);
		pol_deriv[i][j] = 0.;
	}

	derivroots[0] = DBL_MIN;
	derivroots[1] = DBL_MAX;
	n = 0;
	for (deriv = deg - 1; deriv >= 0; deriv--) {
		j = 1;
		d = derivroots[0];
		for (i = 1; i < n + 2; i++) {
			dd = derivroots[i];
			if (find_roots_core(deg - deriv, pol_deriv[deriv], 
						d, dd, derivroots + j))
				j++;
			d = dd;
		}
		if (derivroots[j - 1] < DBL_MAX)
			derivroots[j++] = DBL_MAX;
		n = j - 2;
	}
	for (i = 0; i < n; i++)
		roots[i] = derivroots[i + 1];
	nroots = n;

	i = 0;
	while (i < nroots - 1) {
		if (roots[i] == roots[i + 1]) {
			for (j = i + 1; j < nroots - 1; j++)
				roots[j] = roots[j + 1];
			nroots--;
			continue;
		}
		if (roots[i] > roots[i + 1]) {
			d = roots[i];
			roots[i] = roots[i + 1];
			roots[i + 1] = d;
			if (i > 0)
				i--;
			continue;
		}
		i++;
	}
	return nroots;
}

/*-----------------------------------------------------------------------*/
int
find_poly_optima(int *deg, double **coeff, double skewness, double *optima)
{
	int n, nop, s;

/* zeros */
	nop = find_roots(deg[0], coeff[0], optima);
	nop += find_roots(deg[1], coeff[1], optima + nop);
	for (n = 0; n < nop; n++)
		optima[n] /= skewness;

/* extrema */
	for (s = 0; s < 2; s++) {
		double deriv[MAX_POLY_DEGREE + 1];
		double c[MAX_POLY_DEGREE + 1];
		double sp;
		int d, i;

		d = deg[s];
		sp = 1.;
		for (i = 0; i <= d; i++) {
			c[i] = coeff[s][i] * sp;
			sp *= skewness;
		}
		for (i = 0; i < d; i++)
			deriv[i] = -c[i + 1] * (double) (i + 1);
		deriv[d] = 0.;
		for (i = 0; i < d; i++)
			deriv[i + 1] += c[i] * (double) (d - i);
		while (d) {
			if (deriv[d])
				break;
			else
				d--;
		}
		nop += find_roots(d, deriv, optima + nop);
	}
	optima[nop++] = DBL_MIN;
	optima[nop++] = DBL_MAX;
	return nop;
}
