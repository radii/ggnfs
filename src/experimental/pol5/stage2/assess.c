#include "stage2_impl.h"

#define NSMALLPRIMES      6542

extern const double rho_table[400];

/*----------------------------------------------------------------------*/
static void
pol_mul_lin(unsigned int p, int deg, unsigned int *pol, unsigned int *modpol,
	    unsigned int a)
{
	int i;
	unsigned int coeff; 
	unsigned int pol2[MAX_POLY_DEGREE + 2];

	pol2[0] = 0;
	for (i = 0; i <= deg; i++)
		pol2[i + 1] = pol[i];
	for (i = 0; i <= deg; i++)
		pol2[i] += ((pol[i] * a) % p);
	if (pol2[deg + 1]) {
		coeff = pol2[deg + 1];
		for (i = 0; i <= deg; i++)
			pol2[i] += ((coeff * modpol[i]) % p);
	}
	for (i = 0; i <= deg; i++)
		pol[i] = pol2[i] % p;
}

/*----------------------------------------------------------------------*/
static void
pol_sqr(unsigned int p, int deg, unsigned int *pol, unsigned int *modpol)
{
	int i, j;
	unsigned int coeff; 
	unsigned int prod[2 * (MAX_POLY_DEGREE + 1)];

	for (i = 0; i <= 2 * deg; i++)
		prod[i] = 0;
	for (i = 0; i <= deg; i++)
		prod[2 * i] = (pol[i] * pol[i]) % p;
	for (i = 0; i < deg; i++)
		for (j = i + 1; j <= deg; j++)
			prod[i + j] += (((pol[i] * pol[j]) % p) << 1);
	for (i = deg; i >= 0; i--) {
		if (prod[i + deg]) {
			coeff = prod[i + deg];
			for (j = 0; j <= deg; j++)
				prod[i + j] += ((coeff * modpol[j]) % p);
		}
	}
	for (i = 0; i <= deg; i++)
		pol[i] = prod[i] % p;
}

/*----------------------------------------------------------------------*/
static void
pol_exp(unsigned int p, int deg, unsigned int *pol, unsigned int *modpol,
	unsigned int a, unsigned int ex)
{
	int i, n;
	unsigned int ma, inv;
	unsigned int pol2[MAX_POLY_DEGREE + 1];

	if (modpol[deg] != p - 1) {
		if (modpol[deg] == 0)
			complain("pol_exp\n");
		inv = invert(p - modpol[deg], p);
		for (i = 0; i <= deg; i++) {
			modpol[i] *= inv;
			modpol[i] %= p;
		}
	}
	pol2[0] = 1;
	for (i = 1; i <= deg; i++)
		pol2[i] = 0;
	n = 31;
	ma = 1UL << 31;
	while (n > 0) {
		if (ex & ma)
			break;
		ma >>= 1;
		n--;
	}
	for (; n >= 0; n--) {
		pol_sqr(p, deg, pol2, modpol);
		if (ex & ma)
			pol_mul_lin(p, deg, pol2, modpol, a);
		ma >>= 1;
	}
	for (i = 0; i <= deg; i++)
		pol[i] = pol2[i];
}

/*----------------------------------------------------------------------*/
static int
pol_gcd(unsigned int p, int deg, unsigned int *res, unsigned int *pol1,
	unsigned int *pol)
{
	int i, d2, d3, diff;
	unsigned int h, leadinv; 
	unsigned int pol2[MAX_POLY_DEGREE + 1];
	unsigned int pol3[MAX_POLY_DEGREE + 1];

	for (i = 0; i <= deg; i++) {
		pol2[i] = pol[i];
		pol3[i] = pol1[i];
	}
	d2 = deg;
	d3 = deg;
	while (d2 >= 0)
		if (pol2[d2])
			break;
		else
			d2--;
	while (d3 >= 0)
		if (pol3[d3])
			break;
		else
			d3--;
	if (d2 == -1) {
		for (i = 0; i <= d3; i++)
			res[i] = pol3[i];
		return d3;
	}
	if (d3 == -1) {
		for (i = 0; i <= d2; i++)
			res[i] = pol2[i];
		return d2;
	}
	h = p - pol2[d2];
	leadinv = invert(h, p);
	for (i = 0; i <= d2; i++)
		pol2[i] = ((leadinv * pol2[i]) % p);
	h = p - pol3[d3];
	leadinv = invert(h, p);
	for (i = 0; i <= d3; i++)
		pol3[i] = ((leadinv * pol3[i]) % p);
	while (1) {
		if (d3 == -1) {
			for (i = 0; i <= d2; i++)
				res[i] = pol2[i];
			return d2;
		}
		if (d2 == -1) {
			for (i = 0; i <= d3; i++)
				res[i] = pol3[i];
			return d3;
		}
		if (d3 >= d2) {
			diff = d3 - d2;
			for (i = 0; i <= d2; i++) {
				pol3[i + diff] += p;
				pol3[i + diff] -= pol2[i];
				pol3[i + diff] %= p;
			}
			while (d3 >= 0)
				if (!pol3[d3])
					d3--;
				else
					break;
			if (d3 == -1)
				continue;
			if (pol3[d3] != p - 1) {
				h = p - pol3[d3];
				leadinv = invert(h, p);
				for (i = 0; i <= d3; i++)
					pol3[i] = ((leadinv * pol3[i]) % p);
			}
		}
		else {		/* d3<d2 */
			diff = d2 - d3;
			for (i = 0; i <= d3; i++) {
				pol2[i + diff] += p;
				pol2[i + diff] -= pol3[i];
				pol2[i + diff] %= p;
			}
			while (d2 >= 0)
				if (!pol2[d2])
					d2--;
				else
					break;
			if (d2 == -1)
				continue;
			if (pol2[d2] != p - 1) {
				h = p - pol2[d2];
				leadinv = invert(h, p);
				for (i = 0; i <= d2; i++)
					pol2[i] = ((leadinv * pol2[i]) % p);
			}
		}
	}
	Schlendrian("pol_gcd\n");
	return 0;
}

/*----------------------------------------------------------------------*/
static void
pol_div_lin(unsigned int p, int *deg, unsigned int *coeff, unsigned int r)
{
	int i, d;
	unsigned int h, res;

	d = *deg;
	h = coeff[d];
	for (i = d - 1; i >= 0; i--) {
		res = h;
		h *= r;
		h += coeff[i];
		h %= p;
		coeff[i] = res;
	}
	if (h)
		complain("pol_div: not divisible\n");
	*deg = d - 1;
}

/*----------------------------------------------------------------------*/
static unsigned int
find_root(unsigned int p, int deg, unsigned int *coeff)
{
	unsigned int rpol[MAX_POLY_DEGREE + 1];
	unsigned int rpol1[MAX_POLY_DEGREE + 1]; 
	unsigned int rpol2[MAX_POLY_DEGREE + 1]; 
	unsigned int res, inv, a, p2;
	int i, rd, d;

	for (i = 0; i <= deg; i++)
		rpol[i] = coeff[i];
	if (p == 2) {
		if (!rpol[0])
			return 0;
		a = 0;
		for (i = 0; i <= deg; i++)
			a += rpol[i];
		if ((a & 1) == 0)
			return 1;
		return 2;	/* no root */
	}
	rd = deg;
	p2 = (p - 1) / 2;
	if (rd > 1) {
		pol_exp(p, rd, rpol1, rpol, 0, p);
		if (rpol1[1])
			rpol1[1]--;
		else
			rpol1[1] = p - 1;
		d = 0;
		for (i = 0; i <= rd; i++) {
			if (rpol1[i])
				d = 1;
		}
		if (d) {
			d = pol_gcd(p, rd, rpol, rpol, rpol1);
			rd = d;
		}		/* else rpol has rd simple roots */
		if (rd == 0)
			return p;	/* no root */
	}
	while (rd > 1) {
		a = ((unsigned int) rand()) % p;
		pol_exp(p, rd, rpol1, rpol, a, p2);
		if (rpol1[0])
			rpol1[0]--;
		else
			rpol1[0] = p - 1;
		d = pol_gcd(p, rd, rpol2, rpol1, rpol);
		if (d) {
			rd = d;
			for (i = 0; i <= rd; i++) {
				rpol[i] = rpol2[i];
			}
		}
	}
	if (rd == 0)
		complain("find_root: degree=0?\n");
	inv = invert(rpol[1], p);
	res = inv * rpol[0];
	res %= p;
	if (res)
		res = p - res;
	return res;
}

/*----------------------------------------------------------------------*/
static unsigned int
pol_eval(unsigned int p, int deg, unsigned int *coeff, unsigned int r)
{
	int i;
	unsigned int res;

	res = 0;
	for (i = deg; i >= 0; i--) {
		res *= r;
		res += coeff[i];
		res %= p;
	}
	return res;
}

/*----------------------------------------------------------------------*/
static int
find_pol_roots_homog(unsigned int p, int deg, unsigned int *polmod,
		     unsigned int *roots)
{
	unsigned int r;
	int n, d, i;
	unsigned int mod[MAX_POLY_DEGREE + 1];

	for (i = 0; i <= deg; i++)
		mod[i] = polmod[i] % p;
	d = deg;
	n = 0;
	while ((d) && (mod[d] == 0)) {
		d--;
		roots[n++] = p;
	}
	if (d == 0) {
		if (mod[0] == 0)
			complain("find_pol_roots_homog: zero pol\n");
		return n;
	}
	while (d) {
		r = find_root(p, d, mod);
		if (r == p)
			break;
		if (pol_eval(p, d, mod, r))
			complain("find_pol_roots_homog\n");
		do {
			pol_div_lin(p, &d, mod, r);
			roots[n++] = r;
		} while (pol_eval(p, d, mod, r) == 0);
	}
	return n;
}

/*----------------------------------------------------------------------*/
#define MAX_BF      65535

static double
brute_force(unsigned int p, int deg, mpz_t * gmp_coeff, unsigned int r)
{
	unsigned int qmax, h, q, q0;
	unsigned int lifts[32], powers[32];
	int ind, indmax, j;
	double dp, lo, res;
	unsigned int coeffmod[MAX_POLY_DEGREE + 1];

	qmax = 1;
	indmax = 0;
	while (qmax < MAX_BF / p) {
		qmax *= p;
		powers[indmax++] = qmax;
	}
	if (r == p) {	/* projective root; switch to poly(1/x) */
		r = 0;
		for (j = 0; j <= deg; j++)
			coeffmod[j] = mpz_fdiv_ui(gmp_coeff[deg - j], qmax);
	}
	else {
		for (j = 0; j <= deg; j++)
			coeffmod[j] = mpz_fdiv_ui(gmp_coeff[j], qmax);
	}
	dp = (double) p;
	lo = log(dp) * dp / (dp + 1);
	ind = 0;
	lifts[0] = r;
	res = 0.;
	while (1) {
		q = powers[ind];
		res += lo / ((double) q);
		if (ind < indmax - 1) {	/* lift root */
			h = lifts[ind];
			while (h < q * p) {
				if (pol_eval(q * p, deg, coeffmod, h) == 0)
					break;
				h += q;
			}
			if (h < q * p) {
				ind++;
				lifts[ind] = h;
				continue;
			}
		}		/* find next lift or go back */
		while (ind >= 0) {
			q = powers[ind];
			q0 = q / p;
			h = lifts[ind] + q0;
			while (h < q) {
				if (pol_eval(q, deg, coeffmod, h) == 0)
					break;
				h += q0;
			}
			lifts[ind] = h;
			if (h < q)
				break;	/* found another zero */
			ind--;
		}
		if (ind > 0)
			continue;
		break;		/* finished */
	}
	return res;
}

/*----------------------------------------------------------------------*/
static void
compute_coeffmod(unsigned int p, int deg, 
		mpz_t *coeff, unsigned int *coeffmod)
{
	int i;

	for (i = 0; i <= deg; i++)
		coeffmod[i] = mpz_fdiv_ui(coeff[i], p);
}

/*----------------------------------------------------------------------*/
/* returns 1 if alpha>alpha_targ */
int
murphy_alpha(double *alpha, int deg, assess_t *assess,
		mpz_t * gmp_coeff, double alpha_targ)
{
	double al, dp;
	unsigned int p, r;
	int i, j, n;
	double al_prev, max_rest, rat_rest;
	unsigned int roots[MAX_POLY_DEGREE + 1];
	unsigned int coeffmod[MAX_POLY_DEGREE + 1];

	al = 0.;
	max_rest = assess->alpha_max;
	rat_rest = assess->alpha_random;
	for (i = 0; i < NSMALLPRIMES; i++) {
		al_prev = al;
		p = assess->prime_list[i];
		if (p > assess->prime_bound)
			break;
		compute_coeffmod(p, deg, gmp_coeff, coeffmod);
		n = find_pol_roots_homog(p, deg, coeffmod, roots);
		j = 0;
		dp = (double) p;
		while (j < n) {
			r = roots[j];
			if ((j + 1 < n) && (r == roots[j + 1])) {
				/* multiple zero */
				al -= brute_force(p, deg, gmp_coeff, r);
			}
			else {
				al -= dp * log(dp) / (dp * dp - 1);
			}
			j++;
			while ((j < n) && (roots[j] == r))
				j++;
		}
		al += log(dp) / (dp - 1);
		if (verbose > 3) {
			printf("%u   %.3f %.3f:  ", p, al, al - al_prev);
			for (j = 0; j < n; j++)
				printf("%u ", roots[j]);
			printf("\n");
		}
		if (p > 100) {
			max_rest -= dp * log(dp) / (dp * dp - 1);
			rat_rest -= log(dp) / (dp - 1);
			if (al + rat_rest - 
				max_rest * ((double) deg) > alpha_targ) {
				return 1;
			}
		}
	}
	if (verbose > 3)
		printf("\n");
	*alpha = al;
	return 0;
}

/*----------------------------------------------------------------------*/
void
murphy_alpha_exact(double *alpha, int deg, assess_t *assess,
		    mpz_t * gmp_coeff, unsigned int pb)
{
	double al, dp;
	unsigned int p, r;
	int i, j, n;
	double al_prev;
	unsigned int coeffmod[MAX_POLY_DEGREE + 1];
	unsigned int roots[MAX_POLY_DEGREE + 1];

	al = 0.;
	for (i = 0; i < NSMALLPRIMES; i++) {
		al_prev = al;
		p = assess->prime_list[i];
		if (p > pb)
			break;
		compute_coeffmod(p, deg, gmp_coeff, coeffmod);
		n = find_pol_roots_homog(p, deg, coeffmod, roots);
		j = 0;
		dp = (double) p;
		while (j < n) {
			r = roots[j];
			if ((j + 1 < n) && (r == roots[j + 1])) {
				/* multiple zero */
				al -= brute_force(p, deg, gmp_coeff, r);
			}
			else {
				al -= dp * log(dp) / (dp * dp - 1);
			}
			j++;
			while ((j < n) && (roots[j] == r))
				j++;
		}
		al += log(dp) / (dp - 1);
		if (verbose > 3) {
			printf("%u   %.3f %.3f:  ", p, al, al - al_prev);
			for (j = 0; j < n; j++)
				printf("%u ", roots[j]);
			printf("\n");
		}
	}
	if (verbose > 3)
		printf("\n");
	*alpha = al;
}

/*----------------------------------------------------------------------*/
static int
compare_doubles(const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;

	return (*da > *db) - (*da < *db);
}

/*----------------------------------------------------------------------*/
#define EULER_C   0.577215664901532

static double
dickman(double x)
{
	int i;
	double h;

	if (x <= 1.)
		return 1.;
	if (x >= 20.)
		return 0.;
	i = (int) (20 * x);
	h = 20 * x - (double) i;
	return (h * rho_table[i] + (1 - h) * rho_table[i - 1]);
}

static double
prob(double r, double b)
{
	double u, logr;

	logr = log(fabs(r));
	u = logr / log(b);
	return dickman(u);
}

/*----------------------------------------------------------------------*/
void
murphy_e_score(double *me, int deg0, double *dbl_coeff0, 
		int deg1, double *dbl_coeff1, 
		double alpha0, double alpha1, 
		double skewness, int nsm, assess_t *assess)
{
	int i, j;
	double x, y, sx, sy, theta;
	double v0, v1, xp, e, al0, al1;
	double e0, left, right, theta_left, theta_right, theta_len, theta_inc;
	int interval, nop, nsum;
	int deg[2];
	double *dbl_coeff[2];
	double optima[4 * MAX_POLY_DEGREE + 2];

	sx = sqrt(assess->area * skewness);
	sy = sx / skewness;
	al0 = exp(alpha0);
	al1 = exp(alpha1);
	deg[0] = deg0;
	deg[1] = deg1;
	dbl_coeff[0] = dbl_coeff0;
	dbl_coeff[1] = dbl_coeff1;

	/* find the combined list of roots and turning points
	   of both polynomials; these will be turned into 
	   intervals that separately contribute to the score.
	   The intervals will also include a smallest and largest
	   point */

	nop = find_poly_optima(deg, dbl_coeff, skewness, optima);
	qsort(optima, nop, sizeof(double), compare_doubles);

	e = 0;
	for (interval = 1; interval < nop; interval++) {
		left = optima[interval - 1];
		right = optima[interval];
		theta_left = atan(left);
		theta_right = atan(right);
		theta_len = theta_right - theta_left;
		nsum = (int) (((double) nsm) / M_PI * theta_len);
		if (nsum < 10)
			nsum = 10;
		theta_inc = theta_len / (double) (nsum);
		e0 = 0;
		theta = theta_left + theta_inc / 2.;
		for (i = 0; i < nsum; i++) {
			theta += theta_inc;
			x = sx * sin(theta);
			y = sy * cos(theta);
			v0 = dbl_coeff0[0];
			xp = 1.;
			for (j = 1; j <= deg0; j++) {
				xp *= x;
				v0 *= y;
				v0 += xp * dbl_coeff0[j];
			}
			v1 = dbl_coeff1[0];
			xp = 1.;
			for (j = 1; j <= deg1; j++) {
				xp *= x;
				v1 *= y;
				v1 += xp * dbl_coeff1[j];
			}
			e0 += prob(al0 * v0, assess->bound0) * 
				prob(al1 * v1, assess->bound1);
		}
		e0 /= nsum;
		e += (e0 * theta_len);
	}
	e /= M_PI;
	*me = e;
}

/*----------------------------------------------------------------------*/
void
assess_init(assess_t *a, double b0, double b1, 
		double area, unsigned int pb)
{
	int i;
	double dp;

	a->bound0 = b0;
	a->bound1 = b1;
	a->area = area;
	a->prime_bound = pb;
	prime_table_init();
	a->prime_list = (unsigned int *)xmalloc(NSMALLPRIMES *
						sizeof(unsigned int));
	a->alpha_max = 0.;
	a->alpha_random = 0.;
	for (i = 0; i < NSMALLPRIMES; i++) {
		unsigned int p = get_next_prime();
		a->prime_list[i] = p;
		if (p < 100)
			continue;
		if (p > a->prime_bound)
			break;
		dp = (double)p;
		a->alpha_max += dp * log(dp) / (dp * dp - 1);
		a->alpha_random += log(dp) / (dp - 1);
	}
}

/*----------------------------------------------------------------------*/
void
assess_free(assess_t *a)
{
	free(a->prime_list);
}
