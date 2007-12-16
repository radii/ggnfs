/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	
       				   --jasonp@boo.net 12/13/07
--------------------------------------------------------------------*/

/* Fast Hartley Transform based convolution routines

   See K. Jayasuriya, "Multiprecision Arithmetic Using
   Fast Hartley Transforms", Appl. Math. and Comp. 75:239-251 (1996)

   This paper explains more about the implementation aspects
   of FHTs than any other resource I've seen
*/

#include <fastmult.h>

#define FFT_BASE 65536
#define FFT_WORDS_PER_LIMB 2

#ifndef M_PI
#define M_PI  3.1415926535897932384626433832795029
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.7071067811865475244008443621048490
#endif

/*------------------------------------------------------------------------*/
#ifndef NDEBUG
	/* used to determine the floating point precision
	   level, intended for x86 systems */

#if defined(WIN32) 

#include <float.h>
typedef uint32 precision_t;
static uint32 precision_is_ieee(void) {
	precision_t prec = _control87(0, 0);
	return  ((prec & _MCW_PC) == _PC_53) ? 1 : 0;
}

#elif defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))

#include <float.h>
typedef uint16 precision_t;
static uint32 precision_is_ieee(void) {
	precision_t prec;
	asm volatile ("fnstcw %0":"=m"(prec));
	return ((prec & ~0x0300) == 0x0200) ? 1 : 0;
}

#else

static uint32 precision_is_ieee(void) {
	return 1;
}

#endif
#endif  /* !NDEBUG */

/*------------------------------------------------------------------------*/
double ieee_round = 6755399441055744.0;
double intel_round = 6755399441055744.0 * 2048.0;
double round_test = 2.7;
double round_correct = 3.0;
 
static void init_fast_rounding(fastmult_info_t *info) {

	/* the two paths below are essentially redundant unless
	   this is an x86 processor */

#ifdef NDEBUG
	/* as long as this isn't a debug build, we can
	   determine the FPU precision experimentally */

	double round_me;
 
	info->round_constant[0] = intel_round;
	info->round_constant[1] = intel_round;
	round_me = round_test;
	round_me += info->round_constant[0];
	round_me -= info->round_constant[1];
	if (round_me == round_correct)
		return;

	info->round_constant[0] = ieee_round;
	info->round_constant[1] = ieee_round;
	round_me = round_test;
	round_me += info->round_constant[0];
	round_me -= info->round_constant[1];
	if (round_me != round_correct) {
		printf("error: can't initialize fast rounding\n");
		exit(-1);
	}

#else
	/* for debug builds, quantities are always spilled to
	   registers so the precision would always look like
	   IEEE. Instead, detect the FPU precision directly.
	   We prefer the code path above because the compiler
	   may be using SSE2 registers instead of the FPU stack,
	   and this code would guess wrong in that case */

	if (precision_is_ieee()) {
		info->round_constant[0] = 6755399441055744.0;
		info->round_constant[1] = 6755399441055744.0;
	}
	else {
		info->round_constant[0] = 6755399441055744.0 * 2048.0;
		info->round_constant[1] = 6755399441055744.0 * 2048.0;
	}
#endif

}

/*------------------------------------------------------------------------*/
void fastmult_info_init(fastmult_info_t *info) {

	memset(info, 0, sizeof(fastmult_info_t));
	init_fast_rounding(info);
}

/*------------------------------------------------------------------------*/
void fastmult_info_free(fastmult_info_t *info) {
	
	uint32 i;
	uint32 num_twiddle = MIN(info->log2_runlength, HUGE_TWIDDLE_CUTOFF);
	uint32 num_huge_twiddle = info->log2_runlength - num_twiddle;

	for (i = 0; i <= num_twiddle; i++) {
		free(info->twiddle[i]);
	}
	for (i = 0; i <= num_huge_twiddle; i++) {
		free(info->huge_twiddle[i].small);
		free(info->huge_twiddle[i].large);
	}
	memset(info, 0, sizeof(fastmult_info_t));
}

/*------------------------------------------------------------------------*/
static double *make_twiddle(int32 size) {

	int32 i;
	double arg = (2.0 * M_PI) / size;
	double *w = (double *)xmalloc(2 * (size / 4 - 1) * sizeof(double));

	for (i = 0; i < size / 4 - 1; i++) {
		w[2 * i] = cos((i+1) * arg);
		w[2 * i + 1] = sin((i+1) * arg);
	}

	return w;
}

/*------------------------------------------------------------------------*/
static void make_huge_twiddle(huge_twiddle_t *table, int32 entry) {

	int32 i;
	int32 power = HUGE_TWIDDLE_CUTOFF + entry;
	double arg = (2.0 * M_PI) / (1 << power);
	int32 small_power = (power - 2) / 2;
	int32 large_power = power - 2 - small_power;

	table->large = (double *)xmalloc((2 << large_power) * sizeof(double));
	table->small = (double *)xmalloc((2 << small_power) * sizeof(double));

	for (i = 0; i < (1 << small_power); i++) {
		table->small[2 * i] = cos(i * arg);
		table->small[2 * i + 1] = sin(i * arg);
	}

	for (i = 0; i < (1 << large_power); i++) {
		table->large[2 * i] = cos((i << small_power) * arg);
		table->large[2 * i + 1] = sin((i << small_power) * arg);
	}
}

/*------------------------------------------------------------------------*/
static void create_twiddle(fastmult_info_t *info, int32 power) {

	int32 i;

	if (power <= (int32)(info->log2_runlength))
		return;

	if (power > MAX_FHT_POWER) {
		printf("error: transform size 1 << %d too large\n", power);
		exit(-1);
	}
	for (i = 4; i <= power; i++) {
		if (i < HUGE_TWIDDLE_CUTOFF) {
			if (info->twiddle[i] == NULL)
				info->twiddle[i] = make_twiddle(1 << i);
		}
		else {
			int32 j = i - HUGE_TWIDDLE_CUTOFF;
			huge_twiddle_t *h = info->huge_twiddle + j;
			if (h->large == NULL || h->small == NULL) {
				make_huge_twiddle(h, j);
			}
		}
	}
	info->log2_runlength = power;
}

/*------------------------------------------------------------------------*/
static void fht8_dif(double *x) {

	double t0, t1, t2, t3, t4, t5, t6, t7;
	double t8, t9, t10, t11, t12, t13, t14, t15;

	t0 = x[0];
	t2 = x[2];
	t4 = x[4];
	t6 = x[6];

	t8 = t0 + t4;
	t12 = t0 - t4;
	t10 = t6 + t2;
	t14 = t6 - t2;

	t1 = x[1];
	t3 = x[3];
	t5 = x[5];
	t7 = x[7];

	t9 = t1 + t5;
	t13 = t1 - t5;
	t11 = t7 + t3;
	t15 = t7 - t3;
	t6 = M_SQRT1_2 * (t15 + t13);
	t7 = M_SQRT1_2 * (t15 - t13);

	t0 = t8 + t10;
	t2 = t8 - t10;
	t1 = t11 + t9;
	t3 = t11 - t9;

	x[0] = t0 + t1;
	x[1] = t0 - t1;
	x[2] = t2 + t3;
	x[3] = t2 - t3;

	t0 = t12 + t14;
	t2 = t12 - t14;
	t1 = t7 + t6;
	t3 = t7 - t6;

	x[4] = t0 + t1;
	x[5] = t0 - t1;
	x[6] = t2 + t3;
	x[7] = t2 - t3;
}

/*------------------------------------------------------------------------*/
static void fht_butterfly_dif(double *x, double *w, int32 size) {

	double t0, t1, t2, t3, t4, t5;
	double c, s;
	double *x0, *x1;
	int32 i = size / 4;

	x0 = x + i;
	x1 = x + 3 * i;

	t0 = x0[-i];
	t1 = x0[0];
	t2 = x1[-i];
	t3 = x1[0];
	x0[-i] = t0 + t2;
	x0[0]  = t3 + t1;
	x1[-i] = t0 - t2;
	x1[0]  = t3 - t1;
	i--;

	do {
		t0 = x0[-i];
		t1 = x0[i];
		t2 = x1[-i];
		t3 = x1[i];

		c = w[0];
		s = w[1];

		x0[-i] = t0 + t2;
		x0[i]  = t1 + t3;
		t4 = t0 - t2;
		t5 = t3 - t1;
		x1[-i] = c * t4 + s * t5;
		x1[i] = c * t5 - s * t4;

		w += 2;
	} while(--i);
}

/*------------------------------------------------------------------------*/
static void fht_butterfly_dif_huge(double *x, 
			huge_twiddle_t *h, int32 power) {

	int32 size = 1 << power;
	double t0, t1, t2, t3, t4, t5;
	double c, cs, cl, s, ss, sl;
	double *x0, *x1;
	int32 i = size / 4;
	int32 twiddle_arg;
	int32 shift = (power - 2) / 2;
	int32 mask = (1 << shift) - 1;

	x0 = x + i;
	x1 = x + 3 * i;

	t0 = x0[-i];
	t1 = x0[0];
	t2 = x1[-i];
	t3 = x1[0];
	x0[-i] = t0 + t2;
	x0[0]  = t3 + t1;
	x1[-i] = t0 - t2;
	x1[0]  = t3 - t1;
	i--;
	twiddle_arg = 1;

	do {
		t0 = x0[-i];
		t1 = x0[i];
		t2 = x1[-i];
		t3 = x1[i];

		cl = h->large[2*(twiddle_arg >> shift)];
		sl = h->large[2*(twiddle_arg >> shift) + 1];
		cs = h->small[2*(twiddle_arg & mask)];
		ss = h->small[2*(twiddle_arg & mask) + 1];
		c = cl * cs - sl * ss;
		s = cl * ss + cs * sl;

		x0[-i] = t0 + t2;
		x0[i]  = t1 + t3;
		t4 = t0 - t2;
		t5 = t3 - t1;
		x1[-i] = c * t4 + s * t5;
		x1[i] = c * t5 - s * t4;

		twiddle_arg++;
	} while(--i);
}

/*------------------------------------------------------------------------*/
static void fht_forward(double *x, fastmult_info_t *info, int32 power) {

	int32 size;

	if (power == 3) {
		fht8_dif(x);
		return;
	}

	size = 1 << power;
	if (power < HUGE_TWIDDLE_CUTOFF) {
		fht_butterfly_dif(x, info->twiddle[power], size);
	}
	else {
		fht_butterfly_dif_huge(x, 
			info->huge_twiddle + (power - HUGE_TWIDDLE_CUTOFF), 
			power);
	}

	fht_forward(x, info, power - 1);
	fht_forward(x + size / 2, info, power - 1);
}

/*------------------------------------------------------------------------*/
static void fht8_dit(double *x) {

	double t0, t1, t2, t3, t4, t5, t6, t7;
	double t8, t9, t10, t11;

	t0 = x[0];
	t1 = x[1];
	t2 = x[2];
	t3 = x[3];

	t4 = t0 + t1;
	t5 = t0 - t1;
	t6 = t2 + t3;
	t7 = t2 - t3;

	t0 = t4 + t6;
	t2 = t4 - t6;
	t1 = t5 - t7;
	t3 = t5 + t7;

	t4 = x[4];
	t5 = x[5];
	t6 = x[6];
	t7 = x[7];

	t8 = t4 + t5;
	t9 = t4 - t5;
	t10 = t6 + t7;
	t11 = t6 - t7;

	t4 = t8 + t10;
	t6 = t8 - t10;
	t5 = t9 - t11;
	t7 = t9 + t11;

	x[0] = t0 + t4;
	x[4] = t0 - t4;
	x[2] = t2 - t6;
	x[6] = t2 + t6;

	t0 = M_SQRT1_2 * (t5 - t7);
	t2 = M_SQRT1_2 * (t5 + t7);

	x[1] = t1 + t0;
	x[5] = t1 - t0;
	x[3] = t3 - t2;
	x[7] = t3 + t2;
}

/*------------------------------------------------------------------------*/
static void fht_butterfly_dit(double *x, double *w, int32 size) {

	double t0, t1, t2, t3, t4, t5;
	double c, s;
	double *x0, *x1;
	int32 i = size / 4;

	x0 = x + i;
	x1 = x + 3 * i;

	t0 = x0[-i];
	t1 = x0[0];
	t2 = x1[-i];
	t3 = x1[0];
	x0[-i] = t0 + t2;
	x0[0]  = t1 - t3;
	x1[-i] = t0 - t2;
	x1[0]  = t1 + t3;
	i--;

	do {
		t0 = x0[-i];
		t1 = x0[i];
		t2 = x1[-i];
		t3 = x1[i];

		c = w[0];
		s = w[1];

		t4 = c * t2 - s * t3;
		t5 = s * t2 + c * t3;

		x0[-i] = t0 + t4;
		x0[i]  = t1 - t5;
		x1[-i] = t0 - t4;
		x1[i] = t1 + t5;

		w += 2;
	} while(--i);
}

/*------------------------------------------------------------------------*/
static void fht_butterfly_dit_huge(double *x, 
			huge_twiddle_t *h, int32 power) {

	int32 size = 1 << power;
	double t0, t1, t2, t3, t4, t5;
	double c, cs, cl, s, ss, sl;
	double *x0, *x1;
	int32 i = size / 4;
	int32 twiddle_arg;
	int32 shift = (power - 2) / 2;
	int32 mask = (1 << shift) - 1;

	x0 = x + i;
	x1 = x + 3 * i;

	t0 = x0[-i];
	t1 = x0[0];
	t2 = x1[-i];
	t3 = x1[0];
	x0[-i] = t0 + t2;
	x0[0]  = t1 - t3;
	x1[-i] = t0 - t2;
	x1[0]  = t1 + t3;
	i--;
	twiddle_arg = 1;

	do {
		t0 = x0[-i];
		t1 = x0[i];
		t2 = x1[-i];
		t3 = x1[i];

		cl = h->large[2*(twiddle_arg >> shift)];
		sl = h->large[2*(twiddle_arg >> shift) + 1];
		cs = h->small[2*(twiddle_arg & mask)];
		ss = h->small[2*(twiddle_arg & mask) + 1];
		c = cl * cs - sl * ss;
		s = cl * ss + cs * sl;

		t4 = c * t2 - s * t3;
		t5 = s * t2 + c * t3;

		x0[-i] = t0 + t4;
		x0[i]  = t1 - t5;
		x1[-i] = t0 - t4;
		x1[i] = t1 + t5;

		twiddle_arg++;
	} while(--i);
}

/*------------------------------------------------------------------------*/
static void fht_inverse(double *x, fastmult_info_t *info, int32 power) {

	int32 size;

	if (power == 3) {
		fht8_dit(x);
		return;
	}

	size = 1 << power;
	fht_inverse(x, info, power - 1);
	fht_inverse(x + size / 2, info, power - 1);

	if (power < HUGE_TWIDDLE_CUTOFF) {
		fht_butterfly_dit(x, info->twiddle[power], size);
	}
	else {
		fht_butterfly_dit_huge(x, 
			info->huge_twiddle + (power - HUGE_TWIDDLE_CUTOFF),
			power);
	}
}

/*------------------------------------------------------------------------*/
static void square_hermitian(double *x, int32 size) {

	int32 p;
	double r = 1.0 / (2 * size);

	x[0] = x[0] * x[0] * 2 * r;
	x[1] = x[1] * x[1] * 2 * r;

	for (p = 1; p < size/2; p = p * 2) {

		int32 i;
		int32 lo = 2 * p;
		int32 hi = 4 * p - 1;

		for (i = 0; i < p; i++, lo++, hi--) {
			double d1 = 2 * x[lo] * x[hi];
			double d2 = x[lo] * x[lo] - x[hi] * x[hi];
			x[lo] = (d1 + d2) * r;
			x[hi] = (d1 - d2) * r;
		}
	}
}

/*------------------------------------------------------------------------*/
static void mul_hermitian(double *x, double *y, int32 size) {

	int32 p;
	double r = 1.0 / (2 * size);

	x[0] = x[0] * y[0] * 2 * r;
	x[1] = x[1] * y[1] * 2 * r;

	for (p = 1; p < size/2; p = p * 2) {

		int32 i;
		int32 lo = 2 * p;
		int32 hi = 4 * p - 1;

		for (i = 0; i < p; i++, lo++, hi--) {
			double d1 = x[lo] * y[hi] + x[hi] * y[lo];
			double d2 = x[lo] * y[lo] - x[hi] * y[hi];
			x[lo] = (d1 + d2) * r;
			x[hi] = (d1 - d2) * r;
		}
	}
}

/*------------------------------------------------------------------------*/
static int32 get_log2_fftsize(int32 nwords) {

	int32 i = 1;

	while ((1 << i) < nwords)
		i++;

	return i;
}

/*------------------------------------------------------------------------*/
static void packed_to_double(uint32 *p, int32 pwords,
			 double *d, int32 fftsize) {
	int32 i;
	int32 word;
	int32 borrow;

	for (i = borrow = 0; i < pwords; i++) {
		word = (int32)(p[i] % FFT_BASE) + borrow;
		borrow = 0;
		if (word >= FFT_BASE/2) {
			word -= FFT_BASE;
			borrow = 1;
		}
		d[2 * i] = (double)word;

		word = (int32)(p[i] / FFT_BASE) + borrow;
		borrow = 0;
		if (word >= FFT_BASE/2) {
			word -= FFT_BASE;
			borrow = 1;
		}
		d[2 * i + 1] = (double)word;
	}

	i = 2 * i;
	d[i++] = borrow;
	for (; i < fftsize; i++)
		d[i] = 0.0;
}

/*------------------------------------------------------------------------*/
static void double_to_packed(uint32 *p, int32 max_pwords, int32 pwords, 
				double *d, fastmult_info_t *info) {

	int32 i;
	double recip_base = 1.0 / FFT_BASE;
	double carry = 0.0;
	double round0 = info->round_constant[0];
	double round1 = info->round_constant[1];

	for (i = 0; i < pwords; i++) {
		double digit1, digit2;

		digit1 = d[2 * i] + carry + round0 - round1;
		carry = digit1 * recip_base + round0 - round1;
		digit1 = digit1 - FFT_BASE * carry;
		if (digit1 < 0.0) {
			digit1 += FFT_BASE;
			carry--;
		}

		digit2 = d[2 * i + 1] + carry + round0 - round1;
		carry = digit2 * recip_base + round0 - round1;
		digit2 = digit2 - FFT_BASE * carry;
		if (digit2 < 0.0) {
			digit2 += FFT_BASE;
			carry--;
		}

		if (i < max_pwords)
			p[i] = (int32)digit1 + FFT_BASE * ((int32)digit2);
	}
}

/*------------------------------------------------------------------------*/
static void fht_square(int32 power, uint32 *a, int32 awords,
			uint32 *prod, int32 prod_words,
			fastmult_info_t *info) {

	int32 runlength = 1 << power;
	double *factor1 = (double *)xmalloc(runlength * sizeof(double));

	packed_to_double(a, awords, factor1, runlength);
	fht_forward(factor1, info, power);
	square_hermitian(factor1, runlength);
	fht_inverse(factor1, info, power);
	double_to_packed(prod, prod_words, 2 * awords + 1, factor1, info);
	free(factor1);
}

/*------------------------------------------------------------------------*/
static void fht_mul(int32 power, uint32 *a, int32 awords,
		uint32 *b, int32 bwords, 
		uint32 *prod, int32 prod_words,
		fastmult_info_t *info) {

	int32 runlength = 1 << power;
	double *factor1 = (double *)xmalloc(runlength * sizeof(double));
	double *factor2 = (double *)xmalloc(runlength * sizeof(double));

	packed_to_double(a, awords, factor1, runlength);
	packed_to_double(b, bwords, factor2, runlength);
	fht_forward(factor1, info, power);
	fht_forward(factor2, info, power);
	mul_hermitian(factor1, factor2, runlength);
	fht_inverse(factor1, info, power);
	double_to_packed(prod, prod_words, awords + bwords + 1, factor1, info);
	free(factor1);
	free(factor2);
}

/*------------------------------------------------------------------------*/
void fastmult(uint32 *a, uint32 awords,
		uint32 *b, uint32 bwords,
		uint32 *prod, uint32 prod_words,
		fastmult_info_t *info) {

	/* the product can never be more than awords+bwords
	   in size; however, we use balanced representation,
	   and it's possible for 'a' and 'b' to need one bit
	   more each. If they both need this extra bit, the
	   product will require an extra word. Failing to add
	   the extra word when awords+bwords is a power of 2
	   will mean that the top bit of the product wraps 
	   around to the bottom, corrupting the multiply */

	int32 power = get_log2_fftsize(FFT_WORDS_PER_LIMB * 
					(int32)(awords + bwords + 1));

	create_twiddle(info, power);
	if (a == b && awords == bwords) {
		fht_square(power, a, (int32)awords, 
				prod, (int32)prod_words, info);
	}
	else {
		fht_mul(power, a, (int32)awords, b, (int32)bwords, 
				prod, (int32)prod_words, info);
	}
}
