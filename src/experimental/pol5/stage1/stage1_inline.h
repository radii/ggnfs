#ifndef _STAGE1_INLINE_H_
#define _STAGE1_INLINE_H_

#include "ggnfs.h"

#ifdef __cplusplus
extern "C" {
#endif

/*------------------------------------------------------------------------*/
static INLINE u32 
mp_modinv_1(u32 a, u32 p) {

	/* thanks to the folks at www.mersenneforum.org */

	u32 ps1, ps2, parity, dividend, divisor, rem, q, t;

	q = 1;
	rem = a;
	dividend = p;
	divisor = a;
	ps1 = 1;
	ps2 = 0;
	parity = 0;

	while (divisor > 1) {
		rem = dividend - divisor;
		t = rem - divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) {
			q += ps1; rem = t;
		if (rem >= divisor) {
			q = dividend / divisor;
			rem = dividend % divisor;
			q *= ps1;
		} } } } } } } } }

		q += ps2;
		parity = ~parity;
		dividend = divisor;
		divisor = rem;
		ps2 = ps1;
		ps1 = q;
	}
	
	if (parity == 0)
		return ps1;
	else
		return p - ps1;
}

/*------------------------------------------------------------------------*/
static INLINE u32 
mp_mod64(u64 a, u32 n) {

	u32 ans;

#if (defined(__GNUC__) || defined(__ICL)) && \
	(defined(__i386__) || defined(__x86_64__))
	asm("divl %3  \n\t"
	     : "=d"(ans)
	     : "a"((u32)(a)), "0"((u32)(a >> 32)), "g"(n) : "cc");

#elif defined(_MSC_VER) && !defined(_WIN64)
	__asm
	{
		lea	ecx,a
		mov	eax,[ecx]
		mov	edx,[ecx+4]
		div	n
		mov	ans,edx
	}

#else
	ans = (u32)(a % n);
#endif
	return ans;
}

/*------------------------------------------------------------------------*/
static INLINE u32 
mp_modmul_1(u32 a, u32 b, u32 n) {
	u64 acc = (u64)a * (u64)b;
	return mp_mod64(acc, n);
}

/*------------------------------------------------------------------------*/
static INLINE u32 
mp_expo_1(u32 a, u32 b, u32 n) {

	u32 res = 1;
	while (b) {
		if (b & 1)
			res = mp_modmul_1(res, a, n);
		a = mp_modmul_1(a, a, n);
		b = b >> 1;
	}
	return res;
}

/*------------------------------------------------------------------------*/
static INLINE void
mpz_get_ull_64(uint64_t * targptr, mpz_t src)
{
	if (src->_mp_size == 0)
		*targptr = 0;
	else if (src->_mp_size == 1)
		*targptr = (uint64_t) (src->_mp_d[0]);
	else {
		*targptr = (uint64_t) (src->_mp_d[1]);
		*targptr <<= 32;
		*targptr += (uint64_t) (src->_mp_d[0]);
	}
}

/*------------------------------------------------------------------------*/
#ifdef HAVE_FLOAT64
static INLINE long double
mpz_get_ld(mpz_t src)
{				/* long double has 64bit precision */
	long double res;
	int nb;

	if (mpz_sgn(src) < 0)
		complain("mpz_get_ld: negative\n");
	nb = mpz_size(src) - 1;
	if (nb < 0)
		return 0.L;
	res = (long double) mpz_getlimbn(src, nb);
	nb--;
	if (nb < 0)
		return res;
	res *= 4294967296.L;
	res += (long double) mpz_getlimbn(src, nb);
	nb--;
	if (nb < 0)
		return res;
	res *= 4294967296.L;
	res += (long double) mpz_getlimbn(src, nb);
	while (nb > 0) {
		res *= 4294967296.L;
		nb--;
	}
	return res;
}
#endif

/*------------------------------------------------------------------------*/
static INLINE void 
ulladdmul(uint64_t * resptr, unsigned int ulf, uint64_t * ullfptr)
{
#ifdef HAVE_ASM_INTEL
#if defined(_MSC_VER)
	__asm {
		mov esi, ullfptr 
		mov edi, resptr 
		mov ecx, ulf 
		mov eax,[esi]
		mul ecx 
        add[edi], eax 
		adc[edi + 4], edx 
		mov eax,[esi + 4]
		mul ecx 
		add[edi + 4], eax
	}
#else
	asm volatile("movl (%%esi),%%eax\n"
		     "mull %%ecx\n"
		     "addl %%eax,(%%edi)\n"
		     "adcl %%edx,4(%%edi)\n"
		     "movl 4(%%esi),%%eax\n"
		     "mull %%ecx\n"
		     "addl %%eax,4(%%edi)"
		     :
		     :"S"(ullfptr), "D"(resptr),
		     "c"(ulf):"%edx", "%eax", "cc");
#endif
#else
	uint64_t res, h;

	res = *resptr;
	h = (*ullfptr) * ((uint64_t) ulf);
	res += h;
	*resptr = res;
#endif
}

/*------------------------------------------------------------------------*/
static INLINE void
ull_mulh(uint64_t * resptr, uint64_t * ullf1, uint64_t * ullf2)
{
#ifdef HAVE_ASM_INTEL
#if defined(_MSC_VER)
	__asm {
		mov esi, ullf1 
		mov edi, resptr 
		mov ecx, ullf2 
		mov eax,[esi + 4]
		mov ebx,[ecx + 4]
		mul ebx 
		mov[edi], eax 
		mov[edi + 4], edx 
		mov eax,[esi]
		mul ebx 
		add[edi], edx 
		adc[edi + 4], 0 
		mov ebx,[ecx]
		mov eax,[esi + 4]
		mul ebx 
		add[edi], edx 
		add[edi + 4], 0
	}
#else
	asm volatile("movl 4(%%esi),%%eax\n"
		     "movl 4(%%ecx),%%ebx\n"
		     "mull %%ebx\n"
		     "movl %%eax,(%%edi)\n"
		     "movl %%edx,4(%%edi)\n"
		     "movl (%%esi),%%eax\n"
		     "mull %%ebx\n"
		     "addl %%edx,(%%edi)\n"
		     "adcl $0,4(%%edi)\n"
		     "movl (%%ecx),%%ebx\n"
		     "movl 4(%%esi),%%eax\n"
		     "mull %%ebx\n"
		     "addl %%edx,(%%edi)\n"
		     "addl $0,4(%%edi)"
		     :
		     :"S"(ullf1), "D"(resptr),
		     "c"(ullf2):"%edx", "%eax", "%ebx", "cc");
#endif
#else
	uint64_t res, f1, f2;
	uint64_t f10, f11, f20, f21, h;

	f1 = *ullf1;
	f2 = *ullf2;
	f10 = f1 & 0xffffffffULL;
	f11 = f1 >> 32;
	f20 = f2 & 0xffffffffULL;
	f21 = f2 >> 32;
	res = f11 * f21;
	h = ((f10 * f21) >> 32) + ((f11 * f20) >> 32);
	res += h;
	*resptr = res;
#endif
}

/*------------------------------------------------------------------------*/
#ifdef HAVE_FLOAT64
static INLINE uint64_t
dull(long double d)
{
	long double dh;

	dh = d - (long double) ((int) d);
	dh *= 18446744073709551616.L;
	return (uint64_t) dh;
}
#endif


#ifdef __cplusplus
}
#endif

#endif /* !_STAGE1_INLINE_H_ */
