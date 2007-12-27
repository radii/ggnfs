
#include "profile.h"
#include "if.h"

/*--------------------------------------------------------------------*/
u64 read_clock(void) {

#if (defined(__GNUC__) || defined(__ICL)) && \
	(defined(__i386__) || defined(__x86_64__))
	u32 lo, hi;
	asm("rdtsc":"=d"(hi),"=a"(lo));
	return (u64)hi << 32 | lo;

#elif defined(_MSC_VER)
	LARGE_INTEGER ret;
	QueryPerformanceCounter(&ret);
	return ret.QuadPart;

#else
	struct timeval thistime;   
	gettimeofday(&thistime, NULL);
	return thistime.tv_sec * 1000000 + thistime.tv_usec;
#endif
}

/*--------------------------------------------------------------------*/
void profile_init(profile_t *p, u32 num_events) {

	p->num_events = num_events;
	p->total_time = -read_clock();
	p->events = (u64 *)xcalloc(num_events, sizeof(u64));
}

/*--------------------------------------------------------------------*/
void profile_done(profile_t *p) {
	p->total_time += read_clock();
}

/*--------------------------------------------------------------------*/
void profile_free(profile_t *p) {
	free(p->events);
}

/*--------------------------------------------------------------------*/
void profile_print(profile_t *p, u32 event) {

	printf("%" PRIu64 " (%.2lf%%)", p->events[event],
			100.0 * p->events[event] / p->total_time);
}
