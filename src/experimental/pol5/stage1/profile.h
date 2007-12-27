#ifndef _PROFILE_H_
#define _PROFILE_H_

#include "ggnfs.h"

#if defined (__cplusplus)
extern "C" {
#endif

typedef struct {
	u64 total_time;
	u32 num_events;
	u64 *events;
} profile_t;

u64 read_clock(void);
void profile_init(profile_t *p, u32 num_events);
void profile_done(profile_t *p);
void profile_free(profile_t *p);
void profile_print(profile_t *p, u32 event);

#define PROFILE_START(p,i) (p)->events[i] -= read_clock()
#define PROFILE_STOP(p,i) (p)->events[i] += read_clock()

#if defined (__cplusplus)
}
#endif

#endif
