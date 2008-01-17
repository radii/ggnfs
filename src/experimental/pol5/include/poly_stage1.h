#ifndef _POLY_STAGE1_H_
#define _POLY_STAGE1_H_

#include "ggnfs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	mpz_t gmp_N;
	mpz_t gmp_a5_begin;
	mpz_t gmp_a5_end;
	u32 npr_in_p;
	double norm_max;
	u32 p0_limit;
	FILE *outfile;
} poly_stage1_t;

void poly_stage1_init(poly_stage1_t *data);
void poly_stage1_free(poly_stage1_t *data);
s32 poly_stage1_run(poly_stage1_t *data);

#ifdef __cplusplus
}
#endif

#endif /* !_POLY_STAGE1_H_ */
