#ifndef _POLY_STAGE2_H_
#define _POLY_STAGE2_H_

#include "ggnfs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	mpz_t gmp_N;
	double area;
	double bound0;
	double bound1;
	double min_e;
	double max_norm_1;
	double max_norm_2;
	u32 p_bound;
	FILE *infile;
	FILE *outfile;
} poly_stage2_t;

void poly_stage2_init(poly_stage2_t *data);
void poly_stage2_free(poly_stage2_t *data);
s32 poly_stage2_run(poly_stage2_t *data);

#ifdef __cplusplus
}
#endif

#endif /* !_POLY_STAGE2_H_ */
