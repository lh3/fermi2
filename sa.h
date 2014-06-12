#ifndef FERMI2_SA_H
#define FERMI2_SA_H

#include <stdint.h>
#include "rld0.h"

typedef struct {
	int ms;
	int ss;
	int64_t m, n_ssa;
	uint64_t *r2i; // rank -> index
	uint64_t *ssa; // sampled suffix array
} fmsa_t;

fmsa_t *fm_sa_gen(const rld_t *e, int ssa_shift, int n_threads);
void fm_sa_dump(const fmsa_t *sa, const char *fn);
fmsa_t *fm_sa_restore(const char *fn);
void fm_sa_destroy(fmsa_t *sa);

int64_t fm_sa(const rld_t *e, const fmsa_t *sa, int k, int64_t *si);

#endif
