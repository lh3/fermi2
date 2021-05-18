#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include "fermi2.h"
#include "ketopt.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#define MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))

int fm_extend_to(const rld_t *e, const char *s, int x, int min_ext, int sat_occ, int *occ)
{
	extern unsigned char seq_nt6_table[128];
	int i, k;
	uint64_t l = 0, u = e->mcnt[0];
	*occ = 0;
	if (x + 1 - min_ext < 0) return 0;
	for (i = x, k = 0; i >= 0; --i) {
		int c = (uint8_t)s[i];
		uint64_t l0 = l, u0 = u;
		c = c < 6? c : c < 128? seq_nt6_table[c] : 5;
		l = e->cnt[c] + rld_rank11(e, l, c);
		u = e->cnt[c] + rld_rank11(e, u, c);
		if (l >= u) break; // can't be extended
		++k;
		if (sat_occ > 0) {
			if (k > min_ext && u - l < sat_occ) {
				--k, l = l0, u = u0;
				break;
			}
		} else if (k == min_ext) {
			break;
		}
	}
	if (l < u) *occ = u - l;
	return k;
}

typedef struct {
	int st, en;
	int occ;
} fm_icnt_t;

fm_icnt_t *fm_kprof(const rld_t *e, const char *s, int min_ext, int sat_occ, int *n_)
{
	int i, x, len, k = 0;
	fm_icnt_t *a;
	len = strlen(s);
	MALLOC(a, len);
	for (x = len - 1; x >= min_ext - 1;) {
		int y, occ;
		y = fm_extend_to(e, s, x, min_ext, sat_occ, &occ);
		a[k].st = x + 1 - (y >= min_ext? y : min_ext);
		a[k].en = x + 1;
		a[k++].occ = occ;
		x -= (y >= min_ext? y - min_ext : 0) + 1;
	}
	*n_ = k;
	for (i = 0; i < k>>1; ++i) {
		fm_icnt_t t = a[i];
		a[i] = a[k - i - 1]; a[k - i - 1] = t;
	}
	return a;
}

int main_kprof(int argc, char *argv[])
{
	ketopt_t o = KETOPT_INIT;
	int c, min_ext = 61, sat_occ = 0;
	rld_t *e;
	kseq_t *ks;
	gzFile fp;

	while ((c = ketopt(&o, argc, argv, 1, "k:c:", 0)) >= 0) {
		if (c == 'k') min_ext = atoi(o.arg);
		else if (c == 'c') sat_occ = atoi(o.arg);
	}
	if (argc - o.ind < 2) {
		fprintf(stderr, "Usage: fermi2 kprof [options] <in.fmd> <in.fa>\n");
		fprintf(stderr, "Options:\n");
		fprintf(stderr, "  -k INT    min k-mer size [%d]\n", min_ext);
		fprintf(stderr, "  -c INT    occurrence saturation [%d]\n", sat_occ);
		return 1;
	}

	fp = gzopen(argv[o.ind+1], "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the sequence file\n", __func__);
		return 1;
	}
	ks = kseq_init(fp);
	e = rld_restore(argv[o.ind]);

	while (kseq_read(ks) >= 0) {
		int i, n_a;
		fm_icnt_t *a;
		a = fm_kprof(e, ks->seq.s, min_ext, sat_occ, &n_a);
		for (i = 0; i < n_a; ++i)
			fprintf(stdout, "%s\t%d\t%d\t%d\n", ks->name.s, a[i].st, a[i].en, a[i].occ);
		free(a);
	}
	rld_destroy(e);
	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}
