#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "rld0.h"
#include "fermi2.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

void fm_exact(const rld_t *e, const char *s, int64_t *_l, int64_t *_u)
{
	extern unsigned char seq_nt6_table[128];
	int64_t i, l = 0, u = e->mcnt[0];
	for (i = strlen(s) - 1; i >= 0; --i) {
		int c = (uint8_t)s[i];
		c = c < 6? c : c < 128? seq_nt6_table[c] : 5;
		l = e->cnt[c] + rld_rank11(e, l, c);
		u = e->cnt[c] + rld_rank11(e, u, c);
		if (l >= u) break;
	}
	*_l = l, *_u = u;
}

int main_match(int argc, char *argv[])
{
	int c, use_mmap = 0, max_sa_occ = 10;
	gzFile fp;
	rld_t *e;
	char *fn_sa = 0;
	fmsa_t *sa = 0;
	kseq_t *ks;

	while ((c = getopt(argc, argv, "Ms:m:")) >= 0) {
		if (c == 'M') use_mmap = 1;
		else if (c == 's') fn_sa = optarg;
		else if (c == 'm') max_sa_occ = atoi(optarg);
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: fermi2 match <index.fmd> <seq.fa>\n");
		return 1;
	}

	fp = gzopen(argv[optind+1], "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the sequence file\n", __func__);
		return 1;
	}
	e = rld_restore(argv[optind]);
	if (e == 0) {
		fprintf(stderr, "[E::%s] failed to open the index file\n", __func__);
		gzclose(fp);
		return 1;
	}
	if (fn_sa) sa = fm_sa_restore(fn_sa);
	if (fn_sa && sa == 0) {
		fprintf(stderr, "[E::%s] failed to open the sampled SA file\n", __func__);
		rld_destroy(e);
		gzclose(fp);
		return 1;
	}

	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		int64_t k, l, u;
		printf("SQ\t%s\t%d\n", ks->name.s, ks->seq.l);
		fm_exact(e, ks->seq.s, &l, &u);
		if (l < u) {
			printf("EM\t0\t%d\t%ld", ks->seq.l, (long)(u - l));
			if (sa && u - l <= max_sa_occ) {
				for (k = l; k < u; ++k) {
					int64_t idx, i;
					i = fm_sa(e, sa, k, &idx);
					printf("\t%ld:%ld", (long)idx, (long)i);
				}
			}
			putchar('\n');
		}
		puts("//");
	}
	kseq_destroy(ks);

	if (sa) fm_sa_destroy(sa);
	rld_destroy(e);
	gzclose(fp);
	return 0;
}
