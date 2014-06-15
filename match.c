#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "fermi2.h"
#include "kvec.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

int fmd_smem1_core(const rld_t *e, int min_occ, int len, const uint8_t *q, int x, fmdsmem_v *mem, rldintv_v *curr, rldintv_v *prev)
{ // for more comments, see bwa/bwt.c
	int i, j, c, ret;
	rldintv_t ik, ok[6];
	rldintv_v *swap;
	size_t oldn = mem->n;

	fmd_set_intv(e, q[x], ik);
	ik.info = x + 1;
	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward extension
		c = fmd_comp(q[i]);
		rld_extend(e, &ik, ok, 0);
		if (ok[c].x[2] != ik.x[2]) {
			kv_push(rldintv_t, *curr, ik);
			if (ok[c].x[2] < min_occ) break;
		}
		ik = ok[c]; ik.info = i + 1;
	}
	if (i == len) kv_push(rldintv_t, *curr, ik);
	kv_reverse(rldintv_t, *curr, 0);
	ret = curr->a[0].info;
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) {
		c = i < 0? 0 : q[i];
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			rldintv_t *p = &prev->a[j];
			rld_extend(e, p, ok, 1);
			if (c == 0 || ok[c].x[2] < min_occ) {
				if (curr->n == 0) {
					if (mem->n == oldn || i + 1 < mem->a[mem->n-1].ik.info>>32) {
						fmdsmem_t *q;
						kv_pushp(fmdsmem_t, *mem, &q);
						q->ik = *p; q->ik.info |= (uint64_t)(i + 1)<<32;
						memcpy(q->ok[0], ok, 6 * sizeof(rldintv_t));
					}
				}
			} else if (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]) {
				ok[c].info = p->info;
				kv_push(rldintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}

	kv_reverse(fmdsmem_t, *mem, oldn);
	return ret;
}

int fmd_smem(const rld_t *e, const uint8_t *q, fmdsmem_v *mem, int min_occ, rldintv_v *curr, rldintv_v *prev)
{
	int x = 0, len;
	mem->n = 0;
	len = strlen((char*)q);
	do {
		x = fmd_smem1_core(e, min_occ, len, q, x, mem, curr, prev);
	} while (x < len);
	return mem->n;
}

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

extern void seq_char2nt6(int l, unsigned char *s);

int main_match(int argc, char *argv[])
{
	int c, use_mmap = 0, max_sa_occ = 10, partial = 0, min_occ = 1;
	gzFile fp;
	rld_t *e;
	char *fn_sa = 0;
	fmsa_t *sa = 0;
	kseq_t *ks;
	rldintv_v curr = {0,0,0}, prev = {0,0,0};
	fmdsmem_v smem = {0,0,0};

	while ((c = getopt(argc, argv, "Mps:m:O:")) >= 0) {
		if (c == 'M') use_mmap = 1;
		else if (c == 's') fn_sa = optarg;
		else if (c == 'm') max_sa_occ = atoi(optarg);
		else if (c == 'p') partial = 1;
		else if (c == 'O') min_occ = atoi(optarg);
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
		printf("SQ\t%s\t%d\n", ks->name.s, ks->seq.l);
		seq_char2nt6(ks->seq.l, (uint8_t*)ks->seq.s);
		if (!partial) {
			int64_t k, l, u;
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
		} else {
			size_t i;
			int64_t k;
			fmd_smem(e, (uint8_t*)ks->seq.s, &smem, min_occ, &curr, &prev);
			for (i = 0; i < smem.n; ++i) {
				fmdsmem_t *p = &smem.a[i];
				printf("EM\t%u\t%u\t%ld", (uint32_t)(p->ik.info>>32), (uint32_t)p->ik.info, (long)p->ik.x[2]);
				if (sa && p->ik.x[2] < max_sa_occ) {
					for (k = 0; k < p->ik.x[2]; ++k) {
						int64_t idx, j;
						j = fm_sa(e, sa, p->ik.x[0] + k, &idx);
						printf("\t%ld:%ld", (long)idx, (long)j);
					}
				}
				putchar('\n');
			}
		}
		puts("//");
	}
	kseq_destroy(ks);
	free(curr.a); free(prev.a); free(smem.a);

	if (sa) fm_sa_destroy(sa);
	rld_destroy(e);
	gzclose(fp);
	return 0;
}
