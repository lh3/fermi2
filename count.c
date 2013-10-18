#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "kvec.h"
#include "rld0.h"

typedef struct {
	int64_t k, l;
	int d, c;
} pair64_t;

int main_count(int argc, char *argv[])
{
	rld_t *e;
	int c, min_occ = 100, depth = 51;
	pair64_t *p;
	kvec_t(pair64_t) stack = {0,0,0};
	char *str;

	while ((c = getopt(argc, argv, "d:o:")) >= 0) {
		if (c == 'd') depth = atol(optarg);
		else if (c == 'o') min_occ = atol(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: fermi2 count [-d depth=%d] [-o minOcc=%d] in.fmd\n", depth, min_occ);
		return 1;
	}

	if ((depth&1) == 0) ++depth;
	str = calloc(1, depth + 1);
	e = rld_restore(argv[optind]);
	kv_pushp(pair64_t, stack, &p);
	p->k = 0, p->l = e->mcnt[0], p->d = p->c = 0;
	while (stack.n) {
		uint64_t ok[6], ol[6];
		int a, end;
		pair64_t top = kv_pop(stack);
		if (top.d > 0) str[depth - top.d] = "$ACGTN"[top.c];
		if (top.d == depth) {
			fputs(str, stdout); fputc('\t', stdout);
			printf("%ld\n", (long)(top.l - top.k));
		} else {
			rld_rank2a(e, top.k, top.l, ok, ol);
			end = top.d == depth>>1? 2 : 4;
			for (a = 1; a <= end; ++a) {
				if (ol[a] - ok[a] < min_occ) continue;
				kv_pushp(pair64_t, stack, &p);
				p->k = e->cnt[a] + ok[a];
				p->l = e->cnt[a] + ol[a];
				p->d = top.d + 1;
				p->c = a;
			}
		}
	}
	free(stack.a);
	rld_destroy(e);
	free(str);
	return 0;
}
