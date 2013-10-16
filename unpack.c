#include <stdio.h>
#include "rld0.h"
#include "kstring.h"

int64_t fm_retrieve(const rld_t *e, uint64_t x, kstring_t *s)
{
	uint64_t k = x, *ok;
	ok = alloca(8 * e->asize);
	s->l = 0;
	while (1) {
		int c = rld_rank1a(e, k + 1, ok);
		k = e->cnt[c] + ok[c] - 1;
		if (c == 0) return k;
		kputc(c, s);
	}
}

#include <unistd.h>

int main_unpack(int argc, char *argv[])
{
	int64_t i, k;
	int c, j, tmp, from_stdin = 0;
	rld_t *e;
	kstring_t str = {0,0,0};

	from_stdin = !isatty(fileno(stdin));
	while ((c = getopt(argc, argv, "")) >= 0);
	if (!from_stdin && optind == argc) {
		fprintf(stderr, "Usage: fermi2 unpack [reads.rld]\n");
		return 1;
	}
	e = rld_restore(from_stdin? "-" : argv[optind]);
	for (i = 0; i < e->mcnt[1]; ++i) {
		k = fm_retrieve(e, i, &str);
		for (j = 0; j < str.l; ++j)
			str.s[j] = "$ACGTN"[(int)str.s[j]];
		for (j = 0; j < str.l>>1; ++j)
			tmp = str.s[j], str.s[j] = str.s[str.l-1-j], str.s[str.l-1-j] = tmp;
		fwrite(str.s, 1, str.l, stdout);
		printf("\t%ld\n", (long)k);
	}
	free(str.s);
	rld_destroy(e);
	return 0;
}
