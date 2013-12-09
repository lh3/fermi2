#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include "kvec.h"
#include "rld0.h"

static int dfs_verbose = 3;

/******************
 *** DFS engine ***
 ******************/

#define DFS_SUF_LEN 5

typedef struct {
	int64_t k, l;
	int d, c;
} elem_t;

typedef struct {
	uint64_t c[6];
} fmint6_t;

typedef void (*fmdfs_f)(void *data, int k, char *path, fmint6_t *size, int *cont);

void fm_dfs_core(int n, rld_t **e, int is_half, int max_k, int suf_len, int suf, fmdfs_f func, void *data)
{
	int i, j, c;
	fmint6_t *size, *tk, *tl;
	elem_t *t;
	char *_path, *path;
	uint64_t ok[6], ol[6];
	kvec_t(elem_t) stack = {0,0,0};

	assert(!is_half || (max_k&1));
	t = alloca(sizeof(elem_t) * n);
	size = alloca(sizeof(fmint6_t) * n);
	tk = alloca(sizeof(fmint6_t) * n);
	tl = alloca(sizeof(fmint6_t) * n);
	_path = alloca(max_k + 2);
	path = _path + 1;
	kv_resize(elem_t, stack, n * max_k * 6);
	// descend
	for (i = 0; i < n; ++i) {
		elem_t *p;
		kv_pushp(elem_t, stack, &p);
		p->k = 0, p->l = e[i]->mcnt[0];
		for (j = 0; j < suf_len; ++j) {
			c = (suf>>j*2&3) + 1;
			rld_rank2a(e[i], p->k, p->l, ok, ol);
			p->k = e[i]->cnt[c] + ok[c];
			p->l = e[i]->cnt[c] + ol[c];
		}
		p->d = suf_len, p->c = (suf>>(suf_len-1)*2&3) + 1;
	}
	for (j = 0; j < suf_len; ++j)
		path[max_k - j - 1] = "ACGT"[suf>>j*2&3];
	path[max_k] = 0;
	// traverse
	while (stack.n) {
		int end, cont = 0x1E;
		for (i = n - 1; i >= 0; --i) t[i] = kv_pop(stack);
		if (t->d > max_k) continue;
		path[max_k - t->d] = "\0ACGTN"[t->c];
		for (i = 0; i < n; ++i) {
			assert(t[i].k < e[i]->mcnt[0]);
			rld_rank2a(e[i], t[i].k, t[i].l, tk[i].c, tl[i].c);
			for (c = 0; c < 6; ++c) {
				size[i].c[c] = tl[i].c[c] - tk[i].c[c];
				if (size[i].c[c] == 0) cont &= ~(1<<c);
			}
		}
		func(data, t->d, path + (max_k - t->d), size, &cont);
		end = is_half && t->d == max_k>>1? 2 : 4;
		for (c = 1; c <= end; ++c) {
			if ((cont>>c&1) == 0) continue;
			for (i = 0; i < n; ++i) {
				elem_t *p;
				kv_pushp(elem_t, stack, &p);
				p->k = e[i]->cnt[c] + tk[i].c[c];
				p->l = e[i]->cnt[c] + tl[i].c[c];
				p->d = t->d + 1;
				p->c = c;
			}
		}
	}
	free(stack.a);
}

typedef struct {
	int n, max_k, is_half, suf_len;
	rld_t **e;
	void *data;
	fmdfs_f func;
} shared_t;

void dfs_worker(void *data, long suf, int tid)
{
	shared_t *d = (shared_t*)data;
	fm_dfs_core(d->n, d->e, d->is_half, d->max_k, d->suf_len, suf, d->func, d->data);
	if (dfs_verbose >= 4)
		fprintf(stderr, "[M::%s] processed suffix %d in thread %d\n", __func__, suf, tid);
}

void fm_dfs(int n, rld_t **e, int max_k, int is_half, fmdfs_f func, void *data, int n_threads)
{
	extern void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
	shared_t d;
	int n_suf;
	d.n = n, d.e = e, d.data = data, d.func = func, d.max_k = max_k, d.is_half = is_half;
	d.suf_len = max_k>>1 < DFS_SUF_LEN? max_k>>1 : DFS_SUF_LEN;
	n_suf = 1<<d.suf_len*2;
	n_threads = n_threads < n_suf? n_threads : n_suf;
	kt_for(n_threads, dfs_worker, &d, n_suf);
}

/*************
 *** Count ***
 *************/

typedef struct {
	int max_k, min_occ;
} dfs_count_t;

static void dfs_count(void *data, int k, char *path, fmint6_t *size, int *cont)
{
	dfs_count_t *d = (dfs_count_t*)data;
	int c;
	uint64_t sum = 0;
	for (c = 0; c < 6; ++c) {
		if (size->c[c] < d->min_occ) *cont &= ~(1<<c);
		sum += size->c[c];
	}
	if (k < d->max_k) return;
	printf("%s\t%ld\n", path, (long)sum);
}

int main_count(int argc, char *argv[])
{
	int c, n_threads = 1;
	dfs_count_t d;
	rld_t *e;
	d.max_k = 51, d.min_occ = 10;
	while ((c = getopt(argc, argv, "k:o:t:")) >= 0) {
		if (c == 'k') d.max_k = atoi(optarg);
		else if (c == 'o') d.min_occ = atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: fermi2 count [-k maxK=%d] [-o minOcc=%d] [-t nThreads=1] <in.rld>\n", d.max_k, d.min_occ);
		return 1;
	}
	e = rld_restore(argv[optind]);
	if (!(d.max_k&1)) ++d.max_k;
	fm_dfs(1, &e, d.max_k, 1, dfs_count, &d, n_threads);
	rld_destroy(e);
	return 0;
}

/***********
 *** div ***
 ***********/

#define DIFF_CLASS 12

typedef struct {
	int max_k, min_het_occ, out_k, min_occ[2];
	uint64_t *tot, *cnt[DIFF_CLASS];
} dfs_diff_t;

static void dfs_diff(void *data, int k, char *path, fmint6_t *size, int *cont)
{
	extern double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);
	dfs_diff_t *d = (dfs_diff_t*)data;
	int c, max_c, max2_c, t;
	uint64_t max, max2, sum[2];
	double left, right, two;
	for (c = 0; c < 6; ++c)
		if (size[0].c[c] < d->min_occ[0] || size[1].c[c] < d->min_occ[1]) *cont &= ~(1<<c);
	__sync_fetch_and_add(&d->tot[k], 1);
	max = max2 = 0; max_c = max2_c = 0; sum[0] = sum[1] = 0;
	for (c = 1; c <= 4; ++c) {
		if (max < size[0].c[c] + size[1].c[c])
			max2 = max, max2_c = max_c, max = size[0].c[c] + size[1].c[c], max_c = c;
		else if (max2 < size[0].c[c] + size[1].c[c])
			max2 = size[0].c[c] + size[1].c[c], max2_c = c;
		sum[0] += size[0].c[c], sum[1] += size[1].c[c];
	}
	if (max2_c > 0) {
		int n[4];
		n[0] = size[0].c[max_c],  n[1] = size[1].c[max_c];
		n[2] = size[0].c[max2_c], n[3] = size[1].c[max2_c];
		if (n[0] == 0 || n[1] == 0 || n[2] == 0 || n[3] == 0) {
			kt_fisher_exact(n[0], n[1], n[2], n[3], &left, &right, &two);
			__sync_fetch_and_add(&d->cnt[2][k], 1);
			t = -1;
			if (two < 1e-4) t = 11;
			else if (two < 2e-4) t = 10;
			else if (two < 5e-4) t = 9;
			else if (two < 1e-3) t = 8;
			else if (two < 2e-3) t = 7;
			else if (two < 5e-3) t = 6;
			else if (two < 1e-2) t = 5;
			else if (two < 2e-2) t = 4;
			else if (two < 5e-2) t = 3;
			if (k == d->out_k && two < .01)
				printf("%s\t%.8f\t%ld\t%d\t%d\t%ld\t%d\t%d\n", path, two, (long)sum[0], n[0], n[2], (long)sum[1], n[1], n[3]);
		} else t = 1;
	} else t = 0; // no alternative base
	if (t >= 0) __sync_fetch_and_add(&d->cnt[t][k], 1);
}

static const char *diff2_label[] = { "noAlt", "noZero", "tested", ".05", ".02", ".01", ".005", ".002", ".001", ".0005", ".0002", ".0001" };

int main_diff2(int argc, char *argv[])
{
	dfs_diff_t d;
	int c, i, k, n_threads = 1;
	rld_t *e[2];
	d.max_k = 51, d.min_het_occ = 3, d.out_k = -1, d.min_occ[0] = d.min_occ[1] = 6;
	while ((c = getopt(argc, argv, "k:o:t:h:v:K:")) >= 0) {
		if (c == 'o') {
			char *p;
			d.min_occ[0] = strtol(optarg, &p, 10);
			d.min_occ[1] = (*p && *(p+1))? strtol(p+1, &p, 10) : d.min_occ[1];
		} else if (c == 'k') d.max_k = atoi(optarg);
		else if (c == 'h') d.min_het_occ = atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
		else if (c == 'v') dfs_verbose = atoi(optarg);
		else if (c == 'K') d.out_k = atoi(optarg);
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: fermi2 diff2 <in1.rld> <in2.rld>\n");
		return 1;
	}
	d.tot = calloc(d.max_k, sizeof(uint64_t));
	for (i = 0; i < DIFF_CLASS; ++i)
		d.cnt[i] = calloc(d.max_k, sizeof(uint64_t));
	e[0] = rld_restore(argv[optind+0]);
	e[1] = rld_restore(argv[optind+1]);
	fm_dfs(2, e, d.max_k, 0, dfs_diff, &d, n_threads);
	rld_destroy(e[0]);
	rld_destroy(e[1]);
	printf("#k\ttotal");
	for (i = 0; i < DIFF_CLASS; ++i) printf("\t%s", diff2_label[i]);
	printf("\n");
	for (k = 1; k < d.max_k; ++k) {
		printf("%d\t%ld", k, (long)d.tot[k]);
		for (i = 0; i < DIFF_CLASS; ++i)
			printf("\t%ld", (long)d.cnt[i][k]);
		printf("\n");
	}
	return 0;
}
