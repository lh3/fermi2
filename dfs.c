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
{ // this routine is similar to fmc_collect1()
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
		fprintf(stderr, "[M::%s] processed suffix %ld in thread %d\n", __func__, suf, tid);
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
	int len, min_occ;
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
	if (k < d->len) return;
	printf("%s\t%ld\n", path, (long)sum);
}

int main_count(int argc, char *argv[])
{
	int c, n_threads = 1;
	dfs_count_t d;
	rld_t *e;
	d.len = 51, d.min_occ = 10;
	while ((c = getopt(argc, argv, "k:o:t:")) >= 0) {
		if (c == 'k') d.len = atoi(optarg);
		else if (c == 'o') d.min_occ = atoi(optarg);
		else if (c == 't') n_threads = atoi(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: fermi2 count [-k len=%d] [-o minOcc=%d] [-t nThreads=1] <in.rld>\n", d.len, d.min_occ);
		return 1;
	}
	e = rld_restore(argv[optind]);
	if (!(d.len&1)) {
		++d.len;
		if (dfs_verbose >= 2)
			fprintf(stderr, "[W::%s] %d is an even number; change k to %d\n", __func__, d.len-1, d.len);
	}
	fm_dfs(1, &e, d.len, 1, dfs_count, &d, n_threads);
	rld_destroy(e);
	return 0;
}
