#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "fermi2.h"
#include "kvec.h"
#include "kstring.h"
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
	if (ik.x[2] == 0) return x + 1;
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
extern void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

typedef struct {
	rldintv_v curr, prev;
	fmdsmem_v smem;
	kstring_t str;
} thrmem_t;

typedef struct {
	const rld_t *e;
	const fmsa_t *sa;
	int max_sa_occ, partial, min_occ;

	int n_threads;
	thrmem_t *mem;

	int n_seqs, m_seqs;
	char **name, **seq, **qual, **out;
} global_t;

static void worker(void *data, long jid, int tid)
{
	global_t *g = (global_t*)data;
	thrmem_t *m = &g->mem[tid];
	char *seq = g->seq[jid];
	int l_seq;

	l_seq = strlen(seq);
	seq_char2nt6(l_seq, (uint8_t*)seq);
	m->str.l = 0;
	kputsn("SQ\t", 3, &m->str); kputs(g->name[jid], &m->str); kputc('\t', &m->str); kputw(l_seq, &m->str); kputc('\n', &m->str);
	if (!g->partial) { // full-length match
		int64_t k, l, u;
		fm_exact(g->e, seq, &l, &u);
		if (l < u) {
			kputsn("EM\t0\t", 5, &m->str); kputl(l_seq, &m->str); kputc('\t', &m->str); kputl(u - l, &m->str);
			if (g->sa && u - l <= g->max_sa_occ) {
				for (k = l; k < u; ++k) {
					int64_t idx, i;
					i = fm_sa(g->e, g->sa, k, &idx);
					kputc('\t', &m->str); kputl(idx, &m->str); kputc(':', &m->str); kputl(i, &m->str);
				}
			}
			kputc('\n', &m->str);
		}
	} else { // SMEM
		size_t i;
		int64_t k;
		fmd_smem(g->e, (uint8_t*)seq, &m->smem, g->min_occ, &m->curr, &m->prev);
		for (i = 0; i < m->smem.n; ++i) {
			fmdsmem_t *p = &m->smem.a[i];
			kputsn("EM\t", 3, &m->str); kputw(p->ik.info>>32, &m->str); kputc('\t', &m->str); kputw((uint32_t)p->ik.info, &m->str);
			kputc('\t', &m->str); kputl(p->ik.x[2], &m->str);
			if (g->sa && p->ik.x[2] < g->max_sa_occ) {
				for (k = 0; k < p->ik.x[2]; ++k) {
					int64_t idx, j;
					j = fm_sa(g->e, g->sa, p->ik.x[0] + k, &idx);
					kputc('\t', &m->str); kputl(idx, &m->str); kputc(':', &m->str); kputl(j, &m->str);
				}
			}
			kputc('\n', &m->str);
		}
	}
	kputsn("//", 2, &m->str);
	free(g->qual[jid]); free(g->seq[jid]); free(g->name[jid]);
	g->out[jid] = strdup(m->str.s);
}

int main_match(int argc, char *argv[])
{
	int i, c, use_mmap = 0, batch_size = 10000000, l_seqs;
	gzFile fp;
	char *fn_sa = 0;
	kseq_t *ks;
	global_t g;

	memset(&g, 0, sizeof(global_t));
	g.max_sa_occ = 10, g.min_occ = 1, g.n_threads = 1;
	while ((c = getopt(argc, argv, "Mps:m:O:b:t:")) >= 0) {
		if (c == 'M') use_mmap = 1;
		else if (c == 's') fn_sa = optarg;
		else if (c == 'm') g.max_sa_occ = atoi(optarg);
		else if (c == 'p') g.partial = 1;
		else if (c == 'O') g.min_occ = atoi(optarg);
		else if (c == 't') g.n_threads = atoi(optarg);
		else if (c == 'b') batch_size = atoi(optarg);
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi2 match [options] <index.fmd> <seq.fa>\n\n");
		fprintf(stderr, "Options: -p        find SMEMs (req. both strands in one index)\n");
		fprintf(stderr, "         -t INT    number of threads [%d]\n", g.n_threads);
		fprintf(stderr, "         -b INT    batch size [%d]\n", batch_size);
		fprintf(stderr, "         -s FILE   sampled suffix array [null]\n");
		fprintf(stderr, "         -m INT    show coordinate if the number of hits is no more than INT [%d]\n", g.max_sa_occ);
		fprintf(stderr, "         -s INT    min occurrences [%d]\n", g.min_occ);
		fprintf(stderr, "\n");
		return 1;
	}

	fp = gzopen(argv[optind+1], "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the sequence file\n", __func__);
		return 1;
	}
	g.e = rld_restore(argv[optind]);
	if (g.e == 0) {
		fprintf(stderr, "[E::%s] failed to open the index file\n", __func__);
		gzclose(fp);
		return 1;
	}
	if (g.partial && (g.e->mcnt[2] != g.e->mcnt[5] || g.e->mcnt[3] != g.e->mcnt[4])) {
		fprintf(stderr, "[E::%s] with '-p', the index must include both strands\n", __func__);
		rld_destroy((rld_t*)g.e);
		gzclose(fp);
		return 1;
	}
	if (fn_sa) g.sa = fm_sa_restore(fn_sa);
	if (fn_sa && g.sa == 0) {
		fprintf(stderr, "[E::%s] failed to open the sampled SA file\n", __func__);
		rld_destroy((rld_t*)g.e);
		gzclose(fp);
		return 1;
	}

	g.mem = calloc(g.n_threads, sizeof(thrmem_t));

	batch_size *= g.n_threads;
	ks = kseq_init(fp);
	l_seqs = 0;
	while (kseq_read(ks) >= 0) {
		if (g.n_seqs == g.m_seqs) {
			g.m_seqs = g.m_seqs? g.m_seqs<<1 : 4;
			g.name = realloc(g.name, g.m_seqs * sizeof(char*));
			g.seq  = realloc(g.seq,  g.m_seqs * sizeof(char*));
			g.qual = realloc(g.qual, g.m_seqs * sizeof(char*));
			g.out  = realloc(g.out,  g.m_seqs * sizeof(char*));
		}
		g.name[g.n_seqs] = strdup(ks->name.s);
		g.seq[g.n_seqs]  = strdup(ks->seq.s);
		g.qual[g.n_seqs] = ks->qual.l? strdup(ks->qual.s) : 0; // these will be free'd in worker
		++g.n_seqs;
		l_seqs += ks->seq.l;
		if (l_seqs >= batch_size) {
			kt_for(g.n_threads, worker, &g, g.n_seqs);
			for (i = 0; i < g.n_seqs; ++i) {
				puts(g.out[i]);
				free(g.out[i]);
			}
			g.n_seqs = l_seqs = 0;
		}
	}
	// the last batch
	kt_for(g.n_threads, worker, &g, g.n_seqs);
	for (i = 0; i < g.n_seqs; ++i) {
		puts(g.out[i]);
		free(g.out[i]);
	}
	kseq_destroy(ks);

	for (i = 0; i < g.n_threads; ++i) {
		free(g.mem[i].curr.a); free(g.mem[i].prev.a); free(g.mem[i].smem.a);
		free(g.mem[i].str.s);
	}
	free(g.name); free(g.seq); free(g.qual); free(g.out);
	if (g.sa) fm_sa_destroy((fmsa_t*)g.sa);
	rld_destroy((rld_t*)g.e);
	gzclose(fp);
	return 0;
}
