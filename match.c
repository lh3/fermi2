#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>
#include <zlib.h>
#include "fermi2.h"
#include "kvec.h"
#include "kstring.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

int kvsprintf(kstring_t *s, const char *fmt, va_list ap)
{
	va_list args;
	int l;
	va_copy(args, ap);
	l = vsnprintf(s->s + s->l, s->m - s->l, fmt, args); // This line does not work with glibc 2.0. See `man snprintf'.
	va_end(args);
	if (l + 1 > s->m - s->l) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
		va_copy(args, ap);
		l = vsnprintf(s->s + s->l, s->m - s->l, fmt, args);
		va_end(args);
	}
	s->l += l;
	return l;
}

int ksprintf(kstring_t *s, const char *fmt, ...)
{
	va_list ap;
	int l;
	va_start(ap, fmt);
	l = kvsprintf(s, fmt, ap);
	va_end(ap);
	return l;
}

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
extern void seq_revcomp6(int l, unsigned char *s);
extern void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);

typedef struct {
	rldintv_v curr, prev;
	fmdsmem_v smem;
	kstring_t str, cmp[2];
} thrmem_t;

typedef struct {
	const rld_t *e;
	const fmsa_t *sa;
	int max_sa_occ, min_occ;
	int partial, discovery, kmer;

	int n_threads;
	thrmem_t *mem;

	int n_seqs, m_seqs;
	char **name, **seq, **qual, **out;
} global_t;

static int find_end(const rld_t *e, int l_seq, const char *seq, int64_t size, int start, int is_back)
{
	int i, step = is_back? -1 : 1, end = is_back? -1 : l_seq;
	int64_t l = 0, u = e->mcnt[0];
	for (i = start; i != end; i += step) {
		int c = seq[i];
		if (!is_back) c = fmd_comp(c);
		l = e->cnt[c] + rld_rank11(e, l, c);
		u = e->cnt[c] + rld_rank11(e, u, c);
		if (u - l <= size) break;
	}
	if (u - l == size) i += step;
	return abs(i - start);
}

static void discover(const rld_t *e, const fmdsmem_t *q, const fmdsmem_t *p, int l_seq, const char *seq, const char *qual, kstring_t *s, kstring_t cmp[2])
{
	int start, end, start0, end0, i, occ[2];
	// find the coordinate from which the extension will be applied
	if (q == 0 && p == 0) {
		start = 0; end = l_seq;
		occ[0] = occ[1] = 0;
	} else if (q == 0) {
		if (p->ik.info>>32 == 0) return; // no novel allele
		start = 0; end = p->ik.info>>32;
		occ[0] = 0; occ[1] = p->ik.x[2];
	} else if (p == 0) {
		if ((uint32_t)q->ik.info == l_seq) return; // no novel allele
		start = (uint32_t)q->ik.info; end = l_seq;
		occ[0] = q->ik.x[2]; occ[1] = 0;
	} else {
		start = (uint32_t)q->ik.info, end = p->ik.info>>32;
		if (start >= end && (q->ok[1][0].x[2] == q->ik.x[2] || p->ik.x[2] == p->ok[0][0].x[2])) return; // TODO: is this the desired behavior???
		occ[0] = q->ik.x[2]; occ[1] = p->ik.x[2];
	}
	if (start < end) {
		for (i = start; i < end; ++i)
			if (seq[i] < 5) break;
		if (i == end) return; // the gap is filled with "N"
	}
	// extend
	start0 = start, end0 = end;
	start = start == 0?   0     : start - find_end(e, l_seq, seq, q->ik.x[2], start - 1, 1);
	end   = end == l_seq? l_seq : end   + find_end(e, l_seq, seq, p->ik.x[2], end,       0);
	// find which sequence to output
	cmp[0].l = cmp[1].l = 0;
	kputsn(&seq[start], end - start, &cmp[0]);
	kputsn(&seq[start], end - start, &cmp[1]);
	seq_revcomp6(cmp[1].l, (uint8_t*)cmp[1].s);
	// print
	if (strcmp(cmp[0].s, cmp[1].s) <= 0) {
		ksprintf(s, "NS\t%d\t%d\t%d\t%d\t%ld\t%ld\t+\t", start, start0 - start, end0 - start0, end - end0, (long)occ[0], (long)occ[1]);
		for (i = 0; i < cmp[0].l; ++i)
			kputc("$ACGTN"[(int)cmp[0].s[i]], s);
		kputc('\t', s);
		if (qual) kputsn(&qual[start], end - start, s);
		else kputc('*', s);
	} else {
		ksprintf(s, "NS\t%d\t%d\t%d\t%d\t%ld\t%ld\t-\t", start, end - end0, end0 - start0, start0 - start, (long)occ[1], (long)occ[0]);
		for (i = 0; i < cmp[0].l; ++i)
			kputc("$ACGTN"[(int)cmp[1].s[i]], s);
		kputc('\t', s);
		if (qual) {
			for (i = end - 1; i >= start; --i)
				kputc(qual[i], s);
		} else kputc('*', s);
	}
	kputc('\n', s);
}

static void worker(void *data, long jid, int tid)
{
	global_t *g = (global_t*)data;
	thrmem_t *m = &g->mem[tid];
	char *seq = g->seq[jid], *qual = g->qual? g->qual[jid] : 0;
	int l_seq;

	l_seq = strlen(seq);
	seq_char2nt6(l_seq, (uint8_t*)seq);
	m->str.l = 0;
	ksprintf(&m->str, "SQ\t%s\t%d\n", g->name[jid], l_seq);
	if (!g->partial) { // full-length match
		int64_t k, l, u;
		fm_exact(g->e, seq, &l, &u);
		if (l < u) {
			ksprintf(&m->str, "EM\t0\t%d\t%ld", l_seq, (long)(u - l));
			if (g->sa && u - l <= g->max_sa_occ) {
				for (k = l; k < u; ++k) {
					int64_t idx, i;
					i = fm_sa(g->e, g->sa, k, &idx);
					ksprintf(&m->str, "\t%ld:%ld", (long)idx, (long)i);
				}
			}
			kputc('\n', &m->str);
		}
	} else { // SMEM
		size_t i;
		int64_t k;
		fmd_smem(g->e, (uint8_t*)seq, &m->smem, g->min_occ, &m->curr, &m->prev);
		if (g->discovery) {
			int pre;
			for (i = 0, pre = -1; i < m->smem.n; ++i) {
				fmdsmem_t *p = &m->smem.a[i];
				int start = p->ik.info>>32, end = (uint32_t)p->ik.info;
				if (end - start < g->kmer) continue; // skip short SMEMs
				rld_extend(g->e, &p->ik, p->ok[1], 0);
				discover(g->e, pre < 0? 0 : &m->smem.a[pre], p, l_seq, seq, qual, &m->str, m->cmp);
				pre = i;
			}
			discover(g->e, &m->smem.a[pre], 0, l_seq, seq, qual, &m->str, m->cmp);
		} else {
			for (i = 0; i < m->smem.n; ++i) {
				fmdsmem_t *p = &m->smem.a[i];
				ksprintf(&m->str, "EM\t%u\t%u\t%ld", (uint32_t)(p->ik.info>>32), (uint32_t)p->ik.info, (long)p->ik.x[2]);
				if (g->sa && p->ik.x[2] < g->max_sa_occ) {
					for (k = 0; k < p->ik.x[2]; ++k) {
						int64_t idx, j;
						j = fm_sa(g->e, g->sa, p->ik.x[0] + k, &idx);
						ksprintf(&m->str, "\t%ld:%ld", (long)idx, (long)j);
					}
				}
				kputc('\n', &m->str);
			}
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
	g.max_sa_occ = 10, g.min_occ = 1, g.n_threads = 1, g.kmer = 61;
	while ((c = getopt(argc, argv, "Mdps:m:O:b:t:k:")) >= 0) {
		if (c == 'M') use_mmap = 1;
		else if (c == 's') fn_sa = optarg;
		else if (c == 'm') g.max_sa_occ = atoi(optarg);
		else if (c == 'p') g.partial = 1;
		else if (c == 'd') g.discovery = g.partial = 1;
		else if (c == 'O') g.min_occ = atoi(optarg);
		else if (c == 't') g.n_threads = atoi(optarg);
		else if (c == 'k') g.kmer = atoi(optarg), g.discovery = g.partial = 1;
		else if (c == 'b') batch_size = atoi(optarg);
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   fermi2 match [options] <index.fmd> <seq.fa>\n\n");
		fprintf(stderr, "Options: -p        find SMEMs (req. both strands in one index)\n");
		fprintf(stderr, "         -d        discovery novel alleles (force -p)\n");
		fprintf(stderr, "         -k INT    k-mer length in the discovery mode (force -d) [%d]\n", g.kmer);
		fprintf(stderr, "         -t INT    number of threads [%d]\n", g.n_threads);
		fprintf(stderr, "         -b INT    batch size [%d]\n", batch_size);
		fprintf(stderr, "         -s FILE   sampled suffix array [null]\n");
		fprintf(stderr, "         -m INT    show coordinate if the number of hits is no more than INT [%d]\n", g.max_sa_occ);
		fprintf(stderr, "         -s INT    min occurrences [%d]\n", g.min_occ);
		fprintf(stderr, "\n");
		fprintf(stderr, "Output format:\n\n");
		fprintf(stderr, "    SQ  seqName  seqLen\n");
		fprintf(stderr, "    EM  start  end  occurrence [positions]\n");
		fprintf(stderr, "    NS  start  leftLen  diffLen  rightLen  leftOcc  rightOcc  strand  seq  qual\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "  At an 'NS' line, the length of 'seq' always equals leftLen+diffLen+rightLen.\n");
		fprintf(stderr, "\n");
		return 1;
	}

	fp = gzopen(argv[optind+1], "r");
	if (fp == 0) {
		fprintf(stderr, "[E::%s] failed to open the sequence file\n", __func__);
		return 1;
	}
	g.e = use_mmap? rld_restore_mmap(argv[optind]) : rld_restore(argv[optind]);
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
		free(g.mem[i].cmp[0].s); free(g.mem[i].cmp[1].s); free(g.mem[i].str.s);
	}
	free(g.name); free(g.seq); free(g.qual); free(g.out); free(g.mem);
	if (g.sa) fm_sa_destroy((fmsa_t*)g.sa);
	rld_destroy((rld_t*)g.e);
	gzclose(fp);
	return 0;
}
