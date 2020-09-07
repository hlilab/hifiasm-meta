#include <stdint.h>
#include <zlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "kthread.h"
#include "khashl.h"
#include "kseq.h"
#include "ksort.h"
#include "htab.h"
#include "meta_util.h"
#define __STDC_FORMAT_MACROS 1  // cpp special (ref: https://stackoverflow.com/questions/14535556/why-doesnt-priu64-work-in-this-code)
#include <inttypes.h>  // debug, for printing uint64
#include <time.h>

#define YAK_COUNTER_BITS 12  // note: do not directly modify this, h->pre is not exposed and needs this value to init
#define YAK_COUNTER_BITS1 14  // allow up to 16383
#define YAK_N_COUNTS     (1<<YAK_COUNTER_BITS1)  // used for histogram, so it refers to counter bits, not pre bits
#define YAK_MAX_COUNT    ((1<<YAK_COUNTER_BITS1)-1)

#define HAMT_DIG_KMERRESCUE 50

const unsigned char seq_nt4_table[256] = { // translate ACGT to 0123
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void *ha_flt_tab;
ha_pt_t *ha_idx;

/***************************
 * Yak specific parameters *
 ***************************/

typedef struct {
	int32_t bf_shift, bf_n_hash;
	int32_t k, w, is_HPC;
	int32_t pre;
	int32_t n_thread;
	int64_t chunk_size;
} yak_copt_t;

void yak_copt_init(yak_copt_t *o)
{
	memset(o, 0, sizeof(yak_copt_t));
	o->bf_shift = 0;
	o->bf_n_hash = 4;
	o->k = 31;
	o->w = 1;
	o->pre = YAK_COUNTER_BITS;
	o->n_thread = 4;
	o->chunk_size = 20000000;
}

/************************
 * Blocked bloom filter *
 ************************/

#define YAK_BLK_SHIFT  9 // 64 bytes, the size of a cache line
#define YAK_BLK_MASK   ((1<<(YAK_BLK_SHIFT)) - 1)

typedef struct {
	int n_shift, n_hashes;
	uint8_t *b;
} yak_bf_t;

yak_bf_t *yak_bf_init(int n_shift, int n_hashes)
{
	yak_bf_t *b;
	void *ptr = 0;
	if (n_shift + YAK_BLK_SHIFT > 64 || n_shift < YAK_BLK_SHIFT) return 0;
	CALLOC(b, 1);
	b->n_shift = n_shift;
	b->n_hashes = n_hashes;
	posix_memalign(&ptr, 1<<(YAK_BLK_SHIFT-3), 1ULL<<(n_shift-3));
	b->b = (uint8_t*)ptr;
	bzero(b->b, 1ULL<<(n_shift-3));
	return b;
}

void yak_bf_destroy(yak_bf_t *b)
{
	if (b == 0) return;
	free(b->b); free(b);
}

int yak_bf_insert(yak_bf_t *b, uint64_t hash)
{
	int x = b->n_shift - YAK_BLK_SHIFT;
	uint64_t y = hash & ((1ULL<<x) - 1);
	int h1 = hash >> x & YAK_BLK_MASK;
	int h2 = hash >> b->n_shift & YAK_BLK_MASK;
	uint8_t *p = &b->b[y<<(YAK_BLK_SHIFT-3)];
	int i, z = h1, cnt = 0;
	if ((h2&31) == 0) h2 = (h2 + 1) & YAK_BLK_MASK; // otherwise we may repeatedly use a few bits
	for (i = 0; i < b->n_hashes; z = (z + h2) & YAK_BLK_MASK) {
		uint8_t *q = &p[z>>3], u;
		u = 1<<(z&7);
		cnt += !!(*q & u);
		*q |= u;
		++i;
	}
	return cnt;
}

/********************
 * Count hash table *
 ********************/

#define yak_ct_eq(a, b) ((a)>>YAK_COUNTER_BITS1 == (b)>>YAK_COUNTER_BITS1) // lower bits for counts
#define yak_ct_hash(a) ((a)>>YAK_COUNTER_BITS1)
KHASHL_SET_INIT(static klib_unused, yak_ct_t, yak_ct, uint64_t, yak_ct_hash, yak_ct_eq)

typedef struct {
	yak_ct_t *h;
	yak_bf_t *b;
} ha_ct1_t;

typedef struct {
	int k, pre, n_hash, n_shift;
	uint64_t tot;
	ha_ct1_t *h;
} ha_ct_t;

static ha_ct_t *ha_ct_init(int k, int pre, int n_hash, int n_shift)
{
	ha_ct_t *h;
	int i;
	if (pre < YAK_COUNTER_BITS) return 0;
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1<<h->pre);
	for (i = 0; i < 1<<h->pre; ++i)
		h->h[i].h = yak_ct_init();
	if (n_hash > 0 && n_shift > h->pre) {
		h->n_hash = n_hash, h->n_shift = n_shift;
		for (i = 0; i < 1<<h->pre; ++i)
			h->h[i].b = yak_bf_init(h->n_shift - h->pre, h->n_hash);
	}
	return h;
}

static void ha_ct_destroy_bf(ha_ct_t *h)
{
	int i;
	for (i = 0; i < 1<<h->pre; ++i) {
		if (h->h[i].b)
			yak_bf_destroy(h->h[i].b);
		h->h[i].b = 0;
	}
}

static void ha_ct_destroy(ha_ct_t *h)
{
	int i;
	if (h == 0) return;
	ha_ct_destroy_bf(h);
	for (i = 0; i < 1<<h->pre; ++i)
		yak_ct_destroy(h->h[i].h);
	free(h->h); free(h);
}

static int ha_ct_insert_list(ha_ct_t *h, int create_new, int n, const uint64_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	ha_ct1_t *g;
	if (n == 0) return 0;
	g = &h->h[a[0]&mask];
	for (j = 0; j < n; ++j) {
		int ins = 1, absent;
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		if ((a[j]&mask) != (a[0]&mask)) continue;
		if (create_new) {
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);
			if (ins) {
				k = yak_ct_put(g->h, x << YAK_COUNTER_BITS1 | (g->b? 1 : 0), &absent);
				if (absent) ++n_ins;
				if ((kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
					++kh_key(g->h, k);
			}
		} else {
			k = yak_ct_get(g->h, x<<YAK_COUNTER_BITS1);
			if (k != kh_end(g->h) && (kh_key(g->h, k)&YAK_MAX_COUNT) < YAK_MAX_COUNT)
				++kh_key(g->h, k);
		}
	}
	return n_ins;
}

/*** generate histogram ***/

typedef struct {
	uint64_t c[YAK_N_COUNTS];
} buf_cnt_t;

typedef struct {
	const ha_ct_t *h;
	buf_cnt_t *cnt;
} hist_aux_t;

static void worker_ct_hist(void *data, long i, int tid) // callback for kt_for()
{
	hist_aux_t *a = (hist_aux_t*)data;
	uint64_t *cnt = a->cnt[tid].c;
	yak_ct_t *g = a->h->h[i].h;
	khint_t k;
	for (k = 0; k < kh_end(g); ++k)
		if (kh_exist(g, k))
			++cnt[kh_key(g, k)&YAK_MAX_COUNT];
}

static void ha_ct_hist(const ha_ct_t *h, int64_t cnt[YAK_N_COUNTS], int n_thread)
{
	hist_aux_t a;
	int i, j;
	a.h = h;
	memset(cnt, 0, YAK_N_COUNTS * sizeof(uint64_t));
	CALLOC(a.cnt, n_thread);
	kt_for(n_thread, worker_ct_hist, &a, 1<<h->pre);
	for (i = 0; i < YAK_N_COUNTS; ++i) cnt[i] = 0;
	for (j = 0; j < n_thread; ++j)
		for (i = 0; i < YAK_N_COUNTS; ++i)
			cnt[i] += a.cnt[j].c[i];
	free(a.cnt);
}

/*** shrink a hash table ***/

typedef struct {
	int min, max;
	ha_ct_t *h;
} shrink_aux_t;

typedef struct {
	int min, max;
	float diff_ratio;
	int diff_abs;
	int coverage_cutoff;  // if min(cov1, ...) is small, use diff_abs only since diff_ratio would exaggerate stuff
	ha_ct_t *h;
	All_reads *rs;  // needed to access the read (kmer-based) coverage info
} hamt_shrink_aux_t;

static void worker_ct_shrink(void *data, long i, int tid) // callback for kt_for()
{
	shrink_aux_t *a = (shrink_aux_t*)data;
	ha_ct_t *h = a->h;
	yak_ct_t *g = h->h[i].h, *f;
	khint_t k;
	f = yak_ct_init();
	yak_ct_resize(f, kh_size(g));
	for (k = 0; k < kh_end(g); ++k) {
		if (kh_exist(g, k)) {
			int absent, c = kh_key(g, k) & YAK_MAX_COUNT;
			if (c >= a->min && c <= a->max)
				yak_ct_put(f, kh_key(g, k), &absent);
		}
	}
	yak_ct_destroy(g);
	h->h[i].h = f;
}

static void ha_ct_shrink(ha_ct_t *h, int min, int max, int n_thread)
{
	int i;
	shrink_aux_t a;
	a.h = h, a.min = min, a.max = max;
	kt_for(n_thread, worker_ct_shrink, &a, 1<<h->pre);
	for (i = 0, h->tot = 0; i < 1<<h->pre; ++i)
		h->tot += kh_size(h->h[i].h);
}

// static void worker_hamt_ct_shrink(void *data, long i, int tid) // callback for kt_for()
// {
// 	hamt_shrink_aux_t *a = (hamt_shrink_aux_t*) data;
// 	ha_ct_t *h = a->h;
// 	yak_ct_t *g = h->h[i].h, *f;
// 	khint_t k;
// 	f = yak_ct_init();
// 	yak_ct_resize(f, kh_size(g));
// 	for (k = 0; k < kh_end(g); ++k) {
// 		if (kh_exist(g, k)) {
// 			int absent, c = kh_key(g, k) & YAK_MAX_COUNT;
// 			if (c >= a->min && c <= a->max)  // frequency ok
// 				if ()  // check reads - this needs a count hashtable. Does it worth it???
// 					yak_ct_put(f, kh_key(g, k), &absent);
// 		}
// 	}
// 	yak_ct_destroy(g);
// 	h->h[i].h = f;

// }

// static void hamt_ct_shrink(ha_ct_t *h, int min, int max, float diff_ratio, int diff_abs, int cutoff, All_reads *rs,int n_thread){
// 	// remove high freq kmers, AND the kmers that's share between reads of very different abundances
// 	int i;
// 	hamt_shrink_aux_t a;
// 	a.h = h, a.min = min, a.max = max;
// 	a.diff_abs = diff_abs; a.diff_ratio = diff_ratio; a.coverage_cutoff = cutoff; a.rs = rs;
// 	assert(diff_ratio>=1);
// 	assert(diff_abs>=0);
// 	assert(cutoff>=0);
// 	assert(rs);
// 	kt_for(n_thread, worker_hamt_ct_shrink, &a, 1<<h->pre);
// 	for (i = 0, h->tot = 0; i < 1<<h->pre; ++i)
// 		h->tot += kh_size(h->h[i].h);
// }

/***********************
 * Position hash table *
 ***********************/

KHASHL_MAP_INIT(static klib_unused, yak_pt_t, yak_pt, uint64_t, uint64_t, yak_ct_hash, yak_ct_eq)
#define generic_key(x) (x)
KRADIX_SORT_INIT(ha64, uint64_t, generic_key, 8)

typedef struct {
	yak_pt_t *h;
	uint64_t n;
	ha_idxpos_t *a;
} ha_pt1_t;

struct ha_pt_s {
	int k, pre;
	uint64_t tot, tot_pos;
	ha_pt1_t *h;
};

typedef struct {
	const ha_ct_t *ct;
	ha_pt_t *pt;
} pt_gen_aux_t;

static void worker_pt_gen(void *data, long i, int tid) // callback for kt_for()
{
	pt_gen_aux_t *a = (pt_gen_aux_t*)data;
	ha_pt1_t *b = &a->pt->h[i];
	yak_ct_t *g = a->ct->h[i].h;
	khint_t k;
	for (k = 0, b->n = 0; k != kh_end(g); ++k) {
		if (kh_exist(g, k)) {
			int absent;
			khint_t l;
			l = yak_pt_put(b->h, kh_key(g, k) >> YAK_COUNTER_BITS1 << YAK_COUNTER_BITS1, &absent);  // l = yak_pt_put(b->h, kh_key(g, k) >> a->ct->pre << YAK_COUNTER_BITS, &absent);
			kh_val(b->h, l) = b->n;
			b->n += kh_key(g, k) & YAK_MAX_COUNT;
		}
	}
	yak_ct_destroy(g);
	a->ct->h[i].h = 0;
	CALLOC(b->a, b->n);
}

ha_pt_t *ha_pt_gen(ha_ct_t *ct, int n_thread)
{
	pt_gen_aux_t a;
	int i;
	ha_pt_t *pt;
	ha_ct_destroy_bf(ct);
	CALLOC(pt, 1);
	pt->k = ct->k, pt->pre = ct->pre, pt->tot = ct->tot;
	CALLOC(pt->h, 1<<pt->pre);
	for (i = 0; i < 1<<pt->pre; ++i) {
		pt->h[i].h = yak_pt_init();
		yak_pt_resize(pt->h[i].h, kh_size(ct->h[i].h));
	}
	a.ct = ct, a.pt = pt;
	kt_for(n_thread, worker_pt_gen, &a, 1<<pt->pre);
	free(ct->h); free(ct);
	return pt;
}

int ha_pt_insert_list(ha_pt_t *h, int n, const ha_mz1_t *a)
{
	int j, mask = (1<<h->pre) - 1, n_ins = 0;
	ha_pt1_t *g;
	if (n == 0) return 0;
	g = &h->h[a[0].x&mask];
	for (j = 0; j < n; ++j) {
		uint64_t x = a[j].x >> h->pre;
		khint_t k;
		int n;
		ha_idxpos_t *p;
		if ((a[j].x&mask) != (a[0].x&mask)) continue;
		k = yak_pt_get(g->h, x<<YAK_COUNTER_BITS1);
		if (k == kh_end(g->h)) continue;
		n = kh_key(g->h, k) & YAK_MAX_COUNT;
		assert(n < YAK_MAX_COUNT);
		p = &g->a[kh_val(g->h, k) + n];
		p->rid = a[j].rid, p->rev = a[j].rev, p->pos = a[j].pos, p->span = a[j].span;
		//(uint64_t)a[j].rid<<36 | (uint64_t)a[j].rev<<35 | (uint64_t)a[j].pos<<8 | (uint64_t)a[j].span;
		++kh_key(g->h, k);
		++n_ins;
	}
	return n_ins;
}
/*
static void worker_pt_sort(void *data, long i, int tid)
{
	ha_pt_t *h = (ha_pt_t*)data;
	ha_pt1_t *g = &h->h[i];
	khint_t k;
	for (k = 0; k < kh_end(g->h); ++k) {
		int n;
		uint64_t *p;
		if (!kh_exist(g->h, k)) continue;
		n = kh_key(g->h, k) & YAK_MAX_COUNT;
		p = &g->a[kh_val(g->h, k)];
		radix_sort_ha64(p, p + n);
	}
}

void ha_pt_sort(ha_pt_t *h, int n_thread)
{
	kt_for(n_thread, worker_pt_sort, h, 1<<h->pre);
}
*/
void ha_pt_destroy(ha_pt_t *h)
{
	int i;
	if (h == 0) return;
	for (i = 0; i < 1<<h->pre; ++i) {
		yak_pt_destroy(h->h[i].h);
		free(h->h[i].a);
	}
	free(h->h); free(h);
}

const ha_idxpos_t *ha_pt_get(const ha_pt_t *h, uint64_t hash, int *n)
{
	khint_t k;
	const ha_pt1_t *g = &h->h[hash & ((1ULL<<h->pre) - 1)];
	*n = 0;
	k = yak_pt_get(g->h, hash >> h->pre << YAK_COUNTER_BITS1);
	if (k == kh_end(g->h)) return 0;
	*n = kh_key(g->h, k) & YAK_MAX_COUNT;
	return &g->a[kh_val(g->h, k)];
}

/**********************************
 * Buffer for counting all k-mers *
 **********************************/

typedef struct {
	int n, m;
	uint64_t n_ins;
	uint64_t *a;
	ha_mz1_t *b;
} ch_buf_t;

static inline void ct_insert_buf(ch_buf_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
	int pre = y & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->a, b->m);
	}
	b->a[b->n++] = y;
}

static inline void pt_insert_buf(ch_buf_t *buf, int p, const ha_mz1_t *y)
{
	int pre = y->x & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->b, b->m);
	}
	b->b[b->n++] = *y;
}

static void count_seq_buf(ch_buf_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 1 | (c&1))  & mask;
			x[1] = (x[1] << 1 | (c>>1)) & mask;
			x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
			x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
			if (++l >= k)
				ct_insert_buf(buf, p, yak_hash_long(x));
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
}

static void count_seq_buf_HPC(ch_buf_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
{
	int i, l, last = -1;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			if (c != last) {
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
				if (++l >= k)
					ct_insert_buf(buf, p, yak_hash_long(x));
				last = c;
			}
		} else l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
}

/******************
 * K-mer counting *
 ******************/

KSEQ_INIT(gzFile, gzread)

#define HAF_COUNT_EXACT  0x1
#define HAF_COUNT_ALL    0x2
#define HAF_RS_WRITE_LEN 0x4
#define HAF_RS_WRITE_SEQ 0x8
#define HAF_RS_READ      0x10
#define HAF_CREATE_NEW   0x20
#define HAMTF_FORCE_DONT_INIT 0x40  // meta: override HAF_RS_WRITE_LEN to not allow init_all_reads of ha_count

typedef struct { // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	const void *flt_tab;
	int flag, create_new, is_store;
	uint64_t n_seq;
	kseq_t *ks;
	UC_Read ucr;
	ha_ct_t *ct;
	ha_pt_t *pt;
	const All_reads *rs_in;
	All_reads *rs_out;
	All_reads *rs;  // hamt meta. ALWAYS carry this
} pl_data_t;

typedef struct { // data structure for each step in kt_pipeline()
	pl_data_t *p;
	uint64_t n_seq0;
	int n_seq, m_seq, sum_len, nk;
	int *len;
	char **seq;
	ha_mz1_v *mz_buf;
	ha_mz1_v *mz;
	ch_buf_t *buf;
} st_data_t;

static void worker_for_insert(void *data, long i, int tid) // callback for kt_for()
{
	st_data_t *s = (st_data_t*)data;
	ch_buf_t *b = &s->buf[i];
	if (s->p->pt)
		b->n_ins += ha_pt_insert_list(s->p->pt, b->n, b->b);
	else
		b->n_ins += ha_ct_insert_list(s->p->ct, s->p->create_new, b->n, b->a);
}

static void worker_for_mz(void *data, long i, int tid)
{
	st_data_t *s = (st_data_t*)data;
	ha_mz1_v *b = &s->mz_buf[tid];
	s->mz_buf[tid].n = 0;
	ha_sketch(s->seq[i], s->len[i], s->p->opt->w, s->p->opt->k, s->n_seq0 + i, s->p->opt->is_HPC, b, s->p->flt_tab);
	s->mz[i].n = s->mz[i].m = b->n;
	MALLOC(s->mz[i].a, b->n);
	memcpy(s->mz[i].a, b->a, b->n * sizeof(ha_mz1_t));
}

static void *worker_count(void *data, int step, void *in) // callback for kt_pipeline()
{
	pl_data_t *p = (pl_data_t*)data;
	if (step == 0) { // step 1: read a block of sequences
		// printf("worker_count step0\n"); fflush(stdout);
		int ret;
		st_data_t *s;
		CALLOC(s, 1);
		s->p = p;
		s->n_seq0 = p->n_seq;
		if (p->rs_in && (p->flag & HAF_RS_READ)) {
			while (p->n_seq < p->rs_in->total_reads) {
				int l;
				recover_UC_Read(&p->ucr, p->rs_in, p->n_seq);
				l = p->ucr.length;
				if (s->n_seq == s->m_seq) {
					s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
					REALLOC(s->len, s->m_seq);
					REALLOC(s->seq, s->m_seq);
				}
				MALLOC(s->seq[s->n_seq], l);
				memcpy(s->seq[s->n_seq], p->ucr.seq, l);
				s->len[s->n_seq++] = l;
				++p->n_seq;
				s->sum_len += l;
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
				if (s->sum_len >= p->opt->chunk_size)
					break;
			}
		} else {
			// printf("worker_count step0, will read from file\n"); fflush(stdout);
			while ((ret = kseq_read(p->ks)) >= 0) {
				int l = p->ks->seq.l;
				if (p->n_seq >= 1<<28) {
					fprintf(stderr, "ERROR: this implementation supports no more than %d reads\n", 1<<28);
					exit(1);
				}
				if (p->rs_out) {
					if (p->flag & HAF_RS_WRITE_LEN) {
						assert(p->n_seq == p->rs_out->total_reads);
						ha_insert_read_len(p->rs_out, l, p->ks->name.l);
					} else if (p->flag & HAF_RS_WRITE_SEQ) {
						// printf("woker_count loading seq sancheck: name %s, name length %zu; total name length %" PRIu64 "\n", p->ks->name.s, p->ks->name.l, p->rs_out->total_name_length);
						// printf("n_seq: %" PRIu64 ", name index: %" PRIu64 "\n", p->n_seq, p->rs_out->name_index[p->n_seq]); fflush(stdout);
						int i, n_N;
						assert(l == (int)p->rs_out->read_length[p->n_seq]);
						for (i = n_N = 0; i < l; ++i) // count number of ambiguous bases
							if (seq_nt4_table[(uint8_t)p->ks->seq.s[i]] >= 4)
								++n_N;
						ha_compress_base(Get_READ(*p->rs_out, p->n_seq), p->ks->seq.s, l, &p->rs_out->N_site[p->n_seq], n_N);
						memcpy(&p->rs_out->name[p->rs_out->name_index[p->n_seq]], p->ks->name.s, p->ks->name.l);
					}
				}
				if (s->n_seq == s->m_seq) {
					s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
					REALLOC(s->len, s->m_seq);
					REALLOC(s->seq, s->m_seq);
				}
				MALLOC(s->seq[s->n_seq], l);
				memcpy(s->seq[s->n_seq], p->ks->seq.s, l);
				s->len[s->n_seq++] = l;
				++p->n_seq;
				s->sum_len += l;
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
				if (s->sum_len >= p->opt->chunk_size)
					break;
			}
			// printf("worker_count step0 finished reading from file\n"); fflush(stdout);
		}
		// printf("exit step0\n"); fflush(stdout);
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		st_data_t *s = (st_data_t*)in;
		// printf("worker_count step1, n_seq0=%" PRIu64 ", n_seq=%d, buf size=%" PRIu64 "\n", s->n_seq0, s->n_seq, p->rs->hamt_stat_buf_size); fflush(stdout);
		int i, n_pre = 1<<p->opt->pre, m;
		// allocate the k-mer buffer
		CALLOC(s->buf, n_pre);
		m = (int)(s->nk * 1.2 / n_pre) + 1;
		for (i = 0; i < n_pre; ++i) {
			s->buf[i].m = m;
			if (p->pt) MALLOC(s->buf[i].b, m);
			else MALLOC(s->buf[i].a, m);
		}
		// printf("default buf size:%d", READ_INIT_NUMBER);fflush(stdout);
		// printf("going to fill linear buffer, current hamt buf size%" PRIu64 "\n", s->p->rs_in->index_size); fflush(stdout);
		// fill the buffer
		if (p->opt->w == 1) { // enumerate all k-mers
			for (i = 0; i < s->n_seq; ++i) {
				//////////////////  meta   ////////////////////
				if (s->p->flag & HAMTF_FORCE_DONT_INIT){
					uint64_t rid = s->n_seq0+i;
					assert(rid<p->rs->hamt_stat_buf_size);
					if (p->rs->mask_readnorm[rid]>0) // initial marking hasn't finished; read is marked as dropped
						{continue;}
				}
				//////////////////////////////////////////////
				if (p->opt->is_HPC)
					count_seq_buf_HPC(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i]);
				else
					count_seq_buf(s->buf, p->opt->k, p->opt->pre, s->len[i], s->seq[i]);
				if (!p->is_store) free(s->seq[i]);
			}
		} else { // minimizers only
			uint32_t j;
			// compute minimizers
			CALLOC(s->mz, s->n_seq);
			CALLOC(s->mz_buf, p->opt->n_thread);
			kt_for(p->opt->n_thread, worker_for_mz, s, s->n_seq);
			for (i = 0; i < p->opt->n_thread; ++i)
				free(s->mz_buf[i].a);
			free(s->mz_buf);
			// insert minimizers
			if (p->pt) {
				for (i = 0; i < s->n_seq; ++i){
					//////////////////  meta   ////////////////////
					if (p->flag & HAMTF_FORCE_DONT_INIT){
						uint64_t rid = s->n_seq0+i;
						// printf("assert!\n"); fflush(stdout);
						assert(rid<p->rs->hamt_stat_buf_size);
						if (s->p->rs->mask_readnorm[rid]>0) // initial marking hasn't finished; read is marked as dropped
							{continue;}
					}
					//////////////////////////////////////////////
					for (j = 0; j < s->mz[i].n; ++j)
						pt_insert_buf(s->buf, p->opt->pre, &s->mz[i].a[j]);
				}
			} else {
				for (i = 0; i < s->n_seq; ++i){
					//////////////////  meta   ////////////////////
					if (p->flag & HAMTF_FORCE_DONT_INIT){
						uint64_t rid = s->n_seq0+i;
						if (s->p->rs->mask_readnorm[rid]>0) // initial marking hasn't finished; read is marked as dropped
							{continue;}
					}
					//////////////////////////////////////////////
					for (j = 0; j < s->mz[i].n; ++j)
						ct_insert_buf(s->buf, p->opt->pre, s->mz[i].a[j].x);
				}
			}
			for (i = 0; i < s->n_seq; ++i) {
				free(s->mz[i].a);
				if (!p->is_store) free(s->seq[i]);
			}
			free(s->mz);
		}
		// printf("exit step1\n"); fflush(stdout);
		free(s->seq); free(s->len);
		s->seq = 0, s->len = 0;
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		// printf("worker_count step2\n"); fflush(stdout);
		st_data_t *s = (st_data_t*)in;
		int i, n = 1<<p->opt->pre;
		uint64_t n_ins = 0;
		kt_for(p->opt->n_thread, worker_for_insert, s, n);
		for (i = 0; i < n; ++i) {
			n_ins += s->buf[i].n_ins;
			if (p->pt) free(s->buf[i].b);
			else free(s->buf[i].a);
		}
		if (p->ct) p->ct->tot += n_ins;
		if (p->pt) p->pt->tot_pos += n_ins;
		free(s->buf);
		#if 0
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %ld sequences; %ld %s in the hash table\n", __func__,
				yak_realtime(), yak_cpu_usage(), (long)s->n_seq0 + s->n_seq,
				(long)(p->pt? p->pt->tot_pos : p->ct->tot), p->pt? "positions" : "distinct k-mers");
		#endif
		free(s);
	}
	return 0;
}

static ha_ct_t *yak_count(const yak_copt_t *opt, const char *fn, int flag, ha_pt_t *p0, ha_ct_t *c0, const void *flt_tab, All_reads *rs, int64_t *n_seq)
{
	printf("entered yak_count\n"); fflush(stdout);
	int read_rs = (rs && (flag & HAF_RS_READ));
	pl_data_t pl;
	gzFile fp = 0;
	memset(&pl, 0, sizeof(pl_data_t));
	pl.n_seq = *n_seq;
	if (read_rs) {
		pl.rs_in = rs;
		init_UC_Read(&pl.ucr);
	} else {
		if ((fp = gzopen(fn, "r")) == 0) return 0;
		pl.ks = kseq_init(fp);
	}
	if (rs && (flag & (HAF_RS_WRITE_LEN|HAF_RS_WRITE_SEQ)))
		pl.rs_out = rs;
	pl.rs = rs;
	pl.flt_tab = flt_tab;
	pl.opt = opt;
	pl.flag = flag;
	if (p0) {
		pl.pt = p0, pl.create_new = 0; // never create new elements in a position table
		assert(p0->k == opt->k && p0->pre == opt->pre);
	} else if (c0) {
		pl.ct = c0, pl.create_new = !!(flag&HAF_CREATE_NEW);
		assert(c0->k == opt->k && c0->pre == opt->pre);
	} else {
		pl.create_new = 1; // alware create new elements if the count table is empty
		pl.ct = ha_ct_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	}
	kt_pipeline(3, worker_count, &pl, 3);
	if (read_rs) {
		destory_UC_Read(&pl.ucr);
	} else {
		kseq_destroy(pl.ks);
		gzclose(fp);
	}
	*n_seq = pl.n_seq;
	return pl.ct;
}

ha_ct_t *ha_count(const hifiasm_opt_t *asm_opt, int flag, ha_pt_t *p0, const void *flt_tab, All_reads *rs)
{
	int i;
	int64_t n_seq = 0;
	yak_copt_t opt;
	ha_ct_t *h = 0;
	assert(!(flag & HAF_RS_WRITE_LEN) || !(flag & HAF_RS_WRITE_SEQ)); // not both
	if (rs) {
		if (flag & HAF_RS_WRITE_LEN){
		    if (!(flag & HAMTF_FORCE_DONT_INIT)){
				init_All_reads(rs);
			}
			else{
				printf("reset because flag\n"); fflush(stdout);
				reset_All_reads(rs);
			}
		}
		else if (flag & HAF_RS_WRITE_SEQ){
			printf("ha_count: malloc_all_reads\n");
			malloc_All_reads(rs);
		}
	}
	yak_copt_init(&opt);
	opt.k = asm_opt->k_mer_length;
	opt.is_HPC = !(asm_opt->flag&HA_F_NO_HPC);
	opt.w = flag & HAF_COUNT_ALL? 1 : asm_opt->mz_win;
	opt.bf_shift = flag & HAF_COUNT_EXACT? 0 : asm_opt->bf_shift;
	opt.n_thread = asm_opt->thread_num;
	for (i = 0; i < asm_opt->num_reads; ++i)
		h = yak_count(&opt, asm_opt->read_file_names[i], flag|HAF_CREATE_NEW, p0, h, flt_tab, rs, &n_seq);
	if (h && opt.bf_shift > 0)
		ha_ct_destroy_bf(h);
	return h;
}

/***************************
 * High count filter table *
 ***************************/

KHASHL_SET_INIT(static klib_unused, yak_ft_t, yak_ft, uint64_t, kh_hash_dummy, kh_eq_generic)

static yak_ft_t *gen_hh(const ha_ct_t *h)
{
	int i;
	yak_ft_t *hh;
	hh = yak_ft_init();
	yak_ft_resize(hh, h->tot * 2);
	for (i = 0; i < 1<<h->pre; ++i) {
		yak_ct_t *ht = h->h[i].h;
		khint_t k;
		for (k = 0; k < kh_end(ht); ++k) {
			if (kh_exist(ht, k)) {
				uint64_t y = kh_key(ht, k) >> h->pre << YAK_COUNTER_BITS1 | i;
				int absent;
				yak_ft_put(hh, y, &absent);
			}
		}
	}
	return hh;
}

int ha_ft_isflt(const void *hh, uint64_t y)
{
	yak_ft_t *h = (yak_ft_t*)hh;
	khint_t k;
	k = yak_ft_get(h, y);
	return k == kh_end(h)? 0 : 1;
}

void ha_ft_destroy(void *h)
{
	if (h) yak_ft_destroy((yak_ft_t*)h);
}

/*************************
 * High-level interfaces *
 *************************/

void *ha_ft_gen(const hifiasm_opt_t *asm_opt, All_reads *rs, int *hom_cov)
{
	yak_ft_t *flt_tab;
	int64_t cnt[YAK_N_COUNTS];
	int peak_hom, peak_het, cutoff;
	ha_ct_t *h;
	h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_WRITE_LEN, NULL, NULL, rs);
	ha_ct_hist(h, cnt, asm_opt->thread_num);
	peak_hom = ha_analyze_count(YAK_N_COUNTS, cnt, &peak_het);
	if (hom_cov) *hom_cov = peak_hom;
	if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
	cutoff = (int)(peak_hom * asm_opt->high_factor);
	if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
	ha_ct_shrink(h, cutoff, YAK_MAX_COUNT, asm_opt->thread_num);
	flt_tab = gen_hh(h);
	ha_ct_destroy(h);
	fprintf(stderr, "[M::%s::%.3f*%.2f@%.3fGB] ==> filtered out %ld k-mers occurring %d or more times\n", __func__,
			yak_realtime(), yak_cpu_usage(), yak_peakrss_in_gb(), (long)kh_size(flt_tab), cutoff);
	return (void*)flt_tab;
}

ha_pt_t *ha_pt_gen(const hifiasm_opt_t *asm_opt, const void *flt_tab, int read_from_store, All_reads *rs, int *hom_cov, int *het_cov)
{
	int64_t cnt[YAK_N_COUNTS], tot_cnt;
	int peak_hom, peak_het, i, extra_flag1, extra_flag2;
	ha_ct_t *ct;
	ha_pt_t *pt;
	if (read_from_store) {
		// printf("ha_pt_gen status: all in mem\n"); fflush(stdout);
		extra_flag1 = extra_flag2 = HAF_RS_READ;
	} else if (rs->total_reads == 0) {
		// printf("ha_pt_gen status: none in mem\n"); fflush(stdout);
		extra_flag1 = HAF_RS_WRITE_LEN;
		extra_flag2 = HAF_RS_WRITE_SEQ;
	} else {
		// printf("ha_pt_gen status: knew length, will load seqs in mem\n"); fflush(stdout);
		extra_flag1 = HAF_RS_WRITE_SEQ;
		extra_flag2 = HAF_RS_READ;
	}
	ct = ha_count(asm_opt, HAF_COUNT_EXACT|extra_flag1|HAMTF_FORCE_DONT_INIT, NULL, flt_tab, rs);  // collect minimizer (aware of high freq filter)
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> counted %ld distinct minimizer k-mers\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)ct->tot);
	ha_ct_hist(ct, cnt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s] count[%d] = %ld (for sanity check)\n", __func__, YAK_MAX_COUNT, (long)cnt[YAK_MAX_COUNT]);
	peak_hom = ha_analyze_count(YAK_N_COUNTS, cnt, &peak_het);
	// if (hom_cov) *hom_cov = peak_hom;
	// if (het_cov) *het_cov = peak_het;
	if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
	if (flt_tab == 0) {
		int cutoff = (int)(peak_hom * asm_opt->high_factor);
		if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
		ha_ct_shrink(ct, 2, cutoff, asm_opt->thread_num);
		for (i = 2, tot_cnt = 0; i <= cutoff; ++i) tot_cnt += cnt[i] * i;
	} else {
		ha_ct_shrink(ct, 2, YAK_MAX_COUNT - 1, asm_opt->thread_num);
		for (i = 2, tot_cnt = 0; i <= YAK_MAX_COUNT - 1; ++i) tot_cnt += cnt[i] * i;
	}
	pt = ha_pt_gen(ct, asm_opt->thread_num); // prepare the slots
	ha_count(asm_opt, HAF_COUNT_EXACT|extra_flag2|HAMTF_FORCE_DONT_INIT, pt, flt_tab, rs);  // collect positional info
	fprintf(stderr, "tot_cnt=%" PRIu64 "\n", tot_cnt);
	fprintf(stderr, "tot_pos=%" PRIu64 "\n", pt->tot_pos);
	assert((uint64_t)tot_cnt == pt->tot_pos);
	//ha_pt_sort(pt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> indexed %ld positions\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)pt->tot_pos);
	return pt;
}


//////////////////////////////////////////////////////////////////////
//                      meta
//////////////////////////////////////////////////////////////////////
KRADIX_SORT_INIT(hamt64, uint64_t, uint64_t, 8)  // sort read's kmer count profile
KRADIX_SORT_INIT(hamt16, uint16_t, uint16_t, 2)  // sort read's kmer count profile

typedef struct {  // linear buffer, by reads
	uint64_t n, m;  // note to self: somehow using a[0] and a[1] for this purpose could fail (segfault, realloc is successful but something goes wrong after that?). It shouldn't but I didn't know why. 
	uint64_t *a;  // a linear buffer of hash values
} hamt_ch_buf_t;

typedef struct {  // global data structure for pipeline (pl_data_t)
	yak_copt_t *opt;
	uint64_t *n_seq;
	All_reads *rs_in;
	All_reads *rs_out;
	UC_Read ucr;
	kseq_t *ks;
	ha_ct_t *h;  // counted kmer tables (equals to the h from 1st call of ha_count in vanilla hifiasm)
	uint16_t *buf;  // for round 0; uint16 because yak_ct_t use 12 bits for counting
	// int flag, create_new, is_store;
	int round;
	ha_ct_t *hd;// "hashtables for read drop": the runtime kmer counts
	// uint64_t *san_n_insert;// sancheck for migrating h to hd

	uint8_t flag;
}plmt_data_t;

typedef struct {  // step data structure (st_data_t)
	plmt_data_t* p;
	uint64_t n_seq0;
	int n_seq, m_seq, sum_len, nk;
	int *len;
	char **seq;
	hamt_ch_buf_t *lnbuf;  // $n_seq linear buffers
}pltmt_step_t;

static void worker_insert_lnbuf(void* data, long idx_for, int tid){
	// note to self (bug or?): race of inserting kmers from different thread will segfault for some reason. (i think)
	pltmt_step_t *s = (pltmt_step_t*) data;
	uint64_t *buf = s->lnbuf[idx_for].a;
	uint64_t n = s->lnbuf[idx_for].n;
	yak_ct_t *h = s->p->hd->h[idx_for].h;
	int absent;
	khint_t key;
	uint64_t san = 0;
	for (uint32_t i=0; i<n; i++){
		key = yak_ct_put(h, buf[i]>>s->p->h->pre<<YAK_COUNTER_BITS1, &absent);
		kh_key(h, key)++;
		if (absent) san++;
	}
	// printf("%s::%ld: inserted %" PRIu64 " kmers.\n", __func__, idx_for, san);
}

static void worker_process_one_read(void* data, long idx_for, int tid){
	// !!!!!! LAST UPDATE AUG 20. Based on HPC version but not tested, since non-HPC isn't enabled anywhere !!!!!! //
	pltmt_step_t *s = (pltmt_step_t*) data;
	uint64_t rid = s->n_seq0+idx_for;
	uint16_t *buf = (uint16_t*)malloc(sizeof(uint16_t)*s->len[idx_for]);
	if (!buf) {printf("[error::%s] malloc for buffer failed, thread %d, for_idx %ld.\n", __func__, tid, idx_for); exit(1);}
	int k = s->p->opt->k;
	ha_ct_t *h = s->p->h;

	uint32_t idx = 0;  // index of kmer count in the buffer
	khint_t key;
	int i, l;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->len[idx_for]; ++i) {
		int c = seq_nt4_table[(uint8_t)s->seq[idx_for][i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 1 | (c&1))  & mask;
			x[1] = (x[1] << 1 | (c>>1)) & mask;
			x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
			x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
			if (++l >= k){
					uint64_t hash = yak_hash_long(x);
					yak_ct_t *h_ = h->h[hash & ((1<<h->pre)-1)].h;
					key = yak_ct_get(h_, hash>>h->pre<<YAK_COUNTER_BITS1);
					if (key!=kh_end(h_)){buf[idx] = (kh_key(h_, key) & YAK_MAX_COUNT);}
					else buf[idx] = 0;
					idx++;
				}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
	double mean = meanl(buf, idx);
	double std = stdl(buf, idx, mean);
	radix_sort_hamt16(buf, buf+idx);
	uint16_t median = (uint16_t)buf[idx/2];
	uint8_t code = decide_category(mean, std, buf, idx);
	s->p->rs_out->mean[rid] = mean;
	s->p->rs_out->median[rid] = median;
	s->p->rs_out->std[rid] = std;
	s->p->rs_out->mask_readtype[rid] = code;
	// printf("~%" PRIu64 "\t%f\t%f\t%d\n", rid, mean, std, code);
	free(buf);
	// sequence will be freed by the caller
}

static void worker_process_one_read_HPC(void* data, long idx_for, int tid){
	pltmt_step_t *s = (pltmt_step_t*) data;
	// if (s->p->round) printf("debug\tenter %s (round %d)\n", __func__, s->p->round);fflush(stdout);
	uint64_t rid = s->n_seq0+idx_for;
	uint16_t *buf = (uint16_t*)malloc(sizeof(uint16_t)*s->len[idx_for]);  // bufer used in the 0th round
	uint16_t *buf_norm = (uint16_t*)malloc(sizeof(uint16_t)*s->len[idx_for]);  // bufer used in the 1st+ round
	if (!buf) {printf("[error::%s] malloc for buffer failed, thread %d, for_idx %ld.\n", __func__, tid, idx_for); exit(1);}
	int k = s->p->opt->k;
	ha_ct_t *h = s->p->h;

	uint32_t idx = 0;  // index of kmer count in the buffer
	khint_t key;
	int i, l, last = -1;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	hamt_ch_buf_t *b = 0;// uint64_t *b = 0;

	int has_wanted_kmers = 0;  // flag of rescue

	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->len[idx_for]; ++i) {
		int c = seq_nt4_table[(uint8_t)s->seq[idx_for][i]];
		if (c < 4) { // not an "N" base
			if (c != last) {
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
				if (++l >= k){
					uint64_t hash = yak_hash_long(x);
					if (s->p->round==0){
						yak_ct_t *h_ = h->h[hash & ((1<<h->pre)-1)].h;
						key = yak_ct_get(h_, hash>>h->pre<<YAK_COUNTER_BITS1);
						if (key!=kh_end(h_)){buf[idx] = (kh_key(h_, key) & YAK_MAX_COUNT);}
						else buf[idx] = 0;
						idx++;
					} else{  
						// inert into diginorm linear buffer for kmer counting
						b = &s->lnbuf[hash & ((1<<h->pre)-1)];
						b->a[b->n] = hash;
						// collect runtime count
						yak_ct_t *h_ = s->p->hd->h[hash & ((1<<h->pre)-1)].h;
						key = yak_ct_get(h_, hash>>h->pre<<YAK_COUNTER_BITS1);
						if (key!=kh_end(h_)){buf_norm[idx] = (kh_key(h_, key) & YAK_MAX_COUNT);}
						else buf_norm[idx] = 0;
						idx++;
						b->n++;
						if (b->n==b->m){
							b->m = b->m+((b->m)>>1);
							b->a = (uint64_t*)realloc(b->a, sizeof(uint64_t)*b->m);
						}
						// check whether it's a relatively low-freq (*overall freq) kmer
						yak_ct_t *h2 = s->p->h->h[hash & ((1<<h->pre)-1)].h;
						key = yak_ct_get(h2, hash>>h->pre<<YAK_COUNTER_BITS1);
						if (key!=kh_end(h2)){
							if ((kh_key(h2, key) & YAK_MAX_COUNT) >= HAMT_DIG_KMERRESCUE) has_wanted_kmers++;
						}
					}
				}
				last = c;
			}
		} else l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
	if (s->p->round==0){  // the initial marking
		double mean = meanl(buf, idx);
		double std = stdl(buf, idx, mean);
		radix_sort_hamt16(buf, buf+idx);
		uint16_t median = (uint16_t)buf[idx/2];
		uint8_t code = decide_category(mean, std, buf, idx);
		s->p->rs_out->mean[rid] = mean;
		s->p->rs_out->median[rid] = median;
		s->p->rs_out->std[rid] = std;
		s->p->rs_out->mask_readtype[rid] = code;
		// printf("~%" PRIu64 "\t%f\t%f\t%d\n", rid, mean, std, code);
	}else{  // diginorm
		int code;
		uint16_t cnt = 0;
		uint16_t max_cnt = 0;
		uint16_t median = 0;
		for (uint16_t i=0; i<idx; i++){  // examine if runtime buffer looks like long low-coverage. 
			if (buf_norm[i]<=10)
				cnt+=1;
			else if (buf_norm[i]>30){
				if (cnt>max_cnt) max_cnt = cnt;
				cnt = 0;
			}
		}
		if (((double)max_cnt/idx)>0.3 || (has_wanted_kmers>=3000)) {  // If so, or there's relative low (overall) freq kmers, ignore other hints and keep the read.
		// if (((double)max_cnt/idx)>0.3) {
			code = 0;
		}else{  // get median and decide
			radix_sort_hamt16(buf_norm, buf_norm+idx);
			median = (uint16_t)buf_norm[idx/2];
			code = decide_drop(s->p->rs_out->mean[rid], s->p->rs_out->std[rid], median, s->p->round, s->p->rs_out->mask_readtype[rid]);
		}
		s->p->rs_out->mask_readnorm[rid] = (uint8_t)code;  // todo: bit flag and add debug info?
		// printf("report\truntime median: %" PRIu16 ", coding: %d\n", median, code);
	}
	free(buf);  // sequence will be freed by the caller	
	free(buf_norm);
	// if (s->p->round) printf("debug\tproperly exit %s (round %d)\n", __func__, s->p->round);fflush(stdout);
}

static void *worker_mark_reads(void *data, int step, void *in){  // callback of kt_pipeline
	plmt_data_t *p = (plmt_data_t*)data;
	if (step==0){ // step 1: read a block of sequences
		// printf("debug\tenter step %d, round %d, total seqs= %" PRIu64 "\n", step, p->round, *p->n_seq); fflush(stdout);
		int ret;
		pltmt_step_t *s;
		CALLOC(s, 1);
		s->p = p;
		s->n_seq0 = *p->n_seq;
		s->lnbuf = 0;
		if (p->round!=0) {  // init kmer counting buffers
			s->lnbuf = (hamt_ch_buf_t*)calloc(1<<p->h->pre, sizeof(hamt_ch_buf_t));
			for (int i=0; i<1<<p->h->pre; i++){
				s->lnbuf[i].a = (uint64_t*)calloc(15000, sizeof(uint64_t));
				s->lnbuf[i].n = 0;
				s->lnbuf[i].m = 15000;
			}
		}
		while ((ret = kseq_read(p->ks)) >= 0) {
			int l = p->ks->seq.l;
			if (*p->n_seq >= 1<<28) {
				fprintf(stderr, "ERROR: this implementation supports no more than %d reads\n", 1<<28);
				exit(1);
			}
			if (p->rs_out) {
				assert(*p->n_seq == p->rs_out->total_reads);
				ha_insert_read_len(p->rs_out, l, p->ks->name.l);
			}
			if (s->n_seq == s->m_seq) {
				s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
				REALLOC(s->len, s->m_seq);
				REALLOC(s->seq, s->m_seq);
			}
			MALLOC(s->seq[s->n_seq], l);
			memcpy(s->seq[s->n_seq], p->ks->seq.s, l);
			s->len[s->n_seq++] = l;
			*p->n_seq+=1;
			s->sum_len += l;
			s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
			if (s->n_seq>=2000 || s->sum_len >= p->opt->chunk_size){  // 2000 is because of the linear buffers for diginorm steps
				break;
			}
		}
		// printf("exist 0\n"); fflush(stdout);
		if (s->sum_len == 0) {
			// printf("debug\tterminating step0, round=%d, n_seq=%" PRIu64 "\n", p->round, *p->n_seq);fflush(stdout);
			free(s); 
		}
		else return s;

	}else if (step==1){  // step2: (round 0:)get kmer profile of each read && mark their status (mask_readtype and mask_readnorm) || (round 1:) count kmers into linear buffers, then insert into hashtables

		pltmt_step_t *s = (pltmt_step_t*)in;
		// printf("debug\tstep %d, round %d, n_seq0 %" PRIu64 ", n_seq=%d\n", step, p->round, s->n_seq0, s->n_seq); fflush(stdout);
		if (p->opt->is_HPC)
			kt_for(s->p->opt->n_thread, worker_process_one_read_HPC, s, s->n_seq);
		else
			kt_for(s->p->opt->n_thread, worker_process_one_read, s, s->n_seq);
		if (s->p->round==0){
			// printf("debug\tterminating step1, round=%d, n_seq=%" PRIu64 "\n", p->round, *p->n_seq);fflush(stdout);
			for (int i=0; i<s->n_seq; i++){free(s->seq[i]);}
			free(s->len);
			free(s->seq);
			s->seq = 0; 
			s->len = 0;
			free(s);
		}else{  // insert the linear buffers
			assert(s->p->hd);
			kt_for(s->p->opt->n_thread, worker_insert_lnbuf, s, 1<<p->h->pre);
			// clean up
			for (int i=0; i<s->n_seq; i++){free(s->seq[i]);}
			free(s->len);
			free(s->seq);
			s->seq = 0; 
			s->len = 0;
			// printf("CLEANED SEQS\n"); fflush(stdout);
			for (int i=0; i<1<<p->h->pre; i++){
				free(s->lnbuf[i].a);
				}				
			free(s->lnbuf);
			free(s);
		}
	}
	return 0;
}

// static void gen_hd(void* data, long i, int tid){
// 	plmt_data_t *plmt = (plmt_data_t*)data;
// 	yak_ct_t *g = plmt->hd->h[i].h;
// 	yak_ct_t *H = plmt->h->h[i].h;
// 	khint_t key;
// 	int absent;
// 	for (key=0; key<kh_end(H); key++){
// 		if (kh_exist(H, key)){
// 			uint64_t hash = key>>(plmt->h->pre)<<YAK_COUNTER_BITS1;
// 			yak_ct_put(g, hash, &absent);
// 			if (absent) plmt->san_n_insert[i]++;
// 		}
// 	}
// }

void hamt_mark(const hifiasm_opt_t *asm_opt, All_reads *rs, ha_ct_t *h, int round){
	reset_All_reads(rs); 
	// init data structures
	int i;
	uint64_t n_seq = 0;
	yak_copt_t opt;
	yak_copt_init(&opt);
	opt.k = asm_opt->k_mer_length;
	opt.is_HPC = !(asm_opt->flag&HA_F_NO_HPC);
	opt.n_thread = asm_opt->thread_num;
	plmt_data_t plmt;
	// go
	for (i = 0; i < asm_opt->num_reads; ++i){  // note that num_reads is the number of input files, not reads
		memset(&plmt, 0, sizeof(plmt_data_t));
		plmt.n_seq = &n_seq;
		plmt.opt = &opt;
		plmt.h = h;
		plmt.rs_in = rs;
		plmt.rs_out = rs;
		plmt.round = round;
		plmt.hd = 0;
		if (plmt.round!=0){
			plmt.hd = ha_ct_init(plmt.opt->k, plmt.h->pre, plmt.opt->bf_n_hash, plmt.opt->bf_shift);
			assert(h);
			// plmt.san_n_insert = (uint64_t*)calloc((1<<h->pre), sizeof(uint64_t));
			// kt_for(plmt.opt->n_thread, gen_hd, &plmt, (1<<h->pre));  //gen_hh-like thing
			// uint64_t san = 0;
			// for (int i=0; i<((1<<h->pre)); i++) san += plmt.san_n_insert[i];
			// printf("debug\tmigrated %" PRIu64 "kmers into hd\n", san); fflush(stdout);
		}
		gzFile fp = 0;
		if ((fp = gzopen(asm_opt->read_file_names[i], "r")) == 0) {
			printf("[error::%s] can't open file %s. Abort.\n", __func__, asm_opt->read_file_names[i]);
			exit(1);
		}
		plmt.ks = kseq_init(fp);
		init_UC_Read(&plmt.ucr);
		kt_pipeline(2, worker_mark_reads, &plmt, 2);

		destory_UC_Read(&plmt.ucr);
		if (plmt.hd) ha_ct_destroy(plmt.hd);
		kseq_destroy(plmt.ks);
		gzclose(fp);
	}
	// //debug print
	// for (uint64_t i=0;i<n_seq; i++){
	// 	printf("%" PRIu64 "\t", i);
	// 	printf("%f\t", rs->mean[i]);
	// 	printf("\t%" PRIu8 "\n", rs->mask_readtype[i]);
	// 	fflush(stdout);
	// }
	// //end of debug print

	
}

void hamt_flt(const hifiasm_opt_t *asm_opt, All_reads *rs, int cov, int is_crude){
    // is_crude: be (not) aware to freq distribution on the read; might override cov (if the overall coverage is too low)
	// collect kmers frequencies
	int64_t cnt[YAK_N_COUNTS];
	ha_ct_t *h;
	int peak_het;  // placeholder
	h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_WRITE_LEN, NULL, NULL, rs); 
	ha_ct_hist(h, cnt, asm_opt->thread_num);
	ha_analyze_count(YAK_N_COUNTS, cnt, &peak_het);
	hamt_mark(asm_opt, rs, h, 0);  // initial pass of marking reads with `decide_category`
	hamt_mark(asm_opt, rs, h, 1);  // 1st pass of dropping reads
	reset_All_reads(rs);
	/* debug print (legacy) */
	// for (uint64_t rid=0; rid<rs->total_reads; rid++){
	// 	printf("rid:%" PRIu64 "\t", rid);
	// 	printf("%f\t%f\t%d\t%d\n", rs->mean[rid], rs->std[rid], rs->mask_readtype[rid], rs->mask_readnorm[rid]);
	// }
	/****************/

	// exit(0);

}

void *hamt_ft_gen(const hifiasm_opt_t *asm_opt, All_reads *rs, uint16_t coverage)
{   // use arbitrary coverage because there might be no peaks
	yak_ft_t *flt_tab;
	int64_t cnt[YAK_N_COUNTS];
	int peak_hom, peak_het, cutoff;
	ha_ct_t *h;
	h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_WRITE_LEN|HAMTF_FORCE_DONT_INIT, NULL, NULL, rs);
	ha_ct_hist(h, cnt, asm_opt->thread_num);
	peak_hom = ha_analyze_count(YAK_N_COUNTS, cnt, &peak_het);  // vanilla hifiasm histogram 
	if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
	cutoff = (int)(coverage * asm_opt->high_factor);
	if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
	ha_ct_shrink(h, cutoff, YAK_MAX_COUNT, asm_opt->thread_num);
	flt_tab = gen_hh(h);
	ha_ct_destroy(h);
	fprintf(stderr, "[M::%s::%.3f*%.2f@%.3fGB] ==> filtered out %ld k-mers occurring %d or more times\n", __func__,
			yak_realtime(), yak_cpu_usage(), yak_peakrss_in_gb(), (long)kh_size(flt_tab), cutoff);
	reset_All_reads(rs);
	return (void*)flt_tab;
}