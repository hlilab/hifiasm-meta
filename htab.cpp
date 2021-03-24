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

#define KTPIPE_NB_CPU 32  // threaded step 2

// #define HAMT_DIG_KMERRESCUE 30

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
void *ha_flt_tab_hp;
ha_pt_t *ha_idx_hp;
void *ha_ct_table;

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
///in most cases, n_shift = 25, n_hashes = 4
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
	uint64_t tot;  ///number of distinct k-mers
	ha_ct1_t *h;
} ha_ct_t;

///for 0-th counting, k = 51, pre = 12, n_hash = 4, n_shift = 37
///for 1-th counting, opt.k = 51, opt->pre = 12, opt->bf_n_hash = 4, opt.bf_shift = 0
static ha_ct_t *ha_ct_init(int k, int pre, int n_hash, int n_shift)
{
	ha_ct_t *h;
	int i;
	if (pre < YAK_COUNTER_BITS) return 0;
	CALLOC(h, 1);
	h->k = k, h->pre = pre;
	CALLOC(h->h, 1<<h->pre);
	///i<h->pre = 4096
	///it seems there is a large hash table h, consisting 4096 small hash tables
	for (i = 0; i < 1<<h->pre; ++i)
		h->h[i].h = yak_ct_init();
	///for 0-th counting, enter here; used for bloom filter
	if (n_hash > 0 && n_shift > h->pre) {
		h->n_hash = n_hash, h->n_shift = n_shift;
		for (i = 0; i < 1<<h->pre; ++i)
			h->h[i].b = yak_bf_init(h->n_shift - h->pre, h->n_hash); ///h->n_shift = 37, h->pre = 12, h->n_hash = 4 
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
	///corresponding small hash index
	g = &h->h[a[0]&mask];
	for (j = 0; j < n; ++j) {
		int ins = 1, absent;
		///x is a 64-bit word, h->pre=12
		///all elements at a have the same low 12 bits
		///so low 12 bits are not useful
		uint64_t x = a[j] >> h->pre;
		khint_t k;
		if ((a[j]&mask) != (a[0]&mask)) continue;
		if (create_new) {
			///for 0-th counting, g->b = NULL
			if (g->b)
				ins = (yak_bf_insert(g->b, x) == h->n_hash);
			///for 0-th counting, g->b = NULL
			///x = the high 52 bits of a[j] + low 12 bits 0
			///the low 12 bits are used for counting
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

///YAK_N_COUNTS is also 4096
///used for calculating k-mer histogram
static void ha_ct_hist(const ha_ct_t *h, int64_t cnt[YAK_N_COUNTS], int n_thread)
{
	hist_aux_t a;
	int i, j;
	a.h = h;
	memset(cnt, 0, YAK_N_COUNTS * sizeof(uint64_t));
	CALLOC(a.cnt, n_thread);
	///start 4096 threads
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
	///still start 4096 threads
	kt_for(n_thread, worker_ct_shrink, &a, 1<<h->pre);
	for (i = 0, h->tot = 0; i < 1<<h->pre; ++i)
		h->tot += kh_size(h->h[i].h);
}


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
			///this should be the start index of kh_key's corresponding pos at ha_idxpos_t* a
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

///p = 12
static inline void ct_insert_buf(ch_buf_t *buf, int p, uint64_t y) // insert a k-mer $y to a linear buffer
{
	///assign k-mer to one of the 4096 bins
	///using low 12 bits for assigning
	///so all elements at b have the same low 12 bits
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
	///assign minimizer to one of 4096 bins by low 12 bits
	int pre = y->x & ((1<<p) - 1);
	ch_buf_t *b = &buf[pre];
	if (b->n == b->m) {
		b->m = b->m < 8? 8 : b->m + (b->m>>1);
		REALLOC(b->b, b->m);
	}
	b->b[b->n++] = *y;
}

///buf is the read block, k is the k-mer length, p = 12, len is the read length, seq is the read
static void count_seq_buf(ch_buf_t *buf, int k, int p, int len, const char *seq) // insert k-mers in $seq to linear buffer $buf
{
	int i, l;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		///c = 00, 01, 10, 11
		if (c < 4) { // not an "N" base
			///x[0] & x[1] are the forward k-mer
			///x[2] & x[3] are the reverse complementary k-mer
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
#define HAMTF_FORCE_DONT_INIT 0x100  // meta: override HAF_RS_WRITE_LEN to not allow init_all_reads of ha_count
#define HAMTF_HAS_MARKS 0x80
#define HAF_SKIP_READ    0x40

typedef struct { // global data structure for kt_pipeline()
	const yak_copt_t *opt;
	const void *flt_tab;
	int flag, create_new, is_store;
	uint64_t n_seq; ///number of total reads
	kseq_t *ks;
	UC_Read ucr;
	ha_ct_t *ct;
	ha_pt_t *pt;
	const All_reads *rs_in;
	All_reads *rs_out;
	All_reads *rs;  // hamt meta. ALWAYS carry this

	// debug (profiling)
	double accumulated_time_step1;
	double accumulated_time_step2;
	double accumulated_time_step3;
} pl_data_t;

typedef struct { // data structure for each step in kt_pipeline()
	pl_data_t *p;
	uint64_t n_seq0;  ///the start index of current buffer block at R_INF
	///sum_len = total bases, nk = number of k-mers
	int n_seq, m_seq, sum_len, nk;
	int *len;
	char **seq;
	ha_mz1_v *mz_buf;
	ha_mz1_v *mz;
	ch_buf_t *buf;  // linear buffer
	ch_buf_t **threaded_buf;  // for threaded step2
} st_data_t;

static void worker_for_insert(void *data, long i, int tid) // callback for kt_for()
{
	st_data_t *s = (st_data_t*)data;
	ch_buf_t *b = &s->buf[i];
	if (s->p->pt)
		b->n_ins += ha_pt_insert_list(s->p->pt, b->n, b->b);
	else///for 0-th count, go into here
		b->n_ins += ha_ct_insert_list(s->p->ct, s->p->create_new, b->n, b->a);
}

static void worker_for_mz(void *data, long i, int tid)
{
	st_data_t *s = (st_data_t*)data;
	///get the corresponding minimzer vector of this read
	ha_mz1_v *b = &s->mz_buf[tid];
	s->mz_buf[tid].n = 0;
	///s->p->opt->w = 51, s->p->opt->k
	ha_sketch(s->seq[i], s->len[i], s->p->opt->w, s->p->opt->k, s->n_seq0 + i, s->p->opt->is_HPC, b, s->p->flt_tab);
	s->mz[i].n = s->mz[i].m = b->n;
	MALLOC(s->mz[i].a, b->n);
	memcpy(s->mz[i].a, b->a, b->n * sizeof(ha_mz1_t));
}

typedef struct{
	int is_use_minimizer;
	st_data_t *s;
	pl_data_t *p;
}pl_step2_data_t;
static void worker_count_step2sub_worker(void *data, long i, int tid){  // kt_for() callback
	// Trade memory for speed at kmer counting step 2 (the middle step).
	// The ovec will use a lot of mem anyway, so it doesn't matter too much
	//  if kmer counting gets a few more buffers.
	// A rough profiling by Get_T() suggested step2 is ~2x or slower than the other two steps.
	pl_step2_data_t *pp = (pl_step2_data_t*)data;
	pl_data_t *p = pp->p;
	st_data_t *s = pp->s;

	// handle read selection mask
	if (p->flag & HAMTF_HAS_MARKS){
		uint64_t rid = s->n_seq0+i;
		assert(rid<R_INF.total_reads);
		assert(rid<p->rs->hamt_stat_buf_size);
		if (p->rs->mask_readnorm[rid] & 1) // initial marking hasn't finished; read is marked as dropped
			{return;}
	}

	// count kmer/minimizers
	if (!pp->is_use_minimizer){  // all kmers
		if (p->opt->is_HPC){
			count_seq_buf_HPC(s->threaded_buf[tid], p->opt->k, p->opt->pre, s->len[i], s->seq[i]);
		}else{
			count_seq_buf(s->threaded_buf[tid], p->opt->k, p->opt->pre, s->len[i], s->seq[i]);
		}
	}else{  // minimizers
		if (p->pt){
			for (uint32_t j = 0; j < s->mz[i].n; ++j)
				pt_insert_buf(s->threaded_buf[tid], p->opt->pre, &s->mz[i].a[j]);
		}else{
			for (uint32_t j = 0; j < s->mz[i].n; ++j)
				ct_insert_buf(s->threaded_buf[tid], p->opt->pre, s->mz[i].a[j].x);
		}
	}
	if (!p->is_store) free(s->seq[i]);
}
void worker_count_step2sub(pl_data_t *p, st_data_t *s, int nb_cpu){
	int is_use_minimizer = !(p->opt->w==1);

	pl_step2_data_t pp;
	pp.is_use_minimizer = is_use_minimizer;
	pp.s = s;
	pp.p = p;
	
	kt_for(nb_cpu, worker_count_step2sub_worker, &pp, s->n_seq);
}

static void *worker_count(void *data, int step, void *in) // callback for kt_pipeline()
{
	pl_data_t *p = (pl_data_t*)data;
	double t_profiling = Get_T();
	if (step == 0) { // step 1: read a block of sequences
		int ret;
		st_data_t *s;
		CALLOC(s, 1);
		s->p = p;
		s->n_seq0 = p->n_seq;
		if (p->rs_in && (p->flag & HAF_RS_READ)) {
			while (p->n_seq < p->rs_in->total_reads) {
				if((p->flag & HAF_SKIP_READ) && p->rs_in->trio_flag[p->n_seq] != AMBIGU)
				{
					++p->n_seq;
					continue;
				}
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
			while ((ret = kseq_read(p->ks)) >= 0) {
				int l = p->ks->seq.l;
				if (p->n_seq >= 1<<28) {
					fprintf(stderr, "ERROR: this implementation supports no more than %d reads\n", 1<<28);
					exit(1);
				}
				if (p->rs_out) {
					///for 0-th count, just insert read length to R_INF, instead of read
					if (p->flag & HAF_RS_WRITE_LEN) {
						assert(p->n_seq == p->rs_out->total_reads);
						ha_insert_read_len(p->rs_out, l, p->ks->name.l);
					} else if (p->flag & HAF_RS_WRITE_SEQ) {
						int i, n_N;
						assert(l == (int)p->rs_out->read_length[p->n_seq]);
						for (i = n_N = 0; i < l; ++i) // count number of ambiguous bases
							if (seq_nt4_table[(uint8_t)p->ks->seq.s[i]] >= 4)
								++n_N;
						ha_compress_base(Get_READ(*p->rs_out, p->n_seq), p->ks->seq.s, l, &p->rs_out->N_site[p->n_seq], n_N);
						memcpy(&p->rs_out->name[p->rs_out->name_index[p->n_seq]], p->ks->name.s, p->ks->name.l);
					}
				}
				///for 0-th count, insert both seq and length to local block
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
				///p->opt->chunk_size is the block max size
				if (s->sum_len >= p->opt->chunk_size)
					break;
			}
		}
		p->accumulated_time_step1 += Get_T() - t_profiling;
		if (s->sum_len == 0) free(s);
		else return s;
	} else if (step == 1) { // step 2: extract k-mers
		st_data_t *s = (st_data_t*)in;
		int i, n_pre = 1<<p->opt->pre, m;
		int nb_cpu = p->opt->n_thread<=1? 1 : (p->opt->n_thread>KTPIPE_NB_CPU? KTPIPE_NB_CPU : p->opt->n_thread);

		// allocate the k-mer buffer
		// CALLOC(s->buf, n_pre);
		m = (int)(s->nk * 1.2 / n_pre) + 1;
		s->threaded_buf = (ch_buf_t**)calloc(nb_cpu, sizeof(ch_buf_t*));
		for (int i_cpu=0; i_cpu<nb_cpu; i_cpu++){
			CALLOC(s->threaded_buf[i_cpu], n_pre);
			for (i = 0; i < n_pre; ++i) {
				///for 0-th counting, p->pt = NULL
				s->threaded_buf[i_cpu][i].m = m;
				if (p->pt) MALLOC(s->threaded_buf[i_cpu][i].b, m);  // if (p->pt) MALLOC(s->buf[i].b, m);
				else MALLOC(s->threaded_buf[i_cpu][i].a, m);  // else MALLOC(s->buf[i].a, m);
			}
		}


		// fill the buffer
		if (p->opt->w == 1) { // enumerate all k-mers
			worker_count_step2sub(p, s, nb_cpu);
		} else { // minimizers only
			uint32_t j;
			// compute minimizers
			// s->n_seq is how many reads at this buffer
			// s->mz && s->mz_buf are lists of minimzer vectors
			CALLOC(s->mz, s->n_seq);
			CALLOC(s->mz_buf, p->opt->n_thread);
			///calculate minimzers for each read, each read corresponds to one thread
			kt_for(p->opt->n_thread, worker_for_mz, s, s->n_seq);
			for (i = 0; i < p->opt->n_thread; ++i)
				free(s->mz_buf[i].a);
			free(s->mz_buf);
			
			// insert minimizers
			worker_count_step2sub(p, s, nb_cpu);

			for (i = 0; i < s->n_seq; ++i) {
				free(s->mz[i].a);
			}
			free(s->mz);
		}
		p->accumulated_time_step2 += Get_T() - t_profiling;
		///just clean seq
		free(s->seq); free(s->len);
		s->seq = 0, s->len = 0;
		return s;
	} else if (step == 2) { // step 3: insert k-mers to hash table
		st_data_t *s = (st_data_t*)in;
		int i, n = 1<<p->opt->pre;
		uint64_t n_ins = 0;
		///for 0-th counting, p->pt = NULL

		int nb_cpu = p->opt->n_thread<=1? 1 : (p->opt->n_thread>KTPIPE_NB_CPU? KTPIPE_NB_CPU : p->opt->n_thread);
		for (int i_cpu=0; i_cpu<nb_cpu; i_cpu++){
			s->buf = s->threaded_buf[i_cpu];
			kt_for(p->opt->n_thread, worker_for_insert, s, n);
			for (int i_mer=0; i_mer<n; i_mer++){
				n_ins += s->buf[i_mer].n_ins;
			}
		}
		if (p->ct) p->ct->tot += n_ins;
		if (p->pt) p->pt->tot_pos += n_ins;
		
		///n_ins is number of distinct k-mers
		for (int i_cpu=0; i_cpu<nb_cpu; i_cpu++){
			for (i = 0; i < n; ++i) {
				if (p->pt) free(s->threaded_buf[i_cpu][i].b);
				else free(s->threaded_buf[i_cpu][i].a);
			}
			free(s->threaded_buf[i_cpu]);
		}
		free(s->threaded_buf);
		#if 0
		fprintf(stderr, "[M::%s::%.3f*%.2f] processed %ld sequences; %ld %s in the hash table\n", __func__,
				yak_realtime(), yak_cpu_usage(), (long)s->n_seq0 + s->n_seq,
				(long)(p->pt? p->pt->tot_pos : p->ct->tot), p->pt? "positions" : "distinct k-mers");
		#endif
		free(s);
		p->accumulated_time_step3 += Get_T() - t_profiling;
	}
	return 0;
}

static ha_ct_t *yak_count(const yak_copt_t *opt, const char *fn, int flag, ha_pt_t *p0, ha_ct_t *c0, const void *flt_tab, All_reads *rs, int64_t *n_seq)
{
	///for 0-th counting, flag = HAF_COUNT_ALL|HAF_RS_WRITE_LEN|HAF_CREATE_NEW
	int read_rs = (rs && (flag & HAF_RS_READ));
	pl_data_t pl;
	gzFile fp = 0;
	memset(&pl, 0, sizeof(pl_data_t));
	pl.n_seq = *n_seq;
	if (read_rs) {
		pl.rs_in = rs;
		init_UC_Read(&pl.ucr);
	} else {///for 0-th counting, go into here
		if ((fp = gzopen(fn, "r")) == 0) return 0;
		pl.ks = kseq_init(fp);
	}
	///for 0-th counting, read all reads into pl.rs_out
	if (rs && (flag & (HAF_RS_WRITE_LEN|HAF_RS_WRITE_SEQ)))
		pl.rs_out = rs;
	pl.rs = rs;
	///for 0-th counting, flt_tab = NULL
	///for 1-th counting, flt_tab = NULL
	pl.flt_tab = flt_tab;
	pl.opt = opt;
	pl.flag = flag;
	if (p0) {///for 1-th counting, p0 = NULL
		pl.pt = p0, pl.create_new = 0; // never create new elements in a position table
		assert(p0->k == opt->k && p0->pre == opt->pre);
	} else if (c0) {
		pl.ct = c0, pl.create_new = !!(flag&HAF_CREATE_NEW);
		assert(c0->k == opt->k && c0->pre == opt->pre);
	} else {///for ft-th counting and 1-th counting, go into here
		pl.create_new = 1; // alware create new elements if the count table is empty
		///for 0-th counting, opt.k = 51, opt->pre = 12, opt->bf_n_hash = 4, opt.bf_shift = 37
		///for 1-th counting, opt.k = 51, opt->pre = 12, opt->bf_n_hash = 4, opt.bf_shift = 0
		///building a large hash table consisting of 4096 small hash tables
		pl.ct = ha_ct_init(opt->k, opt->pre, opt->bf_n_hash, opt->bf_shift);
	}
	kt_pipeline(3, worker_count, &pl, 3);

	fprintf(stderr, "[prof::%s] step 1 total %.2f s, step2 %.2f s, step3 %.2f s.\n",
							__func__, pl.accumulated_time_step1, 
									  pl.accumulated_time_step2,
									  pl.accumulated_time_step3);

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
	///for 0-th counting, flag = HAF_COUNT_ALL|HAF_RS_WRITE_LEN
	if (rs) {
		if (flag & HAF_RS_WRITE_LEN){
		    if (!(flag & HAMTF_FORCE_DONT_INIT)){
				init_All_reads(rs);
			}
			else{
				reset_All_reads(rs);
			}
		}
		else if (flag & HAF_RS_WRITE_SEQ){
			malloc_All_reads(rs);
		}
	}
	yak_copt_init(&opt);
	opt.k = asm_opt->k_mer_length;
	///always 0
	opt.is_HPC = !(asm_opt->flag&HA_F_NO_HPC);
	///for ft-counting, shoud be 1
	opt.w = flag & HAF_COUNT_ALL? 1 : asm_opt->mz_win;
	///for ft-counting, shoud be 37
	///for ha_pt_gen, shoud be 0
	opt.bf_shift = flag & HAF_COUNT_EXACT? 0 : asm_opt->bf_shift;
	opt.n_thread = asm_opt->thread_num;
	///asm_opt->num_reads is the number of fastq files
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


void debug_ct_index(void* q_ct_idx, void* r_ct_idx)
{
	ha_ct_t* ct_idx = (ha_ct_t*)q_ct_idx;
	yak_ct_t *g = NULL;
	uint64_t i;
	khint_t k;
	for (i = 0; (int)i < 1<<ct_idx->pre; i++)
	{
		g = ct_idx->h[i].h;
		for (k = 0; k < kh_end(g); ++k) 
		{
			if (kh_exist(g, k)) 
			{
				int c = kh_key(g, k) & YAK_MAX_COUNT;
				uint64_t hash = ((kh_key(g, k) >> ct_idx->pre)<<ct_idx->pre) | i;
				int q = query_ct_index(r_ct_idx, hash);
				if(q!=c)
				{
					fprintf(stderr, "ERROR:c: %d, q: %d\n", c, q);
				}
			}
		}
	}
}

/*************************
 * High-level interfaces *
 *************************/

void *ha_ft_gen(const hifiasm_opt_t *asm_opt, All_reads *rs, int *hom_cov, int is_hp_mode)
{
	yak_ft_t *flt_tab;
	int64_t cnt[YAK_N_COUNTS];
	int peak_hom, peak_het, cutoff = YAK_MAX_COUNT - 1, ex_flag = 0;
	if(is_hp_mode) ex_flag = HAF_RS_READ|HAF_SKIP_READ;
	ha_ct_t *h;
	h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_WRITE_LEN|ex_flag, NULL, NULL, rs);
	if((asm_opt->flag & HA_F_VERBOSE_GFA) && asm_opt->write_new_graph_bins)
	{
		write_ct_index((void*)h, asm_opt->output_file_name);
		// load_ct_index(&ha_ct_table, asm_opt->output_file_name);
		// debug_ct_index((void*)h, ha_ct_table);
		// debug_ct_index(ha_ct_table, (void*)h);
		// ha_ct_destroy((ha_ct_t *)ha_ct_table);
	}
	
	if(!(ex_flag & HAF_SKIP_READ))
	{
		ha_ct_hist(h, cnt, asm_opt->thread_num);
		peak_hom = ha_analyze_count(YAK_N_COUNTS, asm_opt->min_hist_kmer_cnt, cnt, &peak_het);
		if (hom_cov) *hom_cov = peak_hom;
		if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
		///in default, asm_opt->high_factor = 5.0
		cutoff = (int)(peak_hom * asm_opt->high_factor);
		if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
	}
	ha_ct_shrink(h, cutoff, YAK_MAX_COUNT, asm_opt->thread_num);
	flt_tab = gen_hh(h);
	ha_ct_destroy(h);
	fprintf(stderr, "[M::%s::%.3f*%.2f@%.3fGB] ==> filtered out %ld k-mers occurring %d or more times\n", __func__,
			yak_realtime(), yak_cpu_usage(), yak_peakrss_in_gb(), (long)kh_size(flt_tab), cutoff);
	return (void*)flt_tab;
}

ha_pt_t *ha_pt_gen(const hifiasm_opt_t *asm_opt, const void *flt_tab, int read_from_store, int is_hp_mode, All_reads *rs, int *hom_cov, int *het_cov)
{
	fprintf(stderr, "called ha_pt_gen, RS status: has nothing %d, has lengths %d, all in mem %d\n", R_INF.is_has_nothing, R_INF.is_has_lengths, R_INF.is_all_in_mem);
	int64_t cnt[YAK_N_COUNTS], tot_cnt;
	int peak_hom, peak_het, i, extra_flag1, extra_flag2;
	ha_ct_t *ct;
	ha_pt_t *pt;
	if (read_from_store) {///if reads have already been read
		extra_flag1 = extra_flag2 = HAF_RS_READ;
	} else if (rs->total_reads == 0) {///if reads & length have not been scanned
		extra_flag1 = HAF_RS_WRITE_LEN;
		extra_flag2 = HAF_RS_WRITE_SEQ;
	} else {///if length has been loaded but reads have not
		extra_flag1 = HAF_RS_WRITE_SEQ;
		extra_flag2 = HAF_RS_READ;
	}
	if (asm_opt->is_use_exp_graph_cleaning){
		fprintf(stderr, "[M::%s] counting - minimzers\n", __func__);
		ct = ha_count(asm_opt, HAF_COUNT_EXACT|extra_flag1|HAMTF_FORCE_DONT_INIT, NULL, flt_tab, rs);  // collect minimizer (aware of high freq filter)
	}else{
		if(is_hp_mode) extra_flag1 |= HAF_SKIP_READ, extra_flag2 |= HAF_SKIP_READ;

		ct = ha_count(asm_opt, HAF_COUNT_EXACT|extra_flag1, NULL, flt_tab, rs);
	}
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> counted %ld distinct minimizer k-mers\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)ct->tot);
	ha_ct_hist(ct, cnt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s] count[%d] = %ld (for sanity check)\n", __func__, YAK_MAX_COUNT, (long)cnt[YAK_MAX_COUNT]);
	peak_hom = ha_analyze_count(YAK_N_COUNTS, asm_opt->min_hist_kmer_cnt, cnt, &peak_het);
	if (hom_cov) *hom_cov = peak_hom;
	if (het_cov) *het_cov = peak_het;
	if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
	///here ha_ct_shrink is mostly used to remove k-mer appearing only 1 time
	if (flt_tab == 0) {
		int cutoff = (int)(peak_hom * asm_opt->high_factor);
		if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
		if((extra_flag1 & HAF_SKIP_READ) && (extra_flag2 & HAF_SKIP_READ)) cutoff = YAK_MAX_COUNT - 1;
		ha_ct_shrink(ct, 2, cutoff, asm_opt->thread_num);
		for (i = 2, tot_cnt = 0; i <= cutoff; ++i) tot_cnt += cnt[i] * i;
	} else {
		///Note: here is just to remove minimizer appearing YAK_MAX_COUNT times
		///minimizer with YAK_MAX_COUNT occ may apper > YAK_MAX_COUNT times, so it may lead to overflow at ha_pt_gen
		ha_ct_shrink(ct, 2, YAK_MAX_COUNT - 1, asm_opt->thread_num);
		for (i = 2, tot_cnt = 0; i <= YAK_MAX_COUNT - 1; ++i) tot_cnt += cnt[i] * i;
	}
	pt = ha_pt_gen(ct, asm_opt->thread_num); // prepare the slots
	fprintf(stderr, "[M::%s] counting - minimzer positions\n", __func__);
	ha_count(asm_opt, HAF_COUNT_EXACT|extra_flag2|HAMTF_FORCE_DONT_INIT, pt, flt_tab, rs);  // collect positional info
	
	fprintf(stderr, "[debug::%s] tot_cnt is %" PRIu64", pt->tot_pos is %" PRIu64 "\n", __func__,
					tot_cnt, pt->tot_pos);
	assert((uint64_t)tot_cnt == pt->tot_pos);
	
	//ha_pt_sort(pt, asm_opt->thread_num);
	fprintf(stderr, "[M::%s::%.3f*%.2f] ==> indexed %ld positions\n", __func__,
			yak_realtime(), yak_cpu_usage(), (long)pt->tot_pos);

	if (rs->is_has_nothing){
		rs->is_has_nothing = 0;
		rs->is_has_lengths = 1;
	}else if (rs->is_has_lengths && (!rs->is_all_in_mem)){
		rs->is_all_in_mem = 1;
	}
	
	return pt;
}

int query_ct_index(void* ct_idx, uint64_t hash)
{
	ha_ct1_t *g = &(((ha_ct_t*)ct_idx)->h[hash & ((1ULL<<((ha_ct_t*)ct_idx)->pre) - 1)]);
	khint_t k;
	k = yak_ct_get(g->h, hash);
	if (k == kh_end(g->h)) return 0;
	return kh_key(g->h, k)&YAK_MAX_COUNT;
}


int write_ct_index(void *i_ct_idx, char* file_name)
{
	char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.ct_flt", file_name);
    FILE* fp = fopen(gfa_name, "w");
	if (!fp) {
		free(gfa_name);
        return 0;
    }
	ha_ct_t* ct_idx = (ha_ct_t*)i_ct_idx;
	int i;
	ha_ct1_t *g;
	fwrite(&ct_idx->k, sizeof(ct_idx->k), 1, fp);
	fwrite(&ct_idx->pre, sizeof(ct_idx->pre), 1, fp);
	fwrite(&ct_idx->n_hash, sizeof(ct_idx->n_hash), 1, fp);
	fwrite(&ct_idx->n_shift, sizeof(ct_idx->n_shift), 1, fp);
	fwrite(&ct_idx->tot, sizeof(ct_idx->tot), 1, fp);
	for (i = 0; i < 1<<ct_idx->pre; i++)
	{
		g = &(ct_idx->h[i]);
		yak_ct_save(g->h, fp);
	}


	fprintf(stderr, "[M::%s] Index has been written.\n", __func__);
    free(gfa_name);
	fclose(fp);
	return 1;
}

int load_ct_index(void **i_ct_idx, char* file_name)
{
	char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.ct_flt", file_name);
    FILE* fp = fopen(gfa_name, "r");
	if (!fp) {
		free(gfa_name);
        return 0;
    }
	ha_ct_t** ct_idx = (ha_ct_t**)i_ct_idx;
	double index_time = 0;
	int i;
	ha_ct_t *h = 0;
	ha_ct1_t *g;
	CALLOC(h, 1);

	fread(&h->k, sizeof(h->k), 1, fp);
	fread(&h->pre, sizeof(h->pre), 1, fp);
	fread(&h->n_hash, sizeof(h->n_hash), 1, fp);
	fread(&h->n_shift, sizeof(h->n_shift), 1, fp);
	fread(&h->tot, sizeof(h->tot), 1, fp);
	CALLOC(h->h, 1<<h->pre);

	index_time = yak_realtime();
	for (i = 0; i < 1<<h->pre; ++i) 
	{
		g = &(h->h[i]);
		yak_ct_load(&(g->h), fp);
	}

	(*ct_idx) = h;
	fprintf(stderr, "[M::%s::%.3f] ==> Loaded count table\n", __func__, yak_realtime() - index_time);
	fprintf(stderr, "[M::%s] Index has been loaded.\n", __func__);
	free(gfa_name);
	return 1;
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

	uint8_t flag;
	const hifiasm_opt_t *asm_opt;  // hamt experimental, used to toggle things like read selection order. TODO: clean up when ready.

	int nb_reads_kept;  
	int nb_reads_limit;
	int runtime_median_threshold;

	// debug (profiling)
	double accumulated_time_step1;
	double accumulated_time_step2;
	double accumulated_time_step3;

}plmt_data_t;

typedef struct {  // step data structure (st_data_t)
	plmt_data_t* p;
	uint64_t n_seq0;
	int n_seq, m_seq, sum_len, nk;
	int *len;
	char **seq;
	hamt_ch_buf_t *lnbuf;  // $n_seq linear buffers
	hamt_ch_buf_t **threaded_lnbuf;  // step 2 is slow, trade memory for speed
	// exp
	uint64_t *RIDs;
}plmt_step_t;

static void worker_insert_lnbuf(void* data, long idx_for, int tid){
	plmt_step_t *s = (plmt_step_t*) data;
	uint64_t *buf = s->lnbuf[idx_for].a;
	uint64_t n = s->lnbuf[idx_for].n;
	yak_ct_t *h = s->p->hd->h[idx_for].h;
	int absent;
	khint_t key;
	uint64_t san = 0;
	for (uint32_t i=0; i<n; i++){
		key = yak_ct_put(h, buf[i]>>s->p->h->pre<<YAK_COUNTER_BITS1, &absent);
		if ((kh_key(h, key)&YAK_MAX_COUNT)<YAK_MAX_COUNT) kh_key(h, key)++;
		if (absent) san++;
	}
}


#define HAMT_DISCARD 0x1
#define HAMT_VIA_MEDIAN 0x2
#define HAMT_VIA_LONGLOW 0x4
#define HAMT_VIA_KMER 0x8
void worker_process_one_read_HPC(plmt_step_t *s, int idx_seq){
	uint64_t rid;
	if (/*s->p->asm_opt->readselection_sort_order==0 || */s->p->round==0)
		{rid = s->n_seq0+idx_seq;}
	else  // TODO: is legacy
		{rid = s->RIDs[idx_seq];}
	uint16_t *buf = (uint16_t*)malloc(sizeof(uint16_t)*s->len[idx_seq]);  // buffer used in the 0th round
	int k = s->p->opt->k;
	ha_ct_t *h = s->p->h;

	uint32_t idx = 0;  // index of kmer count in the buffer
	khint_t key;
	int i, l, last = -1;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	hamt_ch_buf_t *b = 0;// uint64_t *b = 0;

	int has_wanted_kmers = 0;  // flag of rescue

	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->len[idx_seq]; ++i) {
		int c = seq_nt4_table[(uint8_t)s->seq[idx_seq][i]];
		if (c < 4) { // not an "N" base
			if (c != last) {
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
				if (++l >= k){
					uint64_t hash = yak_hash_long(x);
					// collect for all-reads kmer count hashtable
					yak_ct_t *h_ = h->h[hash & ((1<<h->pre)-1)].h;
					key = yak_ct_get(h_, hash>>h->pre<<YAK_COUNTER_BITS1);
					if (key!=kh_end(h_)){buf[idx] = (kh_key(h_, key) & YAK_MAX_COUNT);}
					else buf[idx] = 1;  // because bf
					idx++;
				}
				last = c;
			}
		} else l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}

	double mean = meanl(buf, idx);
	double std = stdl(buf, idx, mean);
	radix_sort_hamt16(buf, buf+idx);
	uint16_t median = (uint16_t)buf[idx/2];
	uint16_t lowq = (uint16_t)buf[idx/10];
	// uint8_t code = decide_category(mean, std, buf, idx);
	s->p->rs_out->mean[rid] = mean;
	s->p->rs_out->median[rid] = median;
	s->p->rs_out->lowq[rid] = lowq;
	s->p->rs_out->std[rid] = std;
	s->p->rs_out->mask_readtype[rid] = /*code*/ 0;

	free(buf);  // sequence will be freed by the caller	
}
void worker_process_one_read_noHPC(plmt_step_t *s, int idx_seq){
	// mirrored from the HPC version above (Dec 22 2020, ~r20)
	uint64_t rid;
	if (/*s->p->asm_opt->readselection_sort_order==0 || */s->p->round==0)
		{rid = s->n_seq0+idx_seq;}
	else  // TODO: is legacy
		{rid = s->RIDs[idx_seq];}
	uint16_t *buf = (uint16_t*)malloc(sizeof(uint16_t)*s->len[idx_seq]);  // buffer used in the 0th round
	int k = s->p->opt->k;
	ha_ct_t *h = s->p->h;

	uint32_t idx = 0;  // index of kmer count in the buffer
	khint_t key;
	int i, l/*, last = -1*/;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	hamt_ch_buf_t *b = 0;// uint64_t *b = 0;

	int has_wanted_kmers = 0;  // flag of rescue

	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->len[idx_seq]; ++i) {
		int c = seq_nt4_table[(uint8_t)s->seq[idx_seq][i]];
		if (c < 4) { // not an "N" base
			// if (c != last) {
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
				if (++l >= k){
					uint64_t hash = yak_hash_long(x);
					// collect for all-reads kmer count hashtable
					yak_ct_t *h_ = h->h[hash & ((1<<h->pre)-1)].h;
					key = yak_ct_get(h_, hash>>h->pre<<YAK_COUNTER_BITS1);
					if (key!=kh_end(h_)){buf[idx] = (kh_key(h_, key) & YAK_MAX_COUNT);}
					else buf[idx] = 1;  // because bf
					idx++;
				}
				// last = c;
			// }
		} else l = 0, /*last = -1,*/ x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}

	double mean = meanl(buf, idx);
	double std = stdl(buf, idx, mean);
	radix_sort_hamt16(buf, buf+idx);
	uint16_t median = (uint16_t)buf[idx/2];
	uint16_t lowq = (uint16_t)buf[idx/10];
	// uint8_t code = decide_category(mean, std, buf, idx);
	s->p->rs_out->mean[rid] = mean;
	s->p->rs_out->median[rid] = median;
	s->p->rs_out->lowq[rid] = lowq;
	s->p->rs_out->std[rid] = std;
	s->p->rs_out->mask_readtype[rid] = /*code*/ 0;

	free(buf);  // sequence will be freed by the caller	
}

typedef struct{
	plmt_data_t *p;
	plmt_step_t *s;
}pl_flt_withsorting;

static void callback_worker_process_one_read_either(void *data, long i, int tid){
	// callback of kt_for for 
	//   - worker_process_one_read_HPC
	//   - worker_process_one_read_noHPC
	// (used by hamt_flt_withsorting)
	pl_flt_withsorting *d = (pl_flt_withsorting*)data;
	if (d->p->opt->is_HPC)
		worker_process_one_read_HPC(d->s, i);
	else
		worker_process_one_read_noHPC(d->s, i);
}


// single thread version (deprecated)
#if 0
void worker_process_one_read_inclusive(plmt_step_t *s, int idx_seq, int round, int use_HPC){
	// note: this function is not called by kt_for, it's linear
	uint64_t rid = s->RIDs[idx_seq];
	if (s->p->rs_in->mask_readnorm[rid]==0){  // read is already kept
		return;
	}
	// fprintf(stderr, "DEBUG: idx_seq is %d, rid is %d.\n", (int)idx_seq, (int)rid);
	static int flag_round0_said_warning = 0;
	static int flag_round1_said_warning = 0;

	// early termination criteria
	if (round==0){  // keep globally infrequent reads
		// fprintf(stderr, "[process_one_read round0 dbg] read %d, median %d, lowq %d\n", (int)rid, (int) s->p->rs_in->median[rid] , (int)s->p->rs_in->lowq[rid]);
		if (s->p->nb_reads_kept>s->p->nb_reads_limit){
			if (!flag_round0_said_warning){
				fprintf(stderr, "[W::%s] read limit reached while recruiting lowq reads. continue anyway\n", __func__);
				flag_round0_said_warning = 1;
			}
		}
		if (s->p->rs_in->median[rid]>300 || s->p->rs_in->lowq[rid]>50){
			// fprintf(stderr, "[debug::%s] round 0, median %d, lowq %d\n", __func__, 
			// 							(int)s->p->rs_in->median[rid],
			// 							(int)s->p->rs_in->lowq[rid]);
			return ;
		}
		// fprintf(stderr, "[debug::%s] kept a read\n", __func__);
	} else if (round==1){  // recruit remaining reads
		// fprintf(stderr, "round1DEBUG: read %d, round %d, median %d, lowq %d.\n", (int)rid, round, (int)s->p->rs_in->median[rid], (int)s->p->rs_in->lowq[rid]);
		if ((s->p->rs_in->mask_readnorm[rid] & 1)==0){  // this read is already kept
			return ;
		}
		if (s->p->nb_reads_kept>s->p->nb_reads_limit){
			if (!flag_round1_said_warning){
				fprintf(stderr, "[W::%s] read limit reached in round#2.\n", __func__);
				flag_round1_said_warning = 1;
			}
			return ;
		}
	} else {
		fprintf(stderr, "ERROR: %s invalid round.\n", __func__);
		exit(1);
	}

	// collect kmers 
	uint64_t *buf = (uint64_t*)malloc(sizeof(uint64_t)*s->len[idx_seq]);  
	uint16_t *buf_norm = (uint16_t*)malloc(sizeof(uint16_t)*s->len[idx_seq]);  // runtime kmer frequency container
	int idx_buf = 0;
	int k = s->p->opt->k;  
	ha_ct_t *hd = s->p->hd;

	uint32_t idx = 0;  // index of kmer count in the buffer
	khint_t key;
	int i, l, last = -1;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	hamt_ch_buf_t *b = 0;// uint64_t *b = 0;

	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->len[idx_seq]; ++i) {
		int c = seq_nt4_table[(uint8_t)s->seq[idx_seq][i]];
		if (c < 4) { // not an "N" base
			// if (c != last) {
			if (c!=last || (!use_HPC)){
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
				if (++l >= k){
					uint64_t hash = yak_hash_long(x); 
					buf[idx_buf++] = hash;
					yak_ct_t *h_ = hd->h[hash & ((1<<hd->pre)-1)].h;  // runtime count hashtable
					key = yak_ct_get(h_, hash>>hd->pre<<YAK_COUNTER_BITS1);
					if (key!=kh_end(h_)){buf_norm[idx] = (kh_key(h_, key) & YAK_MAX_COUNT);}
					else buf_norm[idx] = 1;  // because bf
					idx++;
				}
				last = c;
			}
		} else l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
	
	int flag_is_keep_read = 0;
	if (round==0) {flag_is_keep_read = 1;}  // median has been checked upfront
	else if (round==1){
		radix_sort_hamt16(buf_norm, buf_norm+idx);
		// if (buf_norm[idx/2]<=s->p->runtime_median_threshold || buf_norm[idx/10]<=150 || buf_norm[idx/20]<=80){
		if (// buf_norm[idx/2]<=s->p->runtime_median_threshold || 
			buf_norm[idx/10]<=s->p->asm_opt->lowq_thre_10 /*|| buf_norm[idx/20]<=s->p->asm_opt->lowq_thre_5*/){
			flag_is_keep_read = 1;
		}
		// else{
		// 	fprintf(stderr, "   debug: runtime median is %d, runtime lowq is %d, discarded. (readID %d) \n", (int)buf_norm[idx/2], (int)buf_norm[idx/10], (int)rid);
		// }
	}

	// flip mask and update linear buffer
	if (flag_is_keep_read){
		s->p->rs_out->mask_readnorm[rid] = 0;
		s->p->nb_reads_kept++;
		for (i=0; i<idx_buf; i++){
			b = &s->lnbuf[buf[i] & ((1<<hd->pre)-1)];// inert into linear buffer for kmer counting
			if ((b->n+3)>=b->m){
				b->m = b->m<16? 16 : b->m+((b->m)>>1);
				REALLOC(b->a, b->m);
				assert(b->a);
			}
			b->a[b->n++] = buf[i];
		}
	}
	free(buf);
	free(buf_norm);
}
#endif

// threaded version of worker_process_one_read_inclusive
void worker_process_one_read_inclusive2(plmt_step_t *s, int idx_seq, int round, int use_HPC, int idx_cpu){
	// no race, each thread has its own linear buffer.
	uint64_t rid = s->RIDs[idx_seq];
	if (s->p->rs_in->mask_readnorm[rid]==0){  // read is already kept
		return;
	}
	static int flag_round0_said_warning = 0;
	static int flag_round1_said_warning = 0;

	// early termination criteria
	if (round==0){  // keep globally infrequent reads
		if (s->p->nb_reads_kept>s->p->nb_reads_limit){
			if (!flag_round0_said_warning){
				fprintf(stderr, "[W::%s] read limit reached while recruiting lowq reads. continue anyway\n", __func__);
				flag_round0_said_warning = 1;
			}
		}
		if (s->p->rs_in->median[rid]>300 || s->p->rs_in->lowq[rid]>50){
			return ;
		}
	} else if (round==1){  // recruit remaining reads
		if ((s->p->rs_in->mask_readnorm[rid] & 1)==0){  // this read is already kept
			return ;
		}
		if (s->p->nb_reads_kept>s->p->nb_reads_limit){
			if (!flag_round1_said_warning){
				fprintf(stderr, "[W::%s] read limit reached in round#2.\n", __func__);
				flag_round1_said_warning = 1;
			}
			return ;
		}
	} else {
		fprintf(stderr, "ERROR: %s invalid round.\n", __func__);
		exit(1);
	}

	// NOTE
	//   (Oct 6 2020) ASan thought there's a heap-use-after-free error right here 
	//		(or a few lines below at the yak_ct_get part).
	//   (Jan 7 2021) This happend again with r024-test. It was runing on a 5% subsampled
	//       readset of the mock dataset. Only ASan complains, no segfaults.
	//   (Mar 21 2021) Maybe not a false positive; the pipeline used to seperate hashtable insertion
	//       into a 3rd step, which races with 2nd step's hashtable querying. It probably didn't
	//       segfault or cause instable results because 2nd step used to be single threaded.
	//       Fixed.

	// collect kmers 
	uint64_t *buf = (uint64_t*)malloc(sizeof(uint64_t)*s->len[idx_seq]);  
	uint16_t *buf_norm = (uint16_t*)malloc(sizeof(uint16_t)*s->len[idx_seq]);  // runtime kmer frequency container
	int idx_buf = 0;
	int k = s->p->opt->k;  
	ha_ct_t *hd = s->p->hd;

	uint32_t idx = 0;  // index of kmer count in the buffer
	khint_t key;
	int i, l, last = -1;
	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	hamt_ch_buf_t *b = 0;// uint64_t *b = 0;

	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < s->len[idx_seq]; ++i) {
		int c = seq_nt4_table[(uint8_t)s->seq[idx_seq][i]];
		if (c < 4) { // not an "N" base
			// if (c != last) {
			if (c!=last || (!use_HPC)){
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
				if (++l >= k){
					uint64_t hash = yak_hash_long(x); 
					buf[idx_buf++] = hash;
					yak_ct_t *h_ = hd->h[hash & ((1<<hd->pre)-1)].h;  // runtime count hashtable
					key = yak_ct_get(h_, hash>>hd->pre<<YAK_COUNTER_BITS1);
					if (key!=kh_end(h_)){buf_norm[idx] = (kh_key(h_, key) & YAK_MAX_COUNT);}
					else buf_norm[idx] = 1;  // because bf
					idx++;
				}
				last = c;
			}
		} else l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
	
	int flag_is_keep_read = 0;
	if (round==0) {flag_is_keep_read = 1;}  // median has been checked upfront
	else if (round==1){
		radix_sort_hamt16(buf_norm, buf_norm+idx);
		if (buf_norm[idx/10]<=s->p->asm_opt->lowq_thre_10 || buf_norm[idx/20]<=s->p->asm_opt->lowq_thre_5||
		    ((s->p->asm_opt->lowq_thre_3>0) && (buf_norm[idx/33]<=s->p->asm_opt->lowq_thre_3)) ){
			flag_is_keep_read = 1;
		}
	}

	// flip mask and update linear buffer
	if (flag_is_keep_read){
		s->p->rs_out->mask_readnorm[rid] = 0;
		s->p->nb_reads_kept++;
		for (i=0; i<idx_buf; i++){
			// b = &s->lnbuf[buf[i] & ((1<<hd->pre)-1)];// inert into linear buffer for kmer counting
			b = &s->threaded_lnbuf[idx_cpu][buf[i] & ((1<<hd->pre)-1)];// inert into linear buffer for kmer counting
			if ((b->n+3)>=b->m){
				b->m = b->m<16? 16 : b->m+((b->m)>>1);
				REALLOC(b->a, b->m);
				assert(b->a);
			}
			b->a[b->n++] = buf[i];
		}
	}
	free(buf);
	free(buf_norm);
}
static void callback_worker_process_one_read_inclusive2(void *data, long i, int tid){
	// callback of kt_for
	plmt_step_t *s = (plmt_step_t*)data;
	worker_process_one_read_inclusive2(s, (int)i, s->p->round, s->p->opt->is_HPC, tid);
}

static void *worker_mark_reads(void *data, int step, void *in){  // callback of kt_pipeline
	plmt_data_t *p = (plmt_data_t*)data;
	double t_profiling = Get_T();
	if (step==0){ // step 1: read a block of sequences
		int ret;
		plmt_step_t *s;
		CALLOC(s, 1);
		s->p = p;
		s->n_seq0 = *p->n_seq;
		s->lnbuf = 0;
		if (p->round!=0) {  // init kmer counting buffers
			s->lnbuf = (hamt_ch_buf_t*)calloc(1<<p->h->pre, sizeof(hamt_ch_buf_t));
			for (int i=0; i<1<<p->h->pre; i++){
				s->lnbuf[i].a = (uint64_t*)calloc(15000, sizeof(uint64_t));
				assert(s->lnbuf[i].a);
				s->lnbuf[i].n = 0;
				s->lnbuf[i].m = 15000;
			}
		}
		if (p->rs_in && (p->rs_in->is_all_in_mem)){
			while (*p->n_seq < p->rs_in->total_reads) {
				int l;
				recover_UC_Read(&p->ucr, p->rs_in, *p->n_seq);
				l = p->ucr.length;
				if (s->n_seq == s->m_seq) {
					s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
					REALLOC(s->len, s->m_seq);
					REALLOC(s->seq, s->m_seq);
				}
				MALLOC(s->seq[s->n_seq], l);
				memcpy(s->seq[s->n_seq], p->ucr.seq, l);
				s->len[s->n_seq++] = l;
				*p->n_seq+=1;
				s->sum_len += l;
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
				if (s->n_seq>=2000 || s->sum_len >= p->opt->chunk_size){  // 2000 becuase linear buffer length
					break;
				}
			}
		}else{
			while ((ret = kseq_read(p->ks)) >= 0) {
				int l = p->ks->seq.l;
				if (*p->n_seq >= (1<<28)) {
					fprintf(stderr, "ERROR: this implementation supports no more than %d reads\n", 1<<28);
					exit(1);
				}
				if (p->rs_out) {
					if (p->rs_in->is_has_nothing){
						assert(*p->n_seq == p->rs_out->total_reads);
						ha_insert_read_len(p->rs_out, l, p->ks->name.l);
					}else{
						assert(p->rs_in->is_has_lengths && (!p->rs_in->is_all_in_mem));
						assert(l == (int)p->rs_out->read_length[*p->n_seq]);
						int i, n_N;
						for (i = n_N = 0; i < l; ++i) // count number of ambiguous bases
							if (seq_nt4_table[(uint8_t)p->ks->seq.s[i]] >= 4)
								++n_N;
						ha_compress_base(Get_READ(*p->rs_out, *p->n_seq), p->ks->seq.s, l, &p->rs_out->N_site[*p->n_seq], n_N);
						memcpy(&p->rs_out->name[p->rs_out->name_index[*p->n_seq]], p->ks->name.s, p->ks->name.l);
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
				*p->n_seq+=1;
				s->sum_len += l;
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
				if (s->n_seq>=2000 || s->sum_len >= p->opt->chunk_size){  // 2000 because linear buffer length
					break;
				}
			}
		}

		p->accumulated_time_step1 += Get_T() - t_profiling;

		if (s->sum_len == 0) {
			free(s); 
		}
		else return s;

	}else if (step==1){  // step2: (round 0:)get kmer profile of each read && mark their status (mask_readtype and mask_readnorm) || (round 1:) count kmers into linear buffers, then insert into hashtables

		plmt_step_t *s = (plmt_step_t*)in;
		int nb_cpu = p->opt->n_thread<=1? 1 : (p->opt->n_thread>KTPIPE_NB_CPU? KTPIPE_NB_CPU : p->opt->n_thread);
		#if 0
		for (int idx_seq=0; idx_seq<s->n_seq; idx_seq++){
			if (p->opt->is_HPC)
				worker_process_one_read_HPC(s, idx_seq);
			else
				worker_process_one_read_noHPC(s, idx_seq);
		}
		#endif
		pl_flt_withsorting pld;
		pld.p = p;
		pld.s = s;
		kt_for(nb_cpu, callback_worker_process_one_read_either, &pld, s->n_seq);
		if (s->p->round==0){
			for (int i=0; i<s->n_seq; i++){free(s->seq[i]);}
			free(s->len);
			free(s->seq);
			s->seq = 0; 
			s->len = 0;
			free(s);
		}else{  // insert the linear buffers
			fprintf(stderr, "warn: shouldnt hit here; legacy read selection code??\n");
			exit(1);
		}
		p->accumulated_time_step2 += Get_T() - t_profiling;
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

void hamt_mark(const hifiasm_opt_t *asm_opt, All_reads *rs, ha_ct_t *h){
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

	for (i = 0; i < asm_opt->num_reads; ++i){  // note that num_reads is the number of input files, not reads
		memset(&plmt, 0, sizeof(plmt_data_t));
		plmt.n_seq = &n_seq;
		plmt.opt = &opt;
		plmt.h = h;
		plmt.rs_in = rs;
		plmt.rs_out = rs;
		plmt.round = 0;
		plmt.hd = 0;
		plmt.asm_opt = asm_opt;
		gzFile fp = 0;
		if ((fp = gzopen(asm_opt->read_file_names[i], "r")) == 0) {
			printf("[E::%s] can't open file %s. Abort.\n", __func__, asm_opt->read_file_names[i]);
			exit(1);
		}
		plmt.ks = kseq_init(fp);
		init_UC_Read(&plmt.ucr);
		kt_pipeline(2, worker_mark_reads, &plmt, 2);

		fprintf(stderr, "[prof::%s] step1 total %.2f s, step2 %.2f s\n", __func__,
								plmt.accumulated_time_step1, plmt.accumulated_time_step2);


		if (rs->is_has_nothing){
			rs->is_has_nothing = 0;
			rs->is_has_lengths = 1;
		}else if (rs->is_has_lengths && (!rs->is_all_in_mem)){
			rs->is_all_in_mem = 1;
		}

		destory_UC_Read(&plmt.ucr);
		if (plmt.hd) ha_ct_destroy(plmt.hd);
		kseq_destroy(plmt.ks);
		gzclose(fp);
	}
}


void *hamt_ft_gen(const hifiasm_opt_t *asm_opt, All_reads *rs, uint16_t coverage, int let_reset)
{   
	// NOTE
	//     Assumes rs* has been inited (it can be blank).
	//     `coverage` is arbitrarily given, it's not calculated from hist and this function will not attempt to do so.
	// 	   `let_reset`: if set, read info except hamt read selection related status will be reset for rs
	yak_ft_t *flt_tab;
	int64_t cnt[YAK_N_COUNTS];
	int peak_hom, peak_het, cutoff;
	ha_ct_t *h;
	if (rs->is_has_nothing){
		h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_WRITE_LEN|HAMTF_FORCE_DONT_INIT|HAMTF_HAS_MARKS, NULL, NULL, rs);  // HAMTF_HAS_MARKS since the array is initialized
		rs->is_has_nothing = 0;
		rs->is_has_lengths = 1;
	}else if (rs->is_has_lengths && (!rs->is_all_in_mem)){
		h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_WRITE_SEQ|HAMTF_HAS_MARKS, NULL, NULL, rs);  // HAMTF_HAS_MARKS bc read global kmer status have been collected / the array is initialized
		rs->is_all_in_mem = 1;
	}else if (rs->is_all_in_mem){
		h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_READ|HAMTF_HAS_MARKS, NULL, NULL, rs);  // HAMTF_HAS_MARKS bc read global kmer status have been collected / the array is initialized
	}else{
		fprintf(stderr, "[E::%s] unexpected, rs status markers are wrong? Traceback:\n", __func__);
		fprintf(stderr, "[E::%s] is_has_nothing:%d, is_has_lengths:%d, is_all_in_mem:%d\n", 
						__func__, rs->is_has_nothing, rs->is_has_lengths, rs->is_all_in_mem);
		exit(1);
	}

	// histogram
	ha_ct_hist(h, cnt, asm_opt->thread_num);
	peak_hom = ha_analyze_count(YAK_N_COUNTS, asm_opt->min_hist_kmer_cnt, cnt, &peak_het);
	if (peak_hom > 0) fprintf(stderr, "[M::%s] peak_hom: %d; peak_het: %d\n", __func__, peak_hom, peak_het);
	
	if (!asm_opt->is_disable_read_selection){
		cutoff = (int)(coverage *2 * asm_opt->high_factor);
	}else{
		cutoff = (int)(coverage * asm_opt->high_factor);
	}
	if (cutoff > YAK_MAX_COUNT - 1) cutoff = YAK_MAX_COUNT - 1;
	ha_ct_shrink(h, cutoff, YAK_MAX_COUNT, asm_opt->thread_num);
	flt_tab = gen_hh(h);
	ha_ct_destroy(h);
	fprintf(stderr, "[M::%s::%.3f*%.2f@%.3fGB] ==> filtered out %ld k-mers occurring %d or more times\n", __func__,
			yak_realtime(), yak_cpu_usage(), yak_peakrss_in_gb(), (long)kh_size(flt_tab), cutoff);
	if (let_reset) reset_All_reads(rs);
	return (void*)flt_tab;
}



/////////////////////////////////////////////////////////////
////            meta debug util                            //
/////////////////////////////////////////////////////////////

void debug_printstat_read_status(All_reads *rs){
	int drp = 0, bcs_median=0, bcs_longlow=0, bcs_kmer=0, san_error = 0;
	for (int i=0; i< (int)rs->total_reads; i++){
		if (rs->mask_readnorm[i] & 1) drp++;
		else{
			if (rs->mask_readnorm[i] & HAMT_VIA_MEDIAN) bcs_median++;
			else if (rs->mask_readnorm[i] & HAMT_VIA_LONGLOW)  bcs_longlow++;
			else if (rs->mask_readnorm[i] & HAMT_VIA_KMER) bcs_kmer++;
			else if (1) san_error++;
		}
		// printf("r#%d\t%" PRIu16 "\t%" PRIu8 "\n",i, rs->median[i], rs->mask_readnorm[i]);
	}
	fprintf(stderr, "[DEBUGhamt::%s] total reads: %d, dropped %d.\n", __func__, (int)rs->total_reads, drp);
	fprintf(stderr, "[DEBUGhamt::%s] retained via median: %d, via longlow: %d, via infrequent kmers %d.\n", __func__, bcs_median, bcs_longlow, bcs_kmer);
	fflush(stderr);
	// exit(0);
}

typedef struct{
	int n, m;
	uint16_t *a;
}debugrkp_linearbuf_t;

typedef struct {  // global data structure for pipeline (pl_data_t)
	yak_copt_t *opt;
	uint64_t *n_seq;
	All_reads *rs_in;
	All_reads *rs_out;
	kseq_t *ks;
	ha_ct_t *h;  // counted kmer tables (equals to the h from 1st call of ha_count in vanilla hifiasm)
	// char *filename_out;
	FILE *fp_out;
}debugrkp_plmt_data_t;

typedef struct {  // step data structure (st_data_t)
	debugrkp_plmt_data_t* p;
	uint64_t n_seq0;
	int n_seq, m_seq, sum_len, nk;
	int *len;
	char **seq;
	char **qname; int *qnamelen;
	debugrkp_linearbuf_t *lnbuf;  // $n_seq linear buffers
}debugrkp_plmt_step_t;


static void worker_get_one_read_kmer_profile_noHPC(void* data, long idx_for, int tid) // insert k-mers in $seq to linear buffer $buf
{
	debugrkp_plmt_step_t *s = (debugrkp_plmt_step_t*) data;
	int i, l;
	int k = s->p->opt->k;
	char *seq = s->seq[idx_for];
	int len = s->len[idx_for];
	ha_ct_t *h = s->p->h;
	debugrkp_linearbuf_t *lnbuf = &s->lnbuf[idx_for];

	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	uint64_t hash;
	uint16_t count;
	khint_t key;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			x[0] = (x[0] << 1 | (c&1))  & mask;
			x[1] = (x[1] << 1 | (c>>1)) & mask;
			x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
			x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
			if (++l >= k){
					if ((lnbuf->n+5)>=(lnbuf->m)){
						lnbuf->m += lnbuf->m>>1;
						REALLOC(lnbuf->a, lnbuf->m);
					}
					hash = yak_hash_long(x);
					yak_ct_t *h_ = h->h[hash&((1<<YAK_COUNTER_BITS)-1)].h;
					key = yak_ct_get(h_, hash>>YAK_COUNTER_BITS<<YAK_COUNTER_BITS1);
					if (key==kh_end(h_)) {  // not in the hashtable; because of bloom filter, count is 1
						count = 1;
					}else{
						count = kh_key(h_, key) & YAK_MAX_COUNT;
						assert(count>=2);
					}
					lnbuf->a[lnbuf->n] = count;
					lnbuf->n++;
				}
		} else l = 0, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
}

static void worker_get_one_read_kmer_profile_HPC(void* data, long idx_for, int tid) // insert k-mers in $seq to linear buffer $buf
{// ch_buf_t *buf, int k, int p, int len, const char *seq
	debugrkp_plmt_step_t *s = (debugrkp_plmt_step_t*) data;
	int i, l, last = -1;
	int k = s->p->opt->k;
	char *seq = s->seq[idx_for];
	int len = s->len[idx_for];
	ha_ct_t *h = s->p->h;
	debugrkp_linearbuf_t *lnbuf = &s->lnbuf[idx_for];
	khint_t key;

	uint64_t x[4], mask = (1ULL<<k) - 1, shift = k - 1;
	uint64_t hash;
	uint16_t count;
	for (i = l = 0, x[0] = x[1] = x[2] = x[3] = 0; i < len; ++i) {
		int c = seq_nt4_table[(uint8_t)seq[i]];
		if (c < 4) { // not an "N" base
			if (c != last) {
				x[0] = (x[0] << 1 | (c&1))  & mask;
				x[1] = (x[1] << 1 | (c>>1)) & mask;
				x[2] = x[2] >> 1 | (uint64_t)(1 - (c&1))  << shift;
				x[3] = x[3] >> 1 | (uint64_t)(1 - (c>>1)) << shift;
				if (++l >= k){
					if ((lnbuf->n+5)>=(lnbuf->m)){
						lnbuf->m += lnbuf->m>>1;
						REALLOC(lnbuf->a, lnbuf->m);
					}
					hash = yak_hash_long(x);
					yak_ct_t *h_ = h->h[hash&((1<<YAK_COUNTER_BITS)-1)].h;
					key = yak_ct_get(h_, hash>>YAK_COUNTER_BITS<<YAK_COUNTER_BITS1);
					if (key==kh_end(h_)) {  // not in the hashtable; because of bloom filter, count is 1
						count = 1;
					}else{
						count = kh_key(h_, key) & YAK_MAX_COUNT;
						assert(count>=2);
					}
					lnbuf->a[lnbuf->n] = count;
					lnbuf->n++;
				}
				last = c;
			}
		} else l = 0, last = -1, x[0] = x[1] = x[2] = x[3] = 0; // if there is an "N", restart
	}
	// printf("sequence length: %d, collected %d.\n", len, lnbuf->n);
}

static void *worker_read_kmer_profile(void *data, int step, void *in){  // callback of kt_pipeline
	debugrkp_plmt_data_t *p = (debugrkp_plmt_data_t*)data;
	if (step==0){ // step 1: read a block of sequences
		int ret;
		debugrkp_plmt_step_t *s;
		CALLOC(s, 1);
		s->p = p;
		s->n_seq0 = *p->n_seq;
		s->lnbuf = 0;

		s->lnbuf = (debugrkp_linearbuf_t*)calloc(1500, sizeof(debugrkp_linearbuf_t));
		for (int i=0; i<1500; i++){
			s->lnbuf[i].a = (uint16_t*)calloc(15000, sizeof(uint16_t));
			s->lnbuf[i].n = 0;
			s->lnbuf[i].m = 15000;
		}

		while ((ret = kseq_read(p->ks)) >= 0) {
			int l = p->ks->seq.l;
			if (*p->n_seq >= 1<<28) {
				fprintf(stderr, "WARNING: too many reads, will continue but mind that the data won't fit in hifiasm\n");
			}
			if (p->rs_out) {  // ignoring flags
				assert(*p->n_seq == p->rs_out->total_reads);
				ha_insert_read_len(p->rs_out, l, p->ks->name.l);
			}
			if (s->n_seq == s->m_seq) {
				s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
				REALLOC(s->len, s->m_seq);
				REALLOC(s->qnamelen, s->m_seq);
				REALLOC(s->seq, s->m_seq);
				REALLOC(s->qname, s->m_seq);
			}
			MALLOC(s->seq[s->n_seq], l);
			memcpy(s->seq[s->n_seq], p->ks->seq.s, l);
			MALLOC(s->qname[s->n_seq], p->ks->name.l+5);
			memcpy(s->qname[s->n_seq], p->ks->name.s, p->ks->name.l);
			s->qname[s->n_seq][p->ks->name.l] = '\0';
			s->qnamelen[s->n_seq] = p->ks->name.l;
			s->len[s->n_seq] = l;
			s->n_seq++;
			*p->n_seq+=1;
			s->sum_len += l;
			s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;

			// if (s->sum_len >= p->opt->chunk_size){
			if (s->n_seq==1500){
				break;
			}
			if (s->n_seq>1500){fprintf(stderr, "?\n"); fflush(stderr); exit(1);}
		}

		if (s->sum_len == 0) {
			free(s->lnbuf);
			for (int i=0; i<s->n_seq; i++){free(s->seq[i]);}
			for (int i=0; i<s->n_seq; i++){free(s->qname[i]);}
			free(s->len); free(s->qnamelen);
			free(s->seq); free(s->qname);
			free(s);
		}
		else return s;

	}else if (step==1){  // step2 get profile
		debugrkp_plmt_step_t *s = (debugrkp_plmt_step_t*)in;
		if (p->opt->is_HPC)
			kt_for(s->p->opt->n_thread, worker_get_one_read_kmer_profile_HPC, s, s->n_seq);
		else
			kt_for(s->p->opt->n_thread, worker_get_one_read_kmer_profile_noHPC, s, s->n_seq);
		return s;
	}else if(step==2){  // write to disk
		debugrkp_plmt_step_t *s = (debugrkp_plmt_step_t*)in;
		int i;
		for (int j=0; j<s->n_seq; j++){
			fprintf(s->p->fp_out, "%s\t", s->qname[j]);
			for (i=0; i<s->lnbuf[j].n; i++){
				fprintf(s->p->fp_out, "%" PRIu16 ",", s->lnbuf[j].a[i]);
			}
			fprintf(s->p->fp_out, "\n");
			free(s->lnbuf[j].a);
		}
		free(s->lnbuf);
		for (int i=0; i<s->n_seq; i++){free(s->seq[i]);}
		for (int i=0; i<s->n_seq; i++){free(s->qname[i]);}
		free(s->len); free(s->qnamelen);
		free(s->seq); free(s->qname);
		free(s);
	}
	return 0;
}


int hamt_read_kmer_profile(hifiasm_opt_t *asm_opt, All_reads *rs){  // debug
	ha_ct_t *h;
	h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_WRITE_LEN, NULL, NULL, rs);
	reset_All_reads(rs);

	debugrkp_plmt_data_t pl;
	gzFile fp = 0;
	pl.h = h;
	pl.rs_in = rs;
	pl.rs_out = rs;
	char *filename_out = (char*)malloc(500);  // just in case there's a super long path name

	yak_copt_t opt;
	yak_copt_init(&opt);
	opt.k = asm_opt->k_mer_length;
	opt.is_HPC = !(asm_opt->flag&HA_F_NO_HPC);
	opt.n_thread = asm_opt->thread_num;
	pl.opt = &opt;
	uint64_t n_seq = 0;
	pl.n_seq = &n_seq;

	sprintf(filename_out, "%s.perread_kmerprofile", asm_opt->output_file_name);
	pl.fp_out = fopen(filename_out, "w");
	for (int i=0; i<asm_opt->num_reads; i++){
		if ((fp = gzopen(asm_opt->read_file_names[i], "r")) == 0) return 0;
		pl.ks = kseq_init(fp);
		kt_pipeline(3, worker_read_kmer_profile, &pl, 3);	
		kseq_destroy(pl.ks);
	}
	fflush(pl.fp_out);
	fclose(pl.fp_out);	
	free(filename_out);
	ha_ct_destroy(h);
	return 0;
}

int hamt_printout_ha_count(hifiasm_opt_t *asm_opt, All_reads *rs){  // debug
	ha_ct_t *h;
	h = ha_count(asm_opt, HAF_COUNT_ALL, NULL, NULL, rs);
	khint_t key;
	char filename_out[500];
	sprintf(filename_out, "%s.ha_count", asm_opt->output_file_name);
	FILE *fp_out = fopen(filename_out, "w");
	if (!fp_out) {fprintf(stderr, "ERROR: specify output filename.\n"); exit(1);}

	uint64_t hash;
	uint16_t count;
	fprintf(stderr, "finished ha_count, writing to disk...\n");

	for (int i=0; i < (1<<h->pre); ++i){
		yak_ct_t *g = h->h[i].h;
		for (key=0; key<kh_end(g); key++){
			if (kh_exist(g, key)){
				hash = kh_key(g, key);
				count = (uint16_t) (hash&YAK_MAX_COUNT);
				if (count<5) continue;
				hash = (uint64_t) (hash>>YAK_COUNTER_BITS1<<YAK_COUNTER_BITS | i);
				fprintf(fp_out, "%" PRIu64 "\t%" PRIu16 "\n", hash, count);
			}
		}
	}
	fflush(fp_out);
	fclose(fp_out);
	ha_ct_destroy(h);
	return 0;
}


///////////////////////////////
//       pre read selection  //
///////////////////////////////

void exp_load_raw_reads(const hifiasm_opt_t *asm_opt, All_reads *rs, int has_len){
	// load all reads into mem (effectively the method used by worker_count with HAF_RS_WRITE_SEQ)
	// exp for experiment. Refactor this into the pipeline if read selection really needs to rely on sorting!
	gzFile fp = 0;
	kseq_t *ks;
	uint32_t idx_seq = 0;
	int ret;  // ks
	if (!has_len) {
		reset_All_reads(rs);

		for (int idx_file =0; idx_file<asm_opt->num_reads; idx_file++){  // 1st pass
			fp = gzopen(asm_opt->read_file_names[idx_file], "r");
			if (fp==0) {fprintf(stderr, "E::%s can't open input file.\n", __func__); exit(1);}
			ks = kseq_init(fp);
			while ((ret = kseq_read(ks))>=0){
				if (rs->total_reads>=(1<<28)){fprintf(stderr, "E::%s too many reads.\n", __func__); exit(1);}
				ha_insert_read_len(rs, ks->seq.l, ks->name.l);
			}
			kseq_destroy(ks); gzclose(fp); 
		}
	}
	
	malloc_All_reads(rs);

	for(int idx_file =0; idx_file<asm_opt->num_reads; idx_file++){  // 2nd pass: get sequences
		int i, n_N, l;
		fp = gzopen(asm_opt->read_file_names[idx_file], "r");
		ks = kseq_init(fp);
		while ((ret = kseq_read(ks))>=0){
			l = ks->seq.l;
			assert(l == (int)rs->read_length[idx_seq]);
			for (i = n_N = 0; i < l; ++i) // count number of ambiguous bases
				if (seq_nt4_table[(uint8_t)ks->seq.s[i]] >= 4)
					++n_N;
			ha_compress_base(Get_READ(*rs, idx_seq), ks->seq.s, l, &rs->N_site[idx_seq], n_N);
			memcpy(&rs->name[rs->name_index[idx_seq]], ks->name.s, ks->name.l);
			idx_seq++;
		}
		kseq_destroy(ks); gzclose(fp); 
	}
	assert(rs->total_reads==idx_seq);
	fprintf(stderr, "M::%s loaded %d sequences.\n", __func__, idx_seq); fflush(stderr);
}

static void *worker_mark_reads_withsorting(void *data, int step, void *in){  // callback of kt_pipeline
    // interate over reads by the specified order, instead using the ordering from input file(s)
	plmt_data_t *p = (plmt_data_t*)data;
	if (step==0){ // collect reads
		plmt_step_t *s;
		// prepare s
		CALLOC(s, 1);
		s->p = p;
		s->n_seq0 = *p->n_seq;
		s->lnbuf = 0;
		s->lnbuf = (hamt_ch_buf_t*)calloc(1<<p->h->pre, sizeof(hamt_ch_buf_t));
		for (int i=0; i<1<<p->h->pre; i++){
			s->lnbuf[i].a = (uint64_t*)calloc(15000, sizeof(uint64_t));
			s->lnbuf[i].n = 0;
			s->lnbuf[i].m = 15000;
			}
		s->RIDs = (uint64_t*)malloc(sizeof(uint64_t)*2000);
		
		// recruit reads
		if (p->rs_in->is_all_in_mem){
			uint64_t readID;
			while (s->n_seq<2000 && ((*p->n_seq) < p->rs_in->total_reads)){
				int l;
				readID = *p->n_seq;  // use the order in the input file
				// if (p->asm_opt->readselection_sort_order==0){
				// 	readID = *p->n_seq;  // use the order in the input file
				// }else if (p->asm_opt->readselection_sort_order==1){  // use sorting order, smallest first
				// 	readID = (uint64_t) ((uint32_t)p->rs_in->statpack[*p->n_seq]);  // touch rarer reads first
				// }else if (p->asm_opt->readselection_sort_order==2){  // use sorting order, largest first
				// 	readID = (uint64_t) ((uint32_t)p->rs_in->statpack[p->rs_in->total_reads - 1 - *p->n_seq]);  // touch more prevalent reads first
				// }else{
				// 	fprintf(stderr, "error at readID, unexpected asm_opt switch?\n"); fflush(stderr);
				// 	exit(1);
				// }
				s->RIDs[s->n_seq] = readID;  // since we can't rely on n_seq0+idx_for
				recover_UC_Read(&p->ucr, p->rs_in, readID);
				l = p->ucr.length;
				if (s->n_seq == s->m_seq) {
						s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
						REALLOC(s->len, s->m_seq);
						REALLOC(s->seq, s->m_seq);
				}
				MALLOC(s->seq[s->n_seq], l);
				memcpy(s->seq[s->n_seq], p->ucr.seq, l);
				s->len[s->n_seq++] = l;
				*p->n_seq = *p->n_seq + 1;
				s->sum_len += l;
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
			}
		}else{
			int ret;
			while ((ret = kseq_read(p->ks)) >= 0) {
				int l = p->ks->seq.l;
				if (*p->n_seq >= 1<<28) {
					fprintf(stderr, "ERROR: this implementation supports no more than %d reads\n", 1<<28);
					exit(1);
				}
				if (p->rs_out) {
					if (p->rs_in->is_has_nothing){
						assert(*p->n_seq == p->rs_out->total_reads);
						ha_insert_read_len(p->rs_out, l, p->ks->name.l);
					} else if (p->flag & HAF_RS_WRITE_SEQ) {
						assert(p->rs_in->is_has_lengths && (!p->rs_in->is_all_in_mem));
						assert(l == (int)p->rs_out->read_length[*p->n_seq]);
						int i, n_N;
						for (i = n_N = 0; i < l; ++i) // count number of ambiguous bases
							if (seq_nt4_table[(uint8_t)p->ks->seq.s[i]] >= 4)
								++n_N;
						ha_compress_base(Get_READ(*p->rs_out, *p->n_seq), p->ks->seq.s, l, &p->rs_out->N_site[*p->n_seq], n_N);
						memcpy(&p->rs_out->name[p->rs_out->name_index[*p->n_seq]], p->ks->name.s, p->ks->name.l);
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
				*p->n_seq = *p->n_seq + 1;
				s->sum_len += l;
				s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
				if (s->n_seq>=2000 || ((*p->n_seq) >= p->rs_in->total_reads))
					break;
			}
		}

		if (s->sum_len == 0) {
			for (int i=0; i<1<<p->h->pre; i++){
				free(s->lnbuf[i].a);
			}
			free(s->lnbuf);
			free(s->RIDs);
			free(s);
		}
		else return s;
	} else if (step==1){  // process reads
		plmt_step_t *s = (plmt_step_t*)in;
		// note: do not use kt_for here, or it's race (that *doesn't always* but occassionally fail)! 
		for (int idx_seq=0; idx_seq<s->n_seq; idx_seq++){
			if (p->opt->is_HPC)
				worker_process_one_read_HPC(s, idx_seq);
			else
				worker_process_one_read_noHPC(s, idx_seq);
		}
		int san = 0;
		for (int i=0; i<1<<p->h->pre; i++){
			san+=s->lnbuf[i].n;
		}
		// update runtime hashtable
		assert(s->p->hd);
		kt_for(s->p->opt->n_thread, worker_insert_lnbuf, s, 1<<p->h->pre);
		// clean up
		for (int i=0; i<s->n_seq; i++){free(s->seq[i]);}
		free(s->len); s->len = 0;
		free(s->seq); s->seq = 0; 
		for (int i=0; i<1<<p->h->pre; i++){
			free(s->lnbuf[i].a);
			}				
		free(s->lnbuf);
		free(s->RIDs);
		free(s);
	} 
	return 0;
}

static void *worker_markinclude_lowq_reads(void *data, int step, void *in){  // callback of kt_pipeline
	// FUNC
	//     This the preovec read selection's worker function.
	//     Any read that has many rare kmers will be retained.
	//      (measurement: lower quantile of global kmer frequency)
	plmt_data_t *p = (plmt_data_t*)data;
	double t_profiling;
	if (step==0){ // collect reads
		t_profiling = Get_T();
		plmt_step_t *s;
		// prepare s
		CALLOC(s, 1);
		s->p = p;
		s->n_seq0 = *p->n_seq;
		s->lnbuf = 0;
		s->lnbuf = (hamt_ch_buf_t*)calloc(1<<p->hd->pre, sizeof(hamt_ch_buf_t));
		for (int i=0; i<1<<p->hd->pre; i++){
			s->lnbuf[i].a = (uint64_t*)calloc(15000, sizeof(uint64_t));
			s->lnbuf[i].n = 0;
			s->lnbuf[i].m = 15000;
			}
		s->RIDs = (uint64_t*)malloc(sizeof(uint64_t)*2000);
		// recruit reads
		uint64_t readID;
		while (s->n_seq<2000 && ((*p->n_seq) < p->rs_in->total_reads)){
			int l;
			readID = (uint64_t) ((uint32_t)p->rs_in->statpack[*p->n_seq]);

			if ((p->rs_in->mask_readnorm[readID] & 1)==0){  // read already retained
				*p->n_seq = *p->n_seq + 1;
				continue;
			}

			s->RIDs[s->n_seq] = readID;  // since we can't rely on n_seq0+idx_for
			recover_UC_Read(&p->ucr, p->rs_in, readID);
			l = p->ucr.length;
			if (s->n_seq == s->m_seq) {
				s->m_seq = s->m_seq < 16? 16 : s->m_seq + (s->m_seq>>1);
				REALLOC(s->len, s->m_seq);
				REALLOC(s->seq, s->m_seq);
			}
			MALLOC(s->seq[s->n_seq], l);
			memcpy(s->seq[s->n_seq], p->ucr.seq, l);
			s->len[s->n_seq++] = l;
			*p->n_seq = *p->n_seq + 1;
			s->sum_len += l;
			s->nk += l >= p->opt->k? l - p->opt->k + 1 : 0;
		}
		p->accumulated_time_step1 += Get_T() - t_profiling;
		if (s->sum_len == 0) {
			for (int i=0; i<1<<p->hd->pre; i++){
				free(s->lnbuf[i].a);
			}
			free(s->lnbuf);
			free(s->RIDs);
			free(s);
		}
		else {return s;}
	} else if (step==1){  // process reads
		// what does this block do:
		//     - retain reads based on runtime kmer frequency
		//     - put kmers into per-thread linear buffers (does not update hastable yet)
		// r032 bug fix note:
		//     Must insert linear buffers to the runtime hashtable within this step,
		//       otherwise the process is not stable (hashtable gets updated while we check lowq values).
		//     This effect was not prominent if we process read with a single thread (before r024),
		//       however it is visible when threaded. 
		t_profiling = Get_T();
		plmt_step_t *s = (plmt_step_t*)in;
		int nb_cpu = s->p->opt->n_thread<=1 ? 1 : (s->p->opt->n_thread>KTPIPE_NB_CPU? KTPIPE_NB_CPU : s->p->opt->n_thread);

		// allocate bufers
		s->threaded_lnbuf = (hamt_ch_buf_t**)calloc(nb_cpu, sizeof(hamt_ch_buf_t*));
		int n_pre = (1<<s->p->hd->pre);
		for (int i_cpu=0; i_cpu<nb_cpu; i_cpu++){
			s->threaded_lnbuf[i_cpu] = (hamt_ch_buf_t*)calloc(n_pre, sizeof(hamt_ch_buf_t));
			for (int i=0; i<n_pre; i++){
				s->threaded_lnbuf[i_cpu][i].a = (uint64_t*)malloc(16*sizeof(uint64_t));
				s->threaded_lnbuf[i_cpu][i].m = 16;
			}
		}
		// process reads
		kt_for(nb_cpu, callback_worker_process_one_read_inclusive2, s, s->n_seq);
		int san = 0;
		for (int i_cpu=0; i_cpu<nb_cpu; i_cpu++){
			for (int i=0; i<n_pre; i++){
				san+= s->threaded_lnbuf[i_cpu][i].n;
			}
		}

		#if 0
		// (single thread approach)
		// note: do not use kt_for here, or it's race
		for (int idx_seq=0; idx_seq<s->n_seq; idx_seq++){
			if (p->opt->is_HPC)
				// worker_process_one_read_HPC_inclusive(s, idx_seq, s->p->round);
				worker_process_one_read_inclusive(s, idx_seq, s->p->round, 1);
			else
				worker_process_one_read_inclusive(s, idx_seq, s->p->round, 0);
		}
		int san = 0;
		for (int i=0; i<1<<p->hd->pre; i++){
			san+=s->lnbuf[i].n;
		}
		#endif

		p->accumulated_time_step2 += Get_T() - t_profiling;

		if (san>0){
			// insert linear buffer into the runtime hashtable
			for (int i_cpu=0; i_cpu<nb_cpu; i_cpu++){
				s->lnbuf = s->threaded_lnbuf[i_cpu];
				kt_for(s->p->opt->n_thread, worker_insert_lnbuf, s, 1<<p->hd->pre);
			}
			p->accumulated_time_step3 += Get_T() - t_profiling;
		}

		// clean up 
		for (int i=0; i<s->n_seq; i++){free(s->seq[i]);}
		free(s->len); s->len = 0;
		free(s->seq); s->seq = 0; 
		for (int i_cpu=0; i_cpu<nb_cpu; i_cpu++){
			for (int i=0; i<n_pre; i++){
				free(s->threaded_lnbuf[i_cpu][i].a);
			}
			free(s->threaded_lnbuf[i_cpu]);
		}
		free(s->threaded_lnbuf);
		#if 0
		for (int i=0; i<1<<p->hd->pre; i++){
			free(s->lnbuf[i].a);
			}				
		free(s->lnbuf);
		#endif
		free(s->RIDs);
		free(s);
	}
	return 0;
}

void hamt_flt_withsorting(const hifiasm_opt_t *asm_opt, All_reads *rs){
	int64_t cnt[YAK_N_COUNTS];
	ha_ct_t *h;
	int peak_het;  // placeholder
	h = ha_count(asm_opt, HAF_COUNT_ALL|HAF_RS_WRITE_LEN, NULL, NULL, rs); 
	ha_ct_hist(h, cnt, asm_opt->thread_num);
	ha_analyze_count(YAK_N_COUNTS, asm_opt->min_hist_kmer_cnt, cnt, &peak_het);

	hamt_mark(asm_opt, rs, h);  // collects mean/median/std and mark global category of each read
	// exp_load_raw_reads(asm_opt, rs, 1);  // already has seq lens, load sequences

	/////// sorting ////////
	// pack: lowFrequency | median | readID
	rs->statpack = (uint64_t*)malloc(sizeof(uint64_t) * rs->total_reads);
	assert(rs->statpack);
	for (uint32_t i=0; i<rs->total_reads; i++){
		// using preovec read selection
		rs->statpack[i] = ((uint64_t) ((YAK_MAX_COUNT -1 - ((uint64_t)rs->lowq[i]))<<48)) | (((uint64_t) rs->median[i])<<32) | i;
	}
	radix_sort_hamt64(rs->statpack, rs->statpack+rs->total_reads);
	fprintf(stderr, "[M::%s] Reads sorted. \n", __func__);
	
	memset(rs->mask_readnorm, 0, rs->total_reads);
	ha_ct_destroy(h);
	
}

void hamt_flt_withsorting_supervised(const hifiasm_opt_t *asm_opt, All_reads *rs, int nb_to_keep){
	// called by hamt_pre_ovec_v2
	// all stats have been collected before this (global median etc + sorting by `hamt_flt_withsorting`)

	memset(R_INF.mask_readnorm, 1, R_INF.total_reads);  // we will be including reads, not excluding
	double t_profiling = Get_T();

	uint64_t n_seq = 0;
	yak_copt_t opt;
	yak_copt_init(&opt);
	opt.k = asm_opt->k_mer_length;
	opt.is_HPC = !(asm_opt->flag&HA_F_NO_HPC);
	opt.n_thread = asm_opt->thread_num;
	plmt_data_t plmt;
	memset(&plmt, 0, sizeof(plmt_data_t));
	plmt.n_seq = &n_seq;
	plmt.opt = &opt;
	plmt.rs_in = rs;
	plmt.rs_out = rs;
	plmt.h = ha_ct_init(plmt.opt->k, plmt.opt->pre, plmt.opt->bf_n_hash, plmt.opt->bf_shift);  // it's a placeholder since some functions may expect h (only to use the pre)
	plmt.hd = ha_ct_init(plmt.opt->k, plmt.opt->pre, plmt.opt->bf_n_hash, plmt.opt->bf_shift);
	plmt.asm_opt = asm_opt;
	plmt.nb_reads_limit = nb_to_keep;
	init_UC_Read(&plmt.ucr);

	// // 1st pass: keep all reads with infrequent kmers and median less than 100
	// plmt.round = 0;
	// plmt.nb_reads_kept = 0;
	// fprintf(stderr, "[M::%s] entering round#1...\n", __func__);
	// kt_pipeline(asm_opt->thread_num, worker_markinclude_lowq_reads, &plmt, 2);
	// fprintf(stderr, "[M::%s] finished round#1, retaining %d reads.\n",__func__, plmt.nb_reads_kept);
	// fprintf(stderr, "[prof::%s] used %.2f s\n", __func__, Get_T() - t_profiling);
	// t_profiling = Get_T();
	plmt.nb_reads_kept = 0;  // just go for runtime lowq selection.

	// 2nd pass: include reads until reaching the limit
	// repack the stat
	if (plmt.nb_reads_kept<plmt.nb_reads_limit){
		plmt.runtime_median_threshold = 50; //  deprecated
		for (uint64_t i=0; i<rs->total_reads; i++){
			rs->statpack[i] = ((uint64_t) rs->median[i]<<48) | (((uint64_t) (YAK_MAX_COUNT -1 -rs->std[i] + .499))<<32) | i;  // std won't exceed 14 bits, therefore safe to cast
		}
		radix_sort_hamt64(rs->statpack, rs->statpack+rs->total_reads);
		plmt.round = 1;
		int rounds_of_dropping = 0, grace=0;
		int prv_nb_reads = plmt.nb_reads_kept;
		{
			fprintf(stderr, "[M::%s] entering round#2, pass#%d...\n", __func__, rounds_of_dropping);
			n_seq = 0;  // reset!
			kt_pipeline(asm_opt->thread_num, worker_markinclude_lowq_reads, &plmt, 2);
			fprintf(stderr, "[prof::%s] step1 %.2f s, step2 %.2f s, step3 %.2f s\n", __func__, 
									plmt.accumulated_time_step1, plmt.accumulated_time_step2, plmt.accumulated_time_step3);
			fprintf(stderr, "[prof::%s] used %.2f s\n", __func__, Get_T() - t_profiling);
			t_profiling = Get_T();

			plmt.runtime_median_threshold +=50;  // deprecated
			rounds_of_dropping++;
		
			fprintf(stderr, "[M::%s] finished selection, retained %d reads (out of %d).\n", __func__, plmt.nb_reads_kept, plmt.nb_reads_limit);	
			prv_nb_reads = plmt.nb_reads_kept;
		
		}
	}else{
		fprintf(stderr, "[M/W::%s] round#1 has collected enough reads, ignoreing lower quantile criteria.\n", __func__);
	}

	destory_UC_Read(&plmt.ucr);
	if (plmt.hd) ha_ct_destroy(plmt.hd);
	if (plmt.h) ha_ct_destroy(plmt.h);
	
}


void hamt_flt_no_read_selection(hifiasm_opt_t *asm_opt, All_reads *rs){
	assert(asm_opt->is_disable_read_selection);
	init_All_reads(rs);
	// exp_load_raw_reads(asm_opt, rs, 0);  // load sequences

	memset(rs->mask_readnorm, 0, rs->hamt_stat_buf_size * sizeof(uint8_t));  
	memset(rs->mask_readtype, 0, rs->hamt_stat_buf_size * sizeof(uint8_t));  
	memset(rs->mean, 0, rs->hamt_stat_buf_size * sizeof(double));  
	memset(rs->median, 0, rs->hamt_stat_buf_size * sizeof(uint16_t));  
	memset(rs->std, 0, rs->hamt_stat_buf_size * sizeof(double));  
}

void hamt_flt_no_read_selection_from_disk_sancheck(hifiasm_opt_t *asm_opt, All_reads *rs){
	// FUNC
	//     Prevent the following case: Command line specifies no read selection, 
	//                                 but bin file did read selection.
	//     Abort if so.
	if (!asm_opt->is_disable_read_selection){return;}
	for (uint64_t i=0; i<rs->total_reads; i++){
		if (rs->mask_readnorm[i] & 1){
			fprintf(stderr, "[E::%s] set to keep all reads, but in bin file read #%d was discarded. Aborting.\n", __func__, (int)i);
			fflush(stderr);
			exit(1);
		}
	}
}

int write_pt_index(void *flt_tab, ha_pt_t *ha_idx, All_reads* r, hifiasm_opt_t* opt, char* file_name)
{
    char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.pt_flt", file_name);
    FILE* fp = fopen(gfa_name, "w");
	if (!fp) {
		free(gfa_name);
        return 0;
    }
	yak_ft_t *ha_flt_tab = (yak_ft_t*)flt_tab;

	if(ha_flt_tab)
	{
		fwrite("f", 1, 1, fp);
		yak_ft_save(ha_flt_tab, fp);
	}


	if(ha_idx)
	{
		int i;
		ha_pt1_t *g;
		fwrite("h", 1, 1, fp);
		fwrite(&ha_idx->k, sizeof(ha_idx->k), 1, fp);
		fwrite(&ha_idx->pre, sizeof(ha_idx->pre), 1, fp);
		fwrite(&ha_idx->tot, sizeof(ha_idx->tot), 1, fp);
		fwrite(&ha_idx->tot_pos, sizeof(ha_idx->tot_pos), 1, fp);

		for (i = 0; i < 1<<ha_idx->pre; ++i) 
		{
			g = &(ha_idx->h[i]);
			yak_pt_save(g->h, fp);
			fwrite(&g->n, sizeof(g->n), 1, fp);
			fwrite(g->a, sizeof(ha_idxpos_t), g->n, fp);
		}
	}

	fwrite(&opt->number_of_round, sizeof(opt->number_of_round), 1, fp);
	fwrite(&opt->hom_cov, sizeof(opt->hom_cov), 1, fp);
	fwrite(&opt->het_cov, sizeof(opt->het_cov), 1, fp);
	fwrite(&opt->max_n_chain, sizeof(opt->max_n_chain), 1, fp);


	write_All_reads(r, gfa_name);

	fprintf(stderr, "[M::%s] Index has been written.\n", __func__);
    free(gfa_name);
	fclose(fp);
	return 1;
}

int load_pt_index(void **r_flt_tab, ha_pt_t **r_ha_idx, All_reads* r, hifiasm_opt_t* opt, char* file_name)
{
	char* gfa_name = (char*)malloc(strlen(file_name)+25);
    sprintf(gfa_name, "%s.pt_flt", file_name);
    FILE* fp = fopen(gfa_name, "r");
	if (!fp) {
		free(gfa_name);
        return 0;
    }

	ha_pt_t *ha_idx = NULL;
	char mode = 0;
	int f_flag, absent, i;
	double index_time, index_s_time, pos_time, pos_s_time;

	
	
	f_flag += fread(&mode, 1, 1, fp);
	if(mode == 'f')
	{
		index_time = yak_realtime();

		yak_ft_load((yak_ft_t **)r_flt_tab, fp);
		
	    f_flag += fread(&mode, 1, 1, fp);

		fprintf(stderr, "[M::%s::%.3f] ==> Loaded flt table\n", __func__, yak_realtime()-index_time);
	}
	///insert using multiple threads???
	if(mode == 'h')
	{
		pos_time = index_time = 0;

		CALLOC(ha_idx, 1);
		ha_pt1_t *g;
		f_flag += fread(&ha_idx->k, sizeof(ha_idx->k), 1, fp);
		f_flag += fread(&ha_idx->pre, sizeof(ha_idx->pre), 1, fp);
		f_flag += fread(&ha_idx->tot, sizeof(ha_idx->tot), 1, fp);
		f_flag += fread(&ha_idx->tot_pos, sizeof(ha_idx->tot_pos), 1, fp);
		CALLOC(ha_idx->h, 1<<ha_idx->pre);
		for (i = 0; i < 1<<ha_idx->pre; ++i) 
		{
			index_s_time = yak_realtime();

			g = &(ha_idx->h[i]);
			yak_pt_load(&(g->h), fp);

			index_time += yak_realtime() - index_s_time;
			
			pos_s_time = yak_realtime();

			f_flag += fread(&g->n, sizeof(g->n), 1, fp);
			MALLOC(g->a, g->n);
			f_flag += fread(g->a, sizeof(ha_idxpos_t), g->n, fp);

			pos_time += yak_realtime() - pos_s_time;
		}
		(*r_ha_idx) = ha_idx;

		fprintf(stderr, "[M::%s::%.3f(index)/%.3f(pos)] ==> Loaded pos table\n", __func__, index_time, pos_time);
	}

	if(mode != 'h' && mode != 'f')
	{
		free(gfa_name);
		fclose(fp);
		return 0;
	}


	f_flag += fread(&absent, sizeof(absent), 1, fp);
	if(absent != opt->number_of_round)
	{
		fprintf(stderr, "ERROR: different number of rounds!\n");
		exit(1);
	}

	f_flag += fread(&opt->hom_cov, sizeof(opt->hom_cov), 1, fp);
	f_flag += fread(&opt->het_cov, sizeof(opt->het_cov), 1, fp);
	f_flag += fread(&opt->max_n_chain, sizeof(opt->max_n_chain), 1, fp);


	fclose(fp);
	
	if(!load_All_reads(r, gfa_name))
	{
		free(gfa_name);
		return 0;
	}


	memset(r->trio_flag, AMBIGU, r->total_reads*sizeof(uint8_t));
	r->paf = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*r->total_reads);
    r->reverse_paf = (ma_hit_t_alloc*)malloc(sizeof(ma_hit_t_alloc)*r->total_reads);
	for (i = 0; i < (long long)r->total_reads; i++)
    {
        init_ma_hit_t_alloc(&(r->paf[i]));
        init_ma_hit_t_alloc(&(r->reverse_paf[i]));
    }

	fprintf(stderr, "[M::%s] Index has been loaded.\n", __func__);

	free(gfa_name);
	return 1;
}
