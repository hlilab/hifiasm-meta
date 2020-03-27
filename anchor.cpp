#include "htab.h"
#include "ksort.h"
#include "Hash_Table.h"

typedef struct {
	uint64_t srt;
	uint32_t self_off;
	uint32_t other_off;
} anchor1_t;

#define an_key1(a) ((a).srt)
#define an_key2(a) ((a).self_off)
KRADIX_SORT_INIT(ha_an1, anchor1_t, an_key1, 8)
KRADIX_SORT_INIT(ha_an2, anchor1_t, an_key2, 4)

typedef struct {
	int n;
	const ha_idxpos_t *a;
} seed1_t;

struct ha_abuf_s {
	uint64_t n_a, m_a;
	uint32_t old_mz_m;
	ha_mz1_v mz;
	seed1_t *seed;
	anchor1_t *a;
};

ha_abuf_t *ha_abuf_init(void)
{
	return (ha_abuf_t*)calloc(1, sizeof(ha_abuf_t));
}

void ha_abuf_destroy(ha_abuf_t *ab)
{
	free(ab->seed); free(ab->a); free(ab->mz.a); free(ab);
}

void ha_get_new_candidates(ha_abuf_t *ab, int64_t rid, UC_Read *ucr, overlap_region_alloc *overlap_list, Candidates_list *cl, double band_width_threshold, int keep_whole_chain)
{
	extern void *ha_flt_tab;
	extern ha_pt_t *ha_idx;
	uint32_t i;
	uint64_t k, l;

    clear_Candidates_list(cl);
    clear_overlap_region_alloc(overlap_list);
	recover_UC_Read(ucr, &R_INF, rid);
	ab->mz.n = 0, ab->n_a = 0;

	ha_sketch(ucr->seq, ucr->length, asm_opt.mz_win, asm_opt.k_mer_length, 0, !asm_opt.no_HPC, &ab->mz, ha_flt_tab);
	if (ab->mz.m > ab->old_mz_m) {
		ab->old_mz_m = ab->mz.m;
		REALLOC(ab->seed, ab->old_mz_m);
	}
	for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
		ab->seed[i].a = ha_pt_get(ha_idx, ab->mz.a[i].x, &ab->seed[i].n);
		ab->n_a += ab->seed[i].n;
	}
	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		ab->m_a = ab->m_a > 16? ab->m_a + (ab->m_a>>1) : 16;
		REALLOC(ab->a, ab->m_a);
	}
	for (i = 0, k = 0; i < ab->mz.n; ++i) {
		uint32_t j;
		ha_mz1_t *z = &ab->mz.a[i];
		seed1_t *s = &ab->seed[i];
		for (j = 0; j < s->n; ++j) {
			const ha_idxpos_t *y = &s->a[j];
			anchor1_t *an = &ab->a[k++];
			uint8_t rev;
			if (z->rev == y->rev) { // forward strand
				rev = 0;
				an->self_off = z->pos;
				an->other_off = y->pos;
			} else { // reverse strand
				rev = 1;
				an->self_off = ucr->length - 1 - (z->pos + 1 - z->span);
				an->other_off = R_INF.read_length[y->rid] - 1 - (y->pos + 1 - y->span);
			}
			an->srt = (uint64_t)rev<<63 | (uint64_t)y->rid << 32 | (0x80000000ULL + (an->other_off - an->self_off));
		}
	}

	radix_sort_ha_an1(ab->a, ab->a + ab->n_a);
	for (k = 1, l = 0; k <= ab->n_a; ++k) {
		if (k == ab->n_a || ab->a[k].srt != ab->a[l].srt) {
			if (k - l > 1)
				radix_sort_ha_an2(ab->a + l, ab->a + k);
			l = k;
		}
	}

	calculate_overlap_region_by_chaining(cl, overlap_list, rid, ucr->length, &R_INF, band_width_threshold, keep_whole_chain);
}
