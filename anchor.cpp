#include <stdio.h>
#include "htab.h"
#include "ksort.h"
#include "Hash_Table.h"
#include "khashl.h"

#define HA_KMER_GOOD_RATIO 0.333

typedef struct { // this struct is not strictly necessary; we can use k_mer_pos instead, with modifications
	uint64_t srt;
	uint32_t self_off:31, good:1;
	uint32_t other_off;
} anchor1_t;

#define an_key1(a) ((a).srt)
#define an_key2(a) ((a).self_off)
KRADIX_SORT_INIT(ha_an1, anchor1_t, an_key1, 8)
KRADIX_SORT_INIT(ha_an2, anchor1_t, an_key2, 4)

#define oreg_xs_lt(a, b) (((uint64_t)(a).x_pos_s<<32|(a).x_pos_e) < ((uint64_t)(b).x_pos_s<<32|(b).x_pos_e))
KSORT_INIT(or_xs, overlap_region, oreg_xs_lt)

#define oreg_ss_lt(a, b) ((a).shared_seed > (b).shared_seed) // in the decending order
KSORT_INIT(or_ss, overlap_region, oreg_ss_lt)

// hamt mitigation of lost containment 
// #define hamt_ov_eq(a, b) ((a) == (b))
// #define hamt_ov_hash(a) ((a))
// KHASHL_SET_INIT(static klib_unused, hamt_ov_t, hamt_ov, uint64_t, hamt_ov_hash, hamt_ov_eq)


typedef struct {
	int n, good;
	const ha_idxpos_t *a;
} seed1_t;

struct ha_abuf_s {
	uint64_t n_a, m_a;///number of anchors (seed positions)
	uint32_t old_mz_m;///number of seeds
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

uint64_t ha_abuf_mem(const ha_abuf_t *ab)
{
	return ab->m_a * sizeof(anchor1_t) + ab->mz.m * (sizeof(ha_mz1_t) + sizeof(seed1_t)) + sizeof(ha_abuf_t);
}

int ha_ov_type(const overlap_region *r, uint32_t len)
{
	if (r->x_pos_s == 0 && r->x_pos_e == len - 1) return 2; // contained in a longer read
	else if (r->x_pos_s > 0 && r->x_pos_e < len - 1) return 3; // containing a shorter read
	else return r->x_pos_s == 0? 0 : 1;
}

void hamt_count_new_candidates(int64_t rid, UC_Read *ucr, All_reads *rs, int sort_mode){
	// work on one read
	// (immitate ha_get_new_candidates, but only count how many candidates are within consideration.)
	// (used for read selection)
	// get & sort seeds, guess how many alignments will be needed

	extern void *ha_flt_tab;
	extern ha_pt_t *ha_idx;
	uint32_t i/*, rlen*/;
	uint64_t k, l;
	double low_occ = 5;
	double high_occ = 500;

	// containers
	ha_mz1_v mz = {0, 0, 0};  // an array of minimizers
	mz.a = (ha_mz1_t*)calloc(16, sizeof(ha_mz1_t));
	mz.m = 16;
	seed1_t *seed; // slots corresponding to the anchors
	anchor1_t *a;

	// get the list of anchors
	ha_sketch(ucr->seq, ucr->length, asm_opt.mz_win, asm_opt.k_mer_length, 0, !(asm_opt.flag & HA_F_NO_HPC), &mz, ha_flt_tab);
	seed = (seed1_t*)malloc(sizeof(seed1_t) * mz.m);
	uint64_t n_a = 0;
	for (i = 0, n_a = 0; i < mz.n; ++i) {
		int n;
		seed[i].a = ha_pt_get(ha_idx, mz.a[i].x, &n);  // start idx of the minimizer in the linear buffer
		seed[i].n = n;  // count of the minimizer
		seed[i].good = (n > low_occ && n < high_occ);
		n_a += n;
	}
	a = (anchor1_t*)malloc(n_a * sizeof(anchor1_t));
	for (i = 0, k = 0; i < mz.n; ++i) {
		int j;
		ha_mz1_t *z = &mz.a[i];
		seed1_t *s = &seed[i];
		for (j = 0; j < s->n; ++j) {  // for each occurence of this minimizer
			const ha_idxpos_t *y = &s->a[j];  // its appearance
			anchor1_t *an = &a[k++];
			uint8_t rev = z->rev == y->rev? 0 : 1;
			an->other_off = y->pos;
			an->self_off = rev? ucr->length - 1 - (z->pos + 1 - z->span) : z->pos;
			an->good = s->good;
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->other_off;
		}
	}

	uint64_t nb_candidates = 0;
	radix_sort_ha_an1(a, a+n_a);  // sort by srt (targetID-strand-posOnTargetRead)
	for (k = 1, l = 0; k <= n_a; ++k) {
		// if (k == n_a || a[k].srt != a[l].srt) {
		if (k==n_a || ((a[k].srt>>33) != (a[l].srt>>33)) ){
			if (k - l > 100)
				nb_candidates++;  // we have 1<<28-1 reads at most, this won't overflow
			l = k;
		}
	}
	// sort anchors + count
	if (sort_mode==0){
		// fprintf(stdout, "debug\t%d\t%d\n", (int)rid, (int)nb_candidates);
		rs->nb_target_reads[rid] = nb_candidates<<32 | ((uint64_t)rid);
	}else if (sort_mode==1){
		nb_candidates = nb_candidates<65535? nb_candidates : 65535;
		rs->nb_target_reads[rid] = nb_candidates<<48 | ((uint64_t)(65535 - rs->lowq[rid]))<<32 | ((uint64_t) rid);  // lowq sorted the reversed order
	}else{
		fprintf(stderr, "[E::%s] unexpected sort_mode, aborting.\n", __func__);
		exit(1);
	}

	// clean up
	free(a);
	free(mz.a);
	free(seed);

}


void ha_get_new_candidates(ha_abuf_t *ab, int64_t rid, UC_Read *ucr, overlap_region_alloc *overlap_list, Candidates_list *cl, double bw_thres, int max_n_chain, int keep_whole_chain, kvec_t_u8_warp* k_flag,
kvec_t_u64_warp* chain_idx, void *ha_flt_tab, ha_pt_t *ha_idx, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct,
int thread_id)
{
	uint32_t i, rlen;
	uint64_t k, l;
	double low_occ, high_occ;
	// hamt_ov_t **hs = (hamt_ov_t**)R_INF.hamt_existed_ov;
	// hamt_ov_t *hamt_existed_ov=NULL; 
	// if (asm_opt.is_final_round)
	// 	hamt_existed_ov= hs[thread_id];  // need per thread hashtable, otherwise expanding the hashtable will be race

	if (asm_opt.is_use_exp_graph_cleaning){
		// might have bugs / TODO
		low_occ = 5;  // hamt, this is and have been arbitrary from the beginning (Dec 25,2020). What's the better way or value?
		high_occ = 100;
	}else{
		low_occ = asm_opt.hom_cov * HA_KMER_GOOD_RATIO;
		high_occ = asm_opt.hom_cov * (2.0 - HA_KMER_GOOD_RATIO);
	}	

	// prepare
    clear_Candidates_list(cl);
    clear_overlap_region_alloc(overlap_list);
	recover_UC_Read(ucr, &R_INF, rid);
	ab->mz.n = 0, ab->n_a = 0;
	rlen = Get_READ_LENGTH(R_INF, rid);

	// get the list of anchors
	ha_sketch_query(ucr->seq, ucr->length, asm_opt.mz_win, asm_opt.k_mer_length, 0, !(asm_opt.flag & HA_F_NO_HPC), &ab->mz, ha_flt_tab, k_flag, dbg_ct);
	// minimizer of queried read
	if (ab->mz.m > ab->old_mz_m) {
		ab->old_mz_m = ab->mz.m;
		REALLOC(ab->seed, ab->old_mz_m);
	}
	for (i = 0, ab->n_a = 0; i < ab->mz.n; ++i) {
		int n;
		ab->seed[i].a = ha_pt_get(ha_idx, ab->mz.a[i].x, &n);  // start idx of the minimizer in the linear buffer
		ab->seed[i].n = n;  // count of the minimizer
		ab->seed[i].good = (n > low_occ && n < high_occ);
		ab->n_a += n;
	}
	if (ab->n_a > ab->m_a) {
		ab->m_a = ab->n_a;
		kroundup64(ab->m_a);
		REALLOC(ab->a, ab->m_a);
	}
	for (i = 0, k = 0; i < ab->mz.n; ++i) {
		int j;
		///z is one of the minimizer
		ha_mz1_t *z = &ab->mz.a[i];
		seed1_t *s = &ab->seed[i];
		for (j = 0; j < s->n; ++j) {  // for each occurence of this minimizer
			const ha_idxpos_t *y = &s->a[j];  // its appearance
			anchor1_t *an = &ab->a[k++];
			uint8_t rev = z->rev == y->rev? 0 : 1;
			an->other_off = y->pos;
			an->self_off = rev? ucr->length - 1 - (z->pos + 1 - z->span) : z->pos;
			an->good = s->good;
			an->srt = (uint64_t)y->rid<<33 | (uint64_t)rev<<32 | an->other_off;
		}
	}

	// sort anchors
	radix_sort_ha_an1(ab->a, ab->a + ab->n_a);  // sort by srt (targetID-strand-posOnTargetRead)
	for (k = 1, l = 0; k <= ab->n_a; ++k) {
		if (k == ab->n_a || ab->a[k].srt != ab->a[l].srt) {
			if (k - l > 1)
				radix_sort_ha_an2(ab->a + l, ab->a + k);  // sort by posOnSelfRead
			l = k;
		}
	}

	// (store the info of "an overlap has existed" before anything)
	// if (asm_opt.is_final_round){
	// 	int absent = 0;
	// 	uint64_t key=0, key_new, xid=rid, yid, yid_old=0;
	// 	for (int i_ov=0; i_ov<ab->n_a; i_ov++){
	// 		overlap_region *handle = &overlap_list->list[i_ov];
	// 		yid = ab->a[i_ov].srt>>33;
	// 		key = xid<<(32+3) | (yid<<3);  // read id is actually at most 28bits. Here using last 6 bits for counting
	// 		yid_old = yid;
	// 		hamt_ov_put(hamt_existed_ov, key, &absent);
	// 	}
	// }


	// copy over to cl
	// (cl for Candidate List)
	if (ab->m_a >= (uint64_t)cl->size) {
		cl->size = ab->m_a;
		REALLOC(cl->list, cl->size);
	}
	for (k = 0; k < ab->n_a; ++k) {
		k_mer_hit *p = &cl->list[k];
		p->readID = ab->a[k].srt >> 33;
		p->strand = ab->a[k].srt >> 32 & 1;
		p->offset = ab->a[k].other_off;
		p->self_offset = ab->a[k].self_off;
		p->good = ab->a[k].good;
	}
	cl->length = ab->n_a;

	calculate_overlap_region_by_chaining(cl, overlap_list, chain_idx, rid, ucr->length, &R_INF, bw_thres, keep_whole_chain, f_cigar);

	#if 0
	if (overlap_list->length > 0) {
		fprintf(stderr, "B\t%ld\t%ld\t%d\n", (long)rid, (long)overlap_list->length, rlen);
		for (int i = 0; i < (int)overlap_list->length; ++i) {
			overlap_region *r = &overlap_list->list[i];
			fprintf(stderr, "C\t%d\t%d\t%d\t%c\t%d\t%ld\t%d\t%d\t%c\t%d\t%d\n", (int)r->x_id, (int)r->x_pos_s, (int)r->x_pos_e, "+-"[r->x_pos_strand],
					(int)r->y_id, (long)Get_READ_LENGTH(R_INF, r->y_id), (int)r->y_pos_s, (int)r->y_pos_e, "+-"[r->y_pos_strand], (int)r->shared_seed, ha_ov_type(r, rlen));
		}
	}
	#endif
	if (overlap_list->length > 0) {
		if (asm_opt.is_dump_relevant_reads && asm_opt.is_final_round){
			char *str_buf = (char*)malloc(1<<10);
			char *str_line_buf = (char*)malloc(1<<9);  // BUG: bug if read name is super long.
			int str_buf_l = 1<<10, str_buf_cnt = 0, tmp;
			sprintf(str_buf, "R\t%.*s\t%d\t%d\n", 
							(int)Get_NAME_LENGTH(R_INF, rid), Get_NAME(R_INF, rid),
							rlen,
							(int)overlap_list->length);
			str_buf_cnt = strlen(str_buf);
			for (int i = 0; i < (int)overlap_list->length; ++i) {
				overlap_region *r = &overlap_list->list[i];
				sprintf(str_line_buf, "C\t%.*s\t%ld\t%d\t%d\t%c\t%.*s\t%ld\t%d\t%d\t%c\t%d\t%d\n", 
								(int)Get_NAME_LENGTH(R_INF, r->x_id), Get_NAME(R_INF, r->x_id), // (int)r->x_id, 
								(long)Get_READ_LENGTH(R_INF, r->x_id), (int)r->x_pos_s, (int)r->x_pos_e, "+-"[r->x_pos_strand],
								(int)Get_NAME_LENGTH(R_INF, r->y_id), Get_NAME(R_INF, r->y_id), // (int)r->y_id, 
								(long)Get_READ_LENGTH(R_INF, r->y_id), (int)r->y_pos_s, (int)r->y_pos_e, "+-"[r->y_pos_strand], 
								(int)r->shared_seed, ha_ov_type(r, rlen));
				tmp = strlen(str_line_buf);
				if (str_buf_cnt+tmp>=str_buf_l){
					str_buf_l = str_buf_l + (str_buf_l<<1);
					str_buf = (char*)realloc(str_buf, str_buf_l);
				}
				sprintf(str_buf+str_buf_cnt, "%s", str_line_buf);
				str_buf_cnt += tmp;
			}
			fprintf(asm_opt.fp_relevant_reads, "%s", str_buf);
			free(str_buf);
			free(str_line_buf);
		}
	}

	if ((int)overlap_list->length > max_n_chain) {
		int32_t w, n[4], s[4];
		n[0] = n[1] = n[2] = n[3] = 0, s[0] = s[1] = s[2] = s[3] = 0;
		ks_introsort_or_ss(overlap_list->length, overlap_list->list);
		for (i = 0; i < (uint32_t)overlap_list->length; ++i) {
			const overlap_region *r = &overlap_list->list[i];
			w = ha_ov_type(r, rlen);
			++n[w];
			if ((int)n[w] == max_n_chain) s[w] = r->shared_seed;
		}
		if (s[0] > 0 || s[1] > 0 || s[2] > 0 || s[3] > 0) {
			for (i = 0, k = 0; i < (uint32_t)overlap_list->length; ++i) {
				overlap_region *r = &overlap_list->list[i];
				w = ha_ov_type(r, rlen);
				if (r->shared_seed >= s[w]) {
					if ((uint32_t)k != i) {
						overlap_region t;
						t = overlap_list->list[k];
						overlap_list->list[k] = overlap_list->list[i];
						overlap_list->list[i] = t;
					}
					++k;
				}
			}
			overlap_list->length = k;
		}
	}

	///ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}


void lable_matched_ovlp(overlap_region_alloc* overlap_list, ma_hit_t_alloc* paf)
{
	uint64_t j = 0, inner_j = 0;
	while (j < overlap_list->length && inner_j < paf->length)
    {
        if(overlap_list->list[j].y_id < paf->buffer[inner_j].tn)
        {
            j++;
        }
        else if(overlap_list->list[j].y_id > paf->buffer[inner_j].tn)
        {
            inner_j++;
        }
        else
        {
            if(overlap_list->list[j].y_pos_strand == paf->buffer[inner_j].rev)
            {
				overlap_list->list[j].is_match = 1;
            }
            j++;
            inner_j++;
        }
    }
}


void ha_get_candidates_interface(ha_abuf_t *ab, int64_t rid, UC_Read *ucr, overlap_region_alloc *overlap_list, overlap_region_alloc *overlap_list_hp, Candidates_list *cl, double bw_thres, 
int max_n_chain, int keep_whole_chain, kvec_t_u8_warp* k_flag, kvec_t_u64_warp* chain_idx, ma_hit_t_alloc* paf, ma_hit_t_alloc* rev_paf, overlap_region* f_cigar, kvec_t_u64_warp* dbg_ct,
int thread_id)
{
	extern void *ha_flt_tab;
	extern ha_pt_t *ha_idx;
	extern void *ha_flt_tab_hp;
	extern ha_pt_t *ha_idx_hp;

	ha_get_new_candidates(ab, rid, ucr, overlap_list, cl, bw_thres, max_n_chain, keep_whole_chain, k_flag, chain_idx, ha_flt_tab, ha_idx, f_cigar, dbg_ct,
						  thread_id);

	if(ha_idx_hp)
	{
		// hamt note: is this block never envoked? ha_idx_hp is only calculated by rescue_hp_reads,
		//            and rescue_hp_reads isn't called anywhere.
		uint32_t i, k, y_id, overlapLen, max_i;
		int shared_seed;
		overlap_region t;
		overlap_region_sort_y_id(overlap_list->list, overlap_list->length);
		ma_hit_sort_tn(paf->buffer, paf->length);
		ma_hit_sort_tn(rev_paf->buffer, rev_paf->length);
		lable_matched_ovlp(overlap_list, paf);
		lable_matched_ovlp(overlap_list, rev_paf);

		for (i = 0, k = 0; i < overlap_list->length; ++i) 
		{
			if(overlap_list->list[i].is_match == 1)
			{
				if(k != i)
				{
					t = overlap_list->list[k];
					overlap_list->list[k] = overlap_list->list[i];
					overlap_list->list[i] = t;
					overlap_list->list[k].is_match = 0;
				}
				k++;
			}
		}
		overlap_list->length = k;


		ha_get_new_candidates(ab, rid, ucr, overlap_list_hp, cl, bw_thres, max_n_chain, keep_whole_chain, k_flag, chain_idx, ha_flt_tab_hp, ha_idx_hp, f_cigar, dbg_ct,
								thread_id);	
		
		if(overlap_list->length + overlap_list_hp->length > overlap_list->size)
		{
			overlap_list->list = (overlap_region*)realloc(overlap_list->list, 
				sizeof(overlap_region)*(overlap_list->length + overlap_list_hp->length));
			memset(overlap_list->list + overlap_list->size, 0, sizeof(overlap_region)*
				(overlap_list->length + overlap_list_hp->length - overlap_list->size));
			overlap_list->size = overlap_list->length + overlap_list_hp->length;
		}
		
		for (i = 0, k = overlap_list->length; i < overlap_list_hp->length; i++, k++)
		{
			t = overlap_list->list[k];
			overlap_list->list[k] = overlap_list_hp->list[i];
			overlap_list_hp->list[i] = t;
		}
		overlap_list->length = k;
		
		overlap_region_sort_y_id(overlap_list->list, overlap_list->length);

		i = k = 0;
		while (i < overlap_list->length)
		{
			y_id = overlap_list->list[i].y_id;
			shared_seed = overlap_list->list[i].shared_seed;
			overlapLen = overlap_list->list[i].overlapLen;
			max_i = i;
			i++;
			while (i < overlap_list->length && overlap_list->list[i].y_id == y_id)
			{
				if((overlap_list->list[i].shared_seed > shared_seed) || 
				  ((overlap_list->list[i].shared_seed == shared_seed) && (overlap_list->list[i].overlapLen <= overlapLen)))
				{
					y_id = overlap_list->list[i].y_id;
					shared_seed = overlap_list->list[i].shared_seed;
					overlapLen = overlap_list->list[i].overlapLen;
					max_i = i;
				}
				i++;
			}

			if(k != max_i)
			{
				t = overlap_list->list[k];
				overlap_list->list[k] = overlap_list->list[max_i];
				overlap_list->list[max_i] = t;
			}
			k++;
		}

		overlap_list->length = k;
	}
	
	ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}

void ha_sort_list_by_anchor(overlap_region_alloc *overlap_list)
{
	ks_introsort_or_xs(overlap_list->length, overlap_list->list);
}
